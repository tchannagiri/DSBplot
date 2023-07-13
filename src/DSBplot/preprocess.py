import os

import shutil
import argparse
import glob

import DSBplot.utils.file_names as file_names
import DSBplot.utils.common_utils as common_utils
import DSBplot.utils.file_utils as file_utils
import DSBplot.utils.log_utils as log_utils
import DSBplot.utils.constants as constants

import DSBplot.get_nhej_data.filter_nhej as filter_nhej
import DSBplot.get_window_data.get_window as get_window
import DSBplot.get_histogram_data.get_variation as get_variation



STAGES = ['0_align', '1_filter', '2_window', '3_variation']

# FIXME RE WRITE THE DOCU MEN TATION FOR THE NEW STAGES
def parse_args():
  parser = argparse.ArgumentParser(
    description = (
      'Perform alignment and preprocessing for raw FASTQ data.' +
      ' This is script is broken in separate stages so that each stage' +
      ' can be run separately. However, the stages must be run in the correct order indicated' +
      ' by their prefix numbers. If running the stages separately, the value of OUTPUT' +
      ' must the same value on each separate invocation. Two experiments should not' +
      ' be given the same OUTPUT directory. See parameter --stages and the README for more information' +
      ' on the individual stages.'
    ),
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
  )
  parser.add_argument(
    '--input',
    nargs = '+',
    type = common_utils.check_file,
    help = (
      'Input files of raw reads.' +
      ' The following file entensions are allowed: ' +
      ' FASTQ: ".fastq", ".fq."; FASTA: "fasta", ".fa", "fna";' +
      ' text: all others. FASTQ files are processed the ' +
      ' Bowtie 2 flag "-q", FASTA files are processed with the Bowtie 2 flag "-f",' +
      ' and text files are processed with the Bowtie 2 flag "-r".' +
      ' Please see the Bowtie 2 manual at http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml.' +
      ' Each file is considered a repeat of the same experiment.' +
      ' Required for stage 0_align. ' +
      ' If 0_align is omitted and a custom alignnment approach is used, ' +
      ' then the aligned SAM files must be placed in OUTPUT directory.'
    ),
  )
  parser.add_argument(
    '--names',
    nargs = '+',
    type = str,
    help = (
      'Identifiers to use for the input libraries.' +
      ' Must be given in the same order as the input files.' +
      ' If omitted, the names will be the input file names without the extension.' +
      ' Required for stages 0_align and 1_filter (but can be omitted because of the default).'
    ),
  )
  parser.add_argument(
    '--output',
    type = common_utils.check_dir_output,
    help = (
      'Working output directory. Required for all stages.' +
      ' If running stages separately, the same directory must be given for all stages.'
    ),
    required = True,
  )
  parser.add_argument(
    '--stages',
    type = str,
    nargs = '+',
    choices = STAGES,
    default = STAGES,
    help = (
      'Stages to run.' +
      ' The stages must be run in the correct order indicated by their prefix numbers' +
      ' Briefly, the stages are:' +
      ' 0_align: Align reads to reference sequence using Bowtie 2 ' +
      ' (may be omitted if a custom alignment approach is used and the resulting SAM files are put in the OUTPUT directory).' +
      ' 1_filter: Filter reads that represent NHEJ repair (all in/dels touch the DSB site and are consecutive).' +
      ' 2_window: Extract the nucleotide sequences and alignment from around the DSB site.' +
      ' 3_variation: Split alignment windows into individual variations (used for plotting variation-position histograms).'
      # ' 3_comparison: Combine the data from two different experiments in preparation for' +
      # ' making comparison graphs (this stage cannot be run from "preprocess.py", only "comparison.py").' +
      # ' 4_graph: Preprocess the data further for graph construction and plotting.' +
      # ' 5_hist: Preprocess the data further for histogram construction and plotting.' +
      # ' For a more detailed explanation of the preprocessing pipeline, see the README.'
    )
  )
  parser.add_argument(
    '--ref_seq_file',
    type = common_utils.check_file,
    help = (
      'Reference sequence file.' +
      ' Should contain a single nucleotide sequence.' +
      ' Must be the same sequence as used for alignment.' +
      f' Must be in FASTA format ({", ".join(constants.FASTA_EXT)}) or' +
      ' text format (all other extensions).' +
      ' Required for stages 0_align and 1_filter.'
    ),
  )
  parser.add_argument(
    '--dsb_pos',
    type = int,
    help = (
      'Position on reference sequence immediately left (5'') of DSB site.' +
      ' I.e., the DSB is between position DSB_POS and DSB_POS + 1.' +
      ' Required only for stage 1_filter.'
    ),
  )
  parser.add_argument(
    '--min_length',
    type = int,
    help = (
      'Minimum length of read sequence to be considered.' +
      ' Reads shorter than this are discarded.' +
      ' Forced to be at least DSB_POS + 1.' +
      ' If not given, will be assigned DSB_POS + 1.' +
      ' Required only for stage 1_filter' +
      ' (but can be omitted because of default).'
    ),
  )
  parser.add_argument(
    '--window_size',
    type = int,
    default = 10,
    help = (
      'Size of window around DSB site to extract.' +
      ' The nucleotides at the positions' +
      ' {DSB_POS - WINDOW_SIZE + 1, ..., DSB_POS + WINDOW_SIZE} are extracted.' +
      ' The actual number of nucleotides extracted from each read may vary depending' +
      ' on the number of insertions/deletions/substitutions in the alignment.' +
      ' Required only for stage 3_window' +
      ' (but can be omitted because of default).'
    ),
  )
  parser.add_argument(
    '--anchor_size',
    type = int,
    default = 20,
    help = (
      'Size of the anchor on the left/right of the extracted window' +
      ' to check for mismatches.' +
      ' Required only for stage 3_window' +
      ' (but can be omitted because of default).'
    ),
  )
  parser.add_argument(
    '--anchor_mismatches',
    type = int,
    default = 1,
    help = (
      'Maximum number of mismatches allowed on the left/right anchor sequences.\n'
      'Reads with more than the allowed number of mismatches on the left/right anchor\n'
      'will be discarded. This limit is applied to the left/right anchors separately.' +
      ' Required only for stage 3_window' +
      ' (but can be omitted because of default).'
    ),
  )
  parser.add_argument(
    '--total_reads',
    type = int,
    help = (
      'Total number reads in each experiment.' +
      ' This may be strictly greater than the number of reads in the input FASTQ' +
      ' files if some reads were discarded during preprocessing.' +
      ' The number of arguments must be the same as the number of ' +
      ' INPUTs. If not provided, the total reads remaining after NHEJ filtering are used.' +
      ' Required only for stage 2_combine' +
      ' (but can be omitted because of default).'
    ),
    nargs = '+',
  )
  parser.add_argument(
    '--label',
    type = str,
    help = (
      'Label of the experiment to be used in plot legends.' +
      ' If not provided, the label is the basename of the OUTPUT directory.' +
      ' Required only for stage 3_window' +
      ' (but can be omitted because of default).'
    ),
  )
  parser.add_argument(
    '--bowtie2_args',
    type = str,
    default = '',
    help = (
      'Additional arguments to pass to Bowtie2.' +
      ' Must be a single string.' +
      ' Required only for stage 0_align' +
      ' (but can be omitted because of default).'
    ),
  )
  parser.add_argument(
    '--quiet',
    action = 'store_true',
    help = 'If present, do no output verbose log message.',
  )
  args = vars(parser.parse_args())
  if args['names'] is not None:
    if len(args['names']) != len(args['input']):
      raise Exception('Number of NAMES must be the same as the number of INPUTs.')
  return args

def get_filter_nhej_file(output, i, rejected=False):
  if rejected:
    return os.path.join(file_names.filter_nhej_dir(output), f'{i}_rejected.tsv')
  else:
    return os.path.join(file_names.filter_nhej_dir(output), f'{i}.tsv')

def do_0_align(
  input,
  names,
  output,
  ref_seq_file,
  bowtie2_args,
):
  if input is None:
    raise Exception('INPUT must be provided for stage 0_align.')
  if ref_seq_file is None:
    raise Exception('REF_SEQ_FILE must be provided for stage 0_align.')
  if names is None:
    names = [file_names.get_file_name(x) for x in input]
  bowtie2_build_file = file_names.bowtie2_build(output)
  file_utils.make_parent_dir(bowtie2_build_file)
  log_utils.log_input('Bowtie2 build file: ' + bowtie2_build_file)
  os.system(f'bowtie2-build-s {ref_seq_file} {bowtie2_build_file} --quiet')

  for i in range(len(input)):
    sam_file = file_names.sam_file(output, names[i])
    file_utils.make_parent_dir(sam_file)
    input_ext = os.path.splitext(input[i])[1]
    if input_ext in ['.fastq', '.fq']:
      flags = '-q'
    elif input_ext in ['.fasta', '.fa', '.fna']:
      flags = '-f'
    else:
      flags = '-r'
    bowtie2_command = f'bowtie2-align-s --no-hd {flags} {bowtie2_args} -x {bowtie2_build_file} {input[i]} -S {sam_file} --quiet'
    log_utils.log('Bowtie 2 command: ' + bowtie2_command)
    os.system(bowtie2_command)
    log_utils.blank_line()

def do_1_filter_nhej(
  names,
  total_reads,
  output,
  ref_seq_file,
  dsb_pos,
  min_length,
  quiet,
):
  if ref_seq_file is None:
    raise Exception('REF_SEQ_FILE must be provided for stage 1_filter.')
  if dsb_pos is None:
    raise Exception('DSB_POS must be provided for stage 1_filter.')
  if min_length is None:
    min_length = dsb_pos + 1

  if names is None:
    # Infer names from the SAM files.
    names = sorted([
      file_names.get_file_name(x) for x in glob.glob(os.path.join(output, '*.sam'))
    ])
  input = [file_names.sam_file(output, x) for x in names]

  min_length = max(dsb_pos + 1, min_length)
  filter_nhej.main(
    input = input,
    names = names,
    total_reads = total_reads,
    ref_seq_file = ref_seq_file,
    output = file_names.filter_nhej(output, 'accepted'),
    output_rejected = file_names.filter_nhej(output, 'rejected'),
    dsb_pos = dsb_pos,
    min_length = min_length,
    quiet = quiet,
    debug_file = file_names.filter_nhej(output, 'debug'),
  )

def do_2_window(
  output,
  ref_seq_file,
  dsb_pos,
  window_size,
  anchor_size,
  anchor_mismatches,
  label,
):
  if ref_seq_file is None:
    raise Exception('REF_SEQ_FILE must be provided for stage 3_window.')
  if dsb_pos is None:
    raise Exception('DSB_POS must be provided for stage 3_window.')
  if window_size is None:
    raise Exception('WINDOW_SIZE must be provided for stage 3_window.')
  if anchor_size is None:
    raise Exception('ANCHOR_SIZE must be provided for stage 3_window.')
  if anchor_mismatches is None:
    raise Exception('ANCHOR_MISMATCHES must be provided for stage 3_window.')
  if label is None:
    label = os.path.basename(output)

  for subst_type in constants.SUBST_TYPES:
    get_window.main(
      input = file_names.filter_nhej(output, 'accepted'),
      ref_seq_file = ref_seq_file,
      output = file_names.window(output, subst_type),
      output_info = file_names.data_info(output),
      dsb_pos = dsb_pos,
      window_size = window_size,
      anchor_size = anchor_size,
      anchor_mismatches = anchor_mismatches,
      subst_type = subst_type,
      label = label,
    )

def do_3_variation(output):
  for subst_type in constants.SUBST_TYPES:
    get_variation.main(
      input = file_names.window(output, subst_type),
      output = file_names.variation(output, subst_type),
    )


def do_stages(
  input,
  names,
  output,
  ref_seq_file,
  dsb_pos,
  min_length,
  window_size,
  anchor_size,
  anchor_mismatches,
  total_reads,
  label,
  bowtie2_args,
  quiet,
  stages = STAGES,
):
  # Copy reference sequence file to output directory.
  if ref_seq_file is not None:
    shutil.copy(ref_seq_file, file_names.ref_seq_file(output))

  if '0_align' in stages:
    do_0_align(
      input = input,
      names = names,
      output = output,
      ref_seq_file = file_names.ref_seq_file(output),
      bowtie2_args = bowtie2_args,
    )

  if '1_filter' in stages:
    do_1_filter_nhej(
      names = names,
      total_reads = total_reads,
      output = output,
      ref_seq_file = ref_seq_file,
      dsb_pos = dsb_pos,
      min_length = min_length,
      quiet = quiet,
    )

  if '2_window' in stages:
    do_2_window(
      output = output,
      ref_seq_file = file_names.ref_seq_file(output),
      dsb_pos = dsb_pos,
      window_size = window_size,
      anchor_size = anchor_size,
      anchor_mismatches = anchor_mismatches,
      label = label,
    )

  if '3_variation' in stages:
    do_3_variation(output = output)

def main(
  input,
  names,
  output,
  ref_seq_file,
  dsb_pos,
  min_length,
  window_size,
  anchor_size,
  anchor_mismatches,
  total_reads,
  label,
  bowtie2_args,
  quiet,
  stages = STAGES,
):
  do_stages(
    input = input,
    names  = names,
    output = output,
    ref_seq_file = ref_seq_file,
    dsb_pos = dsb_pos,
    min_length = min_length,
    window_size = window_size,
    anchor_size = anchor_size,
    anchor_mismatches = anchor_mismatches,
    total_reads = total_reads,
    label = label,
    bowtie2_args = bowtie2_args,
    quiet = quiet,
    stages = stages,
  )

# This allows the "DSBplot-preprocess" command to be run from the command line.
def entry_point():
  main(**parse_args())

if __name__ == '__main__':
  entry_point()
