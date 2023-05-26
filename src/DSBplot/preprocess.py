import os

import argparse
import glob

import DSBplot.utils.file_names as file_names
import DSBplot.utils.common_utils as common_utils
import DSBplot.utils.file_utils as file_utils
import DSBplot.utils.log_utils as log_utils
import DSBplot.utils.constants as constants

import DSBplot.get_nhej_data.filter_nhej as filter_nhej
import DSBplot.get_nhej_data.combine_repeat as combine_repeat
import DSBplot.get_window_data.get_window as get_window
import DSBplot.get_window_data.get_freq as get_freq
import DSBplot.get_window_data.get_freq_comparison as get_freq_comparison
import DSBplot.get_graph_data.get_graph_data as get_graph_data
import DSBplot.get_histogram_data.get_histogram_data as get_histogram_data

STAGES_1 = ['0_align', '1_filter', '2_combine', '3_window']
STAGES_2 = ['3_comparison', '4_graph', '5_histogram']

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
      ' Required only for stages 0_align.'
    ),
  )
  parser.add_argument(
    '--output',
    type = common_utils.check_dir_output,
    help = 'Output directory. Required for all stages.',
    required = True,
  )
  parser.add_argument(
    '--stages',
    type = str,
    nargs = '+',
    choices = STAGES_1 + STAGES_2,
    default = STAGES_1 + [x for x in STAGES_2 if x != '3_comparison'],
    help = (
      'Stages to run.' +
      ' The stages must be run in the correct order indicated by their prefix numbers' +
      ' (3_comparison must come after 3_window).' +
      ' Briefly, the stages are:' +
      ' 0_align: Align reads to reference sequence using Bowtie 2.' +
      ' 1_filter: Filter reads that represent NHEJ repair (all in/dels touch the' +
      ' DSB site and are contiguous).' +
      ' 2_combine: Combine the repeats of from the same experiment into a single file' +
      ' with the mean frequencies of the repeats.' +
      ' 3_window: Extract the nucleotide sequences and alignment from around the DSB site.' +
      ' 3_comparison: Combine the data from two different experiments in preparation for' +
      ' making comparison graphs (this stage cannot be run from "preprocess.py", only "comparison.py").' +
      ' 4_graph: Preprocess the data further for graph construction and plotting.' +
      ' 5_hist: Preprocess the data further for histogram construction and plotting.' +
      ' For a more detailed explanation of the preprocessing pipeline, see the README.'
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
    '--freq_min',
    type = float,
    default = 1e-5,
    help = (
      f'Minimum frequency for output in' +
      f' windows_{constants.FREQ_FILTER_MEAN}.tsv.' +
      f' Sequences with frequencies <= this are discarded.' +
      f' Required only for stage 3_window' +
      f' (but can be omitted because of default).'
    ),
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
    '--quiet',
    action = 'store_true',
    help = 'If present, do no output verbose log message.',
  )
  return vars(parser.parse_args())

def get_sam_file(output, i):
  return os.path.join(file_names.sam_dir(output), f'{i}.sam')

def get_filter_nhej_file(output, i):
  return os.path.join(file_names.filter_nhej_dir(output), f'{i}.tsv')

def do_0_align(
  input,
  output,
  ref_seq_file,
):
  if input is None:
    raise Exception('INPUT must be provided for stage 0_align.')
  if ref_seq_file is None:
    raise Exception('REF_SEQ_FILE must be provided for stage 0_align.')
  bowtie2_build_file = file_names.bowtie2_build(output)
  file_utils.make_parent_dir(bowtie2_build_file)
  log_utils.log_input('Bowtie2 build file: ' + bowtie2_build_file)
  os.system(f'bowtie2-build-s {ref_seq_file} {bowtie2_build_file} --quiet')

  for i, input_1 in enumerate(input, 1):
    sam_file = get_sam_file(output, i)
    file_utils.make_parent_dir(sam_file)
    input_ext = os.path.splitext(input_1)[1]
    if input_ext in ['.fastq', '.fq']:
      flags = '-q'
    elif input_ext in ['.fasta', '.fa', '.fna']:
      flags = '-f'
    else:
      flags = '-r'
    bowtie2_command = f'bowtie2-align-s {flags} -x {bowtie2_build_file} {input_1} -S {sam_file} --quiet'
    log_utils.log('Bowtie 2 command: ' + bowtie2_command)
    os.system(f'bowtie2-align-s {flags} -x {bowtie2_build_file} {input_1} -S {sam_file} --quiet')
    log_utils.blank_line()

def do_1_filter_nhej(
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

  min_length = max(dsb_pos + 1, min_length)
  input_num = len(glob.glob(os.path.join(file_names.sam_dir(output), '*.sam')))

  for i in range(1, input_num + 1):
    debug_file = os.path.join(file_names.filter_nhej_dir(output), f'debug_{i}.txt')
    filter_nhej.main(
      ref_seq_file = ref_seq_file,
      sam_file = get_sam_file(output, i),
      output = get_filter_nhej_file(output, i),
      dsb_pos = dsb_pos,
      min_length = min_length,
      quiet = quiet,
      debug_file = debug_file,
    )

def do_2_combine_repeat(
  output,
  quiet,
):
  input = glob.glob(os.path.join(file_names.filter_nhej_dir(output), '*.tsv'))
  input = sorted(input, key = lambda x: int(os.path.basename(x).split('.')[0]))
  combine_repeat.main(
    input = input,
    column_names = [str(i) for i in range(1, len(input) + 1)],
    output = file_names.combine_repeat_file(output),
    quiet = quiet,
  )

def do_3_window(
  output,
  ref_seq_file,
  dsb_pos,
  window_size,
  anchor_size,
  anchor_mismatches,
  total_reads,
  freq_min,
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
  if freq_min is None:
    raise Exception('FREQ_MIN must be provided for stage 3_window.')
  if label is None:
    label = os.path.basename(output)

  combine_repeat_file = file_names.combine_repeat_file(output)
  window_dir = file_names.window_dir(output)
  for subst_type in constants.SUBST_TYPES:
    get_window.main(
      input = combine_repeat_file,
      output = window_dir,
      ref_seq_file = ref_seq_file,
      dsb_pos = dsb_pos,
      window_size = window_size,
      anchor_size = anchor_size,
      anchor_mismatches = anchor_mismatches,
      subst_type = subst_type,
      label = label,
    )
    get_freq.main(
      input = window_dir,
      output = window_dir,
      subst_type = subst_type,
      total_reads = total_reads,
      freq_min = freq_min,
    )

def do_3_comparison(
  input_1,
  input_2,
  output,
):
  for subst_type in constants.SUBST_TYPES:
    get_freq_comparison.main(
      input = [file_names.window_dir(input_1), file_names.window_dir(input_2)],
      output = file_names.window_dir(output),
      subst_type = subst_type,
    )

def do_4_graph(output):
  for subst_type in constants.SUBST_TYPES:
    get_graph_data.main(
      input = file_names.window_dir(output),
      output = file_names.graph_dir(output),
      subst_type = subst_type,
    )

def do_5_histogram(output):
  for subst_type in constants.SUBST_TYPES:
    get_histogram_data.main(
      input = file_names.graph_dir(output),
      output = file_names.histogram_dir(output),
      subst_type = subst_type,
    )

# The alignment, filtering, and combining stages.
# Everything up to the point of combine samples for comparison.
def do_stages_1(
  input,
  output,
  ref_seq_file,
  dsb_pos,
  min_length,
  window_size,
  anchor_size,
  anchor_mismatches,
  total_reads,
  freq_min,
  label,
  quiet,
  stages = STAGES_1,
):
  if '0_align' in stages:
    do_0_align(
      input = input,
      output = output,
      ref_seq_file = ref_seq_file,
    )

  if '1_filter' in stages:
    do_1_filter_nhej(
      output = output,
      ref_seq_file = ref_seq_file,
      dsb_pos = dsb_pos,
      min_length = min_length,
      quiet = quiet,
    )

  if '2_combine' in stages:
    do_2_combine_repeat(
      output = output,
      quiet = quiet,
    )

  if '3_window' in stages:
    do_3_window(
      output = output,
      ref_seq_file = ref_seq_file,
      dsb_pos = dsb_pos,
      window_size = window_size,
      anchor_size = anchor_size,
      anchor_mismatches = anchor_mismatches,
      total_reads = total_reads,
      freq_min = freq_min,
      label = label,
    )

# The stages from (inclusive) the point of merging samples for comparison.
def do_stages_2(
  output,
  input_comparison = None,
  stages = [x for x in STAGES_2 if x != '3_comparison'],
):
  if ('3_comparison' in stages) and (input_comparison is not None):
    do_3_comparison(
      input_1 = input_comparison[0],
      input_2 = input_comparison[1],
      output = output,
    )

  if '4_graph' in stages:
    do_4_graph(output)

  if '5_histogram' in stages:
    do_5_histogram(output)

def main(
  input,
  output,
  ref_seq_file,
  dsb_pos,
  min_length,
  window_size,
  anchor_size,
  anchor_mismatches,
  total_reads,
  freq_min,
  label,
  quiet,
  stages = STAGES_1 + [x for x in STAGES_2 if x != '3_comparison'],
):
  do_stages_1(
    input = input,
    output = output,
    ref_seq_file = ref_seq_file,
    dsb_pos = dsb_pos,
    min_length = min_length,
    window_size = window_size,
    anchor_size = anchor_size,
    anchor_mismatches = anchor_mismatches,
    total_reads = total_reads,
    freq_min = freq_min,
    label = label,
    quiet = quiet,
    stages = stages,
  )
  do_stages_2(output = output, stages = stages)

# This allows the "DSBplot-preprocess" command to be run from the command line.
def entry_point():
  main(**parse_args())

if __name__ == '__main__':
  entry_point()
