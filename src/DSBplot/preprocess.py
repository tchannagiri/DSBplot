import os

import shutil
import argparse
import glob

import DSBplot.utils.file_names as file_names
import DSBplot.utils.common_utils as common_utils
import DSBplot.utils.file_utils as file_utils
import DSBplot.utils.log_utils as log_utils
import DSBplot.utils.constants as constants

import DSBplot.preprocessing.filter as filter
import DSBplot.preprocessing.get_window as get_window
import DSBplot.preprocessing.get_variation as get_variation

STAGES = ['0_align', '1_filter', '2_window', '3_variation']

def parse_args():
  parser = argparse.ArgumentParser(
    description = (
      'Perform alignment and preprocessing for raw FASTQ data.' +
      ' This is script is broken into separate stages so that each stage' +
      ' can be run separately. However, the stages must be run in the correct order indicated' +
      ' by their prefix numbers. If running the stages separately, the value of OUTPUT' +
      ' must the same value on each separate invocation. Two experiments should not' +
      ' be given the same OUTPUT directory. See parameter --stages and the README for more information' +
      ' on the individual stages. Important: please make sure that the input files for each' +
      ' repeat of the experiment are passed in alphabetical order.'
    ),
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
  )
  group_common = parser.add_argument_group('Common arguments')
  group_align = parser.add_argument_group('Stage 0_align')
  group_filter = parser.add_argument_group('Stage 1_filter')
  group_window = parser.add_argument_group('Stage 2_window')
  group_common.add_argument(
    '--names',
    nargs = '+',
    type = str,
    help = (
      'Identifiers to use for the input libraries.' +
      ' Must be given in the same order as the input files.' +
      ' If omitted, the names will be the input file names without the extension.' +
      ' Required for stages "0_align" and "1_filter" (but can be omitted because of the default).' +
      ' Important: please make sure that the input names are in alphabetical order.'
    ),
  )
  group_common.add_argument(
    '--output',
    type = common_utils.check_dir_output,
    help = (
      'Working output directory. Required for all stages.' +
      ' If running stages separately, the same directory must be given for all stages.' +
      ' The base part of the directory name is used as a unique identifier for the experiment' +
      ' so should not conflict with other experiments.'
    ),
    required = True,
  )
  group_common.add_argument(
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
      ' Will use all SAM files placed in OUTPUT.' +
      ' 2_window: Extract the nucleotide sequences and alignment from around the DSB site.' +
      ' 3_variation: Split alignment windows into individual variations (used for plotting variation-position histograms).'
    )
  )
  group_common.add_argument(
    '--ref',
    type = common_utils.check_file,
    help = (
      'Reference sequence file.' +
      ' Should contain a single nucleotide sequence.' +
      ' Must be the same sequence as used for alignment.' +
      f' Must be in FASTA format ({", ".join(constants.FASTA_EXT)}) or' +
      ' text format (all other extensions).' +
      ' Required for stages "0_align", "1_filter", and "2_window".'
    ),
  )
  group_common.add_argument(
    '--dsb',
    type = int,
    help = (
      'Position on reference sequence immediately left (5\') of DSB site.' +
      ' I.e., the DSB is between 1-based positions DSB_POS and DSB_POS + 1.' +
      ' Required for stage "1_filter" and "2_window".'
    ),
  )
  group_common.add_argument(
    '--quiet',
    action = 'store_true',
    help = 'If present, do no output verbose log message.',
  )
  group_common.add_argument(
    '--no_align',
    action = 'store_true',
    help = 'Shorthand for omitting the "0_align" stage (see --stages).',
  )
  group_align.add_argument(
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
      ' Required for stage "0_align". ' +
      ' If "0_align" is omitted and a custom alignnment approach is used, ' +
      ' then the aligned SAM files must be placed in OUTPUT directory.' +
      ' Important: please make sure that the input files are in alphabetical order.'
    ),
  )
  group_align.add_argument(
    '--bt2',
    type = str,
    default = '',
    help = (
      'Additional arguments to pass to Bowtie 2. Must be a single string.' +
      ' Required for stage "0_align" (but can be omitted because of default).'
    ),
  )
  group_filter.add_argument(
    '--min_len',
    type = int,
    help = (
      'Minimum length of read sequence to be considered. Reads shorter than this are discarded.' +
      ' Forced to be at least DSB_POS + 1 (also the default).' +
      ' Required for stage "1_filter" (but can be omitted because of default).'
    ),
  )
  group_filter.add_argument(
    '--max_sub',
    type = int,
    default = -1,
    help = (
      'Maximum number of substitutions allowed in the alignment.' +
      ' If the alignment has more substitutions, the read is rejected.' +
      ' A large number of substitutions may indicate that an alignment is invalid.' +
      ' Set to -1 to disable this check.' +
      ' Required for stage "1_filter" (but can be omitted because of default).'
    ),
  )
  group_filter.add_argument(
    '--rc',
    action = 'store_true',
    help = (
      'Set if the reads are expected to be reverse-complemented compared to the reference sequence.' +
      ' This is useful if the reads are from the opposite strand as the reference sequence.' +
      ' If this option is used, the alignments in the SAM file must have been aligned' +
      ' against the same reference sequence, and only alignments with the reverse-complement flag (16)' +
      ' will be accepted.' +
      ' Required for stage "1_filter" (but can be omitted because of default).'
    ),
  )
  group_filter.add_argument(
    '--consec',
    type = int,
    choices = [0, 1],
    default = 1,
    help = (
      'Set to 0 to disable the check that all in/dels must be consecutive.' +
      ' If this check is disabled, realignment is not performed.'
    ),
  )
  group_filter.add_argument(
    '--touch',
    type = int,
    choices = [0, 1],
    default = 1,
    help = (
      'Set to 0 to disable the check that some in/dels must touch the DSB.' +
      ' If this check is disabled, the "--dsb" option is ignored' +
      ' and realignment is not performed.'
    ),
  )
  group_filter.add_argument(
    '--reads',
    type = int,
    help = (
      'Total number reads in each experiment.' +
      ' This may be strictly greater than the number of reads in the input FASTQ' +
      ' files if some reads were discarded during preprocessing.' +
      ' The number of arguments must be the same as the number of ' +
      ' INPUTs. If not provided, the total reads remaining after NHEJ filtering are used.' +
      ' Required for stage "1_filter" (but can be omitted because of default).'
    ),
    nargs = '+',
  )
  group_window.add_argument(
    '--window',
    type = int,
    default = 10,
    help = (
      'Size of window around DSB site to obtain.' +
      ' The nucleotides at the positions' +
      ' {DSB_POS - WINDOW_SIZE + 1, ..., DSB_POS + WINDOW_SIZE} are obtained.' +
      ' The actual number of nucleotides obtained from each read may vary depending' +
      ' on the number of insertions/deletions/substitutions in the alignment.' +
      ' Required for stage "2_window" (but can be omitted because of default).'
    ),
  )
  group_window.add_argument(
    '--anchor',
    type = int,
    default = 20,
    help = (
      'Size of the anchor on the left/right of the obtained window to check for substitutions.' +
      ' Required for stage "2_window" (but can be omitted because of default).'
    ),
  )
  group_window.add_argument(
    '--anchor_var',
    type = int,
    nargs = 2,
    default = [1, 0],
    help = (
      'Maximum number of substitutions (arg 1) and indels (arg 2) allowed on the' +
      ' left/right anchor sequences. Reads with more than the allowed number of' +
      ' substitutions or indels on the left/right anchor will be discarded.' +
      ' This limit is applied to the left/right anchors separately.' +
      ' Required for stage "2_window" (but can be omitted because of default).'
    ),
  )
  group_window.add_argument(
    '--label',
    type = str,
    help = (
      'Label of the experiment to be used in plot legends.' +
      ' If not provided, the label is the basename of the OUTPUT directory.' +
      ' Required for stage "2_window" (but can be omitted because of default).'
    ),
  )
  args = vars(parser.parse_args())
  if args['names'] is not None:
    if len(args['names']) != len(args['input']):
      raise Exception('Number of NAMES must be the same as the number of INPUTs.')
    if not all(x == y for x, y in zip(args['names'], sorted(args['names']))):
      raise Exception('NAMES must be in alphabetical order.')
  f_names = [file_names.get_file_name(x) for x in args['input']]
  if not all(x == y for x, y in zip(f_names, sorted(f_names))):
    raise Exception('INPUT must have file names in alphabetical order.')
  if args['no_align']:
    args['stages'] = [x for x in args['stages'] if (x != '0_align')]
  del args['no_align']
  args['consec'] = bool(args['consec'])
  args['touch'] = bool(args['touch'])
  return args

def do_0_align(
  input,
  names,
  output,
  ref,
  bt2,
):
  if input is None:
    raise Exception('INPUT must be provided for stage "0_align".')
  if ref is None:
    raise Exception('REF must be provided for stage "0_align".')
  if names is None:
    names = [file_names.get_file_name(x) for x in input]
  bowtie2_build_file = file_names.bowtie2_build(output)
  file_utils.make_parent_dir(bowtie2_build_file)
  log_utils.log_output('Bowtie2 build file: ' + bowtie2_build_file)
  os.system(f'bowtie2-build-s {ref} {bowtie2_build_file} --quiet')

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
    bowtie2_command = f'bowtie2-align-s --no-hd {flags} {bt2} -x {bowtie2_build_file} {input[i]} -S {sam_file} --quiet'
    log_utils.log('Bowtie 2 command: ' + bowtie2_command)
    os.system(bowtie2_command)
  log_utils.blank_line()

def do_1_filter(
  names,
  output,
  reads,
  ref,
  dsb,
  min_len,
  max_sub,
  rc,
  consec,
  touch,
  quiet,
):
  if ref is None:
    raise Exception('REF must be provided for stage "1_filter".')
  if dsb is None:
    raise Exception('DSB must be provided for stage "1_filter".')
  if min_len is None:
    min_len = dsb + 1

  if names is None:
    # Infer names from the SAM files.
    names = sorted([
      file_names.get_file_name(x)
      for x in glob.glob(os.path.join(output, '*.sam'))
    ])
  input = [file_names.sam_file(output, x) for x in names]

  min_len = max(dsb + 1, min_len)
  filter.main(
    input = input,
    output = [
      file_names.filter(output, 'accepted'),
      file_names.filter(output, 'rejected'),
    ],
    debug = file_names.filter(output, 'debug'),
    names = names,
    reads = reads,
    ref = ref,
    dsb = dsb,
    min_len = min_len,
    max_sub = max_sub,
    rc = rc,
    consec = consec,
    touch = touch,
    quiet = quiet,
  )

def do_2_window(
  output,
  ref,
  dsb,
  window,
  anchor,
  anchor_var,
  label,
):
  if ref is None:
    raise Exception('REF must be provided for stage "2_window".')
  if dsb is None:
    raise Exception('DSB must be provided for stage "2_window".')
  if window is None:
    raise Exception('WINDOW must be provided for stage "2_window".')
  if anchor is None:
    raise Exception('ANCHOR must be provided for stage "2_window".')
  if anchor_var is None:
    raise Exception('ANCHOR_VAR must be provided for stage "2_window".')

  name = file_names.get_file_name(output)
  if label is None:
    label = name

  for subst_type in constants.SUBST_TYPES:
    get_window.main(
      input = file_names.filter(output, 'accepted'),
      ref = ref,
      output = [
        file_names.window(output, subst_type),
        file_names.data_info(output),
      ],
      dsb = dsb,
      window = window,
      anchor = anchor,
      anchor_var = anchor_var,
      sub = subst_type,
      name = name,
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
  ref,
  dsb,
  min_len,
  max_sub,
  rc,
  consec,
  touch,
  window,
  anchor,
  anchor_var,
  reads,
  label,
  bt2,
  quiet,
  stages = STAGES,
):
  # Copy reference sequence file to output directory.
  if ref is not None:
    shutil.copy(ref, file_names.ref_seq_file(output))

  if '0_align' in stages:
    do_0_align(
      input = input,
      names = names,
      output = output,
      ref = ref,
      bt2 = bt2,
    )

  if '1_filter' in stages:
    do_1_filter(
      names = names,
      output = output,
      reads = reads,
      ref = ref,
      dsb = dsb,
      min_len = min_len,
      max_sub = max_sub,
      rc = rc,
      consec = consec,
      touch = touch,
      quiet = quiet,
    )

  if '2_window' in stages:
    do_2_window(
      output = output,
      ref = ref,
      dsb = dsb,
      window = window,
      anchor = anchor,
      anchor_var = anchor_var,
      label = label,
    )

  if '3_variation' in stages:
    do_3_variation(output = output)

def main(
  input,
  names,
  output,
  ref,
  dsb,
  min_len,
  max_sub,
  rc,
  consec,
  touch,
  window,
  anchor,
  anchor_var,
  reads,
  label,
  bt2,
  quiet,
  stages = STAGES,
):
  do_stages(
    input = input,
    names  = names,
    output = output,
    ref = ref,
    dsb = dsb,
    min_len = min_len,
    max_sub = max_sub,
    rc = rc,
    consec = consec,
    touch = touch,
    window = window,
    anchor = anchor,
    anchor_var = anchor_var,
    reads = reads,
    label = label,
    bt2 = bt2,
    quiet = quiet,
    stages = stages,
  )

# This allows the "DSBplot-preprocess" command to be run from the command line.
def entry_point():
  main(**parse_args())

if __name__ == '__main__':
  entry_point()
