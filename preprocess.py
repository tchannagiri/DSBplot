import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'utils'))) # allow importing
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '1_process_nhej'))) # allow importing
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '2_get_window_data'))) # allow importing
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '3_get_graph_data'))) # allow importing
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '4_get_histogram_data'))) # allow importing
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '5_plot_graph'))) # allow importing
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '6_plot_histogram'))) # allow importing

import argparse

import file_names
import common_utils
import file_utils
import log_utils
import constants

import filter_nhej
import combine_repeat
import get_window
import get_freq
import get_freq_comparison
import get_graph_data
import get_histogram_data

STAGES_1 = ['0_align', '1_filter', '2_combine', '3_window']
STAGES_2 = ['3_comparison', '4_graph', '5_histogram']

def parse_args():
  parser = argparse.ArgumentParser(
    description = (
      'Perform alignment and preprocessing for raw FASTQ data.' +
      ' This is script if broken in separate stages so that each stage' +
      ' can be run separately. However, the stages must be run in the correct order indicated' +
      ' by their prefix numbers. See --stages for more information.'
    )
  )
  parser.add_argument(
    '--input',
    nargs = '+',
    type = common_utils.check_file,
    help = (
      'Input FASTQ files of raw reads.' +
      ' Each file is considered a repeat of the same experiment.'
    ),
    required = True,
  )
  parser.add_argument(
    '--output',
    type = common_utils.check_dir_output,
    help = 'Output directory.',
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
      ' making comparison graphs.' +
      ' 4_graph: Preprocess the data further for graph construction and plotting.' +
      ' 5_hist: Preprocess the data further for histogram construction and plotting.' +
      ' For a more detailed explanation of the preprocessing pipeline, see the README.'
    )
  )
  parser.add_argument(
    '--ref_seq_file',
    type = common_utils.check_file,
    help = 'FASTA file with a single nucleotide sequence.',
    required = True,
  )
  parser.add_argument(
    '--dsb_pos',
    type = int,
    required = True,
    help = (
      'Position on reference sequence immediately left (5'') of DSB site.' +
      ' I.e., the DSB is between position DSB_POS and DSB_POS + 1.'
    ),
  )
  parser.add_argument(
    '--min_length',
    type = int,
    default = -1,
    help = (
      'Minimum length of read sequence to be considered.' +
      ' Reads shorter than this are discarded.' +
      ' Forced to be at least DSB_POS + 1.'
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
      ' on the number of insertions/deletions/substitutions in the alignment.'
    ),
  )
  parser.add_argument(
    '--anchor_size',
    type = int,
    default = 20,
    help = (
      'Size of the anchor on the left/right of the extracted window' +
      ' to check for mismatches.'
    ),
  )
  parser.add_argument(
    '--anchor_mismatches',
    type = int,
    default = 1,
    help = (
      'Maximum number of mismatches allowed on the left/right anchor sequences.\n'
      'Reads with more than the allowed number of mismatches on the left/right anchor\n'
      'will be discarded. This limit is applied to the left/right anchors separately.'
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
      ' INPUTs. If not provided, the total reads remaining after NHEJ filtering are used.'
    ),
    nargs = '+',
  )
  parser.add_argument(
    '--freq_min',
    type = float,
    default = 1e-5,
    help = (
      f'Minimum frequency for output in' +
      f' windows_{constants.FREQ_FILTER_MEAN}.' +
      f' Sequences with frequencies <= this are discarded.'
    ),
  )
  parser.add_argument(
    '--label',
    type = str,
    help = 'Label of the experiment to be used in plot legends.',
    required = True,
  )
  parser.add_argument(
    '--quiet',
    action = 'store_true',
    help = 'If present, do no output verbose log message.',
  )
  args = vars(parser.parse_args())
  args['min_length'] = max(args['dsb_pos'] + 1, args['min_length'])
  return args

def get_sam_file(output, i):
  return os.path.join(file_names.sam_dir(output), f'{i}.sam')

def get_filter_nhej_file(output, i):
  return os.path.join(file_names.filter_nhej_dir(output), f'{i}.tsv')

def do_0_align(
  input,
  output,
  ref_seq_file,
):
  bowtie2_build_file = file_names.bowtie2_build(output)
  file_utils.make_parent_dir(bowtie2_build_file)
  log_utils.log_input('Bowtie2 build file: ' + bowtie2_build_file)
  os.system(f'bowtie2-build-s {ref_seq_file} {bowtie2_build_file} --quiet')

  for i, input_1 in enumerate(input, 1):
    sam_file = get_sam_file(output, i)
    file_utils.make_parent_dir(sam_file)
    os.system(f'bowtie2-align-s -x {bowtie2_build_file} {input_1} -S {sam_file} --quiet')
    log_utils.log_output('Bowtie2 SAM file: ' + sam_file)
    log_utils.new_line()

def do_1_filter_nhej(
  input,
  output,
  ref_seq_file,
  dsb_pos,
  min_length,
  quiet,
):
  for i in range(1, len(input) + 1):
    filter_nhej.main(
      ref_seq_file = ref_seq_file,
      sam_file = get_sam_file(output, i),
      output = get_filter_nhej_file(output, i),
      dsb_pos = dsb_pos,
      min_length = min_length,
      quiet = quiet,
    )

def do_2_combine_repeat(
  input,
  output,
  quiet,
):
  combine_repeat.main(
    input = [get_filter_nhej_file(output, i) for i in range(1, len(input) + 1)],
    column_names = [f'r{i}' for i in range(1, len(input) + 1)],
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
      input = [input_1, input_2],
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
      input = input,
      output = output,
      ref_seq_file = ref_seq_file,
      dsb_pos = dsb_pos,
      min_length = min_length,
      quiet = quiet,
    )

  if '2_combine' in stages:
    do_2_combine_repeat(
      input = input,
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
 

if __name__ == '__main__':
  main(**parse_args())
