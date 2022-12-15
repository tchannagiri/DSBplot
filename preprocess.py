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
import library_constants

import filter_nhej
import combine_repeat
import get_window
import get_freq
import get_graph_data
import get_histogram_data

def parse_args():
  parser = argparse.ArgumentParser(
    description = 'Perform alignment and preprocessing for raw FASTQ data.'
  )
  parser.add_argument(
    '--input',
    nargs = '+',
    type = common_utils.check_file,
    help = (
      'Input FASTQ files of raw reads.' +
      ' Each file is considered a repeat of the same experiment.'
    ),
    required = True
  )
  parser.add_argument(
    '--output',
    type = common_utils.check_dir_output,
    help = 'Output directory.',
    required = True
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
      'Position on reference sequence immediately left of DSB site.' +
      ' Ie. the DSB is between position DSB_POS and DSB_POS + 1.'
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
      ' The actual number of nucleotides extracted may vary depending' +
      ' on how many insertions/deletion the alignment of the sequence has.'
    ),
  )
  parser.add_argument(
    '--anchor_size',
    type = int,
    default = 20,
    help = (
      'Size of anchor on left/right of the window to check for mismatches.'
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
      'Total reads for each file.'
      ' Must be the same number of arguments as the number of ' +
      ' Count columns in INPUT.'
    ),
    nargs = '+',
    required = True,
  )
  parser.add_argument(
    '--freq_min',
    type = float,
    default = 1e-5,
    help = (
      f'Minimum frequency for output in' +
      f' windows_{library_constants.FREQ_FILTER_MEAN}.' +
      f' Sequences with frequences <= this are discarded.'
    ),
  )
  parser.add_argument(
    '--label',
    type = common_utils.check_dir_output,
    help = 'Label of the experiment to be used in plot titles.',
    required = True,
  )
  parser.add_argument(
    '--quiet',
    action = 'store_true',
    help = 'If present, do no output verbose log message.',
  )
  return vars(parser.parse_args())

def main(
  input,
  output,
  ref_seq_file,
  dsb_pos,
  window_size,
  anchor_size,
  anchor_mismatches,
  total_reads,
  freq_min,
  label,
  quiet,
):
  bowtie2_build_file = os.path.join(file_names.bowtie2_build(output))
  file_utils.make_parent_dir(bowtie2_build_file)
  log_utils.log_input('Bowtie2 build file: ' + bowtie2_build_file)
  os.system(f'bowtie2-build-s {ref_seq_file} {bowtie2_build_file} --quiet')

  filter_nhej_file_list = []
  for i, input_1 in enumerate(input, 1):
    sam_file = os.path.join(file_names.sam_dir(output), f'{i}.sam')
    file_utils.make_parent_dir(sam_file)
    os.system(f'bowtie2-align-s -x {bowtie2_build_file} {input_1} -S {sam_file} --quiet')
    log_utils.log_output('Bowtie2 SAM file: ' + sam_file)

    filter_nhej_file = os.path.join(file_names.filter_nhej_dir(output), f'{i}.tsv')
    filter_nhej.main(
      ref_seq_file = ref_seq_file,
      sam_file = sam_file,
      output = filter_nhej_file,
      dsb_pos = dsb_pos,
      quiet = quiet,
    )
    filter_nhej_file_list.append(filter_nhej_file)
  
  combine_repeat_file = os.path.join(file_names.combine_repeat_dir(output), 'out.tsv')
  combine_repeat.main(
    input = filter_nhej_file_list,
    column_names = [f'r{i}' for i in range(1, len(input) + 1)],
    output = combine_repeat_file,
    quiet = quiet,
  )

  window_dir = file_names.window_dir(output)
  graph_dir = file_names.graph_dir(output)
  histogram_dir = file_names.histogram_dir(output)
  for subst_type in library_constants.SUBST_TYPES:
    get_window.main(
      input = combine_repeat_file,
      output = file_names.window_dir(output),
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
    get_graph_data.main(
      input = window_dir,
      output = graph_dir,
      subst_type = subst_type,
    )
    get_histogram_data.main(
      input = graph_dir,
      output = histogram_dir,
      subst_type = subst_type,
    )

if __name__ == '__main__':
  sys.argv = [
    '',
    '--input',
    'data/0_fastq/db1_R1.fq',
    'data/0_fastq/db2_R1.fq',
    'data/0_fastq/db3_R1.fq',
    'data/0_fastq/db4_R1.fq',
    '--ref_seq_file',
    'data/0_ref_seq/1DSB_R1_branch.fa',
    '--dsb_pos',
    '50',
    '--output',
    'data/output',
    '--label',
    'db_R1',
    '--total_reads',
    '3000',
    '3000',
    '3000',
    '3000',
  ]
  main(**parse_args())

