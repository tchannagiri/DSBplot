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
import filter_nhej

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
    '--output',
    type = common_utils.check_dir_output,
    help = 'Output directory.',
    required = True
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
  ref_seq_file,
  dsb_pos,
  output,
  label,
  quiet,
):
  bowtie2_build_dir = os.path.join(
    output,
    '0_bowtie2_build',
    'build',
  )
  file_utils.make_parent_dir(bowtie2_build_dir)
  log_utils.log('Bowtie2 build dir: ' + bowtie2_build_dir)
  os.system(f'bowtie2-build-s {ref_seq_file} {bowtie2_build_dir} --quiet')

  for i, input_1 in enumerate(input, 1):
    sam_file = os.path.join(output, 'sam', f'{i}.sam')
    file_utils.make_parent_dir(sam_file)
    log_utils.log('Bowtie2 SAM output file: ' + sam_file)
    os.system(f'bowtie2-align-s -x {bowtie2_build_dir} {input_1} -S {sam_file} --quiet')

    filter_nhej_dir = os.path.join(output, '1_filter_nhej', f'{i}.tsv')
    log_utils.log('Filter NHEJ dir: ' + filter_nhej_dir)
    filter_nhej.main(
      ref_seq_file = ref_seq_file,
      sam_file = sam_file,
      output = filter_nhej_dir,
      dsb_pos = dsb_pos,
      quiet = quiet,
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
  ]
  main(**parse_args())

