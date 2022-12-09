import argparse
import sys
import os
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

import file_utils
import common_utils
import log_utils
import file_names
import library_constants

def main():
  parser = argparse.ArgumentParser(
    description = (
      f'Convert the raw read counts in the input data into frequencies' +
      f' using the input total reads. Outputs 3 files:' +
      f' (1) windows_{library_constants.FREQ}.tsv: contains the all the sequences' +
      f' with the counts converted to frequencies.' +
      f' (2) windows_{library_constants.FREQ_FILTER}.tsv:' +
      f' the previous file with the sequences removed whose frequency is <= FREQ_MIN' +
      f' in any of the repeats.' +
      f' (3) windows_{library_constants.FREQ_FILTER_MEAN}.tsv: ' +
      f' contains the means of the frequencies in the previous file (over all repeats).'
    )
  )
  parser.add_argument(
    '--input',
    type = common_utils.check_dir,
    help = 'Directory with output from get_window.py or get_merged.py.',
    required = True,
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
    '--output',
    type = common_utils.check_dir_output,
    help = 'Output directory.',
    required = True,
  )
  parser.add_argument(
    '--subst_type',
    type = str,
    default = 'without',
    choices = [
      library_constants.SUBST_WITH,
      library_constants.SUBST_WITHOUT,
    ],
    help = 'Whether to process the files with/without substitutions.',
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
  args = parser.parse_args()

  input_file = file_names.window(
    args.input,
    library_constants.COUNT,
    args.subst_type,
  )
  log_utils.log(input_file)
  log_utils.log('------>')

  data = file_utils.read_tsv(input_file)

  count_cols = data.columns[data.columns.str.startswith('count_')]
  if len(count_cols) != len(args.total_reads):
    raise Exception(
      f'Expected {len(args.total_reads)} count columns.' +
      f' Got {len(count_cols)}.'
    )
  freq_cols = count_cols.str.replace('count_', 'freq_')

  data[freq_cols] = data[count_cols].divide(args.total_reads, axis='columns')
  data = data.drop(count_cols, axis='columns')

  output_file = file_names.window(
    args.input,
    library_constants.FREQ,
    args.subst_type,
  )
  file_utils.write_tsv(data, output_file)
  log_utils.log(output_file)

  data = data.loc[data[freq_cols].min(axis='columns') > args.freq_min]
  
  output_file = file_names.window(
    args.input,
    library_constants.FREQ_FILTER,
    args.subst_type,
  )
  file_utils.write_tsv(data, output_file)
  log_utils.log(output_file)

  data['freq_mean'] = data[freq_cols].mean(axis='columns')
  data = data.sort_values('freq_mean', ascending = False)
  output_file = file_names.window(
    args.input,
    library_constants.FREQ_FILTER_MEAN,
    args.subst_type,
  )
  file_utils.write_tsv(
    data[['ref_align', 'read_align', 'freq_mean']],
    output_file,
  )
  log_utils.log(output_file)

  log_utils.new_line()
  
if __name__ == '__main__':
  main()
