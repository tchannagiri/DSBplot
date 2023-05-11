import argparse

import utils.file_utils as file_utils
import utils.common_utils as common_utils
import utils.log_utils as log_utils
import utils.file_names as file_names
import utils.constants as constants

def parse_args():
  parser = argparse.ArgumentParser(
    description = (
      f'Convert the raw read counts in the input data into frequencies' +
      f' using the input total reads. Outputs 3 files:' +
      f' (1) "windows_{constants.FREQ}.tsv": contains the all the sequences' +
      f' with the counts converted to frequencies.' +
      f' (2) "windows_{constants.FREQ_FILTER}.tsv":' +
      f' the previous file with the sequences removed whose frequency is <= FREQ_MIN' +
      f' in any of the repeats.' +
      f' (3) "windows_{constants.FREQ_FILTER_MEAN}.tsv": ' +
      f' contains the means of the frequencies in the previous file (over all repeats).'
    )
  )
  parser.add_argument(
    '--input',
    type = common_utils.check_dir,
    help = 'Directory with output from "get_window.py" or "get_merged.py".',
    required = True,
  )
  parser.add_argument(
    '--total_reads',
    type = int,
    help = (
      'Total reads for each file.'
      ' Must be the same number of arguments as the number of ' +
      ' "Count" columns in INPUT. If not provided, the total reads are' +
      ' are calculated by taking the sum of the "Count" columns in INPUT.'
    ),
    nargs = '+',
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
    default = constants.SUBST_WITHOUT,
    choices = [
      constants.SUBST_WITH,
      constants.SUBST_WITHOUT,
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
      f' windows_{constants.FREQ_FILTER_MEAN}.' +
      f' Sequences with frequencies <= FREQ_MIN are discarded.'
    ),
  )
  return vars(parser.parse_args())

def main(input, output, subst_type, total_reads, freq_min):
  input_file = file_names.window(
    input,
    constants.COUNT,
    subst_type,
  )
  log_utils.log_input(input_file)

  data = file_utils.read_tsv(input_file)

  count_cols = data.columns[data.columns.str.startswith('count_')]

  if total_reads is None:
    total_reads = data[count_cols].sum(axis='index').to_list()
  if len(count_cols) != len(total_reads):
    raise Exception(
      f'Expected {len(total_reads)} count columns.' +
      f' Got {len(count_cols)}.'
    )
  freq_cols = count_cols.str.replace('count_', 'freq_')

  data[freq_cols] = data[count_cols].divide(total_reads, axis='columns')
  data = data.drop(count_cols, axis='columns')

  output_file = file_names.window(
    input,
    constants.FREQ,
    subst_type,
  )
  file_utils.write_tsv(data, output_file)
  log_utils.log_output(output_file)

  data = data.loc[data[freq_cols].min(axis='columns') > freq_min]
  
  output_file = file_names.window(
    output,
    constants.FREQ_FILTER,
    subst_type,
  )
  file_utils.write_tsv(data, output_file)
  log_utils.log_output(output_file)

  data['freq_mean'] = data[freq_cols].mean(axis='columns')
  data = data.sort_values('freq_mean', ascending = False)
  output_file = file_names.window(
    input,
    constants.FREQ_FILTER_MEAN,
    subst_type,
  )
  file_utils.write_tsv(
    data[['ref_align', 'read_align', 'freq_mean']],
    output_file,
  )
  log_utils.log_output(output_file)

  log_utils.new_line()
  
if __name__ == '__main__':
  main(**parse_args())
