import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

import pandas as pd
import argparse
import shutil

import common_utils
import library_constants
import file_names
import file_utils
import log_utils

import get_window

def parse_args():
  parser = argparse.ArgumentParser(
    description = (
      'Combine two individual experiment directories to make a comparison' +
      ' experiment directory for comparison graphs.' +
      ' The experiments must be compatible (have all the same attribute ' +
      ' except for the constructs which must be different).'
    )
  )
  parser.add_argument(
    '--input',
    type = common_utils.check_dir,
    help = (
      'Directory of individual experiment frequency data produced with "get_freq.py".'
    ),
    nargs = 2,
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
    choices = library_constants.SUBST_TYPES,
    help = 'Whether to use with files with/without substitutions.',
    required = True,
  )
  args = parser.parse_args()
  return args

def get_comparison_data(data_1, data_2, join_columns, freq_columns):
  data_1 = data_1[join_columns + freq_columns]
  data_2 = data_2[join_columns + freq_columns]

  data_comparison = pd.merge(
    data_1,
    data_2,
    how = 'outer',
    on = join_columns,
    suffixes = ['_1', '_2'],
  )

  freq_columns_comparison = (
    [x + '_1' for x in freq_columns] +
    [x + '_2' for x in freq_columns]
  )
  data_comparison[freq_columns_comparison] = (
    data_comparison[freq_columns_comparison].fillna(0)
  )
  return data_comparison

def write_comparison_data(
  input_dir_1,
  input_dir_2,
  output_dir,
  subst_type,
):
  for freq_type in [
    library_constants.FREQ_FILTER_MEAN
  ]:
    input_file_1 = file_names.window(input_dir_1, freq_type, subst_type)
    input_file_2 = file_names.window(input_dir_2, freq_type, subst_type)

    data_1 = file_utils.read_tsv(input_file_1)
    data_2 = file_utils.read_tsv(input_file_2)

    data = get_comparison_data(
      data_1,
      data_2,
      ['ref_align', 'read_align'],
      ['freq_mean'],
    )
    output_file = file_names.window(output_dir, freq_type, subst_type)
    file_utils.write_tsv(data, output_file)
    log_utils.log(output_file)

if __name__ == '__main__':
  # Parse args
  args = parse_args()

  # Load data info
  data_info_1 = file_utils.read_tsv_dict(file_names.data_info(args.input[0]))
  data_info_2 = file_utils.read_tsv_dict(file_names.data_info(args.input[1]))

  # Make sure the experiments are compatible
  if not all(
    data_info_1[x] == data_info_2[x]
    for x in [
      'cell_line',
      'dsb_type',
      'guide_rna',
      'strand',
      'control_type',
      'ref_seq_window',
      'version',
    ]
  ):
    raise Exception(f'Incompatible experiments:\n{data_info_1}\n{data_info_2}')

  # Make the data
  write_comparison_data(
    args.input[0],
    args.input[1],
    args.output,
    args.subst_type,
  )

  # Make the comparison info
  get_window.write_data_info(
    dir = args.output,
    format = library_constants.DATA_COMPARISON,
    cell_line = data_info_1['cell_line'],
    dsb_type = data_info_1['dsb_type'],
    guide_rna = data_info_1['guide_rna'],
    strand = data_info_1['strand'],
    constructs = [data_info_1['construct'], data_info_2['construct']],
    control_type = data_info_1['control_type'],
    version = data_info_1['version'],
    ref_seq_window = data_info_1['ref_seq_window'],
    ref_seq = None,
  )
