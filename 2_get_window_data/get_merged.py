import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import argparse

import library_constants
import common_utils
import file_utils
import log_utils
import file_names

import pandas as pd
import argparse

def parse_args():
  parser = argparse.ArgumentParser(
    description = (
      'Merge multiple dsb-sequence-window tables.' +
      ' Currently only used to merge antisense experiments' +
      ' which are nearly (but not totally) technical replicates.'
    )
  )
  parser.add_argument(
    '--input',
    type = common_utils.check_dir,
    help = (
      ' Input directories with "windows_XXX.tsv" output of "get_windows.py".' +
      ' There must be the same number of columns in all' +
      ' input data sets.' +
      ' The count columns must be in the same' +
      ' order they should be merged in.'
    ),
    nargs = 2,
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
      library_constants.SUBST_WITHOUT
    ],
    help = 'Whether to keep or ignore substitutions.',
    required = True,
  )
  parser.add_argument(
    '--new_column_names',
    type = str,
    help = (
      'Names of the new merged columns.' +
      ' Must be the same as the number of count columns in the data.'
    ),
    nargs = '+',
    required = True,
  )
  return vars(parser.parse_args())

def main(
  input,
  output,
  subst_type,
  new_columns_names,
):
  for x in input:
    log_utils.log_input(x)

  data = [
    file_utils.read_tsv(
      file_names.window(x, library_constants.COUNT, subst_type)
    )
    for x in input
  ]
  column_names = []
  num_columns = data[0].shape[1]
  num_count_columns = num_columns - 2
  count_columns_temp = [str(x) for x in range(num_count_columns)]
  for i in range(len(data)):
    if data[i].shape[1] != num_columns:
      raise Exception('Different number of columns in data sets')
    count_columns = [x for x in data[i].columns if x.startswith('count_')]

    if len(count_columns) != len(new_columns_names):
      raise Exception(
        'Number of new names does not match number of count columns in ' +
        input[i]
      )
    
    column_names.append([x.replace('count_', '') for x in count_columns])
    data[i] = data[i].rename(
      dict(zip(count_columns, count_columns_temp), axis='columns'),
      axis = 'columns'
    )
  
  count_columns_new = ['count_' + x for x in new_columns_names]

  data = pd.concat(data, axis='index')
  data = data.rename(dict(zip(count_columns_temp, count_columns_new)), axis='columns')
  data = data.groupby(['ref_align', 'read_align']).sum()
  data['count_min'] = data.min(axis='columns')
  data = data.reset_index()
  data = data.sort_values(
    ['count_min', 'read_align'],
    ascending = [False, True],
  ).reset_index(drop=True).drop('count_min', axis='columns')

  output_file = file_names.window(
    output,
    library_constants.COUNT,
    subst_type,
  )
  file_utils.write_tsv(data, output_file)
  log_utils.log_output(output_file)

  data_info = [file_utils.read_tsv_dict(file_names.data_info(x)) for x in input]
  ref_seq_window = data_info[0]['ref_seq_window']
  if not(all(x['ref_seq_window'] == ref_seq_window for x in data_info)):
    raise Exception('Libraries begin merged have different window reference sequences.')
  data_info = data_info[0]
  data_info['ref_seq'] = None
  output_data_info_file = file_names.data_info(output)
  file_utils.write_tsv(pd.DataFrame(data_info, index=[0]), output_data_info_file)
  log_utils.log_output(output_data_info_file)

  log_utils.new_line()

if __name__ == '__main__':
  main(**parse_args())
