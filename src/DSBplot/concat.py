import os
import argparse

import pandas as pd

import DSBplot.utils.constants as constants
import DSBplot.utils.common_utils as common_utils
import DSBplot.utils.file_utils as file_utils
import DSBplot.utils.log_utils as log_utils
import DSBplot.utils.file_names as file_names
import DSBplot.preprocess as preprocess

import numpy as np

def parse_args():
  parser = argparse.ArgumentParser(
    description = (
      'Concatenate together library files from samples that have been sequenced multiple times.' +
      ' Operates on the output of the preprocessing.'
    ),
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
  )
  parser.add_argument(
    '--input',
    type = common_utils.check_dir,
    help = (
      ' Input directory, which is the output of "preprocess.py".' +
      ' There must be the same number of columns in all input data sets.' +
      ' The count columns must be in the same order they should be merged in.'
    ),
    nargs = '+',
    required = True,
  )
  parser.add_argument(
    '--output',
    type = common_utils.check_dir_output,
    help = (
      'Output directory.' +
      ' The basename of the output directory is used as a unique identifier ' +
      ' for the experiment so it should not conflict with other experiments.'
    ),
    required = True,
  )
  parser.add_argument(
    '--label',
    type = str,
    help = (
      'Label of the output.' +
      ' By defualt use the basename of the input directory.'
    ),
  )
  parser.add_argument(
    '--names',
    type = str,
    help = (
      'Names of the new concatenated columns.' +
      ' Must be the same as the number of "count" columns in the input.'
    ),
    nargs = '+',
    required = True,
  )
  parser.add_argument(
    '--total_reads',
    type = int,
    nargs = '+',
    help = (
      'Total number of reads in each library.' +
      ' This is used to calculate the frequency of each read.' +
      ' Must be arranged in the same order as the input files' +
      ' and the columns in each input file.' +
      ' If not provided, the total number of reads is calculated from the input files.'
    ),
  )
  return vars(parser.parse_args())

def main(
  input,
  output,
  label,
  names,
  total_reads,
):
  count_cols = ['count_' + x for x in names]
  freq_cols = ['freq_' + x for x in names]

  for subst_type in [constants.SUBST_WITH, constants.SUBST_WITHOUT]:
    data_list = []
    for i in range(len(input)):
      in_file = file_names.window(input[i], subst_type)
      data = file_utils.read_csv(in_file)
      log_utils.log_input(in_file)
      count_cols_old = [x for x in data.columns if x.startswith('count_')]
      if len(count_cols_old) != len(names):
        raise Exception('Different number of count and freq columns in ' + input[i])
      data = data.rename(
        dict(zip(count_cols_old, count_cols)),
        axis = 'columns',
      )
      data_list.append(data[['ref_align', 'read_align'] + count_cols])

    if total_reads is None:
      total_reads_1 = (
        np.concatenate([x[count_cols].sum() for x in data_list])
        .tolist()
      )
    else:
      total_reads_1 = total_reads
    if len(total_reads_1) != (len(input) * len(names)):
      raise Exception(
        'Incorrect number of total reads. Expected {}.'
        .format(len(input) * len(names))
      )
    
    total_reads_1 = (
      np.array(total_reads_1).reshape(len(input), len(names))
      .sum(axis=0).tolist()
    )

    data = pd.concat(data_list, axis='index', ignore_index=True)
    data = data.groupby(['ref_align', 'read_align']).sum().reset_index()

    data[freq_cols] = data[count_cols].divide(total_reads_1, axis='columns')
    data['freq_mean'] = data[['freq_' + x for x in names]].mean(axis='columns')
    data = data.sort_values('freq_mean', ascending=False)
    data = data[
      ['ref_align', 'read_align', 'freq_mean'] +
      freq_cols +
      count_cols
    ]

    output_file = file_names.window(output, subst_type)
    file_utils.write_csv(data, output_file)
    log_utils.log_output(output_file)
    log_utils.blank_line()

  ref_seq_window = set(
    file_utils.read_csv_dict(file_names.data_info(x))['ref_seq_window']
    for x in input
  )
  if len(ref_seq_window) > 1:
    raise Exception(
      'Libraries begin merged have different window reference sequences.' +
      ' Got {}.'.format(ref_seq_window)
    )
  
  name = os.path.basename(output)
  if label is None:
    label = name
  data_info = common_utils.make_data_info(
    format = 'individual',
    names = [name],
    labels = [label],
    ref_seqs = [None],
    ref_seq_window = ref_seq_window.pop(),
  )
  output_file = file_names.data_info(output)
  file_utils.write_csv(pd.DataFrame(data_info, index=[0]), output_file)
  log_utils.log_output(output_file)
  log_utils.blank_line()

  preprocess.do_3_variation(output)

# This allows the "DSBplot-concat" command to be run from the command line.
def entry_point():
  main(**parse_args())

if __name__ == '__main__':
  main(**parse_args())
