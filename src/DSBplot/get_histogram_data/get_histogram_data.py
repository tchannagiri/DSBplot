import pandas as pd
import argparse

import DSBplot.utils.file_names as file_names
import DSBplot.utils.file_utils as file_utils
import DSBplot.utils.alignment_utils as alignment_utils
import DSBplot.utils.log_utils as log_utils
import DSBplot.utils.common_utils as common_utils
import DSBplot.utils.constants as constants

def parse_args():
  parser = argparse.ArgumentParser(
    description = (
      'Precomputes the necessary data for the variation-position histograms.' +
      ' Uses as input the output of the "get_graph_data" stage.'
    )
  )
  parser.add_argument(
    '--input',
    type = common_utils.check_dir,
    help = 'Directory where output of "get_graph_data" stage is located.',
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
    choices = constants.SUBST_TYPES,
    help = 'Whether to process the files with/without substitutions.',
  )
  return vars(parser.parse_args())

def split_seqs_into_variations(sequence_data):
  """
    Split sequences into their individual variations.
  """
  variation_data = sequence_data.copy()
  variation_data = variation_data.rename({'id': 'seq_id'}, axis='columns')

  # Explode the different variations
  variation_data = variation_data.rename(
    {
      'ref_align': 'ref_align',
      'read_align': 'read_align',
      'mid_align': 'mid_align',
      'variation_type': 'variation_type_seq',
    },
    axis = 'columns',
  )
  if variation_data.shape[0] > 0:
    variation_data['variation_info'] = (
      variation_data.apply(
        lambda x: [
          dict(zip(['variation_pos', 'variation_type', 'variation_letter'], info_tuple))
          for info_tuple in alignment_utils.get_variation_info(x['ref_align'], x['read_align'])
        ],
        axis = 'columns',
      )
    )
  else:
    variation_data['variation_info'] = []

  variation_data = variation_data.loc[variation_data['variation_info'].apply(len) > 0]
  variation_data = variation_data.explode('variation_info', ignore_index=True)
  if variation_data.shape[0] > 0:
    variation_data[['variation_pos', 'variation_type', 'variation_letter']] = (
      variation_data['variation_info'].apply(
        lambda x: pd.Series([x['variation_pos'], x['variation_type'], x['variation_letter']])
      )
    )
  else:
    variation_data['variation_pos'] = []
    variation_data['variation_type'] = []
    variation_data['variation_letter'] = []
  variation_data = variation_data.drop('variation_info', axis='columns')

  # Add the variation IDs
  variation_data['id'] = (
    variation_data['seq_id'] + '_V' +
    (variation_data.groupby('seq_id')['seq_id'].cumcount() + 1).astype(str)
  )

  # Reorder the columns
  columns = ['id', 'seq_id']
  columns += list(variation_data.columns[~variation_data.columns.isin(columns)])
  variation_data = variation_data[columns]

  return variation_data

def write_variation(input_dir, output_dir, subst_type):
  """
    Make data on individual variations and write to file.
    Sequence data should already be created.
  """
  sequence_data = file_utils.read_tsv(file_names.sequence_data(input_dir, subst_type))
  variation_data = split_seqs_into_variations(sequence_data)

  out_file_name = file_names.variation(output_dir, subst_type)
  file_utils.write_tsv(variation_data, out_file_name)
  log_utils.log_output(out_file_name)

def write_variation_grouped(output_dir, subst_type):
  """
    Groups the variation data by position and number of variations on sequence.
    This data is used for the 3D variation-position histograms.
    Variation data should have already been made.
  """
  variation_data = file_utils.read_tsv(file_names.variation(output_dir, subst_type))
  data_info = file_utils.read_tsv_dict(file_names.data_info(output_dir))
  if data_info['format'] == 'comparison':
    freq_column_list = ['freq_mean_1', 'freq_mean_2']
  elif data_info['format'] == 'individual':
    freq_column_list = ['freq_mean']
  else:
    raise Exception('Invalid format: ' + data_info['format'])
  variation_data = variation_data[[
    'id',
    *freq_column_list,
    'dist_ref',
    'variation_pos',
    'variation_type',
    'variation_letter',
  ]]
  
  aggregate_args = {}
  for freq_column in freq_column_list:
    aggregate_args[freq_column] = (freq_column, 'sum')
  aggregate_args['var_id'] = ('id', common_utils.join_with_comma)
  variation_data = variation_data.groupby([
    'dist_ref',
    'variation_pos',
    'variation_type',
    'variation_letter',
  ]).aggregate(**aggregate_args).reset_index()

  variation_data['freq_mean'] = variation_data[freq_column_list].max(axis='columns')
  variation_data = variation_data.sort_values('freq_mean', ascending=False)
  variation_data['id'] = (
    'GV' + pd.Series(range(1, variation_data.shape[0] + 1), dtype=str)
  )

  variation_data = variation_data[
    ['id', 'var_id'] +
    list(variation_data.columns[~variation_data.columns.isin(['id', 'var_id'])])
  ]
  out_file_name = file_names.variation_grouped(output_dir, subst_type)
  file_utils.write_tsv(variation_data, out_file_name)
  log_utils.log_output(out_file_name)

def main(
  input,
  output,
  subst_type,
):
  log_utils.log_input(input)

  # copy data info files
  if input != output:
    input_data_info_file = file_names.data_info(input)
    output_data_info_file = file_names.data_info(output)
    file_utils.copy(input_data_info_file, output_data_info_file)
    log_utils.log_output(output_data_info_file)

  write_variation(input, output, subst_type)
  write_variation_grouped(output, subst_type)
  log_utils.new_line()

if __name__ == '__main__':
  main(**parse_args())
