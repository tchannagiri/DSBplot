# FIXME: DELETE FILE
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

def split_seqs_into_variations(window_data):
  """
    Split sequences into their individual variations.
  """
  variation_data = window_data.copy()
  variation_data = variation_data.set_index(
    ['ref_align', 'read_align']
  ).mean(axis='columns').rename('freq_mean').reset_index()

  # Separate into individual variations
  new_data = []
  for row in variation_data.to_dict('records'):
    dist_ref = sum(alignment_utils.count_variations(row['ref_align'], row['read_align']))
    if dist_ref > 0:
      for var in alignment_utils.get_variation_info(row['ref_align'], row['read_align']):
        new_data.append({
          'freq_mean': row['freq_mean'],
          'dist_ref': dist_ref,
          'variation_pos': var[0],
          'variation_type': var[1],
          'variation_letter': var[2],
        })
  if len(new_data) > 0:
    variation_data = pd.DataFrame.from_records(new_data)
  else:
    variation_data = pd.DataFrame(
      columns = [
        'freq_mean',
        'dist_ref',
        'variation_pos',
        'variation_type',
        'variation_letter',
      ]
    )
  
  variation_data = variation_data.groupby([
    'dist_ref',
    'variation_pos',
    'variation_type',
    'variation_letter',
  ]).sum().reset_index()

  variation_data = variation_data.sort_values('freq_mean', ascending=False)

  variation_data = variation_data[[
    'freq_mean',
    'dist_ref',
    'variation_pos',
    'variation_type',
    'variation_letter',
  ]]

  return variation_data

def write_variation(input_dir, output_dir, subst_type):
  """
    Make data on individual variations and write to file.
    Sequence data should already be created.
  """
  window_data = file_utils.read_tsv(
    file_names.window(input_dir, constants.FREQ, subst_type)
  )
  data_info = file_utils.read_tsv_dict(file_names.data_info(output_dir))
  if data_info['format'] != constants.DATA_INDIVIDUAL:
    raise Exception('Data format is not individual.')
  variation_data = split_seqs_into_variations(window_data)

  out_file_name = file_names.variation(output_dir, subst_type)
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
  log_utils.blank_line()

if __name__ == '__main__':
  main(**parse_args())
