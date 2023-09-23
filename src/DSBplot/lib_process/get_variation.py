import pandas as pd
import argparse

import DSBplot.utils.file_utils as file_utils
import DSBplot.utils.alignment_utils as alignment_utils
import DSBplot.utils.log_utils as log_utils
import DSBplot.utils.common_utils as common_utils

def parse_args():
  parser = argparse.ArgumentParser(
    description = (
      'Precomputes the necessary data for the variation-position histograms.' +
      ' Splits the alignments windows into individual variations.' +
      ' Uses as input the output of the window extract stage.'
    )
  )
  parser.add_argument(
    '--input',
    type = common_utils.check_file,
    help = 'Input CSV file. Must be output of the window extract stage.',
    required = True,
  )
  parser.add_argument(
    '--output',
    type = common_utils.check_file_output,
    help = 'Output CSV file.',
    required = True,
  )
  return vars(parser.parse_args())

def split_seqs_into_variations(window_data):
  """
    Split sequences into their individual variations.
  """
  value_cols = [
    x for x in window_data.columns
    if (x.startswith('freq_') or x.startswith('count_'))
  ]
  variation_data = window_data.copy()
  variation_data = variation_data[['ref_align', 'read_align'] + value_cols]

  # Separate into individual variations
  new_data = []
  for row in variation_data.to_dict('records'):
    num_var = sum(alignment_utils.count_variations(row['ref_align'], row['read_align']))
    if num_var > 0:
      for var in alignment_utils.get_variation_info(row['ref_align'], row['read_align']):
        new_data.append({
          **{k: row[k] for k in value_cols},
          'num_var': num_var,
          'var_pos': var[0],
          'var_type': var[1],
          'var_letter': var[2],
        })
  if len(new_data) > 0:
    variation_data = pd.DataFrame.from_records(new_data)
  else:
    variation_data = pd.DataFrame(
      columns = (
        value_cols + [
          'freq_mean',
          'num_var',
          'var_pos',
          'var_type',
          'var_letter',
        ]
      )
    )
  
  variation_data = variation_data.groupby([
    'num_var',
    'var_pos',
    'var_type',
    'var_letter',
  ]).sum().reset_index()

  variation_data = common_utils.sort_by_count(
    variation_data,
    [x for x in value_cols if x.startswith('count_')],
    ['num_var', 'var_type', 'var_pos', 'var_letter'],
  )

  variation_data = variation_data[
    value_cols + [
      'num_var',
      'var_pos',
      'var_type',
      'var_letter',
    ]
  ]

  return variation_data

def main(input, output):
  log_utils.log_input(input)
  window_data = file_utils.read_csv(input)
  variation_data = split_seqs_into_variations(window_data)

  file_utils.write_csv(variation_data, output)
  log_utils.log_output(output)

if __name__ == '__main__':
  main(**parse_args())
