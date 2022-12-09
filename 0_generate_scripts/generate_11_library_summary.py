import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import file_utils
import library_constants
import file_names
import generate_constants
import generate_01_filter_nhej
import generate_02_combine_repeat
import generate_03_get_window

import pandas as pd


def get_total_reads(library_info):
  return library_info['total_reads']

def get_filter_nhej_counts(library_info):
  return file_utils.read_tsv(
    generate_01_filter_nhej.get_output_file(library_info['name'])
  )['Count'].sum()

def get_combine_repeat_counts(library_info):
  return file_utils.read_tsv(
    generate_02_combine_repeat.get_output_file(library_info['name_experiment'])
  )['Count_' + library_info['library']].sum()

def get_window_counts(library_info):
  return file_utils.read_tsv(file_names.window(
    generate_03_get_window.get_output_dir(library_info['name_experiment']),
    library_constants.COUNT,
    library_constants.SUBST_WITHOUT
  ))['count_' + library_info['library']].sum()

def get_window_counts_filter(library_info):
  freq = file_utils.read_tsv(file_names.window(
    generate_03_get_window.get_output_dir(library_info['name_experiment']),
    library_constants.FREQ_FILTER,
    library_constants.SUBST_WITHOUT
  ))['freq_' + library_info['library']]
  freq = freq[freq > 1e-5].sum()
  return round(library_info['total_reads'] *freq)


if __name__ == '__main__':
  # Get the raw counts
  count_data_list = []
  for library_info in generate_constants.LIBRARY_INFO.to_dict('records'):
    if (
      (library_info['version'] != library_constants.VERSION_MERGED) and
      (library_info['control_type'] != library_constants.CONTROL_30BPDOWN)
    ):
      log_utils.log(library_info['name'])
      count_data = {'library': library_info['name']}
      count_data['total'] = get_total_reads(library_info)
      count_data['filter_nhej'] = get_filter_nhej_counts(library_info)
      count_data['combine_repeat'] = get_combine_repeat_counts(library_info)
      count_data['window_extract'] = get_window_counts(library_info)
      count_data['window_filter'] = get_window_counts_filter(library_info)
      count_data_list.append(count_data)

  # The raw count data
  data_count = pd.DataFrame.from_records(count_data_list)

  # Sucessive differences in the counts
  freq_columns = data_count.columns[data_count.columns != 'library']
  data_count_diff = data_count.copy()
  for i in range(len(freq_columns) - 1, 0, -1):
    data_count_diff[freq_columns[i]] -= data_count_diff[freq_columns[i - 1]]

  data_dict = {
    'count': data_count,
    'count_diff': data_count_diff
  }

  # The same data as frequencies instead of counts
  for key in ['count', 'count_diff']:
    data = data_dict[key].copy()
    total_reads = data['total'].copy()
    for column in freq_columns:
      data[column] = data[column] / total_reads
    data_dict[key.replace('count', 'freq')] = data
  
  # Write data to files
  for key, data in data_dict.items():
    file_out = os.path.join(
      generate_constants.OUTPUT_DIR['library_summary'],
      'library_summary_' + key + '.xlsx',
    )
    file_utils.make_parent_dir(file_out)
    data.to_excel(file_out)
    log_utils.log(file_out)
