import os
import DSBplot.utils.file_utils as file_utils

def check_file(file_name):
  if os.path.exists(file_name):
    return file_name
  else:
    raise ValueError('File does not exist: ' + str(file_name))

def check_dir(dir_name):
  if os.path.isdir(dir_name):
    return dir_name
  else:
    raise ValueError('Not a directory: ' + str(dir_name))

def check_dir_output(dir_name):
  file_utils.make_dir(dir_name)
  return check_dir(dir_name)

def check_file_output(file_name):
  file_utils.make_parent_dir(file_name)
  return file_name

def join_with_comma(arr):
  return ','.join([str(x) for x in arr])

# Make the data info dict describing an input data set.
def make_data_info(
  format,
  names,
  labels,
  ref_seqs,
  ref_seq_window,
):
  data_info = {
    'format': format,
    'ref_seq_window': ref_seq_window,
  }
  if format == 'individual':
    if len(names) != 1:
      raise Exception(f'Expected 1 name for individual format. Got: {len(names)}')
    if len(labels) != 1:
      raise Exception(f'Expected 1 name for individual format. Got: {len(labels)}')
    if len(ref_seqs) != 1:
      raise Exception(f'Expected 1 reference sequence for individual format. Got: {len(ref_seqs)}')
    data_info['name'] = names[0]
    data_info['label'] = labels[0]
    data_info['ref_seq'] = ref_seqs[0]
  elif format == 'comparison':
    if len(names) != 2:
      raise Exception(f'Expected 2 names for comparison format. Got: {len(names)}')
    if len(labels) != 2:
      raise Exception(f'Expected 2 names for comparison format. Got: {len(labels)}')
    if len(ref_seqs) != 2:
      raise Exception(f'Expected 2 reference sequences for comparison format. Got: {len(ref_seqs)}')
    data_info['name_1'] = names[0]
    data_info['name_2'] = names[1]
    data_info['name'] = names[0] + '_' + names[1]
    data_info['label_1'] = labels[0]
    data_info['label_2'] = labels[1]
    data_info['ref_seq_1'] = ref_seqs[0]
    data_info['ref_seq_2'] = ref_seqs[1]
  else:
    raise Exception('Unknown data format: ' + str(format))
  return data_info

def sort_by_count(data, count_cols, other_cols):
  count_max = data[count_cols].max(axis='columns')
  data['__temp__'] = -count_max # Sort descending
  data = data.sort_values(['__temp__'] + other_cols)
  data = data.drop(columns=['__temp__'])
  data = data.reset_index(drop=True)
  return data
