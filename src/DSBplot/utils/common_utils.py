import os

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
  os.makedirs(dir_name, exist_ok=True)
  return check_dir(dir_name)

def check_file_output(file_name):
  dir_name = os.path.dirname(file_name)
  if dir_name != '':
    os.makedirs(dir_name, exist_ok=True)
  return file_name

def join_with_comma(arr):
  return ','.join([str(x) for x in arr])
