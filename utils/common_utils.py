import pandas as pd
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
  os.makedirs(os.path.dirname(file_name), exist_ok=True)
  return file_name

def check_comma_separated_values(arg_str, type_func=str):
  arg_str = arg_str.split(',')
  return [type_func(x) for x in arg_str]

def check_comma_separated_floats(arg_str):
  return check_comma_separated_values(arg_str, float)

def check_comma_separated_enum(choices):
  def check(arg_str):
    args = check_comma_separated_values(arg_str)
    if not all(x in choices for x in args):
      raise ValueError('Unexpected value: ' + str(arg_str))
  return check

def join_with_comma(arr):
  return ','.join([str(x) for x in arr])
