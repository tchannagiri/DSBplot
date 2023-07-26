import os

import DSBplot.utils.constants as constants

IMAGE_DIR = 'images'

def make_file_name(dir, *args, ext):
  return os.path.join(dir, '_'.join(map(str, args)) + '.' + ext)

def sam_file(dir, name):
  return make_file_name(dir, name, ext='sam')

def ref_seq_file(dir):
  return os.path.join(dir, 'ref_seq.fasta')

def window(dir, subst_type):
  return make_file_name(dir, 'window', subst_type, ext='csv')

def variation(dir, subst_type):
  constants.check_subst_type(subst_type)
  return make_file_name(dir, 'variation', subst_type, ext='csv')

def data_info(dir):
  return make_file_name(dir, 'data_info', ext='json')

def bowtie2_index(dir):
  return os.path.join(dir, 'bowtie2', 'index')

def filter(dir, suffix):
  return os.path.join(dir, 'filter_' + suffix + '.csv')

def get_file_name(full_path):
  return os.path.basename(full_path).split('.')[0]