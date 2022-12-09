import os
import library_constants

IMAGE_DIR = 'images'

def make_file_name(dir, *args, ext):
  return os.path.join(dir, '_'.join(map(str, args)) + os.path.extsep + ext)

def window(dir, freq_type, subst_type):
  library_constants.check_subst_type(subst_type)
  library_constants.check_freq_type(freq_type)
  return make_file_name(dir, 'window', freq_type, subst_type, ext = 'tsv')

def sequence_data(dir, subst_type):
  library_constants.check_subst_type(subst_type)
  return make_file_name(dir, 'sequence_data', subst_type, ext = 'tsv')

def edge_data(dir, subst_type):
  library_constants.check_subst_type(subst_type)
  return make_file_name(dir, 'edge_data', subst_type, ext = 'tsv')

def distance_matrix(dir, subst_type):
  library_constants.check_subst_type(subst_type)
  return make_file_name(dir, 'distance_matrix', subst_type, ext = 'tsv')

def graph_stats(dir, subst_type):
  library_constants.check_subst_type(subst_type)
  return make_file_name(dir, 'graph_stats', subst_type, ext = 'tsv')

def variation(dir, subst_type):
  library_constants.check_subst_type(subst_type)
  return make_file_name(dir, 'variation', subst_type, ext = 'tsv')

def variation_grouped(dir, subst_type):
  library_constants.check_subst_type(subst_type)
  return make_file_name(dir, 'variation_grouped', subst_type, ext = 'tsv')

def data_info(dir):
  return make_file_name(dir, 'data_info', ext = 'tsv')

def ref(dir):
  return make_file_name(dir, 'ref', ext='.fa')

def histogram_3d(data_name, variation_type):
  return '_'.join([data_name, variation_type]) + os.path.extsep + 'png'

def graph_figure(data_name, ext='png'):
  return data_name + os.path.extsep  + ext
