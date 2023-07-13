import os

import DSBplot.utils.constants as constants

IMAGE_DIR = 'images'

def make_file_name(dir, *args, ext):
  return os.path.join(dir, '_'.join(map(str, args)) + '.' + ext)

def bowtie2_build_dir(dir, name):
  return os.path.join(dir, 'bowtie2_build', name)

# def window(dir, freq_type, subst_type):
#   constants.check_subst_type(subst_type)
#   constants.check_freq_type(freq_type)
#   return make_file_name(dir, 'window', freq_type, subst_type, ext = 'tsv')
def window(dir, subst_type):
  return make_file_name(dir, 'window', subst_type, ext='tsv')

def sequence_data(dir, subst_type):
  constants.check_subst_type(subst_type)
  return make_file_name(dir, 'sequence_data', subst_type, ext='tsv')

def edge_data(dir, subst_type):
  constants.check_subst_type(subst_type)
  return make_file_name(dir, 'edge_data', subst_type, ext='tsv')

def graph_stats(dir, subst_type):
  constants.check_subst_type(subst_type)
  return make_file_name(dir, 'graph_stats', subst_type, ext='tsv')

def variation(dir, subst_type):
  constants.check_subst_type(subst_type)
  return make_file_name(dir, 'variation', subst_type, ext='tsv')

def data_info(dir):
  return make_file_name(dir, 'data_info', ext='tsv')

def bowtie2_build(dir):
  return os.path.join(dir, '0_bowtie2_build', 'build')

def sam_dir(dir):
  return os.path.join(dir, '0_sam')

def filter_nhej(dir, suffix):
  return os.path.join(dir, 'filter_nhej_' + suffix + '.csv')

def combine_repeat_file(dir):
  return os.path.join(dir, '2_combine_repeat', 'out.tsv')

def window_dir(dir):
  return os.path.join(dir, '3_window')

def graph_dir(dir):
  return os.path.join(dir, '4_graph')

def histogram_dir(dir):
  return os.path.join(dir, '5_histogram')
