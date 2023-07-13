import os

import DSBplot.utils.constants as constants

IMAGE_DIR = 'images'

def make_file_name(dir, *args, ext):
  return os.path.join(dir, '_'.join(map(str, args)) + '.' + ext)

def bowtie2_build_dir(dir, name):
  return os.path.join(dir, 'bowtie2_build', name)

def sam_file(dir, name):
  return make_file_name(dir, name, ext='sam')

def ref_seq_file(dir):
  return os.path.join(dir, 'ref_seq.fasta')

def window(dir, subst_type):
  return make_file_name(dir, 'window', subst_type, ext='csv')

# FIXME: DELETE
def sequence_data(dir, subst_type):
  constants.check_subst_type(subst_type)
  return make_file_name(dir, 'sequence_data', subst_type, ext='tsv')

# FIXME: DELETE
def edge_data(dir, subst_type):
  constants.check_subst_type(subst_type)
  return make_file_name(dir, 'edge_data', subst_type, ext='tsv')

# FIXME: DELETE
def graph_stats(dir, subst_type):
  constants.check_subst_type(subst_type)
  return make_file_name(dir, 'graph_stats', subst_type, ext='tsv')

def variation(dir, subst_type):
  constants.check_subst_type(subst_type)
  return make_file_name(dir, 'variation', subst_type, ext='csv')

def data_info(dir):
  return make_file_name(dir, 'data_info', ext='csv')

def bowtie2_build(dir):
  return os.path.join(dir, 'bowtie2_build', 'build')

def filter_nhej(dir, suffix):
  return os.path.join(dir, 'filter_nhej_' + suffix + '.csv')

# FIXME: DELETE
def graph_dir(dir):
  return os.path.join(dir, '4_graph')

# FIXME: DELETE
def histogram_dir(dir):
  return os.path.join(dir, '5_histogram')

def get_file_name(full_path):
  return os.path.basename(full_path).split('.')[0]