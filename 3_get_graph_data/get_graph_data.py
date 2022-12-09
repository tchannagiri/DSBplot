import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

import itertools
import pandas as pd
import argparse
import shutil

import file_names
import file_utils
import alignment_utils
import log_utils
import common_utils
import library_constants
import graph_utils

def parse_args():
  parser = argparse.ArgumentParser(
    description = 'Process data for the graph analysis.'
  )
  parser.add_argument(
    '--input',
    type = common_utils.check_dir,
    help = 'Directory where output from "windows" stage is located.',
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
    choices = library_constants.SUBST_TYPES,
    help = 'Whether to process the files with/without substitutions.',
  )
  return parser.parse_args()

def get_sequence_data(data, data_format):
  """
    Get the data for the vertices of the graph.
  """
  dist_ref = []
  variation_type = []
  substitution = []
  insertion = []
  deletion = []
  indel = []

  for row in data.to_dict('records'):
    num_ins, num_del, num_subst = (
      alignment_utils.count_variations(row['ref_align'], row['read_align'])
    )
    if num_ins + num_del + num_subst == 0:
      var_type = 'none'
    elif num_del + num_subst == 0:
      var_type = 'insertion'
    elif num_subst + num_ins == 0:
      var_type = 'deletion'
    elif num_ins + num_del == 0:
      var_type = 'substitution'
    else:
      var_type = 'mixed'
    variation_type.append(var_type)
    dist_ref.append(num_ins + num_del + num_subst)
    substitution.append(num_subst)
    insertion.append(num_ins)
    deletion.append(num_del)
    indel.append(num_ins + num_del)

  all_data = data.to_dict('list')
  all_data.update({
    'is_ref': [x == 0 for x in dist_ref],
    'dist_ref': dist_ref,
    'variation_type': variation_type,
    'substitution': substitution,
    'insertion': insertion,
    'deletion': deletion,
    'indel': indel,
  })

  all_data = pd.DataFrame(all_data)
  all_data['freq_max'] = all_data[library_constants.FREQ_COLUMNS[data_format]].max(axis='columns')
  all_data = all_data.sort_values('freq_max', ascending=False)
  all_data = all_data.drop('freq_max', axis='columns')
  all_data['id'] = 'S' + pd.Series(range(1, all_data.shape[0] + 1), dtype=str)
  all_data = all_data[['id'] + list(all_data.columns[all_data.columns != 'id'])]

  all_data = pd.concat(
    [
      all_data,
      common_utils.get_freq_ranks(
        all_data,
        library_constants.FREQ_COLUMNS[data_format],
        library_constants.FREQ_RANK_COLUMNS[data_format],
      )
    ],
    axis = 'columns',
  )

  return pd.DataFrame(all_data)

def write_sequence_data(input_dir, output_dir, subst_type):
  """
    Make the main node data and write it to a file.
  """

  data = file_utils.read_tsv(
    file_names.window(
      input_dir,
      library_constants.FREQ_FILTER_MEAN,
      subst_type,
    )
  )
  data_info = file_utils.read_tsv_dict(file_names.data_info(output_dir))
  data = get_sequence_data(data, data_info['format'])
  out_file_name = file_names.sequence_data(output_dir, subst_type)
  file_utils.write_tsv(data, out_file_name)
  log_utils.log(out_file_name)

def get_edge_data(sequence_data):
  """
    Make adjacency edge data from sequence data.
  """

  edges = {
    'id_a': [],
    'id_b': [],
    'ref_align_a': [],
    'read_align_a': [],
    'variation_type_a': [],
    'ref_align_b': [],
    'read_align_b': [],
    'variation_type_b': [],
    'edge_type': [],
  }
  for row_a, row_b in itertools.combinations(sequence_data.to_dict('records'), 2):
    ref_align_a = row_a['ref_align']
    read_align_a = row_a['read_align']
    ref_align_b = row_b['ref_align']
    read_align_b = row_b['read_align']
    if graph_utils.is_alignment_adjacent_2(
      read_align_a,
      read_align_b,
    ):
      if (
        (row_a['insertion'] != row_b['insertion']) or
        (row_a['deletion'] != row_b['deletion'])
      ):
        edge_type = 'indel'
      else:
        edge_type = 'substitution'
      
      edges['id_a'].append(row_a['id'])
      edges['id_b'].append(row_b['id'])

      edges['ref_align_a'].append(ref_align_a)
      edges['read_align_a'].append(read_align_a)
      edges['variation_type_a'].append(row_a['variation_type'])

      edges['ref_align_b'].append(ref_align_b)
      edges['read_align_b'].append(read_align_b)
      edges['variation_type_b'].append(row_b['variation_type'])

      edges['edge_type'].append(edge_type)
  return pd.DataFrame(edges)

def write_edge_data(output_dir, subst_type):
  """
    Make adjacency edge data and write to file.
    Sequence data should have been created already.
  """
  in_file_name = file_names.sequence_data(output_dir, subst_type)
  out_file_name = file_names.edge_data(output_dir, subst_type)

  sequence_data = file_utils.read_tsv(in_file_name)
  edge_data = get_edge_data(sequence_data)
  file_utils.write_tsv(edge_data, out_file_name)
  log_utils.log(out_file_name)

def get_distance_matrix(sequence_data):
  """
    Get pairwise distances between vertices.
  """
  distance_matrix = {
    'id_a': [],
    'id_b': [],
    'dist': [],
  }
  for row_a, row_b in itertools.combinations(sequence_data.to_dict('records'), 2):
    read_align_a = row_a['read_align']
    read_align_b = row_b['read_align']
    distance_matrix['id_a'].append(row_a['id'])
    distance_matrix['id_b'].append(row_b['id'])
    distance_matrix['dist'].append(
      graph_utils.get_alignment_distance_2(
        read_align_a,
        read_align_b,
      )
    )
    
  return pd.DataFrame(distance_matrix)

def write_distance_matrix(output_dir, subst_type):
  """
    Get distance matrix and write to file.
    Sequence data should have been created already.
  """

  in_file_name = file_names.sequence_data(output_dir, subst_type)
  out_file_name = file_names.distance_matrix(output_dir, subst_type)

  sequence_data = file_utils.read_tsv(in_file_name)
  distance_matrix = get_distance_matrix(sequence_data)
  file_utils.write_tsv(distance_matrix, out_file_name)
  log_utils.log(out_file_name)

def write_graph_stats(output_dir, subst_type):
  """
    Get graph summary statistics and write to file.
    Sequence data and edge data should have been created already.
  """

  graph = graph_utils.load_graph(output_dir, subst_type)
  data_info = file_utils.read_tsv_dict(file_names.data_info(output_dir))
  graph_stats = graph_utils.get_graph_stats_ref_component(data_info['format'], graph)
  graph_stats = pd.DataFrame.from_records([graph_stats])
  out_file_name = file_names.graph_stats(output_dir, subst_type)
  file_utils.write_tsv(graph_stats, out_file_name)
  log_utils.log(out_file_name)

def main():
  args = parse_args()

  log_utils.log(args.input)
  log_utils.log('------>')

  # copy graphs stats
  input_data_info_file = file_names.data_info(args.input)
  output_data_info_file = file_names.data_info(args.output)
  shutil.copy(input_data_info_file, output_data_info_file)
  log_utils.log(output_data_info_file)

  write_sequence_data(args.input, args.output, args.subst_type)
  write_edge_data(args.output, args.subst_type)
  write_distance_matrix(args.output, args.subst_type)
  write_graph_stats(args.output, args.subst_type)

  log_utils.new_line()

if __name__ == '__main__':
  main()
