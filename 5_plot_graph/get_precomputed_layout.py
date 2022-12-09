import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

import argparse

import pandas as pd
import numpy as np
import networkx as nx
import plotly.graph_objects as go
import plotly.subplots as ps

import library_constants
import common_utils
import log_utils
import kmer_utils
import file_utils
import file_names
import plot_graph

def parse_args():
  parser = argparse.ArgumentParser(
    description = 'Make precomputed layouts for graphs.'
  )
  parser.add_argument(
    '--input',
    type = common_utils.check_dir,
    nargs = '+',
    help = (
      'List of data directories created with make_graph_data.py.\n'
      'All data sets should have the same window reference sequence and\n'
      ' be individual (not comparison) experiments.'
    ),
    required = True,
  )
  parser.add_argument(
    '--reverse_complement',
    choices = ['0', '1'],
    nargs = '+',
    help = (
      'Whether to reverse complement the sequence in the data sets.' +
      ' If present, the number of values must be the same as the number of input directories.' +
      ' "1" mean reverse complement the sequence and "0" means do not.'
      ' Used for making a layout for data sets that have reference sequences'
      ' that are the reverse complements of each other.'
    ),
  )
  parser.add_argument(
    '--output',
    type = common_utils.check_dir_output,
    help = 'Output directory.',
    required = True,
  )
  parser.add_argument(
    '--subst_type',
    choices = library_constants.SUBST_TYPES,
    help = 'Whether to use the data with or without substitutions.',
  )
  parser.add_argument(
    '--layout',
    default = 'radial',
    choices = ['radial', 'mds', 'kamada', 'universal'],
    help = 'Type of layout to use.',
  )
  args = parser.parse_args()
  args.layout += '_layout'
  if args.reverse_complement is None:
    args.reverse_complement = ['0'] * len(args.input)
  if len(args.reverse_complement) != len(args.input):
    raise Exception(
      'Incorrect number of reverse complement flags: '
      f'{len(args.reverse_complement)}. Expected {len(args.input)}.'
    )
  args.reverse_complement = [x == '1' for x in args.reverse_complement]
  return args

def get_precomputed_layout(
  precomputed_layout_dir,
  node_data,
  node_subst_type,
  reverse_complement = False,
):
  """
    Get the precomputed layout coordinates for the node data.
    Assumes that the sequences in node_data are part of the
    precomputed layout in precomputed_layout_dir.If reverse complement is true
    then sequences in node_data will be
    reverse complemented before joining with the common layout.
  """
  layout = file_utils.read_tsv(
    file_names.sequence_data(
      precomputed_layout_dir,
      node_subst_type,
    )
  )

  node_data = node_data.reset_index(drop=True)

  if reverse_complement:
    node_data = node_data.assign(
      ref_align = node_data['ref_align'].apply(kmer_utils.reverse_complement),
      read_align = node_data['read_align'].apply(kmer_utils.reverse_complement),
    )
  
  node_data = pd.merge(
    node_data[['id', 'ref_align', 'read_align']],
    layout[['ref_align', 'read_align', 'x', 'y']],
    on = ['ref_align', 'read_align'],
    how = 'inner',
  )[['id', 'x', 'y']]
  node_data = node_data.set_index('id', drop=True).rename(
    {'x': 0, 'y': 1},
    axis = 'columns',
  ) 
  return node_data


def make_precomputed_layout(
  data_dir_list,
  reverse_complement_list,
  output_dir,
  subst_type,
  layout_type,
):
  """
    Make common layout files by combining the node information in the input data sets.
  """

  for data_dir in data_dir_list:
    log_utils.log(data_dir)
  log_utils.log('------>')

  ### Load node data and edge data ###
  seq_data_list = [
    file_utils.read_tsv(
      file_names.sequence_data(data_dir, subst_type)
    )
    for data_dir in data_dir_list
  ]
  edge_data_list = [
    file_utils.read_tsv(
      file_names.edge_data(data_dir, subst_type),
    )
    for data_dir in data_dir_list
  ]

  ### Reverse complement the reverse strand sequences ###
  for i in range(len(data_dir_list)):
    if reverse_complement_list[i]:
      seq_data_list[i] = seq_data_list[i].assign(
        ref_align = seq_data_list[i]['ref_align'].apply(kmer_utils.reverse_complement),
        read_align = seq_data_list[i]['read_align'].apply(kmer_utils.reverse_complement),
      )

  ### Remove id's from edge data ###
  edge_data_list = [
    edge_data[edge_data.columns[~edge_data.columns.isin(['id_a', 'id_b'])]]
    for edge_data in edge_data_list
  ]

  ### Combine sequence data ###
  seq_data = pd.concat(seq_data_list, axis='index', ignore_index=True)
  seq_data = seq_data.groupby(['ref_align', 'read_align'])
  freq_mean_max = seq_data['freq_mean'].max().reset_index(drop=True)
  seq_data = seq_data.first().reset_index()
  seq_data['freq_mean_max'] = freq_mean_max
  seq_data = seq_data.sort_values('freq_mean_max', ascending=False).reset_index(drop=True)
  seq_data['id'] = 'S' + pd.Series(range(1, seq_data.shape[0] + 1), dtype=str)
  seq_data = seq_data.set_index('id', drop=False)

  ### Combine edge data ###
  edge_data = pd.concat(edge_data_list, axis='index', ignore_index=True)
  edge_data = edge_data.groupby(list(edge_data.columns)).first().reset_index()

  ### Get the new id's for the edges ###
  for suffix in ['_a', '_b']:
    edge_data = pd.merge(
      edge_data,
      seq_data[['id', 'ref_align', 'read_align']].rename(
        {
          'id': 'id' + suffix,
          'ref_align': 'ref_align' + suffix,
          'read_align': 'read_align' + suffix,
        },
        axis = 'columns',
      ),
      on = ['ref_align' + suffix, 'read_align' + suffix],
      how = 'inner',
    )

  ### Make the combined graph ###
  graph = nx.Graph()
  graph.add_nodes_from(seq_data.index)
  nx.set_node_attributes(graph, seq_data.to_dict('index'))

  graph.add_edges_from(zip(
    edge_data['id_a'],
    edge_data['id_b'],
    edge_data.to_dict('records')),
  )

  ### Make the layout ###
  data_info = file_utils.read_tsv_dict(
    file_names.data_info(data_dir_list[0])
  )
  layout = plot_graph.make_graph_layout(
    data_dir = data_dir,
    data_info = data_info,
    node_type = 'sequence_data',
    node_subst_type = subst_type,
    graph = graph,
    layout_type = layout_type,
    node_size_px_dict = None,
    x_size_domain = None,
    y_size_domain = None,
    x_size_px = None,
    y_size_px = None,
    separate_components = False,
    precomputed_layout_dir = None,
  )
  layout.columns = ['x', 'y']

  ### Join layout with sequence data ###
  seq_data = seq_data.join(layout)
  seq_data = seq_data.reset_index(drop=True)

  ### Write to files ###
  file_out = file_names.sequence_data(output_dir, subst_type)
  log_utils.log(file_out)
  file_utils.write_tsv(seq_data, file_out)
  file_out = file_names.edge_data(output_dir, subst_type)
  file_utils.write_tsv(edge_data, file_out)
  log_utils.log(file_out)
  log_utils.new_line()

def main():
  args = parse_args()
  make_precomputed_layout(
    args.input,
    args.reverse_complement,
    args.output,
    args.subst_type,
    args.layout,
  ) 

if __name__ == '__main__':
  main()