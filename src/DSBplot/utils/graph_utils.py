import Levenshtein
import networkx as nx
import pandas as pd
import itertools

import DSBplot.utils.common_utils as common_utils
import DSBplot.utils.alignment_utils as alignment_utils

def get_alignment_distance(read_align_1, read_align_2):
  return Levenshtein.distance(
    alignment_utils.get_orig_seq(read_align_1),
    alignment_utils.get_orig_seq(read_align_2),
  )

def is_alignment_adjacent(read_align_1, read_align_2):
  return get_alignment_distance(read_align_1, read_align_2) == 1

def get_graph(node_data, edge_data):
  graph = nx.Graph()
  graph.add_nodes_from(node_data['id'])
  graph.add_edges_from(zip(
    edge_data['id_a'],
    edge_data['id_b'],
    edge_data.to_dict('records'),
  ))
  nx.set_node_attributes(
    graph,
    node_data.to_dict('index'), # index should be id
  )
  return graph

def get_node_data(data):
  """
    Get the data for the vertices of the graph.
  """
  dist_ref = []
  variation_type = []
  substitution = []
  insertion = []
  deletion = []
  indel = []

  data = data.sort_values(
    ['freq_mean', 'ref_align', 'read_align'],
    ascending = [False, True, True],
  )

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

  data = data[
    ['ref_align', 'read_align'] +
    [x for x in data.columns if x in ['freq_mean', 'freq_mean_1', 'freq_mean_2']]
  ].to_dict('list')
  data.update({
    'is_ref': [x == 0 for x in dist_ref],
    'dist_ref': dist_ref,
    'variation_type': variation_type,
    'substitution': substitution,
    'insertion': insertion,
    'deletion': deletion,
    'indel': indel,
  })

  data = pd.DataFrame(data)
  data = data.sort_values('freq_mean', ascending=False)
  data['id'] = 'S' + pd.Series(range(1, data.shape[0] + 1), dtype=str)
  data = data[['id'] + list(data.columns[data.columns != 'id'])]
  data = data.set_index('id', drop=False)

  return pd.DataFrame(data)

def get_edge_data(node_data):
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
  for row_a, row_b in itertools.combinations(node_data.to_dict('records'), 2):
    ref_align_a = row_a['ref_align']
    read_align_a = row_a['read_align']
    ref_align_b = row_b['ref_align']
    read_align_b = row_b['read_align']
    if is_alignment_adjacent(
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

# Get data combining two experiments for comparison
def get_comparison_data(data_info_1, data_info_2, data_1, data_2):
  # Make sure the experiments are compatible
  if data_info_1['ref_seq_window'] != data_info_2['ref_seq_window']:
    raise Exception(f'Incompatible experiments:\n{data_info_1}\n{data_info_2}')

  # Make the comparison info
  data_info = common_utils.make_data_info(
    format = 'comparison',
    names = [data_info_1['name'], data_info_2['name']],
    labels = [data_info_1['label'], data_info_2['label']],
    ref_seqs = [data_info_1['ref_seq'], data_info_2['ref_seq']],
    ref_seq_window = data_info_1['ref_seq_window'],
  )

  # Make the comparison data
  data_1 = data_1[['ref_align', 'read_align', 'freq_mean']]
  data_2 = data_2[['ref_align', 'read_align', 'freq_mean']]

  data = pd.merge(
    data_1,
    data_2,
    how = 'outer',
    on = ['ref_align', 'read_align'],
    suffixes = ['_1', '_2'],
  )

  data[['freq_mean_1', 'freq_mean_2']] = (
    data[['freq_mean_1', 'freq_mean_2']].fillna(0)
  )
  data['freq_mean'] = (
    data[['freq_mean_1', 'freq_mean_2']].max(axis='columns')
  )
  return data_info, data
