import Levenshtein
import networkx as nx
import pandas as pd
import itertools

import DSBplot.utils.alignment_utils as alignment_utils

# FIXME: DELETE
# def get_extended_align(ref_align, read_align):
#   ord_map = {
#     'A': 0,
#     'C': 1,
#     'G': 2,
#     'T': 3,
#     'N': 4,
#     '-': 5,
#   }
#   extended_alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
#   extended_align = [
#     extended_alphabet[len(ord_map) * ord_map[ref_align[i]] + ord_map[read_align[i]]]
#     for i in range(len(ref_align))
#   ]
#   return ''.join(extended_align)


# FIXME: DELETE
# def get_alignment_distance(ref_align_1, read_align_1, ref_align_2, read_align_2):
#   extended_align_1 = get_extended_align(ref_align_1, read_align_1)
#   extended_align_2 = get_extended_align(ref_align_2, read_align_2)
#   return Levenshtein.distance(extended_align_1, extended_align_2)

# FIXME: DELETE
# def is_alignment_adjacent(ref_align_1, read_align_1, ref_align_2, read_align_2):
#   return get_alignment_distance(ref_align_1, read_align_1, ref_align_2, read_align_2) == 1

def get_alignment_distance(read_align_1, read_align_2):
  return Levenshtein.distance(
    alignment_utils.get_orig_seq(read_align_1),
    alignment_utils.get_orig_seq(read_align_2),
  )

def is_alignment_adjacent(read_align_1, read_align_2):
  return get_alignment_distance(read_align_1, read_align_2) == 1

# FIXME: DELETE
# def load_graph(dir, subst_type):
#   node_data = file_utils.read_tsv(file_names.sequence_data(dir, subst_type))
#   node_data = node_data.set_index('id', drop=False)
#   edge_data = file_utils.read_tsv(file_names.edge_data(dir, subst_type))
#   graph = nx.Graph()
#   graph.add_nodes_from(node_data['id'])
#   graph.add_edges_from(zip(
#     edge_data['id_a'],
#     edge_data['id_b'],
#     edge_data.to_dict('records')),
#   )
#   nx.set_node_attributes(graph, node_data.to_dict('index'))
#   return graph

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

# FIXME: DELETE
# # graph nodes must have the associated variation types!
# def get_graph_stats_ref_component(data_format, graph):
#   ref_id = next(
#     (x[0] for x in graph.nodes(data=True) if x[1]['is_ref']),
#     None,
#   )
#   if ref_id is None:
#     graph = nx.Graph()
#   else:
#     graph = graph.subgraph(nx.node_connected_component(graph, ref_id))

#   num_nodes = len(graph.nodes())
#   num_edges = len(graph.edges())
#   num_edges_indel = sum(
#     1 for x in graph.edges(data=True) if x[2]['edge_type'] == 'indel'
#   )
#   num_edges_substitution = sum(
#     1 for x in graph.edges(data=True) if x[2]['edge_type'] == 'substitution'
#   )
#   if num_nodes == 0:
#     avg_degree = None
#   else:
#     avg_degree = np.mean([x[1] for x in graph.degree()])

#   shorted_path_lengths = dict(nx.all_pairs_shortest_path_length(graph))

#   if ref_id not in shorted_path_lengths:
#     avg_dist_ref = None
#     max_dist_ref = None
#   else:
#     shorted_path_lengths_ref = shorted_path_lengths.get(ref_id, [])
#     if num_nodes == 1:
#       avg_dist_ref = 0
#       max_dist_ref = 0
#     else:
#       total_dist_ref = 0
#       max_dist_ref = 0
#       for length in shorted_path_lengths_ref.values():
#         total_dist_ref += length
#         max_dist_ref = max(max_dist_ref, length)
#       avg_dist_ref = total_dist_ref / (num_nodes - 1)

#   if num_nodes <= 1:
#     max_pairwise_dist = None
#     avg_pairwise_dist = None
#   else:
#     max_pairwise_dist = 0
#     total_pairwise_dist = 0
#     for id_a, shorted_path_lengths_a in shorted_path_lengths.items():
#       for id_b, length in shorted_path_lengths_a.items():
#         if id_a < id_b:
#           total_pairwise_dist += length
#           max_pairwise_dist = max(max_pairwise_dist, length)
#     avg_pairwise_dist = 2 * total_pairwise_dist / num_nodes / (num_nodes - 1)

#   node_view = graph.nodes(data=True)

#   num_seq_substitution = sum(
#     1 for x in node_view if x[1]['variation_type'] == 'substitution'
#   )

#   num_seq_insertion = sum(
#     1 for x in node_view if x[1]['variation_type'] == 'insertion'
#   )

#   num_seq_deletion = sum(
#     1 for x in node_view if x[1]['variation_type'] == 'deletion'
#   )

#   freq_stats = {}
#   if graph.number_of_nodes() > 0:
#     if data_format == 'comparison':
#       freq_columns = ['freq_mean_1', 'freq_mean_2']
#     elif data_format == 'individual':
#       freq_columns = ['freq_mean']
#     else:
#       raise Exception('Unknown data format: {}'.format(data_format))
#     for column in freq_columns:
#       ref_freq = node_view[ref_id][column] if num_nodes > 0 else 0
#       non_ref_freq = sum(x[1][column] for x in node_view if x[0] != ref_id)
#       freq_stats['ref_' + column] = ref_freq
#       freq_stats['non_ref_' + column] = non_ref_freq
#       for var_type in ['insertion', 'deletion']:
#         freq_stats[var_type + '_' + column] = sum(
#           x[1][column] for x in node_view
#           if x[1]['variation_type'] == var_type
#         )

#   return dict(
#     num_nodes = num_nodes,
#     num_edges = num_edges,
#     num_edges_indel = num_edges_indel,
#     num_edges_substitution = num_edges_substitution,
#     avg_degree = avg_degree,
#     avg_dist_ref = avg_dist_ref,
#     max_dist_ref = max_dist_ref,
#     avg_pairwise_dist = avg_pairwise_dist,
#     max_pairwise_dist = max_pairwise_dist,
#     num_seq_substitution = num_seq_substitution,
#     num_seq_insertion = num_seq_insertion,
#     num_seq_deletion = num_seq_deletion,
#     **freq_stats,
#   )

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
  data_info = {
    'format': 'comparison',
    'label_1': data_info_1['label'],
    'label_2': data_info_2['label'],
    'ref_seq_1': data_info_1['ref_seq'],
    'ref_seq_2': data_info_2['ref_seq'],
    'ref_seq_window': data_info_1['ref_seq_window'],
  }

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
