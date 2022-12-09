import file_names
import file_utils
import library_constants
import alignment_utils

import Levenshtein
import networkx as nx
import numpy as np

def get_extended_align(ref_align, read_align):
  ord_map = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3,
    'N': 4,
    '-': 5,
  }
  extended_alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
  extended_align = [
    extended_alphabet[len(ord_map) * ord_map[ref_align[i]] + ord_map[read_align[i]]]
    for i in range(len(ref_align))
  ]
  return ''.join(extended_align)


def get_alignment_distance(ref_align_1, read_align_1, ref_align_2, read_align_2):
  extended_align_1 = get_extended_align(ref_align_1, read_align_1)
  extended_align_2 = get_extended_align(ref_align_2, read_align_2)
  return Levenshtein.distance(extended_align_1, extended_align_2)

def is_alignment_adjacent(ref_align_1, read_align_1, ref_align_2, read_align_2):
  return get_alignment_distance(ref_align_1, read_align_1, ref_align_2, read_align_2) == 1

def get_alignment_distance_2(read_align_1, read_align_2):
  return Levenshtein.distance(
    alignment_utils.get_orig_seq(read_align_1),
    alignment_utils.get_orig_seq(read_align_2),
  )

def is_alignment_adjacent_2(read_align_1, read_align_2):
  return get_alignment_distance_2(read_align_1, read_align_2) == 1

def load_graph(dir, subst_type):
  node_data = file_utils.read_tsv(file_names.sequence_data(dir, subst_type))
  node_data = node_data.set_index('id', drop=False)
  edge_data = file_utils.read_tsv(file_names.edge_data(dir, subst_type))
  graph = nx.Graph()
  graph.add_nodes_from(node_data['id'])
  graph.add_edges_from(zip(
    edge_data['id_a'],
    edge_data['id_b'],
    edge_data.to_dict('records')),
  )
  nx.set_node_attributes(graph, node_data.to_dict('index'))
  return graph

# graph nodes must have the associated variation types!
def get_graph_stats_ref_component(data_format, graph):
  ref_id = next(
    (x[0] for x in graph.nodes(data=True) if x[1]['is_ref']),
    None,
  )
  if ref_id is None:
    graph = nx.Graph()
  else:
    graph = graph.subgraph(nx.node_connected_component(graph, ref_id))

  num_nodes = len(graph.nodes())
  num_edges = len(graph.edges())
  num_edges_indel = sum(
    1 for x in graph.edges(data=True) if x[2]['edge_type'] == 'indel'
  )
  num_edges_substitution = sum(
    1 for x in graph.edges(data=True) if x[2]['edge_type'] == 'substitution'
  )
  if num_nodes == 0:
    avg_degree = None
  else:
    avg_degree = np.mean([x[1] for x in graph.degree()])

  shorted_path_lengths = dict(nx.all_pairs_shortest_path_length(graph))

  if ref_id not in shorted_path_lengths:
    avg_dist_ref = None
    max_dist_ref = None
  else:
    shorted_path_lengths_ref = shorted_path_lengths.get(ref_id, [])
    if num_nodes == 1:
      avg_dist_ref = 0
      max_dist_ref = 0
    else:
      total_dist_ref = 0
      max_dist_ref = 0
      for length in shorted_path_lengths_ref.values():
        total_dist_ref += length
        max_dist_ref = max(max_dist_ref, length)
      avg_dist_ref = total_dist_ref / (num_nodes - 1)

  if num_nodes <= 1:
    max_pairwise_dist = None
    avg_pairwise_dist = None
  else:
    max_pairwise_dist = 0
    total_pairwise_dist = 0
    for id_a, shorted_path_lengths_a in shorted_path_lengths.items():
      for id_b, length in shorted_path_lengths_a.items():
        if id_a < id_b:
          total_pairwise_dist += length
          max_pairwise_dist = max(max_pairwise_dist, length)
    avg_pairwise_dist = 2 * total_pairwise_dist / num_nodes / (num_nodes - 1)

  node_view = graph.nodes(data=True)

  num_seq_substitution = sum(
    1 for x in node_view if x[1]['variation_type'] == 'substitution'
  )

  num_seq_insertion = sum(
    1 for x in node_view if x[1]['variation_type'] == 'insertion'
  )

  num_seq_deletion = sum(
    1 for x in node_view if x[1]['variation_type'] == 'deletion'
  )

  freq_stats = {}
  if graph.number_of_nodes() > 0:
    all_columns = list(node_view[ref_id])
    for column in library_constants.FREQ_COLUMNS[data_format]:
      ref_freq = node_view[ref_id][column] if num_nodes > 0 else 0
      non_ref_freq = sum(x[1][column] for x in node_view if x[0] != ref_id)
      freq_stats['ref_' + column] = ref_freq
      freq_stats['non_ref_' + column] = non_ref_freq
      for var_type in ['insertion', 'deletion']:
        freq_stats[var_type + '_' + column] = sum(
          x[1][column] for x in node_view
          if x[1]['variation_type'] == var_type
        )

  return dict(
    num_nodes = num_nodes,
    num_edges = num_edges,
    num_edges_indel = num_edges_indel,
    num_edges_substitution = num_edges_substitution,
    avg_degree = avg_degree,
    avg_dist_ref = avg_dist_ref,
    max_dist_ref = max_dist_ref,
    avg_pairwise_dist = avg_pairwise_dist,
    max_pairwise_dist = max_pairwise_dist,
    num_seq_substitution = num_seq_substitution,
    num_seq_insertion = num_seq_insertion,
    num_seq_deletion = num_seq_deletion,
    **freq_stats,
  )
