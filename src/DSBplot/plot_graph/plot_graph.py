import os

import argparse

import networkx as nx

import plotly.graph_objects

import pandas as pd
import numpy as np

import sklearn.decomposition
import sklearn.manifold

import PIL.Image

import DSBplot.utils.constants as constants
import DSBplot.utils.common_utils as common_utils
import DSBplot.utils.log_utils as log_utils
import DSBplot.utils.graph_utils as graph_utils
import DSBplot.utils.file_utils as file_utils
import DSBplot.utils.kmer_utils as kmer_utils
import DSBplot.utils.alignment_utils as alignment_utils
import DSBplot.utils.file_names as file_names
import DSBplot.plot_graph.plot_graph_helper as plot_graph_helper

LAYOUT_PROPERTIES = {
 'mds_layout': {
    'only_2d': True,
    'do_pca': True,
    'normalize': True,
    'has_edges': True,
    'distance_matrix': True,
  },
 'radial_layout': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': True,
  },
 'universal_layout': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': True,
  },
 'fractal_layout': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': True,
    'plot_range_x': (-1, 1),
    'plot_range_y': (-1, 1),
  },
 'kamada_layout': {
    'only_2d': False,
    'do_pca': True,
    'normalize': True,
    'has_edges': True,
    'plot_range_x': (0, 1),
    'plot_range_y': (0, 1),
  },
  'spectral_layout': {
    'only_2d': False,
    'do_pca': True,
    'normalize': False,
    'has_edges': True, 
  },
  'spring_layout': {
    'only_2d': False,
    'do_pca': True,
    'normalize': False,
    'has_edges': True, 
  },
  'shell_layout': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': True, 
  },
  'spiral_layout': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': True, 
  },
  'circular_layout': {
    'only_2d': True,
    'do_pca': False,
    'has_edges': True, 
    'normalize': False,
  },
  'multipartite_layout': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': True, 
  },
}

def group_graph_nodes_by(graph, data_name):
  data = pd.Series(dict(graph.nodes(data_name)))
  return list(data.groupby(data).groups.values())

def make_radial_layout(graph):
  node_list = graph.nodes(data=True)

  bucket_dict = {
    'insertion': {},
    'deletion': {},
  }

  ref_nodes = []

  for _, data in node_list:
    if data['dist_ref'] == 0:
      ref_nodes.append(data)
    else:
      dist_ref = data['dist_ref']
      var_type = data['variation_type']
      bucket_dict[var_type].setdefault(dist_ref, [])
      bucket_dict[var_type][dist_ref].append(data)

  xy_dict = {}
  for data in ref_nodes:
    xy_dict[data['id']] = (0, 0)
  for var_type in bucket_dict:
    # Heuristics to make the nodes nicely spaced.
    # Each distance is perturbed slightly so that nodes zig-zag a
    # bit along radii instead od being in a smooth curve. The angle is
    # also offset a bit at each distance level so that the nodes are not
    # collinear. The insertions are placed above and deletion are place below.
    # A bit of offset is added to the y-coord so that the insertions and deletions
    # do not overlap as much.
    if var_type == 'insertion':
      y_sign = 1
      dist_scale = 2
      zig_zag_angle = 15
    elif var_type == 'deletion':
      y_sign = -1
      dist_scale = 1
      zig_zag_angle = 30
    else:
      raise Exception('Impossible: ' + str(var_type))
    for dist_ref in bucket_dict[var_type]:
      bucket = list(sorted(
        bucket_dict[var_type][dist_ref],
        key = lambda x: x['freq_mean'],
        reverse = True,
      ))

      delta_angle = dist_ref + zig_zag_angle * (-1)**dist_ref / dist_ref
      delta_dist = 0
      angle_list = np.linspace(((180 - delta_angle) / 180) * np.pi, (delta_angle / 180) * np.pi, len(bucket))

      for angle, data in zip(angle_list, bucket):
        xy_dict[data['id']] = (
          (dist_ref * dist_scale + delta_dist) * np.cos(angle),
          ((dist_ref * dist_scale + delta_dist) * np.sin(angle) + 2) * y_sign
        )
  return xy_dict

def get_kmer_fractal_x_y(kmer):
  X_COORD_MAP = {'A': 0, 'C': 1, 'G': 0, 'T': 1}
  Y_COORD_MAP = {'A': 1, 'C': 1, 'G': 0, 'T': 0}
  x = 0
  y = 0
  for letter in kmer:
    x = 2 * x + X_COORD_MAP[letter]
    y = 2 * y + Y_COORD_MAP[letter]

  # Transform to (0, 1)
  x = (x + 0.5) / (2 ** len(kmer))
  y = (y + 0.5) / (2 ** len(kmer))

  # Transform to (-1, 1)
  x = 2 * x - 1
  y = 2 * y - 1

  return (x, y) 

def make_fractal_layout(graph, reverse_complement=False):
  node_list = graph.nodes(data=True)

  bucket_dict = {
    'insertion': {},
    'deletion': {},
  }

  ref_nodes = []

  for _, data in node_list:
    if data['dist_ref'] == 0:
      ref_nodes.append(data)
    else:
      dist_ref = data['dist_ref']
      var_type = data['variation_type']
      bucket_dict[var_type].setdefault(dist_ref, [])
      bucket_dict[var_type][dist_ref].append(data)

  xy_dict = {}
  for data in ref_nodes:
    xy_dict[data['id']] = (0, 0)
  for var_type in bucket_dict:
    for dist_ref in bucket_dict[var_type]:
      bucket = list(sorted(
        bucket_dict[var_type][dist_ref],
        key = lambda x: x['freq_mean'],
        reverse = True,
      ))

      for data in bucket:
        ref_align = data['ref_align']
        read_align = data['read_align']

        if reverse_complement:
          ref_align = kmer_utils.reverse_complement(ref_align)
          read_align = kmer_utils.reverse_complement(read_align)

        if var_type == 'insertion':
          xy_dict[data['id']] = get_kmer_fractal_x_y(
            alignment_utils.get_insertion_str(
              ref_align,
              read_align,
            )
          )
        elif var_type == 'deletion':
          # We don't plot the deletion nodes in this layout.
          # We place them all at (0, 0) so they will hopefully be covered by the reference node.
          xy_dict[data['id']] = (0, 0)
        else:
          raise Exception('Impossible.')
  return xy_dict

# Determines how the rows are laid out in the insertion side
# of the universal layout.
def get_universal_layout_insertion_row_spec(dist_ref):
  rows = 2 ** ((dist_ref - 1) // 2)
  cols = 4 ** dist_ref // rows
  row_space = 2 / rows
  return {
    'rows': rows,
    'cols': cols,
    'row_space': row_space,
  }

def get_pos_universal_layout(
  ref_align,
  read_align,
  dist_ref,
  var_type,
  cut_pos_ref,
  x_scale_insertion = constants.GRAPH_UNIVERSAL_LAYOUT_X_SCALE_INSERTION,
  y_scale_insertion = constants.GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_INSERTION,
  x_scale_deletion = constants.GRAPH_UNIVERSAL_LAYOUT_X_SCALE_DELETION,
  y_scale_deletion = constants.GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_DELETION,
):
  if var_type == 'insertion':
    # Place the x-coordinate alphabetically so that A is left-most
    # and T is right-most. This is intended to place the insertions
    # in a tree based on the common prefix of the insertion nucleotides.
    # To prevent overlapping the nodes are placed in multiple rows for
    # higher numbers of insertions.
    kmer_index = kmer_utils.get_kmer_index(
      alignment_utils.get_insertion_str(
        ref_align,
        read_align,
      )
    )
    
    row_spec_curr = get_universal_layout_insertion_row_spec(dist_ref)
    num_rows = row_spec_curr['rows']
    num_cols = row_spec_curr['cols']
    row = kmer_index % num_rows
    col = kmer_index // num_rows
    prev_rows_offset = 0
    for dist_ref_prev in range(1, dist_ref):
      row_spec_prev = get_universal_layout_insertion_row_spec(dist_ref_prev)
      prev_rows_offset += row_spec_prev['rows'] * row_spec_prev['row_space']

    curr_row_offset = row * row_spec_curr['row_space']
    y = 2 + prev_rows_offset + curr_row_offset
    x = (col / num_cols) - 0.5 * (1 - 1 / num_cols)
    return (
      2 * x * x_scale_insertion,
      0.5 * y * y_scale_insertion,
    )
  elif var_type == 'deletion':
    # Place the x coordinate so that the most upstream deletion
    # is the left most, and most downstream deletion is right most.
    # A deletion with equal number of deletions on either side of the
    # cut position should be placed at x = 0.
    first_del_pos = alignment_utils.get_first_deletion_pos(read_align)
    last_del_pos = first_del_pos + dist_ref - 1
    avg_del_pos = (first_del_pos + last_del_pos) / 2
    x = avg_del_pos - (cut_pos_ref + 0.5)
    y = -dist_ref
    return (
      x * x_scale_deletion,
      y * y_scale_deletion,
    )

def make_universal_layout(
  graph,
  cut_pos_ref,
  reverse_complement = False,
  x_scale_insertion = constants.GRAPH_UNIVERSAL_LAYOUT_X_SCALE_INSERTION,
  y_scale_insertion = constants.GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_INSERTION,
  x_scale_deletion = constants.GRAPH_UNIVERSAL_LAYOUT_X_SCALE_DELETION,
  y_scale_deletion = constants.GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_DELETION,
):
  node_list = graph.nodes(data=True)

  bucket_dict = {
    'insertion': {},
    'deletion': {},
  }

  ref_nodes = []

  for _, data in node_list:
    if data['dist_ref'] == 0:
      ref_nodes.append(data)
    else:
      dist_ref = data['dist_ref']
      var_type = data['variation_type']
      if var_type not in bucket_dict:
        raise Exception('Unhandled variation type for universal layout: ' + str(var_type))
      bucket_dict[var_type].setdefault(dist_ref, [])
      bucket_dict[var_type][dist_ref].append(data)

  xy_dict = {}
  for data in ref_nodes:
    xy_dict[data['id']] = (0, 0)
  for var_type in bucket_dict:
    for dist_ref in bucket_dict[var_type]:
      bucket = list(sorted(
        bucket_dict[var_type][dist_ref],
        key = lambda x: x['freq_mean'],
        reverse = True,
      ))
      for data in bucket:
        ref_align = data['ref_align']
        read_align = data['read_align']

        if reverse_complement:
          ref_align = kmer_utils.reverse_complement(ref_align)
          read_align = kmer_utils.reverse_complement(read_align)

        xy_dict[data['id']] = get_pos_universal_layout(
          ref_align = ref_align,
          read_align = read_align,
          dist_ref = dist_ref,
          var_type = var_type,
          cut_pos_ref = cut_pos_ref,
          x_scale_insertion = x_scale_insertion,
          y_scale_insertion = y_scale_insertion,
          x_scale_deletion = x_scale_deletion,
          y_scale_deletion = y_scale_deletion,
        )
  return xy_dict

def make_universal_layout_y_axis(
  figure,
  x_pos,
  ref_length,
  cut_pos_ref, # should be 1 based!
  max_tick_insertion,
  max_tick_deletion,
  y_range = [float('nan'), float('nan')],
  tick_length = 0.25,
  title_font_size = constants.GRAPH_AXES_TITLE_FONT_SIZE,
  tick_font_size = constants.GRAPH_AXES_TICK_FONT_SIZE,
  font_size_scale = constants.GRAPH_FONT_SIZE_SCALE,
  line_width_px = 4,
  x_scale_insertion = constants.GRAPH_UNIVERSAL_LAYOUT_X_SCALE_INSERTION,
  y_scale_insertion = constants.GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_INSERTION,
  x_scale_deletion = constants.GRAPH_UNIVERSAL_LAYOUT_X_SCALE_DELETION,
  y_scale_deletion = constants.GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_DELETION,
):
  tick_list = [{'dist_ref': 0, 'y_pos': 0}]
  for dist_ref in range(1, max(max_tick_insertion, max_tick_deletion) + 1):
    fake_ref_align = (
      ('A' * cut_pos_ref) +
      ('-' * dist_ref) +
      ('A' * (ref_length - cut_pos_ref))
    )
    fake_read_align = 'A' * (ref_length + dist_ref)

    for var_type in ['insertion', 'deletion']:
      if var_type == 'deletion':
        fake_ref_align, fake_read_align = fake_read_align, fake_ref_align

      if (
        ((var_type == 'insertion') and (dist_ref > max_tick_insertion)) or
        ((var_type == 'deletion') and (dist_ref > max_tick_deletion))
      ):
        continue

      y_pos = get_pos_universal_layout(
        ref_align = fake_ref_align,
        read_align = fake_read_align,
        dist_ref = dist_ref,
        var_type = var_type,
        cut_pos_ref = cut_pos_ref,
        x_scale_insertion = x_scale_insertion,
        y_scale_insertion = y_scale_insertion,
        x_scale_deletion = x_scale_deletion,
        y_scale_deletion = y_scale_deletion,
      )[1]

      tick_list.append(
        {
          'dist_ref': dist_ref,
          'y_pos': y_pos,
        }
      )

  for tick in tick_list:
    # tick line
    figure.add_shape(
      type = 'line',
      x0 = x_pos,
      x1 = x_pos + tick_length,
      y0 = tick['y_pos'],
      y1 = tick['y_pos'],
      line_width = line_width_px,
      line_color = 'black',
    )

    # tick label
    figure.add_annotation(
      # set left anchor at tick end
      x = x_pos + tick_length,
      y = tick['y_pos'],
      text = str(tick['dist_ref']),
      showarrow = False,
      font_size = tick_font_size * font_size_scale,
      xanchor = 'left',
      yanchor = 'middle',
      # shift left anchor slighly to the right of tick
      xshift = 0.25 * tick_font_size * font_size_scale
    )

  if np.isnan(y_range[0]):
    y_range[0] = min(tick['y_pos'] for tick in tick_list)
  if np.isnan(y_range[1]):
    y_range[1] = max(tick['y_pos'] for tick in tick_list)

  # axis line
  figure.add_shape(
    type = 'line',
    x0 = x_pos,
    x1 = x_pos,
    y0 = y_range[0],
    y1 = y_range[1],
    line_width = line_width_px,
    line_color = 'black',
  )
  # axis title
  figure.add_annotation(
    # set anchor at bottom of the axis
    x = x_pos,
    y = y_range[0],
    text = str('Number of variations'),
    showarrow = False,
    font_size = title_font_size * font_size_scale,
    xanchor = 'right',
    yanchor = 'bottom',
    # shift left anchor slighly to the right of axis
    xshift = -0.25 * title_font_size * font_size_scale,
    textangle = -90,
  )

# FIXME: should allow line_width_scale parameter
def make_universal_layout_x_axis(
  figure,
  var_type,
  y_pos,
  ref_length,
  cut_pos_ref, # should be 1 based!
  x_range = [float('nan'), float('nan')],
  insertion_axis_type = 'bracket', # tick or bracket
  deletion_label_type = 'relative', # relative or absolute
  deletion_tick_type = 'start', # start or midpoint
  base_tick_length = 0.25,
  title_font_size = constants.GRAPH_AXES_TITLE_FONT_SIZE,
  tick_font_size = constants.GRAPH_AXES_TICK_FONT_SIZE,
  font_size_scale = constants.GRAPH_FONT_SIZE_SCALE,
  line_width_px = 4,
  x_scale_insertion = constants.GRAPH_UNIVERSAL_LAYOUT_X_SCALE_INSERTION,
  y_scale_insertion = constants.GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_INSERTION,
  x_scale_deletion = constants.GRAPH_UNIVERSAL_LAYOUT_X_SCALE_DELETION,
  y_scale_deletion = constants.GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_DELETION,
):
  if insertion_axis_type not in ['tick', 'bracket']:
    raise Exception('Invalid insertion axis type: ' + str(insertion_axis_type))
  if deletion_label_type not in ['absolute', 'relative']:
    raise Exception('Invalid deletion label type: ' + str(deletion_label_type))
  if var_type == 'insertion':
    tick_font_size = 2 * tick_font_size
  tick_list = []
  if var_type == 'insertion':
    for insertion_letter in 'ACGT':
      fake_ref_align = (
        ('A' * cut_pos_ref) +
        '-' +
        ('A' * (ref_length - cut_pos_ref))
      )
      fake_read_align = (
        ('A' * cut_pos_ref) +
        insertion_letter +
        ('A' * (ref_length - cut_pos_ref))
      )
      x_pos = get_pos_universal_layout(
        ref_align = fake_ref_align,
        read_align = fake_read_align,
        dist_ref = 1,
        var_type = 'insertion',
        cut_pos_ref = cut_pos_ref,
        x_scale_insertion = x_scale_insertion,
        y_scale_insertion = y_scale_insertion,
        x_scale_deletion = x_scale_deletion,
        y_scale_deletion = y_scale_deletion,
      )[0]
      tick_list.append({
        'x_pos': x_pos,
        'text': insertion_letter,
        'omit_tick': insertion_axis_type == 'bracket',
        'tick_length': (
          base_tick_length / 2 # used only for spacing in bracket mode
          if insertion_axis_type == 'bracket'
          else base_tick_length
        ),
      })
    if insertion_axis_type == 'bracket':
      tick_space = tick_list[1]['x_pos'] - tick_list[0]['x_pos']
      bracket_x_start = tick_list[0]['x_pos'] - tick_space / 2
      for i in range(5):
        tick_list.append({
          'x_pos': bracket_x_start + tick_space * i,
          'tick_length': 4 * base_tick_length,
        })
  elif var_type == 'deletion':
    tick_list_negative = []
    pos_labels = constants.get_position_labels(deletion_label_type, ref_length)
    for deletion_start in range(1, (ref_length // 2) + 1):
      dist_ref = 1 + (ref_length // 2) - deletion_start
      if (deletion_tick_type == 'midpoint') and ((dist_ref % 2) != 1):
        # only odd deletions so mid point is on an integer
        continue
      deletion_mid = deletion_start + (dist_ref - 1) // 2
      fake_ref_align = 'A' * ref_length
      fake_read_align = (
        'A' * (deletion_start - 1) +
        '-' * dist_ref +
        'A' * (ref_length // 2)
      )
      x_pos = get_pos_universal_layout(
        fake_ref_align,
        fake_read_align,
        dist_ref,
        'deletion',
        cut_pos_ref,
      )[0]

      if deletion_tick_type == 'midpoint':
        final_pos = deletion_mid
      elif deletion_tick_type == 'start':
        final_pos = deletion_start
      else:
        raise Exception(
          'Unsupported deletion tick type: ' +
          str(deletion_tick_type)
        )
      tick_list_negative.append({
        'x_pos': x_pos,
        'text': pos_labels[final_pos - 1],
        'deletion_pos': final_pos,
      })
    
    tick_list_positive = [
      {
        'x_pos': -tick['x_pos'],
        'text': pos_labels[ref_length - tick['deletion_pos']]
      }
      for tick in tick_list_negative
    ]

    tick_list = (
      tick_list_negative +
      [{'x_pos': 0, 'text': '0'}] +
      tick_list_positive
    )

  for tick in tick_list:
    tick_length = tick.get('tick_length', base_tick_length)

    if not tick.get('omit_tick', False):
      # tick line
      figure.add_shape(
        type = 'line',
        x0 = tick['x_pos'],
        x1 = tick['x_pos'],
        y0 = y_pos,
        y1 = y_pos - tick_length,
        line_width = line_width_px,
        line_color = 'black',
      )

    # tick label
    if 'text' in tick:
      figure.add_annotation(
        x = tick['x_pos'],
        # set top anchor at the tick
        y = y_pos - tick_length,
        text = str(tick['text']),
        showarrow = False,
        font_size = tick_font_size * font_size_scale, 
        xanchor = 'center',
        yanchor = 'top',
        # shift top anchor slightly under the tick
        yshift = -0.25 * tick_font_size * font_size_scale,
      )

  if np.isnan(x_range[0]):
    x_range[0] = float('inf')
  if np.isnan(x_range[1]):
    x_range[1] = -float('inf')
  x_range[0] = min(x_range[0], min(tick['x_pos'] for tick in tick_list))
  x_range[1] = max(x_range[1], max(tick['x_pos'] for tick in tick_list))
  # axis line
  figure.add_shape(
    type = 'line',
    x0 = x_range[0],
    x1 = x_range[1],
    y0 = y_pos,
    y1 = y_pos,
    line_width = line_width_px,
    line_color = 'black',
  )
  # axis title
  figure.add_annotation(
    x = x_range[0],
    # set bottom anchor exactly at axis
    y = y_pos,
    text = 'Deletion position' if var_type == 'deletion' else 'Insertion first letter',
    showarrow = False,
    font_size = title_font_size * font_size_scale,
    xanchor = 'left',
    yanchor = 'bottom',
    # shift bottom anchor slightly above axis
    yshift = 0.25 * title_font_size * font_size_scale,
  )

def make_mds_layout(graph, distance_matrix):
  seq_ids = list(graph.nodes())
  seq_ids_set = set(seq_ids)
  distance_matrix = distance_matrix.loc[
    distance_matrix['id_a'].isin(seq_ids_set) &
    distance_matrix['id_b'].isin(seq_ids_set)
  ]
  distance_matrix_transpose = distance_matrix.copy()
  distance_matrix_transpose[['id_b', 'id_a']] = distance_matrix[['id_a', 'id_b']]
  distance_matrix_diag = pd.DataFrame({
    'id_a': seq_ids,
    'id_b': seq_ids,
    'dist': 0,
  })
  distance_matrix = pd.concat(
    [distance_matrix, distance_matrix_transpose, distance_matrix_diag],
    axis = 'index'
  )
  distance_matrix = pd.pivot(distance_matrix, index='id_a', columns='id_b', values='dist')

  seq_ids = list(seq_ids)
  distance_matrix = distance_matrix.reindex(seq_ids)
  distance_matrix = distance_matrix[seq_ids]
  
  mds_out = (
    sklearn.manifold.MDS(n_components=2, dissimilarity='precomputed')
      .fit_transform(distance_matrix)
  )

  xy_dict = {}
  for i in range(len(seq_ids)):
    xy_dict[seq_ids[i]] = (mds_out[i, 0], mds_out[i, 1])
  return xy_dict

def make_graph_layout_single(
  data_info,
  graph,
  layout_type,
  distance_matrix = None,
  reverse_complement = False,
  universal_layout_x_scale_insertion = constants.GRAPH_UNIVERSAL_LAYOUT_X_SCALE_INSERTION,
  universal_layout_y_scale_insertion = constants.GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_INSERTION,
  universal_layout_x_scale_deletion = constants.GRAPH_UNIVERSAL_LAYOUT_X_SCALE_DELETION,
  universal_layout_y_scale_deletion = constants.GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_DELETION,
):
  if layout_type == 'mds_layout':
    layout = make_mds_layout(graph, distance_matrix)
  elif layout_type == 'radial_layout':
    layout = make_radial_layout(graph)
  elif layout_type == 'universal_layout':
    layout = make_universal_layout(
      graph,
      len(data_info['ref_seq_window']) // 2,
      reverse_complement = reverse_complement,
      x_scale_insertion = universal_layout_x_scale_insertion,
      y_scale_insertion = universal_layout_y_scale_insertion,
      x_scale_deletion = universal_layout_x_scale_deletion,
      y_scale_deletion = universal_layout_y_scale_deletion,
    )
  elif layout_type == 'fractal_layout':
    layout = make_fractal_layout(
      graph,
      reverse_complement = reverse_complement,
    )
  elif layout_type == 'kamada_layout':
    layout = nx.kamada_kawai_layout(graph, dim = 2)
  elif layout_type == 'spectral_layout':
    layout = nx.spectral_layout(graph, dim=2)
  elif layout_type == 'spring_layout':
    layout = nx.spring_layout(graph, dim=2)
  elif layout_type == 'shell_layout':
    layout = nx.shell_layout(
      graph,
      dim = 2,
      nlist = group_graph_nodes_by(graph, 'dist_ref'),
    )
  elif layout_type == 'spiral_layout':
    layout = nx.spiral_layout(graph, dim=2)
  elif layout_type == 'circular_layout':
    layout = nx.circular_layout(graph, dim=2)
  elif layout_type == 'multipartite_layout':
    layout = nx.multipartite_layout(graph, subset_key='dist_ref')
  else:
    raise Exception('Unknown layout type: ' + str(layout_type))

  layout = pd.DataFrame.from_dict(layout, orient='index', columns=[0, 1])

  if (layout.shape[0] >= 2) and LAYOUT_PROPERTIES[layout_type]['do_pca']:
    layout = pd.DataFrame(
      data = (
        sklearn.decomposition.PCA(n_components=2)
          .fit_transform(layout.to_numpy())
      ),
      index = layout.index,
      columns = [0, 1],
    )

  if LAYOUT_PROPERTIES[layout_type]['normalize']:
    dim_mins = []
    scales = []
    for i in range(layout.shape[1]):
      dim_min = np.min(layout.loc[:, i])
      dim_max = np.max(layout.loc[:, i])
      dim_mins.append(dim_min)
      if np.isclose(dim_min, dim_max):
        scales.append(0)
      else:
        scales.append(1 / (dim_max - dim_min))
    if LAYOUT_PROPERTIES[layout_type].get('preserve_aspect', False):
      scales = [min(scales)] * layout.shape[1]
    for i in range(layout.shape[1]):
      layout.loc[:, i] = (layout.loc[:, i] - dim_mins[i]) * scales[i]
  return layout


def make_grid_spec(
  num_panels,
  major_panel,
  major_panel_size = 0.85,
  invert_y = True,
  panel_pad = 0.05
):
  grid_spec = []
  if major_panel:
    grid_spec.append({
      'x': 0,
      'y': 0,
      'height': major_panel_size,
      'width': major_panel_size,
    })
    num_minor_panels = num_panels - 1
    top_panels = (num_minor_panels + 1) // 2
    side_panels = num_minor_panels // 2
    
    grid_x = np.linspace(0, 1, top_panels + 1)[:-1]
    for i in range(top_panels):
      grid_spec.append({
        'x': grid_x[i],
        'y': major_panel_size + panel_pad,
        'width': 1 / top_panels - panel_pad,
        'height': 1 - (major_panel_size + panel_pad),
      })
    
    if side_panels > 0:
      grid_y = np.linspace(0, major_panel_size, side_panels + 1)[:-1]
      for i in range(side_panels):
        grid_spec.append({
          'x': major_panel_size + panel_pad,
          'y': grid_y[i],
          'width': 1 - (major_panel_size + panel_pad),
          'height': major_panel_size / side_panels - panel_pad,
        })
  else:
    num_grid_1d = int(np.ceil(np.sqrt(num_panels)))
    grid_1d = np.linspace(0, 1, num_grid_1d + 1)[:-1]
    grid_x, grid_y = np.meshgrid(grid_1d, grid_1d)
    grid_x = grid_x.ravel()
    grid_y = grid_y.ravel()
    grid_spec = [
      {
        'x': x,
        'y': y,
        'height': 1 / num_grid_1d,
        'width': 1 / num_grid_1d,
      }
      for x, y in zip(grid_x, grid_y)
    ]
  if invert_y:
    for grid in grid_spec:
      grid['y'] = 1 - grid['y'] - grid['height']
  grid_spec = grid_spec[:num_panels]
  return grid_spec

def make_graph_layout(
  data_dir,
  data_info,
  node_subst_type,
  graph,
  layout_type,
  graph_layout_common = None,
  separate_components = True,
  reverse_complement = False,
  universal_layout_x_scale_insertion = constants.GRAPH_UNIVERSAL_LAYOUT_X_SCALE_INSERTION,
  universal_layout_y_scale_insertion = constants.GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_INSERTION,
  universal_layout_x_scale_deletion = constants.GRAPH_UNIVERSAL_LAYOUT_X_SCALE_DELETION,
  universal_layout_y_scale_deletion = constants.GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_DELETION,
):
  if graph_layout_common is not None:
    separate_components = False
    node_groups = None
    node_data = pd.DataFrame.from_dict(
      dict(graph.nodes(data=True)),
      orient = 'index',
    ).reset_index(drop=True)
    if reverse_complement:
      node_data = node_data.assign(
        ref_align = node_data['ref_align'].apply(kmer_utils.reverse_complement),
        read_align = node_data['read_align'].apply(kmer_utils.reverse_complement),
      )
    node_data = pd.merge(
      node_data[['id', 'ref_align', 'read_align']],
      graph_layout_common[['ref_align', 'read_align', 'x', 'y']],
      on = ['ref_align', 'read_align'],
      how = 'inner',
    )[['id', 'x', 'y']]
    node_data = node_data.set_index('id', drop=True).rename(
      {'x': 0, 'y': 1},
      axis = 'columns',
    )
    layout_list = [node_data]
  else:
    if separate_components:
      ref_id = next(
        (id for id, data in graph.nodes(data=True) if data['is_ref']),
        None,
      )
      node_groups = list(nx.connected_components(graph))
      if ref_id is not None:
        node_groups = (
          [x for x in node_groups if ref_id in x] +
          [x for x in node_groups if ref_id not in x]
        )
      subgraph_list = [graph.subgraph(group) for group in node_groups]
    else:
      node_groups = None
      subgraph_list = [graph]

    if LAYOUT_PROPERTIES[layout_type].get('distance_matrix', False):
      distance_matrix = file_utils.read_tsv(
        file_names.distance_matrix(data_dir, node_subst_type)
      )
    else:
      distance_matrix = None
    layout_list = [
      make_graph_layout_single(
        data_info = data_info,
        graph = subgraph,
        layout_type = layout_type,
        distance_matrix = distance_matrix,
        reverse_complement = reverse_complement,
        universal_layout_x_scale_insertion = universal_layout_x_scale_insertion,
        universal_layout_y_scale_insertion = universal_layout_y_scale_insertion,
        universal_layout_x_scale_deletion = universal_layout_x_scale_deletion,
        universal_layout_y_scale_deletion = universal_layout_y_scale_deletion,
      )
      for subgraph in subgraph_list
    ]

  if separate_components:
    if ref_id is not None:
      grid_spec = make_grid_spec(len(node_groups), True)
    else:
      grid_spec = make_grid_spec(len(node_groups), False)
    for layout, panel in zip(layout_list, grid_spec):
      layout.loc[:, 0] = layout.loc[:, 0] * panel['width'] + panel['x']
      layout.loc[:, 1] = layout.loc[:, 1] * panel['height'] + panel['y']

  layout = pd.concat(layout_list, axis='index')

  # Center the whole thing a bit
  if LAYOUT_PROPERTIES[layout_type]['normalize']:
    if layout.shape[0] < 10:
      layout = layout.applymap(lambda x: 0.33 + 0.33 * x)
    else:
      layout = layout.applymap(lambda x: 0.1 + 0.8 * x)

  return layout

def make_legend(
  figure,
  legend_title,
  legend_items,
  x_anchor,
  y_anchor,
  x_shift,
  y_shift,
  x_shift_items,
  y_shift_items,
  x_shift_text,
  y_shift_item_step,
  legend_item_scale = 1,
  font_size_scale = 1,
  line_width_scale = 1,
):
  figure.add_annotation(
    text = legend_title,
    xref = 'paper',
    yref = 'paper',
    x = x_anchor,
    y = y_anchor,
    xanchor = 'left',
    yanchor = 'middle',
    xshift = x_shift,
    yshift = y_shift,
    showarrow = False,
    font_size = constants.GRAPH_LEGEND_TITLE_FONT_SIZE * font_size_scale,
  )

  y_shift_step_sign = -1 if y_shift_item_step < 0 else 1
  y_shift_item_step = legend_item_scale * y_shift_item_step
  curr_y_shift = (
    y_shift +
    y_shift_step_sign * constants.GRAPH_LEGEND_TITLE_FONT_SIZE * font_size_scale +
    y_shift_items
  )
  for _, item in enumerate(legend_items):
    if item['type'] == 'circle':
      figure.add_shape(
        type = 'circle',
        xref = 'paper',
        yref = 'paper',
        x0 = x_shift + legend_item_scale * (x_shift_items - item['size'] / 2),
        y0 = curr_y_shift - legend_item_scale * item['size'] / 2,
        x1 = x_shift + legend_item_scale * (x_shift_items + item['size'] / 2),
        y1 = curr_y_shift + legend_item_scale * item['size'] / 2,
        xsizemode = 'pixel',
        ysizemode = 'pixel',
        xanchor = x_anchor,
        yanchor = y_anchor,
        line_color = item.get('line_color', 'black'),
        line_width = item.get('line_width', 1) * font_size_scale,
        fillcolor = item['color'],
      )
    elif item['type'] == 'line':
      figure.add_shape(
        type = 'line',
        xref = 'paper',
        yref = 'paper',
        x0 = x_shift + legend_item_scale * (x_shift_items - item['size'] / 2),
        y0 = curr_y_shift,
        x1 = x_shift + legend_item_scale * (x_shift_items + item['size'] / 2),
        y1 = curr_y_shift,
        xsizemode = 'pixel',
        ysizemode = 'pixel',
        xanchor = x_anchor,
        yanchor = y_anchor,
        line_color = item['color'],
        line_width = item['line_width'] * line_width_scale,
        line_dash = item['line_dash'],
      )
    else:
      raise Exception('Unhandled item type: ' + str(item['type']))
    
    figure.add_annotation(
      text = item['text'],
      xref = 'paper',
      yref = 'paper',
      x = x_anchor,
      y = y_anchor,
      xshift = x_shift + legend_item_scale * (x_shift_items + x_shift_text),
      yshift = curr_y_shift,
      xanchor = 'left',
      yanchor = 'middle',
      showarrow = False,
      font_size = constants.GRAPH_LEGEND_FONT_SIZE * font_size_scale,
    )
    curr_y_shift += y_shift_step_sign * max(
      abs(y_shift_item_step),
      1.5 * constants.GRAPH_LEGEND_FONT_SIZE * font_size_scale,
    )
  return curr_y_shift

def make_edge_legend(
  figure,
  edge_type_list,
  line_size_px,
  line_width_px,
  x_anchor,
  y_anchor,
  x_shift,
  y_shift,
  legend_item_scale = 1,
  font_size_scale = 1,
  line_width_scale = 1,
):
  legend_items = []
  for edge_type in edge_type_list:
    legend_items.append({
      'type': 'line',
      'size': line_size_px,
      'text': constants.EDGE_TYPES[edge_type]['label'],
      'color': constants.EDGE_TYPES[edge_type]['legend_color'],
      'line_dash': constants.EDGE_TYPES[edge_type]['line_dash'],
      'line_width': line_width_px,
    })
  return make_legend(
    figure = figure,
    legend_title = 'Edges',
    legend_items = legend_items,
    x_anchor = x_anchor,
    y_anchor = y_anchor,
    x_shift = x_shift,
    y_shift = y_shift,
    x_shift_items = line_size_px / 2,
    y_shift_items = -50,
    x_shift_text = line_size_px + 10,
    y_shift_item_step = -30,
    legend_item_scale = legend_item_scale,
    font_size_scale = font_size_scale,
    line_width_scale = line_width_scale,
  )

def make_variation_color_legend(
  figure,
  variation_types,
  variation_type_colors,
  node_size_px,
  x_anchor,
  y_anchor,
  x_shift,
  y_shift,
  legend_item_scale = 1,
  font_size_scale = 1,
  line_width_scale = 1,
):
  legend_items = []
  for var_type in variation_types:
    legend_items.append({
      'type': 'circle',
      'size': node_size_px,
      'text': constants.VARIATION_TYPES[var_type]['label'],
      'color': variation_type_colors[var_type],
    })
  return make_legend(
    figure = figure,
    legend_title = 'Variation Types',
    legend_items = legend_items,
    x_anchor = x_anchor,
    y_anchor = y_anchor,
    x_shift = x_shift,
    y_shift = y_shift,
    x_shift_items = node_size_px / 2,
    y_shift_items = -(node_size_px + 10),
    x_shift_text = node_size_px + 10,
    y_shift_item_step = -(node_size_px + 10),
    legend_item_scale = legend_item_scale,
    font_size_scale = font_size_scale,
    line_width_scale = line_width_scale,
  )

def make_outline_legend(
  figure,
  node_size_px,
  x_anchor,
  y_anchor,
  x_shift,
  y_shift,
  node_reference_outline_color,
  node_outline_color,
  node_fill_color,
  legend_item_scale = 1,
  font_size_scale = 1,
  line_width_scale = 1,
):
  legend_items = []
  legend_items.append({
    'type': 'circle',
    'size': node_size_px,
    'text': constants.LABEL_REFERENCE,
    'color': node_fill_color,
    'line_color': node_reference_outline_color,
    'line_width': constants.GRAPH_NODE_REFERENCE_OUTLINE_WIDTH,
  })
  legend_items.append({
    'type': 'circle',
    'size': node_size_px,
    'text': constants.LABEL_NONREFERENCE,
    'color': node_fill_color,
    'line_color': node_outline_color,
    'line_width': constants.GRAPH_NODE_OUTLINE_WIDTH,
  })
  return make_legend(
    figure = figure,
    legend_title = 'Vertex Outline',
    legend_items = legend_items,
    x_anchor = x_anchor,
    y_anchor = y_anchor,
    x_shift = x_shift,
    y_shift = y_shift,
    x_shift_items = node_size_px / 2,
    y_shift_items = -(node_size_px + 10),
    x_shift_text = node_size_px + 10,
    y_shift_item_step = -(node_size_px + 10),
    legend_item_scale = legend_item_scale,
    font_size_scale = font_size_scale,
    line_width_scale = line_width_scale,
  )

def make_size_legend(
  figure,
  node_size_freq_range,
  node_size_px_range,
  x_anchor,
  y_anchor,
  x_shift,
  y_shift,
  legend_item_scale = 1,
  font_size_scale = 1,
  line_width_scale = 1,
):
  node_size_freq_range_log10 = np.round(np.log10(node_size_freq_range)).astype(int)

  num_legend_items = node_size_freq_range_log10[1] - node_size_freq_range_log10[0] + 1

  legend_items = []
  for i in range(num_legend_items):
    freq_log10 = node_size_freq_range_log10[0] + i
    if num_legend_items == 1:
      size = node_size_px_range[0]
    else:
      size = node_size_px_range[0] + (
        i * (node_size_px_range[1] - node_size_px_range[0]) /
        (num_legend_items - 1)
      )
    if freq_log10 == 0:
      text = '1'
    else:
      text = f'10<sup>{freq_log10}</sup>'
    if i == num_legend_items - 1:
      text = '≥' + text
    elif i == 0:
      text = '≤' + text
    legend_items.append({
      'type': 'circle',
      'size': size,
      'text': text,
      'color': 'white',
    })
  legend_items = legend_items[::-1] # Show largest to smallest
  return make_legend(
    figure = figure,
    legend_title = 'Frequency',
    legend_items = legend_items,
    x_anchor = x_anchor,
    y_anchor = y_anchor,
    x_shift = x_shift,
    y_shift = y_shift,
    x_shift_items = node_size_px_range[1] / 2,
    y_shift_items = -(node_size_px_range[1] + 10),
    x_shift_text = legend_item_scale * (node_size_px_range[1] + 10),
    y_shift_item_step = -(node_size_px_range[1] + 10),
    font_size_scale = font_size_scale,
    line_width_scale = line_width_scale,
  )

def make_freq_group_legend(
  label_1,
  label_2,
  color_1,
  color_2,
  freq_ratio_1,
  freq_ratio_2,
  figure,
  node_size_px,
  x_anchor,
  y_anchor,
  x_shift,
  y_shift,
  legend_item_scale = 1,
  font_size_scale = 1,
  line_width_scale = 1,
):
  legend_items = []
  for group, color in [
    ('A', color_1),
    ('B', constants.SIMILAR_FREQ_COLOR),
    ('C', color_2),
  ]:
    legend_items.append({
      'type': 'circle',
      'size': node_size_px,
      'text': constants.get_freq_ratio_label(
        group,
        label_1,
        label_2,
        freq_ratio_1,
        freq_ratio_2,
      ),
      'color': color,
    })
  return make_legend(
    figure = figure,
    legend_title = 'Vertex Color',
    legend_items = legend_items,
    x_anchor = x_anchor,
    y_anchor = y_anchor,
    x_shift = x_shift,
    y_shift = y_shift,
    x_shift_items = node_size_px / 2,
    y_shift_items = -50,
    x_shift_text = node_size_px + 10,
    y_shift_item_step = -(node_size_px + 10),
    legend_item_scale = legend_item_scale,
    font_size_scale = font_size_scale,
    line_width_scale = line_width_scale,
  )

def add_plotly_colorbar(
  figure,
  label_1,
  label_2,
  freq_ratio_1,
  freq_ratio_2,
  content_height_px,
  legend_colorbar_scale = 1,
  legend_x_shift_px = 0,
  legend_y_shift_px = 0,
  line_width_scale = 1,
  font_size_scale = 1,
):
  # Note: Sometimes the entire plot disappears if the colorbar font is too large!
  # Fixes: Increase the colorbar length or make the fonts smaller.
  colorbar_height_px = (
    constants.GRAPH_LEGEND_COLORBAR_HEIGHT_PX *
    legend_colorbar_scale
  )

  colorbar_width_px =  (
    constants.GRAPH_LEGEND_COLORBAR_WIDTH_PX *
    legend_colorbar_scale
  )

  figure.update_traces(
    marker = {
      'colorbar': {    
        'x': 1,
        'y': 1 + legend_y_shift_px / content_height_px,
        'xpad': legend_x_shift_px,
        'ypad': 0,
        'xanchor': 'left',
        'yanchor': 'top',
        'orientation': 'v',
        'lenmode': 'pixels',
        'len': colorbar_height_px,
        'thickness': colorbar_width_px,
        'outlinewidth': 2 * line_width_scale,
        'outlinecolor': 'black',
        'tickmode': 'array',
        'tickvals': [np.log(freq_ratio_1), 0, np.log(freq_ratio_2)],
        'ticktext': [f'{freq_ratio_1:.2f}', '1', f'{freq_ratio_2:.2f}'],
        'title': {
          'text': 'Frequency Ratio<br>' + f'[{label_1} / {label_2}]',
          'font_size': constants.GRAPH_LEGEND_TITLE_FONT_SIZE * font_size_scale,
        },
        'tickfont_size': constants.GRAPH_LEGEND_FONT_SIZE * font_size_scale,
      },
    },
  )
  return legend_y_shift_px - colorbar_height_px

def make_custom_legends(
  figure,
  legend_list,
  content_height_px,
  data_info,
  node_reference_outline_color,
  node_outline_color,
  node_fill_color,
  node_comparison_colors,
  node_variation_type_colors,
  node_filter_variation_types,
  node_size_freq_range,
  node_size_px_range,
  node_freq_ratio_range,
  edge_show_types,
  legend_x_shift_px,
  legend_y_shift_px,
  legend_vertical_space_px,
  legend_item_scale,
  legend_colorbar_scale,
  font_size_scale,
  line_width_scale,
  x_anchor_frac = constants.GRAPH_LEGEND_CUSTOM_X_ANCHOR_FRAC,
  y_anchor_frac = constants.GRAPH_LEGEND_CUSTOM_Y_ANCHOR_FRAC,
):
  y_shift_curr_px = legend_y_shift_px

  for legend in legend_list:
    if legend == 'outline':
      y_shift_curr_px = make_outline_legend(
        figure = figure,
        node_size_px = node_size_px_range[1],
        x_anchor = x_anchor_frac,
        y_anchor = y_anchor_frac,
        x_shift = legend_x_shift_px,
        y_shift = y_shift_curr_px,
        node_reference_outline_color = node_reference_outline_color,
        node_outline_color = node_outline_color,
        node_fill_color = node_fill_color,
        legend_item_scale = legend_item_scale,
        font_size_scale = font_size_scale,
        line_width_scale = line_width_scale,
      )
      y_shift_curr_px -= legend_vertical_space_px

    if (
      (legend == 'variation_type') and
      (data_info['format'] == constants.DATA_INDIVIDUAL)
    ):
      y_shift_curr_px = make_variation_color_legend(
        figure = figure,
        variation_types = node_filter_variation_types,
        variation_type_colors = node_variation_type_colors,
        node_size_px = node_size_px_range[1],
        x_anchor = x_anchor_frac,
        y_anchor = y_anchor_frac,
        x_shift = legend_x_shift_px,
        y_shift = y_shift_curr_px,
        legend_item_scale = legend_item_scale,
        font_size_scale = font_size_scale,
        line_width_scale = line_width_scale,
      )
      y_shift_curr_px -= legend_vertical_space_px

    if (
      (legend == 'freq_ratio_discrete') and
      (data_info['format'] == constants.DATA_COMPARISON)
    ):
      y_shift_curr_px = make_freq_group_legend(
        label_1 = data_info['label_1'],
        label_2 = data_info['label_2'],
        color_1 = node_comparison_colors[0],
        color_2 = node_comparison_colors[1],
        freq_ratio_1 = node_freq_ratio_range[0],
        freq_ratio_2 = node_freq_ratio_range[1],
        figure = figure,
        node_size_px = node_size_px_range[1],
        x_anchor = x_anchor_frac,
        y_anchor = y_anchor_frac,
        x_shift = legend_x_shift_px,
        y_shift = y_shift_curr_px,
        legend_item_scale = legend_item_scale,
        font_size_scale = font_size_scale,
        line_width_scale = line_width_scale,
      )
      y_shift_curr_px -= legend_vertical_space_px
    
    if (
      (legend == 'freq_ratio_continuous') and
      (data_info['format'] == constants.DATA_COMPARISON)
    ):
      y_shift_curr_px = add_plotly_colorbar(
        figure = figure,
        label_1 = data_info['label_1'],
        label_2 = data_info['label_2'],
        freq_ratio_1 = node_freq_ratio_range[0],
        freq_ratio_2 = node_freq_ratio_range[1],
        content_height_px = content_height_px,
        legend_colorbar_scale = legend_colorbar_scale,
        legend_x_shift_px = legend_x_shift_px,
        legend_y_shift_px = y_shift_curr_px,
        line_width_scale = line_width_scale,
        font_size_scale = font_size_scale,
      )
      y_shift_curr_px -= legend_vertical_space_px
    
    if legend == 'edge':
      y_shift_curr_px = make_edge_legend(
        figure = figure,
        edge_type_list = edge_show_types,
        line_size_px = constants.GRAPH_LEGEND_EDGE_ITEM_LINE_SIZE_PX,
        line_width_px = constants.GRAPH_LEGEND_EDGE_ITEM_LINE_WIDTH_PX,
        x_anchor = x_anchor_frac,
        y_anchor = y_anchor_frac,
        x_shift = legend_x_shift_px,
        y_shift = y_shift_curr_px,
        legend_item_scale = legend_item_scale,
        font_size_scale = font_size_scale,
        line_width_scale = line_width_scale,
      )
      y_shift_curr_px -= legend_vertical_space_px

    if legend == 'size':
      y_shift_curr_px = make_size_legend(
        figure = figure,
        node_size_freq_range = node_size_freq_range,
        node_size_px_range = node_size_px_range,
        x_anchor = x_anchor_frac,
        y_anchor = y_anchor_frac,
        x_shift = legend_x_shift_px,
        y_shift = y_shift_curr_px,
        legend_item_scale = legend_item_scale,
        font_size_scale = font_size_scale,
        line_width_scale = line_width_scale,
      )
      y_shift_curr_px -= legend_vertical_space_px
  
  return y_shift_curr_px

def make_graph_stats(
  figure,
  data_dir,
  x,
  y,
  x_shift,
  y_shift,
  x_anchor,
  y_anchor,
  font_size_scale = 1,
):
  graph_stats = file_utils.read_tsv_dict(file_names.graph_stats(data_dir))
  graph_stats = {
    k: 'NA' if pd.isna(v) else
    str(v) if isinstance(v, int) else
    f'{v:.2f}'
    for k, v in graph_stats.items()
  }
  figure.add_annotation(
    xref = 'x domain',
    yref = 'y domain',
    x = x,
    y = y,
    xshift = x_shift,
    yshift = y_shift,
    xanchor = x_anchor,
    yanchor = y_anchor,
    align = 'left',
    font_size = constants.GRAPH_LEGEND_FONT_SIZE * font_size_scale,
    font_family = 'Monospace',
    text = (
      f'Num nodes:            {graph_stats["num_nodes"][0]}<br>'
      f'Num edges:            {graph_stats["num_edges"][0]}<br>'
      f'Avg degree:           {graph_stats["avg_degree"][0]}<br>'
      f'Avg dist from ref:    {graph_stats["avg_dist_ref"][0]}<br>'
      f'Avg pairwise dist:    {graph_stats["avg_pairwise_lev_dist"][0]}<br>'
      f'Max dist from ref:    {graph_stats["max_dist_ref"][0]}<br>'
      f'Max pairwise dist:    {graph_stats["max_pairwise_lev_dist"][0]}<br>'
      f'Num seq substitution: {graph_stats["num_seq_substitution"][0]}<br>'
      f'Num seq insertion:    {graph_stats["num_seq_insertion"][0]}<br>'
      f'Num seq deletion:     {graph_stats["num_seq_deletion"][0]}<br>'
    ),
    showarrow = False,
  )

def make_graph_stats_ref_component(
  figure,
  data_dir,
  data_info,
  subst_type,
  x,
  y,
  x_shift,
  y_shift,
  x_anchor,
  y_anchor,
  font_size_scale = 1,
):
  graph_stats = file_utils.read_tsv_dict(
    file_names.graph_stats(data_dir, subst_type)
  )
  stat_lines = [
    ['Num nodes', graph_stats['num_nodes']],
    ['Num edges', graph_stats['num_edges']],
    ['Avg degree', graph_stats['avg_degree']],
    ['Avg dist from ref', graph_stats['avg_dist_ref']],
    ['Avg pairwise dist', graph_stats['avg_pairwise_dist']],
    ['Max dist from ref', graph_stats['max_dist_ref']],
    ['Max pairwise dist', graph_stats['max_pairwise_dist']],
    ['Num seq insertion', graph_stats['num_seq_insertion']],
    ['Num seq deletion', graph_stats['num_seq_deletion']],
  ]
  if data_info['format'] == 'comparison':
    stat_lines += [
      [
        'Ref seq freq', '{:.3f} & {:.3f}'.format(
          graph_stats['ref_freq_mean_1'],
          graph_stats['ref_freq_mean_2'],
        )
      ],
      [
        'Non-ref seq freq', '{:.3f} & {:.3f}'.format(
          graph_stats['non_ref_freq_mean_1'],
          graph_stats['non_ref_freq_mean_2'],
        )
      ],
      [
        'Insertion freq', '{:.4f} & {:.4f}'.format(
          graph_stats['insertion_freq_mean_1'],
          graph_stats['insertion_freq_mean_2'],
        )
      ],
      [
        'Deletion freq', '{:.4f} & {:.4f}'.format(
          graph_stats['deletion_freq_mean_1'],
          graph_stats['deletion_freq_mean_2'],
        )
      ],
    ]
  elif data_info['format'] == 'individual':
    stat_lines += [
      ['Ref seq freq', '{:.3f}'.format(graph_stats['ref_freq_mean'])],
      ['Non-ref seq freq', '{:.5f}'.format(graph_stats['non_ref_freq_mean'])],
      ['Insertion freq', '{:.5f}'.format(graph_stats['insertion_freq_mean'])],
      ['Deletion freq', '{:.5f}'.format(graph_stats['deletion_freq_mean'])],
    ]
  else:
    raise Exception('Unknown data format: ' + str(data_info['format']))
  for line in stat_lines:
    if not isinstance(line[1], str):
      if pd.isna(line[1]):
        line[1] = 'NA'
      elif line[1] == np.round(line[1]): # integer
        line[1] = str(line[1])
      else:
        line[1] = f'{line[1]:.2f}'
  max_label_len = max(len(line[0]) for line in stat_lines)
  for line in stat_lines:
    line[0] = line[0].ljust(max_label_len) + ': '
  stat_lines = [line[0] + line[1] for line in stat_lines]
  stat_lines = '<br>'.join(stat_lines)
  figure.add_annotation(
    xref = 'x domain',
    yref = 'y domain',
    x = x,
    y = y,
    xshift = x_shift,
    yshift = y_shift,
    xanchor = x_anchor,
    yanchor = y_anchor,
    align = 'left',
    font_size = constants.GRAPH_LEGEND_FONT_SIZE * font_size_scale,
    font_family = 'Monospace',
    text = (
      '<span style="text-decoration: underline;">'
        'Graph Statistics (ref component only)'
      '</span><br>' +
      stat_lines
    ),
    showarrow = False,
  )

def make_graph_figure_helper(
  figure_list,
  data_dir_list,
  data_info_list,
  sequence_reverse_complement_list = None,
  node_subst_type = constants.SUBST_WITHOUT,
  node_filter_freq_range = constants.GRAPH_NODE_FILTER_FREQ_RANGE,
  node_filter_dist_range = constants.GRAPH_NODE_FILTER_DIST_RANGE,
  edge_show = constants.GRAPH_EDGE_SHOW,
  edge_types_show = constants.GRAPH_EDGE_SHOW_TYPES,
  edge_labels_show = constants.GRAPH_EDGE_LABELS_SHOW,
  edge_width_scale = constants.GRAPH_EDGE_WIDTH_SCALE,
  graph_layout_type = constants.GRAPH_LAYOUT_TYPE,
  graph_layout_separate_components = constants.GRAPH_LAYOUT_SEPARATE_COMPONENTS,
  node_labels_show = constants.GRAPH_NODE_LABEL_SHOW,
  node_label_columns = constants.GRAPH_NODE_LABEL_COLUMNS,
  node_label_position = constants.GRAPH_NODE_LABEL_POSITION,
  node_color_type_list = None,
  node_comparison_colors = constants.GRAPH_NODE_COMPARISON_COLORS,
  node_freq_ratio_range = constants.GRAPH_NODE_FREQ_RATIO_RANGE,
  node_reference_outline_color = constants.GRAPH_NODE_REFERENCE_OUTLINE_COLOR,
  node_outline_color = constants.GRAPH_NODE_OUTLINE_COLOR,
  node_fill_color = constants.GRAPH_NODE_FILL_COLOR,
  node_variation_type_colors = constants.GRAPH_NODE_VARIATION_TYPE_COLORS,
  node_size_type = constants.GRAPH_NODE_SIZE_TYPE,
  node_size_px_range = constants.GRAPH_NODE_SIZE_PX_RANGE,
  node_size_freq_range = constants.GRAPH_NODE_SIZE_FREQ_RANGE,
  node_filter_variation_types = constants.GRAPH_NODE_FILTER_VARIATION_TYPES,
  node_outline_width_scale = constants.GRAPH_NODE_OUTLINE_WIDTH_SCALE,
  plot_range_x = constants.GRAPH_PLOT_RANGE_X,
  plot_range_y = constants.GRAPH_PLOT_RANGE_Y,
  legend_plotly_show = constants.GRAPH_LEGEND_CUSTOM_SHOW,
  axes_show = constants.GRAPH_AXES_SHOW,
  font_size_scale = constants.GRAPH_FONT_SIZE_SCALE,
  universal_layout_x_scale_insertion = constants.GRAPH_UNIVERSAL_LAYOUT_X_SCALE_INSERTION,
  universal_layout_y_scale_insertion = constants.GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_INSERTION,
  universal_layout_x_scale_deletion = constants.GRAPH_UNIVERSAL_LAYOUT_X_SCALE_DELETION,
  universal_layout_y_scale_deletion = constants.GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_DELETION,
):
  if sequence_reverse_complement_list is None:
    sequence_reverse_complement_list = (
      [constants.GRAPH_SEQUENCE_REVERSE_COMPLEMENT] * len(data_dir_list)
    )

  ### Load node data ###
  node_data_list = [
    file_utils.read_tsv(file_names.sequence_data(data_dir, node_subst_type))
      .set_index('id', drop=False)
    for data_dir in data_dir_list
  ]

  ### Load graph ###
  graph_list = [
    graph_utils.load_graph(data_dir, node_subst_type)
    for data_dir in data_dir_list
  ]

  ### Node filtering / subgraph ###
  for i in range(len(data_dir_list)):
    if node_filter_variation_types is not None:
      node_data_list[i] = node_data_list[i].loc[
        node_data_list[i]['variation_type'].isin(node_filter_variation_types)
      ]
    node_data_list[i] = node_data_list[i].loc[
      node_data_list[i]['freq_mean']
        .between(node_filter_freq_range[0], node_filter_freq_range[1], inclusive='both')
    ]

    node_data_list[i] = node_data_list[i].loc[
      node_data_list[i]['dist_ref'].between(
        node_filter_dist_range[0],
        node_filter_dist_range[1],
        inclusive = 'both'
      )
    ]

    graph_list[i] = graph_list[i].subgraph(node_data_list[i].index)

  ### Combine the nodes/edges ###

  # Combine node data #
  node_data_list_copy = [node_data.copy() for node_data in node_data_list]
  for i in range(len(data_dir_list)):
    if sequence_reverse_complement_list[i]:
      node_data_list_copy[i] = node_data_list_copy[i].assign(
        ref_align = node_data_list_copy[i]['ref_align'].apply(kmer_utils.reverse_complement),
        read_align = node_data_list_copy[i]['read_align'].apply(kmer_utils.reverse_complement),
      )
  node_data = pd.concat(node_data_list_copy, axis='index', ignore_index=True)
  node_data = node_data.groupby(['ref_align', 'read_align'])
  freq_mean_max = node_data['freq_mean'].max()
  node_data = node_data.first()
  node_data['freq_mean'] = freq_mean_max
  node_data = node_data.sort_values('freq_mean', ascending=False).reset_index()
  node_data['id'] = 'S' + pd.Series(range(1, node_data.shape[0] + 1), dtype=str)
  node_data = node_data.set_index('id', drop=False)

  # Combine edge data
  edge_data_list = [
    file_utils.read_tsv(file_names.edge_data(data_dir, node_subst_type))
      .drop(['id_a', 'id_b'], axis='columns')
    for data_dir in data_dir_list
  ]
  edge_data = pd.concat(edge_data_list, axis='index', ignore_index=True)
  edge_data = edge_data.groupby(list(edge_data.columns)).first().reset_index()

  # Get the new id's for the edges
  for suffix in ['_a', '_b']:
    edge_data = pd.merge(
      edge_data,
      node_data[['id', 'ref_align', 'read_align']].rename(
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
  graph.add_nodes_from(node_data.index)
  nx.set_node_attributes(graph, node_data.to_dict('index'))

  graph.add_edges_from(
    zip(
      edge_data['id_a'],
      edge_data['id_b'],
      edge_data.to_dict('records')
    ),
  )

  ### Make the common graph layout ###
  data_info = file_utils.read_tsv_dict(
    file_names.data_info(data_dir_list[0])
  )
  graph_layout_common = make_graph_layout(
    # FIXME: This mean FAIL for MDS layout! Mention in the args strings!
    data_dir = data_dir_list[0] if (len(data_dir_list) == 0) else None,
    data_info = data_info_list[0],
    node_subst_type = node_subst_type,
    graph = graph,
    layout_type = graph_layout_type,
    separate_components = False,
    graph_layout_common = None,
  )
  graph_layout_common.columns = ['x', 'y']

  ### Join layout with alignment string ###
  graph_layout_common = graph_layout_common.join(node_data[['ref_align', 'read_align']])
  graph_layout_common = graph_layout_common.reset_index(drop=True)

  ### Make graph layout ###
  
  for i in range(len(data_info_list)):
    graph_layout = make_graph_layout(
      data_dir = data_dir_list[i],
      data_info = data_info_list[i],
      node_subst_type = node_subst_type,
      graph = graph_list[i],
      layout_type = graph_layout_type,
      graph_layout_common = graph_layout_common,
      separate_components = graph_layout_separate_components,
      reverse_complement = sequence_reverse_complement_list[i],
      universal_layout_x_scale_insertion = universal_layout_x_scale_insertion,
      universal_layout_y_scale_insertion = universal_layout_y_scale_insertion,
      universal_layout_x_scale_deletion = universal_layout_x_scale_deletion,
      universal_layout_y_scale_deletion = universal_layout_y_scale_deletion,
    )

    ### Plot edges and nodes ###
    edge_traces = []
    if edge_show:
      edge_traces = plot_graph_helper.make_edges_traces(
        data_info = data_info,
        graph = graph_list[i],
        layout = graph_layout,
        show_edge_labels = edge_labels_show,
        show_edge_types = edge_types_show,
        edge_width_scale = edge_width_scale,
        reverse_complement = sequence_reverse_complement_list[i],
      )

    node_traces = plot_graph_helper.make_point_traces(
      data_info = data_info_list[i],
      node_data = node_data_list[i],
      graph_layout = graph_layout,
      show_node_labels = node_labels_show,
      node_label_columns = node_label_columns,
      node_label_position = node_label_position,
      node_label_font_size = constants.GRAPH_LABEL_FONT_SIZE * font_size_scale,
      node_color_type = node_color_type_list[i],
      node_comparison_colors = node_comparison_colors,
      node_freq_ratio_range = node_freq_ratio_range,
      node_reference_outline_color = node_reference_outline_color,
      node_outline_color = node_outline_color,
      node_fill_color = node_fill_color,
      node_variation_type_colors = node_variation_type_colors,
      node_size_type = node_size_type,
      node_size_px_range = node_size_px_range,
      node_size_freq_range = node_size_freq_range,
      node_outline_width_scale = node_outline_width_scale,
      reverse_complement = sequence_reverse_complement_list[i],
    )

    for trace in (edge_traces + node_traces):
      figure_list[i].add_trace(trace)

    ### Format axes ###
    if not axes_show:
      figure_list[i].update_xaxes(visible = False)
      figure_list[i].update_yaxes(visible = False)
    
    if not np.isnan(plot_range_x[0]):
      figure_list[i].update_xaxes(range = plot_range_x)
    if not np.isnan(plot_range_y[0]):
      figure_list[i].update_yaxes(range = plot_range_y)

    ### Enable/disable legend ###
    figure_list[i].update_traces(showlegend = legend_plotly_show)

    ### Format for freq ratio colors ###
    if node_color_type_list[i] == 'freq_ratio_continuous':
      figure_list[i].update_traces(
        marker = {
          'colorscale': constants.get_freq_ratio_color_scale(
            node_comparison_colors[0],
            node_comparison_colors[1],
          ),
          'cmin': np.log(node_freq_ratio_range[0]),
          'cmax': np.log(node_freq_ratio_range[1]),
        }
      )

def get_figure_size_args(
  content_height_px,
  content_width_px,
  margin_top_px,
  margin_bottom_px,
  margin_left_px,
  margin_right_px,
):
  total_height_px = content_height_px + margin_top_px + margin_bottom_px
  total_width_px = content_width_px + margin_left_px + margin_right_px
  return {
    'content_height_px': content_height_px,
    'content_width_px': content_width_px,
    'total_height_px': total_height_px,
    'total_width_px': total_width_px,
  }

def make_graph_figure(
  data_dir_list,
  graph_layout_type = constants.GRAPH_LAYOUT_TYPE,
  graph_layout_separate_components = constants.GRAPH_LAYOUT_SEPARATE_COMPONENTS,
  sequence_reverse_complement_list = None,
  node_subst_type = constants.GRAPH_NODE_SUBST_TYPE,
  node_filter_variation_types = constants.GRAPH_NODE_FILTER_VARIATION_TYPES,
  node_filter_freq_range = constants.GRAPH_NODE_FILTER_FREQ_RANGE,
  node_filter_dist_range = constants.GRAPH_NODE_FILTER_DIST_RANGE,
  node_label_show = constants.GRAPH_NODE_LABEL_SHOW,
  node_label_columns = constants.GRAPH_NODE_LABEL_COLUMNS,
  node_label_position = constants.GRAPH_NODE_LABEL_POSITION,
  node_color_type_list = None,
  node_comparison_colors = constants.GRAPH_NODE_COMPARISON_COLORS,
  node_freq_ratio_range = constants.GRAPH_NODE_FREQ_RATIO_RANGE,
  node_reference_outline_color = constants.GRAPH_NODE_REFERENCE_OUTLINE_COLOR,
  node_outline_color = constants.GRAPH_NODE_OUTLINE_COLOR,
  node_fill_color = constants.GRAPH_NODE_FILL_COLOR,
  node_variation_type_colors = constants.GRAPH_NODE_VARIATION_TYPE_COLORS,
  node_size_type = constants.GRAPH_NODE_SIZE_TYPE,
  node_size_px_range = constants.GRAPH_NODE_SIZE_PX_RANGE,
  node_size_freq_range = constants.GRAPH_NODE_SIZE_FREQ_RANGE,
  edge_show = constants.GRAPH_EDGE_SHOW,
  edge_show_labels = constants.GRAPH_EDGE_SHOW_LABELS,
  edge_show_types = constants.GRAPH_EDGE_SHOW_TYPES,
  edge_width_scale = constants.GRAPH_EDGE_WIDTH_SCALE,
  graph_width_px = constants.GRAPH_WIDTH_PX,
  graph_height_px = constants.GRAPH_HEIGHT_PX,
  title_list = None,
  title_y_shift_px = constants.GRAPH_TITLE_Y_SHIFT_PX,
  legend_plotly_show = constants.GRAPH_LEGEND_PLOTLY_SHOW,
  legend_custom_show = constants.GRAPH_LEGEND_CUSTOM_SHOW,
  legend_custom_list = constants.GRAPH_LEGEND_CUSTOM_LIST,
  legend_x_shift_px = constants.GRAPH_LEGEND_X_SHIFT_PX,
  legend_y_shift_px = constants.GRAPH_LEGEND_X_SHIFT_PX,
  legend_vertical_space_px = constants.GRAPH_LEGEND_VERTICAL_SPACE_PX,
  legend_item_scale = constants.GRAPH_LEGEND_ITEM_SCALE,
  legend_colorbar_scale = constants.GRAPH_LEGEND_COLORBAR_SCALE,
  line_width_scale = constants.GRAPH_LINE_WIDTH_SCALE,
  node_outline_width_scale = constants.GRAPH_NODE_OUTLINE_WIDTH_SCALE,
  plot_range_x = constants.GRAPH_PLOT_RANGE_X,
  plot_range_y = constants.GRAPH_PLOT_RANGE_Y,
  graph_stats_show = constants.GRAPH_STATS_SHOW,
  graph_stats_x_frac = constants.GRAPH_STATS_X_FRAC,
  graph_stats_y_frac = constants.GRAPH_STATS_Y_FRAC,
  graph_stats_x_shift_px = constants.GRAPH_STATS_X_SHIFT_PX,
  graph_stats_y_shift_px = constants.GRAPH_STATS_Y_SHIFT_PX,
  graph_stats_x_anchor = constants.GRAPH_STATS_X_ANCHOR,
  graph_stats_y_anchor = constants.GRAPH_STATS_Y_ANCHOR,
  margin_top_px = constants.GRAPH_MARGIN_TOP_MIN_PX,
  margin_bottom_px = constants.GRAPH_MARGIN_BOTTOM_MIN_PX,
  margin_left_px = constants.GRAPH_MARGIN_LEFT_MIN_PX,
  margin_right_px = constants.GRAPH_MARGIN_RIGHT_MIN_PX,
  font_size_scale = constants.GRAPH_FONT_SIZE_SCALE,
  axes_show = constants.GRAPH_AXES_SHOW,
  universal_layout_x_scale_insertion = constants.GRAPH_UNIVERSAL_LAYOUT_X_SCALE_INSERTION,
  universal_layout_y_scale_insertion = constants.GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_INSERTION,
  universal_layout_x_scale_deletion = constants.GRAPH_UNIVERSAL_LAYOUT_X_SCALE_DELETION,
  universal_layout_y_scale_deletion = constants.GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_DELETION,
):
  if sequence_reverse_complement_list is None:
    sequence_reverse_complement_list = [constants.GRAPH_SEQUENCE_REVERSE_COMPLEMENT] * len(data_dir_list)
  if title_list is None:
    title_list = [constants.GRAPH_TITLE] * len(data_dir_list)
  data_info_list = [
    file_utils.read_tsv_dict(file_names.data_info(data_dir))
    for data_dir in data_dir_list
  ]

  if node_filter_variation_types is None:
    node_filter_variation_types = list(constants.VARIATION_TYPES)
  if node_subst_type == constants.SUBST_WITHOUT:
    node_filter_variation_types = [
      x for x in node_filter_variation_types
      if x not in ['substitution', 'mixed']
    ]

  if LAYOUT_PROPERTIES.get(graph_layout_type, {}).get('plot_range_x', None) is not None:
    plot_range_x_default = LAYOUT_PROPERTIES[graph_layout_type]['plot_range_x']
    plot_range_y_default = LAYOUT_PROPERTIES[graph_layout_type]['plot_range_y']
    if np.isnan(plot_range_x[0]):
      plot_range_x = plot_range_x_default
    if np.isnan(plot_range_y[0]):
      plot_range_y = plot_range_y_default

  edge_show = edge_show and LAYOUT_PROPERTIES[graph_layout_type]['has_edges']

  figure_list = [plotly.graph_objects.Figure() for _ in range(len(data_dir_list))]

  make_graph_figure_helper(
    figure_list = figure_list,
    data_dir_list = data_dir_list,
    data_info_list = data_info_list,
    node_subst_type = node_subst_type,
    node_filter_freq_range = node_filter_freq_range,
    node_filter_dist_range = node_filter_dist_range,
    edge_show = edge_show,
    edge_types_show = edge_show_types,
    edge_labels_show = edge_show_labels,
    edge_width_scale = edge_width_scale,
    graph_layout_type = graph_layout_type,
    graph_layout_separate_components = graph_layout_separate_components,
    sequence_reverse_complement_list = sequence_reverse_complement_list,
    node_labels_show = node_label_show,
    node_label_columns = node_label_columns,
    node_label_position = node_label_position,
    node_color_type_list = node_color_type_list,
    node_comparison_colors = node_comparison_colors,
    node_freq_ratio_range = node_freq_ratio_range,
    node_reference_outline_color = node_reference_outline_color,
    node_outline_color = node_outline_color,
    node_fill_color = node_fill_color,
    node_variation_type_colors = node_variation_type_colors,
    node_size_type = node_size_type,
    node_size_px_range = node_size_px_range,
    node_size_freq_range = node_size_freq_range,
    node_filter_variation_types = node_filter_variation_types,
    node_outline_width_scale = node_outline_width_scale,
    plot_range_x = plot_range_x,
    plot_range_y = plot_range_y,
    legend_plotly_show = legend_plotly_show,
    font_size_scale = font_size_scale,
    axes_show = axes_show,
    universal_layout_x_scale_insertion = universal_layout_x_scale_insertion,
    universal_layout_y_scale_insertion = universal_layout_y_scale_insertion,
    universal_layout_x_scale_deletion = universal_layout_x_scale_deletion,
    universal_layout_y_scale_deletion = universal_layout_y_scale_deletion,
  )

  figure_size_args = get_figure_size_args(
    content_height_px = graph_height_px,
    content_width_px = graph_width_px,
    margin_top_px = margin_top_px,
    margin_bottom_px = margin_bottom_px,
    margin_left_px = margin_left_px,
    margin_right_px = margin_right_px,
  )

  for i in range(len(data_info_list)):
    if graph_stats_show:
      make_graph_stats_ref_component(
        figure = figure_list[i],
        data_dir = data_dir_list[i],
        data_info = data_info_list[i],
        subst_type = node_subst_type,
        x = graph_stats_x_frac,
        y = graph_stats_y_frac,
        x_shift = -margin_left_px + graph_stats_x_shift_px,
        y_shift = graph_stats_y_shift_px,
        x_anchor = graph_stats_x_anchor,
        y_anchor = graph_stats_y_anchor,
        font_size_scale = font_size_scale,
      )

    figure_list[i].update_layout(
      width = figure_size_args['total_width_px'],
      height = figure_size_args['total_height_px'],

      font_color = 'black',

      legend_title_font_size = constants.GRAPH_LEGEND_TITLE_FONT_SIZE * font_size_scale,
      legend_font_size = constants.GRAPH_LEGEND_FONT_SIZE * font_size_scale,
      legend_itemsizing = 'constant',
      legend_itemwidth = constants.GRAPH_LEGEND_PLOTLY_ITEM_WIDTH_PX,
      legend_yanchor = 'top',
      legend_xanchor = 'left',

      margin_t = margin_top_px,
      margin_r = margin_right_px,
      margin_b = margin_bottom_px,
      margin_l = margin_left_px,
      margin_autoexpand = False,

      hovermode = 'closest',
      hoverlabel_font_size = constants.GRAPH_HOVER_LABEL_FONT_SIZE,
      hoverlabel_font_family = constants.GRAPH_HOVER_LABEL_FONT,
      hoverlabel_bgcolor = constants.GRAPH_HOVER_LABEL_BG_COLOR,

      plot_bgcolor = constants.GRAPH_BACKGROUND_COLOR,
    )

    if LAYOUT_PROPERTIES[graph_layout_type].get('preserve_aspect', False):
      figure_list[i].update_yaxes(scaleanchor='x', scaleratio=1)

    if title_list[i] is not None:
      figure_list[i].add_annotation(
        xref = 'paper',
        yref = 'paper',
        text = title_list[i],
        x = 0.5,
        y = 1,
        xanchor = 'center',
        yanchor = 'bottom',
        yshift = title_y_shift_px,
        font_size = constants.GRAPH_TITLE_FONT_SIZE * font_size_scale,
        showarrow = False,
      )

    if legend_custom_show:
      make_custom_legends(
        figure = figure_list[i],
        legend_list = legend_custom_list,
        content_height_px = figure_size_args['content_height_px'],
        data_info = data_info_list[i],
        node_reference_outline_color = node_reference_outline_color,
        node_outline_color = node_outline_color,
        node_fill_color = node_fill_color,
        node_comparison_colors = node_comparison_colors,
        node_variation_type_colors = node_variation_type_colors,
        node_filter_variation_types = node_filter_variation_types,
        node_size_freq_range = node_size_freq_range,
        node_size_px_range = node_size_px_range,
        node_freq_ratio_range = node_freq_ratio_range,
        edge_show_types = edge_show_types,
        legend_x_shift_px = legend_x_shift_px,
        legend_y_shift_px = legend_y_shift_px,
        legend_vertical_space_px = legend_vertical_space_px,
        legend_item_scale = legend_item_scale,
        legend_colorbar_scale = legend_colorbar_scale,
        font_size_scale = font_size_scale,
        line_width_scale = line_width_scale,
      )

  return figure_list

def parse_args():
  parser = argparse.ArgumentParser(
    description = (
      'Layout and plot variation-distance graphs.' +
      ' For more information about the layouts please see the publication FIXME.'
    ),
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
  )
  parser.add_argument(
    '--input',
    type = common_utils.check_dir,
    nargs = '+',
    help = (
      'Input data directories.' +
      ' All libraries specified here must have the same windowed reference sequence' +
      ' (i.e., the 20bp section of the reference around the DSB site should be the same' +
      ' in all libraries). All the libraries will be laid out using common x/y-coordinate' +
      ' assignments to the vertices.'
    ),
    required = True,
  )
  parser.add_argument(
    '--output',
    type = common_utils.check_file_output,
    nargs = '+',
    help = (
      'Output file(s). If not given no output will be written' +
      ' (useful only when using "--interactive").' +
      ' The file extension should be either ".html" for interactive HTML output' +
      ' or any image file format supported by the Plotly package.' +
      ' If not omitted, number of arguments should match' +
      ' the number of input directories.'
    ),
  )
  parser.add_argument(
    '--title',
    type = str,
    nargs = '+',
    help = (
    'If present, adds a title to the plot with this value.' +
    ' Number of arguments should match the number of input' +
    ' files.'
    ),
  )
  parser.add_argument(
    '--layout',
    choices = list(LAYOUT_PROPERTIES),
    default = 'universal_layout',
    help = 'The algorithm to use for laying out the graph.',
  )
  parser.add_argument(
    '--universal_layout_y_axis_x_pos',
    type = float,
    help = (
      'If present, shows a y-axis at the given x position' +
      ' showing the distances to the reference.' +
      ' Univeral layout only.'
    ),
  )
  parser.add_argument(
    '--universal_layout_x_axis_deletion_y_pos',
    type = float,
    help = (
      'If present, shows an x-axis for deletions at the given y position' +
      ' showing the approximate position of the deleted ranges.' +
      ' Univeral layout only.'
    ),
  )
  parser.add_argument(
    '--universal_layout_x_axis_deletion_label_type',
    type = str,
    choices = ['relative', 'absolute'],
    default = 'relative',
    help = (
      'The type of labeling to use for the universal layout deletion x-axis (if present).' +
      ' "relative" labels have 0 in the middle with negative/positive values on the left/right.' +
      ' "absolute" labels have 1 on the left and the length of the reference sequence on the right.'
    ),
  )
  parser.add_argument(
    '--universal_layout_x_axis_insertion_y_pos',
    type = float,
    help = (
      'If present, shows a x-axis for insertions at the given y position' +
      ' showing the first nucleotide of inserted sequences.' +
      ' Univeral layout only.'
    ),
  )
  parser.add_argument(
    '--universal_layout_y_axis_y_range',
    nargs = 2,
    type = float,
    default = [float('nan'), float('nan')],
    help = (
      'If showing an y-axis for the universal layout,' +
      ' the min and max y-position of the line.'
    ),
  )
  parser.add_argument(
    '--universal_layout_x_axis_x_range',
    nargs = 2,
    type = float,
    default = [float('nan'), float('nan')],
    help = (
      'If showing an x-axis for the universal layout,' +
      ' the min and max x-position of the line.'
    ),
  )
  parser.add_argument(
    '--universal_layout_y_axis_deletion_max_tick',
    type = int,
    help = (
      'If showing an y-axis for the universal layout,' +
      ' the max tick value for the deletion side.'
    ),
  )
  parser.add_argument(
    '--universal_layout_y_axis_insertion_max_tick',
    type = int,
    help = (
      'If showing an y-axis for the universal layout,' +
      ' the max tick value for the insertion side.'
    ),
  )
  parser.add_argument(
    '--universal_layout_x_scale_insertion',
    type = float,
    default = constants.GRAPH_UNIVERSAL_LAYOUT_X_SCALE_INSERTION,
    help = (
      'The factor for determining the scale on the universal layout insertion x-axis.' +
      ' Insertion vertex x-coordinates will be multiplied by this value.'
    ),
  )
  parser.add_argument(
    '--universal_layout_y_scale_insertion',
    type = float,
    default = constants.GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_INSERTION,
    help = (
      'The factor for determining the scale on the universal layout insertion y-axis.' +
      ' Insertion vertex y-coordinates will be multiplied by this value.'
    ),
  )
  parser.add_argument(
    '--universal_layout_x_scale_deletion',
    type = float,
    default = constants.GRAPH_UNIVERSAL_LAYOUT_X_SCALE_DELETION,
    help = (
      'The factor for determining the scale on the universal layout deletion x-axis.' +
      ' Deletion vertex x-coordinates will be multiplied by this value.'
    ),
  )
  parser.add_argument(
    '--universal_layout_y_scale_deletion',
    type = float,
    default = constants.GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_DELETION,
    help = (
      'The factor for determining the scale on the universal layout deletion y-axis.' +
      ' Deletion vertex y-coordinates will be multiplied by this value.'
    ),
  )
  parser.add_argument(
    '--subst_type',
    choices = constants.SUBST_TYPES,
    help = 'Whether to plot data with or without substitutions.',
    default = constants.SUBST_WITHOUT,
  )
  parser.add_argument(
    '--node_freq_range',
    nargs = 2,
    type = float,
    default = constants.GRAPH_NODE_SIZE_FREQ_RANGE,
    help = (
      'Min and max frequency to determine node size.' +
      ' Higher frequencies are clipped to this value.'
    ),
  )
  parser.add_argument(
    '--node_px_range',
    nargs = 2,
    type = float,
    default = constants.GRAPH_NODE_SIZE_PX_RANGE,
    help = (
      'Min and max node size in pixels as determined by the frequency.' +
      ' For no size scaling, set both values to the same number.'
    )
  )
  parser.add_argument(
    '--node_outline_scale',
    type = float,
    default = constants.GRAPH_NODE_OUTLINE_WIDTH_SCALE,
    help = (
      'How much to scale the node outline width (thickness).' +
      ' Values > 1 increase the width; values < 1 decrease the width.'
    ),
  )
  parser.add_argument(
    '--node_comparison_colors',
    type = str,
    default = constants.GRAPH_NODE_COMPARISON_COLORS,
    nargs = 2,
    help = (
      'The colors to use in the gradient when the node colors' +
      ' show the frequency ratio of two experiments.' +
      ' May be specified in hex (e.g., "#ff0000" for red) or with' +
      ' recognized keywords such as "red", "blue", "green".'
    ),
  )
  parser.add_argument(
    '--node_comparison_color_type',
    type = str,
    default = constants.GRAPH_NODE_COMPARISON_COLOR_TYPE,
    choices = ['continuous', 'discrete'],
    help = (
      'The type of color scheme to use for coloring nodes in a comparison graph.' +
      ' The "continuous" scheme uses a gradient of colors from the min to max ratio.' +
      ' The "discrete" schemes uses three colors to indicate that the ratio is <' +
      ' the min ratio, between the min and max ratio, or > the max ratio.' +
      ' The min and max ratios are determined by NODE_FREQ_RATIO_RANGE' +
      ' and the corresponding colors are determined by NODE_COMPARISON_COLORS.'
    ),
  )
  parser.add_argument(
    '--node_freq_ratio_range',
    type = float,
    default = constants.GRAPH_NODE_FREQ_RATIO_RANGE,
    nargs = 2,
    help = (
      ' The two frequencies used to determine node colors for comparison graphs.' +
      ' Also controls the range of ratios displayed on the frequency-ratio colorbar legend.' +
      ' Typically, the min value should be < 1 and the max value should be > 1.'
    ),
  )
  parser.add_argument(
    '--node_reference_outline_color',
    type = str,
    default = constants.GRAPH_NODE_REFERENCE_OUTLINE_COLOR,
    help = (
      'Color to make the reference node outline.' +
      ' May be specified in hex (e.g., "#ff0000" for red) or with' +
      ' recognized keywords such as "red", "blue", "green".'
    )
  )
  parser.add_argument(
    '--node_outline_color',
    type = str,
    default = constants.GRAPH_NODE_OUTLINE_COLOR,
    help = (
      'Color to make the default node outline.' +
      ' May be specified in hex (e.g., "#ff0000" for red) or with' +
      ' recognized keywords such as "red", "blue", "green".'
    )
  )
  parser.add_argument(
    '--node_fill_color',
    type = str,
    default = constants.GRAPH_NODE_FILL_COLOR,
    help = (
      'Color to make the default node fill.' +
      ' May be specified in hex (e.g., "#ff0000" for red) or with' +
      ' recognized keywords such as "red", "blue", "green".'
    )
  )
  parser.add_argument(
    '--variation_types',
    nargs = '+',
    default = constants.GRAPH_NODE_FILTER_VARIATION_TYPES,
    choices = list(constants.VARIATION_TYPES),
    help = (
      'The variation types that should be included in the graph.' +
      ' This should be a list of the types:' +
      ' "insertion", "deletion", "substitution", "mixed", "none".' +
      ' Default value: "insertion", "deletion", "none".' +
      ' "mixed" means nodes that have multiples variation types (e.g. insertions and substitutions).' +
      ' "none" means the reference node (no variations).' +
      ' May be specified in hex (e.g., "#ff0000" for red) or with' +
      ' recognized keywords such as "red", "blue", "green".'
    ),
  )
  parser.add_argument(
    '--variation_type_colors',
    type = str,
    nargs = 5,
    default = constants.GRAPH_NODE_VARIATION_TYPE_COLORS,
    help = (
      'The colors for the different variations types.' +
      ' They must be specified in the order: ' +
      ' '.join([x.upper() for x in constants.VARIATION_TYPES]) + '.' +
      ' MIXED is the color for nodes with multiple types of' +
      ' variations (e.g. insertions and substitutions); NONE is the color for' +
      ' the reference node (no variations).' +
      ' May be specified in hex (e.g., "#ff0000" for red) or with' +
      ' recognized keywords such as "red", "blue", "green".'
    ),
  )
  parser.add_argument(
    '--edge_show',
    type = bool,
    choices = [True, False],
    default = constants.GRAPH_EDGE_SHOW,
    help = 'Whether to show edges between nodes.',
  )
  parser.add_argument(
    '--edge_types',
    type = str,
    nargs = '+',
    choices = ['indel', 'substitution'],
    default = constants.GRAPH_EDGE_SHOW_TYPES,
    help = 'The edge types to show.',
  )
  parser.add_argument(
    '--edge_scale',
    type = float,
    default = constants.GRAPH_EDGE_WIDTH_SCALE,
    help = (
      'How much to scale the edges width (thickness).' +
      ' Values > 1 increase the width; values < 1 decrease the width.'
    ),
  )
  parser.add_argument(
    '--stats',
    action = 'store_true',
    help = 'If present, show graph summary statistics in the left margin.',
  )
  parser.add_argument(
    '--stats_x_shift_px',
    type = float,
    default = constants.GRAPH_STATS_X_SHIFT_PX,
    help = 'How much to shift the stats in the x direction (in pixels).',
  )
  parser.add_argument(
    '--stats_y_shift_px',
    type = float,
    default = constants.GRAPH_STATS_Y_SHIFT_PX,
    help = 'How much to shift the stats in the x direction (in pixels).',
  )
  parser.add_argument(
    '--width_px',
    type = int,
    default = constants.GRAPH_WIDTH_PX,
    help = 'The width of the plot in pixels.',
  )
  parser.add_argument(
    '--height_px',
    type = int,
    default = constants.GRAPH_HEIGHT_PX,
    help = 'The height of the plot in pixels.',
  )
  parser.add_argument(
    '--margin_top_px',
    type = int,
    default = constants.GRAPH_MARGIN_TOP_MIN_PX,
    help = 'The size of the top margin in pixels.',
  )
  parser.add_argument(
    '--margin_bottom_px',
    type = int,
    default = constants.GRAPH_MARGIN_BOTTOM_MIN_PX,
    help = 'The size of the bottom margin in pixels.',
  )
  parser.add_argument(
    '--margin_left_px',
    type = int,
    default = constants.GRAPH_MARGIN_LEFT_MIN_PX,
    help = 'The size of the left margin in pixels.',
  )
  parser.add_argument(
    '--margin_right_px',
    type = int,
    default = constants.GRAPH_MARGIN_RIGHT_MIN_PX,
    help = 'The size of the right margin in pixels.',
  )
  parser.add_argument(
    '--line_width_scale',
    type = float,
    default = constants.GRAPH_LINE_WIDTH_SCALE,
    help = (
      'How much to scale the line widths (aka thickness).' +
      ' Values > 1 increase the width; values < 1 decrease the width.'
    ),
  )
  parser.add_argument(
    '--font_size_scale',
    type = float,
    default = constants.GRAPH_FONT_SIZE_SCALE,
    help = (
      'How much to scale the font size.' +
      ' Values > 1 increase the font size; values < 1 decrease it.'
    ),
  )
  parser.add_argument(
    '--crop_x',
    nargs = 2,
    type = float,
    default = constants.GRAPH_CROP_X,
    help = (
      'Range of the horizontal dimension to crop.' +
      ' Specified with normalized coords in range [0, 1].'
    ),
  )
  parser.add_argument(
    '--crop_y',
    nargs = 2,
    type = float,
    default = constants.GRAPH_CROP_Y,
    help = (
      'Range of the vertical dimension to crop.' +
      ' Specified in normalized coords in range [0, 1].'
    ),
  )
  parser.add_argument(
    '--range_x',
    type = float,
    nargs = 2,
    default = constants.GRAPH_PLOT_RANGE_X,
    help = (
      'Range of x-axis for plotting.'
      'If not specified chosen automatically to either show all nodes or a preset value'
      ' for the layout.'
    ),
  )
  parser.add_argument(
    '--range_y',
    type = float,
    nargs = 2,
    default = constants.GRAPH_PLOT_RANGE_Y,
    help = (
      'Range of y-axis for plotting.'
      'If not specified chosen automatically to either show all nodes or a preset value'
      ' for the layout.'
    ),
  )
  parser.add_argument(
    '--legends',
    choices = constants.GRAPH_LEGENDS,
    nargs = '+',
    help = (
      'The types of legends to show.' +
      ' They are drawn from top to bottom on the right margin in the order specified.'
    ),
  )
  parser.add_argument(
    '--legend_x_shift_px',
    type = float,
    default = constants.GRAPH_LEGEND_X_SHIFT_PX,
    help = 'How much to shift the legends in the x direction (in pixels).',
  )
  parser.add_argument(
    '--legend_y_shift_px',
    type = float,
    default = constants.GRAPH_LEGEND_Y_SHIFT_PX,
    help = 'How much to shift the legends in the y direction (in pixels).',
  )
  parser.add_argument(
    '--legend_colorbar_scale',
    type = float,
    default = constants.GRAPH_LEGEND_COLORBAR_SCALE,
    help = 'How much to scale the colorbar legend (for frequency-ratio coloring).',
  )
  parser.add_argument(
    '--legend_spacing_px',
    type = int,
    default = constants.GRAPH_LEGEND_VERTICAL_SPACE_PX,
    help = 'Amount of vertical space in pixels between different legends.',
  )
  parser.add_argument(
    '--separate_components',
    action = 'store_true',
    help = 'If present, separate the connected components of the graph.',
  )
  parser.add_argument(
    '--interactive',
    action = 'store_true',
    help = (
      'If present opens the interactive version in a browser.'
      ' Uses the Ploty library figure.show() function to do so.'
    ),
  )
  parser.add_argument(
    '--reverse_complement',
    choices = ['0', '1'],
    nargs = '+',
    help = (
      'Whether to reverse complement the sequences in the data sets.' +
      ' If present, the number of values must be the same as the number of input directories.' +
      ' "1" mean reverse complement the sequence and "0" means do not.'
      ' Used for making a layout for data sets that have reference sequences'
      ' that are the reverse complements of each other.'
      ' If "1" also uses the reverse complement of sequences when determining the'
      ' display labels and hover text.' +
      ' This affects the universal layout and fractal layout.'
    ),
  )
  parser.add_argument(
    '--quiet',
    action = 'store_true',
    help = 'If present, do not print extra log messages.',
  )
  args = vars(parser.parse_args())

  if len(args['input']) == 0:
    raise Exception('No input arguments.')

  if args['output'] is None:
    args['output'] = [None] * len(args['input'])
  if len(args['output']) != len(args['input']):
    raise Exception(
      'Incorrect number of output args.' +
      f' Got {len(args["output"])}. Expected {len(args["input"])}.'
    )

  if args['title'] is None:
    args['title'] = [None] * len(args['input'])
  if len(args['title']) != len(args['input']):
    raise Exception(
      'Incorrect number of title args.' +
      f' Got {len(args["title"])}. Expected {len(args["input"])}.'
    )

  if args['reverse_complement'] is None:
    args['reverse_complement'] = ['0'] * len(args['input'])
  if len(args['reverse_complement']) != len(args['input']):
    raise Exception(
      'Incorrect number of reverse complement flags.'
      f'Got {len(args["reverse_complement"])}. Expected {len(args["input"])}.'
    )
  args['reverse_complement'] = [x == '1' for x in args['reverse_complement']]

  args['variation_type_colors'] = dict(zip(
    ['insertion', 'deletion', 'substitution', 'mixed', 'none'],
    args['variation_type_colors'],
  ))

  if args['node_comparison_color_type'] == 'continuous':
    args['node_comparison_color_type'] = 'freq_ratio_continuous'
  elif args['node_comparison_color_type'] == 'discrete':
    args['node_comparison_color_type'] = 'freq_ratio_discrete'
  else:
    raise Exception('Impossible.')

  return args

def main(
  input,
  output,
  layout,
  title,
  reverse_complement,
  subst_type,
  node_freq_range,
  node_px_range,
  node_comparison_colors,
  node_comparison_color_type,
  node_freq_ratio_range,
  node_reference_outline_color,
  node_outline_color,
  node_fill_color,
  variation_type_colors,
  variation_types,
  node_outline_scale,
  edge_show,
  edge_types,
  edge_scale,
  width_px,
  height_px,
  margin_top_px,
  margin_bottom_px,
  margin_left_px,
  margin_right_px,
  stats,
  stats_x_shift_px,
  stats_y_shift_px,
  separate_components,
  line_width_scale,
  font_size_scale,
  legends,
  legend_x_shift_px,
  legend_y_shift_px,
  legend_colorbar_scale,
  legend_spacing_px,
  range_x,
  range_y,
  universal_layout_y_axis_insertion_max_tick,
  universal_layout_y_axis_deletion_max_tick,
  universal_layout_y_axis_x_pos,
  universal_layout_y_axis_y_range,
  universal_layout_x_axis_deletion_y_pos,
  universal_layout_x_axis_deletion_label_type,
  universal_layout_x_axis_x_range,
  universal_layout_x_axis_insertion_y_pos,
  universal_layout_x_scale_insertion,
  universal_layout_y_scale_insertion,
  universal_layout_x_scale_deletion,
  universal_layout_y_scale_deletion,
  crop_x,
  crop_y,
  interactive,
  quiet,
):
  data_dir_list = input
  output_list = output
  for data_dir in data_dir_list:
    log_utils.log_input(data_dir)
  data_info_list = [
    file_utils.read_tsv_dict(file_names.data_info(data_dir))
    for data_dir in data_dir_list
  ]

  node_color_type_list = []
  for i in range(len(data_info_list)):
    if data_info_list[i]['format'] == 'individual':
      node_color_type_list.append(constants.GRAPH_NODE_COLOR_TYPE_INDIVIDUAL)
    elif data_info_list[i]['format'] == 'comparison':
      node_color_type_list.append(node_comparison_color_type)
    else:
      # impossible
      raise Exception('Unknown data format: ' + str(data_info_list[i]['format']))

  ref_seq_set = set()
  for i in range(len(data_info_list)):
    if reverse_complement[i]:
      ref_seq_set.add(kmer_utils.reverse_complement(data_info_list[i]['ref_seq_window']))
    else:
      ref_seq_set.add(data_info_list[i]['ref_seq_window'])
  if len(ref_seq_set) > 1:
    raise Exception(
      'Not all reference sequences are identical.' +
      ' Got ' + str(ref_seq_set) + '.'
    )

  figure_list = make_graph_figure(
    data_dir_list = data_dir_list,
    title_list = title,
    sequence_reverse_complement_list = reverse_complement,
    node_subst_type = subst_type,
    node_size_px_range = node_px_range,
    node_size_freq_range = node_freq_range,
    node_filter_variation_types = variation_types,
    node_outline_width_scale = node_outline_scale,
    node_color_type_list = node_color_type_list,
    node_comparison_colors = node_comparison_colors,
    node_freq_ratio_range = node_freq_ratio_range,
    node_reference_outline_color = node_reference_outline_color,
    node_outline_color = node_outline_color,
    node_fill_color = node_fill_color,
    node_variation_type_colors = variation_type_colors,
    graph_width_px = width_px,
    graph_height_px = height_px,
    graph_stats_show = stats,
    graph_stats_x_shift_px = stats_x_shift_px,
    graph_stats_y_shift_px = stats_y_shift_px,
    graph_layout_type = layout,
    graph_layout_separate_components = separate_components,
    margin_top_px = margin_top_px,
    margin_bottom_px = margin_bottom_px,
    margin_left_px = margin_left_px,
    margin_right_px = margin_right_px,
    edge_show = edge_show,
    edge_show_types = edge_types,
    edge_width_scale = edge_scale,
    legend_custom_show = legends is not None,
    legend_custom_list = legends,
    legend_plotly_show = False,
    legend_x_shift_px = legend_x_shift_px,
    legend_y_shift_px = legend_y_shift_px,
    legend_colorbar_scale = legend_colorbar_scale,
    legend_vertical_space_px = legend_spacing_px,
    line_width_scale = line_width_scale,
    font_size_scale = font_size_scale,
    plot_range_x = range_x,
    plot_range_y = range_y,
    universal_layout_x_scale_insertion = universal_layout_x_scale_insertion,
    universal_layout_y_scale_insertion = universal_layout_y_scale_insertion,
    universal_layout_x_scale_deletion = universal_layout_x_scale_deletion,
    universal_layout_y_scale_deletion = universal_layout_y_scale_deletion,
  )

  for i in range(len(figure_list)):
    if not quiet:
      x_min = np.inf
      x_max = -np.inf
      y_min = np.inf
      y_max = -np.inf
      for trace in figure_list[i].data:
        x_min = np.min([x_min] + [x for x in trace.x if x is not None])
        x_max = np.max([x_max] + [x for x in trace.x if x is not None])
        y_min = np.min([y_min] + [y for y in trace.y if y is not None])
        y_max = np.max([y_max] + [y for y in trace.y if y is not None])
      log_utils.log(f'Figure[{i}] x-range: {x_min} to {x_max}')
      log_utils.log(f'Figure[{i}] y-range: {y_min} to {y_max}')
    sequence_data = file_utils.read_tsv(
      file_names.sequence_data(data_dir_list[i], subst_type)
    )
    if universal_layout_y_axis_insertion_max_tick is None:
      try:
        max_tick_insertion = sequence_data.loc[
          sequence_data['variation_type'] == 'insertion',
          'dist_ref'
        ].max()
      except:
        # incase no insertions
        max_tick_insertion = 1
    else:
      max_tick_insertion = universal_layout_y_axis_insertion_max_tick
    
    if universal_layout_y_axis_deletion_max_tick is None:
      try:
        max_tick_deletion = sequence_data.loc[
          sequence_data['variation_type'] == 'deletion',
          'dist_ref'
        ].max()
      except:
        # incase no deletions
        max_tick_deletion = 1
    else:
      max_tick_deletion = universal_layout_y_axis_deletion_max_tick

    if layout == 'universal_layout':
      if universal_layout_y_axis_x_pos is not None:
        make_universal_layout_y_axis(
          figure = figure_list[i],
          x_pos = universal_layout_y_axis_x_pos,
          ref_length = len(data_info_list[i]['ref_seq_window']),
          cut_pos_ref = len(data_info_list[i]['ref_seq_window']) // 2,
          y_range = universal_layout_y_axis_y_range,
          max_tick_deletion = max_tick_deletion,
          max_tick_insertion = max_tick_insertion,
          x_scale_insertion = universal_layout_x_scale_insertion,
          y_scale_insertion = universal_layout_y_scale_insertion,
          x_scale_deletion = universal_layout_x_scale_deletion,
          y_scale_deletion = universal_layout_y_scale_deletion,
          font_size_scale = font_size_scale,
        )
      if universal_layout_x_axis_deletion_y_pos is not None:
        make_universal_layout_x_axis(
          figure = figure_list[i],
          var_type = 'deletion',
          y_pos = universal_layout_x_axis_deletion_y_pos,
          ref_length = len(data_info_list[i]['ref_seq_window']),
          cut_pos_ref = len(data_info_list[i]['ref_seq_window']) // 2,
          x_range = universal_layout_x_axis_x_range,
          deletion_label_type = universal_layout_x_axis_deletion_label_type,
          x_scale_insertion = universal_layout_x_scale_insertion,
          y_scale_insertion = universal_layout_y_scale_insertion,
          x_scale_deletion = universal_layout_x_scale_deletion,
          y_scale_deletion = universal_layout_y_scale_deletion,
          font_size_scale = font_size_scale,
        )
      if universal_layout_x_axis_insertion_y_pos is not None:
        make_universal_layout_x_axis(
          figure = figure_list[i],
          var_type = 'insertion',
          y_pos = universal_layout_x_axis_insertion_y_pos,
          ref_length = len(data_info_list[i]['ref_seq_window']),
          cut_pos_ref = len(data_info_list[i]['ref_seq_window']) // 2,
          x_range = universal_layout_x_axis_x_range,
          x_scale_insertion = universal_layout_x_scale_insertion,
          y_scale_insertion = universal_layout_y_scale_insertion,
          x_scale_deletion = universal_layout_x_scale_deletion,
          y_scale_deletion = universal_layout_y_scale_deletion,
          font_size_scale = font_size_scale,
        )

    if interactive:
      log_utils.log('Opening interactive version in browser.')
      figure_list[i].show()

    if output_list[i] is not None:
      ext = os.path.splitext(output_list[i])[1]
      file_utils.write_plotly(figure_list[i], output_list[i])
      log_utils.log_output(output_list[i])

      crop_x = tuple(crop_x)
      crop_y = tuple(crop_y)
      if (crop_x != (0, 1)) or (crop_y != (0, 1)):
        if ext == '.html':
          raise Exception('Cannot use crop setting with HTML output')

        image = PIL.Image.open(output_list[i])
        width_px, height_px = image.size

        left = crop_x[0] * width_px
        right = crop_x[1] * width_px
        top = crop_y[0] * height_px
        bottom = crop_y[1] * height_px

        image.crop((left, top, right, bottom)).save(output_list[i])
  log_utils.blank_line()

if __name__ == '__main__':
  main(**parse_args())
