import os

import argparse

import networkx as nx

import plotly.graph_objects

import pandas as pd
import numpy as np

import sklearn.decomposition

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
 'radial': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': True,
  },
 'universal': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': True,
  },
 'fractal': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': True,
    'plot_range_x': (-1, 1),
    'plot_range_y': (-1, 1),
  },
 'kamada': {
    'only_2d': False,
    'do_pca': True,
    'normalize': True,
    'has_edges': True,
    'plot_range_x': (0, 1),
    'plot_range_y': (0, 1),
  },
  'spectral': {
    'only_2d': False,
    'do_pca': True,
    'normalize': False,
    'has_edges': True, 
  },
  'spring': {
    'only_2d': False,
    'do_pca': True,
    'normalize': False,
    'has_edges': True, 
  },
  'shell': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': True, 
  },
  'spiral': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': True, 
  },
  'circular': {
    'only_2d': True,
    'do_pca': False,
    'has_edges': True, 
    'normalize': False,
  },
  'multipartite': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': True, 
  },
}

def parse_args():
  parser = argparse.ArgumentParser(
    description = (
      'Layout and plot variation-distance graphs.' +
      ' For more information about the layouts please see the README.'
    ),
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
  )
  io_group = parser.add_argument_group('Input/Output')
  framing_group = parser.add_argument_group('Framing')
  layout_group = parser.add_argument_group('Layouts')
  node_group = parser.add_argument_group('Node Properties')
  node_comparison_group = parser.add_argument_group('Node Comparison Properties')
  edge_group = parser.add_argument_group('Edge Properties')
  legend_group = parser.add_argument_group('Legends')
  universal_group = parser.add_argument_group('Universal Layout')
  misc_group = parser.add_argument_group('Miscellaneous')
  io_group.add_argument(
    '-i',
    type = str,
    nargs = '+',
    help = (
      'Input data directories.' +
      ' All libraries specified here must have the same windowed reference sequence' +
      ' (i.e., the 20bp section of the reference around the DSB site should be the same' +
      ' in all libraries). All the libraries will be laid out using combined x/y-coordinate' +
      ' assignments to the vertices. To plot a comparion graph between two libraries,' +
      f' specify both directories as a single argument separated by "{constants.COMPARISON_SEP}".'
    ),
    required = True,
    metavar = 'INPUT',
    dest = 'input',
  )
  io_group.add_argument(
    '-o',
    type = common_utils.check_file_output,
    nargs = '+',
    help = (
      'Output files. If omitted, no output will be written' +
      ' (useful only when using "--interactive").' +
      ' The file extension should be either ".html" for interactive HTML output' +
      ' or any image file format supported by the Plotly package.' +
      ' If present, number of arguments should match' +
      ' the number of input directories.'
    ),
    metavar = 'OUTPUT',
    dest = 'output',
  )
  io_group.add_argument(
    '--debug',
    type = common_utils.check_dir_output,
    help = (
      'If present, write debug files to the given directory.' +
      ' The debug files contain information about the nodes and edges' +
      ' of the graph.'
    )
  )
  io_group.add_argument(
    '--interactive',
    action = 'store_true',
    help = (
      'If present, opens the interactive version in a browser.'
      ' Uses the Ploty library "figure.show()" function to do so.'
    ),
  )
  layout_group.add_argument(
    '--layout',
    choices = list(LAYOUT_PROPERTIES),
    default = 'universal',
    help = 'The algorithm to use for laying out the graph.',
  )
  layout_group.add_argument(
    '--sep',
    action = 'store_true',
    help = 'If present, separate the connected components of the graph.',
  )
  universal_group.add_argument(
    '--ul_yax_x',
    type = float,
    help = (
      'If present, shows a y-axis at the given x position' +
      ' showing the distances to the reference.' +
      ' Universal layout only.'
      ' To determine appropriate values to set please see the console log, which' +
      ' shows the range of x-values of the nodes.' +
      ' Set to 0 to automatically determine the position.'
    ),
  )
  universal_group.add_argument(
    '--ul_yax_y',
    nargs = 2,
    type = float,
    default = [float('nan'), float('nan')],
    help = (
      'If showing an y-axis for the universal layout,' +
      ' the min and max y-position of the line.' +
      ' Note either or both of these values may be omitted or set to "nan", ' +
      ' and the "ul_yax_*_max_tick" parameters may be used instead.' +
      ' To determine appropriate values to set please see the console log, which' +
      ' shows the range of y-values of the nodes.'
    ),
  )
  universal_group.add_argument(
    '--ul_xax_del_y',
    type = float,
    help = (
      'If present, shows an x-axis for deletions at the given y position' +
      ' showing the approximate position of the deleted ranges.' +
      ' Universal layout only.'
      ' To determine appropriate values to set please see the console log, which' +
      ' shows the range of y-values of the nodes.' +
      ' Set to 0 to automatically determine the position.'
    ),
  )
  universal_group.add_argument(
    '--ul_xax_ins_y',
    type = float,
    help = (
      'If present, shows a x-axis for insertions at the given y position' +
      ' showing the first nucleotide of inserted sequences.' +
      ' Universal layout only.' +
      ' To determine appropriate values to set please see the console log,' +
      ' which shows the range of y-values of the nodes.' +
      ' Set to 0 to automatically determine the position.'
    ),
  )
  universal_group.add_argument(
    '--ul_xax_x',
    nargs = 2,
    type = float,
    default = [float('nan'), float('nan')],
    help = (
      'If showing an x-axis for the universal layout,' +
      ' the min and max x-position of the line.' +
      ' To determine appropriate values to set please see the console log, which' +
      ' shows the range of x-values of the nodes.'
    ),
  )
  universal_group.add_argument(
    '--ul_xax_del_label_type',
    type = str,
    choices = ['rel', 'abs'],
    default = 'rel',
    help = (
      'The type of labeling to use for the universal layout deletion x-axis (if present).' +
      ' "rel" = "relative" labels have 0 in the middle with negative/positive values on the left/right.' +
      ' "abs" = "absolute" labels have 1 on the left and the length of the reference sequence on the right.'
    ),
  )
  universal_group.add_argument(
    '--ul_yax_del_max_tick',
    type = int,
    help = (
      'If showing an y-axis for the universal layout,' +
      ' the max tick value for the deletion side.' +
      ' If omitted, will be caculated as the largest deletion' +
      ' in the data.'
    ),
  )
  universal_group.add_argument(
    '--ul_yax_ins_max_tick',
    type = int,
    help = (
      'If showing an y-axis for the universal layout,' +
      ' the max tick value for the insertion side.' +
      ' If omitted, will be caculated as the lart insertion' +
      ' in the data.'
    ),
  )
  universal_group.add_argument(
    '--ul_x_scale_ins',
    type = float,
    default = constants.GRAPH_UNIVERSAL_X_SCALE_INSERTION,
    help = (
      'The factor for determining the scale on the universal layout insertion x-axis.' +
      ' Insertion vertex x-coordinates will be multiplied by this value.'
    ),
  )
  universal_group.add_argument(
    '--ul_y_scale_ins',
    type = float,
    default = constants.GRAPH_UNIVERSAL_Y_SCALE_INSERTION,
    help = (
      'The factor for determining the scale on the universal layout insertion y-axis.' +
      ' Insertion vertex y-coordinates will be multiplied by this value.'
    ),
  )
  universal_group.add_argument(
    '--ul_x_scale_del',
    type = float,
    default = constants.GRAPH_UNIVERSAL_X_SCALE_DELETION,
    help = (
      'The factor for determining the scale on the universal layout deletion x-axis.' +
      ' Deletion vertex x-coordinates will be multiplied by this value.'
    ),
  )
  universal_group.add_argument(
    '--ul_y_scale_del',
    type = float,
    default = constants.GRAPH_UNIVERSAL_Y_SCALE_DELETION,
    help = (
      'The factor for determining the scale on the universal layout deletion y-axis.' +
      ' Deletion vertex y-coordinates will be multiplied by this value.'
    ),
  )
  node_group.add_argument(
    '--sub',
    type = int,
    choices = [0, 1],
    help = 'Whether to plot data without (0) or with (1) substitutions.',
    default = 0,
  )
  node_group.add_argument(
    '--size',
    nargs = 2,
    type = float,
    default = constants.GRAPH_NODE_SIZE_PX_RANGE,
    help = (
      'Min and max node size in pixels as determined by the frequency.' +
      ' For no size scaling, set both values to the same number.'
    )
  )
  node_group.add_argument(
    '--size_freq',
    nargs = 2,
    type = float,
    default = constants.GRAPH_NODE_SIZE_FREQ_RANGE,
    help = (
      'Min and max frequency to determine node size.' +
      ' Higher frequencies are clipped to this value.'
    ),
  )
  node_group.add_argument(
    '--filter_freq',
    nargs = 2,
    type = float,
    default = constants.GRAPH_NODE_FILTER_FREQ_RANGE,
    help = 'Min and max frequency to filter nodes by.'
  )
  node_group.add_argument(
    '--filter_dist',
    nargs = 2,
    type = float,
    default = constants.GRAPH_NODE_FILTER_DIST_RANGE,
    help = 'Min and max distance to filter nodes by.'
  )
  node_group.add_argument(
    '--outline_scale',
    type = float,
    default = constants.GRAPH_NODE_OUTLINE_WIDTH_SCALE,
    help = (
      'How much to scale the node outline width (thickness).' +
      ' Larger values increase the width and smaller values decrease the width.'
    ),
  )
  node_comparison_group.add_argument(
    '--ratio',
    type = float,
    default = constants.GRAPH_NODE_FREQ_RATIO_RANGE,
    nargs = 2,
    help = (
      ' The two frequency ratio used to determine node colors for comparison graphs.' +
      ' Also controls the range of ratios displayed on the frequency-ratio colorbar legend.' +
      ' Typically, the min value should be < 1 and the max value should be > 1.'
    ),
  )
  node_comparison_group.add_argument(
    '--ratio_colors',
    type = str,
    default = constants.GRAPH_NODE_RATIO_COLORS,
    nargs = 2,
    help = (
      'The colors to use in the gradient when the node colors' +
      ' show the frequency ratio of two experiments.' +
      ' May be specified in hex (e.g., "#ff0000" for red) or with' +
      ' recognized keywords such as "red", "blue", "green".'
    ),
  )
  node_comparison_group.add_argument(
    '--ratio_color_type',
    type = str,
    default = constants.GRAPH_NODE_RATIO_COLOR_TYPE,
    choices = ['cont', 'disc'],
    help = (
      'The type of color scheme to use for coloring nodes in a comparison graph.' +
      ' The "cont" scheme uses a continuous gradient of colors from the min to max ratio.' +
      ' The "disc" schemes uses three colors to indicate that the ratio is <' +
      ' the min ratio, between the min and max ratio, or > the max ratio.' +
      ' The min and max ratios are determined by NODE_FREQ_RATIO_RANGE' +
      ' and the corresponding colors are determined by COMPARISON_COLORS.'
    ),
  )
  node_group.add_argument(
    '--ref_outline_color',
    type = str,
    default = constants.GRAPH_NODE_REFERENCE_OUTLINE_COLOR,
    help = (
      'Color to make the reference node outline.' +
      ' May be specified in hex (e.g., "#ff0000" for red) or with' +
      ' recognized keywords such as "red", "blue", "green".'
    )
  )
  node_group.add_argument(
    '--outline_color',
    type = str,
    default = constants.GRAPH_NODE_OUTLINE_COLOR,
    help = (
      'Color to make the default node outline.' +
      ' May be specified in hex (e.g., "#ff0000" for red) or with' +
      ' recognized keywords such as "red", "blue", "green".'
    )
  )
  node_group.add_argument(
    '--fill_color',
    type = str,
    default = constants.GRAPH_NODE_FILL_COLOR,
    help = (
      'Color to make the default node fill.' +
      ' May be specified in hex (e.g., "#ff0000" for red) or with' +
      ' recognized keywords such as "red", "blue", "green".'
    )
  )
  node_group.add_argument(
    '--var_types',
    nargs = '+',
    default = constants.GRAPH_NODE_FILTER_VARIATION_TYPES,
    choices = constants.VARIATION_TYPES,
    help = (
      'The variation types that should be included in the graph.' +
      ' "ins" means nodes that have only insertions.' +
      ' "del" means nodes that have only deletions.' +
      ' "sub" means nodes that have only substitutions.' +
      ' "mix" means nodes that have multiples variation types (e.g., insertions and substitutions).' +
      ' "none" means the reference node (no variations).'
    ),
  )
  node_group.add_argument(
    '--var_type_colors',
    type = str,
    nargs = len(constants.GRAPH_NODE_VARIATION_TYPE_COLORS),
    default = constants.GRAPH_NODE_VARIATION_TYPE_COLORS,
    help = (
      'The colors for the different variations types.' +
      ' They must be specified in the order: ' +
      ', '.join([x.upper() for x in constants.GRAPH_NODE_VARIATION_TYPE_COLORS]) + '.' +
      ' See the documentation for the "--var_types" argument for more information.' +
      ' May be specified in hex (e.g., "#ff0000" for red) or with' +
      ' recognized keywords such as "red", "blue", "green".'
    ),
  )
  edge_group.add_argument(
    '--edge',
    type = int,
    choices = [0, 1],
    default = 1 if constants.GRAPH_EDGE_SHOW else 0,
    help = 'Whether to show edges between nodes. "0" = "no", "1" = "yes".',
  )
  edge_group.add_argument(
    '--edge_types',
    type = str,
    nargs = '+',
    choices = ['indel', 'sub'],
    default = constants.GRAPH_EDGE_SHOW_TYPES,
    help = 'The edge types to show.',
  )
  edge_group.add_argument(
    '--edge_scale',
    type = float,
    default = constants.GRAPH_EDGE_WIDTH_SCALE,
    help = (
      'How much to scale the edges width (thickness).' +
      ' Larger values increase the width ans smaller values decrease the width.'
    ),
  )
  framing_group.add_argument(
    '--width',
    type = int,
    default = constants.GRAPH_WIDTH_PX,
    help = 'The width of the plot in pixels.',
  )
  framing_group.add_argument(
    '--height',
    type = int,
    default = constants.GRAPH_HEIGHT_PX,
    help = 'The height of the plot in pixels.',
  )
  framing_group.add_argument(
    '--mar_t',
    type = int,
    default = constants.GRAPH_MARGIN_TOP_MIN_PX,
    help = 'The size of the top margin in pixels.',
  )
  framing_group.add_argument(
    '--mar_b',
    type = int,
    default = constants.GRAPH_MARGIN_BOTTOM_MIN_PX,
    help = 'The size of the bottom margin in pixels.',
  )
  framing_group.add_argument(
    '--mar_l',
    type = int,
    default = constants.GRAPH_MARGIN_LEFT_MIN_PX,
    help = 'The size of the left margin in pixels.',
  )
  framing_group.add_argument(
    '--mar_r',
    type = int,
    default = constants.GRAPH_MARGIN_RIGHT_MIN_PX,
    help = 'The size of the right margin in pixels.',
  )
  framing_group.add_argument(
    '--crop_x',
    nargs = 2,
    type = float,
    default = constants.GRAPH_CROP_X,
    help = (
      'Range of the horizontal dimension to crop.' +
      ' Specified with normalized coords in range [0, 1].' +
      ' May only be used with pixel image output formats (e.g., PNG).'
    ),
  )
  framing_group.add_argument(
    '--crop_y',
    nargs = 2,
    type = float,
    default = constants.GRAPH_CROP_Y,
    help = (
      'Range of the vertical dimension to crop.' +
      ' Specified in normalized coords in range [0, 1].' +
      ' May only be used with pixel image output formats (e.g., PNG).'
    ),
  )
  framing_group.add_argument(
    '--range_x',
    type = float,
    nargs = 2,
    default = constants.GRAPH_PLOT_RANGE_X,
    help = (
      'Range of x-axis for plotting.' +
      ' If omitted, chosen automatically to either show all nodes or a preset value for the layout.' +
      ' To determine appropriate values to set please see the console log,' +
      ' which shows the range of x-values of the nodes.'
    ),
  )
  framing_group.add_argument(
    '--range_y',
    type = float,
    nargs = 2,
    default = constants.GRAPH_PLOT_RANGE_Y,
    help = (
      'Range of y-axis for plotting.' +
      ' If omitted, chosen automatically to either show all nodes or a preset value for the layout.' +
      ' To determine appropriate values to set please see the console log,' +
      ' which shows the range of y-values of the nodes.'
    ),
  )
  legend_group.add_argument(
    '--legends',
    choices = constants.GRAPH_LEGENDS,
    nargs = '+',
    help = (
      'The types of legends to show.' +
      ' They are drawn from top to bottom on the right margin in the order specified.'
    ),
  )
  legend_group.add_argument(
    '--legend_x',
    type = float,
    default = constants.GRAPH_LEGEND_X_SHIFT_PX,
    help = 'How much to shift the legends in the x direction (in pixels).',
  )
  legend_group.add_argument(
    '--legend_y',
    type = float,
    default = constants.GRAPH_LEGEND_Y_SHIFT_PX,
    help = 'How much to shift the legends in the y direction (in pixels).',
  )
  legend_group.add_argument(
    '--colorbar_scale',
    type = float,
    default = constants.GRAPH_LEGEND_COLORBAR_SCALE,
    help = 'How much to scale the colorbar legend (for frequency-ratio coloring).',
  )
  legend_group.add_argument(
    '--legend_spacing',
    type = int,
    default = constants.GRAPH_LEGEND_VERTICAL_SPACE_PX,
    help = 'Amount of vertical space in pixels between different legends.',
  )
  misc_group.add_argument(
    '--title',
    type = str,
    nargs = '+',
    help = (
    'If present, adds a title to the plot with this value.' +
    ' Number of arguments should match the number of input' +
    ' files.'
    ),
  )
  misc_group.add_argument(
    '--line_scale',
    type = float,
    default = constants.GRAPH_LINE_WIDTH_SCALE,
    help = (
      'How much to scale the line widths (aka thickness).' +
      ' Larger values increase the width and smaller values decrease the width.'
    ),
  )
  misc_group.add_argument(
    '--font_scale',
    type = float,
    default = constants.GRAPH_FONT_SIZE_SCALE,
    help = (
      'How much to scale the font size.' +
      ' Larger values increase the font size and smaller values decrease it.'
    ),
  )
  misc_group.add_argument(
    '--rc',
    choices = [0, 1],
    type = int,
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
  misc_group.add_argument(
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
      ' Got {}. Expected {}.'.format(len(args['output']), len(args['input']))
    )

  if args['title'] is None:
    args['title'] = [None] * len(args['input'])
  if len(args['title']) != len(args['input']):
    raise Exception(
      'Incorrect number of title args.' +
      ' Got {}. Expected {}.'.format(len(args['title']), len(args['input']))
    )

  if args['rc'] is None:
    args['rc'] = [0] * len(args['input'])
  if len(args['rc']) != len(args['input']):
    raise Exception(
      'Incorrect number of reverse complement flags.'
      'Got {}. Expected {}.'.format(len(args['rc']), len(args['input']))
    )
  args['rc'] = [bool(x) for x in args['rc']]

  args['var_type_colors'] = dict(zip(
    constants.VARIATION_TYPES,
    args['var_type_colors'],
  ))

  if args['ratio_color_type'] == 'cont':
    args['ratio_color_type'] = 'ratio_cont'
  elif args['ratio_color_type'] == 'disc':
    args['ratio_color_type'] = 'ratio_disc'
  else:
    raise Exception('Impossible.')

  args['edge'] = bool(args['edge'])

  args['sub'] = constants.SUBST_TYPES[args['sub']]

  return args

def group_graph_nodes_by(graph, data_name):
  data = pd.Series(dict(graph.nodes(data_name)))
  return list(data.groupby(data).groups.values())

def make_radial_layout(graph):
  node_list = graph.nodes(data=True)

  bucket_dict = {'ins': {}, 'del': {}}

  ref_nodes = []

  for _, data in node_list:
    if data['dist_ref'] == 0:
      ref_nodes.append(data)
    else:
      dist_ref = data['dist_ref']
      var_type = data['var_type']
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
    if var_type == 'ins':
      y_sign = 1
      dist_scale = 2
      zig_zag_angle = 15
    elif var_type == 'del':
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

def make_fractal_layout(graph):
  node_list = graph.nodes(data=True)

  bucket_dict = {'ins': {}, 'del': {}}

  ref_nodes = []

  for _, data in node_list:
    if data['dist_ref'] == 0:
      ref_nodes.append(data)
    else:
      dist_ref = data['dist_ref']
      var_type = data['var_type']
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

        if var_type == 'ins':
          xy_dict[data['id']] = get_kmer_fractal_x_y(
            alignment_utils.get_insertion_str(
              ref_align,
              read_align,
            )
          )
        elif var_type == 'del':
          # We don't plot the deletion nodes in this layout.
          # We place them all at (0, 0) so they will hopefully be covered by the reference node.
          xy_dict[data['id']] = (0, 0)
        else:
          raise Exception('Impossible.')
  return xy_dict

# Determines how the rows are laid out in the insertion side
# of the universal layout.
def get_universal_insertion_row_spec(dist_ref):
  rows = 2 ** ((dist_ref - 1) // 2)
  cols = 4 ** dist_ref // rows
  row_space = 2 / rows
  return {
    'rows': rows,
    'cols': cols,
    'row_space': row_space,
  }

def get_pos_universal(
  ref_align,
  read_align,
  dist_ref,
  var_type,
  cut_pos_ref,
  x_scale_insertion = constants.GRAPH_UNIVERSAL_X_SCALE_INSERTION,
  y_scale_insertion = constants.GRAPH_UNIVERSAL_Y_SCALE_INSERTION,
  x_scale_deletion = constants.GRAPH_UNIVERSAL_X_SCALE_DELETION,
  y_scale_deletion = constants.GRAPH_UNIVERSAL_Y_SCALE_DELETION,
):
  if var_type == 'ins':
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
    
    row_spec_curr = get_universal_insertion_row_spec(dist_ref)
    num_rows = row_spec_curr['rows']
    num_cols = row_spec_curr['cols']
    row = kmer_index % num_rows
    col = kmer_index // num_rows
    prev_rows_offset = 0
    for dist_ref_prev in range(1, dist_ref):
      row_spec_prev = get_universal_insertion_row_spec(dist_ref_prev)
      prev_rows_offset += row_spec_prev['rows'] * row_spec_prev['row_space']

    curr_row_offset = row * row_spec_curr['row_space']
    y = 2 + prev_rows_offset + curr_row_offset
    x = (col / num_cols) - 0.5 * (1 - 1 / num_cols)
    return (
      2 * x * x_scale_insertion,
      0.5 * y * y_scale_insertion,
    )
  elif var_type == 'del':
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
  x_scale_insertion = constants.GRAPH_UNIVERSAL_X_SCALE_INSERTION,
  y_scale_insertion = constants.GRAPH_UNIVERSAL_Y_SCALE_INSERTION,
  x_scale_deletion = constants.GRAPH_UNIVERSAL_X_SCALE_DELETION,
  y_scale_deletion = constants.GRAPH_UNIVERSAL_Y_SCALE_DELETION,
):
  node_list = graph.nodes(data=True)

  bucket_dict = {'ins': {}, 'del': {}}

  ref_nodes = []

  for _, data in node_list:
    if data['dist_ref'] == 0:
      ref_nodes.append(data)
    else:
      dist_ref = data['dist_ref']
      var_type = data['var_type']
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

        xy_dict[data['id']] = get_pos_universal(
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

def make_universal_y_axis(
  figure,
  x_pos,
  ref_length,
  cut_pos_ref, # should be 1 based!
  max_tick_insertion = None,
  max_tick_deletion = None,
  y_range = [float('nan'), float('nan')],
  tick_length = 0.25,
  title_font_size = constants.GRAPH_AXES_TITLE_FONT_SIZE,
  tick_font_size = constants.GRAPH_AXES_TICK_FONT_SIZE,
  font_size_scale = constants.GRAPH_FONT_SIZE_SCALE,
  line_width_px = constants.GRAPH_UNIVERSAL_AXIS_LINE_WIDTH_PX,
  line_width_scale = constants.GRAPH_LINE_WIDTH_SCALE,
  x_scale_insertion = constants.GRAPH_UNIVERSAL_X_SCALE_INSERTION,
  y_scale_insertion = constants.GRAPH_UNIVERSAL_Y_SCALE_INSERTION,
  x_scale_deletion = constants.GRAPH_UNIVERSAL_X_SCALE_DELETION,
  y_scale_deletion = constants.GRAPH_UNIVERSAL_Y_SCALE_DELETION,
):
  tick_list = [{'dist_ref': 0, 'y_pos': 0}]
  dist_ref = 0
  finish_insertion = False
  finish_deletion = False
  while (not finish_insertion) or (not finish_deletion):
    dist_ref += 1
    for var_type in ['ins', 'del']:
      # Make a pair of ref_align and read_align that simulates
      # the correct number of insertions or deletions.
      ref_align = (
        ('A' * cut_pos_ref) +
        ('-' * dist_ref) +
        ('A' * (ref_length - cut_pos_ref))
      )
      read_align = 'A' * (ref_length + dist_ref)
      if var_type == 'del':
        ref_align, read_align = read_align, ref_align

      if var_type == 'ins':
        if (max_tick_insertion is not None) and (dist_ref > max_tick_insertion):
          finish_insertion = True
          continue
      elif var_type == 'del':
        if (max_tick_deletion is not None) and (dist_ref > max_tick_deletion):
          finish_deletion = True
          continue

      y_pos = get_pos_universal(
        ref_align = ref_align,
        read_align = read_align,
        dist_ref = dist_ref,
        var_type = var_type,
        cut_pos_ref = cut_pos_ref,
        x_scale_insertion = x_scale_insertion,
        y_scale_insertion = y_scale_insertion,
        x_scale_deletion = x_scale_deletion,
        y_scale_deletion = y_scale_deletion,
      )[1]

      if var_type == 'ins':
        if (not np.isnan(y_range[1])) and (y_pos > y_range[1]):
          finish_insertion = True
          continue
      elif var_type == 'del':
        if (not np.isnan(y_range[0])) and (y_pos < y_range[0]):
          finish_deletion = True
          continue

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
      line_width = line_width_px * line_width_scale,
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
    line_width = line_width_px * line_width_scale,
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

def make_universal_x_axis(
  figure,
  var_type,
  y_pos,
  ref_length,
  cut_pos_ref, # should be 1 based!
  x_range = [float('nan'), float('nan')],
  insertion_axis_type = 'bracket', # tick or bracket
  deletion_label_type = 'rel', # rel[ative] or abs[olute]
  deletion_tick_type = 'start', # start or midpoint
  base_tick_length = 0.25,
  title_font_size = constants.GRAPH_AXES_TITLE_FONT_SIZE,
  tick_font_size = constants.GRAPH_AXES_TICK_FONT_SIZE,
  font_size_scale = constants.GRAPH_FONT_SIZE_SCALE,
  line_width_px = constants.GRAPH_UNIVERSAL_AXIS_LINE_WIDTH_PX,
  line_width_scale = constants.GRAPH_LINE_WIDTH_SCALE,
  x_scale_insertion = constants.GRAPH_UNIVERSAL_X_SCALE_INSERTION,
  y_scale_insertion = constants.GRAPH_UNIVERSAL_Y_SCALE_INSERTION,
  x_scale_deletion = constants.GRAPH_UNIVERSAL_X_SCALE_DELETION,
  y_scale_deletion = constants.GRAPH_UNIVERSAL_Y_SCALE_DELETION,
):
  if insertion_axis_type not in ['tick', 'bracket']:
    raise Exception('Invalid insertion axis type: ' + str(insertion_axis_type))
  if deletion_label_type not in ['abs', 'rel']:
    raise Exception('Invalid deletion label type: ' + str(deletion_label_type))
  if var_type == 'ins':
    tick_font_size = 2 * tick_font_size
  tick_list = []
  if var_type == 'ins':
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
      x_pos = get_pos_universal(
        ref_align = fake_ref_align,
        read_align = fake_read_align,
        dist_ref = 1,
        var_type = 'ins',
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
  elif var_type == 'del':
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
      x_pos = get_pos_universal(
        fake_ref_align,
        fake_read_align,
        dist_ref,
        'del',
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
        line_width = line_width_px * line_width_scale,
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
    line_width = line_width_px * line_width_scale,
    line_color = 'black',
  )
  # axis title
  figure.add_annotation(
    x = x_range[0],
    # set bottom anchor exactly at axis
    y = y_pos,
    text = 'Deletion position' if var_type == 'del' else 'Insertion first letter',
    showarrow = False,
    font_size = title_font_size * font_size_scale,
    xanchor = 'left',
    yanchor = 'bottom',
    # shift bottom anchor slightly above axis
    yshift = 0.25 * title_font_size * font_size_scale,
  )

def make_graph_layout_single(
  data_info,
  graph,
  layout_type,
  universal_x_scale_insertion = constants.GRAPH_UNIVERSAL_X_SCALE_INSERTION,
  universal_y_scale_insertion = constants.GRAPH_UNIVERSAL_Y_SCALE_INSERTION,
  universal_x_scale_deletion = constants.GRAPH_UNIVERSAL_X_SCALE_DELETION,
  universal_y_scale_deletion = constants.GRAPH_UNIVERSAL_Y_SCALE_DELETION,
):
  if layout_type == 'radial':
    layout = make_radial_layout(graph)
  elif layout_type == 'universal':
    layout = make_universal_layout(
      graph,
      len(data_info['ref_seq_window']) // 2,
      x_scale_insertion = universal_x_scale_insertion,
      y_scale_insertion = universal_y_scale_insertion,
      x_scale_deletion = universal_x_scale_deletion,
      y_scale_deletion = universal_y_scale_deletion,
    )
  elif layout_type == 'fractal':
    layout = make_fractal_layout(graph)
  elif layout_type == 'kamada':
    layout = nx.kamada_kawai_layout(graph, dim = 2)
  elif layout_type == 'spectral':
    layout = nx.spectral_layout(graph, dim=2)
  elif layout_type == 'spring':
    layout = nx.spring_layout(graph, dim=2)
  elif layout_type == 'shell':
    layout = nx.shell_layout(
      graph,
      dim = 2,
      nlist = group_graph_nodes_by(graph, 'dist_ref'),
    )
  elif layout_type == 'spiral':
    layout = nx.spiral_layout(graph, dim=2)
  elif layout_type == 'circular':
    layout = nx.circular_layout(graph, dim=2)
  elif layout_type == 'multipartite':
    layout = nx.multipartite_layout(graph, subset_key='dist_ref')
  else:
    raise Exception('Unknown layout type: ' + str(layout_type))

  layout = pd.DataFrame.from_dict(layout, orient='index', columns=['x', 'y'])

  if (layout.shape[0] >= 2) and LAYOUT_PROPERTIES[layout_type]['do_pca']:
    layout = pd.DataFrame(
      data = (
        sklearn.decomposition.PCA(n_components=2)
        .fit_transform(layout.to_numpy())
      ),
      index = layout.index,
      columns = ['x', 'y'],
    )

  if LAYOUT_PROPERTIES[layout_type]['normalize']:
    dim_mins = {}
    scales = {}
    for i in ['x', 'y']:
      dim_min = np.min(layout[i])
      dim_max = np.max(layout[i])
      dim_mins[i] = dim_min
      if np.isclose(dim_min, dim_max):
        scales[i] = 0
      else:
        scales[i] = 1 / (dim_max - dim_min)
    if LAYOUT_PROPERTIES[layout_type].get('preserve_aspect', False):
      scales = [min(scales)] * layout.shape[1]
    for i in ['x', 'y']:
      layout[i] = (layout[i] - dim_mins[i]) * scales[i]
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
  data_info,
  graph,
  layout_type,
  graph_layout_precomputed = None,
  separate_components = True,
  universal_x_scale_insertion = constants.GRAPH_UNIVERSAL_X_SCALE_INSERTION,
  universal_y_scale_insertion = constants.GRAPH_UNIVERSAL_Y_SCALE_INSERTION,
  universal_x_scale_deletion = constants.GRAPH_UNIVERSAL_X_SCALE_DELETION,
  universal_y_scale_deletion = constants.GRAPH_UNIVERSAL_Y_SCALE_DELETION,
):
  if graph_layout_precomputed is not None:
    separate_components = False
    node_groups = None
    node_data = pd.DataFrame.from_dict(
      dict(graph.nodes(data=True)),
      orient = 'index',
    ).reset_index(drop=True)
    node_data = pd.merge(
      node_data[['id', 'ref_align', 'read_align']],
      graph_layout_precomputed[['ref_align', 'read_align', 'x', 'y']],
      on = ['ref_align', 'read_align'],
      how = 'inner',
    )[['id', 'x', 'y']]
    node_data = node_data.set_index('id', drop=True)
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

    layout_list = [
      make_graph_layout_single(
        data_info = data_info,
        graph = subgraph,
        layout_type = layout_type,
        universal_x_scale_insertion = universal_x_scale_insertion,
        universal_y_scale_insertion = universal_y_scale_insertion,
        universal_x_scale_deletion = universal_x_scale_deletion,
        universal_y_scale_deletion = universal_y_scale_deletion,
      )
      for subgraph in subgraph_list
    ]

  if separate_components:
    if ref_id is not None:
      grid_spec = make_grid_spec(len(node_groups), True)
    else:
      grid_spec = make_grid_spec(len(node_groups), False)
    for layout, panel in zip(layout_list, grid_spec):
      layout['x'] = layout['x'] * panel['width'] + panel['x']
      layout['y'] = layout['y'] * panel['height'] + panel['y']

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
  var_types,
  var_type_colors,
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
  for var_type in var_types:
    legend_items.append({
      'type': 'circle',
      'size': node_size_px,
      'text': constants.VARIATION_TYPES[var_type]['label'],
      'color': var_type_colors[var_type],
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
  node_var_type_colors,
  node_filter_var_types,
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
      (legend == 'var_type') and
      (data_info['format'] == constants.DATA_INDIVIDUAL)
    ):
      y_shift_curr_px = make_variation_color_legend(
        figure = figure,
        var_types = node_filter_var_types,
        var_type_colors = node_var_type_colors,
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
      (legend == 'ratio_disc') and
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
      (legend == 'ratio_cont') and
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

# Load all the data needed for making the graphs
def get_data(
  data_dir_list,
  node_subst_type,
  node_filter_var_types,
  node_filter_freq_range,
  node_filter_dist_range,
  reverse_complement_list,
  debug_dir = None,
):
  data_info_list = []
  node_data_list = []
  graph_list = []

  for data_dir, rc in zip(data_dir_list, reverse_complement_list):
    if constants.COMPARISON_SEP in data_dir:
      data_dir_1, data_dir_2 = data_dir.split(constants.COMPARISON_SEP)
      log_utils.log_input(data_dir_1)
      log_utils.log_input(data_dir_2)
      data_info_1 = file_utils.read_json(file_names.data_info(data_dir_1))
      data_info_2 = file_utils.read_json(file_names.data_info(data_dir_2))
      data_1 = file_utils.read_csv(file_names.window(data_dir_1, node_subst_type))
      data_2 = file_utils.read_csv(file_names.window(data_dir_2, node_subst_type))
      data_info, data = graph_utils.get_comparison_data(
        data_info_1 = data_info_1,
        data_info_2 = data_info_2,
        data_1 = data_1,
        data_2 = data_2,
      )
      node_data = graph_utils.get_node_data(data)
    else:
      log_utils.log_input(data_dir)
      data_info = file_utils.read_json(file_names.data_info(data_dir))
      data = file_utils.read_csv(file_names.window(data_dir, node_subst_type))
      node_data = graph_utils.get_node_data(data)

    node_data = node_data.loc[
      node_data['var_type'].isin(node_filter_var_types)
    ]
    node_data = node_data.loc[
      node_data['freq_mean'].between(
        node_filter_freq_range[0],
        node_filter_freq_range[1],
        inclusive = 'both',
      )
    ]
    node_data = node_data.loc[
      node_data['dist_ref'].between(
        node_filter_dist_range[0],
        node_filter_dist_range[1],
        inclusive = 'both',
      )
    ]
    if rc:
      node_data['ref_align'] = node_data['ref_align'].apply(kmer_utils.reverse_complement)
      node_data['read_align'] = node_data['read_align'].apply(kmer_utils.reverse_complement)
      for x in data_info.keys():
        if x.startswith('ref_seq'):
          if not pd.isna(data_info[x]):
            data_info[x] = kmer_utils.reverse_complement(data_info[x])
          else:
            data_info[x] = None

    edge_data = graph_utils.get_edge_data(node_data)

    graph = graph_utils.get_graph(node_data, edge_data)

    data_info_list.append(data_info)
    node_data_list.append(node_data)
    graph_list.append(graph)

    if debug_dir is not None:
      log_utils.log_output(debug_dir)
      file_utils.write_csv(
        pd.DataFrame(data_info, index=[0]),
        os.path.join(debug_dir, data_info['name'] + '_data_info.csv'),
      )
      file_utils.write_csv(
        node_data,
        os.path.join(debug_dir, data_info['name'] + '_node_data.csv'),
      )
      file_utils.write_csv(
        edge_data,
        os.path.join(debug_dir, data_info['name'] + '_edge_data.csv'),
      )

  # Make the combined data
  node_data_combined = pd.concat(
    [x[['ref_align', 'read_align', 'freq_mean']] for x in node_data_list],
    ignore_index = True,
    axis = 'index',
  )
  node_data_combined = (
    node_data_combined.groupby(['ref_align', 'read_align'])
    .max().reset_index() # max of freq_mean
  )
  node_data_combined = graph_utils.get_node_data(node_data_combined)
 
  edge_data_combined = graph_utils.get_edge_data(node_data_combined)
  graph_combined = graph_utils.get_graph(node_data_combined, edge_data_combined)

  if (debug_dir is not None) and (len(data_dir_list) > 1):
    file_utils.write_csv(
      node_data_combined,
      os.path.join(debug_dir, 'combined_node_data.csv'),
    )
    file_utils.write_csv(
      edge_data_combined,
      os.path.join(debug_dir, 'combined_edge_data.csv'),
    )

  return (
    data_info_list,
    node_data_list,
    graph_list,
    node_data_combined,
    graph_combined,
  )

def make_graph_layout_all(
  data_info_list,
  node_data_combined,
  graph_list,
  graph_combined,
  graph_layout_type,
  graph_layout_separate_components,
  universal_x_scale_insertion,
  universal_x_scale_deletion,
  universal_y_scale_insertion,
  universal_y_scale_deletion,
  debug_dir = None,
):
   # Make the combined graph layout
  graph_layout_combined = make_graph_layout(
    data_info = data_info_list[0],
    graph = graph_combined,
    layout_type = graph_layout_type,
    separate_components = graph_layout_separate_components,
    graph_layout_precomputed = None,
  )

  # Join layout with alignment string
  graph_layout_combined = graph_layout_combined.join(
    node_data_combined[['ref_align', 'read_align']]
  )[['ref_align', 'read_align', 'x', 'y']]
  graph_layout_combined = graph_layout_combined.reset_index(drop=True)

  # Make individual graph layouts
  graph_layout_list = []
  for i in range(len(data_info_list)):
    graph_layout_list.append(make_graph_layout(
      data_info = data_info_list[i],
      graph = graph_list[i],
      layout_type = graph_layout_type,
      graph_layout_precomputed = graph_layout_combined,
      separate_components = graph_layout_separate_components,
      universal_x_scale_insertion = universal_x_scale_insertion,
      universal_y_scale_insertion = universal_y_scale_insertion,
      universal_x_scale_deletion = universal_x_scale_deletion,
      universal_y_scale_deletion = universal_y_scale_deletion,
    ))
  
  if debug_dir is not None:
    if len(graph_layout_list) > 1:
      file_utils.write_csv(
        graph_layout_combined,
        os.path.join(debug_dir, 'combined_graph_layout.csv'),
      )
    for i in range(len(graph_layout_list)):
      file_utils.write_csv(
        graph_layout_list[i].reset_index(),
        os.path.join(debug_dir, data_info_list[i]['name'] + '_graph_layout.csv'),
      )

  return graph_layout_list, graph_layout_combined

def make_graph_figure_helper(
  figure_list,
  data_info_list,
  node_data_list,
  graph_list,
  graph_layout_list,
  edge_show = constants.GRAPH_EDGE_SHOW,
  edge_types_show = constants.GRAPH_EDGE_SHOW_TYPES,
  edge_labels_show = constants.GRAPH_EDGE_LABELS_SHOW,
  edge_width_scale = constants.GRAPH_EDGE_WIDTH_SCALE,
  node_labels_show = constants.GRAPH_NODE_LABEL_SHOW,
  node_label_columns = constants.GRAPH_NODE_LABEL_COLUMNS,
  node_label_position = constants.GRAPH_NODE_LABEL_POSITION,
  node_color_type_list = None, # Should be specified
  node_comparison_colors = constants.GRAPH_NODE_RATIO_COLORS,
  node_freq_ratio_range = constants.GRAPH_NODE_FREQ_RATIO_RANGE,
  node_reference_outline_color = constants.GRAPH_NODE_REFERENCE_OUTLINE_COLOR,
  node_outline_color = constants.GRAPH_NODE_OUTLINE_COLOR,
  node_fill_color = constants.GRAPH_NODE_FILL_COLOR,
  node_var_type_colors = constants.GRAPH_NODE_VARIATION_TYPE_COLORS,
  node_size_type = constants.GRAPH_NODE_SIZE_TYPE,
  node_size_px_range = constants.GRAPH_NODE_SIZE_PX_RANGE,
  node_size_freq_range = constants.GRAPH_NODE_SIZE_FREQ_RANGE,
  node_outline_width_scale = constants.GRAPH_NODE_OUTLINE_WIDTH_SCALE,
  plot_range_x = constants.GRAPH_PLOT_RANGE_X,
  plot_range_y = constants.GRAPH_PLOT_RANGE_Y,
  legend_plotly_show = constants.GRAPH_LEGEND_CUSTOM_SHOW,
  axes_show = constants.GRAPH_AXES_SHOW,
  font_size_scale = constants.GRAPH_FONT_SIZE_SCALE,
):
   # Make individual graph layouts
  for i in range(len(data_info_list)):
    ### Plot edges and nodes ###
    edge_traces = []
    if edge_show:
      edge_traces = plot_graph_helper.make_edges_traces(
        data_info = data_info_list[i],
        graph = graph_list[i],
        layout = graph_layout_list[i],
        show_edge_labels = edge_labels_show,
        show_edge_types = edge_types_show,
        edge_width_scale = edge_width_scale,
      )

    node_traces = plot_graph_helper.make_point_traces(
      data_info = data_info_list[i],
      node_data = node_data_list[i],
      graph_layout = graph_layout_list[i],
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
      node_var_type_colors = node_var_type_colors,
      node_size_type = node_size_type,
      node_size_px_range = node_size_px_range,
      node_size_freq_range = node_size_freq_range,
      node_outline_width_scale = node_outline_width_scale,
    )

    for trace in (edge_traces + node_traces):
      figure_list[i].add_trace(trace)

    # Format axes
    if not axes_show:
      figure_list[i].update_xaxes(visible=False)
      figure_list[i].update_yaxes(visible=False)
    
    if not np.isnan(plot_range_x[0]):
      figure_list[i].update_xaxes(range=plot_range_x)
    if not np.isnan(plot_range_y[0]):
      figure_list[i].update_yaxes(range=plot_range_y)

    # Enable/disable legend
    figure_list[i].update_traces(showlegend=legend_plotly_show)

    # Format for freq ratio colors
    if node_color_type_list[i] == 'ratio_cont':
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
  data_info_list,
  node_data_list,
  graph_list,
  graph_layout_list,
  graph_layout_type = constants.GRAPH_LAYOUT_TYPE,
  node_filter_var_types = constants.GRAPH_NODE_FILTER_VARIATION_TYPES,
  node_label_show = constants.GRAPH_NODE_LABEL_SHOW,
  node_label_columns = constants.GRAPH_NODE_LABEL_COLUMNS,
  node_label_position = constants.GRAPH_NODE_LABEL_POSITION,
  node_color_type_list = None, # Should be specified
  node_comparison_colors = constants.GRAPH_NODE_RATIO_COLORS,
  node_freq_ratio_range = constants.GRAPH_NODE_FREQ_RATIO_RANGE,
  node_reference_outline_color = constants.GRAPH_NODE_REFERENCE_OUTLINE_COLOR,
  node_outline_color = constants.GRAPH_NODE_OUTLINE_COLOR,
  node_fill_color = constants.GRAPH_NODE_FILL_COLOR,
  node_var_type_colors = constants.GRAPH_NODE_VARIATION_TYPE_COLORS,
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
  margin_top_px = constants.GRAPH_MARGIN_TOP_MIN_PX,
  margin_bottom_px = constants.GRAPH_MARGIN_BOTTOM_MIN_PX,
  margin_left_px = constants.GRAPH_MARGIN_LEFT_MIN_PX,
  margin_right_px = constants.GRAPH_MARGIN_RIGHT_MIN_PX,
  font_size_scale = constants.GRAPH_FONT_SIZE_SCALE,
  axes_show = constants.GRAPH_AXES_SHOW,
):
  figure_list = [plotly.graph_objects.Figure() for _ in range(len(data_info_list))]

  make_graph_figure_helper(
    figure_list = figure_list,
    data_info_list = data_info_list,
    node_data_list = node_data_list,
    graph_list = graph_list,
    graph_layout_list = graph_layout_list,
    edge_show = edge_show,
    edge_types_show = edge_show_types,
    edge_labels_show = edge_show_labels,
    edge_width_scale = edge_width_scale,
    node_labels_show = node_label_show,
    node_label_columns = node_label_columns,
    node_label_position = node_label_position,
    node_color_type_list = node_color_type_list,
    node_comparison_colors = node_comparison_colors,
    node_freq_ratio_range = node_freq_ratio_range,
    node_reference_outline_color = node_reference_outline_color,
    node_outline_color = node_outline_color,
    node_fill_color = node_fill_color,
    node_var_type_colors = node_var_type_colors,
    node_size_type = node_size_type,
    node_size_px_range = node_size_px_range,
    node_size_freq_range = node_size_freq_range,
    node_outline_width_scale = node_outline_width_scale,
    plot_range_x = plot_range_x,
    plot_range_y = plot_range_y,
    legend_plotly_show = legend_plotly_show,
    font_size_scale = font_size_scale,
    axes_show = axes_show,
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
        node_var_type_colors = node_var_type_colors,
        node_filter_var_types = node_filter_var_types,
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

def main(
  input,
  output,
  debug,
  layout,
  title,
  rc,
  sub,
  size,
  size_freq,
  ratio_colors,
  ratio_color_type,
  ratio,
  ref_outline_color,
  outline_color,
  fill_color,
  var_type_colors,
  var_types,
  filter_freq,
  filter_dist,
  outline_scale,
  edge,
  edge_types,
  edge_scale,
  width,
  height,
  mar_t,
  mar_b,
  mar_l,
  mar_r,
  sep,
  line_scale,
  font_scale,
  legends,
  legend_x,
  legend_y,
  colorbar_scale,
  legend_spacing,
  range_x,
  range_y,
  ul_yax_x,
  ul_yax_y,
  ul_xax_del_y,
  ul_xax_ins_y,
  ul_xax_x,
  ul_yax_ins_max_tick,
  ul_yax_del_max_tick,
  ul_xax_del_label_type,
  ul_x_scale_ins,
  ul_y_scale_ins,
  ul_x_scale_del,
  ul_y_scale_del,
  crop_x,
  crop_y,
  interactive,
  quiet,
):
  (
    data_info_list,
    node_data_list,
    graph_list,
    node_data_combined,
    graph_combined,
  ) = get_data(
    data_dir_list = input,
    node_subst_type = sub,
    node_filter_var_types = var_types,
    node_filter_freq_range = filter_freq,
    node_filter_dist_range = filter_dist,
    reverse_complement_list = rc,
    debug_dir = debug,
  )

  # Check that all the windowed reference sequences are identical
  ref_seq_set = set(data_info['ref_seq_window'] for data_info in data_info_list)
  if len(ref_seq_set) > 1:
    raise Exception(
      'Not all reference sequences are identical.' +
      ' Got ' + str(ref_seq_set) + '.'
    )

  # Make the layouts
  graph_layout_list, graph_layout_combined = make_graph_layout_all(
    data_info_list = data_info_list,
    node_data_combined = node_data_combined,
    graph_list = graph_list,
    graph_combined = graph_combined,
    graph_layout_type = layout,
    graph_layout_separate_components = sep,
    universal_x_scale_insertion = ul_x_scale_ins,
    universal_x_scale_deletion = ul_x_scale_del,
    universal_y_scale_insertion = ul_y_scale_ins,
    universal_y_scale_deletion = ul_y_scale_del,
    debug_dir = debug,
  )

  # Determine the x/y plot ranges
  if LAYOUT_PROPERTIES.get(layout, {}).get('plot_range_x', None) is not None:
    range_x_default = LAYOUT_PROPERTIES[layout]['plot_range_x']
    range_y_default = LAYOUT_PROPERTIES[layout]['plot_range_y']
    if np.isnan(range_x[0]):
      range_x = range_x_default
    if np.isnan(range_y[0]):
      range_y = range_y_default
  if layout == 'universal':
    padding = 0.1 # Extra padding for of universal axes
  else:
    padding = 0.05
  if np.isnan(range_x[0]):
    range_x = graph_layout_combined['x'].agg(['min', 'max']).to_numpy()
    range_x = (
      range_x +
      (np.array([-padding, padding]) * (range_x[1] - range_x[0]))
     ) # Add padding
  if np.isnan(range_y[0]):
    range_y = graph_layout_combined['y'].agg(['min', 'max']).to_numpy()
    range_y = (
      range_y +
      (np.array([-padding, padding]) * (range_y[1] - range_y[0]))
    ) # Add padding

  # Determine the node color type for each figure
  node_color_type_list = []
  for i in range(len(input)):
    if data_info_list[i]['format'] == 'individual':
      node_color_type_list.append(constants.GRAPH_NODE_COLOR_TYPE_INDIVIDUAL)
    elif data_info_list[i]['format'] == 'comparison':
      node_color_type_list.append(ratio_color_type)
    else:
      # Impossible
      raise Exception('Unknown data format: ' + str(data_info_list[i]['format']))
  
  # Set some default values
  if sub == 'withoutSubst':
    var_types = [
      x for x in var_types
      if x not in ['sub', 'mix']
    ]

  edge = edge and LAYOUT_PROPERTIES[layout]['has_edges']

  if title is None:
    title = [constants.GRAPH_TITLE] * len(input)

  # Make the figures
  figure_list = make_graph_figure(
    data_info_list = data_info_list,
    node_data_list = node_data_list,
    graph_list = graph_list,
    graph_layout_list = graph_layout_list,
    title_list = title,
    node_size_px_range = size,
    node_size_freq_range = size_freq,
    node_filter_var_types = var_types,
    node_outline_width_scale = outline_scale,
    node_color_type_list = node_color_type_list,
    node_comparison_colors = ratio_colors,
    node_freq_ratio_range = ratio,
    node_reference_outline_color = ref_outline_color,
    node_outline_color = outline_color,
    node_fill_color = fill_color,
    node_var_type_colors = var_type_colors,
    graph_width_px = width,
    graph_height_px = height,
    graph_layout_type = layout,
    margin_top_px = mar_t,
    margin_bottom_px = mar_b,
    margin_left_px = mar_l,
    margin_right_px = mar_r,
    edge_show = edge,
    edge_show_types = edge_types,
    edge_width_scale = edge_scale,
    legend_custom_show = legends is not None,
    legend_custom_list = legends,
    legend_plotly_show = False,
    legend_x_shift_px = legend_x,
    legend_y_shift_px = legend_y,
    legend_colorbar_scale = colorbar_scale,
    legend_vertical_space_px = legend_spacing,
    line_width_scale = line_scale,
    font_size_scale = font_scale,
    plot_range_x = range_x,
    plot_range_y = range_y,
  )

  # Log coordinate ranges
  x_min = graph_layout_combined['x'].min()
  x_max = graph_layout_combined['x'].max()
  y_min = graph_layout_combined['y'].min()
  y_max = graph_layout_combined['y'].max()
  if (not quiet) and (len(input) > 1):
    log_utils.log('Combined x-range: {} to {}'.format(x_min, x_max))
    log_utils.log('Combined y-range: {} to {}'.format(y_min, y_max))
  for i in range(len(input)):
    if not quiet:
      log_utils.log(
        'Figure[{}] x-range: {} to {}'.format(
          i,
          graph_layout_list[i]['x'].min(),
          graph_layout_list[i]['x'].max(),
        )
      )
      log_utils.log(
        'Figure[{}] y-range: {} to {}'.format(
          i,
          graph_layout_list[i]['y'].min(),
          graph_layout_list[i]['y'].max(),
        )
      )

  # Make the axes for the universal layout
  if layout == 'universal':
    for i in range(len(input)):
      if ul_yax_ins_max_tick is None:
        max_tick_insertion = node_data_combined['ins'].max()
      else:
        max_tick_insertion = ul_yax_ins_max_tick
      if ul_yax_del_max_tick is None:
        max_tick_deletion = node_data_combined['del'].max()
      else:
        max_tick_deletion = ul_yax_del_max_tick
      if ul_yax_x is not None:
        if ul_yax_x == 0:
          ul_yax_x = x_max + (x_max - x_min) * 0.05
        make_universal_y_axis(
          figure = figure_list[i],
          x_pos = ul_yax_x,
          ref_length = len(data_info_list[i]['ref_seq_window']),
          cut_pos_ref = len(data_info_list[i]['ref_seq_window']) // 2,
          y_range = ul_yax_y,
          max_tick_deletion = max_tick_deletion,
          max_tick_insertion = max_tick_insertion,
          x_scale_insertion = ul_x_scale_ins,
          y_scale_insertion = ul_y_scale_ins,
          x_scale_deletion = ul_x_scale_del,
          y_scale_deletion = ul_y_scale_del,
          font_size_scale = font_scale,
          line_width_px = constants.GRAPH_UNIVERSAL_AXIS_LINE_WIDTH_PX,
          line_width_scale = line_scale,
        )
      if ul_xax_del_y is not None:
        if ul_xax_del_y == 0:
          ul_xax_del_y = y_min - (y_max - y_min) * 0.05
        make_universal_x_axis(
          figure = figure_list[i],
          var_type = 'del',
          y_pos = ul_xax_del_y,
          ref_length = len(data_info_list[i]['ref_seq_window']),
          cut_pos_ref = len(data_info_list[i]['ref_seq_window']) // 2,
          x_range = ul_xax_x,
          deletion_label_type = ul_xax_del_label_type,
          x_scale_insertion = ul_x_scale_ins,
          y_scale_insertion = ul_y_scale_ins,
          x_scale_deletion = ul_x_scale_del,
          y_scale_deletion = ul_y_scale_del,
          font_size_scale = font_scale,
          line_width_px = constants.GRAPH_UNIVERSAL_AXIS_LINE_WIDTH_PX,
          line_width_scale = line_scale,
        )
      if ul_xax_ins_y is not None:
        if ul_xax_ins_y == 0:
          ul_xax_ins_y = y_max + (y_max - y_min) * 0.05
        make_universal_x_axis(
          figure = figure_list[i],
          var_type = 'ins',
          y_pos = ul_xax_ins_y,
          ref_length = len(data_info_list[i]['ref_seq_window']),
          cut_pos_ref = len(data_info_list[i]['ref_seq_window']) // 2,
          x_range = ul_xax_x,
          x_scale_insertion = ul_x_scale_ins,
          y_scale_insertion = ul_y_scale_ins,
          x_scale_deletion = ul_x_scale_del,
          y_scale_deletion = ul_y_scale_del,
          font_size_scale = font_scale,
          line_width_px = constants.GRAPH_UNIVERSAL_AXIS_LINE_WIDTH_PX,
          line_width_scale = line_scale,
        )

  # Do the final output
  for i in range(len(input)):
    if interactive:
      log_utils.log('Opening interactive version in browser.')
      figure_list[i].show()

    if output[i] is not None:
      ext = os.path.splitext(output[i])[1]
      file_utils.write_plotly(figure_list[i], output[i])
      log_utils.log_output(output[i])

      crop_x = tuple(crop_x)
      crop_y = tuple(crop_y)
      if (crop_x != (0, 1)) or (crop_y != (0, 1)):
        if ext == '.html':
          raise Exception('Cannot use crop setting with HTML output')

        log_utils.log(
          'Cropping image to x = ({}, {}) and y = ({}, {})'
          .format(*crop_x, *crop_y)
        )
        image = PIL.Image.open(output[i])
        width, height = image.size

        left = crop_x[0] * width
        right = crop_x[1] * width
        top = crop_y[0] * height
        bottom = crop_y[1] * height

        image.crop((left, top, right, bottom)).save(output[i])
  log_utils.blank_line()

if __name__ == '__main__':
  main(**parse_args())
