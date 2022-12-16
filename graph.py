import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'utils'))) # allow importing
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '5_plot_graph'))) # allow importing

import file_names
import plot_graph

if __name__ == '__main__':
  sys.argv = [
    '',
    '--input',
    'data/output',
    '--output',
    'plot/graph/output.png',
    '--layout',
    'universal',
    '--width',
    '2400',
    '--height',
    '1800',
    '--range_x',
    '-12',
    '13',
    '--range_y',
    '-23',
    '20',
    '--universal_layout_y_axis_x_pos',
    '12',
    '--universal_layout_y_axis_y_range',
    '-20.5',
    '18.5',
    '--universal_layout_y_axis_insertion_max_tick',
    '6',
    '--universal_layout_y_axis_deletion_max_tick',
    '19',
  ]
  args = plot_graph.parse_args()
  args['input'] = file_names.graph_dir(args['input'])
  plot_graph.main(**args)