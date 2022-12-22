import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'utils'))) # allow importing
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '5_plot_graph'))) # allow importing

import file_names
import plot_graph

if __name__ == '__main__':
  sys.argv = "graph.py --input data/db_R1 --output plot/graph/db_R1.png --legend --margin_top_px 1000 --margin_bottom_px 1000 --margin_right_px 1000 --layout universal --width 3600 --height 3600 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19".split(" ")
  args = plot_graph.parse_args()
  args['input'] = file_names.graph_dir(args['input'])
  plot_graph.main(**args)