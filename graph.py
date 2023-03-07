import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'utils'))) # allow importing
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '5_plot_graph'))) # allow importing

import file_names
import plot_graph

if __name__ == '__main__':
  # DO MORE TESTING HERE!
  sys.argv = "graph.py --title Hello-World! --legend_spacing_px 800 --node_comparison_colors black white --node_reference_outline_color #FF0000 --node_outline_color #0000FF --node_fill_color #000000 --font_size_scale 4 --stats --input data/sense_R2 --output plot/graph/sense_R2.png --legend --margin_top_px 1000 --margin_bottom_px 1000 --margin_right_px 1000 --margin_left_px 2000 --layout kamada_layout --width 4000 --height 4000 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19".split(" ")
  # sys.argv = "graph.py --node_comparison_colors #CF191B #33A02C --subst_type withSubst --edge_types indel substitution --title Hello-World! --legend_spacing_px 800 --font_size_scale 4 --stats --input data/sense_db_R2 --output plot/graph/sense_db_R2.png --legend --margin_top_px 1000 --margin_bottom_px 1000 --margin_right_px 1000 --margin_left_px 2000 --layout shell_layout --width 4000 --height 4000 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19".split(" ")
  args = plot_graph.parse_args()
  args['input'] = file_names.graph_dir(args['input'])
  plot_graph.main(**args)