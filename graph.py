import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'utils'))) # allow importing
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '5_plot_graph'))) # allow importing

import file_names
import plot_graph

if __name__ == '__main__':
  # Testing code:
  # sys.argv = "graph.py --variation_type_colors #FFFF00 #FF00FF #FFFF00 #FF0000 #00FF00 --title Sense BranchD pCMVD Sense_BranchD SenseR1 --legend_spacing_px 800 --node_comparison_colors black white --node_reference_outline_color #000000 --node_outline_color #0000FF --node_fill_color #ff0000 --font_size_scale 4 --stats --reverse_complement 0 0 0 0 1 --input data_output/sense_R2 data_output/db_R2 data_output/dcmv_R2 data_output/sense_db_R2 data_output/sense_R1 --output plot/graph/sense_R2.png plot/graph/db_R2.png plot/graph/dcmv_R2.png plot/graph/sense_db_R2.png plot/graph/sense_R1.png --legend --margin_top_px 1000 --margin_bottom_px 1000 --margin_right_px 1000 --margin_left_px 2000 --layout universal_layout --width 4000 --height 4000 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19".split(" ")
  # sys.argv = "graph.py --node_comparison_colors #CF191B #33A02C --node_freq_ratio_range 0.7 1.3 --node_comparison_color_type continuous --subst_type withSubst --edge_types indel substitution --title Hello-World! --legend_spacing_px 800 --font_size_scale 4 --stats --input data_output/sense_db_R2 --output plot/scratch_comparison.png --legend --margin_top_px 1000 --margin_bottom_px 1000 --margin_right_px 1600 --margin_left_px 2000 --layout shell_layout --width 4000 --height 4000 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19".split(" ")
  args = plot_graph.parse_args()
  args['input'] = [file_names.graph_dir(x) for x in args['input']]
  plot_graph.main(**args)