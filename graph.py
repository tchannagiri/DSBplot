import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'utils'))) # allow importing
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '5_plot_graph'))) # allow importing

import file_names
import plot_graph

if __name__ == '__main__':
  args = plot_graph.parse_args()
  args['input'] = file_names.graph_dir(args['input'])
  plot_graph.main(**args)