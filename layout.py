import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'utils'))) # allow importing
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '5_plot_graph'))) # allow importing

import file_names
import get_precomputed_layout

if __name__ == '__main__':
  args = get_precomputed_layout.parse_args()
  for i in range(len(args['input'])):
    args['input'][i] = file_names.graph_dir(args['input'][i])
  get_precomputed_layout.make_precomputed_layout(**args)