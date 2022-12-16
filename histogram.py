import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'utils'))) # allow importing
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '6_plot_histogram'))) # allow importing

import file_names
import plot_histogram

if __name__ == '__main__':
  args = plot_histogram.parse_args()
  args['input'] = file_names.histogram_dir(args['input'])
  plot_histogram.main(**args)