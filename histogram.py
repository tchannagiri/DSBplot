import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'utils'))) # allow importing
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '6_plot_histogram'))) # allow importing

import file_names
import plot_histogram

if __name__ == '__main__':
  # Test code:
  # sys.argv = 'histogram.py --input data_output/db_R1 --output plot/histogram_test.png --color #bfbfbf --freq_range 1e-5 1e-2 --freq_scale linear --variation_type substitution --label_type relative'.split(' ')
  args = plot_histogram.parse_args()
  args['input'] = file_names.histogram_dir(args['input'])
  plot_histogram.main(**args)