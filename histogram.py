import utils.file_names as file_names
import plot_histogram.plot_histogram as plot_histogram

if __name__ == '__main__':
  # Test code:
  # sys.argv = 'histogram.py --input data_output/db_R1 --output plot/histogram_test.png --color #bfbfbf --freq_range 1e-5 1e-2 --freq_scale linear --variation_type substitution --label_type relative'.split(' ')
  args = plot_histogram.parse_args()
  args['input'] = file_names.histogram_dir(args['input'])
  plot_histogram.main(**args)