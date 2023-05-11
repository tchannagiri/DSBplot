import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

import DSBplot.utils.common_utils as common_utils
import DSBplot.utils.file_utils as file_utils
import DSBplot.utils.log_utils as log_utils
import DSBplot.utils.file_names as file_names
import DSBplot.utils.constants as constants

def get_figure_args_pyplot(
  col_widths_px,
  row_heights_px,
  col_space_px,
  row_space_px,
  margin_left_px,
  margin_right_px,
  margin_bottom_px,
  margin_top_px,
  dpi,
):
  figure_width_px = (
    margin_left_px +
    margin_right_px +
    sum(col_widths_px) +
    (len(col_widths_px) - 1) * col_space_px
  )
  figure_height_px = (
    margin_bottom_px +
    margin_top_px +
    sum(row_heights_px) +
    (len(row_heights_px) - 1) * row_space_px
  )
  
  figure_width_inches = figure_width_px / dpi
  figure_height_inches = figure_height_px / dpi

  margin_left_frac = margin_left_px / figure_width_px
  margin_right_frac = 1 - margin_right_px / figure_width_px
  margin_bottom_frac = margin_bottom_px / figure_height_px
  margin_top_frac = 1 - margin_top_px / figure_height_px

  wspace_frac = col_space_px / np.mean(col_widths_px)
  hspace_frac = row_space_px / np.mean(row_heights_px)

  return {
    'figsize': (figure_width_inches, figure_height_inches),
    'dpi': dpi,
    'gridspec_kw': {
      'left': margin_left_frac,
      'right': margin_right_frac,
      'bottom': margin_bottom_frac,
      'top': margin_top_frac,
      'wspace': wspace_frac,
      'hspace': hspace_frac,
      'width_ratios': col_widths_px,
      'height_ratios': row_heights_px,
    }
  }

def get_variation_data(
  data_dir,
  data_info,
  variation_type,
  format,
  y_axis_column = 'dist_ref',
  reverse_pos = False,
):
  ref_length = len(data_info['ref_seq_window'])

  data_long = file_utils.read_tsv(
    file_names.variation_grouped(data_dir, constants.SUBST_WITH)
  )
  data_long = data_long.loc[data_long['variation_type'] == variation_type]
  if reverse_pos:
    data_long['variation_pos'] = ref_length + 1 - data_long['variation_pos']
  data_long = data_long.groupby(['variation_pos', y_axis_column])['freq_mean'].sum()
  data_long = data_long.reindex(pd.MultiIndex.from_product(
    [list(range(ref_length + 1)), list(range(ref_length + 1))],
    names = ['variation_pos', y_axis_column],
  ))
  data_long = data_long.reset_index()
  data_long = data_long.fillna(0)

  if format == 'long':
    return data_long
  elif format == 'grid':
    data_grid = data_long.copy()
    data_grid = data_grid.pivot(
      index = y_axis_column,
      columns = 'variation_pos',
      values = 'freq_mean',
    )
    data_grid = data_grid.rename_axis(
      index = None,
      columns = None,
    )
    data_grid.index = pd.MultiIndex.from_product(
      [
        [y_axis_column],
        data_grid.index,
      ],
    )
    data_grid.columns = pd.MultiIndex.from_product(
      [['variation_pos'], data_grid.columns],
    )
    return data_grid
  else:
    raise Exception('Invalid value for format: ' + str(format))

def plot_histogram_impl(
  data_dir,
  data_info,
  variation_type,
  freq_range,
  freq_log,
  axis,
  label_type,
  color,
  axis_tick_multiplier = constants.HISTOGRAM_AXIS_TICK_MULTIPLIER,
  axis_label_font_size = constants.HISTOGRAM_AXIS_LABEL_FONT_SIZE,
  axis_tick_font_size = constants.HISTOGRAM_AXIS_TICK_FONT_SIZE,
  font_size_scale = constants.HISTOGRAM_FONT_SIZE_SCALE,
  axis_label_pad_px = constants.HISTOGRAM_AXIS_LABEL_PAD_PX,
  reverse_pos = False,
):

  ref_length = len(data_info['ref_seq_window'])
  ref_pos_labels = constants.get_position_labels(label_type, ref_length)

  if freq_log:
    freq_range_axis = np.log10(freq_range)
  else:
    freq_range_axis = freq_range

  data_sub_long = get_variation_data(
    data_dir = data_dir,
    data_info = data_info,
    variation_type = variation_type,
    format = 'long',
    reverse_pos = reverse_pos,
  )

  x = data_sub_long.iloc[:, 0].to_numpy()
  y = data_sub_long.iloc[:, 1].to_numpy()
  z = data_sub_long.iloc[:, 2].to_numpy()

  if freq_log:
    z = np.log10(np.clip(z, freq_range[0], np.inf))
  z -= freq_range_axis[0]

  # The clip with min .0001 is needed to prevent strange plotting
  # artifacts with plotting 0 height bars
  z = np.clip(z, 0.0001, np.inf)
  
  axis.bar3d(
    x = x,
    y = y,
    z = np.full_like(z, freq_range_axis[0]),
    dx = 1,
    dy = 1,
    dz = z,
    shade = True,
    color = color,
  )

  if label_type == 'relative':
    x_label = 'Position (from DSB)'
  elif label_type == 'absolute':
    x_label = 'Position'
  else:
    raise Exception('Impossible')
  x_tick_list = list(zip(range(1, ref_length + 1), ref_pos_labels))
  x_tick_list = [x for x in x_tick_list if (int(x[1]) % axis_tick_multiplier) == 0]
  axis.set_xlabel(
    x_label,
    labelpad = axis_label_pad_px * font_size_scale,
    fontsize = axis_label_font_size * font_size_scale,
  )
  axis.set_xticks(
    ticks = [x[0] for x in x_tick_list],
    labels = [x[1] for x in x_tick_list],
  )
  y_ticks = [
    y for y in
    range(0, ref_length + 1)
    if (y % axis_tick_multiplier) == 0
  ]
  axis.set_yticks(
    ticks = y_ticks,
    labels = [str(y) for y in y_ticks],
  )
  axis.set_ylabel(
    'Variations',
    labelpad = axis_label_pad_px * font_size_scale,
    fontsize = axis_label_font_size * font_size_scale,
  )

  axis.tick_params(
    labelsize = axis_tick_font_size * font_size_scale,
  )

  axis.tick_params(
    axis = 'y',
    pad = 2 * font_size_scale,
  )

  axis.set_zlabel(
    'Frequency',
    labelpad = axis_label_pad_px * 1.5 * font_size_scale,
    fontsize = axis_label_font_size * font_size_scale,
  )
  axis.set_zlim(freq_range_axis[0], freq_range_axis[1])
  axis.tick_params(
    axis = 'z',
    pad = 7 * font_size_scale,
  )
  if freq_log:
    z_ticks = list(range(round(freq_range_axis[0]), round(freq_range_axis[1]) + 1))
    z_labels = []
    for tick in z_ticks:
      if tick == 0:
        z_labels.append(f'1')
      else:
        z_labels.append(f'$10^{{{tick}}}$')
    
    axis.set_zticks(
      ticks = z_ticks,
      labels = z_labels,
    )

def plot_histogram(
  file_out,
  data_dir,
  data_info,
  variation_type,
  freq_range,
  freq_log,
  label_type,
  color,
  reverse_pos = False,
  height_px = constants.HISTOGRAM_HEIGHT_PX,
  width_px = constants.HISTOGRAM_WIDTH_PX,
  margin_left_px = constants.HISTOGRAM_MARGIN_LEFT_PX,
  margin_right_px = constants.HISTOGRAM_MARGIN_RIGHT_PX,
  margin_top_px = constants.HISTOGRAM_MARGIN_TOP_PX,
  margin_bottom_px = constants.HISTOGRAM_MARGIN_BOTTOM_PX,
  axis_tick_multiplier = constants.HISTOGRAM_AXIS_TICK_MULTIPLIER,
  font_size_scale = constants.HISTOGRAM_FONT_SIZE_SCALE,
  title = constants.HISTOGRAM_TITLE,
):
  log_utils.log_input(data_dir)
  if data_info['format'] != 'individual':
    raise Exception('Only applicable for individual data sets')

  figure, axis = plt.subplots(
    nrows = 1,
    ncols = 1,
    subplot_kw = {
      'projection': '3d',
      'proj_type': 'ortho',
    },
    **get_figure_args_pyplot(
      col_widths_px = [width_px],
      row_heights_px = [height_px],
      col_space_px = 0,
      row_space_px = 0,
      margin_left_px = margin_left_px,
      margin_right_px = margin_right_px,
      margin_top_px = margin_top_px,
      margin_bottom_px = margin_bottom_px,
      dpi = constants.HISTOGRAM_DPI,
    ),
  )

  plot_histogram_impl(
    data_dir = data_dir,
    data_info = data_info,
    variation_type = variation_type,
    freq_range = freq_range,
    freq_log = freq_log,
    axis = axis,
    label_type = label_type,
    color = color,
    axis_tick_multiplier = axis_tick_multiplier,
    font_size_scale = font_size_scale,
    reverse_pos = reverse_pos,
  )

  if title is not None:
    figure.suptitle(
      title, 
      fontsize = constants.HISTOGRAM_TITLE_FONT_SIZE * font_size_scale,
    )

  file_utils.write_pyplot(figure, file_out)
  log_utils.log_output(file_out)
  log_utils.new_line()


def parse_args():
  parser = argparse.ArgumentParser(
    description = 'Plot 3d histograms showing variation type/position/frequency.',
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
  )
  parser.add_argument(
    '--input',
    type = common_utils.check_dir,
    help = 'Directory with the data files.',
    required = True,
  )
  parser.add_argument(
    '--output',
    type = common_utils.check_file_output,
    help = (
      'Output image file.' +
      ' Must have a file extension supported by the matplotlib package.'
    ),
    required = True,
  )
  parser.add_argument(
    '--variation_type',
    choices = ['substitution', 'insertion', 'deletion'],
    help = 'Which variation type to show in the histogram.',
    required = True,
  )
  parser.add_argument(
    '--reverse_pos',
    action = 'store_true',
    help = (
      'Whether to reverse the x-axis positions.' +
      ' Useful if comparing reverse strand data with forward strand data.'
    ),
  )
  parser.add_argument(
    '--label_type',
    choices = ['relative', 'absolute'],
    help = (
      'Whether to index the x-axis by "absolute" positions on the' +
      ' reference sequence from 1 to [reference length], or "relative" positions' +
      ' from -[reference length]/2 to [reference length]/2 (skipping 0).' +
      ' [reference length] refers to the length of the reference sequence' +
      ' after extracting the window around the DSB site.'
    ),
    required = True,
  )
  parser.add_argument(
    '--color',
    type = str,
    default = None,
    help = (
      'Color of the bar graph. If not specified,' +
      ' a default color based on VARIATION_TYPE will be chosen.'
    ),
  )
  parser.add_argument(
    '--freq_range',
    type = float,
    nargs = 2,
    default = constants.HISTOGRAM_FREQ_RANGE,
    help = 'Range of the z-axis frequency values to show.',
  )
  parser.add_argument(
    '--freq_scale',
    type = str,
    choices = ['linear', 'log'],
    default = constants.HISTOGRAM_FREQ_SCALE,
    help = 'Whether to use a linear or log scale for the z-axis frequency values.',
  )
  parser.add_argument(
    '--height',
    type = int,
    default = constants.HISTOGRAM_HEIGHT_PX,
    help = 'Height of the output image in pixels.',
  )
  parser.add_argument(
    '--width',
    type = int,
    default = constants.HISTOGRAM_WIDTH_PX,
    help = 'Width of the output image in pixels.',
  )
  parser.add_argument(
    '--margin_left',
    type = int,
    default = constants.HISTOGRAM_MARGIN_LEFT_PX,
    help = 'Left margin of the output image in pixels.',
  )
  parser.add_argument(
    '--margin_right',
    type = int,
    default = constants.HISTOGRAM_MARGIN_RIGHT_PX,
    help = 'Right margin of the output image in pixels.',
  )
  parser.add_argument(
    '--margin_top',
    type = int,
    default = constants.HISTOGRAM_MARGIN_TOP_PX,
    help = 'Top margin of the output image in pixels.',
  )
  parser.add_argument(
    '--margin_bottom',
    type = int,
    default = constants.HISTOGRAM_MARGIN_BOTTOM_PX,
    help = 'Bottom margin of the output image in pixels.',
  )
  parser.add_argument(
    '--tick_multiplier',
    type = int,
    default = constants.HISTOGRAM_AXIS_TICK_MULTIPLIER,
    help = (
      'Multiplier of the axis tick labels.' +
      ' Only multiples of this value will be shown on the axes.'
    ),
  )
  parser.add_argument(
    '--font_size_scale',
    type = float,
    default = constants.HISTOGRAM_FONT_SIZE_SCALE,
    help = (
      'Multiplier of the font sizes for title, axis labels, and axis ticks.' +
      ' Must adjust margins accordingly.'
    ),
  )
  parser.add_argument(
    '--title',
    type = str,
    help = 'Optional title of the plot. Must adjust top margins accordingly.',
  )
  args = vars(parser.parse_args())
  if args['color'] is None:
    args['color'] = constants.VARIATION_TYPES[args['variation_type']]['color_3d']
  return args

def main(
    input,
    output,
    variation_type,
    reverse_pos,
    label_type,
    color,
    freq_range,
    freq_scale,
    height = constants.HISTOGRAM_HEIGHT_PX,
    width = constants.HISTOGRAM_WIDTH_PX,
    margin_left = constants.HISTOGRAM_MARGIN_LEFT_PX,
    margin_right = constants.HISTOGRAM_MARGIN_RIGHT_PX,
    margin_top = constants.HISTOGRAM_MARGIN_TOP_PX,
    margin_bottom = constants.HISTOGRAM_MARGIN_BOTTOM_PX,
    tick_multiplier = constants.HISTOGRAM_AXIS_TICK_MULTIPLIER,
    font_size_scale = constants.HISTOGRAM_FONT_SIZE_SCALE,
    title = constants.HISTOGRAM_TITLE,
  ):
  data_dir = input
  data_info = file_utils.read_tsv_dict(file_names.data_info(input))
  plot_histogram(
    file_out = output,
    data_dir = data_dir,
    data_info = data_info,
    variation_type = variation_type,
    freq_range = freq_range,
    freq_log = freq_scale == 'log',
    label_type = label_type,
    color = color,
    reverse_pos = reverse_pos,
    height_px = height,
    width_px = width,
    margin_left_px = margin_left,
    margin_right_px = margin_right,
    margin_top_px = margin_top,
    margin_bottom_px = margin_bottom,
    axis_tick_multiplier = tick_multiplier,
    font_size_scale = font_size_scale,
    title = title,
  )

if __name__ == '__main__':
  main(**parse_args())