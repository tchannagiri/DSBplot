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
  var_type,
  format,
  y_axis_column = 'num_var',
  reverse_pos = False,
):
  ref_length = len(data_info['ref_seq_window'])

  data_long = file_utils.read_csv(
    file_names.variation(data_dir, 'withSubst')
  )
  data_long = data_long.loc[data_long['var_type'] == var_type]
  if reverse_pos:
    data_long['var_pos'] = ref_length + 1 - data_long['var_pos']
  data_long = data_long.groupby(['var_pos', y_axis_column])['freq_mean'].sum()
  data_long = data_long.reindex(pd.MultiIndex.from_product(
    [list(range(ref_length + 1)), list(range(ref_length + 1))],
    names = ['var_pos', y_axis_column],
  ))
  data_long = data_long.reset_index()
  data_long = data_long.fillna(0)

  if format == 'long':
    return data_long
  elif format == 'grid':
    data_grid = data_long.copy()
    data_grid = data_grid.pivot(
      index = y_axis_column,
      columns = 'var_pos',
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
      [['var_pos'], data_grid.columns],
    )
    return data_grid
  else:
    raise Exception('Invalid value for format: ' + str(format))

def plot_histogram_impl(
  data_dir,
  data_info,
  var_type,
  freq_range,
  freq_log,
  pyplot_axis,
  x_axis_type,
  color,
  axis_tick_multiple = constants.HISTOGRAM_AXIS_TICK_MULTIPLE,
  axis_label_font_size = constants.HISTOGRAM_AXIS_LABEL_FONT_SIZE,
  axis_tick_font_size = constants.HISTOGRAM_AXIS_TICK_FONT_SIZE,
  font_size_scale = constants.HISTOGRAM_FONT_SIZE_SCALE,
  axis_label_pad_px = constants.HISTOGRAM_AXIS_LABEL_PAD_PX,
  reverse_pos = False,
):

  ref_length = len(data_info['ref_seq_window'])
  ref_pos_labels = constants.get_position_labels(x_axis_type, ref_length)

  if freq_log:
    freq_range_axis = np.log10(freq_range)
  else:
    freq_range_axis = freq_range

  data_sub_long = get_variation_data(
    data_dir = data_dir,
    data_info = data_info,
    var_type = var_type,
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
  
  pyplot_axis.bar3d(
    x = x,
    y = y,
    z = np.full_like(z, freq_range_axis[0]),
    dx = 1,
    dy = 1,
    dz = z,
    shade = True,
    color = color,
  )

  if x_axis_type == 'rel':
    x_label = 'Position (from DSB)'
  elif x_axis_type == 'abs':
    x_label = 'Position'
  else:
    raise Exception('Impossible')
  x_tick_list = list(zip(range(1, ref_length + 1), ref_pos_labels))
  x_tick_list = [x for x in x_tick_list if (int(x[1]) % axis_tick_multiple) == 0]
  pyplot_axis.set_xlabel(
    x_label,
    labelpad = axis_label_pad_px * font_size_scale,
    fontsize = axis_label_font_size * font_size_scale,
  )
  pyplot_axis.set_xticks(
    ticks = [x[0] for x in x_tick_list],
    labels = [x[1] for x in x_tick_list],
  )
  y_ticks = [
    y for y in
    range(0, ref_length + 1)
    if (y % axis_tick_multiple) == 0
  ]
  pyplot_axis.set_yticks(
    ticks = y_ticks,
    labels = [str(y) for y in y_ticks],
  )
  pyplot_axis.set_ylabel(
    'Variations',
    labelpad = axis_label_pad_px * font_size_scale,
    fontsize = axis_label_font_size * font_size_scale,
  )

  pyplot_axis.tick_params(
    labelsize = axis_tick_font_size * font_size_scale,
  )

  pyplot_axis.tick_params(
    axis = 'y',
    pad = 2 * font_size_scale,
  )

  pyplot_axis.set_zlabel(
    'Frequency',
    labelpad = axis_label_pad_px * 1.5 * font_size_scale,
    fontsize = axis_label_font_size * font_size_scale,
  )
  pyplot_axis.set_zlim(freq_range_axis[0], freq_range_axis[1])
  pyplot_axis.tick_params(
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
    
    pyplot_axis.set_zticks(
      ticks = z_ticks,
      labels = z_labels,
    )

def plot_histogram(
  file_out,
  data_dir,
  data_info,
  var_type,
  freq_range,
  freq_log,
  x_axis_type,
  color,
  reverse_pos = False,
  height_px = constants.HISTOGRAM_HEIGHT_PX,
  width_px = constants.HISTOGRAM_WIDTH_PX,
  margin_left_px = constants.HISTOGRAM_MARGIN_LEFT_PX,
  margin_right_px = constants.HISTOGRAM_MARGIN_RIGHT_PX,
  margin_top_px = constants.HISTOGRAM_MARGIN_TOP_PX,
  margin_bottom_px = constants.HISTOGRAM_MARGIN_BOTTOM_PX,
  axis_tick_multiple = constants.HISTOGRAM_AXIS_TICK_MULTIPLE,
  font_size_scale = constants.HISTOGRAM_FONT_SIZE_SCALE,
  title = constants.HISTOGRAM_TITLE,
):
  log_utils.log_input(data_dir)
  if data_info['format'] != 'individual':
    raise Exception('Only applicable for individual data sets')

  pyplot_figure, pyplot_axis = plt.subplots(
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
    var_type = var_type,
    freq_range = freq_range,
    freq_log = freq_log,
    pyplot_axis = pyplot_axis,
    x_axis_type = x_axis_type,
    color = color,
    axis_tick_multiple = axis_tick_multiple,
    font_size_scale = font_size_scale,
    reverse_pos = reverse_pos,
  )

  if title is not None:
    pyplot_figure.suptitle(
      title, 
      fontsize = constants.HISTOGRAM_TITLE_FONT_SIZE * font_size_scale,
    )

  file_utils.write_pyplot(pyplot_figure, file_out)
  log_utils.log_output(file_out)
  log_utils.blank_line()


def parse_args():
  parser = argparse.ArgumentParser(
    description = 'Plot 3d histograms showing variation type/position/frequency.',
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
  )
  group_io = parser.add_argument_group('Input/output arguments')
  group_main = parser.add_argument_group('Main aesthestic arguments')
  group_framing = parser.add_argument_group('Framing arguments')
  group_axes = parser.add_argument_group('Axes arguments')
  group_io.add_argument(
    '-i',
    type = common_utils.check_dir,
    help = 'Directory with the data files.',
    required = True,
    metavar = 'INPUT',
    dest = 'input',
  )
  group_io.add_argument(
    '-o',
    type = common_utils.check_file_output,
    help = (
      'Output image file.' +
      ' Must have a file extension supported by the matplotlib package.'
    ),
    required = True,
    metavar = 'OUTPUT',
    dest = 'output',
  )
  group_main.add_argument(
    '--var',
    choices = constants.HISTOGRAM_VARIATION_TYPES,
    help = (
      'Which variation type to show in the histogram.' +
      '"ins" = insertions, "del" = deletions, "sub" = substitutions.'
    ),
    required = True,
    metavar = 'VAR_TYPE',
    dest = 'var_type',
  )
  group_main.add_argument(
    '--color',
    type = str,
    default = None,
    help = (
      'Color of the bar graph. If not specified,' +
      ' a default color based on VAR_TYPE will be chosen.' +
      ' Default colors are: "sub" = "{}", "ins" = "{}", "del" = "{}".'.format(
        constants.VARIATION_TYPES['sub']['color_3d'],
        constants.VARIATION_TYPES['ins']['color_3d'],
        constants.VARIATION_TYPES['del']['color_3d'],
      )
    ),
  )
  group_main.add_argument(
    '--title',
    type = str,
    help = 'Optional title of the plot. Must adjust top margins accordingly.',
  )
  group_main.add_argument(
    '--font_scale',
    type = float,
    default = constants.HISTOGRAM_FONT_SIZE_SCALE,
    help = (
      'Multiplier of the font sizes for title, axis labels, and axis ticks.' +
      ' Must adjust margins accordingly.'
    ),
    dest = 'font_size_scale',
    metavar = 'FONT_SIZE_SCALE',
  )
  group_axes.add_argument(
    '--xax',
    choices = ['rel', 'abs'],
    help = (
      'Whether to index the x-axis by "abs" = "absolute" positions on the' +
      ' reference sequence from 1 to [reference length], or "rel" = "relative" positions' +
      ' from -[reference length]/2 to [reference length]/2 (skipping 0).' +
      ' [reference length] refers to the length of the reference sequence' +
      ' after extracting the window around the DSB site.'
    ),
    required = True,
    metavar = 'X_AXIS_TYPE',
    dest = 'x_axis_type',
  )
  group_axes.add_argument(
    '--zax',
    type = str,
    choices = constants.HISTOGRAM_Z_AXIS_TYPES,
    default = constants.HISTOGRAM_Z_AXIS_TYPE,
    help = 'Whether to use a linear or log scale for the z-axis frequency values.',
    dest = 'z_axis_type',
    metavar = 'Z_AXIS_TYPE',
  )
  group_axes.add_argument(
    '--rev',
    choices = [0, 1],
    default = 0,
    type = int,
    help = (
      'Either reverse (1) the x-axis positions or keep them the same (0).' +
      ' Useful if comparing reverse strand data with forward strand data.'
    ),
    dest = 'reverse_pos',
    metavar = 'REVERSE_POS',
  )
  group_axes.add_argument(
    '--freq',
    type = float,
    nargs = 2,
    default = constants.HISTOGRAM_FREQ_RANGE,
    help = 'Range of the z-axis frequency values to show.',
    dest = 'freq_range',
    metavar = ('FREQ_MIN', 'FREQ_MAX'),
  )
  group_axes.add_argument(
    '--mult',
    type = int,
    default = constants.HISTOGRAM_AXIS_TICK_MULTIPLE,
    help = (
      'Multiple to use on the axis tick labels.' +
      ' Only multiples of this value will be shown on the axes.'
    ),
    dest = 'axis_tick_multiple',
    metavar = 'AXIS_TICK_MULTIPLE',
  )
  group_framing.add_argument(
    '--height',
    type = int,
    default = constants.HISTOGRAM_HEIGHT_PX,
    help = 'Height of the output image in pixels.',
    dest = 'height_px',
    metavar = 'HEIGHT_PX',
  )
  group_framing.add_argument(
    '--width',
    type = int,
    default = constants.HISTOGRAM_WIDTH_PX,
    help = 'Width of the output image in pixels.',
    dest = 'width_px',
    metavar = 'WIDTH_PX',
  )
  group_framing.add_argument(
    '--mar_l',
    type = int,
    default = constants.HISTOGRAM_MARGIN_LEFT_PX,
    help = 'Left margin of the output image in pixels.',
    dest = 'margin_left_px',
    metavar = 'MARGIN_LEFT_PX',
  )
  group_framing.add_argument(
    '--mar_r',
    type = int,
    default = constants.HISTOGRAM_MARGIN_RIGHT_PX,
    help = 'Right margin of the output image in pixels.',
    dest = 'margin_right_px',
    metavar = 'MARGIN_RIGHT_PX',
  )
  group_framing.add_argument(
    '--mar_t',
    type = int,
    default = constants.HISTOGRAM_MARGIN_TOP_PX,
    help = 'Top margin of the output image in pixels.',
    dest = 'margin_top_px',
    metavar = 'MARGIN_TOP_PX',
  )
  group_framing.add_argument(
    '--mar_b',
    type = int,
    default = constants.HISTOGRAM_MARGIN_BOTTOM_PX,
    help = 'Bottom margin of the output image in pixels.',
    dest = 'margin_bottom_px',
    metavar = 'MARGIN_BOTTOM_PX',
  )
  args = vars(parser.parse_args())
  if args['color'] is None:
    args['color'] = constants.VARIATION_TYPES[args['var_type']]['color_3d']
  args['reverse_pos'] = bool(args['reverse_pos'])
  return args

def main(
    input,
    output,
    var_type,
    reverse_pos,
    x_axis_type,
    color,
    freq_range,
    z_axis_type,
    height_px = constants.HISTOGRAM_HEIGHT_PX,
    width_px = constants.HISTOGRAM_WIDTH_PX,
    margin_left_px = constants.HISTOGRAM_MARGIN_LEFT_PX,
    margin_right_px = constants.HISTOGRAM_MARGIN_RIGHT_PX,
    margin_top_px = constants.HISTOGRAM_MARGIN_TOP_PX,
    margin_bottom_px = constants.HISTOGRAM_MARGIN_BOTTOM_PX,
    axis_tick_multiple = constants.HISTOGRAM_AXIS_TICK_MULTIPLE,
    font_size_scale = constants.HISTOGRAM_FONT_SIZE_SCALE,
    title = constants.HISTOGRAM_TITLE,
  ):
  data_dir = input
  data_info = file_utils.read_json(file_names.data_info(input))
  plot_histogram(
    file_out = output,
    data_dir = data_dir,
    data_info = data_info,
    var_type = var_type,
    freq_range = freq_range,
    freq_log = (z_axis_type == 'log'),
    x_axis_type = x_axis_type,
    color = color,
    reverse_pos = reverse_pos,
    height_px = height_px,
    width_px = width_px,
    margin_left_px = margin_left_px,
    margin_right_px = margin_right_px,
    margin_top_px = margin_top_px,
    margin_bottom_px = margin_bottom_px,
    axis_tick_multiple = axis_tick_multiple,
    font_size_scale = font_size_scale,
    title = title,
  )

if __name__ == '__main__':
  main(**parse_args())