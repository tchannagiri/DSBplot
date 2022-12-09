import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

import common_utils
import file_utils
import log_utils
import file_names
import library_constants

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
    file_names.variation_grouped(data_dir, library_constants.SUBST_WITH)
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
  freq_min,
  freq_max,
  freq_log,
  axis,
  show_title,
  label_type,
  tick_modulus = library_constants.HISTOGRAM_AXIS_TICK_MODULUS,
  axis_label_font_size = library_constants.HISTOGRAM_AXIS_LABEL_FONT_SIZE,
  axis_tick_font_size = library_constants.HISTOGRAM_AXIS_TICK_FONT_SIZE,
  font_size_scale = library_constants.HISTOGRAM_TITLE_FONT_SIZE,
  label_pad_px = 10,
  reverse_pos = False,
):

  ref_length = len(data_info['ref_seq_window'])
  ref_pos_labels = library_constants.get_position_labels(label_type, ref_length)

  if freq_log:
    freq_min_axis = np.log10(freq_min)
    freq_max_axis = np.log10(freq_max)
  else:
    freq_min_axis = freq_min
    freq_max_axis = freq_max

  data_sub_long = get_variation_data(
    data_dir = data_dir,
    data_info = data_info,
    variation_type = variation_type,
    format = 'long',
    reverse_pos = reverse_pos,
  )

  color = library_constants.VARIATION_TYPES[variation_type]['color_3d']

  x = data_sub_long.iloc[:, 0].to_numpy()
  y = data_sub_long.iloc[:, 1].to_numpy()
  z = data_sub_long.iloc[:, 2].to_numpy()

  if freq_log:
    z = np.log10(np.clip(z, freq_min, np.inf))
  z -= freq_min_axis

  # The clip with min .0001 is needed to prevent strange plotting
  # artifacts with plotting 0 height bars
  z = np.clip(z, 0.0001, np.inf)
  
  axis.bar3d(
    x = x,
    y = y,
    z = np.full_like(z, freq_min_axis),
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
  x_tick_list = [x for x in x_tick_list if (int(x[1]) % tick_modulus) == 0]
  axis.set_xlabel(
    x_label,
    labelpad = label_pad_px * font_size_scale,
    fontsize = axis_label_font_size * font_size_scale,
  )
  axis.set_xticks(
    ticks = [x[0] for x in x_tick_list],
    labels = [x[1] for x in x_tick_list],
  )
  y_ticks = [
    y for y in
    range(0, ref_length + 1)
    if (y % tick_modulus) == 0
  ]
  axis.set_yticks(
    ticks = y_ticks,
    labels = [str(y) for y in y_ticks],
  )
  axis.set_ylabel(
    'Variations',
    labelpad = label_pad_px * font_size_scale,
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
    labelpad = label_pad_px * 1.5 * font_size_scale,
    fontsize = axis_label_font_size * font_size_scale,
  )
  axis.set_zlim(freq_min_axis, freq_max_axis)
  axis.tick_params(
    axis = 'z',
    pad = 7 * font_size_scale,
  )
  if freq_log:
    z_ticks = list(range(round(freq_min_axis), round(freq_max_axis) + 1))
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

  if show_title:
    axis.set_title(
      library_constants.LABELS[data_info['control_type']] + ' ' + variation_type.capitalize(),
      fontsize = library_constants.HISTOGRAM_TITLE_FONT_SIZE * font_size_scale,
    )

def plot_histogram(
  file_out,
  data_dir,
  data_info,
  variation_type,
  freq_min,
  freq_max,
  freq_log,
  label_type,
  show_title = False,
  reverse_pos = False,
):
  if data_info['format'] != library_constants.DATA_INDIVIDUAL:
    raise Exception('Only applicable for individual data sets')

  figure, axis = plt.subplots(
    nrows = 1,
    ncols = 1,
    subplot_kw = {
      'projection': '3d',
      'proj_type': 'ortho',
    },
    **get_figure_args_pyplot(
      col_widths_px = [library_constants.HISTOGRAM_WIDTH_PX],
      row_heights_px = [library_constants.HISTOGRAM_HEIGHT_PX],
      col_space_px = 0,
      row_space_px = 0,
      margin_left_px = library_constants.HISTOGRAM_MARGIN_LEFT_PX,
      margin_right_px = library_constants.HISTOGRAM_MARGIN_RIGHT_PX,
      margin_top_px = library_constants.HISTOGRAM_MARGIN_TOP_PX,
      margin_bottom_px = library_constants.HISTOGRAM_MARGIN_BOTTOM_PX,
      dpi = library_constants.HISTOGRAM_DPI,
    ),
  )

  font_size_scale = library_constants.HISTOGRAM_FONT_SIZE_SCALE

  plot_histogram_impl(
    data_dir = data_dir,
    data_info = data_info,
    variation_type = variation_type,
    freq_min = freq_min,
    freq_max = freq_max,
    freq_log = freq_log,
    axis = axis,
    show_title = False,
    label_type = label_type,
    tick_modulus = library_constants.HISTOGRAM_AXIS_TICK_MODULUS,
    font_size_scale = font_size_scale,
    reverse_pos = reverse_pos,
  )

  if show_title:
    figure.suptitle(
      library_constants.get_data_label(data_info) + '\n' + variation_type, 
      fontsize = library_constants.HISTOGRAM_TITLE_FONT_SIZE * font_size_scale,
    )

  file_utils.write_pyplot(figure, file_out)
  log_utils.log(file_out)


def parse_args():
  parser = argparse.ArgumentParser(
    description = 'Plot 3d histograms showing variation type/position/frequency.'
  )
  parser.add_argument(
    '--input',
    type = common_utils.check_dir,
    help = 'Directory with the data files.',
    required = True,
  )
  parser.add_argument(
    '--output',
    type = common_utils.check_dir_output,
    help = 'Output directory.',
    required = True,
  )
  parser.add_argument(
    '--reverse_pos',
    action = 'store_true',
    help = (
      'Whether to reverse the x-axis positions.' +
      ' Useful if comparing reverse strand data with forward strand data.'
    )
  )
  parser.add_argument(
    '--label_type',
    choices = ['relative', 'absolute'],
    help = (
      'Whether to index the x-axis by absolute positions on the reference sequence'
      ' from 1 to ref_length, or relative positions from -ref_length / 2 to ref_length / 2'
      ' (skipping 0).'
    ),
    required = True
  )
  return parser.parse_args()

def main():
  args = parse_args()
  data_dir = args.input
  data_info = file_utils.read_tsv_dict(file_names.data_info(args.input))
  data_label = library_constants.get_data_label(data_info)
  for variation_type in ['substitution', 'insertion', 'deletion']:
    plot_histogram(
      file_out = os.path.join(args.output, file_names.histogram_3d(data_label, variation_type)),
      data_dir = data_dir,
      data_info = data_info,
      variation_type = variation_type,
      freq_min = library_constants.HISTOGRAM_FREQ_MIN,
      freq_max = library_constants.HISTOGRAM_FREQ_MAX,
      freq_log = True,
      label_type = args.label_type,
      show_title = False,
      reverse_pos = args.reverse_pos,
    )

if __name__ == '__main__':
  main()