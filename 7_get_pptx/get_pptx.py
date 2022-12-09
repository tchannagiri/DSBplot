import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../2_graph_processing/'))) # allow importing the graphs dir
import argparse

import numpy as np

import common_utils
import file_utils
import file_names
import log_utils
import library_constants
import get_pptx_helpers
import get_pptx_legend

import pptx
import pptx.util
import pptx.enum.shapes
import pptx.enum.text
import pptx.dml.color
import pptx.dml.effect
import PIL

PPTX_TEMPLATE_FILE = os.path.join(os.path.dirname(__file__), 'template.pptx') # make this an arg!

FORMAT_LEGENDS = {
  library_constants.DATA_COMPARISON: [
    'node_size',
    'node_outline',
    'edge_type',
  ],
  library_constants.DATA_INDIVIDUAL: [
    'node_size',
    'node_outline',
    'edge_type',
    'variation_type',
  ],
  'both': [
    'node_size',
    'node_outline',
    'edge_type',
    'variation_type',
  ],
}

FREQ_RATIO_LEGENDS = {
  '1DSB': [
    'freq_ratio_sense_branch',
    'freq_ratio_sense_cmv',
  ],
  '2DSB': [
    'freq_ratio_sense_branch',
    'freq_ratio_sense_cmv',
  ],
  '2DSBanti': [
    'freq_ratio_antisense_splicing',
  ],
}

LEGENDS = {
  'node_size': {'type': 'node_size'},
  'node_outline': {'type': 'node_outline'},
  'edge_type': {'type': 'edge_type'},
  'variation_type': {'type': 'variation_type'},
  'node_type': {'type': 'node_type'},
  'freq_ratio_sense_branch': {
   'type': 'freq_ratio',
   'construct_1': 'sense',
   'construct_2': 'branch',
   'color_bar_file': os.path.join(file_names.IMAGE_DIR, 'freq_ratio_sense_branch.png'),
  },
  'freq_ratio_sense_cmv': {
   'type': 'freq_ratio',
   'construct_1': 'sense',
   'construct_2': 'cmv',
   'color_bar_file': os.path.join(file_names.IMAGE_DIR, 'freq_ratio_sense_cmv.png'),
  },
  'freq_ratio_antisense_splicing': {
   'type': 'freq_ratio',
   'construct_1': 'antisense',
   'construct_2': 'splicing',
   'color_bar_file': os.path.join(file_names.IMAGE_DIR, 'freq_ratio_antisense_splicing.png'),
  },
}

### Font sizes ###
IMAGE_LABEL_FONT_SIZE_PT = 8
TITLE_FONT_SIZE_PT = 18
IMAGE_LABEL_WIDTH_PT = 60
IMAGE_LABEL_HEIGHT_PT = 20
MARGIN_CORNER_LABEL_FONT_SIZE_PT = 12
MARGIN_TOP_LABEL_FONT_SIZE_PT = 9
MARGIN_LEFT_LABEL_FONT_SIZE_PT = 8
LEGEND_TITLE_FONT_SIZE_PT = 10
LEGEND_LABEL_FONT_SIZE_PT = 8

### Height/width of elements ###
TITLE_HEIGHT_PT = 30
MARGIN_TOP_HEIGHT_PT = 20
MARGIN_LEFT_WIDTH_PT = 60
MARGIN_LEFT_SPILL_OVER_PT = 10
MARGIN_RIGHT_SPILL_OVER_PT = 10
MARGIN_RIGHT_WIDTH_PT = 100

LEGEND_TITLE_HEIGHT_PT = 20
LEGEND_ITEM_HEIGHT_PT = 20
LEGEND_SIZE_NODE_OUTLINE_WIDTH_PT = 0.5
LEGEND_NODE_SIZE_PT = 10
LEGEND_NODE_OUTLINE_WIDTH_PT = 0.5
LEGEND_EDGE_LINE_SIZE_PT = 10
LEGEND_EDGE_LINE_WIDTH_PT = 1
LEGEND_FREQ_RATIO_COLOR_BAR_WIDTH_PT = 10
LEGEND_FREQ_RATIO_COLOR_BAR_HEIGHT_PT = 75

LEGEND_HEIGHT_SPACING_PT = 25

HEIGHT_SPACING_PT = 10
WIDTH_SPACING_PT = 5
MULTIPLE_GRIDS_SPACING_PT = 40
CONTENT_LEGEND_SPACING_PT = 20

def get_slide(
  prs,
  image_grid_list,
  total_width_frac_list = None,
  title = None,
  title_height_pt = TITLE_HEIGHT_PT,
  title_font_size_pt = TITLE_FONT_SIZE_PT,
  image_label_grid_list = None,
  image_label_font_size_pt = IMAGE_LABEL_FONT_SIZE_PT,
  image_label_width_pt = IMAGE_LABEL_WIDTH_PT,
  image_label_height_pt = IMAGE_LABEL_HEIGHT_PT,
  node_size_min_freq = None,
  node_size_max_freq = None,
  node_size_min_px = None,
  node_size_max_px = None,
  legend_list = [],
  legend_height_spacing_pt = LEGEND_HEIGHT_SPACING_PT,
  legend_freq_ratio_color_bar_height_pt = LEGEND_FREQ_RATIO_COLOR_BAR_HEIGHT_PT,
  legend_title_font_size_pt = LEGEND_TITLE_FONT_SIZE_PT,
  legend_label_font_size_pt = LEGEND_LABEL_FONT_SIZE_PT,
  margin_label_top_list = None,
  margin_label_left_list = None,
  margin_top_font_size_pt = MARGIN_TOP_LABEL_FONT_SIZE_PT,
  margin_left_font_size_pt = MARGIN_LEFT_LABEL_FONT_SIZE_PT,
  margin_top_height_pt = MARGIN_TOP_HEIGHT_PT,
  margin_left_width_pt = MARGIN_LEFT_WIDTH_PT,
  margin_left_spill_over_pt = MARGIN_LEFT_SPILL_OVER_PT,
  margin_right_spill_over_pt = MARGIN_RIGHT_SPILL_OVER_PT,
):
  num_grids = len(image_grid_list)

  if total_width_frac_list is None:
    total_width_frac_list = [1] * num_grids

  image_label_show = image_label_grid_list is not None
  if image_label_show and (len(image_label_grid_list) != num_grids):
    raise Exception(
      f"Incorrect number of image label grids : {len(image_label_grid_list)}." +
      f" Expected {num_grids}."
    )

  margin_show_top = margin_label_top_list is not None
  margin_show_left = margin_label_left_list is not None
  if margin_show_left:
    if (margin_label_left_list is None) or (len(margin_label_left_list) != num_grids):
      raise Exception(
        f"Incorrect number of left margin lists: {len(margin_label_left_list)}." +
        f" Expected {num_grids}."
      )
  if margin_show_top:
    if (margin_label_left_list is None) or (len(margin_label_left_list) != num_grids):
      raise Exception(
        f"Incorrect number of top margin lists: {len(margin_label_left_list)}." +
        f" Expected {num_grids}."
      )

  title_slide_layout = prs.slide_layouts[6] # should be a blank layout
  slide = prs.slides.add_slide(title_slide_layout)
  
  slide_width_pt = prs.slide_width / pptx.util.Pt(1)

  y_pt = 0
  x_pt = 0
  y_legend_pt = 0
  x_legend_pt = slide_width_pt

  ## Legend constants ###
  legend_const = {
    'v': {
      'node_size': {
        'stride_pt': 20,
        'title_width_pt': 100,
        'title_height_pt': 10,
        'title_x_offset_pt': 0,
        'item_width_pt': 10,
        'label_width_pt': 50,
        'item_height_pt': 20,
        'label_height_pt': 10,
      },
      'node_outline': {
        'stride_pt': 12,
        'title_width_pt': 80,
        'title_height_pt': 10,
        'title_x_offset_pt': -10,
        'item_width_pt': 10,
        'label_width_pt': 70,
        'item_height_pt': 20,
        'label_height_pt': 10,
      },
      'variation_type': {
        'stride_pt': 12,
        'title_width_pt': 80,
        'title_height_pt': 10,
        'title_x_offset_pt': -10,
        'item_width_pt': 10,
        'label_width_pt': 50,
        'item_height_pt': 20,
        'label_height_pt': 20,
      },
      'node_type': {
        'stride_pt': 12,
        'title_width_pt': 80,
        'title_height_pt': 10,
        'title_x_offset_pt': -10,
        'item_width_pt': 10,
        'label_width_pt': 50,
        'item_height_pt': 20,
        'label_height_pt': 20,
      },
      'edge_type': {
        'stride_pt': 10,
        'title_width_pt': 80,
        'title_height_pt': 10,
        'title_x_offset_pt': -10,
        'item_width_pt': 20,
        'label_width_pt': 60,
        'item_height_pt': 20,
        'label_height_pt': 20,
      },
      'freq_ratio': {
        'stride_pt': 20,
        'title_width_pt': 90,
        'title_height_pt': 40,
        'title_x_offset_pt': -10,
        'item_width_pt': 30,
        'label_width_pt': 40,
        'item_height_pt': 20,
        'label_height_pt': 20,
      },
    },
    'h': {
      'node_size': {
        'stride_pt': 20,
        'title_width_pt': 100,
        'title_height_pt': 10,
        'title_x_offset_pt': 10,
        'item_width_pt': 30,
        'label_width_pt': 50,
        'item_height_pt': 20,
        'label_height_pt': 10,
      },
      'node_outline': {
        'stride_pt': 50,
        'title_width_pt': 90,
        'title_height_pt': 10,
        'title_x_offset_pt': 10,
        'item_width_pt': 30,
        'label_width_pt': 70,
        'item_height_pt': 10,
        'label_height_pt': 10,
      },
      'variation_type': {
        'stride_pt': 30,
        'title_width_pt': 90,
        'title_height_pt': 10,
        'title_x_offset_pt': 0,
        'item_width_pt': 30,
        'label_width_pt': 50,
        'item_height_pt': 10,
        'label_height_pt': 10,
      },
      'node_type': {
        'stride_pt': 35,
        'title_width_pt': 90,
        'title_height_pt': 10,
        'title_x_offset_pt': 5,
        'item_width_pt': 30,
        'label_width_pt': 50,
        'item_height_pt': 10,
        'label_height_pt': 10,
      },
      'edge_type': {
        'stride_pt': 40,
        'title_width_pt': 80,
        'title_height_pt': 10,
        'title_x_offset_pt': 0,
        'item_width_pt': 30,
        'label_width_pt': 60,
        'item_height_pt': 5,
        'label_height_pt': 10,
      },
      'freq_ratio': {
        'stride_pt': 40,
        'title_width_pt': 120,
        'title_height_pt': 20,
        'title_x_offset_pt': -20,
        'item_width_pt': 30,
        'label_width_pt': 40,
        'item_height_pt': 20,
        'label_height_pt': 20,
      },
    },
  }
  legend_x_offset_pt = {
    'v': 0,
    'h': 100,
  }

  # Title
  if title is not None:
    get_pptx_helpers.add_textbox_pptx(
      slide = slide,
      text = title,
      x_pt = 0,
      y_pt = y_legend_pt,
      width_pt = slide_width_pt,
      height_pt = title_height_pt,
      font_size_pt =  title_font_size_pt,
    )
    y_pt += title_height_pt

  for i in range(num_grids):
    max_image_width_px = 0
    max_image_height_px = 0

    # Get point/pixel ratio and dimensions
    # All images should be the same width/height
    for row in range(image_grid_list[i].shape[0]):
      for col in range(image_grid_list[i].shape[1]):
        image = PIL.Image.open(image_grid_list[i][row][col])
        max_image_width_px = max(max_image_width_px, image.width)
        max_image_height_px = max(max_image_height_px, image.height)

    content_width_px = image_grid_list[i].shape[1] * max_image_width_px

    content_width_pt = slide_width_pt * total_width_frac_list[i]
    if margin_show_left:
      content_width_pt -= margin_left_width_pt

    ratio_pt_px = content_width_pt / content_width_px
    cell_width_pt = ratio_pt_px * max_image_width_px
    cell_height_pt = ratio_pt_px * max_image_height_px

    content_x_pt = x_pt
    content_y_pt = y_pt
    content_height_pt = image_grid_list[i].shape[0] * cell_height_pt
    content_width_pt = image_grid_list[i].shape[1] * cell_width_pt

    if margin_show_top:
      content_y_pt += margin_top_height_pt
    if margin_show_left:
      content_x_pt += margin_left_width_pt

    # Content images
    get_pptx_helpers.add_picture_grid_pptx(
      slide = slide,
      file_name_grid = image_grid_list[i],
      x_pt = content_x_pt,
      y_pt = content_y_pt,
      cell_width_pt = cell_width_pt,
      cell_height_pt = cell_height_pt,
      cell_width_spacing_pt = 0,
      cell_height_spacing_pt = 0,
    )

    if image_label_show:
      if image_grid_list[i].shape != image_grid_list[i].shape:
        raise Exception(
          f"Incorrect shape of {i}'th label grid: {image_grid_list[i].shape}." +
          f" Expected: {image_grid_list[i].shape}."
        )
      # Image labels
      get_pptx_helpers.add_textbox_grid_pptx(
        slide = slide,
        text_grid = image_label_grid_list[i],
        x_pt = content_x_pt,
        y_pt = content_y_pt,
        cell_width_pt = cell_width_pt,
        cell_height_pt = cell_height_pt,
        cell_width_spacing_pt = 0,
        cell_height_spacing_pt = 0,
        text_width_pt = image_label_width_pt,
        text_height_pt = image_label_height_pt,
        font_size_pt = image_label_font_size_pt,
        text_align = 'center',
      )

    if margin_show_top:
      # Top margin
      get_pptx_helpers.add_textbox_grid_pptx(
        slide = slide,
        text_grid = np.array(margin_label_top_list[i]).reshape(1, -1),
        x_pt = margin_left_width_pt,
        y_pt = y_pt,
        cell_width_pt = cell_width_pt,
        cell_height_pt = margin_top_height_pt,
        cell_width_spacing_pt = 0,
        cell_height_spacing_pt = 0,
        text_height_pt = margin_top_height_pt,
        text_width_pt = cell_width_pt,
        font_size_pt = margin_top_font_size_pt,
      )

    if margin_show_left:
      # Left margin
      cell_width_pt_left_margin = (
        margin_left_width_pt +
        margin_left_spill_over_pt +
        margin_right_spill_over_pt
      )
      get_pptx_helpers.add_textbox_grid_pptx(
        slide = slide,
        text_grid = np.array(margin_label_left_list[i]).reshape(-1, 1),
        x_pt = -margin_left_spill_over_pt,
        y_pt = content_y_pt,
        cell_width_pt = cell_width_pt_left_margin,
        cell_height_pt = cell_height_pt,
        cell_width_spacing_pt = 0,
        cell_height_spacing_pt = 0,
        text_height_pt = cell_height_pt,
        text_width_pt = cell_width_pt_left_margin,
        font_size_pt = margin_left_font_size_pt,
      )

    y_pt = content_x_pt + content_height_pt
    
    # Size legends
    if any(legend['type'] == 'node_size' for legend in legend_list):
      for orientation in ['h', 'v']:
        y_legend_new_pt = get_pptx_legend.get_size_legend_pptx(
          slide = slide,
          x_pt = x_legend_pt + legend_x_offset_pt[orientation],
          y_pt = y_legend_pt,
          stride_pt = legend_const[orientation]['node_size']['stride_pt'],
          title_width_pt = legend_const[orientation]['node_size']['title_width_pt'],
          title_height_pt = legend_const[orientation]['node_size']['title_height_pt'],
          title_x_offset_pt = legend_const[orientation]['node_size']['title_x_offset_pt'],
          item_width_pt = legend_const[orientation]['node_size']['item_width_pt'],
          item_height_pt = legend_const[orientation]['node_size']['item_height_pt'],
          label_width_pt = legend_const[orientation]['node_size']['label_width_pt'],
          label_height_pt = legend_const[orientation]['node_size']['label_height_pt'],
          node_size_min_freq = node_size_min_freq,
          node_size_max_freq = node_size_max_freq,
          node_size_min_pt = node_size_min_px * ratio_pt_px,
          node_size_max_pt = node_size_max_px * ratio_pt_px,
          line_width_pt = LEGEND_SIZE_NODE_OUTLINE_WIDTH_PT,
          legend_title_font_size_pt = legend_title_font_size_pt,
          legend_label_font_size_pt = legend_label_font_size_pt,
          orientation = orientation,
        )
      y_legend_pt = y_legend_new_pt

  x_pt = slide_width_pt # place legends outside the actual slide

  # Legends
  for legend in legend_list:
    for orientation in ['h', 'v']:
      if legend['type'] == 'node_size':
        pass # already handled above
      elif legend['type'] == 'node_outline':
        y_legend_new_pt = get_pptx_legend.get_outline_legend_pptx(
          slide = slide,
          x_pt = x_pt + legend_x_offset_pt[orientation],
          y_pt = y_legend_pt,
          stride_pt = legend_const[orientation]['node_outline']['stride_pt'],
          title_width_pt = legend_const[orientation]['node_outline']['title_width_pt'],
          title_height_pt = legend_const[orientation]['node_outline']['title_height_pt'],
          title_x_offset_pt = legend_const[orientation]['node_outline']['title_x_offset_pt'],
          item_width_pt = legend_const[orientation]['node_outline']['item_width_pt'],
          item_height_pt = legend_const[orientation]['node_outline']['item_height_pt'],
          label_width_pt = legend_const[orientation]['node_outline']['label_width_pt'],
          label_height_pt = legend_const[orientation]['node_outline']['label_height_pt'],
          node_size_pt = LEGEND_NODE_SIZE_PT,
          line_width_pt = LEGEND_NODE_OUTLINE_WIDTH_PT,
          legend_title_font_size_pt = legend_title_font_size_pt,
          legend_label_font_size_pt = legend_label_font_size_pt,
          orientation = orientation,
        )
      elif legend['type'] == 'edge_type':
        y_legend_new_pt = get_pptx_legend.get_edge_legend_pptx(
          slide = slide,
          x_pt = x_pt + legend_x_offset_pt[orientation],
          y_pt = y_legend_pt,
          stride_pt = legend_const[orientation]['edge_type']['stride_pt'],
          title_width_pt = legend_const[orientation]['edge_type']['title_width_pt'],
          title_height_pt = legend_const[orientation]['edge_type']['title_height_pt'],
          title_x_offset_pt = legend_const[orientation]['edge_type']['title_x_offset_pt'],
          item_width_pt = legend_const[orientation]['edge_type']['item_width_pt'],
          item_height_pt = legend_const[orientation]['edge_type']['item_height_pt'],
          label_width_pt = legend_const[orientation]['edge_type']['label_width_pt'],
          label_height_pt = legend_const[orientation]['edge_type']['label_height_pt'],
          line_size_pt = LEGEND_EDGE_LINE_SIZE_PT,
          line_width_pt = LEGEND_EDGE_LINE_WIDTH_PT,
          legend_title_font_size_pt = legend_title_font_size_pt,
          legend_label_font_size_pt = legend_label_font_size_pt,
          orientation = orientation,
        )
      elif legend['type'] == 'variation_type':
        y_legend_new_pt = get_pptx_legend.get_variation_color_legend_pptx(
          slide = slide,
          x_pt = x_pt + legend_x_offset_pt[orientation],
          y_pt = y_legend_pt,
          stride_pt = legend_const[orientation]['variation_type']['stride_pt'],
          title_width_pt = legend_const[orientation]['variation_type']['title_width_pt'],
          title_height_pt = legend_const[orientation]['variation_type']['title_height_pt'],
          title_x_offset_pt = legend_const[orientation]['variation_type']['title_x_offset_pt'],
          item_width_pt = legend_const[orientation]['variation_type']['item_width_pt'],
          item_height_pt = legend_const[orientation]['variation_type']['item_height_pt'],
          label_width_pt = legend_const[orientation]['variation_type']['label_width_pt'],
          label_height_pt = legend_const[orientation]['variation_type']['label_height_pt'],
          variation_types = ['insertion', 'deletion', 'none'],
          node_size_pt = LEGEND_NODE_SIZE_PT,
          line_width_pt = LEGEND_NODE_OUTLINE_WIDTH_PT,
          legend_title_font_size_pt = legend_title_font_size_pt,
          legend_label_font_size_pt = legend_label_font_size_pt,
          orientation = orientation,
        )
      elif legend['type'] == 'node_type':
        y_legend_new_pt = get_pptx_legend.get_node_legend_pptx(
          slide = slide,
          x_pt = x_pt + legend_x_offset_pt[orientation],
          y_pt = y_legend_pt,
          stride_pt = legend_const[orientation]['node_type']['stride_pt'],
          title_width_pt = legend_const[orientation]['node_type']['title_width_pt'],
          title_height_pt = legend_const[orientation]['node_type']['title_height_pt'],
          title_x_offset_pt = legend_const[orientation]['node_type']['title_x_offset_pt'],
          item_width_pt = legend_const[orientation]['node_type']['item_width_pt'],
          item_height_pt = legend_const[orientation]['node_type']['item_height_pt'],
          label_width_pt = legend_const[orientation]['node_type']['label_width_pt'],
          label_height_pt = legend_const[orientation]['node_type']['label_height_pt'],
          node_size_pt = LEGEND_NODE_SIZE_PT,
          line_width_pt = LEGEND_NODE_OUTLINE_WIDTH_PT,
          legend_title_font_size_pt = legend_title_font_size_pt,
          legend_label_font_size_pt = legend_label_font_size_pt,
          orientation = orientation,
        )
      elif legend['type'] == 'freq_ratio':
        construct_1 = legend['construct_1']
        construct_2 = legend['construct_2']
        if orientation == 'v':
          title = f'Ratio\n[{library_constants.LABELS[construct_1]} / {library_constants.LABELS[construct_2]}]'
        elif orientation == 'h':
          title = f'Ratio [{library_constants.LABELS[construct_1]} / {library_constants.LABELS[construct_2]}]'
        else:
          raise Exception('Impossible')
        color_bar_file = legend['color_bar_file']
        y_legend_new_pt = get_pptx_legend.get_freq_ratio_legend_pptx(
          slide = slide,
          x_pt = x_pt + legend_x_offset_pt[orientation],
          y_pt = y_legend_pt,
          title = title,
          title_width_pt = legend_const[orientation]['freq_ratio']['title_width_pt'],
          title_height_pt = legend_const[orientation]['freq_ratio']['title_height_pt'],
          title_x_offset_pt = legend_const[orientation]['freq_ratio']['title_x_offset_pt'],
          label_width_pt = legend_const[orientation]['freq_ratio']['label_width_pt'],
          label_height_pt = legend_const[orientation]['freq_ratio']['label_height_pt'],
          color_bar_minor_axis_pt = LEGEND_FREQ_RATIO_COLOR_BAR_WIDTH_PT,
          color_bar_major_axis_pt = legend_freq_ratio_color_bar_height_pt,
          color_bar_file = color_bar_file,
          legend_title_font_size_pt = legend_title_font_size_pt,
          legend_label_font_size_pt = legend_label_font_size_pt,
          orientation = orientation,
        )
      else:
        raise Exception('Unknown legend type: ' + str(legend))
    y_legend_pt = y_legend_new_pt
    y_legend_pt += legend_height_spacing_pt


def parse_args():
  parser = argparse.ArgumentParser(
    description = 'Create powerpoint figures from images.'
  )
  parser.add_argument(
    '--input',
    nargs = '+',
    type = common_utils.check_file,
    help = 'List of images to include in the grid',
    required = True,
  )
  parser.add_argument(
    '--labels',
    nargs = '+',
    help = 'Labels of images in the grid',
  )
  parser.add_argument(
    '--top_margin_labels',
    nargs = '+',
    help = 'Labels in the top margins of the grids',
  )
  parser.add_argument(
    '--left_margin_labels',
    nargs = '+',
    help = 'Labels in the left margins of the grids',
  )
  parser.add_argument(
    '--output',
    help = 'Output PPTX file.',
    required = True,
  )
  parser.add_argument(
    '--num_grids',
    required = True,
    type = int,
    help = 'Number of separate grids to create',
  )
  parser.add_argument(
    '--num_rows',
    nargs = '+',
    required = True,
    type = int,
    help = (
      'Number of rows in each grid.' +
      ' Number of arguments should match the number of grids.'
    ),
  )
  parser.add_argument(
    '--num_cols',
    nargs = '+',
    required = True,
    type = int,
    help = (
      'Number of columns in each grid.' +
      ' Number of arguments should match the number of grids.'
    ),
  )
  parser.add_argument(
    '--total_width',
    nargs = '+',
    type = float,
    help = (
      'Fraction of the page width to use for each grid .' +
      ' Number of arguments should match the number of grids.'
    ),
  )
  parser.add_argument(
    '--node_max_freq',
    type = float,
    help = (
      'Max frequency to determine node size.' +
      ' Higher frequencies are clipped to this value.'
    ),
    default = library_constants.GRAPH_NODE_SIZE_MAX_FREQ,
  )
  parser.add_argument(
    '--node_min_freq',
    type = float,
    help = (
      'Min frequency to determine node size.' +
      ' Lower frequencies are clipped to this value.'
    ),
    default = library_constants.GRAPH_NODE_SIZE_MIN_FREQ,
  )
  parser.add_argument(
    '--node_max_px',
    type = float,
    help = 'Largest node size as determined by the frequency.',
    default = library_constants.GRAPH_NODE_SIZE_MAX_PX,
  )
  parser.add_argument(
    '--node_min_px',
    type = float,
    help = 'Smallest node size as determined by the frequency.',
    default = library_constants.GRAPH_NODE_SIZE_MIN_PX,
  )
  parser.add_argument(
    '--title',
    default = None,
    help = 'Page title',
  )
  parser.add_argument(
    '--legends',
    nargs = '+',
    choices = list(LEGENDS),
    default = None,
    help = 'Legends to draw outside the page.',
  )
  parser.add_argument(
    '--template',
    default = os.path.join(os.path.dirname(__file__), 'template.pptx'),
    help = 'The PPTX file to use as a template. Controls the page size.',
  )
  args = parser.parse_args()
  return args

if __name__ == '__main__':
  args = parse_args()

  if args.num_grids != len(args.num_rows):
    raise Exception(
      f'Incorrect num rows specification: {args.num_rows}.' +
      f' Expected {args.num_grids} values.' 
    )
  
  if args.num_grids != len(args.num_cols):
    raise Exception(
      f'Incorrect num columns specification: {args.num_cols}.' +
      f' Expected {args.num_grids} values.' 
    )

  if args.total_width is None:
    args.total_width = [1] * args.num_grids
  if len(args.total_width) != args.num_grids:
    raise Exception(
      f'Incorrect number of total widths: {args.total_width}.' +
      f' Expected {args.num_grids} values.' 
    )

  num_images_total = sum(r * c for r, c in zip(args.num_rows, args.num_cols))
  if num_images_total != len(args.input):
    raise Exception(
      f'Incorrect number of input files: {len(args.input)}.' +
      f' Expected {num_images_total} values.' 
    )

  image_grid_list = []
  image_index = 0
  for i in range(args.num_grids):
    num_images = args.num_rows[i] * args.num_cols[i]
    image_grid = np.array(
      [args.input[image_index + j] for j in range(num_images)]
    )
    image_grid = image_grid.reshape((args.num_rows[i], args.num_cols[i]))
    image_grid_list.append(image_grid)
    image_index += num_images

  label_grid_list = None
  if args.labels is not None:
    if len(args.labels) != num_images_total:
      raise Exception(
        f'Incorrect number of input labels: {len(args.input)}.' +
        f' Expected {num_images_total} values.' 
      )
    label_grid_list = []
    index = 0
    for i in range(args.num_grids):
      num_labels = args.num_rows[i] * args.num_cols[i]
      label_grid = np.array(
        [args.labels[index + j] for j in range(num_labels)]
      )
      label_grid = label_grid.reshape((args.num_rows[i], args.num_cols[i]))
      label_grid_list.append(label_grid)
      index += num_labels
  
  margin_label_left_list = None
  if args.left_margin_labels is not None:
    num_left_margin_labels = sum(args.num_rows[i] for i in range(args.num_grids))
    if len(args.left_margin_labels) != num_left_margin_labels:
      raise Exception(
        f'Incorrect number of left margin labels: {len(args.left_margin_labels)}.' +
        f' Expected {num_left_margin_labels} values.' 
      )
    margin_label_left_list = []
    index = 0
    for i in range(args.num_grids):
      num_labels = args.num_rows[i]
      label_list = [args.left_margin_labels[index + j] for j in range(num_labels)]
      margin_label_left_list.append(label_list)
      index += num_labels
  
  margin_label_top_list = None
  if args.top_margin_labels is not None:
    num_top_margin_labels = sum(args.num_cols[i] for i in range(args.num_grids))
    if len(args.top_margin_labels) != num_top_margin_labels:
      raise Exception(
        f'Incorrect number of top margin labels: {len(args.top_margin_labels)}.' +
        f' Expected {num_top_margin_labels} values.' 
      )
    margin_label_top_list = []
    index = 0
    for i in range(args.num_grids):
      num_labels = args.num_cols[i]
      label_list = [args.top_margin_labels[index + j] for j in range(num_labels)]
      margin_label_top_list.append(label_list)
      index += num_labels

  legend_list = []
  if args.legends is not None:
    legend_list = [LEGENDS[x] for x in args.legends]

  # replace \n with newline characters in case using Windows
  for i in range(args.num_grids):
    if label_grid_list is not None:
      for r in range(args.num_rows[i]):
        for c in range(args.num_cols[i]):
          label_grid_list[i][r, c] = label_grid_list[i][r, c].replace('\\n', '\n')
    if margin_label_left_list is not None:
      for r in range(args.num_rows[i]):
        margin_label_left_list[i][r] = margin_label_left_list[i][r].replace('\\n', '\n')
    if margin_label_top_list is not None:
      for c in range(args.num_cols[i]):
        margin_label_top_list[i][c] = margin_label_top_list[i][c].replace('\\n', '\n')

  prs = pptx.Presentation(PPTX_TEMPLATE_FILE)

  get_slide(
    prs,
    title = args.title,
    image_grid_list = image_grid_list,
    image_label_grid_list = label_grid_list,
    margin_label_left_list = margin_label_left_list,
    margin_label_top_list = margin_label_top_list,
    node_size_max_freq = args.node_max_freq,
    node_size_min_freq = args.node_min_freq,
    node_size_max_px = args.node_max_px,
    node_size_min_px = args.node_min_px,
    legend_list = legend_list,
    total_width_frac_list = args.total_width,
  )
  
  log_utils.log(args.output)
  file_utils.make_parent_dir(args.output)
  prs.save(args.output)
