import numpy as np

VERSION_NONE = 'versionNone'
VERSION_OLD = 'old'
VERSION_NEW = 'new'
VERSION_MERGED = 'merged'
VERSIONS = [
  VERSION_NONE,
  VERSION_OLD,
  VERSION_NEW,
  VERSION_MERGED,
]

SUBST_WITH = 'withSubst'
SUBST_WITHOUT = 'withoutSubst'
SUBST_TYPES = [SUBST_WITH, SUBST_WITHOUT]

def check_subst_type(subst_type):
  if subst_type not in SUBST_TYPES:
    raise Exception('Not a valid subst type: ' + str(subst_type))

FREQ = 'freq'
FREQ_FILTER = 'freq_filter'
FREQ_FILTER_MEAN = 'freq_filter_mean'
COUNT = 'count'
FREQ_TYPES = [COUNT, FREQ, FREQ_FILTER, FREQ_FILTER_MEAN]

def check_freq_type(freq_type):
  if freq_type not in FREQ_TYPES:
    raise Exception('Not a valid freq type: ' + str(freq_type))

def get_position_labels(label_type, ref_length):
  if label_type == 'relative':
    return (
      [str(-x) for x in range(1, 1 + ref_length // 2)][::-1] +
      [str(x) for x in range(1, 1 + ref_length // 2)]
    )
  elif label_type == 'absolute':
    return [str(x) for x in range(1, ref_length + 1)]
  else:
    raise Exception('Unknown label type: ' + str(label_type))

POSITION_TITLE = {
  'relative': 'Position (from cut)',
  'absolute': 'Position',
}

EDGE_TYPES = {
  'substitution': {
    'label': 'Substitution',
    'line_dash': 'dash',
    'legend_color': 'black',
    'plot_color': 'rgba(0,0,0,0.5)',
  },
  'indel': {
    'label': 'In/Del',
    'line_dash': 'solid',
    'legend_color': 'black',
    'plot_color': 'rgba(0,0,0,0.5)',
  },
}

VARIATION_TYPES = {
  'insertion': {
    'label': 'Insertion',
    'short_label': 'I',
    'color': '#ffa500',
    'color_3d': '#ffa500',
  },
  'deletion': {
    'label': 'Deletion',
    'short_label': 'D',
    'color': '#8080ff',
    'color_3d': '#8080ff',
  },
  'substitution': {
    'label': 'Substitution',
    'short_label': 'S',
    'color': '#808080',
    'color_3d': '#bfbfbf',
  },
  'mixed': {
    'label': 'Mixed',
    'short_label': 'M',
    'color': '#00ff00',
  },
  'none': {
    'label': 'None',
    'short_label': 'N',
    'color': '#ffffff',
  },
}

LABEL_REFERENCE = 'Reference'
LABEL_NONREFERENCE = 'Nonreference'

SIMILAR_FREQ_COLOR = '#ffffff'

FREQ_GROUP_A = 'A'
FREQ_GROUP_B = 'B'
FREQ_GROUP_C = 'C'
FREQ_RATIO_A = 3/2
FREQ_RATIO_C = 2/3
FREQ_RATIO_LOG_A = np.log(FREQ_RATIO_A)
FREQ_RATIO_LOG_C = np.log(FREQ_RATIO_C)

FREQ_RATIO_COLOR_SCALE_LOG_RANGE = [np.log(2/3), np.log(3/2)]
FREQ_RATIO_COLOR_BAR_TICK_VALS = [
  np.log(2/3),
  np.log(4/5),
  np.log(1),
  np.log(5/4),
  np.log(3/2),
]
FREQ_RATIO_COLOR_BAR_TICK_TEXT = [
  '2/3',
  '4/5',
  '1',
  '5/4',
  '3/2',
]

def get_freq_ratio_color_scale(color_1, color_2):
  return [
    [0, color_2],
    [0.5, SIMILAR_FREQ_COLOR],
    [1, color_1],
  ]

def get_freq_ratio_label(freq_group, label_1, label_2):
  if freq_group == FREQ_GROUP_A:
    return (
      f'Higher in {label_1}<br>' +
      f'(ratio > {FREQ_RATIO_A:0.2f})'
    )
  elif freq_group == FREQ_GROUP_B:
    return (
      f'Similar in both<br>' + 
      f'({FREQ_RATIO_C:0.2f} ≤ ratio ≤ {FREQ_RATIO_A:0.2f})'
    )
  elif freq_group == FREQ_GROUP_C:
    return (
      f'Higher in {label_2}<br>' + 
      f'(ratio < {FREQ_RATIO_C:0.2f})'
    )
  else:
    raise Exception('Unknown freq group: ' + str(freq_group))

### Data formats ###
DATA_INDIVIDUAL = 'individual'
DATA_COMPARISON = 'comparison'

FREQ_COLUMNS = {
  DATA_INDIVIDUAL: ['freq_mean'],
  DATA_COMPARISON: ['freq_mean_1', 'freq_mean_2'],
}

def is_freq_column(x):
  return (
    (x in FREQ_COLUMNS[DATA_INDIVIDUAL]) or
    (x in FREQ_COLUMNS[DATA_COMPARISON])
  )

FREQ_RANK_COLUMNS = {
  DATA_INDIVIDUAL: ['freq_mean_rank'],
  DATA_COMPARISON: ['freq_mean_rank_1', 'freq_mean_rank_2'],
}

def get_data_label(data_info):
  if data_info['format'] == DATA_INDIVIDUAL:
    return data_info['label']
  elif data_info['format'] == DATA_COMPARISON:
    return '_'.join([data_info['label_1'], data_info['label_2']])
  else:
    raise Exception('Unknown format: ' + str(data_info['format']))

### Constants for 3D variation-position histograms ###
HISTOGRAM_TITLE_FONT_SIZE = 16
HISTOGRAM_AXIS_LABEL_FONT_SIZE = 12
HISTOGRAM_AXIS_TICK_FONT_SIZE = 8
HISTOGRAM_AXIS_TICK_MODULUS = 4
HISTOGRAM_FONT_SIZE_SCALE = 6
HISTOGRAM_WIDTH_PX = 1500
HISTOGRAM_HEIGHT_PX = 1500
HISTOGRAM_MARGIN_LEFT_PX = 50
HISTOGRAM_MARGIN_RIGHT_PX = 300
HISTOGRAM_MARGIN_TOP_PX = 0
HISTOGRAM_MARGIN_BOTTOM_PX = 100
HISTOGRAM_DPI = 100
HISTOGRAM_FREQ_RANGE = [1e-5, 1]
HISTOGRAM_VARIATION_TYPE_COLORS = [
  VARIATION_TYPES['insertion']['color_3d'],
  VARIATION_TYPES['deletion']['color_3d'],
  VARIATION_TYPES['substitution']['color_3d'],
]
BASE_FIG_SIZE = 12

### Constants for graphs ###
GRAPH_LAYOUT_TYPE = 'universal_layout'
GRAPH_LAYOUT_SEPARATE_COMPONENTS = False
GRAPH_WIDTH_PX = 2400
GRAPH_HEIGHT_PX = 2400
GRAPH_NODE_SUBST_TYPE = SUBST_WITHOUT
GRAPH_NODE_SIZE_FREQ_RANGE = [1e-5, 1]
GRAPH_NODE_SIZE_PX_RANGE = [10, 120]
GRAPH_NODE_FILTER_VARIATION_TYPES = [
  'insertion',
  'deletion',
  'substitution',
  'mixed',
  'none',
]
GRAPH_NODE_FILTER_FREQ_RANGE = [0, np.inf]
GRAPH_NODE_FILTER_DIST_RANGE = [0, np.inf]
GRAPH_NODE_OUTLINE_WIDTH_SCALE = 4
GRAPH_NODE_LABEL_SHOW = False
GRAPH_NODE_LABEL_COLUMNS = ['id']
GRAPH_NODE_LABEL_POSITION = 'bottom center'
GRAPH_NODE_SIZE_TYPE = 'freq'
GRAPH_NODE_COLOR_TYPE_COMPARISON = 'freq_ratio'
GRAPH_NODE_COLOR_TYPE_INDIVIDUAL = 'variation_type'
GRAPH_NODE_COMPARISON_COLORS = ['#ff0000', '#0000ff']
GRAPH_NODE_REFERENCE_OUTLINE_COLOR = '#32cd32'
GRAPH_NODE_REFERENCE_OUTLINE_WIDTH = 2
GRAPH_NODE_OUTLINE_COLOR = '#000000'
GRAPH_NODE_OUTLINE_WIDTH = 1
GRAPH_NODE_FILL_COLOR = '#ffffff'
GRAPH_NODE_VARIATION_TYPE_COLORS = [
  VARIATION_TYPES['insertion']['color'],
  VARIATION_TYPES['deletion']['color'],
  VARIATION_TYPES['substitution']['color'],
  VARIATION_TYPES['mixed']['color'],
  VARIATION_TYPES['none']['color'],
]
GRAPH_EDGE_SHOW = True
GRAPH_EDGE_SHOW_LABELS = False
GRAPH_EDGE_SHOW_TYPES = ['indel']
GRAPH_EDGE_LABELS_SHOW = True
GRAPH_EDGE_WIDTH_SCALE = 8
GRAPH_LINE_WIDTH_SCALE = 8
GRAPH_FONT_SIZE_SCALE = 2
GRAPH_STATS_SHOW = False
GRAPH_STATS_SEPARATE = False
GRAPH_STATS_WIDTH_PX = 800
GRAPH_STATS_X_FRAC = 0
GRAPH_STATS_Y_FRAC = 1
GRAPH_STATS_X_SHIFT = 20
GRAPH_STATS_Y_SHIFT = -20
GRAPH_STATS_X_ANCHOR = 'left'
GRAPH_STATS_Y_ANCHOR = 'top'
GRAPH_TITLE = None
GRAPH_TITLE_Y_SHIFT_PX = 100
GRAPH_TITLE_FONT_SIZE = 30
GRAPH_AXES_TITLE_FONT_SIZE = 20
GRAPH_AXES_TICK_FONT_SIZE = 16
GRAPH_AXES_SHOW = False
GRAPH_LEGEND_CUSTOM_SHOW = True
GRAPH_LEGEND_CUSTOM_X_ANCHOR_FRAC = 1
GRAPH_LEGEND_CUSTOM_Y_ANCHOR_FRAC = 1
GRAPH_LEGEND_PLOTLY_SHOW = False
GRAPH_LEGEND_X_SHIFT_PX = 0
GRAPH_LEGEND_VERTICAL_SPACE_PX = 200
GRAPH_LEGEND_ITEM_SCALE = 1
GRAPH_LEGEND_TITLE_FONT_SIZE = 24
GRAPH_LEGEND_GROUP_TITLE_FONT_SIZE = 20
GRAPH_LEGEND_FONT_SIZE = 18
GRAPH_LEGEND_COLORBAR_SCALE = 4
GRAPH_LEGEND_COLORBAR_HEIGHT_PX = 500
GRAPH_LEGEND_EDGE_ITEM_LINE_SIZE_PX = 100
GRAPH_LEGEND_EDGE_ITEM_LINE_WIDTH_PX = 2.5
GRAPH_LEGEND_COLORBAR_WIDTH_PX = 50
GRAPH_BACKGROUND_COLOR = 'white'
GRAPH_LABEL_FONT_SIZE = 16
GRAPH_MARGIN_TOP_MIN_PX = 300
GRAPH_MARGIN_BOTTOM_MIN_PX = 300
GRAPH_MARGIN_LEFT_MIN_PX = 300
GRAPH_MARGIN_RIGHT_MIN_PX = 300
GRAPH_HOVER_LABEL_FONT = 'Courier New, monospace'
GRAPH_HOVER_LABEL_FONT_SIZE = 16
GRAPH_HOVER_LABEL_BG_COLOR = 'white'
GRAPH_LEGEND_PLOTLY_ITEM_WIDTH_PX = 100
GRAPH_PLOT_RANGE_X = [float('nan'), float('nan')]
GRAPH_PLOT_RANGE_Y = [float('nan'), float('nan')]
GRAPH_CROP_X = [0, 1]
GRAPH_CROP_Y = [0, 1]
GRAPH_LAYOUT_PRECOMPUTED_DIR = None
GRAPH_SEQUENCE_REVERSE_COMPLEMENT = False

### Universal layout constants ###
GRAPH_UNIVERSAL_LAYOUT_INSERTION_ROW_SPEC = {
  1: {'rows': 1, 'cols': 4, 'row_space': 2},
  2: {'rows': 1, 'cols': 16, 'row_space': 2},
  3: {'rows': 2, 'cols': 32, 'row_space': 1},
  4: {'rows': 2, 'cols': 128, 'row_space': 1},
  5: {'rows': 4, 'cols': 256, 'row_space': 0.5},
  6: {'rows': 8, 'cols': 512, 'row_space': 0.25},
  7: {'rows': 8, 'cols': 2048, 'row_space': 0.25},
  8: {'rows': 8, 'cols': 8192, 'row_space': 0.25},
}
GRAPH_UNIVERSAL_LAYOUT_X_SCALE_INSERTION = 10
GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_INSERTION = 3
GRAPH_UNIVERSAL_LAYOUT_X_SCALE_DELETION = 2
GRAPH_UNIVERSAL_LAYOUT_Y_SCALE_DELETION = 1