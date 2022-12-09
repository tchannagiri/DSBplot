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

# Substitution edges not being shown
# indel edges being show solid
EDGE_TYPES = {
  'substitution': {
    'label': 'Substitution',
    'line_dash': 'solid',
    'legend_color': 'black',
    'plot_color': 'rgba(0,0,0,0.5)',
  },
  'indel': {
    'label': 'In/Del',
    # 'line_dash': 'dash',
    'line_dash': 'solid',
    'legend_color': 'black',
    'plot_color': 'rgba(0,0,0,0.5)',
  },
}

VARIATION_NONE = 'none'
VARIATION_MIXED = 'mixed'
VARIATION_SUBSTITUTION = 'substitution'
VARIATION_DELETION = 'deletion'
VARIATION_INSERTION = 'insertion'
VARIATION_TYPES = {
  VARIATION_NONE: {
    'label': 'None',
    'short_label': 'N',
    'color': '#FFFFFF',
  },
  VARIATION_MIXED: {
      'label': 'Mixed',
      'short_label': 'M',
      'color': '#00FF00',
  },
  VARIATION_SUBSTITUTION: {
    'label': 'Substitution',
    'short_label': 'S',
    'color': '#808080',
    'color_3d': '#BFBFBF',
  },
  VARIATION_DELETION: {
    'label': 'Deletion',
    'short_label': 'D',
    'color': '#8080FF',
    'color_3d': '#8080FF',
  },
  VARIATION_INSERTION: {
    'label': 'Insertion',
    'short_label': 'I',
    # 'color': '#FF8080',
    # 'color_3d': '#FF8080',
    'color': '#FFA500',
    'color_3d': '#FFA500',
  },
}

REFERENCE_OUTLINE_COLOR = '#32cd32'
REFERENCE_OUTLINE_WIDTH = 2

DEFAULT_OUTLINE_COLOR = '#000000'
DEFAULT_OUTLINE_WIDTH = 1
DEFAULT_NODE_COLOR = '#FFFFFF'

REFERENCE_DESCRIPTION = 'Reference'
NON_REFERENCE_DESCRIPTION = 'Non-reference'

### Controls ###
CONTROL_NOT = 'notControl'
CONTROL_NODSB = 'noDSB'
CONTROL_30BPDOWN = '30bpDown'
CONTROLS = [CONTROL_NOT, CONTROL_NODSB, CONTROL_30BPDOWN]

### DSB types ###
DSB_TYPE_1 = '1DSB'
DSB_TYPE_2 = '2DSB'
DSB_TYPE_2anti = '2DSBanti'
DSB_TYPES = [DSB_TYPE_1, DSB_TYPE_2, DSB_TYPE_2anti]

### Hguide types ###
GUIDE_RNA_A = 'sgA'
GUIDE_RNA_B = 'sgB'
GUIDE_RNA_AB = 'sgAB'
GUIDE_RNA_CD = 'sgCD'
GUIDE_RNAS = [GUIDE_RNA_A, GUIDE_RNA_B, GUIDE_RNA_AB, GUIDE_RNA_CD]

### Strand type ###
STRAND_R1 = 'R1'
STRAND_R2 = 'R2'
STRANDS = [STRAND_R1, STRAND_R2]

### Cell line ###
CELL_LINE_WT = 'WT'
CELL_LINE_KO = 'KO'
CELL_LINES = [CELL_LINE_WT, CELL_LINE_KO]

CONSTRUCT_COLOR = {
  'sense': '#CF191B', 
  'branch': '#33A02C',
  'cmv': '#FFE669',
  'antisense': '#CF191B',
  'splicing': '#33A02C',
  'sense_branch': '#ffffff',
  'sense_cmv': '#ffffff',
  'antisense_splicing': '#ffffff',
}

SIMILAR_FREQ_COLOR = '#ffffff'

FREQ_GROUP_A = 'A'
FREQ_GROUP_B = 'B'
FREQ_GROUP_C = 'C'
FREQ_RATIO_A = 3/2
FREQ_RATIO_C = 2/3
FREQ_RATIO_LOG_A = np.log(FREQ_RATIO_A)
FREQ_RATIO_LOG_C = np.log(FREQ_RATIO_C)

FREQ_RATIO_COLOR_SCALE_LOG_MAX = np.log(3/2)
FREQ_RATIO_COLOR_SCALE_LOG_MIN = np.log(2/3)
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

def get_freq_ratio_color_scale(construct_1, construct_2):
  return [
    [0, CONSTRUCT_COLOR[construct_2]],
    [0.5, CONSTRUCT_COLOR[construct_1 + '_' + construct_2]],
    [1, CONSTRUCT_COLOR[construct_1]],
  ]

def get_freq_ratio_label(freq_group, construct_1, construct_2):
  if freq_group == FREQ_GROUP_A:
    return (
      f'Higher in {LABELS[construct_1]}<br>' +
      f'(ratio > {FREQ_RATIO_A:0.2f})'
    )
  elif freq_group == FREQ_GROUP_B:
    return (
      f'Similar in both<br>' + 
      f'({FREQ_RATIO_C:0.2f} ≤ ratio ≤ {FREQ_RATIO_A:0.2f})'
    )
  elif freq_group == FREQ_GROUP_C:
    return (
      f'Higher in {LABELS[construct_2]}<br>' + 
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

FREQ_RANK_COLUMNS = {
  DATA_INDIVIDUAL: ['freq_mean_rank'],
  DATA_COMPARISON: ['freq_mean_rank_1', 'freq_mean_rank_2'],
}

### Constructs ###

## Individual ##
CONSTRUCT_SENSE = 'sense'
CONSTRUCT_BRANCH = 'branch'
CONSTRUCT_CMV = 'cmv'
CONSTRUCT_ANTISENSE = 'antisense'
CONSTRUCT_SPLICING = 'splicing'
CONSTRUCTS_INDIVIDUAL = [
  CONSTRUCT_SENSE,
  CONSTRUCT_BRANCH,
  CONSTRUCT_CMV,
  CONSTRUCT_ANTISENSE,
  CONSTRUCT_SPLICING,
]
CONSTRUCTS_INDIVIDUAL_SENSE = [
  CONSTRUCT_SENSE,
  CONSTRUCT_BRANCH,
  CONSTRUCT_CMV,
]
CONSTRUCTS_INDIVIDUAL_ANTISENSE = [
  CONSTRUCT_ANTISENSE,
  CONSTRUCT_SPLICING,
]

## Comparisons ##
CONSTRUCT_SENSE_BRANCH = CONSTRUCT_SENSE + '_' + CONSTRUCT_BRANCH
CONSTRUCT_SENSE_CMV = CONSTRUCT_SENSE + '_' + CONSTRUCT_CMV
CONSTRUCT_ANTISENSE_SPLICING = CONSTRUCT_ANTISENSE + '_' + CONSTRUCT_SPLICING
CONSTRUCTS_COMPARISON_SENSE = [
  CONSTRUCT_SENSE_BRANCH,
  CONSTRUCT_SENSE_CMV,
]
CONSTRUCTS_COMPARISON_ANTISENSE = [
  CONSTRUCT_ANTISENSE_SPLICING
]
CONSTRUCTS_COMPARISON = [
  CONSTRUCT_SENSE_BRANCH,
  CONSTRUCT_SENSE_CMV,
  CONSTRUCT_ANTISENSE_SPLICING,
]

### Labels ###
LABELS = {
  '1DSB': '1 DSB',
  '2DSB': '2 DSB',
  '2DSBanti': '2 DSB (antisense)',
  'sgA': 'sgRNA A',
  'sgB': 'sgRNA B',
  'sgC': 'sgRNA C/C\'',
  'sgD': 'sgRNA D',
  'sgAB': 'sgRNA A & B',
  'sgCD': 'sgRNA C/C\' & D',
  'KO': 'KO',
  'WT': 'WT',
  'R1': 'Forward strand',
  'R2': 'Reverse strand',
  'NHEJ': 'NHEJ',
  'MMEJ': 'MMEJ',
  'sense': 'Sense',
  'branch': 'Branch∆',
  'cmv': 'pCMV∆',
  'splicing': '5\'-Splicing∆',
  'antisense': 'Antisense',
  'sense_branch': 'Sense & Branch∆',
  'sense_cmv': 'Sense & pCMV∆',
  'antisense_splicing': 'Antisense & 5\'-Splicing∆',
  'noDSB': 'No DSB',
  '30bpDown': '30bp Down',
  'notControl': '',
  'ref_pos': 'Reference sequence position',
  'ref_cut_pos_offset': 'Reference sequence position (from cut)',
  'substitution': 'Substitution',
  'insertion': 'Insertion',
  'deletion': 'Deletion',
}

def get_data_label(data_info):
  if data_info['format'] == DATA_INDIVIDUAL:
    construct_str = data_info['construct']
  elif data_info['format'] == DATA_COMPARISON:
    construct_str = '_'.join([data_info['construct_1'], data_info['construct_2']])
  else:
    raise Exception('Unknown format: ' + str(data_info['format']))
  control_str = None if (data_info['control_type'] == CONTROL_NOT) else data_info['control_type']
  version_str = None if (data_info['version'] == VERSION_NONE) else data_info['version']
  str_list = [
    data_info['cell_line'],
    data_info['guide_rna'],
    data_info['strand'],
    construct_str,
    control_str,
    version_str,
  ]
  return '_'.join(x for x in str_list if x is not None)
  return data_info['name']

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
HISTOGRAM_FREQ_MIN = 1e-5
HISTOGRAM_FREQ_MAX = 1
BASE_FIG_SIZE = 12

### Constants for graphs ###
GRAPH_WIDTH_PX = 2400
GRAPH_HEIGHT_PX = 2400
GRAPH_NODE_SIZE_MIN_FREQ = 1e-5
GRAPH_NODE_SIZE_MAX_FREQ = 1
GRAPH_NODE_FILTER_VARIATION_TYPES = ['insertion', 'deletion', 'none']
GRAPH_NODE_SIZE_MIN_PX = 10
GRAPH_NODE_SIZE_MAX_PX = 120
GRAPH_NODE_OUTLINE_WIDTH_SCALE = 4
GRAPH_EDGE_WIDTH_SCALE = 8
GRAPH_LINE_WIDTH_SCALE = 8
GRAPH_FONT_SIZE_SCALE = 2
GRAPH_SUBPLOT_WIDTH_PX = 1600
GRAPH_SUBPLOT_HEIGHT_PX = 1000
GRAPH_SUBPLOT_ROW_SPACE_PX = 100
GRAPH_SUBPLOT_COL_SPACE_PX = 100
GRAPH_SUBPLOT_TITLE_FONT_SIZE = 24
GRAPH_STATS_SUBPLOT_PX = 800
GRAPH_DESCRIPTION_HEIGHT_PX = 700
GRAPH_TITLE_HEIGHT_PX = 200
GRAPH_TITLE_FONT_SIZE = 30
GRAPH_AXES_TITLE_FONT_SIZE = 20
GRAPH_AXES_TICK_FONT_SIZE = 16
GRAPH_LEGEND_WIDTH_PX = 400
GRAPH_LEGEND_VERTICAL_SPACE_PX = 100
GRAPH_LEGEND_TITLE_FONT_SIZE = 24
GRAPH_LEGEND_GROUP_TITLE_FONT_SIZE = 20
GRAPH_LEGEND_FONT_SIZE = 18
GRAPH_LEGEND_COLORBAR_SCALE = 8
GRAPH_LEGEND_COLORBAR_HEIGHT_PX = 500
GRAPH_EDGE_LEGEND_ITEM_LINE_SIZE_PX = 100
GRAPH_EDGE_LEGEND_ITEM_LINE_WIDTH_PX = 2.5
GRAPH_LEGEND_COLORBAR_WIDTH_PX = 50
GRAPH_BACKGROUND_COLOR = 'white'
GRAPH_LABEL_FONT_SIZE = 16
GRAPH_MARGIN_ROW_HEIGHTS_PX = [500]
GRAPH_MARGIN_COL_WIDTHS_PX = [800, 1500, 1000]
GRAPH_MARGIN_TOP_MIN_PX = 300
GRAPH_MARGIN_BOTTOM_MIN_PX = 300
GRAPH_MARGIN_LEFT_MIN_PX = 300
GRAPH_MARGIN_RIGHT_MIN_PX = 300
GRAPH_MARGIN_FONT_SIZES_TOP = [30, 30]
GRAPH_MARGIN_FONT_SIZES_LEFT = [30, 30, 20]
GRAPH_KAMADA_CUSTOM_INIT = False

# Universal layout constants
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