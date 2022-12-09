import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import pandas as pd
import file_utils
import library_constants

def get_name_library(info):
  name = (
    info['library'] +
    '_' +
    library_constants.get_data_label(info)
  )
  for version in library_constants.VERSIONS:
    name = name.replace('_' + version, '')
  return name

def get_name_experiment(info):
  return library_constants.get_data_label(info)

def get_ref_seq_file(info):
  return (
    info['dsb_type'] +
    '_' + info['strand'] +
    '_' + info['construct'] +
    (('_' + str(info['version'])) if (info['version'] != library_constants.VERSION_NONE) else '') +
    os.path.extsep + 'fa'
  )

DSB_POS = file_utils.read_tsv(os.path.dirname(__file__) + '/dsb_pos.tsv')

def get_dsb_pos(info):
  dsb_pos = DSB_POS.loc[
    (DSB_POS['dsb_type'] == info['dsb_type']) &
    (DSB_POS['strand'] == info['strand']) &
    (DSB_POS['version'] == info['version'])
  ]['dsb_pos']
  if dsb_pos.shape[0] != 1:
    raise Exception(f'Got {dsb_pos.shape[0]} values. Expected 1.')
  dsb_pos = dsb_pos.iloc[0]
  if info['control_type'] == '30bpDown':
    dsb_pos += 30
  return dsb_pos

TOTAL_READS = file_utils.read_tsv(os.path.dirname(__file__) + '/total_reads.tsv')
def get_total_reads(info):
  x = TOTAL_READS.loc[
    (TOTAL_READS['library'] == info['library']),
    info['strand']  
  ]
  if x.shape[0] != 1:
    raise Exception(f'Got {x.shape[0]} values. Expected 1.')
  return x.iloc[0]

MIN_READ_LENGTH = file_utils.read_tsv(os.path.dirname(__file__) + '/min_read_length.tsv')
def get_min_read_length(info):
  x = MIN_READ_LENGTH.loc[
    (MIN_READ_LENGTH['dsb_type'] == info['dsb_type']),
    'min_read_length'
  ]
  if x.shape[0] != 1:
    raise Exception(f'Got {x.shape[0]} values. Expected 1.')
  return x.iloc[0]

LIBRARY_INFO = file_utils.read_tsv(os.path.dirname(__file__) + '/library_info.tsv')
LIBRARY_INFO['format'] = library_constants.DATA_INDIVIDUAL
LIBRARY_INFO['name'] = LIBRARY_INFO.apply(get_name_library, axis='columns')
LIBRARY_INFO['ref_seq_file'] = LIBRARY_INFO.apply(get_ref_seq_file, axis='columns')
LIBRARY_INFO['dsb_pos'] = LIBRARY_INFO.apply(get_dsb_pos, axis='columns')
LIBRARY_INFO['total_reads'] = LIBRARY_INFO.apply(get_total_reads, axis='columns')
LIBRARY_INFO['min_read_length'] = LIBRARY_INFO.apply(get_min_read_length, axis='columns')

def get_library_info(**args):
  library_info = LIBRARY_INFO
  for key in args:
    if args[key] is not None:
      library_info = library_info.loc[library_info[key] == args[key]]
  if library_info.shape[0] != 1:
    raise Exception(f'Got {library_info.shape[0]} rows. Expected 1.')
  return library_info.iloc[0].to_dict()

ANTISENSE_MERGED_PAIRS = {
  ('yjl89', 'yjl349'): 'yjl89n349',
  ('yjl90', 'yjl350'): 'yjl90n350',
  ('yjl91', 'yjl351'): 'yjl91n351',
  ('yjl92', 'yjl352'): 'yjl92n352',
  ('yjl93', 'yjl353'): 'yjl93n353',
  ('yjl94', 'yjl354'): 'yjl94n354',
  ('yjl95', 'yjl355'): 'yjl95n355',
  ('yjl96', 'yjl356'): 'yjl96n356',
}
def get_library_info_antisense_merged():
  info_list = []
  for (library_1, library_2), library_new in ANTISENSE_MERGED_PAIRS.items():
    for strand in library_constants.STRANDS:
      info_1 = get_library_info(library = library_1, strand = strand)
      info_2 = get_library_info(library = library_2, strand = strand)
      info_merged = info_1.copy()
      info_merged['total_reads'] += info_2['total_reads']
      info_merged['library'] = library_new
      info_merged['version'] = library_constants.VERSION_MERGED
      info_merged['name'] = get_name_library(info_merged)
      info_merged['dsb_pos'] = None
      info_merged['ref_seq_file'] = None
      info_list.append(info_merged)
  return pd.DataFrame.from_records(info_list)
LIBRARY_INFO_ANTISENSE_MERGED = get_library_info_antisense_merged()
LIBRARY_INFO = pd.concat(
  [LIBRARY_INFO, LIBRARY_INFO_ANTISENSE_MERGED],
  axis = 'index'
).reset_index(drop=True)
LIBRARY_INFO['name_experiment'] = LIBRARY_INFO.apply(
  get_name_experiment, axis = 'columns'
)

EXPERIMENT_INFO = LIBRARY_INFO.groupby([
  'cell_line',
  'control_type',
  'dsb_type',
  'guide_rna',
  'construct',
  'strand',
  'version',
]).aggregate(
  library_list = ('library', list),
  library_name_list = ('name', list),
  dsb_pos = ('dsb_pos', 'first'),
  ref_seq_file = ('ref_seq_file', 'first'),
  total_reads_list = ('total_reads', list),
  min_read_length = ('min_read_length', 'first'),
).reset_index()
EXPERIMENT_INFO['format'] = library_constants.DATA_INDIVIDUAL
EXPERIMENT_INFO['name'] = EXPERIMENT_INFO.apply(get_name_experiment, axis='columns')

LAYOUT_GROUP_2DSB = '2DSB'
LAYOUT_GROUP_1DSB_A = '1DSB_A'
LAYOUT_GROUP_1DSB_B = '1DSB_B'
LAYOUT_GROUP_2DSBanti = '2DSBanti'
LAYOUT_GROUPS = [
  LAYOUT_GROUP_2DSB,
  LAYOUT_GROUP_1DSB_A,
  LAYOUT_GROUP_1DSB_B,
  LAYOUT_GROUP_2DSBanti,
]

EXPERIMENT_INFO['layout_group'] = None
EXPERIMENT_INFO.loc[
  EXPERIMENT_INFO['dsb_type'] == library_constants.DSB_TYPE_2,
  'layout_group'
] = LAYOUT_GROUP_2DSB
EXPERIMENT_INFO.loc[
  (
    (EXPERIMENT_INFO['dsb_type'] == library_constants.DSB_TYPE_1) &
    (EXPERIMENT_INFO['guide_rna'] == library_constants.GUIDE_RNA_A)
  ),
  'layout_group'
] = LAYOUT_GROUP_1DSB_A
EXPERIMENT_INFO.loc[
  (
    (EXPERIMENT_INFO['dsb_type'] == library_constants.DSB_TYPE_1) &
    (EXPERIMENT_INFO['guide_rna'] == library_constants.GUIDE_RNA_B)
  ),
  'layout_group'
] = LAYOUT_GROUP_1DSB_B
EXPERIMENT_INFO.loc[
  EXPERIMENT_INFO['dsb_type'] == library_constants.DSB_TYPE_2anti,
  'layout_group'
] = LAYOUT_GROUP_2DSBanti
EXPERIMENT_INFO['format'] = library_constants.DATA_INDIVIDUAL

def get_experiment_info(**args):
  info = EXPERIMENT_INFO
  for key in args:
    if args[key] is not None:
      info = info.loc[info[key] == args[key]]
  if info.shape[0] != 1:
    raise Exception(f'Got {info.shape[0]} rows. Expected 1.')
  return info.iloc[0].to_dict()

# Make the comparison experiments
def get_experiment_info_comparison():
  experiments_comparison = []
  experiments_comparison_keys = [
    'cell_line',
    'control_type',
    'dsb_type',
    'guide_rna',
    'strand',
    'version',
    'layout_group',
  ]
  for key, experiments in EXPERIMENT_INFO.groupby(experiments_comparison_keys):
    key_dict = dict(zip(experiments_comparison_keys, key))
    if key_dict['control_type'] == library_constants.CONTROL_NOT:
      if key_dict['dsb_type'] == library_constants.DSB_TYPE_2anti:
        construct_1_list = [library_constants.CONSTRUCT_ANTISENSE]
        construct_2_list = [library_constants.CONSTRUCT_SPLICING]
      else:
        construct_1_list = [library_constants.CONSTRUCT_SENSE]
        construct_2_list = [library_constants.CONSTRUCT_BRANCH, library_constants.CONSTRUCT_CMV]
      for construct_1 in construct_1_list:
        for construct_2 in construct_2_list:
          experiment_1 = experiments.loc[experiments['construct'] == construct_1]
          experiment_2 = experiments.loc[experiments['construct'] == construct_2]
          if experiment_1.shape[0] != 1:
            raise Exception(f'Got {experiment_1.shape[0]} rows. Expected 1.')
          if experiment_2.shape[0] != 1:
            raise Exception(f'Got {experiment_2.shape[0]} rows. Expected 1.')
          experiment_1 = experiment_1.iloc[0].to_dict()
          experiment_2 = experiment_2.iloc[0].to_dict()
          experiment_new = {
            k: v for k, v in experiment_1.items()
            if k in experiments_comparison_keys
          }
          experiment_new['construct_1'] = construct_1
          experiment_new['construct_2'] = construct_2
          experiment_new['construct'] = construct_1 + '_' + construct_2
          experiment_new['name_1'] = experiment_1['name']
          experiment_new['name_2'] = experiment_2['name']
          experiment_new['format'] = library_constants.DATA_COMPARISON
          experiment_new['name'] = get_name_experiment(experiment_new)
          experiments_comparison.append(experiment_new)
  return pd.DataFrame.from_records(experiments_comparison)

EXPERIMENT_INFO_COMPARISON = get_experiment_info_comparison()

def join_path(paths):
  return '/'.join(paths)

def rejoin_path(path):
  return join_path(path.split(os.path.sep))

REF_SEQ_DIR = 'ref_seq'
OUTPUT_DIR = {
  'fastq': 'data_fastq',
  'bowtie2_build': 'data_bowtie2_build',
  'sam': 'data_0_sam',
  'filter_nhej': 'data_1_filter_nhej',
  'combine_repeat': 'data_2_combine_repeat',
  'window': 'data_3_window',
  'graph': 'data_4_graph',
  'histogram': 'data_5_histogram',
  'layout': 'data_6_precomputed_layout',
  'plot_graph': 'plot/graph',
  'plot_histogram': 'plot/histogram',
  'pptx': 'pptx',
  'library_summary': 'data_7_library_summary',
}

def get_output_dir(key):
  return rejoin_path(OUTPUT_DIR[key])

OUTPUT_ENCODING = {
  'sh': 'utf-8',
  'ps1': 'utf_8_sig',
}

ARG_NEWLINE = {
  'sh': '\\n',
  'ps1': '`n',
}

PYTHON_SCRIPTS = {
  'filter_nhej': '1_process_nhej/filter_nhej.py',
  'combine_repeat': '1_process_nhej/combine_repeat.py',
  'get_window': '2_get_window_data/get_window.py',
  'get_merged': '2_get_window_data/get_merged.py',
  'get_freq': '2_get_window_data/get_freq.py',
  'get_freq_comparison': '2_get_window_data/get_freq_comparison.py',
  'get_graph_data': '3_get_graph_data/get_graph_data.py',
  'get_histogram_data': '4_get_histogram_data/get_histogram_data.py',
  'get_precomputed_layout': '5_plot_graph/get_precomputed_layout.py',
  'plot_graph': '5_plot_graph/plot_graph.py',
  'plot_histogram': '6_plot_histogram/plot_histogram.py',
  'get_pptx': '7_get_pptx/get_pptx.py',
}

BOWTIE2_BUILD_COMMAND = {
  'ps1': 'bowtie2-build-s.exe',
  'sh': 'bowtie2-build',
}

BOWTIE2_ALIGN_COMMAND = {
  'ps1': 'bowtie2-align-s.exe',
  'sh': 'bowtie2',
}

def get_python_script(key):
  return rejoin_path(PYTHON_SCRIPTS[key])

def check_scripts_exits():
  for x in PYTHON_SCRIPTS.values():
    if not os.path.exists(x):
      raise Exception('Could not find script: ' + str(x))
check_scripts_exits()

LIBRARY_INFO.to_excel('library_info.xlsx')
EXPERIMENT_INFO.to_excel('experiment_info.xlsx')
EXPERIMENT_INFO_COMPARISON.to_excel('experiment_info_comparison.xlsx')

USE_PRECOMPUTED_LAYOUT = True
LAYOUT_UNIVERSAL = 'universal'
LAYOUT_FRACTAL = 'fractal'
LAYOUT_RADIAL = 'radial'
LAYOUTS = [LAYOUT_UNIVERSAL, LAYOUT_FRACTAL, LAYOUT_RADIAL]
USE_LAYOUT = LAYOUT_UNIVERSAL

GRAPH_HEIGHT_PX = 1800
GRAPH_WIDTH_PX = 2400

