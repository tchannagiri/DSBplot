import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), './utils/'))) # allow importing the utils dir

import library_constants
import kmer_utils

import plotly.graph_objects as go
import plotly.colors as pc
import pandas as pd
import numpy as np

HOVER_TEMPLATE = '%{hovertext}<extra></extra>' # use pre-formatted hovertext field and disable extra margin

def format_hover_html(data_info, the_dict, format_type, reverse_complement=False):
  def format_name_value(name, value): 
    if isinstance(value, float):
      value = f'{value:.2e}'
    return f'<b>{name}:</b> {value}'

  def create_mismatch_htmls(ref_align, read_align):
    if reverse_complement:
      ref_align = kmer_utils.reverse_complement(ref_align)
      read_align = kmer_utils.reverse_complement(read_align)
    ref_html = ''
    read_html = ''
    for i in range(len(ref_align)):
      if ref_align[i] == read_align[i]:
        ref_html += ref_align[i]
        read_html += read_align[i]
      else:
        ref_html += "<span style='color: red;'><b>{}</b></span>".format(ref_align[i])
        read_html += "<span style='color: red;'><b>{}</b></span>".format(read_align[i])
    return ref_html, read_html
  
  title = None
  if format_type in ['sequence_data']:
    column_names = [
      'id',
      *library_constants.FREQ_COLUMNS[data_info['format']],
      *library_constants.FREQ_RANK_COLUMNS[data_info['format']],
      'variation_type',
      'dist_ref',
      'substitution',
      'insertion',
      'deletion',
      'alignments',
    ]
    title = 'Sequence'
  elif format_type == 'variation':
    column_names = [
      'id',
      'sequence',
      *library_constants.FREQ_COLUMNS[data_info['format']],
      *library_constants.FREQ_RANK_COLUMNS[data_info['format']],
      'variation_type',
      'variation_pos',
      'variation_letter',
      'dist_ref',
      'substitutions',
      'insertions',
      'deletions',
      'alignments',
    ]
    title = 'Sequence Variation'
  elif format_type == 'variation_grouped':
    column_names = [
      'id',
      'var_id',
      'dist_ref',
      'variation_pos',
      'variation_type',
      'variation_letter',
      *library_constants.FREQ_COLUMNS[data_info['format']],
      *library_constants.FREQ_RANK_COLUMNS[data_info['format']],
    ]
    title = 'Sequence Variation'
  elif format_type == 'edge':
    column_names = [
      'id_a',
      'id_b',
      'variation_type_a',
      'variation_type_b',
      'edge_type',
    ]
    title = 'Edge'
  else:
    column_names = list(the_dict.keys())

  format = []
  if title:
    format.append(f'<b>{title}</b><br>')
  for name in column_names:
    if name == 'alignments':
      if format_type in ['sequence_data', 'variation']:
        ref_align, read_align = create_mismatch_htmls(
          the_dict['ref_align'],
          the_dict['read_align'],
        )
        align_html = (
          format_name_value('Ref align ', ref_align) + '<br>' +
          format_name_value('Read align', read_align) + '<br>'
        )
        format.append(align_html)
      elif format_type in ['edges']:
        for suffix in ['a', 'b']:
          ref_align, read_align = create_mismatch_htmls(
            the_dict['ref_align_' + suffix],
            the_dict['read_align_' + suffix],
          )
          align_html = (
            format_name_value('Ref align ' + suffix + ' ', ref_align) + '<br>' +
            format_name_value('Read align ' + suffix, read_align) + '<br><br>'
          )
          format.append(align_html)
    elif name in the_dict:
      value = the_dict[name]
      if isinstance(value, str) and (len(value) > 60):
        value = [value[i : (i + 60)] for i in range(0, len(value), 60)]
        value = '<br>    '.join(value)
      format.append(format_name_value(name, value))
  return '<br>'.join(format)

def format_hover_html_all(
  data_info,
  sequence_data,
  format_type,
  reverse_complement = False,
):
  return pd.Series(
    [
      format_hover_html(
        data_info,
        x,
        format_type,
        reverse_complement = reverse_complement
      )
      for x in sequence_data.to_dict('records')
    ],
    index = sequence_data.index
  )

def get_node_label_text(
  data_info,
  node_data,
  node_label_fields,
  reverse_complement = False,
):
  node_labels = {}
  for id, row_data in node_data.to_dict('index').items():
    node_label_row = []
    for label_type in node_label_fields:
      if label_type in row_data:
        if label_type == 'variation_type':
          label = library_constants.VARIATION_TYPES[row_data[label_type]]['short_label']
        elif library_constants.is_freq_column(label_type):
          label = f'{row_data[label_type]: .2e}'
        elif reverse_complement and (label_type in ['ref_align', 'read_align']):
          label = kmer_utils.reverse_complement(str(row_data[label_type]))
        else:
          label = str(row_data[label_type])
      node_label_row.append(label)
    node_labels[id] = '<br>'.join(node_label_row)
  return pd.Series(node_labels)

def color_palette(index):
  colors = pc.qualitative.Plotly
  return colors[index % len(colors)]

def log_transform_scale(
  x,
  min_x,
  max_x,
  min_out,
  max_out,
):
  x = np.clip(x, min_x, max_x)
  x = (np.log10(x) - np.log10(min_x)) / (np.log10(max_x) - np.log10(min_x))
  return min_out + (max_out - min_out) * x

def log_transform_ratio(
  x1,
  x2,
  min_log_ratio,
  max_log_ratio,
):
  both_zero = (x1 == 0) & (x2 == 0)
  x1_only_zero = (x1 == 0) & ~both_zero
  x2_only_zero = (x2 == 0) & ~both_zero
  neither_zero = (x1 != 0) & (x2 != 0)
  
  log_ratios = np.zeros_like(x1)
  log_ratios[both_zero] = 0
  log_ratios[x1_only_zero] = min_log_ratio
  log_ratios[x2_only_zero] = max_log_ratio
  log_ratios[neither_zero] = np.clip(
    np.log(x1[neither_zero]) - np.log(x2[neither_zero]),
    min_log_ratio,
    max_log_ratio,
  )
  return log_ratios

def get_max_freq(data_info, node_data):
  return node_data[library_constants.FREQ_COLUMNS[data_info['format']]].max(axis='columns')

def get_node_size(
  data_info,
  node_data,
  node_type,
  node_size_type,
  node_size_min_px,
  node_size_max_px,
  node_size_min_freq,
  node_size_max_freq,
):
  if node_size_type == 'freq':
    return pd.Series(
      log_transform_scale(
        get_max_freq(data_info, node_data),
        node_size_min_freq,
        node_size_max_freq,
        node_size_min_px,
        node_size_max_px,
      ),
      index = node_data.index,
    ) 
  else:
    return pd.Series(node_size_min_px, index=node_data.index)

def get_node_freq_group(node_data):
  log_ratio = pd.Series(
    log_transform_ratio(
      node_data['freq_mean_1'],
      node_data['freq_mean_2'],
      library_constants.FREQ_RATIO_COLOR_SCALE_LOG_MIN,
      library_constants.FREQ_RATIO_COLOR_SCALE_LOG_MAX,
    ),
    index = node_data.index,
  )
  node_freq_group = pd.Series(
    library_constants.FREQ_GROUP_B,
    index = node_data.index
  )
  node_freq_group.loc[log_ratio > library_constants.FREQ_RATIO_LOG_A] = (
    library_constants.FREQ_GROUP_A
  )
  node_freq_group.loc[log_ratio < library_constants.FREQ_RATIO_LOG_C] = (
    library_constants.FREQ_GROUP_B
  )

def get_node_color(
  data_info,
  node_data,
  node_color_type,
  node_color_min_freq,
  node_color_max_freq,
):
  if node_color_type == 'freq_group':
    if data_info['format'] != library_constants.DATA_COMPARISON:
      raise Exception('Need a comparison data set: ' + data_info['name'])
    node_freq_group = get_node_freq_group(node_data)
    node_color = pd.Series(library_constants.SIMILAR_FREQ_COLOR, index=node_data.index)
    node_color.loc[node_freq_group == library_constants.FREQ_GROUP_A] = (
      library_constants.CONSTRUCT_COLOR[data_info['construct_1']]
    )
    node_color.loc[node_freq_group == library_constants.FREQ_GROUP_C] = (
      library_constants.CONSTRUCT_COLOR[data_info['construct_2']]
    )
    return node_color
  elif node_color_type == 'freq':
     scaled_freq = log_transform_scale(
       get_max_freq(data_info, node_data),
       node_color_min_freq,
       node_color_max_freq,
       0,
       1,
     )
     return pd.Series(
      pc.sample_colorscale(
        pc.get_colorscale('Inferno'),
        scaled_freq,
      ),
      index = node_data.index,
    )
  elif node_color_type == 'freq_ratio':
    if data_info['format'] != library_constants.DATA_COMPARISON:
      raise Exception('Need a comparison data set: ' + data_info['name'])
    return pd.Series(
      log_transform_ratio(
        node_data['freq_mean_1'],
        node_data['freq_mean_2'],
        library_constants.FREQ_RATIO_COLOR_SCALE_LOG_MIN,
        library_constants.FREQ_RATIO_COLOR_SCALE_LOG_MAX,
      ),
      index = node_data.index,
    )
  elif node_color_type == 'variation_type':
    return node_data['variation_type'].apply(
      lambda x: library_constants.VARIATION_TYPES[x]['color']
    )
  else:
    return pd.Series(library_constants.DEFAULT_NODE_COLOR, index=node_data.index)

def get_node_hover_text(
  data_info,
  node_data,
  node_type,
  reverse_complement = False,
):
  return format_hover_html_all(
    data_info,
    node_data,
    node_type,
    reverse_complement =  reverse_complement,
  )

def get_node_trace_group(node_type, node_data, node_group_type):
  group_key_lists = []
  if node_type in ['sequence_data', 'variation']:
    group_key_lists.append(node_data['is_ref'])
  elif node_type in ['variation_grouped']:
    group_key_lists.append(pd.Series(False, index=node_data.index))
  else:
    raise Exception('Unknown node_type: ' + str(node_type))
  if node_group_type == 'freq_group':
    group_key_lists.append(get_node_freq_group(node_data))
  elif node_group_type == 'variation_type':
    group_key_lists.append(node_data['variation_type'])
  group_keys = list(zip(*group_key_lists))
  return pd.Series(group_keys, index=node_data.index)

def make_node_traces(
  data_info,
  node_data,
  layout,
  node_label,
  node_label_font_size,
  node_hover_text,
  node_size,
  node_color,
  node_color_type,
  node_group,
  node_group_type,
  node_label_position,
  line_width_scale = 1,
):
  traces = []

  for group_key, node_data_group in node_data.groupby(node_group):
    is_ref = group_key[0]
    if is_ref:
      line_color = library_constants.REFERENCE_OUTLINE_COLOR
      line_width = library_constants.REFERENCE_OUTLINE_WIDTH
      trace_name = library_constants.REFERENCE_DESCRIPTION
    else:
      if node_group_type == 'variation_type':
        variation_type = group_key[1]
        trace_name = library_constants.VARIATION_TYPES[variation_type]['label']
      elif node_group_type == 'freq_group':
        freq_group = group_key[1]
        trace_name = library_constants.get_freq_ratio_label(
          freq_group,
          data_info['construct_1'],
          data_info['construct_2'],
        )
      else:
        trace_name = library_constants.NON_REFERENCE_DESCRIPTION
      line_color = library_constants.DEFAULT_OUTLINE_COLOR
      line_width = library_constants.DEFAULT_OUTLINE_WIDTH
      
    trace_args = {
      'name': trace_name,
      'x': layout.loc[node_data_group.index, 0],
      'y': layout.loc[node_data_group.index, 1],
      'mode': 'markers+text',
      'text': node_label[node_data_group.index],
      'textposition': node_label_position,
      'textfont_size': node_label_font_size,
      'marker': {
        'color': node_color[node_data_group.index],
        'size': node_size[node_data_group.index],
        'opacity': 1.0,
        'line': {
          'width': line_width * line_width_scale,
          'color': line_color,
        }
      },
      'hovertext': node_hover_text[node_data_group.index],
      'hovertemplate': HOVER_TEMPLATE,
    }
    traces.append(go.Scatter(**trace_args))
  return traces

def make_point_traces(
  data_info,
  node_data,
  graph_layout,
  show_node_labels,
  node_label_columns,
  node_label_position,
  node_label_font_size,
  node_type,
  node_color_type,
  node_size_type,
  node_size_min_px,
  node_size_max_px,
  node_size_min_freq,
  node_size_max_freq,
  node_outline_width_scale = 1,
  reverse_complement = False,
):
  if show_node_labels:
    node_label = get_node_label_text(
      data_info,
      node_data,
      node_type,
      node_label_columns,
      reverse_complement = reverse_complement,
    )
  else:
    node_label = pd.Series('', index=node_data.index)
  hover_text = get_node_hover_text(
    data_info,
    node_data,
    node_type,
    reverse_complement =  reverse_complement,
  )

  node_size = get_node_size(
    data_info = data_info,
    node_data = node_data,
    node_type = node_type,
    node_size_type = node_size_type,
    node_size_min_px = node_size_min_px,
    node_size_max_px = node_size_max_px,
    node_size_min_freq = node_size_min_freq,
    node_size_max_freq = node_size_max_freq,
  )

  node_color = get_node_color(
    data_info = data_info,
    node_data = node_data,
    node_color_type = node_color_type,
    node_color_min_freq = node_size_min_freq,
    node_color_max_freq = node_size_max_freq,
  )

  node_group = get_node_trace_group(
    node_data = node_data,
    node_type = node_type,
    node_group_type = node_color_type,
  )
  
  traces = make_node_traces(
    data_info = data_info,
    node_data = node_data,
    layout = graph_layout,
    node_label = node_label,
    node_label_font_size = node_label_font_size,
    node_hover_text = hover_text,
    node_size = node_size,
    node_color = node_color,
    node_color_type = node_color_type,
    node_label_position = node_label_position,
    node_group = node_group,
    node_group_type = node_color_type,
    line_width_scale = node_outline_width_scale,
  )
  return traces

def make_edges_traces(
  data_info,
  graph,
  layout,
  show_edge_labels,
  show_edge_types,
  edge_width_scale = 1,
  reverse_complement = False,
):
  edge_args = {}
  for edge_type in show_edge_types:
    edge_args[edge_type] = dict(
      x = [],
      y = [],
      name = library_constants.EDGE_TYPES[edge_type]['label'],
      line = dict(
        dash = library_constants.EDGE_TYPES[edge_type]['line_dash'],
        color = library_constants.EDGE_TYPES[edge_type]['plot_color'],
        width = edge_width_scale,
      ),
      mode = 'lines',
      showlegend = True,
    )

  label_args = dict(
    x = [],
    y = [],
    text = [],
    hovertext = [],
    hovertemplate = HOVER_TEMPLATE,
    mode = 'text',
    showlegend = False,
  )
      
  for id_a, id_b in graph.edges():
    edge_type = graph.edges[id_a, id_b]['edge_type']

    if edge_type in show_edge_types:
      x1 = layout.loc[id_a, 0]
      x2 = layout.loc[id_b, 0]
      y1 = layout.loc[id_a, 1]
      y2 = layout.loc[id_b, 1]

      edge_args[edge_type]['x'] += [x1, x2, None]
      edge_args[edge_type]['y'] += [y1, y2, None]

      label_args['x'].append((x1 + x2) / 2)
      label_args['y'].append((y1 + y2) / 2)
      label_args['hovertext'].append(
        format_hover_html(
          data_info,
          graph.edges[id_a, id_b],
          'edge',
          reverse_complement = reverse_complement,
        )
      )
      if show_edge_labels:
        label_args['text'].append(edge_type[0].upper())
  
  traces = []
  for args in edge_args.values():
    traces.append(go.Scatter(**args))
  traces.append(go.Scatter(**label_args))
  return traces
