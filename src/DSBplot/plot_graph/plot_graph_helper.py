import plotly.graph_objects as go
import plotly.colors as pc
import pandas as pd
import numpy as np

import DSBplot.utils.constants as constants
import DSBplot.utils.kmer_utils as kmer_utils

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
  if format_type == 'node':
    if data_info['format'] == 'comparison':
      freq_columns = ['freq_mean_1', 'freq_mean_2']
    elif data_info['format'] == 'individual':
      freq_columns = ['freq_mean']
    else:
      raise Exception(f'Invalid format: {data_info["format"]}')
    column_names = [
      'id',
      *freq_columns,
      'var_type',
      'dist_ref',
      'sub',
      'ins',
      'del',
      'alignments',
    ]
    title = 'Sequence'
  elif format_type == 'edge':
    column_names = [
      'id_a',
      'id_b',
      'var_type_a',
      'var_type_b',
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
      if format_type == 'node':
        ref_align, read_align = create_mismatch_htmls(
          the_dict['ref_align'],
          the_dict['read_align'],
        )
        align_html = (
          format_name_value('ref_align ', ref_align) + '<br>' +
          format_name_value('read_align', read_align) + '<br>'
        )
        format.append(align_html)
      elif format_type == 'edge':
        for suffix in ['a', 'b']:
          ref_align, read_align = create_mismatch_htmls(
            the_dict['ref_align_' + suffix],
            the_dict['read_align_' + suffix],
          )
          align_html = (
            format_name_value('ref_align ' + suffix + ' ', ref_align) + '<br>' +
            format_name_value('read_align ' + suffix, read_align) + '<br><br>'
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
  data,
  format_type,
  reverse_complement = False,
):
  return pd.Series(
    [
      format_hover_html(
        data_info = data_info,
        the_dict = x,
        format_type = format_type,
        reverse_complement = reverse_complement
      )
      for x in data.to_dict('records')
    ],
    index = data.index
  )

def get_node_label_text(
  node_data,
  node_label_columns,
  reverse_complement = False,
):
  node_labels = {}
  for id, row_data in node_data.to_dict('index').items():
    node_label_row = []
    for label_type in node_label_columns:
      if label_type in row_data:
        if label_type == 'var_type':
          label = constants.VARIATION_TYPES[row_data[label_type]]['short_label']
        elif constants.is_freq_column(label_type):
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

def get_node_size(
  data_info,
  node_data,
  node_size_type,
  node_size_px_range,
  node_size_freq_range,
):
  if node_size_type == 'freq':
    return pd.Series(
      log_transform_scale(
        node_data['freq_mean'],
        node_size_freq_range[0],
        node_size_freq_range[1],
        node_size_px_range[0],
        node_size_px_range[1],
      ),
      index = node_data.index,
    ) 
  else:
    return pd.Series(node_size_px_range[0], index=node_data.index)

def get_node_freq_group(node_data, node_freq_ratio_range):
  log_ratio = pd.Series(
    log_transform_ratio(
      node_data['freq_mean_1'],
      node_data['freq_mean_2'],
      node_freq_ratio_range[0],
      node_freq_ratio_range[1],
    ),
    index = node_data.index,
  )
  node_freq_group = pd.Series('B', index=node_data.index)
  node_freq_group.loc[log_ratio > np.log(node_freq_ratio_range[1])] = 'A'
  node_freq_group.loc[log_ratio < np.log(node_freq_ratio_range[0])] = 'C'
  return node_freq_group

def get_node_color(
  data_info,
  node_data,
  node_color_type,
  node_color_freq_range,
  node_comparison_colors,
  node_freq_ratio_range,
  node_fill_color,
  node_var_type_colors,
):
  if node_color_type == 'ratio_disc':
    if data_info['format'] != 'comparison':
      raise Exception('Need a comparison data set: ' + data_info['label'])
    node_freq_group = get_node_freq_group(node_data, node_freq_ratio_range)
    node_color = pd.Series(constants.SIMILAR_FREQ_COLOR, index=node_data.index)
    node_color.loc[node_freq_group == 'A'] = node_comparison_colors[0]
    node_color.loc[node_freq_group == 'C'] = node_comparison_colors[1]
    return node_color
  elif node_color_type == 'freq':
     scaled_freq = log_transform_scale(
       x = node_data['freq_mean'],
       min_x = node_color_freq_range[0],
       max_x = node_color_freq_range[1],
       min_out = 0,
       max_out = 1,
     )
     return pd.Series(
      pc.sample_colorscale(
        pc.get_colorscale('Inferno'),
        scaled_freq,
      ),
      index = node_data.index,
    )
  elif node_color_type == 'ratio_cont':
    if data_info['format'] != 'comparison':
      raise Exception('Need a comparison data set: ' + data_info['label'])
    return pd.Series(
      log_transform_ratio(
        node_data['freq_mean_1'],
        node_data['freq_mean_2'],
        np.log(node_freq_ratio_range[0]),
        np.log(node_freq_ratio_range[1]),
      ),
      index = node_data.index,
    )
  elif node_color_type == 'var_type':
    return node_data['var_type'].apply(lambda x: node_var_type_colors[x])
  else:
    return pd.Series(node_fill_color, index=node_data.index)

def get_node_trace_group(node_data, node_group_type, node_freq_ratio_range):
  group_key_lists = []
  group_key_lists.append(node_data['is_ref'])
  if node_group_type == 'ratio_disc':
    group_key_lists.append(get_node_freq_group(node_data, node_freq_ratio_range))
  elif node_group_type == 'var_type':
    group_key_lists.append(node_data['var_type'])
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
  node_reference_outline_color,
  node_outline_color,
  node_group,
  node_group_type,
  node_label_position,
  node_freq_ratio_range,
  line_width_scale = 1,
):
  traces = []

  for group_key, node_data_group in node_data.groupby(node_group):
    is_ref = group_key[0]
    if is_ref:
      line_color = node_reference_outline_color
      line_width = constants.GRAPH_NODE_REFERENCE_OUTLINE_WIDTH
      trace_name = constants.LABEL_REFERENCE
    else:
      if node_group_type == 'var_type':
        var_type = group_key[1]
        trace_name = constants.VARIATION_TYPES[var_type]['label']
      elif node_group_type == 'ratio_disc':
        freq_group = group_key[1]
        trace_name = constants.get_freq_ratio_label(
          freq_group,
          data_info['label_1'],
          data_info['label_2'],
          node_freq_ratio_range[0],
          node_freq_ratio_range[1],
        )
      else:
        trace_name = constants.LABEL_NONREFERENCE
      line_color = node_outline_color
      line_width = constants.GRAPH_NODE_OUTLINE_WIDTH
      
    trace_args = {
      'name': trace_name,
      'x': layout.loc[node_data_group.index, 'x'],
      'y': layout.loc[node_data_group.index, 'y'],
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
  node_color_type,
  node_comparison_colors,
  node_freq_ratio_range,
  node_reference_outline_color,
  node_outline_color,
  node_fill_color,
  node_var_type_colors,
  node_size_type,
  node_size_px_range,
  node_size_freq_range,
  node_outline_width_scale = 1,
  reverse_complement = False,
):
  if show_node_labels:
    node_label = get_node_label_text(
      node_data = node_data,
      node_label_columns = node_label_columns,
      reverse_complement = reverse_complement,
    )
  else:
    node_label = pd.Series('', index=node_data.index)
  hover_text = format_hover_html_all(
    data_info = data_info,
    data = node_data,
    format_type = 'node',
    reverse_complement =  reverse_complement,
  )

  node_size = get_node_size(
    data_info = data_info,
    node_data = node_data,
    node_size_type = node_size_type,
    node_size_px_range = node_size_px_range,
    node_size_freq_range = node_size_freq_range,
  )

  node_color = get_node_color(
    data_info = data_info,
    node_data = node_data,
    node_color_type = node_color_type,
    node_color_freq_range = node_size_freq_range,
    node_comparison_colors = node_comparison_colors,
    node_freq_ratio_range = node_freq_ratio_range,
    node_fill_color = node_fill_color,
    node_var_type_colors = node_var_type_colors,
  )

  node_group = get_node_trace_group(
    node_data = node_data,
    node_group_type = node_color_type,
    node_freq_ratio_range = node_freq_ratio_range,
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
    node_reference_outline_color = node_reference_outline_color,
    node_outline_color = node_outline_color,
    node_label_position = node_label_position,
    node_group = node_group,
    node_group_type = node_color_type,
    node_freq_ratio_range = node_freq_ratio_range,
    line_width_scale = node_outline_width_scale,
  )

  return traces

def make_edges_traces(
  data_info,
  graph,
  layout,
  show_edge_labels,
  show_edge_types,
  edge_width_scale = constants.GRAPH_EDGE_WIDTH_SCALE,
  reverse_complement = False,
):
  edge_args = {}
  for edge_type in show_edge_types:
    edge_args[edge_type] = dict(
      x = [],
      y = [],
      name = constants.EDGE_TYPES[edge_type]['label'],
      line = dict(
        dash = constants.EDGE_TYPES[edge_type]['line_dash'],
        color = constants.EDGE_TYPES[edge_type]['plot_color'],
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
      x1 = layout.loc[id_a, 'x']
      x2 = layout.loc[id_b, 'x']
      y1 = layout.loc[id_a, 'y']
      y2 = layout.loc[id_b, 'y']

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
