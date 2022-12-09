import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import library_constants
import generate_constants
import generate_07_plot_graph

def get_output_file(cell_line, dsb_type, version):
  version_str = (
    ''
    if (version == library_constants.VERSION_NONE) else
    ('_' + version)
  )
  return generate_constants.join_path(
    [
      generate_constants.get_output_dir('pptx'),
      'graph',
      cell_line + '_' + dsb_type + version_str + os.extsep + 'pptx',
    ]
  )

def get_input_file(experiment_name, layout_name, format, file_ext):
  return generate_constants.join_path(
    [
      generate_07_plot_graph.get_output_dir(
        layout_name,
        format,
        file_ext = file_ext
      ),
      experiment_name + os.path.extsep + file_ext,
    ]
  )

TOTAL_WIDTH = {
  library_constants.DATA_INDIVIDUAL: 1,
  library_constants.DATA_COMPARISON: 0.95,
}
ARG_LEGEND = '--legends node_size freq_ratio_sense_branch freq_ratio_sense_cmv node_type edge_type'

if __name__ == '__main__':
  for ext in ['sh', 'ps1']:
    with open(
      file = 'run_09_pptx_graph' + os.extsep + ext,
      mode = 'w',
      encoding = generate_constants.OUTPUT_ENCODING[ext],
    ) as file_out:
      for dsb_type in library_constants.DSB_TYPES:
        if dsb_type == library_constants.DSB_TYPE_2anti:
          cell_line_list = [library_constants.CELL_LINE_WT]
        else:
          cell_line_list = library_constants.CELL_LINES
        for cell_line in cell_line_list:
          if dsb_type == library_constants.DSB_TYPE_1:
            constructs_individual = library_constants.CONSTRUCTS_INDIVIDUAL_SENSE
            constructs_comparison = library_constants.CONSTRUCTS_COMPARISON_SENSE
            version_list = [library_constants.VERSION_NONE]
          elif dsb_type == library_constants.DSB_TYPE_2:
            constructs_individual = library_constants.CONSTRUCTS_INDIVIDUAL_SENSE
            constructs_comparison = library_constants.CONSTRUCTS_COMPARISON_SENSE
            version_list = [library_constants.VERSION_NONE]
          elif dsb_type == library_constants.DSB_TYPE_2anti:
            constructs_individual = library_constants.CONSTRUCTS_INDIVIDUAL_ANTISENSE
            constructs_comparison = library_constants.CONSTRUCTS_COMPARISON_ANTISENSE
            version_list = [
              library_constants.VERSION_OLD,
              library_constants.VERSION_NEW,
              library_constants.VERSION_MERGED,
            ]
          else:
            raise Exception('Unknown dsb_type: ' + dsb_type)
          
          for version in version_list:
            file_list = []
            label_list = []
            num_grids = 0
            num_rows_list = []
            num_cols_list = []
            total_width_list = []
            for format in [
              library_constants.DATA_INDIVIDUAL,
              library_constants.DATA_COMPARISON,
            ]:
              num_grids += 1
              total_width_list.append(TOTAL_WIDTH[format])
              if format == library_constants.DATA_INDIVIDUAL:
                experiment_info = generate_constants.EXPERIMENT_INFO
                construct_list = constructs_individual
              elif format == library_constants.DATA_COMPARISON:
                experiment_info = generate_constants.EXPERIMENT_INFO_COMPARISON
                construct_list = constructs_comparison
              else:
                raise Exception('Impossible.')
              num_rows_list.append(len(library_constants.STRANDS))
              num_cols_list.append(len(construct_list))
              experiment_info = experiment_info.loc[
                (experiment_info['cell_line'] == cell_line) &
                (experiment_info['dsb_type'] == dsb_type) &
                (experiment_info['version'] == version) &
                (experiment_info['control_type'] == library_constants.CONTROL_NOT)
              ]
              for strand in library_constants.STRANDS:
                for construct in construct_list:
                  info = experiment_info.loc[
                    (experiment_info['strand'] == strand) &
                    (experiment_info['construct'] == construct)
                  ]
                  if info.shape[0] != 1:
                    raise Exception(f'Got {info.shape[0]} rows. Expected 1.')
                  info = info.iloc[0].to_dict()
                  file_list.append(get_input_file(
                    experiment_name = info['name'],
                    layout_name = generate_constants.USE_LAYOUT,
                    format = format,
                    file_ext = 'png',
                  ))
                  label_list.append(
                    library_constants.LABELS[info['guide_rna']] +
                    generate_constants.ARG_NEWLINE[ext] +
                    library_constants.LABELS[strand] +
                    generate_constants.ARG_NEWLINE[ext] +
                    library_constants.LABELS[construct]
                  )
            if dsb_type == library_constants.DSB_TYPE_2anti:
              # transpose the 2'nd grid and make full width
              num_rows_list[-1], num_cols_list[-1] = num_cols_list[-1], num_rows_list[-1]
              total_width_list[-1] = 1
            arg_input = '--input ' + ' '.join(file_list)
            arg_labels = '--labels ' + ' '.join(f'"{x}"' for x in label_list)
            arg_num_grids = '--num_grids ' + str(num_grids)
            arg_num_rows = '--num_rows ' + ' '.join(str(x) for x in num_rows_list)
            arg_num_cols = '--num_cols ' + ' '.join(str(x) for x in num_cols_list)
            arg_total_width = '--total_width ' + ' '.join(str(x) for x in total_width_list)
            arg_output = '--output ' + get_output_file(cell_line, dsb_type, info['version'])
            file_out.write(f"python {generate_constants.get_python_script('get_pptx')} {arg_input} {arg_output} {arg_labels} {arg_num_grids} {arg_num_rows} {arg_num_cols} {arg_total_width} {ARG_LEGEND}\n")
      log_utils.log(file_out.name)

