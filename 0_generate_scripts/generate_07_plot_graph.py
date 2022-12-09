import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import generate_constants
import library_constants
import generate_04_graph_data
import generate_06_precomputed_layout

def get_input_dir(name):
  return generate_04_graph_data.get_output_dir(name)

def get_layout_dir(layout_name, layout_group):
  return generate_06_precomputed_layout.get_output_dir(layout_name, layout_group)

def get_output_dir(layout_name, format, file_ext):
  return generate_constants.join_path(
    [
      generate_constants.get_output_dir('plot_graph'),
      layout_name,
      format,
      file_ext
    ]
  )

if __name__ == '__main__':
  for script_ext in ['sh', 'ps1']:
    for output_ext in ['png', 'html']:
      with open(
        file = os.path.join('run_07_plot_graph_' + output_ext + os.path.extsep + script_ext),
        mode ='w',
        encoding = generate_constants.OUTPUT_ENCODING[script_ext],
      ) as file_out:
        for info in (
          generate_constants.EXPERIMENT_INFO.to_dict('records') +
          generate_constants.EXPERIMENT_INFO_COMPARISON.to_dict('records')
        ):
          if info['control_type'] == library_constants.CONTROL_NOT:
            input_dir = get_input_dir(info['name'])
            output_dir = get_output_dir(
              generate_constants.USE_LAYOUT,
              info['format'],
              file_ext = output_ext
            )
            if generate_constants.USE_LAYOUT == generate_constants.LAYOUT_UNIVERSAL:
              range_x = {
                generate_constants.LAYOUT_GROUP_2DSB: [-12, 13],
                generate_constants.LAYOUT_GROUP_1DSB_A: [-12, 13],
                generate_constants.LAYOUT_GROUP_1DSB_B: [-12, 13],
                generate_constants.LAYOUT_GROUP_2DSBanti: [-12, 13],
              }[info['layout_group']]
              range_y = {
                generate_constants.LAYOUT_GROUP_2DSB: [-22, 22],
                generate_constants.LAYOUT_GROUP_1DSB_A: [-23, 20],
                generate_constants.LAYOUT_GROUP_1DSB_B: [-22, 16],
                generate_constants.LAYOUT_GROUP_2DSBanti: [-22, 27],
              }[info['layout_group']]
              arg_range_x = '--range_x ' + ' '.join(str(x) for x in range_x)
              arg_range_y = '--range_y ' + ' '.join(str(y) for y in range_y)
            else:
              range_x = None
              range_y = None
              arg_range_x = ''
              arg_range_y = ''
            if generate_constants.USE_PRECOMPUTED_LAYOUT:
              arg_precomputed_layout_dir = (
                '--precomputed_layout_dir ' +
                generate_06_precomputed_layout.get_output_dir(
                  generate_constants.USE_LAYOUT,
                  info['layout_group']
                )
              )
            else:
              arg_precomputed_layout_dir = ''
            arg_layout = '--layout ' + generate_constants.USE_LAYOUT
            arg_reverse_complement = (
              (info['strand'] == library_constants.STRAND_R2) and
              (
                (generate_constants.USE_LAYOUT == generate_constants.LAYOUT_UNIVERSAL) or
                (generate_constants.USE_LAYOUT == generate_constants.LAYOUT_FRACTAL) or
                (
                  (generate_constants.USE_LAYOUT == generate_constants.LAYOUT_RADIAL) and
                  (info['dsb_type'] != library_constants.DSB_TYPE_1)
                )
              )
            )
            arg_reverse_complement = '--reverse_complement' if arg_reverse_complement else ''
            arg_width_height = (
              '--width ' +
              str(generate_constants.GRAPH_WIDTH_PX) +
              ' --height ' +
              str(generate_constants.GRAPH_HEIGHT_PX)
            )
            if generate_constants.USE_LAYOUT == generate_constants.LAYOUT_UNIVERSAL:
              arg_universal_layout_axis_pos = (
                f'--universal_layout_y_axis_x_pos {range_x[1] - 1} ' +
                f'--universal_layout_y_axis_y_range {range_y[0] + 2.5} {range_y[1] - 1.5}'
              )
              if info['name'] in [
                'WT_sgAB_R1_sense',
                'WT_sgA_R1_sense',
                'KO_sgAB_R1_sense',
                'KO_sgA_R1_sense',
                'WT_sgCD_R1_antisense_old',
                'WT_sgCD_R1_antisense_new',
                'WT_sgCD_R1_antisense_merged',
              ]:
                arg_universal_layout_axis_pos += (
                  f' --universal_layout_x_axis_deletion_y_pos {range_y[0] + 1.5}' +
                  f' --universal_layout_x_axis_insertion_y_pos {range_y[1] - 0.5}' +
                  f' --universal_layout_x_axis_x_range {range_x[0] + 0.5} {range_x[1] - 1.5}'
                )
            else:
              arg_universal_layout_axis_pos = ''
            if generate_constants.USE_LAYOUT == generate_constants.LAYOUT_UNIVERSAL:
              arg_universal_layout_max_tick_insertion = {
                generate_constants.LAYOUT_GROUP_2DSB: 7,
                generate_constants.LAYOUT_GROUP_1DSB_A: 6,
                generate_constants.LAYOUT_GROUP_1DSB_B: 5,
                generate_constants.LAYOUT_GROUP_2DSBanti: 8,
              }[info['layout_group']]
              arg_universal_layout_max_tick_deletion = {
                generate_constants.LAYOUT_GROUP_2DSB: 17,
                generate_constants.LAYOUT_GROUP_1DSB_A: 19,
                generate_constants.LAYOUT_GROUP_1DSB_B: 18,
                generate_constants.LAYOUT_GROUP_2DSBanti: 17,
              }[info['layout_group']]
              arg_universal_layout_max_tick = (
                '--universal_layout_y_axis_insertion_max_tick ' +
                str(arg_universal_layout_max_tick_insertion) +
                ' --universal_layout_y_axis_deletion_max_tick ' +
                str(arg_universal_layout_max_tick_deletion)
              )
            else:
              arg_universal_layout_max_tick = ''
            arg_ext = '--ext ' + output_ext

            file_out.write(f"python {generate_constants.get_python_script('plot_graph')} --input {input_dir} --output {output_dir} {arg_precomputed_layout_dir} {arg_ext} {arg_layout} {arg_reverse_complement} {arg_width_height} {arg_range_x} {arg_range_y} {arg_universal_layout_axis_pos} {arg_universal_layout_max_tick}\n")
        log_utils.log(file_out.name)