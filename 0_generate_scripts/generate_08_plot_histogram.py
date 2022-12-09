
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import generate_constants
import library_constants
import generate_05_histogram_data

def get_input_dir(name):
  return generate_05_histogram_data.get_output_dir(name)

def get_output_dir(subst_type):
  return generate_constants.join_path(
    [
      generate_constants.get_output_dir('plot_histogram'),
      subst_type,
    ]
  )

if __name__ == '__main__':
  for script_ext in ['sh', 'ps1']:
    with open(
      file = os.path.join('run_08_plot_histogram' + os.path.extsep + script_ext),
      mode = 'w',
      encoding = generate_constants.OUTPUT_ENCODING[script_ext],
    ) as file_out:
      for info in generate_constants.EXPERIMENT_INFO.to_dict('records'):
        output_dir = get_output_dir(library_constants.SUBST_WITH)
        input_dir = get_input_dir(info['name'])

        arg_reverse_pos = (
          ''
          if (info['strand'] == library_constants.STRAND_R2) else
          '--reverse_pos'
        )
        arg_label_type = '--label_type ' + (
          'absolute'
          if (info['control_type'] == library_constants.CONTROL_30BPDOWN) else
          'relative'
        )
        file_out.write(f"python {generate_constants.get_python_script('plot_histogram')} --input {input_dir} --output {output_dir} {arg_reverse_pos} {arg_label_type}\n")
      log_utils.log(file_out.name)