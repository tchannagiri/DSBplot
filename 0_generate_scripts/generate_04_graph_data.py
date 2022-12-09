
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import generate_constants
import library_constants
import generate_03_get_window

def get_input_dir(name):
  return generate_03_get_window.get_output_dir(name)

def get_output_dir(name):
  return generate_constants.join_path(
    [
      generate_constants.get_output_dir('graph'),
      name,
    ]
  )

if __name__ == '__main__':
  for ext in ['sh', 'ps1']:
    with open(
      file = os.path.join('run_04_graph_data' + os.path.extsep + ext),
      mode = 'w',
      encoding = generate_constants.OUTPUT_ENCODING[ext],
    ) as file_out:
      for subst_type in library_constants.SUBST_TYPES:
        for info in (
          generate_constants.EXPERIMENT_INFO.to_dict('records') +
          generate_constants.EXPERIMENT_INFO_COMPARISON.to_dict('records')
        ):
          input_dir = get_input_dir(info['name'])
          output_dir = get_output_dir(info['name'])
          file_out.write(f"python {generate_constants.get_python_script('get_graph_data')} --input {input_dir} --output {output_dir} --subst_type {subst_type}\n")
      log_utils.log(file_out.name)