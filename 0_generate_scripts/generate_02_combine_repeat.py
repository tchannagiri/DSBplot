import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import library_constants
import generate_constants
import generate_01_filter_nhej

def get_input_files(name_list):
  return ' '.join(
    generate_01_filter_nhej.get_output_file(x)
    for x in name_list
  )

def get_output_file(name):
  return generate_constants.join_path(
    [
      generate_constants.get_output_dir('combine_repeat'),
      name + os.extsep + 'tsv',
    ]
  )

if __name__ == '__main__':
  for ext in ['sh', 'ps1']:
    with open(
      file = os.path.join('run_02_combine_repeat' + os.path.extsep + ext),
      mode = 'w',
      encoding = generate_constants.OUTPUT_ENCODING[ext],
    ) as file_out:
      for info in generate_constants.EXPERIMENT_INFO.to_dict('records'):
        if info['version'] != library_constants.VERSION_MERGED:
          input_files = get_input_files(info['library_name_list'])
          output_file = get_output_file(info['name'])
          file_out.write(f"python {generate_constants.get_python_script('combine_repeat')} --input {input_files} --output {output_file} --quiet\n")
      log_utils.log(file_out.name)
