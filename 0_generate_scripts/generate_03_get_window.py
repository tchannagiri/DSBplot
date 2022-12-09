
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import generate_constants
import library_constants
import generate_02_combine_repeat

def get_input_file(name):
  return generate_02_combine_repeat.get_output_file(name)

def get_output_dir(name):
  return generate_constants.join_path(
    [
      generate_constants.get_output_dir('window'),
      name,
    ]
  )

if __name__ == '__main__':
  for ext in ['sh', 'ps1']:
    with open(
      file = os.path.join('run_03_get_window' + os.path.extsep + ext),
      mode = 'w',
      encoding = generate_constants.OUTPUT_ENCODING[ext],
    ) as file_out:
      # get windows
      for info in generate_constants.EXPERIMENT_INFO.to_dict('records'):
        if info['version'] != library_constants.VERSION_MERGED:
          for subst_type in library_constants.SUBST_TYPES:
            input_file = get_input_file(info['name'])
            ref_seq_file = generate_constants.join_path(
              [generate_constants.REF_SEQ_DIR, info['ref_seq_file']]
            )
            dsb_pos = int(info['dsb_pos'])
            dsb_type = info['dsb_type']
            strand = info['strand']
            guide_rna = info['guide_rna']
            cell_line = info['cell_line']
            construct = info['construct']
            control_type = info['control_type']
            version = info['version']
            output_dir = get_output_dir(info['name'])
            file_out.write(f"python {generate_constants.get_python_script('get_window')} --input {input_file} --ref_seq_file {ref_seq_file} --output {output_dir} --dsb_pos {dsb_pos} --dsb_type {dsb_type} --strand {strand} --guide_rna {guide_rna} --cell_line {cell_line} --construct {construct} --subst_type {subst_type} --control_type {control_type} --version {version}\n")
      
      # get merged for old/new antisense libraries
      for info in generate_constants.EXPERIMENT_INFO.to_dict('records'):
        if info['version'] == library_constants.VERSION_MERGED:
          for subst_type in library_constants.SUBST_TYPES:
            input_dirs = ' '.join(
              get_output_dir(info['name'].replace(library_constants.VERSION_MERGED, x))
              for x in [library_constants.VERSION_OLD, library_constants.VERSION_NEW]
            )
            version = library_constants.VERSION_MERGED
            output_dir = get_output_dir(info['name'])
            new_library_names = ' '.join(info['library_list'])
            file_out.write(f"python {generate_constants.get_python_script('get_merged')} --input {input_dirs} --output {output_dir} --subst_type {subst_type} --version {version} --new_library_names {new_library_names}\n")

      # get freqs
      for info in generate_constants.EXPERIMENT_INFO.to_dict('records'):
        total_reads = ' '.join([str(x) for x in info['total_reads_list']])
        for subst_type in library_constants.SUBST_TYPES:
          output_dir = get_output_dir(info['name'])
          file_out.write(f"python {generate_constants.get_python_script('get_freq')} --input {output_dir} --output {output_dir} --subst_type {subst_type} --total_reads {total_reads} --freq_min 1e-5\n")

      # make comparison directories
      for info in generate_constants.EXPERIMENT_INFO_COMPARISON.to_dict('records'):
        input_dir_1 = get_output_dir(info['name_1'])
        input_dir_2 = get_output_dir(info['name_2'])
        output_dir = get_output_dir(info['name'])
        for subst_type in library_constants.SUBST_TYPES:
          file_out.write(f"python {generate_constants.get_python_script('get_freq_comparison')} --input {input_dir_1} {input_dir_2} --output {output_dir} --subst_type {subst_type}\n")

      log_utils.log(file_out.name)