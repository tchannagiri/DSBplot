import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import library_constants
import generate_constants

def get_input_file(library, strand):
  return generate_constants.join_path(
    [
      generate_constants.get_output_dir('fastq'),
      library + '_' + strand + os.extsep + 'fastq'
    ]
  )

def get_output_file(name):
  return generate_constants.join_path(
    [
      generate_constants.get_output_dir('sam'),
      name + os.extsep + 'sam'
    ]
  )

def get_bowtie2_build_file(ref_seq_file):
  return (
    generate_constants.get_output_dir('bowtie2_build') +
    '/' +
    os.path.splitext(os.path.basename(ref_seq_file))[0]
  )

if __name__ == '__main__':
  for ext in ['sh', 'ps1']:
    with open(
      file = os.path.join('run_00_bowtie2_align' + os.path.extsep + ext),
      mode = 'w',
      encoding = generate_constants.OUTPUT_ENCODING[ext],
    ) as file_out:
      for ref_seq_file in generate_constants.LIBRARY_INFO['ref_seq_file'].unique():
        if ref_seq_file is not None:
          ref_seq_bowtie2_build = get_bowtie2_build_file(ref_seq_file)
          file_out.write(f"{generate_constants.BOWTIE2_BUILD_COMMAND[ext]} ref_seq/{ref_seq_file} {ref_seq_bowtie2_build}\n")
      for info in generate_constants.LIBRARY_INFO.to_dict('records'):
        if info['version'] != library_constants.VERSION_MERGED:
          input_file = get_input_file(info['library'], info['strand'])
          output_file = get_output_file(info['name'])
          ref_seq_bowtie2_build = get_bowtie2_build_file(info['ref_seq_file'])
          file_out.write(f"{generate_constants.BOWTIE2_ALIGN_COMMAND[ext]} -x  {ref_seq_bowtie2_build} {input_file} -S {output_file}\n")
      log_utils.log(file_out.name)