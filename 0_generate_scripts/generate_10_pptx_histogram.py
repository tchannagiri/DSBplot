import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import library_constants
import generate_constants
import generate_08_plot_histogram

def get_output_file(cell_line, intron_type, version):
  version_str = '' if (version == library_constants.VERSION_NONE) else ('_' + version)
  return generate_constants.join_path(
    [
      generate_constants.get_output_dir('pptx'),
      'histogram',
      cell_line + '_' + intron_type + version_str + os.extsep + 'pptx',
    ]
  )

def get_input_file(experiment_name, variation_type, file_ext):
  return generate_constants.join_path(
    [
      generate_08_plot_histogram.get_output_dir(
        library_constants.SUBST_WITH
      ),
      experiment_name + '_' + variation_type + os.path.extsep + file_ext,
    ]
  )

VARIATION_TYPES = [
  library_constants.VARIATION_INSERTION,
  library_constants.VARIATION_DELETION,
  library_constants.VARIATION_SUBSTITUTION,
]

if __name__ == '__main__':
  for ext in ['sh', 'ps1']:
    with open(
      file = 'run_10_pptx_histogram' + os.extsep + ext,
      mode = 'w',
      encoding = generate_constants.OUTPUT_ENCODING[ext],
    ) as file_out:
      for cell_line in library_constants.CELL_LINES:
        if cell_line == library_constants.CELL_LINE_WT:
          intron_type_list = ['sense', 'antisense']
        elif cell_line == library_constants.CELL_LINE_KO:
          intron_type_list = ['sense']
        else:
          raise Exception('Unknown cell line: ' + str(cell_line))
        for intron_type in intron_type_list:
          if intron_type == 'sense':
            construct_list = library_constants.CONSTRUCTS_INDIVIDUAL_SENSE
            version_list = [library_constants.VERSION_NONE]
            row_spec_list = [
              {'strand': library_constants.STRAND_R1, 'control_type': library_constants.CONTROL_NOT, 'guide_rna': library_constants.GUIDE_RNA_A},
              {'strand': library_constants.STRAND_R2, 'control_type': library_constants.CONTROL_NOT, 'guide_rna': library_constants.GUIDE_RNA_B},
              {'strand': library_constants.STRAND_R1, 'control_type': library_constants.CONTROL_NOT, 'guide_rna': library_constants.GUIDE_RNA_AB},
              {'strand': library_constants.STRAND_R2, 'control_type': library_constants.CONTROL_NOT, 'guide_rna': library_constants.GUIDE_RNA_AB},
              {'strand': library_constants.STRAND_R1, 'control_type': library_constants.CONTROL_NODSB, 'guide_rna': library_constants.GUIDE_RNA_A},
              {'strand': library_constants.STRAND_R2, 'control_type': library_constants.CONTROL_NODSB, 'guide_rna': library_constants.GUIDE_RNA_B},
              {'strand': library_constants.STRAND_R1, 'control_type': library_constants.CONTROL_30BPDOWN, 'guide_rna': library_constants.GUIDE_RNA_A},
              {'strand': library_constants.STRAND_R2, 'control_type': library_constants.CONTROL_30BPDOWN, 'guide_rna': library_constants.GUIDE_RNA_B},
            ]
          elif intron_type == 'antisense':
            construct_list = library_constants.CONSTRUCTS_INDIVIDUAL_ANTISENSE
            version_list = [library_constants.VERSION_OLD, library_constants.VERSION_NEW, library_constants.VERSION_MERGED]
            row_spec_list = [
              {'strand': library_constants.STRAND_R1, 'control_type': library_constants.CONTROL_NOT, 'guide_rna': library_constants.GUIDE_RNA_CD},
              {'strand': library_constants.STRAND_R2, 'control_type': library_constants.CONTROL_NOT, 'guide_rna': library_constants.GUIDE_RNA_CD},
            ]
          else:
            raise Exception('Impossible.')
          num_rows = len(row_spec_list)
          num_cols = len(construct_list) * len(VARIATION_TYPES)

          # Top labels
          top_labels = []
          for construct in construct_list:
            for variation in VARIATION_TYPES:
              top_labels.append(
                library_constants.LABELS[construct] +
                generate_constants.ARG_NEWLINE[ext] +
                library_constants.LABELS[variation]
              )

          # Left labels
          left_labels = []
          for row_spec in row_spec_list:
            labels = [library_constants.LABELS[row_spec['guide_rna']]]
            if row_spec['guide_rna'] in [library_constants.GUIDE_RNA_AB, library_constants.GUIDE_RNA_CD]:
              labels.append(library_constants.LABELS[row_spec['strand']])
            if row_spec['control_type'] != library_constants.CONTROL_NOT:
              labels.append(library_constants.LABELS[row_spec['control_type']])
            left_labels.append(generate_constants.ARG_NEWLINE[ext].join(labels))

          for version in version_list:
            file_list = []
            for row_spec in row_spec_list:
              for construct in construct_list:
                info = generate_constants.EXPERIMENT_INFO.loc[
                  (generate_constants.EXPERIMENT_INFO['cell_line'] == cell_line) &
                  (generate_constants.EXPERIMENT_INFO['version'] == version) &
                  (generate_constants.EXPERIMENT_INFO['construct'] == construct) &
                  (generate_constants.EXPERIMENT_INFO['strand'] == row_spec['strand']) &
                  (generate_constants.EXPERIMENT_INFO['control_type'] == row_spec['control_type']) &
                  (generate_constants.EXPERIMENT_INFO['guide_rna'] == row_spec['guide_rna'])
                ]
                if info.shape[0] != 1:
                  raise Exception(f'Got {info.shape[0]} rows. Expected 1.')
                info = info.iloc[0].to_dict()
                for variation_type in [
                  library_constants.VARIATION_INSERTION,
                  library_constants.VARIATION_DELETION,
                  library_constants.VARIATION_SUBSTITUTION,
                ]:
                  file_list.append(get_input_file(info['name'], variation_type, file_ext='png'))
            arg_input = '--input ' + ' '.join(file_list)
            arg_top_margin_labels = '--top_margin_labels ' + ' '.join(f'"{x}"' for x in top_labels)
            arg_left_margin_labels = '--left_margin_labels ' + ' '.join(f'"{x}"' for x in left_labels)
            arg_num_grids = '--num_grids 1'
            arg_num_rows = '--num_rows ' + str(num_rows)
            arg_num_cols = '--num_cols ' + str(num_cols)
            arg_output = '--output ' + get_output_file(cell_line, intron_type, version)
            file_out.write(f"python {generate_constants.get_python_script('get_pptx')} {arg_input} {arg_output} {arg_top_margin_labels} {arg_left_margin_labels} {arg_num_grids} {arg_num_rows} {arg_num_cols}\n")
      log_utils.log(file_out.name)

