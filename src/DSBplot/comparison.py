import argparse

import DSBplot.utils.common_utils as common_utils
import DSBplot.preprocess as preprocess

def parse_args():
  parser = argparse.ArgumentParser(
    description = (
      'Combine two individual experiment directories to make a comparison' +
      ' experiment directory for comparison graphs.' +
      ' The experiments must be have the same windowed reference sequence.'
    ),
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
  )
  parser.add_argument(
    '--input',
    nargs = 2,
    type = common_utils.check_dir,
    help = 'Input directories of data created with "preprocess.py".',
    required = True,
  )
  parser.add_argument(
    '--output',
    type = common_utils.check_dir_output,
    help = 'Output directory.',
    required = True,
  )
  parser.add_argument(
    '--stages',
    type = str,
    choices = [x for x in preprocess.STAGES_2 if x != '5_histogram'],
    default = [x for x in preprocess.STAGES_2 if x != '5_histogram'],
    nargs = '+',
    help = 'Stages to run. See the documentation for "preprocess.py".',
  )
  return vars(parser.parse_args())

def main(
  input,
  output,
  stages,
):
  preprocess.do_stages_2(
    output = output,
    input_comparison = input,
    stages = stages,
  )

if __name__ == '__main__':
  main(**parse_args())