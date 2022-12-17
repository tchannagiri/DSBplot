import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'utils'))) # allow importing

import argparse

import common_utils

import preprocess

def parse_args():
  parser = argparse.ArgumentParser(
    description = (
      'Combine two individual experiment directories to make a comparison' +
      ' experiment directory for comparison graphs.' +
      ' The experiments must be have the same windowed reference sequence.'
    )
  )
  parser.add_argument(
    '--input',
    nargs = 2,
    type = common_utils.check_dir,
    help = 'Input directories of data created with "preprocess.py".',
    required = True
  )
  parser.add_argument(
    '--output',
    type = common_utils.check_dir_output,
    help = 'Output directory.',
    required = True,
  )
  return vars(parser.parse_args())

def main(
  input,
  output,
):
  preprocess.do_stage_2(
    output = output,
    input_comparison = input,
  )

if __name__ == '__main__':
  main(**parse_args())