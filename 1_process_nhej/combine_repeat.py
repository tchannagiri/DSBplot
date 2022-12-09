import argparse
import sys
import os
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

import file_utils
import common_utils
import log_utils

def main():
  parser = argparse.ArgumentParser(
    description = (
      'Combine biological repeats of the same experiment.' +
      ' Keeps sequences present in any of the 4 repeats.'
    )
  )
  parser.add_argument(
    '--input',
    type = argparse.FileType('r'),
    help = (
      'TSV files output from script "filter_nhej.py".' +
      ' Must have columns: Sequence, CIGAR, Count, Num_Subst.' +
      ' The files names must be of the form "<name>_XXX",' +
      ' where <name> is the name of the library and will be' +
      ' used to name the columns of the combined output.'
    ),
    nargs = '+',
    required = True,
  )
  parser.add_argument(
    '--output',
    type = common_utils.check_file_output,
    help = 'Output file name',
    required = True,
  )
  parser.add_argument(
    '--quiet',
    help = 'Do not output log messages.',
    action = 'store_true',
  )
  args = parser.parse_args()

  num_repeats = len(args.input)

  for input_file in args.input:
    log_utils.log(input_file.name)

  data = [
    pd.read_csv(args.input[i], sep='\t').set_index(['Sequence', 'CIGAR'])
    for i in range(num_repeats)
  ]
  names = [
    os.path.splitext(os.path.basename(args.input[i].name))[0].split('_')[0]
    for i in range(num_repeats)
  ]

  if not args.quiet:
    for i in range(num_repeats):
      log_utils.log(f"Num sequences {names[i]}: {data[i].shape[0]}")

  data = pd.concat(data, axis='columns', join='outer', keys=names)
  data.columns = data.columns.map(lambda x: '_'.join([x[1], x[0]]))
  data = data.reset_index()
  
  data_combined = pd.DataFrame(
    {
      'Sequence': list(data['Sequence']),
      'Num_Subst': list(data['Num_Subst_' + names[0]]),
      'CIGAR': list(data['CIGAR'])
    },
    index = data.index,
  )
  data_combined['Num_Subst'] = data_combined['Num_Subst'].fillna(0).astype(int)
  for i in range(num_repeats):
    data_combined['Count_' + names[i]] = data['Count_' + names[i]].fillna(0).astype(int)

  if not args.quiet:
    log_utils.log(f"Num sequences combined: {data_combined.shape[0]}\n")
  file_utils.write_tsv(data_combined, args.output)
  log_utils.log('------>')
  log_utils.log(args.output.name)
  log_utils.new_line()
  
if __name__ == '__main__':
  main()
