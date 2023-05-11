import argparse
import pandas as pd

import DSBplot.utils.file_utils as file_utils
import DSBplot.utils.common_utils as common_utils
import DSBplot.utils.log_utils as log_utils

def parse_args():
  parser = argparse.ArgumentParser(
    description = (
      'Combine biological repeats of the same experiment.' +
      ' Keeps sequences present in any of the 4 repeats.'
    )
  )
  parser.add_argument(
    '--input',
    type = common_utils.check_file,
    help = 'TSV files output from script "filter_nhej.py".',
    nargs = '+',
    required = True,
  )
  parser.add_argument(
    '--column_names',
    nargs = '+',
    required = True,
    help = (
      'Names to use as suffixes to the "Count" columns of the output.' +
      ' Number of arguments must match the number of INPUT args.'
    ),
  )
  parser.add_argument(
    '--output',
    type = common_utils.check_file_output,
    help = 'Output TSV file name.',
    required = True,
  )
  parser.add_argument(
    '--quiet',
    help = 'Do not output log messages.',
    action = 'store_true',
  )
  args = vars(parser.parse_args())
  if len(args['input']) != len(args['column_names']):
    raise Exception(
      f"Number of columns names {len(args['column_names'])} does" +
      f" not match number of inputs {len(args['input'])}."
    )
  return args

def main(input, column_names, output, quiet = True):
  num_repeats = len(input)

  for x in input:
    log_utils.log_input(x)

  data = [
    pd.read_csv(input[i], sep='\t').set_index(['Sequence', 'CIGAR'])
    for i in range(num_repeats)
  ]

  if not quiet:
    for i in range(num_repeats):
      log_utils.log(f"Num sequences {column_names[i]}: {data[i].shape[0]}")

  data = pd.concat(data, axis='columns', join='outer', keys=column_names)
  data.columns = data.columns.map(lambda x: '_'.join([x[1], x[0]]))
  data = data.reset_index()
  
  data_combined = pd.DataFrame(
    {
      'Sequence': list(data['Sequence']),
      'Num_Subst': list(data['Num_Subst_' + column_names[0]]),
      'CIGAR': list(data['CIGAR'])
    },
    index = data.index,
  )
  data_combined['Num_Subst'] = data_combined['Num_Subst'].fillna(0).astype(int)
  for i in range(num_repeats):
    data_combined['Count_' + column_names[i]] = data['Count_' + column_names[i]].fillna(0).astype(int)

  if not quiet:
    log_utils.log(f"Num sequences combined: {data_combined.shape[0]}")
  file_utils.write_tsv(data_combined, output)
  log_utils.log_output(output)
  log_utils.new_line()
  
if __name__ == '__main__':
  main(**parse_args())
