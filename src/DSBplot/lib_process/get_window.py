import pandas as pd
import argparse

import DSBplot.utils.common_utils as common_utils
import DSBplot.utils.file_utils as file_utils
import DSBplot.utils.log_utils as log_utils
import DSBplot.utils.alignment_utils as alignment_utils
import DSBplot.utils.constants as constants
import DSBplot.lib_process.filter as filter
import DSBplot.lib_process.alignment_window as alignment_window

PARAMS = {
  '-i': {
    'type': common_utils.check_file,
    'help': 'CSV files output from script "filter.py".',
    'nargs': '+',
    'required': True,
    'metavar': 'INPUT',
    'dest': 'input',
  },  
  '--ref': filter.PARAMS['--ref'],
  '-o': {
    'type': common_utils.check_file_output,
    'help': 'Output CSV file for window data.',
    'required': True,
    'metavar': 'OUTPUT',
    'dest': 'output',
  },
  '--dsb': filter.PARAMS['--dsb'],
  '--window': {
    'type': int,
    'default': 10,
    'help': (
      'Size of window around DSB site to obtain.' +
      ' The nucleotides at the positions' +
      '[DSB - WINDOW + 1, ..., DSB + WINDOW] are obtained.' +
      ' The actual number of nucleotides obtained may vary depending' +
      ' on the number of indels in the alignment.'
    ),
    'dest': 'window_size',
  },
  '--anchor': {
    'type': int,
    'default': 20,
    'help': (
      'Size of anchor on left/right of the window to check for' +
      ' substitutions and indels. See "--anchor_var".'
    ),
    'dest': 'anchor_size',
  },
  '--anchor_vars': {
    'type': int,
    'nargs': 2,
    'default': [1, 0],
    'help': (
      'Maximum number of substitutions (arg 1) and indels (arg 2) allowed on the' +
      ' left/right anchor sequences. Reads with more than the allowed number of' +
      ' substitutions or indels on the left/right anchor will be discarded.' +
      ' This limit is applied to the left/right anchors separately.' +
      ' Set to either/both to -1 to disable the respective limit.'
    ),
    'dest': 'anchor_vars',
  },
  '--sub': {
    'type': int,
    'choices': [0, 1],
    'help': (
      'Whether to ignore (0) or keep (1) alignment substitutions.' +
      ' If ignoring alignment substitutions, the corresponding nucleotide on the' +
      ' read is replaced with the reference sequence nucleotide.'
    ),
    'required': True,
    'dest': 'subst_type',
  },
}

def post_process_args(args):
  args = args.copy()
  if args.get('subst_type') is not None:
    args['subst_type'] = constants.SUBST_TYPES[args['subst_type']]
  if args.get('anchor_vars') is not None:
    args['anchor_substs'] = args['anchor_vars'][0]
    args['anchor_indels'] = args['anchor_vars'][1]
    del args['anchor_vars']
  return args

def parse_args():
  parser = argparse.ArgumentParser(
    description = 'Get windows around DSB from aligned/filtered reads.',
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
  )
  for param, options in PARAMS.items():
    parser.add_argument(param, **options)
  return post_process_args(vars(parser.parse_args()))

def remove_substitutions(ref_align, read_align):
  """
    Remove the substitutions from read_align and replace them with matches.

    Parameters
    ----------
    ref_align : the reference sequence alignment string
    read_align : the read sequence alignment string

    The above strings should be in "alignment matrix" format.
    That is, string over the alphabet ["A", "C", "G", "T", "-"].

    Returns
    -------
    tuple (ref_align, read_align_new)
      ref_align : identical to the parameter
      read_align_new : read_align with the substutitions converted to matches
  """
  read_align_new = ''
  for i in range(len(ref_align)):
    if (
      (ref_align[i] != read_align[i]) and
      (ref_align[i] != '-') and
      (read_align[i] != '-')
    ):
      read_align_new += ref_align[i]
    else:
      read_align_new += read_align[i]
  return read_align_new

def write_window(
  input,
  output,
  ref_seq,
  dsb_pos,
  window_size,
  anchor_size,
  anchor_substs,
  anchor_indels,
  subst_type,
):
  data = file_utils.read_csv(input)

  data_new = {
    'ref_align': [],
    'read_align': [],
  }
  value_cols = [x for x in data.columns if (x.startswith('count_') or x.startswith('freq_'))]
  for col in value_cols:
    data_new[col] = []
  for row in data.to_dict('records'):
    # Convert CIGAR to alignment
    ref_align, read_align = alignment_utils.get_alignment(ref_seq, row['seq'], 1, row['cigar'])

    # Get window around DSB
    ref_align, read_align = alignment_window.get_alignment_window(
      ref_align = ref_align,
      read_align = read_align,
      dsb_pos = dsb_pos,
      window_size = window_size,
      anchor_size = anchor_size,
      anchor_substs = anchor_substs,
      anchor_indels = anchor_indels,
    )
    if ref_align is None:
      continue
  
    # Optionally remove substitutions
    if subst_type == 'withoutSubst':
      read_align = remove_substitutions(
        ref_align,
        read_align,
      )
    
    data_new['ref_align'].append(ref_align)
    data_new['read_align'].append(read_align)

    for col in value_cols:
      data_new[col].append(row[col])

  data = pd.DataFrame(data_new)

  # Sum rows with identical sequences
  # Note: we make an arbitrary choice of which alignment is kept
  data['read_seq'] = data['read_align'].apply(alignment_utils.get_orig_seq)
  data = data.groupby('read_seq').aggregate(
    ref_align = ('ref_align', 'first'),
    read_align = ('read_align', 'first'),
    **{col: (col, 'sum') for col in value_cols}
  ).reset_index()
  data = data.drop(columns='read_seq').reset_index(drop=True)

  data = data.sort_values('freq_mean', ascending=False).reset_index(drop=True)

  file_utils.write_csv(data, output)
  log_utils.log_output(output)

def get_ref_seq_window(ref_seq, dsb_pos, window_size):
  """
    Get the window of the reference sequence around the DSB.
  """
  # Extract the window of a perfect alignment.
  # Do it this way so that it is always consistent
  # with how alignment windows are extracted.
  ref_seq_window, _ = alignment_window.get_alignment_window(
    ref_align = ref_seq,
    read_align = ref_seq,
    dsb_pos = dsb_pos,
    window_size = window_size,
    anchor_size = 0,
    anchor_substs = 0,
    anchor_indels = 0,
  )
  return ref_seq_window

def main(
  input,
  output,
  ref_seq_file,
  dsb_pos,
  window_size,
  anchor_size,
  anchor_substs,
  anchor_indels,
  subst_type,
):
  log_utils.log_input(input)

  write_window(
    input = input,
    output = output,
    ref_seq = file_utils.read_seq(ref_seq_file),
    dsb_pos = dsb_pos,
    window_size = window_size,
    anchor_size = anchor_size,
    anchor_substs = anchor_substs,
    anchor_indels = anchor_indels,
    subst_type = subst_type,
  )

if __name__ == '__main__':
  main(**parse_args())
