import pandas as pd
import argparse

import DSBplot.utils.common_utils as common_utils
import DSBplot.utils.file_utils as file_utils
import DSBplot.utils.log_utils as log_utils
import DSBplot.utils.alignment_utils as alignment_utils
import DSBplot.utils.constants as constants
import DSBplot.preprocessing.alignment_window as alignment_window
import DSBplot.preprocessing.remove_substitution as remove_substitution

def parse_args():
  parser = argparse.ArgumentParser(
    description = 'Get windows from aligned/filtered reads.',
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
  )
  parser.add_argument(
    '--input',
    type = common_utils.check_file,
    help = 'CSV files output from script "filter_nhej.py".',
    nargs = '+',
    required = True,
  )
  parser.add_argument(
    '--ref',
    type = common_utils.check_file,
    help = (
      'Reference sequence FASTA.' +
      ' Should contain a single nucleotide sequence in FASTA format.'
    ),
    required = True,
  )
  parser.add_argument(
    '--output',
    nargs = 2,
    type = common_utils.check_file_output,
    help = 'Output CSV file for window data and metadata, respectively.',
    required = True,
  )
  parser.add_argument(
    '--dsb',
    type = int,
    help = (
      'Position on reference sequence immediately upstream of DSB site.' +
      ' I.e., the DSB is between 1-based positions DSB and DSB + 1.'
    ),
    required = True,
  )
  parser.add_argument(
    '--window',
    type = int,
    default = 10,
    help = (
      'Size of window around DSB site to extract.' +
      ' The nucleotides at the positions' +
      ' {DSB - WINDOW + 1, ..., DSB + WINDOW} are extracted.' +
      ' The actual number of nucleotides extracted may vary depending' +
      ' on how many insertions/deletions in the alignment of each sequence.'
    ),
  )
  parser.add_argument(
    '--anchor',
    type = int,
    default = 20,
    help = 'Size of anchor on left/right of the window to check for substitutions.',
  )
  parser.add_argument(
    '--anchor_sub',
    type = int,
    default = 1,
    help = (
      'Maximum number of substitutions allowed on the left/right anchor sequences.' +
      ' Reads with more than the allowed number of substitutions on the left/right anchor' +
      ' will be discarded. This limit is applied to the left/right anchors separately.'
    ),
  )
  parser.add_argument(
    '--sub',
    type = int,
    choices = [0, 1],
    help = (
      'Whether to ignore (0) or keep (1) alignment substitutions.' +
      ' If ignoring alignment substitutions, the corresponding nucleotide on the' +
      ' read is replaced with the reference sequence nucleotide.'
    ),
    required = True,
  )
  parser.add_argument(
    '--name',
    type = str,
    help = (
      'Name identifying the experiment' +
      ' Should be a legal part of a file name without path separators.'
    ),
    required = True,
  )
  parser.add_argument(
    '--label',
    type = str,
    help = 'Label of the experiment.',
    required = True,
  )
  args = vars(parser.parse_args())
  args['sub'] = constants.SUBST_TYPES[args['sub']]
  return args

def write_window(
  input,
  output,
  ref_seq,
  dsb_pos,
  window_size,
  anchor_size,
  anchor_sub,
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
      anchor_sub = anchor_sub,
    )
    if ref_align is None:
      continue
  
    # Optionally remove substitutions
    if subst_type == 'withoutSubst':
      read_align = remove_substitution.remove_substitutions(
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
    anchor_sub = 0,
  )
  return ref_seq_window

def main(
  input,
  output,
  ref,
  dsb,
  window,
  anchor,
  anchor_sub,
  sub,
  name,
  label,
):
  log_utils.log_input(input)

  ref_seq = file_utils.read_seq(ref)
  write_window(
    input = input,
    output = output[0],
    ref_seq = ref_seq,
    dsb_pos = dsb,
    window_size = window,
    anchor_size = anchor,
    anchor_sub = anchor_sub,
    subst_type = sub,
  )

  data_info = common_utils.make_data_info(
    format = 'individual',
    names = [name],
    labels = [label],
    ref_seqs = [ref_seq],
    ref_seq_window = get_ref_seq_window(
      ref_seq = ref_seq,
      dsb_pos = dsb,
      window_size = window,
    ),
  )
  data_info = pd.DataFrame(data_info, index=[0])
  log_utils.log_output(output[1])
  file_utils.write_csv(data_info, output[1])

  log_utils.blank_line()

if __name__ == '__main__':
  main(**parse_args())
