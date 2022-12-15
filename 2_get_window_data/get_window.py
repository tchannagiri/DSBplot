import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import shutil

import common_utils
import file_utils
import log_utils
import alignment_utils
import fasta_utils
import alignment_window
import remove_substitution
import library_constants
import file_names

import pandas as pd
import argparse

def parse_args():
  parser = argparse.ArgumentParser(
    description = 'Process data for downstream graph and histogram analysis.'
  )
  parser.add_argument(
    '--input',
    type = common_utils.check_file,
    help = (
      'Table of sequences produced with combine_repeat.py.' +
      ' Column format: Sequence, CIGAR, Count_<X1>, Count_<X2>, ..., etc.' +
      ' All the columns after CIGAR should be the counts for each repeat where' +
      ' <Xi> denotes the name of the library.'
    ),
    required = True,
  )
  parser.add_argument(
    '--ref_seq_file',
    type = common_utils.check_file,
    help = 'Reference sequence FASTA. Should contain a single nucleotide sequence in FASTA format.',
    required = True,
  )
  parser.add_argument(
    '--output',
    type = common_utils.check_dir_output,
    help = 'Output directory.',
    required = True,
  )
  parser.add_argument(
    '--dsb_pos',
    type = int,
    help = (
      'Position on reference sequence immediately upstream of DSB site.\n'
      'The DSB is between position DSB_POS and DSB_POS + 1.'
    ),
    required = True,
  )
  parser.add_argument(
    '--window_size',
    type = int,
    default = 10,
    help = (
      'Size of window around DSB site to extract.' +
      ' The nucleotides at the positions' +
      ' {DSB_POS - WINDOW_SIZE + 1, ..., DSB_POS + WINDOW_SIZE} are extracted.' +
      ' The actual number of nucleotides extracted may vary depending' +
      ' on how many insertions/deletion the alignment of the sequence has.'
    ),
  )
  parser.add_argument(
    '--anchor_size',
    type = int,
    default = 20,
    help = (
      'Size of anchor on left/right of the window to check for mismatches.'
    ),
  )
  parser.add_argument(
    '--anchor_mismatches',
    type = int,
    default = 1,
    help = (
      'Maximum number of mismatches allowed on the left/right anchor sequences.\n'
      'Reads with more than the allowed number of mismatches on the left/right anchor\n'
      'will be discarded. This limit is applied to the left/right anchors separately.'
    ),
  )
  parser.add_argument(
    '--subst_type',
    type = str,
    default = library_constants.SUBST_WITHOUT,
    choices = [
      library_constants.SUBST_WITH,
      library_constants.SUBST_WITHOUT
    ],
    help = 'Whether to keep or ignore substitutions.',
    required = True,
  )
  parser.add_argument(
    '--label',
    type = str,
    help = 'Label of the library.',
    required = True,
  )
  return vars(parser.parse_args())

def write_alignment_window(
  input_file,
  output_dir,
  ref_seq,
  dsb_pos,
  window_size,
  anchor_size,
  anchor_mismatches,
  subst_type,
):
  data = file_utils.read_tsv(input_file)
  count_columns_old = [x for x in data.columns if x.startswith('Count_')]
  data = data[['Sequence', 'CIGAR'] + count_columns_old]

  data_new = {
    'ref_align': [],
    'read_align': [],
  }
  count_columns_new = [x.lower() for x in count_columns_old]
  for column in count_columns_new:
    data_new[column] = []
    
  for row in data.to_dict('records'):
    # convert CIGAR to alignment
    ref_align, read_align = alignment_utils.get_alignment(ref_seq, row['Sequence'], 1, row['CIGAR'])

    # extract window around DSB
    ref_align, read_align = alignment_window.get_alignment_window(
      ref_align,
      read_align,
      dsb_pos,
      window_size,
      anchor_size,
      anchor_mismatches,
    )
    if ref_align is None:
      continue
  
    # optionally remove substitutions
    if subst_type == 'withoutSubst':
      read_align = remove_substitution.remove_substitutions(
        ref_align,
        read_align,
      )
    
    data_new['ref_align'].append(ref_align)
    data_new['read_align'].append(read_align)

    for col_old, col_new in zip(count_columns_old, count_columns_new):
      data_new[col_new].append(row[col_old])

  data = pd.DataFrame(data_new)

  # Sum rows with identical sequences
  # Note: we make an arbitrary choice of which alignment is kept
  data['read_seq'] = data['read_align'].apply(alignment_utils.get_orig_seq)
  data = data.groupby('read_seq').aggregate(
    ref_align = ('ref_align', 'first'),
    read_align = ('read_align', 'first'),
    **{col: (col, 'sum') for col in count_columns_new}
  ).reset_index()
  data = data.drop('read_seq', axis='columns').reset_index(drop=True)

  # get the min frequency of the repeats
  data['count_min'] = data[count_columns_new].min(axis='columns')
  data = data.sort_values(
    ['count_min', 'read_align'],
    ascending = [False, True],
  ).reset_index(drop=True)

  # save the unfiltered repeat data
  data = data.drop('count_min', axis='columns')
  output_file = file_names.window(
    output_dir,
    library_constants.COUNT,
    subst_type,
  )
  file_utils.write_tsv(data, output_file)
  log_utils.log_output(output_file)

def write_data_info(
  dir,
  format,
  labels,
  ref_seqs,
  ref_seq_window,
):
  data_info = {
    'format': format,
    'ref_seq_window': ref_seq_window,
  }
  if format == library_constants.DATA_INDIVIDUAL:
    if len(labels) != 1:
      raise Exception(f'Expected 1 name for individual format. Got: {len(labels)}')
    if len(ref_seqs) != 1:
      raise Exception(f'Expected 1 reference sequence for individual format. Got: {len(ref_seqs)}')
    data_info['label'] = labels[0]
    data_info['ref_seq'] = ref_seqs[0]
  elif format == library_constants.DATA_COMPARISON:
    if len(labels) != 2:
      raise Exception(f'Expected 2 names for comparison format. Got: {len(labels)}')
    if len(ref_seqs) != 2:
      raise Exception(f'Expected 2 reference sequences for comparison format. Got: {len(ref_seqs)}')
    data_info['label_1'] = labels[0]
    data_info['label_2'] = labels[1]
    data_info['ref_seq_1'] = ref_seqs[0]
    data_info['ref_seq_2'] = ref_seqs[1]
  else:
    raise Exception('Unknown data format: ' + str(format))
  data_info = pd.DataFrame(data_info, index = [0])
  file_out = file_names.data_info(dir)
  log_utils.log_output(file_out)
  file_utils.write_tsv(data_info, file_out)

def get_ref_seq_window(ref_seq, dsb_pos, window_size):
  """
    Get the window of the reference sequence around the DSB.
  """
  # Extract the window of a perfect alignment.
  # Do it this way so that it is always consistent
  # with how alignment windows are extracted.
  ref_seq_window, _ = alignment_window.get_alignment_window(
    ref_seq,
    ref_seq,
    dsb_pos,
    window_size,
    0,
    0,
  )
  return ref_seq_window

def main(
  input,
  output,
  ref_seq_file,
  dsb_pos,
  window_size,
  anchor_size,
  anchor_mismatches,
  subst_type,
  label,
):
  log_utils.log_input(input)

  ref_seq = fasta_utils.read_fasta_seq(ref_seq_file)
  write_alignment_window(
    input_file = input, 
    output_dir = output,
    ref_seq = ref_seq,
    dsb_pos = dsb_pos,
    window_size = window_size,
    anchor_size = anchor_size,
    anchor_mismatches = anchor_mismatches,
    subst_type = subst_type,
  )
  write_data_info(
    dir = output,
    format = library_constants.DATA_INDIVIDUAL,
    labels = [label],
    ref_seqs = [ref_seq],
    ref_seq_window = get_ref_seq_window(
      ref_seq,
      dsb_pos,
      window_size,
    ),
  )
  log_utils.new_line()

if __name__ == '__main__':
  main(**parse_args())