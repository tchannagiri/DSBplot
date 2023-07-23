import os
import argparse
from collections import defaultdict
import pandas as pd

import DSBplot.utils.constants as constants
import DSBplot.utils.file_names as file_names
import DSBplot.utils.file_utils as file_utils
import DSBplot.utils.common_utils as common_utils
import DSBplot.utils.sam_utils as sam_utils
import DSBplot.utils.alignment_utils as alignment_utils
import DSBplot.utils.log_utils as log_utils


def check_consecutive_indel(ins_pos, del_pos):
  if (len(ins_pos) == 0) and (len(del_pos) == 0): # no in/dels
    return True
  
  ins_pos = list(sorted(set(ins_pos)))
  del_pos = list(sorted(set(del_pos)))
  
  if len(ins_pos) > 1: # all the insertions must be at the same position on the reference
    return False
  
  if len(del_pos) > 0:
    if len(del_pos) != (del_pos[-1] - del_pos[0] + 1): # is del_pos consecutive integers
      return False
    if (len(ins_pos) > 0):
      if not ((ins_pos[0] in del_pos) or ((ins_pos[0] + 1) in del_pos)):
        # the deletions must touch the left or right position of the insertion
        return False
  
  # all checks passed
  return True

def check_dsb_touches_indel(dsb_pos, ins_pos, del_pos):
  return (
    (dsb_pos in ins_pos) or # insertions at DSB
    (dsb_pos in del_pos) or # deletion on left of DSB
    ((dsb_pos + 1) in del_pos) # deletion on right of DSB
  )

def count_mismatch(ref_seq, read_seq, ref_pos):
  num_mismatch = 0
  for i in range(ref_pos - 1, min(len(ref_seq), len(read_seq))):
    if ref_seq[i] != read_seq[i]:
      num_mismatch += 1
  return num_mismatch

def check_insertion_special_case(
  ref_align,
  read_align,
  dsb_pos,
):
  """
    Check if all the insertions can be put at the DSB position.
    Take all the insertions in the alignment and place them
    at the DSB position and check if the number of substitutions does not increase.
    If so, return the new alignment, otherwise return None.
    Can only be used if the alignment has no deletions, otherwise an exception is raised.

    Parameters
    ----------
    ref_align  : the alignment string for the reference sequence.
    read_align : the alignment string for the read sequence.
    dsb_pos    : the position of the DSB on the reference sequence (1-based).

    Returns
    -------
    Tuple (None, None) if the insertions cannot be shifted to the DSB positions.
    Otherwise a tuple of the new reference and read alignments.
  """
  num_ins, num_del, num_subst = alignment_utils.count_variations(ref_align, read_align)

  if num_ins == 0:
    raise Exception('No insertions in alignment')
  if num_del > 0:
    raise Exception('Deletions in alignment')

  ref_seq = alignment_utils.get_orig_seq(ref_align)

  new_ref_align = ref_seq[:dsb_pos] + ('-' * num_ins) + ref_seq[dsb_pos:]

  # check that the number of substitutions has not increased
  _, _, new_num_subst = alignment_utils.count_variations(new_ref_align, read_align)
  if new_num_subst > num_subst:
    return None, None

  return new_ref_align, read_align # read_align remains unchanged

def check_deletion_special_case(
  ref_align,
  read_align,
  dsb_pos,
):
  """
    Check if all the deletions can be made to touch the DSB position.
    Take all the deletions in the alignment and place them
    in all possible ways to touch the DSB position and check if the number of substitutions does not increase.
    If so, return the new alignment, otherwise return None.
    Can only be used if the alignment has no deletions, otherwise an exeption is raised.
  
    Parameters
    ----------
    ref_align  : the alignment string for the reference sequence.
    read_align : the alignment string for the read sequence.
    dsb_pos    : the position of the DSB on the reference sequence (1-based).

    Returns
    -------
    Tuple (None, None) if the deletions cannot be shifted to the DSB positions.
    Otherwise a tuple of the new reference and read alignments.
  """
  num_ins, num_del, num_subst = alignment_utils.count_variations(ref_align, read_align)

  # check if there are insertions
  if num_del == 0:
    raise Exception('No deletions in alignment')
  if num_ins > 0:
    raise Exception('Insertions in alignment')

  read_seq = alignment_utils.get_orig_seq(read_align)

  new_read_align = None
  found_new = False
  # go through all possible ways of placing the deletions
  for del_start in range(dsb_pos - num_del, dsb_pos + 1):
    new_read_align = read_seq[:del_start] + ('-' * num_del) + read_seq[del_start:]
    _, _, new_num_subst = alignment_utils.count_variations(ref_align, new_read_align)
    if new_num_subst <= num_subst: # make sure that the number of substutitions has not increased
      found_new = True
      break

  if not found_new:
    return None, None

  return ref_align, new_read_align # ref_align remains unchanged

def parse_args():
  parser = argparse.ArgumentParser(
    description = 'Filter sequences based on the type and location of their variations.'
  )
  parser.add_argument(
    '--input',
    type = common_utils.check_file,
    help = (
      'Aligned SAM input file.' +
      ' Must be created with Bowtie2 (specific flags from Bowtie2 are used).'
      ' Every read must be aligned with exactly the same reference sequence.'
      ' Multiple files may be specified as long as they all use the same reference sequence.' +
      ' The output will contain a separate count columns for each input file.'
    ),
    nargs = '+',
    required = True,
  )
  parser.add_argument(
    '--output',
    nargs = 2,
    type = common_utils.check_file_output,
    required = True,
    help = 'Output CSV file name for the accepted and rejected reads, respectively.'
  )
  parser.add_argument(
    '--ref',
    type = common_utils.check_file,
    help = (
      'Reference sequence file.' +
      ' Should contain a single nucleotide sequence.' +
      ' Must be the same sequence as used for alignment.' +
      f' Must be in FASTA format ({", ".join(constants.FASTA_EXT)}) or' +
      ' text format (all other extensions).'
    ),
    required = True,
  )
  parser.add_argument(
    '--debug',
    type = common_utils.check_file_output,
    help = (
      'File to output debugging categories.' +
      ' If omitted, the messages are printed to the console (unless --quiet is used).' +
      ' If the extension is ".csv", the output is formatted as a CSV table, otherwise' +
      ' it is formatted as a text file.'
    ),
  )
  parser.add_argument(
    '--names',
    nargs = '+',
    help = (
      'Names to use as suffixes to the value columns of the output.' +
      ' Number of arguments must match the number of INPUT args.' +
      ' If not specified, becomes the names of the INPUT files.'
    ),
  )
  parser.add_argument(
    '--reads',
    type = int,
    nargs = '+',
    help = (
      'Total reads for each file.' +
      ' Must be the same number of arguments as the number of ' +
      ' INPUT files. If not provided, the total reads are' +
      ' are calculated by taking the sum of the "count" columns in each INPUT.'
    ),
  )
  parser.add_argument(
    '--dsb',
    type = int,
    required = True,
    help = (
      'Position on reference sequence immediately upstream of DSB site.' +
      ' I.e., the DSB is between 1-based positions DSB and DSB + 1.'
    ),
  )
  parser.add_argument(
    '--min_len',
    type = int,
    default = -1,
    help = (
      'Minimum length of read sequence to be considered.' +
      ' Reads shorter than this are discarded.' +
      ' Forced to be at least DSB + 1.'
    ),
  )
  parser.add_argument(
    '--rc',
    action = 'store_true',
    help = (
      'Set if the reads are expected to be reverse-complemented compared to the reference sequence.' +
      ' This is useful if the reads are from the opposite strand as the reference sequence.' +
      ' If this option is used, the alignments in the SAM file must have been aligned' +
      ' against the same reference sequence, and only alignments with the reverse-complement flag (16)' +
      ' will be accepted.'
    ),
  )
  parser.add_argument(
    '--quiet',
    help = 'Do not output log messages.',
    action = 'store_true',
  )
  args = vars(parser.parse_args())
  args['min_len'] = max(args['dsb'] + 1, args['min_len'])
  args = vars(parser.parse_args())
  if args['reads'] is not None:
    if len(args['reads']) != len(args['input']):
      raise Exception(
        'Number of total reads must match the number of input files.'
      )
  if args['names'] is not None:
    if len(args['input']) != len(args['names']):
      raise Exception(
        'Number of input files must match the number of names.'
      )
  return args

def do_filter(
  input,
  names,
  total_reads,
  ref_seq_file,
  output,
  output_rejected,
  dsb_pos,
  min_length,
  reverse_complement = False,
  quiet = True,
  debug_file = None,
):
  # read reference sequence from fasta file
  ref_seq = file_utils.read_seq(ref_seq_file)
  log_utils.log_input(ref_seq_file)

  if names is None:
    names = [file_names.get_file_name(x) for x in input]
  
  # For logging
  header = [0] * len(input)
  rejected_repeat = [0] * len(input)
  rejected_no_alignment = [0] * len(input)
  rejected_wrong_rc = [0] * len(input)
  rejected_pos_not_1 = [0] * len(input)
  rejected_not_consecutive = [0] * len(input)
  rejected_not_dsb_touch = [0] * len(input)
  rejected_not_consecutive_and_not_dsb_touch = [0] * len(input)
  rejected_too_short = [0] * len(input)
  accepted_repeat = [0] * len(input)
  accepted_deletion_special = [0] * len(input)
  accepted_insertion_special = [0] * len(input)
  accepted_insertion_and_deletion_special = [0] * len(input)
  accepted_indel_other = [0] * len(input)
  accepted_no_indel = [0] * len(input)

  # Common variables
  read_accepted = set()
  read_rejected = set()
  read_num_subst = {}
  read_aligned = {}
  read_cigar_old = {}
  read_cigar_new = {}
  read_debug = {}

  # Variables for each input file
  read_count_accepted = [defaultdict(lambda: 0) for _ in input]
  read_count_rejected = [defaultdict(lambda: 0) for _ in input]
  read_rank = [None for _ in input]
  total_lines = [0] * len(input)
  total_reads_1 = [0] * len(input)

  expected_rc_flag = 16 if reverse_complement else 0

  for i in range(len(input)): # Loop over input files
    total_lines[i] = file_utils.count_lines(input[i])

    with open(input[i], 'r') as in_h:
      log_utils.log_input(input[i])
      for line_num, line in enumerate(in_h, 1): # Filter reads
        if (line_num % 100000) == 1:
          if not quiet:
            log_utils.log(f"Progress: {i}: {line_num} / {total_lines[i]}")

        if line.startswith('@'): # header line of SAM
          header[i] += 1
          continue

        total_reads_1[i] += 1

        fields = line.rstrip().split('\t')
        mandatory, optional = sam_utils.parse_sam_fields(fields)
        read_seq = mandatory['SEQ']
        cigar = mandatory['CIGAR']

        if read_seq in read_accepted:
          read_count_accepted[i][read_seq] += 1
          accepted_repeat[i] += 1
          continue
        elif read_seq in read_rejected:
          read_count_rejected[i][read_seq] += 1
          rejected_repeat[i] += 1
          continue
        
        read_cigar_old[read_seq] = cigar

        if int(mandatory['FLAG']) & 4: # the read did not align at all
          read_rejected.add(read_seq)
          read_debug[read_seq] = 'no_alignment'
          read_aligned[read_seq] = False
          rejected_no_alignment[i] += 1
          read_count_rejected[i][read_seq] = 1
          continue
        elif (int(mandatory['FLAG']) & 16) != expected_rc_flag: # wrong RC flag
          read_rejected.add(read_seq)
          read_debug[read_seq] = 'wrong_rc'
          read_aligned[read_seq] = True
          rejected_wrong_rc[i] += 1
          read_count_rejected[i][read_seq] = 1
        else:
          read_aligned[read_seq] = True

        if int(mandatory['POS']) != 1:
          read_rejected.add(read_seq)
          read_debug[read_seq] = 'pos_not_1'
          rejected_pos_not_1[i] += 1
          read_count_rejected[i][read_seq] = 1
          continue

        if len(read_seq) < min_length:
          read_rejected.add(read_seq)
          read_debug[read_seq] = 'too_short'
          rejected_too_short[i] += 1
          read_count_rejected[i][read_seq] = 1
          continue

        # XG is the number of gap-extends (aka in/dels). 
        # XM if number of mismatches.
        # Both should always be present for aligned reads.
        num_indel_sam = int(optional['XG']['VALUE'])
        num_subst_sam = int(optional['XM']['VALUE'])

        ref_align, read_align = alignment_utils.get_alignment(ref_seq, read_seq, 1, cigar)

        ins_pos, del_pos, subst_pos = alignment_utils.get_var_pos(ref_align, read_align)
        num_ins = len(ins_pos)
        num_del = len(del_pos)
        num_subst = len(subst_pos)
        if num_indel_sam != (num_ins + num_del):
          raise Exception('Incorrect count of insertions and/or deletions')
        if num_subst_sam != num_subst:
          raise Exception('Incorrect count of substitutions')

        if (num_ins > 0) or (num_del > 0):
          consecutive = check_consecutive_indel(ins_pos, del_pos)
          dsb_touch = check_dsb_touches_indel(dsb_pos, ins_pos, del_pos)
          insertion_special_case = False
          deletion_special_case = False
          if not (consecutive and dsb_touch):
            # Note: The in/del special cases may fail to identify
            # certain edge cases when Bowtie picks an alignment that reduces the
            # number of substitutions by using mixed insertions and deletions, or
            # by using in/dels that are not continguous. This is because the
            # the special cases always try to reduce (or keep equal) the
            # number of substitutions. This could potentially be solved by
            # choosing different scoring parameters in Bowtie, (e.g., penalize
            # gap-open more heavily) but this is not implemented currently.
            if (num_ins > 0) and (num_del == 0):
              new_ref_align, new_read_align = check_insertion_special_case(
                ref_align,
                read_align,
                dsb_pos,
              )
              if new_ref_align is not None:
                ref_align = new_ref_align
                read_align = new_read_align
                insertion_special_case = True

            if (num_del > 0) and (num_ins == 0):
              new_ref_align, new_read_align = check_deletion_special_case(
                ref_align,
                read_align,
                dsb_pos,
              )
              if new_ref_align is not None:
                ref_align = new_ref_align
                read_align = new_read_align
                deletion_special_case = True
            if insertion_special_case or deletion_special_case:
              # if special case used, recompute info and do checks again
              cigar = alignment_utils.get_cigar(ref_align, read_align)
              ins_pos, del_pos, subst_pos = alignment_utils.get_var_pos(ref_align, read_align)
              num_ins = len(ins_pos)
              num_del = len(del_pos)
              num_subst = len(subst_pos)
              consecutive = check_consecutive_indel(ins_pos, del_pos)
              dsb_touch = check_dsb_touches_indel(dsb_pos, ins_pos, del_pos)

          if consecutive and dsb_touch:
            if insertion_special_case and deletion_special_case:
              read_debug[read_seq] = 'insertion_and_deletion_special'
              accepted_insertion_and_deletion_special[i] += 1
            elif insertion_special_case:
              read_debug[read_seq] = 'insertion_special'
              accepted_insertion_special[i] += 1
            elif deletion_special_case:
              read_debug[read_seq] = 'deletion_special'
              accepted_deletion_special[i] += 1
            else:
              read_debug[read_seq] = 'indel_other'
              accepted_indel_other[i] += 1
          else:
            if (not consecutive) and (not dsb_touch): 
              read_debug[read_seq] = 'not_consecutive_and_not_dsb_touch'
              rejected_not_consecutive_and_not_dsb_touch[i] += 1
            elif not consecutive:
              read_debug[read_seq] = 'not_consecutive'
              rejected_not_consecutive[i] += 1
            elif not dsb_touch:
              read_debug[read_seq] = 'not_dsb_touch'
              rejected_not_dsb_touch[i] += 1
            else:
              raise Exception('Impossible to reach.')
            read_rejected.add(read_seq)
            read_count_rejected[i][read_seq] = 1
            continue
        else:
          read_debug[read_seq] = 'no_indel'
          accepted_no_indel[i] += 1

        read_accepted.add(read_seq)
        read_cigar_new[read_seq] = cigar
        read_num_subst[read_seq] = num_subst
        read_count_accepted[i][read_seq] = 1
      # End of loop over reads
    # End with statement
    if len(read_count_accepted[i]) == 0:
      raise Exception('No reads captured. Check input file.')

    read_count = read_count_accepted[i].copy()
    read_count.update(read_count_rejected[i])
    read_count = sorted(read_count.items(), key=lambda x: x[1], reverse=True) # tuples (read_seq, count)
    read_rank[i] = {read_seq : rank for rank, (read_seq, _) in enumerate(read_count, 1)} 
  # End of loop over files

  if total_reads is None:
    total_reads = total_reads_1

  # Makes sure the read data is in the same order for all samples
  read_accepted = list(read_accepted)
  read_rejected = list(read_rejected)

  # Make the accepted read dataframe
  data_list = []
  for i in range(len(input)):
    read_seq_list = list(read_count_accepted[i].keys())
    data = pd.DataFrame({'seq': read_seq_list})
    data['count_' + names[i]] = [read_count_accepted[i][read_seq] for read_seq in read_seq_list]
    data['freq_' + names[i]] = data['count_' + names[i]] / total_reads[i]
    data['rank_' + names[i]] = [read_rank[i][read_seq] for read_seq in read_seq_list]
    data = data.set_index('seq').reindex(read_accepted)
    data['count_' + names[i]] = data['count_' + names[i]].fillna(0).astype(int)
    data['freq_' + names[i]] = data['freq_' + names[i]].fillna(0)
    data['rank_' + names[i]] = data['rank_' + names[i]].fillna(2**31 - 1).astype(int) # Int max
    data_list.append(data)
  data_common = pd.DataFrame({
    'seq': read_accepted,
    'debug': [read_debug[read_seq] for read_seq in read_accepted],
    'num_subst': [read_num_subst[read_seq] for read_seq in read_accepted],
    'cigar': [read_cigar_new[read_seq] for read_seq in read_accepted],
    'cigar_old': [read_cigar_old[read_seq] for read_seq in read_accepted],
  }).set_index('seq')
  data_accepted = (
    pd.concat([data_common] + data_list, join='outer', axis='columns')
    .fillna(0)
    .reset_index()
  )
  for col in data_accepted.columns:
    if col.startswith('count_'):
      data_accepted[col] = data_accepted[col].astype(int)
  data_accepted['num_subst'] = data_accepted['num_subst'].astype(int)
  data_accepted['freq_mean'] = (
    data_accepted[['freq_' + x for x in names]]
    .mean(axis='columns')
  )
  data_accepted = data_accepted[
    ['debug', 'num_subst', 'cigar', 'cigar_old'] +
    ['freq_mean'] +
    ['freq_' + x for x in names] +
    ['count_' + x for x in names] +
    ['rank_' + x for x in names] +
    ['seq']
  ]
  data_accepted = data_accepted.sort_values('freq_mean', ascending=False)

  # Make the rejected read data
  data_list = []
  for i in range(len(input)):
    read_seq_list = list(read_count_rejected[i].keys())
    data = pd.DataFrame({'seq': read_seq_list})
    data['count_' + names[i]] = [read_count_accepted[i][read_seq] for read_seq in read_seq_list]
    data['freq_' + names[i]] = data['count_' + names[i]] / total_reads[i]
    data['rank_' + names[i]] = [read_rank[i][read_seq] for read_seq in read_seq_list]
    data = data.set_index('seq').reindex(read_rejected)
    data['count_' + names[i]] = data['count_' + names[i]].fillna(0).astype(int)
    data['freq_' + names[i]] = data['freq_' + names[i]].fillna(0)
    data['rank_' + names[i]] = data['rank_' + names[i]].fillna(2**31 - 1).astype(int) # Int max
    data_list.append(data)
  data_common = pd.DataFrame({
    'seq': read_rejected,
    'debug': [read_debug[read_seq] for read_seq in read_rejected],
    'cigar_old': [read_cigar_old[read_seq] for read_seq in read_rejected],
    'aligned': [read_aligned[read_seq] for read_seq in read_rejected],
  }).set_index('seq')
  data_rejected = (
    pd.concat([data_common] + data_list, join='outer', axis='columns')
    .fillna(0)
    .sort_values(by='rank_' + names[0])
    .reset_index()
  )
  for col in data_rejected.columns:
    if col.startswith('count_'):
      data_rejected[col] = data_rejected[col].astype(int)
  data_rejected['freq_mean'] = (
    data_rejected[['freq_' + x for x in names]]
    .mean(axis='columns')
  )
  data_rejected = data_rejected[
    ['debug', 'cigar_old', 'aligned'] +
    ['freq_mean'] +
    ['freq_' + x for x in names] +
    ['count_' + x for x in names] +
    ['rank_' + x for x in names] +
    ['seq']
  ]
  data_rejected = data_rejected.sort_values('freq_mean', ascending=False)

  file_utils.write_csv(data_accepted, output)
  file_utils.write_csv(data_rejected, output_rejected)
  log_utils.log_output(output)
  log_utils.log_output(output_rejected)

  accepted_new = [0] * len(input)
  rejected_new = [0] * len(input)
  total_accepted = [0] * len(input)
  total_rejected = [0] * len(input)

  for i in range(len(input)):
    accepted_new[i] = (
      accepted_deletion_special[i] +
      accepted_insertion_special[i] +
      accepted_insertion_and_deletion_special[i] +
      accepted_indel_other[i] +
      accepted_no_indel[i]
    )
    total_accepted[i] = accepted_new[i] + accepted_repeat[i]

    rejected_new[i] = (
      rejected_no_alignment[i] +
      rejected_wrong_rc[i] +
      rejected_pos_not_1[i] +
      rejected_too_short[i] +
      rejected_not_consecutive[i] +
      rejected_not_dsb_touch[i] +
      rejected_not_consecutive_and_not_dsb_touch[i]
    )
    total_rejected[i] = rejected_new[i] + rejected_repeat[i]

    if (total_rejected[i] + total_accepted[i]) != total_reads_1[i]:
      raise Exception("accepted + rejected != total")

    if total_accepted[i] != sum(read_count_accepted[i].values()):
      raise Exception("Total accepted not summing")
    
    if total_rejected[i] != sum(read_count_rejected[i].values()):
      raise Exception("Total rejected not summing")

  debug_data = pd.DataFrame({
    'header': header,
    'total_reads': total_reads_1,
    'total_accepted': total_accepted,
    'accepted_new': accepted_new,
    'accepted_repeat': accepted_repeat,
    'accepted_deletion_special': accepted_deletion_special,
    'accepted_insertion_special': accepted_insertion_special,
    'accepted_insertion_and_deletion_special': accepted_insertion_and_deletion_special,
    'accepted_indel_other': accepted_indel_other,
    'accepted_no_indel': accepted_no_indel,
    'total_rejected': total_rejected,
    'rejected_new': rejected_new,
    'rejected_repeat': rejected_repeat,
    'rejected_no_alignment': rejected_no_alignment,
    'rejected_wrong_rc': rejected_wrong_rc,
    'rejected_pos_not_1': rejected_pos_not_1,
    'rejected_too_short': rejected_too_short,
    'rejected_not_consecutive': rejected_not_consecutive,
    'rejected_not_dsb_touch': rejected_not_dsb_touch,
    'rejected_not_consecutive_and_not_dsb_touch': rejected_not_consecutive_and_not_dsb_touch,
  }).T.rename_axis('debug')
  debug_data.columns = ['count_' + x for x in names]
  debug_data[['freq_' + x for x in names]] = debug_data[['count_' + x for x in names]].divide(total_reads_1, axis='columns')
  debug_data = debug_data.reset_index()
  debug_lines = [
    f'Header lines: ' + ', '.join([str(x) for x in header]),
    f'Total reads: ' + ', '.join([str(x) for x in total_reads_1]),
    f'    Accepted: '  + ', '.join([str(x) for x in total_accepted]),
    f'        Repeat: ' + ', '.join([str(x) for x in accepted_repeat]),
    f'        New: '  + ', '.join([str(x) for x in accepted_new]),
    f'            Insertion special case: ' + ', '.join([str(x) for x in accepted_insertion_special]),
    f'            Deletion special case: ' + ', '.join([str(x) for x in accepted_deletion_special]),
    f'            Insertion and deletion special case: ' + ', '.join([str(x) for x in accepted_insertion_and_deletion_special]),
    f'            In/del other: ' + ', '.join([str(x) for x in accepted_indel_other]),
    f'            No in/del: ' + ', '.join([str(x) for x in accepted_no_indel]),
    f'    Rejected: ' + ', '.join([str(x) for x in total_rejected]),
    f'        Repeat: ' + ', '.join([str(x) for x in rejected_repeat]),
    f'        New: ' + ', '.join([str(x) for x in rejected_new]),
    f'            No alignment: ' + ', '.join([str(x) for x in rejected_no_alignment]),
    f'            Wrong RC flag: ' + ', '.join([str(x) for x in rejected_wrong_rc]),
    f'            POS != 1: ' + ', '.join([str(x) for x in rejected_pos_not_1]),
    f'            Too short: ' + ', '.join([str(x) for x in rejected_too_short]),
    f'            Not consecutive: ' + ', '.join([str(x) for x in rejected_not_consecutive]),
    f'            DSB not touch: ' + ', '.join([str(x) for x in rejected_not_dsb_touch]),
    f'            Not consecutive and DSB not touch: ' + ', '.join([str(x) for x in rejected_not_consecutive_and_not_dsb_touch]),
  ]

  if (debug_file is None) and not quiet:
    for l in debug_lines:
      log_utils.log(l)
  elif debug_file is not None:
    if debug_file.endswith('.csv'):
      file_utils.write_csv(debug_data, debug_file)
    else:
      file_utils.make_parent_dir(debug_file)
      with open(debug_file, 'w') as debug_out:
        for l in debug_lines:
          debug_out.write(l + '\n')
    log_utils.log_output(debug_file)
  log_utils.blank_line()

def main(
  input,
  output,
  ref,
  debug,
  names,
  reads,
  dsb,
  min_len,
  rc,
  quiet,
):
  do_filter(
    input = input,
    output = output[0],
    output_rejected = output[1],
    ref_seq_file = ref,
    debug_file = debug,
    names = names,
    total_reads = reads,
    dsb_pos = dsb,
    min_length = len,
    reverse_complement = rc,
    quiet = quiet,
  )

if __name__ == '__main__':
  main(**parse_args())

# TODO: ADD MAX SUBSTITUTIONS TO ACCEPT
# TODO: ADD CONSECUTIVE IN/DEL TO ACCEPT
# TODO: ADD TOUCH DSB TO ACCEPT