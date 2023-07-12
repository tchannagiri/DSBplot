import os
import argparse

import pandas as pd

import DSBplot.utils.constants as constants
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
    description = 'Filter sequences having mutations near DSB site.'
  )
  parser.add_argument(
    '--ref_seq_file',
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
    '--sam_file',
    type = common_utils.check_file,
    help = (
      'Aligned SAM input file.' +
      ' Must be created with Bowtie2 (specific flags from Bowtie2 are used).'
      ' Every read must be aligned with exactly the same reference sequence.'
    ),
    required = True,
  )
  parser.add_argument(
    '--output',
    type = common_utils.check_file_output,
    required = True,
    help = 'Output TSV file name.'
  )
  parser.add_argument(
    '--output_rejected',
    type = common_utils.check_file_output,
    required = True,
    help = 'Output file for rejected reads.'
  )
  parser.add_argument(
    '--dsb_pos',
    type = int,
    required = True,
    help = (
      'Position on reference sequence immediately upstream of DSB site.' +
      ' I.e. the DSB is between 1-based positions DSB_POS and DSB_POS + 1.'
    ),
  )
  parser.add_argument(
    '--min_length',
    type = int,
    default = -1,
    help = (
      'Minimum length of read sequence to be considered.' +
      ' Reads shorter than this are discarded.' +
      ' Forced to be at least DSB_POS + 1.'
    ),
  )
  parser.add_argument(
    '--quiet',
    help = 'Do not output log messages.',
    action = 'store_true',
  )
  parser.add_argument(
    '--debug_file',
    type = common_utils.check_file_output,
    help = (
      'File to output debugging messages.' +
      'If omitted, the messages are printed to the console (unless --quiet is used).'
    ),
  )
  args = vars(parser.parse_args())
  args['min_length'] = max(args['dsb_pos'] + 1, args['min_length'])
  return args

def main(
  ref_seq_file,
  sam_file,
  output,
  output_rejected,
  dsb_pos,
  min_length,
  quiet = True,
  debug_file = None,
):
  # read reference sequence from fasta file
  ref_seq = file_utils.read_seq(ref_seq_file)
  log_utils.log_input(ref_seq_file)
  
  # For logging
  rejected_header = 0
  rejected_repeat = 0
  rejected_no_alignment = 0
  rejected_pos_not_1 = 0
  rejected_not_consecutive = 0
  rejected_dsb_not_touch = 0
  rejected_not_consecutive_and_dsb_not_touch = 0
  rejected_too_short = 0
  accepted_repeat = 0
  accepted_deletion_special = 0
  accepted_insertion_special = 0
  accepted_insertion_and_deletion_special = 0
  accepted_indel_other = 0
  accepted_no_indel = 0

  # categorize
  read_count_accepted = {}
  read_count_rejected = {}
  read_num_subst = {}
  read_aligned = {}
  read_cigar_old = {}
  read_cigar_new = {}
  total_lines = file_utils.count_lines(sam_file)
  total_reads = 0
  with open(sam_file, 'r') as sam_file_h:
    log_utils.log_input(sam_file)
    for line_num, line in enumerate(sam_file_h, 1):
      if (line_num % 100000) == 1:
        if not quiet:
          log_utils.log(f"Progress: {line_num} / {total_lines}")

      if line.startswith('@'): # header line of SAM
        rejected_header += 1
        continue

      total_reads += 1

      fields = line.rstrip().split('\t')
      mandatory, optional = sam_utils.parse_sam_fields(fields)
      read_seq = mandatory['SEQ']
      cigar = mandatory['CIGAR']
      read_cigar_old[read_seq] = cigar

      if read_seq in read_count_accepted:
        read_count_accepted[read_seq] += 1
        accepted_repeat += 1
        continue
      elif read_seq in read_count_rejected:
        read_count_rejected[read_seq] += 1
        rejected_repeat += 1
        continue

      if int(mandatory['FLAG']) & 4: # the read did not align at all
        rejected_no_alignment += 1
        read_count_rejected[read_seq] = 1
        read_aligned[read_seq] = False
        continue
      else:
        read_aligned[read_seq] = True

      if int(mandatory['POS']) != 1:
        rejected_pos_not_1 += 1
        read_count_rejected[read_seq] = 1
        continue

      if len(read_seq) < min_length:
        rejected_too_short += 1
        read_count_rejected[read_seq] = 1
        continue

      # XG is the number of gap-extends (aka in/dels). 
      # XM if number of mismatches.
      # Both should always be present for aligned reads.
      num_indel_sam = int(optional['XG']['VALUE'])
      num_subst_sam = int(optional['XM']['VALUE'])

      ref_align, read_align = alignment_utils.get_alignment(ref_seq, read_seq, 1, cigar)

      ins_pos, del_pos, subst_pos = alignment_utils.get_variation_pos(ref_align, read_align)
      num_ins = len(ins_pos)
      num_del = len(del_pos)
      num_subst = len(subst_pos)
      if num_indel_sam != (num_ins + num_del):
          raise Exception('Incorrect count of insertions and/or deletions')
      if num_subst_sam != num_subst:
          raise Exception('Incorrect count of substitutions')

      if (num_ins > 0) or (num_del > 0):
        consecutive = check_consecutive_indel(ins_pos, del_pos)
        dsb_touches = check_dsb_touches_indel(dsb_pos, ins_pos, del_pos)
        insertion_special_case = False
        deletion_special_case = False
        if not (consecutive and dsb_touches):
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
            ins_pos, del_pos, subst_pos = alignment_utils.get_variation_pos(ref_align, read_align)
            num_ins = len(ins_pos)
            num_del = len(del_pos)
            num_subst = len(subst_pos)
            consecutive = check_consecutive_indel(ins_pos, del_pos)
            dsb_touches = check_dsb_touches_indel(dsb_pos, ins_pos, del_pos)

        if consecutive and dsb_touches:
          if insertion_special_case and deletion_special_case:
            accepted_insertion_and_deletion_special += 1
          elif insertion_special_case:
            accepted_insertion_special += 1
          elif deletion_special_case:
            accepted_deletion_special += 1
          else:
            accepted_indel_other += 1
        else:
          if (not consecutive) and (not dsb_touches): 
            rejected_not_consecutive_and_dsb_not_touch += 1
          elif not consecutive:
            rejected_not_consecutive += 1
          elif not dsb_touches:
            rejected_dsb_not_touch += 1
          else:
            raise Exception('Impossible to reach.')
          read_count_rejected[read_seq] = 1
          continue
      else:
        accepted_no_indel += 1

      read_count_accepted[read_seq] = 1
      read_cigar_new[read_seq] = cigar
      read_num_subst[read_seq] = num_subst

  if len(read_count_accepted) == 0:
    raise Exception('No reads captured. Check input file.')

  read_count = read_count_accepted.copy()
  read_count.update(read_count_rejected)
  read_count = sorted(read_count.items(), key=lambda x: x[1], reverse=True) # tuples (read_seq, count)
  read_rank = {read_seq : rank for rank, (read_seq, _) in enumerate(read_count, 1)} 

  read_seq_list = sorted(
    read_count_accepted.keys(),
    key = lambda x: read_count_accepted[x],
    reverse = True
  )

  file_utils.make_parent_dir(output)
  with open(output, 'w') as output_h:
    output_h.write('Rank\tCount\tNum_Subst\tCIGAR\tCIGAR_old\tSequence\n')
    for read_seq in read_seq_list:
      rank = read_rank[read_seq]
      cigar = read_cigar_new[read_seq]
      cigar_old = read_cigar_old[read_seq]
      count = read_count_accepted[read_seq]
      num_subst = read_num_subst[read_seq]
      output_h.write(f'{rank}\t{count}\t{num_subst}\t{cigar}\t{cigar_old}\t{read_seq}\n')

  read_seq_list = sorted(
    read_count_rejected.keys(),
    key = lambda x: read_count_rejected[x],
    reverse = True
  )
  with open(output_rejected, 'w') as output_h:
    output_h.write('Rank\tCount\tAligned\tSequence\n')
    for read_seq in read_seq_list:
      rank = read_rank[read_seq]
      count = read_count_rejected[read_seq]
      aligned = int(read_aligned[read_seq])
      output_h.write(f'{rank}\t{count}\t{aligned}\t{read_seq}\n')

  log_utils.log_output(output)
  log_utils.log_output(output_rejected)

  accepted_new = (
    accepted_deletion_special +
    accepted_insertion_special +
    accepted_insertion_and_deletion_special +
    accepted_indel_other +
    accepted_no_indel
  )
  total_accepted = accepted_new + accepted_repeat

  rejected_new = (
    rejected_no_alignment + 
    rejected_pos_not_1 +
    rejected_too_short +
    rejected_not_consecutive +
    rejected_dsb_not_touch +
    rejected_not_consecutive_and_dsb_not_touch
  )
  total_rejected = rejected_new + rejected_repeat

  if (total_rejected + total_accepted) != total_reads:
    raise Exception("accepted + rejected != total")

  if total_accepted != sum(read_count_accepted.values()):
    raise Exception("Total accepted not summing")
  
  if total_rejected != sum(read_count_rejected.values()):
    raise Exception("Total rejected not summing")

  debug_lines = [
    f'Header lines: {rejected_header}',
    f'Total reads: {total_reads}',
    f'    Accepted: {total_accepted}',
    f'        Repeat: {accepted_repeat}',
    f'        New: {accepted_new}',
    f'            Insertion special case: {accepted_insertion_special}',
    f'            Deletion special case: {accepted_deletion_special}',
    f'            Insertion and deletion special case: {accepted_insertion_and_deletion_special}',
    f'            In/del other: {accepted_indel_other}',
    f'            No in/del: {accepted_no_indel}',
    f'    Rejected: {total_rejected}',
    f'        Repeat: {rejected_repeat}',
    f'        New: {rejected_new}',
    f'            No alignment: {rejected_no_alignment}',
    f'            POS != 1: {rejected_pos_not_1}',
    f'            Too short: {rejected_too_short}',
    f'            Not consecutive: {rejected_not_consecutive}',
    f'            DSB not touch: {rejected_dsb_not_touch}',
    f'            Not consecutive and DSB not touch: {rejected_not_consecutive_and_dsb_not_touch}',
  ]

  if (debug_file is None) and not quiet:
    for l in debug_lines:
      log_utils.log(l)
  elif debug_file is not None:
    file_utils.make_parent_dir(debug_file)
    with open(debug_file, 'w') as debug_out:
      for l in debug_lines:
        debug_out.write(l + '\n')
      debug_out.close()
      log_utils.log_output(debug_file)
  log_utils.blank_line()

if __name__ == '__main__':
  main(**parse_args())
