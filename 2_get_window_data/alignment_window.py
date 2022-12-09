import kmer_utils
import numpy as np

def get_alignment_window(
  ref_align,
  read_align,
  dsb_pos,
  window_size,
  anchor_size,
  anchor_mismatch_limit,
):
  """
    Get the part of the alignment in a window around the DSB.

    Parameters
    ----------
    ref_align: the reference sequence alignment.
    read_align: the read sequence alignment.
    dsb_pos: the 1-based DSB position on the reference.
    window_size: size of the window to extract.
    anchor_size: the size of the anchors that must match on the left and right of the DSB.
    anchor_mismatch_limit: the maximum number of mismatches allowed on the left and right anchors.
      The mismatch limit is checked separately on the left or right anchors.
    
    Returns
    -------
    A tuple (ref_align, read_align):
      ref_align: the part of the reference alignment in the window.
      read_align: the part of the read alignment in the window.
    If the extraction fails due to too many anchor mismatches or the read
    is too short to cross dsb_pos, the tuple (None, None) is returned.
  """

  ref_i = 1 # index on the original reference sequence

  ref_align_window = ''
  read_align_window = ''

  window_start = dsb_pos - window_size + 1
  window_end = dsb_pos + window_size
  left_anchor_start = window_start - anchor_size
  left_anchor_end = window_start - 1
  right_anchor_start = window_end + 1
  right_anchor_end = window_end + anchor_size

  left_anchor_mismatches = 0
  right_anchor_mismatches = 0
  for i in range(min(len(ref_align), len(read_align))):
    if ref_i in range(left_anchor_start, left_anchor_end + 1):
      # Check the mismatches/in/dels on the left anchor
      if (ref_align[i] == '-') or (read_align[i] == '-'):
        left_anchor_mismatches = np.inf
      elif ref_align[i] != read_align[i]:
        left_anchor_mismatches += 1
    elif ref_i in range(right_anchor_start, right_anchor_end + 1):
      # Check the anchor/in/dels on the right anchor
      if (ref_align[i] == '-') or (read_align[i] == '-'):
        right_anchor_mismatches = np.inf
      elif ref_align[i] != read_align[i]:
        right_anchor_mismatches += 1
    elif ref_i in range(window_start, window_end + 1):
      # extract the window around the DSB
      ref_align_window += ref_align[i]
      read_align_window += read_align[i]

    # increment counters
    if ref_align[i] != '-':
      ref_i += 1

    if ref_i > right_anchor_end:
      break

  if ref_i <= right_anchor_end:
    # read was not long enough to align across window and anchor
    return None, None

  if (
    (left_anchor_mismatches > anchor_mismatch_limit) or
    (right_anchor_mismatches > anchor_mismatch_limit)
  ):
    return None, None

  return ref_align_window, read_align_window
