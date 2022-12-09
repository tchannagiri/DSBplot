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
