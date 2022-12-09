# get reference sequence
def read_fasta_seq(fasta_file):
  """
    Read the sequence from a FASTA file.
    The FASTA file must have a single sequence or an exception will be raised.

    Parameters
    ----------
    fasta_file : input FASTA file

    Returns
    -------
    seq : string,
      the sequence in the input file
  """
  fasta_file.readline() # skip the first line which should be the seq id
  seq = ''
  for line in fasta_file:
    if line[0] == '>':
      raise Exception(f'There should be only one sequence in reference FASTA file {fasta_file.name}!')
    seq += line.rstrip()
  return seq
