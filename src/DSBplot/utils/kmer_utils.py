def get_kmer_index(kmer):
  nuc_index = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3,
  }
  index = 0
  for x in kmer:
    index *= len(nuc_index)
    index += nuc_index[x]
  return index

def get_num_kmers(kmer_size):
  return 4**kmer_size

def reverse_complement(nucleotides):
  rev_map = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    '-': '-',
    'N': 'N',
  }
  return ''.join(reversed([rev_map[x] for x in nucleotides]))