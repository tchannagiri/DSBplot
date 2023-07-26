import shutil
import csv
import os
import pandas as pd
import re
import json

import DSBplot.utils.constants as constants

def make_parent_dir(file_name):
  dir_name = os.path.dirname(file_name)
  if dir_name != '':
    os.makedirs(os.path.dirname(file_name), exist_ok=True)

def write_csv(data, file, **args):
  if type(file) == str:
    make_parent_dir(file)
  data.to_csv(
    file,
    na_rep = 'NA',
    quoting = csv.QUOTE_NONNUMERIC,
    index = args.get('index', False),
    lineterminator = '\n',
  )

def read_csv(file):
  return pd.read_csv(
    file,
    index_col = False,
    keep_default_na = False,
    na_values = 'NA',
  )

def read_json(file):
  with open(file) as input:
    return json.load(input)

def count_lines(file):
  """Get the number of lines in the file."""
  with open(file) as input:
    return sum(1 for _ in input)

def copy(file_src, file_dst):
  make_parent_dir(file_dst)
  shutil.copy(file_src, file_dst)

def write_pyplot(figure, file):
  if type(file) == str:
    make_parent_dir(file)
  figure.savefig(file)

def write_plotly(figure, file):
  if type(file) == str:
    make_parent_dir(file)
  ext = os.path.splitext(file)[1]
  if ext == '.html':
    figure.write_html(file)
  else:
    figure.write_image(file, engine='kaleido')

def write_json(data, file):
  if type(file) == str:
    make_parent_dir(file)
  with open(file, 'w') as output:
    json.dump(data, output, indent=2)

def read_fasta_seq(fasta_file):
  """
    Read the sequence from a FASTA file.
    The FASTA file must have a single sequence or an exception will be raised.
    The sequence in the FASTA file must only contain letters A, C, G, and T or
    an exception will be raised.

    Parameters
    ----------
    fasta_file : input FASTA file name.

    Returns
    -------
    The sequence in the FASTA file.
  """
  with open(fasta_file) as fasta_h:
    lines = fasta_h.readlines()
    if len(lines) == 0:
      raise Exception(f'The reference FASTA file {fasta_file} is empty.')
    lines = [line.rstrip() for line in lines]
    lines = [line for line in lines if line != '']
    if len(lines) == 0:
      raise Exception(f'The reference FASTA file {fasta_file} contains only empty lines.')
    seq = ''
    i = 0
    if lines[i][0] != '>':
      raise Exception(f'Expected the sequence header ">" in FASTA file {fasta_file}.')
    i += 1
    while (i < len(lines)) and (lines[i][0] != '>'):
      seq += lines[i]
      i += 1
    if (i < len(lines)) and (lines[i][0] == '>'):
      raise Exception(f'The FASTA file {fasta_file} contains multiple sequences.')
    if not(all([x in 'ACGT' for x in seq])):
      raise Exception(f'The sequence in FASTA file {fasta_file} contains invalid characters.')
    if len(seq) == 0:
      raise Exception(f'The sequence in FASTA file {fasta_file} is empty.')
    return seq

def read_text_seq(text_file):
    """
    Read the sequence from a text file.
    The file must contain only A, C, G, T or space and newline characters or
    an exception will be raised. All the space and newline characters will be
    removed.

    Parameters
    ----------
    text_file : input text file name.

    Returns
    -------
    The sequence in the text file.
  """
    with open(text_file, 'r') as input:
      seq = input.read()
      seq = re.sub(r'\s', '', seq)
      seq = seq.replace('\n', '')
      seq = seq.replace(' ', '')
      if not(all([x in 'ACGT' for x in seq])):
          raise Exception(f'The sequence in text file {text_file} contains invalid characters.')
      if len(seq) == 0:
          raise Exception(f'The sequence in text file {text_file} is empty.')
      return seq

def read_seq(seq_file):
  """
    Read the sequence from a FASTA or text file.
    All files with extensions .fa, .fasta, and .fna are considered as FASTA files
    and all other files are considered as text files.
    See the respective documentation in read_fast_seq() and read_text_seq()
    for the requirements of each format.

    Parameters
    ----------
    seq_file : input FASTA or text file name.

    Returns
    -------
    The sequence in the FASTA or text file.
  """
  ext = os.path.splitext(seq_file)[1].replace('.', '')
  if ext in constants.FASTA_EXT:
    return read_fasta_seq(seq_file)
  else:
    return read_text_seq(seq_file)
