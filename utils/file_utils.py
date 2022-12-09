import csv
import os
import pandas as pd

def make_parent_dir(file_name):
  os.makedirs(os.path.dirname(file_name), exist_ok=True)

def write_tsv(data, file, **args):
  if type(file) == str:
    make_parent_dir(file)
  data.to_csv(
    file,
    sep = '\t',
    na_rep = 'NA',
    quoting = csv.QUOTE_NONNUMERIC,
    index = args.get('index', False),
    lineterminator = '\n',
  )

def read_tsv(file):
  return pd.read_csv(
    file,
    index_col = False,
    keep_default_na = False,
    na_values = 'NA',
    sep = '\t',
  )

def read_tsv_dict(file):
  """
    Read a single row tsv as a dict.
  """
  return read_tsv(file).T[0].to_dict()

def count_lines(file):
  """Get the number of lines in the file."""
  with open(file) as input:
    return sum(1 for _ in input)

def write_pyplot(figure, file):
  if type(file) == str:
    make_parent_dir(file)
  figure.savefig(file)


def write_plotly(figure, file):
  if type(file) == str:
    make_parent_dir(file)
  ext = os.path.splitext(file)[1]
  if ext == '.png':
    figure.write_image(file, engine='kaleido')
  elif ext == '.html':
    figure.write_html(file)
  else:
    raise Exception('Unknown extension: ' + str(ext))
