import os

import shutil
import argparse
import glob

import DSBplot.utils.file_names as file_names
import DSBplot.utils.common_utils as common_utils
import DSBplot.utils.file_utils as file_utils
import DSBplot.utils.log_utils as log_utils
import DSBplot.utils.constants as constants

import DSBplot.preprocessing.filter as filter
import DSBplot.preprocessing.get_window as get_window
import DSBplot.preprocessing.get_variation as get_variation

STAGES = ['0_align', '1_filter', '2_window', '3_variation', '4_info']

PARAMS = {
  '--stages': {
    'type': str,
    'nargs': '*',
    'default': STAGES,
    'help': (
      'Stages to run. Choices: "0_align", "1_filter", "2_window", "3_variation", "4_info".' +
      ' Any prefix of a stage name can be used (e.g., "0" for "0_align").' +
      ' The stages must be run in the correct order indicated by their prefix numbers.' +
      ' 0_align: Align reads to reference sequence using Bowtie 2 ' +
      ' (may be omitted if a custom alignment approach is used; resulting SAM files must be put in the OUTPUT directory).' +
      ' 1_filter: Filter reads based on different criteria. Will use all SAM files placed in OUTPUT.' +
      ' 2_window: Extract the nucleotide sequences and alignment from around the DSB site.' +
      ' 3_variation: Split alignment windows into individual variations (used for plotting variation-position histograms).' +
      ' 4_info: Make a data info (metadata) file describing the experiment.'
    ),
  },
  '-o': {
    'type': common_utils.check_dir_output,
    'help': (
      'Working output directory.' +
      ' If running stages separately, the same directory must be given for all stages.' +
      ' The base part of the directory name is used as a unique identifier for the experiment' +
      ' so it should not conflict with other experiments.'
    ),
    'required': True,
    'metavar': 'OUTPUT',
    'dest': 'output',
  },
  '--names': {
    'nargs': '+',
    'type': str,
    'help': (
      'Identifiers to use for the input libraries.' +
      ' Must be given in the same order as the input files.' +
      ' If omitted, the names will be the input file names without the extension.' +
      ' The names must be in alphabetical order and should be legal parts of file names.'
    ),
  },
  '--no_align': {
    'action': 'store_true',
    'help': 'Shorthand for omitting the "0_align" stage (see "--stages").',
  },
  '--ref': filter.PARAMS['--ref'].copy(),
  '--dsb': filter.PARAMS['--dsb'].copy(),
  '-i': {
    'nargs': '+',
    'type': common_utils.check_file,
    'help': (
      'Input files of raw reads.' +
      ' The following file entensions are allowed: ' +
      ' FASTQ: ".fastq", ".fq."; FASTA: "fasta", ".fa", "fna";' +
      ' text: all others. FASTQ files are processed with the ' +
      ' Bowtie 2 flag "-q", FASTA files are processed with the Bowtie 2 flag "-f",' +
      ' and text files are processed with the Bowtie 2 flag "-r".' +
      ' Please see the Bowtie 2 manual at http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml.' +
      ' Each file is considered a repeat of the same experiment.' +
      ' If "0_align" is omitted and a custom alignnment approach is used, ' +
      ' then the aligned SAM files must be placed in OUTPUT directory.' +
      ' The input file basenames must be in alphabetical order and are used to'
      ' name the output SAM files (unless "--names" is given).'
    ),
    'metavar': 'INPUT',
    'dest': 'input',
  },
  '--bt2': {
    'type': str,
    'nargs': '+',
    'default': '',
    'help': 'Additional arguments to pass to Bowtie 2.',
  },
  '--min_len': filter.PARAMS['--min_len'].copy(),
  '--max_sub': filter.PARAMS['--max_sub'].copy(),
  '--rc': filter.PARAMS['--rc'].copy(),
  '--consec': filter.PARAMS['--consec'].copy(),
  '--touch': filter.PARAMS['--touch'].copy(),
  '--reads': filter.PARAMS['--reads'].copy(),
  '--window': get_window.PARAMS['--window'].copy(),
  '--anchor': get_window.PARAMS['--anchor'].copy(),
  '--anchor_var': get_window.PARAMS['--anchor_var'].copy(),
  '--label': {
    'type': str,
    'help': (
      'Label of the experiment to be used in plot legends.' +
      ' If not provided, the label is the basename of the OUTPUT directory.'
    ),
  },
  '--quiet': filter.PARAMS['--quiet'].copy(),
}

PARAMS['--stages']['required'] = False
PARAMS['-o']['required'] = True
PARAMS['--names']['required'] = False
PARAMS['--no_align']['required'] = False
PARAMS['--ref']['required'] = False
PARAMS['--dsb']['required'] = False
PARAMS['-i']['required'] = False
PARAMS['--bt2']['required'] = False
PARAMS['--min_len']['required'] = False
PARAMS['--max_sub']['required'] = False
PARAMS['--rc']['required'] = False
PARAMS['--consec']['required'] = False
PARAMS['--touch']['required'] = False
PARAMS['--reads']['required'] = False
PARAMS['--window']['required'] = False
PARAMS['--anchor']['required'] = False
PARAMS['--anchor_var']['required'] = False
PARAMS['--label']['required'] = False
PARAMS['--quiet']['required'] = False


PARAMS['-o']['help'] += ' Stages: all.'
PARAMS['--names']['help'] += ' Stages: "0_align".'
PARAMS['--ref']['help'] += ' Stages: "0_align", "1_filter",  "2_window".'
PARAMS['--dsb']['help'] += ' Stages: "1_filter", "2_window".'
PARAMS['-i']['help'] += ' Stages: "0_align".'
PARAMS['--bt2']['help'] += ' Stages: "0_align" (may be omitted because of default).'
PARAMS['--min_len']['help'] += ' Stages: "1_filter" (may be omitted because of default).'
PARAMS['--max_sub']['help'] += ' Stages: "1_filter" (may be omitted because of default).'
PARAMS['--rc']['help'] += ' Stages: "1_filter" (may be omitted because of default).'
PARAMS['--consec']['help'] += ' Stages: "1_filter" (may be omitted because of default).'
PARAMS['--touch']['help'] += ' Stages: "1_filter" (may be omitted because of default).'
PARAMS['--reads']['help'] += ' Stages: "1_filter" (may be omitted because of default).'
PARAMS['--window']['help'] += ' Stages: "2_window", "4_info" (may be omitted because of default).'
PARAMS['--anchor']['help'] += ' Stages: "2_window" (may be omitted because of default).'
PARAMS['--anchor_var']['help'] += ' Stages: "2_window" (may be omitted because of default).'
PARAMS['--label']['help'] += ' Stages: "4_info" (may be omitted because of default).'

def post_process_args(args):
  args = args.copy()
  args = filter.post_process_args(args)
  args = get_window.post_process_args(args)
  stages_old = tuple(args['stages'])
  stages = []
  for s in STAGES:
    if s.startswith(stages_old):
      # Some input stage is a prefix of s
      stages.append(s)
      break
  args['stages'] = stages
  if args['names'] is not None:
    if args['input'] is not None:
      if len(args['names']) != len(args['input']):
        raise Exception('Number of NAMES must be the same as the number of INPUTs.')
    if not all(x == y for x, y in zip(args['names'], sorted(args['names']))):
      raise Exception('NAMES must be in alphabetical order.')
  if args['input'] is not None:
    f_names = [file_names.get_file_name(x) for x in args['input']]
    if not all(x == y for x, y in zip(f_names, sorted(f_names))):
      raise Exception('INPUT must have file names in alphabetical order.')
  if args['no_align']:
    args['stages'] = [x for x in args['stages'] if (x != '0_align')]
  del args['no_align']
  return args

def parse_args():
  parser = argparse.ArgumentParser(
    description = (
      'Perform alignment and preprocessing for raw FASTQ data.' +
      ' This is script is broken into separate stages so that each stage' +
      ' can be run separately. However, the stages must be run in the correct order indicated' +
      ' by their prefix numbers. If running the stages separately, the value of OUTPUT' +
      ' must the same value on each separate invocation. Two experiments should not' +
      ' be given the same OUTPUT directory. See parameter "--stages" and the README for more information' +
      ' on the individual stages. The input files for each' +
      ' repeat of the experiment must be passed in alphabetical order.'
    ),
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
  )
  group_multi = parser.add_argument_group('Multiple stages')
  group_align = parser.add_argument_group('Stage 0_align')
  group_filter = parser.add_argument_group('Stage 1_filter')
  group_window = parser.add_argument_group('Stage 2_window')
  group_meta = parser.add_argument_group('Metadata')

  parser.add_argument('--stages', **PARAMS['--stages'])
  group_multi.add_argument('-o', **PARAMS['-o'])
  group_multi.add_argument('--names', **PARAMS['--names'])
  group_multi.add_argument('--ref', **PARAMS['--ref'])
  group_multi.add_argument('--dsb', **PARAMS['--dsb'])
  group_multi.add_argument('--quiet', **PARAMS['--quiet'])
  group_multi.add_argument('--no_align', **PARAMS['--no_align'])
  group_align.add_argument('-i', **PARAMS['-i'])
  group_align.add_argument('--bt2', **PARAMS['--bt2'])
  group_filter.add_argument('--min_len', **PARAMS['--min_len'])
  group_filter.add_argument('--max_sub', **PARAMS['--max_sub'])
  group_filter.add_argument('--rc', **PARAMS['--rc'])
  group_filter.add_argument('--consec', **PARAMS['--consec'])
  group_filter.add_argument('--touch', **PARAMS['--touch'])
  group_filter.add_argument('--reads', **PARAMS['--reads'])
  group_window.add_argument('--window', **PARAMS['--window'])
  group_window.add_argument('--anchor', **PARAMS['--anchor'])
  group_window.add_argument('--anchor_var', **PARAMS['--anchor_var'])
  group_meta.add_argument('--label', **PARAMS['--label'])

  return post_process_args(vars(parser.parse_args()))

def do_0_align(
  input,
  names,
  output,
  ref,
  bt2,
):
  if input is None:
    raise Exception('INPUT must be provided for stage "0_align".')
  if ref is None:
    raise Exception('REF must be provided for stage "0_align".')
  if names is None:
    names = [file_names.get_file_name(x) for x in input]
  bt2 = ' '.join(bt2)
  bowtie2_build_file = file_names.bowtie2_build(output)
  file_utils.make_parent_dir(bowtie2_build_file)
  log_utils.log_output('Bowtie2 build file: ' + bowtie2_build_file)
  os.system(f'bowtie2-build-s {ref} {bowtie2_build_file} --quiet')

  for i in range(len(input)):
    sam_file = file_names.sam_file(output, names[i])
    file_utils.make_parent_dir(sam_file)
    input_ext = os.path.splitext(input[i])[1]
    if input_ext in ['.fastq', '.fq']:
      flags = '-q'
    elif input_ext in ['.fasta', '.fa', '.fna']:
      flags = '-f'
    else:
      flags = '-r'
    bowtie2_command = f'bowtie2-align-s --no-hd {flags} {bt2} -x {bowtie2_build_file} {input[i]} -S {sam_file} --quiet'
    log_utils.log('Bowtie 2 command: ' + bowtie2_command)
    os.system(bowtie2_command)
  log_utils.blank_line()

def do_1_filter(
  names,
  output,
  reads,
  ref,
  dsb,
  min_len,
  max_sub,
  rc,
  consec,
  touch,
  quiet,
):
  if ref is None:
    raise Exception('REF must be provided for stage "1_filter".')
  if dsb is None:
    raise Exception('DSB must be provided for stage "1_filter".')
  if min_len is None:
    min_len = dsb + 1

  if names is None:
    # Infer names from the SAM files.
    names = sorted([
      file_names.get_file_name(x)
      for x in glob.glob(os.path.join(output, '*.sam'))
    ])
  input = [file_names.sam_file(output, x) for x in names]

  min_len = max(dsb + 1, min_len)
  filter.main(
    input = input,
    output = [
      file_names.filter(output, 'accepted'),
      file_names.filter(output, 'rejected'),
    ],
    debug = file_names.filter(output, 'debug'),
    names = names,
    reads = reads,
    ref = ref,
    dsb = dsb,
    min_len = min_len,
    max_sub = max_sub,
    rc = rc,
    consec = consec,
    touch = touch,
    quiet = quiet,
  )

def do_2_window(
  output,
  ref,
  dsb,
  window,
  anchor,
  anchor_var,
):
  if ref is None:
    raise Exception('REF must be provided for stage "2_window".')
  if dsb is None:
    raise Exception('DSB must be provided for stage "2_window".')
  if window is None:
    raise Exception('WINDOW must be provided for stage "2_window".')
  if anchor is None:
    raise Exception('ANCHOR must be provided for stage "2_window".')
  if anchor_var is None:
    raise Exception('ANCHOR_VAR must be provided for stage "2_window".')

  for subst_type in constants.SUBST_TYPES:
    get_window.main(
      input = file_names.filter(output, 'accepted'),
      output = file_names.window(output, subst_type),
      ref = ref,
      dsb = dsb,
      window = window,
      anchor = anchor,
      anchor_var = anchor_var,
      sub = subst_type,
    )

def do_3_variation(output):
  for subst_type in constants.SUBST_TYPES:
    get_variation.main(
      input = file_names.window(output, subst_type),
      output = file_names.variation(output, subst_type),
    )

def do_4_info(
  output,
  ref,
  dsb,
  window,
  label,
):
  if output is None:
    raise Exception('OUTPUT must be provided for stage "4_info".')
  if ref is None:
    raise Exception('REF must be provided for stage "4_info".')
  if dsb is None:
    raise Exception('DSB must be provided for stage "4_info".')
  if window is None:
    raise Exception('WINDOW must be provided for stage "4_info".')
  name = file_names.get_file_name(output)
  if label is None:
    label = name
  ref_seq = file_utils.read_seq(ref)
  data_info = common_utils.make_data_info(
    format = 'individual',
    names = [name],
    labels = [label],
    ref_seqs = [ref_seq],
    ref_seq_window = get_window.get_ref_seq_window(
      ref_seq = ref_seq,
      dsb_pos = dsb,
      window_size = window,
    ),
  )
  out_file = file_names.data_info(output)
  log_utils.log_output(out_file)
  file_utils.write_json(data_info, out_file)

def do_stages(
  input,
  names,
  output,
  ref,
  dsb,
  min_len,
  max_sub,
  rc,
  consec,
  touch,
  window,
  anchor,
  anchor_var,
  reads,
  label,
  bt2,
  quiet,
  stages = STAGES,
):
  # Copy reference sequence file to output directory.
  if ref is not None:
    shutil.copy(ref, file_names.ref_seq_file(output))

  if '0_align' in stages:
    do_0_align(
      input = input,
      names = names,
      output = output,
      ref = ref,
      bt2 = bt2,
    )

  if '1_filter' in stages:
    do_1_filter(
      names = names,
      output = output,
      reads = reads,
      ref = ref,
      dsb = dsb,
      min_len = min_len,
      max_sub = max_sub,
      rc = rc,
      consec = consec,
      touch = touch,
      quiet = quiet,
    )

  if '2_window' in stages:
    do_2_window(
      output = output,
      ref = ref,
      dsb = dsb,
      window = window,
      anchor = anchor,
      anchor_var = anchor_var,
    )

  if '3_variation' in stages:
    do_3_variation(output = output)

  if '4_info' in stages:
    do_4_info(
      output = output,
      ref = ref,
      dsb = dsb,
      window = window,
      label = label,
    )

def main(
  input,
  output,
  names,
  ref,
  dsb,
  min_len,
  max_sub,
  rc,
  consec,
  touch,
  window,
  anchor,
  anchor_var,
  reads,
  label,
  bt2,
  quiet,
  stages = STAGES,
):
  do_stages(
    input = input,
    names  = names,
    output = output,
    ref = ref,
    dsb = dsb,
    min_len = min_len,
    max_sub = max_sub,
    rc = rc,
    consec = consec,
    touch = touch,
    window = window,
    anchor = anchor,
    anchor_var = anchor_var,
    reads = reads,
    label = label,
    bt2 = bt2,
    quiet = quiet,
    stages = stages,
  )

# This allows the "DSBplot-preprocess" command to be run from the command line.
def entry_point():
  main(**parse_args())

if __name__ == '__main__':
  entry_point()
