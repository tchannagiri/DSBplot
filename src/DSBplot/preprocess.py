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
  '--no_align': {
    'action': 'store_true',
    'help': 'Shorthand for omitting the "0_align" stage (see "--stages").',
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
    'dest': 'input_list',
  },
  '--bt2': {
    'type': str,
    'nargs': '+',
    'default': [''],
    'help': 'Additional arguments to pass to Bowtie 2.',
    'dest': 'bowtie2_args',
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
    'dest': 'library_names',
  },
  '--label': {
    'type': str,
    'help': (
      'Label of the experiment to be used in plot legends.' +
      ' If not provided, the label is the basename of the OUTPUT directory.'
    ),
  },
  '--ref': filter.PARAMS['--ref'].copy(),
  '--dsb': filter.PARAMS['--dsb'].copy(),
  '--min_len': filter.PARAMS['--min_len'].copy(),
  '--max_sub': filter.PARAMS['--max_sub'].copy(),
  '--rc': filter.PARAMS['--rc'].copy(),
  '--consec': filter.PARAMS['--consec'].copy(),
  '--touch': filter.PARAMS['--touch'].copy(),
  '--reads': filter.PARAMS['--reads'].copy(),
  '--window': get_window.PARAMS['--window'].copy(),
  '--anchor': get_window.PARAMS['--anchor'].copy(),
  '--anchor_vars': get_window.PARAMS['--anchor_vars'].copy(),
  '--quiet': filter.PARAMS['--quiet'].copy(),
}

PARAMS['--stages']['required'] = False
PARAMS['--no_align']['required'] = False
PARAMS['-o']['required'] = True
PARAMS['-i']['required'] = False
PARAMS['--bt2']['required'] = False
PARAMS['--ref']['required'] = False
PARAMS['--names']['required'] = False
PARAMS['--reads']['required'] = False
PARAMS['--dsb']['required'] = False
PARAMS['--min_len']['required'] = False
PARAMS['--max_sub']['required'] = False
PARAMS['--rc']['required'] = False
PARAMS['--consec']['required'] = False
PARAMS['--touch']['required'] = False
PARAMS['--window']['required'] = False
PARAMS['--anchor']['required'] = False
PARAMS['--anchor_vars']['required'] = False
PARAMS['--label']['required'] = False
PARAMS['--quiet']['required'] = False

PARAMS['-o']['help'] += ' Stages: all.'
PARAMS['-i']['help'] += ' Stages: "0_align".'
PARAMS['--bt2']['help'] += ' Stages: "0_align" (may be omitted because of default).'
PARAMS['--ref']['help'] += ' Stages: "0_align", "1_filter",  "2_window".'
PARAMS['--names']['help'] += ' Stages: "0_align".'
PARAMS['--reads']['help'] += ' Stages: "1_filter" (may be omitted because of default).'
PARAMS['--dsb']['help'] += ' Stages: "1_filter", "2_window".'
PARAMS['--min_len']['help'] += ' Stages: "1_filter" (may be omitted because of default).'
PARAMS['--max_sub']['help'] += ' Stages: "1_filter" (may be omitted because of default).'
PARAMS['--rc']['help'] += ' Stages: "1_filter" (may be omitted because of default).'
PARAMS['--consec']['help'] += ' Stages: "1_filter" (may be omitted because of default).'
PARAMS['--touch']['help'] += ' Stages: "1_filter" (may be omitted because of default).'
PARAMS['--quiet']['help'] += ' Stages: "1_filter" (may be omitted because of default).'
PARAMS['--window']['help'] += ' Stages: "2_window", "4_info" (may be omitted because of default).'
PARAMS['--anchor']['help'] += ' Stages: "2_window" (may be omitted because of default).'
PARAMS['--anchor_vars']['help'] += ' Stages: "2_window" (may be omitted because of default).'
PARAMS['--label']['help'] += ' Stages: "4_info" (may be omitted because of default).'

def post_process_args(args):
  args = args.copy()
  args = filter.post_process_args(args)
  args = get_window.post_process_args(args)
  # Allow any prefix of a stage name.
  args['stages'] = [x for x in STAGES if x.startswith(tuple(args['stages']))]
  if (args['library_names'] is not None) and (args['input_list'] is not None):
    if len(args['library_names']) != len(args['input_list']):
      raise Exception('Number of NAMES must be the same as the number of INPUTs.')
    if not all(x == y for x, y in zip(args['library_names'], sorted(args['library_names']))):
      raise Exception('NAMES must be in alphabetical order.')
  if args['input_list'] is not None:
    f_names = [file_names.get_file_name(x) for x in args['input_list']]
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
      ' repeat of the experiment must be passed in alphabetical order.' +
      ' Certain parameters may be omitted if the corresponding stage is not run,' +
      ' the parameter has a default value, or the previous stage uses the same parameter' +
      ' and the "*_args.json" file for the previous stage is present in the directory' +
      ' (see the help for each parameter).'
    ),
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
  )
  group_multi = parser.add_argument_group('Multiple stages')
  group_align = parser.add_argument_group('Stage 0_align')
  group_filter = parser.add_argument_group('Stage 1_filter')
  group_window = parser.add_argument_group('Stage 2_window')
  group_info = parser.add_argument_group('Metadata')

  parser.add_argument('--stages', **PARAMS['--stages'])
  parser.add_argument('--no_align', **PARAMS['--no_align'])
  parser.add_argument('-o', **PARAMS['-o'])
  group_multi.add_argument('--ref', **PARAMS['--ref'])
  group_multi.add_argument('--names', **PARAMS['--names'])
  group_multi.add_argument('--dsb', **PARAMS['--dsb'])
  group_multi.add_argument('--window', **PARAMS['--window'])
  group_align.add_argument('-i', **PARAMS['-i'])
  group_align.add_argument('--bt2', **PARAMS['--bt2'])
  group_filter.add_argument('--reads', **PARAMS['--reads'])
  group_filter.add_argument('--min_len', **PARAMS['--min_len'])
  group_filter.add_argument('--max_sub', **PARAMS['--max_sub'])
  group_filter.add_argument('--rc', **PARAMS['--rc'])
  group_filter.add_argument('--consec', **PARAMS['--consec'])
  group_filter.add_argument('--touch', **PARAMS['--touch'])
  group_filter.add_argument('--quiet', **PARAMS['--quiet'])
  group_window.add_argument('--anchor', **PARAMS['--anchor'])
  group_window.add_argument('--anchor_vars', **PARAMS['--anchor_vars'])
  group_info.add_argument('--label', **PARAMS['--label'])

  return post_process_args(vars(parser.parse_args()))

def write_args(args, output, prefix):
  out_file = file_names.args_file(output, prefix)
  log_utils.log_output(out_file)
  file_utils.write_json(args, out_file)

def read_args(output, prefix):
  in_file = file_names.args_file(output, prefix)
  log_utils.log_input(in_file)
  return file_utils.read_json(in_file)

def do_0_align(
  output,
  input_list,
  library_names,
  ref_seq_file,
  bowtie2_args,
):
  if input_list is None:
    raise Exception('INPUT must be provided for stage "0_align".')
  if ref_seq_file is None:
    raise Exception('REF must be provided for stage "0_align".')
  if library_names is None:
    library_names = [file_names.get_file_name(x) for x in input_list]
  args = {
    'output': output,
    'input_list': input_list,
    'library_names': library_names,
    'ref_seq_file': ref_seq_file,
    'bowtie2_args': bowtie2_args,
  }
  bowtie2_args = ' '.join(bowtie2_args)
  bowtie2_index_file = file_names.bowtie2_index(output)
  # Existing Bowtie 2 index files must be overwritten.
  # For some reason, Bowtie 2 does not overwrite the index file if it already exists.
  file_utils.make_parent_dir(bowtie2_index_file, overwrite=True)
  log_utils.log_output('Bowtie2 index file: ' + bowtie2_index_file)
  if os.system(f'bowtie2-build-s {ref_seq_file} {bowtie2_index_file} --quiet') != 0:
    raise Exception('Bowtie 2 build failed.')

  for i in range(len(input_list)):
    sam_file = file_names.sam_file(output, library_names[i])
    file_utils.make_parent_dir(sam_file)
    input_ext = os.path.splitext(input_list[i])[1]
    if input_ext in ['.fastq', '.fq']:
      flags = '-q'
    elif input_ext in ['.fasta', '.fa', '.fna']:
      flags = '-f'
    else:
      flags = '-r'
    bowtie2_command = f'bowtie2-align-s --no-hd {flags} {bowtie2_args} -x {bowtie2_index_file} {input_list[i]} -S {sam_file} --quiet'
    log_utils.log('Bowtie 2 command: ' + bowtie2_command)
    if os.system(bowtie2_command) != 0:
      raise Exception('Bowtie 2 alignment failed.')
  write_args(args, output, 'align')

def do_1_filter(
  output,
  library_names,
  total_reads,
  ref_seq_file,
  dsb_pos,
  min_length,
  max_subst,
  reverse_complement,
  consecutive,
  dsb_touch,
  quiet,
):
  if ref_seq_file is None:
    raise Exception('REF must be provided for stage "1_filter".')
  if dsb_pos is None:
    raise Exception('DSB must be provided for stage "1_filter".')
  if min_length is None:
    raise Exception('MIN_LEN must be provided for stage "1_filter".')

  if library_names is None:
    # Infer names from the SAM files.
    library_names = sorted([
      file_names.get_file_name(x)
      for x in glob.glob(os.path.join(output, '*.sam'))
    ])
  input_list = [file_names.sam_file(output, x) for x in library_names]

  args = {
    'input_list': input_list,
    'output': [
      file_names.filter(output, 'accepted'),
      file_names.filter(output, 'rejected'),
    ],
    'debug_file': file_names.filter(output, 'debug'),
    'library_names': library_names,
    'total_reads': total_reads,
    'ref_seq_file': ref_seq_file,
    'dsb_pos': dsb_pos,
    'min_length': min_length,
    'max_subst': max_subst,
    'reverse_complement': reverse_complement,
    'consecutive': consecutive,
    'dsb_touch': dsb_touch,
    'quiet': quiet,
  }

  filter.main(**args)
  write_args(args, output, 'filter')

def do_2_window(
  output,
  ref_seq_file,
  dsb_pos,
  window_size,
  anchor_size,
  anchor_substs,
  anchor_indels,
):
  if ref_seq_file is None:
    raise Exception('REF must be provided for stage "2_window".')
  if dsb_pos is None:
    raise Exception('DSB must be provided for stage "2_window".')
  if window_size is None:
    raise Exception('WINDOW must be provided for stage "2_window".')
  if anchor_size is None:
    raise Exception('ANCHOR must be provided for stage "2_window".')
  if anchor_substs is None:
    raise Exception('ANCHOR_VARS must be provided for stage "2_window".')
  if anchor_indels is None:
    raise Exception('ANCHOR_VARS must be provided for stage "2_window".')

  for subst_type in constants.SUBST_TYPES:
    args = {
      'input': file_names.filter(output, 'accepted'),
      'output': file_names.window(output, subst_type),
      'ref_seq_file': ref_seq_file,
      'dsb_pos': dsb_pos,
      'window_size': window_size,
      'anchor_size': anchor_size,
      'anchor_substs': anchor_substs,
      'anchor_indels': anchor_indels,
      'subst_type': subst_type,
    }
    get_window.main(**args)
    write_args(args, output, 'window_' + subst_type)

def do_3_variation(output):
  for subst_type in constants.SUBST_TYPES:
    args = {
      'input': file_names.window(output, subst_type),
      'output': file_names.variation(output, subst_type),
    }
    get_variation.main(**args)
    write_args(args, output, 'variation_' + subst_type)

def do_4_info(
  output,
  ref_seq_file,
  dsb_pos,
  window_size,
  label,
):
  if output is None:
    raise Exception('OUTPUT must be provided for stage "4_info".')
  if ref_seq_file is None:
    raise Exception('REF must be provided for stage "4_info".')
  if dsb_pos is None:
    raise Exception('DSB must be provided for stage "4_info".')
  if window_size is None:
    raise Exception('WINDOW must be provided for stage "4_info".')
  name = file_names.get_file_name(output)
  if label is None:
    label = name
  ref_seq = file_utils.read_seq(ref_seq_file)
  data_info = common_utils.make_data_info(
    format = 'individual',
    names = [name],
    labels = [label],
    ref_seqs = [ref_seq],
    ref_seq_window = get_window.get_ref_seq_window(
      ref_seq = ref_seq,
      dsb_pos = dsb_pos,
      window_size = window_size,
    ),
  )
  out_file = file_names.data_info(output)
  log_utils.log_output(out_file)
  file_utils.write_json(data_info, out_file)

def main(
  stages,
  output,

  input_list,
  ref_seq_file,
  bowtie2_args,

  library_names,
  total_reads,
  dsb_pos,
  min_length,
  max_subst,
  reverse_complement,
  consecutive,
  dsb_touch,
  quiet,

  window_size,
  anchor_size,
  anchor_substs,
  anchor_indels,

  label,
):
  # Copy reference sequence file to output directory.
  if ref_seq_file is not None:
    shutil.copy(ref_seq_file, file_names.ref_seq_file(output))

  log_utils.blank_line()

  if '0_align' in stages:
    log_utils.log('Running stage "0_align".')
    prev_args = {
      'output': output,
      'input_list': input_list,
      'library_names': library_names,
      'ref_seq_file': ref_seq_file,
      'bowtie2_args': bowtie2_args,
    }
    do_0_align(**prev_args)
    log_utils.blank_line()

  if '1_filter' in stages:
    log_utils.log('Running stage "1_filter".')
    if os.path.exists(file_names.args_file(output, 'align')):
      prev_args = read_args(output, 'align')
      if ref_seq_file is None:
        ref_seq_file = prev_args.get('ref_seq_file')
      if library_names is None:
        library_names = prev_args.get('library_names')

    do_1_filter(
      output = output,
      library_names = library_names,
      total_reads = total_reads,
      ref_seq_file = ref_seq_file,
      dsb_pos = dsb_pos,
      min_length = min_length,
      max_subst = max_subst,
      reverse_complement = reverse_complement,
      consecutive = consecutive,
      dsb_touch = dsb_touch,
      quiet = quiet,
    )
    log_utils.blank_line()

  if '2_window' in stages:
    log_utils.log('Running stage "2_window".')
    prev_args = read_args(output, 'filter')
    if ref_seq_file is None:
      ref_seq_file = prev_args.get('ref_seq_file')
    if dsb_pos is None:
      dsb_pos = prev_args.get('dsb_pos')
    do_2_window(
      output = output,
      ref_seq_file = ref_seq_file,
      dsb_pos = dsb_pos,
      window_size = window_size,
      anchor_size = anchor_size,
      anchor_substs = anchor_substs,
      anchor_indels = anchor_indels,
    )
    log_utils.blank_line()

  if '3_variation' in stages:
    log_utils.log('Running stage "3_variation".')
    do_3_variation(output = output)
    log_utils.blank_line()

  if '4_info' in stages:
    log_utils.log('Running stage "4_info".')
    prev_args = read_args(output, 'window_' + constants.SUBST_TYPES[0])
    if ref_seq_file is None:
      ref_seq_file = prev_args.get('ref_seq_file')
    if dsb_pos is None:
      dsb_pos = prev_args.get('dsb_pos')
    if window_size is None:
      window_size = prev_args.get('window_size')
    do_4_info(
      output = output,
      ref_seq_file = ref_seq_file,
      dsb_pos = dsb_pos,
      window_size = window_size,
      label = label,
    )
    log_utils.blank_line()

# This allows the "DSBplot-preprocess" command to be run from the command line.
def entry_point():
  main(**parse_args())

if __name__ == '__main__':
  entry_point()
