import argparse
import sys

import DSBplot.process as process
import DSBplot.graph as graph
import DSBplot.histogram as histogram
import DSBplot.concat as concat

# This allows the "DSBplot" command to be run from the command line.
def entry_point():
  parser = argparse.ArgumentParser(
    description = (
      'DSBplot main command.' +
      'Use one of the options commands "process", "graph", "histogram", or "concat".' +
      ' Follow by "--help" for more information on individual commands.' +
      ' E.g., "python -m DSBplot process --help".'
    )
  )
  parser.add_argument(
    'command',
    choices = ['process', 'graph', 'histogram', 'concat'],
    help = 'Command to run.',
  )
  parser.add_argument(
    'args',
    nargs = argparse.REMAINDER,
    help = 'Arguments for the command.',
  )
  args = parser.parse_args()
  if sys.argv[1] != args.command:
    raise Exception('Error: "sys.argv[1]" != "args.command".')
  sys.argv.pop(1)
  if args.command == 'process':
    sys.argv[0] = process.__file__
    process.entry_point()
  elif args.command == 'graph':
    sys.argv[0] = graph.__file__
    graph.entry_point()
  elif args.command == 'histogram':
    sys.argv[0] = histogram.__file__
    histogram.entry_point()
  elif args.command == 'concat':
    sys.argv[0] = concat.__file__
    concat.entry_point()
  else:
    raise Exception('Error: Unrecognized command.' + args.command)

if __name__ == '__main__':
  entry_point()