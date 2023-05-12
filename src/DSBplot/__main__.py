import argparse
import sys

import DSBplot.preprocess as preprocess
import DSBplot.comparison as comparison
import DSBplot.graph as graph
import DSBplot.histogram as histogram
import DSBplot.pptx1 as pptx1

# This allows the "DSBplot" command to be run from the command line.
def entry_point():
  parser = argparse.ArgumentParser(
    description = (
      'DSBplot main command.' +
      'Use one of the options commands "preprocess", "comparison", "graph", "histogram", or "pptx".' +
      ' Follow by "--help" for more information on individual commands.' +
      ' E.g., "python -m DSBplot preprocess --help".'
    )
  )
  parser.add_argument(
    'command',
    choices = ['preprocess', 'comparison', 'graph', 'histogram', 'pptx'],
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
  if args.command == 'preprocess':
    sys.argv[0] = preprocess.__file__
    preprocess.entry_point()
  elif args.command == 'comparison':
    sys.argv[0] = comparison.__file__
    comparison.entry_point()
  elif args.command == 'graph':
    sys.argv[0] = graph.__file__
    graph.entry_point()
  elif args.command == 'histogram':
    sys.argv[0] = histogram.__file__
    histogram.entry_point()
  elif args.command == 'pptx':
    sys.argv[0] = pptx1.__file__
    pptx1.entry_point()
  else:
    raise Exception('Error: Unrecognized command.' + args.command)

if __name__ == '__main__':
  entry_point()