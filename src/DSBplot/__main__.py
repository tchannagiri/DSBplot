import argparse
import sys

import preprocess
import comparison
import graph
import histogram
import pptx

if __name__ == '__main__':
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
    preprocess.main(**preprocess.parse_args())
  elif args.command == 'comparison':
    comparison.main(**comparison.parse_args())
  elif args.command == 'graph':
    graph.main(**graph.parse_args())
  elif args.command == 'histogram':
    histogram.main(**histogram.parse_args())
  elif args.command == 'pptx':
    pptx.main(**pptx.parse_args())
  else:
    raise Exception('Error: Unrecognized command.' + args.command)

