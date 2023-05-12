import DSBplot.get_pptx.get_pptx as get_pptx

def parse_args():
  return get_pptx.parse_args()

def main(**args):
  get_pptx.main(**args)

# This allows the "DSBplot-pptx" command to be run from the command line.
def entry_point():
  main(**parse_args())

if __name__ == '__main__':
  entry_point()