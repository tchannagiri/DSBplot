import DSBplot.get_pptx.get_pptx as get_pptx

def parse_args():
  return get_pptx.parse_args()

def main(**args):
  get_pptx.main(**args)

if __name__ == '__main__':
  main(**parse_args())