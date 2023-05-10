import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'utils'))) # allow importing
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '7_get_pptx'))) # allow importing

import get_pptx

if __name__ == '__main__':
  get_pptx.main(**get_pptx.parse_args())