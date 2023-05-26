import os
import unittest
import glob

import DSBplot.utils.file_utils as file_utils

class TestReadSeq(unittest.TestCase):
  def test_read_seq(self):
    """
    Test the file_utils.read_seq() function.
    """
    for file in glob.glob(
      os.path.join(
        os.path.dirname(__file__),
        'data',
        'input',
        'test_seq',
        '*'
      )
    ):
      file_name = os.path.basename(file)
      if file_name == 'seq_1.fa':
        self.assertEqual(file_utils.read_seq(file), 'ACTG', msg=file)
      elif file_name == 'seq_2.fa':
        self.assertEqual(file_utils.read_seq(file), 'ACTG', msg=file)
      elif file_name == 'seq_3.fa':
        with self.assertRaises(Exception, msg=file):
          file_utils.read_seq(file)
      elif file_name == 'seq_4.fa':
        self.assertEqual(file_utils.read_seq(file), 'ACTGACTG', msg=file)
      elif file_name == 'seq_5.fa':
        with self.assertRaises(Exception, msg=file):
          file_utils.read_seq(file)
      elif file_name == 'seq_6.fa':
        with self.assertRaises(Exception, msg=file):
          file_utils.read_seq(file)
      elif file_name == 'seq_7.fa':
        with self.assertRaises(Exception, msg=file):
          file_utils.read_seq(file)
      elif file_name == 'seq_8.fa':
        with self.assertRaises(Exception, msg=file):
          file_utils.read_seq(file)
      elif file_name == 'seq_9.fa':
        with self.assertRaises(Exception, msg=file):
          file_utils.read_seq(file)
      elif file_name == 'seq_1.txt':
        self.assertEqual(file_utils.read_seq(file), 'ACTG', msg=file)
      elif file_name == 'seq_2.txt':
        self.assertEqual(file_utils.read_seq(file), 'ACTG', msg=file)
      elif file_name == 'seq_3.txt':
        self.assertEqual(file_utils.read_seq(file), 'ACTGACTG', msg=file)
      elif file_name == 'seq_4.txt':
        self.assertEqual(file_utils.read_seq(file), 'ACTGACTG', msg=file)
      elif file_name == 'seq_5.txt':
        with self.assertRaises(Exception, msg=file):
          file_utils.read_seq(file)
      elif file_name == 'seq_6.txt':
        with self.assertRaises(Exception, msg=file):
          file_utils.read_seq(file)
      elif file_name == 'seq_7.txt':
        with self.assertRaises(Exception, msg=file):
          file_utils.read_seq(file)
      elif file_name == 'seq_8.txt':
        with self.assertRaises(Exception, msg=file):
          file_utils.read_seq(file)
      elif file_name == 'seq_9.txt':
        with self.assertRaises(Exception, msg=file):
          file_utils.read_seq(file)
      else:
        raise Exception('Unexpected file: ' + file_name)

if __name__ == '__main__':
  unittest.main()
