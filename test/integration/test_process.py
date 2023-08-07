import os
import shutil
import unittest

class TestProcess(unittest.TestCase):
  def get_check_files():
    """
    The files that should be checked for equality between the
    new and old outputs.
    """
    return [
      'data_info.json',
      'ref_seq.fasta',
      'filter_debug.csv',
      'window_withoutSubst.csv',
      'window_withSubst.csv',
      'variation_withoutSubst.csv',
      'variation_withSubst.csv',
    ]
  
  def check_equality(self, file_1, file_2):
    """
    Check if the two files are equal.
    """
    with open(file_1, 'r') as f_1:
      lines_1 = f_1.readlines()
    with open(file_2, 'r') as f_2:
      lines_2 = f_2.readlines()
    self.assertEqual(
      len(lines_1),
      len(lines_2),
      f'The files have different number of lines.\n{file_1}\n{file_2}'
    )
    for i in range(len(lines_1)):
      self.assertEqual(
        lines_1[i],
        lines_2[i],
        f'The files differ at line {i + 1}.\n{file_1}\n{file_2}'
      )

  def test_process(self):
    """
    Test the filter with the given combination of file types for the reads
    and the reference sequence.
    """
    input = [
      os.path.join(os.path.dirname(__file__), 'input', 'fastq', 'Sense_R1_1.fq'),
      os.path.join(os.path.dirname(__file__), 'input', 'fastq', 'Sense_R1_2.fq'),
      os.path.join(os.path.dirname(__file__), 'input', 'fastq', 'Sense_R1_3.fq'),
      os.path.join(os.path.dirname(__file__), 'input', 'fastq', 'Sense_R1_4.fq'),
    ]
    ref = os.path.join(os.path.dirname(__file__), 'input', 'ref_seq', '2DSB_Sense_R1.fa')
    output = os.path.join(os.path.dirname(__file__), 'output', 'Sense_R1')
    output_expected = os.path.join(os.path.dirname(__file__), 'output_expected', 'Sense_R1')
    reads = [3000, 3000, 3000, 3000]

    # Clear output directory (if exists)
    if os.path.exists(output):
      shutil.rmtree(output)
    
    self.assertEqual(
      os.system(
        'DSBplot-process -o {} -i {} {} {} {} --ref {} --dsb 67 --label "Sense (R1)" --reads {} {} {} {}'
        .format(output, *input, ref, *reads)
      ),
      0,
      'DSBplot-process failed.'
    )

    for file in TestProcess.get_check_files():
      print('Testing file: ' + file)
      self.check_equality(os.path.join(output, file), os.path.join(output_expected, file))

if __name__ == '__main__':
  unittest.main()
