# Note: This is a heavily outdated test and will no longer work.

import os
import unittest

import DSBplot.preprocessing as preprocessing

class TestFilterNhej(unittest.TestCase):
  def get_check_files():
    """
    The files that should be checked for equality between the
    new and old outputs.
    """
    return [
      '1_filter_nhej/1.tsv',
      '1_filter_nhej/debug_1.txt',
      '2_combine_repeat/out.tsv',
      '3_window/data_info.tsv',
      '3_window/window_count_withoutSubst.tsv',
      '3_window/window_count_withSubst.tsv',
      '3_window/window_freq_filter_mean_withoutSubst.tsv',
      '3_window/window_freq_filter_mean_withSubst.tsv',
      '3_window/window_freq_filter_withoutSubst.tsv',
      '3_window/window_freq_filter_withSubst.tsv',
      '3_window/window_freq_withoutSubst.tsv',
      '3_window/window_freq_withSubst.tsv',
      '4_graph/data_info.tsv',
      '4_graph/edge_data_withoutSubst.tsv',
      '4_graph/edge_data_withSubst.tsv',
      '4_graph/graph_stats_withoutSubst.tsv',
      '4_graph/graph_stats_withSubst.tsv',
      '4_graph/sequence_data_withoutSubst.tsv',
      '4_graph/sequence_data_withSubst.tsv',
      '5_histogram/data_info.tsv',
      '5_histogram/variation_grouped_withoutSubst.tsv',
      '5_histogram/variation_grouped_withSubst.tsv',
      '5_histogram/variation_withoutSubst.tsv',
      '5_histogram/variation_withSubst.tsv',
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

  def do_test_filter(self, ext_reads, ext_ref_seq):
    """
    Test the NHEJ filter with the given combination of file types for the reads
    and the reference sequence.
    """
    file_reads = os.path.join(
      os.path.dirname(__file__),
      'data/input/reads/',
      'Sense_R1.' + ext_reads
    )
    file_ref_seq = os.path.join(
      os.path.dirname(__file__),
      'data/input/ref_seq/',
      '2DSB_Sense_R1.' + ext_ref_seq
    )
    dir_output = os.path.join(
      os.path.dirname(__file__),
      'data/output/',
      'Sense_R1_' + ext_reads
    )
    dir_output_expected = os.path.join(
      os.path.dirname(__file__),
      'data/output_expected/',
      'Sense_R1_' + ext_reads
    )
    preprocessing.main(
      input = [file_reads],
      output = dir_output,
      ref_seq_file = file_ref_seq,
      dsb_pos = 67,
      min_length = None,
      window_size = 10,
      anchor_size = 20,
      anchor_mismatches = 1,
      total_reads = [3000],
      freq_min = 1e5,
      label = 'Sense',
      quiet = False,
    )
    for file in TestFilterNhej.get_check_files():
      self.check_equality(
        os.path.join(dir_output, file),
        os.path.join(dir_output_expected, file)
      )

  def test_filter(self):
    for ext_read in ['fq', 'fasta', 'txt']:
      for ext_ref_seq in ['fa', 'txt']:
        self.do_test_filter(ext_read, ext_ref_seq)

if __name__ == '__main__':
  unittest.main()
