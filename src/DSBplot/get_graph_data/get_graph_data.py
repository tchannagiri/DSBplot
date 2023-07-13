# FIXME: DELETE

import itertools
import pandas as pd
import argparse

import DSBplot.utils.file_names as file_names
import DSBplot.utils.file_utils as file_utils
import DSBplot.utils.alignment_utils as alignment_utils
import DSBplot.utils.log_utils as log_utils
import DSBplot.utils.common_utils as common_utils
import DSBplot.utils.constants as constants
import DSBplot.utils.graph_utils as graph_utils

# FIXME: This function is not used anywhere. DELETE.
# def parse_args():
#   parser = argparse.ArgumentParser(
#     description = (
#       'Precompute the necessary data for the plotting graphs.' +
#       ' Uses as input the output of the "get_window_data.py" script.'
#     )
#   )
#   parser.add_argument(
#     '--input',
#     type = common_utils.check_dir,
#     help = 'Directory where output from "get_window_data" stage is located.',
#     required = True,
#   )
#   parser.add_argument(
#     '--output',
#     type = common_utils.check_dir_output,
#     help = 'Output directory.',
#     required = True,
#   )
#   parser.add_argument(
#     '--subst_type',
#     type = str,
#     choices = constants.SUBST_TYPES,
#     help = 'Whether to process the files with/without substitutions.',
#   )
#   return vars(parser.parse_args())


# FIXME: This function is not used anywhere. DELETE.
# def write_sequence_data(input_dir, output_dir, subst_type):
#   """
#     Make the main node data and write it to a file.
#   """

#   data = file_utils.read_tsv(
#     file_names.window(
#       input_dir,
#       constants.FREQ_FILTER_MEAN,
#       subst_type,
#     )
#   )
#   data_info = file_utils.read_tsv_dict(file_names.data_info(output_dir))
#   data = get_sequence_data(data, data_info['format'])
#   out_file_name = file_names.sequence_data(output_dir, subst_type)
#   file_utils.write_tsv(data, out_file_name)
#   log_utils.log_output(out_file_name)


# FIXME: This function is not used anywhere. DELETE.
# def write_edge_data(output_dir, subst_type):
#   """
#     Make adjacency edge data and write to file.
#     Sequence data should have been created already.
#   """
#   in_file_name = file_names.sequence_data(output_dir, subst_type)
#   out_file_name = file_names.edge_data(output_dir, subst_type)

#   sequence_data = file_utils.read_tsv(in_file_name)
#   edge_data = get_edge_data(sequence_data)
#   file_utils.write_tsv(edge_data, out_file_name)
#   log_utils.log_output(out_file_name)

# FIXME: This function is not used anywhere. DELETE.
# def write_graph_stats(output_dir, subst_type):
#   """
#     Get graph summary statistics and write to file.
#     Sequence data and edge data should have been created already.
#   """

#   graph = graph_utils.load_graph(output_dir, subst_type)
#   data_info = file_utils.read_tsv_dict(file_names.data_info(output_dir))
#   graph_stats = graph_utils.get_graph_stats_ref_component(data_info['format'], graph)
#   graph_stats = pd.DataFrame.from_records([graph_stats])
#   out_file_name = file_names.graph_stats(output_dir, subst_type)
#   file_utils.write_tsv(graph_stats, out_file_name)
#   log_utils.log_output(out_file_name)

# FIXME: No longer a runnable function. DELETE.
# def main(
#   input,
#   output,
#   subst_type,
# ):
#   log_utils.log_input(input)

#   # copy data info
#   if input != output:
#     input_data_info_file = file_names.data_info(input)
#     output_data_info_file = file_names.data_info(output)
#     file_utils.copy(input_data_info_file, output_data_info_file)

#   write_sequence_data(input, output, subst_type)
#   write_edge_data(output, subst_type)
#   write_graph_stats(output, subst_type)

#   log_utils.blank_line()

# if __name__ == '__main__':
#   main(**parse_args())
