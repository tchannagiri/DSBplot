### Preprocessing ###
python preprocess.py --input data_input/fastq/db1_R1.fq data_input/fastq/db2_R1.fq data_input/fastq/db3_R1.fq data_input/fastq/db4_R1.fq --ref_seq_file data_input/ref_seq/2DSB_R1_branch.fa --dsb_pos 67 --output data_output/db_R1 --label db_R1 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_input/fastq/db1_R2.fq data_input/fastq/db2_R2.fq data_input/fastq/db3_R2.fq data_input/fastq/db4_R2.fq --ref_seq_file data_input/ref_seq/2DSB_R2_branch.fa --dsb_pos 46 --output data_output/db_R2 --label db_R2 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_input/fastq/sense1_R1.fq data_input/fastq/sense2_R1.fq data_input/fastq/sense3_R1.fq data_input/fastq/sense4_R1.fq --ref_seq_file data_input/ref_seq/2DSB_R1_sense.fa --dsb_pos 67 --output data_output/sense_R1 --label sense_R1 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_input/fastq/sense1_R2.fq data_input/fastq/sense2_R2.fq data_input/fastq/sense3_R2.fq data_input/fastq/sense4_R2.fq --ref_seq_file data_input/ref_seq/2DSB_R2_sense.fa --dsb_pos 46 --output data_output/sense_R2 --label sense_R2 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_input/fastq/dcmv1_R1.fq data_input/fastq/dcmv2_R1.fq data_input/fastq/dcmv3_R1.fq data_input/fastq/dcmv4_R1.fq --ref_seq_file data_input/ref_seq/2DSB_R1_cmv.fa --dsb_pos 67 --output data_output/dcmv_R1 --label dcmv_R1 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_input/fastq/dcmv1_R2.fq data_input/fastq/dcmv2_R2.fq data_input/fastq/dcmv3_R2.fq data_input/fastq/dcmv4_R2.fq --ref_seq_file data_input/ref_seq/2DSB_R2_cmv.fa --dsb_pos 46 --output data_output/dcmv_R2 --label dcmv_R2 --total_reads 3000 3000 3000 3000

### Comparison data ###
python comparison.py --input data_output/sense_R1 data_output/db_R1 --output data_output/sense_db_R1
python comparison.py --input data_output/sense_R2 data_output/db_R2 --output data_output/sense_db_R2
python comparison.py --input data_output/sense_R1 data_output/dcmv_R1 --output data_output/sense_dcmv_R1
python comparison.py --input data_output/sense_R2 data_output/dcmv_R2 --output data_output/sense_dcmv_R2

### Plot the graphs ###
python graph.py --input data_output/db_R1 --output plot/graph/universal/db_R1.png --legend --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5  --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/db_R2 --output plot/graph/universal/db_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/sense_R1 --output plot/graph/universal/sense_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/sense_R2 --output plot/graph/universal/sense_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/dcmv_R1 --output plot/graph/universal/dcmv_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/dcmv_R2 --output plot/graph/universal/dcmv_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### PLot the comparison graphs ###
python graph.py --input data_output/sense_db_R1 --node_comparison_colors "#cf191b" "#33a02c" --output plot/graph/universal/sense_db_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/sense_db_R2 --node_comparison_colors "#cf191b" "#33a02c" --output plot/graph/universal/sense_db_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/sense_dcmv_R1 --output plot/graph/universal/sense_dcmv_R1.png --node_comparison_colors "#CF191B" "#FFE669" --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/sense_dcmv_R2 --output plot/graph/universal/sense_dcmv_R2.png --node_comparison_colors "#CF191B" "#FFE669" --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot the 3D histograms ###
python histogram.py --input data_output/db_R1 --output plot/histogram/db_R1_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_output/db_R1 --output plot/histogram/db_R1_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_output/db_R1 --output plot/histogram/db_R1_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_output/db_R2 --output plot/histogram/db_R2_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_output/db_R2 --output plot/histogram/db_R2_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_output/db_R2 --output plot/histogram/db_R2_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_output/sense_R1 --output plot/histogram/sense_R1_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_output/sense_R1 --output plot/histogram/sense_R1_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_output/sense_R1 --output plot/histogram/sense_R1_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_output/sense_R2 --output plot/histogram/sense_R2_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_output/sense_R2 --output plot/histogram/sense_R2_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_output/sense_R2 --output plot/histogram/sense_R2_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_output/dcmv_R1 --output plot/histogram/dcmv_R1_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_output/dcmv_R1 --output plot/histogram/dcmv_R1_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_output/dcmv_R1 --output plot/histogram/dcmv_R1_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_output/dcmv_R2 --output plot/histogram/dcmv_R2_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_output/dcmv_R2 --output plot/histogram/dcmv_R2_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_output/dcmv_R2 --output plot/histogram/dcmv_R2_deletion.png --color "#8080ff" --variation_type deletion --label_type relative

### Plot graphs with common Kamada layout ###
python graph.py `
  --input data_output/db_R1 data_output/db_R2 data_output/sense_R1 data_output/sense_R2 data_output/dcmv_R1 data_output/dcmv_R2 data_output/sense_db_R1 data_output/sense_db_R2 data_output/sense_dcmv_R1 data_output/sense_dcmv_R2 `
  --output plot/graph/common_kamada/db_R1.png plot/graph/kamada_common/db_R2.png plot/graph/common_kamada/sense_R1.png plot/graph/common_kamada/sense_R2.png plot/graph/common_kamada/dcmv_R1.png plot/graph/common_kamada/dcmv_R2.png plot/graph/common_kamada/sense_db_R1.png plot/graph/common_kamada/sense_db_R2.png plot/graph/common_kamada/sense_dcmv_R1.png plot/graph/common_kamada/sense_dcmv_R2.png `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout kamada_layout --width 2400 --height 1800

### Plot graphs with common radial layout ###
python graph.py `
  --input data_output/db_R1 data_output/db_R2 data_output/sense_R1 data_output/sense_R2 data_output/dcmv_R1 data_output/dcmv_R2 data_output/sense_db_R1 data_output/sense_db_R2 data_output/sense_dcmv_R1 data_output/sense_dcmv_R2 `
  --output plot/graph/common_radial/db_R1.png plot/graph/radial_common/db_R2.png plot/graph/common_radial/sense_R1.png plot/graph/common_radial/sense_R2.png plot/graph/common_radial/dcmv_R1.png plot/graph/common_radial/dcmv_R2.png plot/graph/common_radial/sense_db_R1.png plot/graph/common_radial/sense_db_R2.png plot/graph/common_radial/sense_dcmv_R1.png plot/graph/common_radial/sense_dcmv_R2.png `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --range_x -15 15 --range_y -15 15 --layout radial_layout --width 2400 --height 1800
