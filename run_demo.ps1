### Preprocessing ###
python preprocess.py --input data_demo/input/fastq/db1_R1.fq data_demo/input/fastq/db2_R1.fq data_demo/input/fastq/db3_R1.fq data_demo/input/fastq/db4_R1.fq --ref_seq_file data_demo/input/ref_seq/2DSB_R1_branch.fa --dsb_pos 67 --output data_demo/output/db_R1 --label db_R1 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_demo/input/fastq/db1_R2.fq data_demo/input/fastq/db2_R2.fq data_demo/input/fastq/db3_R2.fq data_demo/input/fastq/db4_R2.fq --ref_seq_file data_demo/input/ref_seq/2DSB_R2_branch.fa --dsb_pos 46 --output data_demo/output/db_R2 --label db_R2 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_demo/input/fastq/sense1_R1.fq data_demo/input/fastq/sense2_R1.fq data_demo/input/fastq/sense3_R1.fq data_demo/input/fastq/sense4_R1.fq --ref_seq_file data_demo/input/ref_seq/2DSB_R1_sense.fa --dsb_pos 67 --output data_demo/output/sense_R1 --label sense_R1 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_demo/input/fastq/sense1_R2.fq data_demo/input/fastq/sense2_R2.fq data_demo/input/fastq/sense3_R2.fq data_demo/input/fastq/sense4_R2.fq --ref_seq_file data_demo/input/ref_seq/2DSB_R2_sense.fa --dsb_pos 46 --output data_demo/output/sense_R2 --label sense_R2 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_demo/input/fastq/dcmv1_R1.fq data_demo/input/fastq/dcmv2_R1.fq data_demo/input/fastq/dcmv3_R1.fq data_demo/input/fastq/dcmv4_R1.fq --ref_seq_file data_demo/input/ref_seq/2DSB_R1_cmv.fa --dsb_pos 67 --output data_demo/output/dcmv_R1 --label dcmv_R1 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_demo/input/fastq/dcmv1_R2.fq data_demo/input/fastq/dcmv2_R2.fq data_demo/input/fastq/dcmv3_R2.fq data_demo/input/fastq/dcmv4_R2.fq --ref_seq_file data_demo/input/ref_seq/2DSB_R2_cmv.fa --dsb_pos 46 --output data_demo/output/dcmv_R2 --label dcmv_R2 --total_reads 3000 3000 3000 3000

# Example of running preprocessing stages separately for one experiment
python preprocess.py --input data_demo/input/fastq/db1_R1.fq data_demo/input/fastq/db2_R1.fq data_demo/input/fastq/db3_R1.fq data_demo/input/fastq/db4_R1.fq --ref_seq_file data_demo/input/ref_seq/2DSB_R1_branch.fa --output data_demo/output/db_R1 --stages 0_align
python preprocess.py --ref_seq_file data_demo/input/ref_seq/2DSB_R1_branch.fa --dsb_pos 67 --output data_demo/output/db_R1 --stages 1_filter
python preprocess.py --output data_demo/output/db_R1 --stages 2_combine
python preprocess.py --ref_seq_file data_demo/input/ref_seq/2DSB_R1_branch.fa --dsb_pos 67 --output data_demo/output/db_R1 --label db_R1 --total_reads 3000 3000 3000 3000 --stages 3_window
python preprocess.py --output data_demo/output/db_R1 --stages 4_graph
python preprocess.py --output data_demo/output/db_R1 --stages 5_histogram

### Comparison data ###
python comparison.py --input data_demo/output/sense_R1 data_demo/output/db_R1 --output data_demo/output/sense_db_R1
python comparison.py --input data_demo/output/sense_R2 data_demo/output/db_R2 --output data_demo/output/sense_db_R2
python comparison.py --input data_demo/output/sense_R1 data_demo/output/dcmv_R1 --output data_demo/output/sense_dcmv_R1
python comparison.py --input data_demo/output/sense_R2 data_demo/output/dcmv_R2 --output data_demo/output/sense_dcmv_R2

### Plot the graphs ###
python graph.py --input data_demo/output/db_R1 --output plot_demo/graph/universal/db_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5  --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/db_R2 --output plot_demo/graph/universal/db_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_R1 --output plot_demo/graph/universal/sense_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_R2 --output plot_demo/graph/universal/sense_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/dcmv_R1 --output plot_demo/graph/universal/dcmv_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/dcmv_R2 --output plot_demo/graph/universal/dcmv_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot the comparison graphs ###
python graph.py --input data_demo/output/sense_db_R1 --node_comparison_colors "#cf191b" "#33a02c" --output plot_demo/graph/universal/sense_db_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_db_R2 --node_comparison_colors "#cf191b" "#33a02c" --output plot_demo/graph/universal/sense_db_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_dcmv_R1 --output plot_demo/graph/universal/sense_dcmv_R1.png --node_comparison_colors "#cf191b" "#ffe669" --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_dcmv_R2 --output plot_demo/graph/universal/sense_dcmv_R2.png --node_comparison_colors "#cf191b" "#ffe669" --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot the 3D histograms ###
python histogram.py --input data_demo/output/db_R1 --output plot_demo/histogram/db_R1_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_demo/output/db_R1 --output plot_demo/histogram/db_R1_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_demo/output/db_R1 --output plot_demo/histogram/db_R1_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_demo/output/db_R2 --output plot_demo/histogram/db_R2_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_demo/output/db_R2 --output plot_demo/histogram/db_R2_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_demo/output/db_R2 --output plot_demo/histogram/db_R2_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_demo/output/sense_R1 --output plot_demo/histogram/sense_R1_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_demo/output/sense_R1 --output plot_demo/histogram/sense_R1_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_demo/output/sense_R1 --output plot_demo/histogram/sense_R1_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_demo/output/sense_R2 --output plot_demo/histogram/sense_R2_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_demo/output/sense_R2 --output plot_demo/histogram/sense_R2_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_demo/output/sense_R2 --output plot_demo/histogram/sense_R2_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_demo/output/dcmv_R1 --output plot_demo/histogram/dcmv_R1_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_demo/output/dcmv_R1 --output plot_demo/histogram/dcmv_R1_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_demo/output/dcmv_R1 --output plot_demo/histogram/dcmv_R1_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_demo/output/dcmv_R2 --output plot_demo/histogram/dcmv_R2_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_demo/output/dcmv_R2 --output plot_demo/histogram/dcmv_R2_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_demo/output/dcmv_R2 --output plot_demo/histogram/dcmv_R2_deletion.png --color "#8080ff" --variation_type deletion --label_type relative

### Plot graphs with common Kamada layout ###
python graph.py `
  --input data_demo/output/db_R1 data_demo/output/db_R2 data_demo/output/sense_R1 data_demo/output/sense_R2 data_demo/output/dcmv_R1 data_demo/output/dcmv_R2 data_demo/output/sense_db_R1 data_demo/output/sense_db_R2 data_demo/output/sense_dcmv_R1 data_demo/output/sense_dcmv_R2 `
  --output plot_demo/graph/kamada_common/db_R1.png plot_demo/graph/kamada_common/db_R2.png plot_demo/graph/kamada_common/sense_R1.png plot_demo/graph/kamada_common/sense_R2.png plot_demo/graph/kamada_common/dcmv_R1.png plot_demo/graph/kamada_common/dcmv_R2.png plot_demo/graph/kamada_common/sense_db_R1.png plot_demo/graph/kamada_common/sense_db_R2.png plot_demo/graph/kamada_common/sense_dcmv_R1.png plot_demo/graph/kamada_common/sense_dcmv_R2.png `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout kamada_layout --width 2400 --height 1800

### Plot graphs with common radial layout ###
python graph.py `
  --input data_demo/output/db_R1 data_demo/output/db_R2 data_demo/output/sense_R1 data_demo/output/sense_R2 data_demo/output/dcmv_R1 data_demo/output/dcmv_R2 data_demo/output/sense_db_R1 data_demo/output/sense_db_R2 data_demo/output/sense_dcmv_R1 data_demo/output/sense_dcmv_R2 `
  --output plot_demo/graph/radial_common/db_R1.png plot_demo/graph/radial_common/db_R2.png plot_demo/graph/radial_common/sense_R1.png plot_demo/graph/radial_common/sense_R2.png plot_demo/graph/radial_common/dcmv_R1.png plot_demo/graph/radial_common/dcmv_R2.png plot_demo/graph/radial_common/sense_db_R1.png plot_demo/graph/radial_common/sense_db_R2.png plot_demo/graph/radial_common/sense_dcmv_R1.png plot_demo/graph/radial_common/sense_dcmv_R2.png `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --range_x -20 20 --range_y -20 20 --layout radial_layout --width 2400 --height 1800

### Plot graphs with common fractal layout ###
python graph.py `
  --input data_demo/output/db_R1 data_demo/output/db_R2 data_demo/output/sense_R1 data_demo/output/sense_R2 data_demo/output/dcmv_R1 data_demo/output/dcmv_R2 data_demo/output/sense_db_R1 data_demo/output/sense_db_R2 data_demo/output/sense_dcmv_R1 data_demo/output/sense_dcmv_R2 `
  --output plot_demo/graph/fractal_common/db_R1.png plot_demo/graph/fractal_common/db_R2.png plot_demo/graph/fractal_common/sense_R1.png plot_demo/graph/fractal_common/sense_R2.png plot_demo/graph/fractal_common/dcmv_R1.png plot_demo/graph/fractal_common/dcmv_R2.png plot_demo/graph/fractal_common/sense_db_R1.png plot_demo/graph/fractal_common/sense_db_R2.png plot_demo/graph/fractal_common/sense_dcmv_R1.png plot_demo/graph/fractal_common/sense_dcmv_R2.png `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout fractal_layout --width 2400 --height 1800
