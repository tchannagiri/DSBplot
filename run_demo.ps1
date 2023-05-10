### Preprocessing ###
python preprocess.py --input data_demo/input/fastq/db1_R1.fq data_demo/input/fastq/db2_R1.fq data_demo/input/fastq/db3_R1.fq data_demo/input/fastq/db4_R1.fq --ref_seq_file data_demo/input/ref_seq/2DSB_R1_branch.fa --dsb_pos 67 --output data_demo/output/db_R1 --label "BranchΔ" --total_reads 3000 3000 3000 3000
python preprocess.py --input data_demo/input/fastq/db1_R2.fq data_demo/input/fastq/db2_R2.fq data_demo/input/fastq/db3_R2.fq data_demo/input/fastq/db4_R2.fq --ref_seq_file data_demo/input/ref_seq/2DSB_R2_branch.fa --dsb_pos 46 --output data_demo/output/db_R2 --label "BranchΔ" --total_reads 3000 3000 3000 3000
python preprocess.py --input data_demo/input/fastq/sense1_R1.fq data_demo/input/fastq/sense2_R1.fq data_demo/input/fastq/sense3_R1.fq data_demo/input/fastq/sense4_R1.fq --ref_seq_file data_demo/input/ref_seq/2DSB_R1_sense.fa --dsb_pos 67 --output data_demo/output/sense_R1 --label "Sense" --total_reads 3000 3000 3000 3000
python preprocess.py --input data_demo/input/fastq/sense1_R2.fq data_demo/input/fastq/sense2_R2.fq data_demo/input/fastq/sense3_R2.fq data_demo/input/fastq/sense4_R2.fq --ref_seq_file data_demo/input/ref_seq/2DSB_R2_sense.fa --dsb_pos 46 --output data_demo/output/sense_R2 --label "Sense" --total_reads 3000 3000 3000 3000
python preprocess.py --input data_demo/input/fastq/dcmv1_R1.fq data_demo/input/fastq/dcmv2_R1.fq data_demo/input/fastq/dcmv3_R1.fq data_demo/input/fastq/dcmv4_R1.fq --ref_seq_file data_demo/input/ref_seq/2DSB_R1_cmv.fa --dsb_pos 67 --output data_demo/output/dcmv_R1 --label "pCMVΔ" --total_reads 3000 3000 3000 3000
python preprocess.py --input data_demo/input/fastq/dcmv1_R2.fq data_demo/input/fastq/dcmv2_R2.fq data_demo/input/fastq/dcmv3_R2.fq data_demo/input/fastq/dcmv4_R2.fq --ref_seq_file data_demo/input/ref_seq/2DSB_R2_cmv.fa --dsb_pos 46 --output data_demo/output/dcmv_R2 --label "pCMVΔ" --total_reads 3000 3000 3000 3000

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

### Plot the graphs (PNG) ###
python graph.py --input data_demo/output/db_R1 --output plot_demo/graph/png/universal/db_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5  --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/db_R2 --output plot_demo/graph/png/universal/db_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_R1 --output plot_demo/graph/png/universal/sense_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_R2 --output plot_demo/graph/png/universal/sense_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/dcmv_R1 --output plot_demo/graph/png/universal/dcmv_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/dcmv_R2 --output plot_demo/graph/png/universal/dcmv_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot the comparison graphs (PNG) ###
python graph.py --input data_demo/output/sense_db_R1 --node_comparison_colors "#cf191b" "#33a02c" --output plot_demo/graph/png/universal/sense_db_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_db_R2 --node_comparison_colors "#cf191b" "#33a02c" --output plot_demo/graph/png/universal/sense_db_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_dcmv_R1 --output plot_demo/graph/png/universal/sense_dcmv_R1.png --node_comparison_colors "#cf191b" "#ffe669" --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_dcmv_R2 --output plot_demo/graph/png/universal/sense_dcmv_R2.png --node_comparison_colors "#cf191b" "#ffe669" --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot the 3D histograms (PNG) ###
python histogram.py --input data_demo/output/db_R1 --output plot_demo/histogram/png/db_R1_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_demo/output/db_R1 --output plot_demo/histogram/png/db_R1_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_demo/output/db_R1 --output plot_demo/histogram/png/db_R1_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_demo/output/db_R2 --output plot_demo/histogram/png/db_R2_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_demo/output/db_R2 --output plot_demo/histogram/png/db_R2_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_demo/output/db_R2 --output plot_demo/histogram/png/db_R2_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_demo/output/sense_R1 --output plot_demo/histogram/png/sense_R1_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_demo/output/sense_R1 --output plot_demo/histogram/png/sense_R1_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_demo/output/sense_R1 --output plot_demo/histogram/png/sense_R1_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_demo/output/sense_R2 --output plot_demo/histogram/png/sense_R2_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_demo/output/sense_R2 --output plot_demo/histogram/png/sense_R2_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_demo/output/sense_R2 --output plot_demo/histogram/png/sense_R2_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_demo/output/dcmv_R1 --output plot_demo/histogram/png/dcmv_R1_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_demo/output/dcmv_R1 --output plot_demo/histogram/png/dcmv_R1_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_demo/output/dcmv_R1 --output plot_demo/histogram/png/dcmv_R1_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_demo/output/dcmv_R2 --output plot_demo/histogram/png/dcmv_R2_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_demo/output/dcmv_R2 --output plot_demo/histogram/png/dcmv_R2_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_demo/output/dcmv_R2 --output plot_demo/histogram/png/dcmv_R2_deletion.png --color "#8080ff" --variation_type deletion --label_type relative

### Plot graphs with common Kamada layout (PNG) ###
python graph.py `
  --input data_demo/output/db_R1 data_demo/output/db_R2 data_demo/output/sense_R1 data_demo/output/sense_R2 data_demo/output/dcmv_R1 data_demo/output/dcmv_R2 data_demo/output/sense_db_R1 data_demo/output/sense_db_R2 data_demo/output/sense_dcmv_R1 data_demo/output/sense_dcmv_R2 `
  --output plot_demo/graph/png/kamada_common/db_R1.png plot_demo/graph/png/kamada_common/db_R2.png plot_demo/graph/png/kamada_common/sense_R1.png plot_demo/graph/png/kamada_common/sense_R2.png plot_demo/graph/png/kamada_common/dcmv_R1.png plot_demo/graph/png/kamada_common/dcmv_R2.png plot_demo/graph/png/kamada_common/sense_db_R1.png plot_demo/graph/png/kamada_common/sense_db_R2.png plot_demo/graph/png/kamada_common/sense_dcmv_R1.png plot_demo/graph/png/kamada_common/sense_dcmv_R2.png `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout kamada_layout --width 2400 --height 1800

### Plot graphs with common radial layout (PNG) ###
python graph.py `
  --input data_demo/output/db_R1 data_demo/output/db_R2 data_demo/output/sense_R1 data_demo/output/sense_R2 data_demo/output/dcmv_R1 data_demo/output/dcmv_R2 data_demo/output/sense_db_R1 data_demo/output/sense_db_R2 data_demo/output/sense_dcmv_R1 data_demo/output/sense_dcmv_R2 `
  --output plot_demo/graph/png/radial_common/db_R1.png plot_demo/graph/png/radial_common/db_R2.png plot_demo/graph/png/radial_common/sense_R1.png plot_demo/graph/png/radial_common/sense_R2.png plot_demo/graph/png/radial_common/dcmv_R1.png plot_demo/graph/png/radial_common/dcmv_R2.png plot_demo/graph/png/radial_common/sense_db_R1.png plot_demo/graph/png/radial_common/sense_db_R2.png plot_demo/graph/png/radial_common/sense_dcmv_R1.png plot_demo/graph/png/radial_common/sense_dcmv_R2.png `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --range_x -20 20 --range_y -20 20 --layout radial_layout --width 2400 --height 1800

### Plot the graphs (PDF) ###
python graph.py --input data_demo/output/db_R1 --output plot_demo/graph/pdf/universal/db_R1.pdf --title "BranchΔ (R1)" --legends size variation_type --legend_x_shift_px 100 --legend_y_shift_px 0 --margin_right_px 750 --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5  --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/db_R2 --output plot_demo/graph/pdf/universal/db_R2.pdf --title "BranchΔ (R2)" --legends size variation_type --legend_x_shift_px 100 --legend_y_shift_px 0 --margin_right_px 750 --layout universal_layout --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_R1 --output plot_demo/graph/pdf/universal/sense_R1.pdf --title "Sense (R1)" --legends size variation_type --legend_x_shift_px 100 --legend_y_shift_px 0 --margin_right_px 750 --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_R2 --output plot_demo/graph/pdf/universal/sense_R2.pdf --title "Sense (R2)" --legends size variation_type --legend_x_shift_px 100 --legend_y_shift_px 0 --margin_right_px 750 --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/dcmv_R1 --output plot_demo/graph/pdf/universal/dcmv_R1.pdf --title "pCMVΔ (R1)" --legends size variation_type --legend_x_shift_px 100 --legend_y_shift_px 0 --margin_right_px 750 --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/dcmv_R2 --output plot_demo/graph/pdf/universal/dcmv_R2.pdf --title "pCMVΔ (R2)" --legends size variation_type --legend_x_shift_px 100 --legend_y_shift_px 0 --margin_right_px 750 --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot the comparison graphs (PDF) ###
python graph.py --input data_demo/output/sense_db_R1 --node_comparison_colors "#cf191b" "#33a02c" --output plot_demo/graph/pdf/universal/sense_db_R1.pdf --title "Sense / BranchΔ (R1)"--legend_colorbar_scale 3 --legend_x_shift_px 100 --legend_y_shift_px -100 --margin_right_px 750 --legends freq_ratio_continuous --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_db_R2 --node_comparison_colors "#cf191b" "#33a02c" --output plot_demo/graph/pdf/universal/sense_db_R2.pdf --title "Sense / BranchΔ (R2)"--legend_colorbar_scale 3 --legend_x_shift_px 100 --legend_y_shift_px -100 --margin_right_px 750 --legends freq_ratio_continuous --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_dcmv_R1 --output plot_demo/graph/pdf/universal/sense_dcmv_R1.pdf --node_comparison_colors "#cf191b" "#ffe669" --title "Sense / pCMVΔ (R1)"--legend_colorbar_scale 3 --legend_x_shift_px 100 --legend_y_shift_px -100 --margin_right_px 750 --legends freq_ratio_continuous --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_dcmv_R2 --output plot_demo/graph/pdf/universal/sense_dcmv_R2.pdf --node_comparison_colors "#cf191b" "#ffe669" --title "Sense / pCMVΔ (R2)"--legend_colorbar_scale 3 --legend_x_shift_px 100 --legend_y_shift_px -100 --margin_right_px 750 --legends freq_ratio_continuous --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot the 3D histograms (PDF) ###
python histogram.py --input data_demo/output/db_R1 --output plot_demo/histogram/pdf/db_R1_substitution.pdf --title "BranchΔ (R1)`nSubstitution" --margin_top 250 --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_demo/output/db_R1 --output plot_demo/histogram/pdf/db_R1_insertion.pdf --title "BranchΔ (R1)`nInsertion" --margin_top 250 --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_demo/output/db_R1 --output plot_demo/histogram/pdf/db_R1_deletion.pdf --title "BranchΔ (R1)`nDeletion" --margin_top 250 --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_demo/output/db_R2 --output plot_demo/histogram/pdf/db_R2_substitution.pdf --title "BranchΔ (R2)`nSubstitution" --margin_top 250 --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_demo/output/db_R2 --output plot_demo/histogram/pdf/db_R2_insertion.pdf --title "BranchΔ (R2)`nInsertion" --margin_top 250 --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_demo/output/db_R2 --output plot_demo/histogram/pdf/db_R2_deletion.pdf --title "BranchΔ (R2)`nDeletion" --margin_top 250 --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_demo/output/sense_R1 --output plot_demo/histogram/pdf/sense_R1_substitution.pdf --title "Sense (R1)`nSubstitution" --margin_top 250 --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_demo/output/sense_R1 --output plot_demo/histogram/pdf/sense_R1_insertion.pdf --title "Sense (R1)`nInsertion" --margin_top 250 --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_demo/output/sense_R1 --output plot_demo/histogram/pdf/sense_R1_deletion.pdf --title "Sense (R1)`nDeletion" --margin_top 250 --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_demo/output/sense_R2 --output plot_demo/histogram/pdf/sense_R2_substitution.pdf --title "Sense (R2)`nSubstitution" --margin_top 250 --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_demo/output/sense_R2 --output plot_demo/histogram/pdf/sense_R2_insertion.pdf --title "Sense (R2)`nInsertion" --margin_top 250 --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_demo/output/sense_R2 --output plot_demo/histogram/pdf/sense_R2_deletion.pdf --title "Sense (R2)`nDeletion" --margin_top 250 --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_demo/output/dcmv_R1 --output plot_demo/histogram/pdf/dcmv_R1_substitution.pdf --title "pCMVΔ (R1)`nSubstitution" --margin_top 250 --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_demo/output/dcmv_R1 --output plot_demo/histogram/pdf/dcmv_R1_insertion.pdf --title "pCMVΔ (R1)`nInsertion" --margin_top 250 --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_demo/output/dcmv_R1 --output plot_demo/histogram/pdf/dcmv_R1_deletion.pdf --title "pCMVΔ (R1)`nDeletion" --margin_top 250 --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_demo/output/dcmv_R2 --output plot_demo/histogram/pdf/dcmv_R2_substitution.pdf --title "pCMVΔ (R2)`nSubstitution" --margin_top 250 --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_demo/output/dcmv_R2 --output plot_demo/histogram/pdf/dcmv_R2_insertion.pdf --title "pCMVΔ (R2)`nInsertion" --margin_top 250 --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_demo/output/dcmv_R2 --output plot_demo/histogram/pdf/dcmv_R2_deletion.pdf --title "pCMVΔ (R2)`nDeletion" --margin_top 250 --color "#8080ff" --variation_type deletion --label_type relative

### Plot graphs with common Kamada layout (PDF) ###
python graph.py `
  --input data_demo/output/db_R1 data_demo/output/db_R2 data_demo/output/sense_R1 data_demo/output/sense_R2 data_demo/output/dcmv_R1 data_demo/output/dcmv_R2 data_demo/output/sense_db_R1 data_demo/output/sense_db_R2 data_demo/output/sense_dcmv_R1 data_demo/output/sense_dcmv_R2 `
  --output plot_demo/graph/pdf/kamada_common/db_R1.pdf plot_demo/graph/pdf/kamada_common/db_R2.pdf plot_demo/graph/pdf/kamada_common/sense_R1.pdf plot_demo/graph/pdf/kamada_common/sense_R2.pdf plot_demo/graph/pdf/kamada_common/dcmv_R1.pdf plot_demo/graph/pdf/kamada_common/dcmv_R2.pdf plot_demo/graph/pdf/kamada_common/sense_db_R1.pdf plot_demo/graph/pdf/kamada_common/sense_db_R2.pdf plot_demo/graph/pdf/kamada_common/sense_dcmv_R1.pdf plot_demo/graph/pdf/kamada_common/sense_dcmv_R2.pdf `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout kamada_layout --width 2400 --height 1800

### Plot graphs with common radial layout (PDF) ###
python graph.py `
  --input data_demo/output/db_R1 data_demo/output/db_R2 data_demo/output/sense_R1 data_demo/output/sense_R2 data_demo/output/dcmv_R1 data_demo/output/dcmv_R2 data_demo/output/sense_db_R1 data_demo/output/sense_db_R2 data_demo/output/sense_dcmv_R1 data_demo/output/sense_dcmv_R2 `
  --output plot_demo/graph/pdf/radial_common/db_R1.pdf plot_demo/graph/pdf/radial_common/db_R2.pdf plot_demo/graph/pdf/radial_common/sense_R1.pdf plot_demo/graph/pdf/radial_common/sense_R2.pdf plot_demo/graph/pdf/radial_common/dcmv_R1.pdf plot_demo/graph/pdf/radial_common/dcmv_R2.pdf plot_demo/graph/pdf/radial_common/sense_db_R1.pdf plot_demo/graph/pdf/radial_common/sense_db_R2.pdf plot_demo/graph/pdf/radial_common/sense_dcmv_R1.pdf plot_demo/graph/pdf/radial_common/sense_dcmv_R2.pdf `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --range_x -20 20 --range_y -20 20 --layout radial_layout --width 2400 --height 1800

### Plot the graphs (HTML) ###
python graph.py --input data_demo/output/db_R1 --output plot_demo/graph/html/universal/db_R1.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5  --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/db_R2 --output plot_demo/graph/html/universal/db_R2.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_R1 --output plot_demo/graph/html/universal/sense_R1.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_R2 --output plot_demo/graph/html/universal/sense_R2.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/dcmv_R1 --output plot_demo/graph/html/universal/dcmv_R1.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/dcmv_R2 --output plot_demo/graph/html/universal/dcmv_R2.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot the comparison graphs (HTML) ###
python graph.py --input data_demo/output/sense_db_R1 --node_comparison_colors "#cf191b" "#33a02c" --output plot_demo/graph/html/universal/sense_db_R1.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_db_R2 --node_comparison_colors "#cf191b" "#33a02c" --output plot_demo/graph/html/universal/sense_db_R2.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_dcmv_R1 --output plot_demo/graph/html/universal/sense_dcmv_R1.html --node_comparison_colors "#cf191b" "#ffe669" --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_demo/output/sense_dcmv_R2 --output plot_demo/graph/html/universal/sense_dcmv_R2.html --node_comparison_colors "#cf191b" "#ffe669" --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot graphs with common Kamada layout (HTML) ###
python graph.py `
  --input data_demo/output/db_R1 data_demo/output/db_R2 data_demo/output/sense_R1 data_demo/output/sense_R2 data_demo/output/dcmv_R1 data_demo/output/dcmv_R2 data_demo/output/sense_db_R1 data_demo/output/sense_db_R2 data_demo/output/sense_dcmv_R1 data_demo/output/sense_dcmv_R2 `
  --output plot_demo/graph/html/kamada_common/db_R1.html plot_demo/graph/html/kamada_common/db_R2.html plot_demo/graph/html/kamada_common/sense_R1.html plot_demo/graph/html/kamada_common/sense_R2.html plot_demo/graph/html/kamada_common/dcmv_R1.html plot_demo/graph/html/kamada_common/dcmv_R2.html plot_demo/graph/html/kamada_common/sense_db_R1.html plot_demo/graph/html/kamada_common/sense_db_R2.html plot_demo/graph/html/kamada_common/sense_dcmv_R1.html plot_demo/graph/html/kamada_common/sense_dcmv_R2.html `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout kamada_layout --width 2400 --height 1800

### Plot graphs with common radial layout (HTML) ###
python graph.py `
  --input data_demo/output/db_R1 data_demo/output/db_R2 data_demo/output/sense_R1 data_demo/output/sense_R2 data_demo/output/dcmv_R1 data_demo/output/dcmv_R2 data_demo/output/sense_db_R1 data_demo/output/sense_db_R2 data_demo/output/sense_dcmv_R1 data_demo/output/sense_dcmv_R2 `
  --output plot_demo/graph/html/radial_common/db_R1.html plot_demo/graph/html/radial_common/db_R2.html plot_demo/graph/html/radial_common/sense_R1.html plot_demo/graph/html/radial_common/sense_R2.html plot_demo/graph/html/radial_common/dcmv_R1.html plot_demo/graph/html/radial_common/dcmv_R2.html plot_demo/graph/html/radial_common/sense_db_R1.html plot_demo/graph/html/radial_common/sense_db_R2.html plot_demo/graph/html/radial_common/sense_dcmv_R1.html plot_demo/graph/html/radial_common/sense_dcmv_R2.html `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --range_x -20 20 --range_y -20 20 --layout radial_layout --width 2400 --height 1800

### Arrange the histograms into a grid in a PPTX file ###
python pptx1.py --input `
  plot_demo/histogram/png/db_R1_substitution.png `
  plot_demo/histogram/png/db_R1_insertion.png `
  plot_demo/histogram/png/db_R1_deletion.png `
  plot_demo/histogram/png/sense_R1_substitution.png `
  plot_demo/histogram/png/sense_R1_insertion.png `
  plot_demo/histogram/png/sense_R1_deletion.png `
  plot_demo/histogram/png/dcmv_R1_substitution.png `
  plot_demo/histogram/png/dcmv_R1_insertion.png `
  plot_demo/histogram/png/dcmv_R1_deletion.png `
  plot_demo/histogram/png/db_R2_substitution.png `
  plot_demo/histogram/png/db_R2_insertion.png `
  plot_demo/histogram/png/db_R2_deletion.png `
  plot_demo/histogram/png/sense_R2_substitution.png `
  plot_demo/histogram/png/sense_R2_insertion.png `
  plot_demo/histogram/png/sense_R2_deletion.png `
  plot_demo/histogram/png/dcmv_R2_substitution.png `
  plot_demo/histogram/png/dcmv_R2_insertion.png `
  plot_demo/histogram/png/dcmv_R2_deletion.png `
  --num_grids 2 --num_rows 3 3 --num_cols 3 3 `
  --image_height_spacing_pt 0 --image_width_spacing_pt 0 `
  --grid_height_spacing_pt 200 --total_width_frac 1 1 `
  --image_labels `
  "B/S/R1" "B/I/R1" "B/D/R1" `
  "S/S/R1" "S/I/R1" "S/D/R1" `
  "C/S/R1" "C/I/R1" "C/D/R1" `
  "B/S/R2" "B/I/R2" "B/D/R2" `
  "S/S/R2" "S/I/R2" "S/D/R2" `
  "C/S/R2" "C/I/R2" "C/D/R2" `
  --image_label_font_size_pt 10 --image_label_height_pt 20 `
  --image_label_width_pt 100 --margin_labels_top `
  "Substitution (R1)" "Insertion (R1)" "Deletion (R1)" `
  "Substitution (R2)" "Insertion (R2)" "Deletion (R2)" `
  --margin_top_height_pt 100 --margin_top_font_size_pt 20 `
  --margin_labels_left "BranchΔ" "Sense" "pCMVΔ" "BranchΔ" "Sense" "pCMVΔ" `
  --margin_left_width_pt 100 --margin_left_spill_over_pt 0 0 --margin_left_font_size_pt 20 `
  --title "Demo 2-DSB Histograms" --title_font_size_pt 40 --title_height_pt 40 `
  --output plot_demo/pptx/histograms.pptx