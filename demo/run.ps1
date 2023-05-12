### Preprocessing ###
DSBplot preprocess --input data/input/fastq/db1_R1.fq data/input/fastq/db2_R1.fq data/input/fastq/db3_R1.fq data/input/fastq/db4_R1.fq --ref_seq_file data/input/ref_seq/2DSB_R1_branch.fa --dsb_pos 67 --output data/output/db_R1 --label "BranchΔ" --total_reads 3000 3000 3000 3000
DSBplot-preprocess --input data/input/fastq/db1_R1.fq data/input/fastq/db2_R1.fq data/input/fastq/db3_R1.fq data/input/fastq/db4_R1.fq --ref_seq_file data/input/ref_seq/2DSB_R1_branch.fa --dsb_pos 67 --output data/output/db_R1 --label "BranchΔ" --total_reads 3000 3000 3000 3000
DSBplot-preprocess --input data/input/fastq/db1_R2.fq data/input/fastq/db2_R2.fq data/input/fastq/db3_R2.fq data/input/fastq/db4_R2.fq --ref_seq_file data/input/ref_seq/2DSB_R2_branch.fa --dsb_pos 46 --output data/output/db_R2 --label "BranchΔ" --total_reads 3000 3000 3000 3000
DSBplot-preprocess --input data/input/fastq/sense1_R1.fq data/input/fastq/sense2_R1.fq data/input/fastq/sense3_R1.fq data/input/fastq/sense4_R1.fq --ref_seq_file data/input/ref_seq/2DSB_R1_sense.fa --dsb_pos 67 --output data/output/sense_R1 --label "Sense" --total_reads 3000 3000 3000 3000
DSBplot-preprocess --input data/input/fastq/sense1_R2.fq data/input/fastq/sense2_R2.fq data/input/fastq/sense3_R2.fq data/input/fastq/sense4_R2.fq --ref_seq_file data/input/ref_seq/2DSB_R2_sense.fa --dsb_pos 46 --output data/output/sense_R2 --label "Sense" --total_reads 3000 3000 3000 3000
DSBplot-preprocess --input data/input/fastq/dcmv1_R1.fq data/input/fastq/dcmv2_R1.fq data/input/fastq/dcmv3_R1.fq data/input/fastq/dcmv4_R1.fq --ref_seq_file data/input/ref_seq/2DSB_R1_cmv.fa --dsb_pos 67 --output data/output/dcmv_R1 --label "pCMVΔ" --total_reads 3000 3000 3000 3000
DSBplot-preprocess --input data/input/fastq/dcmv1_R2.fq data/input/fastq/dcmv2_R2.fq data/input/fastq/dcmv3_R2.fq data/input/fastq/dcmv4_R2.fq --ref_seq_file data/input/ref_seq/2DSB_R2_cmv.fa --dsb_pos 46 --output data/output/dcmv_R2 --label "pCMVΔ" --total_reads 3000 3000 3000 3000

# Example of running preprocessing stages separately for one experiment
DSBplot-preprocess --input data/input/fastq/db1_R1.fq data/input/fastq/db2_R1.fq data/input/fastq/db3_R1.fq data/input/fastq/db4_R1.fq --ref_seq_file data/input/ref_seq/2DSB_R1_branch.fa --output data/output/db_R1 --stages 0_align
DSBplot-preprocess --ref_seq_file data/input/ref_seq/2DSB_R1_branch.fa --dsb_pos 67 --output data/output/db_R1 --stages 1_filter
DSBplot-preprocess --output data/output/db_R1 --stages 2_combine
DSBplot-preprocess --ref_seq_file data/input/ref_seq/2DSB_R1_branch.fa --dsb_pos 67 --output data/output/db_R1 --label db_R1 --total_reads 3000 3000 3000 3000 --stages 3_window
DSBplot-preprocess --output data/output/db_R1 --stages 4_graph
DSBplot-preprocess --output data/output/db_R1 --stages 5_histogram

### Comparison data ###
DSBplot-comparison --input data/output/sense_R1 data/output/db_R1 --output data/output/sense_db_R1
DSBplot-comparison --input data/output/sense_R2 data/output/db_R2 --output data/output/sense_db_R2
DSBplot-comparison --input data/output/sense_R1 data/output/dcmv_R1 --output data/output/sense_dcmv_R1
DSBplot-comparison --input data/output/sense_R2 data/output/dcmv_R2 --output data/output/sense_dcmv_R2

### Plot the graphs (PNG) ###
DSBplot-graph --input data/output/db_R1 --output plots/graph/png/universal/db_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5  --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/db_R2 --output plots/graph/png/universal/db_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/sense_R1 --output plots/graph/png/universal/sense_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/sense_R2 --output plots/graph/png/universal/sense_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/dcmv_R1 --output plots/graph/png/universal/dcmv_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/dcmv_R2 --output plots/graph/png/universal/dcmv_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot the comparison graphs (PNG) ###
DSBplot-graph --input data/output/sense_db_R1 --node_comparison_colors "#cf191b" "#33a02c" --output plots/graph/png/universal/sense_db_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/sense_db_R2 --node_comparison_colors "#cf191b" "#33a02c" --output plots/graph/png/universal/sense_db_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/sense_dcmv_R1 --output plots/graph/png/universal/sense_dcmv_R1.png --node_comparison_colors "#cf191b" "#ffe669" --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/sense_dcmv_R2 --output plots/graph/png/universal/sense_dcmv_R2.png --node_comparison_colors "#cf191b" "#ffe669" --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot the 3D histograms (PNG) ###
DSBplot-histogram --input data/output/db_R1 --output plots/histogram/png/db_R1_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/db_R1 --output plots/histogram/png/db_R1_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/db_R1 --output plots/histogram/png/db_R1_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/db_R2 --output plots/histogram/png/db_R2_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/db_R2 --output plots/histogram/png/db_R2_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/db_R2 --output plots/histogram/png/db_R2_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/sense_R1 --output plots/histogram/png/sense_R1_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/sense_R1 --output plots/histogram/png/sense_R1_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/sense_R1 --output plots/histogram/png/sense_R1_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/sense_R2 --output plots/histogram/png/sense_R2_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/sense_R2 --output plots/histogram/png/sense_R2_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/sense_R2 --output plots/histogram/png/sense_R2_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/dcmv_R1 --output plots/histogram/png/dcmv_R1_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/dcmv_R1 --output plots/histogram/png/dcmv_R1_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/dcmv_R1 --output plots/histogram/png/dcmv_R1_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/dcmv_R2 --output plots/histogram/png/dcmv_R2_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/dcmv_R2 --output plots/histogram/png/dcmv_R2_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/dcmv_R2 --output plots/histogram/png/dcmv_R2_deletion.png --color "#8080ff" --variation_type deletion --label_type relative

### Plot graphs with common Universal layout (PNG) ###
DSBplot-graph `
  --input data/output/db_R1 data/output/db_R2 data/output/sense_R1 data/output/sense_R2 data/output/dcmv_R1 data/output/dcmv_R2 data/output/sense_db_R1 data/output/sense_db_R2 data/output/sense_dcmv_R1 data/output/sense_dcmv_R2 `
  --output plots/graph/png/universal_common/db_R1.png plots/graph/png/universal_common/db_R2.png plots/graph/png/universal_common/sense_R1.png plots/graph/png/universal_common/sense_R2.png plots/graph/png/universal_common/dcmv_R1.png plots/graph/png/universal_common/dcmv_R2.png plots/graph/png/universal_common/sense_db_R1.png plots/graph/png/universal_common/sense_db_R2.png plots/graph/png/universal_common/sense_dcmv_R1.png plots/graph/png/universal_common/sense_dcmv_R2.png `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout kamada_layout --width 2400 --height 1800

### Plot graphs with common Kamada layout (PNG) ###
DSBplot-graph `
  --input data/output/db_R1 data/output/db_R2 data/output/sense_R1 data/output/sense_R2 data/output/dcmv_R1 data/output/dcmv_R2 data/output/sense_db_R1 data/output/sense_db_R2 data/output/sense_dcmv_R1 data/output/sense_dcmv_R2 `
  --output plots/graph/png/kamada_common/db_R1.png plots/graph/png/kamada_common/db_R2.png plots/graph/png/kamada_common/sense_R1.png plots/graph/png/kamada_common/sense_R2.png plots/graph/png/kamada_common/dcmv_R1.png plots/graph/png/kamada_common/dcmv_R2.png plots/graph/png/kamada_common/sense_db_R1.png plots/graph/png/kamada_common/sense_db_R2.png plots/graph/png/kamada_common/sense_dcmv_R1.png plots/graph/png/kamada_common/sense_dcmv_R2.png `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout kamada_layout --width 2400 --height 1800

### Plot graphs with common Radial layout (PNG) ###
DSBplot-graph `
  --input data/output/db_R1 data/output/db_R2 data/output/sense_R1 data/output/sense_R2 data/output/dcmv_R1 data/output/dcmv_R2 data/output/sense_db_R1 data/output/sense_db_R2 data/output/sense_dcmv_R1 data/output/sense_dcmv_R2 `
  --output plots/graph/png/radial_common/db_R1.png plots/graph/png/radial_common/db_R2.png plots/graph/png/radial_common/sense_R1.png plots/graph/png/radial_common/sense_R2.png plots/graph/png/radial_common/dcmv_R1.png plots/graph/png/radial_common/dcmv_R2.png plots/graph/png/radial_common/sense_db_R1.png plots/graph/png/radial_common/sense_db_R2.png plots/graph/png/radial_common/sense_dcmv_R1.png plots/graph/png/radial_common/sense_dcmv_R2.png `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --range_x -20 20 --range_y -20 20 --layout radial_layout --width 2400 --height 1800

### Plot the graphs (PDF) ###
DSBplot-graph --input data/output/db_R1 --output plots/graph/pdf/universal/db_R1.pdf --title "BranchΔ (R1)" --legends size variation_type --legend_x_shift_px 100 --legend_y_shift_px 0 --margin_right_px 750 --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5  --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/db_R2 --output plots/graph/pdf/universal/db_R2.pdf --title "BranchΔ (R2)" --legends size variation_type --legend_x_shift_px 100 --legend_y_shift_px 0 --margin_right_px 750 --layout universal_layout --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/sense_R1 --output plots/graph/pdf/universal/sense_R1.pdf --title "Sense (R1)" --legends size variation_type --legend_x_shift_px 100 --legend_y_shift_px 0 --margin_right_px 750 --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/sense_R2 --output plots/graph/pdf/universal/sense_R2.pdf --title "Sense (R2)" --legends size variation_type --legend_x_shift_px 100 --legend_y_shift_px 0 --margin_right_px 750 --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/dcmv_R1 --output plots/graph/pdf/universal/dcmv_R1.pdf --title "pCMVΔ (R1)" --legends size variation_type --legend_x_shift_px 100 --legend_y_shift_px 0 --margin_right_px 750 --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/dcmv_R2 --output plots/graph/pdf/universal/dcmv_R2.pdf --title "pCMVΔ (R2)" --legends size variation_type --legend_x_shift_px 100 --legend_y_shift_px 0 --margin_right_px 750 --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot the comparison graphs (PDF) ###
DSBplot-graph --input data/output/sense_db_R1 --node_comparison_colors "#cf191b" "#33a02c" --output plots/graph/pdf/universal/sense_db_R1.pdf --title "Sense / BranchΔ (R1)" --legend_colorbar_scale 3 --legend_x_shift_px 100 --legend_y_shift_px -100 --margin_right_px 750 --legends freq_ratio_continuous --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/sense_db_R2 --node_comparison_colors "#cf191b" "#33a02c" --output plots/graph/pdf/universal/sense_db_R2.pdf --title "Sense / BranchΔ (R2)" --legend_colorbar_scale 3 --legend_x_shift_px 100 --legend_y_shift_px -100 --margin_right_px 750 --legends freq_ratio_continuous --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/sense_dcmv_R1 --output plots/graph/pdf/universal/sense_dcmv_R1.pdf --node_comparison_colors "#cf191b" "#ffe669" --title "Sense / pCMVΔ (R1)" --legend_colorbar_scale 3 --legend_x_shift_px 100 --legend_y_shift_px -100 --margin_right_px 750 --legends freq_ratio_continuous --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/sense_dcmv_R2 --output plots/graph/pdf/universal/sense_dcmv_R2.pdf --node_comparison_colors "#cf191b" "#ffe669" --title "Sense / pCMVΔ (R2)" --legend_colorbar_scale 3 --legend_x_shift_px 100 --legend_y_shift_px -100 --margin_right_px 750 --legends freq_ratio_continuous --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot the 3D histograms (PDF) ###
DSBplot-histogram --input data/output/db_R1 --output plots/histogram/pdf/db_R1_substitution.pdf --title "BranchΔ (R1)`nSubstitution" --margin_top 250 --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/db_R1 --output plots/histogram/pdf/db_R1_insertion.pdf --title "BranchΔ (R1)`nInsertion" --margin_top 250 --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/db_R1 --output plots/histogram/pdf/db_R1_deletion.pdf --title "BranchΔ (R1)`nDeletion" --margin_top 250 --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/db_R2 --output plots/histogram/pdf/db_R2_substitution.pdf --title "BranchΔ (R2)`nSubstitution" --margin_top 250 --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/db_R2 --output plots/histogram/pdf/db_R2_insertion.pdf --title "BranchΔ (R2)`nInsertion" --margin_top 250 --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/db_R2 --output plots/histogram/pdf/db_R2_deletion.pdf --title "BranchΔ (R2)`nDeletion" --margin_top 250 --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/sense_R1 --output plots/histogram/pdf/sense_R1_substitution.pdf --title "Sense (R1)`nSubstitution" --margin_top 250 --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/sense_R1 --output plots/histogram/pdf/sense_R1_insertion.pdf --title "Sense (R1)`nInsertion" --margin_top 250 --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/sense_R1 --output plots/histogram/pdf/sense_R1_deletion.pdf --title "Sense (R1)`nDeletion" --margin_top 250 --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/sense_R2 --output plots/histogram/pdf/sense_R2_substitution.pdf --title "Sense (R2)`nSubstitution" --margin_top 250 --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/sense_R2 --output plots/histogram/pdf/sense_R2_insertion.pdf --title "Sense (R2)`nInsertion" --margin_top 250 --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/sense_R2 --output plots/histogram/pdf/sense_R2_deletion.pdf --title "Sense (R2)`nDeletion" --margin_top 250 --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/dcmv_R1 --output plots/histogram/pdf/dcmv_R1_substitution.pdf --title "pCMVΔ (R1)`nSubstitution" --margin_top 250 --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/dcmv_R1 --output plots/histogram/pdf/dcmv_R1_insertion.pdf --title "pCMVΔ (R1)`nInsertion" --margin_top 250 --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/dcmv_R1 --output plots/histogram/pdf/dcmv_R1_deletion.pdf --title "pCMVΔ (R1)`nDeletion" --margin_top 250 --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/dcmv_R2 --output plots/histogram/pdf/dcmv_R2_substitution.pdf --title "pCMVΔ (R2)`nSubstitution" --margin_top 250 --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/dcmv_R2 --output plots/histogram/pdf/dcmv_R2_insertion.pdf --title "pCMVΔ (R2)`nInsertion" --margin_top 250 --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/dcmv_R2 --output plots/histogram/pdf/dcmv_R2_deletion.pdf --title "pCMVΔ (R2)`nDeletion" --margin_top 250 --color "#8080ff" --variation_type deletion --label_type relative

### Plot graphs with common Universal layout (PDF) ###
DSBplot-graph `
  --input data/output/db_R1 data/output/db_R2 data/output/sense_R1 data/output/sense_R2 data/output/dcmv_R1 data/output/dcmv_R2 data/output/sense_db_R1 data/output/sense_db_R2 data/output/sense_dcmv_R1 data/output/sense_dcmv_R2 `
  --title "BranchΔ (R1)" "BranchΔ (R2)" "SenseΔ (R1)" "pCMVΔ (R1)" "pCMVΔ (R2)" "pCMVΔ (R2)" "SenseΔ / BranchΔ (R1)" "SenseΔ / BranchΔ (R2)" "SenseΔ / pCMVΔ (R1)" "SenseΔ / pCMVΔ (R2)" `
  --output plots/graph/pdf/universal_common/db_R1.pdf plots/graph/pdf/universal_common/db_R2.pdf plots/graph/pdf/universal_common/sense_R1.pdf plots/graph/pdf/universal_common/sense_R2.pdf plots/graph/pdf/universal_common/dcmv_R1.pdf plots/graph/pdf/universal_common/dcmv_R2.pdf plots/graph/pdf/universal_common/sense_db_R1.pdf plots/graph/pdf/universal_common/sense_db_R2.pdf plots/graph/pdf/universal_common/sense_dcmv_R1.pdf plots/graph/pdf/universal_common/sense_dcmv_R2.pdf `
  --range_x -12 13 --range_y -23 20 `
  --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 `
  --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_deletion_max_tick 19 `
  --legends freq_ratio_continuous variation_type `
  --legend_colorbar_scale 3 --legend_x_shift_px 100 --legend_y_shift_px -100 `
  --margin_right_px 1000 --margin_top_px 400 `
  --font_size_scale 3 --legend_colorbar_scale 3 `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout universal_layout --width 2400 --height 1800

### Plot graphs with common Kamada layout (PDF) ###
DSBplot-graph `
  --input data/output/db_R1 data/output/db_R2 data/output/sense_R1 data/output/sense_R2 data/output/dcmv_R1 data/output/dcmv_R2 data/output/sense_db_R1 data/output/sense_db_R2 data/output/sense_dcmv_R1 data/output/sense_dcmv_R2 `
  --title "BranchΔ (R1)" "BranchΔ (R2)" "SenseΔ (R1)" "pCMVΔ (R1)" "pCMVΔ (R2)" "pCMVΔ (R2)" "SenseΔ / BranchΔ (R1)" "SenseΔ / BranchΔ (R2)" "SenseΔ / pCMVΔ (R1)" "SenseΔ / pCMVΔ (R2)" `
  --output plots/graph/pdf/kamada_common/db_R1.pdf plots/graph/pdf/kamada_common/db_R2.pdf plots/graph/pdf/kamada_common/sense_R1.pdf plots/graph/pdf/kamada_common/sense_R2.pdf plots/graph/pdf/kamada_common/dcmv_R1.pdf plots/graph/pdf/kamada_common/dcmv_R2.pdf plots/graph/pdf/kamada_common/sense_db_R1.pdf plots/graph/pdf/kamada_common/sense_db_R2.pdf plots/graph/pdf/kamada_common/sense_dcmv_R1.pdf plots/graph/pdf/kamada_common/sense_dcmv_R2.pdf `
  --range_x 0.15 0.85 --range_y 0.15 0.85 `
  --legends freq_ratio_continuous variation_type `
  --font_size_scale 4 --legend_colorbar_scale 3 `
  --margin_top_px 400 --margin_right_px 1000 `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout kamada_layout --width 2400 --height 1800

### Plot graphs with common Radial layout (PDF) ###
DSBplot-graph `
  --input data/output/db_R1 data/output/db_R2 data/output/sense_R1 data/output/sense_R2 data/output/dcmv_R1 data/output/dcmv_R2 data/output/sense_db_R1 data/output/sense_db_R2 data/output/sense_dcmv_R1 data/output/sense_dcmv_R2 `
  --title "BranchΔ (R1)" "BranchΔ (R2)" "SenseΔ (R1)" "pCMVΔ (R1)" "pCMVΔ (R2)" "pCMVΔ (R2)" "SenseΔ / BranchΔ (R1)" "SenseΔ / BranchΔ (R2)" "SenseΔ / pCMVΔ (R1)" "SenseΔ / pCMVΔ (R2)" `
  --output plots/graph/pdf/radial_common/db_R1.pdf plots/graph/pdf/radial_common/db_R2.pdf plots/graph/pdf/radial_common/sense_R1.pdf plots/graph/pdf/radial_common/sense_R2.pdf plots/graph/pdf/radial_common/dcmv_R1.pdf plots/graph/pdf/radial_common/dcmv_R2.pdf plots/graph/pdf/radial_common/sense_db_R1.pdf plots/graph/pdf/radial_common/sense_db_R2.pdf plots/graph/pdf/radial_common/sense_dcmv_R1.pdf plots/graph/pdf/radial_common/sense_dcmv_R2.pdf `
  --range_x 0.15 0.85 --range_y 0.15 0.85 `
  --legends freq_ratio_continuous variation_type `
  --font_size_scale 4 --legend_colorbar_scale 3 `
  --margin_top_px 400 --margin_right_px 1000 `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --range_x -20 20 --range_y -20 20 --layout radial_layout --width 2400 --height 1800

### Plot the graphs (HTML) ###
DSBplot-graph --input data/output/db_R1 --output plots/graph/html/universal/db_R1.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5  --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/db_R2 --output plots/graph/html/universal/db_R2.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/sense_R1 --output plots/graph/html/universal/sense_R1.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/sense_R2 --output plots/graph/html/universal/sense_R2.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/dcmv_R1 --output plots/graph/html/universal/dcmv_R1.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/dcmv_R2 --output plots/graph/html/universal/dcmv_R2.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot the comparison graphs (HTML) ###
DSBplot-graph --input data/output/sense_db_R1 --node_comparison_colors "#cf191b" "#33a02c" --output plots/graph/html/universal/sense_db_R1.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/sense_db_R2 --node_comparison_colors "#cf191b" "#33a02c" --output plots/graph/html/universal/sense_db_R2.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/sense_dcmv_R1 --output plots/graph/html/universal/sense_dcmv_R1.html --node_comparison_colors "#cf191b" "#ffe669" --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/sense_dcmv_R2 --output plots/graph/html/universal/sense_dcmv_R2.html --node_comparison_colors "#cf191b" "#ffe669" --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot graphs with common Universal layout (HTML) ###
DSBplot-graph `
  --input data/output/db_R1 data/output/db_R2 data/output/sense_R1 data/output/sense_R2 data/output/dcmv_R1 data/output/dcmv_R2 data/output/sense_db_R1 data/output/sense_db_R2 data/output/sense_dcmv_R1 data/output/sense_dcmv_R2 `
  --output plots/graph/html/universal_common/db_R1.html plots/graph/html/universal_common/db_R2.html plots/graph/html/universal_common/sense_R1.html plots/graph/html/universal_common/sense_R2.html plots/graph/html/universal_common/dcmv_R1.html plots/graph/html/universal_common/dcmv_R2.html plots/graph/html/universal_common/sense_db_R1.html plots/graph/html/universal_common/sense_db_R2.html plots/graph/html/universal_common/sense_dcmv_R1.html plots/graph/html/universal_common/sense_dcmv_R2.html `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout kamada_layout --width 2400 --height 1800

### Plot graphs with common Kamada layout (HTML) ###
DSBplot-graph `
  --input data/output/db_R1 data/output/db_R2 data/output/sense_R1 data/output/sense_R2 data/output/dcmv_R1 data/output/dcmv_R2 data/output/sense_db_R1 data/output/sense_db_R2 data/output/sense_dcmv_R1 data/output/sense_dcmv_R2 `
  --output plots/graph/html/kamada_common/db_R1.html plots/graph/html/kamada_common/db_R2.html plots/graph/html/kamada_common/sense_R1.html plots/graph/html/kamada_common/sense_R2.html plots/graph/html/kamada_common/dcmv_R1.html plots/graph/html/kamada_common/dcmv_R2.html plots/graph/html/kamada_common/sense_db_R1.html plots/graph/html/kamada_common/sense_db_R2.html plots/graph/html/kamada_common/sense_dcmv_R1.html plots/graph/html/kamada_common/sense_dcmv_R2.html `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout kamada_layout --width 2400 --height 1800

### Plot graphs with common Radial layout (HTML) ###
DSBplot-graph `
  --input data/output/db_R1 data/output/db_R2 data/output/sense_R1 data/output/sense_R2 data/output/dcmv_R1 data/output/dcmv_R2 data/output/sense_db_R1 data/output/sense_db_R2 data/output/sense_dcmv_R1 data/output/sense_dcmv_R2 `
  --output plots/graph/html/radial_common/db_R1.html plots/graph/html/radial_common/db_R2.html plots/graph/html/radial_common/sense_R1.html plots/graph/html/radial_common/sense_R2.html plots/graph/html/radial_common/dcmv_R1.html plots/graph/html/radial_common/dcmv_R2.html plots/graph/html/radial_common/sense_db_R1.html plots/graph/html/radial_common/sense_db_R2.html plots/graph/html/radial_common/sense_dcmv_R1.html plots/graph/html/radial_common/sense_dcmv_R2.html `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --range_x -20 20 --range_y -20 20 --layout radial_layout --width 2400 --height 1800

### Arrange the histograms into a grid in a PPTX file ###
DSBplot-pptx --input `
  plots/histogram/png/db_R1_substitution.png `
  plots/histogram/png/db_R1_insertion.png `
  plots/histogram/png/db_R1_deletion.png `
  plots/histogram/png/sense_R1_substitution.png `
  plots/histogram/png/sense_R1_insertion.png `
  plots/histogram/png/sense_R1_deletion.png `
  plots/histogram/png/dcmv_R1_substitution.png `
  plots/histogram/png/dcmv_R1_insertion.png `
  plots/histogram/png/dcmv_R1_deletion.png `
  plots/histogram/png/db_R2_substitution.png `
  plots/histogram/png/db_R2_insertion.png `
  plots/histogram/png/db_R2_deletion.png `
  plots/histogram/png/sense_R2_substitution.png `
  plots/histogram/png/sense_R2_insertion.png `
  plots/histogram/png/sense_R2_deletion.png `
  plots/histogram/png/dcmv_R2_substitution.png `
  plots/histogram/png/dcmv_R2_insertion.png `
  plots/histogram/png/dcmv_R2_deletion.png `
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
  --output plots/pptx/histograms.pptx