### Preprocessing ###
DSBplot-preprocess --input data/input/fastq/BranchD_R1_1.fq data/input/fastq/BranchD_R1_2.fq data/input/fastq/BranchD_R1_3.fq data/input/fastq/BranchD_R1_4.fq --ref_seq_file data/input/ref_seq/2DSB_BranchD_R1.fa --dsb_pos 67 --output data/output/BranchD_R1 --label "BranchΔ" --total_reads 3000 3000 3000 3000
DSBplot-preprocess --input data/input/fastq/BranchD_R2_1.fq data/input/fastq/BranchD_R2_2.fq data/input/fastq/BranchD_R2_3.fq data/input/fastq/BranchD_R2_4.fq --ref_seq_file data/input/ref_seq/2DSB_BranchD_R2.fa --dsb_pos 46 --output data/output/BranchD_R2 --label "BranchΔ" --total_reads 3000 3000 3000 3000
DSBplot-preprocess --input data/input/fastq/Sense_R1_1.fq data/input/fastq/Sense_R1_2.fq data/input/fastq/Sense_R1_3.fq data/input/fastq/Sense_R1_4.fq --ref_seq_file data/input/ref_seq/2DSB_Sense_R1.fa --dsb_pos 67 --output data/output/Sense_R1 --label "Sense" --total_reads 3000 3000 3000 3000
DSBplot-preprocess --input data/input/fastq/Sense_R2_1.fq data/input/fastq/Sense_R2_2.fq data/input/fastq/Sense_R2_3.fq data/input/fastq/Sense_R2_4.fq --ref_seq_file data/input/ref_seq/2DSB_Sense_R2.fa --dsb_pos 46 --output data/output/Sense_R2 --label "Sense" --total_reads 3000 3000 3000 3000
DSBplot-preprocess --input data/input/fastq/pCMVD_R1_1.fq data/input/fastq/pCMVD_R1_2.fq data/input/fastq/pCMVD_R1_3.fq data/input/fastq/pCMVD_R1_4.fq --ref_seq_file data/input/ref_seq/2DSB_pCMVD_R1.fa --dsb_pos 67 --output data/output/pCMVD_R1 --label "pCMVΔ" --total_reads 3000 3000 3000 3000
DSBplot-preprocess --input data/input/fastq/pCMVD_R2_1.fq data/input/fastq/pCMVD_R2_2.fq data/input/fastq/pCMVD_R2_3.fq data/input/fastq/pCMVD_R2_4.fq --ref_seq_file data/input/ref_seq/2DSB_pCMVD_R2.fa --dsb_pos 46 --output data/output/pCMVD_R2 --label "pCMVΔ" --total_reads 3000 3000 3000 3000

# Example of running preprocessing stages separately for one experiment
DSBplot-preprocess --input data/input/fastq/Sense_R1_1.fq data/input/fastq/Sense_R1_2.fq data/input/fastq/Sense_R1_3.fq data/input/fastq/Sense_R1_4.fq --ref_seq_file data/input/ref_seq/2DSB_Sense_R1.fa --output data/output/Sense_R1 --stages 0_align
DSBplot-preprocess --ref_seq_file data/input/ref_seq/2DSB_Sense_R1.fa --dsb_pos 67 --output data/output/Sense_R1 --stages 1_filter
DSBplot-preprocess --output data/output/Sense_R1 --stages 2_combine
DSBplot-preprocess --ref_seq_file data/input/ref_seq/2DSB_Sense_R1.fa --dsb_pos 67 --output data/output/Sense_R1 --label Sense_R1 --total_reads 3000 3000 3000 3000 --stages 3_window
DSBplot-preprocess --output data/output/Sense_R1 --stages 4_graph
DSBplot-preprocess --output data/output/Sense_R1 --stages 5_histogram

### Comparison data ###
DSBplot-comparison --input data/output/Sense_R1 data/output/BranchD_R1 --output data/output/Sense_BranchD_R1
DSBplot-comparison --input data/output/Sense_R2 data/output/BranchD_R2 --output data/output/Sense_BranchD_R2
DSBplot-comparison --input data/output/Sense_R1 data/output/pCMVD_R1 --output data/output/Sense_pCMVD_R1
DSBplot-comparison --input data/output/Sense_R2 data/output/pCMVD_R2 --output data/output/Sense_pCMVD_R2

### Plot the graphs (PNG) ###
DSBplot-graph --input data/output/BranchD_R1 --output plots/graph/png/universal/BranchD_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5  --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/BranchD_R2 --output plots/graph/png/universal/BranchD_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/Sense_R1 --output plots/graph/png/universal/Sense_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/Sense_R2 --output plots/graph/png/universal/Sense_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/pCMVD_R1 --output plots/graph/png/universal/pCMVD_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/pCMVD_R2 --output plots/graph/png/universal/pCMVD_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot the comparison graphs (PNG) ###
DSBplot-graph --input data/output/Sense_BranchD_R1 --node_comparison_colors "#cf191b" "#33a02c" --output plots/graph/png/universal/Sense_BranchD_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/Sense_BranchD_R2 --node_comparison_colors "#cf191b" "#33a02c" --output plots/graph/png/universal/Sense_BranchD_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/Sense_pCMVD_R1 --output plots/graph/png/universal/Sense_pCMVD_R1.png --node_comparison_colors "#cf191b" "#ffe669" --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/Sense_pCMVD_R2 --output plots/graph/png/universal/Sense_pCMVD_R2.png --node_comparison_colors "#cf191b" "#ffe669" --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot the 3D histograms (PNG) ###
DSBplot-histogram --input data/output/BranchD_R1 --output plots/histogram/png/BranchD_R1_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/BranchD_R1 --output plots/histogram/png/BranchD_R1_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/BranchD_R1 --output plots/histogram/png/BranchD_R1_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/BranchD_R2 --output plots/histogram/png/BranchD_R2_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/BranchD_R2 --output plots/histogram/png/BranchD_R2_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/BranchD_R2 --output plots/histogram/png/BranchD_R2_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/Sense_R1 --output plots/histogram/png/Sense_R1_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/Sense_R1 --output plots/histogram/png/Sense_R1_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/Sense_R1 --output plots/histogram/png/Sense_R1_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/Sense_R2 --output plots/histogram/png/Sense_R2_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/Sense_R2 --output plots/histogram/png/Sense_R2_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/Sense_R2 --output plots/histogram/png/Sense_R2_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/pCMVD_R1 --output plots/histogram/png/pCMVD_R1_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/pCMVD_R1 --output plots/histogram/png/pCMVD_R1_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/pCMVD_R1 --output plots/histogram/png/pCMVD_R1_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/pCMVD_R2 --output plots/histogram/png/pCMVD_R2_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/pCMVD_R2 --output plots/histogram/png/pCMVD_R2_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/pCMVD_R2 --output plots/histogram/png/pCMVD_R2_deletion.png --color "#8080ff" --variation_type deletion --label_type relative

### Plot graphs with common Universal layout (PNG) ###
DSBplot-graph `
  --input data/output/BranchD_R1 data/output/BranchD_R2 data/output/Sense_R1 data/output/Sense_R2 data/output/pCMVD_R1 data/output/pCMVD_R2 data/output/Sense_BranchD_R1 data/output/Sense_BranchD_R2 data/output/Sense_pCMVD_R1 data/output/Sense_pCMVD_R2 `
  --output plots/graph/png/universal_common/BranchD_R1.png plots/graph/png/universal_common/BranchD_R2.png plots/graph/png/universal_common/Sense_R1.png plots/graph/png/universal_common/Sense_R2.png plots/graph/png/universal_common/pCMVD_R1.png plots/graph/png/universal_common/pCMVD_R2.png plots/graph/png/universal_common/Sense_BranchD_R1.png plots/graph/png/universal_common/Sense_BranchD_R2.png plots/graph/png/universal_common/Sense_pCMVD_R1.png plots/graph/png/universal_common/Sense_pCMVD_R2.png `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout kamada_layout --width 2400 --height 1800

### Plot graphs with common Kamada layout (PNG) ###
DSBplot-graph `
  --input data/output/BranchD_R1 data/output/BranchD_R2 data/output/Sense_R1 data/output/Sense_R2 data/output/pCMVD_R1 data/output/pCMVD_R2 data/output/Sense_BranchD_R1 data/output/Sense_BranchD_R2 data/output/Sense_pCMVD_R1 data/output/Sense_pCMVD_R2 `
  --output plots/graph/png/kamada_common/BranchD_R1.png plots/graph/png/kamada_common/BranchD_R2.png plots/graph/png/kamada_common/Sense_R1.png plots/graph/png/kamada_common/Sense_R2.png plots/graph/png/kamada_common/pCMVD_R1.png plots/graph/png/kamada_common/pCMVD_R2.png plots/graph/png/kamada_common/Sense_BranchD_R1.png plots/graph/png/kamada_common/Sense_BranchD_R2.png plots/graph/png/kamada_common/Sense_pCMVD_R1.png plots/graph/png/kamada_common/Sense_pCMVD_R2.png `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout kamada_layout --width 2400 --height 1800

### Plot graphs with common Radial layout (PNG) ###
DSBplot-graph `
  --input data/output/BranchD_R1 data/output/BranchD_R2 data/output/Sense_R1 data/output/Sense_R2 data/output/pCMVD_R1 data/output/pCMVD_R2 data/output/Sense_BranchD_R1 data/output/Sense_BranchD_R2 data/output/Sense_pCMVD_R1 data/output/Sense_pCMVD_R2 `
  --output plots/graph/png/radial_common/BranchD_R1.png plots/graph/png/radial_common/BranchD_R2.png plots/graph/png/radial_common/Sense_R1.png plots/graph/png/radial_common/Sense_R2.png plots/graph/png/radial_common/pCMVD_R1.png plots/graph/png/radial_common/pCMVD_R2.png plots/graph/png/radial_common/Sense_BranchD_R1.png plots/graph/png/radial_common/Sense_BranchD_R2.png plots/graph/png/radial_common/Sense_pCMVD_R1.png plots/graph/png/radial_common/Sense_pCMVD_R2.png `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --range_x -20 20 --range_y -20 20 --layout radial_layout --width 2400 --height 1800

### Plot the graphs (PDF) ###
DSBplot-graph --input data/output/BranchD_R1 --output plots/graph/pdf/universal/BranchD_R1.pdf --title "BranchΔ (R1)" --legends size variation_type --font_size_scale 3 --legend_x_shift_px 100 --legend_y_shift_px 0 --margin_right_px 750 --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5  --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/BranchD_R2 --output plots/graph/pdf/universal/BranchD_R2.pdf --title "BranchΔ (R2)" --legends size variation_type --font_size_scale 3 --legend_x_shift_px 100 --legend_y_shift_px 0 --margin_right_px 750 --layout universal_layout --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/Sense_R1 --output plots/graph/pdf/universal/Sense_R1.pdf --title "Sense (R1)" --legends size variation_type --font_size_scale 3 --legend_x_shift_px 100 --legend_y_shift_px 0 --margin_right_px 750 --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/Sense_R2 --output plots/graph/pdf/universal/Sense_R2.pdf --title "Sense (R2)" --legends size variation_type --font_size_scale 3 --legend_x_shift_px 100 --legend_y_shift_px 0 --margin_right_px 750 --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/pCMVD_R1 --output plots/graph/pdf/universal/pCMVD_R1.pdf --title "pCMVΔ (R1)" --legends size variation_type --font_size_scale 3 --legend_x_shift_px 100 --legend_y_shift_px 0 --margin_right_px 750 --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/pCMVD_R2 --output plots/graph/pdf/universal/pCMVD_R2.pdf --title "pCMVΔ (R2)" --legends size variation_type --font_size_scale 3 --legend_x_shift_px 100 --legend_y_shift_px 0 --margin_right_px 750 --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot the comparison graphs (PDF) ###
DSBplot-graph --input data/output/Sense_BranchD_R1 --node_comparison_colors "#cf191b" "#33a02c" --output plots/graph/pdf/universal/Sense_BranchD_R1.pdf --title "Sense / BranchΔ (R1)" --legend_colorbar_scale 3 --font_size_scale 3 --legend_x_shift_px 100 --legend_y_shift_px -100 --margin_right_px 750 --legends freq_ratio_continuous --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/Sense_BranchD_R2 --node_comparison_colors "#cf191b" "#33a02c" --output plots/graph/pdf/universal/Sense_BranchD_R2.pdf --title "Sense / BranchΔ (R2)" --legend_colorbar_scale 3 --font_size_scale 3 --legend_x_shift_px 100 --legend_y_shift_px -100 --margin_right_px 750 --legends freq_ratio_continuous --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/Sense_pCMVD_R1 --output plots/graph/pdf/universal/Sense_pCMVD_R1.pdf --node_comparison_colors "#cf191b" "#ffe669" --title "Sense / pCMVΔ (R1)" --legend_colorbar_scale 3 --font_size_scale 3 --legend_x_shift_px 100 --legend_y_shift_px -100 --margin_right_px 750 --legends freq_ratio_continuous --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/Sense_pCMVD_R2 --output plots/graph/pdf/universal/Sense_pCMVD_R2.pdf --node_comparison_colors "#cf191b" "#ffe669" --title "Sense / pCMVΔ (R2)" --legend_colorbar_scale 3 --font_size_scale 3 --legend_x_shift_px 100 --legend_y_shift_px -100 --margin_right_px 750 --legends freq_ratio_continuous --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot the 3D histograms (PDF) ###
DSBplot-histogram --input data/output/BranchD_R1 --output plots/histogram/pdf/BranchD_R1_substitution.pdf --title "BranchΔ (R1)`nSubstitution" --margin_top 250 --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/BranchD_R1 --output plots/histogram/pdf/BranchD_R1_insertion.pdf --title "BranchΔ (R1)`nInsertion" --margin_top 250 --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/BranchD_R1 --output plots/histogram/pdf/BranchD_R1_deletion.pdf --title "BranchΔ (R1)`nDeletion" --margin_top 250 --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/BranchD_R2 --output plots/histogram/pdf/BranchD_R2_substitution.pdf --title "BranchΔ (R2)`nSubstitution" --margin_top 250 --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/BranchD_R2 --output plots/histogram/pdf/BranchD_R2_insertion.pdf --title "BranchΔ (R2)`nInsertion" --margin_top 250 --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/BranchD_R2 --output plots/histogram/pdf/BranchD_R2_deletion.pdf --title "BranchΔ (R2)`nDeletion" --margin_top 250 --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/Sense_R1 --output plots/histogram/pdf/Sense_R1_substitution.pdf --title "Sense (R1)`nSubstitution" --margin_top 250 --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/Sense_R1 --output plots/histogram/pdf/Sense_R1_insertion.pdf --title "Sense (R1)`nInsertion" --margin_top 250 --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/Sense_R1 --output plots/histogram/pdf/Sense_R1_deletion.pdf --title "Sense (R1)`nDeletion" --margin_top 250 --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/Sense_R2 --output plots/histogram/pdf/Sense_R2_substitution.pdf --title "Sense (R2)`nSubstitution" --margin_top 250 --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/Sense_R2 --output plots/histogram/pdf/Sense_R2_insertion.pdf --title "Sense (R2)`nInsertion" --margin_top 250 --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/Sense_R2 --output plots/histogram/pdf/Sense_R2_deletion.pdf --title "Sense (R2)`nDeletion" --margin_top 250 --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/pCMVD_R1 --output plots/histogram/pdf/pCMVD_R1_substitution.pdf --title "pCMVΔ (R1)`nSubstitution" --margin_top 250 --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/pCMVD_R1 --output plots/histogram/pdf/pCMVD_R1_insertion.pdf --title "pCMVΔ (R1)`nInsertion" --margin_top 250 --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/pCMVD_R1 --output plots/histogram/pdf/pCMVD_R1_deletion.pdf --title "pCMVΔ (R1)`nDeletion" --margin_top 250 --color "#8080ff" --variation_type deletion --label_type relative
DSBplot-histogram --input data/output/pCMVD_R2 --output plots/histogram/pdf/pCMVD_R2_substitution.pdf --title "pCMVΔ (R2)`nSubstitution" --margin_top 250 --color "#bfbfbf" --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/pCMVD_R2 --output plots/histogram/pdf/pCMVD_R2_insertion.pdf --title "pCMVΔ (R2)`nInsertion" --margin_top 250 --color "#ffa500" --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/pCMVD_R2 --output plots/histogram/pdf/pCMVD_R2_deletion.pdf --title "pCMVΔ (R2)`nDeletion" --margin_top 250 --color "#8080ff" --variation_type deletion --label_type relative

### Plot graphs with common Universal layout (PDF) ###
DSBplot-graph `
  --input data/output/BranchD_R1 data/output/BranchD_R2 data/output/Sense_R1 data/output/Sense_R2 data/output/pCMVD_R1 data/output/pCMVD_R2 data/output/Sense_BranchD_R1 data/output/Sense_BranchD_R2 data/output/Sense_pCMVD_R1 data/output/Sense_pCMVD_R2 `
  --title "BranchΔ (R1)" "BranchΔ (R2)" "Sense (R1)" "pCMVΔ (R1)" "pCMVΔ (R2)" "pCMVΔ (R2)" "Sense / BranchΔ (R1)" "Sense / BranchΔ (R2)" "Sense / pCMVΔ (R1)" "Sense / pCMVΔ (R2)" `
  --output plots/graph/pdf/universal_common/BranchD_R1.pdf plots/graph/pdf/universal_common/BranchD_R2.pdf plots/graph/pdf/universal_common/Sense_R1.pdf plots/graph/pdf/universal_common/Sense_R2.pdf plots/graph/pdf/universal_common/pCMVD_R1.pdf plots/graph/pdf/universal_common/pCMVD_R2.pdf plots/graph/pdf/universal_common/Sense_BranchD_R1.pdf plots/graph/pdf/universal_common/Sense_BranchD_R2.pdf plots/graph/pdf/universal_common/Sense_pCMVD_R1.pdf plots/graph/pdf/universal_common/Sense_pCMVD_R2.pdf `
  --range_x -12 13 --range_y -23 20 `
  --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 `
  --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_deletion_max_tick 19 `
  --legends freq_ratio_continuous variation_type `
  --legend_x_shift_px 100 --legend_y_shift_px -100 `
  --margin_top_px 300 --margin_right_px 900 --margin_left_px 0 --margin_bottom_px 0 `
  --font_size_scale 3 --legend_colorbar_scale 3 `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout universal_layout --width 2400 --height 1800

### Plot graphs with common Kamada layout (PDF) ###
DSBplot-graph `
  --input data/output/BranchD_R1 data/output/BranchD_R2 data/output/Sense_R1 data/output/Sense_R2 data/output/pCMVD_R1 data/output/pCMVD_R2 data/output/Sense_BranchD_R1 data/output/Sense_BranchD_R2 data/output/Sense_pCMVD_R1 data/output/Sense_pCMVD_R2 `
  --title "BranchΔ (R1)" "BranchΔ (R2)" "Sense (R1)" "pCMVΔ (R1)" "pCMVΔ (R2)" "pCMVΔ (R2)" "Sense / BranchΔ (R1)" "Sense / BranchΔ (R2)" "Sense / pCMVΔ (R1)" "Sense / pCMVΔ (R2)" `
  --output plots/graph/pdf/kamada_common/BranchD_R1.pdf plots/graph/pdf/kamada_common/BranchD_R2.pdf plots/graph/pdf/kamada_common/Sense_R1.pdf plots/graph/pdf/kamada_common/Sense_R2.pdf plots/graph/pdf/kamada_common/pCMVD_R1.pdf plots/graph/pdf/kamada_common/pCMVD_R2.pdf plots/graph/pdf/kamada_common/Sense_BranchD_R1.pdf plots/graph/pdf/kamada_common/Sense_BranchD_R2.pdf plots/graph/pdf/kamada_common/Sense_pCMVD_R1.pdf plots/graph/pdf/kamada_common/Sense_pCMVD_R2.pdf `
  --range_x 0.15 0.85 --range_y 0.15 0.85 `
  --legends freq_ratio_continuous variation_type `
  --legend_x_shift_px 100 --legend_y_shift_px -100 `
  --font_size_scale 3 --legend_colorbar_scale 3 `
  --margin_top_px 300 --margin_right_px 900 --margin_left_px 0 --margin_bottom_px 0 `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout kamada_layout --width 2400 --height 1800

### Plot graphs with common Radial layout (PDF) ###
DSBplot-graph `
  --input data/output/BranchD_R1 data/output/BranchD_R2 data/output/Sense_R1 data/output/Sense_R2 data/output/pCMVD_R1 data/output/pCMVD_R2 data/output/Sense_BranchD_R1 data/output/Sense_BranchD_R2 data/output/Sense_pCMVD_R1 data/output/Sense_pCMVD_R2 `
  --title "BranchΔ (R1)" "BranchΔ (R2)" "Sense (R1)" "pCMVΔ (R1)" "pCMVΔ (R2)" "pCMVΔ (R2)" "Sense / BranchΔ (R1)" "Sense / BranchΔ (R2)" "Sense / pCMVΔ (R1)" "Sense / pCMVΔ (R2)" `
  --output plots/graph/pdf/radial_common/BranchD_R1.pdf plots/graph/pdf/radial_common/BranchD_R2.pdf plots/graph/pdf/radial_common/Sense_R1.pdf plots/graph/pdf/radial_common/Sense_R2.pdf plots/graph/pdf/radial_common/pCMVD_R1.pdf plots/graph/pdf/radial_common/pCMVD_R2.pdf plots/graph/pdf/radial_common/Sense_BranchD_R1.pdf plots/graph/pdf/radial_common/Sense_BranchD_R2.pdf plots/graph/pdf/radial_common/Sense_pCMVD_R1.pdf plots/graph/pdf/radial_common/Sense_pCMVD_R2.pdf `
  --range_x -18 18 --range_y -16 7 `
  --legends freq_ratio_continuous variation_type `
  --legend_x_shift_px 100 --legend_y_shift_px -100 `
  --font_size_scale 3 --legend_colorbar_scale 3 `
  --margin_top_px 300 --margin_right_px 900 --margin_left_px 0 --margin_bottom_px 0 `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout radial_layout --width 2400 --height 1800

### Plot the graphs (HTML) ###
DSBplot-graph --input data/output/BranchD_R1 --output plots/graph/html/universal/BranchD_R1.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5  --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/BranchD_R2 --output plots/graph/html/universal/BranchD_R2.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/Sense_R1 --output plots/graph/html/universal/Sense_R1.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/Sense_R2 --output plots/graph/html/universal/Sense_R2.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/pCMVD_R1 --output plots/graph/html/universal/pCMVD_R1.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/pCMVD_R2 --output plots/graph/html/universal/pCMVD_R2.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot the comparison graphs (HTML) ###
DSBplot-graph --input data/output/Sense_BranchD_R1 --node_comparison_colors "#cf191b" "#33a02c" --output plots/graph/html/universal/Sense_BranchD_R1.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/Sense_BranchD_R2 --node_comparison_colors "#cf191b" "#33a02c" --output plots/graph/html/universal/Sense_BranchD_R2.html --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/Sense_pCMVD_R1 --output plots/graph/html/universal/Sense_pCMVD_R1.html --node_comparison_colors "#cf191b" "#ffe669" --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/Sense_pCMVD_R2 --output plots/graph/html/universal/Sense_pCMVD_R2.html --node_comparison_colors "#cf191b" "#ffe669" --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot graphs with common Universal layout (HTML) ###
DSBplot-graph `
  --input data/output/BranchD_R1 data/output/BranchD_R2 data/output/Sense_R1 data/output/Sense_R2 data/output/pCMVD_R1 data/output/pCMVD_R2 data/output/Sense_BranchD_R1 data/output/Sense_BranchD_R2 data/output/Sense_pCMVD_R1 data/output/Sense_pCMVD_R2 `
  --output plots/graph/html/universal_common/BranchD_R1.html plots/graph/html/universal_common/BranchD_R2.html plots/graph/html/universal_common/Sense_R1.html plots/graph/html/universal_common/Sense_R2.html plots/graph/html/universal_common/pCMVD_R1.html plots/graph/html/universal_common/pCMVD_R2.html plots/graph/html/universal_common/Sense_BranchD_R1.html plots/graph/html/universal_common/Sense_BranchD_R2.html plots/graph/html/universal_common/Sense_pCMVD_R1.html plots/graph/html/universal_common/Sense_pCMVD_R2.html `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout kamada_layout --width 2400 --height 1800

### Plot graphs with common Kamada layout (HTML) ###
DSBplot-graph `
  --input data/output/BranchD_R1 data/output/BranchD_R2 data/output/Sense_R1 data/output/Sense_R2 data/output/pCMVD_R1 data/output/pCMVD_R2 data/output/Sense_BranchD_R1 data/output/Sense_BranchD_R2 data/output/Sense_pCMVD_R1 data/output/Sense_pCMVD_R2 `
  --output plots/graph/html/kamada_common/BranchD_R1.html plots/graph/html/kamada_common/BranchD_R2.html plots/graph/html/kamada_common/Sense_R1.html plots/graph/html/kamada_common/Sense_R2.html plots/graph/html/kamada_common/pCMVD_R1.html plots/graph/html/kamada_common/pCMVD_R2.html plots/graph/html/kamada_common/Sense_BranchD_R1.html plots/graph/html/kamada_common/Sense_BranchD_R2.html plots/graph/html/kamada_common/Sense_pCMVD_R1.html plots/graph/html/kamada_common/Sense_pCMVD_R2.html `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout kamada_layout --width 2400 --height 1800

### Plot graphs with common Radial layout (HTML) ###
DSBplot-graph `
  --input data/output/BranchD_R1 data/output/BranchD_R2 data/output/Sense_R1 data/output/Sense_R2 data/output/pCMVD_R1 data/output/pCMVD_R2 data/output/Sense_BranchD_R1 data/output/Sense_BranchD_R2 data/output/Sense_pCMVD_R1 data/output/Sense_pCMVD_R2 `
  --output plots/graph/html/radial_common/BranchD_R1.html plots/graph/html/radial_common/BranchD_R2.html plots/graph/html/radial_common/Sense_R1.html plots/graph/html/radial_common/Sense_R2.html plots/graph/html/radial_common/pCMVD_R1.html plots/graph/html/radial_common/pCMVD_R2.html plots/graph/html/radial_common/Sense_BranchD_R1.html plots/graph/html/radial_common/Sense_BranchD_R2.html plots/graph/html/radial_common/Sense_pCMVD_R1.html plots/graph/html/radial_common/Sense_pCMVD_R2.html `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --range_x -20 20 --range_y -20 20 --layout radial_layout --width 2400 --height 1800
