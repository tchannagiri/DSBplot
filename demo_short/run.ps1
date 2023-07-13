### Preprocessing ###
DSBplot-preprocess --input data/input/fastq/Sense_R1_1.fq data/input/fastq/Sense_R1_2.fq data/input/fastq/Sense_R1_3.fq data/input/fastq/Sense_R1_4.fq --ref_seq_file data/input/ref_seq/2DSB_Sense_R1.fa --dsb_pos 67 --output data/output/Sense_R1 --label "Sense" --total_reads 3000 3000 3000 3000

# Example of running preprocessing stages separately
DSBplot-preprocess --input data/input/fastq/Sense_R1_1.fq data/input/fastq/Sense_R1_2.fq data/input/fastq/Sense_R1_3.fq data/input/fastq/Sense_R1_4.fq --ref_seq_file data/input/ref_seq/2DSB_Sense_R1.fa --output data/output/Sense_R1 --stages 0_align
DSBplot-preprocess --ref_seq_file data/input/ref_seq/2DSB_Sense_R1.fa --dsb_pos 67 --output data/output/Sense_R1 --stages 1_filter
DSBplot-preprocess --output data/output/Sense_R1 --stages 2_combine
DSBplot-preprocess --ref_seq_file data/input/ref_seq/2DSB_Sense_R1.fa --dsb_pos 67 --output data/output/Sense_R1 --label Sense_R1 --total_reads 3000 3000 3000 3000 --stages 3_window
DSBplot-preprocess --output data/output/Sense_R1 --stages 4_graph
DSBplot-preprocess --output data/output/Sense_R1 --stages 5_histogram

# TODO: DELETE LINES EVENTUALLY
DSBplot-preprocess --input data/input/fastq/Sense_R1_1.fq data/input/fastq/Sense_R1_2.fq data/input/fastq/Sense_R1_3.fq data/input/fastq/Sense_R1_4.fq --ref_seq_file data/input/ref_seq/2DSB_Sense_R1.fa --output data/output/Sense_R1 --stages 0_align
DSBplot-preprocess --ref_seq_file data/input/ref_seq/2DSB_Sense_R1.fa --dsb_pos 67 --output data/output/Sense_R1 --stages 1_filter --total_reads 3000 3000 3000 3000
DSBplot-preprocess --ref_seq_file data/input/ref_seq/2DSB_Sense_R1.fa --dsb_pos 67 --output data/output/Sense_R1 --label Sense_R1 --total_reads 3000 3000 3000 3000 --stages 2_window
DSBplot-preprocess --output data/output/Sense_R1 --stages 3_variation

### Plot the graphs (PNG) ###
DSBplot-graph --input data/output/Sense_R1 --output plots/graph/Sense_R1_universal.png --layout universal_layout --title "Sense (R1) Universal Layout" --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/Sense_R1 --output plots/graph/Sense_R1_kamada.png --layout kamada_layout --title "Sense (R1) Kamada-Kawaii Layout" --width 2400 --height 1800
DSBplot-graph --input data/output/Sense_R1 --output plots/graph/Sense_R1_radial.png --layout radial_layout --title "Sense (R1) Radial Layout" --width 2400 --height 1800

# TEMP FOR TESTING
DSBplot-graph --input data/output/Sense_R1::data/output/Sense_R1 --output plots/temp.png --layout universal_layout --title "Sense (R1) Universal Layout" --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

### Plot the graphs (HTML) ###
DSBplot-graph --input data/output/Sense_R1 --output plots/graph/Sense_R1_universal.html --layout universal_layout --title "Sense (R1) Universal Layout" --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
DSBplot-graph --input data/output/Sense_R1 --output plots/graph/Sense_R1_kamada.html --layout kamada_layout --title "Sense (R1) Kamada-Kawaii Layout" --width 2400 --height 1800
DSBplot-graph --input data/output/Sense_R1 --output plots/graph/Sense_R1_radial.html --layout radial_layout --title "Sense (R1) Radial Layout" --width 2400 --height 1800

### Plot the 3D histograms (PNG) ###
DSBplot-histogram --input data/output/Sense_R1 --output plots/histogram/Sense_R1_substitution.png --title "Sense (R1) Substitutions" --color "#bfbfbf" --margin_top 200 --variation_type substitution --label_type relative
DSBplot-histogram --input data/output/Sense_R1 --output plots/histogram/Sense_R1_insertion.png --title "Sense (R1) Insertions" --color "#ffa500" --margin_top 200 --variation_type insertion --label_type relative
DSBplot-histogram --input data/output/Sense_R1 --output plots/histogram/Sense_R1_deletion.png --title "Sense (R1) Deletions" --color "#8080ff" --margin_top 200 --variation_type deletion --label_type relative
