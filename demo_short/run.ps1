$input_dir = "input/fastq"
$ref_seq_dir = "input/ref_seq"
$output_dir = "output"
$debug_dir = "debug"

### Preprocessing ###
DSBplot-preprocess --input $input_dir/Sense_R1_1.fq $input_dir/Sense_R1_2.fq $input_dir/Sense_R1_3.fq $input_dir/Sense_R1_4.fq --ref_seq_file $ref_seq_dir/2DSB_Sense_R1.fa --dsb_pos 67 --output $output_dir/Sense_R1 --label "Sense (R1)" --total_reads 3000 3000 3000 3000

# Example of running preprocessing stages separately
DSBplot-preprocess --input $input_dir/Sense_R1_1.fq $input_dir/Sense_R1_2.fq $input_dir/Sense_R1_3.fq $input_dir/Sense_R1_4.fq --ref_seq_file $ref_seq_dir/2DSB_Sense_R1.fa --output $output_dir/Sense_R1 --stages 0_align
DSBplot-preprocess --ref_seq_file $ref_seq_dir/2DSB_Sense_R1.fa --dsb_pos 67 --output $output_dir/Sense_R1 --total_reads 3000 3000 3000 3000 --stages 1_filter
DSBplot-preprocess --ref_seq_file $ref_seq_dir/2DSB_Sense_R1.fa --dsb_pos 67 --output $output_dir/Sense_R1 --label Sense_R1 --total_reads 3000 3000 3000 3000 --stages 2_window
DSBplot-preprocess --output $output_dir/Sense_R1 --stages 3_variation

### Plot the graphs (PNG) ###
DSBplot-graph --input $output_dir/Sense_R1 --output plots/graph/Sense_R1_universal.png --debug $debug_dir/universal --layout universal --title "Sense (R1) Universal Layout" --width 2400 --height 1800 --universal_y_axis_x_pos 0 --universal_x_axis_deletion_y_pos 0 --universal_x_axis_insertion_y_pos 0
DSBplot-graph --input $output_dir/Sense_R1 --output plots/graph/Sense_R1_kamada.png --debug $debug_dir/kamada --layout kamada --title "Sense (R1) Kamada-Kawaii Layout" --width 2400 --height 1800 --separate_components
DSBplot-graph --input $output_dir/Sense_R1 --output plots/graph/Sense_R1_radial.png --debug $debug_dir/radial --layout radial --title "Sense (R1) Radial Layout" --width 2400 --height 1800

### Plot the graphs (HTML) ###
DSBplot-graph --input $output_dir/Sense_R1 --output plots/graph/Sense_R1_universal.html --layout universal --title "Sense (R1) Universal Layout" --width 2400 --height 1800 --universal_y_axis_x_pos 0 --universal_x_axis_deletion_y_pos 0 --universal_x_axis_insertion_y_pos 0
DSBplot-graph --input $output_dir/Sense_R1 --output plots/graph/Sense_R1_kamada.html --layout kamada --title "Sense (R1) Kamada-Kawaii Layout" --width 2400 --height 1800 --separate_components
DSBplot-graph --input $output_dir/Sense_R1 --output plots/graph/Sense_R1_radial.html --layout radial --title "Sense (R1) Radial Layout" --width 2400 --height 1800

### Plot the 3D histograms (PNG) ###
DSBplot-histogram --input $output_dir/Sense_R1 --output plots/histogram/Sense_R1_substitution.png --title "Sense (R1) Substitutions" --color "#bfbfbf" --margin_top 200 --variation_type substitution --label_type relative
DSBplot-histogram --input $output_dir/Sense_R1 --output plots/histogram/Sense_R1_insertion.png --title "Sense (R1) Insertions" --color "#ffa500" --margin_top 200 --variation_type insertion --label_type relative
DSBplot-histogram --input $output_dir/Sense_R1 --output plots/histogram/Sense_R1_deletion.png --title "Sense (R1) Deletions" --color "#8080ff" --margin_top 200 --variation_type deletion --label_type relative
