### Directories. Run this first! ###
$input_dir = "input/fastq"
$ref_seq_dir = "input/ref_seq"
$output_dir = "output"
$debug_dir = "debug"

### Preprocessing ###
DSBplot-preprocess -i $input_dir/Sense_R1_1.fq $input_dir/Sense_R1_2.fq $input_dir/Sense_R1_3.fq $input_dir/Sense_R1_4.fq --ref $ref_seq_dir/2DSB_Sense_R1.fa --dsb 67 -o $output_dir/Sense_R1 --label "Sense (R1)" --reads 3000 3000 3000 3000

# Example of running preprocessing stages separately (R1)
# Some parameters are not required to be specified for a subsequent stage if they are used in the previous stage.
DSBplot-preprocess -i $input_dir/Sense_R1_1.fq $input_dir/Sense_R1_2.fq $input_dir/Sense_R1_3.fq $input_dir/Sense_R1_4.fq --ref $ref_seq_dir/2DSB_Sense_R1.fa -o $output_dir/Sense_R1 --stages 0
DSBplot-preprocess --ref $ref_seq_dir/2DSB_Sense_R1.fa --dsb 67 -o $output_dir/Sense_R1 --reads 3000 3000 3000 3000 --stages 1_filter
DSBplot-preprocess -o $output_dir/Sense_R1 --stages 2_window
DSBplot-preprocess -o $output_dir/Sense_R1 --stages 3_variation
DSBplot-preprocess -o $output_dir/Sense_R1 --stages 4_info

### Plot the graphs (PNG) ###
DSBplot-graph -i $output_dir/Sense_R1 -o plots/graph/Sense_R1_universal.png --debug $debug_dir/universal --layout universal --title "Sense (R1) Universal Layout" --width 2400 --height 1800 --ul_yax_x 0 --ul_xax_del_y 0 --ul_xax_ins_y 0
DSBplot-graph -i $output_dir/Sense_R1 -o plots/graph/Sense_R1_kamada.png --debug $debug_dir/kamada --layout kamada --title "Sense (R1) Kamada-Kawaii Layout" --width 2400 --height 1800 --sep
DSBplot-graph -i $output_dir/Sense_R1 -o plots/graph/Sense_R1_radial.png --debug $debug_dir/radial --layout radial --title "Sense (R1) Radial Layout" --width 2400 --height 1800

### Plot the graphs (HTML) ###
DSBplot-graph -i $output_dir/Sense_R1 -o plots/graph/Sense_R1_universal.html --layout universal --title "Sense (R1) Universal Layout" --width 2400 --height 1800 --ul_yax_x 0 --ul_xax_del_y 0 --ul_xax_ins_y 0
DSBplot-graph -i $output_dir/Sense_R1 -o plots/graph/Sense_R1_kamada.html --layout kamada --title "Sense (R1) Kamada-Kawaii Layout" --width 2400 --height 1800 --sep
DSBplot-graph -i $output_dir/Sense_R1 -o plots/graph/Sense_R1_radial.html --layout radial --title "Sense (R1) Radial Layout" --width 2400 --height 1800

### Plot the 3D histograms (PNG) ###
DSBplot-histogram -i $output_dir/Sense_R1 -o plots/histogram/Sense_R1_sub.png --title "Sense (R1) Substitutions" --color "#bfbfbf" --mar_t 200 --var sub --xax rel
DSBplot-histogram -i $output_dir/Sense_R1 -o plots/histogram/Sense_R1_ins.png --title "Sense (R1) Insertions" --color "#ffa500" --mar_t 200 --var ins --xax rel
DSBplot-histogram -i $output_dir/Sense_R1 -o plots/histogram/Sense_R1_del.png --title "Sense (R1) Deletions" --color "#8080ff" --mar_t 200 --var del --xax rel
