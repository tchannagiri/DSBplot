### Define constants ###
$fastq_dir = "input/fastq"
$process_dir = "output"
$plots_dir = "plots"

$constructs = @(
  "Sense_R1",
  "Sense_R2",
  "BranchD_R1",
  "BranchD_R2",
  "pCMVD_R1",
  "pCMVD_R2"
)

$labels = @{
  "Sense_R1" = "Sense";
  "Sense_R2" = "Sense";
  "BranchD_R1" = "BranchΔ";
  "BranchD_R2" = "BranchΔ";
  "pCMVD_R1" = "pCMVΔ";
  "pCMVD_R2" = "pCMVΔ";
}

$dsb_pos = 67

$reverse_complement = @{
  "Sense_R1" = 0;
  "Sense_R2" = 1;
  "BranchD_R1" = 0;
  "BranchD_R2" = 1;
  "pCMVD_R1" = 0;
  "pCMVD_R2" = 1;
}

$variation_color = @{
  "sub" = "#bfbfbf";
  "ins" = "#ffa500";
  "del" = "#8080ff";
}

$ref_seq_file = @{
  "Sense_R1" = "input/ref_seq/2DSB_Sense_R1.fa";
  "Sense_R2" = "input/ref_seq/2DSB_Sense_R1.fa";
  "BranchD_R1" = "input/ref_seq/2DSB_BranchD_R1.fa";
  "BranchD_R2" = "input/ref_seq/2DSB_BranchD_R1.fa";
  "pCMVD_R1" = "input/ref_seq/2DSB_pCMVD_R1.fa";
  "pCMVD_R2" = "input/ref_seq/2DSB_pCMVD_R1.fa";
}

# PNG and PDF generation for graphs requires kaleido to be installed.
# Please see installation instructions in the README.
$graph_exts = @("png", "pdf", "html")
$graph_layouts = @("universal", "kamada", "radial")
# Whether or not to separate connected components of the graphs.
# Needed for the Kamada-Kawaii because otherwise the structure is muddled.
$sep = @{"universal" = 0; "kamada" = 1; "radial" = 0}

$histogram_exts = @("png", "pdf")
$histograms_var_types = @("sub", "ins", "del")

$constructs_plot = @(
  @("Sense_R1"),
  @("Sense_R2"),
  @("BranchD_R1"),
  @("BranchD_R2"),
  @("pCMVD_R1"),
  @("pCMVD_R2"),
  @("Sense_BranchD_R1", "Sense_R1", "BranchD_R1"),
  @("Sense_BranchD_R2", "Sense_R2", "BranchD_R2"),
  @("Sense_pCMVD_R1", "Sense_R1", "pCMVD_R1"),
  @("Sense_pCMVD_R2", "Sense_R2", "pCMVD_R2")
)

### Processing ###
foreach ($con in $constructs) {
  DSBplot-process -i $fastq_dir/${con}_1.fq $fastq_dir/${con}_2.fq $fastq_dir/${con}_3.fq $fastq_dir/${con}_4.fq -o $process_dir/$con --ref $ref_seq_file[$con] --label $labels[$con] --reads 3000 3000 3000 3000 --dsb $dsb_pos --rc $reverse_complement[$con]
}

# Example of running processing stages separately for one experiment
DSBplot-process -i input/fastq/Sense_R1_1.fq input/fastq/Sense_R1_2.fq input/fastq/Sense_R1_3.fq input/fastq/Sense_R1_4.fq --ref input/ref_seq/2DSB_Sense_R1.fa -o $process_dir/Sense_R1 --stages 0_align
DSBplot-process --ref input/ref_seq/2DSB_Sense_R1.fa --dsb 67 -o $process_dir/Sense_R1 --stages 1_filter
DSBplot-process -o $process_dir/Sense_R1 --stages 2_window
DSBplot-process -o $process_dir/Sense_R1 --stages 3_variation
DSBplot-process -o $process_dir/Sense_R1 --stages 4_info

# Example of concatenating three experiments.
# Note, this does not make biological sense, but is used here for demonstration purposes.
# The three experiments have the same windowed reference sequence and number of repeats,
# which allows them to be concatenated.
DSBplot-concat -i $process_dir/Sense_R1 $process_dir/BranchD_R1 $process_dir/pCMVD_R1 -o $process_dir/concat_R1 --names 1 2 3 4 --reads 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000
DSBplot-concat -i $process_dir/Sense_R2 $process_dir/BranchD_R2 $process_dir/pCMVD_R2 -o $process_dir/concat_R2 --names 1 2 3 4 --reads 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000

### Plot graphs with individual layouts ###
foreach ($ext in $graph_exts) {
  foreach ($layout in $graph_layouts) {
    foreach ($x in $constructs_plot) {
      if ($x.Length -eq 1) {
        $con = $x[0]
        DSBplot-graph -i $process_dir/${con} -o $plots_dir/graph/$ext/$layout/${con}.${ext} --debug debug/$layout --layout $layout --ul_yax_x 0 --ul_xax_del_y 0 --ul_xax_ins_y 0 --width 2400 --height 1800 --font_scale 3 --legend_x 100 --legend_y 0 --mar_r 900 --legends size var_type --sep $sep[$layout]
      } elseif ($x.Length -eq 3) {
        $con = $x[0]
        $con1 = $x[1]
        $con2 = $x[2]
        DSBplot-graph -i $process_dir/${con1}::$process_dir/${con2} -o $plots_dir/graph/$ext/$layout/${con}.${ext} --debug debug/$layout --layout $layout --ul_yax_x 0 --ul_xax_del_y 0 --ul_xax_ins_y 0 --width 2400 --height 1800 --ratio_colors "#cf191b" "#33a02c" --colorbar_scale 3 --font_scale 3 --legend_x 100 --legend_y -100 --mar_r 900 --legends ratio_cont --sep $sep[$layout]
      }
    }
  }
}

### Plot graphs with combined layouts ###
foreach ($ext in $graph_exts) {
  foreach ($layout in $graph_layouts) {
    $output_dir = "$plots_dir/graph/$ext/${layout}_combined"
    if ($layout -eq "universal") {
      # Additional formatting arguments for universal layout
      DSBplot-graph `
        -i $process_dir/BranchD_R1 $process_dir/BranchD_R2 $process_dir/Sense_R1 $process_dir/Sense_R2 $process_dir/pCMVD_R1 $process_dir/pCMVD_R2 $process_dir/Sense_R1::$process_dir/BranchD_R1 $process_dir/Sense_R2::$process_dir/BranchD_R2 $process_dir/Sense_R1::$process_dir/pCMVD_R1 $process_dir/Sense_R2::$process_dir/pCMVD_R2 `
        --title "BranchΔ (R1)" "BranchΔ (R2)" "Sense (R1)" "Sense (R2)" "pCMVΔ (R1)" "pCMVΔ (R2)" "Sense / BranchΔ (R1)" "Sense / BranchΔ (R2)" "Sense / pCMVΔ (R1)" "Sense / pCMVΔ (R2)" `
        -o $output_dir/BranchD_R1.${ext} $output_dir/BranchD_R2.${ext} $output_dir/Sense_R1.${ext} $output_dir/Sense_R2.${ext} $output_dir/pCMVD_R1.${ext} $output_dir/pCMVD_R2.${ext} $output_dir/Sense_BranchD_R1.${ext} $output_dir/Sense_BranchD_R2.${ext} $output_dir/Sense_pCMVD_R1.${ext} $output_dir/Sense_pCMVD_R2.${ext} `
        --debug debug/${layout}_combined `
        --legends ratio_cont var_type size `
        --ul_yax_x 12 --ul_xax_del_y 0 --ul_xax_ins_y 0 `
        --ul_y_scale_ins 10 --ul_y_scale_del 8 `
        --range_x -12 13 --range_y -175 225 `
        --legend_x 100 --legend_y -100 `
        --mar_t 300 --mar_r 1200 --mar_l 0 --mar_b 0 `
        --font_scale 4 --colorbar_scale 3 `
        --layout $layout --width 2400 --height 3000 `
        --sep $sep[$layout]
    } else {
      DSBplot-graph `
        -i $process_dir/BranchD_R1 $process_dir/BranchD_R2 $process_dir/Sense_R1 $process_dir/Sense_R2 $process_dir/pCMVD_R1 $process_dir/pCMVD_R2 $process_dir/Sense_R1::$process_dir/BranchD_R1 $process_dir/Sense_R2::$process_dir/BranchD_R2 $process_dir/Sense_R1::$process_dir/pCMVD_R1 $process_dir/Sense_R2::$process_dir/pCMVD_R2 `
        --title "BranchΔ (R1)" "BranchΔ (R2)" "Sense (R1)" "Sense (R2)" "pCMVΔ (R1)" "pCMVΔ (R2)" "Sense / BranchΔ (R1)" "Sense / BranchΔ (R2)" "Sense / pCMVΔ (R1)" "Sense / pCMVΔ (R2)" `
        -o $output_dir/BranchD_R1.${ext} $output_dir/BranchD_R2.${ext} $output_dir/Sense_R1.${ext} $output_dir/Sense_R2.${ext} $output_dir/pCMVD_R1.${ext} $output_dir/pCMVD_R2.${ext} $output_dir/Sense_BranchD_R1.${ext} $output_dir/Sense_BranchD_R2.${ext} $output_dir/Sense_pCMVD_R1.${ext} $output_dir/Sense_pCMVD_R2.${ext} `
        --debug debug/${layout}_combined `
        --legends ratio_cont var_type `
        --ul_yax_x 0 --ul_xax_del_y 0 --ul_xax_ins_y 0 `
        --legend_x 100 --legend_y -100 `
        --mar_t 300 --mar_r 900 --mar_l 0 --mar_b 0 `
        --font_scale 3 --colorbar_scale 3 `
        --layout $layout --width 2400 --height 1800 `
        --sep $sep[$layout]
    }
  }
}

### Plot histograms ###
foreach ($ext in $histogram_exts) {
  foreach ($var in $histograms_var_types) {
    foreach ($con in $constructs) {
      DSBplot-histogram -i $process_dir/${con} -o $plots_dir/histogram/$ext/${con}_${var}.${ext} --color $variation_color[$var] --var $var --xax rel
    }
  }
}
