### Define constants ###
$fastq_dir = "input/fastq"
$preprocess_dir = "output"
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

$graph_exts = @("png", "pdf", "html")
$graph_layouts = @("universal", "kamada", "radial")

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

### Preprocessing ###
foreach ($con in $constructs) {
  DSBplot-preprocess -i $fastq_dir/${con}_1.fq $fastq_dir/${con}_2.fq $fastq_dir/${con}_3.fq $fastq_dir/${con}_4.fq -o $preprocess_dir/$con --ref $ref_seq_file[$con] --label $labels[$con] --reads 3000 3000 3000 3000 --dsb $dsb_pos --rc $reverse_complement[$con]
}

# Example of running preprocessing stages separately for one experiment
DSBplot-preprocess -i input/fastq/Sense_R1_1.fq input/fastq/Sense_R1_2.fq input/fastq/Sense_R1_3.fq input/fastq/Sense_R1_4.fq --ref input/ref_seq/2DSB_Sense_R1.fa -o $preprocess_dir/Sense_R1 --stages 0_align
DSBplot-preprocess --ref input/ref_seq/2DSB_Sense_R1.fa --dsb 67 -o $preprocess_dir/Sense_R1 --stages 1_filter
DSBplot-preprocess -o $preprocess_dir/Sense_R1 --stages 2_window
DSBplot-preprocess -o $preprocess_dir/Sense_R1 --stages 3_variation
DSBplot-preprocess -o $preprocess_dir/Sense_R1 --stages 4_info

# Example of concatenating three experiments.
# Note, this does not make biological sense, but is used here for demonstration purposes.
# The three experiments have the same windowed reference sequence and number of repeats,
# which allows them to be concatenated.
DSBplot-concat -i $preprocess_dir/Sense_R1 $preprocess_dir/BranchD_R1 $preprocess_dir/pCMVD_R1 -o $preprocess_dir/concat_R1 --names 1 2 3 4 --reads 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000
DSBplot-concat -i $preprocess_dir/Sense_R2 $preprocess_dir/BranchD_R2 $preprocess_dir/pCMVD_R2 -o $preprocess_dir/concat_R2 --names 1 2 3 4 --reads 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000

### Plot graphs with individual layouts ###
foreach ($ext in $graph_exts) {
  foreach ($layout in $graph_layouts) {
    foreach ($x in $constructs_plot) {
      if ($x.Length -eq 1) {
        $con = $x[0]
        DSBplot-graph -i $preprocess_dir/${con} -o $plots_dir/graph/$ext/$layout/${con}.${ext} --debug debug/$layout --layout universal --ul_yax_x 0 --ul_xax_del_y 0 --ul_xax_ins_y 0 --width 2400 --height 1800 --font_scale 3 --legend_x 100 --legend_y 0 --mar_r 900 --legends size var_type 
      } elseif ($x.Length -eq 3) {
        $con = $x[0]
        $con1 = $x[1]
        $con2 = $x[2]
        DSBplot-graph -i $preprocess_dir/${con1}::$preprocess_dir/${con2} -o $plots_dir/graph/$ext/$layout/${con}.${ext} --debug debug/$layout --layout universal --ul_yax_x 0 --ul_xax_del_y 0 --ul_xax_ins_y 0 --width 2400 --height 1800 --ratio_colors "#cf191b" "#33a02c" --colorbar_scale 3 --font_scale 3 --legend_x 100 --legend_y -100 --mar_r 900 --legends ratio_cont
      }
    }
  }
}

### Plot graphs with combined layouts ###
foreach ($ext in $graph_exts) {
  foreach ($layout in $graph_layouts) {
    $output_dir = "$plots_dir/graph/$ext/${layout}_combined"
    DSBplot-graph `
    -i $preprocess_dir/BranchD_R1 $preprocess_dir/BranchD_R2 $preprocess_dir/Sense_R1 $preprocess_dir/Sense_R2 $preprocess_dir/pCMVD_R1 $preprocess_dir/pCMVD_R2 $preprocess_dir/Sense_R1::$preprocess_dir/BranchD_R1 $preprocess_dir/Sense_R2::$preprocess_dir/BranchD_R2 $preprocess_dir/Sense_R1::$preprocess_dir/pCMVD_R1 $preprocess_dir/Sense_R2::$preprocess_dir/pCMVD_R2 `
    --title "BranchΔ (R1)" "BranchΔ (R2)" "Sense (R1)" "Sense (R2)" "pCMVΔ (R1)" "pCMVΔ (R2)" "Sense / BranchΔ (R1)" "Sense / BranchΔ (R2)" "Sense / pCMVΔ (R1)" "Sense / pCMVΔ (R2)" `
    -o $output_dir/BranchD_R1.${ext} $output_dir/BranchD_R2.${ext} $output_dir/Sense_R1.${ext} $output_dir/Sense_R2.${ext} $output_dir/pCMVD_R1.${ext} $output_dir/pCMVD_R2.${ext} $output_dir/Sense_BranchD_R1.${ext} $output_dir/Sense_BranchD_R2.${ext} $output_dir/Sense_pCMVD_R1.${ext} $output_dir/Sense_pCMVD_R2.${ext} `
    --debug debug/${layout}_combined `
    --legends ratio_cont var_type `
    --ul_yax_x 0 --ul_xax_del_y 0 --ul_xax_ins_y 0 `
    --legend_x 100 --legend_y -100 `
    --mar_t 300 --mar_r 900 --mar_l 0 --mar_b 0 `
    --font_scale 3 --colorbar_scale 3 `
    --layout $layout --width 2400 --height 1800
  }
}

### Plot histograms ###
foreach ($ext in $histogram_exts) {
  foreach ($var in $histograms_var_types) {
    foreach ($con in $constructs) {
      DSBplot-histogram -i $preprocess_dir/${con} -o $plots_dir/histogram/$ext/${con}_${var}.${ext} --color $variation_color[$var] --var $var --xax rel
    }
  }
}
