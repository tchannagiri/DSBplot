### Define constants ###
$labels = @{
  "Sense_R1" = "Sense";
  "Sense_R2" = "Sense";
  "BranchD_R1" = "BranchΔ";
  "BranchD_R2" = "BranchΔ";
  "pCMVD_R1" = "pCMVΔ";
  "pCMVD_R2" = "pCMVΔ";
}

$dsb_pos = @{
  "Sense_R1" = 67;
  "Sense_R2" = 46;
  "BranchD_R1" = 67;
  "BranchD_R2" = 46;
  "pCMVD_R1" = 67;
  "pCMVD_R2" = 46;
}

$variation_color = @{
  "sub" = "#bfbfbf";
  "ins" = "#ffa500";
  "del" = "#8080ff";
}

$constructs = @(
  "Sense_R1",
  "Sense_R2",
  "BranchD_R1",
  "BranchD_R2",
  "pCMVD_R1",
  "pCMVD_R2"
)

### Preprocessing ###
foreach ($con in $constructs) {
  DSBplot-preprocess --input input/fastq/${con}_1.fq input/fastq/${con}_2.fq input/fastq/${con}_3.fq input/fastq/${con}_4.fq --output output/$con --ref input/ref_seq/2DSB_${con}.fa --label $labels[$con] --reads 3000 3000 3000 3000 --dsb $dsb_pos[$con]
}

# Example of running preprocessing stages separately for one experiment
DSBplot-preprocess --input input/fastq/Sense_R1_1.fq input/fastq/Sense_R1_2.fq input/fastq/Sense_R1_3.fq input/fastq/Sense_R1_4.fq --ref_seq_file input/ref_seq/2DSB_Sense_R1.fa --output output/Sense_R1 --stages 0_align
DSBplot-preprocess --ref_seq_file input/ref_seq/2DSB_Sense_R1.fa --dsb_pos 67 --output output/Sense_R1 --stages 1_filter
DSBplot-preprocess --ref_seq_file input/ref_seq/2DSB_Sense_R1.fa --dsb_pos 67 --output output/Sense_R1 --label Sense_R1 --total_reads 3000 3000 3000 3000 --stages 2_window
DSBplot-preprocess --output output/Sense_R1 --stages 3_variation

# Example of concatenating three experiments.
# Note, this does not make biological sense, but is used here for demonstration purposes.
# The three experiments have the same windowed reference sequence and number of repeats,
# which allows them to be concatenated.
DSBplot-concat --input output/Sense_R1 output/BranchD_R1 output/pCMVD_R1 --output output/concat_R1 --names 1 2 3 4 --total_reads 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000
DSBplot-concat --input output/Sense_R2 output/BranchD_R2 output/pCMVD_R2 --output output/concat_R2 --names 1 2 3 4 --total_reads 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000 3000

### Plot graphs with individual layouts ###
foreach ($ext in @("png", "pdf", "html")) {
  foreach ($layout in @("universal", "kamada", "radial")) {
    foreach (
      $x in @(
        @("BranchD_R1"),
        @("BranchD_R2"),
        @("Sense_R1"),
        @("Sense_R2"),
        @("pCMVD_R1"),
        @("pCMVD_R2"),
        @("Sense_BranchD_R1", "Sense_R1", "BranchD_R1"),
        @("Sense_BranchD_R2", "Sense_R2", "BranchD_R2"),
        @("Sense_pCMVD_R1", "Sense_R1", "pCMVD_R1"),
        @("Sense_pCMVD_R2", "Sense_R2", "pCMVD_R2")
      )
    ) {
      if ($x.Length -eq 1) {
        $con = $x[0]
        DSBplot-graph --input output/${con} --output plots/graph/$ext/$layout/${con}.${ext} --debug debug/$layout --layout universal --ul_yax_x 0 --ul_xax_del_y 0 --ul_xax_ins_y 0 --width 2400 --height 1800 --font_scale 3 --legend_x 100 --legend_y 0 --mar_r 900 --legends size var_type 
      } elseif ($x.Length -eq 3) {
        $con = $x[0]
        $con1 = $x[1]
        $con2 = $x[2]
        DSBplot-graph --input output/${con1}::output/${con2} --output plots/graph/$ext/$layout/${con}.${ext} --debug debug/$layout --layout universal --ul_yax_x 0 --ul_xax_del_y 0 --ul_xax_ins_y 0 --width 2400 --height 1800 --ratio_colors "#cf191b" "#33a02c" --colorbar_scale 3 --font_scale 3 --legend_x 100 --legend_y -100 --mar_r 900 --legends ratio_cont
      }
    }
  }
}

### Plot graphs with combined layouts ###
foreach ($layout in @("universal", "kamada", "radial")) {
  foreach ($ext in @("png", "pdf", "html")) {
    $output_dir = "plots/graph/$ext/${layout}_combined"
    DSBplot-graph `
    --input output/BranchD_R1 output/BranchD_R2 output/Sense_R1 output/Sense_R2 output/pCMVD_R1 output/pCMVD_R2 output/Sense_R1::output/BranchD_R1 output/Sense_R2::output/BranchD_R2 output/Sense_R1::output/pCMVD_R1 output/Sense_R2::output/pCMVD_R2 `
    --title "BranchΔ (R1)" "BranchΔ (R2)" "Sense (R1)" "Sense (R2)" "pCMVΔ (R1)" "pCMVΔ (R2)" "Sense / BranchΔ (R1)" "Sense / BranchΔ (R2)" "Sense / pCMVΔ (R1)" "Sense / pCMVΔ (R2)" `
    --output $output_dir/BranchD_R1.${ext} $output_dir/BranchD_R2.${ext} $output_dir/Sense_R1.${ext} $output_dir/Sense_R2.${ext} $output_dir/pCMVD_R1.${ext} $output_dir/pCMVD_R2.${ext} $output_dir/Sense_BranchD_R1.${ext} $output_dir/Sense_BranchD_R2.${ext} $output_dir/Sense_pCMVD_R1.${ext} $output_dir/Sense_pCMVD_R2.${ext} `
    --debug debug/${layout}_combined `
    --legends ratio_cont var_type `
    --ul_yax_x 0 --ul_xax_del_y 0 --ul_xax_ins_y 0 `
    --legend_x 100 --legend_y -100 `
    --mar_t 300 --mar_r 900 --mar_l 0 --mar_b 0 `
    --font_scale 3 --colorbar_scale 3 `
    --rc 0 1 0 1 0 1 0 1 0 1 --layout $layout --width 2400 --height 1800
  }
}

### Plot histograms ###
foreach ($ext in @("png", "pdf")) {
  foreach ($var in @("sub", "ins", "del")) {
    foreach ($con in $constructs) {
      DSBplot-histogram --input output/${con} --output plots/histogram/$ext/${con}_${var}.${ext} --color $variation_color[$var] --var $var --xax rel
    }
  }
}
