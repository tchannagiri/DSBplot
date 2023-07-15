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
  "substitution" = "#bfbfbf";
  "insertion" = "#ffa500";
  "deletion" = "#8080ff";
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
  DSBplot-preprocess --input input/fastq/${c}_1.fq input/fastq/${c}_2.fq input/fastq/${c}_3.fq input/fastq/${c}_4.fq --output output/$con --ref_seq_file input/ref_seq/2DSB_${c}.fa --label $labels[$con] --total_reads 3000 3000 3000 3000 --dsb_pos $dsb_pos[$con]
}
  
# Example of running preprocessing stages separately for one experiment
DSBplot-preprocess --input input/fastq/Sense_R1_1.fq input/fastq/Sense_R1_2.fq input/fastq/Sense_R1_3.fq input/fastq/Sense_R1_4.fq --ref_seq_file input/ref_seq/2DSB_Sense_R1.fa --output output/Sense_R1 --stages 0_align
DSBplot-preprocess --ref_seq_file input/ref_seq/2DSB_Sense_R1.fa --dsb_pos 67 --output output/Sense_R1 --stages 1_filter
DSBplot-preprocess --ref_seq_file input/ref_seq/2DSB_Sense_R1.fa --dsb_pos 67 --output output/Sense_R1 --label Sense_R1 --total_reads 3000 3000 3000 3000 --stages 2_window
DSBplot-preprocess --output output/Sense_R1 --stages 3_variation

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
        DSBplot-graph --input output/${con} --output plots/graph/$ext/$layout/${con}.${ext} --debug debug/$layout --layout universal --universal_y_axis_x_pos 0 --universal_x_axis_deletion_y_pos 0 --universal_x_axis_insertion_y_pos 0 --width 2400 --height 1800 --font_size_scale 3 --legend_x_shift_px 100 --legend_y_shift_px 0 --margin_right_px 900 --legends size variation_type 
      } elseif ($x.Length -eq 3) {
        $con = $x[0]
        $con1 = $x[1]
        $con2 = $x[2]
        DSBplot-graph --input output/${con1}::output/${con2} --output plots/graph/$ext/$layout/${con}.${ext} --debug debug/$layout --layout universal --universal_y_axis_x_pos 0 --universal_x_axis_deletion_y_pos 0 --universal_x_axis_insertion_y_pos 0 --width 2400 --height 1800 --node_comparison_colors "#cf191b" "#33a02c" --legend_colorbar_scale 3 --font_size_scale 3 --legend_x_shift_px 100 --legend_y_shift_px -100 --margin_right_px 900 --legends freq_ratio_continuous 
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
    --legends freq_ratio_continuous variation_type `
    --universal_y_axis_x_pos 0 --universal_x_axis_deletion_y_pos 0 --universal_x_axis_insertion_y_pos 0 `
    --legend_x_shift_px 100 --legend_y_shift_px -100 `
    --margin_top_px 300 --margin_right_px 900 --margin_left_px 0 --margin_bottom_px 0 `
    --font_size_scale 3 --legend_colorbar_scale 3 `
    --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout $layout --width 2400 --height 1800
  }
}

### Plot histograms ###
foreach ($ext in @("png", "pdf", "html")) {
  foreach ($var in @("substitution", "insertion", "deletion")) {
    foreach ($con in $constructs) {
      DSBplot-histogram --input output/${con} --output plots/histogram/png/${con}_${var}.png --color $variation_color[$var] --variation_type $var --label_type relative
    }
  }
}
