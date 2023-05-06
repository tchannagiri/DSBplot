$lib_array = (
  ("OX", "BranchD", "R1"),
  ("OX", "BranchD", "R2"),
  ("OX", "Sense", "R1"),
  ("OX", "Sense", "R2"),
  ("OX", "pCMVD", "R1"),
  ("OX", "pCMVD", "R2"),
  ("WT", "BranchD", "R1"),
  ("WT", "BranchD", "R2"),
  ("WT", "Sense", "R1"),
  ("WT", "Sense", "R2"),
  ("WT", "pCMVD", "R1"),
  ("WT", "pCMVD", "R2")
)
foreach ($lib in $lib_array) {
  $cell = $lib[0]
  $construct = $lib[1]
  $strand = $lib[2]
  if ($strand -eq "R1") {
    $dsb_pos = 67
  } else {
    $dsb_pos = 46
  }
  $output = "data_OX/output/${cell}_${construct}_${strand}"
  $ref_seq_file = "data_OX/input/ref_seq/${construct}_${strand}.fa"
  python preprocess.py --output $output --dsb_pos ${dsb_pos} --ref_seq_file ${ref_seq_file} --label ${construct} --stages 2_combine 3_window 4_graph 5_histogram
}

$lib_array_comp = (
  ("OX", "Sense", "BranchD", "R1"),
  ("OX", "Sense", "BranchD", "R2"),
  ("OX", "Sense", "pCMVD", "R1"),
  ("OX", "Sense", "pCMVD", "R2"),
  ("WT", "Sense", "BranchD", "R1"),
  ("WT", "Sense", "BranchD", "R2"),
  ("WT", "Sense", "pCMVD", "R1"),
  ("WT", "Sense", "pCMVD", "R2")
)
foreach ($lib in $lib_array_comp) {
  $cell = $lib[0]
  $construct_1 = $lib[1]
  $construct_2 = $lib[2]
  $strand = $lib[3]
  if ($strand -eq "R1") {
    $dsb_pos = 67
  } else {
    $dsb_pos = 46
  }
  $output_1 = "data_OX/output/${cell}_${construct_1}_${strand}"
  $output_2 = "data_OX/output/${cell}_${construct_2}_${strand}"
  $output = "data_OX/output/${cell}_${construct_1}_${construct_2}_${strand}"
  python comparison.py --output $output --input $output_1 $output_2
}

# Plot graphs for individual data
foreach ($cell in ("OX", "WT")) {
  foreach ($construct in ("Sense", "BranchD", "pCMVD")) {
    foreach ($strand in ("R1", "R2")) {
      python graph.py `
      --input ./data_OX/output/${cell}_${construct}_${strand}/ `
      --output ./plot_OX/graph/${cell}_${construct}_${strand}.png `
      --layout universal_layout --width 2400 --height 1800 `
      --range_x -12 12 --range_y -19 24 `
      --universal_layout_x_axis_deletion_y_pos -18.5 `
      --universal_layout_x_axis_insertion_y_pos 23 `
      --universal_layout_y_axis_x_pos 11 `
      --universal_layout_y_axis_y_range -18 23 `
      --universal_layout_y_axis_deletion_max_tick 18 `
      --universal_layout_y_axis_insertion_max_tick 7 `
      --subst_type withoutSubst
    }
  }
}

# Plot graphs for comparison data
foreach ($cell in ("OX", "WT")) {
  foreach ($construct_pair in (("Sense", "BranchD"), ("Sense", "pCMVD"))) {
    $construct1 = $construct_pair[0]
    $construct2 = $construct_pair[1]
    foreach ($strand in ("R1", "R2")) {
      python graph.py `
      --input ./data_OX/output/${cell}_${construct1}_${construct2}_${strand}/ `
      --output ./plot_OX/graph/${cell}_${construct1}_${construct2}_${strand}.png `
      --node_comparison_colors "#cf191b" "#33a02c" `
      --layout universal_layout --width 2400 --height 1800 `
      --range_x -12 12 --range_y -19 24 `
      --universal_layout_x_axis_deletion_y_pos -18.5 `
      --universal_layout_x_axis_insertion_y_pos 23 `
      --universal_layout_y_axis_x_pos 11 `
      --universal_layout_y_axis_y_range -18 23 `
      --universal_layout_y_axis_deletion_max_tick 18 `
      --universal_layout_y_axis_insertion_max_tick 7 `
      --subst_type withoutSubst
    }
  }
}

# Plot histograms
foreach ($cell in ("OX", "WT")) {
  foreach ($construct in ("Sense", "BranchD", "pCMVD")) {
    foreach ($strand in ("R1", "R2")) {
      foreach ($var in ("deletion", "insertion", "substitution")) {
        $color = @{
          "deletion" = "#8080ff"
          "insertion" = "#ffa500"
          "substitution" = "#bfbfbf"
        }[$var]
        python histogram.py `
        --input ./data_OX/output/${cell}_${construct}_${strand}/ `
        --output ./plot_OX/histogram/${cell}_${construct}_${strand}_${var}.png `
        --variation_type ${var} --color ${color} --label_type relative
      }
    }
  }
}
