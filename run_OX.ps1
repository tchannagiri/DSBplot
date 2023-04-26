$lib_array = (
  ("OX", "BranchD", "R1", ("yjl446", "yjl447", "yjl448", "yjl449")),
  ("OX", "BranchD", "R2", ("yjl446", "yjl447", "yjl448", "yjl449")),
  ("OX", "Sense", "R1", ("yjl442", "yjl443", "yjl444", "yjl445")),
  ("OX", "Sense", "R2", ("yjl442", "yjl443", "yjl444", "yjl445")),
  ("OX", "pCMVD", "R1", ("yjl450", "yjl451", "yjl452", "yjl453")),
  ("OX", "pCMVD", "R2", ("yjl450", "yjl451", "yjl452", "yjl453")),
  ("WT", "BranchD", "R1", ("yjl458", "yjl459", "yjl460", "yjl461")),
  ("WT", "BranchD", "R2", ("yjl458", "yjl459", "yjl460", "yjl461")),
  ("WT", "Sense", "R1", ("yjl454", "yjl455", "yjl456", "yjl457")),
  ("WT", "Sense", "R2", ("yjl454", "yjl455", "yjl456", "yjl457")),
  ("WT", "pCMVD", "R1", ("yjl462", "yjl463", "yjl464", "yjl465")),
  ("WT", "pCMVD", "R2", ("yjl462", "yjl463", "yjl464", "yjl465"))
)
foreach ($lib in $lib_array) {
  $cell = $lib[0]
  $construct = $lib[1]
  $strand = $lib[2]
  $names = $lib[3]
  foreach ($subst in ("withSubst", "withoutSubst")) {
  if (
    (($cell -eq "OX") -and ($construct -eq "Sense") -and ($strand -eq "R1")) -or
    (($cell -eq "WT") -and ($construct -eq "pCMVD") -and ($strand -eq "R1")) -or 
    (($cell -eq "WT") -and ($construct -eq "BranchD") -and ($strand -eq "R1"))
  ) {
    continue
  }
  python ./1_process_nhej/combine_repeat.py --input `
    ./data_input_OX/2DSB_${cell}_${construct}_${strand}/1_filter_nhej/$($names[0])_$strand.tsv `
    ./data_input_OX/2DSB_${cell}_${construct}_${strand}/1_filter_nhej/$($names[1])_$strand.tsv `
    ./data_input_OX/2DSB_${cell}_${construct}_${strand}/1_filter_nhej/$($names[2])_$strand.tsv `
    ./data_input_OX/2DSB_${cell}_${construct}_${strand}/1_filter_nhej/$($names[3])_$strand.tsv `
    --output ./data_input_OX/2DSB_${cell}_${construct}_${strand}/2_combine_repeat/out.tsv `
    --column_names $($names[0]) $($names[1]) $($names[2]) $($names[3])
  }
}

# Extract window data and get graph/histogram data
foreach ($cell in ("OX", "WT")) {
  foreach ($construct in ("Sense", "BranchD", "pCMVD")) {
    foreach ($strand in ("R1", "R2")) {
      foreach ($subst in ("withSubst", "withoutSubst")) {
        $construct2 = if ($construct -eq "BranchD") {
          "branch"
        } else {
          if ($construct -eq "Sense") {
            "sense"
          } else {
            if ($construct -eq "pCMVD") {
              "cmv"
            }
          }
        }
        if (
          (($cell -eq "OX") -and ($construct -eq "Sense") -and ($strand -eq "R1")) -or
          (($cell -eq "WT") -and ($construct -eq "pCMVD") -and ($strand -eq "R1")) -or 
          (($cell -eq "WT") -and ($construct -eq "BranchD") -and ($strand -eq "R1"))
        ) {
          continue
        }
        if ($strand -eq "R1") {
          $dsb_pos = 67
        } else {
          $dsb_pos = 46
        }
        # python ./2_get_window_data/get_window.py `
        #   --input ./data_input_OX/2DSB_${cell}_${construct}_${strand}/2_combine_repeat/out.tsv `
        #   --output ./data_input_OX/2DSB_${cell}_${construct}_${strand}/3_window/ `
        #   --ref_seq_file ./data_input/ref_seq/2DSB_${strand}_${construct2}.fa `
        #   --dsb_pos ${dsb_pos} --subst_type ${subst} --label ${strand}
        # python ./2_get_window_data/get_freq.py `
        #   --input ./data_input_OX/2DSB_${cell}_${construct}_${strand}/3_window/ `
        #   --output ./data_input_OX/2DSB_${cell}_${construct}_${strand}/3_window/ `
        #   --subst_type ${subst}
        python ./3_get_graph_data/get_graph_data.py `
          --input ./data_input_OX/2DSB_${cell}_${construct}_${strand}/3_window/ `
          --output ./data_input_OX/2DSB_${cell}_${construct}_${strand}/4_graph/ `
          --subst_type ${subst}
        python ./4_get_histogram_data/get_histogram_data.py `
          --input ./data_input_OX/2DSB_${cell}_${construct}_${strand}/4_graph/ `
          --output ./data_input_OX/2DSB_${cell}_${construct}_${strand}/5_histogram/ `
          --subst_type ${subst}
      }
    }
  }
}

# Plot graphs
foreach ($cell in ("OX", "WT")) {
  foreach ($construct in ("Sense", "BranchD", "pCMVD")) {
    foreach ($strand in ("R1", "R2")) {
      foreach ($subst in ("withSubst", "withoutSubst")) {
        if (
          (($cell -eq "OX") -and ($construct -eq "Sense") -and ($strand -eq "R1")) -or
          (($cell -eq "WT") -and ($construct -eq "pCMVD") -and ($strand -eq "R1")) -or 
          (($cell -eq "WT") -and ($construct -eq "BranchD") -and ($strand -eq "R1"))
        ) {
          continue
        }
        python graph.py `
        --input ./data_input_OX/2DSB_${cell}_${construct}_${strand}/ `
        --output ./plot/OX/graph/2DSB_${cell}_${construct}_${strand}.png `
        --layout universal_layout --width 2400 --height 1800 `
        --range_x -13 14 --range_y -19 6 `
        --universal_layout_x_axis_deletion_y_pos -18.5 `
        --universal_layout_x_axis_insertion_y_pos 5 `
        --universal_layout_y_axis_x_pos 13 `
        --universal_layout_y_axis_y_range -17 3 `
        --universal_layout_y_axis_deletion_max_tick 17 `
        --universal_layout_y_axis_insertion_max_tick 1 `
        --subst_type withoutSubst --quiet
      }
    }
  }
}

# Plot histograms
foreach ($cell in ("OX", "WT")) {
  foreach ($construct in ("Sense", "BranchD", "pCMVD")) {
    foreach ($strand in ("R1", "R2")) {
      foreach ($var in ("deletion", "insertion", "substitution")) {
        if (
          (($cell -eq "OX") -and ($construct -eq "Sense") -and ($strand -eq "R1")) -or
          (($cell -eq "WT") -and ($construct -eq "pCMVD") -and ($strand -eq "R1")) -or 
          (($cell -eq "WT") -and ($construct -eq "BranchD") -and ($strand -eq "R1"))
        ) {
          continue
        }
        $color = @{
          "deletion" = "#8080ff"
          "insertion" = "#ffa500"
          "substitution" = "#bfbfbf"
        }[$var]
        python histogram.py `
        --input ./data_input_OX/2DSB_${cell}_${construct}_${strand}/ `
        --output ./plot/OX/histogram/2DSB_${cell}_${construct}_${strand}_${var}.png `
        --variation_type ${var} --color ${color} --label_type relative
      }
    }
  }
}
