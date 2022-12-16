python preprocess.py --input data/0_fastq/db1_R1.fq data/0_fastq/db2_R1.fq data/0_fastq/db3_R1.fq data/0_fastq/db4_R1.fq --ref_seq_file data/0_ref_seq/1DSB_R1_branch.fa --dsb_pos 50 --output data/output --label db_R1 --total_reads 3000 3000 3000 3000

python graph.py --input data/output --output plot/graph/output.png --layout universal --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

python histogram.py --input data/output --output plot/histogram/output_substitution.png --variation_type substitution --label_type relative
python histogram.py --input data/output --output plot/histogram/output_insertion.png --variation_type insertion --label_type relative
python histogram.py --input data/output --output plot/histogram/output_deletion.png --variation_type deletion --label_type relative
