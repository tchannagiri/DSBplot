python preprocess.py --input data/0_fastq/db1_R1.fq data/0_fastq/db2_R1.fq data/0_fastq/db3_R1.fq data/0_fastq/db4_R1.fq --ref_seq_file data/0_ref_seq/2DSB_R1_branch.fa --dsb_pos 50 --output data/db_R1 --label db_R1 --total_reads 3000 3000 3000 3000
python preprocess.py --input data/0_fastq/db1_R2.fq data/0_fastq/db2_R2.fq data/0_fastq/db3_R2.fq data/0_fastq/db4_R2.fq --ref_seq_file data/0_ref_seq/2DSB_R2_branch.fa --dsb_pos 50 --output data/db_R2 --label db_R2 --total_reads 3000 3000 3000 3000
python preprocess.py --input data/0_fastq/sense1_R1.fq data/0_fastq/sense2_R1.fq data/0_fastq/sense3_R1.fq data/0_fastq/sense4_R1.fq --ref_seq_file data/0_ref_seq/2DSB_R1_sense.fa --dsb_pos 50 --output data/sense_R1 --label sense_R1 --total_reads 3000 3000 3000 3000
python preprocess.py --input data/0_fastq/sense1_R2.fq data/0_fastq/sense2_R2.fq data/0_fastq/sense3_R2.fq data/0_fastq/sense4_R2.fq --ref_seq_file data/0_ref_seq/2DSB_R2_sense.fa --dsb_pos 50 --output data/sense_R2 --label sense_R2 --total_reads 3000 3000 3000 3000
python preprocess.py --input data/0_fastq/dcmv1_R1.fq data/0_fastq/dcmv2_R1.fq data/0_fastq/dcmv3_R1.fq data/0_fastq/dcmv4_R1.fq --ref_seq_file data/0_ref_seq/2DSB_R1_cmv.fa --dsb_pos 50 --output data/dcmv_R1 --label dcmv_R1 --total_reads 3000 3000 3000 3000
python preprocess.py --input data/0_fastq/dcmv1_R2.fq data/0_fastq/dcmv2_R2.fq data/0_fastq/dcmv3_R2.fq data/0_fastq/dcmv4_R2.fq --ref_seq_file data/0_ref_seq/2DSB_R2_cmv.fa --dsb_pos 50 --output data/dcmv_R2 --label dcmv_R2 --total_reads 3000 3000 3000 3000

python comparison.py --input data/sense_R1 data/db_R1 --output data/sense_db_R1
python comparison.py --input data/sense_R2 data/db_R2 --output data/sense_db_R2
python comparison.py --input data/sense_R1 data/dcmv_R1 --output data/sense_dcmv_R1
python comparison.py --input data/sense_R2 data/dcmv_R2 --output data/sense_dcmv_R2

python graph.py --input data/db_R1 --output plot/graph/db_R1.png --layout universal --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data/db_R2 --output plot/graph/db_R2.png --layout universal --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data/sense_R1 --output plot/graph/sense_R1.png --layout universal --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data/sense_R2 --output plot/graph/sense_R2.png --layout universal --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data/dcmv_R1 --output plot/graph/dcmv_R1.png --layout universal --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data/dcmv_R2 --output plot/graph/dcmv_R2.png --layout universal --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

python graph.py --input data/sense_db_R1 --output plot/graph/sense_db_R1.png --layout universal --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data/sense_db_R2 --output plot/graph/sense_db_R2.png --layout universal --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data/sense_dcmv_R1 --output plot/graph/sense_dcmv_R1.png --layout universal --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data/sense_dcmv_R2 --output plot/graph/sense_dcmv_R2.png --layout universal --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

python histogram.py --input data/db_R1 --output plot/histogram/db_R1_substitution.png --variation_type substitution --label_type relative
python histogram.py --input data/db_R1 --output plot/histogram/db_R1_insertion.png --variation_type insertion --label_type relative
python histogram.py --input data/db_R1 --output plot/histogram/db_R1_deletion.png --variation_type deletion --label_type relative
python histogram.py --input data/db_R2 --output plot/histogram/db_R2_substitution.png --variation_type substitution --label_type relative
python histogram.py --input data/db_R2 --output plot/histogram/db_R2_insertion.png --variation_type insertion --label_type relative
python histogram.py --input data/db_R2 --output plot/histogram/db_R2_deletion.png --variation_type deletion --label_type relative
python histogram.py --input data/sense_R1 --output plot/histogram/sense_R1_substitution.png --variation_type substitution --label_type relative
python histogram.py --input data/sense_R1 --output plot/histogram/sense_R1_insertion.png --variation_type insertion --label_type relative
python histogram.py --input data/sense_R1 --output plot/histogram/sense_R1_deletion.png --variation_type deletion --label_type relative
python histogram.py --input data/sense_R2 --output plot/histogram/sense_R2_substitution.png --variation_type substitution --label_type relative
python histogram.py --input data/sense_R2 --output plot/histogram/sense_R2_insertion.png --variation_type insertion --label_type relative
python histogram.py --input data/sense_R2 --output plot/histogram/sense_R2_deletion.png --variation_type deletion --label_type relative
python histogram.py --input data/dcmv_R1 --output plot/histogram/dcmv_R1_substitution.png --variation_type substitution --label_type relative
python histogram.py --input data/dcmv_R1 --output plot/histogram/dcmv_R1_insertion.png --variation_type insertion --label_type relative
python histogram.py --input data/dcmv_R1 --output plot/histogram/dcmv_R1_deletion.png --variation_type deletion --label_type relative
python histogram.py --input data/dcmv_R2 --output plot/histogram/dcmv_R2_substitution.png --variation_type substitution --label_type relative
python histogram.py --input data/dcmv_R2 --output plot/histogram/dcmv_R2_insertion.png --variation_type insertion --label_type relative
python histogram.py --input data/dcmv_R2 --output plot/histogram/dcmv_R2_deletion.png --variation_type deletion --label_type relative
