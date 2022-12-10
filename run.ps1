# 1. Alignment.

# 1a. Buld Bowtie2 indexes.
bowtie2-build-s ref_seq/2DSB_R1_sense.fa data/0_bowtie2_build/2DSB_R1_sense
bowtie2-build-s ref_seq/2DSB_R2_sense.fa data/0_bowtie2_build/2DSB_R2_sense
bowtie2-build-s ref_seq/2DSB_R1_cmv.fa data/0_bowtie2_build/2DSB_R1_cmv
bowtie2-build-s ref_seq/2DSB_R2_cmv.fa data/0_bowtie2_build/2DSB_R2_cmv
bowtie2-build-s ref_seq/2DSB_R1_branch.fa data/0_bowtie2_build/2DSB_R1_branch
bowtie2-build-s ref_seq/2DSB_R2_branch.fa data/0_bowtie2_build/2DSB_R2_branch

# 1b. Do the alignment.
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R1_branch data/0_fastq/db1_R1.fq -S data/0_sam/db1_R1.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R2_branch data/0_fastq/db1_R2.fq -S data/0_sam/db1_R2.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R1_branch data/0_fastq/db2_R1.fq -S data/0_sam/db2_R1.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R2_branch data/0_fastq/db2_R2.fq -S data/0_sam/db2_R2.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R1_branch data/0_fastq/db3_R1.fq -S data/0_sam/db3_R1.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R2_branch data/0_fastq/db3_R2.fq -S data/0_sam/db3_R2.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R1_branch data/0_fastq/db4_R1.fq -S data/0_sam/db4_R1.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R2_branch data/0_fastq/db4_R2.fq -S data/0_sam/db4_R2.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R1_cmv data/0_fastq/dcmv1_R1.fq -S data/0_sam/dcmv1_R1.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R2_cmv data/0_fastq/dcmv1_R2.fq -S data/0_sam/dcmv1_R2.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R1_cmv data/0_fastq/dcmv2_R1.fq -S data/0_sam/dcmv2_R1.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R2_cmv data/0_fastq/dcmv2_R2.fq -S data/0_sam/dcmv2_R2.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R1_cmv data/0_fastq/dcmv3_R1.fq -S data/0_sam/dcmv3_R1.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R2_cmv data/0_fastq/dcmv3_R2.fq -S data/0_sam/dcmv3_R2.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R1_cmv data/0_fastq/dcmv4_R1.fq -S data/0_sam/dcmv4_R1.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R2_cmv data/0_fastq/dcmv4_R2.fq -S data/0_sam/dcmv4_R2.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R1_sense data/0_fastq/sense1_R1.fq -S data/0_sam/sense1_R1.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R2_sense data/0_fastq/sense1_R2.fq -S data/0_sam/sense1_R2.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R1_sense data/0_fastq/sense2_R1.fq -S data/0_sam/sense2_R1.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R2_sense data/0_fastq/sense2_R2.fq -S data/0_sam/sense2_R2.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R1_sense data/0_fastq/sense3_R1.fq -S data/0_sam/sense3_R1.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R2_sense data/0_fastq/sense3_R2.fq -S data/0_sam/sense3_R2.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R1_sense data/0_fastq/sense4_R1.fq -S data/0_sam/sense4_R1.sam
bowtie2-align-s -x data/0_bowtie2_build/2DSB_R2_sense data/0_fastq/sense4_R2.fq -S data/0_sam/sense4_R2.sam

# 2. NHEJ filtering.
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/db1_R1.sam --ref_seq_file ../ref_seq/2DSB_R1_branch.fa --output data_1_filter_nhej/db1_R1.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/db1_R2.sam --ref_seq_file ../ref_seq/2DSB_R2_branch.fa --output data_1_filter_nhej/db1_R2.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/db2_R1.sam --ref_seq_file ../ref_seq/2DSB_R1_branch.fa --output data_1_filter_nhej/db2_R1.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/db2_R2.sam --ref_seq_file ../ref_seq/2DSB_R2_branch.fa --output data_1_filter_nhej/db2_R2.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/db3_R1.sam --ref_seq_file ../ref_seq/2DSB_R1_branch.fa --output data_1_filter_nhej/db3_R1.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/db3_R2.sam --ref_seq_file ../ref_seq/2DSB_R2_branch.fa --output data_1_filter_nhej/db3_R2.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/db4_R1.sam --ref_seq_file ../ref_seq/2DSB_R1_branch.fa --output data_1_filter_nhej/db4_R1.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/db4_R2.sam --ref_seq_file ../ref_seq/2DSB_R2_branch.fa --output data_1_filter_nhej/db4_R2.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/dcmv1_R1.sam --ref_seq_file ../ref_seq/2DSB_R1_cmv.fa --output data_1_filter_nhej/dcmv1_R1.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/dcmv1_R2.sam --ref_seq_file ../ref_seq/2DSB_R2_cmv.fa --output data_1_filter_nhej/dcmv1_R2.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/dcmv2_R1.sam --ref_seq_file ../ref_seq/2DSB_R1_cmv.fa --output data_1_filter_nhej/dcmv2_R1.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/dcmv2_R2.sam --ref_seq_file ../ref_seq/2DSB_R2_cmv.fa --output data_1_filter_nhej/dcmv2_R2.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/dcmv3_R1.sam --ref_seq_file ../ref_seq/2DSB_R1_cmv.fa --output data_1_filter_nhej/dcmv3_R1.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/dcmv3_R2.sam --ref_seq_file ../ref_seq/2DSB_R2_cmv.fa --output data_1_filter_nhej/dcmv3_R2.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/dcmv4_R1.sam --ref_seq_file ../ref_seq/2DSB_R1_cmv.fa --output data_1_filter_nhej/dcmv4_R1.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/dcmv4_R2.sam --ref_seq_file ../ref_seq/2DSB_R2_cmv.fa --output data_1_filter_nhej/dcmv4_R2.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/sense1_R1.sam --ref_seq_file ../ref_seq/2DSB_R1_sense.fa --output data_1_filter_nhej/sense1_R1.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/sense1_R2.sam --ref_seq_file ../ref_seq/2DSB_R2_sense.fa --output data_1_filter_nhej/sense1_R2.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/sense2_R1.sam --ref_seq_file ../ref_seq/2DSB_R1_sense.fa --output data_1_filter_nhej/sense2_R1.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/sense2_R2.sam --ref_seq_file ../ref_seq/2DSB_R2_sense.fa --output data_1_filter_nhej/sense2_R2.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/sense3_R1.sam --ref_seq_file ../ref_seq/2DSB_R1_sense.fa --output data_1_filter_nhej/sense3_R1.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/sense3_R2.sam --ref_seq_file ../ref_seq/2DSB_R2_sense.fa --output data_1_filter_nhej/sense3_R2.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/sense4_R1.sam --ref_seq_file ../ref_seq/2DSB_R1_sense.fa --output data_1_filter_nhej/sense4_R1.tsv --min_length 50 --dsb_pos 50 --quiet
python ../1_process_nhej/filter_nhej.py --sam_file data_0_sam/sense4_R2.sam --ref_seq_file ../ref_seq/2DSB_R2_sense.fa --output data_1_filter_nhej/sense4_R2.tsv --min_length 50 --dsb_pos 50 --quiet

# 3. Combine repeats.
python ../1_process_nhej/combine_repeat.py --input data_1_filter_nhej/db1_R1.tsv data_1_filter_nhej/db2_R1.tsv data_1_filter_nhej/db3_R1.tsv data_1_filter_nhej/db4_R1.tsv --output data_2_combine_repeat/db_R1.tsv --quiet
python ../1_process_nhej/combine_repeat.py --input data_1_filter_nhej/db1_R2.tsv data_1_filter_nhej/db2_R2.tsv data_1_filter_nhej/db3_R2.tsv data_1_filter_nhej/db4_R2.tsv --output data_2_combine_repeat/db_R2.tsv --quiet
python ../1_process_nhej/combine_repeat.py --input data_1_filter_nhej/dcmv1_R1.tsv data_1_filter_nhej/dcmv2_R1.tsv data_1_filter_nhej/dcmv3_R1.tsv data_1_filter_nhej/dcmv4_R1.tsv  --output data_2_combine_repeat/dcmv_R1.tsv --quiet
python ../1_process_nhej/combine_repeat.py --input data_1_filter_nhej/dcmv1_R2.tsv data_1_filter_nhej/dcmv2_R2.tsv data_1_filter_nhej/dcmv3_R2.tsv data_1_filter_nhej/dcmv4_R2.tsv --output data_2_combine_repeat/dcmv_R2.tsv --quiet
python ../1_process_nhej/combine_repeat.py --input data_1_filter_nhej/sense1_R1.tsv data_1_filter_nhej/sense2_R1.tsv data_1_filter_nhej/sense3_R1.tsv data_1_filter_nhej/sense4_R1.tsv --output data_2_combine_repeat/sense_R1.tsv --quiet
python ../1_process_nhej/combine_repeat.py --input data_1_filter_nhej/sense1_R2.tsv data_1_filter_nhej/sense2_R2.tsv data_1_filter_nhej/sense3_R2.tsv data_1_filter_nhej/sense4_R2.tsv --output data_2_combine_repeat/sense_R2.tsv --quiet

# 4. Extract windows.
python ../2_get_window_data/get_window.py --input data_2_combine_repeat/db_R1.tsv --ref_seq_file ../ref_seq/2DSB_R1_branch.fa --output data_3_window/db_R1 --dsb_pos 50 --dsb_type 2DSB --strand R1 --guide_rna sgAB --cell_line WT --construct branch --subst_type withSubst --control_type notControl --version versionNone
python ../2_get_window_data/get_window.py --input data_2_combine_repeat/db_R1.tsv --ref_seq_file ../ref_seq/2DSB_R1_branch.fa --output data_3_window/db_R1 --dsb_pos 50 --dsb_type 2DSB --strand R1 --guide_rna sgAB --cell_line WT --construct branch --subst_type withoutSubst --control_type notControl --version versionNone
python ../2_get_window_data/get_window.py --input data_2_combine_repeat/db_R2.tsv --ref_seq_file ../ref_seq/2DSB_R2_branch.fa --output data_3_window/db_R2 --dsb_pos 50 --dsb_type 2DSB --strand R2 --guide_rna sgAB --cell_line WT --construct branch --subst_type withSubst --control_type notControl --version versionNone
python ../2_get_window_data/get_window.py --input data_2_combine_repeat/db_R2.tsv --ref_seq_file ../ref_seq/2DSB_R2_branch.fa --output data_3_window/db_R2 --dsb_pos 50 --dsb_type 2DSB --strand R2 --guide_rna sgAB --cell_line WT --construct branch --subst_type withoutSubst --control_type notControl --version versionNone
python ../2_get_window_data/get_window.py --input data_2_combine_repeat/dcmv_R1.tsv --ref_seq_file ../ref_seq/2DSB_R1_cmv.fa --output data_3_window/dcmv_R1 --dsb_pos 50 --dsb_type 2DSB --strand R1 --guide_rna sgAB --cell_line WT --construct cmv --subst_type withSubst --control_type notControl --version versionNone
python ../2_get_window_data/get_window.py --input data_2_combine_repeat/dcmv_R1.tsv --ref_seq_file ../ref_seq/2DSB_R1_cmv.fa --output data_3_window/dcmv_R1 --dsb_pos 50 --dsb_type 2DSB --strand R1 --guide_rna sgAB --cell_line WT --construct cmv --subst_type withoutSubst --control_type notControl --version versionNone
python ../2_get_window_data/get_window.py --input data_2_combine_repeat/dcmv_R2.tsv --ref_seq_file ../ref_seq/2DSB_R2_cmv.fa --output data_3_window/dcmv_R2 --dsb_pos 50 --dsb_type 2DSB --strand R2 --guide_rna sgAB --cell_line WT --construct cmv --subst_type withSubst --control_type notControl --version versionNone
python ../2_get_window_data/get_window.py --input data_2_combine_repeat/dcmv_R2.tsv --ref_seq_file ../ref_seq/2DSB_R2_cmv.fa --output data_3_window/dcmv_R2 --dsb_pos 50 --dsb_type 2DSB --strand R2 --guide_rna sgAB --cell_line WT --construct cmv --subst_type withoutSubst --control_type notControl --version versionNone
python ../2_get_window_data/get_window.py --input data_2_combine_repeat/sense_R1.tsv --ref_seq_file ../ref_seq/2DSB_R1_sense.fa --output data_3_window/sense_R1 --dsb_pos 50 --dsb_type 2DSB --strand R1 --guide_rna sgAB --cell_line WT --construct sense --subst_type withSubst --control_type notControl --version versionNone
python ../2_get_window_data/get_window.py --input data_2_combine_repeat/sense_R1.tsv --ref_seq_file ../ref_seq/2DSB_R1_sense.fa --output data_3_window/sense_R1 --dsb_pos 50 --dsb_type 2DSB --strand R1 --guide_rna sgAB --cell_line WT --construct sense --subst_type withoutSubst --control_type notControl --version versionNone
python ../2_get_window_data/get_window.py --input data_2_combine_repeat/sense_R2.tsv --ref_seq_file ../ref_seq/2DSB_R2_sense.fa --output data_3_window/sense_R2 --dsb_pos 50 --dsb_type 2DSB --strand R2 --guide_rna sgAB --cell_line WT --construct sense --subst_type withSubst --control_type notControl --version versionNone
python ../2_get_window_data/get_window.py --input data_2_combine_repeat/sense_R2.tsv --ref_seq_file ../ref_seq/2DSB_R2_sense.fa --output data_3_window/sense_R2 --dsb_pos 50 --dsb_type 2DSB --strand R2 --guide_rna sgAB --cell_line WT --construct sense --subst_type withoutSubst --control_type notControl --version versionNone

# 5. Get frequencies.
python ../2_get_window_data/get_freq.py --input data_3_window/db_R1 --output data_3_window/db_R1 --subst_type withSubst --total_reads 3000 3000 3000 3000
python ../2_get_window_data/get_freq.py --input data_3_window/db_R1 --output data_3_window/db_R1 --subst_type withoutSubst --total_reads 3000 3000 3000 3000
python ../2_get_window_data/get_freq.py --input data_3_window/db_R2 --output data_3_window/db_R2 --subst_type withSubst --total_reads 3000 3000 3000 3000
python ../2_get_window_data/get_freq.py --input data_3_window/db_R2 --output data_3_window/db_R2 --subst_type withoutSubst --total_reads 3000 3000 3000 3000
python ../2_get_window_data/get_freq.py --input data_3_window/dcmv_R1 --output data_3_window/dcmv_R1 --subst_type withSubst --total_reads 3000 3000 3000 3000
python ../2_get_window_data/get_freq.py --input data_3_window/dcmv_R1 --output data_3_window/dcmv_R1 --subst_type withoutSubst --total_reads 3000 3000 3000 3000
python ../2_get_window_data/get_freq.py --input data_3_window/dcmv_R2 --output data_3_window/dcmv_R2 --subst_type withSubst --total_reads 3000 3000 3000 3000
python ../2_get_window_data/get_freq.py --input data_3_window/dcmv_R2 --output data_3_window/dcmv_R2 --subst_type withoutSubst --total_reads 3000 3000 3000 3000
python ../2_get_window_data/get_freq.py --input data_3_window/sense_R1 --output data_3_window/sense_R1 --subst_type withSubst --total_reads 3000 3000 3000 3000
python ../2_get_window_data/get_freq.py --input data_3_window/sense_R1 --output data_3_window/sense_R1 --subst_type withoutSubst --total_reads 3000 3000 3000 3000
python ../2_get_window_data/get_freq.py --input data_3_window/sense_R2 --output data_3_window/sense_R2 --subst_type withSubst --total_reads 3000 3000 3000 3000
python ../2_get_window_data/get_freq.py --input data_3_window/sense_R2 --output data_3_window/sense_R2 --subst_type withoutSubst --total_reads 3000 3000 3000 3000

# 6. Get graph data.
python ../3_get_graph_data/get_graph_data.py --input data_3_window/db_R1 --output data_4_graph/db_R1 --subst_type withSubst
python ../3_get_graph_data/get_graph_data.py --input data_3_window/db_R1 --output data_4_graph/db_R1 --subst_type withoutSubst
python ../3_get_graph_data/get_graph_data.py --input data_3_window/db_R2 --output data_4_graph/db_R2 --subst_type withSubst
python ../3_get_graph_data/get_graph_data.py --input data_3_window/db_R2 --output data_4_graph/db_R2 --subst_type withoutSubst
python ../3_get_graph_data/get_graph_data.py --input data_3_window/dcmv_R1 --output data_4_graph/dcmv_R1 --subst_type withSubst
python ../3_get_graph_data/get_graph_data.py --input data_3_window/dcmv_R1 --output data_4_graph/dcmv_R1 --subst_type withoutSubst
python ../3_get_graph_data/get_graph_data.py --input data_3_window/dcmv_R2 --output data_4_graph/dcmv_R2 --subst_type withSubst
python ../3_get_graph_data/get_graph_data.py --input data_3_window/dcmv_R2 --output data_4_graph/dcmv_R2 --subst_type withoutSubst
python ../3_get_graph_data/get_graph_data.py --input data_3_window/sense_R1 --output data_4_graph/sense_R1 --subst_type withSubst
python ../3_get_graph_data/get_graph_data.py --input data_3_window/sense_R1 --output data_4_graph/sense_R1 --subst_type withoutSubst
python ../3_get_graph_data/get_graph_data.py --input data_3_window/sense_R2 --output data_4_graph/sense_R2 --subst_type withSubst
python ../3_get_graph_data/get_graph_data.py --input data_3_window/sense_R2 --output data_4_graph/sense_R2 --subst_type withoutSubst

# 7. Get histogram data.
python ../4_get_histogram_data/get_histogram_data.py --input data_4_graph/db_R1 --output data_5_histogram/db_R1 --subst_type withSubst
python ../4_get_histogram_data/get_histogram_data.py --input data_4_graph/db_R1 --output data_5_histogram/db_R1 --subst_type withoutSubst
python ../4_get_histogram_data/get_histogram_data.py --input data_4_graph/db_R2 --output data_5_histogram/db_R2 --subst_type withSubst
python ../4_get_histogram_data/get_histogram_data.py --input data_4_graph/db_R2 --output data_5_histogram/db_R2 --subst_type withoutSubst
python ../4_get_histogram_data/get_histogram_data.py --input data_4_graph/dcmv_R1 --output data_5_histogram/dcmv_R1 --subst_type withSubst
python ../4_get_histogram_data/get_histogram_data.py --input data_4_graph/dcmv_R1 --output data_5_histogram/dcmv_R1 --subst_type withoutSubst
python ../4_get_histogram_data/get_histogram_data.py --input data_4_graph/dcmv_R2 --output data_5_histogram/dcmv_R2 --subst_type withSubst
python ../4_get_histogram_data/get_histogram_data.py --input data_4_graph/dcmv_R2 --output data_5_histogram/dcmv_R2 --subst_type withoutSubst
python ../4_get_histogram_data/get_histogram_data.py --input data_4_graph/sense_R1 --output data_5_histogram/sense_R1 --subst_type withSubst
python ../4_get_histogram_data/get_histogram_data.py --input data_4_graph/sense_R1 --output data_5_histogram/sense_R1 --subst_type withoutSubst
python ../4_get_histogram_data/get_histogram_data.py --input data_4_graph/sense_R2 --output data_5_histogram/sense_R2 --subst_type withSubst
python ../4_get_histogram_data/get_histogram_data.py --input data_4_graph/sense_R2 --output data_5_histogram/sense_R2 --subst_type withoutSubst

# 8. Plot variation-distance graphs.
python ../5_plot_graph/plot_graph.py --input data_4_graph/db_R1 --output plot/graph --ext png --layout universal  --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python ../5_plot_graph/plot_graph.py --input data_4_graph/db_R2 --output plot/graph --ext png --layout universal  --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python ../5_plot_graph/plot_graph.py --input data_4_graph/dcmv_R1 --output plot/graph --ext png --layout universal  --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python ../5_plot_graph/plot_graph.py --input data_4_graph/dcmv_R2 --output plot/graph --ext png --layout universal  --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python ../5_plot_graph/plot_graph.py --input data_4_graph/sense_R1 --output plot/graph --ext png --layout universal  --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python ../5_plot_graph/plot_graph.py --input data_4_graph/sense_R2 --output plot/graph --ext png --layout universal  --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

# 9. Plot variation-position histograms.
python ../6_plot_histogram/plot_histogram.py --input data_5_histogram/db_R1 --output plot/histogram/ --label_type relative
python ../6_plot_histogram/plot_histogram.py --input data_5_histogram/db_R2 --output plot/histogram/ --label_type relative
python ../6_plot_histogram/plot_histogram.py --input data_5_histogram/dcmv_R1 --output plot/histogram/ --label_type relative
python ../6_plot_histogram/plot_histogram.py --input data_5_histogram/dcmv_R2 --output plot/histogram/ --label_type relative
python ../6_plot_histogram/plot_histogram.py --input data_5_histogram/sense_R1 --output plot/histogram/ --label_type relative
python ../6_plot_histogram/plot_histogram.py --input data_5_histogram/sense_R2 --output plot/histogram/ --label_type relative
