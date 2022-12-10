import os

os.system('\n'.join([
  'bowtie2-build-s ref_seq/2DSB_R1_sense.fa data_bowtie2_build/2DSB_R1_sense',
  'bowtie2-build-s ref_seq/2DSB_R2_sense.fa data_bowtie2_build/2DSB_R2_sense',
  'bowtie2-build-s ref_seq/2DSB_R1_cmv.fa data_bowtie2_build/2DSB_R1_cmv',
  'bowtie2-build-s ref_seq/2DSB_R2_cmv.fa data_bowtie2_build/2DSB_R2_cmv',
  'bowtie2-build-s ref_seq/2DSB_R1_branch.fa data_bowtie2_build/2DSB_R1_branch',
  'bowtie2-build-s ref_seq/2DSB_R2_branch.fa data_bowtie2_build/2DSB_R2_branch',
]))

os.system('\n'.join([
  'bowtie2-build-s ref_seq/2DSB_R1_sense.fa data_bowtie2_build/2DSB_R1_sense',
  'bowtie2-build-s ref_seq/2DSB_R2_sense.fa data_bowtie2_build/2DSB_R2_sense',
  'bowtie2-build-s ref_seq/2DSB_R1_cmv.fa data_bowtie2_build/2DSB_R1_cmv',
  'bowtie2-build-s ref_seq/2DSB_R2_cmv.fa data_bowtie2_build/2DSB_R2_cmv',
  'bowtie2-build-s ref_seq/2DSB_R1_branch.fa data_bowtie2_build/2DSB_R1_branch',
  'bowtie2-build-s ref_seq/2DSB_R2_branch.fa data_bowtie2_build/2DSB_R2_branch',
]))