# Variation-Distance Graphs

## Introduction

This tool is intended to allow processing and visualizing high-throughput sequencing data obtained for the purpose of studying double-strand break (DSB) repair due to the nonhomologous end-joining (NHEJ) repair mechanism. For a full description of the protocol, please refer to the manuscript at LINK. This protocol was original used in the study Jeon et al. (LINK) for studying DSB repair in human cells. That publication contains several examples of the graphs in the supplementary figures, as well as a discussion of insights gained from the resulting figures.

## Installation

To install the package, use the command
```
pip install XXX
```
The required dependencies are
* XX
* YY
(TODO)

Bowtie 2 (version >= XX) should be installed and available on the PATH. Particularly, the executables `bowtie2-build-s` and `bowtie2-align-s` should be available as commands.

## Tutorial

Several example input files are available in the `data_input` directory:
* `data_input/ref_seq`: reference sequence FASTA files, representing the perfect repaired sequence for different samples.
* `data_input/fastq`: high-throughput sequencing data for different samples. Note, that the FASTQ samples have been *trimmed*, meaning that we only capture the portion of the read between the primers and low-quality reads have been already filtered out.

The reference sequences must be chose to exactly corresponding with the expected sequence between the primer sites. CONTINUE HERE!