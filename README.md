# Variation-Distance Graphs

## Introduction

This tool is intended to allow processing and visualizing high-throughput sequencing data obtained for the purpose of studying double-strand break (DSB) repair due to the nonhomologous end-joining (NHEJ) repair mechanism. The accompanying article can be found at LINK. This protocol was original used in the study Jeon et al. ([link](https://doi.org/10.1101/2022.11.01.514688)) for studying DSB repair in human cells. That publication contains several examples of the graphs in the supplementary figures, as well as a discussion of insights gained from the resulting figures.

The overall functionality of the package is to (1) preprocess DNA-seq reads that have been obtained from DSB (double-strand break) repair experiments; (2) quantify the variations (insertions, deletions, and substitutions) near the DSB site using sequence alignment; (3) visualize the resulting variations using two types of visualizations: *variation-distance graphs* and *variation-position histograms*.

The expected inputs are DNA-seq libraries that have been treated in a very specific manner to properly allowing infering the variation with alignments. For a more detailed description of the expected input, see section [Commands](#commands) and the publication associated with this software (LINK).

## Terminology

* We use the terms "vertex" and "node" interchageably for the vertices of graphs.
* We sometimes also use the term "sequence" and "alignment" interchangeably with "vertex" and "node", since, as we describe later, there is a one-to-one correspondence between these four things in the constructed graphs. For example all of the following statements may be equivalent when we refer to a single vertex in a graph:
    * The vertex has two insertions.
    * The node has two insertions.
    * The sequence has two insertions.
    * The alignment has two insertions.

## Citation

If you use this software, please use the following citation:

XXX

or in BibTex format:

XXX

## Installation

To install the package, use the command
```
pip install XXX
```
The required dependencies are
* Plotly
* NetworkX
* PIL
* python-pptx
* numpy
* pandas
* SciKit-Learn
* Levenshtein
(TODO)

Bowtie 2 (version >= XX) should be installed and available on the PATH. Particularly, the executables `bowtie2-build-s` and `bowtie2-align-s` should be available as commands.

## Commands

This package is based on the following four commands:

* `preprocess.py`: takes as input the trimmed FASTQ reads files and a FASTA file containing the reference sequence, and creates the intermediate tables needed for plotting the graphs and histograms. The input directories should represent replicate experiments (e.g., biological replicates).
* `comparison.py`: takes as input two directories created by `preprocess.py` and creates a directory that contains the intermediate tables needed for plotting comparison graphs. Both of the directories must have been created with the same reference sequence.
* `graph.py`: takes as input a collection of the output directories of either `preprocess.py` or `comparison.py`, lays out alignment-sequences in all inputs, and plots a separate graph for each input.
* `histogram.py`: takes as input an output directories of `preprocess.py` (outputs of `comparison.py` are not valid) and plots a histogram showing the type and position of variations (insertions, deletion, or substitutions) in the alignment-sequences.

More information about each command is given in the following subsections. The the exposition, we use the notation `NAME` to refer to the value of a command-line parameter set by the user and the notation `--name` to refer to the parameter itself. The notation `file_name.ext` is also used to refer to file names, though the meaning should be clear from context.

### `preprocess.py`

The `preprocess.py` script perform alignment and preprocessing of the raw FASTQ data.

#### Input

We expect that the input FASTQ files to `preprocess.py` are *trimmed*, meaning that the adaptors have been removed. We also expected that the region of DNA between these adaptors is exactly the region of DNA represented by the reference sequence. This mean that if a given read represents a perfectly repaired DNA strand, it should identical with the reference sequence (assuming no substitution errors due to library preparation or sequencing). If multiple input FASTQ files are given, it is assumed that they are repeats of the same treatment condition and are processed identically then combined into a single file (see [preprocessing stages](#prepcocessing-stages) below).

#### Ignoring substitutions

The preprocessing pipeline produces two different versions of most files: one *ignoring* substitutions (suffix "withoutSubst") and another *keeping* substitutions (suffix "withSubst"). The processing for files that ignore substitutions contains an extra step that replaces alginment subsititions (mismatches) with perfect matches with the reference sequence. We chose to ignore substitutions in our analysis because we noticed a consistent distribution of substitutions occurring in both the experiment group (where DNA double-strand breaks were induced) and the control group (where no DSBs were induced). This suggests that the majority of substitutions were likely caused by DNA damage during library preparation or sequencing errors, rather than the process of repairing double-strand breaks. In the command `graph.py`, the `--subst_type` parameter controls whether to use the output with or without substitutions. The `histogram.py` command only uses the output with substitutions, since it is partly used to examine the distribution of substitutions.

#### Preprocessing stages

This `preprocess.py` script is broken in separate stages so that each stage can be run separately (potentially on different machines). However, the stages must be run in the correct order indicated by their prefixes. If running the stages separately, the value of the `OUTPUT` directory must the same value on each separate invocation. Two different experiments should not be given the same `OUTPUT` directory or the data from the second will overwrite the first. The following describes each stage in more detail.

1. **0_align**: Align FASTQ reads against FASTA reference sequence. This stage requires multiple input FASTQ files representing independent repeats or biological/tehchnical replicates of the treatment condition. The alignment is done independently for each of the input FASTQs. The output of this step is a set of [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) files. Please ensure that [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is installed and that the commands `bowtie2-build-s` and `bowtie2-align-s` are available on the PATH (using the Python command `os.system()`). The output directories of this stage are:
    * `0_bowtie2_build`: The Bowtie 2 index files built with `bowtie2-build-s`.
    * `0_sam`: The SAM file output from the alignment by `bowtie2-align-s`.
2. **1_filter**: Filter each SAM file independently using various heuristics to discard alignments that don't represent NHEJ repair. The filtering process involves the following steps:
    1. Discard alignments that represent a failed/invalid aligment in the SAM format.
    2. Discard alignments where the read sequence has length less than `MIN_LENGTH`.
    3. Discard alignments where the left-most (5'-most) position of the read does not align with the left-most position of the reference sequence.
    4. If the in/del positions of the alignment are not adjacent to or around the DSB position, try to shift them towards the DSB position in a way that keeps the number of substitutions roughly the same. If such a modification of the alignment cannot be found, discard the alignment.
    5. Finally, if the resulting alignment has in/dels that are not contiguous, discard the alignment.
   If multiple reads have exactly the same nucleotide sequence but had different alignments with the reference, they will be forced to have the same alignment as the first such read encountered. Note that the left-most (5'-most) position on the read must align with the left-most position on the reference, but the same is not true for the right-most (3'-most) positions. This is because it is possible for reads to be too small and only capture the 5' primer but not the 3' primer. However, the read must have a length of at least `DSB_POS + 1` so that it surrounds the DSB position. The output is a set of tables in TSV format that contain the unique alignments and their read counts. The output is saved in files named `1_filter_nhej/[N].tsv`, where `[N]` is an integer representing the repeat number, in the same order they are specified on the command line.
3. **2_combine**: Combine the tables output from stage 2 into a single table with multiple frequency columns, one for each input. The output will be in `2_combine_repeat/out.tsv`.
4. **3_window**: This stage is broken into the following substages:
    1. Window extraction: For each unique alignment in the output table of stage 2_combine, extract the portion of the read that aligns to the positions `DSB_POS - WINDOW_SIZE + 1` to `DSB_POS + WINDOW_SIZE` on the reference sequence. If different reads become identical after extracting their windows, their read counts will be summed. To ensure that variations near the DSB do not spill outside the window, *anchors* must be present on either side of the extracted windows. Anchors are the parts of the read that align to the `ANCHOR_SIZE` nucleotides on the left (5') and right (3') of the window on the reference. Each anchor must have at most `ANCHOR_MISMATCHES` mismatches and no in/dels, or it is discarded. The output is a table in TSV format. The "withoutSubst" version of the output has substitutions (or mismatches) in the alignment removed by replacing the incorrect nucleotide on the read with the correct nucleotide from the reference. See [here](#substitutions-in-bowtie-2) for the justification of removing substitutions.
    2. Frequency conversion: Convert the raw read counts into frequencies in the range [0,1]. This is done by dividing the read count of each sequence by the total number of reads in the library. Since the total number of reads in each sample might not be exactly the sum of the read counts (e.g., if the trimming process removed some reads), the total reads may be input using `--total_reads`. The output will be in the
    3. Frequency filtering: Sequences whose frequency in any of the repeat is below `FREQ_MIN` will be discarded. Note that this means only sequences present in all libraries are retained. 
    4. Frequency averaging: For each sequence, we obtain a single frequency by taking the mean of the repeat frequencies.
   
   The output will be in the subdirectory `3_window`, in files:
    * `window_count_withoutSubst.tsv` and `window_count_withSubst.tsv`: Count tables with separate columns for the separate repeats.
    * `window_freq_withSubst.tsv`/`window_freq_withoutSubst.tsv`: Frequency tables with separate columns for the separate repeats.
    * `window_freq_filter_withSubst.tsv`/`window_freq_filter_withoutSubst.tsv`: The previous frequency tables with the minimum frequency filter applied.
    * `window_freq_filter_mean_withSubst.tsv`/`window_freq_filter_mean_withoutSubst.tsv`: The previous frequency tables with the separate frequency columns for the repeats collapsed into a single column by taking their mean.
5. **4_graph**: Precompute data needed to plot the graphs, such as adjacency information and summary statistics of the graphs. The output will be in the subdirectory `4_graph`.
6. **5_histogram**: Precompute data needed to plot the histograms, such as the position and frequency of the different types of variations. The output will be in the subdirectory `5_histogram`.

For more information about the parameters, please use the command `python preprocess.py --help`. For more details about the preprocessing, please refer to CITATION.

### `comparison.py`

This stages creates data needed to compare two different sequencing libraries. Specifically, it takes as input two directories that were created using the `preprocess.py` script, and will output a third directory that contains analogous data for creating comparison graphs of the two samples. For meaningful comparisons, it is important that the reference sequences (after restricting to the specified window) are identical between the two libraries. The output data will contain the same columns as that output from `preprocess.py` (i.e., it will have the same information about the variations around the DSBs), but will have two additional columns (with suffixes `_1` and `_2`) for the frequencies of the two different samples being compared. The original frequency column, with no suffix, will contain the maximum of these two frequencies. The output will be in the subdirectories `3_window` and `4_graph`. Note that the histogram data is not recomputed because the comparison data should not be used to create histograms. For more details about the parameters, please use the command `python comparison.py --help`.

### `graph.py`

This command performs the vertex layout and final plotting of the variation-distance graphs. Multiple inputs may be specified to layout the vertices jointly, as long as their windowed reference sequences are identical. See the [Graphs](#graphs) section for a description of the meaning of different visual properties of the graph (e.g., edges, vertex sizes, colors, etc.). For more details about the parameters, please use the command `python graph.py --help`.

### `histogram.py`

This commands plots the 3D histograms, which summarize the distribution of the variations in the windows around the DSB site. Each histogram shows a single type of variation (insertion, deletion, or substitution), which is determined by the `--variation_type` parameter. The axes are described in the following.

* The position of each variation is shown on the $x$-axis of the graphs. The position can be measured in two ways: (1) relative to the position of the DSB, or (2) relative to the 5'-end of the window. Which of these two options is used depends on the value of the `--label_type` parameter. For deletions and substitutions, the position of the variation is unambiguous and corresponds to a specific position on the reference sequence. However, for insertions, the position is assigned to the left (5') position of the two consecutive positions around the DSB on the reference sequence, as the insertion occurs between these two positions.
* The $y$-axis indicates the total number of in/dels on the sequence that the variation originated from.
* The $z$-axis indicates the total frequency of the variations with the given $x$- and $y$-coordinates.

To understand how the frequencies are calculated, consider the following example: Suppose an alignment in the data had one insertion at position 5, two deletions at positions 6 and 7, and an overall frequency of 0.1. Let's assume that the DSB is between positions 5 and 6, and that we are using relative labeling on the $x$-axis. In the insertion histogram, the alignment would contribute +0.1 to the z-value of the bar at the $x$-$y$-coordinate (-1, 3).  In the deletion histogram, the alignment would contribute +0.1 to the z-value of the bars at the $x$-$y$-coordinates (1, 3) and (2, 3). If we used absolute labeling, then the respective $x$-$y$-coordinates would be (5, 3), (6, 3), and (7, 3), since the positions are now labeled according to their absolute positions on the reference sequence.

For more details about the parameters, please use the command `python histogram.py --help`.

## Graphs

In this section, we will explain the visual elements of the output graphs generated by the `graph.py` command. Please note that we use the terms "sequence" and "vertex" interchangeably, as there is a one-to-one correspondence between the windowed sequences produced by the [`preprocess.py`](#preprocesspy) script and the vertices of the graph.

### Layouts

Different *layouts*, specified with the `--layout` parameter, are used to position the vertices. We describe a subset of them in the following.

* Kamada-Kawaii: The Kamada-Kawaii algorithm is a physical simulation algorithm used for laying out graphs (Kamada and Kawai, 1989). In this algorithm, the graph is represented as a spring system where each pair of vertices is connected by a spring whose resting length is proportional to the shortest-path distance between the vertices. The algorithm works by finding the optimal positions of the vertices that minimize the energy of the system. The result is a layout that is both aesthetically pleasing and reveals the graph structure, such as clusters of vertices that are close together. The [NetworkX](#https://networkx.org/) Python package (Hagberg, 2008) implements this algorithm, and as of this writing, it initializes the vertices in a circle. Note that the shortest-path distance between two vertices is at least as large as the Levenshtein distance between them, but it doesn't have to be equal. The final positions of the vertices are always in the range [0, 1], which is important to keep in mind when setting the `--plot_range_x` and `--plot_range_y` parameters.
* Universal layout: We refer to this layout as the "universal" layout because it assigns a predefined position to every possible insertion and deletion, regardless of the subset of insertions and deletions present in the dataset. In contrast to the Kamada-Kawai algorithm, which positions vertices in a data-dependent manner, this layout algorithm is less flexible but highly interpretable and reproducible. This layout also only applies to data *without* substitutions, set with `--subst_type withoutSubst`. The universal layout places the reference vertex at the center of the figure, with insertion vertices above and deletion vertices below. For both insertion and deletion vertices, the $y$-coordinate is determined by the Levenshtein distance between the vertex's sequence and the reference sequence, which is simply the number of nucleotides inserted or deleted in the vertex's sequence. If a vertex's Levenshtein distance is $d$, then deletion vertices are placed $d$ units below the reference, while insertion vertices are placed $d$ units above. The $x$-coordinates of vertices are determined using different rules for insertions and deletions. Insertion vertices are ordered from left to right by alphabetical order, while deletion vertices are ordered from left to right based on the positions of the deleted nucleotides with respect to the DSB site. Figure 1 provides a key to the vertex placement in the universal layout. In this layout, the vertex coordinates typically range between $\pm 2 * window\_size$, where $window\_size$ is a parameter from preprocess.py (this is relevant if you set the `--plot_range_x` and `--plot_range_y` parameters). However, the scaling on the insertion and deletion, $x$ and $y$-coordinates can be controlled by the `--universal_layout_x_scale_insertion`, `--universal_layout_y_scale_insertion`, `--universal_layout_x_scale_deletion`, and `--universal_layout_y_scale_deletion` parameters.
* Radial: This layout arranges the vertices in concentric circles around the reference sequence. The reference vertex is positioned at the center of the figure, while the insertion and deletion vertices are placed above and below the reference, respectively. The Levenshtein distance between a vertex and the reference determines its physical distance from the reference vertex. For instance, vertices with a Levenshtein distance of 1 are placed in the first concentric circle around the reference, and those with a Levenshtein distance of 2 are placed in the second circle, and so on. In each circle, insertion vertices are arranged in a clockwise direction based on their frequency, with the most frequent vertex at the top, while deletion vertices are arranged in a counterclockwise direction, also based on their frequency, with the most frequent vertex at the bottom. To ensure that the vertices do not overlap and edges are not collinear, several heuristics are used to perturb the vertices slightly. The vertex coordinates in this layout usually range between $\pm 4 * window\_size$, where $window\_size$ is the parameter from `preprocess.py` (relevant when setting the `--plot_range_x` and `--plot_range_y` parameters).

![Universal layout vertex placement key](figures/graph_layout_scheme.png)
*Figure 1: Universal layout vertex placement key.*

### Graph aesthetics

In this section, we describe the visual feature of the graphs. Not all parameters are described here; for a full description please use `python graph.py --help`.

* Vertices:
  * Size: For an individual graph, the vertex radius is a function of the log-frequency. This is controlled by the parameters `--node_freq_range` and `--node_px_range`. Vertices with frequencies $\leq$ `NODE_FREQ_RANGE[0]` get `NODE_PX_RANGE[0]` and frequencies $\geq$ `NODE_FREQ_RANGE[1]` get `NODE_PX_RANGE[1]`. Frequencies in between get a radius that varies linearly with the log-frequency.
  * Outline: The outline color is intended to differentiate the reference sequence vertex from the other vertices. The colors are determined by the `--node_reference_outline_color` and `--node_outline_color` parameters.
  * Color: The vertex color is determined by the following modes:
    * Variation type: This mode is only for individual graphs. This mode gives vertices a color depending on the type of variations present in their alignment: insertions only, deletions only, substitutions only, mixed (which is a combinations of insertions, deletions, and substitution), and none (which represents the reference sequence). These colors can be defined by using the `--variation_type_colors` parameter.
    * Frequency ratio (continuous): This mode is only used for comparison graphs, and can be activated by setting the parameter `--node_comparison_color_type continuous`. In this mode, a color gradient is used for each vertex to indicate the frequency ratio of the vertex between the two experiments being compared. The ratio is calculated by putting the vertex's frequency in the first sample in the numerator and the vertex's frequency in the second sample in the denominator (the order is determined during preprocessing). The colors at the two ends of the gradient are specified by the `--node_comparison_colors` parameter, while the colors in between are smoothly interpolated based on the frequency ratio. The range of ratios at the two ends of the gradient is determined by the `--node_freq_ratio_range` parameter.
    * Frequency ratio (discrete): This mode is only used for comparison graphs, and can be activated by setting the parameter `--node_comparison_color_type discrete`. In this mode, three colors are used to indicate the frequency ratio of the vertex between the two experiments being compared. If the frequency ratio is less than `NODE_FREQ_RATIO_RANGE[0]`, the color `NODE_COMPARISON_COLORS[1]` is displayed (meaning higher in sample 2). If the frequency ratio is greater than `NODE_FREQ_RATIO_RANGE[1]`, the color `NODE_FREQ_RATIO_RANGE[1]` is displayed (meaning higher in sample 1). If the frequency ratio is between `NODE_FREQ_RATIO_RANGE[0]` and `NODE_FREQ_RATIO_RANGE[1]` (inclusive), the vertex is colored white.
* Edges: All edges in the graph indicate a 1-nucleotide variation (insertion, deletion, or substitution) between the two vertices that are connected. The display of edges and which types of edges to show can be controlled by using the `--edge_show` and `--edge_types` parameters.
* Legends: The legends describe the vertex size, vertex outline, vertex color, and the edges. These can be all drawn by using the `--legend` parameter. The legends will be drawn in the right margin of the figure. To ensure enough room, use the `--margin_top_px`, `--margin_bottom_px`, `--margin_left_px`, and `--margin_right_px` parameters. The different legends are laid out vertically. To control the spacing between them use the `--legend_spacing_px` parameter.
* Title: A title can be optionally added to the top margin of the figure using the `--title` parameter.
* Universal layout axes: To illustrate the vertex placement for the universal layout, axes can be drawn on the figure. Using the `--universal_layout_x_axis_deletion_y_pos` parameter, a horizontal axis can be drawn on the deletion side (below the reference sequence), which, for each deletion vertex, shows the approximate position of the range of deleted nucleotides (see [Figure 1](#graphs)) for mode details. Using the `--universal_layout_x_axis_insertion_y_pos` parameter, a horizontal axis can be drawn on the insertion side (above the reference sequence) to show the alphabetical order of the nucleotides. Using the `--universal_layout_y_axis_x_pos` parameter, a vertical axis can be drawn to show the Levenshtein distance of vertices from the reference vertex. There are several other parameters to control various aspects of the axes: `--universal_layout_x_axis_deletion_label_type`, `--universal_layout_y_axis_y_range`, `--universal_layout_x_axis_x_range`, `--universal_layout_y_axis_deletion_max_tick`, and `--universal_layout_y_axis_insertion_max_tick`. The following are used to scale the axes by different amount, mostly for the purpose of aesthetics: `--universal_layout_x_scale_insertion`, `--universal_layout_y_scale_insertion`, `--universal_layout_x_scale_deletion`, and `--universal_layout_y_scale_deletion`. To help positioning the axes, a message is printed to the console showing the range of $x$- and $y$-coordinates for each figure, when `graph.py` is run.
* Graph statistics: The `--stats` parameter can be used to display summary statistics of the graph in the left margin of the figure. If the graph has more than one connected component, the component containing the reference sequence is the only one used since the pairwise distance between vertices is not well-defined otherwise. To ensure that the left margin has enough space to display the summary statistics, the `--margin_left_px` parameter can be used.

## Demonstration

Several example input files are available in the `data_input` directory:

* `data_demo/ref_seq`: reference sequence FASTA files, representing the perfect repaired sequence for different experiments.
* `data_demo/fastq`: high-throughput sequencing data in FASTQ format for different samples. Each FASTQ file corresponds to an independent experiment.

Note that the FASTQ samples have been *trimmed*, meaning that we only capture the portion of the read between the adaptors, and low-quality reads have been already filtered out.

### Demonstration: Preprocessing

The [`preprocess.py`](#preprocesspy) command must be run before producing any figures. This stage aligns the input FASTQ files against the reference sequences, extracts the part of the alignment around the DSB position, and creates tables of the unique alignments and their frequencies. The following are examples of using the command:

```
python preprocess.py --input data_input/fastq/sense1_R1.fq data_input/fastq/sense2_R1.fq data_input/fastq/sense3_R1.fq data_input/fastq/sense4_R1.fq --ref_seq_file data_input/ref_seq/2DSB_R1_sense.fa --dsb_pos 67 --output data_output/sense_R1 --label sense_R1 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_input/fastq/sense1_R2.fq data_input/fastq/sense2_R2.fq data_input/fastq/sense3_R2.fq data_input/fastq/sense4_R2.fq --ref_seq_file data_input/ref_seq/2DSB_R2_sense.fa --dsb_pos 46 --output data_output/sense_R2 --label sense_R2 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_input/fastq/db1_R1.fq data_input/fastq/db2_R1.fq data_input/fastq/db3_R1.fq data_input/fastq/db4_R1.fq --ref_seq_file data_input/ref_seq/2DSB_R1_branch.fa --dsb_pos 67 --output data_output/db_R1 --label db_R1 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_input/fastq/db1_R2.fq data_input/fastq/db2_R2.fq data_input/fastq/db3_R2.fq data_input/fastq/db4_R2.fq --ref_seq_file data_input/ref_seq/2DSB_R2_branch.fa --dsb_pos 46 --output data_output/db_R2 --label db_R2 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_input/fastq/dcmv1_R1.fq data_input/fastq/dcmv2_R1.fq data_input/fastq/dcmv3_R1.fq data_input/fastq/dcmv4_R1.fq --ref_seq_file data_input/ref_seq/2DSB_R1_cmv.fa --dsb_pos 67 --output data_output/dcmv_R1 --label dcmv_R1 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_input/fastq/dcmv1_R2.fq data_input/fastq/dcmv2_R2.fq data_input/fastq/dcmv3_R2.fq data_input/fastq/dcmv4_R2.fq --ref_seq_file data_input/ref_seq/2DSB_R2_cmv.fa --dsb_pos 46 --output data_output/dcmv_R2 --label dcmv_R2 --total_reads 3000 3000 3000 3000
```

There are 24 files in the `data_input/fastq` directory, but we only run the preprocessing pipeline six times. This is because we process the biological repeats of each experiment together. For instance, `sense1_R1.fq`, `sense2_R1.fq`, `sense3_R1.fq`, and `sense4_R1.fq` are all biological repeats of the "Sense, 2-DSB" experiment. The `R1` or `R2` in the file names indicate whether the forward or reverse strand of the read was sequenced.

It's important that all the biological repeats for each experiment use the same reference sequence. For example, the four "Sense, 2-DSB" experiments in the examples use the `2DSB_R1_sense.fa` reference sequence.

The output of the preprocessing pipeline will be a collection of tables in TSV (tab-separated value) format. The most significant files are located in the  `4_graph` and `5_histogram` subdirectories, which contain data used for plotting graphs and histograms, respectively. Other intermediate preprocessing files are stored in other subdirectories.

An alternative to running all the stages of preprocessing in a single `preprocess.py` invocation is to run the stages separately using the `--stages` command. Below is an example of preprocessing the db_R1 experiment with only one stage per invocation. Note that not all parameters must be set for all the stages.

```
python preprocess.py --input data_demo/input/fastq/db1_R1.fq data_demo/input/fastq/db2_R1.fq data_demo/input/fastq/db3_R1.fq data_demo/input/fastq/db4_R1.fq --ref_seq_file data_demo/input/ref_seq/2DSB_R1_branch.fa --output data_demo/output/db_R1 --stages 0_align
python preprocess.py --ref_seq_file data_demo/input/ref_seq/2DSB_R1_branch.fa --dsb_pos 67 --output data_demo/output/db_R1 --stages 1_filter
python preprocess.py --output data_demo/output/db_R1 --stages 2_combine
python preprocess.py --ref_seq_file data_demo/input/ref_seq/2DSB_R1_branch.fa --dsb_pos 67 --output data_demo/output/db_R1 --label db_R1 --total_reads 3000 3000 3000 3000 --stages 3_window
python preprocess.py --output data_demo/output/db_R1 --stages 4_graph
python preprocess.py --output data_demo/output/db_R1 --stages 5_histogram
```

If you need more information about the parameters used in the preprocessing pipeline, you can use the command `python preprocess.py --help`. You can also see the [`preprocess.py`](#preprocesspy) section for more details.

### Demonstration: Comparisons

The [comparison](#comparisonpy) script is used to combine two outputs directories of the [preprocessing](#preprocesspy) command for plotting [comparison graphs](#XX). The following are examples:

```
python comparison.py --input data_output/sense_R1 data_output/db_R1 --output data_output/sense_db_R1
python comparison.py --input data_output/sense_R2 data_output/db_R2 --output data_output/sense_db_R2
python comparison.py --input data_output/sense_R1 data_output/dcmv_R1 --output data_output/sense_dcmv_R1
python comparison.py --input data_output/sense_R2 data_output/dcmv_R2 --output data_output/sense_dcmv_R2
```

Note that both experiments must have identical reference sequences after extracting its `2 * WINDOW_SIZE` nucleotides around the DSB position. The whole reference sequences (i.e., the files in the `--ref_seq` parameter to `preprocess.py`) need not be identical so long as this window is.

The output will have a directory structure mirroring the output of `preprocess.py`, except that only the subdirectories `3_window`, `4_graph`, and `5_histogram` will be present (since there is no need to cache the other preprocessing data). The output files contain a combination of the data from both the experiments for creating comparison graphs.

### Demonstration: Graphing

The [`graph.py`](#graphpy) script is used to layout and plot the aligned reads output from the [`preprocess.py`](#preprocesspy) and [`comparison.py`] scripts, into. The input should be a list of directories created by `preprocess.py`/`comparison.py` and the output will be PNG or HTML files that visually represent the graph of sequences extracted from reads around the DSB position. The graphs are represented as a scatter-plot of vertices (representing sequences) and edges (representing single-variation sequence similarity).

Multiple inputs may be specified as long as all the experiments have identical reference sequences (potentially after reverse-complementing; see the `--reverse_complement` parameter) within the DSB window. When multiple inputs are specified, all sequences (potentially after reverse complementing) in all the inputs are combined into a single "supergraph". The supergraph is laid out to determine the coordinates of each vertex. Then the graphs of each input (which are subgraphs of the supergraph) are laid out separately according to the coordinates determined by the supergraph. This allows the output figures for each experiment to have the same coordinates for the same sequence, which aids visual comparison between experiments.

The output may be either PNG or HTML. HTML output may be opened in a browser and allows interactive inspection of the vertices by hovering over them.

To plot individual graphs (graphs with a single experiment), the output directories of `preprocess.py` should be used. To plot comparison graphs (graphs comparing two experiments), the output directories of `comparison.py` should be used.

There are many parameters for customizing the aesthetics of the output graphs, which are covered in more detail in section [Graphs](#graphs).

The following are example commands of using `graph.py` with a single input directory:

```
python graph.py --input data_output/db_R1 --output plot/graph/universal/db_R1.png --legend --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5  --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/db_R2 --output plot/graph/universal/db_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/sense_R1 --output plot/graph/universal/sense_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/sense_R2 --output plot/graph/universal/sense_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/dcmv_R1 --output plot/graph/universal/dcmv_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/dcmv_R2 --output plot/graph/universal/dcmv_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

```

The following is an example of using `graph.py` with multiple input directories:

```
python graph.py `
  --input data_output/db_R1 data_output/db_R2 data_output/sense_R1 data_output/sense_R2 data_output/dcmv_R1 data_output/dcmv_R2 data_output/sense_db_R1 data_output/sense_db_R2 data_output/sense_dcmv_R1 data_output/sense_dcmv_R2 `
  --output plot/graph/kamada_common/db_R1.png plot/graph/kamada_common/db_R2.png plot/graph/kamada_common/sense_R1.png plot/graph/kamada_common/sense_R2.png plot/graph/kamada_common/dcmv_R1.png plot/graph/kamada_common/dcmv_R2.png plot/graph/kamada_common/sense_db_R1.png plot/graph/kamada_common/sense_db_R2.png plot/graph/kamada_common/sense_dcmv_R1.png plot/graph/kamada_common/sense_dcmv_R2.png `
  --reverse_complement 0 1 0 1 0 1 0 1 0 1 --layout kamada_layout --width 2400 --height 1800
```

The following are examples of plotting comparison graphs:

```
python graph.py --input data_output/sense_db_R1 --node_comparison_colors "#cf191b" "#33a02c" --output plot/graph/universal/sense_db_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/sense_db_R2 --node_comparison_colors "#cf191b" "#33a02c" --output plot/graph/universal/sense_db_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/sense_dcmv_R1 --output plot/graph/universal/sense_dcmv_R1.png --node_comparison_colors "#cf191b" "#ffe669" --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/sense_dcmv_R2 --output plot/graph/universal/sense_dcmv_R2.png --node_comparison_colors "#cf191b" "#ffe669" --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
```

### Demonstration: Histograms

The [`histogram.py`](#histogrampy) script allows plotting an alternative representation of the variations extracted from the DSB-window. Instead of displaying the distribution of *sequences* in a graph, the distribution of *variations* (insertion, deletion, or substitution) are shown in a histogram. 

Note that the histogram script can only be run on *individual* data (output from `preprocess.py`) rather than *comparison* data (output from `comparison.py`). Moreover, it always uses the data *with substitutions* (see [here](#substitutions-in-bowtie-2)), so that the substitution distribution can be visualized.

The following are examples of using the `histogram.py` command:

```
python histogram.py --input data_output/db_R1 --output plot/histogram/db_R1_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_output/db_R1 --output plot/histogram/db_R1_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_output/db_R1 --output plot/histogram/db_R1_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_output/db_R2 --output plot/histogram/db_R2_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_output/db_R2 --output plot/histogram/db_R2_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_output/db_R2 --output plot/histogram/db_R2_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_output/sense_R1 --output plot/histogram/sense_R1_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_output/sense_R1 --output plot/histogram/sense_R1_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_output/sense_R1 --output plot/histogram/sense_R1_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_output/sense_R2 --output plot/histogram/sense_R2_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_output/sense_R2 --output plot/histogram/sense_R2_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_output/sense_R2 --output plot/histogram/sense_R2_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_output/dcmv_R1 --output plot/histogram/dcmv_R1_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_output/dcmv_R1 --output plot/histogram/dcmv_R1_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_output/dcmv_R1 --output plot/histogram/dcmv_R1_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
python histogram.py --input data_output/dcmv_R2 --output plot/histogram/dcmv_R2_substitution.png --color "#bfbfbf" --variation_type substitution --label_type relative
python histogram.py --input data_output/dcmv_R2 --output plot/histogram/dcmv_R2_insertion.png --color "#ffa500" --variation_type insertion --label_type relative
python histogram.py --input data_output/dcmv_R2 --output plot/histogram/dcmv_R2_deletion.png --color "#8080ff" --variation_type deletion --label_type relative
```

<!-- SHOULD THE GRAPH DEFINITION BE DESCRIBED EARLIER? -->
<!-- Define the DSB-window term? -->