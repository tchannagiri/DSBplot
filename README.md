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

## Commands

This package is based on the following four commands:

* `preprocess.py`: takes as input the trimmed FASTQ reads files and a FASTA file containing the reference sequence, and creates the intermediate tables needed for plotting the graphs and histograms. The input directories should represent replicate experiments (e.g., biological replicates).
* `comparison.py`: takes as input two directories created by `preprocess.py` and creates a directory that contains the intermediate tables needed for plotting comparison graphs. Both of the directories must have been created with the same reference sequence.
* `graph.py`: takes as input a collection of the output directories of either `preprocess.py` or `comparison.py`, lays out alignment-sequences in all inputs, and plots a separate graph for each input.
* `histogram.py`: takes as input an output directories of `preprocess.py` (outputs of `comparison.py` are not valid) and plots a histogram showing the type and position of variations (insertions, deletion, or substitutions) in the alignment-sequences.

More information about each command is given in the following subsections.

### `preprocess.py`

#### Input

It is expected that the input FASTQ files to `preprocess.py` are *trimmed*, meaning that the adaptors have been removed. It is also expected that the region of DNA between these adaptors is exactly the region of DNA represented by the reference sequence. If a given read ihas been perfectly repaired, it should identical with the reference sequence (assuming no substitution errors due to library preparation or sequencing).

#### Substitutions in Bowtie 2

The preprocessing pipeline produces two different versions of most files: one *ignoring* substitutions (suffix "withoutSubst") and another *keeping* substitutions (suffix "withSubst"). The processing for files that ignore substitutions contains an extra step that takes the alginments output from Bowtie 2 and removes all substitutions (or mismatches) by converting them to matches. The reason for ignoring substitutions was we observed a similar patterns of subtitutions in both the experiment group (where DSBs had been induced in DNA) and control group (where no DSBs had been induced), indicating the subtitutions may largely be due to DNA damage in library preparation or sequencing errors rather than the DSB repair process. In the command `graph.py`, the `--subst_type` parameter controls whether to use the output with or without substitutions. The `histogram.py` command only uses the output with substitutions, since it is used to diagnose the similarity of the substitution profile between control and experiment.

#### Prepcocessing stages

1. Align FASTQ reads against FASTA reference sequence. This stage has multiple inpute FASTQ files representing independent repeats (or biological/tehchnical replicates) of the experimental condition. The alignment is done independently for each of the input FASTQs. The output are [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) files. This stages requires the [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) commands `bowtie2-build-s` and `bowtie2-align-s` to be available on the PATH. The output directories are:
    * `0_bowtie2_build`: The Bowtie 2 index files built with `bowtie2-build-s`.
    * `0_sam`: The SAM file output from the alignment by `bowtie2-align-s`.
2. Filter the resulting SAM files for alignments that represent NHEJ repair using various heuristics. This is done independently for each SAM files from stage 1. The stages of the filtering are:
    1. Discard "alignments" that represent a failed aligment in the SAM format.
    2. Discard alignments where the read sequence has length less than `MIN_LENGTH` (a parameter).
    3. Discard alignment where the left-most (5'-most) position of the read does not align with the left-most position of the reference sequence.
    4. If the in/del positions of the alignment are not adjacent to or around the DSB position, try and shift them towards the DSB position in such a way that remainder of the alignment stays roughly the same (the number of substitutions does not increase). If such a modification of the alignment cannot be found, discard the alignment.
    5. Finally, if the resulting alignment has in/dels that are not at consecutive positions, discard the alignment.
   
   If multiple reads have exactly the same nucleotide sequences but had different alignments with the reference, they will be forced to have the same alignment as the first such read encountered. Note, that while we require the left-most (5'-most) position on the read to align with the left-most on the reference, the same is not true for the right-most (3'-most) positions. This is because it is possible for reads to be too small and only capture the 5' primer but not the 3' primer. However, we require the read to at least length `DSB_POS + 1` (see parameter `--min_length`). The output will be tables (TSV format) that contain the unique alignments and their read counts. The output will be in the files `1_filter_nhej/[N].tsv`, where `[N]` is replaced with an integer representing the repeat number, in the same order they are specified on the command line.
  
3. Combine the (potentially) multiple tables output from stage 2 into a single table (TSV format) with multiple frequency columns, one for each input. The output will be in `2_combine_repeat/out.tsv`.
4. For each unique alignment in the output table of stage 3, extract the alignment-window around the DSB position. That is, the part of the alignment that corresponds to the positions in `DSB_POS - WINDOW_SIZE + 1` to `DSB_POS + WINDOW_SIZE`, where `DSB_POS` and `WINDOW_SIZE` are parameters. If different alignments become identical after extracting their windows, their read counts will be summed. We also require that *anchors* to be present on either side of the extracted windows to make sure that the variations near the DSB do not spill outside the window. Anchors are the parts of the alignment that correspond to the `ANCHOR_SIZE` (a parameter) nucleotides on the left (5') and right (3') of the window region. We require that each anchor have at most `ANCHOR_MISMATCHES` mismatches and no in/dels. The output is a table (TSV format). The output will be in directory `3_window`, in files `window_count_withoutSubst.tsv` and `window_count_withSubst.tsv`. The "withoutSubst" version of the output has substitutions in the alignment removed by replacing the nucleotide on the read sequence with perfect matches with the reference sequence. See [here](#substitutions-in-bowtie-2) for the justification of ignoring substitutions.
5. Covert the raw read counts into frequencies in the range [0,1]. This is done by dividing the read counts by the total number of reads in the library. Since the total number of reads in each sample might not be exactly the sum of the read counts (e.g., if the trimming process removed some reads), the total reads are input as a parameters. In this stage, the separate frequency columns for each repeat are collapsed into a single column by taking their mean. Also, the alignments whose minimum frequency across the repeats is below a user-defined threshold will be discarded. The output will be in the subdirectory `3_window`, in files:
  * `window_freq_withSubst.tsv`/`window_freq_withoutSubst.tsv`: Frequency tables with separate columns for the separate repeats.
  * `window_freq_filter_withSubst.tsv`/`window_freq_filter_withoutSubst.tsv`: The previous frequency tables with the minimum frequency filter applied.
  * `window_freq_filter_mean_withSubst.tsv`/`window_freq_filter_mean_withoutSubst.tsv`: The previous frequency tables with the separate frequency columns for the repeats collapsed into a single columns by taking their mean.
6. Precompute data needed to plot the graphs, such as adjacency information and summary statistics of the graphs. Output in the subdirectory `4_graph`.
7. Precompute data needed to plot the histograms, such as the position and frequency of the different types of variations. Output in the subirectory `5_histogram`.

#### Parameters

FIXME ADDED NEW PARAMETER
Output of `python preprocess.py --help`:
```
usage: preprocess.py [-h] --input INPUT [INPUT ...] --output OUTPUT --ref_seq_file REF_SEQ_FILE --dsb_pos DSB_POS
                     [--window_size WINDOW_SIZE] [--anchor_size ANCHOR_SIZE] [--anchor_mismatches ANCHOR_MISMATCHES] --total_reads    
                     TOTAL_READS [TOTAL_READS ...] [--freq_min FREQ_MIN] --label LABEL [--quiet]

Perform alignment and preprocessing for raw FASTQ data.

options:
  -h, --help            show this help message and exit
  --input INPUT [INPUT ...]
                        Input FASTQ files of raw reads. Each file is considered a repeat of the same experiment.
  --output OUTPUT       Output directory.
  --ref_seq_file REF_SEQ_FILE
                        Size of anchor on left/right of the window to check for mismatches.
  --anchor_mismatches ANCHOR_MISMATCHES
                        Maximum number of mismatches allowed on the left/right anchor sequences. Reads with more than the allowed     
                        number of mismatches on the left/right anchor will be discarded. This limit is applied to the left/right      
                        anchors separately.
  --total_reads TOTAL_READS [TOTAL_READS ...]
                        Total reads for each file. Must be the same number of arguments as the number of Count columns in INPUT.      
  --freq_min FREQ_MIN   Minimum frequency for output in windows_freq_filter_mean. Sequences with frequences <= this are discarded.    
  --label LABEL         Label of the experiment to be used in plot titles.
  --quiet               If present, do no output verbose log message.
```

Please reference to CITATION for more details about the preprocessing.

### `comparison.py`

This stage takes two directories output from the `preprocess.py` stages and creates a third directory that contains analogous data for creating *comparison* graphs of the two samples. The reference sequences (after restricting to the specified window) should be identical between the two library for meaningful comparisons. The output data will contain the same columns as that output from `preprocess.py` but will have two additional columns (with suffixes `_1` and `_2`) for the frequencies of the two different samples (the original frequency column, with no suffix, will contain the maximum of these two). The output will be in the subdirectores `3_window`, `4_graph`, and `5_histogram`.

Output of `python comparson.py --help`:
```
usage: comparison.py [-h] --input INPUT INPUT --output OUTPUT

Combine two individual experiment directories to make a comparison experiment directory for comparison graphs. The experiments must
be have the same windowed reference sequence.

options:
  -h, --help           show this help message and exit
  --input INPUT INPUT  Input directories of data created with "preprocess.py".
  --output OUTPUT      Output directory.
```

### `graph.py`

This command performs the vertex layout and final plotting of the variation-distance graphs. Multiple inputs may be specified to layout the vertices jointly, as long as their windowed reference sequences are identical. See the [Graphs](#graphs) section for a description of the meaning of different visual properties of the graph (e.g., edges, vertex sizes, colors, etc.).

Output of `python graph.py --help`:
```
usage: graph.py [-h] --input INPUT [INPUT ...] [--output OUTPUT [OUTPUT ...]] [--title TITLE [TITLE ...]]
                [--layout {mds_layout,radial_layout,universal_layout,fractal_layout,kamada_layout,spectral_layout,spring_layout,shell_layout,spiral_layout,circular_layout,multipartite_layout}]
                [--universal_layout_y_axis_x_pos UNIVERSAL_LAYOUT_Y_AXIS_X_POS]
                [--universal_layout_x_axis_deletion_y_pos UNIVERSAL_LAYOUT_X_AXIS_DELETION_Y_POS]
                [--universal_layout_x_axis_deletion_label_type {relative,absolute}]
                [--universal_layout_x_axis_insertion_y_pos UNIVERSAL_LAYOUT_X_AXIS_INSERTION_Y_POS]
                [--universal_layout_y_axis_y_range UNIVERSAL_LAYOUT_Y_AXIS_Y_RANGE UNIVERSAL_LAYOUT_Y_AXIS_Y_RANGE]
                [--universal_layout_x_axis_x_range UNIVERSAL_LAYOUT_X_AXIS_X_RANGE UNIVERSAL_LAYOUT_X_AXIS_X_RANGE]
                [--universal_layout_y_axis_deletion_max_tick UNIVERSAL_LAYOUT_Y_AXIS_DELETION_MAX_TICK]
                [--universal_layout_y_axis_insertion_max_tick UNIVERSAL_LAYOUT_Y_AXIS_INSERTION_MAX_TICK]
                [--universal_layout_x_scale_insertion UNIVERSAL_LAYOUT_X_SCALE_INSERTION]
                [--universal_layout_y_scale_insertion UNIVERSAL_LAYOUT_Y_SCALE_INSERTION]
                [--universal_layout_x_scale_deletion UNIVERSAL_LAYOUT_X_SCALE_DELETION]
                [--universal_layout_y_scale_deletion UNIVERSAL_LAYOUT_Y_SCALE_DELETION] [--subst_type {withSubst,withoutSubst}]       
                [--node_freq_range NODE_FREQ_RANGE NODE_FREQ_RANGE] [--node_px_range NODE_PX_RANGE NODE_PX_RANGE]
                [--node_outline_scale NODE_OUTLINE_SCALE] [--node_comparison_colors NODE_COMPARISON_COLORS NODE_COMPARISON_COLORS]    
                [--node_comparison_color_type {continuous,discrete}]
                [--node_freq_ratio_range NODE_FREQ_RATIO_RANGE NODE_FREQ_RATIO_RANGE]
                [--node_reference_outline_color NODE_REFERENCE_OUTLINE_COLOR] [--node_outline_color NODE_OUTLINE_COLOR]
                [--node_fill_color NODE_FILL_COLOR]
                [--variation_types {insertion,deletion,substitution,mixed,none} [{insertion,deletion,substitution,mixed,none} ...]]   
                [--variation_type_colors VARIATION_TYPE_COLORS VARIATION_TYPE_COLORS VARIATION_TYPE_COLORS VARIATION_TYPE_COLORS VARIATION_TYPE_COLORS]
                [--edge_show {True,False}] [--edge_types {indel,substitution} [{indel,substitution} ...]] [--edge_scale EDGE_SCALE]   
                [--stats] [--width_px WIDTH_PX] [--height_px HEIGHT_PX] [--margin_top_px MARGIN_TOP_PX]
                [--margin_bottom_px MARGIN_BOTTOM_PX] [--margin_left_px MARGIN_LEFT_PX] [--margin_right_px MARGIN_RIGHT_PX]
                [--line_width_scale LINE_WIDTH_SCALE] [--font_size_scale FONT_SIZE_SCALE] [--crop_x CROP_X CROP_X]
                [--crop_y CROP_Y CROP_Y] [--range_x RANGE_X RANGE_X] [--range_y RANGE_Y RANGE_Y] [--legend]
                [--legend_colorbar_scale LEGEND_COLORBAR_SCALE] [--legend_spacing_px LEGEND_SPACING_PX] [--separate_components]       
                [--interactive] [--reverse_complement {0,1} [{0,1} ...]]

Lay out and plot variation-distance graphs. Uses the output from "get_graph_data" as input. For more information about the layouts    
please see the publication FIXME.

options:
  -h, --help            show this help message and exit
  --input INPUT [INPUT ...]
                        List of directories with the data files produced with "get_graph_data.py". All libraries specified here must  
                        have the same windowed reference sequence (i.e., the 20bp section of the reference around the DSB site        
                        should be the same in all libraries). All the libraries will be laid out using common x/y-coordinate
                        assignments to the vertices. (default: None)
  --output OUTPUT [OUTPUT ...]
                        Output file. If not given no output will be written (useful only when using "--interactive"). The file        
                        extension should be either ".png" or ".html" for a static PNG or interactive HTML output respectively.        
                        Should be either 0 arguments or the number of arguments should match the number of input directories.
                        (default: None)
  --title TITLE [TITLE ...]
                        If present, adds a title to the plot with this value. Number of arguments should match the number of input    
                        files. (default: None)
  --layout {mds_layout,radial_layout,universal_layout,fractal_layout,kamada_layout,spectral_layout,spring_layout,shell_layout,spiral_layout,circular_layout,multipartite_layout}
                        The algorithm to use for laying out the graph. (default: universal_layout)
  --universal_layout_y_axis_x_pos UNIVERSAL_LAYOUT_Y_AXIS_X_POS
                        If present, shows a y-axis at the given x position showing the distances to the reference. Univeral layout    
                        only. (default: None)
  --universal_layout_x_axis_deletion_y_pos UNIVERSAL_LAYOUT_X_AXIS_DELETION_Y_POS
                        If present, shows an x-axis for deletions at the given y position showing the approximate position of the     
                        deleted ranges. Univeral layout only. (default: None)
  --universal_layout_x_axis_deletion_label_type {relative,absolute}
                        The type of labeling to use for the universal layout deletion x-axis (if present). "relative" labels have 0   
                        in the middle with negative/positive values on the left/right. "absolute" labels have 1 on the left and the   
                        length of the reference sequence on the right. (default: relative)
  --universal_layout_x_axis_insertion_y_pos UNIVERSAL_LAYOUT_X_AXIS_INSERTION_Y_POS
                        If present, shows a x-axis for insertions at the given y position showing the first nucleotide of inserted    
                        sequences. Univeral layout only. (default: None)
  --universal_layout_y_axis_y_range UNIVERSAL_LAYOUT_Y_AXIS_Y_RANGE UNIVERSAL_LAYOUT_Y_AXIS_Y_RANGE
                        If showing an y-axis for the universal layout, the min and max y-position of the line. (default: [nan, nan])  
  --universal_layout_x_axis_x_range UNIVERSAL_LAYOUT_X_AXIS_X_RANGE UNIVERSAL_LAYOUT_X_AXIS_X_RANGE
                        If showing an x-axis for the universal layout, the min and max x-position of the line. (default: [nan, nan])  
  --universal_layout_y_axis_deletion_max_tick UNIVERSAL_LAYOUT_Y_AXIS_DELETION_MAX_TICK
                        If showing an y-axis for the universal layout, the max tick value for the deletion side. (default: None)      
  --universal_layout_y_axis_insertion_max_tick UNIVERSAL_LAYOUT_Y_AXIS_INSERTION_MAX_TICK
                        If showing an y-axis for the universal layout, the max tick value for the insertion side. (default: None)     
  --universal_layout_x_scale_insertion UNIVERSAL_LAYOUT_X_SCALE_INSERTION
                        The factor for determining the scale on the universal layout insertion x-axis. Insertion vertex
                        x-coordinates will be multiplied by this value. (default: 10)
  --universal_layout_y_scale_insertion UNIVERSAL_LAYOUT_Y_SCALE_INSERTION
                        The factor for determining the scale on the universal layout insertion y-axis. Insertion vertex
                        y-coordinates will be multiplied by this value. (default: 3)
  --universal_layout_x_scale_deletion UNIVERSAL_LAYOUT_X_SCALE_DELETION
                        The factor for determining the scale on the universal layout deletion x-axis. Deletion vertex x-coordinates   
                        will be multiplied by this value. (default: 2)
  --universal_layout_y_scale_deletion UNIVERSAL_LAYOUT_Y_SCALE_DELETION
                        The factor for determining the scale on the universal layout deletion y-axis. Deletion vertex y-coordinates   
                        will be multiplied by this value. (default: 1)
  --subst_type {withSubst,withoutSubst}
                        Whether to plot data with or without substitutions. (default: withoutSubst)
  --node_freq_range NODE_FREQ_RANGE NODE_FREQ_RANGE
                        Min and max frequency to determine node size. Higher frequencies are clipped to this value. (default:
                        [1e-05, 1])
  --node_px_range NODE_PX_RANGE NODE_PX_RANGE
                        Largest node size as determined by the frequency. (default: [10, 120])
  --node_outline_scale NODE_OUTLINE_SCALE
                        How much to scale the node outline width (thickness). Values > 1 increase the width; values < 1 decrease the  
                        width. (default: 4)
  --node_comparison_colors NODE_COMPARISON_COLORS NODE_COMPARISON_COLORS
                        The colors to use in the gradient when the node colors show the frequency ratio of two experiments. May be    
                        specified in hex (e.g., "#ff0000" for red) or with recognized keywords such as "red", "blue", "green".        
                        (default: ['#ff0000', '#0000ff'])
  --node_comparison_color_type {continuous,discrete}
                        The type of color scheme to use for coloring nodes in a comparison graph. The "continuous" scheme uses a      
                        gradient of colors from the min to max ratio. The "discrete" schemes uses three colors to indicate that the   
                        ratio is < the min ratio, between the min and max ratio, or > the max ratio. The min and max ratios are       
                        determined by NODE_FREQ_RATIO_RANGE and the corresponding colors are determined by NODE_COMPARISON_COLORS.    
                        (default: continuous)
  --node_freq_ratio_range NODE_FREQ_RATIO_RANGE NODE_FREQ_RATIO_RANGE
                        The two frequencies used to determine node colors for comparison graphs. Also controls the range of ratios    
                        displayed on the frequency-ratio colorbar legend. Typically, the min value should be < 1 and the max value    
                        should be > 1. (default: [0.6666666666666666, 1.5])
  --node_reference_outline_color NODE_REFERENCE_OUTLINE_COLOR
                        Color to make the reference node outline. May be specified in hex (e.g., "#ff0000" for red) or with
                        recognized keywords such as "red", "blue", "green". (default: #32cd32)
  --node_outline_color NODE_OUTLINE_COLOR
                        Color to make the default node outline. May be specified in hex (e.g., "#ff0000" for red) or with recognized  
                        keywords such as "red", "blue", "green". (default: #000000)
  --node_fill_color NODE_FILL_COLOR
                        Color to make the default node fill. May be specified in hex (e.g., "#ff0000" for red) or with recognized     
                        keywords such as "red", "blue", "green". (default: #ffffff)
  --variation_types {insertion,deletion,substitution,mixed,none} [{insertion,deletion,substitution,mixed,none} ...]
                        The variation types that should be included in the graph. This should be a list of the types: "insertion",    
                        "deletion", "substitution", "mixed", "none". Default value: "insertion", "deletion", "none". "mixed" means    
                        nodes that have multiples variation types (e.g. insertions and substitutions). "none" means the reference     
                        node (no variations). May be specified in hex (e.g., "#ff0000" for red) or with recognized keywords such as   
                        "red", "blue", "green". (default: ['insertion', 'deletion', 'substitution', 'mixed', 'none'])
  --variation_type_colors VARIATION_TYPE_COLORS VARIATION_TYPE_COLORS VARIATION_TYPE_COLORS VARIATION_TYPE_COLORS VARIATION_TYPE_COLORS
                        The colors for the different variations types. They must be specified in the order: INSERTION DELETION        
                        SUBSTITUTION MIXED NONE. MIXED is the color for nodes with multiple types of variations (e.g. insertions and  
                        substitutions); NONE is the color for the reference node (no variations). May be specified in hex (e.g.,      
                        "#ff0000" for red) or with recognized keywords such as "red", "blue", "green". (default: ['#ffa500',
                        '#8080ff', '#808080', '#00ff00', '#ffffff'])
  --edge_show {True,False}
                        Whether to show edges between nodes. (default: True)
  --edge_types {indel,substitution} [{indel,substitution} ...]
                        The edge types to show. (default: ['indel'])
  --edge_scale EDGE_SCALE
                        How much to scale the edges width (thickness). Values > 1 increase the width; values < 1 decrease the width.  
                        (default: 8)
  --stats               If present, show graph summary statistics in the left margin. (default: False)
  --width_px WIDTH_PX   The width of the plot in pixels. (default: 2400)
  --height_px HEIGHT_PX
                        The height of the plot in pixels. (default: 2400)
  --margin_top_px MARGIN_TOP_PX
                        The size of the top margin in pixels. (default: 300)
  --margin_bottom_px MARGIN_BOTTOM_PX
                        The size of the bottom margin in pixels. (default: 300)
  --margin_left_px MARGIN_LEFT_PX
                        The size of the left margin in pixels. (default: 300)
  --margin_right_px MARGIN_RIGHT_PX
                        The size of the right margin in pixels. (default: 300)
  --line_width_scale LINE_WIDTH_SCALE
                        How much to scale the line widths (aka thickness). Values > 1 increase the width; values < 1 decrease the     
                        width. (default: 8)
  --font_size_scale FONT_SIZE_SCALE
                        How much to scale the font size. Values > 1 increase the font size; values < 1 decrease it. (default: 2)      
  --crop_x CROP_X CROP_X
                        Range of the horizontal dimension to crop. Specified with normalized coords in range [0, 1]. (default: [0,    
                        1])
  --crop_y CROP_Y CROP_Y
                        Range of the vertical dimension to crop. Specified in normalized coords in range [0, 1]. (default: [0, 1])    
  --range_x RANGE_X RANGE_X
                        Range of x-axis for plotting.If not specified chosen automatically to either show all nodes or a preset       
                        value for the layout. (default: [nan, nan])
  --range_y RANGE_Y RANGE_Y
                        Range of y-axis for plotting.If not specified chosen automatically to either show all nodes or a preset       
                        value for the layout. (default: [nan, nan])
  --legend              Whether to show a legend on the figure. (default: False)
  --legend_colorbar_scale LEGEND_COLORBAR_SCALE
                        How much to scale the colorbar legend (for frequency-ratio coloring). (default: 4)
  --legend_spacing_px LEGEND_SPACING_PX
                        Amount of vertical space in pixels between different legends. (default: 200)
  --separate_components
                        If present, separate the connected components of the graph. (default: False)
  --interactive         If present opens the interactive version in a browser. Uses the Ploty library figure.show() function to do    
                        so. (default: False)
  --reverse_complement {0,1} [{0,1} ...]
                        Whether to reverse complement the sequences in the data sets. If present, the number of values must be the    
                        same as the number of input directories. "1" mean reverse complement the sequence and "0" means do not. Used  
                        for making a layout for data sets that have reference sequences that are the reverse complements of each      
                        other. If "1" also uses the reverse complement of sequences when determining the display labels and hover     
                        text. This affects the universal layout and fractal layout. (default: None)
```

### `histogram.py`

This commands plots the 3D histograms, which summarize the distribution of the variations in the windows around the DSB site. Each histogram shows a single type of variation (insertion, deletion, or substitution), which is determined by the `--variation_type` parameter. The axes are described as:

* The $x$-axis indicates the position of the variations relative to the DSB position (or, alternatively, relative to the 5'-end of the window; see parameter `--label_type`). Deletions and substitutions always have an unambiguous position on the reference sequence. However, insertions are always between the two consecutive positions around the DSB on the reference sequence, and are thus assigned the left (5') position.
* The $y$-axis indicates the total number of in/dels on the sequence that the variation originated from.
* The $z$-axis indicates the total frequency of the variations with the given $x$- and $y$-coordinates.

For example, assume that an alignment had one insertion at position 5,  two deletions at positions 6 and 7, and overall frequency 0.1. Assume that the DSB is between positions 5 and 6, and we use relative labelling on the $x$-axis. Then, in the insertion histogram, the alignment would contribute +0.1 to the $z$-value of the bar at $x$-$y$-coordinate (-1, 3). In the deletion histogram the alignment would contribute +0.1 to the $z$-value of the bars at $x$-$y$-coordinates (1, 3) and (2, 3). If we used absolute labelling, then the respective $x$-$y$-coordinates would be (5, 3), (6, 3), and (7, 3).

Output of `histogram.py --help`:

```
usage: histogram.py [-h] --input INPUT --output OUTPUT --variation_type {substitution,insertion,deletion} [--reverse_pos]
                    --label_type {relative,absolute} [--color COLOR] [--freq_range FREQ_RANGE FREQ_RANGE]
                    [--freq_scale {linear,log}]

Plot 3d histograms showing variation type/position/frequency.

options:
  -h, --help            show this help message and exit
  --input INPUT         Directory with the data files. (default: None)
  --output OUTPUT       Output image file. Must have extension ".png". (default: None)
  --variation_type {substitution,insertion,deletion}
                        Which variation type to show in the histogram. (default: None)
  --reverse_pos         Whether to reverse the x-axis positions. Useful if comparing reverse strand data with forward strand data.    
                        (default: False)
  --label_type {relative,absolute}
                        Whether to index the x-axis by "absolute" positions on the reference sequence from 1 to [reference length],   
                        or "relative" positions from -[reference length]/2 to [reference length]/2 (skipping 0). [reference length]   
                        refers to the length of the reference sequence after extracting the window around the DSB site. (default:     
                        None)
  --color COLOR         Color of the bar graph. If not specified, a default color based on VARIATION_TYPE will be chosen. (default:   
                        None)
  --freq_range FREQ_RANGE FREQ_RANGE
                        Range of the z-axis frequency values to show. (default: [1e-05, 1])
  --freq_scale {linear,log}
                        Whether to use a linear or log scale for the z-axis frequency values. (default: log)
```

## Graphs

Here we describe the visual properties of the output graphs using the `graph.py` command. Note, when describing the graphs, the terms "sequence" and "vertex" are used interchangeably since there is a one-to-one correspondence between the windowed sequences output from [`preprocess.py`](#preprocesspy) and the vertices of the graph. For example, the phrase "the two vertices differ by a single insertion" actually means "the two sequences corresponding to the two vertices differ by a single insertion".

### Layouts

Different *layouts*, specified with the XX parameter, are used to position the vertices. We decribe a subset of them in the following.

* Kamada-Kawaii: a physical simulation algorithm for laying out graphs (Kamada and Kawai, 1989). The graph is modeled as a spring system, where each pair of vertices is connected by a spring of resting length proportional to the shortest-path distance between the two vertices. The algorithms seeks to minimize the energy of the system by finding the optimal positions of the vertices. The result is to bring the vertices to positions such that the Euclidean distance between pairs of vertices is close to the shortest-path distance between them. This is intended to layout the vertices in an aesthetically pleasing way that also reveals graph structure, such as clusters of vertices that are close together. The algorithm is implemented by the NetworkX Python package (Hagberg, 2008). As of this writing, NetworkX places the vertices in a circle for the initialization. Note, the shortest-path distance between two vertices is at least the Levenshtein distance between the two vertices, though equality is not necessary. The final coordinates of the vertices are always in [0, 1] (relevant if setting the `plot_range_x` and `plot_range_y` parameters).
* Universal layout: We call this layout the *universal* layout because every possible insertion and deletion is assigned a predefined position, independent of the subset of in/dels occurring in the dataset. This contrasts with the Kamada-Kawaii algorithm, which has a randomized initialization and positions vertices in a data-dependent manner. While this layout algorithm is less flexible and does not always position vertices based on their Levenshtein distances to one another, it has the advantage that it is reproducible and easily interpretable. The universal layout positions the reference vertex in the center of the figure, the insertion vertices above the reference, and the deletion vertices below the reference. For both the insertion and deletion vertices, the $y$-coordinate is determined by the vertices' Levenshtein distance from the reference, which is simply the number of nucleotides inserted or deleted in the vertices' sequence. If $d$ is a vertex's Levenshetin distance from the reference, the deletion vertices are placed $d$ units below the reference and insertion vertices are place $d$ units above. To determine the $x$-coordinates of vertices, we use different rules for insertions and deletions. Insertion vertices are ordered from left to right by alphabetical order, while deletion vertices are ordered from left to right based on the positions of the deleted vertices with respect to the DSB site. Figure XX gives a key to the vertex placement for the universal layout. In this layout, the vertex coordinates usually range between $\pm 2 * window\_size$, where `window_size` is the parameter from [preprocess.py](#preprocesspy) (relevant if setting the `plot_range_x` and `plot_range_y` parameters).
* Radial: This layout arranges vertices in concentric circles around the reference sequence. The reference vertex is placed at the center of the figure; insertion vertices are placed above the reference; and deletion vertices are placed below the the reference. A vertex's physical distance to the reference is proportional to its Levenshtein distance to the reference. That is, the first concentric circle around the reference is at Levenshtein distance 1, the next is at Levenshtein distance 2, and so on. In each concentric circle, the insertion vertices are arranged clockwise from highest to lowest frequency, and deletion vertices are arranged counterclockwise from highest to lowest frequency. This means that the vertices for both insertions and deletion are placed left to right from highest to lowest frequency. Several heuristics are used to perturb the vertices slightly so that they do not overlap and edges are not collinear. In this layout, the vertex coordinates (with default scaling) usually range between $\pm 4 * window\_size$, where `window_size` is the parameter from [preprocess.py](#preprocesspy) (relevant if setting the `plot_range_x` and `plot_range_y` parameters). However, the scaling on the insertion and deletion, $x$ and $y$-coordinates can be controlled by the FIXME `--universal_layout_x_scale_insertion`, FIXME parameters.
* Fractal. This layout is intended to display the insertion vertices only. The vertices with a single inserted nucleotide (A, C, G, or T) are placed in the center of the corresponding quadrants shown in Figure XX: A in the top left, C in the top right, G in the bottom left, and T in the bottom right. Vertices with two inserted nucleotides, and placed in the center of the corresponding square shown in Figure XX. Note the pattern, the squares for two nucleotides insertions starting with A are obtained by subdividing the A square into four quadrants and arranging AA, AC, AG, AT based on their second nucleotide in the same order as the original A, C, G, T squares: AA in the top left, AC in the top right, AG in the bottom left, and AT in the bottom right. This pattern can be continued for arbitrarily large kmers. The square for kmer $s_1 s_2 \cdots s_k$ will be obtained by subdividing the the square for $s_1 s_2 \cdots s_{k-1}$ into four quadrants and choosing the quadrant in the same manner as before, based on $s_k$. Note, this is layout is very similar to the insertion of the universal layout, except that the both dimensions are used to order the inserted nucleotides alphabetically instead of linearly from left to right. The final coordinates of the vertices are always in [-1, 1] (relevant if setting the `plot_range_x` and `plot_range_y` parameters).

### Graph aesthetics

* Vertices:
  * Size: For an individual graph, the vertex radius is a function of the log-frequency. This is controlled by the parameters FIXME LIST THEM. Vertices with frequencies $\leq$ `MIN_FREQ` get `MIN_RADIUS` and frequencies $\geq$ `MAX_FREQ` get `MAX_RADIUS`. Frequencies in between get a radius that varies linearly with the log-frequency.
  * Outline: the outline color is intended to differentiate the reference sequence vertex from the other vertices. The colors are determined by the XXX and XXX parameters.
  * Color: There color depends on whether the graph being plotted is an individual experiment or a compairison of experiments:
    * Variation type (`variation_type`): Only for individual graphs. Gives vertices a color depending on the type of variations present in their alignment: insertions only, deletions only, substitutions only, mixed (a combinations of insertions, deletions, and substitution), and none (the reference sequence). These colors may be specified in the `XXX` variable.
    * Frequency ratio (continuous): Only for comparison graphs. Indicate "continuous" in the XX parameter. Uses a color gradient that indicates the ratio of frequencies in the two experiments, for each vertex. The ratio is computed by putting the first sample in the numerator and second sample in the denominator (this order is determined in the [preprocessing](#preprocesspy) stage). The colors at the extreme ends of the gradient are determined by the XX parameter, and values in the middle are smoothly interpolated. The ratios at the extreme ends of the gradient are determined by the XX parameter.
    * Frequency ratio (discrete): Only for comparison graphs. Indicate "discrete" in the XX parameter. Uses three colors to indicate the ratio of frequencies in the two experiments, for each vertex. If the frequency ratio is less than the first value in XX, the second color in XX is displayed (meaning higher in sample 2); If the frequency ratio is greater than the second value in XX, the first color in XX is displayed (meaning higher in sample 1); if the frequency ratio is between the values in XX (inclusive), the vertex is colored white.
* Edges: All edges indicate a 1-nucleotide variation (insertion, deletion, or substitution) between the two vertices. Whether to show and which types of edges to show can be controlled by the FIXME parameters.
* Legends: The legends describe the vertex size, vertex outline, vertex color, and the edges. These can be all drawn by using the `--legend` parameter. The legends will be drawn in the right margin of the figure. To ensure enough room, use the `--margin_*` parameters. The different legend are laid out vertically, so to control the spacing between then use the `--legend_vertical_spacing` parameter.
* Title: A title can be optionally added to the top margin of the figure using the `--title` parameter.
* Universal layout axes: To illustrate the vertex placement for the universal layout, axes can be drawn on the figure. Using the `FIXME` parameter, a horizontal axis can be drawn on the deletion side (below the reference sequence)to show the midpoint position of the range of deletions that lies at that x-coordinate. Using the `FIXME` parameter, a horizontal axis can be drawn on the insertion side (above the reference sequence) to show the alphabetical order of the nucleotides. Using the `FIXME` parameter, a vertical axis can be drawn to show the Levenshtein distance of vertices from the reference vertex. To control the $x$ or $y$-coordinate of the axes, use the `FIXME` parameters.
* Graph statistics: Use the `--stats` parameter to display summary statistics of the graph in the left margin of the figure. If the graph has more than one connected component, only the component with the reference sequence is used, since otherwise the pairwise distance between vertices is not well-defined. To ensure the margin has enough space, use the `FIXME` parameters.

## Tutorial

Several example input files are available in the `data_input` directory:

* `data_input/ref_seq`: reference sequence FASTA files, representing the perfect repaired sequence for different samples.
* `data_input/fastq`: high-throughput sequencing data in FASTQ format for different samples. Each FASTQ file corresponds to an independent experiment.

Note, the FASTQ samples have been *trimmed*, meaning that we only capture the portion of the read between the adaptors, and low-quality reads have been already filtered out.

### Tutorial: Preprocessing

The [`preprocess.py`](#preprocesspy) command must be run before producing any figures. This stage aligns the input FASTQ files against the reference sequences, extracts the part of the alignment around the DSB position, and precomputes data tables describing the unique alignments and their frequencies. The following are examples of using the command:

```
python preprocess.py --input data_input/fastq/sense1_R1.fq data_input/fastq/sense2_R1.fq data_input/fastq/sense3_R1.fq data_input/fastq/sense4_R1.fq --ref_seq_file data_input/ref_seq/2DSB_R1_sense.fa --dsb_pos 67 --output data_output/sense_R1 --label sense_R1 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_input/fastq/sense1_R2.fq data_input/fastq/sense2_R2.fq data_input/fastq/sense3_R2.fq data_input/fastq/sense4_R2.fq --ref_seq_file data_input/ref_seq/2DSB_R2_sense.fa --dsb_pos 46 --output data_output/sense_R2 --label sense_R2 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_input/fastq/db1_R1.fq data_input/fastq/db2_R1.fq data_input/fastq/db3_R1.fq data_input/fastq/db4_R1.fq --ref_seq_file data_input/ref_seq/2DSB_R1_branch.fa --dsb_pos 67 --output data_output/db_R1 --label db_R1 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_input/fastq/db1_R2.fq data_input/fastq/db2_R2.fq data_input/fastq/db3_R2.fq data_input/fastq/db4_R2.fq --ref_seq_file data_input/ref_seq/2DSB_R2_branch.fa --dsb_pos 46 --output data_output/db_R2 --label db_R2 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_input/fastq/dcmv1_R1.fq data_input/fastq/dcmv2_R1.fq data_input/fastq/dcmv3_R1.fq data_input/fastq/dcmv4_R1.fq --ref_seq_file data_input/ref_seq/2DSB_R1_cmv.fa --dsb_pos 67 --output data_output/dcmv_R1 --label dcmv_R1 --total_reads 3000 3000 3000 3000
python preprocess.py --input data_input/fastq/dcmv1_R2.fq data_input/fastq/dcmv2_R2.fq data_input/fastq/dcmv3_R2.fq data_input/fastq/dcmv4_R2.fq --ref_seq_file data_input/ref_seq/2DSB_R2_cmv.fa --dsb_pos 46 --output data_output/dcmv_R2 --label dcmv_R2 --total_reads 3000 3000 3000 3000
```

Although there are 24 files in `data_input/fastq`, we only run the preprocessing pipeline six times because the biological repeats of each experiment areprocessed together. For example, `sense1_R1.fq`, `sense2_R1.fq`, `sense3_R1.fq`, and `sense4_R1.fq` are biological repeats of the "Sense, 2-DSB" experiment (the `R1` or `R2` indicates whether the forward or reverse strand of the read was sequenced).

All the biological repeats must, of course, use the same reference sequences. For example the four "Sense, 2-DSB" experiments use the `2DSB_R1_sense.fa` reference sequence.

The output will be a collection of tables in TSV (tab-separated value) format. The most important files are located in the subdirectories `4_graph` and `5_histogram`, which are the data used for plotting the graphs and histograms.  Other intermediate preprocessing files are stored in other files.

For more information about the parameters, use `python preprocess.py --help` and see section [`preprocess.py`](#preprocesspy).

### Tutorial: Comparisons

The [comparison](#comparisonpy) pipline is used to combine two outputs of the [preprocessing](#preprocesspy) command for plotting [comparison graphs](#XX). The following are examples:

```
python comparison.py --input data_output/sense_R1 data_output/db_R1 --output data_output/sense_db_R1
python comparison.py --input data_output/sense_R2 data_output/db_R2 --output data_output/sense_db_R2
python comparison.py --input data_output/sense_R1 data_output/dcmv_R1 --output data_output/sense_dcmv_R1
python comparison.py --input data_output/sense_R2 data_output/dcmv_R2 --output data_output/sense_dcmv_R2
```

Note, both experiments must have identical reference sequences after extracting its `2 * WINDOW_SIZE` nucleotides around the DSB position. The whole reference sequence of the experiments (i.e., the files in the `--ref_seq` parameter to `preprocess.py`) need not be identical so long as this window is.

The output will have a directory structure mirroring the output of `preprocess.py`, except that only the subdirectories `3_window`, `4_graph`, and `5_histogram` will be present (since there is no need to cache the other preprocessing data). The output files contain a combination of the data from both the experiments for creating comparison graphs.

### Tutorial: Graphing

The [`graph.py`](#graphpy) script is used to layout the aligned reads output from the [`preprocess.py`](#preprocesspy) and [`comparison.py`] scripts. The input should be a list of directories created by `preprocess.py`/`comparison.py` and the output will be PNG or HTML files that visually represent the graph of sequences extracted from reads around the DSB position.

Multiple inputs may be specified as long as all the experiments have identical reference sequences within the DSB window (potentially, after reverse complementing; see the `--reverse_complement` parameter). When multiple inputs are specified all sequences in all the inputs (potentially, after reverse complementing) are combined into a single "supergraph", the graph is laid out to determine the coordinates of each vertex, and then each graph is laid out separately according to the coordinates determined by the supergraph. This allows the output figures for each experiment to have the same coordinates for the same sequence, which facilitates visual comparison betwen experiments when when using data-dependent layouts such as the Kamada-Kawaii.

The output may be either PNG or HTML. HTML output may be opened in a browser and allows interactive inspection of the vertices by hovering over them.

To plot individual graphs (graphs with a single experiment), the output directories of `preprocess.py` should be used. To plot comparison graphs (graphs comparing two experiments), the output directories of `comparison.py` should be used.

There are many parameters for customizing the aesthetics of the output graphs, which are covered in more detail in section [Graphs](#graphs).

The following are example commands of using `graph.py` with single inputs with PNG output:
```
python graph.py --input data_output/db_R1 --output plot/graph/universal/db_R1.png --legend --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5  --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/db_R2 --output plot/graph/universal/db_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/sense_R1 --output plot/graph/universal/sense_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/sense_R2 --output plot/graph/universal/sense_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/dcmv_R1 --output plot/graph/universal/dcmv_R1.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python graph.py --input data_output/dcmv_R2 --output plot/graph/universal/dcmv_R2.png --layout universal_layout --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

```

The following is an example of using multiple inputs with PNG output:
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

### Tutorial: Histograms

The [`histogram.py`](#histogrampy) script allows plotting an alternative representation of the variations extracted from the DSB-window. Instead of displaying the distribution of *sequences* in a graph, the distribution of *variations* (insertion, deletion, or substitution) are shown in a histogram. 

Note, the histogram script can only be run on *individual* data (output from `preprocess.py`) rather than *comparison* data (output from `comparison.py`). Moreover, it always uses the data *without substitution* (see [here](#substitutions-in-bowtie-2)), so that the substitution distribution can be visualized also.

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

<!-- The alignment done between the reads in the input FASTQ and the reference sequence expects the first base of each read to align with the first nucleotide of the reference sequence. Therefore, the reads must be trimmed and the reference sequence selected in such a way that their left-most (5'-most) base pairs align. In the experiments performed in the study by Jeon et al. (LINK), the primers were designed to be about 50-150 base pairs away from the induced DSB site. In principle, this would mean reads repaired by NHEJ would have variations near the DSB site (say, within $\pm$ nucleotides), and allow the remainder of the read to otherwise match the reference perfectly (not counting substitution errors due to sequencing or library preparation). CONTINUE HERE. -->

<!-- SHOULD THE GRAPH DEFINITION BE DESCRIBED EARLIER? -->
<!-- Define the DSB-window term? -->