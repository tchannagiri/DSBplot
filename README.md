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

## Usage

This package is based on the following four commands:

* `preprocess.py`: takes as input the trimmed FASTQ reads files and a FASTA file containing the reference sequence, and creates the intermediate tables needed for plotting the graphs and histograms. The input directories should represent replicate experiments (e.g., biological replicates).
* `comparison.py`: takes as input two directories created by `preprocess.py` and creates a directory that contains the intermediate tables needed for plotting comparison graphs. Both of the directories must have been created with the same reference sequence.
* `graph.py`: takes as input a collection of the output directories of either `preprocess.py` or `comparison.py`, lays out alignment-sequences in all inputs, and plots a separate graph for each input.
* `histogram.py`: takes as input an output directories of `preprocess.py` (outputs of `comparison.py` are not valid) and plots a histogram showing the type and position of variations (insertions, deletion, or substitutions) in the alignment-sequences.

More information about each command is given in the following subsections.

### `preprocess.py`

The full output of `python preprocess.py --help` is given below.

```
usage: preprocess.py [-h] --input INPUT [INPUT ...] --output OUTPUT --ref_seq_file REF_SEQ_FILE --dsb_pos DSB_POS [--window_size WINDOW_SIZE]                     [--anchor_size ANCHOR_SIZE] [--anchor_mismatches ANCHOR_MISMATCHES] --total_reads TOTAL_READS [TOTAL_READS ...]
                     [--freq_min FREQ_MIN] --label LABEL [--quiet]

Perform alignment and preprocessing for raw FASTQ data.

options:
  -h, --help            show this help message and exit
  --input INPUT [INPUT ...]
                        Input FASTQ files of raw reads. Each file is considered a repeat of the same experiment.
  --output OUTPUT       Output directory.
  --ref_seq_file REF_SEQ_FILE
                        FASTA file with a single nucleotide sequence.
  --dsb_pos DSB_POS     Position on reference sequence immediately left of DSB site. Ie. the DSB is between position DSB_POS and DSB_POS + 1.  
  --window_size WINDOW_SIZE
                        Size of window around DSB site to extract. The nucleotides at the positions {DSB_POS - WINDOW_SIZE + 1, ..., DSB_POS   
                        + WINDOW_SIZE} are extracted. The actual number of nucleotides extracted may vary depending on how many
                        insertions/deletion the alignment of the sequence has.
  --anchor_size ANCHOR_SIZE
                        Size of anchor on left/right of the window to check for mismatches.
  --anchor_mismatches ANCHOR_MISMATCHES
                        Maximum number of mismatches allowed on the left/right anchor sequences. Reads with more than the allowed number of    
                        mismatches on the left/right anchor will be discarded. This limit is applied to the left/right anchors separately.     
  --total_reads TOTAL_READS [TOTAL_READS ...]
                        Total reads for each file. Must be the same number of arguments as the number of Count columns in INPUT.
  --freq_min FREQ_MIN   Minimum frequency for output in windows_freq_filter_mean. Sequences with frequences <= this are discarded.
  --label LABEL         Label of the experiment to be used in plot titles.
  --quiet               If present, do no output verbose log message.
```

It is expected that the input FASTQ files to `preprocess.py` are *trimmed*, meaning that the adaptors have been removed. It is also expected that the region of DNA between these adaptors is exactly the region of DNA represented by the reference sequence. If a given read ihas been perfectly repaired, it should identical with the reference sequence (assuming no substitution errors due to library preparation or sequencing).

(MENTON THE WITH/WITHOUT SUBST OPTION HERE!)

The stages of the preprocessing are:
1. Align FASTQ reads against FASTA reference sequence. This stage has multiple inpute FASTQ files representing independent "repeats" of the experimental condition. The alignment is done independently for each of the input FASTQs. The output are [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) files.  MENTION EXACT FILE NAMES. MENTION THAT THE BOWTIE TABLES WILL BE BUILT ALSO.
2. Filter the resulting SAM files for alignments that represent NHEJ repair using various heuristics. This is done independently for each SAM files from stage 1. The output are tables (TSV format) that contain the unique alignments and their read counts. SHOULD THE FILTERING BET DESCRIBED IN MORE DETAIL HERE?  MENTION EXACT FILE NAMES.
3. Combine the (potentially) multiple tables output from stage 2 into a single table (TSV file) with multiple frequency columns, one for each input.
4. For each unique alignment in the output table of stage 3, extract the alignment-window around the DSB position. The length of this window and the DSB position is an input parmeter. If different alignents have the same alignment-window, their read counts will be summed. The output is a table (TSV format). MENTION EXACT FILE NAMES.
5. Covert the raw read counts into frequencies in the range [0,1]. This is done by dividin the read counts by the total number of reads in the library. Since the total number of reads in each sample might not be exactly the sum of the read counts (e.g., if the trimming process removed some reads), the total reads are input as a parameters. In this stage, the multiple separate frequency columns are collapsed into a single column by taking the mean. Also, the alignments that have a minimum frequecy across all the repeats is under a user-defined threshold will be discarded. The output is three tables:
  * `window_freq_withSubst.tsv`/`window_freq_withoutSubst.tsv`: frequency tables with separate columns for the separate repeats.
  * `window_freq_filter_withSubst.tsv`/`window_freq_filter_withoutSubst.tsv`: the previous frequency tables with the low frequency alignments removed.
  * `window_freq_filter_mean_withSubst.tsv`/`window_freq_filter_mean_withoutSubst.tsv`: the previous frequency tables with the multiples frequency columns for the repeats collapsed into a single columns by taking the mean.
The output will be in the `3_window_freq` directory.
6. Precompute data needed to plot the graphs, such as adjacency information and summary statistics of the graphs. Output in the subdirectory `4_graph`.
7. Precompute data needed to plot the histograms, such as the position and frequency of the different types of variations. Output in the subirectory `5_histogram`.

Please reference to CITATION for more details about the preprocessing.

### `comparison.py`

This stage takes two directories output from the `preprocess.py` stages and creates a third directory that contains analogous data for creating *comparison* graphs of the two samples. The reference sequences (after restricting to the specified window) should be identical between the two library for meaningful comparisons. The output data will contain the same columns as that output from `preprocess.py` but will have two additional columns (with suffixes `_1` and `_2`) for the frequencies of the two different samples (the original frequency column, with no suffix, will contain the maximum of these two). The output will be in the subdirectores `3_window`, `4_graph`, and `5_histogram`.

SHOW THE OUTPUT OF THE --HELP.

### `graph.py`

This stage performs the vertex layout and final plotting of the variation-distance graphs. Multiple inputs may be specified to layout the vertices jointly, as long as their windowed reference sequences are identical. See the [Figures](#figures) section for a complete description of the meaning of different visual properties of the graph (e.g., edges, vertex sizes, colors, etc.).

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
                [--node_reference_outline_color NODE_REFERENCE_OUTLINE_COLOR] [--node_outline_color NODE_OUTLINE_COLOR]
                [--node_fill_color NODE_FILL_COLOR]
                [--variation_types {insertion,deletion,substitution,mixed,none} [{insertion,deletion,substitution,mixed,none} ...]]
                [--variation_type_colors VARIATION_TYPE_COLORS VARIATION_TYPE_COLORS VARIATION_TYPE_COLORS VARIATION_TYPE_COLORS VARIATION_TYPE_COLORS]  
                [--edge_show EDGE_SHOW [EDGE_SHOW ...]] [--edge_types EDGE_TYPES [EDGE_TYPES ...]] [--edge_scale EDGE_SCALE] [--stats]
                [--width_px WIDTH_PX] [--height_px HEIGHT_PX] [--margin_top_px MARGIN_TOP_PX] [--margin_bottom_px MARGIN_BOTTOM_PX]
                [--margin_left_px MARGIN_LEFT_PX] [--margin_right_px MARGIN_RIGHT_PX] [--line_width_scale LINE_WIDTH_SCALE]
                [--font_size_scale FONT_SIZE_SCALE] [--crop_x CROP_X CROP_X] [--crop_y CROP_Y CROP_Y] [--range_x RANGE_X RANGE_X]
                [--range_y RANGE_Y RANGE_Y] [--legend] [--legend_colorbar_scale LEGEND_COLORBAR_SCALE] [--legend_spacing_px LEGEND_SPACING_PX]
                [--separate_components] [--interactive] [--reverse_complement {0,1} [{0,1} ...]]

Lay out and plot variation-distance graphs. Uses the output from "get_graph_data" as input. For more information about the layouts please see the        
publication FIXME.

options:
  -h, --help            show this help message and exit
  --input INPUT [INPUT ...]
                        List of directories with the data files produced with "get_graph_data.py". All libraries specified here must have the same       
                        windowed reference sequence (i.e., the 20bp section of the reference around the DSB site should be the same in all libraries).   
                        All the libraries will be laid out using common x/y-coordinate assignments to the vertices.
  --output OUTPUT [OUTPUT ...]
                        Output file. If not given no output will be written (useful only when using "--interactive"). The file extension should be       
                        either ".png" or ".html" for a static PNG or interactive HTML output respectively. Should be either 0 arguments or the number    
                        of arguments should match the number of input directories.
  --title TITLE [TITLE ...]
                        If present, adds a title to the plot with this value. Number of arguments should match the number of input files.
  --layout {mds_layout,radial_layout,universal_layout,fractal_layout,kamada_layout,spectral_layout,spring_layout,shell_layout,spiral_layout,circular_layout,multipartite_layout}
                        The algorithm to use for laying out the graph.
  --universal_layout_y_axis_x_pos UNIVERSAL_LAYOUT_Y_AXIS_X_POS
                        If present, shows a y-axis at the given x position showing the distances to the reference. Univeral layout only.
  --universal_layout_x_axis_deletion_y_pos UNIVERSAL_LAYOUT_X_AXIS_DELETION_Y_POS
                        If present, shows an x-axis for deletions at the given y position showing the approximate position of the deleted ranges.        
                        Univeral layout only.
  --universal_layout_x_axis_deletion_label_type {relative,absolute}
                        The type of labeling to use for the universal layout deletion x-axis (if present). "relative" labels have 0 in the middle with   
                        negative/positive values on the left/right. "absolute" labels have 1 on the left and the length of the reference sequence on     
                        the right.
  --universal_layout_x_axis_insertion_y_pos UNIVERSAL_LAYOUT_X_AXIS_INSERTION_Y_POS
                        If present, shows a x-axis for insertions at the given y position showing the first nucleotide of inserted sequences. Univeral   
                        layout only.
  --universal_layout_y_axis_y_range UNIVERSAL_LAYOUT_Y_AXIS_Y_RANGE UNIVERSAL_LAYOUT_Y_AXIS_Y_RANGE
                        If showing an y-axis for the universal layout, the min and max y-position of the line.
  --universal_layout_x_axis_x_range UNIVERSAL_LAYOUT_X_AXIS_X_RANGE UNIVERSAL_LAYOUT_X_AXIS_X_RANGE
                        If showing an x-axis for the universal layout, the min and max x-position of the line.
  --universal_layout_y_axis_deletion_max_tick UNIVERSAL_LAYOUT_Y_AXIS_DELETION_MAX_TICK
                        If showing an y-axis for the universal layout, the max tick value for the deletion side.
  --universal_layout_y_axis_insertion_max_tick UNIVERSAL_LAYOUT_Y_AXIS_INSERTION_MAX_TICK
                        If showing an y-axis for the universal layout, the max tick value for the insertion side.
  --universal_layout_x_scale_insertion UNIVERSAL_LAYOUT_X_SCALE_INSERTION
                        The factor for determining the scale on the universal layout insertion x-axis. X-axis values will be between +/- this value.     
  --universal_layout_y_scale_insertion UNIVERSAL_LAYOUT_Y_SCALE_INSERTION
                        The factor for determining the scale on the universal layout insertion y-axis. Each level of the y-axis (each level has
                        vertices with the same number of insertions) will be this large.
  --universal_layout_x_scale_deletion UNIVERSAL_LAYOUT_X_SCALE_DELETION
                        The factor for determining the scale on the universal layout deletion x-axis. Shifting a deletion left/right by 1 nucleotide     
                        will shift the corresponding vertex by this much.
  --universal_layout_y_scale_deletion UNIVERSAL_LAYOUT_Y_SCALE_DELETION
                        The factor for determining the scale on the universal layout deletion y-axis. Each level of the y-axis (each level has vertices  
                        with the same number of deletions) will be this large.
  --subst_type {withSubst,withoutSubst}
                        Whether to plot data with or without substitutions.
  --node_freq_range NODE_FREQ_RANGE NODE_FREQ_RANGE
                        Min and max frequency to determine node size.Higher frequencies are clipped to this value.
  --node_px_range NODE_PX_RANGE NODE_PX_RANGE
                        Largest node size as determined by the frequency.
  --node_outline_scale NODE_OUTLINE_SCALE
                        How much to scale the node outline width (thickness). Values > 1 increase the width; values < 1 decrease the width.
  --node_comparison_colors NODE_COMPARISON_COLORS NODE_COMPARISON_COLORS
                        The colors to use in the gradient when the node colors show the frequency ratio of two experiments. May be specified in hex      
                        (e.g., "#FF0000" for red) or with recognized keywords such as "red", "blue", "green".
  --node_reference_outline_color NODE_REFERENCE_OUTLINE_COLOR
                        Color to make the reference node outline. May be specified in hex (e.g., "#FF0000" for red) or with recognized keywords such as  
                        "red", "blue", "green".
  --node_outline_color NODE_OUTLINE_COLOR
                        Color to make the default node outline. May be specified in hex (e.g., "#FF0000" for red) or with recognized keywords such as    
                        "red", "blue", "green".
  --node_fill_color NODE_FILL_COLOR
                        Color to make the default node fill. May be specified in hex (e.g., "#FF0000" for red) or with recognized keywords such as       
                        "red", "blue", "green".
  --variation_types {insertion,deletion,substitution,mixed,none} [{insertion,deletion,substitution,mixed,none} ...]
                        The variation types that should be included in the graph. This should be a list of the types: "insertion", "deletion",
                        "substitution", "mixed", "none". Default value: "insertion", "deletion", "none". "mixed" means nodes that have multiples
                        variation types (e.g. insertions and substitutions). "none" means the reference node (no variations). May be specified in hex    
                        (e.g., "#FF0000" for red) or with recognized keywords such as "red", "blue", "green".
  --variation_type_colors VARIATION_TYPE_COLORS VARIATION_TYPE_COLORS VARIATION_TYPE_COLORS VARIATION_TYPE_COLORS VARIATION_TYPE_COLORS
                        The colors for the different variations types. They must be specified in the order: INSERTION DELETION SUBSTITUTION MIXED NONE.  
                        MIXED is the color for nodes with multiple types of variations (e.g. insertions and substitutions); NONE is the color for the    
                        reference node (no variations). May be specified in hex (e.g., "#FF0000" for red) or with recognized keywords such as "red",     
                        "blue", "green".
  --edge_show EDGE_SHOW [EDGE_SHOW ...]
                        Whether to show edges between nodes.
  --edge_types EDGE_TYPES [EDGE_TYPES ...]
                        The edge types to show.
  --edge_scale EDGE_SCALE
                        How much to scale the edges width (thickness). Values > 1 increase the width; values < 1 decrease the width.
  --stats               If present, show graph summary statistics in the left margin.
  --width_px WIDTH_PX   The width of the plot in pixels.
  --height_px HEIGHT_PX
                        The height of the plot in pixels.
  --margin_top_px MARGIN_TOP_PX
                        The size of the top margin in pixels.
  --margin_bottom_px MARGIN_BOTTOM_PX
                        The size of the bottom margin in pixels.
  --margin_left_px MARGIN_LEFT_PX
                        The size of the left margin in pixels.
  --margin_right_px MARGIN_RIGHT_PX
                        The size of the right margin in pixels.
  --line_width_scale LINE_WIDTH_SCALE
                        How much to scale the line widths (aka thickness). Values > 1 increase the width; values < 1 decrease the width.
  --font_size_scale FONT_SIZE_SCALE
                        How much to scale the font size. Values > 1 increase the font size; values < 1 decrease it.
  --crop_x CROP_X CROP_X
                        Range of the horizontal dimension to crop. Specified with normalized coords in range [0, 1].
  --crop_y CROP_Y CROP_Y
                        Range of the vertical dimension to crop. Specified in normalized coords in range [0, 1].
  --range_x RANGE_X RANGE_X
                        Range of x-axis for plotting.If not specified chosen automatically to either show all nodes or a preset value for the layout.    
  --range_y RANGE_Y RANGE_Y
                        Range of y-axis for plotting.If not specified chosen automatically to either show all nodes or a preset value for the layout.    
  --legend              Whether to show a legend on the figure.
  --legend_colorbar_scale LEGEND_COLORBAR_SCALE
                        How much to scale the legend color bar (for freq ratio coloring).
  --legend_spacing_px LEGEND_SPACING_PX
                        Amount of vertical space in pixels between different legends.
  --separate_components
                        If present, separate the connected components of the graph.
  --interactive         If present opens the interactive version in a browser. Uses the Ploty library figure.show() function to do so.
  --reverse_complement {0,1} [{0,1} ...]
                        Whether to reverse complement the sequences in the data sets. If present, the number of values must be the same as the number    
                        of input directories. "1" mean reverse complement the sequence and "0" means do not. Used for making a layout for data sets      
                        that have reference sequences that are the reverse complements of each other. If "1" also uses the reverse complement of
                        sequences when determining the display labels and hover text. This affects the universal layout and fractal layout.
```

### `histogram.py`

## Figures

Here we describe the visual properties of the output graphs.

### Layouts

Different *layouts*, specified with the XX parameter, are used to position the vertices. We decribe a subset of them in the following.

* Kamada-Kawaii: a physical simulation algorithm for laying out graphs (Kamada and Kawai, 1989). The graph is modeled as a spring system, where each pair of vertices is connected by a spring of resting length proportional to the shortest-path distance between the two vertices. The algorithms seeks to minimize the energy of the system by finding the optimal positions of the vertices. The result is to bring the vertices to positions such that the Euclidean distance between pairs of vertices is close to the shortest-path distance between them. This is intended to layout the vertices in an aesthetically pleasing way that also reveals graph structure, such as clusters of vertices that are close together. The algorithm is implemented by the NetworkX Python package (Hagberg, 2008). As of this writing, NetworkX places the vertices in a circle for the initialization. Note, the shortest-path distance between two vertices is at least the Levenshtein distance between the two vertices, though equality is not necessary. The final coordinates of the vertices are always in [0, 1] (relevant if setting the `plot_range_x` and `plot_range_y` parameters).
* Universal layout: We call this layout the *universal* layout because every possible insertion and deletion is assigned a predefined position, independent of the subset of in/dels occurring in the dataset. This contrasts with the Kamada-Kawaii algorithm, which has a randomized initialization and positions vertices in a data-dependent manner. While this layout algorithm is less flexible and does not always position vertices based on their Levenshtein distances to one another, it has the advantage that it is reproducible and easily interpretable. The universal layout positions the reference vertex in the center of the figure, the insertion vertices above the reference, and the deletion vertices below the reference. For both the insertion and deletion vertices, the $y$-coordinate is determined by the vertices' Levenshtein distance from the reference, which is simply the number of nucleotides inserted or deleted in the vertices' sequence. If $d$ is a vertex's Levenshetin distance from the reference, the deletion vertices are placed $d$ units below the reference and insertion vertices are place $d$ units above. To determine the $x$-coordinates of vertices, we use different rules for insertions and deletions. Insertion vertices are ordered from left to right by alphabetical order, while deletion vertices are ordered from left to right based on the positions of the deleted vertices with respect to the DSB site. Figure XX gives a key to the vertex placement for the universal layout. In this layout, the vertex coordinates usually range between $\pm 2 * window\_size$, where `window_size` is the parameter from [preprocess.py](#preprocesspy) (relevant if setting the `plot_range_x` and `plot_range_y` parameters).
* Radial: This layout arranges vertices in concentric circles around the reference sequence. The reference vertex is placed at the center of the figure; insertion vertices are placed above the reference; and deletion vertices are placed below the the reference. A vertex's physical distance to the reference is proportional to its Levenshtein distance to the reference. That is, the first concentric circle around the reference is at Levenshtein distance 1, the next is at Levenshtein distance 2, and so on. In each concentric circle, the insertion vertices are arranged clockwise from highest to lowest frequency, and deletion vertices are arranged counterclockwise from highest to lowest frequency. This means that the vertices for both insertions and deletion are placed left to right from highest to lowest frequency. Several heuristics are used to perturb the vertices slightly so that they do not overlap and edges are not collinear. In this layout, the vertex coordinates (with default scaling) usually range between $\pm 4 * window\_size$, where `window_size` is the parameter from [preprocess.py](#preprocesspy) (relevant if setting the `plot_range_x` and `plot_range_y` parameters). However, the scaling on the insertion and deletion, $x$ and $y$-coordinates can be controlled by the FIXME `--universal_layout_x_scale_insertion`, FIXME parameters.
* Fractal. This layout is intended to display the insertion vertices only. The vertices with a single inserted nucleotide (A, C, G, or T) are placed in the center of the corresponding quadrants shown in Figure XX: A in the top left, C in the top right, G in the bottom left, and T in the bottom right. Vertices with two inserted nucleotides, and placed in the center of the corresponding square shown in Figure XX. Note the pattern, the squares for two nucleotides insertions starting with A are obtained by subdividing the A square into four quadrants and arranging AA, AC, AG, AT based on their second nucleotide in the same order as the original A, C, G, T squares: AA in the top left, AC in the top right, AG in the bottom left, and AT in the bottom right. This pattern can be continued for arbitrarily large kmers. The square for kmer $s_1 s_2 \cdots s_k$ will be obtained by subdividing the the square for $s_1 s_2 \cdots s_{k-1}$ into four quadrants and choosing the quadrant in the same manner as before, based on $s_k$. Note, this is layout is very similar to the insertion of the universal layout, except that the both dimensions are used to order the inserted nucleotides alphabetically instead of linearly from left to right. The final coordinates of the vertices are always in [-1, 1] (relevant if setting the `plot_range_x` and `plot_range_y` parameters).

### Aesthetics

SOME OF THIS MIGHT BE UNNECESSARY! ANYTHING DESCRIBED IN THE ARGPARASE HELP SHOULD NOT BE REPEATED!

* Figure size: The height and width of the figure in pixels is controlled by the FIXME parameters. The margins of the figure in pixels is controlled by the FIXME parameteres.
* Scaling: Different elements of the figure can be independently scaled using the following parameters.
  * `--font_size_scale`: scales the defaut font sizes.
  * `Line width`: FIXME Scales the width of edges between vertices, among other things.
  * `Node outline width`: FIXME Scales the outline width of vertices.
  * FIXME: SHOULD INCLUDE THE UNIVERSAL LAYOUT ARGS?
* Coordinates: In the Kamada layout, the $x$- and $y$-coordinates are not interpretable, though the distance between vertices is intended to approximate their shortest-path distance in the graph. For the radial, universal, and fractal layouts please see the description of the [layouts](#layouts) for the interpretation of the $x$- and $y$-coordinates.
* Vertices:
  * Size: For an individual graph, the vertex radius is a function of the log-frequency. This is controlled by the parameters FIXME LIST THEM. Vertices with frequencies $\leq$ min_freq get min_radius and frequencies $\geq$ `MAX_FREQ` get `MAX_RADIUS`. Frequencies in between get a radius that varies linearly with the log-frequency.
  * Outline: the outline color is intended to differentiate the reference sequence vertex from the other vertices. The colors are determined by the XXX and XXX parameters.
  * Color: There color depends on whether the graph being plotted is an individual experiment or a compairison of experiments:
    * Variation type (`variation_type`): Only for individual graphs. Gives vertices a color depending on the type of variations present in their alignment: insertions only, deletions only, substitutions only, mixed (a combinations of insertions, deletions, and substitution), and none (the reference sequence). These colors may be specified in the `XXX` variable.
    * Frequency ratio (continuous): Only for comparison graphs. Indicate "continuous" in the XX parameter. Uses a color gradient that indicates the ratio of frequencies in the two experiments, for each vertex. The ratio is computed by putting the first sample in the numerator and second sample in the denominator (this order is determined in the [preprocessing](#preprocesspy) stage). The colors at the extreme ends of the gradient are determined by the XX parameter, and values in the middle are smoothly interpolated. The ratios at the extreme ends of the gradient are determined by the XX parameter.
    * Frequency ratio (discrete): Only for comparison graphs. Indicate "discrete" in the XX parameter. Uses three colors to indicate the ratio of frequencies in the two experiments, for each vertex. If the frequency ratio is less than the first value in XX, the second color in XX is displayed (meaning higher in sample 2); If the frequency ratio is greater than the second value in XX, the first color in XX is displayed (meaning higher in sample 1); if the frequency ratio is between the values in XX (inclusive), the vertex is colored white.
* Edges: All edges indicate a 1-nucleotide variation (insertion, deletion, or substitution) between the two vertices. Whether to show and which types of edges to show can be controlled by the FIXME parameters.
* Legends: The legends describe the vertex size, vertex outline, vertex color, and the edges. These can be all drawn by using the `--legend` parameter. The legends will be drawn in the right margin of the figure. To ensure enough room, use the `--margin_*` parameters.
* Title: A title can be optionally added to the top margin of the figure using the `--title` parameter.
* Universal layout axes: To illustrate the vertex placement for the universal layout, axes can be drawn on the figure. Using the `FIXME` parameter, a horizontal axis can be drawn on the deletion side (below the reference sequence)to show the midpoint position of the range of deletions that lies at that x-coordinate. Using the `FIXME` parameter, a horizontal axis can be drawn on the insertion side (above the reference sequence) to show the alphabetical order of the nucleotides. Using the `FIXME` parameter, a vertical axis can be drawn to show the Levenshtein distance of vertices from the reference vertex. To control the $x$ or $y$-coordinate of the axes, use the `FIXME` parameters.
* Graph statistics: Use the `--stats` parameter to display summary statistics of the graph in the left margin of the figure. If the graph has more than one connected component, only the component with the reference sequence is used, since otherwise the pairwise distance between vertices is not well-defined. To ensure the margin has enough space, use the `FIXME` parameters.

NEED TO MENTION SOMEWHERE THAT "SEQUENCE" AND "VERTEX" CAN BE USED INTERCHANGEABLY?
NEED TO DESCRIBE THE WITH AND WITHOUT SUBST PARAMETER.

## Tutorial

Several example input files are available in the `data_input` directory:
* `data_input/ref_seq`: reference sequence FASTA files, representing the perfect repaired sequence for different samples.
* `data_input/fastq`: high-throughput sequencing data for different samples. Note, the FASTQ samples have been *trimmed*, meaning that we only capture the portion of the read between the primers and low-quality reads have been already filtered out.

The alignment done between the reads in the input FASTQ and the reference sequence expects the first base of each read to align with the first base of the reference sequence. Therefore, the reads must be trimmed and the reference sequence selected in such a way that their first base pairs align. In the experiments performed in the study by Jeon et al. (LINK), the primers were designed to be about 50-150 base pairs away from the induced DSB site. In principle, this would allow reads repaired by NHEJ to have variations nears the DSB (within +/- 10 base pairs), and allow them remaining sequence to the primers to otherwise match the reference perfectly (not counting substitution errors due to sequencing or library preparation). See diagram XX for more clarification.