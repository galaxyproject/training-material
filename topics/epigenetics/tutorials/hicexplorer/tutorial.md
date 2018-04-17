---
layout: tutorial_hands_on
topic_name: epigenetics
tutorial_name: hicexplorer
---

# Introduction
{:.no_toc}

In this HiCExplorer tutorial we will generate and plot a Hi-C contact matrix.
For this the following steps are necessary to be performed:
1. Map the Hi-C reads to the reference genome
2. Creation of a Hi-C matrix
3. Plotting the Hi-C matrix
4. Correction of Hi-C matrix
5. TAD Calling

After a corrected Hi-C matrix is created other tools can be used to visualize it, call TADS or compare it with other matrices.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. Mapping Hi-C data
> 2. Create a contact matrix
> 3. Learn how to correct a contact matrix
> 4. Visualize a contact matrix
> 5. Compute and visualize topological associating domains (TADs)
> {:toc}
>
{: .agenda}

## Data upload

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history
> 2. Import from [Zenodo](https://doi.org/10.5281/zenodo.1183661).
>
>    > ### {% icon tip %} Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > Below are the links to the read files that can be copied and pasted in the upload manager.
>    >
>    > ```
>    > https://zenodo.org/record/1183661/files/HiC_S2_1p_10min_lowU_R1.fastq.gz
>    > https://zenodo.org/record/1183661/files/HiC_S2_1p_10min_lowU_R2.fastq.gz
>    > ```
>    >
>    > * Paste the links into the text field
>    > * Press **Start**
>    {: .tip}
>
> 3. Rename the data set to something meaningful, e.g. `HiC_S2_1p_10min_lowU_R1` and `HiC_S2_1p_10min_lowU_R2`.
> By default, when data is imported via its link, Galaxy names it with its URL.
>
{: .hands_on}

# Reads mapping

Mates have to be mapped individually to avoid mapper specific heuristics designed for standard paired-end libraries.

We have used the HiCExplorer successfully with bwa, bowtie2 and hisat2. In this tutorial we will be using **Map with BWA-MEM** tool. It is important to remember to:
- use local mapping, in contrast to end-to-end. A fraction of Hi-C reads are chimeric and will not map end-to-end thus, local mapping is important to increase the number of mapped reads
- tune the aligner parameters to penalize deletions and insertions. This is important to avoid aligned reads with gaps if they happen to be chimeric.
- If bowtie2 or hisat2 are used, `--reorder` option and as a file format `bam_native` needs to be used. Regular `bam` files are sorted by Galaxy and can not be used as an input for HiCExplorer.

> ### {% icon hands_on %} Hands-on: Mapping reads
>
> 1. **Map with BWA-MEM 0.8.0** {% icon tool %}: Run Map with BWA-MEM on both strands `HiC_S2_1p_10min_lowU_R1` and `HiC_S2_1p_10min_lowU_R2` with:
>    - "Will you select a reference genome from your history or use a built-in index?" to `Use a built-in index`
>    - "Select a reference genome" to `dm3`
>    - "Is this library mate-paired?" to `Single-end or interleaved paired-end`
>    - Set multiple data sets
>    - "FASTQ file" to `HiC_S2_1p_10min_lowU_R1`and `HiC_S2_1p_10min_lowU_R2`
>    - "BWA settings to use" to `Full parameter List`
>    - "Gap extension penalty (-E)" to `50`
>    - "Penalty for clipping (-L)" to `0`
>
> 2. Rename the output of the tool according to the corresponding files: `R1.sam` and `R2.sam`
>
{: .hands_on}

# Creation of a Hi-C matrix

Once the reads have been mapped the Hi-C matrix can be built.

For this step we will use [hicBuildMatrix](http://hicexplorer.readthedocs.io/en/latest/content/tools/hicBuildMatrix.html#hicbuildmatrix) tool, which builds the matrix of read counts over the bins in the genome, considering the sites around the given restriction site.

> ### {% icon hands_on %} Hands-on: hicBuildMatrix
>
> 1. **hicBuildMatrix** {% icon tool %}: Run hicBuildMatrix on the `R1.sam` and `R2.sam` from previous step with modifying the following parameters:
>    - "1: Sam/Bam files to process" to `R1.sam`
>    - "2: Sam/Bam files to process" to `R2.sam`
>    - "Choose to use a restriction cut file or a bin size" to `Bin size`
>    - "Bin size in bp" to `10000`
>    - "Sequence of the restriction site" to `GATC`
>
>       > ### {% icon comment %} Comment
>       >
>       > *hicBuildMatrix* creates two files, a bam file containing only the valid Hi-C read pairs and a matrix containing the Hi-C contacts at the given resolution. The bam file is useful to check the quality of the Hi-C library on the genome browser. A good Hi-C library should contain piles of reads near the restriction fragment sites. In the QC folder a html file is saved with plots containing useful information for the quality control of the Hi-C sample like the number of valid pairs, duplicated pairs, self-ligations etc. Usually, only 25%-40% of the reads are valid and used to build the Hi-C matrix mostly because of the reads that are on repetitive regions that need to be discarded.
>       {: .comment}
>
>       > ### {% icon comment %} Comment
>       >
>       > Normally 25% of the total reads are selected. The output matrices have counts for the genomic regions. The extension of output matrix files is .h5.
>       > A quality report is created in e.g. `hicMatrix/R1_10kb_QC`, have a look at the report hicQC.html.
>       {: .comment}
>
> 2. Rename the output to `10 kb contact matrix`.
>
{: .hands_on}

# Plotting the Hi-C matrix

A 10kb bin matrix is too large to plot, it's better to reduce the resolution. We usually run out of memory for a 1 kb or a 10 kb bin matrix and the time to plot it is very long (minutes instead of seconds). In order to reduce the resolution we use the tool [hicMergeMatrixBins](http://hicexplorer.readthedocs.io/en/latest/content/tools/hicMergeMatrixBins.html#hicmergematrixbins).

[hicMergeMatrixBins](http://hicexplorer.readthedocs.io/en/latest/content/tools/hicMergeMatrixBins.html#hicmergematrixbins) merges the bins into larger bins of given number (specified by –numBins). We will merge 100 bins in the original (uncorrected) matrix and then correct it. The new bin size is going to be 10.000 bp * 100 = 1.000.000 bp = 1 Mb

> ### {% icon hands_on %} Hands-on: hicMergeMatrixBins
>
> 1. **hicMergeMatrixBins** {% icon tool %}: Run hicMergeMatrixBins on the output from previous step setting the following parameter:
>    - "Number of bins to merge" to `100`
>
> 2. Rename the output to `1 MB contact matrix`.
>
> 3. **hicPlotMatrix** {% icon tool %}: Run hicPlotMatrix on the output from hicMergeMatrixBins `1 MB contact matrix` adjusting the parameters:
>    - "Plot title" to `Hi-C matrix for dm3`
>    - "Remove masked bins from the matrix" to `True`
>    - "Plot the log1p of the matrix values: `True`
>    - "Chromosomes to include (and order to plot in)" to `chr2L`
>    - "+ Insert Chromosomes to include (and order to plot in):" to `chr2R`
>    - "+ Insert Chromosomes to include (and order to plot in):" to `chr3L`
>    - "+ Insert Chromosomes to include (and order to plot in):" to `chr3R`
>    - "+ Insert Chromosomes to include (and order to plot in):" to `chrX`
>
>    > ### {% icon tip %} Tip: log1p
>    >
>    > Because of the large differences in counts found in the matrix, it is better to plot the counts using the *–log1p* option.
>    {: .tip}
>
{: .hands_on}

The resulting plot of the 1 Mb contact matrix should look like: 
![TAD plot](../../images/plotMatrix.png)

# Correction of Hi-C matrix


[hicCorrectMatrix](http://hicexplorer.readthedocs.io/en/latest/content/tools/hicCorrectMatrix.html#hiccorrectmatrix) corrects the matrix counts in an iterative manner. For correcting the matrix, it’s important to remove the unassembled scaffolds (e.g. NT_) and keep only chromosomes, as scaffolds create problems with matrix correction. Therefore we use the chromosome names (chr2R, chr2L, chr3R, chr3L, chrX) here.


Matrix correction works in two steps: first a histogram containing the sum of contact per bin (row sum) is produced. This plot needs to be inspected to decide the best threshold for removing bins with lower number of reads. The second steps removes the low scoring bins and does the correction.

> ### {% icon hands_on %} Hands-on: Matrix diagnostic
> 
> 1. **hicCorrectMatrix** {% icon tool %}: Run hicCorrectMatrix on the output from hicBuildMatrix `10 kb contact matrix` adjusting the parameters:
>    - "Range restriction (in bp)" to `Diagnostic plot`
>    - "Chromosomes to include (and order to plot in)" to `chr2L`
>    - "+ Insert Chromosomes to include (and order to plot in):" to `chr2R`
>    - "+ Insert Chromosomes to include (and order to plot in):" to `chr3L`
>    - "+ Insert Chromosomes to include (and order to plot in):" to `chr3R`
>    - "+ Insert Chromosomes to include (and order to plot in):" to `chrX`
>
{: .hands_on}

The output of the program prints a threshold suggestion that is usually accurate but is better to revise the histogram plot. The threshold is visualized in the plot as a black vertical line.

In our case the distribution describes the counts per bin of a genomic distance. To remove all bins with a z-score threshold less / more than X means to remove all bins which have less / more counts than X of mean of their specific distribution in units of the standard deviation. Looking at the distribution, we can select the value of -1.6 (lower end) and 1.8 (upper end) to remove. This is given by the –filterThreshold option in hicCorrectMatrix set to 'correct matrix' mode.

> ### {% icon hands_on %} Hands-on: Matrix correction
>
> 1. **hicCorrectMatrix** {% icon tool %}: Run hicCorrectMatrix on the original matrix `10 kb contact matrix` adjusting the parameters:
>    - "Range restriction (in bp)" to `Correct matrix`
>    - "Normalize each chromosome separately" to `True`
>    - "Remove bins of low coverage" to `-1.6`
>    - "Remove bins of large coverage" to `1.8`
>    - "Chromosomes to include (and order to plot in)" to `chr2L`
>    - "+ Insert Chromosomes to include (and order to plot in):" to `chr2R`
>    - "+ Insert Chromosomes to include (and order to plot in):" to `chr3L`
>    - "+ Insert Chromosomes to include (and order to plot in):" to `chr3R`
>    - "+ Insert Chromosomes to include (and order to plot in):" to `chrX`
>
> 2. Rename the corrected matrix to `10 kb corrected contact matrix`.
>
{: .hands_on}

It can happen that the correction stops with:

```
ERROR:iterative correction:*Error* matrix correction produced extremely large values.
This is often caused by bins of low counts. Use a more stringent filtering of bins.
```

This can be solved by a more stringent z-score values for the filter threshold or by a look at the plotted matrix. For example, chromosomes with 0 reads in its bins can be excluded from the correction by not defining it for the set of chromosomes that should be corrected (parameter 'Include chromosomes').

### Plotting the corrected Hi-C matrix

We can now plot chromosome 2L with the corrected matrix.

> ### {% icon hands_on %} Hands-on: Plotting the corrected Hi-C matrix
>
> 1. **hicPlotMatrix** {% icon tool %}: Run hicPlotMatrix on `10 kb corrected contact matrix` adjusting the parameters:
>    - "Plot title" to `Hi-C matrix for dm3`
>    - "Plot per chromosome" to `False`
>    - "Plot only this region" to `chr2L`
>    - "Plot the log1p of the matrix values" to `True`
>
{: .hands_on}


### Load new data

The steps so far would have led to long run times if real data would have been used. We therefore prepared a new matrix for you, `corrected contact matrix dm3 large`. Please load it into your history.

# TAD calling

“The partitioning of chromosomes into topologically associating domains (TADs) is an emerging concept that is reshaping our understanding of gene regulation in the context of physical organization of the genome” [Ramirez et al. 2017](https://doi.org/10.1101/115063).

TAD calling works in two steps: First HiCExplorer computes a TAD-separation score based on a z-score matrix for all bins. Then those bins having a local minimum of the TAD-separation score are evaluated with respect to the surrounding bins to assign a p-value. Then a cutoff is applied to select the bins more likely to be TAD boundaries.

[hicFindTADs](http://hicexplorer.readthedocs.io/en/latest/content/tools/hicFindTADs.html#hicfindtads) tries to identify sensible parameters but those can be change to identify more stringent set of boundaries.

> ### {% icon hands_on %} Hands-on: Finding TADs
>
> 1. **hicFindTADs** {% icon tool %}: Run hicFindTADs on `corrected contact matrix dm3 large` adjusting the parameters:
>    - "Minimum window length (in bp) to be considered to the left and to the right of each Hi-C bin." to `30000`
>    - "Maximum window length (in bp) to be considered to the left and to the right of each Hi-C bin." to `100000`
>    - "Step size when moving from minDepth to maxDepth" to `10000`
>    - "Multiple Testing Corrections" to `False discovery rate`
>    - "q-value" to `0.05` 
>    - "Minimum threshold of the difference between the TAD-separation score of a putative boundary and the mean of the TAD-sep. score of surrounding bins." to `0.001`
>
{: .hands_on}

As an output we get the boundaries, domains and scores separated files. We will use in the plot later only the TAD-score file.


> ### {% icon hands_on %} Hands-on: Computing A / B compartments
>
> 1. **hicPCA** {% icon tool %}: Run hicPCA adjusting the parameters:
>    - "Matrix to compute on" to `corrected contact matrix dm3 large`
>    - "Output file format" to `bigwig`
>
{: .hands_on}

# Integrating Hi-C and other data

We can plot the TADs for a given chromosomal region. For this we will use [hicPlotTADs](http://hicexplorer.readthedocs.io/en/latest/content/tools/hicPlotTADs.html). But before make sure to import [gene track file](https://zenodo.org/record/1176070/files/dm6_genes.bed) in .bed format from [Zenodo](https://doi.org/10.5281/zenodo.1176070).

For the next step we need additional data tracks. Please load `dm3_genes.bed`, `H3K27me3.bw`, `H3K36me3.bw` and `H4K16ac.bw` to your history.

> ### {% icon hands_on %} Hands-on: Plotting TADs
>
> 1. **hicPlotTADs** {% icon tool %}: Run hicPlotTADs adjusting the parameters:
>    - "Region of the genome to limit the operation" to `chr2L:14500000-16500000`
>    - "Choose style of the track" to `TAD visualization`
>         - "Plot title" to `HiC dm3 chr2L:14500000-16500000`
>         - "Matrix to compute on." to the corrected matrix from hicCorrectMatrix step
>         - "Depth" to `750000`
>         - "Width" to `4`
>         - "Boundaries file" to `hicFindTads on data XX: TAD domains`
>         - "Show x labels" to `Yes`
>
>    - "+Insert Include tracks in your plot"
>        - "Choose style of the track" to `Bedgraph matrix track`
>        - "Plot title" to `TAD separation score`
>        - "Track file bedgraph format" to `hicFindTads on data XX: TAD information in bm file`
>        - "Width" to `4`
>        - "Set type to lines" to 'True'
>
>
>    - "+Insert Include tracks in your plot"
>        - "Choose style of the track" to `Bigwig track`
>        - "Plot title" to `PC1`
>        - "Track file bigwig format" the first computed `hicPCA` result
>        - "Minimum value" to `-0.03`
>        - "Maximum value" to `0.03`
>        - "Width" to `1.5`
>        - "Color of track" to a color of your choice
>
>    - "+Insert Include tracks in your plot"
>        - "Choose style of the track" to `Bigwig track`
>        - "Plot title" to `PC2`
>        - "Track file bigwig format" the second computed `hicPCA` result
>        - "Minimum value" to `-0.03`
>        - "Maximum value" to `0.03`
>        - "Width" to `1.5`
>        - "Color of track" to a color of your choice
>
>    - "+Insert Include tracks in your plot"
>        - "Choose style of the track" to `Bigwig track`
>        - "Plot title" to `H3K36me3`
>        - "Track file bigwig format" to `H3K36me3`
>        - "Width" to `1.5`
>        - "Color of track" to a color of your choice
>
>    - "+Insert Include tracks in your plot"
>        - "Choose style of the track" to `Bigwig track`
>        - "Plot title" to `H3K27me3`
>        - "Track file bigwig format" to `H3K27me3`
>        - "Width" to `1.5`
>        - "Color of track" to a color of your choice
>
>    - "+Insert Include tracks in your plot"
>        - "Choose style of the track" to `Bigwig track`
>        - "Plot title" to `H4K16ac`
>        - "Track file bigwig format" to `H4K16ac`
>        - "Width" to `1.5`
>        - "Color of track" to a color of your choice
>
>    - "+Insert Include tracks in your plot"
>        - "Choose style of the track" to `Gene track`
>        - "Plot title" to `dm3 genes`
>        - "Track file bedgraph format" the imported .bed file
>        - "Width" to `3`
>        - "Type" to `genes`
>        - "Gene rows" to `15`
>        - "Color of track" to a color of your choice
>
{: .hands_on}

The resulting image should look like this one: 
![TAD plot](../../images/plotTADs.png)
# Conclusion
{:.no_toc}

In this tutorial we used HiCExplorer to analyze drosophila melanogaster cells. We mapped the chimeric reads and created a contact matrix, to reduce noise this contact matrix was normalized. We showed how to visualize a contact matrix and how we can investigate topological associating domains and relate them to additional data like gene tracks.


 To improve your learned skills we offer an additional tutorial based on mouse stem cells: [following work](http://hicexplorer.readthedocs.io/en/latest/content/example_usage.html).
