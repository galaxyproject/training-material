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
This tutorial is based on the [following work](http://hicexplorer.readthedocs.io/en/latest/content/example_usage.html).

We will be using a Hi-C dataset from [Marks et. al. 2015](http://www.genomebiology.com/2015/16/1/149), on mouse ESCs.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Data upload

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history
> 2. Import the [SRR1956527_1.fastq](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/007/SRR1956527/SRR1956527_1.fastq.gz), [SRR1956527_2.fastq](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/007/SRR1956527/SRR1956527_2.fastq.gz), [SRR1956528_1](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/008/SRR1956528/SRR1956528_1.fastq.gz), [SRR1956528_2](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/008/SRR1956528/SRR1956528_2.fastq.gz), [SRR1956529_1](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/009/SRR1956529/SRR1956529_1.fastq.gz) and [SRR1956529_2](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/009/SRR1956529/SRR1956529_2.fastq.gz).
>
>    > ### {% icon tip %} Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > Below are the links to the read files that can be copied and pasted in the upload manager.
>    >
>    > ``` ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/007/SRR1956527/SRR1956527_1.fastq.gz
>    > ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/007/SRR1956527/SRR1956527_2.fastq.gz
>    > ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/008/SRR1956528/SRR1956528_1.fastq.gz
>    > ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/008/SRR1956528/SRR1956528_2.fastq.gz
>    > ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/009/SRR1956529/SRR1956529_1.fastq.gz
>    > ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/009/SRR1956529/SRR1956529_2.fastq.gz
>    >```
>    >
>    > * Paste the links into the text field
>    > * Press **Start**
>    {: .tip}
>
>    > ### {% icon tip %} Tip: Change the file type `fastq` to `fastqsanger` once the data file is in your history
>    >
>    > * Click on the pencil button displayed in your data file in the history
>    > * Choose **Datatype** on the top
>    > * Select `fastqsanger`
>    > * Press **save**
>    {: .tip}
>
> 3. Rename the dataset to something meaningful
> By default, when data is imported via its link, Galaxy names it with its URL.
>
{: .hands_on}

# Reads mapping

Mates have to be mapped individually to avoid mapper specific heuristics designed for standard paired-end libraries.

We have used the HiCExplorer sucessfuly with bwa, bowtie2 and hisat2. In this tutorial we will be using **Map with BWA-MEM** tool. It is important to remember to:
- use local mapping, in contrast to end-to-end. A fraction of Hi-C reads are chimeric and will not map end-to-end thus, local mapping is important to increase the number of mapped reads
- tune the aligner parameters to penalize deletions and insertions. This is important to avoid aligned reads with gaps if they happen to be chimeric.

> ### {% icon hands_on %} Hands-on: Mapping reads
>
> 1. **Map with BWA-MEM** {% icon tool %}: Run Map with BWA-MEM on `SRR1956527_1` with:
>    - "Will you select a reference genome from your history or use a built-in index?" to `Use a built-in genome index`
>    - "Using reference genome" to `mm10`
>    - "Single or Paired-end reads" to `Single`
>    - "Select fastq dataset" to `SRR1956527_1`
>    - "Select analysis mode" to `Full list of options`
>    - "Set algorithmic options?" to `Set`
>    - "Skip seeds with more than that many occurrences" to `10000`
>    - "Set scoring options?" to `Set`
>    - "Gap extension penalties; a gap of size k cost '-O + -E*k'. If two numbers are specified, the first is the penalty of extending a deletion and the second for extending an insertion" to `50,50`
>    - "Penalties for 5'-end and 3'-end clipping" to `0`
>
> 2. Rename the output of the tool according to the corresponding file: `SRR1956527_1.sam`
>
> 3. **Map with BWA-MEM** {% icon tool %}: Run Map with BWA-MEM on the remaining imported files individually(`SRR1956527_2`, `SRR1956528_1`, `SRR1956528_2`, `SRR1956529_1` and `SRR1956529_2`) *individually* with the same settings as during step 1.
>
{: .hands_on}

# Creation of a Hi-C matrix

Once the reads have been mapped the Hi-C matrix can be built. In the following we will create three Hi-C matrices and merge them to one.

For this step we will use [hicBuildMatrix](http://hicexplorer.readthedocs.io/en/latest/content/tools/hicBuildMatrix.html#hicbuildmatrix) tool, which builds the matrix of read counts over the bins in the genome, considering the sites around the given restriction site.
In order to increase the depth of reads we will merge the counts from these three replicates. We will be using *hicSumMatrices* tool.

> ### {% icon hands_on %} Hands-on: hicBuildMatrix
>
> 1. **hicBuildMatrix** {% icon tool %}: Run hicBuildMatrix on the `SRR1956527_1.sam` and `SRR1956527_2.sam` from previous step with modifying the following parameters:
>    - "1: Sam/Bam files to process" to `SRR1956527_1.sam`
>    - "2: Sam/Bam files to process" to `SRR1956527_2.sam`
>    - "Choose to use a restriction cut file or a bin size" to `Bin size`
>    - "Bin size in bp" to `1000`
>    - "Sequence of the restriction site" to `GATC`
>
>       > ### {% icon comment %} Comment
>       >
>       > *hicBuildMatrix* creates two files, a bam file containing only the valid Hi-C read pairs and a matrix containing the Hi-C contacts at the given resolution. The bam file is useful to check the quality of the Hi-C library on the genome browser. A good Hi-C library should contain piles of reads near the restriction fragment sites. In the QCfolder a html file is saved with plots containing useful information for the quality control of the Hi-C sample like the number of valid pairs, duplicated pairs, self-ligations etc. Usually, only 25%-40% of the reads are valid and used to build the Hi-C matrix mostly because of the reads that are on repetitive regions that need to be discarded.
>       {: .comment}
> 2. **hicBuildMatrix** {% icon tool %}: Re-run hicBuildMatrix on `SRR1956528_1.sam` and `SRR1956528_2.sam`, `SRR1956529_1.sam` and `SRR1956529_2.sam` respectively.
>
>    > ### {% icon question %} Question
>    >
>    > How many selected reads do the output bam files have?
>    >
>    > <details>
>    > <summary>Click to view answers</summary>
>    > The output bam files show that we have around 34M, 54M and 58M selected reads for `SRR1956527`, `SRR1956528` & `SRR1956529`, respectively.
>    > </details>
>    {: .question}
>
>       > ### {% icon comment %} Comment
>       >
>       > Normally 25% of the total reads are selected. The output matrices have counts for the genomic regions. The extension of output matrix files is .h5.
>       > A quality report is created in e.g. `hicMatrix/SRR1956527_10kb_QC`, have a look at the report hicQC.html.
>       {: .comment}
>
> 3. **hicSumMatrices** {% icon tool %}: Run hicSumMatrices on `.h5` outputs of hicBuildMatrix tool with default settings.
>
{: .hands_on}

# Plotting the Hi-C matrix

A 10kb bin matrix is quite large to plot and is better to reduce the resolution (to know the size of a Hi-C matrix use the tool [hicInfo](http://hicexplorer.readthedocs.io/en/latest/content/tools/hicInfo.html#hicinfo)), i.e. we usually run out of memory for a 1 kb or a 10 kb matrix and second, the time to plot is very long (minutes instead of seconds). For this we use the tool [hicMergeMatrixBins](http://hicexplorer.readthedocs.io/en/latest/content/tools/hicMergeMatrixBins.html#hicmergematrixbins).

[hicMergeMatrixBins](http://hicexplorer.readthedocs.io/en/latest/content/tools/hicMergeMatrixBins.html#hicmergematrixbins) merges the bins into larger bins of given number (specified by –numBins). We will merge 1000 bins in the original (uncorrected) matrix and then correct it. The new bin size is going to be 10.000 bp * 100 = 1.000.000 bp = 1 Mb

> ### {% icon hands_on %} Hands-on: hicMergeMatrixBins
>
> 1. **hicMergeMatrixBins** {% icon tool %}: Run hicMergeMatrixBins on the output from previous step setting the following parameter:
>    - "Number of bins to merge" to `100`
>
> 2. **hicPlotMatrix** {% icon tool %}: Run hicPlotMatrix on the output from hicMergeMatrixBins adjusting the parameters:
>    - "Plot title" to `Hi-C matrix for mESC`
>    - "Plot per chromosome" to `True`
>    - "Chromosomes to include (and order to plot in)" to `chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY`
>    - "Remove masked bins from the matrix" to `True`
>
>    > ### {% icon tip %} Tip: log1p
>    >
>    > Because of the large differences in counts found int he matrix, it is better to plot the counts using the *–log1p* option.
>    {: .tip}
>
>    - **Plot the log1p of the matrix values**: True
>
{: .hands_on}

# Correction of Hi-C matrix

[hicCorrectMatrix](http://hicexplorer.readthedocs.io/en/latest/content/tools/hicCorrectMatrix.html#hiccorrectmatrix) corrects the matrix counts in an iterative manner. For correcting the matrix, it’s important to remove the unassembled scaffolds (e.g. NT_) and keep only chromosomes, as scaffolds create problems with matrix correction. Therefore we use the chromosome names (1-19, X, Y) here. **Important**: Use ‘chr1 chr2 chr3 etc.’ if your genome index uses chromosome names with the ‘chr’ prefix.

Matrix correction works in two steps: first a histogram containing the sum of contact per bin (row sum) is produced. This plot needs to be inspected to decide the best threshold for removing bins with lower number of reads. The second steps removes the low scoring bins and does the correction.

First we will re-run hicMergeMatrixBins on the hicSumMatrices output and set the bin size to 20.

> ### {% icon hands_on %} Hands-on: Matrix diagnostic
>
> 1. **hicMergeMatrixBins** {% icon tool %}: Run hicMergeMatrixBins on the icSumMatrices output setting the following parameter:
>    - "Number of bins to merge" to `2`
>
> 2. **hicCorrectMatrix** {% icon tool %}: Run hicCorrectMatrix on the output from previous step adjusting the parameters:
>    - "Range restriction (in bp)" to `Diagnostic plot`
>    - "Include chromosomes" to `chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY`
> 
> 3. **hicCorrectMatrix** {% icon tool %}: Now try again but setting chromosomes as:
>    - "Range restriction (in bp)" to `Diagnostic plot`
>    - "Include chromosomes" to `chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY`
>
{: .hands_on}

The output of the program prints a threshold suggestion that is usually accurate but is better to revise the histogram plot. The threshold is visualized in the plot as a black vertical line.

In our case the distribution describes the counts per bin of a genomic distance. To remove all bins with a z-score threshold less / more than X means to remove all bins which have less / more counts than X of mean of their specific distribution in units of the standard deviation. Looking at the distribution, we can select the value of -2 (lower end) and 3 (upper end) to remove. This is given by the –filterThreshold option in hicCorrectMatrix set to 'correct matrix' mode.

> ### {% icon hands_on %} Hands-on: Matrix correction
>
> 1. **hicCorrectMatrix** {% icon tool %}: Run hicCorrectMatrix on the output from previous step adjusting the parameters:
>    - "Range restriction (in bp)" to `Correct matrix plot`
>    - "Normalize each chromosome separately" to `True`
>    - "Remove bins of low coverage" to `-2`
>    - "Remove bins of large coverage" to `3`
>    - "Include chromosomes" to `chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY`
>
{: .hands_on}

It can happen that the correction stops with:

```
ERROR:iterative correction:*Error* matrix correction produced extremely large values.
This is often caused by bins of low counts. Use a more stringent filtering of bins.
```

This can be solved by a more stringent z-score values for the filter threshold or by a look at the plotted matrix. In our case we see that chromosome Y is having more or less 0 counts in its bins. This chromosome can be excluded from the correction by not defining it for the set of chromosomes that should be corrected (parameter 'Include chromosomes').

### Plotting the corrected Hi-C matrix

We can now plot the one of the chromosomes (e.g. chromosome X) , with the corrected matrix.

> ### {% icon hands_on %} Hands-on: Plotting the corrected Hi-C matrix
>
> 1. **hicPlotMatrix** {% icon tool %}: Run hicPlotMatrix on the corrected matrix adjusting the parameters:
>    - "Plot title" to `Hi-C matrix for mESC`
>    - "Plot per chromosome" to `False`
>    - "Chromosomes to include (and order to plot in)" to `chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY`
>    - "Plot only this region" to `chrX` 'or' `X`
>    - "Plot the log1p of the matrix values" to `True`
>
{: .hands_on}

# TAD calling

“The partitioning of chromosomes into topologically associating domains (TADs) is an emerging concept that is reshaping our understanding of gene regulation in the context of physical organization of the genome” [Ramirez et al. 2017](https://doi.org/10.1101/115063).

TAD calling works in two steps: First HiCExplorer computes a TAD-separation score based on a z-score matrix for all bins. Then those bins having a local minimum of the TAD-separation score are evaluated with respect to the surrounding bins to decide assign a p-value. Then a cutoff is applied to select the bins more likely to be TAD boundaries.

[hicFindTADs](http://hicexplorer.readthedocs.io/en/latest/content/tools/hicFindTADs.html#hicfindtads) tries to identify sensible parameters but those can be change to identify more stringent set of boundaries.

> ### {% icon hands_on %} Hands-on: Finding TADs
>
> 1. **hicFindTADs** {% icon tool %}: Run hicFindTADs on the corrected matrix adjusting the parameters:
>    - "Minimum window length (in bp) to be considered to the left and to the right of each Hi-C bin." to `60000`
>    - "Maximum window length (in bp) to be considered to the left and to the right of each Hi-C bin." to `120000`
>    - "Step size when moving from minDepth to maxDepth" to `20000`
>    - "Multiple Testing Corrections" to `False discovery rate`
>    - "q-value" to `0.05` 
>    - "Minimum distance between boundaries (in bp)." to `80000`
>
{: .hands_on}

As an output we get the boundaries, domains and scores separated files. We will use in the plot later only the TAD-score file.

We can plot the TADs for a given chromosomal region. For this we need to create a track file containing the instructions to build the plot. The [hicPlotTADs](http://hicexplorer.readthedocs.io/en/latest/content/tools/hicPlotTADs.html) documentation contains the instructions to build the track file.

In following plot we will use the listed track file. Please store it as track.ini with Upload Manager.

```[hic]
file = hicMatrix/replicateMerged.Corrected_20kb.h5
title = HiC mESC chrX:99974316-101359967
colormap = RdYlBu_r
depth = 2000000
width = 7
transform = log1p
x labels = yes
type = interaction
file_type = hic_matrix
boundaries_file =  TADs/marks_et-al_TADs_20kb-Bins_domains.bed

[x-axis]
fontsize=16
where=top

[tad score]
file = TADs/marks_et-al_TADs_20kb-Bins_score.bedgraph
title = "TAD separation score"
width = 2
type = lines
color = blue
file_type = bedgraph

[spacer]

[gene track]
file = mm10_genes_sorted.bed
width = 10
title = "mm10 genes"
width = 5
labels = off
```

> ### {% icon hands_on %} Hands-on: Plotting TADs
>
> 1. **hicPlotTADs** {% icon tool %}: Run hicPlotTADs adjusting the parameters:
>    - "Region of the genome to limit the operation" to `chrX:98000000-105000000`
>
{: .hands_on}

# Conclusion
{:.no_toc}

In this tutorial we used HiCExplorer to analyze a published dataset. 
