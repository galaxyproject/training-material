---
layout: tutorial_hands_on
topic_name: sRNA-Seq
tutorial_name: small RNA sequencing analysis
---

# Introduction

**Give a short introduction to small RNAs and why we care about them.** The data used in this Galaxy tutorial are from polyphosphatase-treated small RNA sequencing experiments in *Drosophila*. The goal of this study was to determine how piRNA and piRNA target expression changes in flies mutant for Kinesin-like protein at 10A (*klp10A*). To that end, mRNA-seq experiments were performed in parallel to determine whether targets of differentially expressed piRNAs were also differentially expressed and can be analyzed by following the *de novo* transcriptome reconstruction tutorial **(add link)**. Because of the long processing time for the large original files - which contained 7-22 million reads - we have downsampled the fastq reads to include only those that align to **something interesting (flamenco, gypsy (retro-elements), diminutive (*dm/Myc*), ZAM element (LTR-retrotransposon), cluster 20A, maybe just all of the X-chromosome)**.

# Analysis strategy

The goal of this exercise is to identify what small RNAs, specifically piRNAs, are present in wild-type (WT) flies and flies treated with *klp10A* RNAi (*klp10A* KD). In this study, biological triplicate small RNA- and mRNA-seq samples for both WT and *klp10A* KD flies. We will quantify piRNA expressed from aligned reads as well as identify differentially expressed piRNAs. We will generally follow a popular piRNA analysis pipeline developed by the Zamore Lab and ZLab at UMass Med School called [PiPipes](https://github.com/bowhan/piPipes). Although PiPipes was developed for analysis of piRNAs, many of the basical principles can be applied to other classes of small RNAs. It is of note that this tutorial uses datasets that have been de-multiplexed so that the input data files are a single FASTQ formatted file for each sample. This tutorial also uses datasets for which the quality scores are encoded using the Sanger/Illumina 1.9 encoding scheme (**Check with FASTQC tool, and if not, use the FASTQ Groomer tool to convert FASTQ files to Sanger/Illumina 1.9+ encoding)**. Because small RNAs are, well, small, single-end sequencing is almost always used for sRNA-seq libraries. This tutorial uses the *Collections* feature of Galaxy to orgainze each set of replicates into a single group, making tool form submission easier.

> ### Agenda
>
> In this tutorial, we will address:
>
> 1. Data upload and organization
> 1. Read quality checking and trimming
> 1. Read alignment
> 1. Small RNA annotation
> 1. Small RNA abundance estimation
> 1. Small RNA differential expression testing
> 1. Small RNA and mRNA integration
> 1. Visualization

## Data upload and organization

Due to the large size of the original sRNA-seq datasets, we have downsampled them to only inlcude reads mapping to **something interesting, probably chromosome X**. These datasets are avaialble at [`Zenodo`](https://zenodo.org/record/####), where you can find the FASTQ files corresponding to replicate sRNA-seq and mRNA-seq libraries and an annotation file of known RefSeq transcripts for the *Drosophila melanogaster* genome version dm3.

> ### :pencil2: Hands-on: Data upload and organization
>
> 1. Create a new history and name it something meaningful (*e.g.* sRNA-seq tutorial)
> 1. Open the data upload manager by selecting *Get Data* from the Tool Panel and clicking *Upload File*
> 1. Select *Paste/Fetch Data*
> 1. Copy and paste each link for the 6 read (.fq) and 1 annotation (.gtf) files into a separate text field
>    - Set the datatype of the read (.fq) files to **fastqsanger**
>    - Set the datatype of the annotation (.gtf) file to **gtf** and assign the Genome as **dm3**
> 1. Click *Start*
> 1. Rename the files in your history to something meaningful (*e.g.* control_sRNA_rep1.fq)
> 1. Build a *Dataset list* for each set of replicates
>    - Click the *Operations on multiple datasets* check box at the top of the history panel
>    - Check the three boxes next to the control RNAi (control) sRNA-seq samples
>    - Click *For all selected...* and choose *Build dataset list*
>    - Ensure the three control samples are the only ones selected, and enter a name for the new collection (*e.g.* control sRNA-seq)
>    - Click *Create list*
>    - Repeat steps ii. through v. for the *klp10A* RNAi samples.
>
> {: .hands_on}

## Read quality checking and trimming

Small RNA sequencing library preparations involve adding an artificial adaptor sequence to both the 5' and 3' ends of the small RNAs. While the 5' adaptor anchors reads to the sequencing surface and thus are not sequenced, the 3' adaptor is typically sequence immediately follow the RNA sequence. These artificial sequences need to be removed before attemping to align the reads to a reference. We know that the datasets used here contain the Illumina universal small RNA-seq 3' adaptor, but let's confirm this using the [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) tool (further described in the [NGS-QC tutorial](../../NGS-QC/tutorials/dive_into_qc)).

> ### :pencil2: Hands-on: Quality control
>
> 1. **FastQC** :wrench:: Run `FastQC` on the FASTQ read files to identify adaptors and assess the quality of the reads. Select the `FASTQC` tool. Under the **Short read data from your current history** option select the *Dataset collections* tab and then choose one of the dataset collections to analyze. Repeat for the second dataset collection.
>
>    > ### :question: Questions
>    >
>    > 1. What is the read length for each sample?
>    > 1. What does the base/read quality look like?
>    > 1. Are there any adaptors present in these reads? Which one(s)?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>All samples have a read length of 51 nt. </li>
>    >    <li>The read quality across the entire length of the reads is good (phred score > 28 for the most part). </li>
>    >    <li>Yes, "Illumina Small RNA 3' Adapters" are present. </li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 1. **Trim Galore!** :wrench:: Trim off Illumina adaptors from the 3' ends of the reads by running `Trim Galore!` on every FASTQ file with the following parameters:
>    - **Is this library paired- or single-end?**: Single-end
>    - **Reads in FASTQ format**: Select "Dataset collections" and then select the control sRNA-seq dataset
>    - **Trimming reads?**: Illumina small RNA adapters
>    - **Trim Galore! advanced settings**: Full parameter list
>    - **Overlap with adapter sequence required to trim a sequence**: 6
>    - **Discard reads that became shorter than length INT**: 12
>
>    ![](../images/image.png)
>
> 1. **FastQC** :wrench:: Re-run `FastQC` on trimmed reads and inspect the differences.
>
>    > ### :question: Questions
>    >
>    > 1. What is the read length?
>    > 1. Are there any adaptors present in these reads? Which one(s)?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>The read lengths range from 12 to 51 nt after trimming.</li>
>    >    <li>No, Illumina Small RNA 3' Adaptors are no longer present. No other adaptoers are present.</li>
>    >    </ol>
>    >    </details>
>    {: .question}
> ![](../images/image.png)
> {: .hands_on}

Now that we have trimmed our reads of the Illumina Small RNA 3' adaptors, we will align our trimmed reads to the reference *Drosophila* genome (dm3). For miRNA analyses, it is useful to align to the reference set of known miRNAs first, and then re-align any unaligned reads to the reference genome.

## Read alignment

To quantify small RNA abundance and identify their putative targets, we need to know where the sequenced reads align to a reference genome. In the case of a eukaryotes, some small RNAs are transcribed from mRNA templates, which means that some small RNAs can originate from an exon-exon (spliced) boundary. Therefore, a splice-aware aligner must be used to account for this possibility. [`HISAT2`](https://ccb.jhu.edu/software/hisat2/index.shtml) is an accurate and fast tool for aligning spliced reads to a genome, and we will be using `HISAT2` in this tutorial.

> ### :pencil2: Hands-on: Splice-aware alignment
>
> 1. **HISAT2** :wrench:: Run `HISAT2` on one collection of trimmed reads with the following parameters:
>    - **Single end or paired reads?**: Individual unpaired reads
>    - **Reads**: Select "Dataset collection" and choose one collection of trimmed FASTQ files
>    - **Source for the reference genome to align against**: Use a built-in genome
>    - **Select a reference genome**: Fruit fly (Drosophila melanogaster): dm3
>    - **Spliced alignment parameters**: Specify spliced alignment parameters
>    - **Specify strand-specific information**: First Strand (R/RF)
>
>       ![](../images/image.png)
>
> 1. **HISAT2** :wrench:: Run `HISAT` on the second collection of trimmed reads with the same parameters.
>

**UPDATES STOPPED HERE**

## Small RNA annotation
**TODO Describe what is meant by small RNA annotation. Describe sense v. antisense meaning. Describe that different classes of small RNAs are in a small RNA-seq library. Talk about size distribution and nt composition biases in fly piRNAs and in other classes of small RNAs.**

> ### :pencil2: Hands-on: Small RNA annotation
>
> 1. **Tool** :wrench:: Run `Tool` on one collection of `HISAT` alignments using the default parameters.
>    - Use batch mode to run all four samples from one tool form.
>
> ![](../images/image.png)

## Small RNA abundance estimation

**TODO RPM for small RNA counts.**

> ### :pencil2: Hands-on: Small RNA abundance estimation
>
> 1. **Tool** :wrench:: Run `Tool` on the `Tool`-annotated small RNAs.
>    - Use batch mode to inlcude all four `Stringtie` assemblies.
>    - **Use Reference Annotation**: Yes, then select the "RefSeq GTF mm10" file.
> ![](../images/image.png)
>

## Small RNA differential expression testing

### Analysis of the differential gene expression

We just generated a transriptome database that represents the transcripts present in the G1E and megakaryocytes samples. This database provides the location of our transcripts with non-redundant identifiers, as well as information regarding the origin of the transcript.

We now want to identify which transcripts are differentially expressed between the G1E and megakaryocyte cellular states. To do this we will implement a counting approach using `FeatureCounts` to count reads per transcript. Then we will provide this information to `DESeq2` to generate normalized transcript counts (abundance estimates) and significance testing for differential expression.

### Count the number of reads per transcript

To compare the abundance of transcripts between different cellular states, the first essential step is to quantify the number of reads per transcript. [`FeatureCounts`](http://bioinf.wehi.edu.au/featureCounts/) is one of the most popular tools for counting reads in genomic features. In our case, we'll be using `FeatureCounts` to count reads aligning in exons of our `Cuffmerge` generated transcriptome database.

The recommended mode is "union", which counts overlaps even if a read only shares parts of its sequence with a genomic feature and disregards reads that overlap more than one feature.

> ### :pencil2: Hands-on: Counting the number of reads per transcript
>
> 1. **FeatureCounts** :wrench:: Run `FeatureCounts` on the aligned reads (`HISAT` output) using the `Cuffmerge` transcriptome database as the annotation file.
>
>    - Using the batch mode for input selection, choose the four `HISAT` aligned read files
>    - **Gene annotation file**:  in your history, then select the GTF file output by Cuffmerge (this specifies the "union" mode)
>    - Expand **Options for paired end reads**
>    - **Orientation of the two read from the same pair**: Forward, Reverse (fr)
>    - Expand **Advanced options**
>    - **GFF gene identifier**: enter "transcript_id"
>    - **Strand specificity of the protocol**: select "Stranded (forwards)"
> ![](../images/image.png)
>
> {: .hands_on}

### Perform differential gene expression testing

Transcript expression is estimated from read counts, and attempts are made to correct for variability in measurements using replicates. This is absolutely essential to obtaining accurate results. We recommend having at least two biological replicates.

[`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) is a great tool for differential gene expression analysis. It takes read counts produced by `FeatureCounts` and applies size factor normalization:

- Computation for each gene of the geometric mean of read counts across all samples
- Division of every gene count by the geometric mean
- Use of the median of these ratios as sample's size factor for normalization

> ### :pencil2: Hands-on:
>
> 1. **DESeq2** :wrench:: Run `DESeq2` with the following parameters:
>    - Specify "G1E" as the first factor level (condition) and select the count files corresponding to the two replicates
>    - Specify "Mega" as the second factor level (condition) and select the count files corresponding to the two replicates
>
>       > ### :nut_and_bolt: Comment
>       >
>       > You can select several files by holding down the CTRL (or COMMAND) key and clicking on the desired files
>       {: .comment}
>    - **Visualising the analysis results**: Yes
>    - **Output normalized counts table**: Yes
>
{: .hands_on}

The first output of `DESeq2` is a tabular file. The columns are:

1.	Gene identifiers
2.	Mean normalized counts, averaged over all samples from both conditions
3.	Logarithm (base 2) of the fold change (the values correspond to up- or downregulation relative to the condition listed as Factor level 1)
4.	Standard error estimate for the log2 fold change estimate
5.	[Wald](https://en.wikipedia.org/wiki/Wald_test) statistic
6.	*p*-value for the statistical significance of this change
7.	*p*-value adjusted for multiple testing with the Benjamini-Hochberg procedure which controls false discovery rate ([FDR](https://en.wikipedia.org/wiki/False_discovery_rate))


> ### :pencil2: Hands-on:
>
>1. **Filter** :wrench:: Run `Filter` to extract genes with a significant change in gene expression (adjusted *p*-value less than 0.05) between treated and untreated samples
>
>    > ### :question: Question
>    >
>    > How many transcripts have a significant change in expression between these conditions?
>    >
>    > <details>
>    > <summary>Click to view answers</summary>
>    > To filter, use "c7 lessthan 0.05". And we get 278 transcripts with a significant change in gene expression between the G1E and megakaryocyte cellular states.
>    > </details>
>    {: .question}
>
> 2. **Filter** :wrench:: Determine how many transcripts are up or down regulated in the G1E state.
>
>    > ### :nut_and_bolt: Comments
>    > Rename your datasets for the downstream analyses
>    {: .comment}
>
>    > ### :question: Question
>    >
>    > Are there more upregulated or downregulated genes in the treated samples?
>    >
>    > <details>
>    > <summary>Click to view answers</summary>
>    > To obtain the up-regulated genes in the G1E state, we filter the previously generated file (with the significant change in transcript expression) with the expression "c3>0" (the log2 fold changes must be greater than 0). We obtain 131  genes (47.1% of the genes with a significant change in gene expression). For the down-regulated genes in the G1E state, we did the inverse and we find 147 transcripts (52.9% of the genes with a significant change in transcript expression)
>    > </details>
>    {: .question}
{: .hands_on}

In addition to the list of genes, `DESeq2` outputs a graphical summary of the results, useful to evaluate the quality of the experiment:

1. Histogram of *p*-values for all tests

    ![](../images/image.png)

2. [MA plot](https://en.wikipedia.org/wiki/MA_plot): global view of the relationship between the expression change of conditions (log ratios, M), the average expression strength of the genes (average mean, A), and the ability of the algorithm to detect differential gene expression. The genes that passed the significance threshold (adjusted p-value < 0.1) are colored in red.

    ![](../images/image.png)

3. Principal Component Analysis ([PCA](https://en.wikipedia.org/wiki/Principal_component_analysis)) and the first two axes

    ![](../images/image.png)

    Each replicate is plotted as an individual data point. This type of plot is useful for visualizing the overall effect of experimental covariates and batch effects.

4. Heatmap of sample-to-sample distance matrix: overview over similarities and dissimilarities between samples

    ![](../images/image.png)

5. Dispersion estimates: gene-wise estimates (black), the fitted values (red), and the final maximum a posteriori estimates used in testing (blue)

    ![](../images/image.png)

    This dispersion plot is typical, with the final estimates shrunk from the gene-wise estimates towards the fitted estimates. Some gene-wise estimates are flagged as outliers and not shrunk towards the fitted value. The amount of shrinkage can be more or less than seen here, depending on the sample size, the number of coefficients, the row mean and the variability of the gene-wise estimates.


For more information about `DESeq2` and its outputs, you can have a look at [`DESeq2` documentation](https://www.bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf).

## Small RNA and mRNA integration

**TODO**

## Visualization
Now that we have a list of transcript expression levels and their differential expression levels, it is time to visually inspect our transcript structures and the reads they were predicted from. It is a good practice to visually inspect (and present) loci with transcripts of interest. Fortuantely, there is a built-in genome browser in Galaxy, **Trackster**, that make this task simple (and even fun!).

In this last section, we will convert our aligned read data from BAM format to bigWig format to simplify observing where our stranded RNA-seq data aligned to. We will then initiate a session on Trackster, load it with our data, and visually inspect our interesting loci.

> ### :pencil2: Hands-on: Converting aligned read files to bigWig format
>
> 1. **bamCoverage** :wrench:: Run `bamCoverage` on all four aligned read files (`HISAT` output) with the following parameters:
>    - **Bin size in bases**: 1
>    - **Effective genome size**: mm9 (2150570000)
>    - Expand the **Advanced options**
>    - **Only include reads originating from fragments from the forward or reverse strand**: forward
> 2. **Rename** :wrench:: Rename the outputs to reflect the origin of the reads and that they represent the reads mapping to the PLUS strand
>![](../images/image.png)
>
> 3. **bamCoverage** :wrench:: Repeat Step 1 except changing the following parameter:
>    - **Only include reads originating from fragments from the forward or reverse strand**: reverse
> 4. **Rename** :wrench:: Rename the outputs to reflect the origin of the reads and that they represent the reads mapping to the MINUS strand
> ![](../images/image.png)

> ### :pencil2: Hands-on: Trackster based visualization
>
> 1. **Viz** :wrench:: On the center console at the top of the Galaxy interface, choose " Visualization" -> "New track browser"
>    - Name your visualization someting descriptive under "Browser name:"
>    - Choose "Mouse Dec. 2011 (GRCm38/mm10) (mm10)" as the "Reference genome build (dbkey)
>    - Click "Create" to initiate your Trackster session
> ![](../images/image.png)
>
> 2. **Viz** :wrench:: Click "Add datasets to visualization"
>    - Select the "RefSeq GTF mm10" file
>    - Select the output files from `Stringtie`
>    - Select the output file from `Cuffmerge`
>    - Select the output files from `bamCoverage`
>
> 3. :wrench:: Using the grey labels on the left side of each track, drag and arrange the track order to your preference
>
> 4. :wrench:: Hover over the grey label on the left side of the "RefSeq GTF mm10" track and click the "Edit settings" icon.
>    - Adjust the block color to blue (#0000ff) and antisense strand color to red (#ff0000)
>
> 5. :wrench:: Repeat the previous step on the output files from `StringTie` and `Cuffmerge`
>
> 6. :wrench:: Hover over the grey label on the left side of the "G1E R1 plus" track and click the "Edit settings" icon.
>    - Adjust the color to blue (#0000ff)
>
> 7. :wrench:: Repeat the previous step on the other three bigWig files representing the plus strand
>
> 8. :wrench:: Hover over the grey label on the left side of the "G1E R1 minus" track and click the "Edit settings" icon.
>    - Adjust the color to red (#ff0000)
>
> 9. :wrench:: Repeat the previous step on the other three bigWig files representing the minus strand
>
> 10. :wrench:: Adjust the track height of the bigWig files to be consistant for each set of plus strand and minus strand tracks
> ![](../images/image.png)
> 11. :wrench:: Direct Trackster to the coordinates: chr11:96193539-96206376, what do you see?
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>There are two clusters of transcripts that are exclusively expressed in the G1E background</li>
>    >    <li>The left-most transcript is the Hoxb13 transcript</li>
>    >    <li>The center cluster of transcripts are not present in the RefSeq annotation and are determined by `Cuffmerge` to be "u" and "x"</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>

## Conclusion

In this tutorial, we have analyzed real RNA sequencing data to extract useful information, such as which genes are up- or down-regulated by depletion of the Pasilla gene and which genes are regulated by the Pasilla gene. To answer these questions, we analyzed RNA sequence datasets using a reference-based RNA-seq data analysis approach. This approach can be sum up with the following scheme:


![](../images/schematic_for_sRNAseq_tutorial.png)
