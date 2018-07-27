---
layout: tutorial_hands_on
topic_name: chip-seq
tutorial_name: formation_of_super-structures_on_xi
---

# Introduction
{:.no_toc}

This exercise uses the published data from [Chen-Yu Wang et al., 2018](https://www.cell.com/cell/fulltext/S0092-8674(18)30584-1).
The goal of this research was to investigate the mechanism by which SMCHD1 gene shapes the Xi and represses gene expression.

To this end, several ChIP-seq experiments were performed on both wild-type and SMCHD1 gene knockdown cells including ChIP data for H3K27me3, H3K4me3 histone marks and CTCF transcription factor which will be used in the upcoming tutorial to present a step by step ChIP-seq data analysis.

During the following steps, the corresponding 'input' samples, for which the same treatment as the ChIP-seq samples was done except for the immunoprecipitation step, are also used along with the 'ChIP-seq' samples to identify the potential sequencing bias and help for differential analysis.

Because of the long processing time for the large original files, we have extracted the data from the chromosome X and provide you with the already processed data for the subsequent steps.

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Step 1: Quality control and treatment of the sequences

The first step of any ChIP-Seq data analysis is quality control of the raw sequencing data.

> ### {% icon hands_on %} Hands-on: Quality control
>
> 1. Create a new history for this tutorial and give it a proper name
> 2. Import `wt_H3K4me3_read1.fastq` (link in [Zenodo](https://zenodo.org/record/1321974/files/wt_H3K4me3_read1.fastq)) and `wt_H3K4me3_read2.fastq` (link in [Zenodo](https://zenodo.org/record/1321974/files/wt_H3K4me3_read2.fastq)) from
 [Zenodo](https://zenodo.org/) or from the data library into the history
>
>    > ### {% icon tip %} Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * Press **Start** and **Close** the window
>    > * Click on the `pencil` icon once the file is imported
>    > * If format is not `fastqsanger`:
>    > > * Click on **Datatype** in the middle panel
>    > > * Select `fastqsanger` as **New Type**
>    {: .tip}
>
>    > ### {% icon tip %} Tip: Importing data from a data library
>    >
>    > * Go into "Shared data" (top panel) then "Data libraries"
>    > * Click on "Training data" and then "Analyses of ChIP-Seq data"
>    > * Select interesting files
>    > * Click on "Import selected datasets into history"
>    > * Import in the history
>    {: .tip}
>
>    As default, Galaxy takes the link as name, so rename them.
>
>
> 3. Inspect the file by clicking on the `eye` icon
>
>    > ### {% icon question %} Questions
>    >
>    > 1. How are the DNA sequences stored?
>    > 2. What are the other entries?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. The DNA sequences are stored in the second line of every 4-line group
>    > > 2. This file is called a [FastQ file](https://en.wikipedia.org/wiki/FASTQ_format). It stores sequence information and quality information. Each sequence is represented by a group of 4 lines with the 1st line being the sequence id, the second the sequence of nucleotides, the third a transition line and the last one a sequence of quality score for each nucleotide.
>    > {: .solution }
>    {: .question}
>
> 4. Run **FastQC** {% icon tool %} with
>    - "Short read data from your current history" to the imported files
>    -  "Execute"
>
>    Inspect the generated files
{: .hands_on}

It is often necessary to trim sequenced read, for example, to get rid of bases that were sequenced with high uncertainty (= low quality bases).

> ### {% icon hands_on %} Hands-on: Quality control
>
> 1. Run **Trim Galore!** {% icon tool %} with
>    - "Is this library paired- or single-end?" to `Paired-end`
>    - "Reads in FASTQ format" to the imported files
>    - "Trim Galore! advanced settings" to `Full parameter list`
>    - "Trim low-quality ends from reads" to `15`
>    - "Overlap with adapter sequence required to trim a sequence" to `3`
>    - "Execute"
>
>    > ### {% icon tip %} Tip: Importing data from a data library
>    >
>    > If your FASTQ files cannot be selected, you might check whether their format is FASTQ with Sanger-scaled quality values (`fastqsanger`). You can edit the data type by clicking on the `pencil` symbol.
>    {: .tip}
>
{: .hands_on}

# Step 2: Mapping of the reads

In order to figure where the sequenced DNA fragments originated from in the genome, the short reads must be aligned to the reference genome. This is equivalent to solving a jigsaw puzzles, but unfortunately, not all pieces are unique. In principle, you could do a BLAST analysis to figure out where the sequenced pieces fit best in the known genome. Aligning millions of short sequences this way may, however, take a couple of weeks.

### Running Bowtie2

Nowadays, there are many read alignment programs for shotgun sequenced DNA, Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is one of them.

> ### {% icon hands_on %} Hands-on: Mapping
>
> 1. **Bowtie2** {% icon tool %} with
>    - "Is this single or paired library" to `Paired-end`
>    - "FASTA/Q file" to the Trim Galore! output with the trimmed reads
>    - "Will you select a reference genome from your history or use a built-in index?" to `Use a built-in genome index`
>    - "Select reference genome" to `Mouse (Mus musculus): mm10`
>    - "Execute"
>
>    > ### {% icon question %} Questions
>    >
>    > 1. How many reads where mapped?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. This information can be accessed by clicking on the resulting history entry. You can see some basic mapping statistics once the alignment is completed. 100220 (81.55%) aligned concordantly exactly 1 time
>    22350 (18.19%) aligned concordantly >1 times.
>    > {: .solution }
>    {: .question}
>
{: .hands_on}

The read alignment step with bowtie2 resulted in a compressed, binary file (BAM) that is not human-readable. It's like the zipped version of a text file.

We will show you two ways to inspect the file:

1. Visualization using a Genome Browser
2. Converting the binary format into its text file equivalent


## Visualization using a Genome Browser

> ### {% icon hands_on %} Hands-on: Visualization of the reads in IGV
>
> 1. Click on the `display with IGV local` to load the reads into the IGV browser
> 2. Zoom at the start of chromosome X (or `chrX:106,435,807-106,439,377`)
{: .hands_on}

The reads have a direction: they are mapped to the forward or reverse strand, respectively. When hovering over a read, extra information is displayed

> ### {% icon question %} Questions
>
> 1. Some reads have colored lines included. What is this?
>
> > ### {% icon solution %} Solution
> > 1. Try to zoom in in one of those lines and you will see the answer!
> {: .solution }
{: .question}

> ### {% icon comment %} Comments
> Because the number of reads over a region can be quite large, the IGV browser by default only allows to see the reads that fall into a small window. This behaviour can be changed in the IGV from `view > Preferences > Alignments`.
{: .comment}


## Inspection of the SAM format

As mentioned above, you can convert the binary BAM file into a simple (but large!) text file, which is called a [SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) (Sequence Alignment Map) file.

> ### {% icon hands_on %} Hands-on: Conversion into a SAM file
>
> 1. **BAM-to-SAM** {% icon tool %} with
>    - "BAM File to Convert" to the file generated by Bowtie2
>    - "Header options" to `Include header in SAM output`
>    - "Execute"
>
> 2. Inspect the file by clicking on `eye` icon
>
{: .hands_on}

A SAM file is a file with

- A header with the chromosome names and lengths
- A file content as a tabular file with the location and other information of each read found in the FASTQ file and the mapping information

> ### {% icon question %} Questions
>
> 1. Which information do you find in a SAM/BAM file? What is the additional information compared to a FASTQ file.
>
> > ### {% icon solution %} Solution
> > 1. Sequences and Quality information, like FASTQ
> > 2. Mapping information; Location of the read on the chromosome; Mapping quality ...
> > 3. For more information please see [https://en.wikipedia.org/wiki/SAM_(file_format)](https://en.wikipedia.org/wiki/SAM_(file_format))
> {: .solution }
{: .question}

# Step 3: ChIP-seq Quality Control

We already checked the quality of the raw sequencing reads in the first step.
Now we would like to test the quality of the ChIP-seq preparation, to know if your ChIP-seq samples are more enriched than the control (input) samples.

## Correlation between samples

To assess the similarity between the replicates of the ChIP and the input, respectively, it is a common technique to calculate the correlation of
read counts on different regions for all different samples.
We expect that the replicates of the ChIP-seq experiments should be clustered more closely to each other than the replicates of the input sample.
That is, because the input samples should not have enriched regions included - remember the immuno-precipitation step was skiped during the sample preparation.

To compute the correlation between the samples we are going to to use the QC modules of deepTools (http://deeptools.readthedocs.io/), a software package for the QC, processing and analysis of NGS data. Before computing the correlation a time consuming step is required, which is to compute the read coverage over a large number of regions from each of the inputed BAM files. For this we will use the tool **multiBamSummary** {% icon tool %}.

Since in this tutorial we are interested in assessing H3K4me3, H3K27me3 and CTCF ChIP samples, at first we need to catch up with all the replicates of ChIP samples as well as the input samples and re-run the previous steps (quality control and mapping) on each sample.
To save time, we already did that and you can now work directly on the BAM files of the 8 samples. For simplicity, the files include only ChrX.

> ### {% icon hands_on %} Hands-on: Correlation between samples
>
> 1. Create a new history
> 2. Import the 8 BAM files from [Zenodo](https://zenodo.org/record/1321974) or from the data library into the history
>    - [`wt_CTCF_rep1.bam`](https://zenodo.org/record/1321974/files/wt_CTCF_rep1.bam)
>    - [`wt_CTCF_rep2.bam`](https://zenodo.org/record/1321974/files/wt_CTCF_rep2.bam)
>    - [`wt_H3K4me3_rep1.bam`](https://zenodo.org/record/1321974/files/wt_H3K4me3_rep1.bam)
>    - [`wt_H3K4me3_rep2.bam`](https://zenodo.org/record/1321974/files/wt_H3K4me3_rep2.bam)
>    - [`wt_H3K27me3_rep1.bam`](https://zenodo.org/record/1321974/files/wt_H3K27me3_rep1.bam)
>    - [`wt_H3K27me3_rep2.bam`](https://zenodo.org/record/1321974/files/wt_H3K27me3_rep2.bam)
>    - [`wt_input_rep1.bam`](https://zenodo.org/record/1321974/files/wt_input_rep1.bam)
>    - [`wt_input_rep2.bam`](https://zenodo.org/record/1321974/files/wt_input_rep2.bam)
>
> 3. **multiBamSummary** {% icon tool %} with
>    - "Sample order matters" to `No`
>    - "Bam file" to the 8 imported BAM files
>    - "Choose computation mode" to `Bins`
>    - "Bin size in bp" to `1000`
>       
>       This corresponds to the length of the fragments that were sequenced; it is not the read length!
>
>    - "Distance between bins" to `500` (to reduce the computation time for the tutorial)
>    - "Region of the genome to limit the operation to" to `chrX`
>    - "Execute"
>  
>    Using these parameters, the tool will take bins of 1000 bp separated by 500. For each bin the overlapping reads in each sample will be computed
>
> 4. **plotCorrelation** {% icon tool %} with
>    - "Matrix file from the multiBamSummary tool" to the generated multiBamSummary output
>    - Correlation method `Pearson`
>    
>    To compute and visualize the sample correlation we use plotCorrelation from deepTools. This is a fast process that allows the user to quickly try different color combinations and outputs. Feel free to try different parameters.
>    
{: .hands_on}

> ### {% icon question %} Questions
>
> ![Output for plotCorrelation with the correlation scores between the 8 samples](../../images/formation_of_super-structures_on_xi/plotCorrelation.png "Correlation scores between the 8 samples")
>
> 1. How are your samples clustered? Does that correspond to your expectations?
>
> > ### {% icon solution %} Solution
> >  As one could expect, the input replicates cluster together and the ChIP replicates cluster together.
> {: .solution }
{: .question}

> ### {% icon comment %} Comments
> More information on these two tools can be found at the [deepTools documentation page](https://deeptools.readthedocs.io/en/latest/content/list_of_tools.html).
{: .comment}

## IP strength estimation

To evaluate the quality of the immuno-precipitation step, we can compute the IP strength. It determines how well the signal in the ChIP-seq sample can be differentiated from the background distribution of reads in the control sample. To do that we take the data from the `rep1` of the `wt_H3K4me3` ChIP-seq sample and compare it with its corresponding input sample.

> ### {% icon comment %} Comments
> For more information on how the IP strength is estimated, you can check the [deepTools documentation page](https://deeptools.readthedocs.io/en/latest/content/list_of_tools.html).
{: .comment}

> ### {% icon hands_on %} Hands-on: IP strength estimation
>
> 1. **plotFingerprint** {% icon tool %} with
>    - "Sample order matters" to `No`
>    - "Bam file" to `wt_input_rep1`and `wt_H3K4me3_rep1`
>    - "Region of the genome to limit the operation to" to `ChrX`
>    - "Execute"
{: .hands_on}

The plotFingerprint tool generates a fingerprint plot. You need to intepret it to know the IP strength. The [deepTools documentation](https://deeptools.readthedocs.io/en/latest/content/list_of_tools.html) explains it clearly:

![A guide to interpret a fingerprint plot](../../images/estrogen-receptor-binding-site-identification/QC_fingerprint.png "How to interpret a fingerprint plot? Image extracted from the deepTools documentation")

> ### {% icon question %} Questions
>
> ![Output for plotFingerprint with the fingerprint plot to estimate the IP strength](../../images/formation_of_super-structures_on_xi/plotFingerPrint.png "Fingerprint plot for the first replicates to estimate the IP strength")
>
> 1. What do you think about the quality of the IP for this experiment?
>
> > ### {% icon solution %} Solution
> >  There is a obvious signal between H3K4me3 and input.
> >
> >  Almost 20% of chromosome 1 are not sequenced at all.
> {: .solution }
{: .question}

> ### {% icon hands_on %} (Optional) Hands-on: IP strength estimation (other samples)
>
> 1. Run the same analysis on data of the 3 other patients
{: .hands_on}


# Step 4: Normalization
One of the goals in ChIP-seq data analysis is finding regions on the genome which are enriched for the ChIP data of interest (regions with significantly higher read coverage for the ChIP data comparing to its corresponding input). However, to reach a reliable comparison the data needs to be normalized to remove any technical bias.
In the following exercise we would like to know where the H3K4me3 binding sites are. For this we need
to extract which parts of the genome have been enriched (more reads mapped) within the samples that underwent immunoprecipitation.

For the normalization we have two options.

1. Normalization by sequencing depth
2. Normalization by input file


## Generation of coverage files normalized by sequencing depth

We first need to make the samples comparable. Indeed, the different samples have usually a different sequencing depth, i.e. a different number of reads.
These differences can bias the interpretation of the number of reads mapped to a specific genome region.

> ### {% icon hands_on %} Hands-on: Coverage file normalization
>
> 1. **IdxStats** {% icon tool %} with
>    - "BAM file" to "Multiple datasets": `wt_H3K4me3_rep1.bam` and `wt_input_rep1.bam`
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What is the output of this tool?
>    > 2. How many reads has been mapped on chrX for the input and for the ChIP-seq samples?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. This tool estimates how many reads mapped to which chromosome. Furthermore, it tells the chromosome lengths and naming convention (with or without 'chr' in the beginning)
>    > > 2. 1,204,821 for ChIP-seq samples and 1,893,595 for the input
>    > {: .solution }
>    {: .question}
>
> 2. **deeptools-bamCoverage** {% icon tool %} with
>    - "BAM file" to "Multiple datasets": `wt_H3K4me3_rep1.bam` and `wt_input_rep1.bam`
>    - "Bin size in bases" to `25`
>    - "Scaling/Normalization method" to `Normalize coverage to 1x`
>    - "Effective genome size" to `GRCm38/mm10 (2308125349)`
>    - "Coverage file format" to `bedgraph`
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What are the different columns of a `bedgraph` file?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. chrom, chromStart, chromEnd and a data value
>    > {: .solution }
>    {: .question}
>
> 3. **bamCoverage** {% icon tool %} with the same parameters but to generate a `bigWig` output file
> 4. **IGV** {% icon tool %} to inspect both signal coverages (input and ChIP samples) in IGV
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What is a bigWig file?
>    >
>    > > ### {% icon solution %} Solution
>    > > A bigWig file is a compressed bedgraph file. Similar in relation as BAm to SAM, but this time just for coverage data. This means bigWig and bedgraph
>    > > files are much smaller than BAM or SAM files.
>    > {: .solution }
>    {: .question}
>
{: .hands_on}

## Generation of input-normalized coverage files and their visualization

To extract only the information induced by the immunoprecipitation, we normalize for each patient the coverage file for the sample that underwent immunoprecipitation by the coverage file for the input sample. Here we use the tool bamCompare which compare 2 BAM files while caring for sequencing depth normalization.

> ### {% icon hands_on %} Hands-on: Generation of input-normalized coverage files
>
> 1. **bamCompare** {% icon tool %} with
>    - "First BAM file (e.g. treated sample)" to `wt_H3K4me3_rep1.bam`
>    - "Second BAM file (e.g. control sample)" to `wt_input_rep1.bam`
>    - "Bin size in bases" to `50`
>    - "How to compare the two files" to `Compute log2 of the number of reads ratio`
>    - "Coverage file format" to `bedgraph`
>    - "Region of the genome to limit the operation to" to `chrX` (to reduce the computation time for the tutorial)
>
>      > ### {% icon question %} Questions
>      >
>      > 1. What does mean a positive or a negative value in the 4th column?
>      >
>      > > ### {% icon solution %} Solution
>      > > 1. The 4th column contains the log2 of the number of reads ratio between the ChIP-seq sample and the input sample. A positive value means that the coverage on the portion is more important in the ChIP-seq sample than in the input sample
>      > {: .solution }
>      {: .question}
>
> 2. **bamCompare** {% icon tool %} with the same parameters but to generate a `bigWig` output file
> 3. **IGV** {% icon tool %} to inspect the log2 ratio
>
>    Remember that the bigWig file contains only the signal on chromosome X!
>
{: .hands_on}

# Step 5: Detecting enriched regions (peak calling)

We can also call the enriched regions, or peaks, found in the ChIP-seq samples.

> ### {% icon hands_on %} Hands-on: Peak calling
>
> 1. **MACS2 callpeak** {% icon tool %} with
>    - "ChIP-Seq Treatment File" to `wt_H3K4me3_rep1.bam`
>    - "ChIP-Seq Control File" to `wt_input_rep1.bam`
>    - "Format of Input Files" to `Paired-end BAM`
>    - "Effective genome size" to `M.musculus(1.87e9)`
>    - "Outputs" to `Summary page (html)`
>
>    > ### {% icon comment %} Comments
>    > The advanced options may be adjusted, depending of the samples.
>    > If your ChIP-seq experiment targets regions of broad enrichment, *e.g.* non-punctuate histone modifications, select calling of broad regions.
>    > If your sample has a low duplication rate (*e.g.* below 10%), you might keep all duplicate reads (tags). Otherwise, you might use the 'auto' option to estimate the maximal allowed number of duplicated reads per genomic location.
>    {: .comment}
>
> 2. **IGV** {% icon tool %} to inspect with the signal coverage and log2 ratio tracks
>
{: .hands_on}

The called peak regions can be filtered by, *e.g.* fold change, FDR and region length for further downstream analysis.

# Step 6: Plot the signal on the peaks between samples

Plotting your region of interest will involve using two tools from the **deepTools** suite.
+ computeMatrix : Computes the signal on given regions, using the bigwig coverage files from different samples.
+ plotHeatmap : Plots heatMap of the signals using the computeMatrix output.

Optionally, you can also use `plotProfile`to create a profile plot using to computeMatrix output.

## computeMatrix

> ### {% icon hands_on %} Hands-on: Visualization of the coverage
>
> 1. **UCSC Main** {% icon tool %} with
>    - "assembly" to `mm10`
>    - "track" to `RefSeq genes`
>    - "region" to `position` with `chrX`
>    - "output format" to `BED`
>    - "Send output to" to `Galaxy`
>    - "Get output"
>    - "Send query to Galaxy"
>
> 2. **computeMatrix** {% icon tool %} with
>    - "Regions to plot" to the imported UCSC file
>    - "Score file" to the bigwig file generated by bamCompare
>    - "computeMatrix has two main output options" to `scale-regions`
>
>       This option stretches or shrinks all regions in the BED file (here: genes) to the same length (bp) as indicated by the user
>
>    - "Show advanced options" to `yes`
>    - "Convert missing values to 0?" to `Yes`
>    - "Execute"
>
>    This tool prepares a file with scores per genomic region, which is required as input for the next tool.
>
{: .hands_on}

## plotHeatmap

> ### {% icon hands_on %} Hands-on: Visualization of the coverage
> 3. **plotHeatmap** {% icon tool %} with
>    - "Matrix file from the computeMatrix tool" to the generated matrix
>    - "Show advanced options" to `yes`
>    - "Did you compute the matrix with more than one groups of regions?" to the correct setting
>    - "Execute"
>
{: .hands_on}

# Conclusion
{:.no_toc}

![Summary of the different steps of the tutorial and the generated files](../../images/estrogen-receptor-binding-site-identification/tutorial-scheme.png "Different steps of the tutorials with the generated files")


# Additional exercise (if you have finished all above)

## Plotting heatmap from multiple samples with clustering

> ### {% icon hands_on %} Hands-on: plotting multiple samples
> The goal is to visualize more than a sample on the heatmap. To this end, we ask you to generate a matrix which includes 2 different samples (H3K4me3, CTCF). Since You have already generated the required files for the H3K4me3 sample, let's make them only for the CTCF sample.
> 1. Run **bamCompare** {% icon tool %} with  same parameters as above, for the first replicate of CTCF:
>   - "First BAM file (e.g. treated sample)" to `wt_CTCF_rep1.bam`
>   - "Second BAM file (e.g. control sample)" to `wt_input_rep1.bam`
>   - save as a `bigwig` file.
>
> 2. Perform peak calling again using treatment file : `wt_CTCF_rep1.bam` and control `wt_input_rep1.bam`, using macs2 parameters same as above.
>
> 3. Concatenate the MACS2 outputs (summits in BED) from `CTCF` and `H3K4me3` using `Operate on Genomic Intervals` --> `Concatenate`
>
> 4. Sort the output `Operate on Genomic Intervals` --> `sortBED`
>
> 5. Merge the overlapping intervals using `Operate on Genomic Intervals` --> `MergeBED`
>
> 6. **computeMatrix** {% icon tool %} with the same parameters but:
>    - Regions to plot : select the merged bed from above
>    - Score file: bigwig files from bamCompare for both H3K4me3 and CTCF samples.
>    - Output option : `reference-point`
>    - The reference point for the plotting: `center of region`
>    - Distance upstream of the start site of the regions defined in the region file : `3000`
>    - Distance downstream of the end site of the given regions: `3000`
>
>    With this option, it considers only those genomic positions before (downstream) and/or after (upstream) a reference point (*e.g.* TSS, which corresponds to the annotated gene start in our case)
>
> 7. **plotHeatmap** {% icon tool %} with
>    - "Matrix file from the computeMatrix tool" to the generated matrix
>    - "Show advanced options" to `yes`
>    - Choose the right label for `Reference point label`.
>    - "Did you compute the matrix with more than one groups of regions?" to `No, I used only one group`
>    - "Clustering algorithm" to `Kmeans clustering`
>    - "Number of clusters to compute" to `2`
> Inspect the output
>
{: .hands_on}
