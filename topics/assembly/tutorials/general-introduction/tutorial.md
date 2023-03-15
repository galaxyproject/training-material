---
layout: tutorial_hands_on

title: "An Introduction to Genome Assembly"
zenodo_link: "https://doi.org/10.5281/zenodo.582600"
questions:
  - "How do we perform a very basic genome assembly from short read data?"
objectives:
  - "assemble some paired end reads using Velvet"
  - "examine the output of the assembly."
time_estimation: "30m"
level: Introductory
key_points:
  - "We assembled some Illumina fastq reads into contigs using a short read assembler called Velvet"
  - "We showed what effect one of the key assembly parameters, the k-mer size, has on the assembly"
  - "It looks as though there are some exploitable patterns in the metric data vs the k-mer size."
contributors:
  - slugger70
---

# Genome assembly with Velvet: Background


Velvet is one of a number of *de novo* assemblers that use short read sets as input (*e.g.* Illumina Reads). The assembly method is based on the manipulation of de Bruijn graphs, via the removal of errors and the simplification of repeated regions.

> <comment-title></comment-title>
>
> For information about Velvet, you can check its (nice) [Wikipedia page](https://en.wikipedia.org/wiki/Velvet_assembler).
{: .comment}

For this tutorial, we have a set of reads from an imaginary *Staphylococcus aureus* bacterium with a miniature genome (197,394 bp). Our mutant strain read set was sequenced with the whole genome shotgun method, using an Illumina DNA sequencing instrument. From these reads, we would like to rebuild our imaginary *Staphylococcus aureus* bacterium via a *de novo* assembly of a short read set using the Velvet assembler.

We also have a sequence for a reference genome that we will use later in the tutorial.

> <agenda-title></agenda-title>
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get the data

> <hands-on-title>Getting the data</hands-on-title>
>
> 1. Create and name a new history for this tutorial.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import from [Zenodo](https://doi.org/10.5281/zenodo.582600) or from the data library the files:
>
>    ```
>    https://zenodo.org/record/582600/files/mutant_R1.fastq
>    https://zenodo.org/record/582600/files/mutant_R2.fastq
>    https://zenodo.org/record/582600/files/wildtype.fna
>    ```
>
>    > <tip-title>Importing data via links</tip-title>
>    >
>    > * Copy the link location (Right-click on the filename then "Copy Link Address")
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * For the read files, change the data-type to **fastqsanger**
>    > * Press **Start**
>    {: .tip}
>
> 3. Change the name of the files to `mutant_R1`, `mutant_R2` and `wildtype.fna`.
>
>    As a default, Galaxy uses the link as the name of the new dataset. It also does not link the dataset to a database or a reference genome.
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
> 4. Inspect the content of a file.
>
>    > <tip-title>Inspecting the content of a dataset</tip-title>
>    >
>    > * Click on the {% icon galaxy-eye %} (eye) icon next to the relevant history entry
>    > * View the content of the file in the central panel
>    {: .tip}
>
>    > <question-title></question-title>
>    >
>    > 1. What are four key features of a FASTQ file?
>    > 2. What is the main difference between a FASTQ and a FASTA file?
>    >
>    > > <solution-title></solution-title>
>    > > 1. Each sequence in a FASTQ file is represented by 4 lines: 1st line is the id, 2nd line is the sequence, 3rd line is not used, and 4th line is the quality of sequencing per nucleotide
>    > > 2. In a FASTQ file, not only are the sequences present, but information about the quality of sequencing is also included.
>    > {: .solution }
>    {: .question}
>
{: .hands_on}

The reads have been sequenced from an imaginary *Staphylococcus aureus* bacterium using an Illumina DNA sequencing instrument. We obtained the 2 files we imported (`mutant_R1` and `mutant_R2`)

> <question-title></question-title>
>
> Why do we have 2 files here if we only sequenced the bacteria once?
>
> > <solution-title></solution-title>
> > 1. The bacteria has been sequenced using paired-end sequencing. The first file corresponds to forward reads and the second file to reverse reads.
> {: .solution }
{: .question}

# Evaluate the input reads

Before doing any assembly, the first questions you should ask about your input reads include:

- What is the coverage of my genome?
- How good is my read set?
- Do I need to ask for a new sequencing run?
- Is it suitable for the analysis I need to do?

We will evaluate the input reads using the FastQC tool. This tool runs a standard series of tests on your read set and returns a relatively easy-to-interpret report. We will use it to evaluate the quality of our FASTQ files and combine the results with MultiQC.

> <hands-on-title>FastQC on a fastq file</hands-on-title>
>
> 1. **FastQC** {% icon tool %} with the following parameters
>    - {% icon param-files %} *"Raw read data from your current history"*:  `mutant_R1.fastq` and `mutant_R2.fastq`
>
>    {% snippet faqs/galaxy/tools_select_multiple_datasets.md %}
>  
> 2. **MultiQC** {% icon tool %} with the following parameters
>    - *"Software name"*: `FastQC`
>    - *"Result file"*:  All raw data files generated by FastQC.
>
{: .hands_on}

MultiQC generates a webpage combining reports for FastQC on both datasets. It includes these graphs and tables:

- General statistics

    This is important in setting maximum k-mer size for an assembly.

    > <comment-title>Getting the length of sequences</comment-title>
    >
    > * Click on **Configure Columns**
    > * Check **Length**
    > * Close the window
    {: .comment}

    > <question-title></question-title>
    >
    > 1. How long are the sequences?
    > 2. What is the average coverage of the genome, given our imaginary *Staphylococcus aureus* bacterium has a genome of 197,394 bp?
    >
    > > <solution-title></solution-title>
    > > 1. The sequences are 150 bp long
    > > 2. We have 2 x 12,480 sequences of 150 bp, so the average genome coverage is: 2 * 12480 * 150 / 197394, or approximately 19 X coverage.
    > {: .solution }
    {: .question}

- Sequence Quality Histograms

    Dips in quality near the beginning, middle or end of the reads may determine the trimming/cleanup methods and parameters to be used, or may indicate technical problems with the sequencing process/machine run.

    ![Sequence Quality Histograms with the mean quality value across each base position in the read](../../images/fastqc_per_base_sequence_quality_plot.png "The mean quality value across each base position in the read")

    > <question-title></question-title>
    >
    > 1. What does the y-axis represent?
    > 2. Why is the quality score decreasing across the length of the reads?
    >
    > > <solution-title></solution-title>
    > > 1. The y-axis represents the quality score for each base (an estimate of the error during sequencing).
    > > 2. The quality score is decreasing accross the length of the reads because the sequencing become less and less reliable at the end of the reads.
    > {: .solution }
    {: .question}

- Per Sequence GC Content

    High GC organisms tend not to assemble well and may have an uneven read coverage distribution.

- Per Base N Content

    The presence of large numbers of Ns in reads may point to a poor quality sequencing run. You will need to trim these reads to remove Ns.

- k-mer content

    The presence of highly recurring k-mers may point to contamination of reads with barcodes or adapter sequences.


> <comment-title></comment-title>
>
> For a fuller discussion of FastQC outputs and warnings, see the [FastQC website link](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), including the section on each of the output [reports](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/), and examples of ["good"](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html) and ["bad"](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html) Illumina data.
{: .comment}

We won't be doing anything to these data to clean it up as there isn't much need. Therefore we will get on with the assembly!


# Assemble reads with Velvet

Now, we want to assemble our reads to find the sequence of our imaginary *Staphylococcus aureus* bacterium. We will perform a *de novo* assembly of the reads into long contiguous sequences using the Velvet short read assembler.

The first step of the assembler is to build a de Bruijn graph. For that, it will break our reads into k-mers, *i.e.* fragments of length *k*. Velvet requires the user to input a value of *k* (k-mer size) for the assembly process. Small k-mers will give greater connectivity, but large k-mers will give better specificity.

> <hands-on-title>Assemble the reads</hands-on-title>
>
> 1. **FASTQ interlacer** {% icon tool %} with the following parameters
>    - *"Type of paired-end datasets"*: `2 separate datasets`
>    - *"Left-hand mates"*: `mutant_R1.fastq`
>    - *"Right-hand mates"*: `mutant_R2.fastq`
>
>    Currently our paired-end reads are in 2 files (one with the forward reads and one with the reverse reads), but Velvet requires only one file, where each read is next to its mate read. In other words, if the reads are indexed from 0, then reads 0 and 1 are paired, 2 and 3, 4 and 5, etc. Before doing the assembly *per se*, we need to prepare the files by combining them.
>
> 2. **velveth** {% icon tool %} with the following parameters
>    - *"Hash Length"*: `29`
>    - *"Input Files"*
>       - Click on {% icon param-repeat %} *"Input Files"*
>       - In "1: Input Files"
>          - *"file format"*: `fastq`
>          - *"read type"*: `shortPaired reads`
>          - {% icon param-files %} *"Dataset"*:  outputs of **FASTQ interlacer** 
>
>    The tool takes our reads and break them into k-mers.
>
> 3. **velvetg** {% icon tool %} with the following parameters
>    - {% icon param-files %} *"Velvet Dataset"*:  outputs of **velveth**
>    - *"Using Paired Reads"*: `Yes`
>
>    This last tool actually does the assembly.
{: .hands_on}

Two files are generated:

- A **Contigs** file

    This file contains the sequences of the contigs longer than 2k. In the header of each contig, a bit of information is added:
    - the k-mer length (called "length"): For the value of k chosen in the assembly, a measure of how many k-mers overlap (by 1 bp each overlap) to give this length
    - the k-mer coverage (called "coverage"): For the value of k chosen in the assembly, a measure of how many k-mers overlap each base position (in the assembly).

    ![Contigs output](../../images/image10.png)

- A **Stats** file

    This is a tabular file giving for each contig the k-mer lengths, k-mer coverages and other measures.

    ![Contigs stats output](../../images/image11.png)

# Collect some statistics on the contigs

> <question-title></question-title>
>
> 1. How many contigs have been built?
> 2. What is the mean, min and max length of the contigs?
>
> > <solution-title></solution-title>
> > 1. 190
> > 2. To compute this information, we can use the Datamash tool on the 2nd columns (length). Be careful with the first line, the header. As a result, we obtain: 597.82 as mean, 1 as min and 12904 as max. It would mean that the smallest contig has a length of 1 bp, even smaller than k. The length on the 2nd column corresponds to length of the contig in k-mers. This means that the smallest contig has a length of 1k = 29. So to obtain the real length, we need to add k-1 to the length. We then obtain a mean contig length of 625.82 bp, a min contig of 29 bp and a max contig of 12,932 bp.
> {: .solution }
{: .question}

This table is limitted, but we will now collect more basic statistics on our assembly.

> <hands-on-title>Collect fasta statistics on our contigs</hands-on-title>
>
> 1. **Quast** {% icon tool %} with
>    - *"Contigs/scaffolds output file"*: outputs of **velvetg**
>    - *"Type of assembly"*: `Genome`
>    - *"Use a reference genome?"*: `Yes`
>    - *"Reference genome"*: `wildtype.fna`
>    - *"Type of organism"*: `Prokaryotes`
>    - *"Lower Threshold"*: `500`
>    - *"Comma-separated list of contig length thresholds"*: `0,1000`
{: .hands_on}

This tool generates 5 output files, but we will focus on the HTML report and the Icarus viewer.

> <question-title></question-title>
>
> 1. What is represented in the Icarus viewer?
>
> > <solution-title></solution-title>
> > 1. Icarus is a novel genome visualizer for accurate assessment and analysis of genomic draft assemblies. It draws contigs ordered from longest to shortest, highlights N50, N75 (NG50, NG75) and long contigs larger than a user-specified threshold
> {: .solution }
{: .question}

The HTML report reports many statistics computed by QUAST to assess the quality of the assembly:

- Statistics about the quality of the assembly when compared to the reference (fraction of the genome, duplication ratio, etc)
- Misassembly statistics, including the number of misassemblies

    A misassembly is a position in the contigs (breakpoints) that satisfy one of the following criteria:
    - the left flanking sequence aligns over 1 kbp away from the right flanking sequence on the reference;
    - flanking sequences overlap on more than 1 kbp
    - flanking sequences align to different strands or different chromosomes

- Unaligned regions in the assembly
- Mismatches compared to the reference genomes
- Statistics about the assembly *per se*, such as the number of contigs and the length of the largest contig

> <question-title></question-title>
>
> 1. How many contigs have been constructed?
> 2. Which proportion of the reference genome do they represent?
> 3. How many misassemblies have been found?
> 4. Has the assembly introduced mismatches and indels?
> 5. What are N50 and L50?
> 6. Is there a bias in GC percentage induced by the assembly?
>
> > <solution-title></solution-title>
> > 1. 190 contigs have been constructed, but only 47 have a length > 500 bp.
> > 2. The contigs represents 87.965% of the reference genome.
> > 3. 1 misassembly has been found: it corresponds to a relocation, *i.e.* a misassembly event (breakpoint) where the left flanking sequence aligns over 1 kbp away from the right flanking sequence on the reference genome.
> > 4. 8.06 mismatches per 100 kbp and 4.03 indels per 100 kbp are found.
> > 5. N50 is the length for which the collection of all contigs of that length or longer covers at least half an assembly. In other words, if contigs were ordered from small to large, half of all the nucleotides will be in contigs this size or larger. And L50 is the number of contigs equal to or longer than N50: L50 is the minimal number of contigs that cover half the assembly.
> > 6. The GC % in the assembly is 33.64%, really similar to the one of the reference genome (33.43%).
> {: .solution }
{: .question}

# Discussion

> <hands-on-title>(Optional) Rerun for values <i>k</i> ranging from 31 to 101</hands-on-title>
>
> 1. **velveth** {% icon tool %} with the same parameters as before except
>    - *"Hash Length"*: a value between 31 and 101
> 2. **velvetg** {% icon tool %} with the same parameters as before
> 3. **Quast** {% icon tool %} with the same parameters as before
{: .hands_on}

We have completed an assembly on this data set for a number of k values ranging from 29 to 101. A few of the assembly metrics appear below.

![contigs](../../images/number_of_contigs.png "Number of contigs in the assembly for various k-mer sizes")

![largest_contig](../../images/largest_contig.png "Largest contig in each of the assemblies by k-mer size")

![total_bp](../../images/total_bp.png "Total number of base pairs in all the contigs for each assembly by k-mer size")

![n50](../../images/n50.png "N50 metric for each of the assemblies by k-mer size")

> <question-title></question-title>
>
> 1. Are there any distinct features in the charts?
> 2. Does it look like one assembly might be better than some of the others?
>
{: .question}

The reasons for these patterns will be discussed in detail in the De Bruijn graph assembly slides and tutorial.
