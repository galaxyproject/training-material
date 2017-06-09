---
layout: tutorial_hands_on
topic_name: sequence-analysis
tutorial_name: srna
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
> 1. Read quality checking
> 1. Adaptor trimming
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
> 1. Open the Data Upload Manager by selecting *Get Data* from the Tool Panel and clicking *Upload File*
> 1. Select *Paste/Fetch Data*
> 1. Copy each link for the 6 read (.fq), 1 annotation (.gtf), and 2 reference sequence (.fa) files, and paste each link into a separate text field
>    - Set the datatype of the read (.fq) files to **fastq**
>    - Set the datatype of the annotation (.gtf) file to **gtf** and assign the Genome as **dm3**
>    - Set the datatype of the reference (.fa) files to **fasta** and assign the Genome as **dm3**
> 1. Click *Start*
> 1. Rename the files in your history to something meaningful (*e.g.* control_sRNA_rep1.fq)
> 1. Build a *Dataset list* for each set of replicate FASTQ files
>    - Click the *Operations on multiple datasets* check box at the top of the history panel
>    - Check the three boxes next to the control RNAi (control) sRNA-seq samples
>    - Click *For all selected...* and choose *Build dataset list*
>    - Ensure the three control samples are the only ones selected, and enter a name for the new collection (*e.g.* control sRNA-seq)
>    - Click *Create list*
>    - Repeat for the three *klp10A* RNAi samples
>
> {: .hands_on}

## Read quality checking

Read quality scores (phred scores) in FASTQ-formatted data can be encoded by one of a few different encoding schemes. Most Galaxy tools assume that input FASTQ files are using the Sanger/Illumina 1.9 encoding scheme, and if the files are using another scheme, the tools will not interpret the quality score correctly. It is good practice to confirm the quality encoding scheme of your data and then convert to Sanger/Illumina 1.9 if necessary. We can check the quality encoding scheme using the [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) tool (further described in the [NGS-QC tutorial](../../NGS-QC/tutorials/dive_into_qc)).

> ### :pencil2: Hands-on: Quality checking
>
> 1. **FastQC** :wrench:: Run `FastQC` on the FASTQ read files to assess the quality of the reads using the following parameters:
>    - **Short read data from your current history**: Click the "Dataset collection" tab and then select the control sRNA-seq dataset collection
>    - Repeat for the *klp10A* RNAi dataset collection
>
>    > ### :question: Questions
>    >
>    > 1. What quality score encoding scheme is being used for each sample?
>    > 1. What is the read length for each sample?
>    > 1. What does the base/read quality look like for each sample?
>    > 1. Are there any adaptors present in these reads? Which one(s)?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>All samples use the Illumina 1.5 quality encoding scheme. We will need to convert to Sanger/Illumina 1.9. </li>
>    >    <li>All samples have a read length of 51 nt. </li>
>    >    <li>The base quality across the entire length of the reads is good (phred score > 28 for the most part). </li>
>    >    <li>Yes, "Illumina Small RNA 3' Adapters" are present. </li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 1. **FASTQ Groomer** :wrench:: Run `FASTQ Groomer` on the FASTQ read files to convert the quality scores from Illumina 1.5 to Sanger/Illumina 1.9 encoding using the following parameters:
>    - **File to groom**: Click the "Dataset collection" tab and then select the control sRNA-seq dataset collection
>    - **Input FASTQ quality scores type**: Illumina 1.3-1.7
>    - Repeat for the *klp10A* RNAi dataset collection
>
>    ![](../images/image.png)
>
> {: .hands_on}

After `FASTQ Groomer` finishes, click on the groomed control sRNA-seq dataset collection and then click on the name of one of the datasets. You should see that the format is **fastqsanger** instead of **fastq**, meaning we have successfully converted the quality score encoding scheme.

    ![](../images/image.png)

If we go back to the FASTQC output and scroll down to the "Adapter Content" section, we can see that Illumina Small RNA adapters are present in ~80% of our reads. The next step is to remove these artificial adaptors because they will not map to the reference genome. If your reads contain a different adapter, update the **Adapter sequence to be trimmed off** in the `Trim Galore!` step.

    ![](../images/image.png)

## Adaptor trimming

Small RNA sequencing library preparations involve adding an artificial adaptor sequence to both the 5' and 3' ends of the small RNAs. While the 5' adaptor anchors reads to the sequencing surface and thus are not sequenced, the 3' adaptor is typically sequenced immediately following the small RNA sequence. In the example datasets here, the 3' adaptor sequence is `TGGAATTCTCGGGTG`, and needs to be removed from each read before aligning to a reference. We will be using the Galaxy tool `Trim Galore!` which implements the [`cutadapt`](https://cutadapt.readthedocs.io/en/stable/) tool for adapter trimming.

> ### :pencil2: Hands-on: Adaptor trimming
>
> 1. **Trim Galore!** :wrench:: Remove Illumina adapters from the 3' ends of reads by running `Trim Galore!` with the following parameters:
>    - **Is this library paired- or single-end?**: Single-end
>    - **Reads in FASTQ format**: Click the "Dataset collection" tab and then select the control sRNA-seq dataset
>    - **Trimming reads?**: User defined adapter trimming
>    - **Adapter sequence to be trimmed off**: TGGAATTCTCGGGTG
>    - **Trim Galore! advanced settings**: Full parameter list
>    - **Trim low-quality ends from reads in addition to adapter removal**: 0
>    - **Overlap with adapter sequence required to trim a sequence**: 6
>    - **Discard reads that became shorter than length INT**: 12
>    - **Generate a report file**: Yes
>
>    ![](../images/image.png)
>
>    We don't want to trim for quality because the adapter-trimmed sequences represent a full small RNA molecule, and we want to maintain the integrity of the entire molecule. We increase the minimum read length required to keep a read because small RNAs can potentially be shorter than 20 nt (the default value). We can check out the report file for any sample and see the command for the tool, a summary of the total reads processed and number of reads with an adapter identified, and histogram data of the length of adaptor trimmed. We also see that a very small percentage of low-quality bases have been trimmed
>
> 1. **FastQC** :wrench:: Re-run `FastQC` on trimmed reads and inspect the differences.
>
>    > ### :question: Questions
>    >
>    > 1. What is the quality encoding scheme?
>    > 1. What is the read length?
>    > 1. Are there any adaptors present in these reads? Which one(s)?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>All samples are using the Sanger/Illumina 1.9 quality encoding scheme.</li>
>    >    <li>The read lengths range from 12 to 51 nt after trimming.</li>
>    >    <li>No, Illumina Small RNA 3' Adaptors are no longer present. No other adapters are present.</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> ![](../images/image.png)
>
> {: .hands_on}

Now that we have trimmed our reads of the Illumina Small RNA 3' adaptors and converted, we will align our trimmed reads to the reference *Drosophila* genome (dm3). For miRNA analyses, it is useful to align to the reference set of known miRNAs first, and then re-align any unaligned reads to the reference genome.

## Read alignment

To quantify small RNA abundance and identify their putative targets, we need to know where the sequenced reads align to a reference genome. In the case of a eukaryotes, some small RNAs are transcribed from mRNA templates, which means that some small RNAs can originate from an exon-exon (spliced) boundary. Therefore, a splice-aware aligner must be used to account for this possibility. [`HISAT2`](https://ccb.jhu.edu/software/hisat2/index.shtml) is an accurate and fast tool for aligning spliced reads to a genome, and we will be using `HISAT2` in this tutorial.

> ### :pencil2: Hands-on: Heirarchical alignment to rRNA, miRNA, and genome reference sequences
>
> 1. **HISAT2** :wrench:: Run `HISAT2` to align one collection of trimmed reads to reference rRNA sequences with the following parameters:
>    - **Single end or paired reads?**: Individual unpaired reads
>    - **Reads**: Click the "Dataset collection" tab and then select the control sRNA-seq dataset of trimmed FASTQ files
>    - **Source for the reference genome to align against**: Use a genome from history
>    - **Select the reference genome**: Dmel_rRNA_sequences.fa
>    - **Spliced alignment parameters**: Specify spliced alignment parameters
>    - **Specify strand-specific information**: First Strand (R/RF)
>
>       ![](../images/image.png)
>
>    We now need to extract the *unaligned* reads from the output BAM file for aligning to reference miRNA sequences. We can do this by using the `Filter SAM or BAM, output SAM or BAM` tool to obtain reads with the bit flag 4 (meaning the read is unaligned) and then converting the filtered BAM file to FASTQ format with the `Convert from BAM to FastQ` tool.
>
> 1. **Filter SAM or BAM, output SAM or BAM** :wrench:: Run `Filter SAM or BAM, output SAM or BAM` on one collection of HISAT output BAM files with the following parameters:
>    - **SAM or BAM file to filter**: Click the "Dataset collection" tab and then select the control sRNA-seq dataset of aligned HISAT BAM files
>    - **Filter on bitwise flag**: Yes
>    - **Only output alignments with all of these flag bits set**: Check the box next to "The read in unmapped"
>
> 1. **Convert from BAM to FastQ** :wrench:: Run `Convert from BAM to FastQM` on one collection of filtered HISAT output BAM files with the following parameters:
>    - **Convert the following BAM file to FASTQ**: Click the "Dataset collection" tab and then select the control sRNA-seq dataset of filtered HISAT BAM files
>
>    Next we will align the non-rRNA reads to a known set of miRNA sequences to remove miRNAs from the whole genome alignment step.
>
> 1. **HISAT2** :wrench:: Run `HISAT2` to align one collection of non-rRNA reads to reference miRNA sequences using the following parameters:
>    - **Single end or paired reads?**: Individual unpaired reads
>    - **Reads**: Click the "Dataset collection" tab and then select the control sRNA-seq dataset of non-rRNA FASTQ files
>    - **Source for the reference genome to align against**: Use a genome from history
>    - **Select the reference genome**: Dmel_miRNA_sequences.fa
>    - **Spliced alignment parameters**: Specify spliced alignment parameters
>    - **Specify strand-specific information**: First Strand (R/RF)
>
>       ![](../images/image.png)
>
>    We now need to extract *unaligned* reads from the output BAM file for aligning to reference genome sequence. Repeat the `Filter SAM or BAM, output SAM or BAM` and `Convert from BAM to FastQ` steps to do this. Rename the converted FASTQ files something meaningful (*e.g.* "non-rRNA/miRNA control RNAi sRNA").
>
> 1. **HISAT2** :wrench:: Run `HISAT` to align a collection of non-rRNA/non-miRNA reads to the reference genome using the following parameters.
>    - **Single end or paired reads?**: Individual unpaired reads
>    - **Reads**: Click the "Dataset collection" tab and then select the control sRNA-seq dataset of non-rRNA/non-miRNA FASTQ files
>    - **Source for the reference genome to align against**: Use a built-in genome
>    - **Select a reference genome**: Fruit Fly (D. melanogaster): dm3
>    - **Primary alignments**: 1000000
>    - **Spliced alignment parameters**: Specify spliced alignment parameters
>    - **Specify strand-specific information**: First Strand (R/RF)
>
>    We use a ridiculously high value for "Primary alignments" (-k) here to ensure that we retain piRNAs. In flies, piRNA often align to transposable and repeat elements in the genome. If a piRNA read aligns to more than -k loci, it is removed from the output. We want to make sure we don't lose these reads, so we set this parameter to something extremely high. This slows down the alignment step, but is necessary.
>
> 1. Repeat the alignment steps for the *klp10A* RNAi sRNA-seq dataset.
>
> {: .hands_on}

**UPDATES STOPPED HERE**

## Small RNA annotation
**TODO Describe what is meant by small RNA annotation. Describe sense v. antisense meaning. Describe that different classes of small RNAs are in a small RNA-seq library. Talk about size distribution and nt composition biases in fly piRNAs and in other classes of small RNAs.**

> ### :pencil2: Hands-on: Small RNA annotation
>
> 1. **Tool** :wrench:: Run `Tool` on one collection of `HISAT` alignments using the default parameters.
>    - Use batch mode to run all four samples from one tool form.
>
> ![](../images/image.png)
>
> {: .hands_on}

## Small RNA abundance estimation

**TODO RPM for small RNA counts.**

> ### :pencil2: Hands-on: Small RNA abundance estimation
>
> 1. **Tool** :wrench:: Run `Tool` on the `Tool`-annotated small RNAs.
>    - Use batch mode to inlcude all four `Stringtie` assemblies.
>    - **Use Reference Annotation**: Yes, then select the "RefSeq GTF mm10" file.
> ![](../images/image.png)
>
> {: .hands_on}
>

## Small RNA differential expression testing

### Analysis of the differential gene expression

We want to identify which small RNAs are differentially expressed between the control and *klp10A* RNAi conditions. To do this we will implement a counting approach using `FeatureCounts` to count small RNAs per genome feature. Then we will provide this information to `DESeq2` to generate normalized counts and significance testing for differential expression.

### Count the number of reads per genome feature

To compare the abundance of small RNAs mapping to genomic features between different RNAi conditions, the first essential step is to quantify the number of small RNAs per feature. [`FeatureCounts`](http://bioinf.wehi.edu.au/featureCounts/) is one of the most popular tools for counting reads in genomic features. In our case, we'll be using `FeatureCounts` to count small RNA reads aligning to dm3 genomic features in a custom GTF file that contains transposon and repeat elements, which are common targets of the piRNA subclass of small RNAs.

The recommended mode is "union", which counts overlaps even if a read only shares parts of its sequence with a genomic feature and disregards reads that overlap more than one feature.

> ### :pencil2: Hands-on: Counting the number of small RNA reads per feature
>
> 1. **FeatureCounts** :wrench:: Run `FeatureCounts` on the aligned reads (`HISAT` output) using the gene+transposon+repBase GTF as the annotation file.
>
>    - **TODO**
>
> ![](../images/image.png)
>
> {: .hands_on}
>

### Perform differential "expression" testing

TODO

[`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) is a great tool for differential expression analysis. It takes read counts produced by `FeatureCounts` and applies size factor normalization.

> ### :pencil2: Hands-on: Differential expression testing
>
> 1. **DESeq2** :wrench:: Run `DESeq2` with the following parameters:
>
>    - **TODO**
>
> 1. **Filter** :wrench:: Run `Filter` to extract feature with a significant differences in small RNA mappings (adjusted *p*-value less than 0.05) between control and *klp10A* RNAi conditions.
>
>    - **TODO**
>
>    > ### :question: Question
>    >
>    > How many features have a significant change in mapped small RNAs between these conditions?
>    >
>    > <details>
>    > <summary>Click to view answers</summary>
>    > To filter, use "c7 lessthan 0.05". And we get ## features with a significant change in mapped small RNAs.
>    > </details>
>    {: .question}
>
> 1. **Filter** :wrench:: Determine how many features have increased or decreased mapped small RNA counts.
>
>    - **TODO**
>
> {: .hands_on}

For more information about `DESeq2` and its outputs, you can have a look at [`DESeq2` documentation](https://www.bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf).

## Small RNA and mRNA integration

**TODO**

## Visualization

Now that we have a list of transcript expression levels and their differential expression levels, it is time to visually inspect our transcript structures and the reads they were predicted from. It is a good practice to visually inspect (and present) loci with transcripts of interest. Fortuantely, there is a built-in genome browser in Galaxy, **Trackster**, that make this task simple (and even fun!).

In this last section, we will convert our aligned read data from BAM format to bigWig format to simplify observing where our stranded RNA-seq data aligned to. We will then initiate a session on Trackster, load it with our data, and visually inspect our interesting loci.

> ### :pencil2: Hands-on: Visualizing data on genome browser
>
> 1. **bamCoverage** :wrench:: Run `bamCoverage` on all four aligned read files (`HISAT` output) with the following parameters:
>
>    - **TODO**
>
> 1. **Viz** :wrench:: On the center console at the top of the Galaxy interface, choose " Visualization" -> "New track browser"
>
>    - **TODO**
>
> {: .hands_on}

## Conclusion

**TODO**

![](../images/schematic_for_sRNAseq_tutorial.png)
