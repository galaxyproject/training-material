Reference based RNA-seq data analysis
=====================================

:grey_question: ***Questions***

- *What are the effects of Pasilla (PS) gene depletion on splicing events?*
- *Second question*
- *Third question*
- *...*

:dart: ***Objectives***

- *Analysis of RNA sequencing data using a reference genome*
- *Analysis of the differential gene expression*
- *Third objective*
- *...*

:heavy_check_mark: ***Requirements***

- *Galaxy introduction*
- *Quality control*
- *Mapping*
- *...*

:hourglass: ***Time estimation*** *1d/3h/6h*

[:book: **Associated slide deck**](http://bgruening.github.io/training-material/RNA-Seq/slides/ref_based.html)

# Introduction

In this exercise, we use RNA-seq data from the study by [Brooks et al. 2011](http://genome.cshlp.org/content/21/2/193.long). In this study, the Pasilla (PS) gene was depleted in *Drosophila melanogaster* by RNAi. The authors wanted to analyze the effects of Pasilla gene depletion on splicing events.

Total RNA was isolated and used for preparing either single-end or paired-end RNA-seq libraries used for sequencing. RNA sequencing data are then available. The genome of *Drosophila melanogaster* is known and can be used for this analysis. The effects of Pasilla gene depletion on splicing events are then analyzed with reference based RNA-seq data analysis.

Reference based RNA-seq data analysis is ...
*add a small introduction about reference based RNA-seq data analysis*

In this tutorial, we will analyze the data with:

1. Pretreatments
2. Mapping
3. ...

# Pretreatments

## Data upload

The data is available at NCBI Gene Expression Omnibus (GEO) under accession number [GSE18508](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18508).
Here, we have 7 samples:
- 3 treated samples with Pasilla (PS) gene depletion (called "pasilla-depleted")
- 4 untreated samples (called "wt")

Each sample constitutes a separate biological replicate of the corresponding condition (treated or untreated). Moreover, two of the treated and two of the untreated samples are from a paired-end sequencing assay, while the remaining samples are from a single-end sequencing experiment.

:pencil2: ***Hands on!*** Data upload

1. Create a new history for this RNA-seq exercise
2. Import a FASTQ file pair (*e.g.* [`GSM461177_untreat_paired_subset_1`](https://zenodo.org/record/61771/files/GSM461177_untreat_paired_subset_1.fastq) and [`GSM461177_untreat_paired_subset_2`](https://zenodo.org/record/61771/files/GSM461177_untreat_paired_subset_2.fastq)) from [Zenodo](http://dx.doi.org/10.5281/zenodo.61771)

    > :bulb: **Importing data via links**
    > * Copy the link location
    > * Open the Galaxy Upload Manager
    > * Select **Paste/Fetch Data**
    > * Paste the link into the text field
    > * Press **Start**    
    <br>

    > :bulb: **Changing the file type `fastq` to `fastqsanger` once data is in history**
    > * Click on the pencil button displayed in your dataset in the history
    > * Choose **Datatype** on the top
    > * Select `fastqsanger`
    <br>
    <br>

    > :+1: **Edit the "Database/Build" to select "dm3"**
    <br>

    > :+1: **Rename the datasets according to the samples**
    > As default, Galaxy takes the link as name.
    <br>

Both files contain the first 100.000 paired-end reads of one sample. The sequences are raw sequences from sequencing. They needs to be controlled for their quality.

## Quality control

For quality control, we use similar tools as described in [NGS-QC tutorial](https://github.com/bgruening/training-material/blob/master/NGS-QC/tutorials/dive_into_qc.md).

:pencil2: ***Hands on!*** Quality control

1. **FastQC** :wrench:: Run on one of the two FASTQ files to control the quality of the reads as

    :question: What is the read length? Is there anything what you find striking?

2. **Trim Galore** :wrench:: Trim low quality bases from the 3' end using  on both paired-end datasets.
3. **FastQC** :wrench:: Re-run and inspect the differences

As the genome of *Drosophila melanogaster* is known, we can use this information and map the sequences on this genome to identify the effects of Pasilla gene depletion on splicing events.

# Mapping

To make sense of the reads, their positions within *Drosophila melanogaster* genome must be determined. This process is known as aligning or 'mapping' the read to the reference genome.

> Want to learning more about mapping? Follow our [training](http://bgruening.github.io/training-material/NGS-mapping/slides)

For this tutorial, we will use [TopHat](https://ccb.jhu.edu/software/tophat/index.shtml), a mapping tool developed for RNA-Seq reads. It aligns RNA-Seq reads to mammalian-sized genomes using [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml), and then analyzes the mapping results to identify splice junctions between exons.

To help annotations of RNA sequences, we can take advantage from already known reference gene annotations.

:pencil2: ***Hands on!***

1. Load the Ensembl gene annotation for *Drosophila melanogaster* ([`Drosophila_melanogaster.BDGP5.78.gtf`](https://zenodo.org/record/61771/files/Drosophila_melanogaster.BDGP5.78.gtf)) from [Zenodo](http://dx.doi.org/10.5281/zenodo.61771) into your current Galaxy history

TopHat also needs to know two important parameters about the sequencing library

- The *strandedness*: *unstranded* or *stranded*?  and if *stranded* there are many types)
- The *inner distance* between the two reads for paired-end data

**These information should usually come with your FASTQ files!!!** If not, try to find them on the site where you downloaded the data or in the corresponding publication. Another option is to estimate these parameters with a *preliminary mapping* of a *downsampled* file and some analysis programs. Afterward, the actual mapping can be redo on the original files with the optimized parameters.

## Preliminary mapping

*This step is not necessary if you don't need to estimate parameters*

:pencil2: ***Hands on!***

1. **Select first** :wrench:: Downsample the FASTQ file to 200k to 1M reads

    :warning: For the provided files downsampling is not necessary as they only contain 100k reads

2. **TopHat** :wrench:: Run **TopHat** with:
    - "Paired-end (as individual datasets)" instead of "Single-end"
    - "Drosophila melanogaster: dm3" as reference genome
    - the defaults for *strandedness* and *insert size*
    
3. **Inner Distance** :wrench:: Run **Inner Distance** on the BAM file using the `Drosophila_melanogaster.BDGP5.78.gtf` reference gene model to estimate the *inner distance*
4. Inspect the resulting PDF

    :question: What is the mean value for the inner distance?

    If you already have read the corresponding paper carefully you might know that the fragment size is ~200bp. With read lengths of 2x37bp an educated guess could also be `125` for the inner distance. It's up to you decision, which value you prefer...

5. **Infer Experiment** :wrench:: Run **Infer Experiment** with the same files
6. Check the results and search the tool's documentation for help on the meaning.

As it is sometimes quite difficult to find out which settings correspond to those of other programs, the following table might be helpful to identify the library type:

Library type | **Infer Experiment** | **TopHat** | **HISAT2** | **htseq-count** | **featureCounts**
--- | --- | --- | --- | --- | ---
PE | "1++,1--,2+-,2-+" | "FR Second Strand" | "FR" | "yes" | "1"
PE | "1+-,1-+,2++,2--" | "FR First Strand" | "RF" | "reverse" | "2"
SE | "++,--" | "FR Second Strand" | "F" | "yes" | "1"
SE | "+-,-+" | "FR First Strand" | "R" | "reverse" | "2"
SE,PE | undecided | "FR Unstranded" | default | "no" | "0"

## Actual mapping

With the sequencing library parameters, the full RNA sequences can be mapped on the *Drosophila melanogaster* genome.

:pencil2: ***Hands on!***

1. **TopHat** :wrench:: Run **TopHat** with the full parameter set to get the best mapping results:
    - "Paired-end (as individual datasets)" instead of "Single-end"
    - "Mean Inner Distance" to "112" (or "125"?)
    - "Drosophila melanogaster: dm3" as reference genome
    - "Full parameter list" for "TopHat settings to use"
    - "FR First Strand" as "Library type"
    - "18" for the "Minimum length of read segments"

        By default, TopHat proposes to fix the minimum length of read segments to 25, but a value of `18` seems to be a more appropriate value for this input data.

        :question: Why?

    - "Yes" for use of own junction data
    - "Yes" for use of Gene Annotation Model
    - `Drosophila_melanogaster.BDGP5.78.gtf` as Gene Model Annotations (to enable transcriptome alignment)
    - "Yes (--coverage-search)" to use coverage-based search for junctions

        The TopHat algorithm splits reads into segments to map the reads across splice junctions. Coverage-based search for junctions increases the sensitivity

**TopHat** generated a BAM file with the mapped reads and three BED files containing splice junctions, insertions and deletions.

The mapping exercise worked for you? Great! :tada:

## Inspection of TopHat results

However, the datasets are too small to get you a good impression of how real data looks like. So we run TopHat for you on a real dataset. We extract only the reads mapped to Chromosome 4 of *Drosophila*.

:pencil2: ***Hands on!***

1. Import from [Zenodo](http://dx.doi.org/10.5281/zenodo.61771) the following files:
    - [`GSM461177_untreat_paired_chr4.bam`](https://zenodo.org/record/61771/files/GSM461177_untreat_paired_chr4.bam)
    - [`GSM461177_untreat_paired_deletions_chr4.bed`](https://zenodo.org/record/61771/files/GSM461177_untreat_paired_deletions_chr4.bed)
    - [`GSM461177_untreat_paired_insertions_chr4.bed`](https://zenodo.org/record/61771/files/GSM461177_untreat_paired_insertions_chr4.bed)
    - [`GSM461177_untreat_paired_junctions_chr4.bed`](https://zenodo.org/record/61771/files/GSM461177_untreat_paired_junctions_chr4.bed)
2. **IGV** :wrench:: Visualize this BAM file and the three BED files

    You might for example inspect the region on Chromosome 4 between 560 kb to 600 kb (`chr4:560,000-600,000`)

    :question: Which information does each of the BED files contain?

    > :+1: **Change the data type from "tabular" to "bed"**
    <br>

3. **IGV** :wrench:: Inspect the results using a **Sashimi plot**

    > :bulb: **Creation of a Sashimi plot**
    > * Right click on the bam file
    > * Select **Sashimi Plot** from the context menu
    <br>

4. **IGV** :wrench:: Look around to find other regions with in interesting junctions, *e.g.* `chr4:870,000-940,000`.

# Analysis of the differential gene expression

## Count the number of reads per annotated gene

## Analysis of the differential gene expression

## Analysis of the functional enrichment among differentially expressed genes

# Conclusion

Conclusion about the technical key points. And then relation between the technics and the biological question to end with a global view.

:grey_exclamation: ***Key Points***

- *Simple sentence to sum up the first key point of the tutorial (Take home message)*
- *Second key point*
- *Third key point*
- *...*

# :clap: Thank you
