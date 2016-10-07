Reference based RNA-seq data analysis
=====================================

:grey_question: ***Questions***

- *What are the effects of Pasilla (PS) gene depletion on splicing events?*
- *Second question*
- *Third question*
- *...*

:dart: ***Objectives***

- *Identification of *
- *Second objective*
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

    What is the read length? Is there anything what you find striking?

2. **Trim Galore** :wrench:: Trim low quality bases from the 3' end using  on both paired-end datasets.
3. **FastQC** :wrench:: Re-run and inspect the differences

As the genome of *Drosophila melanogaster* is known, we can use this information and map the sequences on this genome to identify the effects of Pasilla gene depletion on splicing events.

# Mapping

What is mapping? link to mapping tutorial

TopHat introduction....

If you want TopHat to take advantage from already known reference gene annotations, load a reference annotation file into your current Galaxy history. For this exercise please import the Ensembl gene annotation for Drosophila melanogaster (`Drosophila_melanogaster.BDGP5.78.gtf`). From [Zenodo](http://dx.doi.org/10.5281/zenodo.61771)

## Preliminary mapping

*This step is not necessary if you don't need to estimate parameters*

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

Some blabla

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

# Part 2

Short introduction about this subpart.

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

## Subpart 2

Short introduction about this subpart.

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

Some blabla

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

# Conclusion

Conclusion about the technical key points. And then relation between the technics and the biological question to end with a global view.

:grey_exclamation: ***Key Points***

- *Simple sentence to sum up the first key point of the tutorial (Take home message)*
- *Second key point*
- *Third key point*
- *...*

# :clap: Thank you
