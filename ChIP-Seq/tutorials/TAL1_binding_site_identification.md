Identification of the binding sites of the T-cell acute lymphocytic leukemia protein 1 (TAL1) with ChIP-sequencing
=============

:grey_question: ***Questions***

- *How is raw ChIP-seq data processed and analyzed?*
- *What are the binding sites of TAL1?*
- *Which genes are regulated by TAL1?*
- *...*

:dart: ***Objectives***

- *Analysis of ChIP-seq data*
- *Mastering the ChIP-seq data analysis workflow*
- *Identifying TAL1 binding sites*
- *Determining genes regulated by TAL1*
- *...*

:heavy_check_mark: ***Requirements***

- *Galaxy introduction*
- *NGS-QC*
- *NGS-mapping*
- *IGV introduction*

:hourglass: ***Time estimation*** *3h*

# Introduction

This tutorial uses ChIP-seq datasets from a study published by [Wu et al., 2012](http://genome.cshlp.org/content/24/12/1945.full.pdf+html).
The goal of this study was to investigate "the dynamics of occupancy and the role in gene regulation of the transcription factor TAL1, a critical regulator of hematopoiesis, at multiple stages of hematopoietic differentiation."
To this end, ChIP-seq was performed in the G1E cell line - a GATA-null immortalized cell line derived from targeted disruption of GATA-1 in mouse embryonic stem cells - and megakaryocytes.
This dataset (GEO Accession [GSE51338](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51338)) consists of biological replicate TAL1 ChIP-seq experiments and input control experiments for which the same treatment as the ChIP-seq samples was done except for the immunoprecipitation step.
Input control experiments are used to identify and remove sampling bias, for example open/accessible chromatin or GC bias.

Because of the long processing time for the large original files, we have downsampled the original raw data files to include only a subset of genomic loci.

# Pretreatments

As for any NGS data analysis, ChIP-seq data must be [quality controlled](../../NGS-QC/slides/dive_into_qc.html) before being aligned to a reference genome.

:pencil2: ***Hands on!***

1. Import the datasets
2. Control quality of sequences using `FastQC` and `TrimGalore` as explained in [NGS-QC tutorial](./../NGS-QC/tutorials/dive_into_qc.md)

    We recommend you trim low quality bases with a cut-off of 15 and adapter sequences with an overlap of at least 3 and also to clip a possibly trailing adapter sequence

3. Align the reads to the mouse genome (version mm10) using `bowtie2`

    How many reads were mapped?

4. Visualize the reads in IGV and zoom in to position ##-## on chromosome 11 (e.g. `chr11:##-##`)

    Notice that the reads have a _direction_ (i.e., they are mapped to the forward or reverse strand, respectively). When hovering over a read, extra information is displayed. Some reads have lines over them. Try to zoom in in one of those lines to identify the reason for these.

    Because the number of reads over a region can be quite large, the IGV browser by default only allows to see the reads that fall into a small window. This behaviour can, in principle, be changed in the preferences panel.

# Inspection of the SAM format

The binary BAM file can be converted into a simple (but large!) text file, which is called a SAM (Sequence Alignment Map) file

:pencil2: ***Hands on!***

1. Run the tool `BAM-to-SAM` in Galaxy using the BAM file that was created  after mapping
2. Click 'include header in output'
3. Click on the view icon ('eye')

    The first part of the file, the header, contains the chromosome names and lengths. After the header, the location and other information of each read found in the FASTQ file is given.

# Correlation between samples

# Conclusion

Conclusion about the technical key points. And then relation between the technics and the biological question to end with a global view.

:grey_exclamation: ***Key Points***

- *Simple sentence to sum up the first key point of the tutorial (Take home message)*
- *Second key point*
- *Third key point*
- *...*

# :clap: Thank you
