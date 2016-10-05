Identification of the binding sites of the Estrogen receptor with ChIP-sequencing
=============

:grey_question: ***Questions***

- *How can we analyze raw ChIP-seq data?*
- *What are the binding sites of the Estrogen receptor?*
- *Third question*
- *...*

:dart: ***Objectives***

- *Analysis of ChIP-seq data*
- *Mastering the ChIP-seq data analysis workflow*
- *Identification of the binding sites of the Estrogen receptor*
- *Third objective*
- *...*

:heavy_check_mark: ***Requirements***

- *Galaxy introduction*
- *NGS-QC*
- *NGS-mapping*
- *IGV introduction*

:hourglass: ***Time estimation*** *1d/3h/6h*

# Introduction

This exercise uses the dataset from the Nature publication by [Ross-Inness et al., 2012](http://www.ncbi.nlm.nih.gov/pubmed/22217937).
The goal was to identify the binding sites of the Estrogen receptor, a transcription factor known to be associated with different types of breast cancer.
To this end, ChIP-seq was performed in breast cancer cells from 4 patients of different outcomes (good and poor).
For each ChIP-seq experiment there is a matching technical control, i.e., there are 8 samples in total, half of which are the so-called 'input' samples for which the same treatment as the ChIP-seq samples was done except for the immunoprecipitation step.
The input files are used to identify sequencing bias like open chromatin or GC bias.

Because of the long processing time for the large original files, we have selected small samples for practice and provide already processed data for subsequent steps.

# Pretreatments

As for each NGS data analysis, ChIP-seq data must be [quality controlled](../../NGS-QC/slides/dive_into_qc.html) and mapped on reference genomes.

:pencil2: ***Hands on!***

1. Import the datasets
2. Control quality of sequences using `FastQC` and `TrimGalore` as explained in [NGS-QC tutorial](./../NGS-QC/tutorials/dive_into_qc.md)

    We recommend you to trim low quality bases with a cut-off of 15 and adapter sequences with an overlap of at least 3 and also to clip a possibly trailing adapter sequence

3. Map the reads on human genome (version GRCh38) using `bowtie2`

    How many reads where mapped?

4. Visualize the reads in IGV and zoom in the start of chromosome 11 (`chr11:1,562,200-1,591,483` for example)

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
