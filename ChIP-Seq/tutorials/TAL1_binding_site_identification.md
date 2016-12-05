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

# Analysis

### Step 1: Quality control

As for any NGS data analysis, ChIP-seq data must be [quality controlled](../../NGS-QC/slides/dive_into_qc.html) before being aligned to a reference genome.

:pencil2: ***Hands on!***

1. Create and name a new history for this tutorial

1. Import the ChIP-seq raw data from Zenodo

    - Import all files from Zenodo. Load them into Galaxy by right-clicking → copy link location and paste the link in Galaxy → Get Data → Upload File from your computer → paste/fetch data → Start.

2. Examine the data by clicking on the 'eye' icon. 

    - Can you spot where the DNA sequence is stored?
    - Can you guess what the other entries mean?        

2. Run the tool `FastQC` on each FASTQ file to determine the quality of the raw data. An explanation of the results can be found on the [FastQC web page](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

### Step 2: Trimming and clipping reads

It is often necessary to trim a sequenced read, for example, to remove bases sequenced with high uncertainty (= 'low quality bases'). In addition, artificial adaptor sequences used in the library preparation protocol need to be removed before attempting to align the reads to a reference genome. More explanation of quality trimming can be found in the [NGS-QC tutorial](./../NGS-QC/tutorials/dive_into_qc.md).

:pencil2: ***Hands on!***

1. Run the tool `TrimGalore` on each FASTQ file to trim low-quality bases and remove any adaptor sequences. Explore the full parameter list for `TrimGalore` in the Tool Form and tell `TrimGalore` to:
    - **trim low-quality bases** using a cut-off of 15
    - **clip possible adapter sequence** using the default sequence that appears (this is the generic Illumina adaptor) with an overlap of at least 3.

2. Rerun the tool `FastQC` on each trimmed/clipped FASTQ file to determine whether low-quality and adaptor sequences were correctly removed.

**Note**: If your FASTQ files cannot be selected, you might check whether their format is FASTQ with Sanger-scaled quality values (*fastqsanger*). You can edit the data type by clicking on the 'pencil' symbol.

### Step 3: Aligning reads to a reference

In order to figure where the sequenced DNA fragments originated from in the genome, the short reads must be aligned to the reference genome. This is equivalent to solving a jigsaw puzzles, but unfortunately, not all pieces are unique. In principle, you could do a BLAST analysis to figure out where the sequenced pieces fit best in the known genome. Aligning millions of short sequences this way may, however, take a couple of weeks.
Nowadays, there are many read alignment programs for shot-gun sequenced DNA, `bowtie2` being one of them.

:pencil2: ***Hands on!***

1. Run the tool `bwa` to map the single-end reads to the mouse genome version mm10.

    - By clicking on the resulting history entry, you can see some basic mapping statistics once the alignment is completed. 
    - How many reads were mapped?

### Step 4: Visualizing aligned reads in IGV

The read alignment step with `bowtie2` resulted in a compressed, binary file (BAM) that is not human-readable (it's like the zipped version of a text file). We will show you two ways to inspect the file, (1) by visualization using a Genome Browser and (2) by converting the binary format into its text file equivalent.

:pencil2: ***Hands on!***

1. Load the reads into the IGV browser by clicking the option 'display with IGV local'. To see the reads in IGV you will need to zoom in the start of chromosome 11.

    - Try to zoom in to position ##-## on chromosome 11 (e.g. `chr11:##-##`).

2. Notice that the reads have a _direction_ (i.e., they are mapped to the forward or reverse strand, respectively). When hovering over a read, extra information is displayed. Some reads have lines over them. Try to zoom in in one of those lines to identify the reason for these.

**Note**: Because the number of reads over a region can be quite large, the IGV browser by default only allows to see the reads that fall into a small window. This behaviour can, in principle, be changed in the preferences panel.

### Step 5: Inspection of SAM format

The binary BAM file can be converted into a simple (but large!) text file, which is called a SAM (Sequence Alignment Map) file

:pencil2: ***Hands on!***

1. Run the tool `BAM-to-SAM` in Galaxy using the BAM file that was created  after mapping

2. Click 'include header in output'

3. Click on the view icon ('eye')

    - The first part of the file - the header - contains the chromosome names and lengths. After the header, the location and other information of each read found in the FASTQ file is given.

### Step 6: Correlation between samples

TODO

# Conclusion

Conclusion about the technical key points. And then relation between the technics and the biological question to end with a global view.

:grey_exclamation: ***Key Points***

- *Simple sentence to sum up the first key point of the tutorial (Take home message)*
- *Second key point*
- *Third key point*
- *...*

# :clap: Thank you
