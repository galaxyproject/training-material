Identification of the binding sites of the T-cell acute lymphocytic leukemia protein 1 (TAL1) with ChIP-sequencing
=============

:grey_question: ***Questions***

- *How is raw chIP-seq data processed and analyzed?*
- *What are the binding sites of Tal1?*
- *Which genes are regulated by Tal1?*
- *...*

:dart: ***Objectives***

- *Inspect read quality with FastQC*
- *Perform read trimming with Trimmomatic*
- *Align trimmed reads with BWA*
- *Identify Tal1 "peaks" with MACS2 (and generate bedgraph files for visualization)*
- *Intersect Tal1 peaks from G1E and Megakaryocytes to identify unique and common Tal1 binding sites*
- *Identify unique/common Tal1 peaks occupying gene promoters*
- *Visually inspect Tal1 peaks with Trackster*

:heavy_check_mark: ***Requirements***

- *Galaxy introduction* (add link)
- *NGS-QC* (add link)
- *NGS-mapping* (add link)
- *Trackster* (add link)

:hourglass: ***Time estimation*** *3h*

# Introduction
 
This tutorial uses ChIP-seq datasets from a study published by [Wu et al., 2012](http://genome.cshlp.org/content/24/12/1945.full.pdf+html).
The goal of this study was to investigate "the dynamics of occupancy and the role in gene regulation of the transcription factor Tal1, a critical regulator of hematopoiesis, at multiple stages of hematopoietic differentiation."
To this end, chIP-seq was performed in the G1E cell line - a GATA-null immortalized cell line derived from targeted disruption of GATA-1 in mouse embryonic stem cells - and megakaryocytes.
This dataset (GEO Accession [GSE51338](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51338)) consists of biological replicate Tal1 ChIP-seq experiments and input control experiments for which the same treatment as the ChIP-seq samples was done except for the immunoprecipitation step.
Input control experiments are used to identify and remove sampling bias, for example open/accessible chromatin or GC bias.

Because of the long processing time for the large original files, we have downsampled the original raw data files to include only a subset of genomic loci.

INSERT TABLE WITH SAMPLE META-DATA

# Analysis

### Step 1: Quality control

As for any NGS data analysis, ChIP-seq data must be [quality controlled](../../NGS-QC/slides/dive_into_qc.html) before being aligned to a reference genome.

:pencil2: ***Hands on!***

1. Create and name a new history for this tutorial

2. Import the ChIP-seq raw data from Zenodo

    - Download all files from [Zenodo](link).
    - Import files directly into Galaxy by providing links: for each file right-click → copy link location and paste the link in Galaxy → Get Data → Upload File from your computer → paste/fetch data → set datatype to 'fastqsanger' → Start.

3. Examine the data by clicking on the 'eye' icon. 

    - This is what raw Illuimna sequencing reads look like.
    
4. Run the tool `FastQC` on each FASTQ file to determine the quality of the raw data. An explanation of the results can be found on the [FastQC web page](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

### Step 2: Trimming and clipping reads

It is often necessary to trim a sequenced read, for example, to remove bases sequenced with high uncertainty (= 'low quality bases'). In addition, artificial adaptor sequences used in the library preparation protocol need to be removed before attempting to align the reads to a reference genome. More explanation of quality trimming can be found in the [NGS-QC tutorial](./../NGS-QC/tutorials/dive_into_qc.md).

:pencil2: ***Hands on!***

1. Run the tool `Trimmomatic` on each FASTQ file to trim low-quality bases. Explore the full parameter list for `Trimmomatic` in the Tool Form and set `Trimmomatic` parameters to:
    - **Paired end data?** No
    - **Perform initial ILLUMINACLIP?** NO
    - **Select Trimmomatic operation to perform** Sliding window trimming (SLIDINGWINDOW)
    - **Number of bases to average across** 4
    - **Average quality required** 20
    
2. Rerun the tool `FastQC` on each trimmed/clipped FASTQ file to determine whether low-quality and adaptor sequences were correctly removed.

**Note**: If your FASTQ files cannot be selected, you might check whether their format is FASTQ with Sanger-scaled quality values (*fastqsanger*). If you didn't set this datatype when importing the data, you can edit the data type by clicking on the 'pencil' symbol.

### Step 3: Aligning reads to a reference genome

In order to figure where the sequenced DNA fragments originated from in the genome, the short reads must be aligned to the reference genome. This is equivalent to solving a jigsaw puzzle, but unfortunately, not all pieces are unique. In principle, you could do a BLAST analysis to figure out where the sequenced pieces fit best in the known genome. Aligning millions of short sequences this way may, however, take a couple of weeks.
Nowadays, there are many read alignment programs for shot-gun sequenced DNA, `BWA` being one of them.

:pencil2: ***Hands on!***

1. Run the tool `Map with BWA` to map the single-end reads to the mouse genome version mm10. Change the following parameters:
    - **Using reference genome** Mouse (mus musculus): mm10
    - **Single or Paired-end reads** Single

    - By clicking on the resulting history entry, you can see some basic mapping statistics once the alignment is completed. 
    - How many reads where mapped?

### Step 4: Determining Tal1 binding sites 

Now that 'BWA' has aligned the reads to the genome, we will use the tool 'MACS2' to identify regions of Tal1 occupancy, which are called "peaks". 
'MACS2' will perform two tasks: 1) identify regions of Tal1 occupancy ("peaks") and 2) generate bedgraph files for visual inspection of the data on a genome browser. 

:pencil2: ***Hands on!***

1. Run the tool 'MACS2 callpeak' with the aligned read files from the previous step as Treatment (Tal1) and Control (input). 
    - **Effective genome size** Mouse (2,150,570,000)
2. Rename your files to reflect the origin and contents.

### Step 5: Inspection of peaks and aligned data

It is critical to visualize your NGS data on a genome browser after alignemnt. Evaluation criteria will differ for the various NGS experiment types, but for chIP-seq data we want to ensure reads from a Treatment sample are enriched at "peaks" and do not localize non-specifically (like the Control condition).
'MACS2' generated a bedgraph and a bed file that we'll use to visualize read abundance and peaks, respectively, at regions 'MACS2' determined to be Tal1 peaks using the genome browser Trackster. We'll first need to tidy up the peak file before we send it to Trackster. We'll also import a gene annotation file so we can visualize aligned reads and Tal1 peaks relative to gene features and positions.

:pencil2: ***Hands on!***

1. Run the tool `Cut` on the peak file choosing columns 'c1,c2,c3,c4' and rename this file to reflect the origin and contents. 

2. Import the gene annotations file here (COMPLETE THIS MO!!!!)

3. Click 'Visualizations' on the page header and select 'New Track Browser'. Name this session something descriptive and choose 'mm10' as the 'Reference genome build (dbkey)'.  

3. Click 'Add Datasets to visualization' and select the history containing the data from this analysis. Select the bedgraph files and the peak files (that you renamed).  

4. Navigate to the Gata1 locus (chr16:92501466-92926074) to inspect the aligned reads and Tal1 peak calls. (Mo:This region should have Tal1 peaks in both cellular states)
    - What do you see?
    
### Step 6: Identifying unique and common Tal1 peaks between the G1E and megakaryocyte states

We've just processed chIP-seq data from two stages of hematopoiesis and have lists of Tal1 occupied sites (peaks) in both cellular states. Now lets identify Tal1 peaks that are shared between the two cellular states and also those that are specific to one cellular state.

:pencil2: ***Hands on!***

1.  


# Conclusion

In this exercise you imported raw Illumina sequencing data, evaluated the quality before and after you trimmed reads with low confidence scores, algined the trimmed reads, identified Tal1 peaks relative to the negative control (background), and visualized the aligned reads and Tal1 peaks relative to gene structures and positions. 

:grey_exclamation: ***Key Points***

- *Simple sentence to sum up the first key point of the tutorial (Take home message)*
- *Second key point*
- *Third key point*
- *...*

# :clap: Thank you
