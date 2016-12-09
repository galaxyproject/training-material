Identification of the binding sites of the T-cell acute lymphocytic leukemia protein 1 (TAL1) with ChIP-sequencing
=============

:grey_question: ***Questions***

- *How is raw ChIP-seq data processed and analyzed?*
- *What are the binding sites of Tal1?*
- *Which genes are regulated by Tal1?*

:dart: ***Objectives***

- *Inspect read quality with FastQC*
- *Perform read trimming with Trimmomatic*
- *Align trimmed reads with BWA*
- *Assess quality and reproducibility of experiments*
- *Identify Tal1 binding sites with MACS2*
- *Determine unique/common Tal1 binding sites from G1E and Megakaryocytes*
- *Identify unique/common Tal1 peaks occupying gene promoters*
- *Visually inspect Tal1 peaks with Trackster*

:heavy_check_mark: ***Requirements***

- *Galaxy introduction* ([link](https://github.com/nekrut/galaxy/wiki/Galaxy101-1))
- *NGS-QC* ([link](../../NGS-QC/slides/dive_into_qc.html))
- *NGS-mapping* ([link](../../NGS-mapping/slides/dive_into_mapping.md))
- *Trackster* ([link](https://wiki.galaxyproject.org/Learn/Visualization))

:hourglass: ***Time estimation*** *3h*

# Content

* [Introduction](#introduction)
* [Hands-on example of ChIP-seq analysis](#analysis)
  - [Step 1: Quality control](#step-1-quality-control)
  - [Step 2: Trimming/clipping reads](#step-2-trimming-and-clipping-reads)
  - [Step 3: Aligning reads to a genome](#step-3-aligning-reads-to-a-reference-genome)
  - [Step 4: Assessing correlation between samples](#step-4-assessing-correlation-between-samples)
  - [Step 5: Assessing IP strength](#step-5-assessing-ip-strength)
  - [Step 6: Determining Tal1 binding sites](#step-6-determining-tal1-binding-sites)
  - [Step 7: Inspection of Tal1 peaks](#step-7-inspection-of-peaks-and-aligned-data)
  - [Step 8: Identifying unique/common Tal1 peaks](#step-8-identifying-unique-and-common-tal1-peaks-between-states)
  - [Additional optional analyses](#additional-optional-analyses)
* [Concluding remarks](#conclusion)
* [Useful literature](#useful-literature)
  
# Introduction
 
This tutorial uses ChIP-seq datasets from a study published by [Wu et al., 2014](http://genome.cshlp.org/content/24/12/1945.full.pdf+html).
The goal of this study was to investigate "the dynamics of occupancy and the role in gene regulation of the transcription factor Tal1, a critical regulator of hematopoiesis, at multiple stages of hematopoietic differentiation."
To this end, ChIP-seq experiments were performed in multiple mouse cell types including G1E - a GATA-null immortalized cell line derived from targeted disruption of GATA-1 in mouse embryonic stem cells - and megakaryocytes.
This dataset (GEO Accession: [GSE51338](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51338)) consists of biological replicate Tal1 ChIP-seq and input control experiments.
Input control experiments are used to identify and remove sampling bias, for example open/accessible chromatin or GC bias.

Because of the long processing time for the large original files, we have downsampled the original raw data files to include only reads that align to chromosome 19 and a subset of interesting genomic loci identified in the Wu et al. (2014) publication.

**Table 1**: Metadata for ChIP-seq experiments in this tutorial. SE: single-end.

| Cellular state | Datatype | ChIP Ab | Replicate | SRA Accession | Library type | Read length | Stranded? | Data size (MB) |
|---|---|:-:|:-:|---|:-:|:-:|:-:|---|
| G1E | ChIP-seq | input | 1 | SRR507859 | SE | 50 | No | 35.8 |
| G1E | ChIP-seq | input | 2 | SRR507860 | SE | 50 | No | 427.1 |
| G1E | ChIP-seq | Tal1 | 1 | SRR492444 | SE | 50 | No | 32.3 |
| G1E | ChIP-seq | Tal1 | 2 | SRR492445 | SE | 50 | No | 62.7 |
| Megakaryocyte | ChIP-seq | input | 1	| SRR492453 | SE | 50 | No | 57.2 |
| Megakaryocyte | ChIP-seq | input | 2 | SRR492454 | SE | 50 | No | 403.8 |
| Megakaryocyte | ChIP-seq | Tal1 | 1 | SRR549006 | SE | 50 | No | 340.3 |
| Megakaryocyte | ChIP-seq | Tal1 | 2 | SRR549007 | SE | 50 | No | 356.9 |

# Analysis

### Step 1: Quality control

As for any NGS data analysis, ChIP-seq data must be quality controlled before being aligned to a reference genome.

:pencil2: ***Hands on!***

1. Create and name a new history for this tutorial.

2. Import the ChIP-seq raw data (\*.fastqsanger) from Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.197100.svg)](https://doi.org/10.5281/zenodo.197100)

    - In Zenodo: Right-click a FASTQ filename → Copy Link Address
    - In Galaxy: Get Data → Upload File from your computer → Paste/Fetch data
    - Paste the copied link into the dialog box and set the datatype to 'fastqsanger'
    - Click Start
    - Repeat for all FASTQ files

    ![upload](../images/upload_data_page.png)
    <figcaption><b>Figure X:</b> Data can be imported directly with links.</figcaption>
    
    ![data](../images/data_uploaded.png)
    <figcaption><b>Figure X:</b> Imported datasets will appear in the history panel.</figcaption>
    
3. In Galaxy, examine the data in a FASTQ file by clicking on the eye icon.

    | :grey_question: Questions |
    |:---|
    | <ul><li>What are four key features of a FASTQ file?</li><li>What is the main difference between a FASTQ and a FASTA file?</li></ul> |
    
4. Run the tool `FastQC` on each FASTQ file to assess the quality of the raw data. An explanation of the results can be found on the [FastQC web page](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

    **Tip**: You can run this tool - and many other tools - on all the FASTQ files at once! To do this, first select the "Multiple datasets" icon (two stacked pages) under the "Input FASTQ file" heading in the `FASTQC` Tool Form, then shift+click to select multiple FASTQ files.

    ![fastqbefore](../images/fastqc_before.png)
    <figcaption><b>Figure X:</b> Sequence quality per base generated by FastQC <b>before</b> end trimming.</figcaption>
    
    | :grey_question: Questions |
    |:---|
    | <ul><li>What does the y-axis represent in Figure X?</li><li>Why is the quality score decreasing across the length of the reads?</li></ul> |

### Step 2: Trimming and clipping reads

It is often necessary to trim a sequenced read to remove bases sequenced with high uncertainty (*i.e.* low-quality bases). In addition, artificial adaptor sequences used in library preparation protocols need to be removed before attempting to align the reads to a reference genome.

:pencil2: ***Hands on!***

1. Run the tool `Trimmomatic` on each FASTQ file to trim low-quality bases (remember how to run tools on all files at once?). Explore the full parameter list for `Trimmomatic` in the Tool Form and set the following `Trimmomatic` parameters:

    - **Paired end data?**: No
    - **Perform initial ILLUMINACLIP?**: No
    - **Select Trimmomatic operation to perform**: Sliding window trimming (SLIDINGWINDOW)
    - **Number of bases to average across**: 4
    - **Average quality required**: 20
    
    **Tip**: If the FASTQ files cannot be selected, check whether their format is FASTQ with Sanger-scaled quality values (*fastqsanger*). If not, you can edit the data type by clicking on the pencil symbol next to a file in the history, clicking the 'Datatype' tab, and choosing *fastqsanger* as the 'New Type'.

2. Rerun the tool `FastQC` on each trimmed/clipped FASTQ file to determine whether low-quality and adaptor sequences were correctly removed.

    | :grey_question: Questions |
    |:---|
    | <ul><li>How did the range of read lengths change after trimming/clipping?</li><li>What do you think could account for the enriched k-mers (**Kmer Content** heading in `FASTQC` output) observed in the Megakaryocytes Tal1 R2 ChIP-seq experiment?</li></ul> |
    
    ![fastqafter](../images/fastqc_after.png)
    <figcaption><b>Figure X:</b> Sequence quality per base generated by FastQC <b>after</b> end trimming.</figcaption>
    
### Step 3: Aligning reads to a reference genome

To determine where DNA fragments originated from in the genome, the sequenced reads must be aligned to a reference genome. This is equivalent to solving a jigsaw puzzle, but unfortunately, not all pieces are unique. In principle, you could do a BLAST analysis to figure out where the sequenced pieces fit best in the known genome. Aligning millions of short sequences this way, however, can take a couple of weeks.
Nowadays, there are many read alignment programs for sequenced DNA, `BWA` being one of them. You can read more about the BWA algorithm and tool [here](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btp324).

:pencil2: ***Hands on!***

1. Run the tool `Map with BWA` to map the trimmed/clipped reads to the mouse genome. Set the following `Map with BWA` parameters:

    - **Using reference genome**: Mouse (mus musculus) mm10
    - **Single or Paired-end reads**: Single

    ![bwa](../images/BWA_tool_form.png)
    <figcaption><b>Figure X:</b> Select the trimmed read files to be aligned.</figcaption>
    
2. After `BWA` finishes, rename your files to reflect the origin and contents. Then, click on a file produced by running `BWA`.

    | :grey_question: Questions |
    |:---|
    | <ul><li>What datatype is the `BWA` output file?</li><li>How many reads were mapped from each file?</li></ul> |

3. Run the tool `IdxStats` and look at the output (poke it in the eye!).

    | :grey_question: Questions |
    |:---|
    | <ul><li>What does each column in the output represent (tip: look at the Tool Form)?</li><li>How many reads were mapped to chromosome 19 in each experiment?</li><li>If the mouse genome has 21 pairs of chromosomes, what are the other reference chromosomes (e.g. chr1_GL456210_random)?</li></ul> |

### Step 4: Assessing correlation between samples

To assess the similarity between the replicates sequencing datasets, it is a common technique to calculate the correlation of read counts for the different samples.
We expect that the replicates of the ChIP-seq experiments should be clustered more closely to each other than the replicates of the input samples.
We will be use tools from the package `deepTools` for the next few steps. More information on `deepTools` can be found [here](http://deeptools.readthedocs.io/en/latest/content/list_of_tools.html).

1. Run the tool `multiBamSummary` from the `deepTools` package. This tool splits the reference genome into bins of equal size and counts the number of reads in each bin from each sample. We set a small **Bin size** here because we are working with a subset of reads that align to only a fraction of the genome. Set the following `multiBamSummary` parameters:

    - Select all of the aligned BAM files
    - **Bin size in bp**: 1000

    ![multiBamSummary](../images/multiBamSummary_1000bin.png)
    <figcaption><b>Figure X:</b> Select all of the aligned BAM files and change Bin size.</figcaption>

2. Run the tool `plotCorrelation` from the `deepTools` package to visualize the results from the previous step. Feel free to try different parameters. To start, set the following `plotCorrelation` parameters:

    - **Correlation method**: Pearson
    - **Plotting type**: Heatmap
    - **Plot the correlation value**: Yes
    - **Skip zeros**: Yes
    - **Remove regions with very large counts**: Yes

    ![corr](../images/plotCorrelation.png)
    <figcaption><b>Figure X:</b> Select the newly generation correlation matrix file from the previous step.</figcaption>

| :grey_question: Questions |
|:---|
| <ul><li>Why do we want to skip zeros in `plotCorrelation`?</li><li>What happens if the Spearman correlation method is used instead of the Pearson method?</li><li>What does the output of making a Scatterplot instead of a Heatmap look like?</li></ul> |

![heatmap](../images/plotCorrelation_heatmap_pearson.png)
<figcaption><b>Figure X:</b> Heatmap of correlation matrix generated by `plotCorrelation`.</figcaption>

For additional informaton on how to interpret `plotCorrelation` plots, read the information [here](http://deeptools.readthedocs.io/en/latest/content/tools/plotCorrelation.html#background)

### Step 5: Assessing IP strength

We will now evaluate the quality of the immuno-precipitation step in the ChIP-seq protocol.

1. Run the tool `plotFingerprint` from the `deepTools` package.

    - Select all of the aligned BAM files
    - **Show advanced options**: yes
    - **Bin size in bases**: 50
    - **Skip zeros**: Yes

    ![fingerprint1](../images/plotFingerprint1.png)
    ![fingerprint2](../images/plotFingerprint2.png)
    <figcaption><b>Figure X:</b> Select all of the aligned BAM file to assess IP strength.</figcaption>

2. View the output image.

    | :grey_question: Questions |
    |:---|
    | <ul><li>What does this graph represent?</li><li>How do input datasets differ from IP datasets?</li><li>What do you think about the quality of the IP for this experiment?</li></ul> |
    
For additional informaton on how to interpret `plotFingerprint` plots, read the information [here](http://deeptools.readthedocs.io/en/latest/content/tools/plotFingerprint.html#background)

### Step 6: Determining Tal1 binding sites 

Now that `BWA` has aligned the reads to the genome, we will use the tool `MACS2` to identify regions of Tal1 occupancy, which are called "peaks". Peaks are determined from pileups of sequenced reads across the genome that correspond to where Tal1 binds.
`MACS2` will perform two tasks: 1) identify regions of Tal1 occupancy (peaks) and 2) generate bedgraph files for visual inspection of the data on a genome browser. More information about `MACS2` can be found [here](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137).

:pencil2: ***Hands on!***

1. Run the tool `MACS2 callpeak` with the aligned read files from the previous step as Treatment (Tal1) and Control (input). 

    - Select replicate **ChIP-Seq Treatment Files** for one cell type
    - Select replicate **ChIP-Seq Control Files** for the same cell type

    ![macs2](../images/MACS2_tool_form.png)
    <figcaption><b>Figure X:</b> Select the appropriate control and treatment files.</figcaption>
    
2. After `MACS2 callpeak` finishes, rename your files to reflect the origin and contents.

### Step 7: Inspection of peaks and aligned data

It is critical to visualize your NGS data on a genome browser after alignemnt. Evaluation criteria will differ for the various NGS experiment types, but for ChIP-seq data we want to ensure reads from a Treatment sample are enriched at "peaks" and do not localize non-specifically (like the Control condition).
`MACS2` generated a bedgraph and a BED file that we'll use to visualize read abundance and peaks, respectively, at regions `MACS2` determined to be Tal1 peaks using Galaxy's in-house genome browser, Trackster. We'll first need to tidy up the peak file before we send it to Trackster. We'll also import a gene annotation file so we can visualize aligned reads and Tal1 peaks relative to gene features and positions.

:pencil2: ***Hands on!***

1. Run the tool `Cut` on the peak file choosing columns 'c1,c2,c3,c4' and rename this file to reflect the origin and contents. 

2. Import the gene annotations file from [Zenodo] [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.197100.svg)](https://doi.org/10.5281/zenodo.197100)

3. Click 'Visualizations' on the page header and select 'New Track Browser'.

    ![vizbutton](../images/Highlighting_viz_button.png)
    <figcaption><b>Figure X:</b> Trackster can be accessed from the Visualizations button at the top of the screen.</figcaption>
    
4. Give this session a descriptive name and choose 'mm10' as the 'Reference genome build (dbkey)'.  
    
    ![tracksterdb](../images/Trackster_session_and_db.png)
    <figcaption><b>Figure X:</b> Session name and assigning reference genome.</figcaption>

5. Click 'Add Datasets to visualization' and select the history containing the data from this analysis. Select the bedgraph files and the peak files that you renamed.

    ![tracksteradd](../images/Trackster_add_datasets.png)
    <figcaption><b>Figure X:</b> Load your data into Trackster with the 'Add Datasets to visualization' feature. </figcaption>

    ![tracksterselect](../images/Trackster_selecting_datasets.png)
    <figcaption><b>Figure X:</b> Select data from your histories to view in Trackster. </figcaption>

6. Navigate to the Runx1 locus (chr16:92501466-92926074) to inspect the aligned reads and Tal1 peak calls. 

    | :grey_question: Questions |
    |:---|
    | <ul><li>What do you see?</li><li>What gene(s) other than Runx1 could be regulated by Tal1?</li></ul> |

    ![runx1](../images/Trackster_Runx1_locus.png)
    <figcaption><b>Figure X:</b> The Runx1 locus.</figcaption>

### Step 8: Identifying unique and common Tal1 peaks between states

We've just processed ChIP-seq data from two stages of hematopoiesis and have lists of Tal1 occupied sites (peaks) in both cellular states. Now let's identify Tal1 peaks that are shared between the two cellular states and also those that are specific to one cellular state.

:pencil2: ***Hands on!***

1. Run the tool `bedIntersect` to find peaks that exist both in G1E and megakaryocytes. 

  - Select the "Tal1 G1E peaks" and "Tal1 Mega peaks" files as the inputs. Running this tool with the default settings will return overlapping peaks of both files. 
  
2. Run the tool `bedIntersect` to find peaks that exist only in G1E.

  - Select "Tal1 G1E peaks" as the first input and "Tal1 Mega peaks" as the second input file.
  - **Report only those alignments that \*\*do not\*\* overlap the BED file**: Yes
  
3. Run the tool `bedIntersect` to find peaks that exist only in megakaryocytes.

  - Select "Tal1 Mega peaks" as the first input and "Tal1 G1E peaks" as the second input file.
  - **Report only those alignments that \*\*do not\*\* overlap the BED file**: Yes  
  
4. Re-name the three files we generated to reflect their contents. 

    | :grey_question: Questions |
    |:---|
    | <ul><li>How many Tal1 binding sites are common to both G1E cells and megakaryocytes?</li><li>How many are unique to G1E cells?</li><li>How many are unique to megakaryocytes?</li></ul> |

# Additional optional analyses

### Assessing GC bias

A common problem of PCR-based protocols is the observation that GC-rich regions tend to be amplified more readily than GC-poor regions.
We will now check whether the samples have more reads from regions of the genome with high GC.
 
1. Run the tool `computeGCbias` from the `deepTools` package.

    - First, select an aligned BAM files for an input dataset
    - **Reference genome**: locally cached
    - **Using reference genome**: mm10
    - **Effective genome size**: user specified
    - **Effective genome size**: 10000000
    - **Fragment length used for the sequencing**: 50
    
    ![gcbias](../images/computeGCBias.png)
    <figcaption><b>Figure X:</b> Select a single aligned BAM file to check GC bias for that dataset.</figcaption>

    | :grey_question: Questions |
    |:---|
    | <ul><li>Why would we worry more about checking for GC bias in an input file?</li><li>Does this dataset have a GC bias?</li></ul> |

2. Explore the tool `correctGCbias` from the `deepTools` package.

    | :grey_question: Questions |
    |:---|
    | <ul><li>What does this tool do?</li><li>What is the output of this tool?</li><li>What are some caveats to be aware of if using the output of this tool in downstream analyses?</li></ul> |

For additional informaton on how to interpret `computeGCbias` plots, read the information [here](http://deeptools.readthedocs.io/en/latest/content/tools/computeGCBias.html#background)


# Conclusion

In this exercise you imported raw Illumina sequencing data, evaluated the quality before and after you trimmed reads with low confidence scores, algined the trimmed reads, identified Tal1 peaks relative to the negative control (background), and visualized the aligned reads and Tal1 peaks relative to gene structures and positions. Additional, you assessed the "goodness" of the experiments by looking at metrics such as GC bias and IP enrichment. 

:grey_exclamation: ***Key Points***

- *Sophisticated analysis of ChIP-seq data is possible using tools hosted by Galaxy.*
- *Genomic dataset analyses require multiple methods of quality assessment to ensure that the data are appropriate for answering the biology question of interest.*
- *By using the sharable and transparent Galaxy platform, data analyses can easily be shared and reproduced.*

## :clap: Congratulations on successfully completing this tutorial!

# Useful literature

Useful information regarding ChIP-seq analysis in general, descriptions and paper references for the tools used in this tutorial, and literature for ChIP-seq-related analysis techniques and interpretations can be found [here](../README.md).
