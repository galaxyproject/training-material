---
layout: tutorial_hands_on

title: "Pre-processing of 10X Single-Cell ATAC-seq Datasets"
subtopic: end-to-end
priority: 3
redirect_from:
  - /topics/transcriptomics/tutorials/satac-preprocessing-tenx/tutorial
zenodo_link: "https://zenodo.org/record/3457880"
tags:
  - single-cell
  - 10x
questions:
  - What is 10X?
  - What are single-cell ATAC fragments and MTX files?
  - What is an HDF5 file, and why is it important?
objectives:
  - Demultiplex single-cell FASTQ data from 10X Genomics
  - Learn about transparent matrix formats
  - Create a count matrix starting from scATAC-seq FASTQ files
time_estimation: 1h
key_points:
  - Barcode FASTQ Reads are used to parse cDNA sequencing Reads
  - Creating a scATAC-seq count matrix requires positions of open chromatin regions
requirements:
  -
    type: "internal"
    topic_name: single-cell
    tutorials:
      - scrna-preprocessing-tenx

follow_up_training:
  -
    type: "internal"
    topic_name: single-cell
    tutorials:
      - scrna-scanpy-pbmc3k

contributors:
  - pavanvidem

gitter: Galaxy-Training-Network/galaxy-single-cell

---



# Introduction

Similar to bulk ATAC-Seq, single-cell ATAC-Seq (scATAC-seq) leverages the hyperactive Tn5 Transposase to profile open chromatin regions but at single-cell resolution. Thus helps in understanding cell type-specific chromatin accessibility from a heterogeneous cell population.

### Era of 10x Genomics

10x genomics has provided not only a cost-effective high-throughput solution to understanding sample heterogeneity at the individual cell level but has defined the standards of the field that many downstream analysis packages are now scrambling to accommodate. The gain in resolution reduces the granularity and noise issues that plagued the field of single-cell omics not long ago, where now individual clusters are much easier to decipher due to the added stability added by this gain in information.

### Library Preparation

![Library Preparation]({% link topics/single-cell/images/scrna-pre-processing/tenx_libprep_scatac.png %} "An overview of the 10x single-nuclei ATAC-seq library preparation")

The 10X barcoded gel beads consist of a pool of barcodes which are used to separately index each cell. The individual gel barcodes are delivered to each cell via flow cytometry, where each cell is fed single-file along a liquid tube and tagged with a 10X gel bead. The cells are then isolated from one another within thousands of nanoliter droplets, where each droplet is described by a unique 10x barcode that all reads in that droplet are associated with. The oil is then removed and all (now barcoded) DNA reads are pooled together to be sequenced.

More information can be found in the [reagent kit documentation](https://www.10xgenomics.com/support/single-cell-atac/documentation/steps/library-prep/chromium-single-cell-atac-reagent-kits-user-guide-v-2-chemistry).


### 10x Chemistries

There are two main reagent kits used during the library preparation. Below we can see the layout of the primers used in both chemistries. We can ignore most of these as they are not relevant, namely: the P5 and P7 Illumina primers are used in the Illumina [bridge amplification](https://en.wikipedia.org/wiki/Illumina_dye_sequencing#Bridge_amplification) process; the Sample Index is an 8bp primer which is related to the Chromium system that balances nucleotide bias and ensures that there is no sample overlap during the multiplexed sequencing. The primer of interest to us is the Cell Barcode (CB). Unlike scRNA-seq, Unique Molecular Identifiers (UMIs) are not used in scATAC-seq.

![chem]({% link topics/single-cell/images/scrna-pre-processing/tenx_primers_scatac.svg %} "10x scATAC-Seq sequencing library")

There are two major chemistries (v1 and v2) used in 10x scATAC-seq library construction. Technically, v2 chemistry is more sensitive and has better data yield. But the final sequencing library is the same. Hence, same the data analysis applies to both chemistries. Both chemistries generate the 4 FASTQ files containing:
* 8bp sample index sequence (I1)
* 50bp forward read (R1)
* 16bp barcode sequence (R2)
* 49bp reverse read (R3)

The major differences to 10x scRNA-seq sequencing libraries are: (a) barcode sequences get a separate FASTQ file and (b) there are no UMIs. Note that as opposed to the traditional FASTQ naming, R3 contains reverse mate-pair and R2 contains barcode sequence. 


# Analysis Strategy

> <agenda-title></agenda-title>
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

The tutorial is structured into three parts. The first part deals with the preprocessing of FASTQ files, quality control and mapping. Then in the second part, we identify the open chromatin regions and finally create a count matrix using the cell barcode and open chromatin regions information. 

Here we will use the 5k PBMCs from a Healthy Donor (v1 chemistry) from 10x genomics, consisting of 5k Peripheral blood mononuclear cells (PBMCs) extracted from a healthy donor, where PBMCs are primary cells. Analyzing these full datasets requires several hours of computation time. Hence, we subsampled the data so that the tutorial can be finished in a reasonable amount of time with sensible results. In the sub-sampled data, we have the reads that were mostly generated from chromosome Y. The original data has FASTQ files from two sequencing lanes *L001* and *L002*. We merged the sequences from both lanes. We also ignored the sample index FASTQ file as it is not useful for this analysis. Finally, we are left with three FASTQ files: **R1** (forward reads), **R2** (barcodes) and **R3** (reverse reads).

  * atac\_v1\_pbmc\_5k\_S1\_**R1**\_001.chrY.fastq.gz
  * atac\_v1\_pbmc\_5k\_S1\_**R2**\_001.chrY.fastq.gz
  * atac\_v1\_pbmc\_5k\_S1\_**R3**\_001.chrY.fastq.gz

These files are provided in the [Zenodo](https://zenodo.org/record/7844000) data repository. If you want to analyze the full datasets, the original FASTQ files are available at [10x genomics website](https://www.10xgenomics.com/resources/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-1-standard-1-2-0) 

# FASTQ processing and mapping

## Data upload and organization

We import all three FASTQ files into a history.

> <hands-on-title>Data upload and organization</hands-on-title>
>
> 1. Create a new history and rename it (e.g. scATAC-seq 10X preprocessing tutorial)
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 1. Import the sub-sampled FASTQ data from [`Zenodo`](https://zenodo.org/record/7844000) or from the data library (ask your instructor)
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
>    ```
>    https://zenodo.org/record/7844000/files/atac_v1_pbmc_5k_S1_R1_001.chrY.fastq.gz
>    https://zenodo.org/record/7844000/files/atac_v1_pbmc_5k_S1_R2_001.chrY.fastq.gz
>    https://zenodo.org/record/7844000/files/atac_v1_pbmc_5k_S1_R3_001.chrY.fastq.gz
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
{: .hands_on}

## Attach cell barcodes to FASTQ headers

First, we will attach the barcodes from the barcodes FASTQ file to the ids of forward and reverse reads. Hence, we will be left with forward and reverse reads only after this step.

> <comment-title></comment-title>
>
> {% tool [Sinto barcode](toolshed.g2.bx.psu.edu/repos/iuc/sinto/sinto_barcode/0.9.0) %}  consumes the barcode sequence from the `atac_v1_pbmc_5k_S1_R2_001.chrY.fastq.gz` file and uses `:` delimiter to prepend it to both the forward and reverse read ids. The tool expects the same number of reads in all three FASTQ files and should be present in the same order.
{: .comment}

> <question-title></question-title>
>
> 1. Which file contains the barcode sequence?
> 1. What is the length of the barcode?
>
> > <solution-title></solution-title>
> > 1. In this case, R2 contains the cell barcodes
> > 1. 16bp
> >
> {: .solution}
{: .question}

> <hands-on-title></hands-on-title>
>
> {% tool [Sinto barcode](toolshed.g2.bx.psu.edu/repos/iuc/sinto/sinto_barcode/0.9.0) %}  with the following parameters:
>    - *"FASTQ file containing cell barcode sequences"*: `atac_v1_pbmc_5k_S1_R2_001.chrY.fastq.gz`
>    - *"Single or Paired-end data"*: `Paired`
>    - *"Forward reads FASTQ file"*: `atac_v1_pbmc_5k_S1_R1_001.chrY.fastq.gz`
>    - *"Reverse reads FASTQ file"*: `atac_v1_pbmc_5k_S1_R3_001.chrY.fastq.gz`
>    - *"Number of bases to extract from barcode-containing FASTQ"*: `16`
>
>    > <comment-title>Sinto barcode output</comment-title>
>    >
>    > After a successful run, the barcode sequences from the barcode FASTQ file are prepended to the read names. We can inspect the output FASTQ files by clicking on the {% icon galaxy-eye %} symbol of the *barcoded read 1* file. 
>    {: .comment}
>
{: .hands_on}


Now let's proceed with the next steps QC, mapping and peak calling as we do with the bulk ATAC-Seq data. For a thorough description of bulk ATAC-Seq analysis, please follow the bulk [ATAC-Seq data analysis tutorial]({% link topics/epigenetics/tutorials/atac-seq/tutorial.md %}). Here we use a slightly different workflow. Here we map reads using **BWA-MEM** instead of **Bowtie2**.

## FASTQ Quality control
First things first. Let's do a basic FASTQ quality control using **FastQC**

> <hands-on-title></hands-on-title>
>
> 1. {% tool [FastQC](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.73+galaxy0) %} with the following parameters:
>       - *"Short read data from your current history"*: Use {% icon param-files %} **Multiple datasets** to choose both `barcoded read 1` and `barcoded read 2` (outputs of **Sinto barcode** {% icon tool %}).
> 2. Inspect the web page output of **FastQC** {% icon tool %} for the `barcoded read 1` sample.
>
> > <question-title></question-title>
> >
> > 1. How many reads are in the FASTQ?
> > 2. Which sections have a warning?
> >
> > > <solution-title></solution-title>
> > >
> > > 1. There are 685686 reads.
> > > 2. The 2 steps below have warnings:
> > >
> > >    1. **Per base sequence content**
> > >
> > >       It is well known that the Tn5 has a strong sequence bias at the insertion site. You can read more about it in {% cite Green2012 %}.
> > >
> > >    2. **Sequence Duplication Levels**
> > >
> > >       The read library quite often has PCR duplicates that are introduced
> > >       simply by the PCR itself. We will remove these duplicates later on.
> > >
> > {: .solution}
> >
>    {: .question}
{: .hands_on}

> <comment-title></comment-title>
>
> There is only a little amount of Nextera Transposase Sequence in the reads (see Adapter Content section of **FastQC** output). In this case, we skip the adapter trimming step and proceed to the mapping step. If your data has a significant amount of Nextera Transposase Sequences, then proceed to adapter trimming using either **cutadapt** or **Trim Galore!**. **Trim Galore** can automatically detect the adapter sequence. If you use **Cutadapt**, follow the bulk ATAC-Seq tutorial [Trimming Reads]({% link topics/epigenetics/tutorials/atac-seq/tutorial.md#trimming-reads %}) step.
{: .comment}

## Mapping reads to a reference genome
Now we map the reads to a reference genome using {% tool [BWA-MEM](toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.17.2) %}.

> <hands-on-title>Mapping</hands-on-title>
>
> {% tool [Map with BWA-MEM](toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.17.2) %} with the following parameters:
>    - *"Using reference genome"*: `hg38`
>    - In *"Single or Paired-end reads"*: `Paired`
>      - In *"Select first set of reads "*: `barcoded read 1` (output of **Sinto barcode** {% icon tool %})
>      - In *"Select second set of reads "*: `barcoded read 2` (output of **Sinto barcode** {% icon tool %})
>    - *"Select analysis mode"*: `1.Simple Illumina mode`
>    - *"BAM sorting mode"*: `Sort by chromosomal coordinates`
{: .hands_on}

# Peak calling
In scRNA-seq we always have the standard set of genes (usually downloaded from public databases) to quantify the expression levels. For scATAC-seq, there is no such reference open chromatin regions because the regions of chromatin accessibility are tissue dependent. Hence, we first have to detect open chromatin regions. There are several ways of doing this. The most common ways are the following:
* Find a published bulk ATAC-seq data that is the closest to your scATAC-seq data and use the regions from the bulk ATAC-seq data
* Chumk the genome into equal-sized bins and use these bins as reference locations for quantification
* Combine the data from all the cells of the scATAC-seq data together and detect open chromatin regions from the data. Later use these detected regions for quantification

In this tutorial, we opt for the 3rd option. 

## Create scATAC-seq fragments file
An ATAC-seq fragment file is a BED file with Tn5 integration sites, the cell barcode associated with the fragment, and the frequency of the sequenced fragment. PCR duplicates are collapsed. Here we filter out the alignments with low mapping quality. 

> <comment-title>on Tn5 insertion</comment-title>
>
> When Tn5 cuts an accessible chromatin locus it inserts adapters separated by 9bp ({% cite Kia2017 %}):
> ![Nextera_Library_Construction]({% link topics/epigenetics/images/atac-seq/NexteraLibraryConstruction.jpg %} "Nextera Library Construction")
>
> This means, to have the read start site reflecting the center of where Tn5 bound, the reads on the positive strand should be shifted 4 bp to the right and reads on the negative strands should be shifted 5 bp to the left as in [Buenrostro et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3959825). **Genrich** can apply these shifts when ATAC-seq mode is selected. In most cases, we do not have 9bp resolution so we don't take it into account but if you are interested in the footprint, this is important.
{: .comment}

> <hands-on-title>Create scATAC fragments</hands-on-title>
> 1. {% tool [Sinto fragments](toolshed.g2.bx.psu.edu/repos/iuc/sinto/sinto_fragments/0.9.0) %} with the following parameters:
>    - *"Input BAM file"*: `mapped reads in BAM format` (output of **Map with BWA-MEM** {% icon tool%})`
>    - *"Minimum MAPQ required to retain fragment"*: `30`
>    - *"Regular expression used to extract cell barcode from read name"*: `[^:]*` (matches all characters up to the first colon)
>    - *"Number of bases to shift Tn5 insertion position by on the forward strand"*: `4`
>    - *"Number of bases to shift Tn5 insertion position by on the reverse strand"*: `-5`
>
> 1. {% tool [bedtools SortBED](toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_sortbed/2.30.0+galaxy2) %} with the following parameters:
>    - *"Sort the following BED/bedGraph/GFF/VCF/EncodePeak file *"*: `fragments BED` (output of **Sinto fragments** {% icon tool%})`
>    - *"Sort by"*: `chromosome, then by start position (asc)`
> 1. Rename the datasets `sorted fragments`
{: .hands_on}

## Call Peaks

We have now finished the data preprocessing. Next, to find regions corresponding to potential open chromatin regions, we want to identify regions where reads have piled up (peaks) greater than the background read coverage. The tools which are currently used are [Genrich](https://github.com/jsh58/Genrich) and [MACS2](https://github.com/taoliu/MACS). MACS2 is more widely used.
At this step, two approaches exist:

- The first one is to select only paired whose fragment length is below 100bp corresponding to nucleosome-free regions and to use a peak calling like you would do for a ChIP-seq, joining signal between mates. The disadvantage of this approach is that you can only use it if you have paired-end data and you will miss small open regions where only one Tn5 is bound.
- The second one chosen here is to use all reads to be more exhaustive. In this approach, it is very important to re-center the signal of each read on the 5' extremity (read start site) as this is where Tn5 cuts. Indeed, you want your peaks around the nucleosomes and not directly on the nucleosome:

![ATAC-Seq_reads_relative_to_nucleosomes]({% link topics/epigenetics/images/atac-seq/schemeWithLegend.jpg %} "Scheme of ATAC-Seq reads relative to nucleosomes")

If we only assess the coverage of the 5' extremity of the reads, the data would be too sparse and it would be impossible to call peaks. Thus, we will extend the start sites of the reads by 200bp (100bp in each direction) to assess coverage. We call peaks with MACS2. In order to get the coverage centered on the 5' extended 100bp on each side we will use `--shift -100` and `--extend 200`:
![MACS2_peak_centering_options]({% link topics/epigenetics/images/atac-seq/macs2Options.jpg %} "MACS2 options to get 100bp each side")

> <hands-on-title>Call peaks with MACS2</hands-on-title>
>
> 1. {% tool [MACS2 callpeak](toolshed.g2.bx.psu.edu/repos/iuc/macs2/macs2_callpeak/2.1.1.20160309.6) %} with the following parameters:
>    - *"Are you pooling Treatment Files?"*: `No`
>    - *"ChIP-Seq Treatment File"*: `sorted fragments`
>    - *"Do you have a Control File?"*: `No`
>    - *"Format of Input Files"*: `Single-end BED`
>    - *"Effective genome size"*: `H. sapiens (2.7e9)`
>    - *"Build Model"*: `Do not build the shifting model (--nomodel)`
>        - *"Set extension size"*: `200`
>        - *"Set shift size"*: `-100`. It needs to be - half the extension size to be centered on the 5'.
>    - *"Additional Outputs"*:
>        - Check `Peaks as tabular file (compatible with MultiQC)`
>        - Check `Peak summits`
>        - Check `Scores in bedGraph files`
>    - In *"Advanced Options"*:
>        - *"Composite broad regions"*: `No broad regions`
>            - *"Use a more sophisticated signal processing approach to find subpeak summits in each enriched peak region"*: `Yes`
>        - *"How many duplicate tags at the exact same location are allowed?"*: `all`
>
>    > <comment-title>Why keeping all duplicates is important</comment-title>
>    >
>    > We previously removed duplicates using **Sinto fragment** {% icon tool %} using paired-end information. If two pairs had identical R1 but different R2, we knew it was not a PCR duplicate. Because we converted the BAM to BED we lost the pair information. If we keep the default (removing duplicates) one of the 2 identical R1 would be filtered out as duplicate.
>    {: .comment}
>    > <comment-title>BED / encode narrowPeak format</comment-title>
>    > If you are not familiar with BED format or encode narrowPeak format, see the [BED Format](https://genome.ucsc.edu/FAQ/FAQformat.html)
>    {: .comment}
>
{: .hands_on}

# Count matrix creation

Now we have ATAC fragments with cell barcode information as well as the peaks to create a scATAC count matrix. In the context of ATAC-seq, the peaks are also often called features. But the term *feature* is a generic term that denotes a *peak* in the case of scATAC-seq data and a *gene* in case of scRNA-seq data. So the count matrix we build is similar to that of scRNA-seq data with the exception that here we have open chromatin regions instead of expressed genes.

For count matrix creation, we can use **Build count matrix** from **EpiScanpy** tool suite. It takes fragments file and the peaks in either BED format or directly a narrowPeak file of **MACS2** output.

> <hands-on-title>Build count matrix with EpiScanpy</hands-on-title>
>
> 1. {% tool [Build count matrix with EpiScanpy](toolshed.g2.bx.psu.edu/repos/iuc/episcanpy/episcanpy_build_matrix/0.3.2+galaxy0) %} with the following parameters:
>    - *"ATAC fragments file"*: `fragments BED` (output of **Sinto fragments** {% icon tool %})
>    - *"Features file"*: `narrow Peaks` (output of **MACS2** {% icon tool %})
>    - *"Number of bases to extend both sides of peaks"*: `0`
>    - In *"Select the chromosomes of the species you are considering"*: `Custom list of chromosomes`
>        - *"Enter comma seperated list of chromosome ids (without chr prefix)"*: `Y`
>
>    > <tip-title>Chromosome selection</tip-title>
>    >
>    > In this example data, our reads are from the *chrY*. Hence, we input `Y` for chromosome ids. For human or mouse genomes, directly select from the dropdown list. Make sure that all the canonical chromosomes are present in the fragments file.
>    {: .tip}
>
> > <question-title></question-title>
> >
> > How many barcodes and features are there?
> >
> > > <solution-title></solution-title>
> > >
> > > There are 30328 barcodes and 204 features
> > >
> > {: .solution}
> >
>    {: .question}
{: .hands_on}

# Conclusion


In this workflow, we have learned to quickly perform mapping, peak calling and quantification of scATAC-seq FASTQ data.

A full pipeline that produces both an AnnData and tabular file for inspection is provided [in this workflow](workflows/scatac_tenx.ga).

This tutorial is part of the https://singlecell.usegalaxy.eu portal ({% cite tekman2020single %}).