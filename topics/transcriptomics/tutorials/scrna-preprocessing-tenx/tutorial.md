---
layout: tutorial_hands_on
title: "Pre-processing of 10X Single-Cell RNA Datasets"
zenodo_link: "https://zenodo.org/record/3383115"
tags:
  - single-cell
  - 10x
questions:
  - "What is 10X?"
  - "What is STARSolo and what is Cell Ranger?"
  - "What are BCL and MTX files?"
  - "What is an HDF5 file, and why is it important?"
objectives:
  - "Demultiplex single-cell FASTQ data from 10X Genomics"
  - "Learn about transparent matrix formats"
  - "Understand the importance of high and low quality cells"
time_estimation: "1h"
key_points:
  - "Filtering a raw matrix into a high quality one for downstream analysis"
  - "Removing unwanted barcodes"
requirements:
  -
    type: "internal"
    topic_name: transcriptomics
    tutorials:
      - scrna-preprocessing
#
# follow_up_training:
#   -
#     type: "internal"
#     topic_name: transcriptomics

contributors:
  - mtekman

---



# Introduction
{:.no_toc}

* Split-out single-cell into it's own intro section
* (optional) Check out the original pre-processing tutorial.

## Era of 10x Genomics

10x genomics has provided not only a cost-effective high-throughput solution to understanding sample heterogeneity at the individual cell level, but has defined the standards of the field that many downstream analysis packages are now scrambling to accommodate.

![clusters]({{ site.baseurl }}{% link topics/transcriptomics/images/tenx_clusters_intro.png %} "From less than 1K to over 10K with 10x genomics: Analyses of two separate scRNA datasets using the (left) CelSEQ2 protocol, and the (right) 10x Chromium system.")

The gain in resolution reduces the granularity and noise issues that plagued the field of scRNA-seq not long ago, where now individual clusters are much easier to decipher due to the added stability added by this gain in information.

## Library Preparation

![Library Preparation]({{ site.baseurl }}{% link topics/transcriptomics/images/tenx_libprep.png %} "An overview of the library preparation")

The 10X barcoded gel beads consist of a pool barcodes which are used to separately index each cell's transcriptome. The individual gel barcodes are delivered to each cell via flow-cytometry, where each cell is fed single-file along a liquid tube and tagged with a 10X gel bead. The cells are then isolated from one another within thousands of nanoliter droplets, where each droplet described by a unique 10x barcode that all reads in that droplet are associated with once they undergo reverse-transcription (RT) which reconstructs the mRNA into a cDNA counterpart. The oil is then removed and all (now barcoded) cDNA reads are pooled together to be sequenced.

Though there are ~750,000 10X gel barcodes used, the amount actually qualitatively profiled in a sample is ~10,000 due to majority of droplets (>90%) being empty in order to ensure that the remainder contains only one cell.

More information can be found in the [reagent kit documentation](https://support.10xgenomics.com/single-cell-gene-expression/library-prep/doc/user-guide-chromium-single-cell-3-reagent-kits-user-guide-v2-chemistry).

TODO:
 * Discuss droplet-seq based workflows
 * What is Chromium, and do users need to know what it is?


# Analysis Strategy
{:.no_toc}

The tutorial is structured into two parts:

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

![Overview of workflow]({{ site.baseurl }}{% link topics/transcriptomics/images/tenx_workflow.png %} "An overview of the workflow")

The first part of this tutorial is essentially a one-click "fire and forget" solution to demultiplexing and quantifying scRNA-seq data, where much of the intricacy and complexity required in this extremely crucial stage is obscured behind several layers of safety glass.

However, Science favours the brave, and for those more interested in the gritty details of how FASTQ files are transformed into a count matrix, please see the [Pre-processing of Single-Cell RNA Data]({{ site.baseurl }}{% link topics/transcriptomics/tutorials/scrna-preprocessing/tutorial.md %}) tutorial.


Since STARSolo is a drop-in solution to the CellRanger pipeline, the first part of the tutorial is a one-click solution where users are encouraged to launch their **RNA STARSolo** jobs and then spend the time deciding on what type of caffeine-inspired beverage they can make to maximise their time during the waiting period.

The second part of this tutorial also has a one-click solution if an identical CellRanger matrix is the desired result, but the more interesting aspects of the pipeline are explored in the flexible method part of the tutorial.


# Producing a Count Matrix from FASTQ

Here we will use the [1k PBMCs from a Healthy Donor (v3 chemistry)](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3) from 10x genomics, consisting of 1000 Peripheral blood mononuclear cells (PBMCs) extracted from a healthy donor, where PBMCs are primary cells with relatively small amounts of RNA (~1pg RNA/cell).

The source file consists of 6 FASTQ files split into two sequencing lanes of three reads: R1 (barcodes), R2 (cDNA sequence), I1 (illumina lane info).

The CellRanger pipeline requires all three files to perform the demultiplexing and quantification, but **RNA STAR Solo** does [not require](https://github.com/alexdobin/STAR/issues/640) the I1 lane file to perform the analysis. These source files are provided in the [Zenodo](https://zenodo.org/record/3383115) data repository, but they require approximately 2 hours to process. For this tutorial, we will use datasets sub-sampled from the source files to contain approximately 300 cells instead of 1000. Details of this sub-sampling process can be viewed at the [Zenodo link](https://zenodo.org/api/files/858fcfb3-aca8-408d-afa2-eade47bd623a/subsetting_data.txt).

### Data upload and organization

For the mapping, we require the sub-sampled source files, as well as a list of (737,000) known cell barcodes. The barcodes in the R1 FASTQ data are checked against these known cell barcodes in order assign a specific read to a specific known cell. The barcodes are designed in such a manner that there is virtually no chance that they will align to a place in the reference genome. In this tutorial we will be using hg19 (GRCh37) version of the human genome, and will therefore also need to use a hg19 GTF file to annotate our reads.


> ### {% icon hands_on %} Hands-on: Data upload and organization
>
> 1. Create a new history and rename it (e.g. scRNA-seq 10X dataset tutorial)
>
>    {% include snippets/create_new_history.md %}
>
> 1. Import the sub-sampled FASTQ data and the Cell Barcodes from [`Zenodo`](https://zenodo.org/record/3383115) or from the data library (ask your instructor)
>
>    ```
>    https://zenodo.org/record/3383115/737K-august-2016.txt
>    https://zenodo.org/record/3383115/subset_pbmc_1k_v3_S1_L001_R1_001.fastq.gz
>    https://zenodo.org/record/3383115/subset_pbmc_1k_v3_S1_L001_R2_001.fastq.gz
>    https://zenodo.org/record/3383115/subset_pbmc_1k_v3_S1_L002_R1_001.fastq.gz
>    https://zenodo.org/record/3383115/subset_pbmc_1k_v3_S1_L002_R2_001.fastq.gz
>    ```
>
>    {% include snippets/import_via_link.md %}
>
> 3. Import the Gene Annotations and Cell Barcodes from [`Zenodo`](https://zenodo.org/record/3383115) or from the data library (ask your instructor)
>
>    ```
>    https://zenodo.org/record/3383115/files/Homo_sapiens.GRCh37.75.gtf
>    https://zenodo.org/record/3383115/files/737K-august-2016.txt
>    ```
>
{: .hands_on}


## 10x Chemistries

There are two main reagent kits used during the library preparation, and the choice of one will influence the size of the sequences we work with. Below we can see the layout of the primers used in both chemistries. We can ignore most of these as they are not relevant, namely: the P5 and P7 illumina primers are used in the illumina [bridge amplification](https://en.wikipedia.org/wiki/Illumina_dye_sequencing#Bridge_amplification) process; the Sample Index is an 8bp primer which is related to the Chromium system that balances nucleotide bias and ensures that there is no sample overlap during the multiplexed sequencing; and the Poly(dT) VN primer used to capture RNA sequences with poly-A tails (i.e. mRNA).

![chem]({{ site.baseurl }}{% link topics/transcriptomics/images/tenx_primers.svg %} "10x Chromiumv2 and Chromiumv3 Chemistries")

The primers of interest to us are the Cell Barcode (CB) and the Unique Molecular Identifiers (UMI) used in the Read 1 sequencing primer, as they describe to us how to demultiplex and deduplicate our reads. It is highly advised that the [Plates, Batches, and Barcodes]({{ site.baseurl }}{% link topics/transcriptomics/tutorials/scrna-preprocessing/tutorial.md %}) slides are revisited to refresh your mind on these concepts.


| Chemistry | Read 2 | Read 1 (CB + UMI) | Insert (Read 2 + Read 1) |
|:---------:|:------:|:-----------------:|:------------------------:|
| v2        |  98    | 26 (16 + 10)      | 124 |
| v3        |  91    | 28 (16 + 12)      | 119 |

The table above gives a summary of the primers used in the image and the number of basepairs occupied by each.

> ### {% icon question %} Questions
>
> 1. What has stayed constant between the chemistry versions?
> 1. What advantage does this constant factor give?
> 1. What do UMIs do?
> 1. What advantage does the 2 extra bp in the v3 UMIs have over v2 UMIs?
>
> > ### {% icon solution %} Solution
> >
> > 1. The Cell Barcode (CB) has remained at 16bp for both chemistries.
> > 1. This has the advantage that the same set of barcodes can be used in both chemistries, which is important because barcodes are *very* hard to design.
> >    * They need to be designed in such a way to minimise accidentally aligning to the reference they were prepared to be used for.
> >    * Longer barcodes tend to be more unique, so this is a problem that is being solved as the barcodes increase in size, allowing for barcodes that can be used on more than one reference to be more common, as seen above.
> > 1. UMIs (or Unique Molecular identifiers) do not delineate cells as Cell Barcodes do, but instead serve as random 'salt' that tag molecules randomly and are used to mitigate amplification bias by deduplicating any two reads that map to the same position with the same UMI, where the chance of this happening will be astronomically small unless one read is a direct amplicon of the other.
> > 1. $$4^10 = 1,048,576$$ unique molecules tagged, vs. $$4^12 = 16,777,216$$ unique molecules tagged. The reality is much much smaller due to edit distances being used that would reduce both these numbers substantially (as seen in the [*Plates, Batches, and Barcodes*](({{ site.baseurl }}{% link topics/transcriptomics/tutorials/scrna-plates-batches-barcodes/slides.html %})) slides), but the scale factor of 16 times more molecules ($$4^{12-10} = 16$$) can be uniquely tagged is true
> {: .solution}
{: .question}

The differences in the chemistries is a slight change in the library size, where the v2 aims to capture *on average* 50,000 reads per cell, whereas the v3 aims to capture *at minimum* 20,0000 reads per cell. This greatly reduces the lower-tail of the library size compared to the previous version.


### Determining what Chemistry our Data Contains

To perform the demultiplexing, we need to tell **RNA STAR Solo** where to look in the R1 FASTQ to find the cell barcodes. We can do this by simply counting the number of basepairs in any read of the R1 files.

> ### {% icon question %} Question
> Peek at one of the R1 FASTQ files using the {% icon galaxy-eye %} symbol below the dataset name.
>
> 1. How many basepairs are there in any given read?
> 2. Which library preparation chemistry version was this read generated from?
>
> > ### {% icon solution %} Solution
> > 1. There are 28 basepairs
> > 2. The v2 has 26 basepairs, but the v3 has 28 basepairs. Therefore the reads we have here use the Chromium v3 chemistry.
> {: .solution}
{: .question}


## Performing the Demultiplexing and Quantification

We will now proceed to demultiplex, map, and quantify both sets of reads using the correct chemistry discovered in the previous sub-section.


> ### {% icon hands_on %} Hands-on
>
> **RNA STARSolo** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"RNA-Seq FASTQ/FASTA file, cDNA reads"*: `pbmc_1k_v3_S1_L001_R1_001.fastq.gz`
>    - {% icon param-file %} *"RNA-Seq FASTQ/FASTA file, Barcode reads"*: `pbmc_1k_v3_S1_L001_R2_001.fastq.gz`
>    - {% icon param-repeat %} *Insert Input Pairs*
>    - {% icon param-file %} *"RNA-Seq FASTQ/FASTA file, cDNA reads"*: `pbmc_1k_v3_S1_L002_R1_001.fastq.gz`
>    - {% icon param-file %} *"RNA-Seq FASTQ/FASTA file, Barcode reads"*: `pbmc_1k_v3_S1_L002_R2_001.fastq.gz`
>      (*pay attention to the* **L001** *and* **L002** *names*)
>    - {% icon param-file %} *"RNA-Seq Cell Barcode Whitelist"*: `737K-august-2016.txt`
>    - *"Custom or built-in reference genome"*: `Use a built-in index`
>        - *"Reference genome with or without an annotation"*: `use genome reference without builtin gene-model`
>            - *"Select reference genome"*: `Human (Homo Sapiens): hg19 Full`
>            - *"Gene model (gff3,gtf) file for splice junctions"*: `Homo_sapiens.GRCh37.75.gtf`
>    - In *"Advanced Settings"*:
>        - *"Configure Chemistry Options"*: `Cell Ranger v3`
>
>    > ### {% icon comment %} Comment
>    >
>    > We leave the *Genomic features to collect UMI counts upon* at `Gene` and *UMI deduplication (collapsing) algorithm* at `All`, as these are the options that emulate the CellRanger pipeline.
>    >
>    {: .comment}
>
{: .hands_on}


## Inspecting the Output Files

At this stage **RNA STARSolo** has output 5 files, 2 mapping quality files and 3 matrix files:
 1. Log
 1. Feature Statistic Summaries
 1. Matrix Gene Counts
 1. Barcodes
 1. Genes


### Mapping Quality

Let us investigate the output log. This type of quality control is essential in any RNA-based analysis and it is strongly recommended that you familiarise yourself with the [Quality Control]({{site.baseurl}}{% link topics/sequence-analysis/tutorials/quality-control/tutorial.html}) tutorial.


> ### {% icon hands_on %} Hands-on
>
> **MultiQC** {% icon tool %} with the following parameters:
>    - In *"Results"*:
>        - Under *"1:Results"*:
>            - *"Which tool was used generate logs?"*: `STAR`
>            - In *"STAR output"*:
>                - Under *"1:STAR output"*:
>                    - *"Type of STAR output?"*: `Log`
>                    - {% icon param-file %} *"STAR log output"*: `log output of RNA STARSolo`


> ### {% icon question %} Question
> 
> What percentage of reads are uniquely mapped?
> 
> > ### {% icon solution %} Solution
> > 87.5%
> > - This is good, and is expected of 10x datasets
> >
> {: .solution}
{: .question}

### Quantification Quality

Let us investigate the STARSolo specific log. We can look at this directly by clicking on the {% icon galaxy-eye %} symbol of the *Feature Statistic Summaries* file. 

```
             Barcodes:
           nNinBarcode            358
       nUMIhomopolymer            707
              nTooMany              0
              nNoMatch        7491731
                 Gene:
             nUnmapped            196
            nNoFeature          74087
         nAmbigFeature           9537
 nAmbigFeatureMultimap           7473
              nTooMany            691
         nNoExactMatch           1750
           nExactMatch         106272
                nMatch         106868
         nCellBarcodes            272
                 nUMIs          52748
```

The explanation of these parameters can be seen in the [RNA STAR Manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) under the *STARsolo* section, with most of the information given in the details box below.

> ### {% icon details %} Details
> 
> * nNinBarcode: number of reads with more than 2 Ns in cell barcode (CB)
> * nUMIhomopolymer: number of reads with homopolymer in CB
> * nTooMany: not used at the moment
> * nNoMatch: number of reads with CBs that do not match whitelist even with one mismatch
> 
> All of the above reads are discarded from Solo output.  Remaining reads are checked for overlap with features (e.g. genes):
> 
> * nUnmapped: number of reads unmapped to the genome
> * nNoFeature: number of reads that map to the genome but do not belong to a feature
> * nAmbigFeature: number of reads that belong to more than one feature
> * nAmbigFeatureMultimap: number of reads that belong to more than one feature and are also multimapping to the genome (this is a subset of the `nAmbigFeature`)
> * nTooMany: number of reads with ambiguous CB (i.e. CB matches whitelist with one mismatch but with posterior probability 0.95)
> * nNoExactMatch: number of reads with CB that matches a whitelist barcode with 1 mismatch, but this whitelist barcode does not get any other reads with exact matches of CB
>
> These metrics can be grouped into more broad categories:
> * nNinBarcode+ nUMIhomopolymer + nNoMatch + nTooMany + nNoExactMatch = number of reads with CBs that do not match whitelist.
> * nUnmapped + nAmbigFeature = number of reads without defined feature (gene)
> * nMatch = number of reads that are output as solo counts
>
> The three categoties above summed together should be equal to the total number of reads.
> 
{: .details}

The main information to gather at this stage is that the `nCellBarcodes` tell us how many cells were detected in our sample, where we see 272 which is expected of our sub-sampled data. Another metric to take into account is that the number of matches (`nMatch`) has the largest value, and that the number of reads that map to the genome but not to a feature/gene given in the GTF (`nNoFeature`) is not too large.

# Producing a Quality Count Matrix

The matrix files produced by **RNA STAR Solo** are in the *bundled* format, meaning that the information to create a tabular matrix of Genes vs Cells are separated into different files. These files are already 10x analysis datasets, compatible with any downstream single-cell RNA analysis pipeline, however the number of cells represented here are greatly over-represented, as the they have not yet been filtered for high quality cells, and therefore the matrix represents *any* cells that were unambiguously detected in the sample. 

To get a high quality count matrix we must apply the **DropletUtils** tool, which will produce a filtered set which will produce a dataset more representative of the CellRanger pipeline. 

## CellRanger Method

> ### {% icon hands_on %} Hands-on: Task description
>
> **DropletUtils** {% icon tool %} with the following parameters:
>    - *"Format for the input matrix"*: `Bundled (barcodes.tsv, genes.tsv, matrix.mtx)`
>        - {% icon param-file %} *"Count Data"*: `output_matrix` (output of **RNA STARSolo** {% icon tool %})
>        - {% icon param-file %} *"Genes List"*: `output_genes` (output of **RNA STARSolo** {% icon tool %})
>        - {% icon param-file %} *"Barcodes List"*: `output_barcodes` (output of **RNA STARSolo** {% icon tool %})
>    - *"Operation"*: `Filter for Barcodes`
>        - *"Method"*: `DefaultDrops`
>        - *"Format for output matrices"*: `Tabular`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

## Flexible Method

> ### {% icon hands_on %} Hands-on: Task description
>
> **DropletUtils** {% icon tool %} with the following parameters:
>    - *"Format for the input matrix"*: `Bundled (barcodes.tsv, genes.tsv, matrix.mtx)`
>        - {% icon param-file %} *"Count Data"*: `output_matrix` (output of **RNA STARSolo** {% icon tool %})
>        - {% icon param-file %} *"Genes List"*: `output_genes` (output of **RNA STARSolo** {% icon tool %})
>        - {% icon param-file %} *"Barcodes List"*: `output_barcodes` (output of **RNA STARSolo** {% icon tool %})
>    - *"Operation"*: `Filter for Barcodes`
>        - *"Method"*: `EmptyDrops`
>        - *"Format for output matrices"*: `Tabular`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}



# Conclusion
{:.no_toc}

In this tutorial we have learned the importance of barcoding; namely how to define, extract, and annotate them from our reads and into the read headers, in order to preserve them during mapping. We have discovered how these barcoded reads are transformed into counts, where the cell barcode and UMI barcode are used to denote individual cells and to correct against reads that are PCR duplicates. Finally, we have learned how to combine separate batch data as well as being able to check and correct against cross-contamination.

A full pipeline which produces both an hda5 and tabular file for inspection is provided [here](workflows/scrna_tenx.ga).

<!-- ![Recap of workflow]({{site.baseurl}}{% link topics/transcriptomics/images/scrna_workflow.svg %} "A recap of the entire workflow") -->
