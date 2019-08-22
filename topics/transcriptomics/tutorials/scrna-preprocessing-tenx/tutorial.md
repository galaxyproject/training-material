---
layout: tutorial_hands_on
title: "Pre-processing of 10X Single-Cell RNA Datasets"
zenodo_link: "https://zenodo.org/record/UPDATETHIS"
tags:
  - single-cell
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

### Data upload and organization

> ### {% icon hands_on %} Hands-on: Data upload and organization
>
> 1. Create a new history and rename it (e.g. scRNA-seq 10X dataset tutorial)
>
>    {% include snippets/create_new_history.md %}
>
> 1. Import the subset FASTQ paired data from [`Zenodo`](https://zenodo.org/record/UPDATETHIS) or from the data library (ask your instructor)
>
>    ```
>    https://zenodo.org/record/UPDATETHIS/files/
>    https://zenodo.org/record/UPDATETHIS/files/
>    ```
>
>    {% include snippets/import_via_link.md %}
>
> 3. Import the Gene Annotations and Barcodes from [`Zenodo`](https://zenodo.org/record/UPDATETHIS) or from the data library (ask your instructor)
>
>    ```
>    https://zenodo.org/record/UPDATETHIS/files/
>    https://zenodo.org/record/UPDATETHIS/files/
>    ```
>
> 4. Set the datatype of the `something` to `tabular`
>
{: .hands_on}


> ### {% icon comment %} BCL files and Lane Information
> STARSolo doesn't use this, link to the [original issue](https://github.com/alexdobin/STAR/issues/640)
{: .comment}



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
> 1. What has stayed constant between the chemistries?
> 1. What advantage does this have?
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

The differences in the chemistries is slight, when the v2 aims to capture *on average* 50,000 reads per cell, whereas the v3 aims to capture *at minimum* 20,0000 reads per cell.


### Determining what Chemistry our Data Contains

We inspect the lengths

> ### {% icon comment %} Note
> 
{: .comment}


## Performing the Demultiplexing and Quantification


> ### {% icon hands_on %} Hands-on
>
> 1. **RNA STARSolo** {% icon tool %} with the following parameters:
>
{: .hands_on}


## Inspecting the Unfiltered Output Matrix


# Producing a Quality Count Matrix

## CellRanger Method


## Flexible Method

# Conclusion
{:.no_toc}

In this tutorial we have learned the importance of barcoding; namely how to define, extract, and annotate them from our reads and into the read headers, in order to preserve them during mapping. We have discovered how these barcoded reads are transformed into counts, where the cell barcode and UMI barcode are used to denote individual cells and to correct against reads that are PCR duplicates. Finally, we have learned how to combine separate batch data as well as being able to check and correct against cross-contamination.

<!-- ![Recap of workflow]({{site.baseurl}}{% link topics/transcriptomics/images/scrna_workflow.svg %} "A recap of the entire workflow") -->
