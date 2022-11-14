---
layout: tutorial_hands_on
subtopic: deconvo
priority: 4
title: Comparing inferred cell compositions using MuSiC deconvolution
zenodo_link: https://zenodo.org/record/7319925
tags:
  - single-cell
  - mouse
  - human
  - deconvolution
  - bulk
  - transcriptomics
questions:
- How do the cell type distributions vary in bulk RNA samples across my variable of interest?
- For example, are beta cell proportions different in the pancreas data from diabetes and healthy patients?
objectives:
- Apply the MuSiC deconvolution to samples and compare the cell type distributions
- Compare the results from analysing different types of input, for example, whether combining disease references or not yields better results
time_estimation: 2H
key_points:
- Deconvolution can be used to compare cell type distributions from bulk RNA-seq datasets
contributors:
- nomadscientist
- mtekman
requirements:
  -
    type: "internal"
    topic_name: single-cell
    tutorials:
      - bulk-music
      - bulk-music-2-preparescref
      - bulk-music-3-preparebulk

gitter: Galaxy-Training-Network/galaxy-single-cell
---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

The goal of this tutorial is to apply bulk RNA deconvolution techniques to a problem with multiple variables - in this case, a model of diabetes is compared with its healthy counterparts. All you need to compare inferred cell compositions are well-annotated, high quality reference scRNA-seq datasets, transformed into MuSiC-friendly Expression Set objects, and your bulk RNA-samples of choice (also transformed into MuSiC-friendly Expression Set objects). For more information on how MuSiC works, you can check out their github site [MuSiC](https://xuranw.github.io/MuSiC/articles/MuSiC.html) or published article {% cite wang2019bulk %}


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Data

In the standard MuSiC tutorial, we used human pancreas data. We will now use the same single cell reference dataset {%cite segerstolpe2016single %} withits 10 samples of 6 healthy subjects and 4 with Type-II diabetes (T2D), as well as the bulk RNA-samples from the same lab. Both of these datasets were accessed from the public EMBL-EBI repositories and transformed into Expression Set objects in the previous two tutorials. For both the single cell reference and the bulk samples of interest, you have generated Expression Set objects with only T2D samples, only healthy samples, and a final everything-combined sample. The plan is to analyse this data in three ways: using a combined reference and combined bulk-files; using only the healthy single cell reference; and separating to infer healthy cells from a healthy reference and diseased cells from a diseased reference.

![Three colours of arrows connect bulk combined and single cell combined; bulk healthy and single cell healthy & bulk diseased with single cell diseased; and bulk diseased and healthy with the single cell healthy reference. These are labelled altogether, like4like, and healthyscref, respectively.](../../images/bulk-music/comparisons.png "Plan of analysis")

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial *"Deconvolution: Compare"*
> 2. Import the files from [Zenodo]({{ page.zenodo_link }})
>
>    * Human single cell RNA ESet objects (tag: `#singlecell`)
>
>      ```
>    {{ page.zenodo_link }}/files/ESet_object_sc_combined.rdata
>    {{ page.zenodo_link }}/files/ESet_object_sc_T2D.rdata
>    {{ page.zenodo_link }}/files/ESet_object_sc_healthy.rdata
>      ```
>
>    * Human bulk RNA ESet objects (tag: `#bulk`)
>      ```
>    {{ page.zenodo_link }}/files/ESet_object_bulk_combined.rdata
>    {{ page.zenodo_link }}/files/ESet_object_bulk_healthy.rdata
>    {{ page.zenodo_link }}/files/ESet_object_bulk_T2D.rdata
>      ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Rename the datasets
>
> 5. Add to each file a tag corresponding to `#bulk` and `#scrna`
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

# Infer cellular composition & compare

It's finally time!

# Comparing combined ESet objects

## Sub-step with **MuSiC Compare**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [MuSiC Compare](toolshed.g2.bx.psu.edu/repos/bgruening/music_compare/music_compare/0.1.1+galaxy4) %} with the following parameters:
>    - In *"New scRNA Group"*:
>        - {% icon param-repeat %} *"Insert New scRNA Group"*
>            - *"Name of scRNA Dataset"*: `scRNA_set`
>            - In *"Advanced scRNA Parameters"*:
>                - *"Cell Types Label from scRNA dataset"*: `Inferred cell type - author labels`
>                - *"Samples Identifier from scRNA dataset"*: `Individual`
>                - *"Comma list of cell types to use from scRNA dataset"*: `alpha cell,beta cell,delta cell,gamma cell,acinar cell,ductal cell`
>            - In *"Bulk Datasets in scRNA Group"*:
>                - {% icon param-repeat %} *"Insert Bulk Datasets in scRNA Group"*
>                    - *"Name of Bulk Dataset"*: `Bulk_set`
>                    - {% icon param-file %} *"Bulk RNA Dataset"*: `output` (Input dataset)
>                    - *"Factor Name"*: `Disease`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}
