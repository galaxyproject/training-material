---
layout: tutorial_hands_on
title: "Downstream Single-cell RNA Plant analysis with ScanPy"
subtopic: single-cell
priority:

zenodo_link: 'https://zenodo.org/record/4597857'
tags:
  - single-cell
questions:
  - Can we reclaim cell markers using a different analysis method?
  - Are highly variable genes paramount to the analysis?
objectives:
  - Fill in later
requirements:
  -
    type: "internal"
    topic_name: transcriptomics
    tutorials:
        - scrna-scanpy-pbmc3k

time_estimation: 1H
key_points:
  - Fill in later
contributors:
  - mtekman

gitter: Galaxy-Training-Network/galaxy-single-cell


---


# Introduction
{:.no_toc}

> ### {% icon comment %} Comment
>
> Please familiarise yourself with the ["Clustering 3K PBMCs with ScanPy"]() tutorial first, as much of the concepts and processes are the same.
>
{: .comment}

This tutorial replicates the paper ["Spatiotemporal Developmental Trajectories in the Arabidopsis Root Revealed Using High-Throughput Single-Cell RNA Sequencing"](https://doi.org/10.1016/j.devcel.2019.02.022).


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data

In this paper we replicate the analysis of the plant dataset.

Some intro on why plant datasets are the new thing, why it's worth studying.

We have two datasets, wildtype and mutant


* GSE123818_Root_single_cell_shr_datamatrix.fixednames.transposed.csv.gz
* GSE123818_Root_single_cell_wt_datamatrix.fixednames.transposed.csv.gz

As explained in the Zenodo link, the datasets have been modified to use more common gene names and the datasets have been transposed to better suit AnnData.


## Data upload

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the `GSE123818_Root_single_cell_shr_datamatrix.fixednames.transposed.csv.gz` and `GSE123818_Root_single_cell_wt_datamatrix.fixednames.transposed.csv.gz` from [Zenodo]({{ page.zenodo_link }}) or from the shared data library
>
>    ```
>    {{ page.zenodo_link }}/files/GSE123818_Root_single_cell_shr_datamatrix.fixednames.transposed.csv.gz
>    {{ page.zenodo_link }}/files/GSE123818_Root_single_cell_wt_datamatrix.fixednames.transposed.csv.gz
>    ```
>
>    {% snippet snippets/import_via_link.md %}
>    {% snippet snippets/import_from_data_library.md %}
>
> 3. Rename the datasets to `cells_shr` and `cells_wt`.
{: .hands_on}

> ### {% icon question %} Questions
>
> We can peek at the `cells_shr` file, and determine the dimensionality and naming scheme of the data. The rows and the columns depict different variables.
>
> 1. Which are the cell barcodes? Rows or Columns?
> 1. Which are the gene names? Rows or Columns?
> 1. How many cells in the dataset?
> 1. How many genes in the dataset?
>
> > ### {% icon solution %} Solution
> >
> > 1. Rows, notice the long A-C-T-G based names.
> > 1. Columns, "ATXG" are the prefixes commonly used by XXXX.
> > 1. There are 1 099 lines and therefore 1 099 cells
> > 1. If we scroll the preview window all the way to the right we see 27630 as the last number. Given that the first column are the cell names, we must subtract by 1 to get 27 629 cells.
> >
> {: .solution}
>
{: .question}

If the above feels like a convoluted way to get the dimensionality, that's because we haven't imported the data into the right format. For this we need to import both datasets into a single AnnData object (see the [AnnData]() section in the "ScanPy Tutorial" tutorial).


## AnnData

> ### {% icon hands_on %} Hands-on: Transform matrix and all into AnnData object
>
> 1. {% tool [Import Anndata and loom](toolshed.g2.bx.psu.edu/repos/iuc/anndata_import/anndata_import/0.7.5+galaxy0) %} with the following parameters:
>    - *"hd5 format to be created"*: `Anndata file`
>    - *"Format for the annotated data matrix"*: `Tabular, CSV, TSV`
>        - {% icon param-file %} *"Annotated data matrix"*: (Select Multiple Datasets) `cells_shr, cells_wt`
>        - *"Does the first column store the row names?"*: `Yes`
>
{: .hands_on}

Inspect and verify that the AnnData datasets reflect the dimensionality that we see before.

(Peeking exercise)


We now concatenate the datasets into the same file using the manipulate tool


> 9. {% tool [Manipulate Anndata](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.7.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Output of AnnData import (WT)`
>    - *"Function to manipulate the object"*: `Concatenate along the observations axis`
>      - {% icon param-file %} *"Annotate data matrix to add"*: `Output of AnnData import (SHR)`
>      - *"Join method"*: `Intersection of variables`
>      - *"Key to add the batch annotation to obs"*: `batch`
>      - *"Separator to join the existing index names with the batch category"*: `-`

We do the rest when the tools are updated