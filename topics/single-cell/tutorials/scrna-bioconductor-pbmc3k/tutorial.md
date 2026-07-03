---
layout: tutorial_hands_on
title: Clustering 3K PBMCs with Bioconductor
level: Introductory
subtopic: firstsc
priority: 3
zenodo_link: https://zenodo.org/record/3581213
questions:
- How can we identify cell types in single cell RNA-Seq data?
- What are the steps for clustering single cell data with Bioconductor packages?
objectives:
- Explain the steps involved in clustering single cell data
- Evaluate the quality of single cell data and filter out low quality cells
- Prepare single cell data for analysis with Bioconductor packages
- Perform clustering with Bioconductor packages
- Be ready to apply apply Bioconductor tools to new datasets
time_estimation: 8H
key_points:
- Bioconductor is a large repository of software packages, some of them for single cell data analysis
- Clustering makes single cell datasets easier for us to understand
- Different tools and parameters should be considered when analysing different datasets
requirements:
- type: internal
  topic_name: single-cell
  tutorials:
  - scrna-preprocessing
  - scrna-preprocessing-tenx
follow_up_training:
- type: internal
  topic_name: single-cell
  tutorials:
  - EBI-retrieval
tags:
- 10x
contributions:
  authorship:
  - kevinrue
recordings:
- youtube_id: -tIOSBXeXPs
  length: 2H11M
  galaxy_version: 24.2.4.dev0
  date: '2025-04-15'
  speakers:
  - MarisaJL
  captioners:
  - MarisaJL
  bot-timestamp: 1744754154


---


Single cell RNA-seq analysis enables us to explore differences in gene expression between cells. It can reveal the heterogenity within cell populations and help us to identify cell types that could play roles in development, disease, or other processes. Single cell omics is a relatively young field, but there are a few commonly-used analysis pipelines that you will often see in the literature. In this tutorial, we will use one of these pipelines, Seurat, to cluster single cell data from a 10X Genomics experiment ({% cite seurat2023v5 %}). You can follow the same analysis using the Scanpy pipeline in the [Clustering 3K PBMCs with Scanpy]({% link topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.md %}) tutorial.

Clustering is typically the first type of analysis we will perform on a single cell dataset. It groups together cells that are expressing similar genes, which makes the data easier to understand and often helps us to identify specific cell types.

> <comment-title></comment-title>
>
> This tutorial is based on the [Seurat - Guided Clustering Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial). The SCTransform sections also draw from the [Using sctransform in Seurat](https://satijalab.org/seurat/articles/sctransform_vignette.html) tutorial.
>
{: .comment}

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Important tips for easier analysis

{% snippet topics/single-cell/faqs/single_cell_omics.md %}

{% snippet faqs/galaxy/tutorial_mode.md %}

{% snippet faqs/galaxy/analysis_troubleshooting.md sc=true %}

# Data

For this tutorial, we will analyze a dataset of Peripheral Blood Mononuclear Cells (PBMC) extracted from a healthy donor, which is freely available from 10X Genomics. The dataset contains 2700 single cells sequenced using Illumina NextSeq 500. The raw sequences have been processed by the [**cellranger**](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) pipeline from 10X to extract a unique molecular identifier (UMI) count matrix, in a similar way to that explained in the [Pre-processing of 10X Single-Cell RNA Datasets]({% link topics/single-cell/tutorials/scrna-preprocessing-tenx/tutorial.md %}) tutorial.

In this matrix, the values represent the number of each feature (i.e. gene; row) detected in each cell (column). Single cell matrices can be quite large: here there are 2700 columns with 32,738 rows, so for each of our 2700 cells we know how many times we found RNAs matching each of the 32,738 genes. Since most of these genes weren't detected in most of the cells, the matrix is largely filled with zeros, i.e. it is an extremely sparse matrix. To optimize the storage of such a table and the information about the genes and cells, **cellranger** creates 3 files:

- `genes.tsv`: a tabular file with information about the 32,738 genes in 2 columns (Ensembl gene id and the gene symbol)
- `barcodes.tsv`: a tabular file with the barcode for each of the 2700 cells
- `matrix.mtx`: a condensed version of the count matrix (including the non-zero values only)

    The count matrix is represented by its non-zero values - we don't need to store all of those zeroes as long as we know where our non-zero values are in the matrix. Each non-zero value is represented by its line number (1st column), its column number (2nd column) and its value (3rd column). The first row gives the total number of rows (genes), columns (cells) and non-zero values. More information on the Matrix Market Exchange (mtx) format can be found [in this documentation](https://math.nist.gov/MatrixMarket/formats.html)
