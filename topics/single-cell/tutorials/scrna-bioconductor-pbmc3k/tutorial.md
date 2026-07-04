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
---


Single cell RNA-seq analysis enables us to explore differences in gene expression between cells.
It can reveal the heterogenity within cell populations and help us to identify cell types that could play roles in development, disease, or other processes.

In this tutorial, we showcase packages from the Bioconductor repository.

> <comment-title></comment-title>
>
> This tutorial is significantly based on ["Orchestrating Single-Cell Analysis with Bioconductor"](https://bioconductor.org/books/release/OSCA/) {% cite amezquita2019orchestrating %} and draws from the [Clustering 3K PBMCs with Seurat](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-seurat-pbmc3k/tutorial.html) and [Clustering 3K PBMCs with Scanpy](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.html) tutorials.
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

# Data

For this tutorial, we will analyze a dataset of Peripheral Blood Mononuclear Cells (PBMC) extracted from a healthy donor, which is freely available from 10X Genomics. The dataset contains 2700 single cells sequenced using Illumina NextSeq 500. The raw sequences have been processed by the [**cellranger**](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) pipeline from 10X to extract a unique molecular identifier (UMI) count matrix, in a similar way to that explained in the [Pre-processing of 10X Single-Cell RNA Datasets]({% link topics/single-cell/tutorials/scrna-preprocessing-tenx/tutorial.md %}) tutorial.

In this matrix, the values represent the number of each feature (i.e. gene; row) detected in each cell (column). Single cell matrices can be quite large: here there are 2700 columns with 32,738 rows, so for each of our 2700 cells we know how many times we found RNAs matching each of the 32,738 genes. Since most of these genes weren't detected in most of the cells, the matrix is largely filled with zeros, i.e. it is an extremely sparse matrix. To optimize the storage of such a table and the information about the genes and cells, **cellranger** creates 3 files:

- `genes.tsv`: a tabular file with information about the 32,738 genes in 2 columns (Ensembl gene id and the gene symbol)
- `barcodes.tsv`: a tabular file with the barcode for each of the 2700 cells
- `matrix.mtx`: a condensed version of the count matrix (including the non-zero values only)

    The count matrix is represented by its non-zero values - we don't need to store all of those zeroes as long as we know where our non-zero values are in the matrix. Each non-zero value is represented by its line number (1st column), its column number (2nd column) and its value (3rd column). The first row gives the total number of rows (genes), columns (cells) and non-zero values. More information on the Matrix Market Exchange (mtx) format can be found [in this documentation](https://math.nist.gov/MatrixMarket/formats.html)

# Libraries

This tutorial showcases the following Bioconductor packages:

- [bluster](https://bioconductor.org/packages/bluster/) provides clustering algorithms
- [DropletUtils](https://bioconductor.org/packages/DropletUtils/) provides utilities for handling single-cell droplet data
- [scater](https://bioconductor.org/packages/scater/) provides a toolkit for the analysis of single-cell gene expression data
- [scran](https://bioconductor.org/packages/scran/) provides additional methods for single-cell RNA-Seq data analysis
- [scuttle](https://bioconductor.org/packages/scuttle/) provides additional utility functions for performing single-cell analyses.

We also use the following R packages from the CRAN repository:

- [cowplot](https://cran.r-project.org/package=cowplot) provides various features that help with creating publication-quality figures with [ggplot2](https://cran.r-project.org/package=ggplot2)

## Data upload

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial
>
> 2. Import the `genes.tsv`, `barcodes.tsv` and `matrix.mtx` from [Zenodo]({{ page.zenodo_link }}) or from the shared data library
>
>    ```
>    {{ page.zenodo_link }}/files/genes.tsv
>    {{ page.zenodo_link }}/files/barcodes.tsv
>    {{ page.zenodo_link }}/files/matrix.mtx
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets as `genes`, `barcodes`, and `matrix` if necessary
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
> 4. Check the datatypes are correct - the `genes` and `barcodes` files should be tsv or tabular while the `matrix` should be an mtx file
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md %}
>
> 5. Inspect the `matrix` file by clicking on the {% icon galaxy-eye %} icon
>
>    {% snippet faqs/galaxy/histories_dataset_item.md %}
>
{: .hands_on}

The beginning of the file should look like this:

> <question-title></question-title>
>
> ```
> 32738	2700	2286884
> 32709	1	4
> 32707	1	1
> 32706	1	10
> 32704	1	1
> ```
>
> 1. How many non-zero values are in the matrix?
> 2. How many counts were found for the 32,706th gene in the 1st cell?
>
> > <solution-title></solution-title>
> >
> > 1. The first row tells us there are 2,286,884 non-zero values for the 32,738 genes (rows) and 2,700 cells (columns) - so only 2.6% of the 88,392,600 potential values we could have in this matrix are non-zero. Getting rid of all those zeros has made the matrix much more compact.
> > 2. 10 counts were found for the 32,706th row (gene) and 1st column (cell), so we collected 10 RNAs that the first cell had produced from this particular gene.
> >
> {: .solution}
>
{: .question}

Representing the matrix with these three files is convenient for sharing the data, but not for processing them. Different single cell analysis packages have attempted to solve the problem of storage and analysis by inventing their own formats, which has led to the proliferation of many different 'standards' in the scRNA-seq package ecosystem.

## SingleCellExperiment objects



```{r}
sce <- DropletUtils::read10xCounts(
    samples = "data",
    row.names = "symbol"
)
sce
```

For this task, the tool [DropletUtils Read10x](https://usegalaxy.eu/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fdropletutils_read_10x%2Fdropletutils_read_10x%2F1.0.4%2Bgalaxy0&version=latest) can be used.

Notes:

- This tool produces an `rdata` file.
  RData objects are deemed insecure as discussed in this [GitHub issue](https://github.com/galaxyproject/tools-iuc/issues/3921).
  The tool should be updated to produce a `loom` file containing a `LoomExperiment` object.
- This tool does not offer the option to use gene symbols as rownames.
  The default Ensembl gene identifiers complicates the interpretation of results later in the workflow.
  The tool should be updated to offer the possibility of using gene symbols.
  