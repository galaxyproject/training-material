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

In this tutorial, we showcase packages from the Bioconductor repository, combined into a workflow from data import to the identification of cell clusters and associated marker genes.

> <comment-title></comment-title>
>
> This tutorial is significantly based on ["Orchestrating Single-Cell Analysis with Bioconductor"](https://bioconductor.org/books/release/OSCA/) {% cite amezquita2019orchestrating %} and draws from the [Clustering 3K PBMCs with Seurat]({% link topics/single-cell/tutorials/scrna-seurat-pbmc3k/tutorial.md %}) and [Clustering 3K PBMCs with Scanpy]({% link topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.md %}) tutorials.
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

# Libraries

This tutorial showcases the following Bioconductor packages:

- [bluster](https://bioconductor.org/packages/bluster/) provides clustering algorithms.
- [DropletUtils](https://bioconductor.org/packages/DropletUtils/) provides utilities for handling single-cell droplet data.
- [LoomExperiment](https://bioconductor.org/packages/LoomExperiment/) provides a means to easily convert the Bioconductor "Experiment" classes to loom files.
- [scater](https://bioconductor.org/packages/scater/) provides a toolkit for the analysis of single-cell gene expression data.
- [scran](https://bioconductor.org/packages/scran/) provides additional methods for single-cell RNA-Seq data analysis.
- [scuttle](https://bioconductor.org/packages/scuttle/) provides additional utility functions for performing single-cell analyses.

We also use the following R packages from the CRAN repository:

- [cowplot](https://cran.r-project.org/package=cowplot) provides various features that help with creating publication-quality figures with [ggplot2](https://cran.r-project.org/package=ggplot2).

# Data

To create a parallel tutorial to the tutorial [Clustering 3K PBMCs with Seurat]({% link topics/single-cell/tutorials/scrna-seurat-pbmc3k/tutorial.md %}), input data files are downloaded from these Zenodo links:

```
https://zenodo.org/record/3581213/files/genes.tsv
https://zenodo.org/record/3581213/files/barcodes.tsv
https://zenodo.org/record/3581213/files/matrix.mtx
```

## Data upload

This section was copied from the tutorial [Clustering 3K PBMCs with Seurat]({% link topics/single-cell/tutorials/scrna-seurat-pbmc3k/tutorial.md %}) to set up input files in the same way.

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

The common data structure adopted by Bioconductor packages to store the matrix as well as gene and cell annotations is called `SingleCellExperiment`.

We can import the matrix and annotations of genes and cells into an `SingleCellExperiment` object using the function `read10xCounts()` of the package [DropletUtils](https://bioconductor.org/packages/DropletUtils/).

We set the `samples=` argument to the directory that contains the three 10x files.
Changing `row.names=` from the default `"id"` to `"symbol"` ensures that human-readable gene symbols are used instead of the database-friendlier Ensembl gene identifiers.
Either way, both Ensembl gene identifiers and gene symbols will be stored in the `rowData()` component of the `SingleCellExperiment` and available at all times.

In the context of a Galaxy workflow, we want to save the R object to a file that can be re-used by downstream tasks.
R objects are traditionally saved to `.RData` or `.rds` files.
However those file formats are specific to R and [deemed insecure](https://github.com/galaxyproject/tools-iuc/issues/3921).

To be saved to the safer `loom` file format, `SingleCellExperiment` objects must be converted to the `SingleCellLoomExperiment` class before they can be saved to disk using the `export()` function from the `BiocIO` package.

Altogether, the first tool in our Galaxy workflow should execute the following code:

```{r}
sce <- DropletUtils::read10xCounts(
	samples = "data/",
	row.names = "symbol"
)
dir.create("outputs")
scle <- as(sce, "SingleCellLoomExperiment")
if (!file.exists("outputs/sce.loom")) {
	export(object = scle, con = "outputs/sce.loom", format = "loom")
}
```

For this task, the tool [DropletUtils Read10x](https://usegalaxy.eu/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fdropletutils_read_10x%2Fdropletutils_read_10x%2F1.0.4%2Bgalaxy0&version=latest) can be used.

However:

- This tool produces an `rdata` file.
  RData objects are deemed insecure as discussed in this [GitHub issue](https://github.com/galaxyproject/tools-iuc/issues/3921).
  The tool should be updated to produce a `loom` file containing a `LoomExperiment` object.
- This tool does not offer the option to use gene symbols as `rownames`.
  The default Ensembl gene identifiers are not immediately interpretable and complicate the interpretation of results later in the workflow (e.g., list of highly variable genes, cluster marker genes).
  The tool should be updated to offer the possibility of using gene symbols as `rownames` for the `SingleCellExperiment`.
  It is worth pointing out that both Ensembl gene identifiers and gene symbols are stored in the `rowData` component of the `SingleCellExperiment` object, meaning they remain available throughout the analysis

### Inspect the `SingleCellExperiment` object.

Given that R sessions do not persist between steps of a Galaxy workflow, a Galaxy tool would need to execute the following code to re-import the `SingleCellLoomExperiment` from the file written in the previous step before displaying a summary view of it.

In the context of a Galaxy workflow, the output would then be saved to a `.txt` file that the user could inspect after the job completes.

```{r}
sce <- import('outputs/sce.loom', format = "loom", type = "SingleCellLoomExperiment")
print(sce)
```

For this task, a tool similar to the [Seurat Data Management](https://usegalaxy.eu/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fseurat_data%2Fseurat_data%2F5.4.0%2Bgalaxy2&version=latest) tool method "Inspect Seurat Object" is needed.
