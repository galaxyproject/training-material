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

In this tutorial, we showcase packages from the Bioconductor repository, combined into a workflow guiding users from data import to the identification of cell clusters and associated marker genes.

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
Either way, both Ensembl gene identifiers and gene symbols will be stored in the `rowData` component of the `SingleCellExperiment` object and available at all times.

In the context of a Galaxy workflow, we want to save the R object to a file that can be re-used by downstream tasks.
R objects are traditionally saved to `.RData` or `.rds` files.
However those file formats are specific to R and [deemed insecure](https://github.com/galaxyproject/tools-iuc/issues/3921).

To be saved to the safer `loom` file format, `SingleCellExperiment` objects must be converted to `SingleCellLoomExperiment` before they can be saved to disk using the `export()` function from the [BiocIO](https://bioconductor.org/packages/BiocIO/) package.

> <hands-on-title>R code</hands-on-title>
>
> For this task, a Galaxy tool should execute the following code:
>
> ```{r}
> library(DropletUtils)
> library(LoomExperiment)
> sce <- DropletUtils::read10xCounts(
> 	samples = "data/",
> 	row.names = "symbol"
> )
> scle <- as(sce, "SingleCellLoomExperiment")
> export(object = scle, con = "outputs/sce.loom", format = "loom")
> ```
{: .hands_on}

> <comment-title></comment-title>
>
> For this task, the tool [DropletUtils Read10x](https://usegalaxy.eu/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fdropletutils_read_10x%2Fdropletutils_read_10x%2F1.0.4%2Bgalaxy0&version=latest) exists, but:
> 
> - It produces an `rdata` file.
>   RData objects are deemed insecure as discussed in this [GitHub issue](https://github.com/galaxyproject/tools-iuc/issues/3921).
>   The tool should be updated to produce a `loom` file containing a `LoomExperiment` object.
> - It does not offer the option to use gene symbols as `rownames`.
>   The Ensembl gene identifiers used by default are not immediately interpretable to human eyes and complicate the interpretation of results later in the workflow (e.g., list of highly variable genes, cluster marker genes).
>   The tool should be updated to offer the possibility of using gene symbols as `rownames` for the `SingleCellExperiment` object.
>   It is worth pointing out that both Ensembl gene identifiers and gene symbols are automatically stored in the `rowData` component of the `SingleCellExperiment` object, meaning they remain available for conversion and interpretation throughout the analysis.
{: .comment}

### Inspect the SingleCellLoomExperiment object

Having saved the `SingleCellLoomExperiment` object to a `loom` file, a Galaxy tool would need to execute the following code to re-import the `SingleCellLoomExperiment` and display a summary view of it.

In the context of a Galaxy workflow, the output would then be saved to a `.txt` file that the user could inspect after the job completes.

> <hands-on-title>R code</hands-on-title>
>
> For this task, a Galaxy tool should execute the following code:
>
> ```{r}
> library(LoomExperiment)
> sce <- import('outputs/sce.loom', format = "loom", type = "SingleCellLoomExperiment")
> capture.output(print(sce), file = "outputs/sce_summary.txt")
> ```
{: .hands_on}

> <comment-title></comment-title>
>
> For this task, a new tool similar to the [Seurat Data Management](https://usegalaxy.eu/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fseurat_data%2Fseurat_data%2F5.4.0%2Bgalaxy2&version=latest) tool method "Inspect Seurat Object" is needed.
> It seems worth starting a `SingleCellLoomExperiment Data Management` tool suite with that sort of miscellaneous utilities.
{: .comment}

In this instance, the summary view of the object looks as follows:

```
class: SingleCellLoomExperiment 
dim: 32738 2700 
metadata(5): CreatedWith LOOM_SPEC_VERSION LoomExperiment-class MatrixName Samples
assays(1): counts
rownames(32738): MIR1302-10 FAM138A ... AC002321.2 AC002321.1
rowData names(2): ID Symbol
colnames: NULL
colData names(2): Barcode Sample
reducedDimNames(0):
mainExpName: NULL
altExpNames(0):
rowGraphs(0): NULL
colGraphs(0): NULL
```

## Preprocessing

### Computation of QC metrics

#### Identify mitochondrial genes.

In preparation for the calculation of quality control metrics, we need to identify mitochondrial genes in our dataset.

Conveniently, human genes that are encoded in the mitochondrial DNA (rather than in the cell nucleus) have names beginning with `MT-`.

In the context of a Galaxy workflow, that list of rownames should be stored in a `.txt` file that can be passed as input to downstream tools.

> <hands-on-title>R code</hands-on-title>
>
> For this task, a Galaxy tool should execute the following code:
>
> ```{r}
> library(LoomExperiment)
> sce <- import('outputs/sce.loom', format = "loom", type = "SingleCellLoomExperiment")
> mt_gene_ids <- rownames(sce)[grep('^MT-', rownames(sce))]
> cat(mt_gene_ids, file = "outputs/mt_gene_ids.txt", sep = "\n")
> ```
{: .hands_on}

In this instance, the output file - listing rownames corresponding to mitochondrial genes - looks as follows:

```
MT-ND1
MT-ND2
MT-CO1
MT-CO2
MT-ATP8
MT-ATP6
MT-CO3
MT-ND3
MT-ND4L
MT-ND4
MT-ND5
MT-ND6
MT-CYB
```

> <comment-title></comment-title>
>
> For this task, it seems like a new tool could be added to the `SingleCellLoomExperiment Data Management` tool suite proposed above.
{: .comment}

It is worth mentioning that the tutorial [Single-cell quality control with scater]({% link topics/single-cell/tutorials/scrna-scater-qc/tutorial.md %}) provides a precomputed file of mitochondrial genes, which could be a convenient workaround that could be used here to save a step in a workflow and a new tool wrapper in Galaxy.

#### Calculate the Proportion of Mitochondrial Reads

Equipped with a `SingleCellLoomExperiment` object and a file that contains the list of rownames that correspond to the mitochondrial genes, we can use the `addPerCellQCMetrics()` function from the [scuttle](https://bioconductor.org/packages/scuttle/) package to compute some generic quality control metrics alongside metrics focusing on mitochondrial genes.

In particular, the percentage of UMI assigned to mitochondrial genes is frequently used as a measure of cell integrity allowing us to remove droplets that contain damaged or otherwise suspicious cells.

> <hands-on-title>R code</hands-on-title>
>
> For this task, a Galaxy tool should execute the following code:
>
> ```{r}
> library(LoomExperiment)
> library(scuttle)
> sce <- import('outputs/sce.loom', format = "loom", type = "SingleCellLoomExperiment")
> mt_gene_ids <- scan("outputs/mt_gene_ids.txt", what = "character", quiet = TRUE)
> gene_subsets <- list()
> gene_subsets[["MT"]] <- mt_gene_ids
> sce <- scuttle::addPerCellQCMetrics(
> 	x = sce,
> 	subsets = gene_subsets
> )
> export(object = sce, con = "outputs/sce_after_qc.loom", format = "loom")
> ```
{: .hands_on}

> <comment-title></comment-title>
>
> For this task, the tool [Scater: calculate QC metrics](https://usegalaxy.eu/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fscater_create_qcmetric_ready_sce%2Fscater_create_qcmetric_ready_sce%2F1.22.0&version=latest) exists, but:
> 
> - It requires the expression matrix in tabular format.
> - It masks the fact that [scater](https://bioconductor.org/packages/scater/) simply re-exports a functionality that is actually implemented in [scuttle](https://bioconductor.org/packages/scuttle/).
> - A new tool 'scuttle addPerCellQCMetrics' should be created as a replacement that takes a `loom` file containing a `SingleCellLoomExperiment` object and produces a new, updated object.
{: .comment}

### Inspect the quality control metrics

The newly computed per-cell quality control metrics are stored in the `colData` component of the object.
This a tabular `DataFrame` structure which can be inspected by re-importing the `SingleCellLoomExperiment` object and printing the `colData()` component.

In the context of a Galaxy workflow, the output would then be saved to a `.txt` file that the user could inspect after the job completes.

> <hands-on-title>R code</hands-on-title>
>
> For this task, a Galaxy tool should execute the following code:
>
> ```{r}
> library(LoomExperiment)
> sce <- import('outputs/sce_after_qc.loom', format = "loom", type = "SingleCellLoomExperiment")
> write.table(
>	as.data.frame(colData(sce)), file = "outputs/sce_coldata.tsv", sep = "\t", 
>    quote = FALSE, row.names = FALSE)
> ```
{: .hands_on}

> <comment-title></comment-title>
>
> For this task, a new tool similar to the [Seurat Data Management](https://usegalaxy.eu/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fseurat_data%2Fseurat_data%2F5.4.0%2Bgalaxy2&version=latest) tool method "Inspect Seurat Object" is needed (as mentioned earlier), with the option to display only the `colData` component of the object.
{: .comment}

In this instance, the first few lines of the output file would look as follows:

```
Barcode	Sample	detected	subsets_MT_detected	subsets_MT_percent	subsets_MT_sum	sum	total
AAACATACAACCAC-1	data/	781	10	3.0152829409335	73	2421	2421
AAACATTGAGCTAC-1	data/	1352	10	3.79359575769937	186	4903	4903
AAACATTGATCAGC-1	data/	1131	8	0.889171165449349	28	3149	3149
AAACCGTGCTTCCG-1	data/	960	10	1.74308450170519	46	2639	2639
AAACCGTGTATGCG-1	data/	522	5	1.22324159021407	12	981	981
AAACGCACTGGTAC-1	data/	782	7	1.66358595194085	36	2164	2164
AAACGCTGACCAGT-1	data/	783	10	3.81433823529412	83	2176	2176
AAACGCTGGTTCTT-1	data/	790	9	3.09734513274336	70	2260	2260
```

### Filtering of low-quality cells

#### Visualise QC Metrics.

Having computed the quality control metrics and stored them in the `colData` component of the `SingleCellLoomExperiment` object, we can reload that object and use the `plotColData()` function of the [scater](https://bioconductor.org/packages/scater) package to visualise the quality control metrics to determine appropriate thresholds for filtering out low-quality cells.

Which metrics are available can be determined by inspecting the `colData()` component of the object as shown above.

Given that `plotColData()` can only plot one QC metric at a time on the y-axis, we use the `cowplot` package to combine them into a single figure.

> <hands-on-title>R code</hands-on-title>
>
> For this task, a Galaxy tool should execute the following code:
>
> ```{r}
> library(LoomExperiment)
> library(scater)
> sce <- import('outputs/sce_after_qc.loom', format = "loom", type = "SingleCellLoomExperiment")
> plot_list <- list()
> plot_list[["detected"]] <- scater::plotColData(
> 	object = sce,
> 	y = "detected"
> )
> plot_list[["sum"]] <- scater::plotColData(
> 	object = sce,
> 	y = "sum"
> )
> plot_list[["subsets_MT_percent"]] <- scater::plotColData(
> 	object = sce,
> 	y = "subsets_MT_percent"
> )
> cowplot::plot_grid(plotlist = plot_list, nrow = 1)
> ```
{: .hands_on}

> <comment-title></comment-title>
>
> For this task, a new tool similar to the tools [Seurat Visualize](https://usegalaxy.eu/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fseurat_plot%2Fseurat_plot%2F5.4.0%2Bgalaxy2&version=latest) or [scanpy plot](https://usegalaxy.eu/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fscanpy_plot%2Fscanpy_plot%2F1.11.5%2Bgalaxy0&version=latest) is needed, with the `plotColData()` function being one of multiple choices available through the tool.
> The `plotColData()` function is quite versatile and will need careful wrapping to manage required and optional parameters.
> In this case, taking inspiration from the [Seurat Visualize](https://usegalaxy.eu/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fseurat_plot%2Fseurat_plot%2F5.4.0%2Bgalaxy2&version=latest) or [scanpy plot](https://usegalaxy.eu/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fscanpy_plot%2Fscanpy_plot%2F1.11.5%2Bgalaxy0&version=latest) tool, users could provide a comma-separated list of columns names in `colData`, and the tool wrapper would call `scater::plotColData()` for each column before using [cowplot](https://cran.r-project.org/package=cowplot) package to combine the multiple plots into a single image.
{: .comment}

![Violin Plots showing the unique features (detected), total counts (sum) and the proportion of reads coming from mitochondial genes (subsets_MT_percent) for all cells.](../../images/scrna-bioconductor-pbmc3k/scater-plotcoldata-qc-violin.png "Violin Plots showing the unique features (detected), total counts (sum) and the proportion of reads coming from mitochondial genes (subsets_MT_percent) for all cells.")

The `plotColData()` function can also plot one QC metric on the y-axis against another QC metric on the x-axis.
This can be additional insights into the relationship between QC metrics and help determine appropriate thresholds for filtering out low-quality cells.

> <hands-on-title>R code</hands-on-title>
>
> For this task, a Galaxy tool should execute the following code (once for each combination of features):
>
> ```{r}
> library(LoomExperiment)
> library(scater)
> sce <- import('outputs/sce_after_qc.loom', format = "loom", type = "SingleCellLoomExperiment")
> scater::plotColData(
> 	object = sce,
> 	y = "detected",
> 	x = "sum"
> ```
>
> ```{r}
> library(LoomExperiment)
> library(scater)
> sce <- import('outputs/sce_after_qc.loom', format = "loom", type = "SingleCellLoomExperiment")
> scater::plotColData(
> 	object = sce,
> 	y = "subsets_MT_percent",
> 	x = "sum"
> )
> ```
{: .hands_on}

> <comment-title></comment-title>
>
> In contrast to the 'violin plot' use case above, where users can give a list of column names that the wrapper can simply loop over, this 'scatter plot' use case requires two features per plot - one of the y-axis and one for the x-axis - which immediately complicates the idea of producing multiple plots for multiple combinations of features on the X-Y axes.
> Again, taking inspiration from the [Seurat Visualize](https://usegalaxy.eu/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fseurat_plot%2Fseurat_plot%2F5.4.0%2Bgalaxy2&version=latest), the tool should only allow users to specify one feature for the x-axis and one feature for the y-axis, forcing them to call the tool once per combination of features on the X-Y axes.
{: .comment}

![Scatter plots showing the relationships between the total counts (sum) and A. the number of unique features (detected) and B. the proportion of mitochondrial reads (subsets_MT_percent).](../../images/scrna-bioconductor-pbmc3k/scater-plotcoldata-qc-scatter.png "Scatter plots showing the relationships between the total counts (sum) and A. the number of unique features (detected) and B. the proportion of mitochondrial reads (subsets_MT_percent).")

#### Filter Out Low Quality Cells

Having visualised the quality control metrics, we can now filter out low-quality cells based on thresholds for the number of detected genes and the percentage of mitochondrial reads.

> <hands-on-title>R code</hands-on-title>
>
> For this task, a Galaxy tool should execute the following code (once for each combination of features):
>
> ```{r}
> library(LoomExperiment)
> sce <- import('outputs/sce_after_qc.loom', format = "loom", type = "SingleCellLoomExperiment")
> keep_cells <- sce$detected >= 200 & sce$detected <= 2500 & sce$subsets_MT_percent <= 5
> sce <- sce[, keep_cells]
> export(object = sce, con = "outputs/sce_after_qc_filter.loom", format = "loom")
> ```
{: .hands_on}

> <comment-title></comment-title>
>
> Not every user will want to apply filters on all the same metrics.
> Taking inspiration from the [Seurat Create](https://usegalaxy.eu/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fseurat_create%2Fseurat_create%2F5.4.0%2Bgalaxy3&version=latest) tool, the tool should offer two fields for each possible metric - one for the minimum, one for the maximum values allowed - and the logic of the tool should only apply filters where values were given by the user, skipping over the others.
> The resulting script should look something like this:
>
> ```{r}
> keep_cells <- rep(TRUE, ncol(sce))
> if (!is.null(detected_min)) {
>   keep_cells <- keep_cells & sce$detected <= detected_max
> }
> if (!is.null(detected_max)) {
>   keep_cells <- keep_cells & sce$detected >= detected_max
> }
> sce_filtered <- sce[, keep_cells]
> ```
{: .comment}

We can then produce the same plots as earlier to visualise the new distribution of quality control metrics in our data.

> <hands-on-title>R code</hands-on-title>
>
> Same code as earlier, using the filtered `SingleCellLoomExperiment` object.
>
> ```{r}
> library(LoomExperiment)
> library(scater)
> sce <- import('outputs/sce_after_qc.loom', format = "loom", type = "SingleCellLoomExperiment")
> plot_list <- list()
> plot_list[["detected"]] <- scater::plotColData(
> 	object = sce,
> 	y = "detected"
> )
> plot_list[["sum"]] <- scater::plotColData(
> 	object = sce,
> 	y = "sum"
> )
> plot_list[["subsets_MT_percent"]] <- scater::plotColData(
> 	object = sce,
> 	y = "subsets_MT_percent"
> )
> cowplot::plot_grid(plotlist = plot_list, nrow = 1)
> ```
{: .hands_on}

![Violin Plots showing the unique features (detected), total counts (sum) and the proportion of reads coming from mitochondial genes (subsets_MT_percent) for cells that remain after filtering.](../../images/scrna-bioconductor-pbmc3k/scater-plotcoldata-qc-violin-filtered.png "Violin Plots showing the unique features (detected), total counts (sum) and the proportion of reads coming from mitochondial genes (subsets_MT_percent) for cells that remain after filtering")

![Scatter plots showing the relationships between the total counts (sum) and A. the number of unique features (detected) and B. the proportion of mitochondrial reads (subsets_MT_percent) after filtering.](../../images/scrna-bioconductor-pbmc3k/scater-plotcoldata-qc-scatter-filtered.png "Scatter plots showing the relationships between the total counts (sum) and A. the number of unique features (detected) and B. the proportion of mitochondrial reads (subsets_MT_percent) after filtering.")

### Further Preprocessing

#### Log-normalisation

Having filtered droplets that do not pass our chosen quality control filters, the next step in a standard single-cell RNA-sequencing analysis workflow is normalisation.

In Bioconductor, the simple log-normalisation method is implemented in the function `logNormCounts()` of the package [scuttle](https://bioconductor.org/packages/scuttle/).

> <hands-on-title>R code</hands-on-title>
>
> For this task, a Galaxy tool should execute the following code (once for each combination of features):
>
> ```{r}
> library(LoomExperiment)
> library(scuttle)
> sce <- import('outputs/sce_after_qc_filter.loom', format = "loom", type = "SingleCellLoomExperiment")
> sce <- scuttle::logNormCounts(x = sce)
> export(object = sce, con = "outputs/sce_after_lognorm.loom", format = "loom")
> ```
{: .hands_on}

> <comment-title></comment-title>
>
> The [Galaxy ToolShed](https://toolshed.g2.bx.psu.edu/) does not seem to contain any tool offering access to this functionality.
> It seems that a whole new tool is needed.
> It is worth identifying the full set of [scuttle](https://bioconductor.org/packages/scuttle/) functions needed for this tutorial to decide how to design this new tool,
> as we may host all those functions in the same tool suite.
{: .comment}
