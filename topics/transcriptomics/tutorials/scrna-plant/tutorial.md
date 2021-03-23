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


## CSV to AnnData

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Import Anndata and loom](toolshed.g2.bx.psu.edu/repos/iuc/anndata_import/anndata_import/0.7.5+galaxy0) %} with the following parameters:
>    - *"hd5 format to be created"*: `Anndata file`
>        - *"Format for the annotated data matrix"*: `Tabular, CSV, TSV`
>            - {% icon param-file %} *"Annotated data matrix"*: Multi-select both `SHR` and `WT` datasets
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


## Merge Batches and Relabel

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.7.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `SHR dataset` (output of **Import Anndata and loom** {% icon tool %})
>    - *"Function to manipulate the object"*: `Concatenate along the observations axis`
>        - {% icon param-file %} *"Annotated data matrix to add"*: `WT dataset` (output of **Import Anndata and loom** {% icon tool %})
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
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.7.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata` (output of **Manipulate AnnData** {% icon tool %})
>    - *"Function to manipulate the object"*: `Rename categories of annotation`
>        - *"Key for observations or variables annotation"*: `batch`
>        - *"Comma-separated list of new categories"*: `shr, wt`
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

## Quality Control

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata` (output of **Manipulate AnnData** {% icon tool %})
>    - *"Method used for inspecting"*: `Calculate quality control metrics, using 'pp.calculate_qc_metrics'`
>
> 1. **TODO** {% tool [Scatter Plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata` (output of **Manipulate AnnData** {% icon tool %})
>    - *"Method used for inspecting"*: `Calculate quality control metrics, using 'pp.calculate_qc_metrics'`
> 
{: .hands_on}


### Filter

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Inspect and manipulate** {% icon tool %})
>    - *"Method used for filtering"*: `Filter cell outliers based on counts and numbers of genes expressed, using 'pp.filter_cells'`
>        - *"Filter"*: `Minimum number of genes expressed`
>            - *"Minimum number of genes expressed required for a cell to pass filtering"*: `TBD`
>
> 1. {% tool [Filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Filter** {% icon tool %})
>    - *"Method used for filtering"*: `Filter genes based on number of cells or counts, using 'pp.filter_genes'`
>        - *"Filter"*: `Minimum number of cells expressed`
>            - *"Minimum number of cells expressed required for a gene to pass filtering"*: `TBD`
>
> 1. {% tool [Filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Filter** {% icon tool %})
>    - *"Method used for filtering"*: `Filter genes based on number of cells or counts, using 'pp.filter_genes'`
>        - *"Filter"*: `Maximum number of counts`
>            - *"Maximum number of counts required for a gene to pass filtering"*: `TBD`
>
> 1. {% tool [Filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Filter** {% icon tool %})
>    - *"Method used for filtering"*: `Filter cell outliers based on counts and numbers of genes expressed, using 'pp.filter_cells'`
>        - *"Filter"*: `Maximum number of counts`
>            - *"Maximum number of counts required for a cell to pass filtering"*: `TBD`
>
{: .hands_on}

## Freeze Raw Data

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.7.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Filter** {% icon tool %})
>    - *"Function to manipulate the object"*: `Freeze the current state into the 'raw' attribute`
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

## Confounder Removal

Normalizing, scaling and regressing out library size

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Normalize](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_normalize/scanpy_normalize/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata` (output of **Manipulate AnnData** {% icon tool %})
>    - *"Method used for normalization"*: `Normalize counts per cell, using 'pp.normalize_total'`
>        - *"Target sum"*: `10000.0`
>        - *"Exclude (very) highly expressed genes for the computation of the normalization factor (size factor) for each cell"*: `No`
>
> 1. {% tool [Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Normalize** {% icon tool %})
>    - *"Method used for inspecting"*: `Logarithmize the data matrix, using 'pp.log1p'`
>
> 1. {% tool [Remove confounders](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_remove_confounders/scanpy_remove_confounders/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Inspect and manipulate** {% icon tool %})
>    - *"Method used for plotting"*: `Regress out unwanted sources of variation, using 'pp.regress_out'`
>        - *"Keys for observation annotation on which to regress on"*: `total_counts`
>
> 1. {% tool [Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Remove confounders** {% icon tool %})
>    - *"Method used for inspecting"*: `Scale data to unit variance and zero mean, using 'pp.scale'`
>        - *"Maximum value"*: `10.0`
>
{: .hands_on}

## Dimension Reduction

PCA and UMAP

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Cluster, infer trajectories and embed](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_cluster_reduce_dimension/scanpy_cluster_reduce_dimension/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Inspect and manipulate** {% icon tool %})
>    - *"Method used"*: `Computes PCA (principal component analysis) coordinates, loadings and variance decomposition, using 'tl.pca'`
>        - *"Type of PCA?"*: `Full PCA`
>            - *"SVD solver to use"*: `ARPACK wrapper in SciPy`
>
> 1. {% tool [Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Cluster, infer trajectories and embed** {% icon tool %})
>    - *"Method used for inspecting"*: `Compute a neighborhood graph of observations, using 'pp.neighbors'`
>        - *"The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation"*: `10`
>        - *"Number of PCs to use"*: `40`
>
> 1. {% tool [Cluster, infer trajectories and embed](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_cluster_reduce_dimension/scanpy_cluster_reduce_dimension/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Inspect and manipulate** {% icon tool %})
>    - *"Method used"*: `Embed the neighborhood graph using UMAP, using 'tl.umap'`
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

## Clustering

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Cluster, infer trajectories and embed](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_cluster_reduce_dimension/scanpy_cluster_reduce_dimension/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Cluster, infer trajectories and embed** {% icon tool %})
>    - *"Method used"*: `Cluster cells into subgroups, using 'tl.leiden'`
>        - *"Coarseness of the clusterin"*: `0.3`
>
> 1. {% tool [Plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Cluster, infer trajectories and embed** {% icon tool %})
>    - *"Method used for plotting"*: `Embeddings: Scatter plot in UMAP basis, using 'pl.umap'`
>        - *"Keys for annotations of observations/cells or variables/genes"*: `leiden, batch`
>        - *"Show edges?"*: `No`
>        - In *"Plot attributes"*:
>            - *"Location of legend"*: `on data`
>            - *"Legend font size"*: `14`
>            - *"Colors to use for plotting categorical annotation groups"*: `rainbow (Miscellaneous)`
>
{: .hands_on}

## DotPlot and Validating Cell Types

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Cluster, infer trajectories and embed** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Makes a dot plot of the expression values, using 'pl.dotplot'`
>        - *"Variables to plot (columns of the heatmaps)"*: `Subset of variables in 'adata.var_names'`
>            - *"List of variables to plot"*: `AT1G26680,COBL2,AT3G55180,NCED2,AT5G22550,ATL63,RGF8,TEL1,DUF9,AGL42,PLT1,PLT2,PER3,PER9,PBL15,TBL30,AT2G48130,AT3G22620,AED3,T3P18-7,MBK21-8,NPF6-4,CYP71B34,IAMT1,AT1G31950,F14G6-22,GL2,ANL2,PMEI4,MUD21-4,EXPA7,MES15,ADF8,COBL9,TPR14,BHLH144,MYB83,SCPL35,MYB46,COBL4,PIN7,PIN3   `
>        - *"The key of the observation grouping to consider"*: `leiden`
>        - *"Number of categories"*: `9`
>        - *"Use 'raw' attribute of input if present"*: `Yes`
>        - In *"Group of variables to highlight"*:
>            - {% icon param-repeat %} *"Insert Group of variables to highlight"*
>                - *"Start"*: `0`
>                - *"End"*: `5`
>                - *"Label"*: `Columella`
>            - {% icon param-repeat %} *"Insert Group of variables to highlight"*
>                - *"Start"*: `6`
>                - *"End"*: `9`
>                - *"Label"*: `QC`
>            - {% icon param-repeat %} *"Insert Group of variables to highlight"*
>                - *"Start"*: `10`
>                - *"End"*: `11`
>                - *"Label"*: `NC`
>            - {% icon param-repeat %} *"Insert Group of variables to highlight"*
>                - *"Start"*: `12`
>                - *"End"*: `17`
>                - *"Label"*: `Endodermis`
>            - {% icon param-repeat %} *"Insert Group of variables to highlight"*
>                - *"Start"*: `18`
>                - *"End"*: `23`
>                - *"Label"*: `Cortex`
>            - {% icon param-repeat %} *"Insert Group of variables to highlight"*
>                - *"Start"*: `24`
>                - *"End"*: `29`
>                - *"Label"*: `Atrichoblast`
>            - {% icon param-repeat %} *"Insert Group of variables to highlight"*
>                - *"Start"*: `30`
>                - *"End"*: `34`
>                - *"Label"*: `Trichoblast`
>            - {% icon param-repeat %} *"Insert Group of variables to highlight"*
>                - *"Start"*: `35`
>                - *"End"*: `39`
>                - *"Label"*: `Xylem`
>            - {% icon param-repeat %} *"Insert Group of variables to highlight"*
>                - *"Start"*: `40`
>                - *"End"*: `41`
>                - *"Label"*: `VC`
>        - *"Custom figure size"*: `No: the figure width is set based on the number of variable names and the height is set to 10.`
>
{: .hands_on}

## Relabel clusters

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.7.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Cluster, infer trajectories and embed** {% icon tool %})
>    - *"Function to manipulate the object"*: `Rename categories of annotation`
>        - *"Key for observations or variables annotation"*: `leiden`
>        - *"Comma-separated list of new categories"*: `0, 1, 2, 3, trichoblast, 5, 6, cortex, 8, 9, 10, columella+QC+NC, xylem, 13`
>
> 1. {% tool [Plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata` (output of **Manipulate AnnData** {% icon tool %})
>    - *"Method used for plotting"*: `Embeddings: Scatter plot in UMAP basis, using 'pl.umap'`
>        - *"Keys for annotations of observations/cells or variables/genes"*: `leiden, batch`
>        - *"Show edges?"*: `No`
>        - In *"Plot attributes"*:
>            - *"Location of legend"*: `on data`
>            - *"Legend font size"*: `14`
>            - *"Colors to use for plotting categorical annotation groups"*: `rainbow (Miscellaneous)`
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}


## Conclusions

None yet. We need to know whether or not we do lineage, maybe using Seurat?