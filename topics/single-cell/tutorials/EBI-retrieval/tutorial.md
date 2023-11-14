---
layout: tutorial_hands_on
subtopic: datamanipulation
priority: 3
title: EBI Single Cell Expression Atlas files to AnnData | Creating preprocessed dataset for sc-RNA Filter, Plot, Explore tutorial
questions:
- How do I use the EBI Single Cell Expression Atlas?
- How can I reformat and manipulate the downloaded files to create the correct input for downstream analysis?
objectives:
- You will retrieve raw data from the EBI Single Cell Expression Atlas.
- You will manipulate the metadata and matrix files.
- You will combine the metadata and matrix files into an AnnData object for downstream analysis.
  
time_estimation: "15m"
key_points:
- The EMBL-EBI Single-cell Expression Atlas contains high quality datasets.
- Metadata manipulation is key for generating the correctly formatted files.
- To use Scanpy tools, you have to transform your metadata into an AnnData object.
contributions:
  authorship:
    - wee-snufkin
  testing:
    - nomadscientist
  funding: 
  - elixir-fair-data

requirements:
  -
    type: "internal"
    topic_name: single-cell
    tutorials:
      - scrna-case_alevin
      - scrna-case_alevin-combine-datasets

follow_up_training:
  -
    type: "internal"
    topic_name: single-cell
    tutorials:
        - scrna-case_basic-pipeline

tags:
  - transcriptomics
  - data management
---

# Introduction 

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Getting data from Single Cell Expression Atlas

If you happen to be interested in analysing publicly available data, particularly from the [Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/home), you may be interested in the following tool {% cite Moreno2020.04.08.032698 %} which combines all the preprocessing steps shown in [the previous tutorial]({% link topics/single-cell/tutorials/scrna-case_alevin/tutorial.md %}) into one! For this tutorial, the dataset can be seen [at the EBI](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-6945/downloads) with experiment id of `E-MTAB-6945`.

> <hands-on-title>Retrieving data from Single Cell Expression Atlas</hands-on-title>
>
> 1. {% tool [EBI SCXA Data Retrieval](toolshed.g2.bx.psu.edu/repos/ebi-gxa/retrieve_scxa/retrieve_scxa/v0.0.2+galaxy2) %} with the following parameters:
>      - *"SC-Atlas experiment accession"*: `E-MTAB-6945`
>      - *"Choose the type of matrix to download"*: `Raw filtered counts`
>
{: .hands_on}

It's important to note that this matrix is processed somewhat through the SCXA pipeline, which is quite similar to the pre-processing that has been shown in this case study tutorial series, and it contains any and all metadata provided by their pipeline as well as the authors (for instance, more cell or gene annotations). So don't worry if the plots generated using this input method are slightly different! 

# Metadata manipulation

Before creating an AnnData object, we need to make a small modification in experimental design table. The dataset contains information about 7 samples N701 – N707), however in the experimental design table (cell metadata) they are just numbered from 1 to 7. The plotting tool that we will going to use later will fail if the entries are integers and not categoricals, so we will change "1" to "N01" and so on. You can simply preview the experimental design dataset and move to the column "Sample Characteristic[individual]" (that's where the information about batch is - don't worry, we will rename the column header later!). Make a note of the number of that column - number 12 - we will need it to change the batch number to batch name. 

> <hands-on-title> Change batch numbers into names </hands-on-title>
>
> 1. Change the datatype of `EBI SCXA Data Retrieval on E-MTAB-6945 exp_design.tsv` to `tabular`:
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular" %}
>
> 2. {% tool [Column Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regexColumn1/1.0.3) %} with the following parameters:
>    - *"Select cells from"*: `EBI SCXA Data Retrieval on E-MTAB-6945 exp_design.tsv`
>    - *"using column"*: `Column: 12`
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `1`
>            - *"Replacement"*: `N01`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `2`
>            - *"Replacement"*: `N02`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `3`
>            - *"Replacement"*: `N03`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `4`
>            - *"Replacement"*: `N04`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `5`
>            - *"Replacement"*: `N05`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `6`
>            - *"Replacement"*: `N06`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `7`
>            - *"Replacement"*: `N07`
>
> 4. Rename {% icon galaxy-pencil %} output `Cell metadata`
> 
{: .hands_on}

Now we can create an AnnData object!

# Creating AnnData object

> <hands-on-title> Task description </hands-on-title>
>
> 
> 1. {% tool [Scanpy Read10x](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_read_10x/scanpy_read_10x/1.8.1+galaxy9) %}
> 2. Make sure you are using version **1.8.1+galaxy9** of the tool (change by clicking on {% icon tool-versions %} Versions button):
>   ![List of available tool versions shown when clicking on the 'Versions' button on the top of the page.](../../images/scrna-casestudy/version.png "How to change the version of the tool")
> 3. Set the following parameters:
>    - *"Expression matrix in sparse matrix format (.mtx)"*: `EBI SCXA Data Retrieval on E-MTAB-6945 matrix.mtx (Raw filtered counts)`
>    - *"Gene table"*:  `EBI SCXA Data Retrieval on E-MTAB-6945 genes.tsv (Raw filtered counts)`
>    - *"Barcode/cell table"*: `EBI SCXA Data Retrieval on E-MTAB-6945 barcodes.tsv (Raw filtered counts)`
>    - *"Cell metadata table"*: `Cell metadata`
>
> 4. Rename {% icon galaxy-pencil %} output `AnnData object`
> 
{: .hands_on}

# AnnData manipulation

Now we will do several modifications within the AnnData object so that you can follow this tutorial despite the other way of getting data! 
We would like to flag mitochondrial genes. They can be identified quite easily since they names start with mt. Since the tool for flagging the mitochondrial genes is case-sensitive, it might be a good idea to check what is the formatting of mitochondrial genes in our dataset.

> <hands-on-title> Check the format of mitochondrial genes names </hands-on-title>
>
> 1. {% tool [Search in textfiles](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_grep_tool/1.1.1) %} with the following parameters:
>    - *"Select lines from"*: `EBI SCXA Data Retrieval on E-MTAB-6945 genes.tsv (Raw filtered counts)`
>    - *"that"*: `Match`
>    - *"Regular Expression"*: `mt`
>    - *"Match type"*: `case insensitive`
>    - *"Output"*: `Highlighted HTML (for easier viewing)`
> 
> 3. Rename {% icon galaxy-pencil %} output `Mito genes check`
>
{: .hands_on}

If you click on that dataset, you will see all the genes containing 'mt' in their name. We can now clearly see that mitochondrial genes in our dataset start with 'mt-'. Keep that in mind, we will use it in a moment!

Speaking about gene names, we will also change the header of the column containing those names from `gene_symbols` to `Symbol`. This edit is only needed to make our AnnData object compatible with this tutorial's workflow. 

As I mentioned at the beginning, we will also change the header of the column storing information about batch. Actually, we will change several other headers as well.

And the good news is that we can do all those steps using only one tool!

> <hands-on-title> Modify AnnData object </hands-on-title>
>
> 1. {% tool [AnnData Operations](toolshed.g2.bx.psu.edu/repos/ebi-gxa/anndata_ops/anndata_ops/1.8.1+galaxy92) %}
> 2. Make sure you are using version **1.8.1+galaxy92** of the tool (change by clicking on {% icon tool-versions %} Versions button)
> 3. Set the following parameters:
>    - In *"Input object in hdf5 AnnData format"*: `AnnData object`
>    - In *"Change field names in AnnData observations"*:
>        - {% icon param-repeat %} *"Insert Change field names in AnnData observations"*
>            - *"Original name"*: `Sample Characteristic[genotype]`
>            - *"New name"*: `genotype`
>        - {% icon param-repeat %} *"Insert Change field names in AnnData observations"*
>            - *"Original name"*: `Sample Characteristic[individual]`
>            - *"New name"*: `batch`
>        - {% icon param-repeat %} *"Insert Change field names in AnnData observations"*
>            - *"Original name"*: `Sample Characteristic[sex]`
>            - *"New name"*: `sex`
>        - {% icon param-repeat %} *"Insert Change field names in AnnData observations"*
>            - *"Original name"*: `Sample Characteristic[cell type]`
>            - *"New name"*: `cell_type`
>    - In *"Change field names in AnnData var"*:
>        - {% icon param-repeat %} *"Insert Change field names in AnnData var"*
>            - *"Original name"*: `gene_symbols`
>            - *"New name"*: `Symbol`
>    - *"Gene symbols field in AnnData"*: `Symbol`
>    - In *"Flag genes that start with these names"*:
>        - {% icon param-repeat %} *"Insert Flag genes that start with these names"*
>            - *"Starts with"*: `mt-`
>            - *"Var name"*: `mito`
>
> 4. Rename {% icon galaxy-pencil %} output `Mito-counted AnnData for downstream analysis`
>
{: .hands_on}

And that's all! What's even more exciting about AnnData Operations tool is that it automatically calculates a bunch of metrics, such as log1p_mean_counts, log1p_total_counts, mean_counts, n_cells, n_cells_by_counts, n_counts, pct_dropout_by_counts, total_counts. Amazing, isn't it?

# Conclusion
Now you can use this object as input for the [Filter, Plot, Explore tutorial]({% link topics/single-cell/tutorials/scrna-case_basic-pipeline/tutorial.md %}) and its associated workflow! 

Even though this tutorial was designed specifically to modify the AnnData object to be compatible with the subsequent tutorial, it also shows useful tools that you can use for your own, independent data analysis. You can find the [workflow](https://singlecell.usegalaxy.eu/u/j.jakiela/w/ebi-single-cell-expression-atlas-files-to-anndata) and the [answer key history](https://singlecell.usegalaxy.eu/u/j.jakiela/h/ebi-single-cell-expression-atlas-files-to-anndata-1). However, if you want to use the workflow from this tutorial, you have to keep in mind that different datasets may have different column names. So you have to check them first, and only then you can modify them. 
