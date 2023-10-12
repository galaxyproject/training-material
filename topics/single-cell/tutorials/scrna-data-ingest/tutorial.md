---
layout: tutorial_hands_on
subtopic: datamanipulation
priority: 3
title: Single cell data ingestion and formats conversion
questions:
- What are the most popular single cell datatypes?
- What if the format of my files is different than that used in the tutorial I want to follow?
- Where should I start the analysis depending on the format of my data?
- How to ingest data into Galaxy?
- How to convert between the formats?
objectives:
- You will get to know single cell files formats.
- You will import single cell data to Galaxy using different methods.
- You will manipulate the metadata and matrix files.
- You will perform conversions between the most common single cell formats.
time_estimation: 1H
key_points:
- Single cell data from different sources may have unfamiliar formatting and thus may require different way of ingesting it into Galaxy.
- There are many ways of importing single cell files into Galaxy and converting between single cell formats. 
contributions:
  authorship:
    - hexhowells
    - wee-snufkin
    - nomadscientist

follow_up_training:
  -
    type: "internal"
    topic_name: single-cell

tags:
  - single-cell
  - data-management
  - transcriptomics
---

# Introduction 

You finally decided to analyse some single cell data, you got your files either from the lab or publicly available sources, you opened the first tutorial available on Galaxy Training Network and... you hit the wall - the format of your files is not compatible with the one used in tutorial! Have you been there? 
This tutorial was created to help you overcome that problem and ensure data interoperability in single cell analysis. Once you get your data into Galaxy in the right format, that's already 50% of success. Additionally, by using format conversion, you will be able to use different packages presented in tutorials that may require different datatypes. 


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Single cell datatypes

To start with, here are the most common formats and datatypes that you might come across if you work with single cell data:
- HDF5 - Hierarchical Data Format - can store datasets and groups. A dataset is a  a multidimensional array of data elements, together with supporting metadata, and a group is a structure for organizing objects in an HDF5 file. It allows for storing both the count matrices and all metadata in a single file rather than having separate features, barcodes and matrix files. 
- AnnData objects - [anndata](link https://anndata.readthedocs.io/en/latest/) is a Python package for handling annotated data matrices. In Galaxy, you'll see AnnData objects in **h5ad** format, which is based on the standard HDF5 (h5) format. There are lots of Python tools that work with this format, such as Scanpy, MUON, Cell Oracle, SquidPy, etc. 
- Loom - it is simply an HDF5 file that contains specific groups containing the main matrix as well as row and column attributes and can be read by any language supporting HDF5. [Loompy](https://linnarssonlab.org/loompy/) has been released as a Python API to interact with loom files, and [loomR](https://github.com/mojaveazure/loomR) is its implementation in R. 
- Tabular - simply using TSV, CSV or TXT formats to store expression matrix as well as cells and genes metadata. 
- MTX - it's just a sparse matrix format with genes on the rows and cells on the columns as output by Cell Ranger.
- Zarr - a Python package providing an implementation of compressed, chunked, N-dimensional arrays, designed for use in parallel computing. The Zarr file format offers powerful compression options, supports multiple data store backends, and can read/write your NumPy arrays.
- Seurat objects - a representation of single-cell expression data for R, in Galaxy you might see them in **rdata** format.
- Single Cell Experiment (SCE) object - defines a S4 class for storing data from single-cell experiments and provides a more formalized approach towards construction and accession of data. The S4 system is one of R's systems for object oriented programing. In Galaxy you might see SCE objects in **rdata** format.
- CellDataSet (CDS) object - the main class used by Monocle to hold single cell expression data. In Galaxy you might see CDS objects in **rdata** format.

<!---
TO FURTHER IMPROVE THE TUTORIAL:
include images showing the structure of those files
-->

{% snippet faqs/galaxy/tutorial_mode.md %}


# Data ingestion
As you can see above, there are multiple ways to store single cell data. Therefore, there are also many ways how you can get that data! 

## EBI SCXA Data Retrieval

If you want to use publicly available data, then EBI's [Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/home) is a great place to get resources from. You can search datasets according to various criteria either using search box in **Home** tab or choosing kingdom, experiment collection, technology type (and others) in **Browse experiments** tab. When you find the experiment you are interested in, just click on it and the experiment ID will be displayed in the website URL, as shown below.

![Arrow pointing to the website URL where you can find experiment ID.](../../images/path/exp_id.jpg "Where to find experiment ID on the EBI Single Cell Expression Atlas website.")

Once you know the experiment ID, you can use EBI SCXA Data Retrieval tool in Galaxy! 

> <hands-on-title>Retrieving data from Single Cell Expression Atlas</hands-on-title>
>
> 1. {% tool [EBI SCXA Data Retrieval](toolshed.g2.bx.psu.edu/repos/ebi-gxa/retrieve_scxa/retrieve_scxa/v0.0.2+galaxy2) %} with the following parameters:
>      - *"SC-Atlas experiment accession"*: `E-MTAB-6945`
>      - *"Choose the type of matrix to download"*: `Raw filtered counts`
>
>    Now we need to transform this into an AnnData object.
>
> 2. {% tool [Scanpy Read10x](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_read_10x/scanpy_read_10x/1.8.1+galaxy0) %} with the following parameters:
>    - *"Expression matrix in sparse matrix format (.mtx)"*: `EBI SCXA Data Retrieval on E-MTAB-6945 matrix.mtx (Raw filtered counts)`
>    - *"Gene table"*:  `EBI SCXA Data Retrieval on E-MTAB-6945 genes.tsv (Raw filtered counts)`
>    - *"Barcode/cell table"*: `EBI SCXA Data Retrieval on E-MTAB-6945 barcodes.tsv (Raw filtered counts)`
>    - *"Cell metadata table"*: `EBI SCXA Data Retrieval on E-MTAB-6945 exp_design.tsv`
{: .hands_on}

<!---
CAN FINISH HERE 
or go further with Anndata operations and input into Filter Plot Explore
-->

It's important to note that this matrix is processed somewhat through the SCXA pipeline, which is quite similar to the pre-processing that has been shown in this case study tutorial series, and it contains any and all metadata provided by their pipeline as well as the authors (for instance, more cell or gene annotations). So don't worry if the plots generated using this input method are slightly different! 

Before creating an AnnData object, we need to make a small modification in experimental design table. The dataset contains information about 7 samples N701 – N707), however in the experimental design table (cell metadata) they are just numbered from 1 to 7. The plotting tool that we will going to use later will fail if the entries are integers and not categoricals, so we will change "1" to "N701" and so on. You can simply preview the experimental design dataset and move to the column "Sample Characteristic[individual]" (that's where the information about batch is - don't worry, we will rename the column header later!). Make a note of the number of that column - number 12 - we will need it to change the batch number to batch name. 

> <hands-on-title> Change batch numbers into names </hands-on-title>
>
> 1. Change the datatype of `EBI SCXA Data Retrieval on E-MTAB-6945 exp_design.tsv` to *tabular*:
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular" %}
>
> 2. {% tool [Column Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regexColumn1/1.0.3) %} with the following parameters:
>    - *"using column"*: `c12`
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `1`
>            - *"Replacement"*: `N701`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `2`
>            - *"Replacement"*: `N702`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `3`
>            - *"Replacement"*: `N703`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `4`
>            - *"Replacement"*: `N704`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `5`
>            - *"Replacement"*: `N705`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `6`
>            - *"Replacement"*: `N706`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `7`
>            - *"Replacement"*: `N707`
>
> 3. Rename {% icon galaxy-pencil %} output `Cell metadata`
> 
{: .hands_on}

Now we can create an AnnData object!

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy Read10x](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_read_10x/scanpy_read_10x/1.8.1+galaxy9) %} with the following parameters:
>    - *"Expression matrix in sparse matrix format (.mtx)"*: `EBI SCXA Data Retrieval on E-MTAB-6945 matrix.mtx (Raw filtered counts)`
>    - *"Gene table"*:  `EBI SCXA Data Retrieval on E-MTAB-6945 genes.tsv (Raw filtered counts)`
>    - *"Barcode/cell table"*: `EBI SCXA Data Retrieval on E-MTAB-6945 barcodes.tsv (Raw filtered counts)`
>    - *"Cell metadata table"*: `Cell metadata`
>
> 2. Rename {% icon galaxy-pencil %} output `AnnData object`
> 
{: .hands_on}

Now we will do several modifications within the AnnData object so that you can follow this tutorial despite the other way of getting data! 
We would like to flag mitochondrial genes. They can be identified quite easily since they names start with mt. Since the tool for flagging the mitochondrial genes is case-sensitive, it might be a good idea to check what is the formatting of mitochondrial genes in our dataset.

> <hands-on-title> Check the format of mitochondrial genes names </hands-on-title>
>
> 1. {% tool [Search in textfiles](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_grep_tool/1.1.1) %} with the following parameters:
>    - *"Regular Expression"*: `mt`
>    - *"Output"*: `Highlighted HTML (for easier viewing)`
> 
> 2. Rename {% icon galaxy-pencil %} output `Mito genes check`
>
{: .hands_on}

If you click on that dataset, you will see all the genes containing 'mt' in their name. We can now clearly see that mitochondrial genes in our dataset start with 'mt-'. Keep that in mind, we will use it in a moment!

Speaking about gene names, we will also change the header of the column containing those names from `gene_symbols` to `Symbol`. This edit is only needed to make our AnnData object compatible with this tutorial's workflow. 

As I mentioned at the beginning, we will also change the header of the column storing information about batch. Actually, we will change several other headers as well.

And the good news is that we can do all those steps using only one tool!

> <hands-on-title> Modify AnnData object </hands-on-title>
>
> 1. {% tool [AnnData Operations](toolshed.g2.bx.psu.edu/repos/ebi-gxa/anndata_ops/anndata_ops/1.8.1+galaxy92) %} with the following parameters:
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
> 2. Rename {% icon galaxy-pencil %} output `Mito-counted AnnData for downstream analysis`
>
{: .hands_on}

And that's all! What's even more exciting about AnnData Operations tool is that it automatically calculates a bunch of metrics, such as log1p_mean_counts, log1p_total_counts, mean_counts, n_cells, n_cells_by_counts, n_counts, pct_dropout_by_counts, total_counts. Amazing, isn't it?




## Downsampling 
Sometimes it is useful to work on smaller subset of data (especially for teaching / learning purposes). 

# Format conversion

