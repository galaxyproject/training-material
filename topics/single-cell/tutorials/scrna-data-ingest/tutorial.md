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
    - wee-snufkin
    - hexhowells
    - nomadscientist
  
funding: 
  - elixir-fair-data

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
{: .hands_on}

At this point you might want to do some modifications in the files before downstream analysis. That can include re-formating the cell metadata or changing the names of the column headers, it all depends on your dataset and how you want to perfrom your analysis. It's also fine to transform those files straight away. Now you have the choice to create AnnData object or Seurat object. 

{% include _includes/cyoa-choices.html option1="Scanpy" option2="Seurat" default="Scanpy"
       text="You can choose whether you want to create an AnnData object for Scanpy Analysis or an RDS object for Seurat Analysis. Galaxy has more resources for Scanpy analysis, but sometimes Seurat might have what you want." %}

 <div class="AnnData object" markdown="1">

> <hands-on-title>Create AnnData object</hands-on-title>
>
> {% tool [Scanpy Read10x](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_read_10x/scanpy_read_10x/1.8.1+galaxy0) %} with the following parameters:
>    - *"Expression matrix in sparse matrix format (.mtx)"*: `EBI SCXA Data Retrieval on E-MTAB-6945 matrix.mtx (Raw filtered counts)`
>    - *"Gene table"*:  `EBI SCXA Data Retrieval on E-MTAB-6945 genes.tsv (Raw filtered counts)`
>    - *"Barcode/cell table"*: `EBI SCXA Data Retrieval on E-MTAB-6945 barcodes.tsv (Raw filtered counts)`
>    - *"Cell metadata table"*: `EBI SCXA Data Retrieval on E-MTAB-6945 exp_design.tsv`
{: .hands_on}

</div>

<div class="Seurat object" markdown="1">

> <hands-on-title>Create Seurat object / Loom / SCE </hands-on-title>
>
> {% tool [Seurat Read10x](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_read10x/seurat_read10x/3.2.3+galaxy0) %} with the following parameters:
>    - *"Choose the format of the input"*: `10X-type MTX`
>    - *"Expression matrix in sparse matrix format (.mtx)"*: `EBI SCXA Data Retrieval on E-MTAB-6945 matrix.mtx (Raw filtered counts)`
>    - *"Gene table"*:  `EBI SCXA Data Retrieval on E-MTAB-6945 genes.tsv (Raw filtered counts)`
>    - *"Barcode/cell table"*: `EBI SCXA Data Retrieval on E-MTAB-6945 barcodes.tsv (Raw filtered counts)`
>    - *"Cell Metadata"*: `EBI SCXA Data Retrieval on E-MTAB-6945 exp_design.tsv`
>
> You can now choose if you want to get Seurat object, Loom or Single Cell Experiment by selecting your option in *"Choose the format of the output"*.
{: .hands_on}

</div>

<!---
HCA doesn't work well for other datasets...
-->
## Human Cell Atlas Matrix Downloader

Matrix market format:  matrix mtx, genes tsv, barcodes tsv, exp design tsv
Scnapy read10x to transform those to AnnData object
Flagging by using AnnData Operations works well (only need to check name of the column with gene symbols):
Case sensitive
No parentheses 
Dash important 


# Downsampling FASTQ files

Sometimes it is useful to work on smaller subset of data (especially for teaching / learning purposes). Here is an example of how you can downsample your FASTQ files.
First, let's get some toy data. We just need two FASTQ files - one containing barcodes, the other with transcripts. 

> <hands-on-title>Get toy data</hands-on-title>
>
> You can simply download the files by pasting the links below into "Upload Data" searchbox.
>
>    ```
>    https://zenodo.org/record/4574153/files/SLX-7632.TAAGGCGA.N701.s_1.r_1.fq-400k.fastq
>    https://zenodo.org/record/4574153/files/SLX-7632.TAAGGCGA.N701.s_1.r_2.fq-400k.fastq
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
{: .hands_on}

Funnily enough, those files are already downsampled so that you don't have to wait for too long to download them. We are not going to analyse that data today anyway, it's just for demonstration purposes.

Quick check now which file contains barcodes and which transcripts. If you click on the two datasets, you will see that one has shorter sequences, while the other has longer. It's quite straight-forward to deduce that shorter sequences are barcodes, so let's rename the file `s_1.r_1` as `Barcodes` and file `s_1.r_2` as `Transcripts`.

> <hands-on-title>Rename the files</hands-on-title>
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
{: .hands_on}

Now we will convert the FASTQ files to tabular:

> <hands-on-title> FASTQ to tabular </hands-on-title>
>
> 1. {% tool [FASTQ to Tabular](toolshed.g2.bx.psu.edu/repos/devteam/fastq_to_tabular/fastq_to_tabular/1.1.5) %} with the following parameters:
>    - {% icon param-file %} *"FASTQ file to convert"*: `Barcodes` 
>
> 3. Rename {% icon galaxy-pencil %} the dataset `Barcodes tabular`.
>    
> 4. Repeat the same with `Transcripts` file and rename it as `Transcripts tabular`.
> 
{: .hands_on}

Now let's select the number of the reads we would like to keep. It's totally up to you, we choose 100000 here.

> <hands-on-title> Downsampling </hands-on-title>
>
> 1. {% tool [Select last](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_tail_tool/1.1.0) %} with the following parameters:
>    - {% icon param-file %} *"Text file"*:  `Barcodes tabular` (output of **FASTQ to Tabular** {% icon tool %})
>    - *"Operation"*: `Keep last lines`
>    - *"Number of lines"*: `100000`
>
> 2. Rename {% icon galaxy-pencil %} the dataset `Barcodes cut`.
>    
> 3. Repeat the same with `Transcripts tabular` file and rename it as `Transcripts cut`
> 
{: .hands_on}

All done, now we just need to go back to FASTQ from Tabular again!


> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Tabular to FASTQ](toolshed.g2.bx.psu.edu/repos/devteam/tabular_to_fastq/tabular_to_fastq/1.1.5) %} with the following parameters:
>    - {% icon param-file %} *"Tabular file to convert"*: `Barcodes cut` (output of **Select last** {% icon tool %})
>    - *"Identifier column"*: `c1`
>    - *"Sequence column"*: `c2`
>    - *"Quality column"*: `c3`
>
> 2. Rename {% icon galaxy-pencil %} the dataset `Downsampled barcode/UMI read`.
>    
> 3. Repeat the same with `Transcripts cut` file and rename it as `Downsampled transcript read`
> 
{: .hands_on}

And that's all! Your downsampled data is ready to use. You can check your answers in this [example history](https://usegalaxy.eu/u/j.jakiela/h/how-to-downsample-fastq-files) or if you want to accelerate this process, feel free to use the [workflow](https://singlecell.usegalaxy.eu/u/j.jakiela/w/workflow-constructed-from-history-copy-of-cs1generating-a-single-cell-matrix-using-alevin---how-to-downsample) next time!


# Format conversion

In Galaxy Toolshed there is a wonderful tool called {% tool [SCEasy](toolshed.g2.bx.psu.edu/repos/iuc/sceasy_convert/sceasy_convert/0.0.7+galaxy1) %} which allows to convert between common single cell formats, such as:
- AnnData to CellDataSet (CDS)
- AnnData to Seurat
- Loom to AnnData
- Loom to SingleCellExperiment (SCE)
- SingleCellExperiment (SCE) to AnnData
- SingleCellExperiment (SCE) to Loom
- Seurat to AnnData
- Seurat to SingleCellExperiment (SCE)

> <warning-title>Two SCEasy tools</warning-title>
>  
> The updated SCEasy tool is called [**SCEasy Converter**](toolshed.g2.bx.psu.edu/repos/iuc/sceasy_convert/sceasy_convert/0.0.7+galaxy1) and it's only available on usegalaxy.eu. The second tool is called [**SCEasy convert**](toolshed.g2.bx.psu.edu/repos/ebi-gxa/sceasy_convert/sceasy_convert/0.0.5+galaxy1) and it works on usegalaxy.org, however has limited conversion options.
> 
{: .warning}

However, sometimes it is useful to know how to do this conversion manually or at least to know how it all works. Therefore, below are some examples showing how to do it. 

## AnnData -> Seurat

Let's get an AnnData object that we can further work on. It's the object used in [previous tutorials]({% link topics/single-cell/tutorials/scrna-case_basic-pipeline/tutorial.md %}), so check it out if you're curious. 

> <hands-on-title>Get toy data</hands-on-title>
>
> You can simply download the files by pasting the links below into "Upload Data" searchbox.
>
>    ```
>    https://zenodo.org/record/7053673/files/Mito-counted_AnnData
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
{: .hands_on}


First, we will extract observations and the full matrix from our AnnData.

> <hands-on-title> Inspect AnnData </hands-on-title>
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy1) %} with the following parameters:
>    - *"Annotated data matrix"*: `Mito-counted_AnnData`
>    - *"What to inspect?"*: `Key-indexed observations annotation (obs)`
>   
> 2. Rename {% icon galaxy-pencil %} the output `Observations`.
>
> 3. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy1) %} with the following parameters:
>    - *"Annotated data matrix"*: `Mito-counted_AnnData`
>    - *"What to inspect?"*: `The full data matrix`
>
> 4. Rename {% icon galaxy-pencil %} the output `Matrix`.
>    
{: .hands_on}

> <question-title></question-title>
>
> What are the rows, and what are the columns in the retrieved Matrix?
>
> > <solution-title></solution-title>
> >
> > If you just click on the `Matrix` dataset, you will see a preview, showing barcodes in the first column, while genes in the first row.
> > 
> {: .solution}
>
{: .question}

To proceed with the conversion, we must have the matrix where the genes are listed in the first column while all the barcodes should be in the first row. Therefore, we need to transpose the current matrix.

> <hands-on-title> Transpose the matrix </hands-on-title>
>
> 1. {% tool [Transpose](toolshed.g2.bx.psu.edu/repos/iuc/datamash_transpose/datamash_transpose/1.8+galaxy0) %} with the following parameters:
>    - *"Input tabular dataset"*: `Matrix`
>
{: .hands_on}

And now we are ready to input that data to **DropletUtils** tool.

> <hands-on-title> DropletUtils </hands-on-title>
>
> 1. {% tool [DropletUtils](toolshed.g2.bx.psu.edu/repos/iuc/dropletutils/dropletutils/1.10.0+galaxy2) %} with the following parameters:
>    - *"Format for the input matrix"*: `Tabular`
>    - *"Count Data"*: output of **Transpose** {% icon tool %}
>    - *"Operation"*: `Filter for Barcodes`
>        - *"Method"*: `DefaultDrops`
>            - *"Expected Number of Cells"*: `31178`
>            - *"Upper Quantile"*: `1.0`
>            - *"Lower Proportion"*: `0.0`
>        - *"Format for output matrices"*: `Bundled (barcodes.tsv, genes.tsv, matrix.mtx)`
>        - *"Random Seed"*: `100`
>
{: .hands_on}

Finally, let's combine those files that we have just generated and turn them into the Seurat object!

> <hands-on-title> Create Seurat object </hands-on-title>
>
> 1. {% tool [Seurat Read10x](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_read10x/seurat_read10x/3.2.3+galaxy0) %} with the following parameters:
>    - *"Choose the format of the input"*: `10X-type MTX`
>    - *"Expression matrix in sparse matrix format (.mtx)"*: `DropletUtils 10X Matrices`
>    - *"Gene table"*: `DropletUtils 10X Genes`
>    - *"Barcode/cell table"*: `DropletUtils 10X Barcodes`
>    - *"Cell Metadata"*: `Observations`
>    - *"Choose the format of the output"*: `RDS with a Seurat object`
>
> 2. Rename {% icon galaxy-pencil %} the output `Converted Seurat object`.
>    
{: .hands_on}

As usual, you can check the [example history](https://usegalaxy.eu/u/j.jakiela/h/anndata---seurat) and the dedicated [workflow](https://usegalaxy.eu/u/j.jakiela/w/anndata---seurat-conversion).


## AnnData -> SingleCellExperiment (SCE)

We will work on the same AnnData object so if you create a new history for this exercise, you can either get this file from Zenodo again or just copy this dataset from the previous history. 

> <hands-on-title>Get toy data</hands-on-title>
>
> You can simply download the files by pasting the links below into "Upload Data" searchbox.
>
>    ```
>    https://zenodo.org/record/7053673/files/Mito-counted_AnnData
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
{: .hands_on}

{% snippet faqs/galaxy/histories_copy_dataset.md %}

First, we will extract observations and the full matrix from our AnnData.

> <hands-on-title> Inspect AnnData </hands-on-title>
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy1) %} with the following parameters:
>    - *"Annotated data matrix"*: `Mito-counted_AnnData`
>    - *"What to inspect?"*: `Key-indexed observations annotation (obs)`
>   
> 2. Rename {% icon galaxy-pencil %} the output `Observations`.
>
> 3. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy1) %} with the following parameters:
>    - *"Annotated data matrix"*: `Mito-counted_AnnData`
>    - *"What to inspect?"*: `The full data matrix`
>
> 4. Rename {% icon galaxy-pencil %} the output `Matrix`.
>    
{: .hands_on}

> <question-title></question-title>
>
> What are the rows, and what are the columns in the retrieved Matrix?
>
> > <solution-title></solution-title>
> >
> > If you just click on the `Matrix` dataset, you will see a preview, showing barcodes in the first column, while genes in the first row.
> > 
> {: .solution}
>
{: .question}

To proceed with the conversion, we must have the matrix where the genes are listed in the first column while all the barcodes should be in the first row. Therefore, we need to transpose the current matrix.

> <hands-on-title> Transpose the matrix </hands-on-title>
>
> 1. {% tool [Transpose](toolshed.g2.bx.psu.edu/repos/iuc/datamash_transpose/datamash_transpose/1.8+galaxy0) %} with the following parameters:
>    - *"Input tabular dataset"*: `Matrix`
>
{: .hands_on}

And now we are ready to input that data to **DropletUtils** tool.

> <hands-on-title> DropletUtils </hands-on-title>
>
> 1. {% tool [DropletUtils](toolshed.g2.bx.psu.edu/repos/iuc/dropletutils/dropletutils/1.10.0+galaxy2) %} with the following parameters:
>    - *"Format for the input matrix"*: `Tabular`
>    - *"Count Data"*: output of **Transpose** {% icon tool %}
>    - *"Operation"*: `Filter for Barcodes`
>        - *"Method"*: `DefaultDrops`
>            - *"Expected Number of Cells"*: `31178`
>            - *"Upper Quantile"*: `1.0`
>            - *"Lower Proportion"*: `0.0`
>        - *"Format for output matrices"*: `Bundled (barcodes.tsv, genes.tsv, matrix.mtx)`
>        - *"Random Seed"*: `100`
>
{: .hands_on}

Finally, let's combine those files that we have just generated and turn them into the SingleCellExperiment!

> <hands-on-title> Create SCE object </hands-on-title>
>
> 1. {% tool [DropletUtils Read10x](toolshed.g2.bx.psu.edu/repos/ebi-gxa/dropletutils_read_10x/dropletutils_read_10x/1.0.4+galaxy0) %} with the following parameters:
>    - *"Expression matrix in sparse matrix format (.mtx)"*: `DropletUtils 10X Matrices`
>    - *"Gene table"*: `DropletUtils 10X Genes`
>    - *"Barcode/cell table"*: `DropletUtils 10X Barcodes`
>    - *"Should metadata file be added?"*: {% icon param-toggle %} `Yes`
>        - *"Metadata file"*: `Observations`
>        - *"Cell ID column"*: `index`
>     
> 2. Rename {% icon galaxy-pencil %} the output `Converted SCE object`.
>
{: .hands_on}

As usual, you can check the [example history](https://usegalaxy.eu/u/j.jakiela/h/anndata---singlecellexperiment-sce) and the dedicated [workflow](https://usegalaxy.eu/u/j.jakiela/w/anndata-to-singlecellexperiment-sce-conversion).


## Anndata -> Cell Data Set (CDS)

Cell Data Set (CDS) format is usually used when working with a package called Monocle3 and in the [corresponding tutorial]({% link /topics/single-cell/tutorials/scrna-case_monocle3-trajectories/tutorial.md %}) it was shown how to transform annotated AnnData to CDS object. Please note that depending on your dataset and pre-processing, you might need to refer to the method used in the mentioned tutorial which uses both annotated and unprocessed matrices. However, here we will just show the main principle of that conversion, so we will continue working on previously used dataset, so you can copy it from your history.


> <hands-on-title>Get toy data</hands-on-title>
>
> You can simply download the files by pasting the links below into "Upload Data" searchbox.
>
>    ```
>    https://zenodo.org/record/7053673/files/Mito-counted_AnnData
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
{: .hands_on}

{% snippet faqs/galaxy/histories_copy_dataset.md %}

