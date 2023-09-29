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
This tutorial was created to help you overcome that problem. Once you get your data into Galaxy in the right format, that's already 50% of success. Additionally, by using format conversion, you will be able to use different packages presented in tutorials that may require different datatypes. 


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
- Seurat objects 
- Zarr
- Single Cell Experiment
- CDS 

# Data ingestion

{% snippet faqs/galaxy/tutorial_mode.md %}


# Format conversion

