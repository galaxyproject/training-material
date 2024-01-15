---
layout: tutorial_hands_on
subtopic: datamanipulation
priority: 3
title: EBI Single Cell Expression Atlas files to AnnData (Scanpy) or Seurat Object | Creating preprocessed dataset for sc-RNA Filter, Plot, Explore tutorial
questions:
- How do I use the EBI Single Cell Expression Atlas?
- How can I reformat and manipulate the downloaded files to create the correct input for downstream analysis?
objectives:
- You will retrieve raw data from the EBI Single Cell Expression Atlas.
- You will manipulate the metadata and matrix files.
- You will combine the metadata and matrix files into an AnnData or Seurat object for downstream analysis.

time_estimation: "15m"
key_points:
- The EMBL-EBI Single-cell Expression Atlas contains high quality datasets.
- Metadata manipulation is key for generating the correctly formatted files.
- To use Scanpy tools, you have to transform your metadata into an AnnData object.
- To use Seurat tools, you have to transform your metadata into a Seurat object.
contributions:
  authorship:
    - wee-snufkin
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
  - data import
  - data management
---

# Introduction

Public single cell datasets seem to accumulate by the second. Well annotated, quality datasets are slightly trickier to find. which is why projects like the [Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/home) (SCXA) exist - to curate datasets for public use. Here, we will guide you through transforming data imported from the SCXA repository into the input file required for the [Filter, Plot, Explore tutorial]({% link topics/single-cell/tutorials/scrna-case_basic-pipeline/tutorial.md %}).

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Getting data from the Single Cell Expression Atlas

Galaxy has a specific tool for importing data from the SCXA ({% cite Moreno2020.04.08.032698 %}), which combines all the preprocessing steps shown in [the previous tutorial]({% link topics/single-cell/tutorials/scrna-case_alevin/tutorial.md %}) into one! For this tutorial, the dataset can be seen [at the EBI](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-6945/downloads) with experiment ID of `E-MTAB-6945`.

> <hands-on-title>Retrieving data from Single Cell Expression Atlas</hands-on-title>
>
> 1. {% tool [EBI SCXA Data Retrieval](toolshed.g2.bx.psu.edu/repos/ebi-gxa/retrieve_scxa/retrieve_scxa/v0.0.2+galaxy2) %} with the following parameters:
>      - *"SC-Atlas experiment accession"*: `E-MTAB-6945`
>      - *"Choose the type of matrix to download"*: `Raw filtered counts`
>
{: .hands_on}

It's important to note that this matrix is processed somewhat through the SCXA pipeline, which is quite similar to the pre-processing that has been shown in this case study tutorial series. The resultant datasets contain any and all metadata provided by the SCXA pipeline as well as the metadata contributed by the original authors (for instance, more cell or gene annotations). So while the AnnData object generated at the end of this tutorial will be similar to that generated using the [Alevin workflows]({% link topics/single-cell/tutorials/scrna-case_alevin/tutorial.md %}) on the original FASTQ files, some of the metadata will be slightly different. Relevant results and interpretation will not change, however!

## Examine the imported files

> <question-title></question-title>
>
> 1. What format has this tool imported?
>
> > <tip-title></tip-title>
> >
> > Selecting the title of each resultant dataset will expand the dataset in the Galaxy history.
> {: .tip}
> > <solution-title></solution-title>
> >
> > Matrix Market Format! We can tell this because our first file helpfully says `MatrixMarket` in the first line.
> > ![Green box containing first output, the matrix.mtx file. Columns are labelled 14458, 5218, and 5308559](../../tutorials/EBI-retrieval/images/matrix-output.png "Matrix Market Output")
> >
> > This {% icon param-file %} `matrix.mtx` file, in Matrix Market format, contains a column referring to each gene (column 1), to each cell (column 2), and the expression values themselves in the final column. To be useful, then, we need to know which genes and cells the numbers are referring to. That's why this format comes with two more files.
> >
> > ![Green box containing second output, the genes.tsv file. The first column contains EnsemblIDs such as ENSMUSG######, while the second column contains gene names. There are 14,457 lines.](../../tutorials/EBI-retrieval/images/genes-output.png "Genes Output")
> >
> > The {% icon param-file %} `genes.tsv` file lists each EnsemblID and its gene name. The lines (14,457) corresponds with the 14458 in the Matrix file...but the 14458 contains a header, so that's why it has one more than the genes file!
> >
> > ![Green box containing third output, the barcodes.tsv file. The file consists of 5,217 lines and a single column containing the cell barcode, variations of ERR2704656-AAAACACTCTGA.](../../tutorials/EBI-retrieval/images/cells-output.png "Cells Output")
> >
> > The {% icon param-file %} `barcodes.tsv` file lists each barcode. The lines (5,217) again correspond with the 5,218 lines in the Matrix file...which adds in the header again!
> >
> > ![Green box containing fourth output, the exp_design.tsv file. The file consists of 5,218 lines and numerous columns starting with 'Assay' and 'Sample Characteristic'.](../../tutorials/EBI-retrieval/images/cell-metadata.png "Experimental Design")
> >
> > Finally, and helpfully, the tool also includes cell metadata where the `Assay` column corresponds with the barcodes in the {% icon param-file %} `barcodes.tsv` file. While this is not a required file to create an AnnData object from the three matrix market files, it is extremely necessary for actually interpreting the data. Imagine not knowing which barcodes came from which sample!
> >
> {: .solution}
{: .question}


# Metadata manipulation

Before creating an AnnData object, we need to make a small modification in experimental design table. The dataset contains information about the 7 experimental samples (N701 – N707). However, in the {% icon param-file %} `exp_design.tsv` dataset, which contains the cell metadata, these samples are just numbered from 1 to 7.

> <details-title>How did we know this?</details-title>
>
> You can preview this column in the the {% icon param-file %} `exp_design.tsv` dataset by selecting the {% icon galaxy-eye %} in the {% icon galaxy-history %} Galaxy history. If you scroll to the right, and move to the column `Sample Characteristic[individual]`, you will find the batch information. Don't worry, we're about to rename and reformat this whole dataset to more useful titles. Make a note of the number of that column - number 12 - as we will need it to change the batch number to a batch name shortly.
>
{: .details}

The plotting tool that we will use later will fail if the entries are integers and not categorical values, so we will change `1` to `N01` and so on.

> <hands-on-title> Change batch numbers into names </hands-on-title>
>
> 1. Change the datatype of {% icon param-file %} `EBI SCXA Data Retrieval on E-MTAB-6945 exp_design.tsv` to `tabular`:
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular" %}
>
> 2. {% tool [Column Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regexColumn1/1.0.3) %} with the following parameters:
>    - {% icon param-file %} *"Select cells from"*: `EBI SCXA Data Retrieval on E-MTAB-6945 exp_design.tsv`
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
{: .hands_on}

While we're renaming things, let's also fix our titles.

> <hands-on-title> Change cell metadata titles </hands-on-title>
>
> 1. {% tool [Replace parts of text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_find_and_replace/1.1.4) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: output from **Column Regex and Replace** {% icon tool %}
>    - In *"Find and Replace"*:
>        - *"Find pattern"*: `"Sample Characteristic[genotype]"`
>        - *"Replace with"*: `genotype`
>    - {% icon param-repeat %} *"Insert Find and Replace"*
>        - *"Find pattern"*: `"Sample Characteristic[individual]"`
>        - *"Replace with"*: `batch`
>    - {% icon param-repeat %} *"Insert Find and Replace"*
>        - *"Find pattern"*: `"Sample Characteristic[sex]"`
>        - *"Replace with"*: `sex`
>    - {% icon param-repeat %} *"Insert Find and Replace"*
>        - *"Find pattern"*: `"Sample Characteristic[cell type]"`
>        - *"Replace with"*: `cell_type`
> 2. Rename {% icon galaxy-pencil %} output `Cell metadata`
{: .hands_on}

# Check mitochondrial gene name format

We might like to flag mitochondrial genes. They can be identified quite easily since - depending on the species and formatting convention - their names often start with `mt`. Since tools for flagging mitochondrial genes are often case-sensitive, it might be a good idea to check the formatting of the mitochondrial genes in our dataset.

> <hands-on-title> Check the format of mitochondrial genes names </hands-on-title>
>
> 1. {% tool [Search in textfiles](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_grep_tool/1.1.1) %} with the following parameters:
>    - {% icon param-file %}  *"Select lines from"*: `EBI SCXA Data Retrieval on E-MTAB-6945 genes.tsv (Raw filtered counts)`
>    - *"that"*: `Match`
>    - *"Regular Expression"*: `mt`
>    - *"Match type"*: `case insensitive`
>    - *"Output"*: `Highlighted HTML (for easier viewing)`
> 3. Rename {% icon galaxy-pencil %} output `Mito genes check`
>
{: .hands_on}

If you click on that dataset, you will see all the genes containing `mt` in their name. We can now clearly see that mitochondrial genes in our dataset start with `mt-`. Keep that in mind, we might use it in a moment!

Now we can create our single cell object!

{% include _includes/cyoa-choices.html option1="Scanpy" option2="Seurat" default="Scanpy"
       text="You can choose whether you want to create an AnnData object for Scanpy Analysis or an RDS object for Seurat Analysis. Galaxy has more resources for Scanpy analysis, but sometimes Seurat might have what you want. The two packages are constantly trying to outdo the other! It often depends on what is more 'standard' in your work environment!" %}

 <div class="Scanpy" markdown="1">


# Creating the AnnData object

We will do several modifications within the AnnData object so that you can follow [the next tutorial]({% link topics/single-cell/tutorials/scrna-case_basic-pipeline/tutorial.md %}).

> <hands-on-title> Create the AnnData Object </hands-on-title>
>
> 1. {% tool [Scanpy Read10x](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_read_10x/scanpy_read_10x/1.8.1+galaxy9) %}
> 2. Make sure you are using version **1.8.1+galaxy9** of the tool (change by clicking on {% icon tool-versions %} Versions button):
>   ![List of available tool versions shown when clicking on the 'Versions' button on the top of the page.](../../images/scrna-casestudy/version.png "How to change the version of the tool")
>
> 3. Set the following parameters:
>    - {% icon param-file %} *"Expression matrix in sparse matrix format (.mtx)"*: `EBI SCXA Data Retrieval on E-MTAB-6945 matrix.mtx (Raw filtered counts)`
>    - *"Gene table"*:  `EBI SCXA Data Retrieval on E-MTAB-6945 genes.tsv (Raw filtered counts)`
>    - *"Barcode/cell table"*: `EBI SCXA Data Retrieval on E-MTAB-6945 barcodes.tsv (Raw filtered counts)`
>    - *"Cell metadata table"*: `Cell metadata`
> 4. Rename {% icon galaxy-pencil %} output `AnnData object`
>
{: .hands_on}

# AnnData manipulation

We will now change the header of the column containing gene names from `gene_symbols` to `Symbol`. This edit is only needed to make our AnnData object compatible with this tutorial's workflow. We will also flag the mitochondrial genes.

And the good news is that we can do both those steps using only one tool!

> <hands-on-title> Modify AnnData object </hands-on-title>
>
> 1. {% tool [AnnData Operations](toolshed.g2.bx.psu.edu/repos/ebi-gxa/anndata_ops/anndata_ops/1.8.1+galaxy92) %}
> 2. Make sure you are using version **1.8.1+galaxy92** of the tool (change by clicking on {% icon tool-versions %} Versions button)
> 3. Set the following parameters:
>    - {% icon param-file %} In *"Input object in hdf5 AnnData format"*: `AnnData object`
>    - In *"Change field names in AnnData var"*:
>        - {% icon param-repeat %} *"Insert Change field names in AnnData var"*
>            - *"Original name"*: `gene_symbols`
>            - *"New name"*: `Symbol`
>    - *"Gene symbols field in AnnData"*: `Symbol`
>    - In *"Flag genes that start with these names"*:
>        - {% icon param-repeat %} *"Insert Flag genes that start with these names"*
>            - *"Starts with"*: `mt-`
>            - *"Var name"*: `mito`
> 4. Rename {% icon galaxy-pencil %} output `Mito-counted AnnData for downstream analysis`
>
{: .hands_on}

And that's all! What's even more exciting about the {% icon tool %}) **AnnData Operations** tool is that it automatically calculates a bunch of metrics, such as `log1p_mean_counts`, `log1p_total_counts`, `mean_counts`, `n_cells`, `n_cells_by_counts`, `n_counts`, `pct_dropout_by_counts`, and `total_counts`. Amazing, isn't it?

# Conclusion
Now you can use this object as input for the [Filter, Plot, Explore tutorial]({% link topics/single-cell/tutorials/scrna-case_basic-pipeline/tutorial.md %}) and its associated workflow!

Even though this tutorial was designed specifically to modify the AnnData object to be compatible with the subsequent tutorial, it also shows useful tools that you can use for your own, independent data analysis. You can find the [workflow](https://singlecell.usegalaxy.eu/u/j.jakiela/w/workflow-constructed-from-history-ebi-workflow) and the [answer key history](https://singlecell.usegalaxy.eu/u/j.jakiela/h/ebi-scxa---anndata-scanpy-or-seurat-object-1). However, if you want to use the workflow from this tutorial, you have to keep in mind that different datasets may have different column names. So you have to check them first, and only then you can modify them.

</div>

<div class="Seurat" markdown="1">

# Creating the Seurat Object

> <hands-on-title> Create the Seurat Object </hands-on-title>
>
> 1. {% tool [Seurat Read10x](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_read10x/seurat_read10x/3.2.3+galaxy0) %}
> 2. Set the following parameters:
>    - {% icon param-file %} *"Expression matrix in sparse matrix format (.mtx)"*: `EBI SCXA Data Retrieval on E-MTAB-6945 matrix.mtx (Raw filtered counts)`
>    - *"Gene table"*:  `EBI SCXA Data Retrieval on E-MTAB-6945 genes.tsv (Raw filtered counts)`
>    - *"Barcode/cell table"*: `EBI SCXA Data Retrieval on E-MTAB-6945 barcodes.tsv (Raw filtered counts)`
>    - *"Cell metadata"*: `Cell metadata`
> 3. Rename {% icon galaxy-pencil %} output `Seurat object`
>
{: .hands_on}

# Conclusion
And you're there! You now have a usable Seurat object for analysis with Seurat tools in your history! {% icon congratulations %} Congrats!

</div>
