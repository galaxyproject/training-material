---
layout: tutorial_hands_on
subtopic: deconvo
priority: 2
title: Creating the single-cell RNA-seq reference dataset for deconvolution
questions:
- Where can I find good quality scRNA-seq reference datasets?
- How can I reformat and manipulate these downloads to create the right format for MuSiC?
objectives:
- You will retrieve raw data from the EMBL-EBI Single cell expression atlas.
- You will manipulate the metadata and matrix files.
- You will combine the metadata and matrix files into an ESet object for MuSiC deconvolution.
- You will create multiple ESet objects - both combined and separated out by disease phenotype for your single cell reference.
time_estimation: 1H
key_points:
- The EMBL-EBI Single-cell expression atlas contains high quality datasets.
- Metadata manipulation is key for generating the correctly formatted resource.
contributions:
  authorship:
    - nomadscientist
    - mtekman
  testing:
    - MarisaJL

requirements:
  -
    type: "internal"
    topic_name: single-cell
    tutorials:
      - bulk-music

follow_up_training:
  -
    type: "internal"
    topic_name: single-cell
    tutorials:
        - bulk-music-3-preparebulk

tags:
  - single-cell
  - human
  - deconvolution
  - bulk
  - transcriptomics
---


# Introduction


After completing the [MuSiC](https://xuranw.github.io/MuSiC/articles/MuSiC.html) {% cite wang2019bulk %} deconvolution tutorial, you are hopefully excited to apply this analysis to data of your choice. Annoyingly, getting data in the right format is often what prevents us from being able to successfully apply analyses. This tutorial is all about reformatting a raw scRNA-seq dataset pulled from a public resource (the EMBL-EBI single cell expression atlas {% cite Moreno2021 %}. Let's get started!


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Metadata Manipulation

First, we will tackle the metadata. We are roughly following the same concept as in the previous bulk deconvolution tutorial, by comparing human pancreas data across a disease variable (type II diabetes vs healthy), but using public datasets to do it.

## Find the data
We explored the [single cell expression atlas](https://www.ebi.ac.uk/gxa/sc/experiments), browsing experiments in order to find a pancreas dataset ({% cite segerstolpe2016single %}). You can [explore this dataset](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-5061/results/tsne) using their browser. These cells come from 6 healthy individuals and 4 individuals with Type II diabetes, so we will create reference Expression Set objects for the total as well as separating out by phenotype, as you may have reason to do this in your analysis (or you may not!).

{% snippet faqs/galaxy/tutorial_mode.md %}

Galaxy has a specific tool for ingesting data from the Single cell expression atlas, so there are no uploads for this tutorial.

> <hands-on-title>Data retrieval</hands-on-title>
>
> 1. {% tool [EBI SCXA Data Retrieval](toolshed.g2.bx.psu.edu/repos/ebi-gxa/retrieve_scxa/retrieve_scxa/v0.0.2+galaxy2) %} with the following parameters:
>    - *"SC-Atlas experiment accession"*: `E-MTAB-5061`
>
> Data management is going to be key in this analysis, so trust me now to start adding tags.
> 5. Add to the **EBI SCXA Data Retrieval on E-MTAB-5061 exp_design.tsv** file the following tags: `#ebi #metadata #singlecell`
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
{: .hands_on}

This tool will retrieve four files: a barcodes list, a genes list, an experimental design file, and a matrix market format (where columns refer to genes, cells, and quantities). We (mostly) only need the experimental design file, but keep in mind this will have data on all the cells reported by the authors.

> <question-title></question-title>
>
> 1. How many cells are in the sample?
> 2. How many cells were submitted by the authors?
>
> > <div id="solution-1" class="box-title"><button type="button" aria-controls="solution-1-contents" aria-expanded="true" aria-label="Toggle solution box: "><i class="far fa-eye" aria-hidden="true"></i><span class="visually-hidden"></span> Solution<span role="button" class="fold-unfold fa fa-minus-square"></span></button></div>
> >
> > 1. If you select the {% icon param-file %} **barcodes.tsv** file, you'll find that it contains 2914 lines - this corresponds to the 2914 cells, because each cell is given a barcode.
> > 2. The nature of public repositories is that they ingest data from many places, which means they usually apply a uniform analysis to samples. This rarely means they yield the same cell numbers as the original authors. If you check the {% icon param-file %} **exp_design.tsv** file, which refers to the data submitted by the authors, you'll find it contains 3514 lines - referring to 3514 cells submitted by authors. It's important to know that (currently) these files differ.
> >
> {: .solution}
{: .question}

## Prepare the experimental design file

Let's get rid of a bunch of repetitive columns in the metadata we don't need. You can find out what each column is by inspecting the dataset {% icon galaxy-eye %} in the history window.

> <hands-on-title>Cutting necessary metadata columns</hands-on-title>
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1,c4,c6,c8,c10,c14,c20,c24,c26,c30,c32,c34`
>    - {% icon param-file %} *"From"*: `design_tsv` (output of **EBI SCXA Data Retrieval** {% icon tool %})
>
{: .hands_on}

You can inspect the dataset {% icon galaxy-eye %} to see that it's full of annoying "" everywhere, and overly long descriptions of each columns.

![Columns in Galaxy history window contain "Assay" and "" around every word or ID](../../images/bulk-music/annoying_ebimetadata.png "Annoying metadata")

Now, there might be a better way to do this in Galaxy (or you might consider downloading the file locally and changing it in a spreadsheet application or something), but this is what will work to reformat all that annoying text.

> <hands-on-title>Reformatting the metadata</hands-on-title>
>
> 1. {% tool [Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.2) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `out_file1` (output of **Cut** {% icon tool %})
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `"Sample Characteristic\[individual\]"`
>            - *"Replacement"*: `Individual`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `"Sample Characteristic\[sex\]"`
>            - *"Replacement"*: `Sex`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `"Sample Characteristic\[age\]"`
>            - *"Replacement"*: `Age`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `"Sample Characteristic\[body mass index\]"`
>            - *"Replacement"*: `BMI`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `kilogram per square meter`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `HbA1c `
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `"Sample Characteristic\[clinical information\]"`
>            - *"Replacement"*: `HbA1c`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `%`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `"Sample Characteristic\[disease\]"`
>            - *"Replacement"*: `Disease`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `"Sample Characteristic\[single cell quality\]"`
>            - *"Replacement"*: `Single cell quality`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `"Sample Characteristic\[submitted single cell quality\]"`
>            - *"Replacement"*: `"Submitted single cell quality"`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `"Factor Value\[inferred cell type - ontology labels\]"`
>            - *"Replacement"*: `Inferred cell type - ontology label`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `"Factor Value\[inferred cell type - authors labels\]"`
>            - *"Replacement"*: `Inferred cell type - author labels`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `""`
>            - *"Replacement"*: 
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `"`
>            - *"Replacement"*: 
>
>    > <comment-title></comment-title>
>    >
>    > What's with the `\` everywhere? That's because the `[]` symbols usually call the code to do something, rather than just read it as a normal character. the `\` prevents this.
>    {: .comment}
>
> 2. Change the datatype to tabular.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular" %}
>
{: .hands_on}

Great, this file is now ready to go! But, it contains all those extra cells that didn't pass filtration with the EBI pipeline and therefore won't exist in the matrix. We need to remove them for future steps to work. We can use our barcodes list to remove the extra cells.

## Prepare the barcodes file

> <hands-on-title>Adding a header</hands-on-title>
>
> 1. {% tool [Add line to file](toolshed.g2.bx.psu.edu/repos/bgruening/add_line_to_file/add_line_to_file/0.1.0) %} with the following parameters:
>    - *"text to add"*: `Cell`
>    - {% icon param-file %} *"input file"*: `barcode_tsv` (output of **EBI SCXA Data Retrieval** {% icon tool %})
> 
> 2. Change the datatype to tabular.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular" %}
> 
>    > <comment-title></comment-title>
>    >
>    > This is an annoying step we have to do to get the right format, otherwise future steps won't work.
>    {: .comment}
>
{: .hands_on}

## Use the barcodes list to filter out cells in the experimental design file

> <hands-on-title>Joining datasets</hands-on-title>
>
> 1. {% tool [Join two Datasets](join1) %} with the following parameters:
>    - {% icon param-file %} *"Join"*: `outfile` (output of **Add line to file** {% icon tool %})
>    - *"using column"*: `c1`
>    - {% icon param-file %} *"with"*: `out_file1` (output of **Regex Find And Replace** {% icon tool %})
>    - *"and column"*: `c1`
>    - *"Fill empty columns"*: `No`
>    - *"Keep the header lines"*: `Yes`
>    > <comment-title></comment-title>
>    >
>    > Make sure that you join the files in the same order as above - put the output of Add line to file in first - otherwise your columns will be in a different order for the next step. Everything will still work, but you would need to change the number of the column you remove using Advanced Cut. 
>    {: .comment}
{: .hands_on}

> <question-title></question-title>
>
> 1. How many cells are now in your table?
> 2. Is your table ready to go?
>
> > <div id="solution-1" class="box-title"><button type="button" aria-controls="solution-1-contents" aria-expanded="true" aria-label="Toggle solution box: "><i class="far fa-eye" aria-hidden="true"></i><span class="visually-hidden"></span> Solution<span role="button" class="fold-unfold fa fa-minus-square"></span></button></div>
> >
> > 1. If you select the output dataset in your history, you will find `2915` lines, corresponding to 2914 cells and a header. Success!
> > 2. Not quite - notice how you have two identical columns `Cell` and `Assay`? Let's get rid of one.
> >
> {: .solution}
{: .question}

> <hands-on-title>Remove duplicate columns</hands-on-title>
>
> 1. {% tool [Advanced Cut](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/1.1.0) %} with the following parameters:
>    - {% icon param-file %} *"File to cut"*: `out_file1` (output of **Join two Datasets** {% icon tool %})
>    - *"Operation"*: `Discard`
>    - *"Cut by"*: `fields`
>        - *"List of Fields"*: `c1`
>        
>    > <comment-title></comment-title>
>    >
>    > Advanced cut works slightly differently in a workflow versus running the tool independently. Independently, there is a list and you can click through the list to note your columns, while in a workflow it appears as a text option and you put each column on a different line. The point is, each number above represents a column, so remove them!
>    {: .comment}
{: .hands_on}

Fantastic! You've completed part 1 - making the single cell metadata file. It should now look like this:

![Columns in the history window of a dataset contain words without any extra symbols or ""](../../images/bulk-music/corrected_ebimetadata.png "Pretty scRNA metadata")

You can use the [workflow for this portion of the tutorial](https://usegalaxy.eu/u/wendi.bacon.training/w/music-deconvolution-data-generation--sc--metadata), and access an [example history](https://usegalaxy.eu/u/wendi.bacon.training/h/music-deconvolution-data-generation--sc--metadata).


# Manipulate the expression matrix

Currently, the matrix data is in a 3-column format common in 10x outputs, where you need the barcodes and the genes files to interpret the matrix. What you actually need is an expression matrix with cells on one axis and genes on another. While we aren't running a Scanpy analysis, we can still use our Scanpy tools to get this format.

## Reformat the matrix

> <hands-on-title>Task description</hands-on-title>
>
> 1. {% tool [Scanpy Read10x](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_read_10x/scanpy_read_10x/1.8.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Expression matrix in sparse matrix format (.mtx)"*: `matrix_mtx` (output of **EBI SCXA Data Retrieval** {% icon tool %})
>    - {% icon param-file %} *"Gene table"*: `genes_tsv` (output of **EBI SCXA Data Retrieval** {% icon tool %})
>    - {% icon param-file %} *"Barcode/cell table"*: `barcode_tsv` (output of **EBI SCXA Data Retrieval** {% icon tool %})
>    - *"Format of output object"*: `AnnData format (h5 for older versions)`
>
> 2. Change the datatype to `h5ad`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="h5ad" %}
>
{: .hands_on}

Now your precious matrix is stored in the 10x AnnData object. Let's retrieve it!

> <hands-on-title>Inspect the matrix</hands-on-title>
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output_h5` (output of **Scanpy Read10x** {% icon tool %})
>    - *"What to inspect?"*: `The full data matrix`
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Which are currently the rows in your matrix, cells or genes?
>
> > <div id="solution-1" class="box-title"><button type="button" aria-controls="solution-1-contents" aria-expanded="true" aria-label="Toggle solution box: "><i class="far fa-eye" aria-hidden="true"></i><span class="visually-hidden"></span> Solution<span role="button" class="fold-unfold fa fa-minus-square"></span></button></div>
> >
> > 1. You may remember from earlier that the sample should have `2914 cells` in it. If you inspect the dataset in your history, you will find that it contains `2915` lines (1 for the header), which means that rows correspond to cells. Unfortunately... that's not what you need.
> >
> {: .solution}
{: .question}

> <hands-on-title>Transpose the matrix</hands-on-title>
>
> 1. {% tool [Transpose](toolshed.g2.bx.psu.edu/repos/iuc/datamash_transpose/datamash_transpose/1.1.0+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Input tabular dataset"*: `X` (output of **Inspect AnnData** {% icon tool %})
>
>
{: .hands_on}

> <question-title></question-title>
>
> 1. How many genes are in your sample?
>
> > <div id="solution-1" class="box-title"><button type="button" aria-controls="solution-1-contents" aria-expanded="true" aria-label="Toggle solution box: "><i class="far fa-eye" aria-hidden="true"></i><span class="visually-hidden"></span> Solution<span role="button" class="fold-unfold fa fa-minus-square"></span></button></div>
> >
> > 1. You should have `30,416` lines in it, meaning your sample has `30,415` genes.
> >
> {: .solution}
{: .question}

## Collapse EnsemblIDs

Ok, real talk here. Technically, the best way of analysing anything is by using the EnsemblIDs for any given RNA transcript, because they are more specific than gene names and also cover more of the transcriptome than our gene names...
But...
As biologists, it's very difficult to interpret ENSIDs. And it's an awful shame to get to the end of the MuSiC deconvolution and have all our plots show sad ENS IDs. So, courtesy of the excellent @mtekman, we steal his workflow to collapse the ENS IDs into gene names.

![First table shows 4 rows with different ENS IDs and a second column with geney symbols and some overlap. Arrow pointing to second table shows 3 rows having collapsed the overlap.](../../images/bulk-music/ensid_collapse.png "Collapsing ENS IDs")

> <hands-on-title>Convert from Ensembl to GeneSymbol using workflow</hands-on-title>
>
> 1. Import this [workflow](https://usegalaxy.eu/u/wendi.bacon.training/w/convert-from-ensembl-to-genesymbol-summing-duplicate-genes).
>
>    {% snippet faqs/galaxy/workflows_import.md %}
>
> 2. Run the workflow on your sample with the following parameters:
>    - *"Organism"*: `Human`
>    - {% icon param-file %} *"Expression Matrix (Gene Rows)"*: `output_h5` (output of **Transpose** {% icon tool %})
>
>    {% snippet faqs/galaxy/workflows_run.md %}
>
{: .hands_on}

The output will likely be called **Text transformation** and will look like this:

![Alphabetised gene symbols appear in column one with decimal pointed integers in the following columns corresponding to cells](../../images/bulk-music/ensid_output.png "Output of the ENS ID collapsing workflow")

# Construct Expression Set Objects

We're nearly there! We have three more tasks to do: first, we need to create the expression set object with all the phenotypes combined. Then, we also want to create two separate objects - one for healthy and one for diseased as references.

> <hands-on-title>Creating the combined object</hands-on-title>
>
> 1. {% tool [Construct Expression Set Object](toolshed.g2.bx.psu.edu/repos/bgruening/music_construct_eset/music_construct_eset/0.1.1+galaxy4) %} with the following parameters:
>    - {% icon param-file %} *"Assay Data"*: `out_file` #matrix (output of **Text transformation** {% icon tool %})
>    - {% icon param-file %} *"Phenotype Data"*: `output` (output of **Advanced Cut** {% icon tool %})
>
> 2. Remove the `#metadata #matrix` tags from the output **RData ESet Object**
>
> 3. Add the tag `#combined` to the output **RData ESet Object**
>
>
{: .hands_on}

> <question-title></question-title>
>
> 1. How many genes are in your sample now?
>
> > <div id="solution-1" class="box-title"><button type="button" aria-controls="solution-1-contents" aria-expanded="true" aria-label="Toggle solution box: "><i class="far fa-eye" aria-hidden="true"></i><span class="visually-hidden"></span> Solution<span role="button" class="fold-unfold fa fa-minus-square"></span></button></div>
> >
> > 1. If you select the {% icon galaxy-eye %} of the output **General Info** dataset in the history, you will find it contains 21671 features and 2914 samples, or rather, `21671` genes and `2914` cells. That's a huge reduction in genes thanks to the ENS ID collapsing!
> >
> {: .solution}
{: .question}

> <hands-on-title>Creating the disease-only object</hands-on-title>
>
> 1. {% tool [Manipulate Expression Set Object](toolshed.g2.bx.psu.edu/repos/bgruening/music_manipulate_eset/music_manipulate_eset/0.1.1+galaxy4) %} with the following parameters:
>    - {% icon param-file %} *"Expression Set Dataset"*: `out_rds` (output of **Construct Expression Set Object** {% icon tool %})
>    - *"Concatenate other Expression Set objects?"*: `No`
>    - *"Subset the dataset?"*: `Yes`
>        - *"By"*: `Filter Samples and Genes by Phenotype Values`
>            - In *"Filter Samples by Condition"*:
>                - {% icon param-repeat %} *"Insert Filter Samples by Condition"*
>                    - *"Name of phenotype column"*: `Disease`
>                    - *"List of values in this column to filter for, comma-delimited"*: `type II diabetes mellitus`
>
> 2. Remove the `#combined` tag from the output **RData ESet Object**
>
> 3. Add the tag `#T2D` to the output **RData ESet Object**
>
>
>
{: .hands_on}

You can either re-run this tool or set it up again to create the healthy-only object.

> <hands-on-title>Creating the healthy-only object</hands-on-title>
>
> 1. {% tool [Manipulate Expression Set Object](toolshed.g2.bx.psu.edu/repos/bgruening/music_manipulate_eset/music_manipulate_eset/0.1.1+galaxy4) %} with the following parameters:
>    - {% icon param-file %} *"Expression Set Dataset"*: `out_rds` (output of **Construct Expression Set Object** {% icon tool %})
>    - *"Concatenate other Expression Set objects?"*: `No`
>    - *"Subset the dataset?"*: `Yes`
>        - *"By"*: `Filter Samples and Genes by Phenotype Values`
>            - In *"Filter Samples by Condition"*:
>                - {% icon param-repeat %} *"Insert Filter Samples by Condition"*
>                    - *"Name of phenotype column"*: `Disease`
>                    - *"List of values in this column to filter for, comma-delimited"*: `normal`
>
> 2. Remove the `#combined` tag from the output **RData ESet Object**
>
> 3. Add the tag `#healthy` to the output **RData ESet Object**
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Why are you making a healthy-only and diseased-only reference objects?
>
> > <div id="solution-1" class="box-title"><button type="button" aria-controls="solution-1-contents" aria-expanded="true" aria-label="Toggle solution box: "><i class="far fa-eye" aria-hidden="true"></i><span class="visually-hidden"></span> Solution<span role="button" class="fold-unfold fa fa-minus-square"></span></button></div>
> >
> > 1. We could imagine that the cells will express different transcript levels, but that the deconvolution tools will have to take some sort of average. Perhaps it might be more accurate to infer like from like, i.e. healthy from healthy? Or perhaps that is skewing the data through a more 'supervised' approach. We're not   quite sure, and it likely depends on the biology, so we're covering all our bases by making sure you can do this every way. (We've tested it on our dataset in all the ways and got the same results, so it doesn't make much of a difference as far as we can tell!)
> >
> {: .solution}
{: .question}

# Conclusion

You have successfully performed, essentially, three workflows. You can find the [workflows for generating the ESet object](https://usegalaxy.eu/u/wendi.bacon.training/w/music-deconvolution-data-generation--sc--matrix--eset) and the [answer key history for this entire tutorial](https://usegalaxy.eu/u/wendi.bacon.training/h/music-deconvolution-data-generation--sc--matrix--eset).

![6 boxes in the workflow editor](../../images/bulk-music/scref-metadata.png "Workflow: Manipulating metadata")

![9 boxes including a subworkflow in the workflow editor](../../images/bulk-music/scref-eset.png "Workflow: Creating the ESet single-cell object")

![8 boxes in the workflow editor](../../images/bulk-music/ens-id-collapse.png "Subworkflow: Collapsing the Ensembl IDs into gene names")

With these workflows, you've created three Expression Set objects, capable of running in the MuSiC Compare tutorial. Now you just need the bulk RNA-seq Expression Set objects!

This tutorial is part of the [https://singlecell.usegalaxy.eu](https://singlecell.usegalaxy.eu) portal ({% cite tekman2020single %}).
