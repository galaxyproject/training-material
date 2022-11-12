---
layout: tutorial_hands_on

title: Creating the single-cell RNA-seq reference dataset
questions:
- Where can I find good quality scRNA-seq reference datasets?
- How can I reformat and manipulate these downloads to create the right format for MuSiC?
objectives:
- You will retrieve raw data from the EMBL-EBI Single cell expression atlas.
- You will manipulate the metadata and matrix files.
- You will combine the metadata and matrix files into an ESet object for MuSiC deconvolution.
- You will create multiple ESet objects - both combined and separated out by disease phenotype for your single cell reference.
time_estimation: 2H
key_points:
- The EMBL-EBI Single-cell expression atlas contains high quality datasets.
- Metadata manipulation is key for generating the correctly formatted resource.
contributors:
- nomadscientist
- mtekman

follow_up_training:
  -
    type: "internal"
    topic_name: single-cell
    tutorials:
        - bulk-music-3-preparebulk

requirements:
  -
    type: "internal"
    topic_name: single-cell
    tutorials:
      - bulk-music

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

After completing the MuSiC {% cite wang2019bulk %} deconvolution tutorial, you are hopefully excited to apply this analysis to data of your choice. Annoyingly, getting data in the right format is often what prevents us from being able to successfully apply analyses. This tutorial is all about reformatting a raw dataset pulled from a public resource (the EMBL-EBI single cell expression atlas {% cite Moreno2021 %}.  [MuSiC](https://xuranw.github.io/MuSiC/articles/MuSiC.html) or published article {% cite wang2019bulk %}. Let's get started!


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Metadata Manipulation

First, we will tackle the metadata. We are roughly following the same concept as in the previous bulk deconvolution tutorial, by comparing human pancreas data across a disease variable (type II diabetes vs healthy), but using public datasets to do it. 

## Finding data 
We explored the [single cell expression atlas](https://www.ebi.ac.uk/gxa/sc/experiments), browsing experiments in order to find a pancreas dataset: {% cite Segerstolpe2016 %}. You can [explore this dataset here](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-5061/results/tsne) using their browser. These cells come from 6 healthy individuals and 4 individuals with Type II diabetes, so we will create reference Expression Set objects for the total as well as separating out by phenotype, as you may have reason to do this in your analysis (or you may not!).

## Get data

Galaxy has a specific tool for ingesting data from the Single cell expression atlas, so there are no uploads for this tutorial.

> ### {% icon hands_on %} Hands-on: Data retrieval
>
> 1. {% tool [EBI SCXA Data Retrieval](retrieve_scxa/v0.0.2+galaxy2) %} with the following parameters:
>    - *"SC-Atlas experiment accession"*: `E-MTAB-5061`
>
> Data management is going to be key in this analysis, so trust me now to start adding tags.
> 5. Add to the **EBI SCXA Data Retrieval on E-MTAB-5061 exp_design.tsv** file the following tags: #ebi #metadata #singlecell
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

Let's get rid of a bunch of repetitive columns in the metadata we don't need. You can find out what each columns is by inspecting the dataset {% icon galaxy-eye %} in the history window.

> ### {% icon hands_on %} Hands-on: Cutting necessary metadata columns
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1,c4,c6,c8,c10,c14,c20,c24,c26,c30,c32,c34`
>    - {% icon param-file %} *"From"*: `design_tsv` (output of **EBI SCXA Data Retrieval** {% icon tool %})
>
{: .hands_on}

You can inspect the dataset {% icon galaxy-eye %} to see that it's full of annoying "" everywhere, and overly long descriptions of each columns.

![Columns in Galaxy history window contain "Assay" and "" around every word or ID](../../images/bulk-music/annoying_ebimetadata.png History window "Annoying metadata")

Now, there might be a better way to do this in Galaxy (or you might consider downloading the file locally and changing it in a spreadsheet application or something), but this is what will work to reformat all that annoying text.

> ### {% icon hands_on %} Hands-on: Reformatting the metadata
>
> 1. {% tool [Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.2) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `out_file1` (output of **Cut** {% icon tool %})
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Sample Characteristic\[individual\]`
>            - *"Replacement"*: `Individual`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Sample Characteristic\[sex\]`
>            - *"Replacement"*: `Sex`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Sample Characteristic\[age\]`
>            - *"Replacement"*: `Age`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Sample Characteristic\[body mass index\]`
>            - *"Replacement"*: `BMI`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: ` kilogram per square meter`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `HbA1c `
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Sample Characteristic\[clinical information\]`
>            - *"Replacement"*: `HbA1c`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `%`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Sample Characteristic\[disease\]`
>            - *"Replacement"*: `Disease`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Sample Characteristic\[single cell quality\]`
>            - *"Replacement"*: `Single cell quality`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Sample Characteristic\[submitted single cell quality\]`
>            - *"Replacement"*: `Submitted single cell quality`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Factor Value\[inferred cell type - ontology labels\]`
>            - *"Replacement"*: `Inferred cell type - ontology label`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Factor Value\[inferred cell type - authors labels\]`
>            - *"Replacement"*: `Inferred cell type - author labels`
>
> 2. Change the datatype to tabular.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular" %}
>
{: .hands_on}

Great, this file is now ready to go! But, it contains all those extra cells that didn't pass filtration with the EBI pipeline and therefore won't exist in the matrix. We need to remove them for future steps to work. We can use our barcodes list to remove the extra cells.

## Prepare the barcodes file

> ### {% icon hands_on %} Hands-on: Adding a header
>
> 1. {% tool [Add line to file](toolshed.g2.bx.psu.edu/repos/bgruening/add_line_to_file/add_line_to_file/0.1.0) %} with the following parameters:
>    - *"text to add"*: `Cell`
>    - {% icon param-file %} *"input file"*: `barcode_tsv` (output of **EBI SCXA Data Retrieval** {% icon tool %})
>
>    > ### {% icon comment %} Comment
>    >
>    > This is an annoying step we have to do to get the right format, otherwise future steps won't work.
>    {: .comment}
>
{: .hands_on}

## Use the barcodes list to filter out cells in the experimental design file

> ### {% icon hands_on %} Hands-on: Joining datasets
>
> 1. {% tool [Join two Datasets](join1) %} with the following parameters:
>    - {% icon param-file %} *"Join"*: `outfile` (output of **Add line to file** {% icon tool %})
>    - *"using column"*: `c1`
>    - {% icon param-file %} *"with"*: `out_file1` (output of **Regex Find And Replace** {% icon tool %})
>    - *"and column"*: `c1`
>    - *"Fill empty columns"*: `No`
>    - *"Keep the header lines"*: `Yes`
>
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

> ### {% icon hands_on %} Hands-on: Remove duplicate columns
>
> 1. {% tool [Advanced Cut](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/1.1.0) %} with the following parameters:
>    - {% icon param-file %} *"File to cut"*: `out_file1` (output of **Join two Datasets** {% icon tool %})
>    - *"Operation"*: `Discard`
>    - *"Cut by"*: `fields`
>        - *"List of Fields"*: `c1`
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


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.