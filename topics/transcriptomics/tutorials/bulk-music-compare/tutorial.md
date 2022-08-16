---
layout: tutorial_hands_on
subtopic: single-cell
priority: 11
title: Bulk deconvolution with MuSiC across multiple variables
zenodo_link: https://zenodo.org/record/5719228
tags:
  - single-cell
  - mouse
  - human
  - deconvolution
  - bulk
questions:
- How do the cell type distributions vary in bulk RNA samples across my variable of interest?
- For example, are beta cell proportions different in the pancreas data from diabetes and healthy patients?
objectives:
- Apply the MuSiC deconvolution to samples and compare the cell type distributions
- Use statistical outputs to identify significant differences
  what you should focus on during the course
- Describe methods of input - i.e. using a single scRNA-seq datatype or multiple scRNA-seq files
time_estimation: 2H
key_points:
- Deconvolution can be used to compare cell type distributions from bulk RNA-seq datasets
contributors:
- nomadscientist
- mtekman
requirements:
-
    topic_name: transcriptomics
    tutorials:
        - bulk-music
---


# Introduction
{:.no_toc}

<!-- This is a comment. -->
#TODO
- TRY A WORKFLOW WITH 3 DIFFERENT VARIABLES (does the stats handle this?)
- Fix main tutorial dataset descriptions
- Remake workflow with separated scRNA-seq references
- Future: include link to 'how we made the datasets'

The goal of this tutorial is to apply bulk RNA deconvolution techniques to a problem with multiple variables - in this case, a model of diabetes is compared with its healthy counterparts. All you need is well-annotated, high quality reference scRNA-seq dataset (or multiple!) and your bulk RNA-samples of choice. For more information on how MuSiC works, you can check out their github site [MuSiC](https://xuranw.github.io/MuSiC/articles/MuSiC.html) or published article {% cite wang2019bulk %}


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Data

We will use the same data from the Deconvolution with MuSiC tutorial. As a reminder, this means we will extract cell proportions from a bulk data of human pancreas data from {%cite fadista2014global %} concerning 56 638 genes across 89 samples, using a single cell human pancreas dataset from {%cite segerstolpe2016single %} containing 25 453 genes across 2209 cells, clustered into 14 cell types, from 6 healthy subjects and 4 with Type-II diabetes (T2D). One of our first tasks will be separating these datasets into healthy and T2D so we can analyse them separately.

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial *"Deconvolution: Cell Type inference of Human Pancreas Data"*
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    * Human pancreas bulk RNA datasets (tag: `#bulk`)
>
>      ```
>      https://zenodo.org/record/5719228/files/GSE50244bulkeset.expression.tabular
>      https://zenodo.org/record/5719228/files/GSE50244bulkeset.phenotype.tabular
>      ```
>    * Human pancreas single-cell RNA datasets (tag: `#scrna`)
>      ```
>      https://zenodo.org/record/5719228/files/EMTABesethealthy.expression.tabular
>      https://zenodo.org/record/5719228/files/EMTABesethealthy.phenotype.tabular
>      ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Rename the datasets
>
> 4. Check the datatype
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular" %}
>
> 5. Add to each `expression` file a tag corresponding to `#bulk` and `#scrna`
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

{% snippet faqs/galaxy/tutorial_mode.md %}

### Exploring the Datasets


The data consists of:
  * The bulk human pancreas dataset is 89 samples across 56 638 genes

  * The single cell human pancreas datasets is 2209 cells across 24 453 genes.

If you examine {% icon galaxy-eye %} the bulk dataset phenotype, you will remember that it appears as follows:

![peek_tabular_bulk_pheno.png](../../images/bulk-music/peek_tabular_bulk_pheno.png "Peeking at the tabular bulk RNA-seq phenotype dataset")

The expression meanwhile looks like this:
![bulk_tab](../../images/bulk-music/peek_tabular_bulk_expr.png "Peeking at the tabular bulk RNA-seq expression dataset")

Importantly, the first six columns are healthy controls, while the final four represent a disease model. We have our two variables to compare! Let's separate out these datasets.

# Title of the section usually corresponding to a big step in the analysis

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


## Sub-step with **Cut**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1,c2,c3,c4,c5,c6,c7`
>    - {% icon param-file %} *"From"*: `output` (Input dataset)
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

## Sub-step with **Split file**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Split file](toolshed.g2.bx.psu.edu/repos/bgruening/split_file_on_column/tp_split_on_column/0.4) %} with the following parameters:
>    - {% icon param-file %} *"File to select"*: `output` (Input dataset)
>    - *"on column"*: `c4`
>    - *"Include the header in all splitted files?"*: `Yes`
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

## Sub-step with **Construct Expression Set Object**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Construct Expression Set Object](toolshed.g2.bx.psu.edu/repos/bgruening/music_construct_eset/music_construct_eset/0.1.1+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Assay Data"*: `out_file1` (output of **Cut** {% icon tool %})
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

## Sub-step with **MuSiC Compare**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [MuSiC Compare](toolshed.g2.bx.psu.edu/repos/bgruening/music_compare/music_compare/0.1.1+galaxy4) %} with the following parameters:
>    - In *"New scRNA Group"*:
>        - {% icon param-repeat %} *"Insert New scRNA Group"*
>            - *"Name of scRNA Dataset"*: `scrna-total`
>            - In *"Advanced scRNA Parameters"*:
>                - *"Comma list of cell types to use from scRNA dataset"*: `Endo,Podo,PT,LOH,DCT,CD-PC,CD-IC,Fib,Macro,Neutro,B lymph,T lymph,NK`
>            - In *"Bulk Datasets in scRNA Group"*:
>                - {% icon param-repeat %} *"Insert Bulk Datasets in scRNA Group"*
>                    - *"Name of Bulk Dataset"*: `bulk`
>                    - {% icon param-file %} *"Bulk RNA Dataset"*: `out_rds` (output of **Construct Expression Set Object** {% icon tool %})
>                    - *"Factor Name"*: `Control`
>                - {% icon param-repeat %} *"Insert Bulk Datasets in scRNA Group"*
>                    - *"Name of Bulk Dataset"*: `APOL`
>                    - *"Factor Name"*: `Control`
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

## Sub-step with **MuSiC Compare**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [MuSiC Compare](toolshed.g2.bx.psu.edu/repos/bgruening/music_compare/music_compare/0.1.1+galaxy4) %} with the following parameters:
>    - In *"New scRNA Group"*:
>        - {% icon param-repeat %} *"Insert New scRNA Group"*
>            - *"Name of scRNA Dataset"*: `scrna-total`
>            - In *"Advanced scRNA Parameters"*:
>                - *"Comma list of cell types to use from scRNA dataset"*: `Endo,Podo,PT,LOH,DCT,CD-PC,CD-IC,Fib,Macro,Neutro,B lymph,T lymph,NK`
>            - In *"Bulk Datasets in scRNA Group"*:
>                - {% icon param-repeat %} *"Insert Bulk Datasets in scRNA Group"*
>                    - *"Name of Bulk Dataset"*: `bulk-control`
>                    - {% icon param-file %} *"Bulk RNA Dataset"*: `out_rds` (output of **Construct Expression Set Object** {% icon tool %})
>                    - *"Factor Name"*: `Control`
>                    - In *"Advanced Bulk Parameters"*:
>                        - *"Phenotype factors"*: `Control`
>                - {% icon param-repeat %} *"Insert Bulk Datasets in scRNA Group"*
>                    - *"Name of Bulk Dataset"*: `APOL-bulk`
>                    - *"Factor Name"*: `Control`
>                    - In *"Advanced Bulk Parameters"*:
>                        - *"Phenotype factors"*: `Control`
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
