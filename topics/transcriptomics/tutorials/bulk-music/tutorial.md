---
layout: tutorial_hands_on

title: Bulk RNA Deconvolution with MuSiC
zenodo_link: ''
questions:
- How do we infer cell type proportions from bulk RNA-seq data?
- How are these cell types grouped together?
objectives:
- Construct Bulk and scRNA Expression Set Objects
- Inspect these objects for various properties
- Measure the abundance of certain cell type cluster markers compared to another
time_estimation: '2H'
key_points:
- Deconvolution tools show individual cell type proportions in bulk RNA-seq data
- Bulk RNA-seq can be complimented by scRNA-seq data
contributors:
- mtekman
- nomadscientist

---


# Introduction
{:.no_toc}

<!-- @Wendi - I'm using info from here: https://xuranw.github.io/MuSiC/articles/MuSiC.html -->

Bulk RNA-seq expression data obtained from RNA-sequencing contains a mixture of the expression of several types of cells. We wish to deconvolve this data to obtain more precise estimates of the proportions of this type of data.

By combining bulk data with multi-subject single cell expression data obtained from single-cell RNA-sequencing, we can use this as a reference for estimating the cell type proportions in the bulk data.

In this tutorial we will be using the [MuSiC](https://xuranw.github.io/MuSiC/articles/MuSiC.html) tool suite {% cite wang2019bulk %}.


**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Bulk RNA-seq Cell Type Deconvolution

Explanation of some terms


## Background

What is deconvoltion and why is it needed.


## MusiC

Music utilizes blah blah blah

<!-- Nice image here (https://xuranw.github.io/MuSiC/articles/MuSiC.html) of how it works -->

## Expression Set Objects

Expressions Sets are Rdata sets invented by Bioconductor for reasons X Y Z


# Workflow Overview 

Here is a nice image of the workflow that will be used in this tutorial

![workflow1](../../images/bulk-music/workflow1.png "Workflow of Steps")


## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}



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


## Sub-step with **Construct Expression Set Object**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Construct Expression Set Object](music_construct_eset) %} with the following parameters:
>    - {% icon param-file %} *"Assay Data"*: `output` (Input dataset)
>    - {% icon param-file %} *"Phenotype Data"*: `output` (Input dataset)
>    - In *"Meta Data"*:
>        - {% icon param-repeat %} *"Insert Meta Data"*
>            - *"Label"*: `{'id': 2, 'output_name': 'output'}`
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
> 1. {% tool [Construct Expression Set Object](music_construct_eset) %} with the following parameters:
>    - {% icon param-file %} *"Assay Data"*: `output` (Input dataset)
>    - *"Annotation"*: `{'id': 4, 'output_name': 'output'}`
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

## Sub-step with **Inspect Expression Set Object**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Inspect Expression Set Object](music_inspect_eset) %} with the following parameters:
>    - {% icon param-file %} *"ESet Dataset"*: `out_rds` (output of **Construct Expression Set Object** {% icon tool %})
>    - *"Inspect"*: `Feature Data Table`
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

## Sub-step with **Inspect Expression Set Object**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Inspect Expression Set Object](music_inspect_eset) %} with the following parameters:
>    - {% icon param-file %} *"ESet Dataset"*: `out_rds` (output of **Construct Expression Set Object** {% icon tool %})
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

## Sub-step with **Inspect Expression Set Object**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Inspect Expression Set Object](music_inspect_eset) %} with the following parameters:
>    - {% icon param-file %} *"ESet Dataset"*: `out_rds` (output of **Construct Expression Set Object** {% icon tool %})
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

## Sub-step with **MuSiC**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [MuSiC](music_deconvolution) %} with the following parameters:
>    - {% icon param-file %} *"scRNA Dataset"*: `out_rds` (output of **Construct Expression Set Object** {% icon tool %})
>    - {% icon param-file %} *"Bulk RNA Dataset"*: `out_rds` (output of **Construct Expression Set Object** {% icon tool %})
>    - *"Purpose"*: `Estimate Proportions`
>        - *"Comma list of cell types to use from scRNA dataset"*: `{'id': 5, 'output_name': 'output'}`
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

## Sub-step with **MuSiC**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [MuSiC](music_deconvolution) %} with the following parameters:
>    - {% icon param-file %} *"scRNA Dataset"*: `out_rds` (output of **Construct Expression Set Object** {% icon tool %})
>    - {% icon param-file %} *"Bulk RNA Dataset"*: `out_rds` (output of **Construct Expression Set Object** {% icon tool %})
>    - *"Purpose"*: `Compute Dendrogram`
>        - In *"Cluster Groups"*:
>            - {% icon param-repeat %} *"Insert Cluster Groups"*
>                - *"Cluster ID"*: `C1`
>            - {% icon param-repeat %} *"Insert Cluster Groups"*
>                - *"Cluster ID"*: `C2`
>                - {% icon param-file %} *"List of Gene Markers"*: `output` (Input dataset)
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