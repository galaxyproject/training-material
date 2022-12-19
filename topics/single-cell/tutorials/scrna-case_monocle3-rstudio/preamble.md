# Introduction

This tutorial is the next one in the [Single-cell RNA-seq: Case Study]({% link topics/single-cell/index.md %}) series. This tutorial also focuses on trajectory analysis using [monocle3](https://cole-trapnell-lab.github.io/monocle3/), similarly to the [previous one]({% link topics/single-cell/tutorials/scrna-case_monocle3-trajectories/tutorial.md %}), but instead using Galaxy buttons, we will have a look what’s happening behind, in the code – we will be using R programming language. Sometimes you might encounter some limitations when working with Galaxy tools or you might want to make a wee modification that has to be done manually – it is useful then to be able to switch between RStudio and Galaxy smoothly. If you are not feeling confident enough with using R, [this tutorial]({% link topics/data-science/tutorials/r-basics/tutorial.md %}) is a good place to start. However, our tutorial is quite straightforward to follow and at the end you will feel like a programmer! On the other hand, if you are not confident with the biological or statistical theory behind trajectory analysis, check out the [slide deck]({% link topics/single-cell/tutorials/scrna-case_monocle3-trajectories/slides.html %}). With those resources (including the previous case study tutorials) you are well-equipped to go through this tutorial with ease. Let’s get started! 


## Get data
In the [previous tutorial]({% link topics/single-cell/tutorials/scrna-case_monocle3-trajectories/tutorial.md %}), we showed that Monocle3 works great with annotated data, but what if your data is not annotated yet? Is it still possible to use Monocle? The answer is yes, Monocle also allows annotating cells according to their type and it will be shown in this tutorial. First, we need to get appropriate data to work with. We will continue to work on the case study data from a mouse model of fetal growth restriction {% cite Bacon2018 %} (see [the study in Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-6945/results/tsne) and [the project submission](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6945/)). We will use the filtered AnnData object, before normalisation and annotation, generated in the [filtering tutorial]({% link topics/single-cell/tutorials/scrna-case_basic-pipeline/tutorial.md %}). You can simply go to the history of this tutorial, find step 20: Filtered Object and download it. For ease of use, that was already done for you and you can import the file from Zenodo below. 

{% include topics/single-cell/tutorials/scrna-case_monocle3-rstudio/get_data.md %} 

> <hands-on-title>Data upload </hands-on-title>
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



## Preparing the files
Monocle uses cell_data_set class to hold expression data; it requires three input files: expression_matrix, cell_metadata and gene_metadata. We will extract that information from our AnnData object. 


It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> <details-title> More details about the theory </details-title>
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output` (Input dataset)
>    - *"What to inspect?"*: `Key-indexed observations annotation (obs)`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}


> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output` (Input dataset)
>    - *"What to inspect?"*: `Key-indexed annotation of variables/features (var)`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}


> <hands-on-title>Task description</hands-on-title>
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output` (Input dataset)
>    - *"What to inspect?"*: `The full data matrix`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

We are now ready to open RStudio and start the coding part!  


## Launching RStudio
Thanks to available interactive tools, you can easily launch RStudio in Galaxy. 

{% snippet faqs/galaxy/interactive_tools_rstudio_launch.md %}


From now on, you can either view this tutorial in the RMarkdown, which will allow you to read the material and simultaneously execute the code cells, or write/paste the lines of code in RStudio. 

{% snippet faqs/galaxy/interactive_tools_rstudio_rmarkdown.md %}
