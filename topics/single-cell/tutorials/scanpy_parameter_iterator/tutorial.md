---
layout: tutorial_hands_on

title: Scanpy Parameter Iterator
zenodo_link: 'https://zenodo.org/record/8011681'
subtopic: tricks
priority: 3
questions:
- 
objectives:
- 
requirements:
-
    type: "internal"
    topic_name: single-cell
    tutorials:
        - scrna-case_alevin
        - scrna-case_alevin-combine-datasets
        - scrna-case_basic-pipeline

time_estimation: 1H
key_points:
- 

tags:
- single-cell
- tips&tricks
- parameter-iterator

contributions:
  authorship:
    - wee-snufkin
  funding:
    - eosc-life


---


# Introduction

The magic of bioinformatic analysis is that we use maths, statistics and complicated algorithms to deal with huge amounts of data to help us interpret the biology behind. However, it’s not always so straightforward – each tool has various parameters so eventually, we might end up with very different outcomes depending on the values we choose. With analysing scRNA-seq data, it’s almost like you need to know about 75% of your data and make sure your analysis shows that, for you to then identify the 25% new information. 
Since there is a vast number of values that we can specify in the tools, how can we know if the values we choose are the most optimal ones or at least good enough? Well, we can use different values and then compare the outputs to see which is consistent with our understanding of the underlying biology. 
And here the Parameter Iterator comes in – it allows to run the analysis using different variables quickly and easily. Now you don’t have to execute your workflow every time you change one value to compare the output. This tutorial will show you how to use Parameter Iterator to generate multiple outputs with different parameter values at one go. 

{% snippet faqs/galaxy/tutorial_mode.md %}


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get Data
The data used in this tutorial is from a mouse dataset of fetal growth restriction ({% cite Bacon2018 %}). You can download the dataset below or import the history with the starting data.

> <comment-title></comment-title> 
> If you've been working through the Single-cell RNA-seq: Case Study then you can use your dataset from the [Filter, Plot and Explore Single-cell RNA-seq Data]({% link topics/single-cell/tutorials/scrna-case_basic-pipeline/tutorial.md %}) tutorial here. We’ll be working on the output of **Scanpy RunPCA** {% icon tool %}, so you can run the analysis using this dataset. If you haven’t completed that tutorial but you’re interested in how we get to this point, feel free to have a look.
{: .comment}

Here are several ways of getting our toy dataset – choose whichever you like! 

> <hands-on-title>Option 1: Data upload - Import history</hands-on-title>
>
> 1. Import history from: [example input history](https://usegalaxy.eu/u/j.jakiela/h/scanpy-parameter-iterator)
>
>
>    {% snippet faqs/galaxy/histories_import.md %}
>
> 2. **Rename** {% icon galaxy-pencil %} the the history to your name of choice.
>
{: .hands_on}

> <hands-on-title>Option 2: Data upload - Add to history</hands-on-title>
>
> 1. Create a new history for this tutorial
>
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    {{ page.zenodo_link }}
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>    
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the dataset if you wish: 
>
>    {% snippet  faqs/galaxy/datasets_rename.md %}
>
> 4. Check that the datatype is `h5ad`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}



# Inputs
Scanpy ParameterIterator tool currently works only for the following parameters:
1.	Number of neighbours to derive kNN graph (for **Scanpy ComputeGraph** {% icon tool %})
2.	Perplexity (for **Scanpy RunTSNE** {% icon tool %})
3.	Resolution (for **Scanpy FindCluster** {% icon tool %})

There are two formats of the input values:
1.	List of all parameter values to be iterated
2.	Step increase values to be iterated

> <comment-title></comment-title>
> Iterating the parameters within one tool will give you a list with X datasets: each dataset is the output with the given parameter value. However, if you want to use Parameter Iterator again within another tool, specifying Y parameter values, you **will not** get X x Y datasets as you might expect. Therefore you have to choose **just one** output file to be passed on to the next tool which will use Parameter Iterator again. Alternatively, you can use Parameter Iterator once and run the rest of the tools on dataset collection with just one parameter value. 
{: .comment}

# Number of neighbours to derive kNN graph (for **Scanpy ComputeGraph** {% icon tool %})
Our dataset is right after PCA. Therefore we will now use Scanpy ComputeGraph to derive kNN graph. We can use Parameter Iterator to check how different values of the number of neighbors affect the final outcome. It is important that n-neighbours is an integer. 

Warning
Using ‘Step increase values to be iterated’ as the format of the input values automatically generates float values instead of integers. Therefore in this case you have to use ‘List of all parameter values to be iterated’ with your chosen values.

k-nearest neighbor (kNN) graph will be needed for plotting a UMAP. From [UMAP developers](https://github.com/lmcinnes/umap): “Larger neighbor values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50, with a choice of 10 to 15 being a sensible default”. Therefore, let’s pick some values bigger and smaller than 15 to check how it changes the final UMAP.

HANDS ON

Now you have two options: either pick one of the generated output files and proceed to the next tool with another parameter iteration or continue with the current collection of datasets. We choose the second option as only then you will be able to see the effect of using different n-neighbours values. However, the disadvantage of this attitude is that you have to come up with one value for the subsequent parameters in the workflow to see the changes in the final plots.

The part of the workflow we are using in this tutorial is as follows:
Scanpy 
For the detailed explanation of the tools used below, check out [this tutorial]({% link topics/single-cell/tutorials/scrna-case_basic-pipeline/tutorial.md %}).

HANDS ON

If you compare the UMAP graphs, you can see the differences that were caused by changing the value of n-neighbors. Relying on your biological knowledge, you can now choose which parameter value works best and use it for further analysis. 
We will go forward with n-neighbor value equal to 15. 
HANDS ON 
Unhide n-neighbours equal to 15


# Title of the section usually corresponding to a big step in the analysis

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


## Sub-step with **Scanpy ComputeGraph**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy ComputeGraph](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_compute_graph/scanpy_compute_graph/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy RunPCA** {% icon tool %})
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

## Sub-step with **Scanpy RunTSNE**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy RunTSNE](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_run_tsne/scanpy_run_tsne/1.8.1+galaxy9) %} with the following parameters:
>    - *"Use the indicated representation"*: `X_pca`
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

## Sub-step with **Scanpy FindCluster**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy FindCluster](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_find_cluster/scanpy_find_cluster/1.8.1+galaxy9) %} with the following parameters:
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

## Sub-step with **Scanpy RunTSNE**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy RunTSNE](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_run_tsne/scanpy_run_tsne/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy ComputeGraph** {% icon tool %})
>    - *"Use the indicated representation"*: `X_pca`
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

## Sub-step with **Scanpy RunUMAP**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy RunUMAP](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_run_umap/scanpy_run_umap/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy RunTSNE** {% icon tool %})
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

## Sub-step with **Scanpy FindMarkers**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy FindMarkers](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_find_markers/scanpy_find_markers/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy FindCluster** {% icon tool %})
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

## Sub-step with **Scanpy RunUMAP**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy RunUMAP](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_run_umap/scanpy_run_umap/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy RunTSNE** {% icon tool %})
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

## Sub-step with **Scanpy FindCluster**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy FindCluster](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_find_cluster/scanpy_find_cluster/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy RunUMAP** {% icon tool %})
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

## Sub-step with **Scanpy PlotEmbed**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy PlotEmbed](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_plot_embed/scanpy_plot_embed/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy FindMarkers** {% icon tool %})
>    - *"name of the embedding to plot"*: `pca`
>    - *"color by attributes, comma separated texts"*: `louvain`
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

## Sub-step with **Scanpy PlotEmbed**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy PlotEmbed](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_plot_embed/scanpy_plot_embed/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy FindMarkers** {% icon tool %})
>    - *"color by attributes, comma separated texts"*: `louvain`
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

## Sub-step with **Scanpy PlotEmbed**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy PlotEmbed](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_plot_embed/scanpy_plot_embed/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy FindMarkers** {% icon tool %})
>    - *"name of the embedding to plot"*: `tsne`
>    - *"color by attributes, comma separated texts"*: `louvain`
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

## Sub-step with **Scanpy FindCluster**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy FindCluster](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_find_cluster/scanpy_find_cluster/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy RunUMAP** {% icon tool %})
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

## Sub-step with **Scanpy FindMarkers**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy FindMarkers](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_find_markers/scanpy_find_markers/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy FindCluster** {% icon tool %})
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

## Sub-step with **Scanpy FindMarkers**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy FindMarkers](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_find_markers/scanpy_find_markers/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy FindCluster** {% icon tool %})
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

## Sub-step with **Scanpy PlotEmbed**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy PlotEmbed](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_plot_embed/scanpy_plot_embed/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy FindMarkers** {% icon tool %})
>    - *"name of the embedding to plot"*: `tsne`
>    - *"color by attributes, comma separated texts"*: `louvain`
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

## Sub-step with **Scanpy PlotEmbed**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy PlotEmbed](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_plot_embed/scanpy_plot_embed/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy FindMarkers** {% icon tool %})
>    - *"color by attributes, comma separated texts"*: `louvain`
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


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
