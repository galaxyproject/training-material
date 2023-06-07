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
> 3. Rename the dataset if you wish: `Scanpy RunPCA: AnnData object`
>
>    {% snippet  faqs/galaxy/datasets_rename.md %}
>
> 4. Check that the datatype is `h5ad`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}


# Workflow

This tutorial is an extension of the full analysis shown in [Filter, Plot and Explore Single-cell RNA-seq Data]({% link topics/single-cell/tutorials/scrna-case_basic-pipeline/tutorial.md %}) tutorial in the Single-cell RNA-seq: Case Study series. So if you've been working through it, you can use your dataset from that tutorial here. If you haven’t completed it but you’re interested in how we get to this point, feel free to have a look at the mentioned tutorial.

Our starting data will be the output of **Scanpy RunPCA** {% icon tool %}. It is part of the full analysis tutorial, but we will only focus on a smaller and shortened bit of the [workflow]() to show the application of the Parameter Iterator. Our workflow consists of the following steps:

![Workflow that we are going to use in this tutorial: Scanpy RunPCA, Scanpy ComputeGraph, Scanpy RunTSNE, Scanpy RunUMAP, Scanpy FindCluster, Scanpy PlotEmbed](../../images/scrna-casestudy_parameter-iterator/workflow.png "Shortened workflow that we are going to use in this tutorial.")

For the detailed explanation of the tools presented above, check out [this tutorial]({% link topics/single-cell/tutorials/scrna-case_basic-pipeline/tutorial.md %}).


# Inputs
Scanpy ParameterIterator tool currently works only for the following parameters:
1.	Number of neighbours to derive kNN graph (for **Scanpy ComputeGraph** {% icon tool %})
2.	Perplexity (for **Scanpy RunTSNE** {% icon tool %})
3.	Resolution (for **Scanpy FindCluster** {% icon tool %})

There are two formats of the input values:
1.	List of all parameter values to be iterated
2.	Step increase values to be iterated


# Number of neighbours to derive kNN graph (for **Scanpy ComputeGraph** {% icon tool %})

Our dataset is right after PCA. Therefore we will now use **Scanpy ComputeGraph** {% icon tool %} to derive kNN graph. We can use Parameter Iterator to check how different values of the number of neighbours affect the final outcome. It is important that **n-neighbours is an integer**. 

> <warning-title>Float vs integer</warning-title>
> Using ‘Step increase values to be iterated’ as the format of the input values automatically generates float values instead of integers. Therefore in this case you have to use ‘List of all parameter values to be iterated’ with your chosen values.
{: .warning}


k-nearest neighbor (kNN) graph will be needed for plotting a UMAP. From [UMAP developers](https://github.com/lmcinnes/umap): “Larger neighbor values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50, with a choice of 10 to 15 being a sensible default”. Therefore, let’s pick some values bigger and smaller than 15 to check how it changes the final UMAP. This is where the Parameter Iterator comes in!

> <hands-on-title> Set your values in Parameter Iterator </hands-on-title>
>
> 1. {% tool [Scanpy ParameterIterator](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_parameter_iterator/scanpy_parameter_iterator/0.0.1+galaxy9) %} with the following parameters:
>    - *"Parameter type"*: `n-neighbours`
>    - *"Choose the format of the input values"*: `List of all parameter values to be iterated
>    - *"User input values"*: `5,10,15,20,25,30,35,40`
> 
> 2. **Rename** {% icon galaxy-pencil %} the resulting list of datasets: `Parameter iterated - n-neighbours` (you have to first click on the collection so that you see the datasets, and then raname it)
> 
> 3. **Tag** {% icon galaxy-tags %} each dataset with its corresponding value: 
>    - *n-neighbours_10*: `#n-neighbours_10` etc.
>    If you want to refresh your memory on how to add tags to datasets, have a look here: 
>    {% snippet faqs/galaxy/datasets_add_tag.md %} 
>    
{: .hands_on}

The output of the Parameter Iterator is the list of datasets. We will be working on dataset collections quite a lot, so if you want to gain more understanding of collection operations, visit the [corresponding tutorial](% link topics/galaxy-interface/tutorials/collections/tutorial.md %)

> <hands-on-title> Derive kNN graph with iterated parameter </hands-on-title>

> 1. {% tool [Scanpy ComputeGraph](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_compute_graph/scanpy_compute_graph/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `Scanpy RunPCA: AnnData object`
>    - *"Use programme defaults"*: {% icon history-share %} `No`
>    - *"File with n_neighbours, use with parameter iterator. Overrides the n_neighbors setting"*:
>       - Click on {% icon param-collection %} (*Dataset collection*)
>       - Choose `Parameter iterated - n-neighbours`
>    - *"Use the indicated representation"*: `X_pca`
>    - *"Number of PCs to use"*: `20`
>
{: .hands_on}

You should now see the output `Scanpy ComputeGraph on collection 2: Graph object AnnData` - it’s also a collection. If you click on that, you will see Anndata files, differing only by the n-neighbour value. You can access those files separately if you go back to your history and click on {% icon eye-slash %} *Show hidden*. You can bring each individual dataset to the visible and active datasets by clicking {% icon eye-slash %} *Unhide*.

Now you have two options: either pick one of the generated output files and proceed to the next tool with another parameter iteration or continue with the current collection of datasets. We choose the second option as only then you will be able to see the effect of using different n-neighbours values. However, the disadvantage of this option is that you have to come up with one value for the subsequent parameters in the workflow to see the changes in the final plots.

> <comment-title>Why only one Parameter Iteration per workflow?</comment-title>
> Iterating the parameters within one tool will give you a list with X datasets: each dataset is the output with the given parameter value. However, if you want to use Parameter Iterator again within another tool, specifying Y parameter values, you **will not** get X x Y datasets as you might expect. Therefore you have to choose **just one** output file to be passed on to the next tool which will use Parameter Iterator again. Alternatively, you can use Parameter Iterator once and run the rest of the tools on dataset collection with just one parameter value. 
{: .comment}

Where are we now in our workflow?
![Image showing the step we are at: after Scanpy RunPCA, already run Scanpy ComputeGraph, and before Scanpy RunTSNE, Scanpy RunUMAP, Scanpy FindCluster, Scanpy PlotEmbed](../../images/scrna-casestudy_parameter-iterator/workflow_n-neighbours.png "We used Parameter Iterator for the n-neighbours to derive kNN graph, now we’ll complete our small workflow to see the differences at the end.")

> <hands-on-title> Complete the workflow </hands-on-title>
>
> 1. {% tool [Scanpy RunTSNE](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_run_tsne/scanpy_run_tsne/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-collection %} *"Input object in AnnData/Loom format"* (make sure you choose *Dataset collection*): `Scanpy ComputeGraph on collection X: Graph object AnnData` 
>    - *"Use the indicated representation"*: `X_pca`
>    - *"Use programme defaults"*: {% icon history-share %} `No`
>    - *"The perplexity is related to the number of nearest neighbours, select a value between 5 and 50"*: `30`
>
> 2. {% tool [Scanpy RunUMAP](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_run_umap/scanpy_run_umap/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-collection %} *"Input object in AnnData/Loom format"* (make sure you choose *Dataset collection*): `Scanpy RunTSNE on collection X: tSNE object AnnData` 
>    - *"Use programme defaults"*: {% icon history-share %} `Yes`>
>
> 3. {% tool [Scanpy FindCluster](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_find_cluster/scanpy_find_cluster/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-collection %} *"Input object in AnnData/Loom format"* (make sure you choose *Dataset collection*): `Scanpy RunUMAP on collection X: UMAP object AnnData` 
>    - *"Use programme defaults"*: {% icon history-share %} `No`
>    - *"Resolution, high value for more and smaller clusters"*: `0.6`
{: .hands_on}

In fact, the differences will be only seen in the UMAP embedding, so we’ll plot only them. However, when you run your own analysis, you might want to check if there aren’t any changes in other embeddings as well. 

> <hands-on-title> Plot UMAP embedding </hands-on-title>
>
> 1. {% tool [Scanpy PlotEmbed](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_plot_embed/scanpy_plot_embed/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-collection %} *"Input object in AnnData/Loom format"* (make sure you choose *Dataset collection*): `Scanpy FindCluster on collection X: Clusters AnnData`
>    - *"name of the embedding to plot"*: `umap`
>    - *"color by attributes, comma separated texts"*: `louvain`
>    - *"Use raw attributes if present"*: `No`
>
{: .hands_on}

If you click on the resulting collection you will see several plots. Click on {% icon galaxy-eye %} to see how they differ. Galaxy's {% icon galaxy-scratchbook %} Window Manager, which you can enable (and disable again) from the menu bar can be very helpful for comparing multiple datasets.

![Eight graphs showing the differences between UMAP embeddings caused by different values of n-neighbours.](../../images/scrna-casestudy_parameter-iterator/n_neighbours.png "Comparison of UMAP embedding with different values of n-neighbours, perplexity set to 30 and resolution to 0.6.")

If you compare the UMAP graphs, you can see the differences that were caused by changing the value of n-neighbours. Relying on your biological knowledge, you can now choose which parameter value works best and use it for further analysis. 
We will go forward with n-neighbour value equal to 15. 

HANDS ON 
Unhide n-neighbours equal to 15

## Sub-step with **Scanpy PlotEmbed**


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
