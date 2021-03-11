---
layout: tutorial_hands_on

title: Filtering and Plotting Single-cell RNA-seq Data
zenodo_link: ''
questions:
- Which biological questions are addressed by the tutorial?
- Which bioinformatics techniques are important to know for this type of data?
objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentences describing what a learner should be able to do once they
  have completed the tutorial
- You can use Bloom's Taxonomy to write effective learning objectives
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- contributor1
- contributor2

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

There are many packages, which have collection of functions for running this workflow. Such as Seurat, Scanpy, Monocle3, Scater. They offer streamlined workflow and tutorials that can be easily followed.

Note - this tutorial is similar to another fantastic tutorial: [Clustering 3k PBMC](
https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/scrna-scanpy-pbmc3k/tutorial.html). That tutorial will go into much further depth on the analysis, in particular the visualisation and science behind identifying marker genes. Plus, their experimental data is clean and well annotated which illustrates the steps beautifully. Here, we work more as a case study with slightly messier (but real!) data. Furthermore, we have set up the steps in this tutorial to allow for easy access to associated tutorials for plotting trajectories, batch correction, decision-making, and plotting. We highly recommend you work through all the galaxy single cell tutorials to build confidence and expertise!

 from the  from . This tutorial will guide you to how to check the quality of these matrix files and analyse them. This includes: checking the quality of the generated matrix files, preprecessing the data for analysis, correcting batch effect and visualising them by clustering. This is most basic workflow that you will run in any single-cell data analysis.


Today, we will follow the scanpy workflow. Scanpy is python package based on object called anndata, which is data storage format made easy to store metadata and processed data from each stages.

General introduction about the topic and then an introduction of the
tutorial (the questions and the objectives). It is nice also to have a
scheme to sum up the pipeline used during the tutorial. The idea is to
give to trainees insight into the content of the tutorial and the (theoretical
and technical) key concepts they will learn.

You may want to cite some publications; this can be done by adding citations to the
bibliography file (`tutorial.bib` file next to your `tutorial.md` file). These citations
must be in bibtex format. If you have the DOI for the paper you wish to cite, you can
get the corresponding bibtex entry using [doi2bib.org](https://doi2bib.org).

With the example you will find in the `tutorial.bib` file, you can add a citation to
this article here in your tutorial like this:
{% raw %} `{% cite Batut2018 %}`{% endraw %}.
This will be rendered like this: {% cite Batut2018 %}, and links to a
[bibliography section](#bibliography) which will automatically be created at the end of the
tutorial.


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

# Title for your first section

Give some background about what the trainees will be doing in the section.
Remember that many people reading your materials will likely be novices,
so make sure to explain all the relevant concepts.

## Title for a subsection
Section and subsection titles will be displayed in the tutorial index on the left side of
the page, so try to make them informative and concise!

# Hands-on Sections
Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% include snippets/add_tag.md %}
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

#Filtering

You have generated an annotated AnnData object from your raw scRNA-seq fastq files. However, you have only completed a 'rough' filter of your dataset - there will still be a number of 'cells' that are actually just background from empty droplets, as well as genes that could be sequencing artifacts or appear so rarely, that they will be difficult to perform any usable statistics on. This background garbage of both cells and genes not only makes it harder to distinguish real biological information from the background, as well as makes it computationally heavy to analyse. These spurious reads take a lot of computational power for each step to analyse! So, first on our agenda is to filter this matrix to give us cleaner data to extract meaningful insight from, and to allow faster analysis.

## Calculating Mitochondrial Content **AnnData Operations**

In our annotated AnnData object, we flagged mitochondrial-associated genes with a 'true' and 'false' column. It's time to do some calculating!

> ### {% icon hands_on %} Hands-on: Counting mitochondrial genes
>
> 1. **AnnData Operations** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input object in hdf5 AnnData format"*: `Input dataset` (Input dataset)
>    - *"Gene symbols field in AnnData"*: `NA.`
>    - In *"Flag genes that start with these names"*:
>        - {% icon param-repeat %} *"Insert Flag genes that start with these names"*
>            - *"Starts with"*: `True`
>            - *"Var name"*: `mito`
>    - *"Copy observations (such as clusters)"*: `Yes`
>    - *"Copy embeddings (such as UMAP, tSNE)"*: `Yes`
>    - *"Copy uns"*: `Yes`
>
> 2. Rename {% icon galaxy-pencil %} output `Mito-counted AnnData`
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
> 1. What has this tool done? How can you figure that out?
> 2. While you are figuring that out, how many genes and cells are in your object?
>
>   > ### {% icon tip %} Hint
>   > You want to use the same tool you used in the previous tutorial to examine your AnnData, since it's not necessarily as simple as clicking on it!
>   >
>   >   > ### {% icon hands_on %} Hands-on: Inspecting AnnData Objects
>   >   >
>   >   > 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy0) %} with the following parameters:
>   >   >    - {% icon param-file %} *"Annotated data matrix"*: Mito-counted AnnData`
>   >   >    - *"What to inspect?"*: `General information about the object`
>   >   > 2. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy0) %} with the following parameters:
>   >   >    - {% icon param-file %} *"Annotated data matrix"*: Mito-counted AnnData`
>   >   >    - *"What to inspect?"*: `Key-indexed observations annotation (obs)`
>   >   > 3. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy0) %} with the following parameters:
>   >   >    - {% icon param-file %} *"Annotated data matrix"*: Mito-counted AnnData`
>   >   >    - *"What to inspect?"*: `Key-indexed annotation of variables/features (var)`
>   >   {: .hands_on}
>   {: .tip}
> > ### {% icon solution %} Solution
> >
> > 1. If you examine your AnnData object, you'll find the addition of a number of different quality control metrics for boths cells (`obs`) and genes (`var`). For instance, you can see a `n_cells` under `var`, which counts the number of cells that gene appears in. In the `obs`, you have both discrete and log-based metrics for `n_genes`, how many genes are counted in a cell, and `n_counts`, how many UMIs are counted per cell, so for instance you might count multiple GAPDHs so your `n_counts` will be higher than `n_cells`. But what about the mitochondria?? Within the cells information (`obs`), the `total_counts_mito`,  `log1p_total_counts_mito`, and `pct_counts_mito` has been calculated for each cell.
> > 2. You can see in the {% icon tool %} **General information about the object** output that the matrix is `25281 x 35734`. This is `obs x vars`, or rather, `cells x genes`, so there are `25281 cells` and `35734 genes` in the matrix.
> >
> {: .solution}
>
{: .question}

> ### {% icon tip %} Tip: Turn the 3 **Inspect AnnData** outputs above into a workflow for quick access!
>
> {% include snippets/create_new_workflow.md %}
>
{: .tip}

## Generate QC Plots

We want to filter our cells, but first we need to know what our data looks like. There are a number of 'decision' steps within scRNA-seq analysis, so we need to make our best informed decisions about where to set our thresholds. We're going to plot our data a few different ways. Different bioinformaticians might prefer to see the data in different ways, we are only generating some of the myriad of plots you can use, and ultimately you need to go with what makes the most sense to you. I strongly recommend turning the following QC plots into a workflow so you can re-run it easily as you filter your data, when you're working with your real samples!

{% include snippets/create_new_workflow.md %}

Add the following tools to that workflow, and then run {% icon workflow-run %} the workflow on {% icon param-file %} `Mito-counted AnnData`.

### Creating the plots

> ### {% icon comment %} I don't know how to make a workflow?
>
> You can simply run these tools instead, as normal, without making a workflow. Long-term, however, please see the **How to make a workflow subsection** [here](https://training.galaxyproject.org/training-material/topics/introduction/tutorials/galaxy-intro-101-everyone/tutorial.html) - so that you can make your Galaxy life a lot easier in the future!
{: .comment}

> ### {% icon hands_on %} Hands-on: Making QC plots
>
> 1. **Plot** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Mito-counted AnnData`
>    - *"Method used for plotting"*: `Generic: Violin plot, using 'pl.violin'`
>        - *"Keys for accessing variables"*: `Subset of variables in 'adata.var_names' or fields of '.obs'`
>            - *"Keys for accessing variables"*: `log1p_total_counts,log1p_n_genes_by_counts,pct_counts_mito`
>        - *"The key of the observation grouping to consider"*: `genotype`
>        - In *"Violin plot attributes"*:
>            - *"Add a stripplot on top of the violin plot"*: `Yes`
>                - *"Add a jitter to the stripplot"*: `Yes`
>            - *"Display keys in multiple panels"*: `No`
> 2. Rename {% icon galaxy-pencil %} output `Violin - genotype - log`
>
> 3. **Plot** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Mito-counted AnnData`
>    - *"Method used for plotting"*: `Generic: Violin plot, using 'pl.violin'`
>        - *"Keys for accessing variables"*: `Subset of variables in 'adata.var_names' or fields of '.obs'`
>            - *"Keys for accessing variables"*: `log1p_total_counts,log1p_n_genes_by_counts,pct_counts_mito`
>        - *"The key of the observation grouping to consider"*: `sex`
>        - In *"Violin plot attributes"*:
>            - *"Add a stripplot on top of the violin plot"*: `Yes`
>                - *"Add a jitter to the stripplot"*: `Yes`
>            - *"Display keys in multiple panels"*: `No`
>
> 4. Rename {% icon galaxy-pencil %} output `Violin - sex - log`
>
> 5. **Plot** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Mito-counted AnnData`
>    - *"Method used for plotting"*: `Generic: Violin plot, using 'pl.violin'`
>        - *"Keys for accessing variables"*: `Subset of variables in 'adata.var_names' or fields of '.obs'`
>            - *"Keys for accessing variables"*: `log1p_total_counts,log1p_n_genes_by_counts,pct_counts_mito`
>        - *"The key of the observation grouping to consider"*: `batch`
>        - In *"Violin plot attributes"*:
>            - *"Add a stripplot on top of the violin plot"*: `Yes`
>                - *"Add a jitter to the stripplot"*: `Yes`
>            - *"Display keys in multiple panels"*: `No`>
>
> 6. Rename {% icon galaxy-pencil %} output `Violin - batch - log`
>
> 7. **Plot** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Mito-counted AnnData`
>    - *"Method used for plotting"*: `Generic: Scatter plot along observations or variables axes, using 'pl.scatter'`
>        - *"Plotting tool that computed coordinates"*: `Using coordinates`
>            - *"x coordinate"*: `log1p_total_counts`
>            - *"y coordinate"*: `pct_counts_mito`
>            - *"Use the layers attribute?"*: `No`
>
> 6. Rename {% icon galaxy-pencil %} output `Scatter - mito x UMIs`
>
> 7. **Plot** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Mito-counted AnnData`
>    - *"Method used for plotting"*: `Generic: Scatter plot along observations or variables axes, using 'pl.scatter'`
>        - *"Plotting tool that computed coordinates"*: `Using coordinates`
>            - *"x coordinate"*: `log1p_n_genes_by_counts`
>            - *"y coordinate"*: `pct_counts_mito`
>            - *"Use the layers attribute?"*: `No`
>
> 8. Rename {% icon galaxy-pencil %} output `Scatter - mito x genes`
>
> 9. **Plot** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Mito-counted AnnData`
>    - *"Method used for plotting"*: `Generic: Scatter plot along observations or variables axes, using 'pl.scatter'`
>        - *"Plotting tool that computed coordinates"*: `Using coordinates`
>            - *"x coordinate"*: `log1p_total_counts`
>            - *"y coordinate"*: `log1p_n_genes_by_counts`
>            - *"Color by"*: `pct_counts_mito`
>            - *"Use the layers attribute?"*: `No`
>
> 10. Rename {% icon galaxy-pencil %} output `Scatter - genes x UMIs`
>
{: .hands_on}

### Analysing the plots

That's a lot of information! Let's attack this in sections and see what questions these plots can help us answer. The scratchbook {% icon galaxy-scratchbook %} may help here to look at the different plots at the same time!

{% include snippets/use_scratchbook.md %}

> ### {% icon question %} Question - Batch Variation
>
> Are their differences in sequencing depth across the samples?
> 1. Which plot(s) addresses this?
> 2. How do you interpret it?
>
> > ### {% icon solution %} Solution
> >
> > 1. The plot `violin - batch - log` will have what you're looking for!
> > ![Violin - batch - log](../../images/wab-violin-batch-log.png "Violin - batch - log (Raw)")
> >
> > 2. Keeping in mind that this is a log scale, so small differences can mean large differences, the violin plots probably look pretty similar. N707 and N703 might be a bit lower on genes and counts (or UMIs), but the differences aren't catastrophic. Also, the pct_counts_mito looks pretty similar across the batches, so this also looks good. Nothing here would cause me to eliminate a sample from my analysis, but if you see a sample looking completely different from the rest, you would need to question why that is and consider eliminating it from your experiment!
> >
> {: .solution}
>
{: .question}

> ### {% icon question %} Question - Biological Variables
>
> Are their differences in sequencing depth across sex? Genotype?
> 1. Which plot(s) addresses this?
> 2. How do you interpret the `sex` differences?
> 3. How do you interpret the `genotype` differences?
>
> > ### {% icon solution %} Solution
> >
> > 1. Similar to above, the plots `violin - genotype - log` and `violin - sex - log` will have what you're looking for!
> > ![Violin - genotype - log](../../images/wab-violin-genotype-log.png "Violin - genotype - log (Raw)")
> >
> > ![Violin - sex - log](../../images/wab-violin-sex-log.png "Violin - sex - log (Raw)")
> >
> > 2. There isn't a major difference, I would say - though you are welcome to disagree! It is clear there are far fewer female cells, which makes sense given that only one sample was female. Note - that was an unfortunate discovery made long after generating libraries. It's quite hard to identify the sex of a neonate in the lab! In practice, try hard to not let such a confounding factor into your data! You could consider re-running all the following analysis without that female sample, if you wish.  
> >
> > 3. In `genotype`, however, we can see there is a difference. The `knockout` samples have clearly a lower amount of both genes and counts.
> > - From an experimental point of view, we can consider, does this make sense?
> >    * Would we biologically expect that those cells would be smaller or having fewer transcripts? Possibly, in this case, given that these cells were generated by growth restricted neonatal mice, and in which case we don't need to worry.
> >    * On the other hand, it may be that those cells didn't survive dissociation as well as the healthy ones (in which case we'd expect higher mitochondrial-associated genes, which we don't see, so we can rule that out!).
> > - What do we do about it?
> >    * Ideally, we consider re-sequencing all the samples but with a higher concentration of the knockout samples in the library.
> > Any bioinformatician will tell you that the best way to get clean data is in the lab, not the computer! Sadly, absolute best practice isn't necessarily always a realistic option in the lab (as it can get quite time-consuming and expensive to generate 'perfect' data), so sometimes, we have to make the best of it. There are options to try and address such discrepancy in sequencing depth. Thus, we're going to take these samples forward and see if we can find biological insight despite the technical differences.

Now that we've assessed the differences in our samples, we will look at the libraries overall to identify appropriate thresholds for our analysis.

> ### {% icon question %} Question - Filter Thresholds
>
> What threshold should you set for `log1p_n_genes_by_counts`?
> 1. Which plot(s) addresses this?
> 2. What number would you pick?
>
> > ### {% icon solution %} Solution
> >
> > 1. Any plot with `log1p_n_genes_by_counts` would do here, actually! Some people prefer scatterplots.
> > ![Scatter-genesxmito](../../images/wab-scatter-genesxmito.png "Scatterplot - log1p_n_genes_by_counts x pct_counts_mito (Raw)")
> >
> > 2. In `Scatter - mito x genes` you can see how cells with `log1p_n_genes_by_counts` up to around, perhaps, `5.5` (around 250 genes) often have high `pct_counts_mito`. You can plot this as just `n_counts` and see this same trend at around 250 genes, but with this data the log format was clearer so that's how we're presenting it. You could also use the violin plots to come up with the threshold, and thus keeping batch into account. It's good to look at that as well, because you don't want to accidentally cut out an entire sample (i.e. N703 and N707). Some bioinformaticians would recommend filtering each sample individually, but this is difficult in larger scale and in this case (you're welcome to give it a go! You'd have to filter separately and then concatenate) won't make a notable difference in the final interpretation.
> >
> {: .solution}
>
> What threshold should you set for `log1p_total_counts`?
> 1. Which plot(s) addresses this?
> 2. What number would you pick?
>
> > ### {% icon solution %} Solution
> >
> > 1. As before, any plot with `log1p_n_total_counts` will do! Again, we'll use a scatterplot here, but you can use a violin plot if you wish!
> > ![Scatter-countsxmito](../../images/wab-scatter-countsxmito.png "Scatterplot - log1p_total_counts (Raw)")
> >
> > 2. We can see that we will need to set a higher threshold (which makes sense, as you'd expect more unique counts per cell rather than unique genes!). Again, perhaps being a bit aggressive, we might choose `5.8`, for instance (which amounts to around 300 counts/cell). In an ideal world, you'll see a clear population of real cells separated from a clear population of debris. Many samples, like this one, are under-sequenced, and such separation would likely be seen after deeper sequencing!
> >
> {: .solution}
>
> What threshold should you set for `pct_counts_mito`?
> 1. Which plot(s) addresses this?
> 2. What number would you pick?
>
> > ### {% icon solution %} Solution
> >
> > 1. Any plot with `pct_counts_mito` would do here, however the scatterplots are likely the easiest to interpret. We'll use the same as last time.
> > ![Scatter-countsxmito](../../images/wab-scatter-countsxmito.png "Scatterplot - log1p_total_counts (Raw)")
> >
> > 2. We can see a clear trend wherein cells that have around `5%` mito counts or higher also have far fewer total counts. These cells are low quality, will muddy our data, and are likely stressed or ruptured prior to encapsulation in a droplet. We will use quite a common cut-off of 5% here, however you must adapt all cut-offs to your data - metabolically active cells might have higher mitochondrial RNA in general, and you don't want to lose a cell population because of a cut-off.
> >
> {: .solution}
{: .question}

## Apply the Thresholds

It's now time to apply these thresholds to our data! First, a reminder of how many cells and genes are in your object: `25281 cells` and `35734 genes`. Let's see how that changes each time!

> ### {% icon details %} Working in a group? Decision-time!
> If you are working in a group, you can now divvy up a decision here with one *control* and the rest varied numbers so that you can compare results throughout the tutorials.
> - Control
>      - **log1p_n_genes_by_counts** > `5.5`
>      - **log1p_total_counts** > `5.8`
>      - **pct_counts_mito** > `5%`
> - Everyone else: Choose your own thresholds and compare results!
{: .details}

> >
> {: .solution}
>
{: .question}

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Scanpy FilterCells** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **AnnData Operations** {% icon tool %})
>    - In *"Parameters to select cells to keep"*:
>        - {% icon param-repeat %} *"Insert Parameters to select cells to keep"*
>            - *"Name of parameter to filter on"*: `pct_counts_mito`
>            - *"Max value"*: `4.5`
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

You could run this entire QC step using `n_counts` and `n_genes` instead! Either way works, so pick whichever one makes it easiest to interpret your data.

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

## Sub-step with **Inspect AnnData**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Inspect AnnData** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output_h5ad` (output of **AnnData Operations** {% icon tool %})
>    - *"What to inspect?"*: `General information about the object`
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



## Sub-step with **Inspect AnnData**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Inspect AnnData** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output_h5ad` (output of **AnnData Operations** {% icon tool %})
>    - *"What to inspect?"*: `Key-indexed observations annotation (obs)`
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

## Sub-step with **Plot**

> ### {% icon hands_on %} Hands-on: Task description


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

## Sub-step with **Inspect AnnData**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Inspect AnnData** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output_h5ad` (output of **AnnData Operations** {% icon tool %})
>    - *"What to inspect?"*: `Key-indexed annotation of variables/features (var)`
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

## Sub-step with **Scanpy FilterCells**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Scanpy FilterCells** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy FilterCells** {% icon tool %})
>    - In *"Parameters to select cells to keep"*:
>        - {% icon param-repeat %} *"Insert Parameters to select cells to keep"*
>            - *"Name of parameter to filter on"*: `log1p_n_genes_by_counts`
>            - *"Min value"*: `5.7`
>            - *"Max value"*: `20.0`
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

## Sub-step with **Plot**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Plot** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output_h5ad` (output of **Scanpy FilterCells** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Violin plot, using 'pl.violin'`
>        - *"Keys for accessing variables"*: `Subset of variables in 'adata.var_names' or fields of '.obs'`
>            - *"Keys for accessing variables"*: `log1p_total_counts,log1p_n_genes_by_counts,pct_counts_mito`
>        - *"The key of the observation grouping to consider"*: `genotype`
>        - In *"Violin plot attributes"*:
>            - *"Add a stripplot on top of the violin plot"*: `Yes`
>                - *"Add a jitter to the stripplot"*: `Yes`
>            - *"Display keys in multiple panels"*: `No`
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

## Sub-step with **Scanpy FilterCells**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Scanpy FilterCells** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy FilterCells** {% icon tool %})
>    - In *"Parameters to select cells to keep"*:
>        - {% icon param-repeat %} *"Insert Parameters to select cells to keep"*
>            - *"Name of parameter to filter on"*: `log1p_total_counts`
>            - *"Min value"*: `6.3`
>            - *"Max value"*: `20.0`
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

## Sub-step with **Plot**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Plot** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output_h5ad` (output of **Scanpy FilterCells** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Violin plot, using 'pl.violin'`
>        - *"Keys for accessing variables"*: `Subset of variables in 'adata.var_names' or fields of '.obs'`
>            - *"Keys for accessing variables"*: `log1p_total_counts,log1p_n_genes_by_counts,pct_counts_mito`
>        - *"The key of the observation grouping to consider"*: `genotype`
>        - In *"Violin plot attributes"*:
>            - *"Add a stripplot on top of the violin plot"*: `Yes`
>                - *"Add a jitter to the stripplot"*: `Yes`
>            - *"Display keys in multiple panels"*: `No`
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

## Sub-step with **Scanpy NormaliseData**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Scanpy NormaliseData** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy FilterCells** {% icon tool %})
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

## Sub-step with **Plot**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Plot** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output_h5ad` (output of **Scanpy NormaliseData** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Violin plot, using 'pl.violin'`
>        - *"Keys for accessing variables"*: `Subset of variables in 'adata.var_names' or fields of '.obs'`
>            - *"Keys for accessing variables"*: `log1p_total_counts,log1p_n_genes_by_counts,pct_counts_mito`
>        - *"The key of the observation grouping to consider"*: `genotype`
>        - In *"Violin plot attributes"*:
>            - *"Add a stripplot on top of the violin plot"*: `Yes`
>                - *"Add a jitter to the stripplot"*: `Yes`
>            - *"Display keys in multiple panels"*: `No`
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

## Sub-step with **Plot**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Plot** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output_h5ad` (output of **Scanpy NormaliseData** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Scatter plot along observations or variables axes, using 'pl.scatter'`
>        - *"Plotting tool that computed coordinates"*: `Using coordinates`
>            - *"x coordinate"*: `log1p_total_counts`
>            - *"y coordinate"*: `log1p_n_genes_by_counts`
>            - *"Color by"*: `pct_counts_mito`
>            - *"Use the layers attribute?"*: `No`
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

## Sub-step with **Plot**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Plot** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output_h5ad` (output of **Scanpy NormaliseData** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Scatter plot along observations or variables axes, using 'pl.scatter'`
>        - *"Plotting tool that computed coordinates"*: `Using coordinates`
>            - *"x coordinate"*: `log1p_n_genes_by_counts`
>            - *"y coordinate"*: `pct_counts_mito`
>            - *"Use the layers attribute?"*: `No`
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

## Sub-step with **Plot**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Plot** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output_h5ad` (output of **Scanpy NormaliseData** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Violin plot, using 'pl.violin'`
>        - *"Keys for accessing variables"*: `Subset of variables in 'adata.var_names' or fields of '.obs'`
>            - *"Keys for accessing variables"*: `log1p_total_counts,log1p_n_genes_by_counts,pct_counts_mito`
>        - *"The key of the observation grouping to consider"*: `sex`
>        - In *"Violin plot attributes"*:
>            - *"Add a stripplot on top of the violin plot"*: `Yes`
>                - *"Add a jitter to the stripplot"*: `Yes`
>            - *"Display keys in multiple panels"*: `No`
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

## Sub-step with **Plot**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Plot** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output_h5ad` (output of **Scanpy NormaliseData** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Violin plot, using 'pl.violin'`
>        - *"Keys for accessing variables"*: `Subset of variables in 'adata.var_names' or fields of '.obs'`
>            - *"Keys for accessing variables"*: `log1p_total_counts,log1p_n_genes_by_counts,pct_counts_mito`
>        - *"The key of the observation grouping to consider"*: `batch`
>        - In *"Violin plot attributes"*:
>            - *"Add a stripplot on top of the violin plot"*: `Yes`
>                - *"Add a jitter to the stripplot"*: `Yes`
>            - *"Display keys in multiple panels"*: `No`
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

## Sub-step with **Plot**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Plot** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output_h5ad` (output of **Scanpy NormaliseData** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Scatter plot along observations or variables axes, using 'pl.scatter'`
>        - *"Plotting tool that computed coordinates"*: `Using coordinates`
>            - *"x coordinate"*: `log1p_total_counts`
>            - *"y coordinate"*: `pct_counts_mito`
>            - *"Use the layers attribute?"*: `No`
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

## Sub-step with **Scanpy FindVariableGenes**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Scanpy FindVariableGenes** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy NormaliseData** {% icon tool %})
>    - *"Flavor of computing normalised dispersion"*: `Seurat`
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

## Sub-step with **Scanpy ScaleData**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Scanpy ScaleData** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy FindVariableGenes** {% icon tool %})
>    - *"Truncate to this value after scaling"*: `10.0`
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

## Sub-step with **Scanpy RunPCA**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Scanpy RunPCA** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy ScaleData** {% icon tool %})
>    - *"Perform incremental PCA by chunks"*: `Yes`
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

## Sub-step with **Plot**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Plot** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output_h5ad` (output of **Scanpy RunPCA** {% icon tool %})
>    - *"Method used for plotting"*: `PCA: Scatter plot in PCA coordinates, using 'pl.pca_variance_ratio'`
>        - *"Number of PCs to show"*: `50`
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

## Sub-step with **Scanpy ComputeGraph**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Scanpy ComputeGraph** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy RunPCA** {% icon tool %})
>    - *"Use programme defaults"*: `Yes`
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

## Sub-step with **Scanpy FindCluster**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Scanpy FindCluster** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy ComputeGraph** {% icon tool %})
>    - *"Use programme defaults"*: `Yes`
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

## Sub-step with **Scanpy RunUMAP**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Scanpy RunUMAP** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy FindCluster** {% icon tool %})
>    - *"Use programme defaults"*: `Yes`
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

## Sub-step with **Scanpy RunTSNE**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Scanpy RunTSNE** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy RunUMAP** {% icon tool %})
>    - *"Use programme defaults"*: `Yes`
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

## Sub-step with **Scanpy FindMarkers**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Scanpy FindMarkers** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy RunTSNE** {% icon tool %})
>    - *"Use programme defaults"*: `Yes`
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

## Sub-step with **Scanpy PlotEmbed**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Scanpy PlotEmbed** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy FindMarkers** {% icon tool %})
>    - *"name of the embedding to plot"*: `pca`
>    - *"color by attributes, comma separated texts"*: `louvain,genotype,batch,sex,Nusap1`
>    - *"Field for gene symbols"*: `Symbol`
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

## Sub-step with **Scanpy PlotEmbed**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Scanpy PlotEmbed** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy FindMarkers** {% icon tool %})
>    - *"color by attributes, comma separated texts"*: `louvain,genotype,batch,sex,Nusap1`
>    - *"Field for gene symbols"*: `Symbol`
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

## Sub-step with **Scanpy PlotEmbed**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Scanpy PlotEmbed** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `output_h5ad` (output of **Scanpy FindMarkers** {% icon tool %})
>    - *"name of the embedding to plot"*: `tsne`
>    - *"color by attributes, comma separated texts"*: `louvain,genotype,batch,sex,Nusap1`
>    - *"Field for gene symbols"*: `Symbol`
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

## Sub-step with **Join two Datasets**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Join two Datasets** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Join"*: `output_tsv` (output of **Scanpy FindMarkers** {% icon tool %})
>    - *"using column"*: `cColumn: 4`
>    - {% icon param-file %} *"with"*: `var` (output of **Inspect AnnData** {% icon tool %})
>    - *"and column"*: `cColumn: 2`
>    - *"Keep lines of first input that do not join with second input"*: `Yes`
>    - *"Keep lines of first input that are incomplete"*: `Yes`
>    - *"Fill empty columns"*: `No`
>    - *"Keep the header lines"*: `Yes`
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

## Sub-step with **Cut**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Cut** {% icon tool %} with the following parameters:
>    - *"Cut columns"*: `c1,c2,c3,c4,c11,c5,c6,c7,c8`
>    - {% icon param-file %} *"From"*: `out_file1` (output of **Join two Datasets** {% icon tool %})
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
