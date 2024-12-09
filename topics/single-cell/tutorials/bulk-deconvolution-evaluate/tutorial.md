---
layout: tutorial_hands_on

title: Evaluating Reference Data for Bulk RNA Deconvolution
subtopic: deconvo
priority: 3
zenodo_link: ''

questions:
- How do I evaluate my reference data?
- How do I compare different deconvolution tools?
- What are the best metrics for determining tool accuracy?
objectives:
- Generate psuedo-bulk data from single-cell RNA data
- Process the single-cell and psuedo-bulk data using various deconvolution tools
- Evaluate and visualse the results of the different deconvolution methods
time_estimation: 2H
key_points:
- It is important to validate the accuracy of both deconvolution tools and reference data
- There are various visualisation and quantative methods of analysing results
- Different deconvolution tools have varying accuracy, comparing them against the same reference is a useful test to determine the best tool for your data

tags:
- transcriptomics

contributions:
  authorship:
    - hexhowells
  funding:
    - elixir-fair-data

follow_up_training:
  -
    type: "internal"
    topic_name: single-cell

---

There are various methods to estimate the proportions of cell types in bulk RNA data. Since the actual cell proportions of the data are unknown, how do we know if our tools are producing accurate results?

In this tutorial we will be using single-cell data with known cell-type proportions in order to create pseudo-bulk RNA data. We will then estimate the cell-type proportions of this pseudo-bulk data using the currently available deconvolution tools within Galaxy. Since we know the true proportion values, we will be able to measure and compare the accuracy of the tools' predictions.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get the single-cell data

First we need to create a new history in Galaxy and load in our single-cell data. We are going to use the single-cell dataset from a previous deconvolution tutorial found here: [https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/bulk-music/tutorial.html]({% link topics/single-cell/tutorials/bulk-music/tutorial.md %}).

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial *"Deconvolution: Evaluating Reference Data"*
> 2. Import the files from [Zenodo](https://zenodo.org/records/5719228) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    * Human pancreas single-cell RNA datasets (tag: `#scrna`)
>      ```
>      https://zenodo.org/record/5719228/files/EMTABesethealthy.expression.tabular
>      https://zenodo.org/record/5719228/files/EMTABesethealthy.phenotype.tabular
>      ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 4. Check the datatypes are `tabular`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular" %}
>
> 5. Add a `#metadata` tag to `EMTABesethealthy.phenotype.tabular` and a `#expression` tag to `EMTABesethealthy.expression.tabular`
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

## Inspecting the single-cell data

Before continuing lets quickly inspect our single-cell data. We can find all of the cell types present in the data alongside their proportions by using the count tool to count the occurrence of each cell type category in the metadata file.

> <hands-on-title>Get cell counts</hands-on-title>
>
> 1. {% tool [Count](Count1) %} with the following parameters:
>    - {% icon param-file %} *"from dataset"*: `EMTABesethealthy.phenotype.tabular`
>    - *"Count occurrences of values in column(s)"*: `Column 5`
>    - *"Delimited by"*: `Tab`
>    - *"How should the results be sorted?"*: `With the most common value first`
>
> 2. **Rename** {% icon galaxy-pencil %} output `Cell type counts`
>
{: .hands_on}

We can see from the output table below, there are various cell types present in the data. Note that many of the cell types have very low proportion values, this should be kept in mind later on as cell types that appear only a hand full of times (or even just once!) in the data may not be very useful and only add noise. 

| Cell Type               | Count |
|-------------------------|-------|
| alpha                   | 443   |
| beta                    | 171   |
| ductal                  | 135   |
| acinar                  | 112   |
| gamma                   | 75    |
| delta                   | 59    |
| unclassified endocrine  | 29    |
| co-expression           | 26    |
| PSC                     | 23    |
| endothelial             | 13    |
| epsilon                 | 5     |
| mast                    | 4     |
| unclassified            | 1     |
| MHC class II            | 1     |

> <question-title>Inspecting the single-cell data</question-title>
>
> 1. How many cells are in the single-cell data?
> 2. How many cell types are present in the data?
>
> > <solution-title></solution-title>
> >
> > 1. Inspecting the general information of `EMTABesethealthy.expression.tabular` we can see that there are **1,097** cells in the data as there are 1,098 columns (we need to subtract 1 for the header).
> > 2. Looking at the output of the {% tool [Count](Count1) %} tool (or the above table), there are **14** distinct cell types in the data.
> >
> {: .solution}
>
{: .question}

# Process the single-cell data

In order to get a good understanding of the accuracy of our deconvolution tools, we are going to run our evaluations multiple times. This approach ensures that a single good or bad evaluation does not disproportionately represent the tool's overall performance.

However, instead of running all of our tools multiple times for each evaluation (which would be quite time consuming!), we will leverage "batch computation" in Galaxy. By storing our data in collections, any tools or workflows that use those collections will automatically run multiple times (once for each element in the collection). We will now perform some pre-processing of our data to get it into the right format.

## Transpose expression matrix

If we inspect the expression data file downloaded earlier, we can see that currently the rows represent genes and columns represent cells. However, this needs to be swapped for the later workflows. To fix this we will transpose the expression matrix.

> <hands-on-title>Transpose expression matrix</hands-on-title>
>
> 1. {% tool [Transpose](toolshed.g2.bx.psu.edu/repos/iuc/datamash_transpose/datamash_transpose/1.8+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input tabular dataset"*: `EMTABesethealthy.expression.tabular`
>
> 2. **Rename** {% icon galaxy-pencil %} output `Transposed expression matrix`
>
{: .hands_on}

## Generate batch mode collections

For this tutorial we will run the evaluations **20** times, this will both help improve the sample size and allow us to determine the consistency of the tools, whilst being small enough to run in a reasonable amount of time.

We will now duplicate our single-cell data 20 times and store it in a collection. This will be done for both the expression data and metadata files.

> <hands-on-title>Generate collections from data</hands-on-title>
>
> 1. {% tool [Duplicate file to collection](__DUPLICATE_FILE_TO_COLLECTION__) %} with the following parameters:
>    - {% icon param-file %} *"Input Dataset"*: `EMTABesethealthy.phenotype.tabular`
>    - *"Size of output collection"*: `20`
>
> 2. **Rename** {% icon galaxy-pencil %} output `Metadata`
>
> 3. {% tool [Duplicate file to collection](__DUPLICATE_FILE_TO_COLLECTION__) %} with the following parameters:
>    - {% icon param-file %} *"Input Dataset"*: `Transposed expression matrix`
>    - *"Size of output collection"*: `20`
>
> 4. **Rename** {% icon galaxy-pencil %} output `Expression data`
>
{: .hands_on}

## Generate expression set objects

Next we will need to use the single-cell data to build an expression set object, this will be used later in the evaluation when we perform the actual deconvolution. 

**Note: We are using the original imported data here, not the transposed data or collections.**

> <hands-on-title>Build the Expression Set object</hands-on-title>
>
> 1. {% tool [Construct Expression Set Object](toolshed.g2.bx.psu.edu/repos/bgruening/music_construct_eset/music_construct_eset/0.1.1+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Assay Data"*: `EMTABesethealthy.expression.tabular` (Input dataset)
>    - {% icon param-file %} *"Phenotype Data"*: `EMTABesethealthy.phenotype.tabular` (Input dataset)
>
>    > <comment-title></comment-title>
>    >
>    > An ExpressionSet object has many data slots, the principle of which are the experiment data (*exprs*), the phenotype data (*pData*), as well metadata pertaining to experiment information and additional annotations (*fData*).
>    {: .comment}
>
{: .hands_on}

Similar to the expression data, this ExpressionSet object needs to be duplicated 20 times into a collection for later batch processing.

> <hands-on-title>Generate ESet collection</hands-on-title>
>
> 1. {% tool [Duplicate file to collection](__DUPLICATE_FILE_TO_COLLECTION__) %} with the following parameters:
>    - {% icon param-file %} *"Input Dataset"*: `RData ESet Object` (output of **Construct Expression Set Object** {% icon tool %})
>    - *"Size of output colection"*: `20`
>
> 2. **Rename** {% icon galaxy-pencil %} output `ESet Object`
>
{: .hands_on}


# Create pseudo-bulk and actual cell proportions

We are now going to run our first workflow! This workflow will extract a subsample from the data containing 200 random cells. The workflow will then perform two things with this subsample:

1. Count the cell types and proportions of the data in order to be used as reference later against the predicted proportion values
2. Remove the cell types and convert the single-cell data into pseudo-bulk data to be later inputted into the deconvolution tools.

The above will be done twice to emulate multiple "subjects". Since the deconvolution tools will be expecting the bulk-RNA data to comprise of at least 2 subjects (each with their own bulk data). For this tutorial our subjects will simply be called **A** and **B**. However, in the real world these subjects could be different patients, tissue samples, diseased/healthy, etc. 

> <comment-title>Different Results</comment-title>
> Note that since we are selecting 20 samples, each containing 200 randomly selected cells. The plots and results presented in this tutorial will differ from your own. There will be some similarities such as certain cells being in higher proportion to others but the exact values with differ!
{: .comment}

**Remember** since we have a collection of 20 inputs, the output of this workflow will be a collection of 20 elements, each corresponding to the input elements. Each output will have its own random selection of 200 cells.

> <comment-title>Inputting multiple datasets</comment-title>
>
> In order to upload the input collections into the workflow, you first need to set the input type to **Multiple datasets** in the input file selection.
> 
> ![Multiple Datasets](../../images/bulk-deconvolution-evaluate/batch-mode.png "Multiple Datasets button in Galaxy")
{: .comment}

> <hands-on-title>Run pseudo-bulk and actual proportions workflow</hands-on-title>
>
> 1. **Import the workflow** into Galaxy
>    - Copy the URL (e.g. via right-click) of [this workflow](https://usegalaxy.eu/u/hexhowells/w/deconv-eval-stage-1) or download it to your computer.
>    - Import the workflow into Galaxy
>
>    {% snippet faqs/galaxy/workflows_run_trs.md path="topics/transcriptomics/tutorials/rna-seq-reads-to-counts/workflows/qc_report.ga" title="QC Report" %}
>
> 2. Run **Workflow pseudobulk and actual proportions** {% icon workflow %} using the following parameters:
>    - {% icon param-collection %} *"Metadata"*: `Metadata`
>    - {% icon param-collection %} *"Expression Data"*: `Expression Data`
>
>    {% snippet faqs/galaxy/workflows_run.md %}
> 3. Add a tag labelled `#A` to the first "Actual cell proportions" and "Pseudobulk" collections
>
> 4. Add a tag labelled `#B` to the second "Actual cell proportions" and "Pseudobulk" collections
{: .hands_on}

<iframe title="Galaxy Workflow Embed" style="width: 100%; height: 700px; border: none;" src="https://usegalaxy.eu/published/workflow?id=cb27f805d076ee9f&embed=true&buttons=true&about=false&heading=false&minimap=true&zoom_controls=true&initialX=-20&initialY=-20&zoom=0.5"></iframe>

The output of this workflow will be the psuedo-bulk and actual cell proportions for both samples A and B. If you inspect one of the elements in the `Actual Cell Proportions` collection, you should see a table similar to the following:

|                         | A_actual   |
|-------------------------|------------|
| acinar                  | 0.090000   |
| alpha                   | 0.415000   |
| beta                    | 0.170000   |
| co-expression           | 0.050000   |
| delta                   | 0.070000   |
| ductal                  | 0.105000   |
| endothelial             | 0.015000   |
| gamma                   | 0.050000   |
| mast                    | 0.010000   |
| unclassified endocrine  | 0.025000   |

Comparing the above table with the cell-type counts of the original single-cell data, does this look correct? Well the top 3 cell-types with the highest proportion in the single-cell data are: alpha, beta, ductal. Which aligns with the proportion values of the above data! There may be some variance due to the randomly selected cells. Also note that some of the lesser common cell types (like `MHC class II`) aren't present in the above table, again this is due to the 200 randomly selected cells for this specific sample and isn't of concern.

# Perform deconvolution on the pseudo-bulk data

Now that we have our pseudo-bulk data alongside the actual proportion values. Our next step is to run deconvolution to get predicted cell-type proportions! Currently, Galaxy contains two tools for performing deconvolution: **MuSiC** and **NNLS**. We will use both of these tools in this tutorial and compare their results together.

The following workflow will take the two pseudo-bulk samples (A and B), as well as the original single-cell data as reference and output the deconvolution results for both samples and deconvolution methods. Thus producing 4 output collections. The pdf results of the deconvolution tools will also be outputted from the workflow but won't be needed for the tutorial.

> <hands-on-title>Run inferring cellular proportions workflow</hands-on-title>
>
> 1. **Import the workflow** into Galaxy
>    - Copy the URL (e.g. via right-click) of [this workflow](https://usegalaxy.eu/u/hexhowells/w/deconv-eval-stage-2) or download it to your computer.
>    - Import the workflow into Galaxy
>
>    {% snippet faqs/galaxy/workflows_run_trs.md path="topics/transcriptomics/tutorials/rna-seq-reads-to-counts/workflows/qc_report.ga" title="QC Report" %}
>
> 2. Run **Workflow inferring cellular proportions** {% icon workflow %} using the following parameters:
>    - {% icon param-collection %} *"Pseudobulk - A"*: `expression data - A`
>    - {% icon param-collection %} *"Pseudobulk - B"*: `expression data - B`
>    - {% icon param-collection %} *"ESet Reference scRNA-seq"*: `ESet Object`
>    - *"Cell Types Label from scRNA dataset"*: `cellType`
>    - *"Samples Identifier from scRNA dataset"*: `sampleID`
>    - *"Cell types to use from scRNA dataset"*:`alpha,beta,ductal,acinar,gamma,delta,unclassified endocrine,co-expression,PSC,endothelial,epsilon,mast,unclassified,MHC class II`
>    - {% icon param-collection %} *"Actual - B"*: `actual - B`
>    - {% icon param-collection %} *"Actual - A"*: `actual - A`
>
>    {% snippet faqs/galaxy/workflows_run.md %}
> 3. Rename output collections / add tags?
{: .hands_on}

<iframe title="Galaxy Workflow Embed" style="width: 100%; height: 700px; border: none;" src="https://usegalaxy.eu/published/workflow?id=6a47b42bf753aa3a&embed=true&buttons=true&about=false&heading=false&minimap=true&zoom_controls=true&initialX=-20&initialY=-20&zoom=0.5"></iframe>

# Visualise results

Now that we have our deconvolution results, the next step is to analyse the predictions and determine how accurate our tools are given our reference data. Since our pseudo-subjects **A** and **B** come from the same data, there isn't much point inspecting them both. So for the rest of the tutorial we will just focus our analysis on subject **A**.

In order to determine if our tools have produced accurate results, we will create various plots and compute different metrics to visualise and quantify the outputs of our tools.

## Pre-process the output results

Before visualising or inspecting the outputs of the deconvolution tools, we first need to perform some pre-processing. Up until now we have been working with collections in order to perform our evaluations multiple times in parallel. However, for analysing our data, collections will be a bit messy and are no longer needed. The following workflow will combine all the collections of the MuSiC and NNLS outputs into two tables:

1. A results table presenting the predicted and actual proportion values of each cell-type of each subsample
2. An error table showing the difference between the actual and predicted values. Which will be needed for a later plot.

> <hands-on-title>Run visualisation pre-processing workflow</hands-on-title>
>
> 1. **Import the workflow** into Galaxy
>    - Copy the URL (e.g. via right-click) of [this workflow](https://usegalaxy.eu/u/hexhowells/w/deconv-eval-stage-3-process) or download it to your computer.
>    - Import the workflow into Galaxy
>
>    {% snippet faqs/galaxy/workflows_run_trs.md path="topics/transcriptomics/tutorials/rna-seq-reads-to-counts/workflows/qc_report.ga" title="QC Report" %}
>
> 2. Run **Workflow preprocess visualisations** {% icon workflow %} using the following parameters:
>    - {% icon param-collection %} *"Cell Proportions"*: `A - Music Results`
>
> 3. Run **Workflow preprocess visualisations** {% icon workflow %} using the following parameters:
>    - {% icon param-collection %} *"Cell Proportions"*: `B - NNLS Results`
>
>    {% snippet faqs/galaxy/workflows_run.md %}
{: .hands_on}

<iframe title="Galaxy Workflow Embed" style="width: 100%; height: 700px; border: none;" src="https://usegalaxy.eu/published/workflow?id=76d3408d0d22ad05&embed=true&buttons=true&about=false&heading=false&minimap=true&zoom_controls=true&initialX=-20&initialY=-20&zoom=0.5"></iframe>

The following table shows a snippet of the `Results Table` for the MuSiC tool. A header has been added for better reading but has been omitted in the workflow output as it will interfere with the visualisation tools.

| Cell Type      | Actual Proportion   | Predicted Proportion   |
|----------------|---------------------|------------------------|
| acinar         | 0.090000            | 0.0814442584577275     |
| alpha          | 0.415000            | 0.427718807911522      |
| beta           | 0.170000            | 0.256954867012044      |
| co-expression  | 0.050000            | 0                      |
| delta          | 0.070000            | 0.0929840465107452     |
| ...            | ...                 | ...                    |

Already at first glance we can see some interesting results! Firstly we can see that the tool is able to make predictions close to the actual values such as with `acinar, alpha, delta`. We also see the tool failing to make any type of prediction for `co-expression` cells with a predicted proportion value of 0. This however isn't a compete surprise since `co-expression` cells are of small proportion in the bulk and reference data. 

But this is only a small sample of the results. Lets create some visualisations to see the whole picture!

## Plot scatter plots of the results

The first type of visualisation we will do is a scatter plot. This plot will compare the actual and predicted proportion values for each cell across each subsample. We will also colour each point on the plot to indicate which cell type it belongs to. Let's do that now for both the MuSiC and NNLS results.

{% snippet faqs/galaxy/tools_rerun.md %}

> <hands-on-title>Plot the actual and inferred data</hands-on-title>
>
> 1. {% tool [Scatterplot with ggplot2](toolshed.g2.bx.psu.edu/repos/iuc/ggplot2_point/ggplot2_point/3.4.0+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input in tabular format"*: `Results Table (Music)`
>    - *"Column to plot on x-axis"*: `2`
>    - *"Column to plot on y-axis"*: `3`
>    - *"Plot title"*: `Correlation between inferred and actual cell-type proportions`
>    - *"Label for x axis"*: `Actual proportions`
>    - *"Label for y axis"*: `Inferred proportions`
>    - In *"Advanced options"*:
>       - *"Plotting multiple groups"*: `Plot multiple groups of data on one plot`
>           - *"column differentiating the different groups"*: `1`
>           - *"Color schemes to differentiate your groups"*: `Paired - predefined color pallete (discrete, max=12 colors)`
>           - *"Reverse color scheme"*: `Default order of color scheme`
>    - In *"Output options"*:
>       - *"width of output"*: `6.0`
>       - *"height of output"*: `4.0`
>
> 2. **Rename** {% icon galaxy-pencil %} output `Scatter plot - Music`
>
> 3. Add a `#plot` tag to `Scatter plot - Music`
>
> 4. {% tool [Scatterplot with ggplot2](toolshed.g2.bx.psu.edu/repos/iuc/ggplot2_point/ggplot2_point/3.4.0+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input in tabular format"*: `Results Table (NNLS)`
>    - *"Column to plot on x-axis"*: `2`
>    - *"Column to plot on y-axis"*: `3`
>    - *"Plot title"*: `Correlation between inferred and actual cell-type proportions`
>    - *"Label for x axis"*: `Actual proportions`
>    - *"Label for y axis"*: `Inferred proportions`
>    - In *"Advanced options"*:
>       - *"Plotting multiple groups"*: `Plot multiple groups of data on one plot`
>           - *"column differentiating the different groups"*: `1`
>           - *"Color schemes to differentiate your groups"*: `Paired - predefined color pallete (discrete, max=12 colors)`
>           - *"Reverse color scheme"*: `Default order of color scheme`
>    - In *"Output options"*:
>       - *"width of output"*: `6.0`
>       - *"height of output"*: `4.0`
>
> 5. **Rename** {% icon galaxy-pencil %} output `Scatter plot - NNLS`
>
> 6. Add a `#plot` tag to `Scatter plot - NNLS`
>
{: .hands_on}

The output of this tool should produce two scatter plots that looks like the image below. Each point on the plot represents a cell-type for a specific subsample, so there should be 20 points of each colour (one for each subsample created earlier). Since we are comparing the actual and inferred proportions, the ideal scatter plot would have all of the points be at the `y=x` line. The further the deviations are from this ideal line, the less accurate the tool is. We can also use this plot to determine if the tool is under or over predicting proportion values for each cell-type, or if the tool is struggling to predict certain cell types.

![Scatter plot MuSiC](../../images/bulk-deconvolution-evaluate/scatterplot-music.png "Scatter plot of Music results")

> <question-title>Interpreting the Scatter Plots</question-title>
>
> 1. Which method has the most accurate results?
> 2. Which cell type has the biggest proportion in the dataset?
> 3. Do either of the tools struggle with any cell types?
>
> {% snippet faqs/galaxy/features_scratchbook.md datatype="tabular" %}
>
> > <solution-title></solution-title>
> >
> > ![Scatter plot comparison](../../images/bulk-deconvolution-evaluate/scatterplot-compare.png "Scatter plot comparison between Music and NNLS")
> >
> > 1. Comparing scatter plots, the MuSiC tool has the most accurate results since the points fall closer onto the x=y line
> > 2. Both scatter plots show `alpha` cells having the highest proportion by a large margin
> > 3. The MuSiC tool seems to handle all cell types well. However, NNLS appears to struggle predicting the proportions of beta cells, with many of the samples being predicted as having a proportion of 0
> >
> {: .solution}
>
{: .question}


## Plot violin plots of the errors

Next we will plot the distribution of errors between the predicted and actual cellular proportions for a select number of cell types. We could plot all cell types in the output, however too many will cause the visualisations to be messy and difficult to interpret.

We can use the cell-type counts we computed at the beginning of the tutorial to determine the best cell types to use. We will use the top 5 most abundant cell types in the single-cell data being: `alpha, beta, gamma, ductal, acinar`. Before plotting we will extract only these cell types from our table of errors.

{% snippet faqs/galaxy/tools_rerun.md %}

> <hands-on-title>Extract Cell Types</hands-on-title>
>
> 1. {% tool [Advanced Cut](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/9.3+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"File to cut"*: `Error Table (Music)`
>    - *"Operation"*: `Keep`
>    - *"Cut by"*: `fields`
>       - *"Delimited by"*: `Tab`
>       - *"Is there a header for the data's columns ?"*: `Yes`
>           - *"List of Fields"*: `Select the columns containing: alpha, beta, gamma, ductal, acinar`
>
> 2. **Rename** {% icon galaxy-pencil %} output `Music Errors`
>
> 3. {% tool [Advanced Cut](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/9.3+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"File to cut"*: `Error Table (NNLS)`
>    - *"Operation"*: `Keep`
>    - *"Cut by"*: `fields`
>       - *"Delimited by"*: `Tab`
>       - *"Is there a header for the data's columns ?"*: `Yes`
>           - *"List of Fields"*: `Select the columns containing: alpha, beta, gamma, ductal, acinar`
>
> 4. **Rename** {% icon galaxy-pencil %} output `NNLS Errors`
>
{: .hands_on}

Now we have our table of errors consisting of only the top 5 cell-types, we can plot the violin plots.

{% snippet faqs/galaxy/tools_rerun.md %}

> <hands-on-title>Plot violin plots</hands-on-title>
>
> 1. {% tool [Violin plot w ggplot2](toolshed.g2.bx.psu.edu/repos/iuc/ggplot2_violin/ggplot2_violin/3.4.0+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input in tabular format"*: `Music Errors`
>    - *"Plot title"*: `Error Distribution`
>    - *"Label for x axis"*: `Cell Type`
>    - *"Label for y axis"*: `Difference Error`
>    - In *"Advanced Options"*:
>       - *"Violin border options"*: `Purple`
>    - In *"Output Options"*:
>       - *"width of output"*: `3.0`
>       - *"height of output"*: `2.0`
>
> 2. **Rename** {% icon galaxy-pencil %} output `Violin Plot - Music`
>
> 3. Add a `#plot` tag to `Violin Plot - Music`
>
> 4. {% tool [Violin plot w ggplot2](toolshed.g2.bx.psu.edu/repos/iuc/ggplot2_violin/ggplot2_violin/3.4.0+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input in tabular format"*: `NNLS Errors`
>    - *"Plot title"*: `Error Distribution`
>    - *"Label for x axis"*: `Cell Type`
>    - *"Label for y axis"*: `Difference Error`
>    - In *"Advanced Options"*:
>       - *"Violin border options"*: `Purple`
>    - In *"Output Options"*:
>       - *"width of output"*: `3.0`
>       - *"height of output"*: `2.0`
>
> 5. **Rename** {% icon galaxy-pencil %} output `Violin Plot - NNLS`
>
> 6. Add a `#plot` tag to `Violin Plot - NNLS`
>
{: .hands_on}

The output of this tool will be two violin plots that will look similar to the below image. Here we can see the distribution of errors for each cell type. Since we are using normal errors and not absolute or squared errors, we are also able to see whether the tool has under or over estimated the cell type. An ideal plot would have all the violin plots being short in height and close to 0 indicating that the estimated and actual values are close together (resulting in an error close to 0).

![Violin plot MuSiC](../../images/bulk-deconvolution-evaluate/violin-music.png "Violin plot of Music results")

> <question-title>Interpreting the Violin Plots</question-title>
>
> 1. Which method has the least errors?
> 2. Which method is the most balanced when to over and under estimating proportions?
> 3. Which is the most overestimated cell type in NNLS?
>
> {% snippet faqs/galaxy/features_scratchbook.md datatype="tabular" %}
>
> > <solution-title></solution-title>
> >
> > ![Scatter plot comparison](../../images/bulk-deconvolution-evaluate/violin-compare.png "Scatter plot comparison between Music and NNLS")
> >
> > 1. Comparing the two violin plots, MuSiC has the better error results, with more samples closer to zero. Inspecting the y-axis of the plots also show that the MuSiC errors span a smaller range compared to NNLS.
> > 2. MuSiC can be seen as having the most balanced results with the bulk of the estimates being around 0. Whereas the NNLS results show large amounts of both overestimation and underestimation of various cell types.
> > 3. From the NNLS violin plot it can be seen that ductal cells are greatly overestimated.
> >
> {: .solution}
>
{: .question}


# Compute accuracy metrics

Visualisations are a great tool for getting an intuitive overview of the data. However, some of the interpretations from visualisations can be subjective. Having quantitative results alongside visualisations can offer concrete and precise values about the data that can more easily be compared. We will use two different quantitative metrics in this tutorial; Pearson correlation and RMSE.

## Pearson Correlation

The Pearson correlation coefficient is a statistical value that represents the direction and correlation between two variables, the value of this metric ranges between -1 and 1, where:

- -1 = negative correlation
- 0 = no correlation
- 1 = positive correlation

The equation for calculating the Pearson correlation can be seen below, the workflow to compute this metric breaks down this formula into smaller steps.

![Pearson Correlation Equation](../../images/bulk-deconvolution-evaluate/pearson.png "Pearson Correlation Equation")

Where
- `x` = actual proportion values
- `x̄` = mean of actual proportion values
- `y` = predicted proportion values
- `ȳ` = mean of predicted proportion values

## Root Mean Squared Error (RMSE)

Root Mean Squared Error or RMSE is a common metric for measuring a tools prediction error. This metric calculates the average error between the predicted and actual values for each prediction then takes the mean and square root of the error to produce a final value. Lower RMSE values (close to 0) indicate accurate predictions similar to the actual value, as the value increases the accuracy score worsens.

The equation for calculating this metric is seen below, the implementation of this calculation is in the workflow alongside the Pearson correlation.

![Root Mean Squared Error Equation](../../images/bulk-deconvolution-evaluate/rmse.png "Root Mean Squared Error Equation")

Where
- `n` = number of samples
- `y` = actual proportion
- `ŷ` = predicted proportion


## Compute metrics

With a basic understanding of some useful metrics, we will now compute these to get quantitative values alongside our visualisation results. The following workflow needs to be run for both the MuSiC and NNLS results table.

> <hands-on-title>Run metrics workflow</hands-on-title>
>
> 1. **Import the workflow** into Galaxy
>    - Copy the URL (e.g. via right-click) of [this workflow](https://usegalaxy.eu/u/hexhowells/w/deconv-eval-stage-3-metrics) or download it to your computer.
>    - Import the workflow into Galaxy
>
>    {% snippet faqs/galaxy/workflows_run_trs.md path="topics/transcriptomics/tutorials/rna-seq-reads-to-counts/workflows/qc_report.ga" title="QC Report" %}
>
> 2. Run **Workflow compute metrics** {% icon workflow %} using the following parameters:
>    - {% icon param-collection %} *"Cell Proportions"*: `Results Table (Music)`
>
> 3. Run **Workflow compute metrics** {% icon workflow %} using the following parameters:
>    - {% icon param-collection %} *"Cell Proportions"*: `Results Table (NNLS)`
>
>    {% snippet faqs/galaxy/workflows_run.md %}
{: .hands_on}

<iframe title="Galaxy Workflow Embed" style="width: 100%; height: 700px; border: none;" src="https://usegalaxy.eu/published/workflow?id=7caa58e36df6ad03&embed=true&buttons=true&about=false&heading=false&minimap=true&zoom_controls=true&initialX=-20&initialY=-20&zoom=0.5"></iframe>

After running the workflow on both the MuSiC and NNLS results we should have the Pearson and RMSE metrics for both tools in various outputs. Below combines these metrics into a single summary table.

| Tool  | Pearson Correlation | RMSE  |
|-------|---------------------|-------|
| MuSiC | 0.982               | 0.022 |
| NNLS  | 0.778               | 0.678 |

From the table we can now see concrete values representing the error and correlation between the predictions and actual proportion values. We can see from the table that the MuSiC tool has a much better accuracy with a higher correlation score and lower error compared to NNLS.

The conclusions to draw from this analysis, is that our reference data is effective for use in deconvolution analysis since both tools were able to have high accuracy and low error scores. We also determined that (for at least this data) the MuSiC tool was the more effective/accurate tool and thus would likely be the more trustworthy when performing deconvolution with this single-cell reference data.

# Conclusion

{% icon congratulations %} Congratulations! You made it to the end of the tutorial!

In this tutorial we took some single-cell data with known cell-type proportions, subsampled the data, and converted them to pseudo-bulk data. We then used this pseudo-bulk data to perform deconvolution using the two tools available in Galaxy: MuSiC and NNLS. Using the known cell-type proportions we were able to analyse the predicted proportions to the ground truth in order to determine if the reference data can be used and which tool is the most effective. We used various visualisation and statistical techniques to analyse and quantify the tools accuracy, reliability, and error.
