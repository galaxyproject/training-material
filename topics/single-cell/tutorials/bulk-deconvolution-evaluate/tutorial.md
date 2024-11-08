---
layout: tutorial_hands_on

title: Evaluating Reference Data for Bulk RNA Deconvolution
subtopic: deconvo
priority: 3
zenodo_link: ''

questions:
- How do I evaluate my reference data
- How do I compare different deconvolution tools
- What are the best metrics for determining tool accuracy
objectives:
- Generate psuedo-bulk data from single-cell RNA data.
- Process the single-cell and psuedo-bulk data using various deconvolution tools
- Evaluate and visualse the results of the different deconvolution methods
time_estimation: 2H
key_points:
- Something about deconv tools having different accuracies
- Its important to validate the accuracy of tools 

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

There are various methods to estimate the proportions of cell types in bulk RNA data. Since the actual proportions of the data are unknown, how do we know if our tools are producing accurate results?

In this tutorial we will be using single-cell data with known cell proportions in order to create' pseudo' bulk RNA data. We will then estimate the proportions of this this pseudo-bulk data using the currently available deconvolution tools within Galaxy. Since we know the true proportions values, we will be able to measure and compare the accuracy of the tools.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get single-cell data

First we need to create a new history in Galaxy and load our single-cell data We are going to use the single-cell dataset found in a previous deconvolution tutorial [https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/bulk-music/tutorial.html](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/bulk-music/tutorial.html).

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial *"Deconvolution: Evaluating Reference Data"*
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
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
> 3. Rename the datasets
>
> 4. Check the datatype
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular" %}
>
> 5. Add a `#metadata` tag to `EMTABesethealthy.phenotype.tabular` and a `#expression` tag to `EMTABesethealthy.expression.tabular`
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}



# Process single-cell data

In order to get a good understanding of the accuracy of our deconvolution tools, we are going to run our evaluations multiple times. This is to ensure that a single good or bad evaluation doesn't become indicative of the overall tool's performance.

However, instead of running all of our tools multiple times for each evaluation (which would be quite time consuming!), we will leverage "batch computation" in Galaxy. By storing our data in collections, any tools or workflows that use with that data will run individually for each element in the collection. We will now perform some pre-processing of our data to get it into the right format.

## Transpose Expression Matrix

If we inspect the expression data file downloaded earlier, we can see that currently the rows represent genes and columns represent cells. However, this needs to be swapped for the later workflows. To fix this we will transpose the expression matrix.

> <hands-on-title>Transpose expression matrix</hands-on-title>
>
> 1. {% tool [Transpose](toolshed.g2.bx.psu.edu/repos/iuc/datamash_transpose/datamash_transpose/1.8+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input tabular dataset"*: `EMTABesethealthy.expression.tabular`
>
> 2. **Rename** {% icon galaxy-pencil %} output `Transposed expression matrix`
>
{: .hands_on}

## Generate Batch Mode Collections

For this tutorial we will run the evaluations **20** times, this will both help improve the sample size and allow us to determine the consistency of the tools, whilst being small enough to run in a reasonable amount of time!

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

## Generate Expression Set Objects

Next we will need to use the single-cell data to build and expression set object, this will be used later in the evaluation when we perform the actual deconvolution. 

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
>    - {% icon param-file %} *"Input Dataset"*: `ESet Object` (output of **Construct Expression Set Object** {% icon tool %})
>    - *"Size of output colection"*: `20`
>
> 2. **Rename** {% icon galaxy-pencil %} output `ESet Object`
>
{: .hands_on}


# Create pseudo-bulk and actual cell proportions

We are now going to run our first workflow! This workflow will extract a subsample from the data containing 200 cells. This data will be used to generate our pseudo-bulk data along with the actual cell proportions used to evaluate/compare with the output of the deconvolutional tools.

> <comment-title>Inputting Multiple Datasets</comment-title>
>
> In order to upload the input collections into the workflow, you first need to set the input type to **Multiple datasets** in the input file selection.
{: .comment}

> <hands-on-title>Run pseudobulk and actual proportions workflow</hands-on-title>
>
> 1. **Import the workflow** into Galaxy
>    - Copy the URL (e.g. via right-click) of [this workflow]({{ site.baseurl }}{{ page.dir }}workflows/qc_report.ga) or download it to your computer.
>    - Import the workflow into Galaxy
>
>    {% snippet faqs/galaxy/workflows_run_trs.md path="topics/transcriptomics/tutorials/rna-seq-reads-to-counts/workflows/qc_report.ga" title="QC Report" %}
>
> 2. Run **Workflow pseudobulk and actual proportions** {% icon workflow %} using the following parameters:
>    - {% icon param-collection %} *"Metadata"*: `metadata collection`
>    - {% icon param-collection %} *"Expression Data"*: `expression data collection`
>
>    {% snippet faqs/galaxy/workflows_run.md %}
> 3. Rename output collections / add tags?
>
> 4. Inspect `cell type counts`
{: .hands_on}

Need to identify the list of cells first before running the second workflow, this can be done by looking at the cell type counts output. For this tutorial we will use all of the cell types as this will be the most realistic scenario. However, cell types with counts below a threshold value (e.g. 10) can be removed.

# Perform Deconvolution on the Pseudo-Bulk Data

> <hands-on-title>Run inferring cellular proportions workflow</hands-on-title>
>
> 1. **Import the workflow** into Galaxy
>    - Copy the URL (e.g. via right-click) of [this workflow]({{ site.baseurl }}{{ page.dir }}workflows/qc_report.ga) or download it to your computer.
>    - Import the workflow into Galaxy
>
>    {% snippet faqs/galaxy/workflows_run_trs.md path="topics/transcriptomics/tutorials/rna-seq-reads-to-counts/workflows/qc_report.ga" title="QC Report" %}
>
> 2. Run **Workflow inferring cellular proportions** {% icon workflow %} using the following parameters:
>    - {% icon param-collection %} *"Subsample_A - matrix"*: `expression data - A`
>    - {% icon param-collection %} *"Subsample_B - matrix"*: `expression data - B`
>    - {% icon param-collection %} *"ESet Reference scRNA-seq"*: `ESet Object`
>    - *"Cell Types Label from scRNA dataset"*: `cellType`
>    - *"Samples Identifier from scRNA dataset"*: `sampleID`
>    - *"Cell types to use from scRNA dataset"*:`acinar,alpha,beta,delta,ductal,gamma`
>    - {% icon param-collection %} *"B_actuals"*: `actual - B`
>    - {% icon param-collection %} *"A_actuals"*: `actual - A`
>
>    {% snippet faqs/galaxy/workflows_run.md %}
> 3. Rename output collections / add tags?
>
> 4. Inspect `cell type counts`
{: .hands_on}

# Visualise Results

## Pre-process the output results

> <hands-on-title>Run visualisation pre-processing workflow</hands-on-title>
>
> 1. **Import the workflow** into Galaxy
>    - Copy the URL (e.g. via right-click) of [this workflow]({{ site.baseurl }}{{ page.dir }}workflows/qc_report.ga) or download it to your computer.
>    - Import the workflow into Galaxy
>
>    {% snippet faqs/galaxy/workflows_run_trs.md path="topics/transcriptomics/tutorials/rna-seq-reads-to-counts/workflows/qc_report.ga" title="QC Report" %}
>
> 2. Run **Workflow preprocess visualisations** {% icon workflow %} using the following parameters:
>    - {% icon param-collection %} *"Cell Proportions"*: `Music Results`
>
> 3. Run **Workflow preprocess visualisations** {% icon workflow %} using the following parameters:
>    - {% icon param-collection %} *"Cell Proportions"*: `NNLS Results`
>
>    {% snippet faqs/galaxy/workflows_run.md %}
{: .hands_on}

## Plot scatter plots of the results

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
>       - *"width of output"*: `5.0`
>       - *"height of output"*: `3.0`
>
> 2. **Rename** {% icon galaxy-pencil %} output `Scatter plot - Music`
>
> 3. {% tool [Scatterplot with ggplot2](toolshed.g2.bx.psu.edu/repos/iuc/ggplot2_point/ggplot2_point/3.4.0+galaxy1) %} with the following parameters:
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
>       - *"width of output"*: `5.0`
>       - *"height of output"*: `3.0`
>
> 4. **Rename** {% icon galaxy-pencil %} output `Scatter plot - NNLS`
>
{: .hands_on}

> <question-title>Interpreting the Scatter Plots</question-title>
>
> 1. Which method has the most accurate results?
> 2. Which cell type has the biggest proportion in the dataset?
> 3. Are there any cell types the tools fail to predict?
>
> {% snippet faqs/galaxy/features_scratchbook.md datatype="tabular" %}
>
> > <solution-title></solution-title>
> >
> > ![Scatter plot comparison](../../images/bulk-deconvolution-evaluate/scatterplot-compare.png "Scatter plot comparison between Music and NNLS")
> >
> > 1. Comparing scatter plots, the MuSiC tool has the most accurate results since the points fall closer onto the x=y line
> > 2. Both scatter plots show alpha cells having the highest proportion by a large margin
> > 3. In both plots, only the top 6 cell types are predicted accurately, the rest are predicted as 0. This is likely due to how small the proportions of these cell types are
> >
> {: .solution}
>
{: .question}


## Plot Violin plots of the error

Next we will plot the distribution of errors between the predicted and actual cellular proportions for a select number of cell types. We could plot all cell types in the output, however too many will cause the visualisations to be messy and difficult to interpret.

> <hands-on-title>Get cell counts</hands-on-title>
>
> 1. {% tool [Count](Count1) %} with the following parameters:
>    - {% icon param-file %} *"from dataset"*: `EMTABesethealthy.phenotype.tabular`
>    - *"Count occurrences of values in column(s)"*: `Column 5`
>    - *"Delimited by"*: `Tab`
>    - *"How should the results be sorted?"*: `With the most common value first`
>
> 2. **Rename** {% icon galaxy-pencil %} output `Cell Type Counts`
>
{: .hands_on}

Run Counts on the original metadata file and sort by most frequent:

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

We can see that most cell types have very low proportions, so for this visualisation we will only look at 5 cell types with the highest proportion values. For the above table these cell types are: `alpha, beta, gamma, ductal, acinar`. Before we visualise the data we first need to extract only these cell types from the error table.

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
> 3. {% tool [Violin plot w ggplot2](toolshed.g2.bx.psu.edu/repos/iuc/ggplot2_violin/ggplot2_violin/3.4.0+galaxy1) %} with the following parameters:
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
> 4. **Rename** {% icon galaxy-pencil %} output `Violin Plot - NNLS`
>
{: .hands_on}


# Compute Accuracy Metrics

Visualisations are a great tool for getting an intuitive overview of the data. However, some of the interpretations from visualisations can be subjective. Having quantitative results alongside visualisations can offer concrete and precise values about the data that can more easily be compared. We will use two different quantitative metrics in this tutorial; Pearson correlation and RMSE.

## Pearson Correlation

The Pearson correlation coefficient is a statistical value that represents the direction and correlation between two variables, the value of this metric ranges between -1 and 1, where:

- -1 = negative correlation
- 0 = no correlation
- 1 = positive correlation

The equation for calculating the Pearson correlation can be seen below, the workflow to compute this metric breaks down this formula into smaller steps.

![Pearson Correlation Equation](../../images/bulk-deconvolution-evaluate/pearson.png "Pearson Correlation Equation")


## Root Mean Squared Error (RMSE)

Root Mean Squared Error or RMSE is a common error metric for measuring a tools prediction accuracy. This metric calculates the average error between the predicted and actual values for each prediction then takes the mean and square root of the error to produce a final value. Lower RMSE values (close to 0) indicate accurate predictions similar to the actual value, as the value increases the accuracy score is worse.

The equation for calculating this metric is seen below, the implementation of this calculation is in the workflow alongside the Pearson correlation.

![Root Mean Squared Error Equation](../../images/bulk-deconvolution-evaluate/rmse.png "Root Mean Squared Error Equation")


## Compute Metrics

With a basic understanding of some useful metrics, we will now compute these to get quantitative values alongside our visualisation results. The following workflow needs to be run for both the MuSiC and NNLS results table.

> <hands-on-title>Run visualisation workflow</hands-on-title>
>
> 1. **Import the workflow** into Galaxy
>    - Copy the URL (e.g. via right-click) of [this workflow]({{ site.baseurl }}{{ page.dir }}workflows/qc_report.ga) or download it to your computer.
>    - Import the workflow into Galaxy
>
>    {% snippet faqs/galaxy/workflows_run_trs.md path="topics/transcriptomics/tutorials/rna-seq-reads-to-counts/workflows/qc_report.ga" title="QC Report" %}
>
> 2. Run **Workflow preprocess visualisations** {% icon workflow %} using the following parameters:
>    - {% icon param-collection %} *"Cell Proportions"*: `Music Results`
>
> 3. Run **Workflow preprocess visualisations** {% icon workflow %} using the following parameters:
>    - {% icon param-collection %} *"Cell Proportions"*: `NNLS Results`
>
>    {% snippet faqs/galaxy/workflows_run.md %}
{: .hands_on}

After running the workflow on both the MuSiC and NNLS results we should have the Pearson and RMSE metrics for both tools. The below table summarises these scores.

| Tool  | Pearson Correlation | RMSE  |
|-------|---------------------|-------|
| MuSiC | 0.982               | 0.030 |
| NNLS  | 0.954               | 0.037 |

From the table we can now see concrete values representing the error and correlation between the predictions and actual proportion values. We can see from the table that the MuSiC tool has a better accuracy with a higher correlation score and lower error compared to NNLS.
