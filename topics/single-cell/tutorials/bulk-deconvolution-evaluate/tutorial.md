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

## Transpose Expression Matrix

Inspecting the expression data file, we can see that currently the rows represent genes and columns represent cells. However, this needs to be swapped for the later workflows. So first we will transpose the expression matrix.

> <hands-on-title>Transpose expression matrix</hands-on-title>
>
> 1. {% tool [Transpose](toolshed.g2.bx.psu.edu/repos/iuc/datamash_transpose/datamash_transpose/1.8+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input tabular dataset"*: `EMTABesethealthy.expression.tabular`
>
> 2. **Rename** {% icon galaxy-pencil %} output `Transposed expression matrix`
>
{: .hands_on}

## Generate Batch Mode Collections

In order to run our evaluations 20 times we need to provide our workflow with a collection, doing this will run the workflow in "batch mode" which will run the workflow individually for each element in the collection.

We are going to duplicate our single-cell data 20 times and store it in a collection. This will be done for both the expression data and metadata files.

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

Next we will need to use the single-cell data to build and expression set object, this will be used later in the evaluation when we perform the actual deconvolution. **Note: We are using the original imported data here, not the transposed data or collections.**

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

Similar to the expression data, this object needs to be duplicated 20 times into a collection for later batch processing.

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

Here we are going to run our first workflow, this workflow will extract a subsample from the data containing 200 cells. This data will be used to generate our pseudo-bulk data along with the actual cell proportions used to evaluate/compare with the output of the deconvolutional tools.

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

# Compute Accuracy and Error Metrics

Combine the collection of tables together

> <hands-on-title>Combine output tables</hands-on-title>
>
> 1. {% tool [Column join](toolshed.g2.bx.psu.edu/repos/iuc/collection_column_join/collection_column_join/0.0.3) %} with the following parameters:
>    - {% icon param-file %} *"Tabular files"*: `A proportions actual-infer`
>    - *"Identifier column"*: `1`
>    - *"Number of header lines in each input line"*: `1`
>    - *"Fill character"*: `0.0`
>
> 2. **Rename** {% icon galaxy-pencil %} output `Sample A output table`
>
> 3. {% tool [Column join](toolshed.g2.bx.psu.edu/repos/iuc/collection_column_join/collection_column_join/0.0.3) %} with the following parameters:
>    - {% icon param-file %} *"Tabular files"*: `B proportions actual-infer`
>    - *"Identifier column"*: `1`
>    - *"Number of header lines in each input line"*: `1`
>    - *"Fill character"*: `0.0`
>
> 4. **Rename** {% icon galaxy-pencil %} output `Sample B output table`
>
{: .hands_on}

# Visualise outputs
