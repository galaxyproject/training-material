---
layout: tutorial_hands_on
topic_name: visualization
tutorial_name: rna-seq-viz-with-cummerbund
---

# Introduction
{:.no_toc}

RNA-Seq analysis helps researchers finding new genes and spliced variants, and provides a quantification of cell- and context- specific gene expression. Its data is however complex to handle, and requires both computer science and mathematical knowledge to be managed and interpreted.

Visualization techniques are key to overcome the complexity of RNA-Seq data, and represent valuable tools to gather information and insights.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Reasons for visualizing RNA-Seq results

To make sense of the available RNA-Seq data, and overview the condition-specific gene expression levels of the provided samples, we need to visualize our results. Here we will use [CummeRbund](http://compbio.mit.edu/cummeRbund/).

CummeRbund is an open-source tool that simplifies the analysis of a CuffDiff RNA-Seq output. In particular, it helps researchers:
- managing, integrating, and visualizing the data produced by CuffDiff
- simplifying data exploration
- providing a bird's-eye view of the expression analysis by describing relationships betweeen genes, transcripts, transcription start sites, and CDS regions
- exploring subfeatures of individual genes or gene-sets
- creating publication-ready plots

A typical workflow for the visualization of RNA-Seq data involving CummeRbund:

![workflow](../../images/cummerbund-rna-seq-viz-with-cummerbund.png)

CummeRbund reads your RNA-Seq results from a [SQLite](https://www.sqlite.org/) database. This database has to be created using CuffDiff's SQLite output option.

> ### {% icon tip %} Tip: SQLite output with CuffDiff
>
> Instruct CuffDiff to organize its output in a SQLite database for later read with CummeRbund
>
> ![SQLite output](../../images/cummerbund-cuffdiff-set-sqlite.png)
{: .tip}

# Importing RNA-Seq result data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history
> 2. Import the [`CuffDiff SQLite`](https://zenodo.org/record/1001880/files/CuffDiff_SQLite_database.sqlite) dataset
>
>    > ### {% icon tip %} Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * Press **Start**
>    {: .tip}
>
>    > ### {% icon comment %} Comments
>    >
>    > Rename the dataset to "RNA-Seq SQLite result data"
>    {: .comment}
> By default, when data is imported via its link, Galaxy names it with its URL.
{: .hands_on}

CuffDiff's output data is organized in a SQLite database, so we need to extract it to be able to see how it looks like.

For this tutorial, we are interested in CuffDiff's tested transcripts for differential expression.

> ### {% icon hands_on %} Hands-on: Extract CuffDiff results
>
> 1. **Extract CuffDiff** {% icon tool %}: Lookup for this tool in the search bar, and select it
> 2. Click on "Select tables to output", and select only the table called "Transcript differential expression testing"
> 3. Click execute. Extract CuffDiff will extract the selected table
> 4. Inspect the table
>
>    > ### {% icon tip %} Tip: Inspecting the content of a file in Galaxy
>    >
>    > * Click on the eye ("View data") on the right of the file name in the history
>    > * Inspect the content of the file on the middle
>    {: .tip}
> 5. Each entry represents a differentially expressed gene, but not all are significant. We want to keep only those that are reported as significant differentially expressed. 
>
>    > ### {% icon question %} Questions
>    >
>    > 1. How to retain only the significant differentially expressed genes?
>    > 2. Which column stores this information?
>    >
>    > <details>
>    > <summary>Click to view the answers</summary>
>    > <ol type="1">
>    > <li>We need to filter on the column storing the record's significance</li>
>    > <li>Column 14</li>
>    > </ol>
>    > </details>
>    {: .question}
{: .hands_on}

# Filtering and sorting

We now want to first highlight the most significant differentially expressed genes in our analysis, and then obtain informative visualizations.

> ### {% icon hands_on %} Hands-on: Extract CuffDiff's most significant differentially expressed genes
>
> 1. **Filter** {% icon tool %}: Look up for this tool in the search bar, and select it
> 2. Select the extracted table as its input, and filter over the target column
>
>    > ### {% icon question %} Questions
>    > 1. What column stores the information of significance for each record?
>    > 2. Which conditional expression has to be set to filter all records on the selected column?
>    > 3. What happened to the records in the original table?
>    >
>    > <details>
>    > <summary>Click to view the answers</summary>
>    > <ol type="1">
>    > <li>column 14</li>
>    > <li>c14=='yes'</li>
>    > <li>All records whose "significant" field was set to "yes" have been retained, while the others filtered out</li>
>    > </ol>
>    > </details>
>    {: .question}
> 3. **Sort** {% icon tool %}: Look up for this tool in the search bar, and select it
> 4. Look at your data. The differential expression values are stored on column 10, we will sort (descending) all records on the basis of their value at the 10th column
> 4. Select the filtered table as the input, and provide the column on which the records have to be sorted, the sorting flavor, and the order
>
>    > ### {% icon question %} Questions
>    > 1. Since the start of our filtering process, how many records now represent the significant subset for extracting informations?
>    > 2. What does this shrinking of the number of lines represent?
>    >
>    > <details>
>    > <summary>Click to view the answers</summary>
>    > <ol type="1">
>    > <li>Click on the boxes in your history, their small preview higlights the number of lines: from ~140,000 to 219</li>
>    > <li>This process represents a necessary step to gather insights on the biological meaning of our samples in our analyses: putting the original raw RNA-Seq result data into context, cutting down the less-meaningful records to focus on what is needed to go from data to information</li>
>    > </ol>
>    > </details>
>    {: .question}
{: .hands_on}

# CummeRbund

With CummeRbund we can visualize our RNA-Seq results of interest.

CummeRbund generates always two outputs:
- the plot
- the R script responsible for generating the plot

We are interested in visualizing all expression values of all transcripts relative to the most significant differentially expressed gene we found in the previous section.

> ### {% icon hands_on %} Hands-on: Visualization
>
> 1. **CummeRbund** {% icon tool %}: Lookup for this tool in the search bar, and select it
> 2. Click on "Insert plot"
> 3. Select 800 x 600 as the width and height of the resulting image
> 4. Select "Expression Plot" as the plot type
> 5. Select "Isoforms" as the expression levels to plot
> 6. Select the gene "NDUFV1" as the gene to plot
> 7. Execute
{: .hands_on}

![Expression plot](../../images/cummerbund-expression-plot.png)

The Expression Plot represents the expression of all isoforms of a single gene (NDUFV1) with replicate FPKMs exposed.

Our plot has a modest number of isoforms, and is therefore already readable. However, in case of 5 or 6 isoforms, the plot can look very busy. We can therefore change the visualization type by selecting another type of plot.

> ### {% icon hands_on %} Hands-on: Visualization
>
> 1. **CummeRbund** {% icon tool %}: Lookup for this tool in the search bar, and select it
> 2. Click on "Insert plot"
> 3. Select 800 x 600 as the width and height of the resulting image
> 4. Select "Expression Bar Plot" as the plot type
> 5. Select "Isoforms" as the expression levels to plot
> 6. Select the gene "NDUFV1" as the gene to plot
> 7. Execute
{: .hands_on}

![Expression bar plot](../../images/cummerbund-expression-bar-plot.png)

Expression Bar Plot of a single gene (NDUFV1) with replicate FPKMs exposed.

> ### {% icon comment %} Comment
>
> These plots are shown also in [this](https://vimeo.com/channels/884356/128265982) Galaxy video tutorial.
>
> Would you like to obtain more sophisticated visualization of your RNA-Seq analysis results? Select different CummeRbund plot options, and look at their parametrizations according to the filtering and sorting operations we performed
> {: .comment}

# Conclusion
{:.no_toc}

Visualization tools help researchers making sense of data, providing a bird's-eye view of the underlying analysis results.
In this tutorial we overviewed the advantages of visualizing RNA-Seq results with CummeRbund, and gained insights on CuffDiff's big-data output by plotting informations relative to the most significant differentially expressed genes in our RNA-Seq analysis.

