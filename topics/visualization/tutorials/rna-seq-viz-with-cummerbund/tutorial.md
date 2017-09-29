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

# Importing RNA-Seq data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history
> 2. Import the SQLite database [`cuffdiff-sqlite`](https:zenodo.org/record/XXX)
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
>    > Rename the dataset to "RNA-Seq SQLite data"
>    {: .comment}
> By default, when data is imported via its link, Galaxy names it with its URL.
{: .hands_on}

# Reasons for visualizations

To make sense of the available RNA-Seq data, and overview the condition-specific gene expression levels of the provided samples, we need to visualize the imported dataset.

[CummeRbund](http://compbio.mit.edu/cummeRbund/) is an open-source tool that simplifies the analysis of a Cuffdiff RNA-Seq output. In particular, it helps researchers
- managing, integrating, and visualizing the data produced by Cuffdiff
- simplifying data exploration
- providing a bird's-eye view of the expresion analysis by describing relationships betweeen genes, transcripts, transcription start sites, and CDS regions
- exploring subfeatures of individual genes or gene-sets
- creating publication-ready plots

> ### {% icon hands_on %} Hands-on: Visualization
>
> 1. **CummeRbund** {% icon tool %}: Run CummeRbund on the imported SQLite dataset with default parameters
> 2. Inspect CummerBund's report on its webpage output
>
>    > ### {% icon tip %} Tip: Inspecting the content of a file in Galaxy
>    >
>    > * Click on the eye ("View data") on the right of the file name in the history
>    > * Inspect the content of the file on the middle
>    {: .tip}
>
>    > ### {% icon question %} Questions
>    >
>    > 1. How XXX?
>    > 2. How YYY?
>    > 3. How ZZZ?
>    >
>    > <details>
>    > <summary>Click to view the answers</summary>
>    > <ol type="1">
>    > <li></li>
>    > <li></li>
>    > <li></li>
>    > </ol>
>    > </details>
>    {: .question}
{: .hands_on}

# CASE 1

Short introduction about this subpart.

> ### {% icon comment %} Comment
>
> Do you want to learn more about the principles behind mapping? Follow our [training](../../NGS-mapping)
> {: .comment}

# CASE 2

# CASE 3

# Conclusion
{:.no_toc}

Visualization tools help researchers making sense of data, providing a bird's-eye view of the underlying analysis.
In this tutorial we overviewed the advantages of visualizing RNA-Seq results with CummeRbund, and gained insights on Cuffdiff's big-data output by plotting informations relative to CASE1, CASE2, CASE3.

