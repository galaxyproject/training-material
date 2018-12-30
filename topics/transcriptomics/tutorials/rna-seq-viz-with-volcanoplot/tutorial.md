---
layout: tutorial_hands_on
title: Visualization of RNA-Seq results with Volcano Plot
zenodo_link: "https://zenodo.org/record/2491305"
questions:
  - "How to generate a volcano plot from RNA-seq data?"
objectives:
  - "Create a volcano plot of RNA-seq data to visualize significant genes"
time_estimation: "1h"
key_points:
  - "A volcano plot can be used to quickly visualize significant genes in RNA-seq results"
contributors:
  - mblue9
---


# Introduction
{:.no_toc}

Volcano plots are commonly used to display the results of RNA-seq or other omics experiments. A volcano plot is a type of scatterplot that shows significance versus fold-change. It enables quick visual identification of genes with large fold changes that are also statistically significant. In a volcano plot, the most upregulated genes are towards the right, the most downregulated genes are towards the left, and the most statistically significant genes are towards the top.

To generate a volcano plot of RNA-seq results, we need a file of differentially expressed results which is provided for you here. To generate this file yourself, see the [RNA-seq counts to genes]({{ site.baseurl }}/topics/transcriptomics/tutorials/rna-seq-counts-to-genes/tutorial.html) tutorial. The file used here was generated from limma-voom but you could use a file from any RNA-seq differential expression tool, such as edgeR or DESeq2, as long as it has the required columns (see below).

The data for this tutorial comes from a Nature Cell Biology paper, [EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival](https://www.ncbi.nlm.nih.gov/pubmed/25730472)), Fu et al. 2015. This study examined the expression profiles of basal and luminal cells in the mammary gland of virgin, pregnant and lactating mice. Here we will visualise the results of the luminal pregnant vs lactate comparison.


> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Preparing the inputs

We will use one file for this analysis:

 * **Differentially expressed results file** (genes in rows, and 4 required columns: raw P values, adjusted P values (FDR), log fold change and gene labels)

## Import data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this RNA-seq exercise e.g. `RNA-seq volcano plot`
> 2. Import the differentially results table.
>
>     To import the file, there are two options:
>     - Option 1: From a shared data library if available (ask your instructor)
>     - Option 2: From [Zenodo](https://zenodo.org/record/2491305)
>
>         > ### {% icon tip %} Tip: Importing data via links
>         >
>         > * Copy the link location
>         > * Open the Galaxy Upload Manager
>         > * Select **Paste/Fetch Data**
>         > * Paste the link into the text field
>         > * Press **Start**
>         {: .tip}
>
>         - You can paste the link below into the **Paste/Fetch** box:
>
>           ```
>       https://zenodo.org/record/2491305/files/limma-voom_luminalpregnant-luminallactate
>           ```
>
>         - Select *"Genome"*: `mm10`
>
> 2. Rename the counts dataset as `DE results` using the {% icon galaxy-pencil %} (pencil) icon.
> 3. Check that the datatype is `tabular`.
>    If the datatype is not `tabular`, please change the file type to `tabular`.
>
>    > ### {% icon tip %} Tip: Changing the datatype
>    > * Click on the {% icon galaxy-pencil %} (pencil) icon displayed in your dataset in the history
>    > * Choose **Datatype** on the top
>    > * Select `tabular`
>    > * Press **Save**
>    {: .tip}
{: .hands_on}


## Create volcano plot

We will create a volcano plot highlighting genes with FDR < 0.01 and a log fold change of 0.58 (equivalent to a fold-change of 1.5), as these were the values used in the original paper for this dataset.

> ### {% icon hands_on %} Hands-on: Create a Volcano plot
>
> 1. **Volcano Plot** {% icon tool %} to create a volcano plot
>    - {% icon param-file %} *"Specify an input file"*: `DE results`
>    - *"FDR (adjusted P value)"*: `Column 8`
>    - *"P value (raw)"*: `Column 7`
>    - *"Log Fold Change"*: `Column 4`
>    - *"Labels"*: `Column 2`
>    - *"Significance threshold"*: `0.01`
>    - *"LogFC threshold to colour"*: `0.58`
{: .hands_on}

![Volcano plot highlighting significant genes](../../images/rna-seq-viz-with-volcanoplot/volcanoplot.png)

In the plot above the genes are coloured if they pass the thresholds for FDR and Log Fold Change., red if they are upregulated and blue if they are downregulated.

## Create volcano plot highlighting genes of interest

We can also label one or more genes of interest in a volcano plot. This enables us to visualize where these genes are in terms of significance and in comparison to the other genes. In the original paper using this dataset, there is a heatmap of 31 genes (Fig. 6b). These genes are a set of 30 cytokines/growth factor identified as differentially expressed, and the authors' main gene of interest, Mcl1. We will label these genes in the volcano plot.


```
GeneID
Mcl1
Hbegf
Tgfb2
Cxcl16
Csf1
Pdgfb
Edn1
Lif
Kitl
Bmp1
Pdgfa
Cmtm3
Cx3cl1
Ctgf
Wnt5a
Ptn
Spp1
Bmp3
Cmtm8
Gmfg
Cxcl2
Cxcl3
Il15
Egf
Cmtm7
Il34
Pdgfd
Nov
Cmtm6
Ccl28
Cxcl1
```

> ### {% icon hands_on %} Hands-on: Create a Volcano plot labelling genes of interest
> 1. Create a file of the gene symbols of interest
>    - Paste the information above (the 31 gene symbols and header) into the Galaxy Data Uploader Paste/Fetch box
>    - Set File Type to `tabular`
>    - Use the {% icon galaxy-pencil %} (pencil) icon to rename the file to `volcano genes`
> 2. **Volcano Plot** {% icon tool %} to create a volcano plot
>    - {% icon param-file %} *"Specify an input file"*: `DE results`
>    - *"FDR (adjusted P value)"*: `Column 8`
>    - *"P value (raw)"*: `Column 7`
>    - *"Log Fold Change"*: `Column 4`
>    - *"Labels"*: `Column 2`
>    - *"Significance threshold"*: `0.01`
>    - *"LogFC threshold to colour"*: `0.58`
>    - *"Points to label"*: `Input from file`
>        - *"File of labels"*: `volcano genes`
{: .hands_on}

![Volcano plot highlighting genes of interest](../../images/rna-seq-viz-with-volcanoplot/volcanoplot_custom_genes.png)

As in the previous plot, genes are coloured if they pass the thresholds for FDR and Log Fold Change. The genes of interest in the file we supplied are labelled, and also coloured red or blue if they pass the thresholds. Here all 31 labelled genes are significant (red or blue) except for two genes. One is the authors' gene of interest, Mcl1, and this result is expected, as they showed it's expression did change, but it was not significant at the transcription level. The other gene Gmfg, has an FDR very slightly outside the significance threshold we used of 0.01 (0.0105).

# Conclusion
{:.no_toc}

In this tutorial we have seen how a volcano plot can be generated from RNA-seq data and used to quickly visualize significant genes.