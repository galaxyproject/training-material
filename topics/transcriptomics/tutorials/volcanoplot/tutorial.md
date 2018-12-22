---
layout: tutorial_hands_on
title: Visualization of RNA-Seq results with Volcano Plot
zenodo_link: "https://zenodo.org/record/2491305"
questions:
  - "How to generate a volcano plot from RNA-seq data?"
objectives:
  - "Create a volcano plot of RNA-seq data"
time_estimation: "1h"
key_points:
  - "A volcano plot can be used to quickly visualize significant genes" 
contributors:
  - mblue9
---


# Introduction
{:.no_toc}

Volcano plots are commonly used to display the results of RNA-seq or other omics experiments. A volcano plot is a type of scatterplot that plots significance versus fold-change. A volcano plot enables quick visual identification of the genes that show large fold changes that are also statistically significant. The most significant genes are towards the top of the plot, the most upregulated genes are towards the right, nad the most downregulated genes towards the left.

To generate a volcano plot of RNA-seq results, we need a file of differetially expressed results, which is provided for you here. To generate this file yourself, see the [RNA-seq counts to genes]({{ site.baseurl }}/topics/transcriptomics/tutorials/limma-voom/tutorial.html) tutorial. You could also use a file generated from other RNA-seq differential expression tools such as edgeR or DESeq2.

The data for this tutorial comes from a Nature Cell Biology paper, [EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival](https://www.ncbi.nlm.nih.gov/pubmed/25730472)), Fu et al. 2015. This study examined the expression profiles of basal and luminal cells in the mammary gland of virgin, pregnant and lactating mice. Six groups are present, with one for each combination of cell type and mouse status. Here we will visualise the results of the luminal pregnant vs lactate comparison.


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

 * **Differentially expressed results file** (with genes in rows, and columns for raw P values, adjusted P values (FDR) and log fold change)

## Import data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this RNA-seq exercise e.g. `RNA-seq volcano plot`
> 2. Import the differentially results table.
>
>     To import the file, there are two options:
>     - Option 1: From a shared data library if available (ask your instructor)
>     - Option 2: From [Zenodo](add link)
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
>    - *"Points to label"*: `Significant`
>        - *"Only label top most significant"*: `10`
{: .hands_on}


## Create volcano plot of custom genes

We can also label genes of interest in a volcano plot, to visualize where they lie in relation to other genes. In the original paper using this dataset, there is a heatmap (Fig. 6b below) of 31 genes. These genes are the authors' main gene of interest in the paper, Mcl1, and a set of cytokines/growth factors, identified as differentially expressed in the luminal lactate comparison. We will label these genes in the volcano plot.


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

## Create volcano plot highlighting genes of interest

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


# Conclusion
{:.no_toc}

In this tutorial we have seen how a volcano plot can be generated from RNA-seq data and used to quickly visualize significant genes. This uses the same dataset from the tutorials, [RNA-seq reads to counts]({{ site.baseurl }}/topics/transcriptomics/tutorials/limma-voom_fastqs_to_counts/tutorial.html), [RNA-seq counts to genes]({{ site.baseurl }}/topics/transcriptomics/tutorials/limma-voom/tutorial.html), and [RNA-seq genes to pathways]({{ site.baseurl }}/topics/transcriptomics/tutorials/limma-voom_gene_set_testing/tutorial.html).