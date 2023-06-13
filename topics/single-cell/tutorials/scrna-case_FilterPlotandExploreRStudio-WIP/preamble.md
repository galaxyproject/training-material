# Introduction

You’ve previously done all the work to make a single cell matrix. Now it’s time to fully process our data using Seurat: remove low quality cells, reduce the many dimensions of data that make it difficult to work with, and ultimately try to define clusters and find some biological meaning and insights! There are many packages for analysing single cell data - Seurat (Satija et al. 2015), Scanpy (Wolf et al. 2018), Monocle (Trapnell et al. 2014), Scater (McCarthy et al. 2017), and many more. We’re working with Seurat in RStudio because it is well updated, broadly used and highly trusted within the field of bioinformatics.

> <comment-title></comment-title>
> This tutorial is significantly based on the [Seurat documentation](https://satijalab.com/seurat) as well as [Seurat's Guided Clustering Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html).
{: .comment}

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get Data onto Galaxy 
> Let's get our dataset loaded into Galaxy, first. 

> <hands-on-title>GetData</hands-on-title>
{% tool [EBI SCXA Data Retrieval](toolshed.g2.bx.psu.edu/repos/ebi-gxa/retrieve_scxa/retrieve_scxa/v0.0.2+galaxy2) %} with the following parameters: 
> - *"SC-Atlas experiment accession"*: 'E-MTAB-6945'
> - *"Choose the type of matrix to download"*: 'Raw filtered counts'
{: .hands_on}

# Open RStudio in Galaxy 
You should now see a matrix.mtx, genes.tsv, barcodes.tsv, and exp_design.tsv files in your Galaxy history. For the rest of the workflow, let's move onto RStudio and get coding!
> <hands-on-title>Open RStudio in Galaxy</hands-on-title>
{% tool [RStudio](interactive_tool_rstudio)}
{: .hands_on}

The interactive RStudio tool should begin to load now. Make your way over to your Active Interactive Tools page (User (in the top bar)> Active Interactive Tools > RStudio)