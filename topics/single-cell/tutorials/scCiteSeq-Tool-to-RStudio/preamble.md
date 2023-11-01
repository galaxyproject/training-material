# Introduction
Multiomic analyses are a new and exciting way to understand the world of biology through bioinformatics! Cite-Seq (% cite Satija&Smibert2017 %) is merely one multiomic technology which enables us to measure both transcriptomes and cell surface proteins simultaneously, from the same cell! Transcriptomic measurements are achieved via RNA sequencing techniques and the surface protein abundanc emeasurements are quantified via DNA barcoded antibodies.

Seurat has kept up to date with the capacities of multimodal technologies such as Cite-Seq, which means once you've familiarized yourself with Seurat (check out [Filter, Plot, and Explore with Seurat] ({% link topics/single-cell/tutorials/scrna-case_FilterPlotandExploreRStudio/tutorial.md %}) to start doing so in RStudio with an scRNA-seq dataset!), you can seamlessly continue to use the package to analyze and explore many other types multimodal single-cell datasets.

> <comment-title></comment-title>
> This tutorial is significantly based on the Seurat documentation({% cite Satija2015 %}) as well as [Seurat's Guided Clustering Tutorial](../scrna-case_FilterPlotandExploreRStudio/tutorial.bib).
{: .comment}

# Get Your Data
For this tutorial, we'll use a publicly available dataset of 8,617 cord blood mononuclear cells (CBMCs) which have been sequenced for transcriptomic measurements as well as for 11 surface proteins. 

First on the to-do list is importing our dataset. You can do this in a number of ways: 

1. but to start let's use the Upload Data button on the upper left of your screen, above the : 

<hands-on-title>GetData</hands-on-title>
>
> Run{% tool [EBI SCXA Data Retrieval](toolshed.g2.bx.psu.edu/repos/ebi-gxa/retrieve_scxa/retrieve_scxa/v0.0.2+galaxy2) %} with the following parameters:
> - *"SC-Atlas experiment accession"*: `E-MTAB-6945`
> - *"Choose the type of matrix to download"*: `Raw filtered counts`
{: .hands_on}
