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


#ONCE WE'RE IN R:
look at pre analysis RNAseq matrix
Note that this dataset also contains ~5% of mouse cells, which we can use as negative controls for the protein measurements. For this reason, the gene expression matrix has HUMAN_ or MOUSE_ appended to the beginning of each gene
```{r}
gx_get(1)
RNA<-read.csv('/import/1')
```

look at pre analysis ADT matrix 
```{r}
gx_get(2)
ADT<-read.csv('/import/2')
```

run log from the Seurat tool 
would be nice to be able to look at in R (html format though) 
can also just look at this on Galaxy 
```{r}
gx_get(3)
```

markers from ADT Assay 
AKA protein markers

```{r}
gx_get(4)
protein_markers<-read.table('/import/4', header = T)
```

let's subset this huge list down to just those that are statistically relevant
```{r}
protein_markers<-subset(protein_markers, p_val_adj < 0.045)
```

markers from RNA assay 
aka rna markers 
we'll also subset these down to only include statistically significant ones
```{r}
gx_get(5)
rna_markers<-read.table('/import/5', header = T)
rna_markers<-subset(rna_markers, p_val_adj < 0.045)
```

Cite-seq graphs (pdf) 
```{r}
gx_get(6)
```

seurat object post marker analysis 
this is the most processed object coming out of the Seurat tool
```{r}
library(Seurat)
library(SeuratObject)
gx_get(7)
srt<-readRDS('/import/7')
```

another marker list (I presume combined protein + rna?)
```{r}
gx_get(8)
markers<-read.table('/import/8', header = T)
markers<-subset(markers, p_val_adj < 0.045)
```

normalize ADT assay
would be cool to add to the Seurat tool 
```{r}
srt <- NormalizeData(srt, normalization.method = "CLR", margin = 2, assay = "ADT")
```

visualize CD19 protein and rna side by side 
```{r}
library(ggplot2)
DefaultAssay(srt)<-"ADT"
adt_cd19<-FeaturePlot(srt, features = "CD19") + ggtitle("CD19 Protein")
DefaultAssay(srt)<-"RNA"
rna_cd19<-FeaturePlot(srt, features = "CD19") + ggtitle("CD19 RNA")
adt_cd19|rna_cd19
```

As always, there are a couple of ways we can get the same output:
instead of swapping the default assay back and forth, we can use specific assay keys to do this
let's find out waht the RNA key is: 
```{r}
Key(srt[["RNA"]])
```
what about the ADT 
```{r}
Key(srt[["ADT"]])
```

now that we know what the key is, we can just include that in the featurename, overriding the default assay of the object
```{r}
cd19_adt<-FeaturePlot(srt, features = "adt_CD19") + ggtitle("CD19 Protein")
cd19_rna<-FeaturePlot(srt, features = "rna_CD19") + ggtitle("CD19 RNA")
cd19_adt|cd19_rna
```
We can also make ADT scatter plots (functioning similarly to biaxial plots for FACS). 
```{r}
FeatureScatter(srt, feature1 = "adt_CD19", feature2 = "adt_CD3")
```
we can also see the relationship between protein and rna 
```{r}
FeatureScatter(srt, feature1 = "adt_CD3", feature2 = "rna_CD3E")

FeatureScatter(srt, feature1 = "adt_CD4", feature2 = "adt_CD8")
```

MOVING FORWARD: 
Seurat v4 also includes additional functionality for the analysis, visualization, and integration of multimodal datasets. For more information, please explore the resources below:
- Defining cellular identity from multimodal data using WNN analysis in Seurat v4 vignette
- Mapping scRNA-seq data onto CITE-seq references [vignette]
- Introduction to the analysis of spatial transcriptomics analysis [vignette]
- Analysis of 10x multiome (paired scRNA-seq + ATAC) using WNN analysis [vignette]
- Signac: Analysis, interpretation, and exploration of single-cell chromatin datasets [package]
- Mixscape: an analytical toolkit for pooled single-cell genetic screens [vignette]