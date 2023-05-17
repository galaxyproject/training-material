---
layout: tutorial_hands_on

title: 'Filter, Plot, and Explore with Seurat in RStudio'
subtopic: single-cell-CS
priority: 6
zenodo_link: ''

questions:
- Is my single cell dataset a quality dataset?
- How do I pick thresholds and parameters in my analysis? What’s a “reasonable” number, and will the world collapse if I pick the wrong one?
- How do I generate and annotate cell clusters?

objectives:
- Interpret quality control plots to direct parameter decisions
- Repeat analysis from matrix to clustering to labelling clusters
- Identify decision-making points
- Appraise data outputs and decisions
- Explain why single cell analysis is an iterative process (i.e. the first plots you generate are not final, but rather you go back and re-analyse your data repeatedly)

time_estimation: 3H

key_points:
- Being able to switch between Galaxy and RStudio when analyzing datasets can be useful when looking to adjust default parameters within Seurat's functions and workflow. 
- Seurat in RStudio gives more flexibility and specificity of analyses, but Galaxy offers great reproducibility and ease of analysis.
- Beginning to learn the syntax and use of R will expand your 

requirements:
-
    type: "internal"
    topic_name: single-cell
    tutorials:
        - scrna-case_alevin-combine-datasets
        - scrna-case_basic-pipeline


tags:
- single-cell
- seurat
- rstudio

contributions:
  authorship:
    - camila-goclowski
  editing:
    - nomadscientist
  funding:
    - eosc-life

notebook:
  language: r
  snippet: topics/single-cell/tutorials/scrna-case_FilterPlotandExploreRStudio-WIP/preamble.md
---

# Introduction

You’ve previously done all the work to make a single cell matrix. Now it’s time to fully process our data using Seurat: remove low quality cells, reduce the many dimensions of data that make it difficult to work with, and ultimately try to define clusters and find some biological meaning and insights! There are many packages for analysing single cell data - Seurat (Satija et al. 2015), Scanpy (Wolf et al. 2018), Monocle (Trapnell et al. 2014), Scater (McCarthy et al. 2017), and many more. We’re working with Seurat in RStudio because it is well updated, broadly used and highly trusted within the field of bioinformatics.

> <comment-title></comment-title>
> This tutorial is significantly based on the [Seurat documentation](https://satijalab.com/seurat) as well as [Seurat's Guided Clustering Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) with implementations from [clustree]() and []()
{: .comment}

## Get Data onto Galaxy 
> Let's get our dataset loaded into Galaxy, first. 

> <hands-on-title>GetData</hands-on-title>
{% tool [EBI SCXA Data Retrieval](toolshed.g2.bx.psu.edu/repos/ebi-gxa/retrieve_scxa/retrieve_scxa/v0.0.2+galaxy2) %} with the following parameters: 
> - *"SC-Atlas experiment accession"*: 'E-MTAB-6945'
> - *"Choose the type of matrix to download"*: 'Raw filtered counts'
{: .hands_on}

## Open RStudio in Galaxy 
You should now see a matrix.mtx, genes.tsv, barcodes.tsv, and exp_design.tsv files in your Galaxy history. For the rest of the workflow, let's move onto RStudio and get coding!
> <hands-on-title>Open RStudio in Galaxy</hands-on-title>
{% tool [RStudio](interactive_tool_rstudio)}
{: .hands_on}

The interactive RStudio tool should begin to load now. Make your way over to your Active Interactive Tools page (User > Active Interactive Tools > RStudio)

## Setting the environment and files 
> <hands-on-title>Setting the environment</hands-on-title>
> First thing's first, we should load the packages we will need into our environment. In your console (likely in the lower left corner of your RStudio window) run the following line of code. 

```r
library(Matrix)
library(Seurat)
library(dplyr) 
```
{: .hands_on}

> Alright, the packages are called - now let's get our data files moved from the Galaxy history and into our RStudio enviornment so that we can create a Seurat object.

## Upload, view and modify the files
> Now that we have made it into RStudio and called the packages we'll use, let's begin loading the datasets we retrieved in Galaxy into RStudio. Galaxy expedites this process by providing us the gx_get() function, which will output the file path to locate each dataset within our history. 

> <hands-on-title>Get Data</hands-on-title>
>> So, for example, our matrix was the first to be saved in our history. As such, we will ask for the filepath for the first slot in our history with the following line:  

```{r}
gx_get(1) #matrix.mtx 
```

> Now we have the file path! We can use the Matrix package to read our counts matrix into our environment, using our file path to let it know where to find the matrix. 

```{r}
library(Matrix)
matrix.mtx<-readMM("/import/1") 
```
> Now we will do the same thing with the feature, barcode, and experimental design files. Don't try to skip ahead and fill in the position of the dataset, reading in the files will not work without having used the gx_get() function. 

```{r}
gx_get(2) #genes.tsv
genes.tsv<-read.delim("/import/2", header = FALSE) 

gx_get(3) #barcodes.tsv
barcodes.tsv<-read.delim("/import/3", header = FALSE) 

gx_get(4) #exp_design.tsv
exp_design.tsv<-read.delim("/import/4")
```
> The formatting of the experimental design dataset has cell barcodes as the first column of data as opposed to the row names. In order to Seurat to properly use this dataset, we will need to make the cell barcodes  the row names. This can be accomplished by doing the following: 

```{r}
rownames(exp_design.tsv)<-exp_design.tsv$Assay 
```
{: .hands_on}

> Now, in our RStudio environment, we should have all of the data sets necessary to create a Seurat Object: the matrix, a file with feature (gene) names, a file with cell barcodes, and an optional, but highly useful, experimental design file containing sample (cell-level) metadata. 

## Generating Seurat object
> Next we will add our dimension names to our matrix. In the end, this will provide us with a matrix whose rows are gene names, columns are cell barcodes, and values are expression of a given gene in a given cell. The first dimension name will be assigned to the genes (rows), and the second dimension name will be assigned to the cells (columns):

><hands-on-title>Add Dimension Names</hands-on-title>

```{r}
matrix.mtx@Dimnames[[1]]<-genes.tsv$V2
matrix.mtx@Dimnames[[2]]<-barcodes.tsv$V1
```
{: .hands_on}

> In a more typical Seurat pipeline, or on a local version of RStudio, this step would be replaced with Read10x. Read10x is Seurat's automated function to add in feature and barcode names. However, due to the nature of how Galaxy histories and RStudio interact, we'll use this manual method. 

> Now we will create a Seurat object using our newly labelled counts matrix! Make sure you have called the Seurat library, first, or RStudio will not recognize the function. 

><hands-on-title>Create Seurat Object</hands-on-title>

```{r}
library(Seurat)
srt<-CreateSeuratObject(counts = matrix.mtx)
```
{: .hands_on}

> You've created a Seurat object, congratulations!

# Adding Cell Level Metadata
> Now that we have an object, we can add in our metadata from our experimental design dataframe (table). This will be useful to us shortly as we begin to visualize our data! 

><comment-title> </comment-title>
>The code preceding the left pointing arrow will indicate where to put your metadata (the name of your new metadata column: object@metadata$newcolumnname), and the code following the arrow will denote where to find that metadata information (metadatatable$columnname) 
{: .comment}

><hands-on-title>Add Cell Metadata</hands-on-title>

```{r}
srt@meta.data$Sex<-exp_design.tsv$Sample.Characteristic.sex.
srt@meta.data$Organism<-exp_design.tsv$Sample.Characteristic.organism.
srt@meta.data$Strain<-exp_design.tsv$Sample.Characteristic.strain.
srt@meta.data$Developmental.Stage<-exp_design.tsv$Sample.Characteristic.developmental.stage.
srt@meta.data$Age<-exp_design.tsv$Sample.Characteristic.age.
srt@meta.data$Individual<-exp_design.tsv$Sample.Characteristic.individual.
srt@meta.data$Disease<-exp_design.tsv$Sample.Characteristic.disease.
srt@meta.data$Genotype<-exp_design.tsv$Sample.Characteristic.genotype.
srt@meta.data$Organism.Part<-exp_design.tsv$Sample.Characteristic.organism.part.
srt@meta.data$Cell.Type<-exp_design.tsv$Sample.Characteristic.cell.type.
srt@meta.data$Factor.Value.Genotype<-exp_design.tsv$Factor.Value.genotype.
```
{: .hands_on}

> Now that we have our almost fully annotated object, we will add one more metadata column: percent mitochondrial (perc.mt). This metadata column will denote what percentage of a cell's feature (gene) expression is mitochondrial. 

> <hands-on-title>Add Percent Mitochondrial</hands-on-title>

```{r}
srt <- PercentageFeatureSet(srt, 
                               pattern = "^mt-", 
                               col.name = "perc.mt")
```
{: .hands_on}

><comment-title> </comment-title>
> For the sake of this data set, and many others, the mitochondrial genes will all be marked with an "mt" as the prefix, so that is how we have asked Seurat's PercentageFeatureSet function to search for mitochondrial genes. With that being said, once you are analyzing your own data, it is highly recommended that you figure out how your data set has labelled mitochondrial genes to ensure that you are calculating the correct percentage--the mt prefix may not always include all mitochondrial genes depending on how your dataset was labelled. 
{: .comment}

## QC Plots
> Now that we have a complete Seurat object, we can begin the filtering process. There will  be a number of ‘cells’ that are actually just empty droplets or low-quality. 
>There will also be genes that may be sequencing artifacts or that appear with such low frequency that statistical tools will fail to accurately analyse them. 
>This background noise of both cells and genes not only makes it harder to distinguish real biological information from sequencing artifacts, but also makes it computationally difficult to analyse. 
>First on our agenda is to filter the matrix to give us cleaner data off which to extract meaningful insight from. It will also allow for faster analysis of your data.

> We want to filter our cells, but first we need to know what our data looks like. There are a number of subjective choices to make within scRNA-seq analysis, for instance we now need to make our best informed decisions about where to set our thresholds (more on that soon!). 
>We’re going to plot our data a few different ways. Different bioinformaticians might prefer to see the data in different ways, and here we are only generating a few of the plots you can use. Ultimately you need to go with what makes the most sense to you.

> So let's generate some QC plots. First off, let's check our dataset for batch effect:

><hands-on-title>Check for Batch Effect</hands-on-title>
```{r}
VlnPlot(srt,
        group.by = "Individual",
        features = "nCount_RNA", 
        log = TRUE)
```
{: .hands_on}

> This plot shows us the number of cells split by the individual (mouse) from which the cells came from. Now, depending on your experimental design, batch may be represented by something other than individual--like timepoint or even the wet lab researcher who isolated the cells. 
> Ideally, we would like to see a relatively even distribution of counts for each individual (or batch) but if there isn’t, fear not, we can regress this variable out in a later step.

><comment-title> </comment-title>
>In order to accurately assess potential batch effects, use the "group.by" argument to indicate the variable which differed across experiments.   
{: .comment}

> Now let's get an idea of how different variables, like the sex or genotype of the mice, might be represented across our dataset. 

><hands-on-title>Biologicsl Variables</hands-on-title>
> 1. Sex?
 ```{r}
 VlnPlot(srt, group.by = "Sex",features = "nCount_RNA",log = TRUE)
> ```
> 2. Genotype?
  ```{r}
  VlnPlot(srt,
        group.by = "Genotype",
        features = "nCount_RNA", 
        log = TRUE)
 ```
{: .hands_on}

## Finding Our Filtering Parameters
> Now that we have a better understanding of what our data looks like, we can begin identifying those spurious reads and low quality cells for removal. First we'll plot the percent mito (perc.mt) against the cell count (nCount_RNA) to get an idea of what threshold we should set for nCount:

><hands-on-title>Find nCount Thresholds</hands-on-title>
```{r}
plot(x = srt$nCount_RNA, 
     y = srt$perc.mt, 
     main = "UMI Counts x Percent Mito", 
     xlab = "UMI_count", 
     ylab="% mito")
```
> We are looking for cell counts with high mitochondrial percentages in their feature expression. 
> High mito expression typically indicates stressed out cells (typically due to the extraction, sorting, or sample prep protocols). These cells won't tell us much biologically, rather, they will contribute noise that we will aim to filter out of our data. With that being said, there is a level of metabolic activity that is expected but will be specific to your samples/tissue/organism--so it is worth looking into what that might look like when it comes time to analyze your own data.  

> We can also zoom in on the x-axis to get a better idea of what threshold to set by adjusting the xlim argument:
>```{r}
>plot(x = srt$nCount_RNA, 
>     y = srt$perc.mt, 
>     main = "UMI Counts x Percent Mito", 
>     xlab = "UMI_count", 
>     ylab="% mito", 
>     xlim = c(0,1750))
>```
> It looks like just before nCount_RNA = 1750, the perc.mito peaks above 2%--a conservative threshold that still encompasses the majority of other cells.  
{: .hands_on}

><hands-on-title>Find mito Thresholds</hands-on-title>
> Now we can take a closer look at the y-axis to decide on a mito threshold to set. Once more, we want to get rid of as few cells as possible while still removing those with unexpectedly high mito percentages. 

>```{r}
>plot(x = srt$nCount_RNA, 
>     y = srt$perc.mt, 
>     main = "UMI Counts x Percent Mito", 
>     xlab = "UMI_count", 
>     ylab="% mito", 
>     ylim = c(0,3))
>```
> We can see a clear trend wherein cells that have around 3% mito counts or higher also have far fewer total counts. These cells are low quality, will muddy our data, and are likely stressed or ruptured prior to encapsulation in a droplet.

> Take a look at what setting those thresholds will include and disclude from our dataset: 

```{r}
prop.table(table(srt@meta.data$nCount_RNA > 1750)) 
prop.table(table(srt@meta.data$perc.mt > 3))
```
> If we are happy with those thresholds for cells and percent mito, we can look at the the gene count threshold next. If not, repeat the preceding steps to hone in on a threshold more suited for your dataset.

> To do so, let's plot the gene counts (nFeature_RNA) against the percent mito (perc.mt):

```{r}
plot(x = srt$nFeature_RNA, 
     y = srt$perc.mt, 
     main = "Gene Counts x Percent Mito", 
     xlab= "gene_count", 
     ylab="% mito")
```

> Once again, let's zoom in on the x-axis but this time to get an idea of which nFeature_RNA threshold to set:

```{r}
plot(x = srt$nFeature_RNA, 
     y = srt$perc.mt, 
     main = "Gene Counts x Percent Mito", 
     xlab= "gene_count", 
     ylab="% mito", 
     xlim = c(0,1275))
```
> You can see how cells with nFeature_RNA up to around, perhaps 575 genes, often have high perc.mt. The same can be said for cells with nFeature_RNA above 1275. We could also use the violin plots to come up with these thresholds, and thus also take batch into account. It’s good to look at the violins as well, because you don’t want to accidentally cut out an entire sample (i.e. N703 and N707 which both have cell counts on the lower side).

> Now let's take a look at what those nFeature_RNA thresholds will include and disclude from our data. 

```{r}
prop.table(table(srt@meta.data$nFeature_RNA > 1275 | srt@meta.data$nFeature_RNA < 575))
```
## Applying our Thresholds 
> Once we are happy with our filtering thresholds, it’s now time to apply them to our data!

><hands-on-title>Filter</hands-on-title>

```{r}
subset_srt<-subset(srt,
                      nCount_RNA > 1750 & nFeature_RNA > 1275 & perc.mt < 3 | nFeature_RNA < 600) 
```
{: .hands_on}

><comment-title> </comment-title>
> In this step we are also creating a new object (notice the new object name preceding the subset function you just ran) so that we may compare back and forth between our unfiltered and filtered data set if we please. 
{: .comment}

> Next, we want to filter out genes that no longer show nay expression in the cells remaining in our filtered dataset. In order to do so we will extract that filtered matrix from our filtered object. 

><hands-on-title>Filter Genes: Extract filtered matrix</hands-on-title>

```{r}
subset_matrix<-GetAssayData(subset_srt)
```
{: .hands_on}

> Now that you’ve removed a whole heap of cells, and since the captured genes are sporadic (i.e. a small percentage of the overall transcriptome per cell) this means there are a number of genes in your matrix that are currently not in any of the remaining cells. Genes that do not appear in any cell, or even in only 1 or 2 cells, will make some analytical tools break and overall will not be biologically informative. So let’s remove them! 
>We can take the filtered matrix we just extracted and create a new Seurat object, this time including an additional min.cells argument. This will remove any genes from our matrix that have less than 3 cells expressing them. Note that 3 is not necessarily the best number, rather it is a fairly conservative threshold. You could go as high as 10 or more.

><hands-on-title>Filter Genes: Create New Object</hands-on-title>

```{r}
filtered_srt <- CreateSeuratObject(counts = subset_matrix, 
                                      meta.data = subset_srt@meta.data, 
                                      min.cells = 3)
```
{: .hands_on}

> Now that we have filtered out both noisy "cells" and genes from our dataset, let's clean up our environment. Remove objects that we no longer need to ensure that we stay organized and RStudio has enough memory capacity to perform downstream analyses. This likely will not be an issue while doing this tutorial, but in practice it will help things run smoothly. 

><hands-on-title>Remove Old Objects</hands-on-title>

```{r}
rm(subset_matrix, subset_srt)
```
{: .hands_on}

## Processing
> Currently, we still have quite big data. We have two issues here:
> 1. We already saw in our filtering plots that there are differences in how many transcripts and genes have been counted per cell. This technical variable can obscure biological differences. 
> 2. We like to plot things on x/y plots, so for instance Gapdh could be on one axis, and Actin can be on another, and you plot cells on that 2-dimensional axis based on how many of each transcript they possess. While that would be fine, adding in a 3rd dimension (or, indeed, in our case, many more dimensions), is a bit trickier. So, our next steps will be to transform our big data object into something that is easy to analyse and easy to visualize.

> We will run SCTransform, a combinatorial function by Seurat that normalizes the data, finds variable features, and then scales the data. In their initial workflow, and in the Scanpy version of this tutorial, these steps are run individually. However, with the second version of SCTransform comes time efficiency and optimization for downstream analyses. 

><hands-on-title>SCTransform</hands-on-title>

```{r}
filtered_srt<- SCTransform(filtered_srt, 
                          vars.to.regress = c("perc.mt", "nFeature_RNA", "nCount_RNA"),
                          verbose = TRUE, 
                          return.only.var.genes = FALSE, 
                          seed.use = 1448145)
```
{: .hands_on}

> Normalisation helps reduce the differences between gene and UMI counts by fitting total counts to 10,000 per cell. The inherent log-transform (by log(count+1)) aligns the gene expression level better with a normal distribution. This is fairly standard to prepare for any future dimensionality reductions.

> We also have loads of genes, but not all of them vary in expression from cell to cell. For instance, housekeeping genes are defined as not changing much from cell to cell, so we could remove these from our data to simplify our analyses. The find variable genes step within SCTransform (and Seurat's FindVariableFeatures function) will flag the genes that do vary across the cells for future analysis.

> Then, SCTransform (or Seurat's ScaleData function) will scale the data so that all genes have the same variance and a zero mean. 
>This is an important step to set up our data for further dimensionality reduction. It also helps negate sequencing depth differences between samples, since the gene levels across the cells become comparable. 

><comment-title> </comment-title>
>Note, that the differences from scaling etc. are not the values you have at the end - i.e. if your cell has average GAPDH levels, it will not appear as a ‘0’ when you calculate gene differences between clusters.
{: .comment}

> Although we've made our expression values comparable to one another and our overall dataset less computationally demanding, we still have too many dimensions (n cells x n genes!). 

> Transcript changes are not usually singular - which is to say, genes were in pathways and in groups. It would be easier to analyse our data if we could more easily group these changes.To address this we will run principal component analysis (PCA). 
>Principal components are calculated from highly dimensional data to find the most representative spread in the dataset. So in our highly variable gene dimensions, there will be one line (axis) that yields the most spread and variation across the cells. That will be our first principal component. We can calculate the first handful of principal components in our data to drastically reduce the number of dimensions:

><hands-on-title>Principal Component Analysis</hands-on-title>

```{r}
filtered_srt <- RunPCA(filtered_srt, 
              features = VariableFeatures(object = filtered_srt))
```
{: .hands_on}

> To visualize how our principal components (PCs) represent our data, let's create an elbow plot: 

><hands-on-title> Visualize Principal Components</hands-on-title>

```{r}
ElbowPlot(filtered_srt, ndims=50)
```
{: .hands_on}

> We can see that there is really not much variation explained past the 9th PC. So we might save ourselves a great deal of time and muddied data by focusing on the top 10 PCs. 
> You can also think about it like choosing a threshold of variance explained. Conservatively, 2.5 standard deviations are explained by about 10 of the PCs. 

> We’re still looking at around 10 dimensions at this point--likely not the easiest to visualize. To make our lives easier, we need to identify how similar a cell is to another cell, across every cell across each of these dimensions. 
>For this, we will use the k-nearest neighbor (kNN) graph, to identify which cells are close together and which are not. 
>The kNN graph plots connections between cells if their distance (when plotted in this 10 dimensional space!) is amongst the k-th smallest distances from that cell to other cells. This will be crucial for identifying clusters, and is necessary for plotting a UMAP--which is what will ultimately allow us to visualize our data in 2 dimensions. 
>From UMAP developers: “Larger neighbor values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50, with a choice of 10 to 15 being a sensible default”.

>Let's now use the 10 PC threshold we chose from the Elbowplot and apply it to find neighbors:

><hands-on-title>Find Neighbor Cells</hands-on-title>

```{r}
filtered_srt <- FindNeighbors(filtered_srt, dims = 1:10)
```
{: .hands_on}

> Now we can uses the neighborhood graph to identify clusters of cells whose transcriptional profiles look similar to one another.

><hands-on-title>Find Clusters</hands-on-title>

```{r}
filtered_srt <- FindClusters(filtered_srt, resolution = 0.5)
```
{: .hands_on}

>Unfortunately, identifying clusters is not as majestic as biologists often think - the math doesn’t necessarily identify true cell clusters. Every algorithm for identifying cell clusters falls short of a biologist knowing their data, knowing what cells should be there, and proving it in the lab. Sigh. 
> So, we’re going to make the best of it as a starting point and see what happens! We will define clusters from the kNN graph, based on how many connections cells have with one another. Roughly, this will depend on a resolution parameter for how granular you want to be.

> Now that we have made note within our object of which cells cluster together, we can start to visualize our data! Two major visualizations for this data are tSNE and UMAP. We must calculate the coordinates for both prior to visualization. For tSNE, the parameter perplexity can be changed to best represent the data, while for UMAP the main change would be to change the kNN graph above itself, by changing the neighbors.

><hands-on-title>UMAP</hands-on-title>

```{r}
filtered_srt <- RunUMAP(filtered_srt, dims = 1:10)
```
{: .hands_on}

## Take a Look
> Now that we have run dimensionality reduction on our dataset, it is ready for visualization. Let's take a look at what our cells look like in a UMAP projection: 

><hands-on-title>DimPlots</hands-on-title>

```{r}
DimPlot(filtered_srt, reduction = "umap", label = TRUE, label.box = TRUE)+ NoLegend()
```
{: .hands_on}

> We can also look for expression of particular genes and see how those map to our UMAP projection. This is often useful in getting a quick and initial understanding of which clusters might be representing which cell types.

><hands-on-title>FeaturePlots</hands-on-title>
```{r}
FeaturePlot(filtered_srt, features = "Gapdh")
```
{: .hands_on}

> We just plotted a housekeeping gene, Gapdh, so the broad expression is expected. 

>In practice, it is helpful to plot known markers of cell types you expect to be in your dataset. This will give you a first look at how your cells are clustered. 
>For example, we can plot macrophage marker Aif1 and get an idea of which cells and/or clusters might resemble macrophages: 

><hands-on-title>FeaturePlot</hands-on-title>
```{r}
FeaturePlot(filtered_srt, features = "Aif1")
```
{: .hands_on}

It is a good idea, when analyzing your own data, to plot some markers of cell types you expect to be present. Later on we can also use these FeaturePlots to visualize manual annotation of clusters. 

## Differential Expression Testing: Finding Markers
> Because each cluster of cells was grouped based on similar transcriptome profiles, each cluster will inherently differ from one another based on a set of "marker" genes. 
>Following an initial look at the DimPlots and FeaturePlots, we can take an even closer look at which genes are driving the clustering. 
>In order to do so we can run cluster level differential expression tests. First, we will need to set our object's active identity to be the clusters. This will ensure that when Seurat's differential expression function is run, the groupings of cells across which it will compare are the clusters. 
> Then, we'll run Seurat's FindAllMarkers function, which will compare each identity (in this case cluster) against every other identity within its class (all the other clusters). This function of marker finding is partoicularly useful in identifying up, or down, regulated genes that drive the differences in identity. 

><hands-on-title>Find Cluster Markers</hands-on-title>
```{r}
Idents(filtered_srt)<- filtered_srt$seurat_clusters
cluster_markers<-FindAllMarkers(object = filtered_srt)
```
{: .hands_on}

>We'll use these marker lists later on to label our cell types. 

> We can also see which genes are differentially expressed across other variables in our metadata. For example, you can see which genes are up or down regulated across the different genotypes present in our dataset. To do so, let's first get a list of all the identity classes in our data. This information is kept in the metadata column, and any categorical variable will do. Here, let's pick genotype. 

><hands-on-title>See your options and Change your object's identity</hands-on-title>
```{r}
metadata<-as.data.frame(filtered_srt@meta.data)
Idents(filtered_srt)<- filtered_srt$Genotype
```
{: .hands_on}

> The "metadata" object now in your environment is a dataframe with column names representing the different identities you may choose to group your cells by when running differential expression. The second line above sets our object's identity class to the genotype from which the cell came from. 

> Now, let's see what genes differentiate our wildtype from our mutant cells. First, we can identify how many different genotypes are in our data: 

><hands-on-title>List the Genotypes in your Data</hands-on-title>
```{r}
unique(filtered_srt$Genotype) 
```
{: .hands_on}

> This output helpfully shows us what the genotypes are, and how they are labelled in our metadata. The small details, like capitalization, are important for referencing metadata information--our references must perfectly match the labelling in the object, otherwise it will not be recognized by the functions. 

> Now that we know how our wildtype and mutant cells are labelled, we can use that information to directly compare the two. This time we will use a pairwise comparison method by using Seurat's FindMarkers function (not to be confused with FindAllMarkers which has a comprehensive comparison approach):

><hands-on-title>Genotype Markers</hands-on-title>
```{r}
markers<-FindMarkers(object = filtered_srt, ident.1 = "wild type genotype", ident.2 = "Igf2-p0 heterozygous knockout", test.use = "wilcox")
```
{: .hands_on}

> The above function will find all of the differentially expressed genes between ident.1 (wildtype) and ident.2 (mutant) using the Wilcoxon test. The resulting output will show genes with positive fold changes (denoting a higher expression in the first identity--wildtype) and negative fold changes (denoting a higher expression value in the second identity--mutant). 

> This same test of differential expression can be run using any identity class and any two identities within the same class. As this is a more fine tuned comparison than FindAllMarkers, it can be useful to uncover differences across specific samples. 

## Biological Interpretations
>Now it’s the fun bit! We can see where genes are expressed, and start considering and interpreting the biology of it. At this point, it’s really about what information you want to get from your data - the following is only the tip of the iceberg. However, a brief exploration is good, because it may help give you ideas going forward with for your own data. Let us start interrogating our data!

Let's take a look at what our clusters look like: 
><hands-on-title> Clusters</hands-on-title>
```{r}
DimPlot(object = filtered_srt, reduction = "umap", group.by = "seurat_clusters")
```
{: .hands_on}

> Note that the cluster numbering is based on size alone - clusters 0 and 1 are not necessarily related, they are just the clusters containing the most cells. 
>It would be nice to know what these cells are. This analysis (googling all of the marker genes, both checking where the ones you know are as well as going through the marker tables you generated!) is a fun task for any individual experiment, so we’re going to speed past that and nab the assessment from the original paper!

| Clusters | Markers                 | Cell Type                           |
|----------|-------------------------|-------------------------------------|
| 4        | Il2ra                   | Double negative (early T-cell)      |
| 0,1,2,6  | Cd8b1, Cd8a, Cd4        | Double positive (middle T-cell)     |
| 5        | Cd8b1, Cd8a, Cd4 - high | Double positive (late middle T-cell)|
| 3        | Itm2a                   | Mature T-cell                       |
| 7        | Aif1                    | Macrophages                         |

>