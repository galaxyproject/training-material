# Introduction
> <comment-title></comment-title>
> This tutorial is significantly based on the Seurat documentation({% cite Satija2015 %}) as well as Seurat's vignette entitled [Using Seurat with multimodal data](https://satijalab.org/seurat/articles/multimodal_vignette).
{: .comment}
 
Multiomic analyses are a new and exciting way to understand the world of biology through bioinformatics! Cite-Seq (% cite Satija&Smibert2017 %) is one of many multiomic technologies. Cite-Seq enables us to measure both transcriptomes and cell surface proteins from the same cell. Transcriptomic measurements are achieved via RNA sequencing techniques and the surface protein abundance measurements are quantified via DNA barcoded antibodies. As of current, it is only possible to tag a small number of surface proteins--typically around 10-15. 

Seurat has kept up to date with the capacities of multimodal technologies such as Cite-Seq, which means once you've familiarized yourself with Seurat, you can seamlessly continue to use the package to analyze and explore many other types multimodal single-cell datasets. 

<comment-title></comment-title>
Check out [Filter, Plot, and Explore with Seurat] ({% link topics/single-cell/tutorials/scrna-case_FilterPlotandExploreRStudio/tutorial.md %}) to start doing so in RStudio with an scRNA-seq dataset!
{: .comment}

Before we can start exploring, you'll process our transcriptomic and surface protein measdurements into a Seurat object. Galaxy has kindly provided us a Seurat tool with Cite-Seq functionality that can take our raw csv files and output Seurat objects which are easy to explore! 

<comment-title></comment-title>
If you're interested in what the tool is doing behind the scenes, check out Seurat's [Using Seurat with multimodal data](https://satijalab.org/seurat/articles/multimodal_vignette) vignette. The first portion of the tutorial is what the
{: .comment}

# Get Your Data
For this tutorial, we'll use a publicly available dataset of 8,617 cord blood mononuclear cells (CBMCs) which have been sequenced for transcriptomic measurements as well as for 11 surface proteins. 

First on the to-do list is importing our csv files. You can do this in a couple of ways: 

1.  Use the Upload Data button on the upper left of your screen:
    [plot1]()

    Select the 
    [plot2]()
    
    Copy the following links into the box
    ADT: ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz
    RNA : ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz
    
    Select "Start" and then close once the files have been 100% imported. 


2. Import [my history]() 
    {% snippet faqs/galaxy/histories_import.md %}

3. Import from Zenodo


Next we'll run those csv files through the Cite-Seq tool with the following parameters:
> Run{% tool [Seurat](toolshed.g2.bx.psu.edu/repos/iuc/seurat/seurat/4.3.0.1+galaxy1) %} with the following parameters:
> - *"Which Seuray method should be run"*: `Cite-seq`
> - *"RNA counts file"*: `1: GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz`
> - *"Protein counts file"*: `2: GSE100866_CBMC_8K_10X-ADT_umi.csv.gz`
> - *"Minimum cells"*: `5`
> - *"Minimum genes"*: `10`
> - *"Low threshold for filtering cells"*: `1`
> - *"High threshold for filtering cells"*: ``
> - *"Include violin plot and scatter plot of cell features"*: `Yes`
> - *"Output seurat object after data normalization"*: `No`
> - *"Include plot of variable features"*: `Yes`
> - *"Output seurat object after data scaling"*: `No`
> - *"Number of PCs to use in plots"*: `15`
> - *"Include PCA plots"*: `Yes`
> - *"Output seurat object after PCA analysis"*: `No`
> - *"Perplexity parameter"*: ``
> - *"Resolution parameter"*: `0.8`
> - *"Include UMAP and TSNE plots"*: `Yes`
> - *"Output seurat object after TSNE and UMAP analysis"*: `No`
> - *"Minimum percent cells"*: `0.1`
> - *"Log fold change threshold"*: `0.25`
> - *"Include heatmaps of markers"*: `Yes`
> - *"Output marker data"*: `Yes`
> - *"Output list of cite-seq markers"*: `Yes`
> - *"Compare specific feature's effect on protein and rna expression?"*: `No`
> - *"Compare top RNA and protein features graphicaly against themselves and one another"*: `Yes`
> - *"How many of the top featyre should be shown"*: `5`
{: .hands_on}

