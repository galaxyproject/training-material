---
layout: tutorial_hands_on

title: 'Filter, plot, and explore single cell RNA-seq data with Seurat'
subtopic: single-cell-CS
priority: 3
zenodo_link: 'https://zenodo.org/record/7053673'

input_histories:
  - label: "UseGalaxy.eu"
    history: https://singlecell.usegalaxy.eu/u/camila-goclowski/h/tool-based-seurat-fpe-input-data


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

requirements:
-
    type: "internal"
    topic_name: single-cell
    tutorials:
        - scrna-intro
        - scrna-preprocessing-tenx

tags:
- 10x
- paper-replication
- MIGHTS

contributions:
  authorship:
    - Camila-goclowski
    - pcm32
  editing:
    - pavanvidem
    - hexylena
---


You’ve previously done all the work to make a single cell matrix. Now it’s time to fully process our data using Seurat: remove low quality cells, reduce the many dimensions of data that make it difficult to work with, and ultimately try to define clusters and find some biological meaning and insights! There are many packages for analysing single cell data - Seurat ({% cite Satija2015 %}), Scanpy ({% cite Wolf2018 %}), Monocle ({% cite Trapnell2014 %}), Scater ({% cite McCarthy2017 %}), and many more. We’re working with Seurat because it is well updated, broadly used, and highly trusted within the field of bioinformatics.

> <comment-title></comment-title>
> This tutorial is significantly based on the Seurat documentation({% cite Satija2015 %}) as well as [Seurat's Guided Clustering Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html).
{: .comment}

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

We’ll provided you with experimental data to analyse from a mouse dataset of fetal growth restriction ({% cite Bacon2018 %}). This is the full dataset generated from [this tutorial]({% link topics/single-cell/tutorials/scrna-case_alevin/tutorial.md %}).

# Get Data onto Galaxy and generate a Seurat object

## Get Data onto Galaxy
To start, let's get our dataset loaded into Galaxy.

{% include _includes/cyoa-choices.html option1="EBI Data Retrieval" option2="Importing from a history" option3="Uploading from Zenodo" default="EBI-Data-Retrieval" text="There are multiple ways in which to collect the data for this tutorial. I find it easiest to do so via the EBI Data Retrieval." %}

<div class='EBI-Data-Retrieval' markdown='1'>
> <hands-on-title>EBI Data Retrieval</hands-on-title>
>
> You may retrieve that files necessary to construct a Seurat Object in this way. Doing to will alleviate the necessity to convert AnnData (Python) objects into Seurat (R) objects:
>
> Run{% tool [EBI SCXA Data Retrieval](toolshed.g2.bx.psu.edu/repos/ebi-gxa/retrieve_scxa/retrieve_scxa/v0.0.2+galaxy2) %} with the following parameters:
> - *"SC-Atlas experiment accession"*: `E-MTAB-6945`
> - *"Choose the type of matrix to download"*: `Raw filtered counts`
{: .hands_on}
</div>

<div class='Importing-from-a-history' markdown='1'>
> <hands-on-title>Importing from a history</hands-on-title>
> You can import [this history](https://singlecell.usegalaxy.eu/u/camila-goclowski/h/tool-based-seurat-fpe-input-data)
>
> {% snippet faqs/galaxy/histories_import.md %}
>
> This also alleviates the necessity to convert the AnnData object into a Seurat one, which is an additional step you must complete if you choose to use the next method.
{:.hands_on}
</div>

<div class='EBI-Data-Retrieval Importing-from-a-history' markdown='1'>

## Generating a Seurat object
You now should have imported the `matrix.mtx`, `genes.tsv`, `barcodes.tsv`, and `exp_design.tsv` files into your Galaxy history. In order for Seurat tools to work, we will have to convert the data into a format that Seurat recognizes. To do so, we will add row and column names to our matrix. In the end, this will leave us with a matrix whose rows are gene names, columns are cell barcodes, and each value in the matrix represent the expression value of a given gene in a given cell.

This can be accomplished via the Read10x step. **Read10x** tool implements Seurat's function to create a matrix and add in feature and barcode names simultaneously:

> <hands-on-title>Read10X</hands-on-title>
>
> Run{% tool [Seurat Read10x](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_read10x/seurat_read10x/4.0.4+galaxy0) %} with the following parameters:
> - *"Expression matrix in sparse matrix format (.mtx)"*: `EBI SCXA Data Retrieval on E-MTAB-6945 matrix.mtx (Raw filtered counts)`
> - *"Gene table"*: `EBI SCXA Data Retrieval on EMTAB-6945 genes.tsv (Raw filtered counts)`
> - *"Barcode/cell table"*: `EBI SCXA Data Retrieval on E-MTAB-6945 barcodes.tsv (Raw filtered counts)`
> - *"Cell Metadata"*: `EBI SCXA Data Retrieval on E-MTAB-6945 exp_design.tsv`
> - *"Minimum cells to include features"*: `5`
> - *"Choose the format of the output"*: `RDS with a Seurat object`
>
> **Rename** {% icon galaxy-pencil %} output `Initial Seurat Object`
{: .hands_on}

The output of this tool will result in a Seurat object with row/column names as described above.

</div>

<div class='Uploading-from-Zenodo' markdown='1'>

## Get Data onto Galaxy
> <hands-on-title>Uploading from Zenodo</hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the AnnData object from [Zenodo]({{ page.zenodo_link }})
>
>    ```
>    {{ page.zenodo_link }}/files/Mito-counted_AnnData
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. **Rename** {% icon galaxy-pencil %} the datasets `Mito-counted AnnData`
>
>    {% snippet faqs/galaxy/datasets_rename.md name="Mito-counted AnnData" %}
>
> 4. Check that the datatype is `h5ad`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="h5ad" %}
>
{: .hands_on}

## Generating a Seurat object
When you uploaded your data from Zenodo, it came in AnnData format, thus we will need to convert this to a Seurat Object. This can be easily accomplished using the Seurat FilterCells tool.
Simply run the tool without any actual filtering thresholds and with the following parameters:

> <hands-on-title>Filter Cells</hands-on-title>
>
> Run{% tool [Seurat FilterCells](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_filter_cells/seurat_filter_cells/4.0.4+galaxy0) %} with the following parameters:
> - *"Choose the format of the input"*: `AnnData`
> - *"AnnData file"*: `Mito-counted AnnData`
> - *"Choose the format of the output"*: `RDS with a Seurat object`
{: .hands_on}

</div>

You've created a Seurat object, congratulations!

# QC Plots
Now that we have a complete Seurat object, we can begin the filtering process.

There will  be a number of ‘cells’ that are actually just empty droplets or low-quality. There will also be genes that could be sequencing artifacts or that appear with such low frequency that statistical tools will fail to accurately analyse them.

This background noise of both cells and genes not only makes it harder to distinguish real biological information from artifacts, but also makes it computationally demanding to analyze.

We want to filter our cells, but first we need to know what our data looks like. There are a number of subjective decisions to make within scRNA-seq analysis, for instance we now need to make informed decisions about where to set our thresholds (more on that soon!).

We’re going to plot our data a few different ways. Different bioinformaticians might prefer to see the data in different ways, and here we are only generating a few of the plots you can use. Ultimately you need to go with what makes the most sense to you.

So let's generate some QC plots. First off, let's visualize the spread of our data:

> <hands-on-title>Visualize Counts</hands-on-title>
>
> Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Initial Seurat Object` (output of **Seurat Read10x** {% icon tool %})
> - *"Plot_type_selector"*: `VlnPlot`
> - *"Features"*: `nCount_RNA`
> - *"Log"*: `Yes`
{: .hands_on}

![Violin Plot of Counts](../../images/scrna-case_FPE_SeuratTools/nCount_RNA_vln_plot.png "Violin Plot of counts.")

This plot will show us the spread of cells in our data containing a given number of counts (or transcripts) observed in a given cell. We can use this plot, and others like it in a moment, to help filter out the uninformative cells.

In a similar fashion we can visualize the spread of cells in our data expressing a given number of features (or genes):

> <hands-on-title>Visualize Features</hands-on-title>
>
> Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Initial Seurat Object` (output of **Seurat Read10x** {% icon tool %})
> - *"Plot_type_selector"*: `VlnPlot`
> - *"Features"*: `nFeature_RNA`
> - *"Log"*: `Yes`
{: .hands_on}

![Violin Plot of Features](../../images/scrna-case_FPE_SeuratTools/nFeature_RNA_vln_plot.png "Violin Plot of features.")

Now, we could just pick filtering thresholds based on these plots, and in a typical pipeline we would also plot the proportion of features that map to the mitochondrial genome (a tool coming soon to do so!). In the meantime, let's do some QC checks.

We can, and should, ask a number of questions about the quality of our data before conducting any actual analyses. Batch effect, for example carries the potential to alter the conclusions we make. Let's take a look at whether this may be the case here:  

> <hands-on-title>Visualize Counts Split by Individual</hands-on-title>
>
> Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Initial Seurat Object` (output of **Seurat Read10x** {% icon tool %})
> - *"Plot_type_selector"*: `VlnPlot`
> - *"Features"*: `nCount_RNA`
> - *"Group by"*: `Sample.Characteristic.individual.`
> - *"Split by"*: `Sample.Characteristic.individual.`
> - *"Log"*: `Yes`
> - *"Fill by"*: `ident`
{: .hands_on}

![Violin Plot of Counts Split by Individual](../../images/scrna-case_FPE_SeuratTools/nCount_split_by_Individual_vln_plot.png "Violin Plot of counts split by individual.")

This plot shows us the number of cells split by the individual (mouse) from which the cells came from. Now, depending on your experimental design, batch may be represented by something other than individual--like timepoint or even the wet lab researcher who isolated the cells.

Ideally, we would like to see a relatively even distribution of counts for each individual (or batch) but if there isn’t, fear not, we can regress this variable out in a later step.

><tip-title>Plotting Lesson</tip-title>
>In order to accurately assess potential batch effects, use the "group.by" parameter to indicate the variable which differed across experiments.
>If you are analyzing your own data, try plotting another variable which you know differed across experiments or even just samples.
>
{: .tip}

Now let's get an idea of how other variables, like  sex or genotype of the mice, might be represented across our dataset.

## 1. Sex Differences in Counts?

   > <hands-on-title>Visualize Counts Split by Sex</hands-on-title>
  >
  > Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
  > - *"RDS file"*: `Initial Seurat Object` (output of **Seurat Read10x** {% icon tool %})
  > - *"Plot_type_selector"*: `VlnPlot`
  > - *"Features"*: `nCount_RNA`
  > - *"Group by"*: `Sample.Characteristic.sex.`
  > - *"Split by"*: `Sample.Characteristic.sex.`
  > - *"Log"*: `Yes`
  > - *"Fill by"*: `ident`
  {: .hands_on}

![Violin Plot of Counts Split by Sex with males representing more of the counts than females but distributions are about the same across both sexes](../../images/scrna-case_FPE_SeuratTools/nCount_split_by_Sex_vln_plot.png "Violin Plot of counts split by sex.")

## 2. Genotype Differences in Counts?

  > <hands-on-title>Visualize Counts Split by Genotype</hands-on-title>
  >
  > Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
  > - *"RDS file"*: `Initial Seurat Object` (output of **Seurat Read10x** {% icon tool %})
  > - *"Plot_type_selector"*: `VlnPlot`
  > - *"Features"*: `nCount_RNA`
  > - *"Group by"*: `Sample.Characteristic.genotype.`
  > - *"Split by"*: `Sample.Characteristic.genotype.`
  > - *"Log"*: `Yes`
  > - *"Fill by"*: `ident`
  {: .hands_on}

  ![Violin Plot split by Genotype--Mutant versus Control](../../images/scrna-case_FPE_SeuratTools/nCount_split_by_Genotype_vln_plot.png "Violin Plot of counts split by Genotype--Mutant versus Control.")

# Finding Our Filtering Parameters
Now that we have a better understanding of what our data looks like, we can begin identifying those spurious reads and low quality cells and then remove them.

In a standard workflow, we would plot the percentage of mitochondrial reads (perc.mt) against the total number of reads per cell (nCount_RNA) and number of genes per cell (nFeature_RNA) to get an idea of our thresholds.

><comment-title>High Mitochondrial Reads</comment-title>
>High perc.mt will typically indicate stressed out cells (often due to the extraction, sorting, or sample prep protocols).
>
{: .comment}

However, this is not always necessary, and in fact, filtering based on counts and features (indeed, even just counts alone) will often remove the cells with spuriously high mitochondrial transcripts. These tools (and this tutorial) will soon be updated to allow us to do so--in the meantime, please see [Filter, Plot, and Explore single cell RNA-seq data (Seurat, R)]({% link topics/single-cell/tutorials/scrna-case_FilterPlotandExploreRStudio/tutorial.md %}) [Filter, plot and explore single-cell RNA-seq (Scanpy)]({% link topics/single-cell/tutorials/scrna-case_basic-pipeline/tutorial.md %}), or [Filter, plot and explore single-cell RNA-seq data (Scanpy, Python)]({% link topics/single-cell/tutorials/scrna-case-jupyter_basic-pipeline/tutorial.md %}) if you hope to include perc.mt in your own data analysis adventures.

For now, we will use just transcript and gene counts to filter our data. Let's take a look back at our nFeature Violin Plot to pick our gene threshold:

![Violin Plot of Features](../../images/scrna-case_FPE_SeuratTools/nFeature_RNA_vln_plot.png "Violin Plot of features.")

><comment-title>Interpretations</comment-title>
>You can see that very few cells in the dataset contain fewer than ~500 genes. Biologically, this makes sense, and the cells appear to be outliers in the data. As such, we will set our lower threshold of genes (nFeature_RNA) at 500.
{: .comment}

Now, what about transcripts (nCount_RNA)? Let's take a look:

![Violin Plot of Counts](../../images/scrna-case_FPE_SeuratTools/nCount_RNA_vln_plot.png "Violin Plot of counts.")

><comment-title>Interpretations</comment-title>
>This one is a bit more difficult to visualize but we see a severe drop off in the number of cells that contain fewer than 500 and more than 10,000 transcripts. These will be our nCount thresholds that we filter based on.
{: .comment}

These cells won't tell us much biologically, rather, they will contribute noise that we'll want to filter out of the data. With that being said, filtering scRNA-seq data will always be an iterative process--so label your work well and be ready to revisit these thresholds if your analyses seem strange down the line.

# Applying our Thresholds
It’s time to filter our cells by applying the above thresholds!

In order to include more than one parameter by which to filter, use the "Insert Subsets used to filter cells" button below the first parameter box.

> <hands-on-title>Filter Cells</hands-on-title>
>
> Run{% tool [Seurat FilterCells](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_filter_cells/seurat_filter_cells/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Initial Seurat Object` (output of **Seurat Read10x** {% icon tool %})
> - In *"Subsets used to filter cells"*:
>    - *"Name of Parameter to filter on"*: `nCount_RNA`
> - *"Min value"*: `500.0`
> - *"Max value"*: `10000`
>
> - {% icon param-repeat %} *"Insert Subsets used to filter cells"*
>     - *"Name of Parameter to filter on"*: `nFeature_RNA`
> - *"Min value"*: `500.0`
> - *"Max value"*: `1000000000.0`
>
> - *"Choose the format of the output"*: `RDS with a Seurat object`
>
> **Rename** {% icon galaxy-pencil %} output `Filtered Seurat Object`
{: .hands_on}

In this step we are creating a new Seurat object (notice that the selected output of this tool will be an RDS file as opposed to the png plots we have thus far been creating).

Now, genes that do not appear in any cell, or even in only 1 or 2 cells, may break some analytical tools and will generally not be biologically informative.

Since you’ve removed a whole heap of cells, and the captured genes are sporadic (i.e. a small percentage of the overall transcriptome per cell) this means there are a number of genes still present in your matrix that are not expressed in any of the cells.

The removal of these genes is by no means necessary, but will speed up your analyses. The developers are currently working to enable a means of doing this through the Seurat Tools, but, in the meantime if you are analyzing your own data and would like to filter genes--please see [Filter, Plot, and Explore single cell RNA-seq data (Seurat, R)]({% link topics/single-cell/tutorials/scrna-case_FilterPlotandExploreRStudio/tutorial.md %}) [Filter, plot and explore single-cell RNA-seq (Scanpy)]({% link topics/single-cell/tutorials/scrna-case_basic-pipeline/tutorial.md %}), or [Filter, plot and explore single-cell RNA-seq data (Scanpy, Python)]({% link topics/single-cell/tutorials/scrna-case-jupyter_basic-pipeline/tutorial.md %}).

# Processing
Currently, we still have quite big data. We have two issues here
 1. We already saw in our filtering plots that there are differences in how many transcripts and genes have been counted per cell. This technical variable could, and likely will, distract us from identifying true biological differences.
 2. We like to plot things on 2-dimensional X/Y plots. So, for instance, Gapdh could be on one axis, and Actin could be on another, and then each cell is plotted onto that 2D axis based on how many of each transcript they possess.

Although that would be fine, adding in a 3rd dimension (or, indeed, in our case, a dimension for each of the thousands of genes), is a bit trickier.

So, our next steps will be to transform our big data object into something that is easy to analyse and easy to visualize: this is commonly referred to as preprocessing of the data and a typical scRNA-seq preprocessing pipeline will include the following steps:

## 1. Normalization

What is Normalization?

Normalisation helps reduce the differences between gene and UMI counts by fitting total counts across cells in our data to be comparable to one another. SCTransform regularizes the gene expression profiles via a negative binomial regression while also controlling for overfitting of the data.

> <hands-on-title>Normalize Data</hands-on-title>
>
> Run{% tool [Seurat NormaliseData](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_normalise_data/seurat_normalise_data/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Filtered Seurat Object` (output of **Seurat FilterCells** {% icon tool %})
> - *"Normalisation method"*: `Log Normalise`
>
> **Rename** {% icon galaxy-pencil %} output `Normalised Seurat Object`
{: .hands_on}

## 2. Identifying Variable Genes

What are variable genes?

The datasets have loads of genes, but not all of them vary in expression from cell to cell. For instance, housekeeping genes are defined as not changing much from cell to cell, so we could remove these from our data to simplify our analyses.

The find variable genes step flags genes that *do* vary across cells to expedite future analyses and ensure that we, and Seurat, don't waste time looking for meaningful differences where they don't exist.

> <hands-on-title>Find Variable Genes</hands-on-title>
>
> Run{% tool [Seurat FindVariableGenes](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_find_variable_genes/seurat_find_variable_genes/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Normalised Seurat Object` (output of **Seurat NormaliseData** {% icon tool %})
> - *"Choose the format of the output"*: `RDS with a Seurat object`
>
> **Rename** {% icon galaxy-pencil %} output `Normalised Seurat Object with Variable Features`
{: .hands_on}

This tool will output two new pieces of data into our Galaxy history:
  1. a new Seurat object with variable features identified and flagged
  2. a tabular file with a list of these variable genes.

This gene list may be used as a sneak peak into understanding what the dataset will look like! We can begin to understand which genes are going to be driving downstream clustering of our cells and maybe even make some decisions about whether we are happy with our filtering based on this list.

## 3. Scale Data

Now we will scale the data.

What is scaling?

This is an important step to set up our data for further dimensionality reduction. It will transform the dataset such that all genes have the same variance and a zero mean. It helps negate sequencing depth differences between samples, since the gene levels across the cells become comparable.

><comment-title>Don't Worry!</comment-title>
> Note, that the differences from scaling etc. are not the values you have at the end - i.e. if your cell has average GAPDH levels, it will not appear as a ‘0’ when you calculate gene differences between clusters.
>
{: .comment}

> <hands-on-title>Scale Data </hands-on-title>
>
> Run{% tool [Seurat ScaleData](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_scale_data/seurat_scale_data/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Normalised Seurat Object with Variable Features` (output of **SeuratFindVariableGenes** {% icon tool %})
> - *"Choose the format of the output"*: `RDS with a Seurat object`
> - *"Genes to use"*: `Seurat FindVariableGenes on data 12: Variable genes tabular file`
> - *"Vars to regress"*: `nCount_RNA`
> - *"Statistical model"*: `Linear model`
>
> **Rename** {% icon galaxy-pencil %} output `Preprocessed Seurat Object`
{: .hands_on}

You now have a preprocessed Seurat object!

><comment-title>Regressing Variables</comment-title>
> Take note of the "Vars to regress" argument in the above tool. This function allow us to mitigate the effects of confounding factors in our dataset.
> In true research practice, I often regress out multiple variables including but not limited to perc.mt, cell cycle scoring, and feature count.
> As currently written, this tool only allows us to regress out a single variable: so feel free to pick another to regress and see how it changes the downstream analyses!
{: .comment}

# Dimensionality Reduction
Although we've made our expression values comparable to one another and our overall dataset less computationally demanding, we still have way too many dimensions (n cells x n genes!).

Transcript changes are not usually singular--which is to say, genes function and exist in pathways and groups. It would be easier to analyse our data if we could group these differences. To address this we will run principal component analysis (PCA).

><comment-title>What is PCA?</comment-title>
>Principal components (PCs) are calculated from highly dimensional data to find the most representative spread in the dataset. So in our highly variable gene dimensions, there will be one line (axis) that yields the most spread and variation across the cells. That will be our first principal component.
{: .comment}

We can calculate the first handful of principal components in our data to drastically reduce the number of dimensions:

><tip-title>Running Computationally Demanding Steps on Variable Features </tip-title>
>You'll notice that the **Seurat RunPCA** tool is run using the variable features from the previous step. This signficantly decreases the number of genes, and their expression changes, that must be grouped into principal components by this step.
{: .tip}

> <hands-on-title>Run PCA </hands-on-title>
>
> Run{% tool [Seurat RunPCA](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_run_pca/seurat_run_pca/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Preprocessed Seurat Object` (output of **Seurat RunPCA** {% icon tool%})
> - *"Choose the format of the output"*: `RDS with a Seurat object`
> - *"Genes to scale"*: `Seurat FindVariableGenes on data 12: Variable genes tabular file`
{: .hands_on}

This tool will output you with four new datasets into your history:
  1. Seurat RDS which includes all of the following PCA metadata
  2. Embeddings: Principal component values for each of the cells in your dataset
  3. Loadings: Prinicial component values for each of the genes in your dataset  
  4. Standard deviations of each principal component coordinates

><tip-title>Visualizing PCA</tip-title>
>In order to use the PCA information which was just calculated, we must visualize it. The currently available tools unfortunately do not [YET] carry the capacity to do so, but I will provide the plot so that we may make an informed decision together:
>
>![Elbow Plot of Principal Components (PCs) against the standard deviation (variability in the dataset) they each account for](../../images/scrna-case_FPE_SeuratTools/ElbowPlot.png "Elbow Plot")
>
>We can see that there is really not much variation explained past the 9th PC. So we might save ourselves a great deal of time and muddied data by focusing on the top 15 PCs to be conservative.
>You can also think about it like choosing a threshold of variance explained. Conservatively, 2.5 standard deviations are explained by about 10 of the PCs.
{: .tip}

We’re still looking at around 15 dimensions at this point--likely not the easiest to visualize. To make our lives easier, we need to identify how similar a cell is to another cell, across every cell across each of these dimensions.

For this, we will use the k-nearest neighbor (kNN) graph, to identify which cells are close together and which are not.

The kNN graph plots connections between cells if their distance (when plotted in this 10 dimensional space!) is amongst the k-th smallest distances from that cell to other cells. This will be crucial for identifying clusters, and is necessary for plotting a UMAP--which is what will ultimately allow us to visualize our data in 2 dimensions.

><comment-title>From UMAP developers:</comment-title>
>“Larger neighbor values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50, with a choice of 10 to 15 being a sensible default”.
{: .comment}

Let's now use the 15 PC threshold we chose from the Elbowplot and apply it to find neighbors:

> <hands-on-title>Find Neighbors </hands-on-title>
>
> Run{% tool [Seurat FindNeighbours](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_find_neighbours/seurat_find_neighbours/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Preprocessed Seurat Object` (output of **Seurat RunPCA** {% icon tool %})
> - *"Reduction"*: `pca`
> - *"Dimensions"*: `1,2,3,4,5,6,7,8,9,10,11,12,13,14,15`
> - *"Assay"*: `RNA`
{: .hands_on}

Now we can use the neighborhood graph to identify clusters of cells whose transcriptional profiles appear most similar to one another: we can identify and label clusters:

> <hands-on-title>Find Clusters </hands-on-title>
>
> Run{% tool [Seurat FindClusters](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_find_clusters/seurat_find_clusters/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Preprocessed Seurat Object` (output of **Seurat FindNeighbors** {% icon tool %})
> - In *"Advanced Options "*
>   - *"Resolution"*: `0.5`
{: .hands_on}

This tool will output two new datasets: as usual, a new Seurat object which includes a metadata column denoting which cluster each cell was assigned to, and a csv file of the same information.

Unfortunately, identifying clusters is not as majestic as biologists often think - the math doesn’t necessarily identify true cell clusters. Every algorithm for identifying cell clusters falls short of a biologist knowing their data, knowing what cells should be there, and proving it in the lab.

So, we’re going to make the best of it as a starting point and see what happens! We will define clusters from the kNN graph, based on how many connections cells have with one another. Roughly, this will depend on a resolution parameter for how granular you want to be.

><tip-title>On Clustering Resolution</tip-title>
>The resolution parameter available in the advanced options of this tool allows for you, the bioinformatician, to dictate the granularity of the clusters.
>
>For example, a higher clustering resolution dictates increased granularity, and more stringent clusters. That is--cells must more closely resemble one another in order to be grouped into the same cluster than at a lower clustering resolution.
>
>In general, I find it easiest to think of a higher resolution producing more clusters and conversely, a lower resolution will produce less clusters. This parameter is a useful one that you will use often to help decipher how many true populations of cells are present in your data!
{: .tip}

Now that we have made note within our object of which cells cluster together, we can start to visualize our data! Two major visualizations for this data type are tSNE and UMAP. We can calculate the coordinates for both prior to visualization. For tSNE, the parameter perplexity can be changed to best represent the data, while for UMAP the main change would be to change the kNN graph above itself, via the **Seurat FindNeighbors** tool.

><tip-title>On UMAP</tip-title>
>UMAP is the most recently developed, and most widely used dimensionality reduction for visualization of principal component data. It has been optimized since tSNE to better preserve global structure and is less computationally demanding.
>
{: .tip}

> <hands-on-title>Run UMAP </hands-on-title>
>
> Run{% tool [Seurat UMAP](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_run_umap/seurat_run_umap/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Preprocessed Seurat Object` (output of **Seurat FIndClusters** {% icon tool %})
> - *"Choose the format of the output"*: `RDS with a Seurat object`
> - *"Dims"*: `1:15`
>
> **Rename** {% icon galaxy-pencil %} output `Final Preprocessed Seurat Object`
{: .hands_on}

You now have a completely preprocessed and ready to be analyzed Seurat object--congratulations!

# Let's Take a Look
Now that we have run dimensionality reduction on our dataset, it is ready for visualization. Let's take a look at what our cells look like in a UMAP projection:

> <hands-on-title>Plot UMAP </hands-on-title>
>
> Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Final Preprocessed Seurat Object` (output of **Seurat UMAP** {% icon tool %})
> - *"Plot_type_selector"*: `DimPlot`
> - *"Group by"*: `RNA_nn_res.0.5`
{: .hands_on}

![DimPlot colored by 0.5 resolution clusters with a total of 8 clusters--two noticeably different "island" populations and six clustering more closely to one another](../../images/scrna-case_FPE_SeuratTools/DimPlot_GroupBy_pt5Res.png "DimPlot colored by 0.5 resolution cluster.")

Good work! It looks like with a clustering resolution of 0.5, we are able to identify 8 clusters of cells in our data.

We can also look for expression of particular genes and see how those map to our UMAP projection. This is often useful in getting an initial understanding of which clusters might be representative of which cell types.

> <hands-on-title>Plot Gapdh </hands-on-title>
>
> Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Final Preprocessed Seurat Object` (output of **Seurat UMAP** {% icon tool %})
> - *"Plot_type_selector"*: `FeaturePlot`
> - *"Features"*: `Gapdh`
{: .hands_on}

![FeaturePlot of Gapdh expression shown boradly across the whole dataset](../../images/scrna-case_FPE_SeuratTools/FeaturePlot_Gapdh.png "FeaturePlot: Gapdh")

We just plotted a housekeeping gene, Gapdh, so the broad expression we observe is expected.

In practice, it is helpful to plot known markers of cell types you expect to be in your dataset. This will give you a first look at how your cells are clustered.

For example, we can plot early T-cell marker Il2ra and get an idea of which cells and/or clusters might resemble the early T-cells:

> <hands-on-title>Plot Il2ra </hands-on-title>
>
> Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Final Preprocessed Seurat Object` (output of **Seurat UMAP** {% icon tool %})
> - *"Plot_type_selector"*: `FeaturePlot`
> - *"Features"*: `Il2ra`
{: .hands_on}

![FeaturePlot of Il2ra expression most highlky expressed in the lower island population (cluster 2)](../../images/scrna-case_FPE_SeuratTools/FeaturePlot_Il2ra.png "FeaturePlot: Il2ra")

It is a good idea, when analyzing your own data, to plot some markers of cell types you expect to be present. Later on we can also use these FeaturePlots to visualize manual annotation of clusters.

# Differential Expression Testing: Finding Markers
Because each cluster of cells was grouped based on similar transcriptome profiles, each cluster will inherently differ from one another based on a set of "marker" genes.

Following an initial look at the DimPlots and FeaturePlots, we can take an even closer look at which genes are driving the clustering.

To do so, we'll run Seurat's FindMarkers function, which will compare each identity (in this case cluster) against every other identity within its class (all the other clusters). This function of marker finding is particularly useful in identifying up, or down, regulated genes that drive differences in identity/cluster.

> <hands-on-title>Find Markers </hands-on-title>
>
> Run{% tool [Seurat FindMarkers](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_find_markers/seurat_find_markers/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Final Preprocessed Seurat Object` (output of **Seurat UMAP** {% icon tool %})
{: .hands_on}

The marker list that has been output by this tool will be useful to us shortly for identifying which cells types are represented by the various clusters.\

><comment-title>On Finding Markers</comment-title>
> Differential expression can be run using any identity class and any two identities within the same class. The tools as of current do not allow for such direct comparisons (but will soon!) and I implore you to try [Filter, Plot, and Explore single cell RNA-seq data (Seurat, R)]({% link topics/single-cell/tutorials/scrna-case_FilterPlotandExploreRStudio/tutorial.md %}) if you hope to conduct marker identification across, say, genotypes.
{: .comment}

# Biological Interpretations
Now it’s the fun bit! We can see where genes are expressed, and start considering and interpreting the biology of it. At this point, it’s really about what information you want to get from your data--the following is only the tip of the iceberg. However, a brief exploration is good, because it may help give you ideas going forward for your own data. Let's start interrogating our data!

Let's take another look at what our clusters look like:

> <hands-on-title>Plot UMAP </hands-on-title>
>
> Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Final Preprocessed Seurat Object` (output of **Seurat UMAP** {% icon tool %})
> - *"Plot_type_selector"*: `DimPlot`
> - *"Group by"*: `RNA_nn_res.0.5`
{: .hands_on}

![DimPlot colored by 0.5 resolution cluster showing 8 total clusters](../../images/scrna-case_FPE_SeuratTools/DimPlot_GroupBy_pt5Res.png "DimPlot colored by 0.5 resolution cluster.")

><comment-title>On Cluster Numbering</comment-title>
>Note that Seurat's cluster numbering is based on size alone, so clusters 0 and 1 are not necessarily related, they are just the clusters containing the most cells.
{: .comment}

It would be nice to know what these cells are. This analysis (googling all of the marker genes, both checking where the ones you know are and then going through marker tables we generated) is a fun task for any individual experiment, so we’re going to speed past that and nab the assessment from the original paper!

Here are the markers per cell type that the paper uses to classify the avrious cell types which are expected to be present in the data:

| Markers                 | Cell Type                           |
|-------------------------|-------------------------------------|
| Il2ra                   | Double negative (early T-cell)      |
| Cd8b1, Cd8a, Cd4        | Double positive (middle T-cell)     |
| Cd8b1, Cd8a, Cd4 - high | Double positive (late middle T-cell)|
| Itm2a                   | Mature T-cell                       |

We can plot these markers as a means of discerning which cluster might be representing a given cell type. Let's start with the Mature T-cell marker Itm2a:

> <hands-on-title>Plot Itm2a </hands-on-title>
>
> Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Final Preprocessed Seurat Object` (output of **Seurat UMAP** {% icon tool %})
> - *"Plot_type_selector"*: `FeaturePlot`
> - *"Features"*: `Itm2a`
{: .hands_on}

![FeaturePlot of Itm2a expression most highly localized to cluster 3](../../images/scrna-case_FPE_SeuratTools/FeaturePlot_Itm2a.png "FeaturePlot of Itm2a")

We can see that Cluster 3 seems to be most highly expressing this gene, and when we look back at the marker list we created earlier, Itm2a appears as a marker for cluster 3 with an average log fold change of 3! Therefore, we can quite conficently say that cluster 3 is likely to be representing Mature T-cells!

Now what about the opposite end of the spectrum: the double negative early T-cells? We can see from our table about that Il2ra is a known marker of this cell type, let's see where it's most strongly expressed in the data:

> <hands-on-title>Plot Il2ra </hands-on-title>
>
> Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Final Preprocessed Seurat Object` (output of **Seurat UMAP** {% icon tool %})
> - *"Plot_type_selector"*: `FeaturePlot`
> - *"Features"*: `Il2ra`
{: .hands_on}

![FeaturePlot of Il2ra expression mainly in cluster 2](../../images/scrna-case_FPE_SeuratTools/FeaturePlot_Il2ra.png "FeaturePlot of Il2ra")

Il2ra expression seems to be most prominent in Cluster 2 and as seen in our marker list, has a fold change of 2.79!

It can, and should, act as a sanity check that the most (and least) differentiated cell types expected in the data appear as two "island" clusters.

Now for the intermediate populations--which may be a bit more tricky to deconvolute. Based on our known markers, we see that both the double-positive populations express many of the same genes: Cd8b1, Cd8a, and Cd4. Let's make sure that the remaining clusters (0, 1, 3, 4, 5, 6, and 7) all express these:

> <hands-on-title>Plot Cd8b1 </hands-on-title>
>
> Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Final Preprocessed Seurat Object` (output of **Seurat UMAP** {% icon tool %})
> - *"Plot_type_selector"*: `FeaturePlot`
> - *"Features"*: `Cd8b1`
{: .hands_on}

![FeaturePlot of Cd8b1 expression expressed across each of the "body" clusters (0, 1, 4, 5, 6, and 7)](../../images/scrna-case_FPE_SeuratTools/FeaturePlot_Cd8b1.png "FeaturePlot of Cd8b1")

It looks like clusters 1, 4, 5, and 6 pretty strongly express Cd8b1, now what about Cd8a?

> <hands-on-title>Plot Cd8a </hands-on-title>
>
> Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Final Preprocessed Seurat Object` (output of **Seurat UMAP** {% icon tool %})
> - *"Plot_type_selector"*: `FeaturePlot`
> - *"Features"*: `Cd8a`
{: .hands_on}

![FeaturePlot of Cd8a expression in the same clsuters as Cd8b1 (0, 1, 4, 5, 6, and 7)](../../images/scrna-case_FPE_SeuratTools/FeaturePlot_Cd8a.png "FeaturePlot of Cd8a")

This looks pretty consistent with the Cd8b1 plot, which is expected as these are markers of a double positive population! The true discernment between the two populations will be the level of Cd4 expression:

> <hands-on-title>Plot Cd4 </hands-on-title>
>
> Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Final Preprocessed Seurat Object` (output of **Seurat UMAP** {% icon tool %})
> - *"Plot_type_selector"*: `FeaturePlot`
> - *"Features"*: `Cd4`
{: .hands_on}

![FeaturePlot of Cd4 expression most highly expressed in cluster 1 and lowly expressed in clusters 4, 5, 6, and 7](../../images/scrna-case_FPE_SeuratTools/FeatuerePlot_Cd4.png "FeaturePlot of Cd4")

Looks like the top portion of these clusters, and mainly cluster 1 hold the vast majority of Cd4 expression. This, coupled with Cd4 being identified as a marker of cluster 1 in our marker table above tells us pretty confidently that the upper part of this "body" cluster are likely the late middle t-cells!

Here is a breakdown of the cell type labeling we just accomplished:

| Clusters     | Markers                 | Cell Type                           |
|--------------|-------------------------|-------------------------------------|
| 2            | Il2ra                   | Double negative (early T-cell)      |
| 0            | Il2ra, Cd8b1, Cd8a      | [maybe] Early transitioning T-cell  |
| 4, 5, 6, 7   | Cd8b1, Cd8a, Cd4        | Double positive (middle T-cell)     |
| 1            | Cd8b1, Cd8a, Cd4 - high | Double positive (late middle T-cell)|
| 3            | Itm2a                   | Mature T-cell                       |

Now we can begin to feel a bit more oriented in exploring our data. The clusters are labelled with cell types, and our object has been processed enough such that we may now begin to answer some real biological questions! Now that we know what we’re dealing with, let’s examine the effect of our variable, real science!

## Keep Digging
Are there any differences in genotype? Or in biological terms, is there an impact of growth restriction on T-cell development in the thymus? We can begin to answer this question visually by using the "split.by" parameter in Seurat's plot functions.

> <hands-on-title>Visualize Counts Split by Genotype</hands-on-title>
>
> Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Final Preprocessed Seurat Object` (output of **Seurat UMAP** {% icon tool %})
> - *"Plot_type_selector"*: `VlnPlot`
> - *"Features"*: `nCount_RNA`
> - *"Group by"*: `RNA_nn_res.0.5`
> - *"Split by"*: `Sample.Characteristic.genotype.`
> - *"Log"*: `Yes`
> - *"Fill by"*: `ident`
{: .hands_on}

![VlnPlot colored by cluster and split by genotype showing the mutant clusters appearing more sparse than the control](../../images/scrna-case_FPE_SeuratTools/VlnPlot_GroupedbyCluster_SplitbyGenotype.png "Violin Plot colored by assigned clusters & split by genotype")

We can see that there seems to be a decrease in cellcounts across the celltypes in the het mutant... INTERESTING! What next? We might look further at the transcripts present in both those populations, and perhaps also look at the genotype marker table… So much to investigate! But before we set you off to explore to your heart’s delight, let’s also look at this a bit more technically.

# Technical Assessment
Is our analysis real? Is it right? Well, we can assess that a little bit.

First thing's first, is there a batch effect?

> <hands-on-title>Visualize Counts Split by Individual</hands-on-title>
>
> Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Final Preprocessed Seurat Object` (output of **Seurat UMAP** {% icon tool %})
> - *"Plot_type_selector"*: `DimPlot`
> - *"Group by"*: `RNA_nn_res.0.5`
> - *"Split by"*: `Sample.Characteristic.individual.`
{: .hands_on}

![DimPlot colored by labelled celltype split by individual/batch showing that individual 3 and 4 contribute the vast majority of immature T-cells](../../images/scrna-case_FPE_SeuratTools/DimPlotColorbyClusters_SplitbyIndividual.png "DimPlot colored by assigned clusters split by individual/batch")

While some differences across batch are expected and nothing to be concerned about, the immature T-cells looks to be mainly comprised of Individuals 3 and 4. There might be a bit of batch effect, so you could consider using batch correction on this dataset. However, if we focus our attention on the other cluster - mature T-cells - where there is batch mixing, we can still assess this biologically even without batch correction.

Additionally, we will also look at the confounding effect of sex:

> <hands-on-title>Visualize Counts Split by Sex</hands-on-title>
>
> Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Final Preprocessed Seurat Object` (output of **Seurat UMAP** {% icon tool %})
> - *"Plot_type_selector"*: `DimPlot`
> - *"Group by"*: `RNA_nn_res.0.5`
> - *"Split by"*: `Sample.Characteristic.sex.`
{: .hands_on}

![DimPlot colored by cluster and split by sex showing that there are fewere female cells (a product of fewer female samples being present in the data as a whole)](../../images/scrna-case_FPE_SeuratTools/DimPlot_SplitbySex.png "DimPlot colored by assigned clusters and split by Sex")

We note that the one female sample - unfortunately one of merely three knockout samples - seems to be distributed in the same areas as the knockout samples at large, so luckily, this doesn’t seem to be a confounding factor and we can still learn from our data. Ideally, this experiment would be re-run with either more female samples all around or swapping out this female from the male sample.

Are there any clusters or differences being driven by sequencing depth, a technical factor?

> <hands-on-title>Plot Counts </hands-on-title>
>
> Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Final Preprocessed Seurat Object` (output of **Seurat UMAP** {% icon tool %})
> - *"Plot_type_selector"*: `FeaturePlot`
> - *"Features"*: `nCount_RNA`
{: .hands_on}

![FeaturePlot colored by counts showing no differences in sequencing depth across the dataset](../../images/scrna-case_FPE_SeuratTools/FeaturePlot_nCount.png "FeaturePlot colored by counts")

There doesn't visually appear to be any differences in sequencing depth across the clusters, but let's check out some of those other variables we grouped by:

> <hands-on-title>Plot Features </hands-on-title>
>
> Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Final Preprocessed Seurat Object` (output of **Seurat UMAP** {% icon tool %})
> - *"Plot_type_selector"*: `FeaturePlot`
> - *"Features"*: `nCount_RNA`
> - *"split.by"*: `Sample.Characteristic.individual.`
{: .hands_on}

![FeaturePlot colored by counts & Split by Individual showing that Individual 3 seems to have been most deeply sequenced](../../images/scrna-case_FPE_SeuratTools/FeaturePlot_Counts_SplitbyIndividual.png "FeaturePlot colored by counts & split by individual")

There we go! This might explain the dramatic shift in early to middle T-Cell between wildtype and knockout cells--the leftmost early to middle T-cells simply have a higher sequencing depth represented by Individual 3 (UMIs/cell) than the ones on the right side. Well, that explains some of the sub-cluster that we’re seeing in that splurge (specifically this likely accounts for the discernment between clusters 0, 4, 5, 6, and 7).

Luckily, and importantly, we don’t see the double negative or mature T-cells being similarly affected. So, although, this variable of sequencing depth, or moreso, Individual, might be something to regress out somehow, it doesn’t seem to be impacting our dataset such that we cannot draw meaningful insights.

><tip-title>Overprocessing</tip-title>
>The less you can regress/modify your data, in general, the better--you want to stay as true as you can to the raw data, and only use maths to correct your data when you really need to (and not to create insights where there are none!).
{: .tip}

Do you think we processed these samples well enough? We have seen in the previous images that these clusters are not very tight or distinct, so we could consider stronger filtering. Let's take a look at gene expression of a gene we know should not be expressed in tCells as a sanity check:

> <hands-on-title>Plot Hba-a1 </hands-on-title>
>
> Run{% tool [Plot with Seurat](toolshed.g2.bx.psu.edu/repos/ebi-gxa/seurat_plot/seurat_plot/4.0.4+galaxy0) %} with the following parameters:
> - *"RDS file"*: `Final Preprocessed Seurat Object` (output of **Seurat UMAP** {% icon tool %})
> - *"Plot_type_selector"*: `FeaturePlot`
> - *"Features"*: `Hba-a1`
{: .hands_on}

![FeaturePlot of Hba-a1 showing sparse expression across the whole dataset](../../images/scrna-case_FPE_SeuratTools/FeaturePlot_Hba-a1.png "FeaturePlot of Hba-a1")

Hemoglobin--a red blood cell marker that should NOT be found in T-cells--appears throughout the entire dataset in low numbers and as a likely marker of Cluster 6. This suggests that some background noise may have been introduced by the media the cells were in. We might consider in the wet lab trying to get a purer, happier sample, with less background or in the dry lab, we can take advantage of techniques such as SoupX or others to remove this technical noise.

><tip-title>Removing Noise</tip-title>
>Adjusting the filtering settings (increasing minimum counts/cell, etc.) is often the place to start in these scenarios.
{: .tip}

Do you think the clustering is appropriate? i.e. are there single clusters that you think should be separate, and multiple clusters that could be combined?

Important to note, lest all bioinformaticians combine forces to attack the biologists: just because a cluster doesn’t look like a cluster by eye is NOT enough to say it’s not a cluster! But looking at the biology here, we struggled to find marker genes to distinguish the double positive populations, which we know are also affected by depth of sequencing. That’s a reasonable argument that Clusters 0, 4, 5, 6, and 7 might not be all that different. Maybe we need more depth of sequencing across all those cells, or to compare these explicitly to each other (consider variations on FindMarkers!).

However, the late double positive cluster is both seemingly leaving the larger body of clusters and also has fewer knockout cells, so we might go and look at what those cells are expressing in the marker genes. If we look at the mature T-cells further, we can see that their marker gene--Itm2a--is only expressed in half of the cluster. You might consider sub-clustering this to investigate further, either through changing the resolution or through analysing this cluster alone.

If we look at the differences between genotypes alone (so the pseudo-bulk), we can see that many, if not most, of the genes in that list are actually ribosomal. This could be a housekeeping background, it might be cell cycle related, it may be biological, or some combination of all three. You might consider investigating the cycling status of the cells, or even regressing this out (which is what the authors did).

Ultimately, there are quite a lot ways to analyse your single-cell data, both within the confines of this tutorial (the many parameters that could be changed throughout) and outside of it (batch correction, sub-clustering, cell-cycle scoring, inferred trajectories, etc.) Most analyses will still yield the same general output, though: there are fewer knockout cells in the mature T-cell population, suggesting some sort of abberant development of T-cells in the Igf2-p0 hets.

Congratulations! You have interpreted your plots in several important ways!
