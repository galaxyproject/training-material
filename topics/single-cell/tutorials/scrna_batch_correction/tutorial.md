---
layout: tutorial_hands_on

title: "Batch Correction and Integration with Seurat or Scanpy"
subtopic: tricks
zenodo_link: 'https://zenodo.org/records/20574474'
priority: 4

answer_histories:
  - label: Using Seurat
    history: https://singlecell.usegalaxy.eu/u/marisa_jl/h/batch-correction--integration-with-seurat---answer-key
    date: 2025-01-21

questions:
- What is the difference between batch correction and integration?
- How can we perform batch correction or integration using the Scanpy and Seurat pipelines?

objectives:
- Understand what batch correction and integration are and how they are different
- Know when to perform batch correction or integration on single cell data
- Perform batch correction or integration using either the Scanpy or Seurat pipelines

time_estimation: 2H

key_points:
- Batch correction can remove some of the technical differences between experimental groups
- Integration enables us to combine different datasets
- Both the Scanpy and Seurat pipelines include tools that can be used for batch correction or integration
requirements:
-
    type: "internal"
    topic_name: single-cell
    tutorials:
        - scrna-scanpy-pbmc3k
        - scrna-seurat-pbmc3k

tags:
- 10x

contributions:
  authorship:
    - MarisaJL
  editing:
    - nomadscientist
    - dianichj
  funding:
    - biofair

---

Single cell analyses can be complex. We may have data from different experimental batches, perhaps because we ran our experiments at different times, in different labs, or using different sequencing platforms. Sometimes we might want to combine multiple datasets, for example if we want to compare our own experimental data to a similar public dataset.

Running different batches or combined datasets through a clustering pipeline like Scanpy or Seurat without any pre-processing may not yield useful results. This is because clustering works by identifying genes with the largest differences in expression and grouping cells that share similar expression patterns. When data comes from multiple experimental batches or studies, however, the largest sources of variation often reflect technical differences between batches rather than the biology of interest like cell type. As a result, clusters may end up capturing batch or dataset identity rather than anything biologically meaningful.

To look beyond these technical differences, we can perform batch correction or integration. Both Scanpy and Seurat include tools for correcting differences between experimental batches and integrating datasets, and in practice we use the same tools for both. In this tutorial, you will learn how to use these tools in either pipeline. The choice of Scanpy or Seurat is yours.

> <comment-title></comment-title>
>
> This tutorial is based on the [Introduction to scRNA-seq integration](https://satijalab.org/seurat/articles/integration_introduction) and [Integrative analysis in Seurat v5](https://satijalab.org/seurat/articles/seurat5_integration) tutorials.
>
{: .comment}


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Important tips for easier analysis

{% snippet topics/single-cell/faqs/single_cell_omics.md %}

{% snippet faqs/galaxy/tutorial_mode.md %}

{% snippet faqs/galaxy/analysis_troubleshooting.md sc=true %}

# Batch Correction or Integration?

We will often need to perform batch correction or integration when working with different experimental batches, donors, conditions, or datasets. We need to look beyond the technical differences between them and batch correction or integration are the techniques we use to do this. Both work by identifying cell subpopulations that are shared across groups, effectively matching cells of similar types or states.

The terms batch correction and integration are closely related and are often used interchangeably, since they refer to the same underlying process and use the same tools in the same way. The distinction is mainly one of context: batch correction typically refers to removing unwanted technical variation between groups from a single study, such as different experimental batches, while integration refers to aligning cell populations across separate datasets from multiple studies.

In this tutorial, you can choose whether you want to use the Scanpy or Seurat pipelines for clustering and batch correction. Scanpy and Seurat's integration tools will create a dimensional reduction that captures the shared sources of variation across the batches or datasets. The dimensional reduction can be used to find clusters or produce visualisations such as UMAP.

# Scanpy or Seurat?

Scanpy and Seurat are two of the most commonly used pipelines for analysing single cell data. Both include tools for preprocessing, dimensional reduction (such as PCA), neighbourhood graph construction, and clustering. Clustering is often a key goal in single cell analysis, as it groups cells with similar expression profiles. These groups frequently correspond to specific cell types or states, making the data easier to interpret and understand.

Although both pipelines follow the same basic steps, small differences in how those steps are implemented mean results can vary slightly depending on your choice. Broadly speaking, though, you should reach the same conclusions either way.

The main difference between these two pipelines is that Scanpy is written for Python while Seurat is written for R. If we were working in a Python or R environment, then we would need to choose the appropriate pipeline. However, since we're working on Galaxy, we're free to choose either set of tools.

## Get data

{% include _includes/cyoa-choices.html option1="Scanpy" option2="Seurat" default="Scanpy" text="Choose the single cell pipeline you want to use." %}

<div class='Scanpy' markdown='1'>

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
>
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from the shared data library
>
>    ```
>    {{ page.zenodo_link }}/files/Input_Anndata.h5ad
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>    
> 3. Check that the datatype is `h5ad`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

</div>

<div class='Seurat' markdown='1'>

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
>
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from the shared data library
>
>    ```
>    {{ page.zenodo_link }}/files/Input_SeuratObject%20.rds
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Check that the datatype is `rds`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

</div>

You should now have either an AnnData or SeuratObject dataset in your history. Both datasets contain the same information, but in the different formats required by the Scanpy and Seurat pipelines.

We are using a single cell dataset of human Peripheral Blood Mononuclear Cells (PBMCs) that was also used in Seurat's [Integrative analysis in Seurat v5](https://satijalab.org/seurat/articles/seurat5_integration) tutorial. The original study compared the results from seven different single cell and single nuclear techniques {% cite Ding2020 %}.

Let's take a look at our data before we begin the analysis to see whether we might need to perform batch correction or integration. While batches or combined datasets often do require correction, we should examine the data first to check whether it is necessary.

Correcting batch effects is important, but we also need to be careful not to overcorrect. Applying batch correction or integration when it isn't needed, or overcorrecting when it is, risks removing the very biological differences we are interested in. Either way, it is important to interpret results carefully to determine whether they reflect true biological variation or technical artefacts.

<div class='Scanpy' markdown='1'>

> <hands-on-title> Inspect the Data </hands-on-title>
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.11.4+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output` (Input dataset)
>    - *"What to inspect?"*: `General information about the object`
>
{: .hands_on}

> <question-title></question-title>
>
> 1. How many cells and genes are in this dataset?
>
> > <solution-title></solution-title>
> >
> > 1. If we click on the {% icon galaxy-eye %} of the new output in our history, we can see that this dataset contains information about the expression of 33,694 vars (genes) in 10,434 obs (cells).
> >
> {: .solution}
>
{: .question}

Now we'll take a closer look at the metadata describing how the dataset was produced. It can tell us whether the dataset is made up of different batches that might require correction.

> <hands-on-title> Inspect the Cell Metadata </hands-on-title>
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.11.4+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output` (Input dataset)
>    - *"What to inspect?"*: `Key-indexed observations annotation (obs)`
> 2. {% tool [Count](Count1) %} with the following parameters:
>    - {% icon param-file %} *"from dataset"*: `obs` (output of **Inspect AnnData** {% icon tool %})
>    - *"Count occurrences of values in column(s)"*: `c['11']`
{: .hands_on}

</div>

<div class='Seurat' markdown='1'>

> <hands-on-title> Inspect the Data </hands-on-title>
>
> 1. {% tool [Seurat Data Management](toolshed.g2.bx.psu.edu/repos/iuc/seurat_data/seurat_data/5.0+galaxy0) %} with the following parameters:
>    - *"Method used"*: `Inspect Seurat Object`
>        - *"Display information about"*: `General`
>
{: .hands_on}

> <question-title></question-title>
>
> 1. How many cells and genes are in this dataset?
> 2. How many layers are in this dataset?
>
> > <solution-title></solution-title>
> >
> > 1. If we click on the {% icon galaxy-eye %} of the new output in our history, we can see that this dataset contains information about the expression of 33,694 features (genes) in 10,434 samples (cells).
> > 2. Data in SeuratObjects are stored in layers. In this case, we only have one layer called `counts`. The `counts` layer is the raw data that we'll be using in this tutorial.
> >
> {: .solution}
>
{: .question}

Now we'll take a closer look at the metadata describing how the dataset was produced. It can tell us whether the dataset is made up of different batches that might require correction.

> <hands-on-title> Inspect the Cell Metadata </hands-on-title>
>
> 1. {% tool [Seurat Data Management](toolshed.g2.bx.psu.edu/repos/iuc/seurat_data/seurat_data/5.0+galaxy0) %} with the following parameters:
>    - *"Method used"*: `Inspect Seurat Object`
>        - *"Display information about"*: `Cell Metadata`
> 2. {% tool [**Count** occurrences of each record](Count1) %} with the following parameters:
>    - *"Count occurrences of values in column(s)"*: `c11` This is the 11th column in your table, which contains the `Method` metadata
{: .hands_on}

</div>

> <question-title></question-title>
>
> 1. What does the `Method` column represent in the cell metadata?
> 2. Do you think batch correction or integration is needed for this analysis?
>
> > <solution-title></solution-title>
> >
> > 1. The dataset that we're using comes from a study that compared different single cell techniques. The `Method` column tells us which technique was used on each cell. Use the {% icon galaxy-eye %} to look at both outputs. The first one shows the metadata for all cells, with `Method` in column 11. The second output shows how many times each method appeared in this column. We can see there are nine different `Method` batches (as well as the `Method` heading which the tool has counted too!). Three of the batches used the same 10X Chromium (v2) method, but they appear to have been processed separately as they have been placed in different batches. We have nine entires in the `Method` column that represent nine batches using seven different experimental techniuques.
> > 2. Each experimental technique can be considered as its own experimental batch, as can the three different batches using the 10X Chromium (v2) method. Each of these batches was processed independently, which by itself can be enough to require batch correction, even if the same experimental protocol is used. Batches can vary simply because they were processed at different times or by different people in the same lab! In this case, we have an even stronger reason to believe that these batches will differ - we know that these batches were produced using different techniques. It seems likely that we'll need to perform batch correction, but we'll check what happen when we cluster without correction first. Batches often require correction, but we should always examine the data first to be sure. If we decide that correction is needed, we would consider this to be batch correction rather than integration because these data all came from the same original study.
> >
> {: .solution}
>
{: .question}

> <comment-title></comment-title>
>
> The cell metadata is any information about the cells that the original authors have included with the dataset. As well as the cell barcode or identifier for each individual cell, the metadata will usually include information such as which donor or sample the cell came from, or which experimental group it was in. Sometimes, this metadata will include a lot of useful details, such as demographic information about human donors. This information can help us to better understand our results.
{: .comment}

# Clustering without Batch Correction

We suspect that batch correction will be needed because of the different technologies used to construct this dataset, but we'll try clustering without any correction first. It is considered good practice to perform this uncorrected clustering to confirm whether batch correction is truly needed. We don't want to perform a correction unless there are technical differences between batches that need to be removed, otherwise we risk overcorrecting our data and eliminating the biological differences we're interested in. We will check whether we need to correct on the basis of `Method`. Comparing the results we get now with those we'll get after batch correction should also help us to understand what batch correction is doing to our single cell data.

Since our focus is batch correction/integration, we won't go into too much detail on the clustering process. We just want to see how the integration steps fit into the main clustering pipeline and understand the impact it has on our data. If you aren't already familiar with this process, then you can learn more about clustering using the [Scanpy]({% link topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.html %}) or [Seurat]({% link topics/single-cell/tutorials/scrna-seurat-pbmc3k/tutorial.html %}) pipelines from the other single cell tutorials available on the GTN.

<div class='Scanpy' markdown='1'>

We'll follow the default Scanpy pipeline here, except that we'll use `30` PCs to build the neighborhood graph and cluster with a resolution of `2` as these were the parameters used in [the original Seurat version of this tutorial](https://satijalab.org/seurat/articles/seurat5_integration).

> <hands-on-title> Cluster without Batch Correction </hands-on-title>
>
> 1. {% tool [Scanpy normalize](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_normalize/scanpy_normalize/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output` (Input dataset)
>    - *"Method used for normalization"*: `Normalize counts per cell, using 'pp.normalize_total'`
>        - *"Exclude (very) highly expressed genes for the computation of the normalization factor (size factor) for each cell"*: `No`
>
>    > <comment-title></comment-title>
>    >
>    > We will use the output from `Scanpy normalize` in the  following section when we perform batch correction.
>    >
>    > If you're already familiar with the Scanpy clustering pipeline and you just want to try using the  {% icon tool %} **Scanpy remove confounders** tools, then you can skip ahead to the **Clustering after Integration** step now. In a real analysis, it would be best to complete the clustering without batch correction first to check if it is needed.
>    {: .comment}
>
> 2. {% tool [Scanpy Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy normalize** {% icon tool %})
>    - *"Method used for inspecting"*: `Logarithmize the data matrix, using 'pp.log1p'`
>
> 3. {% tool [Scanpy filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used for filtering"*: `Annotate (and filter) highly variable genes, using 'pp.highly_variable_genes'`
>        - *"Choose the flavor for identifying highly variable genes"*: `Seurat`
>
> 4. {% tool [Scanpy Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy filter** {% icon tool %})
>    - *"Method used for inspecting"*: `Scale data to unit variance and zero mean, using 'pp.scale'`
>        - *"Maximum value"*: `10.0`
>
> 5. {% tool [Scanpy cluster, embed](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_cluster_reduce_dimension/scanpy_cluster_reduce_dimension/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used"*: `Computes PCA (principal component analysis) coordinates, loadings and variance decomposition, using 'pp.pca'`
>        - *"Type of PCA?"*: `Full PCA`
>
> 6. {% tool [Scanpy Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy cluster, embed** {% icon tool %})
>    - *"Method used for inspecting"*: `Compute a neighborhood graph of observations, using 'pp.neighbors'`
>        - *"Number of PCs to use"*: `30`
>
> 7. {% tool [Scanpy cluster, embed](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_cluster_reduce_dimension/scanpy_cluster_reduce_dimension/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used"*: `Cluster cells into subgroups, using 'tl.louvain'`
>        - *"Flavor for the clustering"*: `vtraag (much more powerful than igraph)`
>            - *"Resolution"*: `2.0`
>
> 8. {% tool [Scanpy cluster, embed](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_cluster_reduce_dimension/scanpy_cluster_reduce_dimension/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy cluster, embed** {% icon tool %})
>    - *"Method used"*: `Embed the neighborhood graph using UMAP, using 'tl.umap'`
>
{: .hands_on}

Now let's take a look at our results. We'll plot one version of the UMAP showing the clusters we've just identified and another coloured by `Method` to see if that might be influencing our results.

> <tip-title>Interactive Visualisations</tip-title>
>
> * Select `Yes` for `Make an interactive plot?` if you want to explore the data further. You'll be able to colour the interactive plot by `Method` or any of the other metadata categories.
{: .tip}

> <hands-on-title> Visualise the Results </hands-on-title>
> 1. {% tool [Scanpy plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy cluster, embed** {% icon tool %})
>    - *"Method used for plotting"*: `Embeddings: Scatter plot in UMAP basis, using 'pl.umap'`
>        - *"Keys for annotations of observations/cells or variables/genes"*: `louvain,Method`
>        - *"Show edges?"*: `No`
>
{: .hands_on}

![Two UMAP plots showing many small and fragmented clusters of cells. The first UMAP is coloured into 47 clusters. The second UMAP shows many clusters as a single colour of cells analysed with the same method.](../../images/scrna_batch_correction/UMAP_before_Scanpy.png "UMAPs before batch correction integration coloured by A. louvain (cluster) and B. Method")

> <question-title></question-title>
>
> 1. How many clusters did we identify?
> 2. Are the batches well mixed?
>
> > <solution-title></solution-title>
> >
> > 1. The first plot is coloured by cluster. We can see there are 47 clusters (Scanpy numbers the clusters starting from 0). That's a lot - this is partly because we used a relatively high clustering resolution, but these fragmented clusters could also be a sign that something has gone wrong with our analysis.
> > 2. The second plot shows the UMAP coloured by `Method`. Each colour here represents cells that were sequenced by a different experimental technique. We can see lots of clusters and patches of cells that are only made up of one colour - this suggests that our cells are grouping together by batch rather than a biologically relevant characteristic such as cell type. This is a problem as it means we're not learning anything new about our cells since we already knew which batches they were in! If we want to find out something more interesting, we'll need to get rid of the technical differences between the batches.
> >
> {: .solution}
>
{: .question}

</div>

<div class='Seurat' markdown='1'>

In the Seurat pipeline, we will split each of our batches into its own 'layer' within the SeuratObject before we begin the analysis. We could do the same with each dataset if we were integrating multiple datasets together. Splitting will affect some aspects of preprocessing but it also sets up the dataset for the integration tools, which expect each batch or dataset to be in its own layer.

Splitting our data into layers means that the Seurat preprocessing tools can work on each layer separately. Seurat can treat each layer as if it were a separate dataset during preprocessing. Each layer (in this case, each of our batches) will be normalised independently. We'll also identify the highly variable genes within each batch, rather than across the whole dataset. Seurat will then create a single consensus list of highly variable genes to use for the whole dataset.

The other tools in the Seurat pipeline, such as `RunPCA` and `FindClusters` will still work on the entire dataset.

Splitting the batches into separate layers within our SeuratObject can act as a very mild form of batch correction. The separate preprocessing of each layer can reduce some of the technical differences between the batches. We have to wait for the results to see if this has been enough to eliminate these differences or if full batch correction is still needed.

> <hands-on-title> Split the Batches into Layers </hands-on-title>
>
> 1. {% tool [Seurat Integrate](toolshed.g2.bx.psu.edu/repos/iuc/seurat_integrate/seurat_integrate/5.0+galaxy0) %} with the following parameters:
>    - *"Method used"*: `Split data into layers using 'split'`
>        - *"Factor or group to use to split data"*: `Method`
>
>    > <comment-title></comment-title>
>    >
>    > We are splitting our data on `Method` as this is the column in our metadata that represents our batches. Each of the methods listed in this column will be split into its own layer.
>    {: .comment}
>
{: .hands_on}

Let's take a look to see what we've done to our data.

> <hands-on-title> Inspect the Data </hands-on-title>
>
> 2. {% tool [Seurat Data Management](toolshed.g2.bx.psu.edu/repos/iuc/seurat_data/seurat_data/5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file with the Seurat object"*: `rds_out` (output of **Seurat Integrate** {% icon tool %})
>    - *"Method used"*: `Inspect Seurat Object`
>        - *"Display information about"*: `General`
>
{: .hands_on}

> <question-title></question-title>
>
> 1. How many layers do we now have in our dataset?
> 2. What do these layers represent?
>
> > <solution-title></solution-title>
> >
> > 1. We can see that there are now 9 layers in our SeuratObject.
> > 2. We started out with one layer of raw data, called `counts`. That layer has now been split up according to `Method`. We now have nine `counts` layers. Each layer represents one of the batches named in the `Method` column of the cell metadata. We can see the names of the methods in the layer names. For example, the `counts.Drop-seq` layer contains the raw counts produced using the Drop-seq technique. Seven different methods were used in this study, but one of them was applied to three different batches - you should be able to see three layers with `Chromium_v2` in their names.
> >
> {: .solution}
>
{: .question}

Now that we have split our data so that each batch is in its own layer, we will cluster it. We won't perform any batch correction, so we'll see if the differences in `Method` are causing any problems that might require correction.

We'll follow the default Seurat pipeline here, except that we'll use `30` PCs to build the neighborhood graph and cluster with a resolution of `2` as these were the parameters used in [the original Seurat version of this tutorial](https://satijalab.org/seurat/articles/seurat5_integration). We'll also give our clusters and UMAP more recognisable names as we'll be running these tools again later, after batch correction.

> <comment-title></comment-title>
> Seurat has another option for preprocessing - rather than use the three separate functions we're using in this tutorial, you can use a single function called `SCTransform` to preform normalisation, identification of variable genes, and scaling all in one go. You will find this option on Galaxy's {% icon tool %} **Seurat Preprocessing** tool.
>
> If you use `SCTransform` for preprocessing, then you'll need to choose `Yes` for `Use SCT as Normalization Method` when you run `IntegrateLayers`. The `SCTransform` normalises the data in its own way, so we just need to let the tool know what to expect!
>
> The next step after identifying clusters would usually be to look for marker genes that are differentially expressed between clusters. If you perform integration/batch correction after using `SCTransform`, then you will need to run the `PrepSCTFindMarkers` function before using tools such as `FindMarkers`. You'll find this in the {% icon tool %} **Seurat Integrate** tool.
>
> The rest of the workflow will be the same as shown in this tutorial, but you will end up with slightly different results because `SCTransform` handles preprocessing in a different way than the three separate tools. If you want to learn more about these differences then you can follow the SCTransform route in the [Clustering 3k PBMCs with Seurat]({% link topics/single-cell/tutorials/scrna-seurat-pbmc3k/tutorial.html %}) tutorial.
{: .comment}

> <hands-on-title> Cluster without Batch Correction </hands-on-title>
>
> 1. {% tool [Seurat Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/seurat_preprocessing/seurat_preprocessing/5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file with the Seurat object"*: `rds_out` (output of **Seurat Integrate** {% icon tool %})
>    - *"Method used"*: `Normalize with 'NormalizeData'`
>        - *"Method for normalization"*: `LogNormalize`
>
> 2. {% tool [Seurat Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/seurat_preprocessing/seurat_preprocessing/5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file with the Seurat object"*: `rds_out` (output of **Seurat Preprocessing** {% icon tool %})
>    - *"Method used"*: `Identify highly variable genes with 'FindVariableFeatures'`
>        - *"Method to select variable features"*: `vst`
>        - *"Output list of most variable features"*: `No`
>
> 3. {% tool [Seurat Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/seurat_preprocessing/seurat_preprocessing/5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file with the Seurat object"*: `rds_out` (output of **Seurat Preprocessing** {% icon tool %})
>    - *"Method used"*: `Scale and regress with 'ScaleData'`
>        - *"Regress out a variable"*: `No`
>        - *"Features to scale"*: `Variable Features`
>
> 4. {% tool [Seurat Run Dimensional Reduction](toolshed.g2.bx.psu.edu/repos/iuc/seurat_reduce_dimension/seurat_reduce_dimension/5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file with the Seurat object"*: `rds_out` (output of **Seurat Preprocessing** {% icon tool %})
>    - *"Method used"*: `Run a PCA dimensionality reduction using 'RunPCA'`
>
>    > <comment-title></comment-title>
>    >
>    > We will use the output from `RunPCA` in the  following section when we perform batch correction.
>    >
>    > If you're already familiar with the Seurat clustering pipeline and you just want to try using the  {% icon tool %} **Seurat Integrate** tools, then you can skip ahead to the **Clustering after Integration** step now. In a real analysis, it would be best to complete the clustering without batch correction first to check if it is needed.
>    {: .comment}
>
> 5. {% tool [Seurat Find Clusters](toolshed.g2.bx.psu.edu/repos/iuc/seurat_clustering/seurat_clustering/5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file with the Seurat object"*: `rds_out` (output of **Seurat Run Dimensional Reduction** {% icon tool %})
>    - *"Method used"*: `Compute nearest neighbors with 'FindNeighbors'`
>        - *"Number of dimensions from reduction to use as input"*: `30`
>
> 6. {% tool [Seurat Find Clusters](toolshed.g2.bx.psu.edu/repos/iuc/seurat_clustering/seurat_clustering/5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file with the Seurat object"*: `rds_out` (output of **Seurat Find Clusters** {% icon tool %})
>    - *"Method used"*: `Identify cell clusters with 'FindClusters'`
>        - *"Resolution"*: `2.0`
>        - *"Algorithm for modularity optimization"*: `1. Original Louvain`
>        - *"Name for output clusters"*: `unintegrated_clusters`
>
>    > <warning-title></warning-title>
>    >
>    > Make sure that you change the default name for the clusters to `unintegrated_clusters`!
>    {: .warning}
>
> 7. {% tool [Seurat Run Dimensional Reduction](toolshed.g2.bx.psu.edu/repos/iuc/seurat_reduce_dimension/seurat_reduce_dimension/5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file with the Seurat object"*: `rds_out` (output of **Seurat Find Clusters** {% icon tool %})
>    - *"Method used"*: `Run a UMAP dimensional reduction using 'RunUMAP'`
>        - *"UMAP implementation to run"*: `uwot`
>        - *"Run UMAP on dimensions, features, graph or KNN output"*: `dims`
>            - *"Number of dimensions from reduction to use as input"*: `30`
>        - In *"Advanced Options"*:
>            - *"Name for dimensional reduction"*: `umap.unintegrated`
>
>    > <warning-title></warning-title>
>    >
>    > Make sure that you change the default name for the UMAP results to `umap.unintegrated`!
>    {: .warning}
>
{: .hands_on}

Now let's take a look at our results. We'll first plot a UMAP showing the clusters we've just identified. Then, we will colour this plot in by `Method` to see if that might be influencing our results.

> <hands-on-title> Visualise the Results </hands-on-title>
>
> 1. {% tool [Seurat Visualize](toolshed.g2.bx.psu.edu/repos/iuc/seurat_plot/seurat_plot/5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file with the Seurat object"*: `rds_out` (output of **Seurat Run Dimensional Reduction** {% icon tool %})
>    - *"Method used"*: `Visualize Dimensional Reduction with 'DimPlot'`
>        - *"Name of reduction to use"*: `umap.unintegrated`
>
> 2. {% tool [Seurat Visualize](toolshed.g2.bx.psu.edu/repos/iuc/seurat_plot/seurat_plot/5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file with the Seurat object"*: `rds_out` (output of **Seurat Run Dimensional Reduction** {% icon tool %})
>    - *"Method used"*: `Visualize Dimensional Reduction with 'DimPlot'`
>        - *"Name of reduction to use"*: `umap.unintegrated`
>        - In *"Advanced Options"*:
>            - *"Factor to group cells by"*: `Method`
>
{: .hands_on}

![Two UMAP plots showing many small and fragmented clusters of cells. Image A is coloured into 48 clusters. Image B shows many clusters as a single colour of cells analysed with the same method.](../../images/scrna_batch_correction/UMAP_Before_Seurat.png "UMAP before batch correction integration coloured by A: cluster, and B: Method")

> <question-title></question-title>
>
> 1. How many clusters did we identify?
> 2. Are the batches well mixed?
>
> > <solution-title></solution-title>
> >
> > 1. The first plot is coloured by cluster. We can see there are 48 clusters (Seurat numbers the clusters starting from 0). That's a lot - this is partly because we used a relatively high clustering resolution, but these fragmented clusters could also be a sign that something has gone wrong with our analysis.
> > 2. The second plot shows the UMAP coloured by `Method`. Each colour here represents cells that were sequenced by a different experimental technique. We can see lots of clusters and patches of cells that are only made up of one colour - this suggests that our cells are grouping together by batch rather than a biologically relevant characteristic such as cell type. This is a problem as it means we're not learning anything new about our cells since we already knew which batches they were in! If we want to find out something more interesting, we'll need to get rid of the technical differences between the batches.
> >
> {: .solution}
>
{: .question}

</div>

# Clustering with Batch Correction

It looks like we do need to perform batch correction on our dataset. The Scanpy and Seurat pipelines both provide tools that can reduce the technical differences between batches. If we were combining multiple datasets, we could use the same tools in the same way to perform integration (i.e. to correct for technical differences between the datasets).

<div class='Scanpy' markdown='1'>

Two simple changes will enable us to perform batch corrections within the Scanpy pipeline.

First, we will go back to the step where we identified Highly Variable Genes (HVGs). This time, we will add in `Method` as the batch key. Now, the tool will select the most variable genes within each batch before merging them into a shared list. Doing this can prevent selection of batch-specific genes and acts as a lightweight form of batch correction.

Secondly, we will add in one more step using a tool called Harmony. We will use Harmony in between performing the PCA and constructing the neighborhood graph. Harmony will take the principal components and adjust them to remove batch effects. It will create a corrected low-dimensional representation that we can use instead of the uncorrected PCA reduction. We will then use `X_pca_harmony` instead of the PCA when we construct the neighbourhood graph.

> <comment-title></comment-title>
>
> {% tool Scanpy remove confounders %} inclludes several methods for batch correction/integration, which all work in different ways. You might want to experiment by using different methods to see how they affect the results. When you are working on your own data, it can be a good idea to try a few different integration methods to see which one produces the best results. The best integration or batch correction would be the one that eliminates the most of the technical differences between datasets or batches while producing biologically meaningful results. If we end up with completely unexpected results rather than clusters that match up well with known cell types, then we know that something has gone wrong!
{: .comment}

> <hands-on-title> Recluster with Batch Correction </hands-on-title>
>
> 1. {% tool [Scanpy filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used for filtering"*: `Annotate (and filter) highly variable genes, using 'pp.highly_variable_genes'`
>        - *"Choose the flavor for identifying highly variable genes"*: `Seurat`
>        - *"Specify the batch key"*: `Method`
>
>    > <comment-title> short description </comment-title>
>    >
>    > We will specify `Method` as the batch key. This means that Highly Variable Genes will be identified within each batch and then combined.
>    {: .comment}
>
> 2. {% tool [Scanpy Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy filter** {% icon tool %})
>    - *"Method used for inspecting"*: `Scale data to unit variance and zero mean, using 'pp.scale'`
>        - *"Maximum value"*: `10.0`
>
> 3. {% tool [Scanpy cluster, embed](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_cluster_reduce_dimension/scanpy_cluster_reduce_dimension/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used"*: `Computes PCA (principal component analysis) coordinates, loadings and variance decomposition, using 'pp.pca'`
>        - *"Type of PCA?"*: `Full PCA`
>
> 4. {% tool [Scanpy remove confounders](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_remove_confounders/scanpy_remove_confounders/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy cluster, embed** {% icon tool %})
>    - *"Method used for plotting"*: `Integrate multiple single-cell experiments with Harmony, using 'external.pp.harmony_integrate'`
>        - *"The name of the column in adata.obs that differentiates among experiments/batches"*: `Method`
>
> 5. {% tool [Scanpy Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy remove confounders** {% icon tool %})
>    - *"Method used for inspecting"*: `Compute a neighborhood graph of observations, using 'pp.neighbors'`
>        - *"Number of PCs to use"*: `30`
>        - *"Use the indicated representation"*: `X_pca_harmony`
>
>    > <comment-title></comment-title>
>    >
>    > We will use the corrected embedding, 'X_pca_harmony' to calculate the neighborhood graph. Make sure to enter this as the representation to use.
>    {: .comment}
>
> 6. {% tool [Scanpy cluster, embed](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_cluster_reduce_dimension/scanpy_cluster_reduce_dimension/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used"*: `Cluster cells into subgroups, using 'tl.louvain'`
>        - *"Flavor for the clustering"*: `vtraag (much more powerful than igraph)`
>            - *"Resolution"*: `2.0`
>        - *"Key under which to add the cluster labels"*: `louvain_integrated`
>
>    > <comment-title></comment-title>
>    >
>    > We'll use a different name for this clustering so that we don't get confused. Enter 'louvain_integrated' as the key to add the cluster labels under. We'll use this when we plot our integrated clusters.
>    {: .comment}
>
> 7. {% tool [Scanpy cluster, embed](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_cluster_reduce_dimension/scanpy_cluster_reduce_dimension/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy cluster, embed** {% icon tool %})
>    - *"Method used"*: `Embed the neighborhood graph using UMAP, using 'tl.umap'`
>
{: .hands_on}

Now let's visualise the results again to see if the batch correction has worked. As before, we'll make one plot coloured by cluster and another coloured by batch (`Method`). We're hoping that the batches in that second plot will be more mixed together instead of forming separate groups like they did before batch correction.

> <tip-title>Interactive Visualisations</tip-title>
>
> * Select `Yes` for `Make an interactive plot?` if you want to explore the data further. You'll be able to colour the interactive plot by `Method` or any of the other metadata categories.
{: .tip}

> <hands-on-title> Visualise the Results </hands-on-title>
> 1. {% tool [Scanpy plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy cluster, embed** {% icon tool %})
>    - *"Method used for plotting"*: `Embeddings: Scatter plot in UMAP basis, using 'pl.umap'`
>        - *"Keys for annotations of observations/cells or variables/genes"*: `louvain_integrated,Method`
>        - *"Show edges?"*: `No`
>
>    > <comment-title> short description </comment-title>
>    >
>    > Make sure to use 'louvain_integrated' rather than 'louvain' for the cluster annotation to plot. This is the name we used for our integrated clusters.
>    {: .comment}
>
{: .hands_on}

![Two UMAP plots showing three large groups of cells, each made up of multiple clusters. The first UMAP is coloured as 25 clusters. The second UMAP shows a mix of colours across all the clusters.](../../images/scrna_batch_correction/UMAP_after_Scanpy.png "UMAP after batch correction coloured by A. louvain_integrated (cluster) and B. Method")

> <question-title></question-title>
>
> 1. How many clusters did we identify?
> 2. How well mixed are the batches?
>
> > <solution-title></solution-title>
> >
> > 1. The first plot shows 25 clusters (remember that Scanpy starts from cluster 0!). Although the high resolution means we still have plenty of clusters, the batch correction has reduced the number. The clusters also look less fragmented than they did before. The reduced number of clusters doesn't necessarily mean the analysis is better, but removing the batch-specific structure has allowed biologically similar cells from different batches to cluster together.
> > 2. When we colour in the plot by `Method`, we can see that all the colours are mixed together across all of the clusters. We don't have any clusters that contain only one colour. The batch correction has successfully removed the differences between the batches so that they're no longer dominating the results.
> >
> {: .solution}
>
{: .question}

</div>

<div class='Seurat' markdown='1'>

We will now run Seurat's batch correction tool - it's called `IntegrateLayers`, but despite the name we can use the same tool to address differences between batches as we would for integrating datasets.

> <comment-title></comment-title>
>
> {% tool Seurat Integrate %} provides several integration methods, which all perform the integration or batch correction in their own way. You might want to experiment by using different methods to see how they affect the results. When you are working on your own data, it can be a good idea to try a few different integration methods to see which one produces the best results. The best integration or batch correction would be the one that eliminates the most of the technical differences between datasets or batches while producing biologically meaningful results. If we end up with completely unexpected results rather than clusters that match up well with known cell types, then we know that something has gone wrong!
{: .comment}

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Seurat Integrate](toolshed.g2.bx.psu.edu/repos/iuc/seurat_integrate/seurat_integrate/5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file with the Seurat object"*: `rds_out` (output of **Seurat Run Dimensional Reduction** {% icon tool %})
>    - *"Method used"*: `Apply integration methods with 'IntegrateLayers'`
>        - *"Integration method to use"*: `CCA Integration`
>        - *"Name for new dimensional reduction"*: `integrated.cca`
>
>    > <comment-title> Remember the name </comment-title>
>    >
>    > Make sure you remember the name you've used for the new dimensional reduction - we'll be using this later instead of the PCA we produced previously.
>    {: .comment}
>
{: .hands_on}

It's good practice to rejoin our layers now, so that those separate layers/batches will end up together. We don't actually need to do this now (as it won't affect the clustering results), but it is important if we want to perform downstream analyses such as Differential Expression analysis.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Seurat Integrate](toolshed.g2.bx.psu.edu/repos/iuc/seurat_integrate/seurat_integrate/5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file with the Seurat object"*: `rds_out` (output of **Seurat Integrate** {% icon tool %})
>    - *"Method used"*: `Join layers with 'JoinLayers'`
>
{: .hands_on}

> <question-title></question-title>
>
> 1. How many layers are now in our dataset?
>
> > <solution-title></solution-title>
> >
> > 1. You might think that we should only have one layer in our dataset now, because we split it into nine layers that we have now rejoined. However, if you use {% tool Seurat Data Management %} to check, you'll see that we actually have three layers now! This is because the preprocessing functions we ran (normalisation and scaling) created their own layers of data. We still have the original raw `counts` layer, but we now have a normalised layer called `data` and a scaled one called `scale.data` as well.
> >
> >    In fact, if you run {% tool Seurat Data Management %} on the previous dataset in your history, from before we rejoined the layers, you'll see that it actually had 19 layers in it - each of the nine `counts` layers we split the dataset was normalised into its own `data` layer. We then had the `scale.data` layer too.
> >
> {: .solution}
>
{: .question}

Now let's try clustering our integrated data. We'll repeat the steps we performed earlier, but this time we'll be using the dimensional reduction produced by `Integrate Layers` instead of the PCA. The clustering will be based on the integrated embedding rather than the original PCA embedding. Let's also give our clusters and UMAP results some new names to differentiate them from the uncorrected results.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Seurat Find Clusters](toolshed.g2.bx.psu.edu/repos/iuc/seurat_clustering/seurat_clustering/5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file with the Seurat object"*: `rds_out` (output of **Seurat Integrate** {% icon tool %})
>    - *"Method used"*: `Compute nearest neighbors with 'FindNeighbors'`
>        - *"Name of reduction to use"*: `integrated.cca`
>        - *"Number of dimensions from reduction to use as input"*: `30`
>
>    > <comment-title></comment-title>
>    >
>    > Make sure to use `integrated.cca` as the reduction, not the `pca` we made previously.
>    {: .comment}
>
> 2. {% tool [Seurat Find Clusters](toolshed.g2.bx.psu.edu/repos/iuc/seurat_clustering/seurat_clustering/5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file with the Seurat object"*: `rds_out` (output of **Seurat Find Clusters** {% icon tool %})
>    - *"Method used"*: `Identify cell clusters with 'FindClusters'`
>        - *"Resolution"*: `2.0`
>        - *"Algorithm for modularity optimization"*: `1. Original Louvain`
>        - *"Name for output clusters"*: `cca_clusters`
>
>    > <comment-title></comment-title>
>    >
>    > Make sure that you know what name you used for your clusters as we'll use this for the UMAP!
>    {: .comment}
>
> 3. {% tool [Seurat Run Dimensional Reduction](toolshed.g2.bx.psu.edu/repos/iuc/seurat_reduce_dimension/seurat_reduce_dimension/5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file with the Seurat object"*: `rds_out` (output of **Seurat Find Clusters** {% icon tool %})
>    - *"Method used"*: `Run a UMAP dimensional reduction using 'RunUMAP'`
>        - *"Name of reduction to use"*: `integrated.cca`
>        - *"UMAP implementation to run"*: `uwot`
>        - *"Run UMAP on dimensions, features, graph or KNN output"*: `dims`
>            - *"Number of dimensions from reduction to use as input"*: `30`
>        - In *"Advanced Options"*:
>            - *"Name for dimensional reduction"*: `umap.cca`
>
>    > <comment-title></comment-title>
>    >
>    > Make sure that you know what name you used for your UMAP results as we'll use this for the plots!
>    {: .comment}
>
{: .hands_on}

Let's see how the batch correction has changed our results. As before, we'll make one plot coloured by cluster and then another coloured by batch (`Method`). We're hoping that the batches in that second plot will be more mixed together instead of forming separate groups like they did before batch correction.

> <hands-on-title> Visualise the Results </hands-on-title>
>
> 1. {% tool [Seurat Visualize](toolshed.g2.bx.psu.edu/repos/iuc/seurat_plot/seurat_plot/5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file with the Seurat object"*: `rds_out` (output of **Seurat Run Dimensional Reduction** {% icon tool %})
>    - *"Method used"*: `Visualize Dimensional Reduction with 'DimPlot'`
>        - *"Name of reduction to use"*: `umap.cca`
>
> 2. {% tool [Seurat Visualize](toolshed.g2.bx.psu.edu/repos/iuc/seurat_plot/seurat_plot/5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file with the Seurat object"*: `rds_out` (output of **Seurat Run Dimensional Reduction** {% icon tool %})
>    - *"Method used"*: `Visualize Dimensional Reduction with 'DimPlot'`
>        - *"Name of reduction to use"*: `umap.cca`
>        - In *"Advanced Options"*:
>            - *"Factor to group cells by"*: `Method`
>
{: .hands_on}

![Two UMAP plots showing three large groups of cells, each made up of multiple clusters. Image A is coloured as 25 clusters. Image B shows a mix of colours across all the clusters.](../../images/scrna_batch_correction/UMAP_After_Seurat.png "UMAP after batch correction coloured by A. cluster B. Method")

> <question-title></question-title>
>
> 1. How many clusters did we identify?
> 2. How well mixed are the batches?
>
> > <solution-title></solution-title>
> >
> > 1. The first plot shows 25 clusters (remember that Seurat starts from cluster 0!). Although the high resolution means we still have plenty of clusters, the batch correction has reduced the number. The clusters also look less fragmented than they did before. The reduced number of clusters doesn't necessarily mean the analysis is better, but removing the batch-specific structure has allowed biologically similar cells from different batches to cluster together.
> > 2. When we colour in the plot by `Method`, we can see that all the colours are mixed together across all of the clusters. We don't have any clusters that contain only one colour. The batch correction has successfully removed the differences between the batches so that they're no longer dominating the results.
> >
> {: .solution}
>
{: .question}

</div>

# Comparing the Results

Let's take another look at our UMAPs coloured by `Method` to see what the batch correction process has done to our data. In the 'before' picture, we can see that the different batches are forming their own clusters or patches in the UMAP plot, with very little mixing between colours. The differences between batches or methods are having a big impact on the clustering, which means that the biological differences we're interested in are being missed. In the 'after' picture, we can see that the colours are all mixed up and the clusters are no longer separating out based on method. We have removed the technical differences between batches, so hopefully these clusters are now based on the biological differences we're interested in.

<div class='Scanpy' markdown='1'>

![Two UMAP plots showing three large groups of cells, each made up of multiple clusters. Each UMAP is coloured in based on Method. The first shows patches and clusters that are all the same colour or method, with little mixing between colours. The second UMAP shows the colours are well mixed and all clusters include a mixture of multiple colours](../../images/scrna_batch_correction/Scanpy_before_after.png "UMAP coloured by Method A. before and B. after batch correction")

</div>

<div class='Seurat' markdown='1'>

![Two UMAP plots showing three large groups of cells, each made up of multiple clusters. Each UMAP is coloured in based on Method. The first shows patches and clusters that are all the same colour or method, with little mixing between colours. The second UMAP shows the colours are well mixed and all clusters include a mixture of multiple colours](../../images/scrna_batch_correction/Seurat_before_after.png "UMAP coloured by Method A. before and B. after batch correction")

</div>

# Checking the Clusters are Biologically Meaningful

Our plots suggest that the batch correction has successfully brought the different methods together, but this alone is not enough to confirm that it has worked. As always in single cell analysis, we also need to verify that the clusters we have found are biologically meaningful. Scanpy and Seurat will always produce clusters, but it is up to us to evaluate whether those results actually make sense.

In order to do this, we would usually take a closer look at the clusters to work out what they represent, for example by looking for clusters expressing genes that are known to be present in specific cell types. If you've worked through the [Scanpy]({% link topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.md %}) or [Seurat]({% link topics/single-cell/tutorials/scrna-seurat-pbmc3k/tutorial.md %}) clustering tutorials then you'll already have seen how this can be done using the top differentially expressed genes or known markers of gene types. If you haven't already completed these tutorials then they can tell you more about identifying cell types.

We don't need to go through this process again now, because we have the annotations provided by the researchers who created this dataset. If you look back at the cell metadata table we created at the beginning of this tutorial, you'll see there is an annotation called `CellType`. We can colour in our UMAPs using this annotation instead of the `Method`. If our clusters make biological sense, we should see that these cell types are clumped together because cells of the same type should be close to each other. If the cell types are all blended together across the entire UMAP (as with the methods in our integrated plots) then this would be a sign that something has gone wrong - we want the different methods to be mixed up together, but we'd like the biologically meaningful differences between cell types to be preserved. When we are performing batch correction or integration, there is a risk that we could over-integrate the data, eliminating the biological differences we're interested in alongside the technical differences we wanted to remove.

You can rerun the UMAP plots yourself if you like, or just take a look at the plots below to see how the integration has grouped together the cells in a biologically meaningful way. The `CellType` annotation won't match up exactly with our clusters (remember we used a high resolution to make lots of clusters!) but they certainly shouldn't be scattered across the whole plot.

<div class='Scanpy' markdown='1'>

![Two UMAP plots showing three large groups of cells, each made up of multiple clusters. Each UMAP is coloured in based on cell types. The first shows patches of the same cell types appearing in different parts of the plot. The second UMAP shows each colour or cell type is centred in just one cluster or area of the plot](../../images/scrna_batch_correction/CellType_Scanpy.png "UMAP coloured by cell type A. before and B. after batch correction")

</div>

<div class='Seurat' markdown='1'>

![Two UMAP plots showing three large groups of cells, each made up of multiple clusters. Each UMAP is coloured in based on cell types. The first shows patches of the same cell types appearing in different parts of the plot. The second UMAP shows each colour or cell type is centred in just one cluster or area of the plot](../../images/scrna_batch_correction/CellType_Seurat.png "UMAP coloured by cell type A. before and B. after batch correction")

</div>

# Conclusion

{% icon congratulations %} Well done, you've successfully performed batch correction to remove technical effects while clustering single cell data with Scanpy or Seurat. You might want to check your results against the example histories for the [Scanpy](https://singlecell.usegalaxy.eu/u/marisa_jl/h/batch-correction-integration-with-scanpy-answer-key) or [Seurat](https://singlecell.usegalaxy.eu/u/marisa_jl/h/batch-correction--integration-with-seurat---answer-key) pipelines. You can also take a look at the whole workflow for [Scanpy](https://singlecell.usegalaxy.eu/u/marisa_jl/w/batch-correction-scanpy) or [Seurat](https://singlecell.usegalaxy.eu/u/marisa_jl/w/batch-correction---seurat).

In this tutorial, we've learned how to perform batch correction or integration when analysing single cell data with either the Scanpy or Seurat pipelines. If you want to learn more about these pipelines then you might want to try analysing a slightly trickier dataset in the [Scanpy]({% link topics/single-cell/tutorials/scrna-case_basic-pipeline/tutorial.md %}) or [Seurat]({% link topics/single-cell/tutorials/scrna-case_FilterPlotandExplore_SeuratTools/tutorial.md %}) case study tutorials.
