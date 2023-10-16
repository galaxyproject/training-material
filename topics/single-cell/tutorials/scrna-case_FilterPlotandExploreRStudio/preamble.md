# Introduction

You’ve previously done all the work to make a single cell matrix. Now it’s time to fully process our data using Seurat: remove low quality cells, reduce the many dimensions of data that make it difficult to work with, and ultimately try to define clusters and find some biological meaning and insights! There are many packages for analysing single cell data - Seurat ({% cite Satija2015 %}), Scanpy ({% cite Wolf2018 %}), Monocle ({% cite Trapnell2014 %}), Scater ({% cite McCarthy2017 %}), and many more. We’re working with Seurat in RStudio because it is well updated, broadly used, and highly trusted within the field of bioinformatics.

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

# Get Data onto Galaxy
To start, let's get our dataset loaded into Galaxy.

You can access the data for this tutorial in multiple ways:
1. **EBI Data Retrieval** - You may retrieve that files necessary to construct a Seurat Object in this way.Doing to will alleviate the necessity to convert AnnData (Python) objects into Seurat (R) objects:

> <hands-on-title>GetData</hands-on-title>
>
> Run{% tool [EBI SCXA Data Retrieval](toolshed.g2.bx.psu.edu/repos/ebi-gxa/retrieve_scxa/retrieve_scxa/v0.0.2+galaxy2) %} with the following parameters:
> - *"SC-Atlas experiment accession"*: `E-MTAB-6945`
> - *"Choose the type of matrix to download"*: `Raw filtered counts`
{: .hands_on}

2. **Importing from a history** - You can import [this history](https://usegalaxy.eu/u/camila-goclowski/h/fpe)

   {% snippet faqs/galaxy/histories_import.md %}
This also alleviates the necessity to convert the AnnData object into a Seurat one.

3. **Uploading from Zenodo** (see below)

> <hands-on-title>Option 3: Uploading from Zenodo</hands-on-title>
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
> 4. Check that the datatype is `h5ad`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="h5ad" %}
>
{: .hands_on}

# Important tips for easier analysis

{% snippet faqs/galaxy/tutorial_mode.md %}

> <comment-title></comment-title>
> - The Galaxy tool search panel sometimes doesn't find the tools we need from the thousands available.
> - You'll have a much easier time selecting tools from the panel (if you aren't using tutorial mode!) if you are on the [https://humancellatlas.usegalaxy.eu](https://humancellatlas.usegalaxy.eu)
{: .comment}

# Open RStudio in Galaxy
You now should have imported the matrix.mtx, genes.tsv, barcodes.tsv, and exp_design.tsv files into your Galaxy history. For the rest of the workflow, let's move onto RStudio and get coding!
> <hands-on-title>Open RStudio in Galaxy</hands-on-title>
> Run {% tool [RStudio](interactive_tool_rstudio)%}
{: .hands_on}


><comment-title>Next Step</comment-title>
> The interactive RStudio tool should begin to load now. Make your way over to your Active Interactive Tools page (User (in the top bar of the usegalaxy page)> Active Interactive Tools > RStudio)
>
>Alternatively, you may use the view (eye) button in your Galaxy History to open the interactive RStudio environment.
{: .comment}
