# Introduction

You’ve done all the work to make a single cell matrix, with gene counts and mitochondrial counts and buckets of cell metadata from all your variables of interest. Now it’s time to fully process our data, to remove low quality cells, to reduce the many dimensions of data that make it difficult to work with, and ultimately to try to define our clusters and to find our biological meaning and insights! There are many packages for analysing single cell data - Seurat Satija et al. 2015, Scanpy Wolf et al. 2018, Monocle Trapnell et al. 2014, Scater McCarthy et al. 2017, and so forth. We’re working with Scanpy, the python iteration of the most widely used single cell toolkit.

> <details-title>Python Version</details-title>
>
> This tutorial is an adaptation of [Filter, Plot and Explore]({% link topics/single-cell/tutorials/scrna-case_basic-pipeline/tutorial.md %}). The workflow has been converted into a Jupyter notebook that can be ran in Galaxy through `JupyterLab`. The notebook runs in Python and primarily relies on the Scanpy library for performing most tasks. Running through this notebook will allow you to see and modify the code being run at each step, so feel free to experiment with the different code cells to gain a deeper understanding of the analysis process.
>
{: .details}


## Get data

We've provided you with experimental data to analyse from a mouse dataset of fetal growth restriction {% cite Bacon2018 %}. This is the full dataset generated from [this tutorial]({% link topics/single-cell/tutorials/scrna-case_alevin-combine-datasets/tutorial.md %}) if you used the full FASTQ files rather than the subsampled ones (see the [study in Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-6945/results/tsne) and the [project submission](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6945/)). You can find this dataset in this [input history](https://usegalaxy.eu/u/wendi.bacon.training/h/cs3-answerkey) or download from Zenodo below.

You can access the data for this tutorial in multiple ways:

1. **Your own history** - If you're feeling confident that you successfully ran a workflow on all 7 samples from the previous tutorial, and that your resulting 7 AnnData objects look right (you can compare with the [answer key history](https://usegalaxy.eu/u/wendi.bacon.training/h/cs2combining-datasets-after-pre-processing---input-1)), then you can use those! To avoid a million-line history, I recommend dragging the resultant datasets into a fresh history

   {% snippet faqs/galaxy/histories_copy_dataset.md %}

2. **Importing from a history** - You can import [this history](https://usegalaxy.eu/u/wendi.bacon.training/h/cs3-answerkey)

   {% snippet faqs/galaxy/histories_import.md %}

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

## Launching JupyterLab

Notebook's can be run in Galaxy using the JupyterLab tool. JupyterLab is an interactive development environment that allows users to create interactive workflows combining code, text descriptions, and visualisations.

> {% snippet faqs/galaxy/interactive_tools_jupyter_launch.md %}

> <warning-title>Save your data!</warning-title>
> Data saved in interactive notebooks only persist with the instance currently running, stopping JupyterLab will result in loosing any stored data. Regularly download any important files/data!
{: .warning}

## Prepare the notebook

You have two options for running the code in this tutorial, you can either download a completed version of the notebook from the overview box at the start of this tutorial or follow along and create the notebook yourself as you go.

> <hands-on-title>Option 1: Importing notebook from file</hands-on-title>
>
> 1. Download the notebook from the {% icon external-link %} **Supporting Materials** section in the overview box
>
> 2. Click the {% icon galaxy-upload %} **Upload Files** button in JupyterLab
>
> 3. Select the downloaded notebook ```filter_plot_and_explore.ipynb```
>
> 4. The notebook should appear on the left hand side, click on the file to open it (if prompted to select a kernel select ```Python```)
>
{: .hands_on}

> <hands-on-title>Option 2: Creating a new notebook</hands-on-title>
>
> 1. Under the **Notebook** section in the JupyterLab select ```Python 3```
>
> 2. Rename the new notebook from the left-hand side panel
>
{: .hands_on}

{% icon congratulations %} Congratulations! you now have JupyterLab running and can start the tutorial!
