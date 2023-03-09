# Introduction


You've done all the hard work of: preparing a single cell matrix, processing it, plotting it, interpreting it, and finding lots of lovely genes. Now you want to infer trajectories, or relationships between cells... and you've been threatened with learning Python to do so! Well, fear not. If you can have a run-through of a basic python coding introduction such as [this one](https://www.w3schools.com/python/), then that will help you make more sense of this tutorial, however you'll be able to make and interpret glorious plots even without understanding the Python coding language. This is the beauty of Galaxy - all the 'set-up' is identical across computers, because it's browser based. So fear not!

Traditionally, we thought that differentiating or changing cells jumped between discrete states, so 'Cell A' became 'Cell B' as part of its maturation. However, most data shows otherwise. Generally, there is a spectrum (a 'trajectory', if you will...) of small, subtle changes along a pathway of that differentiation. Trying to analyse cells every 10 seconds can be pretty tricky, so 'pseudotime' analysis takes a single sample and assumes that those cells are all on slightly different points along a path of differentiation. Some cells might be slightly more mature and others slightly less, all captured at the same 'time'. We 'assume' or 'infer' relationships between cells.

We will use the same sample from the previous three tutorials, which contains largely T-cells in the thymus. We know T-cells differentiate in the thymus, so we would assume that we would capture cells at slightly different time points within the same sample. Furthermore, our cluster analysis alone showed different states of T-cell. Now it's time to look further!

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Get data

We've provided you with experimental data to analyse from a mouse dataset of fetal growth restriction {% cite Bacon2018 %}. This is the full dataset generated from [this tutorial]({% link topics/single-cell/tutorials/scrna-case_basic-pipeline/tutorial.md %}) (see the [study in Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-6945/results/tsne) and the [project submission](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6945/)). You can find the final dataset in this [input history](https://usegalaxy.eu/u/wendi.bacon.training/h/cs4inferred-trajectory-analysis-using-python-jupyter-notebook-in-galaxy---input) or download from Zenodo below.

> <hands-on-title>Option 1: Data upload - Import history</hands-on-title>
>
> 1. Import history from: [input history](https://usegalaxy.eu/u/wendi.bacon.training/h/cs4inferred-trajectory-analysis-using-python-jupyter-notebook-in-galaxy---input)
>
>
>    {% snippet faqs/galaxy/histories_import.md %}
>
> 2. **Rename** {% icon galaxy-pencil %} the the history to your name of choice.
>
{: .hands_on}

> <hands-on-title>Option 2: Data upload - Add to history</hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the AnnData object from [Zenodo]({{ page.zenodo_link }})
>
>    ```
>    {{ page.zenodo_link }}/files/Final_cell_annotated_object.h5ad
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. **Rename** {% icon galaxy-pencil %} the .h5ad object as `Final cell annotated object`
>
>    {% snippet faqs/galaxy/datasets_rename.md name="Final cell annotated object" %}
>
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

## Filtering for T-cells

One problem with our current dataset is that it's not just T-cells: we found in the previous tutorial that it also contains macrophages. This is a problem, because trajectory analysis will generally try to find relationships between all the cells in the sample. We need to remove those cell types to analyse the trajectory.

> <hands-on-title>Removing macrophages</hands-on-title>
>
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Final cell annotated object`
>    - *"Function to manipulate the object"*: `Filter observations or variables`
>    - *"What to filter?"*: `Observations (obs)`
>    - *"Type of filtering?"*: `By key (column) values`
>    - *"Key to filter"*: `cell_type`
>    - *"Type of value to filter"*: `Text`
>    - *"Filter"*: `not equal to`
>    - *"Value"*: `Macrophages`
>
> 3. **Rename** {% icon galaxy-pencil %} output h5ad `T-cell_object.h5ad`
>
{: .hands_on}

You should now have `8569` cells, as opposed to the `8605` you started with. You've only removed a few cells (the contaminants!), but it makes a big difference in the next steps.

Take note of what `#` this dataset is in your history, as you will need that shortly!

## Launching Jupyter

> <warning-title>Data uploads and Jupyter</warning-title>
> There are a few ways of importing and uploading data in Jupyter. You might find yourself accidentally doing this differently than the tutorial, and that's ok. There are a few key steps where you will call files from a location - if these don't work from you, check that the file location is correct and change accordingly!
{: .warning}

JupyterLab is a bit like RStudio but for other coding languages. What, you've never heard of [RStudio](https://www.rstudio.com/products/rstudio/features/)? Then don't worry, just follow the instructions!

{% icon warning %} Please note: this is only currently available on the [usegalaxy.eu](https://usegalaxy.eu) and [usegalaxy.org](https://usegalaxy.org) sites.

> <hands-on-title>Launching JupyterLab</hands-on-title>
>
> 1. {% tool [Interactive JupyTool and Notebook](interactive_tool_jupyter_notebook) %} with the following parameters:
>    - *"Do you already have a notebook?"*: `Start with a fresh notebook`
>
>    This may take a moment, but once the `Executed notebook` in your dataset is **orange**, you are up and running!
>
> 2. Either click on the blue `User menu` that appeared when you executed the **Interactive JupyTool and Notebook**, or go to the top of the screen and choose `User` and then `Active InteractiveTools`.
>
> 3. Click on the newest `JupyTool interactive tool`.
>
{: .hands_on}

Welcome!

> <warning-title>Danger: You can lose data!</warning-title>
> Do NOT delete or close this notebook dataset in your history. YOU WILL LOSE IT!
{: .warning}

You have two options for how to proceed with this tutorial - either you download the tutorial notebook and run the lines in that notebook, or you can copy and paste the code for each step into a fresh notebook and run it yourself. The initial instructions for both options are below.

> <hands-on-title>Option 1A: Open the notebook directly in the JupyterLab</hands-on-title>
>
> 1. Open a Terminal in JupyterLab with File -> New -> Terminal
>
> 2. Run 
>    ```
>    wget {{ ipynbpath }}
>    ```
>
> 3. Select the notebook that appears in the list of files on the left.
>
{: .hands_on}

> <hands-on-title>Option 1B: Downloading the tutorial notebook & Upload</hands-on-title>
>
> 1. You will need to download the tutorial notebook locally to your own computer. Do this by clicking on {% icon notebook %}`Jupyter notebook` in the `Supporting Materials` section at the very beginning of the tutorial, in the Overview box.
>
> 2. In the folder window, {% icon galaxy-upload %} Upload the downloaded notebook from your computer. It should appear in the file window.
>
> 3. Open it by double clicking it in the file window.
>
{: .hands_on}

> <hands-on-title>Option 2: Creating a notebook</hands-on-title>
>
> 1. Click the **Python 3** icon under **Notebook**
>
>   ![Python 3 icon](../../images/scrna-casestudy/wab-python3logo.png "Python 3 Button")
>
> 2. Save your file (**File**: **Save**, or click the {% icon galaxy-save %} Save icon at the top left)
>
> 3. If you right click on the file in the folder window at the left, you can rename your file `whateveryoulike.ipynb`
>
{: .hands_on}

> <warning-title>You should <b>Save</b> frequently!</warning-title>
> This is both for good practice and to protect you in case you accidentally close the browser. Your environment will still run, so it will contain the last saved notebook you have. You might eventually stop your environment after this tutorial, but ONLY once you have saved and exported your notebook (more on that at the end!) Note that you can have multiple notebooks going at the same time within this JupyterLab, so if you do, you will need to save and export each individual notebook. You can also download them at any time.
{: .warning}
