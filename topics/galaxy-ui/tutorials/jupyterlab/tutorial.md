---
layout: tutorial_hands_on

title: JupyterLab in Galaxy
zenodo_link: ""
questions:
- How can I manipulate data using JupyterLab in Galaxy?
- How can I create a new Python or R Jupyter notebook in JupyterLab?
- How can I import/export dataset from/to my history to/from the notebook?
- How can I save my notebook to my history?
objectives:
- Launch JupyterLab in Galaxy
- Start a Python or R Jupyter notebook
- "Install Libraries with pip or Conda"
- "Use get() to import datasets from your history to the notebook"
- "Use put() to export datasets from the notebook to your history"
- "Save your notebook into your history"
follow_up_training:
-
    type: "internal"
    topic_name: introduction
    tutorials:
        - r-basics
time_estimation: 3H
key_points:
- Why it's helpful to be able to work with JupyterLab interactively within Galaxy
contributors:
  - annefou
---


# Introduction
{:.no_toc}

{% include topics/galaxy-ui/tutorials/jupyterlab/tutorial_origin.md %}

[JupyterLab](https://jupyterlab.readthedocs.io/en/stable) is an [Integrated Development Environment (IDE)](https://en.wikipedia.org/wiki/Integrated_development_environment). 
Like most IDEs, it provides a graphical interface for R/Python, making it more user-friendly, and providing dozens of useful features. 
We will introduce additional benefits of using JupyterLab as you cover the lessons. 

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# JupyterLab

## Opening up JupyterLab

{% include snippets/launch_jupyterlab.md %}

You should now be looking at a page with the JupyterLab interface:

![rstudio default session](../../images/jupyterlab/jupyterlab_session_default.png)

## Creating your first R script

Now that we are ready to start exploring R, we will want to keep a record of the commands we are using. To do this we can create an R script.

> ### {% icon hands_on %} Hands-on: Create a R notebook
>
> 1. Click on the  **R icon**  in the **Notebook** section
{: .hands_on}

A new R notebook appears in the centre panel. Before we go any further, you should save your script.

> ### {% icon hands_on %} Hands-on: Save a R notebook
>
> 1. Click the **File** menu and select **Save Notebook As...
>    Alternatively, you can also:
>    - Click the {% icon galaxy-save %} icon (**Save the notebook contents and create checkpoint**) in the bar above the first line in the script editor
>    - Click the **File** menu and select **Save Notebook**
>    - Type <kbd>CTRL</kbd>+<kbd>S</kbd> (<kbd>CMD</kbd>+<kbd>S</kbd> on OSX)
>
> 2. In the **Save Notebook As** window that opens, name your file `jupyterlab_r_basics`
>    Alternatively, you can also rename your Jupyter Notbook:
>    - Right click on the name (`Untitled.ipynb`) in the bar above the first line in the script editor and select **Rename Notebook**
{: .hands_on}

The new script `jupyterlab_r_basics.ipynb` should appear in the **File Browser*** in the left panel. By convention, Jupyter notebooks  end with the file extension `.ipynb` independently of the programming language (R, Python, etc.).

## Overview and customization of the JupyterLab layout

Here are the major windows (or panels) of the JupyterLab environment:

![jupyterlab default session](../../images/jupyterlab/jupyterlab_session_default_layout.png)


> ### {% icon comment %} Working with the terminal
> Although we won't be working with R or Pyrhon at the terminal, there are lots of reasons to.
>
>
> For more on running an R or Python Script at the terminal see the dedicated [Software Carpentry lesson](https://swcarpentry.github.io/r-novice-inflammation/05-cmdline/).
{: .comment}

# Interaction between JupyterLab and Galaxy

Getting data in and out from Galaxy


# Conclusion
{:.no_toc}
