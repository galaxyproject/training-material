---
layout: tutorial_hands_on
title: Functionally Assembled Terrestrial Ecosystem Simulator (FATES) with Galaxy Climate JupyterLab
zenodo_link: ''
requirements:
  -
    type: "internal"
    topic_name: galaxy-ui
    tutorials:
        - jupyterlab
  -
    type: "external"
    title: "Programming with Python"
    link: "https://swcarpentry.github.io/python-novice-inflammation/"

questions:
- What is the Community Terrestrial System Model (CTSM)?
- What is FATES?
- Why do we need a "host" Land Surface Model for running FATES?
- Why running single point CLM-FATES simulations?
- Which configurations of CLM-FATES can I use?
- What type of data is available to compare the model outputs with observations?
objectives:
- Understanding ecological, biogeochemical, biogeophysical, and hydrologic theory underpinning the CTSM.
- Running CLM-FATES in Galaxy for single-point locations where in-situ measurements are available.
- Analyzing CLM-FATES output and comparing it with observations.
- Sharing data (results from the simulations, in-situ data) using FAIR principles.
- Composing, executing and publishing the corresponding Galaxy workflow. 
time_estimation: 12H
key_points:
- CTSM
- FATES
- Model analysis
- Comparison with observations.
contributors:
- annefou

---


# Introduction
{:.no_toc}


The practical aims at familiarzing you with CLM-FATES. 

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

> ### {% icon comment %} Background
>
> Background about ESMs, CTSM and FATES. We can have slides too (TO DO separately).
>
> ##### Earth System Modelling (ESM)
>
> ##### The Community Terrestrial Systems Model
>
> ##### The Community Land Model
>
> ##### Functionally Assembled Terrestrial Ecosystem Simulator (FATES)
>
{:  .comment}

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial. If you are not inspired, you can name it *fates*.
>    {% include snippets/create_new_history.md %}
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>     TODO: input data for running FATES (so it can be run anywhere if not in data library).
>    https://zenodo.org/record/
>    ```
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Check the datatype is **tar**
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 4. Rename Datasets
>
>    As "https://zenodo.org/record/?????/files/fates_emerald_inputdata.tar" is not a beautiful name and can give errors for some tools, it is a good practice to change the dataset name by something more meaningfull. For example by removing `https://zenodo.org/record/?????/files/` to obtain `fates_emerald_inputdata.tar`, respectively.
>
>    {% include snippets/rename_dataset.md %}
>
> 5. Add a tag to the dataset corresponding to `fates`
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

## Opening up Climate JupyterLab

> ### {% icon hands_on %} Hands-on: Launch JupyterLab for Ocean / Atmosphere / Land / Climate Python ecosystem in Galaxy
>
> Currently JupyterLab for Ocean / Atmosphere / Land / Climate Python ecosystem in Galaxy is available on [Live.useGalaxy.eu](https://live.usegalaxy.eu) only. JupyterLab for Ocean / Atmosphere / Land / Climate Python ecosystem and not the default JupyterLab in Galaxy contains all the python packages and additional software we need for running Earth System Model, including Functionally Assembled Terrestrial Ecosystem Simulator (FATES). The default JupyterLab in Galaxy would not be sufficient for executing all the tasks in th
is tutorial.
>
> 1. Open the JupyterLab tool {% icon tool %} by clicking [here](https://live.usegalaxy.eu/?tool_id=interactive_tool_climate_notebook){:target="_blank"}​ with the following par
ameters:
>   - *Include data into the environment*: `fates_inputdata.tar`
> 2. Click Execute
> 3. The tool will start running and will stay running permanently
> 4. Click on the "User" menu at the top and go to "Active Interactive Tools" and locate the JupyterLab instance you started.
> 5. Click on your JupyterLab instance (please not that it may take a few minutes before you can click on the link to your jupyterLab instance).
>
{: .hands_on}


You should now be looking at a page with the JupyterLab interface:

![Jupyterlab climate session](../../images/jupyterlab_climate_session.png)

The input dataset is located in the `data` folder. We will explain later how to use it for running FATES.


# CLM-FATES single point simulations

# Analysis

## Analyzing FATES-CLM model outputs

## Comparisons with observations

# Customize your runs

# Create a workflow

# Share your workflow

# Conclusion

{:.no_toc}

We have learnt to run single-point simulations with FATES-CLM through the Galaxy Climate JupyterLab.
