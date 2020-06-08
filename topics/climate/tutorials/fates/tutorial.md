---
layout: tutorial_hands_on
title: Functionally Assembled Terrestrial Ecosystem Simulator (FATES) Modelling with Galaxy
zenodo_link: ''
questions:
- What is the Community Terrestrial System Model (CTSM)?
- What is FATES?
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

> ### {% icon comment %} Comment
>
>
{: .comment}

The practical aims at familiarzing you with FATES. 

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
{:  .comment}

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial. If you are not inspired, you can name it *fates*.
>    {% include snippets/create_new_history.md %}
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>     TODO: input data for running FATES (so it can be run anywhere).
>    https://zenodo.org/record/
>    ```
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Check the datatype
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 4. Add a tag to the dataset corresponding to `fates`
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

# What is CTSM?

## Ecological, biogeochemical, biogeophysical, and hydrologic theory underpinning the CTSM

# Starting Interactive JupyterLab

# CLM-FATES single point simulations

# Analysis

## Analyzing FATES-CLM model outputs

## Comparisons with observations

# Create a workflow

# Share your workflow

# Conclusion

{:.no_toc}

We have learnt to run single-point simulations with FATES-CLM through the Galaxy Climate JupyterLab.
