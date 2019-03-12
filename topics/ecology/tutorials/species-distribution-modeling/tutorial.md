---
layout: tutorial_hands_on
title: "Species distribution modeling"
zenodo_link: 
questions:
    - "How to build & evaluate a niche model, visualize and project it ?"
objectives:
    - "Find and download occurences data from GBIF"
    - "Find and download environmental data"
    - "Process both occurences and environmental"
    - "Partition occurence data"
    - "Model a theorical ecological niche and predict species distribution in a future climate scneario by using SDM"
requirements:
time_estimation: "1h"
key_points:
    - ""
contributors:
    - emichn
    - yvanlebras
    - bebatut
---

# Introduction
{:.no_toc}

Species Distribution Modeling (SDM) can help understand the distribution of a species depending on its environment. It can also attempt to quantify the impact of climate change on the species habitat, direct conservation efforts and predict invasive species distributions. This is done by associating data of species occurences (observations) with a set of environmental data (such as temperature and precipitation).

The goal of this tutorial is to model a theorical ecological niche and predict species distribution in a future climate scenario by using SDM with the Wallace interactive environment on Galaxy. We'll use the data occurrences of US *Chrysemys picta* ([Painted turtle](https://en.wikipedia.org/wiki/Painted_turtle)) from the North America region.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Step 1: Import data from public databank

In this study the datasets are all imported from the [GBIF](https://www.gbif.org/) databanks. It is also possible to import from [bison](https://www.gbif.org/), [iNaturalist](https://www.inaturalist.org/) and others.


> ### {% icon hands_on %} Hands-on: Import the data
> 1. Create a new history for this tutorial and give it a proper name
>
>    {% include snippets/create_new_history.md %}
>    {% include snippets/rename_history.md %}
> 
> 2. **Get species occurrences data** {% icon tool %} with the following parameters
>    - *"Scientific name"*: `Chrysemys picta`
>    - *"Data source"*: `gbif`
>    - *"Number of records to return"*: `10000`
>
> 3. Inspect the generated file with the Scratchbook
{: .hands_on}

We have now a tabular file with about 10,001 lines and many columns.

> ### {% icon question %} Questions
> 1. What do the lines in the files represent?
> 2. What are the columns 1, 2, 3 and 50 of the generated file?
>
> > ### {% icon solution %} Solution
> > 1. Each line represent an observation of Chrysemys picta, with its location, and some other metadata
> > 2. Column 1 is the name, column 2 the longitude of the observation of Chrysemys picta, column 3 the latitude and column 50 the country code of the location
> {: .solution }

We would to extraction the 4 columns previously cited and keep only occurence records from US.

> ### {% icon hands_on %} Hands-on: Import the data
> 1. **Cut columns** {% icon tool %} with the following parameters
>    - *"Cut columns"*: `c1,c2,c3,c50`
>    - *"Delimited by"*: `Tab`
>    - {% icon param-file %} *"From"*: output of **Get species occurrences data**
>
> 2. Inspect the generated file and check the 4 columns
> 3. **Filter data** {% icon tool %} with the following parameters
>    - {% icon param-file %} *"Filter"*: output of **Cut columns**
>    - *"With following condition"*: `c4=='US'`
>    - *"Number of header lines to skip"*: `1`
> 3. Inspect the generated files
>
>    > ### {% icon question %} Questions
>    >
>    > How many occurences have been conserved?
>    >
>    > ### {% icon solution %} Solution
>    > 71.18% of the occurences (7,118) are conserved
>    > {: .solution }
>    {: .question}
>
> 4. **Tabular to CSV** {% icon tool %}
>    - {% icon param-file %} *"tabular file"*: output of **Filter**
>    - *"output csv Separator"*: `,`
>    - *"Header in file"*: `Yes`
{: .hands_on}

# Step 2: Development and evaluation of SDM using Wallace

[Wallace](https://wallaceecomod.github.io/) is a R Shiny app integrated into Galaxy. It provides an interactive environment for the rapid and effective development and evaluation of SDM, including data download, cleaning, partitioning, modeling, visualisation and predictions.

## Obtain occurrence data

With this you can either upload file you've loaded earlier from Galaxy data or you can download data directly from Wallace. Let's use the data from your Galaxy history:

> ### {% icon hands_on %} Hands-on: Launch Wallace
> 1. Click on **Visualization** (top panel) and then **Interactive Environments**
> 2. Select
>    - *"GIE"*: `interactiveShiny`
>    - *"Image"*: `valentinchdock/wallace-galaxy-ie`
>    - *"Datasets"*: output of **Tabular to CSV**
> 3. In Wallace
>    1. Go to **1 Occ Data**
>    2. In **Obtain Occurrence Data**
>       - *"Modules Available"*: `Galaxy History User`
>       - *"Select from your Galaxy History User csv file"*: output of **Tabular to CSV**
>       - Click on **Load Occurrences**
{: .hands_on}

In the main panel of Wallace, the data will load. They can now be inspected using the map, the table, etc

You now have your occurence records on Wallace!

> ### {% icon question %} Questions
>
> *TODO: Add a question related to the data there*
> 
> > ### {% icon solution %} Solution
> > 
> >    {: .solution }
{: .question}


## Process occurrence data

For the modelling, we need to the choose the occurences we want to use. There is 4 ways to do it in Wallace in the **2 Process Occ** panel:

1. Selecting occurrences on the map by delimiting an interesting geographic area
2. Removing occurrences by ID
3. Delimiting a spacial thin

   This allow to select occurrences by setting a minimum distance (in km) between the different occurrences. For exemple, by typing 30km, we will get all the occurrences on the map which are at minimum 30km from each other.

4. Considering all the occurences

We want to work on the data from the US, we will select the occurrences with the 1st option

> ### {% icon hands_on %} Hands-on: Select interesting occurrences
> 1. Go to **2 Process Occ**
> 2. In **Process Occurrence Data**
>    - *"Modules Available"*: ` Select Occurrences On Map`
> 3. In the middle panel,
>    1. Click on the polygon icon on the map
>    2. Draw a polygon to delimite the area
> 4. On the left panel, click on **Select Occurrences**
{: .hands_on}


## Obtain Environmental Data

We would like now to provide environmental variables, that are considered relevant for biological purposes. The [Bioclimatic variables](http://www.worldclim.org/bioclim) consist of such original and derived variables. We will associate them with the occurrence data, using the WorldClim Bioclims module. 

This module provides a raster with environmental variables from online sources. The raster is composed of environmental information. Each layer of the raster contains a climatic variable: starting from BIO1 = Annual mean temperature, to BIO19 = Precipitation of Coldest Quarter.

> ### {% icon hands_on %} Hands-on: Obtain Environmental Data
> 1. Go to **3 Env Data**
> 2. In **Obtain Environmental Data**
>    - *"Modules Available"*: `WorldClim Bioclims`
>
>      > ### {% icon comment %} Using your own raster
>      > If you prefer to load your own raster:
>      > 1. Import it into Galaxy in the GeoTIFF format
>      > 2. In **Obtain Environmental Data**
>      >    - *"Modules Available"*: `Galaxy History User`
>      >    - *"Select from your Galaxy History User csv file"*: imported GeoTIFF file
>      {: .comment}
>
>    - *"Select WorldClim bioclimatic variable resolution"*: `10 arcmin`
>
>      This is the lowest resolution.
>
> 3. Click on the **Load Env Data**
{: .hands_on}

## Process environmental Data

Wallace will now associate environmental data and occurrences data to train a model.

1. It first creates a buffer zone around the occurrences

    We can choose the size of the buffer zone to control the area we will be working with and on which a map of suitability will be made.

    Several background extent can be used:
    - Bounding box will define an area where the occurrences centered
    - Minimum convex polygon will make an area considering the repartition of the occurrences
    - Point buffers will use occurrences localities to build a buffer zone around each occurences

2. It samples Background Points

> ### {% icon hands_on %} Hands-on: Obtain Environmental Data
> 1. Go to **4 Process Envs**
> 2. In **Process Environmental Data**
>    - *"Modules Available"*: `Select Study Region`
> 3. In **Step 1: Choose Background Extent**
>    - *"Background Extents"*: `Minimum convex polygon`
>    - *"Study region buffer distance (degree)"*: `1`
> 4. Click on **Select** 
> 5. In **Step 2: Sample Background Points**
>    - *"No. of background points"*: `100000`
> 6. Click on **Sample**
{: .hands_on}

## Partition Occurrence Data

By partitioning data, a dataset is divided into subsets (i.e. bins). A model is then built on each of subsets but one and test it on the last one (assuming that all the groups are independent). 

There is 2 possible way to partition data:

1. Non-spatial Partition: a partition used when there is no bias due to space, time or sampling method
    1. Jakknife (k=n): each occurrence in the dataset is considered as equal to a bin. This is usually used when you have a small dataset with no known bias
    2. Random k-fold: partition the data randomly in a number of bins set by the user with the option `Number of Folds`
2. Spatial Partition: a partition used when there could be biais due to time, space or sampling method
    1. Block (k=4): the area is divided in four and the different occurrences are put equally into four bins
    2. Checkerboard (k=2): two bins are used according to the position of the occurence on the grid
    3. Checkerboard (k=4): four bins are used according to the position of the occurence on the grid. This require an aggregation factor, which is the size of a second grid put on a first one. For exemple, with a factor 4, the grids size will be 4x4

        ![Checkerboard 2](https://raw.githubusercontent.com/emichn/training-material/patch-1/topics/ecology/tutorials/species-distribution-modeling/Images/Checkerboard.png "Example of checkerboard (k=4) with a factor 4. Source: doi:10.0.4.87/2041-210X.12261")

For both of these technics the number of occurrences into each bin may vary.


`Use spatial partition` on these data choose `Checkerboard 2 (k=4)` and an aggregation factor of `6`


## Build and Evaluate Niche Model

Wallace can build different models using either: 1) the presence-only approach BIOCLIM (Module BIOCLIM); or 2) the presence-background (presense-pseudo absence) algorithm Maxent (Module Maxent). To evaluate these models, Wallace computes the performance on a hold-out dataset (data not used for training) and provide evaluation metrics as the AUC (Area Under the Curve) mean. As a rule of thumb, an AUC of 0.75 and above is considered good, and closer to 1 is better.

## Visualize Model Results

After the precedent step you can now model your theoretical niche.

First with `BIOCLIM Envelope Plots`you can make a chart and choose the parameters of interest and see how the data responds adapting the threshold for more accuracy.

With our data we can see that when we make a chart to simutate an ecological niche using `bio1` as y axis and `bio12` as x axis, with a threshol of `0,75`, the optimum environment parameter for this species is between `5°C` and `15°C`(on the graph, values are *10), for an annual precipitation between aproximatively `700mm` and `1250mm`.

Then with `Map Prediction`, you can either select `no threshold` to have a gradient of predicted presence or use, like in this study, `minimum training presence` to have a map with the predicted presence and predicted absence.


## Project Model

Wallace can also help you use the model you just trained to predict possible species distributions in a different area, outside of the sampled one.

We'll apply the concept on a part of the Canada. Select `Project to New Extent` then draw a polygone around a part of Canada and use `Minimum Training Presence`.

You can also project not only in space, but time, and a different climate. Select `Project to New Time` if you want to do this.
You need to choose a `GCM` (global circulation model). They are model made to predict atmospheric fluctuation and then study climate change. Each one is different and use parameters like ocean atmosphere and others.

Here using different model you can see the evolution of the predicted presence of `Chrysemys Picta` in Canada in 2050.

For exemple we can use the model [CCSM4](http://www.cesm.ucar.edu/models/ccsm4.0/), a US model based on earth circulation. We can try it with differents RCP scenarios, wich are scernarios about the amout of greenhouse gases emitted in the near futur. This allows us to have differents predicted presence models. Once again, we will use the `Minimum Training Presence` threshold.

![CCM4 with a 2.6 RCP](https://raw.githubusercontent.com/emichn/training-material/patch-1/topics/ecology/tutorials/species-distribution-modeling/Images/GCM_CCSM4_RCP_2.6.png "CCM4 with a 2.6 RCP")

![CCM4 with a 8.5 RCP](https://raw.githubusercontent.com/emichn/training-material/patch-1/topics/ecology/tutorials/species-distribution-modeling/Images/GCM_CCSM4_RPC_8.5.png "CCM4 with a 8.5 RCP")

You can see here the list of the different [Global Circulation Model](https://en.wikipedia.org/wiki/General_circulation_model#Atmospheric_and_oceanic_models)

# Conclusion

We’ve been able through Galaxy-E, to load a dataset of occurrences used in the shiny app Wallace and model the repartition of *Chrysemys picta* (Painted turtle) with the Species Distribution Modeling (SDM) method. It allowed us to visualize it’s ecological niche and how climate change can influence it’s futur repartition on North America. The project saved can help for futur similar studies.

# References

About [Wallace](https://wallaceecomod.github.io/) 

Guissan, A. et *al.*, Predicting species distributions for conservation decisions. Ecology Letters, 16, 1424–1435  (2013).

Booth, T. H. et *al.*, BIOCLIM: the first species distribution modelling package, its early applicationsand relevance to most current MAXENT studies. Diversity and Distributions, 20, 1–9 (2014).

Muscarella, R. et *al.*, ENMeval: An R package for conducting spatially independent evaluations and estimating optimal model complexity for MAXENT ecological niche models. British Ecological Society, Methods in Ecology and Evolution, 5, 1198–1205 (2014).

[Here](http://www.ipcc-data.org/guidelines/pages/gcm_guide.html) for informations on Global Circulation Model(GCM), how it's done what is taken in consideration and more.

