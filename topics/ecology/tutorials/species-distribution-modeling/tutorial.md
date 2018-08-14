---
layout: tutorial_hands_on
topic_name: sdm
tutorial_name: species distribution modeling
---

# Introduction
{:.no_toc}

Species Distribution Modeling (SDM) can help understand the distribution of a species depending on its environment. It can also attempt to quantify the impact of climate change on the species habitat, direct conservation efforts and predict invasive species distributions. This is done by associating data of species occurences (observations) with a set of environmental data (such as temperature and precipitation).

The goal of this tutorial is to model a theorical ecological niche and predict species distribution in a future climate scneario by using SDM with the Wallace interactive environment on Galaxy. We'll use the data occurrences of US *Chrysemys picta* ([Painted turtle](https://en.wikipedia.org/wiki/Painted_turtle)) from the North America region.

# Step 1: Loading a dataset

In this study the datasets are all imported from the tool `Get species occurrences data and taxref informations` in Galaxy-E. With this tool, data are available from differents databanks like [GBIF](https://www.gbif.org/), [bison](https://www.gbif.org/), [iNaturalist](https://www.inaturalist.org/) and others.

>    > ### {% icon tip %} Tip: Importing data set from a data bank

>    > * Go into the "upload files" tool section (top left panel) then select the "Get species occurrences data and taxref informations" tool
>    > * Fill the "Scientific name" with `Chrysemys picta`
>    > * Choose the data source `gbif` and set the number of occurrences on `10000`
>    > * Set "Get species information from taxref database" to "No"
>    > * Click on "Execute"

You should now have a file with about `9562` occurrences.

Because you only need informations about occurrences and their location:
> use the tool `Couper des colonnes d'un jeu de données tabulé`
> then in "Couper les colonnes" type `c1,c2,c3,c44`

>    > ### {% icon question %} Questions
>    >
>    > 1. For what stand c1, c2, c3, c44 ?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>c1 is the species name and c2 &c3 are respectively longitude an latitude corresponding to each occurrences of the file. The fourth column contain country code to have a possibility to easily filter occurences by countries.</li>
>    >    </ol>
>    >    </details>
>    {: .question}

Then, as we want to keep only occurence records from US, we will use a tool to filter the data based on the fourth column. Can you find a way to do that alone?

>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>Using the "Filtrer des données sur une colonne en utilisant des expressions simples", you can enter the following condition ```c4=='US'```. You should specify the following field as 1 in order to skip the first line (it is a header line).</li>
>    >    </ol>
>    >    </details>
>    {: .question}

Finally, use the tool `Tabular to CSV` (use the default options) on this new file to convert it to a format suitable for further processing with Wallace.

# Step 2: Using Wallace

[Wallace](https://wallaceecomod.github.io/) ([source code](https://github.com/wallaceEcoMod/wallace), [CRAN page](https://cran.r-project.org/web/packages/wallace/index.html)) is a R Shiny app integrated into Galaxy, providing an interactive environment for the development and evaluation of SDM. It covers all the main parts of such a workflow, including data download, cleaning, partitioning, modeling, visualisation and predictions. It allows for a rapid and effective development of SDM.

## Obtain occurrence data

With this you can either upload file you've loaded earlier from Galaxy-E data or you can download data directly from Wallace. Let's use the data from your Galaxy history:

> 1. Upload data from Galaxy-E
>    > * Check `Galaxy History User`
>    > * Select the correct csv file with your *Chrysemys Picta* occurrences informations

You now have your occurence records on Wallace!

>    > ### {% icon question %} Questions
>    >
>    > 1. One point appears near the African continent. Can you propose a reason ? Is this an error? Which ID is it?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>It's likely an error where the coordinates were not filled. It proprely have a country code and an ID: 6783.</li>
>    >    </ol>
>    >    </details>
>    {: .question}


## Process occurrence data

Here you'll have to chose the occurrences you want to use for the rest of your model. To do so, you have 4 options

> * Option 1: `Select Occurrences On Map`
> With this, you have to select your occurrence on the map by delimiting a geographic area you want to use. Note that first you should draw the polygon (the icon is on the map), and then click on "Select Occurences". If you make a mistake you can click Reset.

> * Option 2: `Remove Occurrences By ID`
> You'll be removing occurrences you don't need, or want, by their ID.

> * Option 3: `Spacial thin`
> This can allow you to select occurrences by setting a minimum distance (in km) between the different occurrences. For exemple:
>    > If you type 30 km, you'll end up with all the occurrences on the map which are at minimum 30km from each other.

> * Option 4: skip this step and consider all records.

Because we want to work on the data from the US we'll select all the occurrences there with the first option: `Select Occurrences On Map`.

## Obtain Environmental Data

The `WorldClim Bioclims` module will provide a raster with environmental variables from online sources. The [Bioclimatic variables](http://www.worldclim.org/bioclim) consist of original and derived variables that are considered relevant for biological purposes. Those will later be associated with the occurrence data.
The raster is composed of environmental information. Each layer of the raster contains a climatic variable; starting from BIO1 = Anual mean temperature, to BIO19 = Precipitation of Coldest Quarter.

> To load these data you can either:
>    > * Use `WorldClim Bioclims` module
>    > * Load your own raster from your Galaxy-E history by checking the box `Galaxy History User`, the file must be in GeoTIFF format.

Here we use the `WorldClim Bioclims` module with the lowest resolution `10 arcmin`.
>    > nb: the occurrences situated where there is no environmental data are removed

After loading the environmental data, you can go to the next point.

## Process environmental Data

Wallace will now associate environmental data and occurrences data to train a model.
> * First: `Choose Background Extent` creates a buffer zone around the occurrences. You can choose the size of the buffer zone by choosing the distance on `Study region buffer distance` parameter. This allows you to control the area you'll be working with and on which a map of suitability will be made.

This is why you have to know what type of background extent you want to use.
>    > `Bounding box` will define an area whit occurrences centered
>    >
>    > `Minimum convex polygon` will make an area considering the repartition of your occurrences
>    >
>    > `Point buffers` will use occurrences localities to build a buffer zone aroud each of them

Then, to associate your occurrences to the environmental data, you'll have to choose the number of points to sample.

For our study we made a `Minimum convex polygon` background extent with a `buffer area` of `1 degree` and use `100000` backgound points.

## Partition Occurrence Data

Partitioning data allows to divide a data set into subsets (ie bins), then make a model on each of all the groups but one and test it on the last one (assuming that all the groups are independent). You'll have two options:

* `Non-spatial Partition`, is a type of partition used when there is no bias due to space, time or sampling method
>    > 1. `Jakknife (k=n)` consider that each occurrence in the dataset is equal to a bin. This is usually used when you have a small dataset with no known biais.
>    > 2. `Random k-fold` partition the data randomly in a number of bins set by the user with the option `Number of Folds`

* `Spatial Partition`, when there could be biais due to time, space or sampling method
>    > 1. `Block (k=4)` divide the area in four and put equally into four bins, the different occurrences.
>    > 2. `Checkerboard (k=2)` uses two bins according to the position of the occurence on the grid.
>    > 3. `Checkerboard (k=4)` uses four bins according to the position of the occurence on the grid. This require an aggregation factor, which is the size of a second grid put on a first one.
>    >    >
>    >    > * For exemple: if you use a factor 4, the grids size will be 4x4

![Checkerboard 2](https://raw.githubusercontent.com/emichn/training-material/patch-1/topics/ecology/tutorials/species-distribution-modeling/Images/Checkerboard.png "doi:10.0.4.87/2041-210X.12261")

> For both of these technics the number of occurrences into each bin may vary.

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

