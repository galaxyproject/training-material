---
layout: tutorial_hands_on
topic_name: sdm
tutorial_name: species distribution modeling
---

# Introduction
{:.no_toc}

Species Distribution Modeling can help understand the distribution of a species depending of environmental parameters such as temperature and precipitation. It can also help understand the impact of climate change on the repartition of some species. This is done by associating data occurrences of a species with environmental data.

The goal of this study is to model a theorical ecological niche and predict future repartition using Species Distribution Modeling through the use of Wallace interactive environment on Galaxy. We'll use the data occurrences of US *Chrysemys Picta* ([Painted turtle](https://fr.wikipedia.org/wiki/Tortue_peinte)) from the North America region.   

# Step 1: Loading a dataset

In this study the datasets are all imported from the tool `Get species occurrences data and taxref informations` in Galaxy-E. With this tool, data are available from differents databanks like [GBIF](https://www.gbif.org/), [bison](https://www.gbif.org/), [iNaturalist](https://www.inaturalist.org/) and others.

>    > ### {% icon tip %} Tip: Importing data set from a data bank

>    > * Go into the "upload files" tool section (top left panel) then select the "Get species occurrences data and taxref informations" tool
>    > * Fill the "Scientific name" with `Chrysemys Picta`
>    > * Choose the data source `Gbif` and set the number of occurrences on `10000` 
>    > * Click on "Execute"
You now have a file with about `9508` occurrences

Because you only need informations about occurrences and their location: 
> use the tool `Couper des colonnes d'un jeu de données tabulé` 
> then in "Couper les colonnes" type `c1,c2,c3,c44`

TODO question: For what stand c1, c2, c3, c44
answer: c1 is the species name and c2 &c3 are respectively longitude an latitude corresponding to each occurrences of the file. The fourth column contain country code to have a possibility to easily filter occurences by countries.


Then using the tool `Tabular to CSV` on this new file to have the right format to use in Wallace.

Finally, we want to keep only occurence records from US. To do that, we will use a tool to filter dataset on the fourth column. Can you find a way to do that alone ?

TODO tips/answer: Using the "Filtrer des données sur une colonne en utilisant des expressions simples", you can enter the following condition ```c4=='US'``` specifying that we need to don't consider the first line as it's a header.

# Step 2: Using Wallace

[Wallace](https://wallaceecomod.github.io/) ([source code](https://github.com/wallaceEcoMod/wallace), [CRAN page](https://cran.r-project.org/web/packages/wallace/index.html)) is a R Shiny app integrated into Galaxy as an interactive environment which can model a species distribution.

## Obtain occurrence data

With this you can either upload file you've loaded earlier from Galaxy-E data or you can upload  data directly from Wallace. Here you will select data previously imported and filtered on your Galaxy history.

> 1. Upload data from Galaxy-E
>    > * Check `Galaxy History User`
>    > * Select the correct csv file with your *Chrysemys Picta* occurrences informations 

You now have your occurence records on Wallace!

## Process occurrence data

Here you'll have to chose the occurrences you want to use for the rest of your model. To do so, you have three options

> * Option 1: `Select Occurrences On Map`
> With this, you have to select your occurrence on the map by delimiting a geographic area you want to use.


> * Option 2: `Remove Occurrences By ID`
> You'll be removing occurrences you don't need, or want, by their ID.

> * Option 3: `Spacial thin`
> This can allow you to select occurrences by setting a minimum distance (in km) between the different occurrences. For exemple:
>    > If you type 30 km, you'll end up with all the occurrences on the map which are at minimum 30km from each other.

Because we want to work on the data from the US we'll select all the occurrences there with the first option:`Select Occurrences On Map`  

## Obtain Environmental Data

The `WorldClim Bioclims` module will provide a raster with environmental variations from online sources. The [Bioclimatic variables](http://www.worldclim.org/bioclim) are describing temperature and precipitations variations. This will later be associated with the occurrences data.
The raster is composed of environmental information. Each layer of the raster contains a climatic variable; going from BIO1 = Anual mean temperature, to BIO19 = Precipitation of Coldest Quarter.

> To load these data you can either:
>    > * Use `WorldClim Bioclims` module
>    > * Load your own raster from your Galaxy-E history by checking the box `Galaxy History User`

Here we use the `WorldClim Bioclims` module with the lowest resolution `10 arcmin`.
>    > nb: the occurrences situated where there is no environmental data are removed

After loading the environmental data, you can go to the next point.

## Process environnemental Data

Wallace will now associate environnmental data and occurrences data to make an area for your model.
> * First: `Choose Background Extent` creates a buffer zone around the occurrences. You can choose the size of the buffer zone by choosing the distance on `Study region buffer distance` parameter. This allows you to control the area you'll be working with and on which a map of suitability will be made.

This is why you have to know what type of background extent you want to use.
>    > `Bounding box` will define an area whit occurrences centered

>    > `Minimum convex polygon` will make an area considering the repartition of your occurrences

>    > `Point buffers` will use occurrences localities to build a buffer zone aroud each of them

Then, to associate your occurrences to the environmental data, you'll have to choose the number of points to sample. This will cross

For our study we made a `Minimum convex polygon` background extent with a `buffer area` of `1 degree` and use `100000` backgound points. 

## Partition Occurrence Data

Partitioning data allows to divide a data set into subsets (ie bins), then make a model on each of all the groups but one and test it on the last one (assuming that all the groups are independent). You'll have two options: 

* `Non-spatial Partition`, is a type of partition used when there is no biais due to space, time or sampling method
>    > 1. `Jakknife (k=n)` consider that each occurrence in the dataset is equal to a bin. This is usually used when you have a small dataset with no known biais.
>    > 2. `Random k-fold` partition the data randomly in a number of bins set by the user with the option `Number of Folds`

* `Spatial Partition`, when there could be biais due to time, space or sampling method
>    > 1. `Block (k=4)` divide the area in four and put equally into four bins, the different occurrences.
>    > 2. `Checkerboard (k=2)` uses two bins according to the position of the occurence on the grid.
>    > 3. `Checkerboard (k=4)` uses four bins according to the position of the occurence on the grid. This require an aggregation factor, which is the size of a second grid put on a first one.

>    > * For exemple: if you use a factor 4, the grids size will be 4x4 "insert image from muscarella & al"

> For both of these technics the number of occurrences into each bin may vary.
 
flag : `Use spatial partition` on these data choose `Checkerboard 2 (k=4)` and an aggregation factor of `6`

## Build and Evaluate Niche Model

Wallace can now build different models using either: 1) the presence-only approach BIOCLIM (Module BIOCLIM); or 2) the presence-background algorithm Maxent (Module Maxent). To evaluate these models, Wallace determines the ability to predict localities that are held out of the training ones (the ones who are used to construct models) and provide evaluation metrics as the AUC mean. One often used AUC threshold is 0,75 to consider a model accurate. The closer to 1 the better.

## Visualize Model Results

After the precedent step you can now model your theoretical niche.

First with `BIOCLIM Envelope Plots`you can make a chart and choose the parameters of interest and see how the data responds adapting the threshold for more accuracy.

With our data we can see that when we make a chart to simutate an ecological niche using `bio1` as y axis and `bio12` as x axis, with a threshol of `0,75`, the optimum environement parameter for this species is between `5°C` and `15°C` for an annual precipitation between aproximatively `700mm` and `1250mm`.

You can either select `no threshold` to have a gradient of predicted presence or use, like in this study, `minimum training presence` to have a map with the predicted presence and predicted absence.
flag: insert photo of map prediction.

## Project Model

Wallace can also apply the model you've juste created on a specified zone, here the US, to another zone and create a map of suitability.

We'll aply the concept on a part of the Canada.
Flag: insert photo

To go further, you can test your model on a specified zone, years from now.
You need to choose a `GCM` (global circulation model). They are model made to predict atmospheric fluctuation and then study climate change. Each one is different and use parameters like ocean atmosphere and other.

Here using different model you can see the evolution of the predicted presence of `Chrysemys Picta` in Canada in 2050.

For exemple we can use the model [HadGEM2-ES](https://portal.enes.org/models/earthsystem-models/metoffice-hadley-centre/hadgem2-es) (Earth System), a model based on earth circulation, and try it with differents RCP scenarios. RCP are scernarios about the amout of greenhouse gases emitted in the near futur, this allows us to have differents predicted presence models.
Flag: insert two photos of 4,5 & 8,5

You can see here a brief list of some of the [Global Circulation Model](http://www.ipcc-data.org/sim/gcm_monthly/SRES_AR4/index.html)

# Conclusion

