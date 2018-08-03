---
layout: tutorial_hands_on
topic_name: sdm
tutorial_name: species distribution modeling
---

# Introduction
{:.no_toc}

Species Modeling Distribution can help understand the distribution of a species depending of environmental parameters such as temperature and precipitation. It can also help understand the impact of climate change on the repartition of some species. This is done by associating data occurrences of a species with environmental data.

The goal of this study is to model a theorical ecological niche and predict futur repartition using Species Modeling Distribution through the use of Wallace interactive environment on Galaxy. We'll use the data occurrences of Chrysemys Picta (Painted turtle) from the tool `Get species occurrences data and taxref informations` in Galaxy-E applied on the North America region.   

# Step 1: Loading a dataset

In this study the datasets are all imported from the tool `Get species occurrences data and taxref informations` in Galaxy-E. With this tool, data are available from differents databanks like [GBIF](https://www.gbif.org/), [bison](https://www.gbif.org/), [iNaturalist](https://www.inaturalist.org/) and others.

>    > ### {% icon tip %} Tip: Importing data set from a data bank

>    > * Go into "upload files" (top left panel) then "Get species occurrences data and taxref informations"
>    > * Fill the "Scientific name" with `"Chrysemys Picta"`
>    > * Choose the data source `"Gbif"` and set the number of occurrences on `"10000"` 
>    > * Click on "Execute"
You now have a file with about 9508 occurrences

Because you only need informations about occurrences and their location: 
> use the tool `Couper des colonnes d'un jeu de données tabulé` 
> then in "Couper les colonnes" type `"c1,c2,c3"`

Then using the tool `Tabular to CSV` on this new file to have the right format to use in Wallace.

TODO question: For what stand c1, c2,c3
answer: c1 is the species name and c2 &c3 are respectively longitude an latitude correspoinding to each occurrences of the file

# Step 2: Using Wallace

[Wallace](https://wallaceecomod.github.io/) ([source code](https://github.com/wallaceEcoMod/wallace), [CRAN page](https://cran.r-project.org/web/packages/wallace/index.html)) is a R Shiny app integrated into Galaxy as an interactive environment which can simulate a species modeling distribution.

## Obtain occurrence data

With this you can either upload file you've loaded earlier from Galaxy-E data or you can upload  data directly from Wallace

> 1. Upload data from Galaxy-E
>    > * Check `Galaxy History User`
>    > * Select the correct csv file
> ### or
> 2. If you want to upload data from Wallace
>    > * Check `Query Database`
>    > * Select the databank of your interest between Gbif, VertNet or BISON 
>    > * Type the name of the wanted species
>    > * And then set the number of occurrences

You now have your data for the next step.

## Process occurrence data

Here you'll have to chose the occurrences you want to use for the rest of your model. To do so, you have three options

> * Option 1: `Select Occurrences On Map`
> With this, you have to select your occurrence on the map by delimiting a geographic area you want to use.


> * Option 2: `Remove Occurrences By ID`
> You'll be removing occurrences you don't need, or want, by their ID.

> * Option 3: `Spacial thin`
> This can allow you to select occurrences by setting a minimum distance (in km) between the different occurrences. For exemple:
>    > If you type 30 km, you'll end up with all the occurrences on the map which are at minimum 30km from each other.

To have better idea about the robustness of you species distrbution model, it would be better at this step, to don't select the entire set of data occurrences, but repeat the complete Wallace workflow on occurences subgroups, one by one.
## Obtain Environmental Data

The `WorldClim Bioclims` module will provide a raster with environmental variations from online sources. The [Bioclimatic variables](http://www.worldclim.org/bioclim) are describing temperature and precipitations variations. This will later be associated with the occurrences data.
The raster is composed of environmental information. Each layer of the raster contain a climatic variable; going from BIO1 = Anual mean temperature, to BIO19 = Precipitation of Coldest Quarter.

> To load these data you can either:
>    > * Use `WorldClim Bioclims` module
>    > * Load your own raster from your Galaxy-E history by checking the box `Galaxy History User`

Here we use the `WorldClim Bioclims` module.
After loading the environmental data, you can go to the next point.

## Process environnemental Data

Wallace will now associate environnmental data and occurrences data to make an area for your model.
> * First: `Choose Background Extent` make a buffer zone around the occurrences. You can chosse the size of the buffer zone by choosing the distance in `Study region buffer distance`. This allows you to control the area you'll be working with and on which a map of suitability will be made.

This is why you have to know what type of background extent you want to use.
>    > `Bounding box` will define an area whit occurrences centered

>    > `Minimum convex polygon` will make an area considering the repartition of your occurrences

>    > `Point buffers` will use occurrences localities to build a buffer zone aroud each of them

Then, to associate your occurrences to the environmental data, you'll have to choose the number of points to sample. This will cross

## Partition Occurrence Data

Partitioning data allows to divide a data set into subsets (ie bins), then make a model on each of all the groups but one and test it on the last one (assuming that all the groups are independent). You'll have two options: 

* `Non-spatial Partition`, is a type of partition used when there is no biais due to space, time or sampling method
>    > 1. `Jakknife (k=n)` consider that each occurrence in the dataset is equal to a bin. This is usually used when you have a dataset with no known biais.
>    > 2. `Random k-fold` partition the data randomly in a number of bins set by the user with the option `Number of Folds`

* `Spatial Partition`, when there could be biais due to time, space or sampling method
>    > 1. `Block (k=4)` divide the area in four an put equally into four bins, the different occurrences.
>    > 2. `Checkboard (k=2)` uses two bins according to the position of the occurence on the grid.
>    > 3. `Checkboard (k=4)` uses four bins according to the position of the occurence on the grid.
>
> For both of these technics the number of occurrences into each bin may vary.
> ⚠️ 

## Build and Evaluate Niche Model

Wallace can now build different models using either: 1) the presence-only approach BIOCLIM (Module BIOCLIM); or 2) the presence-background algorithm Maxent (Module Maxent). To evaluate these models, Wallace determines the ability to predict localities that are held out of the training ones (the ones who are used to construct models) and provide evaluation metrics as the AUC mean. One often used AUC threshold is 0,75 to consider a model accurate. The closer to 1 the better.

## Visualize Model Results

After the precedent step you can now model your theoretical niche.
first with `BIOCLIM Envelope Plots`you can make a chart and choose the parameters of interest and see how the data responds adapting the threshold for more accuracy 
