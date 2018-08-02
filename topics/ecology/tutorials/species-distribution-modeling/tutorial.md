---
layout: tutorial_hands_on
topic_name: sdm
tutorial_name: species distribution modeling
---

# Introduction
{:.no_toc}

Species Modeling Distribution can help understand the distribution of a species depending of environmental parameters such as temperature and precipitation. It can also help understand the impact of climate change on the repartition of some species. This is done by associating data occurrences of a species with environmental data.

The goal of this study is to model a theorical ecological niche using Species Modeling Distribution through the interactive environment Wallace.    

# Step 1: Loading a dataset

In this study the data set are all imported from the tool `Get species occurrences data and taxref informations` in Galaxy-E. With this tool, data are available from differents databanks like gbif, bison, inat and other.

>    > ### {% icon tip %} Tip: Importing data set from a data bank

>    > * Go into "upload files" (top panel) then "Get species occurrences data and taxref informations"
>    > * Fill the scientific name of the species wanted 
>    > * Choose the data source and the number of occurrences needed
>    > * Click on "Execute"

Then using `Tabular to CSV` on the new file to have the right format to use in Wallace.

# Step 2: Using Wallace

Wallace is an interactive interface which can simulate a species modeling distribution.

## Obtain occurrence data

With this you can either upload data directly from Wallace or you can upload file you've loaded earlier from Galaxy-E

> 1. Upload data from Galaxy-E
>    > * Check `Galaxy History User`
>    > * Select the corect csv file
> ### or
> 2. If you want to upload data from Wallace
>    > * Check `Query Database`
>    > * Select the databank of your interstest between Gbif, VertNet and BISON 
>    > * Type the name of the wanted species
>    > * And then set the number of occurrences

You now have your data for the next step.

## Process occurrence data

Here you'll have to chose the occurrences you want to use for the rest of your model. To do so, you have three options

> * Option 1: `Select Occurrences On Map`
> With this, you have to select your occurrence on the map by delimiting a geografic area you want to use.

> * Option 2: `Remove Occurrences By ID`
> You'll be remooving occurrences you don't need, or want, by their ID.

> * Option 3: `Spacial thin`
> This can allow you to select occurrences by setting a minimum distance (in km) beetween the different occurrences. For exemple:
>    > If you type 30 km, you'll end up with all the occurrences on the map which are at minimum 30km from each other.

For this part you have to know where you want your study to take place with the occurrence informations you have.
Best would be not to do it on data occurrences, but to create multiple sub-habitats to compare the parametres after.

## Obtain Environmental Data

The `WorldClim Bioclims` module will provide a raster with environmental variations from online sources. The Variables are mostly about temperature and precipitations. This will later be associated with the occurrences data.
The raster is composed of environmental info. Each layer of the raster contain a climatic variable; going from BIO1 = Anual mean temperature, to BIO19 = Precipitation of Coldest Quarter.

> To load these data you can either:
>    > * Use `WorldClim Bioclims` module
>    > * Load your own raster from your Galaxy-E history by checking the box `Galaxy History User`

Here we use the `WorldClim Bioclims` module.
After loading the environmental data, you can go to the next point.

## Process environnemental Data

Wallace will now associate environnmental data and occurrences data to make an area for your model.
> * First: `Choose Background Extent` make a buffer zone around the occurrences. You can chosse the size of the buffer zone around your occurrences by choosing the distance in `Study region buffer distance`. This allows you to control the area you'll be working with and on which a map of suitability will be made.

This is why you have to know what type of background extent you want to use.
>    > `Bounding box` will difine an area whit the occurrence centered

>    > `Minimum convex polygon` will make an area considering the repartition of your occurrences

>    > `Point buffers` will use occurrences to build a buffer zone aroud them

Then, to associate your occurrences to the environmental data, you'll have to choose the number of point to sample. This will cross

## Partition Occurrence Data

Partitioning data allows to divide a data set into subset, then make a model on each of all the groups but one and test it on the last one (assuming that all the groups are independent. You'll have two options: 

* `Non-spatial Partition`,is a type of partition used when sure that no biais due to space, time or sampling methode
>    > 1. `Jakknife (k=n)` consider that each occurrence in the dataset is equal to a bin. This is usually used when you have a dataset with no known biais.
>    > 2. `Random k-fold` partition de data randomly in a nuber of bin set by the user with the option `Number of Folds`

* `Spatial Partition`, when there could be biais due to time, space or sampling methode
>    > 1.`Block (k=4)`
>    > 2.
>    > 3.

## Build and Evaluate Niche Model

## Visualize Model Results


{: .hands_on}
