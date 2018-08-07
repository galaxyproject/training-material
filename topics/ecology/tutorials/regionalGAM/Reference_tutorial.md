---
layout: tutorial_hands_on
topic_name: ecology
tutorial_name: Regional GAM
---

# Introduction
{:.no_toc}


This tutorial will show how to study species phenology through the computation of abundance index and trends. It will explain you how to use different [regionalGAM](https://github.com/RetoSchmucki/regionalGAM) tools on Galaxy-E allowing you to deal with datasets containing occurences informations for various species per site and per date through a couple of years.
After a certain numbers of steps, you will be able to extract single species data and study related phenology through the years. The goal of this exercise is to be able to create abundance trend over time and biodiversity indicators. Following these indicators allow to follow trends in terms of population dynamics. You could for example try to predict the occurences of one specific species in a certain type of environnement using the prediction model of climate evolution. Based on charts that you will generate, you could try to explain the evolution of a species with environmental data (temperatures variations, modifications of the environmental conditions).
You will basically learn how to create a file on the basis of which you can create a visual material that can be quite easily understood and therefore be efficient for a large audience.


> ### Agenda
> In this tutorial, we will cover:
1. Pre-processing
> {:pre-processing}
2. Displaying the occurence of the chosen species through the years
> {:Displaying the occurence of the chosen species through the years} 
3. Conclusion 
> {:conclusion}

# Step 1: Pre-processing

The goal of the first step is to upload and prepare the file so that it will be usable for the *regional GAM* analysis (See [this warning](#inputdatawarning) for more information about the input file.
First of all, you will have to upload the files on Galaxy-E and then you might have to use some data handling tools to be able to use *regional GAM* tools.

>  ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial and give it a proper name.
> 2. Import the following files from [Zenodo](https://zenodo.org/record/1324204#.W2BmRn7fNE4) or from a data
>    library named `regional GAM data tutorial`
>
>    ```
>    CSV dataset with only one species:
>    https://zenodo.org/record/1324204/files/regional%20GAM%20data.csv?download=1
>    ```
>   
> ### {% icon tip %} Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch data**
>    > * Paste the link into the text field
>    > * Press **Start** and **Close** the window
>    {: .tip}
>
> ### {% icon comment %} Comment
>
> ⚠️ <a name="inputdatawarning"></a>Please note that the file must contain a header corresponding to: ```"SITES","SPECIES","YEAR","MONTH","DAY","COUNT"```, and that all the non numeric content must be between double quotes as "x" and that separators have to be ",". 

>    > ## <a name="resampling"></a>Re-sampling.  
When the dataset contains many details, it lengthens the file processing time therefore it can be very useful to learn how to hide the informations you don't need. For example, the list of SITE (look at the column with header `SITE)`  of the dataset you are using is really long and the SITES are classified into sub-sites. Here, we will assume that your file doesn't really need be as precise and this is the reason why you have to specify you don't want the sub-sites. To create a new "down-sampled" file, you can follow these steps:   

> ### {% icon hands_on %} Hands-on: hiding some informations

>    > 1. Use the `CSV to tabular` tool to first create a tabular file from your csv one (with only one species). This is a mandatory step as further tools are only working on tabular files!
>    > 2. Search for the tool `trouver et remplacer des patterns dans des colonnes` on the file on CSV with the following  parameters.
>    >    >  * Select the input file & the column with the SITE header.
>    >    >  * Click on`"insert checks"`
>    
>    >    >  * "Trouver l'expression suivante": `"(\.[0-9]+)"` which specifies that you don't want the sub-sites (all suites of digits following a "." character) to be taken into account.
>    >    >  * "Remplacement": leave it empty.

>   > ### {% icon question %} Questions
>   >
>    > 1. After having successfully deleted the sub-sites informations, can you look at the original dataset and this new one and say how many sites you had, and you have now? You will maybe need to use tools like `Compter le nombre d'occurrences de chaque enregistrement`?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>The dataset contains 6 sites now against 1144 before "down-sampling". </li>
>    >    </ol>
>    >    </details>
>    {: .question}

>    > 3. Now you can regenerate a CSV file using the `tabular to CSV` tool on the ouptut from **trouver et remplacer des patterns dans des colonnes**.

{: .hands_on}


 # Step 2: Displaying the occurence of the chosen species through the years 
 
Now you have a file containing all the data on the species of interest. The main goal of this step is basically to create a material that can be used to generate charts. What you could also do, for example, would be to compare the evolution of various species through the years in the same site. You would have to superpose the different graphs on one another.
>

> ### {% icon hands_on %} Hands-on: Phenology

This step will allow you the show the phenology of a species and then to create charts representing it. In the second part, you will learn that it is possible to show the phenology of various species on a single chart allowing to compare them and analyse them more easily. 


> 1. Search for the tool `flight curve` with the following parameters: 
> * "Fichier de comptage": `output` from **flight curve**.
>
> 🔹 Based on the flight curve, you can create a line chart which shows the occurence of the species through the years on a very visual material 

![Phenology chart](https://github.com/Claraurf/training-material/blob/ecology/topics/ecology/tutorials/regionalGAM/Images/Phenology%20regional%20GAM.png "This shows the occurrence of Ardea alba and Bombycilla garrulus")

If you want to access the chart on an interactive interface, you can click on the following link: [Chart on Galaxy](http://openstack-192-168-100-19.genouest.org/plugins/visualizations/charts/saved?id=417e33144b294c21)
>
>    > ### {% icon tip %} Visualiser 
>    > 1. Select `Charts`
>    > 2. Give it a proper name
>    > 3. Select a visualization: "line chart (NVD 3) 
>    > 4. Select data:
>    > * "Provide a label": `The name of the species`
>    > * "Pick a series color": Choose a color for the line 
>    > * "Data point labels": `Column 1`
>    > * "Values for x-axis": `Column 2`
>    > * "Values for y-axis": `Column 6`
>    > * "X-Axis label": `nm values`  
>    > * Y-Axis label: `Year`
>    > 8. Click on {% icon tip %} `Visualize`
>    > 9. Click on {% icon tip %} `save this visualization`if you are willing
>    > 5. Visualize
>    > 6. Click on {% icon tip %} `save this visualization`if you are willing to keep it
>
> ⚠️ Please note that it is possible to show the occurences of more than one species on a single chart. If you are interested in doing so, you should follow the tip below.

 > ### {% icon tip %} Visualiser 
 
>    > 1. Select `Charts`
>    > 2. Give it a proper name
>    > 3. Select a visualization: "line chart (NVD 3) 
>    > 4. Select data 
>    > * "Provide a label": `The name of the species 1` 
>    > * "Pick a series color": Choose a color
>    > * "Data point labels": `Column corresponding to the name of the species 1` 
>    > * "Values for x-axis": `Column corresponding to the year of the species 1`
>    > * "Values for y-axis": `Column corresponding to nm of the species 1`
>    > 5. Insert data series:
>    > * "Provide a label": `The name of the species 2` 
>    > * "Pick a series color": Choose a different color
>    > * "Data point labels": `Column corresponding to the name of the species 2` 
>    > * "Values for x-axis": `Column corresponding to the year of the species 2`
>    > * "Values for y-axis": `Column nm of the species 2`
>    > 6. You may repeat "Insert data series" as many times as needed depending on the number of different species you want represent
>    > 7. Click on {% icon tip %} `Customize`
>    > * "X-Axis label": `nm values`  
>    > * Y-Axis label: `Year`
>    > * "Use multi-panels": click on `No`(or you will have three different charts)
>    > 8. Click on {% icon tip %} `Visualize`
>    > 9. Click on {% icon tip %} `save this visualization`if you are willing to keep it

{: .hands_on}
>
>
> ### Abundance per year and per site

This will allow you to create a file showing the abundance per year of a chosen species in a certain site. Based on this file you will then learn how to represent this abundance on a chart. 
>
> 1. Look for the tool `Abundance index` with the following parameters:
> * "Fichier de comptage": `output` from **tabular to CSV**.  
> * "Flight curve output": `output` from **flight curve**.


> 🔹 Based on the abundance index, we can create a chart showing the annual abundance trend of a certain species per site. 
>    > 1. Select `Charts`
>    > 2. Give it a proper name
>    > 3. Select a visualization: "Bar diagram (NVD 3)" 
>    > 4. Select data 
>    > * "Data point labels": `Column 1` 
>    > * "Values for x-axis": `Column 3`
>    > * "Values for y-axis": `Column 4`
>    > 5. Visualize
>    > 6. Click on {% icon tip %} `save this visualization`if you are willing to keep it

{: .hands_on}

> ### Expected temporal trend

The expected temporal trend allows you to have an overview of the evolution of a species in a certain type of environment in the futur.

![Expected temporal trend](https://github.com/Claraurf/training-material/blob/ecology/topics/ecology/tutorials/regionalGAM/Images/Expected%20temporal%20trend.png "This shows the expected evolution of Ardea alba")

> ### {% icon hands_on %} Hands-on: Expected temporal trend
>    > 1. Look for the tool `Expected temporal trend` with the the following parameters: 
>    > * "Fichier tabulé, produit par l'outil ab_index": `output` from **abundance index**.
>    
> ⚠️ Please note that sometimes the expected temporal trend can't be done on dataset. If you want this action to work, the occurences on your dataset must lie between the month of April and the end of the month of September.

Note also that you will obtain two files resulting of the action above. The first will be the graph and the second will contains the values of "x".

> ### Linear regression 

The point of doing a linear regression is to determinate if the year has an influence on the abundance of a species. 

>    > 1. Look for the tool `linear regression` with the following parameters.
>    > * "Fichier produit par l'outil glmmpql/Expected temporal trend": `output 2` from **temporal trend**. 
>    > * "Fichier produit par l'outil ab_index": `output` from **abundance index**.

{: .hands_on}
 
# Conclusions  

{:.no_toc}

In this tutorial, you have analyzed regional GAM data to extract useful informations in order to be able to show different tendencies of a chosen species. Therefore, you are now able to treat the dataset so that it shows only the data concerning one specific species of your choice. From there, you can show the occurrence of this species through the years first on a dataset and then on a visual chart. You have also learned how to represent on a single chart the occurences of various species. Afterwards, we have shown you how to create a dataset containing the informations on the abundance of a species per year and per site. Based on which you can henceforth visually represent the annual abundance trend on a chart. Thereafter, you have the possibility of showing the expected temporal trend, based on which you will be able to try predicting the future evolution a given species. The last part of this tutorial has shown you how to calculate the linear regression allowing you to determinate wether the year has an influence on the abundance of a species or not. 
