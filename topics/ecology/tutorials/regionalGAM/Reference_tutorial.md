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

> ### {% icon comment %} Comment
>
> ⚠️ Please note that there are two other tutorials on regionalGAM, one is specific to [the use of RData input dataset](https://github.com/Claraurf/training-material/blob/ecology/topics/ecology/tutorials/regionalGAM/Rdata_tutorial.md) and the other explains you how to deal with [a multispecies dataset](https://github.com/Claraurf/training-material/blob/ecology/topics/ecology/tutorials/regionalGAM/Multispecies_tutorial.md). If you are interested in learning more about these specific cases, do not hesitate to follow these two tutorials in complement. 

> {: .comment}

> ### Agenda
> In this tutorial, we will cover:
1. Pre-processing
> {:pre-processing}
2. Making sure the dataset concerns only one species
> {:Making sure the dataset concerns only one species} 
3. Displaying the occurrence of the chosen species through the years
> {:Displaying the occurence of the chosen species through the years}

# Step 1: Pre-processing

The goal of the first step is to upload and prepare the file so that it will be usable for the *regional GAM* analysis (See [this warning](#inputdatawarning) for more information about the input file).
First of all, you will have to upload the files on Galaxy-E and then you might have to use some data handling tools to be able to use *regional GAM* tools.

>  ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial and give it a proper name.
> 2. Import the following file from [Zenodo](https://zenodo.org/record/1324204#.W2BmRn7fNE4) or from a data
>    library named `regional GAM data tutorial`
>
>    ```
>    CSV dataset with only one species:
>    https://zenodo.org/record/1324204/files/regional%20GAM%20data.csv?download=1
>    ```
>   
> ### {% icon tip %} Tip: Importing data via links
>    > 1. Search for the tool `Téléverser un ou plusieurs fichiers de votre ordinateur ou d'un serveur distant`
>    > 2. To import the dataset, you have two options:
>    > * If you have uploaded the file on your computer, drop it in the box provided for that purpose.
>    >    * Press **Start** and **Close** the window
>    > * If you have copied the link location:
>    >    * Select **Paste/Fetch data**
>    >    * Paste the link into the text field
>    >    * Choose the type: CSV
>    >    * Press **Start** and **Close** the window

>    {: .tip}

> ### {% icon comment %} Comment
>
> ⚠️ <a name="inputdatawarning"></a>Please note that the file must contain a header corresponding to: ```"SITES","SPECIES","YEAR","MONTH","DAY","COUNT"```, and that all the non numeric content must be between double quotes as "x" and that separators have to be ",". 

> {: .comment}

>    > ## <a name="resampling"></a> Re-sampling.  
When the dataset contains many details, it lengthens the file processing time therefore it can be very useful to learn how to hide the informations you don't need. For example, the list of SITE (look at the column with header `SITE)`  of the dataset you are using is really long and the SITES are classified into sub-sites. Here, we will assume that your file doesn't really need be as precise and this is the reason why you have to specify you don't want the sub-sites. To create a new "down-sampled" file, you can follow these steps:   

> ### {% icon hands_on %} Hands-on: hiding some informations

>    > 1. Use the `CSV to tabular` tool to first create a tabular file from your csv one (with only one species). This is a mandatory step as further tools are only working on tabular files!
>    > 2. Search for the tool `trouver et remplacer des patterns dans des colonnes` on the file on CSV with the following  parameters.
>    >  * Select the input file & the column with the SITE header.
>    >  * Click on `insert checks`
>    >  * "Trouver l'expression suivante": `(\.[0-9]+)` which specifies that you don't want the sub-sites (all suites of digits following a "." character) to be taken into account.
>    >  * "Remplacement": leave it empty.

>   > ### {% icon question %} Questions
>   >
>    > 1. After having successfully deleted the sub-sites informations, can you look at the original dataset and this new one and say how many sites you had, and you have now? You will maybe need to use tools like `Compter le nombre d'occurrences de chaque enregistrement`?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>The dataset contains 5 sites now against 1143 before "down-sampling". </li>
>    >    </ol>
>    >    </details>
>    {: .question}

# Step 2: Making sure the dataset concerns only one species 
 
The second step of any Regional GAM data analysis is making sure to have a dataset of only one specific species that you will then be able to use. If you want to create a graph showing abundance evolution by years of several species, you will have to superimpose the graphs on one another. 

> ### {% icon hands_on %} Hands-on: How many species are taken into account in this dataset
>
> As the dataset is quite big and may countain heterogeneous informations, you want to know wether the data are about one species or more. 
> 1. Search for the tool `compter le nombre d'occurence de chaque enrégistrement`with the following parameters:
> * "Sur le jeu de données": `output`from **tabular to CSV**
> * "Compter les occurrences des valeurs présentes dans la(les) colonne(s)": `column 1`
> * "Délimité par": `tabulation`.
> * "Comment les résultats doivent t'ils être triés ?": `Avec la valeur la plus présente en premier`.
> 2. Inspect the file by clicking on the `eye` icon to check that the dataset is on one species only.
> 3. Now you can regenerate a CSV file using the `tabular to CSV` tool on the ouptut from **trouver et remplacer des patterns dans des colonnes**.

> ### {% icon comment %} Comment

❗In case the dataset contains informations on more than on species, you can follow the multispecies tutorial on regionalGAM which is a complement to this one from the second step. If you are interested in doing so, please click on the following link: [Selectionning one specific species](https://github.com/Claraurf/training-material/blob/ecology/topics/ecology/tutorials/regionalGAM/Multispecies_tutorial.md#selectionningonespecificspecies)

{: .hands_on}


# Step 3: <a name="displayingtheoccurrenceofthespecies"></a> Displaying the occurrence of the chosen species through the years  
 
 
Now you have a file containing all the data on the species of interest. The main goal of this step is basically to create a material that can be used to generate charts. What you could also do, for example, would be to compare the evolution of various species through the years in the same site. You would have to superpose the different graphs on one another.
>

> ### {% icon hands_on %} Hands-on: Phenology

This step will allow you the show the phenology of a species and then to create charts representing it. In the second part, you will learn that it is possible to show the phenology of various species on a single chart allowing to compare them and analyse them more easily. 


> 1. Search for the tool `flight curve` and execute it specifying the following parameters: 
> * "Fichier de comptage": `output` from **tabular to CSV**.
>
> 🔹 Based on the `output` from **flight curve**, you can create a line chart which shows the occurence of the species through the years on a very visual material 

![Phenology chart](https://raw.githubusercontent.com/Claraurf/training-material/ecology/topics/ecology/tutorials/regionalGAM/Images/Phenology%20chart%20.png "This shows the occurrence of Pyronia tithonus")

If you want to access the chart on an interactive interface, you can click on the following link: [Chart on Galaxy](http://openstack-192-168-100-19.genouest.org/plugins/visualizations/charts/saved?id=d413a19dec13d11e)
>
>    > ### {% icon tip %} Visualiser. 
>    > 1. Click on the `Output` dataset from **flight curve** 
>    > 2. Click on the "Visualize" button then select `Charts` 
>    > 3. Give it a proper name
>    > 4. Select a visualization: "line chart (NVD 3)" 
>    > 5. Click on {% icon tip %} Select data:
>    > * "Provide a label": `The name of the species`
>    > * "Pick a series color": Choose a color for the line 
>    > * "Data point labels": `Column 1`
>    > * "Values for x-axis": `Column 2`
>    > * "Values for y-axis": `Column 6`
>    > 6. Click on {% icon tip %} Customize:
>    > * "X-Axis label": `Year`
>    > * "Y-Axis label": `nm values`
>    > 7. Click on {% icon tip %} `Visualize`
>    > 8. Click on {% icon tip %} `save this visualization`if you are willing
>
>
> ### {% icon comment %} Comment
>
> ⚠️ Please note, that if you want your chart to be more precise and to specify that the x-axis coresponds to "Week and year", it is possible. In order to do so, follow the tip below:
>    > ### {% icon tip %} Tip: Creating a new column of the dataset containing the week and the year 
First of all, you have to know how many years are taken into account in your dataset.
>    > 1. Search for the tool `Compter le nombre d'occurrences de chaque enregistrement` with the following parameters 
>    > * "Sur le jeu de données": `output` from **trouver et Remplacer des patterns dans des colonnes en utilisant des expressions régulières**.
>    > * "Select": `Column 3` (the on headed with `SITE`)
>    > * "Délimité par": `Tabulation`.
>    > * "Comment les résultats doivent t'ils être triés ?": `Par les valeurs comptées`.
>    > 2. Inspect the file by clicking on the `eye` icon to check how many years are taken into account.
>    > 3. Search for the tool`Trouver et Remplacer des patterns dans des colonnes en utilisant des expressions régulières` with the following parameters:
>    > * "Selectionner les cellules à partir de": `output` from **flight curve**.
>    > * "la colonne": `Column 2` (corresponding to the one headed with `YEAR`)
>    > * Click on `Insert check`:
>    >   * "Trouver l'expression suivante": `(20[0-9][0-9])`
>    >   * "Remplacement": `-\1` 
>    > 5. Inspect the file by clicking on the `eye` icon to check if all the years are now written with a "-" before the digits. 
>    > 6. Search for the tool `Merger des colonnes ensemble` with the following parameters:
>    > * "Selection du jeu de données": `output` from the last **Trouver et Remplacer des patterns dans des colonnes en utilisant des expressions régulières**.
>    > * "Merger les colonnes": `Column 3`(corresponding to the one headed with `WEEK`)
>    > * "avec les colonnes": `Column 2`(corresponding to the one headed with `YEAR`)

With the `output` from **Merger des colonnes ensemble** you can now generate a new chart which will have a x-axis corresponding to your column `Column "week""year"`.
>    > ### {% icon tip %} Visualiser
>    > 1. With the  `output` from **Merger des colonnes ensemble**.
>    > 2. Select `Charts`
>    > 3. Give it a proper name
>    > 4. Select a visualization: "line chart (NVD 3) 
>    > 5. Click on {% icon tip %} Select data:
>    > * "Provide a label": `The name of the species`
>    > * "Pick a series color": Choose a color for the line 
>    > * "Data point labels": `Column 1`
>    > * "Values for x-axis": `Column 7`
>    > * "Values for y-axis": `Column 6`
>    > 6. Click on {% icon tip %} Customize:
>    > * "X-Axis label": `Week-Year`
>    > * Y-Axis label: `nm values`
>    > 7. Click on {% icon tip %} `Visualize`
>    > 8. Click on {% icon tip %} `save this visualization`if you are willing to save to chart.
>   
> {: .comment}

> ⚠️ Please note that it is possible to show the occurrences of more than one species on a single chart. If you are interested in doing so, you should click here : [Various species chart explanations](https://github.com/Claraurf/training-material/blob/ecology/topics/ecology/tutorials/regionalGAM/Multispecies_tutorial.md#variousoccurencesonasinglechartexplanations)

 
> ###  <a name="Abundanceindex"></a>Abundance per year and per site

This will allow you to create a file showing the abundance per year of a chosen species in a certain site. Based on this file you will then learn how to represent this abundance on a chart. 
>
> 1. Look for the tool `Abundance index` with the following parameters:
> * "Fichier de comptage": `output` from **tabular to CSV**.  
> * "Flight curve output": `output` from **flight curve**.


> 🔹 Based on the  `output` from **abundance index**, we can create a chart showing the annual abundance trend of a certain species per site. 
>    > 1. Select `Charts`
>    > 2. Give it a proper name
>    > 3. Select a visualization: "Bar diagram (NVD 3)" 
>    > 4. Select data 
>    > * "Data point labels": `Column 1` 
>    > * "Values for x-axis": `Column 3`
>    > * "Values for y-axis": `Column 4`
>    > 5. Visualize
>    > 6. Click on {% icon tip %} `save this visualization`if you are willing to keep it

>   > ### {% icon question %} Questions
>   >
>    > 1. What do you think about this visualization? Maybe not so good? Search a way to display the content of the file using charts in a more accurate manner... To do so, you can use tools like **Trouver et Remplacer des patterns dans des colonnes en utilisant des expressions régulières**, **Merger des colonnes ensemble** and **Trier les données dans un ordre ascendant ou descendant**
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>You can use the **Trouver et Remplacer des patterns dans des colonnes en utilisant des expressions régulières** tool to first replace `(20[0-9][0-9])` on the column 3 by `-\1` then on the result of this tool execution, replace `"` by nothing on the column 1. Furthermore, you can merge column 1 and column 3 of the resulting dataset. Finally, you can sort the new dataset by column 1 (alphabetical/ascending) and column 3 (alphabetical/ascending). </li>
>    >    </ol>
>    >    </details>
>    {: .question}

{: .hands_on}

> ### Expected temporal trend

The expected temporal trend allows you to have an overview of the evolution of a species in a certain type of environment in the futur.

> ### {% icon hands_on %} Hands-on: Expected temporal trend
>    > 1. Look for the tool `Expected temporal trend` with the the following parameters: 
>    > * "Fichier tabulé, produit par l'outil ab_index": `output` from **abundance index**.
>    
![Expected temporal trend](https://raw.githubusercontent.com/Claraurf/training-material/ecology/topics/ecology/tutorials/regionalGAM/Images/Expected%20temporal%20trend.png "This shows the expected evolution of Aglais io")

> ⚠️ Please note that sometimes the expected temporal trend can't be done on dataset. If you want this action to work, the occurences on your dataset must lie between the month of April and the end of the month of September.

Note also that you will obtain two files resulting of the action above. The first one will be the graph and the second one will contains the values of "x".

> ### Linear regression 

The point of doing a linear regression is to determinate if the year has an influence on the abundance of a species. 

>    > 1. Look for the tool `linear regression` with the following parameters.
>    > * "Fichier produit par l'outil glmmpql/Expected temporal trend": `output 2` from **temporal trend**. 
>    > * "Fichier produit par l'outil ab_index": `output` from **abundance index**.

{: .hands_on}
 
# Conclusions  

{:.no_toc}

In this tutorial, you have analyzed regional GAM data to extract useful informations in order to be able to show different tendencies of a chosen species. Therefore, you are now able to treat the dataset so that it shows only the data concerning one specific species of your choice. From there, you can show the occurrence of this species through the years first on a dataset and then on a visual chart. You have also learned how to represent on a single chart the occurences of various species. Afterwards, we have shown you how to create a dataset containing the informations on the abundance of a species per year and per site. Based on which you can henceforth visually represent the annual abundance trend on a chart. Thereafter, you have the possibility of showing the expected temporal trend, based on which you will be able to try predicting the future evolution a given species. The last part of this tutorial has shown you how to calculate the linear regression allowing you to determinate wether the year has an influence on the abundance of a species or not. 
