---
layout: tutorial hands on
topic_name: ecology
tutorial_name: Regional GAM
---

# Introduction
{:.no_toc}


This tutorial will show how to study species phenology through the computation of abundance index and trends. It will explain you how to use different [regionalGAM](https://github.com/RetoSchmucki/regionalGAM) tools on Galaxy-E allowing you to deal with datasets containing occurences informations for various species per site and per date through a couple of years.
After a certain numbers of steps, you will be able to extract single species data and study related phenology through the years. The goal of this exercise is to be able to create abundance trend over time and biodiversity indicators. You could for example try to predict the occurences of one specific species in a certain type of environnement using the prediction model of climate evolution. Based on charts that you will generate, you could try to explain the evolution of a species with environmental data (temperatures variations, modifications of the environmental conditions).
You will basically learn how to create a file on the basis of which you can create a visual material that can be quite easily understood and therefore be efficient for a large audience.


> ### Agenda
> In this tutorial, we will cover:
1. Pre-processing
> {:pre-processing}
2. Selectionning one specific species and show all corresponding data
> {:selectionning one specific species and show all corresponding data}
3. Displaying the occurence of the chosen species through the years
> {:Displaying the occurence of the chosen species through the years}
4. Conclusion 
> {:conclusion}

# Step 1: Pre-processing

The goal of the first step is to upload and prepare the file so that it will be usable for the *regional GAM* analysis.
First of all, you will have to upload the files on Galaxy-E and then you maybe will have to use some specific formating tools to be able to use *regional GAM* tools.

>  ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial and give it a proper name.
> 2. Import the following files from [Zenodo](https://zenodo.org/record/1324204#.W2BmRn7fNE4) or from a data
>    library named `regional GAM data tutorial`
>
>    ```
>    CSV dataset with only one species:
>    https://zenodo.org/record/1324204/files/regional%20GAM%20data.csv?download=1
>
>    RData dataset with only one species:
>    https://zenodo.org/record/1324204/files/gatekeeper_CM%20.RData?download=1csv
>
>    CSV dataset with several species: 
>    https://zenodo.org/record/1324204/files/Dataset%20multispecies%20Regional%20GAM.csv?download=1
>    ```
>
> ### {% icon tip %} Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * Press **Start** and **Close** the window
>    {: .tip}
>
> ### {% icon comment %} Comment
>
> ⚠️ Please note that the file must contain a header corresponding to: ```"SITES","SPECIES","YEAR","MONTH","DAY","COUNT"```, and that all the non numeric content must be between double quotes as "x" and that separators have to be ",". 
> 
>  Note also that the first file mentionned above is on the RData format and if you choose to use it, you will have to process this binary file to obtain an appropriate CSV dataset. To do so, you can use the following tools:
>   > * Search for the tool `RData binary file reader`with the following parameters:
>   >      * "Rdata binary file to explore": "dataset on RData" 
>   > * Search for the tool `RData parser` with the following parameters:
>   >      * "Rdata file to explore": "dataset on RData"
>   >      * "File with .Rdata content details": file of **`RData binary file reader`**
>   >      * "Select which attribute(s) you want to extract": select everything but "trend"
>   >      * ⚠️ Please note that the tool `RData parser` creates separated files, each of them containing one column. The file with the "TREND" header can be let aside as we don't need it for what will follow.
>   > * Search for the tool `Coller deux jeux de données l'un à côté de l'autre` to create a file comporting all the data     required with the following parameters:
>   >      * "Coller":  file of **RData parser** headed with "SPECIES"
>   >      * "et": file of **RData parser** with headed with "SITE"
>   >      * Repeat `Coller deux jeux de données l'un à côté de l'autre` as many times as there are separated files in order to create a final dataset with all the columns. First you must paste 2 columns together, then you must paste this last file with a third column and do this action again and again until your final file countains all the columns. 
>   >      * Repeat `Coller deux jeux de données l'un à côté de l'autre` pasting the file containing 2 columns with the one headed by "YEAR".
>   >      * Repeat `Coller deux jeux de données l'un à côté de l'autre`pasting the file containing 3 columns with the one headed by "MONTH". 
>   >      * Repeat `Coller deux jeux de données l'un à côté de l'autre`pasting the file containing 4 columns with the one headed by "DAY".
>   >      * Repeat `Coller deux jeux de données l'un à côté de l'autre`pasting the file you containing 5 with the one headed by "COUNT". 
>   >  
>   > {: .comment}
>
>   > ### {% icon question %} Questions
>   >
>   > 1. Why do you need to use `Coller deux jeux de données l'un à côté de l'autre` 
>   >
>   >    <details>
>   >    <summary>Click to view answers</summary>
>   >    <ol type="1">
>   >    1. Because you want to create a single dataset which countains all the data on a species. You decided to upload a dataset on RData and therefore you had to use `RData binary file reader` and `RData parser`. This last tool treats the file and allows you to open it on Galaxy-E but it creates as many files as there are columns (when RData object is composed from a unique data table). This is the reason why you had to carry out on a set of actions ending by the creation of one complete file. 
>    >    </details>


>    > ## Re-sampling. 
When the dataset contains many details, it lengthens the file processing time therefore it can be very useful to learn how to hide the informations you don't need. For example, the list of SITE of the dataset you are using are is really long and the SITES are classified into sub-sites. Your file doesn't really need be as precised and this is the reason why you have to specify you don't want the sub-sites to be considered in order to create a new  file that you will be able to use.   

>    > 1. Search for the tool `trouver et remplacer des patterns dans des colonnes` on the file on CSV with the following  parameters.
>    >  * Click on`"insert checks"`
>    >  * "Trouver l'expression suivante": `"(\.[0-9]+)"` which specifies that you don't want the sub-sites to be taken into                                    >          account.
>    >  * "Remplacement":`"leave it empty"`.
>    > 3. Search for the tool `tabular to CSV`and select the file of **trouver et remplacer des patterns dans des colonnes**.
>    
>
>

# Step 2: Selectionning one specific species and showing all the data corresponding to it


The second step of any Regional GAM data analysis is making sure to have one dataset of only one specific species that you will then be able to use. Because a graph is occurent only if it shows the occurence of one species trough the years. If you want to compare this evolution with the one of another species, you will have to superimpose the graphs on one another. 



> ## How many species are taken into account in this dataset
> As the dataset is important and countains manw informations, you want to know wether the data are about only one species or  more.
> 1. Search for the tool `compter le nombre d'occurence de chaque enrégistrement`with the following parameters.
> * "Compter les occurrences des valeurs présentes dans la(les) colonne(s)": `column 1`
> * "Délimité par": `tabulation`.
> * "Comment les résultats doivent t'ils être triés ?": `Avec la valeur la plus présente en premier`.
> 2. Inspect the file by clicking on the `eye` icon to check how many species are taken into account.
> > If there is ony one species you can't skip the following steps ang go directly convert your file using the tool `tabular  to CSV`
>
>    > ### Create a new file concerning only the data of one species
>    > 1. Copy the name of the species you are interested in from the file on CSV (for example: "Pyronia tithonus").
>    > 2. Search for the tool`filtrer des données dur une colonne en utilisant des expressions simples`with the following   parameters.
>    > * En utilisant la condition suivante: `"c1=='habitat2'"` replacing 'habitat2' with the name of the species (for example: `"c1=='"Pyronia tithonus"'"`  instead of c1=='habitat2')  
>    > * Nombre de lignes d'en-tête à passer: `"0"`.
>    > * You can repeat this set of actions as much as is necessary changing only the name of the species taken into account.
>    > 3. Search for the tool`tabular to CSV` with the following parameters 
>    > * Choose the file you've just created 
>    > * Separators: `","`.
>    > 4. Repeat this last operation on all the differents files on different species. 
>    
>



 # Step 3: Displaying the occurence of the chosen species through the years
 
Now that you have a file containing all the data on the species you chose, this step is goint to explain to you what you can do with it and how useful it can be. The main goal of this step is basically to create a material that can be used to generate charts.  What you could also do, for example, would be to compare the evolution of various species trough the years in the same site. You would have to superpose the different graphs on one another. The charts have the benefit of being very visual and to be easier to interpret than a dataset.
>
> ### Create a flight curve 

> 1. Search for the tool `flight curve`
> * Select the file on CSV with the data on one species
>
> Based on the flight curve, you can create a line chart which shows the occurence of the species through the years on a very visual material 

![Flight curve chart](https://github.com/Claraurf/training-material/blob/ecology/topics/ecology/tutorials/regionalGAM/Images/Phenology%20Regional%20GAM.png)

If you want to access the chart on an interactive interface, you can use click on the following link: [Chart on Galaxy](http://openstack-192-168-100-19.genouest.org/plugins/visualizations/charts/saved?id=e85a3be143d5905b)
>
>    > ### {% icon tip %} Visualiser 
>    > 1. Select `Charts`
>    > 2. Give it a proper name
>    > 3. Select a visualization: "line chart (NVD 3) 
>    > 4. Select data 
>    > * "Data point labels": `"Column 1"`
>    > * "Values for x-axis": `"Column 2"`
>    > * "Values for y-axis": `"Column 6"`
>    > 5. Visualize
>    > 6. Click on {% icon tip %} `save this visualization`if you are willing to keep it
>
>
> ### Abundance per year and per site 
>
> 1. Look for the tool `Abundance index` with the following parameters:
> * "Fichier de comptage": onespecies dataset on CSV.  
> * "Flight curve output": flight curved obtained .


> Based on the abundance index, we can create a chart showing the anual trend abundance of a certain species per site. 
>    > 1. Select `"Charts"`
>    > 2. Give it a proper name
>    > 3. Select a visualization: "Bar diagram (NVD 3)" 
>    > 4. Select data 
>    > * "Data point labels": `"Column 1"` 
>    > * "Values for x-axis": `"Column 3"`
>    > * "Values for y-axis": `"Column 4"`
>    > 5. Visualize
>    > 6. Click on {% icon tip %} `save this visualization`if you are willing to keep it


> ### Expected temporal trend

The expected temporal trend allows you to have an overview of the evolution of a species in a certain type of environment in the futur.

![Expected temporal trend](https://github.com/Claraurf/training-material/blob/ecology/topics/ecology/tutorials/regionalGAM/Images/Expected%20temporal%20trend.png)

>    > 1. Look for the tool `Expected temporal trend` with the the following parameters: 
>    > * "Fichier tabulé, produit par l'outil ab_index": `output` of **abundance index**.
>    
> ⚠️ Please note that sometimes the expected temporal trend can't be done on dataset. If you want this action to work, the occurences on your dataset must lie between the month of April and the end of the month of September.

> ### Linear regression 

The point of doing a linear regression is to determinate if the year has an influence on the abundnce of a species. 

>    > 1. Look for the tool `linear regression` with the following parameters.
>    > * "Fichier produit par l'outil glmmpql/Expected temporal trend": `output 2` of **temporal trend**. 
>    > * "Fichier produit par l'outil ab_index": `output` of **abundance index**.
 
# Conclusions  

{:.no_toc}

In this tutorial, we have analyzed regional GAM data to extract useful informations in order to be able to show different tendencies of a chosen species. 
