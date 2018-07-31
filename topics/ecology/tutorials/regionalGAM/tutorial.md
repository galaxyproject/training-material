---
layout: tutorial_hands_on
topic_name: ecology
tutorial_name: Regional GAM
---

# Introduction
{:.no_toc}

<!-- This is a comment. -->


This exercise uses the public dataset of Reto Schmucki [regional GAM](https://zenodo.org/record/1321885/files/gatekeeper_CM%20.RData?download=1). It is a file recording the presence of various species per site and per days through a couple of years. 
The goal of this exercise is to be able to create biodiversity indicators and abundance trend over time based on this dataset. Using different tools we will show the occurrence of one specie through the years. 
You could for example try to predict the occurence of one specific species in a certain type of environnement using the prediction model of climate evolution.
You will basically learn how to create a file with which you can create a visual material that can be quite easily understood and therefore be efficient for a large audience.

> ### Agenda
>
1. Pre processing.
2. Selectionning one specific species and show all the data corresponding to it.
3. Show the occurence of the chosen species through the years.
> {:toc

>
>
# Step 1: Pre-processing

The goal of the first step is to upload and prepare the file so that it will be usable for the regional GAM analysis.
First of all the dataset must be uploaded.

> ### {% icon hands_on %} Hands-on: Upoading dataset
>
> 1. Create a new history for this tutorial and give it a proper name.
> 2. Then you can either import `gatekeeper_CM .RData` from [Zenodo](https://zenodo.org/record/1321885/files/gatekeeper_CM%20.RData?download=1) or either import `regional GAM data.csv`from [Zenodo](https://zenodo.org/record/1321885/files/regional%20GAM%20data.csv?download=1).

> ### {% icon tip %} Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * Press **Start** and **Close** the window
>    > * Click on the `pencil` icon once the file is imported
>    > * Click on **Datatype** in the middle panel
>    > * Select `fastqsanger` as **New Type**
>    {: .tip}
>
>    > ### {% icon tip %} Tip: Importing data from a data library
>    >
>    > * Go into "Shared data" (top panel) then "Data libraries"
>    > * Click on "Training data" and then "Analyses of ChIP-Seq data"
>    > * Select interesting file
>    > * Click on "Import selected datasets into history"
>    > * Import in the history
>    {: .tip}

>   >⚠️ Please note that the file must contain the headers "SITES","SPECIES","YEAR", "MONTH","DAY","COUNT", that all the non numerum content must have "x" and that separators have to be ",".
>   > * 1. Note also that he first file mentionned above is on RData and if you choose to use it, you will have to use the following tools:
>   > * `RData binary file reader`. 
>   > * `RData parser`.
>   > * Search for the tool `Coller deux jeux de données l'un à côté de l'autre` to create a file comporting all the data     required.
>   >  > * Repeat `Coller deux jeux de données l'un à côté de l'autre`and select the file on "SPECIES" and the file on "SITE".   >   >  > * Repeat `Coller deux jeux de données l'un à côté de l'autre`  and select the file obtained above and the one on                 "YEAR". 
>   >  > * Repeat `Coller deux jeux de données l'un à côté de l'autre`select the file obtained above and the one on "MONTH".
>   >  > * Repeat `Coller deux jeux de données l'un à côté de l'autre`select the file obtained above and the one on "DAY".
>   >  > * Repeat `Coller deux jeux de données l'un à côté de l'autre`select the file obtained above and the one on "COUNT".
> 
>   >  > ⚠️ Please note that we left the file on "TREND" aside on purpose because we will be able to show it afterwards using some tools.
>
>
>
>
> ### {% icon question %} Questions
>
> 1. Why do you need to use `Coller deux jeux de données l'un à côté de l'autre` 
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Because you want to create a single file which countains all the data. You decided to upload a dataset on RData and therefore you had to use `RData binary file reader` and `RData parser`. This last tool treats the file and allows you to open it on Galaxy-E but it creates as many files as there are columns.
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}




>
>
>    > ## Re-sampling. 
>    > The list of SITES of the dataset you are using are is really long and the SITES are classified into sub-categories.    >    > Your file doesn't really need  to contain all these sub-sites because it lengthens the file processing time on Galaxy-     >    > E.This is the reason why you have to specify you don't want the sub-sites to be considered in order to create a new     >    > file that you will be able to use.   
>    >
>    >
>    > > 1. Search for the tool `trouver et remplacer des patterns dans des colonnes` on the file on CSV with the following   >        parameters.
>    > > * Click on`"insert checks"`
>    > > * Trouver l'expression suivante: `"(\.[0-9]+)"` which specifies that you don't want the sub-sites to be taken into                                    >          account.
>    > > * Remplacement:`"leave it empty"`.
>    > > 3. Search for the tool `tabular to CSV`and select the file obtained above.
>    
>
>

# Step 2: Selectionning one specific species and showing all the data corresponding to it

The second step of any Regional GAM data analysis is making sure to have one dataset of only one specific specie that you will then be able to use. Because a graph is occurent only if it shows the occurence of one species trough the years. If you want to compare this evolution with the one of another species, you will have to superimpose the graphs onto the others. 
>
>    
> ## Count on many species are taken into account in the dataset your are using 
> As the dataset is important and countains manw informations, you want to know if the data are about only one species or  >more.
> 1. Search for the tool `compter le nombre d'occurence de chaque enrégistrement`with the following parameters.
> * Compter les occurrences des valeurs présentes dans la(les) colonne(s): select: `"column 1"`
> * Délimité par: `"tabulation"`.
> * Comment les résultats doivent t'ils être triés ?: `"Avec la valeur la plus présente en premier"`.
> 2. Inspect the file by clicking on the `eye` icon to check how many species are taken into account.
> > If there is ony one species you can't skip the following steps ang go directly convert your file using the tool `tabular t    to CSV`
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
>
>
>
 # Step 3: Show the occurrence of the chosen species through the years
>
> ### Create a flight curve 
>
> 1. Search for the tool `flight curve`
> * Select the file on CSV with the data on one species
>
> Based on the flight curve, you can create a line chart which shows the occurence of the species through the years on a very visual material 
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
> ### Regional GAM per year and per site 
>
> 1. Look for the tool `Abundance index`
> * "Fichier de comptage": tabular on CSV corresponding to the data of one species. 
> * "Flight curve output": flight curved obtained above.
>
>
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
>    > 1. Look for the tool `Expected temporal trend`
>    > * Select the file obtained afer the abundance index.
>    
> ⚠️ Please note that sometimes the expected temporal trend can't be done. It happens that the data are to random or that to many parameters must be taken into account and therefore it isn't possible for the software to estimate the expected temporal trend.



# Step 4: Show the occurrence of the chosen species through the years



