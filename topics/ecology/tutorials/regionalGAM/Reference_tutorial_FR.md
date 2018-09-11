---
layout: tutorial_hands_on
topic_name: ecology
tutorial_name: Regional GAM
---

# Introduction
{:.no_toc}

érentes espèces par site et par date. années.
Après un certain nombre d'étapes, vous pourrez extraire des données sur une seule espèce et étudier la phénologie associée au fil des ans. Le but de cet exercice est de pouvoir créer une tendance de l’abondance au fil du temps et des indicateurs de biodiversité. Ces indicateurs permettent de suivre les tendances en termes de dynamique de la population. Vous pourriez par exemple essayer de prédire les occurrences d'une espèce spécifique dans un certain type d'environnement en utilisant le modèle de prédiction de l'évolution du climat. Sur la base de graphiques que vous allez générer, vous pouvez essayer d’expliquer l’évolution d’une espèce avec des données environnementales (variations de température, modifications des conditions environnementales).
Vous apprendrez essentiellement comment créer un fichier sur la base duquel vous pourrez créer un matériel visuel facilement compréhensible et donc efficace pour un large public.
> ### {% icon comment %} Comment
>
> ⚠️ Veuillez noter qu'il existe deux autres didacticiels sur regionalGAM, l'un spécifique à [l'utilisation du jeu de données d'entrée RData] (https://github.com/Claraurf/training-material/blob/ecology/topics/ecology/tutorials/regionalGAM/ Rdata_tutorial.md) et l’autre vous explique comment gérer un [ensemble de données multispécifiques] (https://github.com/Claraurf/training-material/blob/ecology/topics/ecology/tutorials/regionalGAM/Multispecies_tutorial.md). Si vous souhaitez en savoir plus sur ces cas spécifiques, n'hésitez pas à suivre ces deux didacticiels en complément.
> {: .comment}

> ### Agenda
> Dans ce tutoriel, nous couvrirons:
1. le "Pré-processing" des données
> {:pre-processing}
2. Vérifier que les données ne concernent bien qu'une seule espèce
> {:Making sure the dataset concerns only one species} 
3. Traiter les occurences d'une espèce à travers les années et les sites
> {:Displaying the occurence of the chosen species through the years}

# Step 1: Pré-processing

Le but de la première étape est de télécharger et de préparer le fichier afin qu'il soit utilisable pour l'analyse * regionalGAM * (Voir [cet avertissement] (# inputdatawarning) pour plus d'informations sur le fichier d'entrée).
Tout d'abord, vous devrez télécharger les fichiers sur [Galaxy-E] (https://openstack-192-168-100-19.genouest.org) et vous devrez peut-être utiliser des outils de traitement des données pour pouvoir utiliser * les outils regionalGAM *.
>  ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Créez un nouvel historique pour ce tutoriel et donnez-lui un nom propre comme «Tuto training regionalGAM».
> 2. Importez le fichier suivant à partir de [Zenodo] (https://zenodo.org/record/1324204#.W2BmRn7fNE4) ou à partir d'une donnée
> bibliothèque nommée `regional GAM data tutorial`
>
>    ```
>    Fichier au format CSV avec une seule espèce :
>    https://zenodo.org/record/1324204/files/regional%20GAM%20data.csv?download=1
>    ```
>   
> ### {% icon tip %} Tip: Importaion via hyperliens
>    > 1. Cherchez l'outil `Téléverser un ou plusieurs fichiers de votre ordinateur ou d'un serveur distant`
>    > 2. Pour importer le jeu de données :
>    > * Sélectionnez **Paste/Fetch data**
>    > * Coler l'hyperlien dans le champ vide
>    > * Spécifier le format : CSV
>    > * Presser le bouton **Start** puis **Close** pour lancer le téléversement et fermer la boîte de dialogue.

>    {: .tip}

> ### {% icon comment %} Comment
>
> ⚠️ <a name="inputdatawarning"></a>Veuillez noter que le fichier doit contenir un en-tête correspondant à : ```"SITES","SPECIES","YEAR","MONTH","DAY","COUNT"```, and that all the non numeric content must be between double quotes as "x" and that separators have to be ",". 

> {: .comment}

>    > ## <a name="resampling"></a> Re-sampling.  
Lorsque le jeu de données contient de nombreux détails, cela rallonge le temps de traitement des fichiers. Il peut donc être très utile d'apprendre à masquer les informations inutiles. Par exemple, la liste de SITE (consultez la colonne avec l'en-tête `SITE) du jeu de données que vous utilisez est très longue et les sites SITE sont classés en sous-sites (comme ESBMS.12, ESBMS.28, ESBMS.55, ...). Ici, nous supposerons que votre fichier n'a pas vraiment besoin d'être aussi précis et c'est la raison pour laquelle vous devez spécifier que vous ne voulez pas les sous-sites. Pour créer un nouveau fichier "échantillonné" (supprimez donc les mentions ---. 12, ---. 28), vous pouvez suivre ces étapes:
> ### {% icon hands_on %} Hands-on: Dégrader le niveau d'information du jeu de données.

>    > 1. Utiliser l'outil `CSV to tabular` pour créer un fichier tabulé à partir du CSV (contenant les information sur une seule espèce). Cette étape est obligatoire car d'autres outils ne fonctionnent que sur des fichiers tabulés.!
>    > 2. Chercher l'outil `trouver et remplacer des patterns dans des colonnes` et l'exécuter sur le fichier CSV en spécifiant les paramètres suivants :
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
> 1. Search for the tool `compter le nombre d'occurence de chaque enregistrement` with the following parameters:
> * "Sur le jeu de données": `output`from **CSV to Tabular**
> * "Compter les occurrences des valeurs présentes dans la(les) colonne(s)": `column 1`
> * "Délimité par": `tabulation`.
> * "Comment les résultats doivent t'ils être triés ?": `Avec la valeur la plus présente en premier`.
> 2. Inspect the file by clicking on the `eye` icon to check that the dataset is on one species only.
> 3. Now, as regionalGAM tools use CSV files as input, you can regenerate a CSV file using the `tabular to CSV` tool on the output from **trouver et remplacer des patterns dans des colonnes**. Please, tag your new dataset with an explicit tag as "Count and/or rename this dataset like "Counting file".

> ### {% icon comment %} Comment

❗In case the dataset contains informations on more than on species, you can follow the multispecies tutorial on regionalGAM which is a complement to this one from the second step. If you are interested in doing so, please click on the following link: [Selectionning one specific species](https://github.com/Claraurf/training-material/blob/ecology/topics/ecology/tutorials/regionalGAM/Multispecies_tutorial.md#selectionningonespecificspecies)

{: .hands_on}


# Step 3: <a name="displayingtheoccurrenceofthespecies"></a> Displaying the occurrence of the chosen species through the years  
 
 
Now you have a file containing all the data on the species of interest. The main goal of this step is basically to create a material that can be used to generate charts. What you could also do, for example, would be to compare the evolution of various species through the years in the same site. You would have to superpose the different graphs on one another.
>

> ### {% icon hands_on %} Hands-on: Phenology

This step will allow you to compute and display the phenology of a species. In the second part, you will learn that it is possible to show the phenology of various species on a single chart allowing to compare them and analyse them more easily. 


> 1. Search for the tool `flight curve` and execute it specifying the following parameters: 
> * "Fichier de comptage": `output` from **tabular to CSV**.
>
> 🔹 Based on the `output` from **flight curve**, you can create a line chart which shows species occurences through the years on a very visual material 


![Phenology chart](https://github.com/yvanlebras/training-material/blob/ecology/topics/ecology/tutorials/regionalGAM/Images/Phenology%20chart%20.png "This shows the occurrence of Pyronia tithonus")

If you want to access the chart on an interactive interface, you can click on the following link: [Chart on Galaxy-E](https://openstack-192-168-100-19.genouest.org/u/ylebras/v/pyronia-tithonus-phenology-2)
>
>    > ### {% icon tip %} Visualization. 
>    > 1. Click on the name of output dataset from **flight curve** (something like `Flight curve on data 6`) to expand dataset related details and options
>    > 2. Click on the "Visualize" button then select `Charts` 
>    > 3. Give it a proper name like `Pyronia tithonus phenology`
>    > 4. Select a visualization type: "line chart (NVD 3)" 
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
> ### {% icon comment %} Comments
>
> ⚠️ Please note, that if you want to create a "stacked" visualization, overlapping each year, you can use several executions (one execution by year) of the `Filtrer des données sur une colonne en utilisant des expressions simples` tool specifying the year you want on the `condition suivante` parameter, `c2==2004` for 2004, then `c2==2005` for 2005, etc... then you can paste all resulting files side by side using one or several executions of the `Coller deux jeux de données l'un à côté de l'autre` tool so you can specify on the "Select Data" tab of visualization, several Data series (one by year).
> ⚠️ Please note, that if you want your chart to be more precise and to specify that the x-axis corresponds to "Week and year", it is possible. In order to do so, follow the tip below:
>    > ### {% icon tip %} Tip: Creating a new column of the dataset containing the week and the year 
First of all, you have to know how many years are taken into account in your dataset.
>    > 1. Search for the tool `Compter le nombre d'occurrences de chaque enregistrement` with the following parameters 
>    > * "Sur le jeu de données": `output` from **trouver et Remplacer des patterns dans des colonnes en utilisant des expressions régulières**.
>    > * "Select": `Column 3` (the on headed with `YEAR`)
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
>    > 4. Select a visualization type: "line chart (NVD 3) 
>    > 5. Click on {% icon tip %} Select data:
>    > * "Provide a label": The name of the species like `Pyronia tithonus`
>    > * "Pick a series color": Choose a color for the line 
>    > * "Data point labels": `Column 6` (the nm column) or another one
>    > * "Values for x-axis": `Column 7` (the "week-year" column)
>    > * "Values for y-axis": `Column 6` (the nm column)
>    > 6. Click on {% icon tip %} Customize:
>    > * "X-Axis label": `Week-Year`
>    > * Y-Axis label: `nm values`
>    > 7. Click on {% icon tip %} `Visualize`
>    > 8. Click on {% icon tip %} `save this visualization` if you are willing to save to chart. You will obtain something like [that](https://openstack-192-168-100-19.genouest.org/u/ylebras/v/pyronia-tithonus-phenology-week-year)
>   
> {: .comment}

> ⚠️ Please note, that if you want to create a "stacked" visualization, overlapping each year, you cn use the same "tip" than at the previous step.
> ⚠️ Please note that it is possible to show the occurrences of more than one species on a single chart. If you are interested in doing so, you should click here : [Various species chart explanations](https://github.com/Claraurf/training-material/blob/ecology/topics/ecology/tutorials/regionalGAM/Multispecies_tutorial.md#variousoccurencesonasinglechartexplanations)

 
> ###  <a name="Abundanceindex"></a>Abundance per year and per site

This will allow you to create a file showing the abundance per year of a chosen species in a certain site. Based on this file you will then learn how to represent this abundance on a chart. 
>
> 1. Look for the tool `Abundance index` with the following parameters:
> * "Fichier de comptage": `output` from **tabular to CSV** (normally renamed "Counting file" and/or tagged "Count").  
> * "Flight curve output": `output` from **flight curve**.


> 🔹 Based on the  `output` from **abundance index**, we can create a chart showing the annual abundance trend of a certain species per site. 
>    > 1. Select `Charts`
>    > 2. Give it a proper name (`Pyronia tithonus abundance index ` for example)
>    > 3. Select a visualization type: "Bar diagram (NVD 3)" 
>    > 4. Select data 
>    > * "Data point labels": `Column 1` 
>    > * "Values for x-axis": `Column 3`
>    > * "Values for y-axis": `Column 4`
>    > 5. Customize 
>    > * "X-Axis label": `Year`
>    > * "Y-Axis label": `regional_gam`
>    > 5. Visualize
>    > 6. Click on {% icon tip %} `save this visualization`if you are willing to keep it

>   > ### {% icon question %} Questions
>   >
>    > 1. What do you think about this visualization? Maybe not so good? Search a way to display the content of the file using charts in a more accurate manner... To do so, you can use tools like **Trouver et Remplacer des patterns dans des colonnes en utilisant des expressions régulières**, **Merger des colonnes ensemble**, **Supprimer les premières lignes d'un fichier** , and **Trier les données dans un ordre ascendant ou descendant**
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>You can use the **Trouver et Remplacer des patterns dans des colonnes en utilisant des expressions régulières** tool to first replace `(20[0-9][0-9])` on the column 3 (the "YEAR" one) by `-\1` then on the result of this tool execution, replace `"` by nothing on the column 1. Furthermore, you can merge column 1 and column 3 of the resulting dataset. Finally, after deleting the first line (the header) with **Supprimer les premières lignes d'un fichier**, you can sort the new dataset by column 1 (alphabetical/ascending) and column 3 (alphabetical/ascending). </li>
>    >    </ol>
>    >    </details>
>    {: .question}
> 🔹 Based on the  new dataset, we can display a better chart showing the annual abundance trend of a certain species per site. 
>    > 1. Select `Charts`
>    > 2. Give it a proper name (`Pyronia tithonus abundance index ` for example)
>    > 3. Select a visualization type: "Bar diagram (NVD 3)" 
>    > 4. Select data 
>    > * "Data point labels": `Column 6` 
>    > * "Values for x-axis": `Column 6`
>    > * "Values for y-axis": `Column 4`
>    > 5. Customize 
>    > * "X-Axis label": `Site-Year`
>    > * "Y-Axis label": `regional_gam`
>    > 5. Visualize
>    > 6. Click on {% icon tip %} `save this visualization`if you are willing to keep it

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

You can then execute the **Linear regression ajusted** tool to take into account autocorrelation of model residuals .

Finally, a global trend (over years) can be displayed using the **Trend** tool.

{: .hands_on}
 
# Conclusions  

{:.no_toc}

In this tutorial, you have analyzed regional GAM data to extract useful informations in order to be able to show different tendencies of a chosen species. Therefore, you are now able to treat the dataset so that it shows only the data concerning one specific species of your choice. From there, you can show the occurrence of this species through the years first on a dataset and then on a visual chart. You have also learned how to represent on a single chart the occurences of various species. Afterwards, we have shown you how to create a dataset containing the informations on the abundance of a species per year and per site. Based on which you can henceforth visually represent the annual abundance trend on a chart. Thereafter, you have the possibility of showing the expected temporal trend, based on which you will be able to try predicting the future evolution a given species. The last part of this tutorial has shown you how to calculate the linear regression allowing you to determinate wether the year has an influence on the abundance of a species or not. 
