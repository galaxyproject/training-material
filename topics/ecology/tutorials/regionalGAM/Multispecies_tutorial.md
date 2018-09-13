---
layout: tutorial_hands_on
topic_name: ecology
tutorial_name: Regional GAM
---

# Introduction
{:.no_toc}

⚠️ You might be willing to follow this tutorial if you are interested in working on a multispecies file.

❗Please be aware that this tutorial is only a complement to [refence_tutorial on regionalGAM](training-material/topics/ecology/tutorials/regionalGAM/Reference_tutorial.md) and that therefore there are some missing parts. 
Follow the steps bellow and then when indicated, you will be redirected to the complete tutorial. 

This tutorial will show how to study species phenology through the computation of abundance index and trends. It will explain you how to use different [regionalGAM](https://github.com/RetoSchmucki/regionalGAM) tools on Galaxy-E allowing you to deal with datasets containing occurences informations for various species per site and per date through a couple of years.
After a certain numbers of steps, you will be able to extract single species data and study related phenology through the years. The goal of this exercise is to be able to create abundance trend over time and biodiversity indicators. Following these indicators allow to follow trends in terms of population dynamics. You could for example try to predict the occurences of one specific species in a certain type of environnement using the prediction model of climate evolution. Based on charts that you will generate, you could try to explain the evolution of a species with environmental data (temperatures variations, modifications of the environmental conditions).
You will basically learn how to create a file on the basis of which you can create a visual material that can be quite easily understood and therefore be efficient for a large audience.


> ### Agenda
> In this tutorial, we will cover:
1. Pre-processing
> {:pre-processing}
2. Selectionning one specific species and show all corresponding data
> {:selectionning one specific species and show all corresponding data}
3. Creating a chart showing more than one occurrence
> {:creating a chart showing more than one occurrence}

# Step 1: Pre-processing

The goal of the first step is to upload and prepare the file so that it will be usable for the *regional GAM* analysis (See [this warning](#inputdatawarning) for more information about the input file).
First of all, you need to use a Galaxy instance with related regionalGAM tools. You can deploy your own local instance through Docker as a  Galaxy flavour or use our [Galaxy-E test instance](https://openstack-192-168-100-96.genouest.org/). 
After uploading input files, you might have to use some data handling tools to be able to use *regional GAM* tools.

>  ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial and give it a proper name as `Tuto training regionalGAM Multispecies`.
> 2. Import the following file from [Zenodo](https://zenodo.org/record/1324204#.W2BmRn7fNE4).
>
>    ```
>    CSV dataset with several species: 
>    https://zenodo.org/record/1324204/files/Dataset%20multispecies%20Regional%20GAM.csv?download=1
>    ```
>   
> ### {% icon tip %} Tip: Importing data via links
>    > 1. Click on the Upload button: 
>    > 
>    > ![upload button](Images/upload.PNG "Galaxy upload access button")
>    > 
>    > 2. To import the dataset:
>    > * Select **Paste/Fetch data**
>    > * Paste the link into the text field
>    > * Galaxy will normally automatically find the right format, here CSV, so you don't have to specify it
>    > * Press **Start** and **Close** the window

>    {: .tip}

> ### {% icon comment %} Comment
>
> ⚠️ <a name="inputdatawarning"></a>Please note that the file must contain a header corresponding to: ```"SITES","SPECIES","YEAR","MONTH","DAY","COUNT"```, and that all the non numeric content must be between double quotes as "x" and that separators have to be ",". 

> {: .comment}

>    > ## Re-sampling. 
When the dataset contains many details, it lengthens the file processing time therefore it can be very useful to learn how to hide the informations you don't need. For example, the list of SITE of the dataset you are using is really long and the SITES are classified into sub-sites. Here, we will assume that your file doesn't really need be as precise and this is the reason why you have to specify you don't want the sub-sites. To create a new "down-sampled" file, you can follow these steps:   

> ### {% icon hands_on %} Hands-on: hiding some informations
>    > 1. Use the `CSV to tabular` tool to first create a tabular file from your csv one (with only one species). This is a mandatory step as further tools are only working on tabular files!
>    > 2. Search for the tool `Column Regex Find And Replace` on the tabular file with the following  parameters.
>    >  * Select the input file & the column with the `SITE` header.
>    >  * "Trouver l'expression suivante": `(\.[0-9]+)` which specifies that you don't want the sub-sites (all suites of digits following a "." character) to be taken into account.
>    >  * "Remplacement":`leave it empty`.
>    > 3. Now, as regionalGAM tools use CSV files as input, you can regenerate a CSV file using the `tabular to CSV` tool on the output from **Column Regex Find And Replace**. Please, tag your new dataset with an explicit tag as "Count" and/or rename this dataset like "Count file".
>    
{: .hands_on}

# Step 2: <a name="selectionningonespecificspecies"></a> Selectionning one specific species and showing all the data corresponding to it

The second step of any Regional GAM data analysis is making sure to have one dataset of only one specific species that you will then be able to use. If you want to create a graph showing abundance evolution by years of several species, you will have to superimpose the graphs on one another. 

> ### {% icon hands_on %} Hands-on: How many species are taken into account in this dataset
>
> As the dataset is quite big and countains heterogeneous informations, you want to know wether the data are about one species or more.
> 1. Search for the tool `Count occurrences of each record`with the following parameters:
> * "from dataset": `output`from **tabular to CSV**
> * "Count occurrences of values in column(s)": Specify the `SPECIES` column, normally `column 1`
> * "Delimited by": `Tab`.
> * "How should the results be sorted?": `With the most common values first`.
> 2. Inspect the file by clicking on the `eye` icon to check how many species are taken into account.

>    > ### Creating a new file concerning only the data of one species
>    > 1. Copy the name of the species you are interested in from the CSV file (for example: "Aglais io").
>    > 2. Search for the tool`Filter data on any column using simple expressions` with the following   parameters.
>    > * With following condition: `c1=='habitat2'` replacing 'habitat2' with the name of the species (for example: `c1='"Aglais io"''`)  
>    > * Number of header lines to skip: `1`.
>    > * You can repeat this set of actions as much as necessary, changing only the name of the species taken into account.  By doing this, you will obtain separated dataset, each of them concerning a different species.
>    > 3. Search for the tool `tabular to CSV` with the following parameters 
>    > * Select the `output` from **Filter data on any column using simple expressions**
>    > * Separators: `","`.
>    > 4. Repeat `tabular to CSV` on all the different `outputs`from **Filter data on any column using simple expressions** that you have.
>    
>   > {: .comment}
>
>   > ### {% icon question %} Questions
>   >
>    > 1. How many species does your dataset take into account ?
>    > 2. What are their names ? 
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>The dataset contains informations on 2 different species. </li>
>    >    <li>Their names are "Pyronia tithonus" and	"Aglais io".</li>
>    >    </ol>
>    >    </details>
>    {: .question}

{: .hands_on}

❗❗ Now that you have done the steps specific to a multispecies dataset, you will be redirected to the reference tutorial on regionalGAM. In order to go further on regionalGAM analysis, start the tutorial from the second step. To do so, you can click here:
[Step 2 - Displaying the occurrence of a chosen species through the years of the reference tutorial on RegionalGAM](Reference_tutorial.md#displayingtheoccurrenceofthespecies).

# Step 3:  <a name="variousoccurrencesonasinglechartexplanations"></a>Creating a chart showing more than one occurrence 

It can sometimes be interesting to have the occurrences of various species reprensented on the same chart.
The example below shows you what a unispecies phenology chart looks like. 
![Phenology chart](Images/Phenology%20chart%20.png "This shows the occurrence of Pyronia tithonus"). 

> ### {% icon hands_on %} Hands-on: Add various occurrences on a single chart

First of all, you will have to combine the different `outputs` from  each **flight curve** (each of them containing the data on one species) in order to have a single dataset with all the species you wish to appear on your chart. If you want to do so, please follow the tip below: 

> ### {% icon tip %} Tip: Creating a file comporting all the data on various species
>    > 1. Search for the tool `Paste two files side by side` with the following parameters:
>    > * "Paste": `the output` from **merger des colonnes** (with the dataset concerning species 1)
>    > * "and": `the output` from **merger des colonnes** (with the dataset concerning species 2)
>    > * "Delimited by": tabulation 

> ### {% icon comment %} Comment
❗ Note that you can repeat `Coller deux jeux de données l'un à côté de l'autre` with as many files as you want (each of them concerning a different species). In order to do so you will have to do as explained below:
>    > * Search for the tool `Coller deux jeux de données l'un à côté de l'autre` with the following parameters:
>    >    * "Coller": the `output` from **Paste two files side by side** (with the dataset concerning species 1 and 2)
>    >    * "et": `the output` from **Merge Columns together** (with the dataset concerning species 3)
>    >    * "Délimité par": tabulation 
>    > * Repeat `Paste two files side by side` with `the output` from **Paste two files side by side** (with the data concerning species 1, 2 and 3) and with `the output` from **Merge Columns together** (with the dataset concerning species 4).

>   > {: .comment}

> ### {% icon tip %} Tip: Generating a multispecies chart
>    > Click on the last `output` from **Paste two files side by side**
>    > 1. Click on: {% icon tip %} Visualiser  
>    > 2. Select `Charts`
>    > 3. Give it a proper name
>    > 4. Select a visualization: "line chart (NVD 3) 
>    > 5. Select data 
>    > * "Provide a label": `The name of the species 1` 
>    > * "Pick a series color": Choose a color
>    > * "Data point labels": `Column corresponding to the name of the species 1` 
>    > * "Values for x-axis": `Column corresponding to the "week and year" of the species 1`
>    > * "Values for y-axis": `Column corresponding to nm of the species 1`
>    > 6. Insert data series:
>    > * "Provide a label": `The name of the species 2` 
>    > * "Pick a series color": Choose a different color
>    > * "Data point labels": `Column corresponding to the name of the species 2` 
>    > * "Values for x-axis": `Column corresponding to the "week and year" of the species 2`
>    > * "Values for y-axis": `Column corresponding to nm of the species 2`
>    > 7. You may repeat "Insert data series" as many times as needed depending on the number of different species you want to represent on your chart.
>    > 8. Click on {% icon tip %} `Customize`
>    > * "X-Axis label": `Week and Year`  
>    > * Y-Axis label: `nm values`
>    > * "Use multi-panels": click on `No`(or you will have separated charts, one for each species)
>    > 9. Click on {% icon tip %} `Visualize`
>    > 10. Click on {% icon tip %} `save this visualization`if you are willing to keep it

{: .hands_on}

❗❗Now that you have accomplished this part, please go back to tu reference tutorial in order to accomplish the final steps. To do so, please click here: [Abundance per year and per site from reference tutorial on regionalGAM](Reference_tutorial.md#Abundanceindex) 
