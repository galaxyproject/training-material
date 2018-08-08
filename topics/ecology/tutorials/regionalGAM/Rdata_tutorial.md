---
layout: tutorial_hands_on
topic_name: ecology
tutorial_name: Regional GAM
---

# Introduction
{:.no_toc}

⚠️ You might be willing to follow this tutorial if you want to learn how to deal with a dataset which is on RData format.

❗Please be aware that this tutorial is only a complement to [the refence_tutorial on regionalGAM](training-material/topics/ecology/tutorials/regionalGAM/Reference_tutorial.md) and that therefore there are some missing parts. 
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
3. Displaying the occurence of the chosen species through the years
> {:Displaying the occurence of the chosen species through the years}
4. Conclusion 
> {:conclusion}

# Step 1: Pre-processing

The goal of the first step is to upload and prepare the file so that it will be usable for the *regional GAM* analysis (See [this warning](#inputdatawarning) for more information about the input file.
First of all, you will have to upload the file on Galaxy-E and then you might have to use some data handling tools to be able to use *regional GAM* tools.

>  ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial and give it a proper name.
> 2. Import the following files from [Zenodo](https://zenodo.org/record/1324204#.W2BmRn7fNE4) or from a data
>    library named `regional GAM data tutorial`
>
>    ```
>    RData dataset with only one species:
>    https://zenodo.org/record/1324204/files/gatekeeper_CM%20.RData?download=1csv
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
>    >    * Press **Start** and **Close** the window
>    {: .tip}
>
> ### {% icon comment %} Comment
>
> ⚠️ <a name="inputdatawarning"></a>Please note that the file must contain a header corresponding to: ```"SITES","SPECIES","YEAR","MONTH","DAY","COUNT"```, and that all the non numeric content must be between double quotes as "x" and that separators have to be ",". 
> {: .comment}

> ❗As the dataset you have just uploaded on the RData format, you will have to process this binary file to obtain an appropriate CSV dataset. To do so, you can use the following tools:
>   > * Search for the tool `RData binary file reader`with the following parameters:
>   >      * "Rdata binary file to explore": "dataset on RData" 
>   > * Search for the tool `RData parser` with the following parameters:
>   >      * "Rdata file to explore": `"dataset on RData"`
>   >      * "File with .Rdata content details": file of **`RData binary file reader`**
>   >      * "Select which attribute(s) you want to extract": select everything but "trend"
>   >      * ⚠️ Please note that the tool `RData parser` creates separated files, each of them containing one column. The file with the "TREND" header can be let aside as we don't need it for what will follow.
>   > * Search for the tool `Coller deux jeux de données l'un à côté de l'autre` to create a file comporting all the data     required with the following parameters:
>   >      * "Coller":  outut from **RData parser** headed with "SPECIES"
>   >      * "et": output from **RData parser** with headed with "SITE"
>   >      * Repeat `Coller deux jeux de données l'un à côté de l'autre` as many times as there are separated files in order to create a final dataset with all the columns. First you must paste 2 columns together, then you must paste this last file with a third column and do this action again and again until your final file countains all the columns. 
>   >      * Repeat `Coller deux jeux de données l'un à côté de l'autre` pasting the file containing 2 columns with the one headed by "YEAR".
>   >      * Repeat `Coller deux jeux de données l'un à côté de l'autre` pasting the file containing 3 columns with the one headed by "MONTH". 
>   >      * Repeat `Coller deux jeux de données l'un à côté de l'autre` pasting the file containing 4 columns with the one headed by "DAY".
>   >      * Repeat `Coller deux jeux de données l'un à côté de l'autre` pasting the file containing 5 columns with the one headed by "COUNT". 
>   >  
>   > {: .comment}

>    > ### {% icon question %} Questions
>    >
>    > 1. In which specific case do you have to proceed to a particular set of actions on your dataset in order to be able to use it ?
>    > 2. Why do you need to use `Coller deux jeux de données l'un à côté de l'autre`? 
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>You only have to do these actions when you are using a dataset on the Rdata format. </li>
>    >    <li>Because you want to create a single dataset which countains all the data on a chosen species. You decided to upload a dataset on RData format and therefore you had to use the tools `RData binary file reader` and `RData parser`. This last tool treats the file and allows you to open it on Galaxy-E but it creates as many files as there are columns (when RData object is composed from a unique data table). This is the reason why you had to carry out on a set of actions ending by the creation of one complete file.</li>
>    >    </ol>
>    >    </details>
>    {: .question}

{: .hands_on}

❗❗ Now that you have done the steps specific to a RData dataset, you will be redirected to the reference tutorial on regionalGAM. In order to go further on regionalGAM analysis, start the tutorial from #Re-sampling. To do so, you can click here:
[Re-sampling of the reference tutorial on RegionalGAM](https://github.com/Claraurf/training-material/blob/ecology/topics/ecology/tutorials/regionalGAM/Reference_tutorial.md#resampling) 
