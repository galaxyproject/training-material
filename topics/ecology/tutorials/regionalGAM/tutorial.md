---
layout: tutorial_hands_on
title: Regional GAM
zenodo_link: "https://zenodo.org/record/1324204#.W2BmRn7fNE4"
edam_ontology: "topic_0610"
questions:
    - "What are abundance and trends of a butterfly species?"
objectives:
    - "Obtain and filter/manipulate occurence data"
    - "Compute and visualize phenology of a species through the years"
    - "Compute temporal abundance trends"
time_estimation: "1h"
key_points:
    - ""
contributors:
    - Claraurf
    - emichn
    - yvanlebras
---

# Introduction
{:.no_toc}

This tutorial will show how to study species phenology through the computation of abundance index and trends. It will explain you how to use different [regionalGAM](https://github.com/RetoSchmucki/regionalGAM) tools on Galaxy-E allowing you to deal with datasets containing occurences informations for various species per site and per date through a couple of years.
After a certain numbers of steps, you will be able to extract single species data and study related phenology through the years. The goal of this exercise is to be able to create abundance trend over time and biodiversity indicators. Following these indicators allow to follow trends in terms of population dynamics. You could for example try to predict the occurences of one specific species in a certain type of environnement using the prediction model of climate evolution. Based on charts that you will generate, you could try to explain the evolution of a species with environmental data (temperatures variations, modifications of the environmental conditions).
You will basically learn how to create a file on the basis of which you can create a visual material that can be quite easily understood and therefore be efficient for a large audience.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Step 1: Pre-processing

The goal of the first step is to upload and prepare the file so that it will be usable for the *regional GAM* analysis.
First of all, you need to use a Galaxy instance with related regionalGAM tools. You can deploy your own local instance through Docker as a  Galaxy flavour or use our [Galaxy-E test instance](https://openstack-192-168-100-96.genouest.org/). 
After uploading input files, you might have to use some data handling tools to be able to use *regional GAM* tools.

>  ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial and give it a proper name as `Tuto training regionalGAM`.
> 2. Import the CSV dataset file with only one species  from [Zenodo](https://zenodo.org/record/1324204#.W2BmRn7fNE4).
>
>    ```
>    https://zenodo.org/record/1324204/files/regional%20GAM%20data.csv
>    ```
>   
>    {% include snippets/import_via_link.md %}
>
> 2. Check that the file contains a header corresponding to: ```"SITES","SPECIES","YEAR","MONTH","DAY","COUNT"```, and that all the non numeric content are between double quotes as "x" and that separators have to be ","
> {: .hands_on}

When the dataset contains many details, it lengthens the file processing time therefore it can be very useful to learn how to hide the informations you don't need. For example, the list of sites (look at the column with header `SITE)` of the dataset you are using is really long and the SITES are classified into sub-sites (like ESBMS.12, ESBMS.28,ESBMS.55,...). Here, we will assume that your file doesn't really need be as precise and this is the reason why you have to specify you don't want the sub-sites. To create a new "down-sampled" file (so delete the ---.12, ---.28 mentions), you can follow these steps:

1. This tool creates a tabular file from your csv one (with only one species). This is a mandatory step as further tools are only working on tabular files!

> ### {% icon hands_on %} Hands-on: hiding some informations
> 1. Run **CSV to tabular** {% icon tool %} with the following parameters:
>       - {% icon param-files %} *"CSV file"*: imported dataset
>       - {% icon param-files %} *"Separator"*: ","
>       - {% icon param-files %} *"Header in file"*: Yes
>   
> 2. **Column Regex Find And Replace** {% icon tool %} with the following parameters:
>       - *"Select cells from"*: Select the input file
>       - *"using column"*: the column with the `SITE` header.
>       - *"Find Regex"*: `(\.[0-9]+)`
>
>           It specifies that you don't want the sub-sites (all suites of digits following a "." character) to be taken into account.
>
>       - *"Replacement"*: leave it empty
> {: .hands_on}


> ### {% icon question %} Questions
>
> After having successfully deleted the sub-sites informations, can you look at the original dataset and this new one and say how many sites you had, and you have now? You will maybe need to use tools like `Count occurrences of each record` (don't forget that if you want to run the same tool with same parameters to several input files, you can directly specify the "multiple dataset" option on the tool form for the *"from dataset"* parameter).
>
> > ### {% icon solution %} Solution
> > The dataset contains 5 sites now against 1143 before "down-sampling"
> > {: .solution}
{: .question}

> ### {% icon details %} IF your original data is on RData format
> 
> > ### {% icon hands_on %} Hands-on: Data upload.
> > 1. Import the RData
> >
> >    For example, you can upload: 
> >
> >    ```
> >    https://zenodo.org/record/1324204/files/gatekeeper_CM%20.RData
> >    ```
> >
> > 2. **RData binary file reader** {% icon tool %} with the following parameters:
> >    -  *"Rdata binary file to explore"*: imported RData
> > 
> > 2. **RData parser** {% icon tool %} with the following parameters
> >    -  *"Rdata file to explore"*: imported RData
> >    -  *"File with .Rdata content details"*: output of **RData binary file reader** {% icon tool %}
> >    -  *"Select which attribute(s) you want to extract"*: select everything but "trend"
> >    -  *"Bind variables in a single tabular when its possible"*: `Yes`
> {: .hands_on}
> 
> Please note that if the tool **RData parser** {% icon tool %} don't succeed to create a single tabular file, it will generate separate files, each of them containing one column. The file with the "TREND" header can be let aside as we don't need it for what will follow.
> 
> > ### {% icon question %} Questions
> >
> > If Rdata parser fails to generate a single unified tabular file, can you propose a way to regenerate such a dataset ?
> >
> > > ### {% icon solution %} Solution
> > > You can do that using:
> > > 1. **Paste two files side by side tool** {% icon tool %} with the following parameters:
> >     - *"paste"*: output from **RData parser** {% icon tool %} headed with "SPECIES"
> >     - *"and"*: output from **RData parser** {% icon tool %} with headed with "SITE"
> > > 2. Repeat **Paste two files side by side** {% icon tool %} executions as many times as there are separated files in order to create a final dataset with all the columns:
> > >     1. Repeat **Paste two files side by side tool** {% icon tool %} to paste the file containing 2 columns with the one headed by "YEAR".
> > >     1. Repeat **Paste two files side by side tool** {% icon tool %} to paste the file containing 3 columns with the one headed by "MONTH". 
> > >     1. Repeat **Paste two files side by side tool** {% icon tool %} to paste the file containing 4 columns with the one headed by "DAY".
> > >     1. Repeat **Paste two files side by side tool** {% icon tool %} to paste the file containing 5 columns with the one headed by "COUNT". 
> > {: .solution}
> {: .question}
{: .details}



### Making sure the dataset concerns only one species 
 
The second step of any Regional GAM data analysis is making sure to have a dataset of only one specific species that you will then be able to use. If you want to create a graph showing abundance evolution by years of several species, you will have to superimpose the graphs on one another. 

> ### {% icon hands_on %} Hands-on: How many species are taken into account in this dataset
>
> As the dataset is quite big and may countain heterogeneous informations, you want to know wether the data are about one species or more. 
> 1. **Count occurrences of each record** {% icon tool %} with the following parameters:
>    - "from dataset": `output` from **CSV to Tabular**
> * "Count occurrences of values in column(s)": Specify the `SPECIES` column, normally `column 1`
> * "Delimited by": `Tab`.
> * "How should the results be sorted?": `With the most common values first`.
> 2. Inspect the file by clicking on the `eye` icon to check that the dataset is on one species only.
> 3. Now, as regionalGAM tools use CSV files as input, you can regenerate a CSV file using the `tabular to CSV` tool on the output from **Column Regex Find And Replace**. Please, tag your new dataset with an explicit tag as "Count" and/or rename this dataset like "Count file".
{: .hands_on}


> ### {% icon details %} Datasets containings informations about more than one species
>
> If your dataset contains informations about more than one species, you can apply the previous steps and then run an extra-step to select one specific species and show all the data corresponding to it.
>
> As the dataset is quite big and contains heterogeneous informations, you want to know wether the data are about one species or more. So the first step consists to count how many species are taken into account in this dataset.
>
> > ### {% icon hands_on %} Hands-on: How many species are taken into account
> > 1. **Count occurrences of each record** {% icon tool %} with the following parameters:
> >      - *"from dataset"*: `output` from **Column Regex Find And Replace**
> >      - *"Count occurrences of values in column(s)"*: the `SPECIES` column (normally `column 1`)
> >      - *"Delimited by"*: `Tab`
> >      - *"How should the results be sorted?": `With the most common values first`
> 2. Inspect the file by clicking on the {% icon galaxy-eye %} icon to check how many species are taken into account.
> {: .hands_on}
>
> To test these steps, you can use the following dataset: 
> 
>   ```
>   https://zenodo.org/record/1324204/files/Dataset%20multispecies%20Regional%20GAM.csv
>   ```
> 
> > ### {% icon question %} Questions
> >
> > 1. How many species does your initial dataset take into account ?
> > 2. What are their names ? 
> >
> > ### {% icon solution %} Solutions
> > 1. The dataset contains informations on 2 different species
> > 2. Their names are "Pyronia tithonus" and "Aglais io".
> > {: .solution}
> {: .question}
> 
> We now need to create a new file concerning only the data of one species
> 
> > ### {% icon hands_on %} Hands-on: Creating a new file concerning only the data of one species
> > 1. Copy the name of the species you are interested in (for example: "Aglais io").
> > 2. **Filter data on any column using simple expressions** {% icon tool %}
> >      - *"Filter"*: output of **Column Regex Find And Replace** {% icon tool %}
> >      - *"With following condition"*: `c1=='"Aglais io"'` or (another species name)  
> >      - *"Number of header lines to skip"*: `1`
> >
> >    You can repeat this set of actions as much as necessary, changing only the name of the species taken into account. By doing this, you will obtain separated dataset, each of them concerning a different species.
> >
> > 3. **tabular to CSV** {% icon tool %}
> >      - *"tabular file"*: output of **Filter data on any column using simple expressions** {% icon tool %}
> >      - *"output csv Separator"*: `,`
> >      - *"Header in file"*: `Yes`
> >
> >    Repeat this step on all the different `outputs` from **Filter data on any column using simple expressions** {% icon tool %} that you have, one by species. Please, tag your new dataset with an explicit tags as "Count" and "Aglais io" and/or rename this dataset like "Aglais io count file".
> {: .hands_on}
>
> If you want to create a graph showing abundance evolution by years of several species, you will have to superimpose the graphs on one another. 
{: .details}

# Step 2: Analyze phenology of a species through the years  

 
Now you have a file containing all the data on the species of interest. The main goal of this step is to treat phenology related informations and create a material that can be used to generate charts. What you could also do, for example, would be to compare the phenology through the years and sites.

This step will allow you to compute and display the phenology of a species. In the second part, you will learn that it is possible to show the phenology of various species on a single chart allowing to compare them and analyse them more easily.

> ## {% icon hands_on %} Hands-on: Phenology
> 1. **flight curve** {% icon tool %} with the following parameters: 
>    - *"Count file"*: output file you just generated with the **tabular to CSV** {% icon tool %}
>
> 2. Generate the chart using the visualization
>    1. Inspect and expand the output data from **flight curve** {% icon tool %}
>    2. Click on the {% icon galaxy-barchart %} (**Visualize**) icon
>    3. Select `Charts` 
>    3. Select a visualization type: "line chart (NVD 3)"
>    4. Give it a proper name, i.e. `Pyronia tithonus phenology raw simple vizu` 
>    5. On **Select data** area, specify:
>       - *"Provide a label"*: `Pyronia tithonus phenology from 2003 to 2012` for example
>       - *"Pick a series color"*: Choose a color for the line 
>       - *"Data point labels"*: `Column 1`
>       - *"Values for x-axis"*: `Column 2`
>       - *"Values for y-axis"*: `Column 6`
>    6. On **Change settings**, specify:
>       - *"X-Axis label"*: `Year`
>       - *"Y-Axis label"*: `nm values`
>    7. Click on **Visualize**
>    8. Click on **save this visualization**
>

> ## {% icon hands_on %} (Optional) Hands-on: Creating a new column of the dataset containing the week and the year 
> 1. **Count occurrences of each record** {% icon tool %} with the following parameters 
>    - *"from dataset"*: output from **Flight curve**
>    - *"Select"*: `Column 2` (the on headed with `YEAR`)
>    - *"Delimited by"*: `Tab`.
>    - *"How should the results be sorted?"*: `By the values being counted`.
> 2. Inspect and expand the output data from **Count occurrences of each record** {% icon tool %}
> 3. **Column Regex Find And Replace** {% icon tool %} with the following parameters:
>    - "File to process": output file from **flight curve**.
>    - "in column": `Column 2` (corresponding to the one headed with `YEAR`)
>    - Click on `Insert check`:
>       - "Find pattern": `(20[0-9][0-9])`
>       - "Replace with": `-\1` 
> 4. Inspect the file by clicking on the `eye` icon to check if all the years are now written with a "-" before the digits. 
> 5. **Merge Columns together** {% icon tool %} with the following parameters:
>    - "Select data": output from the last **Column Regex Find And Replace**.
>    - "Merge column": `Column 3`(corresponding to the one headed with `WEEK`)
>    - "with column": `Column 2`(corresponding to the one headed with `YEAR`)
> 7. Use **Remove beginning of a file** {% icon tool %} to remove first line
> 8. Use one more **Column Regex Find And Replace** {% icon tool %} to recreate the original content for the `YEAR` column with the following parameters:
>    - "File to process": output file from **Remove beginning of a file**.
>    - "in column": `Column 2` (corresponding to the one headed with `YEAR`)
>    - Click on `Insert check`:
>       - "Find pattern": `-(20[0-9][0-9])`
>       - "Replace with": `\1` 
> 
>    It is mandatory step to avoid header to be part of the visualization
>    If your dataset contains informations about more than one species, you can apply the previous steps and then run an extra-step to select one specific species and show all the data corresponding to it.
> 
> 2. Generate the chart using the visualization with the x-axis corresponding to your column `"week""year"`.
>    1. Inspect and expand the output data from **Remove beginning of a file** {% icon tool %}
>    2. Click on the {% icon galaxy-barchart %} (**Visualize**) icon
>    3. Select `Charts` 
>    3. Select a visualization type: "line chart (NVD 3)"
>    4. Give it a proper name, i.e. `Pyronia tithonus phenology simple vizu` 
>    5. On **Select data** area, specify:
>       - *"Provide a label"*: `Pyronia tithonus phenology from 2003 to 2012` for example
>       - *"Pick a series color"*: Choose a color for the line 
>       - *"Data point labels"*: `Column 6` (the nm column) or another one
>       - *"Values for x-axis"*: `Column 7` (the "week-year" column)
>       - *"Values for y-axis"*: `Column 6` (the nm column)
>    6. On **Change settings**, specify:
>       - *"X-Axis label"*: `Week-Year`
>       - *"Y-Axis label"*: `nm values`
>    7. Click on **Visualize**
>    8. Click on **save this visualization**
> 
{: .hands_on}

![Phenology chart](../../Images/Pyronia%20tithonus%20phenology%20explicit%20ID.png "This shows the occurrence of Pyronia tithonus")

Please note, that if you want to create a "stacked" visualization, overlapping each year, you can use several executions (one execution by year) of the `Filter data on any column using simple expressions` tool specifying the year you want on the `With following condition` parameter, `c2==2003` for 2003, then `c2==2004` for 2004, etc... then you can paste all resulting files side by side using one or several executions of the `Paste two files side by side` tool so you can specify on the "Select Data" tab of visualization, several Data series (one by year). WARNING The use of this `Paste two files side by side` tool must be done carefully as in case of differences in term of number of lines between datasets to paste, it will mix informations from columns. Here, as we are working on temporal series over years, and as some years have 365 days, others 366, and as in this case, phenology is not centered on winter, we can delete informations from the 366 days of some years without any problems so we will have the same number of lines between datasets to paste side by side. To do so, use the `Select first lines from a dataset` tool to select first 365 lines from filtered datasets. As it will be of interest to reuse this combination of tools in a next tutorial step, you can create a workflow that you can named something like `Phenology "stacked" visualization creation`.

With the `output` from **Remove beginning of a file** you can now generate a new 'stacked' chart which will have a x-axis corresponding to your column `"week"`.

> ### {% icon hands_on %} Generate a 'stacked' chart
>
> 1. Inspect and expand the output data from **Remove beginning of a file** {% icon tool %}
> 2. Click on the {% icon galaxy-barchart %} (**Visualize**) icon
> 3. Select `Charts` 
> 3. Select a visualization type: "line chart (NVD 3)"
> 4. Give it a proper name, i.e. `Pyronia tithonus phenology` 
> 5. On **Select data** area, specify:
>    - *"Provide a label"*: `2003` for example
>    - *"Pick a series color"*: Choose a color for the line 
>    - *"Data point labels"*: `Column 6` (the 2003 nm column) or another one
>    - *"Values for x-axis"*: `Column 7` (the 2003 "week" column)
>    - *"Values for y-axis"*: `Column 6` (the 2003 nm column)
> 6. Insert a new Data Series, choose a different color:
>    - *"Provide a label"*: `2004` for example
>    - *"Pick a series color"*: Choose a color for the line 
>    - *"Data point labels"*: `Column 12` (the 2004 nm column) or another one
>    - *"Values for x-axis"*: `Column 9` (the 2004 "week" column)
>    - *"Values for y-axis"*: `Column 12` (the 2004 nm column)
> 6. Insert a new Data Series, choose a different color:
>    - *"Provide a label"*: `2005` for example
>    - *"Pick a series color"*: Choose a color for the line 
>    - *"Data point labels"*: `Column 18` (the 2005 nm column) or another one
>    - *"Values for x-axis"*: `Column 15` (the 2005 "week" column)
>    - *"Values for y-axis"*: `Column 18` (the 2005 nm column)
> 6. Insert a new Data Series, choose a different color:
>    - *"Provide a label"*: `2006` for example
>    - *"Pick a series color"*: Choose a color for the line 
>    - *"Data point labels"*: `Column 24` (the 2006 nm column) or another one
>    - *"Values for x-axis"*: `Column 21` (the 2006 "week" column)
>    - *"Values for y-axis"*: `Column 24` (the 2006 nm column)
> 6. Insert a new Data Series, choose a different color:
>    - *"Provide a label"*: `2007` for example
>    - *"Pick a series color"*: Choose a color for the line 
>    - *"Data point labels"*: `Column 30` (the 2007 nm column) or another one
>    - *"Values for x-axis"*: `Column 27` (the 2007 "week" column)
>    - *"Values for y-axis"*: `Column 30` (the 2007 nm column)
> 6. Insert a new Data Series, choose a different color:
>    - *"Provide a label"*: `2008` for example
>    - *"Pick a series color"*: Choose a color for the line 
>    - *"Data point labels"*: `Column 36` (the 2008 nm column) or another one
>    - *"Values for x-axis"*: `Column 33` (the 2008 "week" column)
>    - *"Values for y-axis"*: `Column 36` (the 2008 nm column)
> 6. Insert a new Data Series, choose a different color:
>    - *"Provide a label"*: `2009` for example
>    - *"Pick a series color"*: Choose a color for the line 
>    - *"Data point labels"*: `Column 42` (the 2009 nm column) or another one
>    - *"Values for x-axis"*: `Column 39` (the 2009 "week" column)
>    - *"Values for y-axis"*: `Column 42` (the 2009 nm column)
> 6. Insert a new Data Series, choose a different color:
>    - *"Provide a label"*: `2010` for example
>    - *"Pick a series color"*: Choose a color for the line 
>    - *"Data point labels"*: `Column 48` (the 2010 nm column) or another one
>    - *"Values for x-axis"*: `Column 45` (the 2010 "week" column)
>    - *"Values for y-axis"*: `Column 48` (the 2010 nm column)
> 6. Insert a new Data Series, choose a different color:
>    - *"Provide a label"*: `2011` for example
>    - *"Pick a series color"*: Choose a color for the line 
>    - *"Data point labels"*: `Column 54` (the 2011 nm column) or another one
>    - *"Values for x-axis"*: `Column 51` (the 2011 "week" column)
>    - *"Values for y-axis"*: `Column 54` (the 2011 nm column)
> 6. Insert a new Data Series, choose a different color:
>    - *"Provide a label"*: `2012` for example
>    - *"Pick a series color"*: Choose a color for the line 
>    - *"Data point labels"*: `Column 60` (the 2012 nm column) or another one
>    - *"Values for x-axis"*: `Column 57` (the 2012 "week" column)
>    - *"Values for y-axis"*: `Column 60` (the 2012 nm column)
> 6. On **Change settings**, specify:
>    - *"X-Axis label"*: `Week`
>    - *"Y-Axis label"*: `nm values`
>    - *"Use multi-panels"*: `No`
> 7. Click on **Visualize**
> 8. Click on **save this visualization**
>    
{: .hands_on}

![Stacked Phenology chart](../../Images/Pyronia_tithonus_phenology_stacked_explicit_ID.png "This shows the occurrence of Pyronia tithonus")

> ### {% icon tip %} Tip: Creating a file comporting all the data on various species
>    > 1. Search for the tool `Paste two files side by side` with the following parameters:
>    > * "Paste": `the output` from **merger des colonnes** (with the dataset concerning species 1)
>    > * "and": `the output` from **merger des colonnes** (with the dataset concerning species 2)
>    > * "Delimited by": tabulation 
>    >WARNING The use of this `Paste two files side by side` tool must be done carefully as in case of differences in term of number of lines between datasets to paste, it will mix informations from columns. Here, datasets have the same number of lines.

> ### {% icon comment %} Comment
> Concerning a different species,. In order to do so you will have to do as explained below:
>    > * Search for the tool `Paste two files side by side` with the following parameters:
>    >    * "Coller": the `output` from **Paste two files side by side** (with the dataset concerning species 1 and 2)
>    >    * "et": `the output` from **Merge Columns together** (with the dataset concerning species 3)
>    >    * "Délimité par": tabulation 
>    > * Repeat `Paste two files side by side` with `the output` from **Paste two files side by side** (with the data concerning species 1, 2 and 3) and with `the output` from **Merge Columns together** (with the dataset concerning species 4).

> ### {% icon details %} Generating a multispecies chart
> 
> If your input dataset contains informations about more than one species, you can now generate char for the multispecies:
> 
> > ## {% icon hands_on %}
> > 1. Inspect and expand the output data from **flight curve** {% icon tool %}
> > 2. Click on the {% icon galaxy-barchart %} (**Visualize**) icon
> > 3. Select `Charts`
> > 3. Select a visualization: "line chart (NVD 3) 
> > 4. Give it a proper name like `Aglais io & Pyronia tithonus phenology`
> > 5. Select data 
> >     * "Provide a label": The name of the first species, for example `Aglais io`
> >     * "Pick a series color": Choose a color
> >     * "Data point labels": `Column corresponding to the name of the species 1` 
> >     * "Values for x-axis": `Column corresponding to the "week and year" of the species 1`
> >     * "Values for y-axis": `Column corresponding to nm of the species 1`
> > 6. Insert data series:
> >     * "Provide a label": he name of the second species, for example `Pyronia tithonus`
> >     * "Pick a series color": Choose a different color
> >     * "Data point labels": `Column corresponding to the name of the species 2` 
> >     * "Values for x-axis": `Column corresponding to the "week and year" of the species 2`
> >     * "Values for y-axis": `Column corresponding to nm of the species 2`
> > 7. You may repeat "Insert data series" as many times as needed depending on the number of different species you want to represent on your chart.
> > 8. Click on {% icon tip %} `Customize`
> >     * "X-Axis label": `Week and Year`  
> >     * Y-Axis label: `nm values`
> >     * "Use multi-panels": click on `No`(or you will have separated charts, one for each species)
> > 9. Click on {% icon tip %} `Visualize`
> > 10. Click on {% icon tip %} `save this visualization`if you are willing to keep it
> {: .hands_on}
{: .details}

## Compute Abundance Index across sites and years

This will allow you to create a file showing the abundance index per year of a chosen species in a certain site. Based on this file you will then learn how to represent this abundance on a chart. 
>
> 1. Look for the tool `Abundance index` with the following parameters:
> * "Count file": `output` from **tabular to CSV** (normally renamed "Counting file" and/or tagged "Count").  
> * "Flight curve output": `output` from **flight curve**.


>  Based on the  output from **abundance index**, we can create a chart showing the annual abundance trend of a certain species per site. 
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
>    > 1. What do you think about this visualization? Maybe not so good? Search a way to display the content of the file using charts in a more accurate manner... To do so, you can apply approaches used before on this tutotrial using tools like **Column Regex Find And Replace**, **Merge Columns together**, **Remove beginning of a file**, and **Sort data in ascending or descending order** on one hand to create a new column of more explicit identifiers (ie `Site-Year`) and/or  **Filter data on any column using simple expressions**, **Paste two files side by side**, **Remove beginning of a file** on another hand to create a stacked visualization.
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>You can use the **Column Regex Find And Replace** tool to first replace `(20[0-9][0-9])` on the column 3 (the "YEAR" one) by `-\1` then on the result of this tool execution, replace `"` by nothing on the column 1. Furthermore, you can merge column 1 and column 3 of the resulting dataset. Finally, after deleting the first line (the header) with **Remove beginning of a file**, you can sort the new dataset by column 1 (alphabetical/ascending) and column 3 (alphabetical/ascending). </li>
>    >    </ol>
>    >    </details>
>    {: .question}
>  If you choose to create a new column of more explicit identifiers (`Site-Year`), we can now display a better chart showing the annual abundance trend of a certain species per site. 
>    > 1. Select `Charts` from the **Sort data in ascending or descending order** execution
>    > 2. Select a visualization type: "Bar diagram (NVD 3)" 
>    > 3. Give it a proper name (`Pyronia tithonus abundance index` for example)
>    > 4. Select data 
>    > * "Data point labels": `Column 6` 
>    > * "Values for x-axis": `Column 6`
>    > * "Values for y-axis": `Column 4`
>    > 5. Customize 
>    > * "X-Axis label": `Site-Year`
>    > * "Y-Axis label": `regional_gam`
>    > 5. Visualize
>    > 6. Click on {% icon tip %} `save this visualization`if you are willing to keep it

![Abundance index chart](../../Images/Pyronia%20tithonus%20Abundance%20index%20explicit%20ID.png "This shows the occurrence of Pyronia tithonus")


>  If you choose to create a stacked visualization, we can now display a better chart showing the annual abundance trend of a certain species per site. 
>    > 1. Select `Charts` from the last execution of **Remove beginning of a file**
>    > 2. Select a visualization type: "Bar diagram (NVD 3)" 
>    > 3. Give it a proper name (`Pyronia tithonus abundance index` for example)
>    > 4. Select data
>    > * "Provide a label": This can be here the site, `UKBMS`
>    > * "Data point labels": `Column 5` (the UKBMS prop_pheno_sampled column)
>    > * "Values for x-axis": `Column 3` (the UKBMS YEAR column)
>    > * "Values for y-axis": `Column 4` (the UKBMS regional_gam column)
>    > 5. Insert a new Data Series, choose a different color and specify:
>    > * "Provide a label": This can be here the site, `NLBMS`
>    > * "Pick a series color": Choose a color for the line 
>    > * "Data point labels": `Column 10` (the NLBMS prop_pheno_sampled column)
>    > * "Values for x-axis": `Column 8` (the NLBMS YEAR column)
>    > * "Values for y-axis": `Column 9` (the NLBMS regional_gam column)
>    > 6. Insert a new Data Series, choose a different color and specify:
>    > * "Provide a label": This can be here the site, `ESBMS`
>    > * "Pick a series color": Choose a color for the line 
>    > * "Data point labels": `Column 15` (the ESBMS prop_pheno_sampled column)
>    > * "Values for x-axis": `Column 13` (the ESBMS YEAR column)
>    > * "Values for y-axis": `Column 14` (the ESBMS regional_gam column)
>    > 7. Insert a new Data Series, choose a different color and specify:
>    > * "Provide a label": This can be here the site, `FRBMS`
>    > * "Pick a series color": Choose a color for the line 
>    > * "Data point labels": `Column 20` (the FRBMS prop_pheno_sampled column)
>    > * "Values for x-axis": `Column 18` (the FRBMS YEAR column)
>    > * "Values for y-axis": `Column 19` (the FRBMS regional_gam column)
>    > 8. Insert a new Data Series, choose a different color and specify:
>    > * "Provide a label": This can be here the site, `DEBMS`
>    > * "Pick a series color": Choose a color for the line 
>    > * "Data point labels": `Column 25` (the 2004 prop_pheno_sampled column)
>    > * "Values for x-axis": `Column 23` (the 2004 YEAR column)
>    > * "Values for y-axis": `Column 24` (the 2004 regional_gam column)

![Abundance index chart](../../Images/Pyronia%20tithonus%20Abundance%20index%20stacked.png "This shows the occurrence of Pyronia tithonus")


>    > 5. Customize 
>    > * "X-Axis label": `Year`
>    > * "Y-Axis label": `regional_gam`
>    > 5. Visualize
>    > 6. Click on {% icon tip %} `save this visualization`if you are willing to keep it
{: .hands_on}

> ## Compute a collated index for each year and estimates the temporal trend

The expected temporal trend allows you to have an overview of the evolution of a species in a certain type of environment in the futur.

> ### {% icon hands_on %} Hands-on: Expected temporal trend
>    > 1. Look for the tool `Expected temporal trend` with the the following parameters: 
>    > * "Tabular file generated by the ab_index tool": output from **abundance index**.
>    
![Expected temporal trend](../../Images/Expected%20temporal%20trend.png "This shows the expected evolution of Abundance")

> {% icon warning %} Please note that sometimes the expected temporal trend can't be done on dataset. If you want this action to work, the occurences on your dataset must lie between the month of April and the end of the month of September.

Note also that you will obtain two files resulting of the action above. The first one will be the graph and the second one will contains the values of "x".

> ## Model temporal trend with a simple linear regression 

The point of doing a linear regression is to determinate if the year has an influence on the abundance of a species. 

>    > 1. Look for the tool `Model temporal trend with a simple linear regression` with the following parameters.
>    > * "File generated by the glmmpql/Expected temporal trend tool": `output 2` from **temporal trend**. 
>    > * "File generated by the ab_index tool": `output` from **abundance index**.

> ## Model temporal trend taking into account autocorrelation of residuals

Here we apply the same approach than at the previous step with addition of a correslation structure (ARMA(2,0)) to adjust the model.

>    > 1. Look for the tool `Linear regression ajusted for autocorrelation in the residuals` with the following parameters.
>    > * "File generated by the glmmpql/Expected temporal trend tool": `output 2` from **temporal trend**. 
>    > * "File generated by the ab_index tool": `output` from **abundance index**.


> ## Plot collated abundance index with trend line

Finally, a global trend (over years) can be computed and displayed using the **Plot abundance with trend line** tool selecting outputs from **abundance index** and **Linear regression ajusted for autocorrelation in the residuals** or **Model temporal trend with a simple linear regression**.

![Simple expected temporal trend](../../Images/trend_simple.png "This shows the expected evolution of Abundance")

![Adjusted expected temporal trend](../../Images/trend_adjusted.png "This shows the expected evolution of adjsuted Abundance")


{: .hands_on}
 
# Conclusions  

{:.no_toc}

In this tutorial, you have analyzed regional GAM data to extract useful informations in order to be able to show different tendencies of a chosen species. Therefore, you are now able to treat the dataset so that it shows only the data concerning one specific species of your choice. From there, you can show the occurrence of this species through the years first on a dataset and then on a visual chart. You have also learned how to represent on a single chart the occurences of various species. Afterwards, we have shown you how to create a dataset containing the informations on the abundance of a species per year and per site. Based on which you can henceforth visually represent the annual abundance trend on a chart. Thereafter, you have the possibility of showing the expected temporal trend, based on which you will be able to try predicting the future evolution a given species. The last part of this tutorial has shown you how to calculate the linear regression allowing you to determinate wether the year has an influence on the abundance of a species or not. 
