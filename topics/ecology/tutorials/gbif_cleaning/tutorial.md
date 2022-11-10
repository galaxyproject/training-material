---
layout: tutorial_hands_on

title: Cleaning GBIF data for the use in Ecology
zenodo_link: ''
questions:
- How can I get ecological data from GBIF?
- How do I check and clean the data from GBIF?
- Which ecoinformatics techniques are important to know for this type of data?

objectives:
- Get occurrence data on a species
- Visualize the data to understand them
- Clean GBIF dataset for further analyses
time_estimation: 0H30M
key_points:
- Take the time to look at your data first, manipulate it before analyzing it
tags:
  - gbif
  - data management
  - data cleaning
contributors:
- yvanlebras
- sbenateau

---


# Introduction


GBIF (Global Biodiversity Information Facility, www.gbif.org) is for sure THE most remarkable biodiversity data aggregator worldwide giving access to more than 1 billion records across all taxonomic groups. The data provided via these sources are highly valuable for research. However, some issues exist concerning data heterogeneity, as they are obtained from various collection methods and sources.

In this tutorial we will propose a way to clean occurrence records retrieved from GBIF.

This tutorial is based on the Ropensci {% cite zizka2018 %} tutorial.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. Download data from GBIF
> 2. Clean data
> 3. Convert text data into GIS format
> 4. Visualize spatialized data
> {:toc}
>
{: .agenda}

# Retrive data from GBIF

## Get data

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the files from GBIF: **Get species occurrences data** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Scientific name of the species"*: write the scientific name of something you are interested on, for example `Loligo vulgaris`
>    - *"Data source to get data from"*: `Global Biodiversity Information Facility : GBIF`
>    - *"Number of records to return"*: `999999` is a minimum value
>
>    > <comment-title></comment-title>
>    >
>    > The spocc Galaxy tool allows you to search species occurrences across a single or many data sources (GBIF, eBird, iNaturalist, EcoEngine, VertNet, BISON). Changing the number of records to return allows you to have all or limited numbers of occurrences. Specifying more than one data source will change the manner the output dataset is formatted.
>    >
>    {: .comment}
>
>
> 3. **Check the datatype** {% icon galaxy-pencil %}, it should be `tabular`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular" %}
>
> 4. **Add tags** {% icon galaxy-tags %} to the dataset
>    - make them propagating tags (tags starting with `#`)
>    - make a tag corresponding to the species (`#LoligoVulgaris` for example here)
>    - and another tag mentioning the data source (#GBIF for example here).
>
>    Tagging dataset like this is good practice in Galaxy, and will help you 1/ finding content of particular interest (using the filtering option on the history search form for example) and 2/ visualizing rapidly (notably thanks to the propagated tags) which dataset is associated to which content.
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}


## Where do the records come from?

Here we propose to investigate the content of the dataset looking notably at the "basisOfRecord" attribute to know more about heterogeneity related to the data collection origin.

> <hands-on-title>"basisOfRecord" filtering</hands-on-title>
>
> 1. **Count** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"from dataset"*: `output` (output of **Get species occurrences data** {% icon tool %})
>    - *"Count occurrences of values in column(s)"*: `c[17]`
>
>
>    > <comment-title></comment-title>
>    >
>    > This tool is one of the important "classical" Galaxy tool who allows you to better synthesize information content of your data. Here we apply this tool to the 17th column (corresponding to the basisOfRecord attribute) but don't hesitate to investigate others attributes!
>    {: .comment}
>
{: .hands_on}

> <question-title></question-title>
>
> 1. How many different types of data collection origin are there?
> 2. What is your assumption regarding this heterogeneity?
>
> > <solution-title></solution-title>
> >
> > 1. 5
> > 2. each basisOfRecord type is related to different collection method so different data quality
> >
> {: .solution}
>
{: .question}



## Filtering data based on the data origin

> <hands-on-title>Filter data on basisOfRecord GBIF attribute</hands-on-title>
>
> 1. **Filter** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `output` (output of **Get species occurrences data** {% icon tool %})
>    - *"With following condition"*: `c17=='HUMAN_OBSERVATION' or c17=='OBSERVATION' or c17=='PRESERVED_SPECIMEN'`
>    - *"Number of header lines to skip"*: `1`
>
>    > <comment-title></comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
>
>    > <question-title></question-title>
>    >
>    > 1. How many records are kept and what is the percentage of filtered data?
>    > 2. Why are we keeping only these 3 types of data collection origin?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > 1. 470 and 8.79% of records were drop out
>    > > 2. These data collection methods are the most relevant
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
> 3. Add to the output dataset a propagating tag corresponding to the filtering criteria adding `#basisOfRecord` string for example
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}


## Have a look at the number of counts per record

Here we propose to have a look at the number of counts by record to know if there is some possible records with errors.

> <hands-on-title>Summary statistics of count</hands-on-title>
>
> 1. **Summary Statistics** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Summary statistics on"*: `out_file1` (output of **Filter** {% icon tool %})
>    - *"Column or expression"*: `c72`
>
> 2. Add to the output dataset a propagating tag corresponding to the filtering criteria adding `#individualCount` string for example
>
{: .hands_on}

> <question-title></question-title>
>
> 1. What is the min and max of individual counts?
>
> > <solution-title></solution-title>
> >
> > 1. From 1 to 100
> >
> {: .solution}
>
{: .question}



## **Filtering** data on individual counts

> <hands-on-title>Filter data on individualCount GBIF attribute</hands-on-title>
>
> 1. **Filter** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `out_file1` (output of **Filter** {% icon tool %})
>    - *"With following condition"*: `c72>0 and c72<99`
>    - *"Number of header lines to skip"*: `1`
>
>
> > <question-title></question-title>
> >
> > 1. How many records are kept and what is the percentage of filtered data?
> > 2. How can you explain this result?
> > 3. Which propagated tag you can propose to add here?
> >
> > > <solution-title></solution-title>
> > >
> > > 1. 50 and 89.29% o records were drop out
> > > 2. An important percentage of data were drop out because of many records whithout any value for this individual count field
> > > 3. As for the previous "count" step you are dealing with the `individualCount` column, you can add a to the output dataset a  `#individualCount` tag for example
> > >
> > {: .solution}
> >
> {: .question}
{: .hands_on}


## Have a look at the age of records

> <hands-on-title>Here we propose to have a look at the age of records, through the `year` GBIF attribute to know if there is some ancient data to maybe not consider.</hands-on-title>
>
> 1. **Summary Statistics** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Summary statistics on"*: `out_file1` (output of **Filter** {% icon tool %})
>    - *"Column or expression"*: `c41`
>
> 2. Add to the output dataset a propagating tag corresponding to the filtering criteria adding `#ageOfRecord` string for example
>
{: .hands_on}

> <question-title></question-title>
>
> 1. What is the year of the older and younger records?
> 2. Why do you think of interest to treat differently ancient and recent records?
>
> > <solution-title></solution-title>
> >
> > 1. From 1903 to 2018
> > 2. We can assume ancient records are not made in the same way than recent one so keeping ancient records can enhance heterogeneity of our dataset.
> >
> {: .solution}
>
{: .question}


## Filtering data based on the age of records

> <hands-on-title>Filter data on ageOfRecord GBIF attribute</hands-on-title>
>
> 1. **Filter** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `out_file1` (output of **Get species occurrences data** {% icon tool %})
>    - *"With following condition"*: `c41>1945`
>    - *"Number of header lines to skip"*: `1`
>
>    > <comment-title></comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
>
>    > <question-title></question-title>
>    >
>    > 1. How many records are kept and what is the percentage of filtered data?
>    > 2. Why are we keeping only data from 1945?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > 1. 44 and 11.76% of records were drop out
>    > > 2. This arbitrary date allow to have only quite recent records, but you can specify another year.
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
> 3. Add to the output dataset a propagating tag corresponding to the filtering criteria adding `#ageOfRecord` string for example
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

## Taxonomic investigation

> <hands-on-title>Investigate the taxonomic coverage, at the family level</hands-on-title>
>
> 1. **Count** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"from dataset"*: `out_file1` (output of **Filter** {% icon tool %})
>    - *"Count occurrences of values in column(s)"*: `c[31]`
>
>
>    > <comment-title></comment-title>
>    >
>    > This column allows us to look at the different families associated to records. Normally, looking at a unique species, we will obtain only one family
>    {: .comment}
>
{: .hands_on}

## Filtering

> <hands-on-title>Filter data on family attribute</hands-on-title>
>
> 1. **Filter** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `out_file1` (output of **Filter** {% icon tool %})
>    - *"With following condition"*: `c31=='Loliginidae'`
>    - *"Number of header lines to skip"*: `1`
>
>    > <comment-title></comment-title>
>    >
>    > We here select only records with the family of interest, Loliginidae
>    {: .comment}
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Is the filtering here of interest ?
> 2. Why keeping this step can be of interest?
>
> > <solution-title></solution-title>
> >
> > 1. No, because 100% of records are kept
> > 2. Because this is an important step we have to take into account in such a GBIF data treatment, and if your goal is to create your own workflow you plan to use on others species, this can be of interest to keep this step
> >
> {: .solution}
>
{: .question}

## Sub-step with **OGR2ogr**

> <hands-on-title>Convert occurrence dataset to GIS one for visualization</hands-on-title>
>
> 1. **OGR2ogr** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Gdal supported input file"*: `out_file1` (output of **Filter** {% icon tool %})
>    - *"Conversion format"*: `GEOJSON`
>    - *"Specify advanced parameters"*: `Yes, see full parameter list.`
>        - In *"Add an input dataset open option"*:
>            - {% icon param-repeat %} *"Insert Add an input dataset open option"*
>                - *"Input dataset open option"*: `X_POSSIBLE_NAMES=longitude`
>            - {% icon param-repeat %} *"Insert Add an input dataset open option"*
>                - *"Input dataset open option"*: `Y_POSSIBLE_NAMES=latitude`
>
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Did you have access to standard output and error of the original R script?
> 2. What kind of information you can retrieve here in the standard output and/or error?
>
> > <solution-title></solution-title>
> >
> > 1. Yes, of course ;) A previsualization of stdout is visible when clicking on the history output dataset and full report accessible through the information button, then stdout or stderr (here you can see warnings on the stderr)
> > 2. The stderr is showing several warning related to automatic variable name mapping from GBIF to OGR plus information about application of a truncate process on a particularly long GeoJSON value
> >
> {: .solution}
>
{: .question}


## Visualize your data on a GIS oriented visualization

From your GeoJSON Galaxy history dataset, you can launch GIS visualization.

> <hands-on-title>Launch OpenLayers to visualize a map with your filtered records</hands-on-title>
>
> 1. Click on the *Visualize* tab on the upper menu and select `Create Visualization`
> 2. Click on the OpenLayers icon
> 3. Select the GeoJSON file from your history
> 4. Click on `Create Visualization`
> 5. Select Openlayers
>
> > <question-title></question-title>
> >
> > 1. You don't see Opebnlayers? Did you know why?
> >
> > > <solution-title></solution-title>
> > >
> > > 1.If you don't see Openlayers but others visualization types like Cytoscape, this means your datatype is JSON, not geojson. You have to change the datafile manually before visualizing it
> > >
> > {: .solution}
> >
> {: .question}
>
{: .hands_on}


# Conclusion

In this tutorial we learned how to get occurrence records from GBIF and several steps to filter these data to be ready to analyze it! So now, let's go for the show!
