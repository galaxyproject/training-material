---
layout: tutorial_hands_on
title:  Making clinical datasets FAIR
abbreviations:
  FAIR: Findable, Accessible, Interoperable, Reusable
zenodo_link: ''
questions:
- Why make clinical data FAIR?
- How to make clinical data FAIR?
- test
objectives:
- Learn how to make clinical datasets FAIR
- Recognise why FAIR datasets are important
time_estimation: "10M"
key_points:
- FAIR data are data which meet principles of findability, accessibility, interoperability, and reusability (FAIR).
tags:
- fair
- open
- data stewardship
priority: 1
contributions:
  authorship:
    - SNG888
    - kkamieniecka
    
subtopic: fair-data

requirements:
  - type: "internal"
    topic_name: fair
    tutorials:
      - fair-intro

---


# Introduction

The life science community is generally very good at sharing omics data on platforms such as GEO and ArrayExpress.  However, the metadata and clinical data associated with the omics datasets are often incredibly sparse.  Tracking down the meta- and clinical data of omics sets can be time-consuming for both data owner and researcher.  This also means exploration of enriched omics is overlooked in a time when there are so many tools in data analyses, machine learning and AI.

In this tutorial, you will learn the benefits of making clinical datasets FAIR and how to make them FAIR.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Clinical datasets

Metadata that accompanies transcriptomic data
In GEO and ArrayExpress the elements that we are describing under clinical data is the Characteristics or Source Characteristics.

Hands-on: explore GEO and ArrayExpress for characteristics
 Click on the link below for ArrayExpress and open the section for Source Characteristics.  Also check out the Experimental Factors
https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12829 (few characteristics)
Click on the link below for GEO.  The first link describes the dataset.  The second link drills down to one of the samples 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65391
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1594476 (this dataset has a lot of characteristics)
Explore ArrayExpress and/or GEO for datasets within the domain(s) that you are familiar with
Take note the availability of characteristics data for the different datasets

Data collected during clinical trials
[to write up]

Data collected during observational studies or registries
[to write up]

Health records?
[to write up]

Hands-on: explore HDRUK Innovation Gateway https://www.healthdatagateway.org/

## The importance of making clinical datasets FAIR

In some recent work in rheumatoid arthritis, it was demonstrated that simple patient demographics such as sex, ethnicity and age are all drivers of expression variation in addition to disease activity.  Principal component driver plots highlighted critical associations between diverse clinical features and omics.  This shows how rich clinical information may be key to analysis in some, if not many diseases.  

![RA-MAP.](../../images/ra_map.png "")

For many diseases (or even pan-diseases), there is  need for well-powered cohorts in the hundreds to thousands in number.  FAIR principles have been widely promoted for omics datasets, but as shown, there are key challenges to meta-analysis.  This is illustrated by an example of collating 17 public rheumatoid arthritis studies.  The clinical data in most cases were sparse and incomplete.  For example, of the curated 17 studies, 3 did not include sex, 5 did not include age and only one study included ethnicity or race.  The limitations in the availability of clinical data, substantially diminishes the value of public data sets and violates FAIR principles.

Hands-on:  For your domain, look at datasets on ArrayExpress and GEO and see how many of them have sex, age, ethnicity/race, and some measure of disease activity.
Do the datasets have all these characteristics?

Usually, after finding potentially interesting datasets on GEO or ArrayExpress, the next step to trying to enrich these public datasets is to contact the PI of the study.  From experience, the PI is usually happy to share additional data, especially if there have been previous collaborations.  However, even given willingness to share, there is the contacting of the person who actually has the data and agreeing what can and will be shared.  This at best will take weeks, but most likely months, especially if there’s been a lapse in communications between parties.


# Methods to make clinical data FAIR
There are certain concerns when making clinical data freely available through public libraries such as GEO and ArrayExpress.  Data needs to be anonymised or at least pseudo-anoymised.  

Hands-on:  what data (screenshot of some data) would be problematic to share?
[need to create a dataset with items such as DOB, Names, Address/PostCode]

## Some practical techniques:

Converting dates to time periods 
Avoid including dates, and convert these into time periods.  For example, use date of birth and date of assessment to calculate Age (at assessment or baseline).  Similarly, can calculate Age at Diagnosis, and Age at Onset, if dates are given.

If there are more than one visits, and if these are not within roughly defined time periods, for example 4 weeks, 12 weeks, calculate the number of days between the visits.

Hands-on: from the dataset shown, which columns would you change to time periods from dates? For those columns, convert to time periods
[need to create a dataset that can be used for the next few hands-on]

When publishing omics datasets to libraries such as GEO and ArrayExpress, consider at least including Age, Sex, Race/Ethnicity, and if possible some disease activity measure.  Publishing as much as possible is great for researchers to access enriched datasets, and also reduces fielding queries and requests to the PI.

There are cases when PIs would rather not have the clinical datasets on free to access platforms.  In these situations, the clinical datasets can still be made FAIR.  By having clear documentation of the datasets, these documents can be published and attributed to the project or consortia.  For each individual clinical dataset, it is useful to have a data dictionary

## Data dictionaries
A data dictionary is the metadata of the dataset.  It should include details such as data type, value range or allowed values, any relationship to other data elements and their meaning or purpose.  The details of the data dictionary often originate from how the data was input either at source or into the database.

Here’s an example of a data dictionary used on a project


![data dictionary](../../images/data_dictionary.png "")

## Key things to note are:
The variable name is often the short form name of the data item.
The definition describes the data item more fully.
In this example data items are encoded, and so their format is numeric, with constraints.
Validation rule, usually refer to instructions given when data was originally input, but could be the calculation formula for the data item.
Codes and Labels here are the encoding definitions

Hands-on: create a data dictionary for the dataset
[need to create data dictionary for the hands-on dataset]

## Data catalogues
[to write up, include links to the rdmbites]

## Publishing on Zenodo
[to write up, show example of CLUSTER datasets on Zenodo, look for links that describe how to upload to Zenodo]

## Data access request forms
Data management committee
[these two can be combined, but are the mechanisms for sharing data based]
