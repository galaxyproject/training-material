---
layout: tutorial_hands_on
title:  Making clinical datasets FAIR
abbreviations:
  FAIR: Findable, Accessible, Interoperable, Reusable
zenodo_link: ''
questions:
- Why make clinical data FAIR?
- How to make clinical data FAIR?
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
Hands-on: describing clinical datasets - time 3 minutes - silent reflection
 How would you describe clinical datasets? In general terms, and specifically

In general, clincial data are clinical reports (from trials and observational studies) and individual patient data.  In this section we will explore some of the different types of clinical data.

## Metadata that accompanies transcriptomic data
In GEO and ArrayExpress the elements that we are describing as clinical data is the Characteristics or Source Characteristics.

Hands-on: explore GEO and ArrayExpress for characteristics
 Click on the link below for ArrayExpress and open the section for Source Characteristics.  Also check out the Experimental Factors
https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12829 (few characteristics)
Click on the link below for GEO.  The first link describes the dataset.  The second link drills down to one of the samples 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65391
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1594476 (this dataset has a lot of characteristics)
Explore ArrayExpress and/or GEO for datasets within the domain(s) that you are familiar with
Take note the availability of characteristics data for the different datasets

## Interventional clinical trials
In these trials, usually the participants have a condition or disease, and the research question is related to the effectiveness of a treatment.  Some of the participants are given the treatment(s) under investigation, and the other participants are given either a placebo or non-treatment drug.  Usually these are under double-blind conditions, both participant and doctor do not know whether they are receiving/dispensing treatment or placebo.  In addition, there may also be a set of control participants that do not have the condition or disease.  Data is usually collected at baseline or screening, before any treatment or placebo is given, and then at least one time point in the future.  The breadth and type of data collected can vary across trials, and this can be limiting when trying to pool data from different trials.  The data collected could be minimal, basic metadata to accompany sample or imagining data.  But also could be quite extensive including quality of life questionnaires filled out by the patient.

An example of an interventional clinical trial is TRACTISS, where the aim of the study was to establish whether rituximab improves clinical outcomes for Primary Sjögren's Syndrome patients
https://bmcmusculoskeletdisord.biomedcentral.com/articles/10.1186/1471-2474-15-21

## Observational clinical trials
In these trials, data is collected at least once, there may not be a follow-up collection.  The research question is usually related directly about the disease or condition, for example understanding biomarkers of a disease.  As with interventional clinical trials, the breadth of data collected varies on the trial.  There is often a control group but not always.  

PRECISEADS is an example of an observational clinical trial where the study's aim is molecular reclassification to find clinically useful biomarkers for systemic autoimmune diseases
https://www.imi.europa.eu/projects-results/project-factsheets/precisesads

## Registries
Registries usually record (over a period of time), disease progression of patients with the same disease or conditiion.  Often referred to a hub/collection point once diagnosis confirmed.  Patients are continually added to the registries and therefore registries can provide a large pool of patients with relatively similar data collected.  As these are registries of patients with a particular condition, there are no control participants.

Here are links to two registry projects.  Both websites describe the projects, the aims and the outputs from the studies.
UK Sjögren's registry: https://www.sjogrensregistry.org/index.php
UK JIA Biologics Registers: https://www.sjogrensregistry.org/index.php

## Electronic health records (EHR)
These are generally the data collected when a person has contact with the health services.  In the UK, EHR data is either at the Trust level or at country level, e.g. NHS England, and it is usually necessary to specify primary care (GP, pharmacy, dental and optometry) and/or secondary care (hospital and specialists) data.  This data is usually received as annoymised data, to avoid the ability to indentify individuals from the data.  

In the UK, Clinical Practice Research Datalink (CPRD) collects anonymised patient data from a network of GP practices, and links this data to a range of other health related data and provides a longitudinal, representative UK population health dataset.  For more information about CPRD go to https://cprd.com

The Health Data Research Innovation Gateway (https://www.healthdatagateway.org/) can help find and discover UK health datasets, tools and other resources to further health research.

Hands-on: explore HDRUK Innovation Gateway for datasets
Click https://web.www.healthdatagateway.org/search?search=&datasetSort=latest&tab=Datasets
In the filters, expand Publisher and in Search Filter type CPRD and enter.  Select CPRD
Expand Phenotype and search for rheumatoid, and select Rheumatoid Arthritis
Take a look at the information about the two datasets

Clear the filters, and search for datasets related to your areas of interest


# The importance of making clinical datasets FAIR

## Associations between clinical features and omics
In some recent work in rheumatoid arthritis, it was demonstrated that simple patient demographics such as sex, ethnicity and age are all drivers of expression variation in addition to disease activity.  Principal component driver plots highlighted critical associations between diverse clinical features and omics.  This shows how rich clinical information may be key to analysis in some, if not, many diseases.  

![RA-MAP.](../../images/RAMAP_graphic.png "")

Hands-on:  For your domain, look at datasets on ArrayExpress and GEO and see how many of them have sex, age, ethnicity/race, and some measure of disease activity.
Do the datasets have all these characteristics?

## Powering up cohorts for analyses and machine learning
For many diseases (or even pan-diseases), there is need for well-powered cohorts in the hundreds to thousands in number, for analyses purposes.  FAIR principles have been widely promoted for omics datasets, but as shown, there are key challenges to meta-analysis.  This is illustrated by an example of collating 17 public rheumatoid arthritis studies.  The clinical data in most cases were sparse and incomplete.  For example, of the curated 17 studies, 3 did not include sex, 5 did not include age and only one study included ethnicity or race.  The limitations in the availability of clinical data, substantially diminishes the value of public data sets and violates FAIR principles.

Similarly, in order to make the most of machine learning techniques, datasets need to be reasonably large.  

## Time-saving
Usually, after finding potentially interesting datasets on GEO or ArrayExpress, the next step to trying to enrich these public datasets is to contact the PI of the study.  From experience, the PI is usually happy to share additional data, especially if there have been previous collaborations.  However, even given willingness to share, there is the contacting of the person who actually has the data and agreeing what can and will be shared.  This at best will take weeks, but most likely months, especially if there’s been a lapse in communications between parties.

# Practical techniques to make clinical data FAIR
There are certain concerns when making clinical data freely available through public libraries such as GEO and ArrayExpress.  Data needs to be anonymised or at least pseudo-anoymised.  

Hands-on:  In the data below, what data would be problematic to share?

![PatientData.](../../images/Patient_data.png "")

Solution: In the above dataset, remove columns
C - FIRST
D - LAST
H - ADDRESS
L - ZIP
as these are clear identifiers to the person

In terms of location, preferable to keep this as broad as possible.  In this case keep column J - STATE, but drop column I - CITY and column K - COUNTY is debatable to keep or not.

Column B - BIRTHDATE is another clear identifier but could be replaced by age.

It is preferable not to keep dates (Columns M - DATEOFASSESSMENT and N - DATEOFDISEASEONSET).  We will discuss date handling in the next section

## Converting dates to time periods 
Avoid including dates, and convert these into time periods.  For example, use date of birth and date of assessment to calculate Age (at assessment or baseline).  Similarly, can calculate Age at Diagnosis, and Age at Onset, if dates are given.

If there are more than one visits, and if these are not within roughly defined time periods, for example 4 weeks, 12 weeks, calculate the number of days between the visits.

Hands-on: from the dataset shown, which columns would you change to time periods from dates? For those columns, convert to time periods

![Patient_Dates.](../../images/Patient_dates.png "")

Solution: Replace BIRTHDATE with AGE (at DATEOFASSESSMENT) and DATEOFDISEASEONSET with AGEATONSET.  Remove BIRTHDATE, DATEOFASSESSMENT and DATEOFDISEASEONSET

![DatesRemoved.](../../images/Datesremoved.png "")

## Comprehensive metadata for omics datasets
When publishing omics datasets to libraries such as GEO and ArrayExpress, consider at least including Age, Sex, Race/Ethnicity, and if possible some disease activity measure.  Publishing as much as possible is great for researchers to access enriched datasets, and also reduces fielding queries and requests to the PI.

## Data dictionaries
There are cases when PIs would rather not have the clinical datasets on free to access platforms.  In these situations, the clinical datasets can still be made FAIR.  By having clear documentation of the datasets, these documents can be published and attributed to the project or consortia.  For each individual clinical dataset, it is useful to have a data dictionary

A data dictionary is the metadata of the dataset.  It should include details such as data type, value range or allowed values, any relationship to other data elements and their meaning or purpose.  The details of the data dictionary often originate from how the data was input either at source or into the database.

Here’s an example of a data dictionary used on a project

![datadictionary](../../images/ExampleDataDict1.png "")

Key things to note are:
The variable name is often the short form name of the data item.
The definition describes the data item more fully.
In this example data items are encoded, and so their format is numeric, with constraints.
Validation rule, usually refer to instructions given when data was originally input, but could be the calculation formula for the data item.
Codes and Labels here are the encoding definitions

Hands-on: Create a data dictionary for the dataset below

![dataDictHandsOn](../../images/DataDictionaryHandsOn.png "")

Solution: Here is one example of a suitable data dictionary for the dataset

![dataDictSoln](../../images/DataDictSolution.png "")

A key benefit of data dictionaries, is that the information captured is generally not confidential or patient-sensitive and therefore, there are few if any restrictions to sharing them.

## Data catalogues
Data catalogues are useful for teams, projects or consortia where there are many datasets being used and generated by the group.  Data catalogues are the metadata for the individual datasets, and they can provide context and provenance for each dataset within the group of datasets.

Here is an example of a data catalogue

![DataCat](../../images/DataCat.png "")

You may not want to capture all these metadata, and you may want to include other items such as phenotype, links to any publications related to the dataset, and links to the published dataset.

Additionally, if there are complex data sharing agreements, this information could be added to the catalogue (for internal use).

When projects and consortia are live, the corresponding data catalogue is a work in progress alongside the project and needs to be updated as required.  Although, data catalogues can be somewhat time-consuming to set up, this is more than made up with the time saved looking for details about datasets, especially as people switch out of the team.

As with data dictionaries, data catalogues rarely capture any information that need to be restricted access and usually can be shared in full.  The exception potentially maybe the contact person's email address.  

[To add: More information about data catalogues can be found in these RDMbites (add links)]

## Publishing on Zenodo
Data dictionaries and data catalogues include a wealth of information about the clinical datasets they reference, and in general, contain no sensitive or confidential data.  Zenodo (https://www.zenodo.org) is an open reposititory for research outcomes hosted by CERN.  It is free to use and each upload is citeable through a Digital Object Identifier (DOI).  Uploads are instantly available, with version control features.  Publishing on platforms such as Zenodo can help make your data findable.

Here you can see published data dictionaries for a consortium
![Zenodo](../../images/CLUSTER_Zenodo.png "")

To upload, create a log in using GitHub ID, ORCID, or email.  Click on Upload (at the top of the page, left of centre), and follow the instructions.

## Managing Data Access Requests
Regardless of whether or not you publish data dictionaries and data catalogues, it is worth considering the mechanisms when receiving data access requests.  

When a project is on-going, if there is a Data Management Committee (DMC), the responsibility of data access requests would be with the DMC.  It is really helpful to have a centralised point to provide data goverance including guidelines for data sharing.

It is really worth while clarifying the process of handling data access requests to minimise effort and delays for both the data holder and the data requester.  This process may need to be revised as the project ends.

Points to consider when handling data access requests:
- who will be the point of contact and the email address (e.g. create a general project email that is handled by the project administration team)
- expected length of time before requester receives response
- involve the institution's legal team 
  - terms and conditions of use of data and/or samples
  - publication guidelines
  - 
- standardise the process with a Data Access Request Form

An example of a data access process is shown below
![DataAccessProcess](../../images/DataAccessProcess.png "")

[To do: upload example Data access request form to Zenodo and add link (check Sharepoint for authorship guidelines)]
This is an example of a data and sample access request form and contracts associated with data and samples


In the UK, the HDRUK Innovation Gateway may also be a way to manage the data access request process. 

