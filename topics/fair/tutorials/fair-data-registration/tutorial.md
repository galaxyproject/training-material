---
layout: tutorial_hands_on
title: Data Registration
abbreviations:
  FAIR: Findable, Accessible, Interoperable, Reusable
zenodo_link: ''
questions:
- What is data registration?
- Why should you upload your data to a data repository?
- What types of data repositories are there?
- How to choose the right repository for your dataset?
objectives:
- Describe why indexed data repositories are important.
- Summarise resources enabling you to choose a searchable repository.
time_estimation: "40M"
key_points:
- A good way to FAIRify your (meta)data is through submission to a public repository if it indexes and exposes the appropriate level of metadata to serve your specific use case or serve your envisaged users.
- Use Repositories that  support controlled access to data if necessary.
- FAIRsharing is a useful resource to locate relevant public repositories.
tags:
- fair
- dmp
- data stewardship
priority: 3
contributions:
  authorship:
    - robertmand
    - nsjuty
    - smza
    - nsoranzo
    - saramorsy
    - kellsnow
    - khens
    - proccaserra
    - lcooper
    - sitjart
    - asmasonomics
    - bfranicevic
    - saskia-lawson-tovey
    - kkamieniecka
    - khaled196
    - poterlowicz-lab
  funding:
    - elixir-uk-dash
subtopic: pointers
follow_up_training:
  - type: "internal"
    topic_name: fair
    tutorials:
      - fair-access
requirements:
  - type: "internal"
    topic_name: fair
    tutorials:
      - fair-intro
      - fair-origin
      - fair-metadata
---

The concept of data registration is defined as well as ways in which data registration can be achieved.  Learners will be able to describe why indexed data repositories are important as well as resources enabling you to choose a searchable repository.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}




# Data registration and the FAIR Principles
Data registration relates to the following 3 FAIR Principles (Table 3.1).  
We will discuss and signpost these in this Episode.

| The FAIR Guiding Principles |                                                                                                                                                                                                                                                                                                                                       |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| To be Findable:             | F1. **(meta)data are assigned a globally unique and persistent identifier**<br>F2. data are described with rich metadata (defined by R1 below)<br>F3. metadata clearly and explicitly include the identifier of the data it describes<br>F4. **(meta)data are registered or indexed in a searchable resource**                                   |
| To be Accessible:           | A1. (meta)data are retrievable by their identifier using a standardized communications protocol <br>A1.1 the protocol is open, free, and universally implementable<br>A1.2 the protocol allows for an authentication and authorization procedure, where necessary<br> A2. metadata are accessible, even when the data are no longer available |
| To be Interoperable:        | I1. (meta)data use a formal, accessible, shared, and broadly applicable language for knowledge representation. <br>I2. (meta)data use vocabularies that follow FAIR principles<br>I3. (meta)data include qualified references to other (meta)data                                                                                         |
| To be Reusable:             | R1. meta(data) are richly described with a plurality of accurate and relevant attributes <br>R1.1. **(meta)data are released with a clear and accessible data usage license**<br>R1.2. (meta)data are associated with detailed provenance<br>R1.3. (meta)data meet domain-relevant community standards                                        

Table 3.1: The 15 FAIR Guiding Principles. Principles relating to data registration in **black**.

# What is data deposition and registration?

Data deposition and registration refer to the process of uploading data to a searchable resource and providing appropriate metadata to facilitate its discoverability.
For example, a data repository, where data and metadata can be uploaded, may enable it to be discovered, preserved and accessed.  Here we use the general term data repository to describe any online storage location that can host deposited (meta)data.

In the context of FAIR, data deposition relates to a number of the Guiding Principles. Firstly,  _“(meta)data are registered or indexed in a searchable resource, Indexed in a searchable resource: a resource where (meta)data are organised so that they can be queried based on defined fields, ”_ (FAIR Principle F4).  Searchable (indexed) metadata enables humans and computers to query and discover data of interest, though this depends on what is indexed.  Here, indexing refers to a process that occurs within the architecture of the data repository (local indexing) where metadata are organised so that they can be queried based on a defined field.  It is worth noting that community resources focused on a particular domain (for example, the human database in [Ensembl](https://www.ensembl.org/Homo_sapiens/Info/Index)) are better indexed for a particular community, rather than generic repositories (for example, [Zenodo](https://zenodo.org/)) which may not index the community specific components, and may focus on higher level metadata.  Indexing by an internet search engine is another example of this.  Google (and other search engines, such as Yahoo and Yandex) have an agreed vocabulary ([schema.org](https://schema.org/)), within web pages, that are ‘scraped’ and indexed.  While the focus of this vocabulary was originally intended for commercial products, community-specific efforts to facilitate discipline-specific indexing are underway (for example, [Bioschemas](https://faircookbook.elixir-europe.org/content/recipes/findability/seo/bioschemas-data-page.html)). 

# Why should I upload my data to a data repository?

Data repositories are generally preferred to file storage systems (such as Dropbox) or sharing data on an ad hoc basis since they often better support FAIR best practice.  Repositories will assign citable, _“globally unique and persistent identifiers”_ (FAIR Principle F1) to data, and in some cases enable a data submitter to apply a data usage licence through association with the resource (FAIR Principle R1.1).  

Although not exclusively, data repositories support the creation of metadata through curation interfaces providing drop-downs and text fields for metadata entry and validation.  Often in the case of a domain or data-specific data repository, such as BioStudies shown in the previous Episode, drop-downs for metadata curation will link community-endorsed vocabularies (FAIR Principle R1.3).

# Types of data repository

General public data repositories, such as [Zenodo](https://zenodo.org/), are multidisciplinary and permit registration and upload of open and closed access (meta)data.  Metadata curation is relatively high level and made searchable via indexing.  Relating to data in the Life Sciences, Zenodo is often used to publish and provide citable URLs to supplementary data within articles, usually in instances where a domain repository does not exist.

Institutional repositories work similarly and provide an online archive for hosting, indexing and preserving research output specific to an institution.  Typically these house more than data, providing a repository often for documents and articles.  Institutions will have their own systems supported locally or buy into company solutions. 

Discipline-specific repositories cater for communities and datatypes, and typically provide web interfaces to annotate rich metadata at the point when data are submitted.  Examples of these belong to the suite of data repositories at the [European Bioinformatics Institute (EBI)](https://www.ebi.ac.uk/submission/) where rich metadata creation is supported by teams of curators.


> <question-title></question-title>
>
> An example of a discipline-specific repository is [ArrayExpress](https://www.ebi.ac.uk/biostudies/arrayexpress) database.  ArrayExpress stores data from high-through functional genomics assays, such as RNAseq, ChIPseq and expression microarrays.
> The data submission interface of ArrayExpress is called [Annotare](https://www.ebi.ac.uk/fg/annotare/login/).  Without creating a login, what help is given to a person looking to submit a dataset for the first time?
>
> > <solution-title></solution-title>
> >
> > Both a [submission guide](https://www.ebi.ac.uk/fg/annotare/help/index.html) and [YouTube video](https://www.youtube.com/watch?v=ANr2PTVy7JA) is provided.
> >
> {: .solution}
{: .question}


> <question-title></question-title>
>
> Finding more help on how to upload data to specific repositories
> The [FAIR Cookbook](https://faircookbook.elixir-europe.org/content/home.html) is an online open resource housing specific ‘how to’ guides or recipes.  Use the FAIR Cookbook to find two recipes for “depositing data to Zenodo” and “registering datasets with Wikidata”, respectively.
>
> > <solution-title></solution-title>
> >
> > Open the **Findability** pulldown on the left hand banner to find recipes for the following:
> > [Depositing to generic repositories - Zenodo use case](https://faircookbook.elixir-europe.org/content/recipes/findability/zenodo-deposition.html) and [Registering Datasets in Wikidata](https://faircookbook.elixir-europe.org/content/recipes/findability/registeringDatasets.html).
> >
> {: .solution}
{: .question}




> <question-title></question-title>
>
> Choosing the right data repository for your data
> [FAIRsharing](https://fairsharing.org/) helps researchers identify suitable data repositories, standards and policies relating to their data.  Use this resource to identify data repositories for proteomic data.
>
> > <solution-title></solution-title>
> >
> > Access the search bar for the [FAIRsharing database registry](https://fairsharing.org/search?fairsharingRegistry=Database).  Search for proteomics and select “repository” under “Record Type”.
> > ![A screenshot of FAIRsharing showing the results for a search for all proteomics repositories.](../../images/figure3-1-fairsharing.png "FAIRsharing allows you to search for specific record types that are relevant for your area.")
> >
> {: .solution}
{: .question}


# Useful Resources
- Registries and lists of public repositories: [FAIR Cookbook](https://fairsharing.org/search?page=1&recordType=repository) and [nature](https://www.nature.com/sdata/policies/repositories) journal
- Publishing your data: [RDMkit](https://rdmkit.elixir-europe.org/data_publication)
- Using Bioschemas to embed metadata into webpages: [FAIR Cookbook Bioschemas](https://faircookbook.elixir-europe.org/content/recipes/findability/seo/bioschemas-datacatalog.html)
