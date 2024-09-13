---
layout: tutorial_hands_on
title: Introduction to Data Management Plan (DMP) for Peatland Research and PeatDataHub
level: Introductory
abbreviations:
  DMP: Data Management Plan
  FAIR: Findable, Accessible, Interoperable, Reusable
zenodo_link: ''
questions:
- Why create a DMP for peatland research?
- How to apply FAIR data management solutions for peatland research?
objectives:
- Understand the importance of data management in peatland research.
- Learn how to create a data management plan (DMP) for peatland research
- Create a first draft of a data management plan for peatland research
time_estimation: "30M"
key_points:
- A Data Management Plan (DMP) is a live document that provides a framework for managing data material during and after the research project.
- FAIR data are data which meet principles of findability, accessibility, interoperability, and reusability (FAIR).

tags:
- fair
- open
- data management plan
priority: 1
contributions:
  authorship:
    - GabrielaLopez
    - kkamieniecka
    - robertmand
    - angiolini
    - poterlowicz-lab
  editing:
    - hexylena
subtopic: fair-data

recordings:
  - speakers:
    - GabrielaLopez
    captioners:
    - kkamieniecka
    date: '2024-09-12'
    length: 1H1M41S
    youtube_id: "MX-4AWyKP8E"

---


Peatland monitoring and research generates data that is used for management, conservation and policy related to these ecosystems. To ensure that your project data is collected systematically and can be analysed, shared, and re-used it is important to have a data management plan.

This tutorial is a short introduction to writing a Data Management Plan (DMP) with examples for peatland monitoring projects and initiatives.

You can find more information on DMPs in the [Research Data Management Toolkit for Life Sciences]( https://rdmkit.elixir-europe.org/)

> <agenda-title></agenda-title>
>
> In this tutorial we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# What is a Data Management Plan?

A Data Management Plan (DMP) is a live document that provides a framework for managing data material during and after the research project.

A data management plan will facilitate the process of collecting and managing your peatland data and will allow you to discuss with your team how to share data, considering the requirements of your organisation and funders.

# Key components of a DMP
DMP templates vary by institution and funding agencies. 

| DMP component  | What to include|
| ------------- | ------------- |
| Descriptions of the data  |    &bull; Project details <br> &bull; Types of data <br> collected and generated <br> &bull; Origin of data <br> &bull; Format and size of data |
|Data collection/generation | &bull; Data collection methods <br> &bull; Data quality standards| 
|Data management, documentation, and curation | &bull; Storage and accessibility <br> &bull; Metadata standards and documentation <br> &bull; Curation and preservation |  
|Data sharing and access | &bull; Where the data will be shared <br> &bull; When will the data be available |
|Data security (where relevant)|&bull; Risks to data security<br> &bull; GDPR considerations|
|Capabilities| &bull; Evaluation of institutional capabilities to preserve and manage data|
|Maintaining and implementing the Data Management Plan|&bull; Actions to maintain and deliver the DMP |
|Environmental considerations|&bull; Energy cost considerations related to data management|
|Responsibilities|&bull; Agreement on the responsibilities of named individuals for data management|
|Relevant institutional, departmental or study policies on data sharing and data security| &bull; Review of institutional policies for data management|

## Hands-on: Create a DMP for peatland research
* Download the DMP template for [Peatland Research and Monitoring](https://zenodo.org/records/11185774).
* Fill in the DMP template using the prompt questions.
* Review your DMP with your team and collaborators.

# FAIR data management solutions

After you have filled in the DMP template, review how you can make your DMP FAIR.

There are several ways to set up FAIR (Findable, Accessible, Interoperable, Reusable) data management plans (DMPs) 

* Findable (F): Data description and collection or reuse of existing data.  Using metadata and community recognised standards for curation that enable someone to find the data.
* Accessible (A): Standardised authentication or authorisation (e.g. HTTP, HTTPS) for open and restricted data access.
* Interoperable (I): Documentation and data quality, and community standards for data exchange (e.g. tab-delimited text).
* Reusable (R): Storage and backup supported by legal and ethical requirements, as well as data release policies.

More information is available in the [FAIR data management tutorial]({% link topics/fair/tutorials/data-management/tutorial.md %})

# Conclusion

Creating a DMP will allow you and your institution to make the most of your peatland data. This introduction to DMP aims to highlight the usefulness of data management for peatland research and monitoring.



##  Suggested reading
* [FAIR in a nutshell]({% link topics/fair/tutorials/fair-intro/tutorial.md %}) 
* [The published FAIR Guiding Principles: Wilkinson, M., Dumontier, M., Aalbersberg, I. et al. The FAIR Guiding Principles for scientific data management and stewardship. Sci Data 3, 160018 (2016)]({% cite Wilkinson2016 %})
* [Recipes for data FAIRification, written by domain experts giving real-world examples: FAIR Cookbook](https://faircookbook.elixir-europe.org/content/home.html)
* [Documentation and frameworks for data FAIRification. Each of the 15 FAIR principles is put into context with real data examples: GO FAIR](https://www.go-fair.org/fair-principles)
* [FAIR walkthrough using examples from across all academic disciplines: How to FAIR](https://howtofair.dk/)

## Resources
* [How to write a data-sharing agreement between collaborators](https://www.youtube.com/watch?v=iaZInoaHa04)
* [DMP online](https://dmponline.dcc.ac.uk)
* [DSWizard](https://bio.tools/tool/Data_Stewardship_Wizard)
* [WaterLands](https://waterlands.eu/)
* [PeatDataHub](https://peatdatahub.net/)
