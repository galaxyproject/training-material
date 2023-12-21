---
layout: tutorial_hands_on
title: FAIR data management solutions
abbreviations:
  FAIR: Findable, Accessible, Interoperable, Reusable
  DMP: Data Management Plan
  PID: Persistent Identifier

zenodo_link: ''
questions:
- Is there a reproducibility crisis?
- What can go wrong with data analysis?
objectives:
- Learn best practices in data management
- Learn how to introduce computational reproducibility in your research
time_estimation: "10M"
key_points:
- FAIR data management allows machines to automatically find and use the data accordingly.
tags:
- fair
- dmp
- data stewardship
priority: 2
contributions:
  authorship:
    - kkamieniecka
    - poterlowicz-lab
  editing:
    - hexylena
  funding:
      - ELIXIR-UK-DaSH
subtopic: fair-data

requirements:
  - type: "internal"
    topic_name: fair
    tutorials:
      - fair-intro

---


The FAIR (Findable, Accessible, Interoperable, Reusable)  data stewardship created the foundation for sharing and publishing digital assets, especially data. This apply to machine accessibility and emphasize that all digital assets should share data in a way that will enable maximum use and reuse.

This tutorial is a short introduction to the FAIR data management framework. You can find out more at the [FAIR Pointers](https://elixir-uk-dash.github.io/FAIR-Pointers/ep1/index.html) online course.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data management planning
In recent years we have notice a data explosion. Number of sequence records in each release of [GenBank](https://www.ncbi.nlm.nih.gov/genbank/statistics/), from 1982 to the present, doubled in size approximately every 18 months. Great amounts of data available followed by expanding range of tools and computational solutions result in reproducibility crisis. Having **data management plan (DMP)** in place is essential to achieve FAIR data management. DMPs are often described as living documents and should be updated according to changing circumstances.

There are several ways to set up FAIR (Findable, Accessible, Interoperable, Reusable) data management plans (DMPs) :
  - Findable (F): Data description and collection or reuse of existing data
  - Accessible (A): Standardised authentication or authorisation (e.g. HTTP, HTTPS)
  - Interoperable (I): Documentation and data quality
  - Reusable (R): Storage and backup supported by legal and ethical requirements

## Data description and collection or reuse of existing data
Reusing legacy datasets from institutional repositories or the digital libraries data collections can be FAIRified retrospectively. Support for data collection and development, throughout the life cycle can be provided and followed by change management and capacity improvement.

Multi-part FAIR research need a way of wrapping up, describing and sharing to promote the reuse of data. **Data sharing agreements** define the purpose of data sharing. Reference roles and responsibilities; specifies the purpose and legal requirements, e.g. for data security.

An institutional aim should be to create an integrated view and context over fragmented resources using their **persistent identifiers (PIDs)** and **metadata**. To make datasets findable, these metadata need to be as widely available as possible.

Enhancing reproducibility, quality and transparency by ensuring information flow and showcasing secondary use is also a part of data management. Promoting hands-on data experience and events activities built an collaborative environment for reproducible science.


## Documentation and data quality
Having access to local knowledge and encouraging best practises at the departmental level is a smart way to offer direction on a variety of standards and methods. In order to implement FAIR data practises within an institution, resources and infrastructure are needed. To increase the possibility of data reuse, several FAIR requirements can be satisfied using freely available guidelines e.g. [RDMkit](https://rdmkit.elixir-europe.org/), [FAIR Cookbook](https://faircookbook.elixir-europe.org/content/home.html), [ELIXIR-UK DaSH Fellowship](https://sites.google.com/view/navigation-portal-fellowship/home?authuser=0) initiative and repositories e.g. [Zenodo](https://zenodo.org/), [Harvard Dataverse](https://dataverse.harvard.edu/) and [figshare](https://figshare.com/).

## Storage and backup
Systems for storage, backup and collaboration depend upon technical infrastructure.The **'3-2-1 rule'**, a recommendation for saving three copies of the research data—two locally and one off-site—is a standard backup strategy for research data. Data, metadata, and other research artefacts, such as ontologies, software, documentation, and papers, must all be kept in locations where they are adequately safeguarded, backed up, and accessible to maximise their potential for reuse. Appropriate access management is essential, in addition to backup and restoration services that protect researchers against data loss, theft, malfunctioning computers or storage media, and accidental deletion or inadvertent alterations to the data.

The fundamental component of infrastructure required for the FAIR research data lifecycle are repository services. They allow access to the data, a persistent identifier, and the descriptive metadata that support interoperability. Repositories can include basic data storage, resource finding, managing access and use of private information, facilitating peer review of information related to publications or services requiring digital preservation, and more.

The [OpenAIRE](https://www.openaire.eu/opendatapilot-repository-guide) repository guide advises users to check the availability of a suitable repository in following order:

1. The most effective option (if available) is to maintain the data in accordance with acknowledged discipline-specific criteria using an established, dedicated (external) data archive or repository that caters specifically to the study domain.
2. Making use of institutional data repositories is the second-best option.
3. If none of those options is practical, a free data repository should be used.

Up-to-date lists of available registered data repositories can be found at [re3data](https://www.re3data.org/) and [FAIRsharing](https://fairsharing.org/search?fairsharingRegistry=Database).

## Legal and ethical requirements
Institutional support network (data stewards, ethics boards, IP, legal and financial offices) need to guide researchers in safeguarding data management responsibilities and resources.

# Conclusion
You will have the advantage of saving time and resources by planning how to FAIRify your data in the early phases of your research endeavour. To put this into action, a data management strategy, or DMP, must be written. A DMP is also where you outline your data collection, storage, processing, sharing, and disposal procedures. Planning the management and FAIRification of your data reduces the possibility of issues down the road, whether they be practical, legal, or technical.

Keep in mind that creating FAIR data is a complex process. Consider how you can make your data FAIR one step at a time at each stage of the creation, collection, documentation, storage, sharing, archiving, and preservation processes. The framework for the rest of your study planning is laid by incorporating your data documentation. Imagine you want to use a dataset that was generated by another researcher and how you would like it to be found. Hope this quick introduction to FAIR data management solutions will help you improve not only your experience with data but also influence others by using your guidance and FAIRified resources.
