---
layout: tutorial_hands_on
title: FAIR and its Origins
abbreviations:
  FAIR: Findable, Accessible, Interoperable, Reusable
zenodo_link: ''
questions:
- What is FAIR and the FAIR Guiding Principles?
- Where does FAIR come from?
objectives:
- Identify the FAIR principles and their origin.
- Explain the difference between FAIR and open data.
- Contextualise the main principles of FAIR around the common characteristics of identifiers, access, metadata and registration.
time_estimation: "40M"
key_points:
- FAIR stands for Findable, Accessible, Interoperable and Reusable.
- Metadata, Identifiers, Registration and Access are 4 key components in the process of FAIRification.
- FAIR data is as open as possible, and as closed as necessary.
tags:
- fair
- dmp
- data stewardship
priority: 1
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
      - fair-metadata
requirements:
  - type: "internal"
    topic_name: fair
    tutorials:
      - fair-intro
---

FAIR as an acronym, its origins and its guiding principles are introduced.  Learners will be able to explain the difference between FAIR and open data and contextualise the main principles of FAIR around the common characteristics of identifiers, access, metadata and registration.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}



# FAIR

The concepts underlying the FAIR principles are intuitively grounded in good scientific practice, though formalising them as easily measurable is not straightforward.  Being completely ‘FAIR’ is somewhat of an ideal, where it is often more pragmatic to be ’FAIR enough’ for a particular purpose or use case.

FAIR’s primary goal is to maximise data reuse by researchers.  The FAIR principles enable reuse by helping researchers share and manage their data, though FAIR is not limited to data alone and can be applied to services, software, training and workflows.  The word FAIR is an acronym, derived from its major components: ‘F’indable, ‘A’ccessible, ‘I’nteroperable, and ‘R’eusable which form the foundation of the FAIR Guiding Principles. 

**Findable** means that data and its metadata can be found/discovered by humans and computers.  Part of this is making rich metadata and keywords available to search engines and data repositories, so that companion data can be discovered.

**Accessible** means that once discovered, data and metadata can be accessed/downloaded by humans and computers.  Typically this means the commitment of the resource to its long term hosting and availability, with a suitable licence, and in appropriate format.

**Interoperable** means that data and metadata are supplied in formats that can be easily used and interpreted by humans and computers. The file formats and terms (vocabularies) used can be integrated easily with other datasets and software.

**Reusable** means that metadata is rich, enabling appropriate reuse.  Commonly FAIR will encourage the use of community standards for data curation.

These 4 major components form headings for the 15 FAIR Guiding Principles shown in **Table 1.1**

| The FAIR Guiding Principles |                                                                                                                                                                                                                                                                                                                                       |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| To be Findable:             | F1. (meta)data are assigned a globally unique and persistent identifier<br>F2. data are described with rich metadata (defined by R1 below)<br>F3. metadata clearly and explicitly include the identifier of the data it describes <br>F4. (meta)data are registered or indexed in a searchable resource                                   |
| To be Accessible:           | A1. (meta)data are retrievable by their identifier using a standardized communications protocol <br>A1.1 the protocol is open, free, and universally implementable<br>A1.2 the protocol allows for an authentication and authorization procedure, where necessary <br>A2. metadata are accessible, even when the data are no longer available |
| To be Interoperable:        | I1. (meta)data use a formal, accessible, shared, and broadly applicable language for knowledge representation. <br>I2. (meta)data use vocabularies that follow FAIR principles<br>I3. (meta)data include qualified references to other (meta)data                                                                                         |
| To be Reusable:             | R1. meta(data) are richly described with a plurality of accurate and relevant attributes <br>R1.1. (meta)data are released with a clear and accessible data usage license<br>R1.2. (meta)data are associated with detailed provenance<br>R1.3. (meta)data meet domain-relevant community standards                                        |

Table 1.1: The FAIR guiding principles as described in Wilkinson, M., Dumontier, M., Aalbersberg, I. et al. The FAIR Guiding Principles for scientific data management and stewardship. Sci Data 3, 160018 (2016) {% cite Wilkinson2016 %}. 


> <question-title></question-title>
>
> Look at the wording of the FAIR principles in **Table 1.1**. Which terms are used more than once?  Which terms are you seeing for the first time?
>
> > <solution-title></solution-title>
> >
> > Shared terms are: **data** (F1, F2, F4, A1, A2, I1, I2, I3, R1, R1.1, R1.2, R1,3), **metadata** (F1, F2, F3, F4, A1, A2, I1, I2, I3,  R1, R1.1, R1.2, R1,3), **identifier** (F1, F3, A1), **protocol** (A1, A1.1, A1.2). <br> Not so obviously shared terms are: **vocabularies** (I1, I2), **access** (A1, A1.1, A1.2, R1.1) <br> Terms which you may be seeing for the first time are: metadata, persistent identifier, searchable resource, standardised communications protocol, authentication and authorisation procedure, knowledge representation, vocabularies, qualified references, data usage license, provenance, domain-relevant community standards.
> >
> {: .solution}
{: .question}


The aim of this course is to put these principles into context and familiarise learners with  terms used commonly around FAIR.  A glossary of terms appears at the bottom of this episode as a point of reference.

# Open data and FAIR 

![FAIR_is_not_the_same_as_Open_data](../../images/figure1-2_fair-vs-open.png "FAIR is not the same as Open data")
 

Since FAIR promotes data sharing, it is often misunderstood as Open data.
  
The [Open Data handbook](https://opendatahandbook.org/) defines Open data as **“data that can be freely used, reused and redistributed by anyone.”**  Commonly, FAIR data is open however FAIR compliance does not mandate access without restriction. Instead, for instances where FAIR data is subject to restricted access, the conditions of access need to be stated to be compliant with FAIR, for example access around sensitive data.  Often in FAIR data, people adhere to the philosophy of **“as open as possible, and as closed as necessary”** thereby maximising opportunities to reuse data.

> <question-title></question-title>
>
> Under which circumstances would restricted or closed  access to FAIR data be advisable?
>
> > <solution-title></solution-title>
> >
> > Where data is sensitive or subject to intellectual property.  Protecting sensitive data overrules mandating research data should be open access. <br> Note though that in most cases, people following FAIR principles will be looking to share their data openly. Also note that sensitive data can be released through anonymisation and in many cases subject to controlled access by the authority of the principal investigator or data access committee. 
> >
> {: .solution}
{: .question}


# What is meant by FAIRification and FAIRness of data? 
  
**FAIRification** is the process of making your data FAIR compliant by applying the 15 Guiding Principles shown in **Table 1.1**.  The extent to which you apply these principles defines the **FAIRness** of your data.  In other words, FAIRness refers to the extent by which your data is FAIR and implies some implicit means of measuring its compliance.

# FAIR’s origins

A [report from the European Commission Expert Group on FAIR data](https://zenodo.org/record/1285272#.ZD_9KhXMIqt) describes the origins of FAIR and its development in 2014-2015 by a FORCE11 Working Group.  The following exercise dips into this report and asks you to investigate some of FAIR’s history and foundation. 


> <question-title></question-title>
>
> Read page 11 of [the European Commission report](https://zenodo.org/record/1285272#.ZD_9QxXMIqu), under the heading “Origins and definitions of FAIR”. What benefit did the FORCE11 Working Group see to coining the word FAIR?
>
> > <solution-title></solution-title>
> >
> > The report states: “a FORCE11 Working Group coined the FAIR data definition, latching onto an arresting and rhetorically useful acronym. The wordplay with fairness, in the sense of equity and justice, has also been eloquent in communicating the idea that FAIR data serves the best interests of the research community, and the advancement of science as a public enterprise that benefits society.”
> >
> {: .solution}
{: .question}


# Course structure
  
During the remainder of this course we will put the FAIR Guiding Principles into context using the 4 FAIR characteristics of **metadata**, **data registration**, **access** and **identifiers**, and devote an episode to each of these.  Whilst we teach, we will map content to the appropriate FAIR principle and define all relevant terms identified.
 
# Glossary of FAIR terms
  
All terms in this glossary are mentioned in the FAIR Guiding Principles (Table 1.2) and referenced in the following episodes of this course.
  
  
  | The FAIR Guiding Principles |                                                                                                                                                                                                                                                                                                                                       |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| To be Findable:             | F1. (meta)data are assigned a **globally unique and persistent identifier**<br>F2. data are described with rich metadata (defined by R1 below)<br>F3. metadata clearly and explicitly include the **identifier** of the data it describes <br>F4. (meta)data are registered or **indexed in a searchable resource**                                   |
| To be Accessible:           | A1. (meta)data are retrievable by their **identifier** using a **standardized communications protocol** <br>A1.1 the protocol is open, free, and **universally implementable**<br>A1.2 the protocol allows for an **authentication and authorization procedure**, where necessary <br>A2. metadata are accessible, even when the data are no longer available |
| To be Interoperable:        | I1. (meta)data use a formal, accessible, shared, and broadly applicable **language for knowledge representation**. <br>I2. (meta)data use **vocabularies** that follow FAIR principles<br>I3. (meta)data include **qualified references** to other (meta)data                                                                                         |
| To be Reusable:             | R1. meta(data) are richly described with a plurality of accurate and relevant attributes <br> R1.1. (meta)data are released with a clear and accessible **data usage license**<br>R1.2. (meta)data are associated with detailed **provenance**<br>R1.3. (meta)data meet domain-relevant **community standards**            

Table 1.2: Terms used by the FAIR Principles that appear in this glossary are highlighted in **black**.

**Globally unique and persistent identifier**: a reference to a digital resource such as a dataset, a document, a database, etc, usually given as a URL, that takes the user to that resource.  The persistent identifier (PID) is unique, globally, in the sense that it is not used to identify any other digital resource.  The PID is persistent, in that it enables a resource to be located long term, preferably permanently.  For more see: [How to FAIR](https://howtofair.dk/how-to-fair/persistent-identifiers/) and [FAIR Cookbook](https://faircookbook.elixir-europe.org/content/recipes/findability/identifiers.html), episode 5 (identifiers).

**Indexed in a searchable resource**: repositories and catalogues that can be queried (search box) are examples of searchable resources.  Indexing refers to processes within the architecture of the data repository where (meta)data are organised so that they can be queried based on defined fields.  Indexing by an internet search engine (for example, Google) is another example of an indexed and searchable resource.  For more see: [GO FAIR](https://www.go-fair.org/fair-principles/f4-metadata-registered-indexed-searchable-resource/).
  
**Standardised communications protocol**: a method that connects two computers and ensures secure data transfer. Examples of this include the hypertext transfer protocol ([http](https://en.wikipedia.org/wiki/HTTP)(s)) and the file transfer protocol ([ftp](https://en.wikipedia.org/wiki/File_Transfer_Protocol)) that permit data to be requested and downloaded by selecting a link on a webpage, or launched from within a script, for example using an (application programming interface) API. For more see the [Australian Research Data Commons](https://ardc.edu.au/resource/standardised-communications-protocols/).
  
**Universally implementable**: in the context of a standardised communications protocol such as http(s), universally implementable means it can be used by a number of resources and software, for example Firefox, Chrome and Unix commands such as wget.  For more see the [Australian Research Data Commons](https://ardc.edu.au/resource/standardised-communications-protocols/).
  
**Authentication and authorisation procedure**: data repositories with authentication and authorisation procedures generally require a login (and password) to access (meta)data.  For more see the [Australian Research Data Commons](https://ardc.edu.au/resource/standardised-communications-protocols/).

**Language for knowledge representation**: in the context of (meta)data exchange between computers, (meta)data should be in formats that are universally recognised (interoperable standards).  For more see [GO FAIR](https://www.go-fair.org/fair-principles/i1-metadata-use-formal-accessible-shared-broadly-applicable-language-knowledge-representation/) and [FAIR Cookbook](https://faircookbook.elixir-europe.org/content/recipes/introduction/FAIR-and-knowledge-graphs.html).
  
**Vocabularies**: (or **controlled vocabulary**) is a dictionary of terms you can use when producing (meta)data.  Controlled vocabularies are often shared between databases and communities so by using them you can allow data from different sources to be merged (interoperable), based on a shared understanding of the concepts.   Ontologies are related to vocabularies, where terms in the vocabulary are organised by relations between them.  A commonly used ontology is the [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) where the term ‘Homo sapiens’ belongs to a hierarchy of parent terms such as ‘Primates’ and ‘Mammalia’.   The ontology defines the vocabularies and the parent/child relationships.   For more see Ten simple rules for making a vocabulary FAIR {% cite Cox2021 %}.
  
**Qualified references**: are terms used to describe relationships to pieces of (meta)data.  For more see [GO FAIR](https://www.go-fair.org/fair-principles/i3-metadata-include-qualified-references-metadata/).
  
**Data usage license**: describes the legal rights on how others use your data.  For more details on considerations in selecting a licence, or to find out more about types of licence that are available, see [RDMkit](https://rdmkit.elixir-europe.org/licensing). 
  
**Data provenance**: refers to metadata describing the origin of a piece of data, including information such as version, original location of the data, and usually an audit trail up to the current version.  For more see [RDMkit](https://rdmkit.elixir-europe.org/data_provenance#how-to-document-and-track-your-data-provenance).
  
**Community standards**: are standard guidelines used to structure and exchange data, usually supported by community-developed resources and/or software.  In the context of (meta)data, community standards relate often to the standardised ontologies used by a domain of research, and minimum information guidelines allowing data to be interoperable.  For more see [FAIRDOM](https://fair-dom.org/community_standards).
  
**Machine-readable**: though this term is not referenced in the FAIR principles, it is often discussed within the FAIR context. Machine-readable (meta)data is supplied in a structured format that can be read by a computer.  For more see [Open Data Handbook](https://opendatahandbook.org/glossary/en/terms/machine-readable/) and [RDMkit](https://rdmkit.elixir-europe.org/machine_actionability).
  
**(Meta)data**: is shorthand for ‘metadata and data’.

# Useful Resources
- The published FAIR Guiding Principles: Wilkinson, M., Dumontier, M., Aalbersberg, I. et al. The FAIR Guiding Principles for scientific data management and stewardship. Sci Data 3, 160018 (2016) {% cite Wilkinson2016 %}.
- Recipes for data FAIRification, written by domain experts giving real-world examples: [FAIR Cookbook](https://faircookbook.elixir-europe.org/content/home.html)
- Documentation and frameworks for data FAIRification.  Each of the 15 FAIR principles is put into context with real data examples: [GO FAIR](https://www.go-fair.org/fair-principles/)
- FAIR walkthrough using examples from across all academic disciplines: [How to FAIR](https://howtofair.dk/)