---
layout: tutorial_hands_on
title: Access
abbreviations:
  FAIR: Findable, Accessible, Interoperable, Reusable
zenodo_link: ''
questions:
- What is data access in the context of FAIR
- What are the different types of data access?
- What is a data usage licence?
- How can you share sensitive data?
objectives:
- To illustrate data access in terms of the FAIR Principles using companion terms including communications protocol and authentication.
- To interpret the data usage licence associated with different data sets.
time_estimation: "50M"
key_points:
- Data access is supported by standardised communication protocols allowing for authentication where appropriate.
- Metadata must be available even if data are deleted.
- Data usage licences are applied to (meta)data detailing terms of use.
- Sensitive data can be subject to restricted access and/or de-identification.
tags:
- fair
- dmp
- data stewardship
priority: 4
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
      - fair-persistent-identifiers
requirements:
  - type: "internal"
    topic_name: fair
    tutorials:
      - fair-intro
      - fair-origin
      - fair-metadata
      - fair-data-registration
---
Data access including levels of access are described.  Learners will be able to illustrate data access in terms of the FAIR Principles using companion terms including communications protocol and authentication.  Learners will also be able to interpret the data usage licence associated with different data sets

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data access and the FAIR Principles
Data access relates to the following 5 FAIR Principles (Table 4.1).  We will discuss and signpost these in this episode.

| The FAIR Guiding Principles |                                                                                                                                                                                                                                                                                                                                       |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| To be Findable:             | F1. (meta)data are assigned a globally unique and persistent identifier<br>F2. data are described with rich metadata (defined by R1 below)<br>F3. metadata clearly and explicitly include the identifier of the data it describes <br>F4. (meta)data are registered or indexed in a searchable resource                                   |
| To be Accessible:           | **A1. (meta)data are retrievable by their identifier using a standardized communications protocol <br>A1.1 the protocol is open, free, and universally implementable<br>A1.2 the protocol allows for an authentication and authorization procedure, where necessary <br>A2. metadata are accessible, even when the data are no longer available** |
| To be Interoperable:        | I1. (meta)data use a formal, accessible, shared, and broadly applicable language for knowledge representation. <br>I2. (meta)data use vocabularies that follow FAIR principles<br>I3. (meta)data include qualified references to other (meta)data                                                                                         |
| To be Reusable:             | R1. meta(data) are richly described with a plurality of accurate and relevant attributes <br>**R1.1. (meta)data are released with a clear and accessible data usage license**<br>R1.2. (meta)data are associated with detailed provenance<br>R1.3. (meta)data meet domain-relevant community standards        


Table 4.1: The 15 FAIR Guiding Principles.  Principles relating to data access are highlighted in **black**.

# What is data access?
Making data accessible means it can be made available to use by both humans and computers, though FAIR data does not necessarily mandate that all data is openly accessible and we discuss this in a minute.

As outlined in the Accessibility Principles (meta)data must be retrievable by their identifier using a standardised communication protocol (FAIR Principle A1) which is open, free and universally implementable (FAIR Principle A1.1).

A standardised communication protocol is something like http(s) or FTP that allows data to be requested and downloaded, for example, by clicking on a link on a webpage.  Simply put, **a protocol** is a method that connects two computers and ensures secure data transfer.  Web browsers such as Firefox and Chrome can use the http(s) communication protocol since it is universally implementable, open and free.  

It is important to note that (meta)data access is not limited to humans clicking links on webpages. For a machine ‘user’, examples of accessing data include the use of an ‘application programming interface’ (API), and Unix command line, using  wget and [curl](https://curl.se/docs/manpage.html).

We mentioned above that FAIR data does not necessarily mandate open data.  Commonly, this relates to controlled access, which is discussed later in this Episode, but in the exercise below we investigate a dataset that has been deleted where it has been previously accessible.  FAIR states: _“metadata are accessible even when the data are no longer available”_ (FAIR Principle A2) which is one of the few principles that relate solely to metadata (and not data).  This means data can be deleted from its original online location at any time, but its original metadata must remain accessible.  PIDs (usually connected through URLs) to the original data must remain live, though the record may, for example, change to display metadata only. This is useful where originally referenced data changes over time, or becomes obsolete or deprecated; a record of the original metadata, and where appropriate reasons for its removal, or redirection to updated records, provide a provenance trail from the original data, which may have been referenced, for example, in publications, where it is not possible to update the PID.  


> <question-title></question-title>
>
> The [Protein Data Bank in Europe](https://www.ebi.ac.uk/pdbe/) (PDBe) is a searchable repository of biological macromolecular structures.  Please take a look at the following record that has been retired: [1ins](https://www.ebi.ac.uk/pdbe/entry/pdb/1ins).  You will see that the crystal structure is no longer available, but what metadata is available?
>
> > <solution-title></solution-title>
> >
> > Metadata includes the original citation, an explanation of why the record is no longer available and a redirection to the replacement  entry.
> >
> {: .solution}
{: .question}

# Types of data access

Where restricted access is required (for example, sensitive data or data subject to intellectual property) generally only the recommended parts are restricted, whereas associated metadata and non-sensitive data are openly accessible.  As part of the FAIR Principles, terms of access need to be stated, usually as part of a data licence (FAIR Principle R1.1).  Remember, FAIR data is **as open as possible, and as closed as necessary**.
Resources for working with sensitive data are given at the end of this Episode.

There are four types of data access as described in [RDMkit](https://rdmkit.elixir-europe.org/sharing#what-should-be-considered-for-data-sharing):

**Open access**: Anyone can access the data, and use it for any purpose.

**Registered access or authentication procedure**: researchers are required to register and authenticate to have the right to access the data (login and password)

**Controlled access or Data Access Committees (DACs)**: researchers will apply for access, and their application reviewed by a data access committee

**Access upon request (not recommended)**: a researcher provides their contact details for access. Contact details should be provided in the metadata which will be publicly available.

Any access requiring login and password makes use of _“[a protocol] allowing for an authentication and authorisation procedure, where necessary”_ (FAIR Principle A1.2).  Commonly for authentication, a researcher will be assigned a unique ID, or the system may support sign-in with an ORCID ID, which is the case of many data repositories including [Zenodo](https://zenodo.org/).  In some cases, there may be an option to use a Google account or institutional email to sign in, and many infrastructures also support ‘single sign in’.

# Data usage licence

A data usage licence describes the legal rights on how others use your data.  When you publish your data, you should describe clearly in what capacity your data can be used. 

There are many types of licences that can be used, including the  [MIT licence](https://opensource.org/license/mit/) (for software) or the commonly used [Creative Commons licences](https://creativecommons.org/licenses/) (a selectable collection of licences). These licences provide precise descriptions of the rights to data use, where the latter defines rights for sharing, adapting and commercialisation. Open access data usually carries the [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) or [CC0](https://creativecommons.org/share-your-work/public-domain/cc0/) licence permitting open sharing and adaptation, even for commercial purposes.  The licence is applied by adding the licence declaration to the data similar to this training page.  Take a look at the banner at the bottom of this page.  It states: Licensed under [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) 2018–2023 by [The Carpentries](https://carpentries.org/). 

# Making sensitive data accessible

Controlled access is often afforded to sensitive data or commonly any data that could potentially harm, see [RDMkit](https://rdmkit.elixir-europe.org/sensitive_data#is-your-data-sensitive) for a full definition.

Where access is granted, sensitive data is often de-identified, meaning that identifying (meta) data is removed or reassigned, leaving the analytics-appropriate component.

For example, a participant’s name can be removed from a questionnaire, and their home address can be substituted for the name of the town they live in.  This anonymises the data since the participant can no longer be located.  Note though that other information in the questionnaire could compromise this.  If other data reveal that the participant won the town’s 10K road race in 2023, we could potentially identify the individual using the name of the town and an online search.  If more information in the questionnaire states that the participant has a rare disease, we are broaching disclosure of sensitive, personal data.

Though people refer to anonymisation when de-identifying (meta)data, often they mean pseudonymisation.  Data anonymisation and pseudonymisation are slightly different.

**Data anonymisation** is the process of irreversibly de-identifying personal data such that an individual cannot be identified by anyone, including the study team and the individual themselves.  If data are anonymised, no one can link data back to the subject. 

**Pseudonymisation** is more commonly used.  It is a process where identifying-fields are replaced by artificial identifiers called pseudonyms or pseudonymised IDs.
Commonly, a person’s name or medical ID will be replaced with a unique participant ID within the study.
Pseudonymisation ensures no one can link data back to the individual, apart from nominated members of the study team who will be able to link pseudonyms to identifying information, such as medical records.


> <question-title></question-title>
>
> Use RDMkit's guidelines on [sensitive data](https://rdmkit.elixir-europe.org/sensitive_data) to familiarise yourself further on de-identification of data. What further training can you identify?
>
> > <solution-title></solution-title>
> >
> > At the bottom of the page, under “Training”, useful resources are given.  The [TeSS training portal](https://tess.elixir-europe.org/) permits users to search for courses, events, videos and other learning material for data in the life sciences.
> >
> {: .solution}
{: .question}


# Useful resources
- More about Data Access and Sharing: [RDMkit](https://rdmkit.elixir-europe.org/sharing), [FAIR Cookbook](https://faircookbook.elixir-europe.org/content/recipes/accessibility.html)
- Definition and examples of Standardised Communication protocol using authentication: [Australian Research Data Commons](https://ardc.edu.au/resource/standardised-communications-protocols/)
- More about Data licensing: [RDMkit](https://rdmkit.elixir-europe.org/licensing), [FAIR Cookbook](https://faircookbook.elixir-europe.org/content/recipes/reusability/ATI-licensing.html)
- The [Creative Commons licences](https://creativecommons.org/choose/)
- Working with sensitive data: [EBI training](https://www.ebi.ac.uk/training/events/working-sensitive-data/)
