---
layout: tutorial_hands_on
logo: "Galrollup"
title: "5. Stories: Galaxy project community development in practice"
questions:
  - "Who controls new activities being added to the project?"
  - "Are there examples of new project initiatives being developed by participants?"
  - "What structures support community development initiatives"
  - "How are technical decisions made about the source code?"
  - "How are big decisions really made?"
  - "Who's really in charge of the project?"
  - "Why is contribution from the community so important?"
  - "How can people participate in the project as contributors?"
objectives:
  - "Understand some Galaxy community development stories as examples of how things are done"
  - "Understand how different project communities work together to make things happen"
  - "Understand the importance of community development and contributor engagement to project sustainability"
 
time_estimation: 10m
subtopic: anatomy
key_points:
  - "Galaxy is a global collaboration"
  - "It provides great value to scientists as a platform for transparent analysis"
  - "That value is created by coordinating the work of hundreds of individuals"
  - "The collaboration welcomes new ideas to expand and improve the project"
  - "Stories about how people found their way in starting new ideas help explain how the project functions"
  - "They may offer guidance and ideas for initiatives"
  - "New data rich fields, contributors and investigators are all invited to work with us"
 
contributors:
  - ggsc
  - fubar2
  - timothygriffin
---

> <comment-title>Note to contributors</comment-title>
> - WARNING! Work in Progress! 
> - This is far from complete - please add what's missing and fix what's broken
>
> - Idea here is to gather some success stories introducing new project initiatives.
> - Tim's proteomics story is here as a start
> - There are some other stubs to invite some text :)
> - Not meant to be complete - not useful to duplicate or replace the Hub.
> - Tell the story briefly of who/what/why/how it began and what opportunities it offers now.
> - Looking for stories illustrating how to get new things started
> - Please feel free to contribute your own, to make this more useful to future readers.
>
{: .comment}


> <agenda-title>Field Guide Part 4. Stories of engagement and community development</agenda-title>
>
> The stories here recount how participants created some existing project initiatives. All open projects encourage this, and this training material is designed to make it easier to navigate the process in Galaxy.
> 
> Community initiatives and new collaborations make the project grow organically, rather than by design, and the project now has active communities in a wide range of scientific fields. Project source code supports a highly adaptable, scalable open science analysis framework, readily "flavoured" with open source command line analysis packages as tools, to suit the complex analysis needs of researchers in any data rich field.
> Each of the initiatives and experiences described below have lives of their own in the Galaxy ecosystem. The kinds of leadership needed and the processes that were followed are useful examples, providing a view from inside some successful initiatives. Individual experiences show the varied ways it is possible for interested participants to make a difference.
> 
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Many paths to making new things happen in Galaxy

- Galaxy is already big and complicated, with many moving parts. Starting a new activity is hard enough, made even harder when it's not clear where to begin. Stories offer insight for planning new project activities.
- Stories were contributed by the people who helped organise the work. Individual leadership, responding to genuine need, supported by participant engagement, are the common themes. They illustrate a core shared value in open projects, that encourages everyone to engage in activities that interest them, where they can help improve the project's value to scientists.
- This part of the field guide has community initiatives, team initiatives and individual stories in separate sections, so readers can focus on their own interests. 

# Communities of Practice

## Proteomics communities: a case-study

Galay's unique features of open source tool flexibility and availability, provenance tracking, workflow development, scalability and convenient user interface, offered a very competitive low cost alternative to the predominantly Windows-based, standalone software options in the field.  As an example, the Galaxy for proteomics project at the University of Minnesota around 2010 was borne of the need for a scalable, enterprise workflow solution to support analyses of the large and increasingly numerous datasets being produced on campus.  The Minnesota Supercomputing Institute (MSI) supported these computational needs and sought out solutions.  Because of connections to the genomics community, MSI personnel were aware of Galaxy.  The forward-looking design of Galaxy offered potential for applications beyond genomics, and exploratory work was initiated to investigate suitability for proteomics needs.

These explorations resulted in revelations about new possibilities for using Galaxy beyond domain-centric analyses (proteomics), and further developing multi-omic applications such as proteogenomics which integrates genomic, transcriptomic and proteomic data for novel molecular characterization of systems.  The innovation offered by multi-omics, coupled with the numerous unique advantages of Galaxy (access and training in use of complex workflows by novice users in particular), proved a competitive advantage in seeking funding to support development work.  These qualities have carried on since early days of development to gained continued support of Galaxy-P, as plausible arguments can be made to reviewers that this is not just another software tool but rather an ecosystem with unique, beneficial features.  A number of other groups across the world (Europe, Australia) have also recognized these benefits, and via coordination through the Galaxy-P team, a community of developers and users has been built and sustained.


## Climate science communities

Anne?


## Small scale system administrators

Hans-Rudolf?

## Micro

## Ecology


# Team initiated specialised communities

These are open to community participation, but their activities are very technical in nature, so technical decisions are made by a small group of members. They are open to anyone with an interest in the maintenance of technical and related standards for reliable and heavily used tools and workflows. 

## IWC

## IUC

*Ross placeholder. Delete me. Not a good example here because it was a team initiative rather than from the community?*

- The IUC was created in 2012.
- An initiative by Greg Von Kuster.
- Invited most of the active toolshed contributors at the time.
- Unlike most other stories here, it was started by the Galaxy team.
- A group responsible for identifying reliable and trustworthy Toolshed tools
- The Toolshed has always been an open resource to which any author could contribute.
- Some automated checks at upload, but no human review of contributions.
- Contributions are of variable reliability.
- Users had to do their own due diligence for security and reliability.
- Demand for curated and reliable Galaxy tools was growing.
- Response was creating the IUC to help identify and mark reliable tools. 
- Like everything else, it has adapted to changing needs over time.
- Now provides: 
    - "Best practice" tool wrapper implementation standards.
    - Tool wrapper Github source code repository automation scripts.
    - Tool wrapper development utilities and support.
    - Wrappers for many popular and important new packages.
    - Training topic content contributions and support.
- Open to participation from community participants interested in helping out.
- IUC Gitter channel provides support and advice for contributors developing new tool wrappers.

> <details-title>Formation of the IUC</details-title>
> Greg Von Kuster <greg@bx.psu.edu>
> 
> Tue, 11 Sept 2012, 06:40
> 
> to Edward, Peter, jj@umn.edu, Brad, Ross, Ira, Björn, Dave, Anton, Daniel, Dave, Nate, Greg
> 
> Hello everyone,
> 
> Thanks to each of you for agreeing to be a charter member of the new Intergalactic Utilities Commission (IUC).  In agreeing to membership, you are now deemed to be an Intergalactic Utilities Commissioner!
> 
> The goals of the IUC include enabling the pervasive use of the main Galaxy tool shed by:
> 
> a. Ensuring the repositories available in the tool shed include contents that are functionally correct and optimized for installation into local Galaxies.  New features will soon be introduced to the tool shed enabling repositories to be flagged as "IUC Approved".  When appropriate, guidance should be provided to repository owners so that they can understand what changes are necessary for their repository to be approved.
> 
> b. Providing prioritization guidance to the Galaxy Tool Shed developers for the implementation of upcoming new features. 
> 
> A detailed understanding of the tool shed and it's current and upcoming features and behaviors is necessary in order for the IUC to be successful with these goals.  We would like to schedule a teleconference towards the end of September using Google Hangout to start guiding you through this process.  We'll work with you to schedule the conference so it is as convenient as possible for the various time zones, and we'll send out an agenda a few days prior to the scheduled discussion.
> 
> In preparation for this first teleconference, if you have not already done so, please make sure to read the tool shed wiki at http://wiki.g2.bx.psu.edu/Tool%20Shed.  I know that the wiki is a bit monolithic, and I am planning to break it up as soon as I get a chance, but it contains valuable details that are essential to understand as we start this process.
> 
> In addition, I am willing to answer some questions from you via email to help prepare for this if needed.  Let us know if you feel that yet another email list is appropriate for this (e.g., iuc@galaxyproject.org or something similar).
> 
> Thanks very much,
> 
> Greg Von Kuster <greg@bx.psu.edu>
>
{: .details}

# Individual contributor stories



## [Interesting user story](https://static.sched.com/hosted_files/gcc2022/13/NEWCOMER%20TO%20CONTRIBUTOR%20Awan%20Collins%20Savage.pdf)

## Insert your story here or email text to Ross please?

### Related field guide content for further reading

- [Introduction and definitions](../introduction/tutorial.html)
- [People and their interactions](../people/tutorial.html)
- [Resources used in project activities](../resources/tutorial.html)
- [Project outputs and deliverables](../outputs/tutorial.html)





