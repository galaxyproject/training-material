---
layout: tutorial_hands_on
logo: "Galrollup"
title: "1. Introduction to the field guide"
summary: "Introduction to the people, resource inputs and project outputs of Galaxy. Work in progress. Please help make it better?"
questions:
  - "What are the main identifiable project components?"
  - "How do they all fit together?"
  - "Who pays for the 'free' analysis computational resources?"
  - "How are decisions really made?"
  - "Who's in charge?"
  - "Why is contribution from the community so important?"
  - "How can I join in?"
objectives:
  - "Understand Galaxy community development stories as examples of how things get done"
  - "Understand how different communities work together to make things happen"
  - "Understand opportunities for engagement and contributing your skills"
 
time_estimation:
subtopic: anatomy
key_points:
  - "Galaxy is a complicated, participatory self-governing collaborative community project"
  - "Participants engage with the project on their own individual terms"
  - "Multiple perspectives and opinions are needed to see the whole interacting structure"
  - "Vision: Build hard features like reproducibility into an analysis framework with pluggable tools"
  - "Running code as a free service provides effective stress testing for software defects and useability"
  - "Project depends on and collaborates as part of the open science ecosystem"
  - "Started in genomics but now adapted to suit many data rich disciplines"
  - "New contributors and investigators are highly valued and warmly welcomed"
 
contributors:
  - ggsc
  - fubar2
 
---

> <comment-title>Note to contributors</comment-title>
> - Work in progress!! First draft to try to get a structure to make sense.
> - Needs many contributors to make it useful. What would you like to have known, when you first tried getting things done in Galaxy?
> Please add what's missing and fix what's broken. Headings are mostly stubs waiting to be edited and extended 
> - Trying to describe the big picture will necessarily be big. Will probably need to break this already very wordy module into separate parts.
> - Add your story or stories to the stories tutorial too please!
> - Ross has strong opinions. 
> - Many of them are probably wrong but he doesn't know which ones yet. 
> - Please feel free to contribute your own, to make this more useful to future readers.
>
{: .comment}


> <comment-title>Note to readers</comment-title>
> - This module introduces the most important parts of the Galaxy project. They are grouped arbitrarily and there are many overlaps, because Galaxy grows organically through collaborations, rather than by design.
> - The project continues to expand rapidly, so this module will need updating regularly
> - The Hub provides much more detail about many of the same structures and their activities, but this material is designed to provide simplified views of the project, so the Hub becomes easier to navigate.
> - This is an attempt at a kind of field guide to the ecosystem generating those Hub activities, for the use of participants trying to navigate it.
> 
{: .comment}

> <agenda-title>Contents</agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Motivation, goals and definitions

## Why does Galaxy need a field guide?

Galaxy is a complicated global enterprise, with many separate components that have a substantial impact on open science analysis practice, by making every analysis transparent and reproducible, without requiring any special effort from the user. Project deliverables address a much wider range of open science research than most other open source projects, because Galaxy has proven to be easily adopted for use in many distinct scientific fields. New kinds of analyses are easily accomodated, by packaging relevant existing open source command line applications, into interoperable and sharable Galaxy workflow tools. 

Many project activities build on this generic, adaptable, analysis infrastructure source code core, adding value in the form of free analysis and training resources and services for open science researchers. Topic based communities, such as for Proteomics and Climate science, maintain "best practice" tools and workflows for researchers in those fields, improving productivity, analysis reliability and trustworthiness for results, for each community of researchers. Those communities also contribute specialised training resources through the Galaxy Training Network. Material includes short videos and slides, together with hands-on exercises that are integrated into the Galaxy user interface, giving step by step guidance, for new users learning to use the Galaxy platform in general, and for useful, complicated specialised analysis workflows using real data.

In short, there is a lot going on. This makes it seem more confusing, but provides a wide range of different ways to engage with the project, depending on individual interests. For example, it serves as

- A convenient and reliable way to share complex open science analyses, for researchers,
- Inter-dependent open source code repositories, for software developers,
- A global open science collaboration, for investigators.

Each of these perspectives is valid, but woefully incomplete. Seeing the big picture requires many perspectives and opinions, because the project has grown so large, that few individuals can devote enough time to remain engaged in, and fully informed, about every project activity.

This guide is intended to allow many perspectives to be combined, to provide a more comprehensive view of the project. For those engaged in the project, it may provide useful insights into how new active communities and new collaborations can be started, and how individuals have become engaged on their own terms.

## How to use the field guide.

The [Hub](https://galaxyproject.org/) provides much more detail and current activities. This training material refers to the Hub where possible, but aims to present an orderly overview of the main structures. 

This introduction is recommended reading, because it describes how the guide is structured to help make it easier to use. It also provides a restricted definition for the term "community" to simplify the guide. Related topic sections are linked at the end of this module, or from the main Topic page.

The field guide includes descriptions by those involved, of how activities were started. These stories can provide information and strategies for participants thinking about initiating new activities, based on existing successful initiatives.

> <tip-title>Classifying components in a global open science project is a hard problem</tip-title>
> Every open science project functions in, and depends on the context of the global open science ecosystem. That context is an essential part of a complete project description, but is too complicated to address here. It is assumed that the reader is familiar with it, because they have found this guide, so it is not addressed further in this material.
> Given the complexities, descriptive categories and divisions are necessarily arbitrary. The project has grown organically, and of its own accord, not neatly from a blueprint. It is driven by shared vision and values, and constantly adapting to rapidly changing externalities and project growth.
> A useful guide breaks all the complexities into smaller, more manageable chunks, since the material is potentially overwhelming. Large open projects are complex, and each has particular complexities, so there are many alternative ways to describe them. 
> - Many participants, structures and processes could fit into more than one category, with many users engaging with more than one community for example. It seems that there is no obvious, “best” way to describe an open science collaboration. It is a kind of virtual ecosystem with many interdependent and dynamic components with complex interactions and synergies. Descriptive terminology for the components and interactions needs to be invented, since it is not yet in wide use.
> 
{:  .tip}
 
## Structure for the field guide

The arbitrary division used for this guide gives 3 high level categories, related in the following trivial model:

### *People + Resource inputs = Project outputs*

This simplifies the challenge of seeing it all at once a little, and allows the guide to be split into corresponding sections. "Details" below provide more information about each component in overview. This division does not change the fact that in practice, the project depends on them all working together for success. Galaxy flourishes, because all these components govern themselves, in an efficient and productive global collaboration.
 
> <details-title>People working together make the project possible</details-title>
> 
> Open source is a very productive way of delivering software. Shared values and participatory self-governance help people get things done. Galaxy participants must govern the collaboration for themselves, because there are many independent institutions and investigators, so no single institution has control
> 
> Individuals engage with the project according to their interestests, such as:
>  - Analysing and sharing experimental data
>  - Working as a member of the core professional team
>  - Contributing skills, support, code and ideas while working in a related field.
>  - Building best practice tool kits and workflows for specific fields
>  - Leading new activities in the collaboration by providing community and project leadership.
>  - Contributing to project governance
> 
> - A detailed module on [people and their interactions is available](../people/tutorial.html)
{: .details}


> <details-title>Resources from grants sustain the communities, team, collaborations and services</details-title>
> 
> - Computational resource allocations used for free services are provided by collaborating institutions, adding substantial value to the source code in the form of computing power for free analysis and training services.
> - Core "corporate" project services like outreach, communication and administration are needed to keep the project on track. Professional, dedicated staff are needed for source code curation, user and contributor support, usegalaxy.* services, ToolShed maintenance, GTN material and training services, and many other related resources. These also depend on community contributions, and there are many project management and administrative tasks in coordinating such a large and complicated global enterprise, that require dedicated effort. Professional staff are needed to support project related core functions:
>   - outreach, communication, project management.
>   - user support, system administration and software engineers.
> 
> - Without these, the source code repository would still be a valuable resource. However, the open science value added by these services is likely to be a multiple of the total grant investment in terms of return. Substantial highly skilled community effort, adds greatly to total project impact and value, greatly exceeding the specifics of all the individual grants.
> 
> - A detailed guide to [Resources used in project activities is available.](../resources/tutorial.html)
{: .details}

> <details-title>Open source resources are essential parts of project code</details-title>
> 
> Galaxy source code depends on thousands of other open source projects
> Flavours for Galaxy servers depend on thousands of open source command line packages
> 
> Galaxy builds on and adds value to all those resources.
> 
> It makes command line packages interoperable and reproducible, through a uniform interface.
> Provides a GUI and sharing for users, without requiring effort from the package developer.
{:  .details}

> <details-title>Project outputs: Project outputs for open science and scientists</details-title>
> 
> Galaxy source code is the core project deliverable. It is widely deployed in public and many more private settings, but it is used to support many other important project outputs. Those other project activities build on the source code, to provide a range of open science benefits such as:
> 
> - The usegalaxy.* public deployments are large, free analysis services that support tens of thousands of researchers each day.
> - The Galaxy Training Network supports training to improve researcher productivity, integrated into Galaxy's user interface.
> - Specialised training, toolkits and workflows for specific fields are supported and distributed by communities of practice
> 
> These and other outputs are [described in more detail here](../outputs/tutorial.html) 
{:  .details}

     

### Communities: Definition and importance in this guide.

Active, self-governing communities are an important part of the project ecology, because their activities extend the project and produce added value for open science. When participants with a shared interest are sufficiently motivated to organise themselves to work on project activities related to that shared interest, a new community is created. For example, participants with research interests in some particular fields, such as climate science, have created their own open, active, collaborative communities of practice, to share ideas, tool wrappers and workflows.

#### Definition

For clarity, the term is used here in a restrictive way, referring to communities that organise publicised activities for interested participants. It is hard for participants to engage with a community that does not organise any activities, although it may be important in other ways. The *functional definition* used in this guide, is that *a project community arises when those interactions become sufficiently frequent and distinct from existing communities, to require their own pages on the Hub*.

#### Origins and sustainability

Communities form when participants come together to work on a project related initiative. Participants who share a common interest are more likely to interact when there are project activities related to that interest. That special interest might be ecology, proteomics, muon science, India or almost anything else that generates enough activity and participation. 

Communities organise and publicise activities for interested participants through the Hub. Leadership and participation are responsive to the real needs of the community because they are motivated and organised by participants. The project provides support and publicity, but the community must be able to sustain itself, with minimal project resources.

Some communities have been initiated by the project team, but many have come from participants organising themselves around an interest. The are all open and welcome participants. 

#### Some clarification about what "open" means

In some highly specialised groups, such as the IUC, formal decision making on technical issues is restricted to designated "members". They are open, like other groups, because regular communication channels, code and documentation are public and visible. Their work depends on specific skills. They welcome contribution and ideas from anyone with an interest. In practice, developers who contribute regularly, are typically recruited by consensus of the existing members.

#### Importance

Communities are a core element of the success of any large open project. They represent self-selected individuals, who choose to work on a project activity in which they have a particular interest. Communities are an important and pervasive idea in the main "People" part of the field guide.

Ironically, the largest "community" of all, comprising the tens of thousands of researchers who use Galaxy for open science analysis each day, does not organise its own activities, so fails the restrictive definition. That does not diminish their importance, such as demonstrating productivity for grant renewals, and serving as the major recruitment source of individuals choosing to engage in active project communities. The sharable, reproducible analysis results generated using Galaxy, are important downstream open science outputs, so they appear in the "Products" part of this guide.

#### Shared values help communities flourish

Communities are open for participation, such as:

   - Developers contributing to the code repository are an example of an important community.
   - Communities of practice maintain "best practice" tools, workflows and training for specific kinds of open science analysis.
   - Regional communities support local training and other distributed activities.

Newly organised communities related to the project are welcomed. Open source shared values ensure that community activities are safe, productive and enjoyable:

- inclusive, participatory, professional behaviours are encouraged
- The project explicitly strives to welcome and engage contributors and users
- Community members are encouraged to help make the project better.

Communities are a core resource, adding substantial value to project grant resources. Project success is the result of efficient and coherent self-governing collaboration, involving many specialised active communities of contributors.  


# A field guide to the Galaxy collaboration in 4 parts. 

## 2. [People and their interactions](../people/tutorial.html)

## 3. [Resources used in project activities](../resources/tutorial.html)

## 4. [Project outputs and deliverables](../outputs/tutorial.html)

## 5. [Community development success stories](../stories/tutorial.html)

