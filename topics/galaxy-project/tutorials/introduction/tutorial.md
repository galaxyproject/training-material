---
layout: tutorial_hands_on
logo: "Galrollup"
title: "1. Introduction and definitions for the field guide"
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

# Introduction to a Galaxy field guide

## Motivation: Why maintain a field guide?

- Galaxy is a complicated global enterprise, so there are many different ways to describe all the project activities.
- It has a much wider range of deliverables than most open source projects. Project activities build on the code to support many kinds of analysis and training, across many diverse areas of open science.
- This provides many different opportunities for participants, according to their interests. For example, it is seen as:
    - A convenient and reliable way to share complex open science analyses by users
    - Inter-dependent open source code repositories by developers.
    - A global open science collaboration by investigators
- All these views are true, but incomplete. Seeing the big picture requires many perspectives and opinions.
- A separate Lesson gives examples of how activities were started. These stories can provide information and strategies for participants thinking about initiating new activities. 
- This section introduces the field guide for participants. The [Hub](https://galaxyproject.org/) can provide much more detail and current activities, particularly once you know what you are looking for. Where possible, links to the relevant Hub pages are provided in the detailed sections of the guide. The introduction shows how the guide is structured and gives a definition for the term "community" because it is used widely. Sections are linked at the end of this module, or from the main Topic page.

> <tip-title>Classifying components in a global open science project is a hard problem</tip-title>
> - Categories and divisions are necessarily arbitrary. 
> - The project has grown of its own accord, not neatly from a blueprint.
> - Driven by shared vision and values, adapting to rapid external developments
> - Many participants, structures and processes that fit into more than one category. 
>    - Most users identify with multiple different communities for example.
> - No obvious “best” way to describe an open science collaborative ecosystem.
> - No known terminology for the components or interactions.
> - Need to invent a descriptive approach.
{:  .tip}
 
## Structure of this guide

- The project exists in the context of a global open science ecosystem so that context is an essential part of a complete description. The reader is assumed to have some experience of these external issues.
- Successful projects are complex, and each has particular complexities, so there are many alternative ways to describe them. Here we adopt one arbitrary but useful starting point, with three categories in the form of a simple model:

### People + resource inputs = project outputs

- This simple equation summarise the complicated project, allowing a logical division of material for the field guide, to make it easier to break down some of the complexity. "Details" tabs below show more overview information about them.
- The field guide is divided into three corresponding sections. In practice, the project depends on them all, efficiently working together. Galaxy flourishes, because all these parts govern themselves in an efficient and productive collaboration.

 
> <details-title>People working together make the project possible</details-title>
> - Open source is a very productive way of delivering software. Shared values and participatory self-governance help people get things done
> - Galaxy participants must governed the collaboration themselves, because there are many independent institutions and investigators, so no single institution has control
> - Individuals engage with the project in their own ways, such as:
>     - Analysing and sharing experimental data
>     - Working as a member of the core professional team
>     - Contributing skills, support, code and ideas while working in a related field.
>     - Helping build best practice tool kits and workflows for specific fields
>     - Leading new activities in the collaboration by providing community and project leadership.
> - Contributing to project governance
{: .details}


> <details-title>Resources from grants sustain the communities, team, collaborations and services</details-title>
> - Core "corporate" project services like outreach, communication and administration are needed to keep the project on track
> - Source code management, user support, usegalaxy.*, GTN and other services.
> - Computational resource allocations for free services
>     - adding substantial value to the source
>     - Free analysis and training services
> - Professional staff
>     - outreach, communication, project management.
>     - user support, system administration and software engineers.
> - Without these, the source code repository would still be a valuable resource
> - These additional activites add value, in the form of highly valued resources for scientists.
{: .details}

> <details-title>Open source resources are essential parts of project code</details-title>
> - Galaxy source code depends on thousands of other open source projects
> - Flavours for Galaxy servers depend on thousands of open source command line packages
> - Galaxy builds on and adds value to all those resources.
> - It makes command line packages interoperable and reproducible, through a uniform interface.
> - Provides a GUI and sharing for users, without requiring effort from the package developer.
{:  .details}

> <details-title>Project outputs: Project outputs for open science and scientists</details-title>
> - Source code
> - Online analysis services
> - ToolShed automated tools for servers
> - Best practice standards and automation for writing tool wrappers and workflows - IUC
> - Communities of practice specialty toolkits such as in proteomics or ecology.
> - GTN provides training material access 24x7, and regular training intensives.
>     - Crucial capacity building resources help Galaxy grow:
>         - Introductory Galaxy training for users
>         - Specialised training 
>             - for users doing common analyses in a range of specialties 
>             - for tool wrapper and software developers
>             - for trainers
>             - for server system administrators   
{:  .details}

     

### Communities: Definition for this guide.

- Communities sustain Galaxy, although they are not always easy to define. They appear when participants organise themselves to get things done in the project. 
- They are a core element of the success of any large open project and feature heavily in the main "People" part of the field guide, but the term is used here in a restrictive way, to refer to communities that organise activities for interested participants. If a community exists but does not organise activities, it is hard to engage so not so useful in the "People" part of the guide.
- One major exception is arguably the largest "community", comprising all the researchers who use Galaxy for open science analysis. They do not organise their own activities, but they are a very important group, for demonstrating productivity for grants, and as a major source of individuals who join in and contribute to active communities. That group and the research they do using Galaxy are probably the most important open science outputs of the entire project, so they appear in the "Products" part of the guide.
- Communities form themselves through regular interactions within the project.             Participants are more likely to interact with others who share a distinguishing interest. That interest might be proteomics, muon science, India or almost anything else that attracts participants to activities. A functional definition is that a community is important when those interactions become sufficiently frequent and distinct from existing communities to have their own Hub page.
- Some communities have been initiated by the project team, but many have come from participants organising themselves around an interest. 
- As a definition for this guide, communities organise and publicise activities for interested participants through the Hub. Leadership and participation are responsive to the real needs of the community because they are motivated and organised by participants. The project provides support and publicity, but the community must sustain itself with minimal external resources.
- Communities are open, for participation or to initiate new community activities, such as:
            Developers contributing to the code repository are an example of an important community.
            Communities of practice maintain "best practice" tools, workflows and training for specific kinds of open science analysis.
            Regional communities support local training and other distributed activities.
- Open source shared values ensure that community activities are safe, productive and enjoyable:
            inclusive, participatory, professional behaviours are encouraged
            The project explicitly strives to welcome and engage contributors and users
            Community members are encouraged to help make the project better.
- Communities are a core resource, adding substantial value to project grant resources. The project succeeds because it is a self governing collaboration between many specialised communities of contributors and investigators. 

# A field guide to the project ecosystem 

## 2. [People and their interactions](../people/tutorial.html)

## 3. [Resources used in project activities](../resources/tutorial.html)

## 4. [Project outputs and deliverables](../outputs/tutorial.html)

## 5. [Community development success stories](../stories/tutorial.html)

