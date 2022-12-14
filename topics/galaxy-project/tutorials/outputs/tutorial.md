---
layout: tutorial_hands_on
logo: "Galrollup"
title: "4: Galaxy project outputs"
summary: "Deliverables part of the field guide. Work in progress. Please help make it better?"
questions:
  - "What are the main identifiable project outputs?"
  - "How can Galaxy offer 'free' analysis computational resources?"
  - "How does the open source community support Galaxy"
  - "How are new project outputs or resources created?"
objectives:
  - "Understand Galaxy outputs and open science services"
  - "Understand how different communities work together to make things happen"
  - "Understand opportunities for engagement and contributing your skills"
 
time_estimation:
subtopic: anatomy
key_points:
  - "The Galaxy project creates many highly valued open science resources"
  - "Unlike most open source projects, Galaxy offers important services and community activities in addition to source code"
  - "Participants engage with, and add value to the project outputs on their own individual terms"
 
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
> - The most important outputs of the Galaxy project are grouped arbitrarily here, and there are many overlaps, because Galaxy grows organically through collaborations, rather than by design.
> - The project continues to expand rapidly, so this module will need updating regularly
> - The Hub provides much more detail about many of the same structures and their activities, but this material is designed to provide simplified views of the project, so the Hub becomes easier to navigate.
> - This is an attempt at a kind of field guide to the ecosystem generating those Hub activities, for the use of participants trying to navigate it.
> 
{: .comment}

> <agenda-title>Field Guide Part 3. Project outputs and impact</agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

# 3. Project outputs

Source code is the primary deliverable for the project, because most other project activities build on it. It is also probably the most important resource for the project in a sense, but for simplicity, it is an output in this guide. The most widely used outputs of the project for researchers, are free services for training and analyses, made possible by the source code, and all the people and the resources described in this guide.

The downstream impact of Galaxy is important but hard to quantify, such as:
    - increased research outputs from analyses in open science represent increased productivity for researchers. Access to efficient and reliable analysis methods for large, complex data resources, is likely to greatly increase their actual use in research, across many scientific fields.
    - Improved trustworthiness of sharable, replicable computation for analyses.
    

## Source code

The core framework source code is supported by many other project repositories. For example, providing tools for system administrators and developers, and ToolShed code for tool distribution services. These are not very useful outside the project, since they are specific to Galaxy. Their impact on open science is through their support for Galaxy. 

### Do these need to be listed here? Will most readers care about detail?

- Generic analysis framework source code
- ToolShed source code
- Developer and system administrator utilities
    -  planemo, ephemeris, ansible...
- Tool wrappers to "flavour" framework instances
    - +8k variable quality tools in public toolshed
    - Click to install from ToolShed, in any framework server
    - Many well maintained tools from IUC and communities of practice. 
- Your ideas here please?

## Open science analysis services

- "Free" project supported services
    - usegalaxy.*
    - Large Australian, European and US research infrastructure allocations
    - Tens of thousands of users
    - Professional user support
    - Stress testing framework code and tools
- 100+ specialised public instances
- Unknown number of private installations
- Your ideas here please?

## Capacity building: training resources and services

Providing training to build community capacity is an essential activity for the project, to ensure wide, well managed deployment and long term sustainability. 

 - GTN user training integrated directly into the Galaxy user interface, helps new users to gain the skills they need to be productive and efficient. 
 - Training system administrators helps support the public usegalaxy.* and the many private servers that operate in academic and commercial laboratories. 
 - Training for external developers makes it easier for them to contribute efficiently, improving Galaxy code and wrapping new tools.
- The Galaxy Trainng Network (GTN) is central to building community capability.
- Offers free training to enhance global open science research capacity.
    - Generic aspects of using Galaxy for new users
    - Specific kinds of analyses with common types of open data.
    - System administrators are key to running reliable framework services
    - Software developers can contribute more efficiently with appropriate training
    - Trainers can learn how to prepare material for new GTN topics
- Your ideas here please?

## Downstream: repeatable open science analyses

The "black box" of a closed source software package hides the source code, assumptions and methods used in an analyses, so they cannot be readily scrutinised. The software is not readily accessible for independent replication without fees, so the scientific trustworthiness of the results cannot be tested readily. It is said that *many eyes make bugs shallow*.

Closed source code may be perfect, but unfortunately, experience suggests that all complex software contains errors, many of which can only be found after widespread and thorough independent scrutiny. This is equally true of expensive commercial software, and of open source software. Open source package assumptions, methods and code are readily accessible for review, testing and improvement. Open projects encourage and facilitate scrutiny and replication, in order to decrease the risks to scientific integrity and trustworthiness from hidden coding or methodological errors. 

- Analysis of provable scientific integrity are arguably the most important project deliverable
- 10k+ publications
- Tens of thousands of scientists trained
- Millions of jobs run.
- Your ideas here please?

### Related field guide content for further reading

- [Introduction and definitions](../introduction/tutorial.html)
- [People and their interactions](../people/tutorial.html)
- [Resources used in project activities](../resources/tutorial.html)
- [Community development success stories](../stories/tutorial.html)
