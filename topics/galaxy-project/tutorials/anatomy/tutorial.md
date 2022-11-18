---
layout: tutorial_hands_on
logo: "Galrollup"
title: "The people, values, resources and source code of Galaxy"
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
  - "Galaxy is a global collaboration"
  - "Write once, use often for framework features like analysis reproducibility"
  - "Running code as a free service provides effective stress testing for defects"
  - "Embeded in the much bigger open source and open science ecosystem"
  - "Started in genomics but adapted for use by many data rich disciplines"
  - "New contributors and investigators are highly valued"
 
contributors:
  - fubar2
---

> <comment-title></comment-title>
>
> Ross has strong opinions. 
> Many of them are probably wrong but he doesn't know which ones yet. 
> Please feel free to contribute your own, to make this more useful to future readers.
>
{: .comment}


> <agenda-title></agenda-title>
>
> - This module introduces some of the most important components of the Galaxy project
> - These are necessarily somewhat arbitrary and overlapping.
> - They have grown and continue to grow organically over time rather than by design.
> - The project continues to expand rapidly, so this module will need updating regularly
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Alternative ways of seeing Galaxy

- Galaxy is a big, complicated project with many moving parts making it hard to comprehend completely
- What you see depends on how you look at it, so it is:
    - A convenient and reliable way to get complex analyses done for users
    - Inter-dependent open source code repositories for developers.
    - A global collaboration for investigators
- All models of Galaxy are wrong, but some may be useful.
- This module outlines the major components, and how they fit together
- Drafted by a human, it represents only one of many possible explanations
- Most participants will be associated with more than one of the parts described here
- They all need to be considered together but it's easier to list them individually

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
 
# Three related top level categories for this guide

People + resource inputs = project outputs

- An equation that summarises the Galaxy project
- More explanation in the optional material below
- People are the key.
- They are considered in detail in the first section of the field guide

> <details-title>People working together make the project possible</details-title>
> - That might seem trite, but it's how open source succeeds.
> - Shared values help people get things done
> - Self governed collaboration because no single "owner"
>     - Many independent institutions and investigators 
>     - No single institution has control
> -  Individuals engage with the project on their own terms
>     - Leading and participating in community activities
>     - Communities of practice
>     - Collaborating investigators with grant resources
>     - Contributing training material, documentation, source code 
>     - Participation and leadership in governance.
{:  .details}


> <details-title>Grant resource inputs help sustain the communities, team, collaborations and services</details-title>
> - Core "corporate" project services like outreach, communication and administration.
> - Source code management, user support, usegalaxy.*, GTN and other services.
> - Computational resource allocations for free services
>     - adding substantial value to the source
>     - Free analysis and training services
> - Professional staff
>     - outreach, communication, project management.
>     - user support, system administration and software engineers.
{:  .details}

> <details-title>Project code and services depend on resources from open source community</details-title>
> - Galaxy source code depends on interacting resources provided by thousands of other open source projects
> - Flavours for Galaxy servers depend on command line software packages wrapped as tools
{:  .details}

> <details-title>Project outputs - things the project delivers for open science and scientists</details-title>
> - Source code
> - Online analysis services
> - ToolShed automated tools for servers
> - Best practice standards and automation for writing tool wrappers and workflows - IUC
> - Communities of practice specialty toolkits such as in proteomics or ecology.
> - GTN open access to training material any time and to periodic training intensives.
>     - Provides capacity building resources:
>         - Introductory Galaxy useage training
>         - Specialised training 
>             - for users doing common analyses in a range of specialties 
>             - for tool wrapper and software developers
>             - for trainers
>             - for server system administrators  
{:  .details}

# Field guide to the project ecosystem 

## People: Individuals, collaborations, communities and governance
- Starts with the major ways people have organised themselves in the project
- Many individuals are part of many of the organisations and communities below
- How people engage and work together is central to understanding how Galaxy works
- People bring the resource inputs and project outputs are what people produce.
- For this module, how they get things done is of more interest than exactly what they do.

### Communities 

- Participants more likely to interact when they share a particular interest
- Informal and formal communities arise from these interactions.
- They are an important part of any open project 
- Open source shared values drive participation
    - inclusive, participatory, professional behaviours are encouraged
    - Project strives to welcome and engage contributors and users
    - Community members are encouraged to help improve the project
- A community might define itself as participants sharing a particular interest.
    - Recognisable communities organise activities for interested participants
    - Code contributors collaborate as one of many interacting open source communities
    - Communities of practice contribute consensus best practice, tool kits and training material.
- Project growth and sustainability are dependent on providing value through community engagement
- Communities are a core resource, adding substantial value to project grant resources.
- Some important communities are described individually below

### Communities of Practice

- Communities developing and supporting specialised fieldd
    - proteomics
    - climate science
    - public halth 
    - small scale admins
        - for people responsible for keeping Galaxy services running

### Regions

### Developers


### Collaborating investigators

- Independent principal investigators
- Multiple independent grants that support interacting project deliverables
- Support source code, hub, usegalaxy.* analysis services, Galaxy Training Network.
- "Free" Galaxy analysis services are sustained by large research grants


> <details-title>Collaborative resources in more details</details-title>
> - Free analysis consumes substantial research infrastructure and computational resources 
> - Project activities are supported by highly skilled professionals
> - Independent grants operate efficiently for their own goals and can pool resources
> - Efficient collaboration magnifies their individual impact.
> - No single grant could encompass the complexity of institutions and resources involved
> - Some are code focussed, others provide resources such as the GTN and usegalaxy.* services
> - Community engagement and contribution adds import additional value
> - These must work together efficiently for the project to succeed.
{:  .details}
 


### Contributors
 
- Galaxy is an open project, where contributors are welcomed and supported
- Project success and sustainability depends on contributions from community members.
- Project complexity offers many opportunities for contribution
- Effort is needed to make it easier to navigate.
- In the past, many community members have become highly valued contributors by finding their own particular way to join in.
- More recently, effort is being devoted to documenting and clarifying the process, to make it far easier to navigate, for all kinds of contributors. 
- This training module is part of that effort and there is a separate module on contributor pathways and stories to help show how many useful contributions have been made.

> <details-title>Examples for different kinds of skills and interests</details-title>
> - Issues can be raised and pull request submitted for review at the github repository. 
> - The GTN provides integrated tutorials on using Galaxy, and technical capacity building training material, such as on server administration, and on making GTN tutorials. Infrastructure for slides narrated in multiple languages, and hands-on tutorials applying real Galaxy tools to real data are supported.
> - The IUC provides best practice guidance and GTN training modules to help developers build new tool wrappers for Galaxy, and automated infrastructure to help maintain them efficiently.
> - GGSC offers opportunities for motivated community members to initiate and lead community development collaborations that can help the project to expand into activities that help support and expand subcommunities of users.
> - Working groups provide opportunities for the contribution of technical and other skills in focussed project areas ranging from GUI development through to training and outreach.
{:  .details}

### Users

- GUI tool users who do not want to write code
- Bioinformaticians and other users who want to write code


### Institutional adopters

- Institutions: organisations deploying and using private Galaxy services


### Developers
- Source code. Team and many independent contributors
- Tool wrappers. IUC. Communities of Practice
- Analysis packages - external authors.

### Educators

### Governance

## Resource inputs from grants



## Project deliverables

### Code for the generic analysis framework

### Tools to "flavour" framework instances

### Free analysis services

### Capacity building training resources and services

### Downstream: repeatable open science analyses
