---
layout: tutorial_hands_on
logo: "Galrollup"
title: "A field guide to the Galaxy project for participants"
summary: "People, resource inputs and project outputs. Work in progress. Please help make it better?"
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

> <comment-title>Work in progress</comment-title>
> - First draft to try to get the structure to make sense.
> - Needs many contributors to make it useful
> - Please add what's missing and fix what's broken
> - Most of the headings are stubs waiting to be edited and extended 
> - Hoping this way of splitting the material makes for useful reading?
> - Add your story or stories to the stories tutorial too please!
> - Ross has strong opinions. 
> - Many of them are probably wrong but he doesn't know which ones yet. 
> - Please feel free to contribute your own, to make this more useful to future readers.
>
{: .comment}


> <comment-title>Note to readers</comment-title>
> - This module introduces the most important parts of the Galaxy project
> - They are grouped arbitrarily and there are many overlaps.
> - Galaxy grows organically through collaborations, rather than by design.
> - The project continues to expand rapidly, so this module will need updating regularly
> - The Hub provides details about many of the same structures
> - This material is designed to provide simplified views of the project
> - A kind of field guide to the ecosystem generating those Hub activities.
> - For the use of participants trying to navigate it.
> 
{: .comment}

> <agenda-title>Contents</agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Many ways of seeing Galaxy

- Galaxy is a complicated global enterprise.
- Individuals engage in many different ways, according to their interests.
- For example, it is seen as:
    - A convenient and reliable way to share complex open science analyses by users
    - Inter-dependent open source code repositories by developers.
    - A global open science collaboration by investigators
- The project is all of these at the same time and more
- Seeing the big picture requires multiple perspectives.
- A separate tutorial on how activities were started, helps illustrate how the project grows.
- This tutorial outlines the major components as a field guide.
- The [Hub](https://galaxyproject.org/) can provide much more detail, particularly once you know what you are looking for.

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
 
# Top level categories used in this guide

- The project exists in the context of a global open science ecosystem
- That context is an essential part of any description. 
- Many possible top level divisions for describing those complexities
- Three top level categories are related in a simple model:

### People + resource inputs = project outputs

- Interacting high level project components
- Click on the "Details" tabs for more information about these.
- The field guide separates those three components artificially.
- The project depends on them all, efficiently working together.
- Galaxy flourishes, because all these parts govern themselves in an efficient collaboration.

 
> <details-title>People working together make the project possible</details-title>
> - Open source is a very productive way of delivering software
> - Shared values and participatory self governance help people get things done
> - Galaxy has to be a self-governed collaboration
>     - Many independent institutions and investigators 
>     - No single institution has control
> - Individuals engage with the project in their own way, such as:
>     - Analysing and sharing experimental data
>     - Working as professional member of the core team
>     - Contributing skills, support, code and ideas
>     - Maintaining best practice tool kits and workflows for specific fields
>     - Bringing new activities to the collaboration
>     - Providing community and project leadership
>     - Contributing to project governance
{: .details}


> <details-title>Resources from grants sustain the communities, team, collaborations and services</details-title>
> - Core "corporate" project services like outreach, communication and administration.
> - Source code management, user support, usegalaxy.*, GTN and other services.
> - Computational resource allocations for free services
>     - adding substantial value to the source
>     - Free analysis and training services
> - Professional staff
>     - outreach, communication, project management.
>     - user support, system administration and software engineers.
> - Without these, the source code repository would still be a valuable resource
> - They add value, in the form of highly valued resources for scientists.
{: .details}

> <details-title>Open source community resources are essential parts of project code</details-title>
> - Galaxy source code depends on thousands of other open source projects
> - Flavours for Galaxy servers depend on thousands of open source command line packages
> - Galaxy builds on and adds value to all those resources.
> - It makes command line packages interoperable and reproducible, through a uniform interface.
> - Provides a GUI and sharing for users, without requiring effort from the package developer.
{:  .details}

> <details-title>Project outputs: things the project delivers for open science and scientists</details-title>
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

## 1. People: Individuals, collaborations, communities and governance
- This field guide begins with people and their interactions in communities and structures
- Seeing how they have organised themselves provides one useful perspective on the project.  
- Many individuals are part of more than one of the organisations and communities below
- People provide the project inputs
- Investigators bring resources.
- Contributors volunteer their skills and time to make things even better
- The sum of their effort produces the project outputs.
- Focus here is on showing *how* they get things done, rather than exactly *what* they do.
- For those details, see the [Hub](https://galaxyproject.org/)
- If any activity of the project looks interesting to you, please get in touch!

> <tip-title>Communities in more detail</tip-title>
> - Communities sustain Galaxy, although the idea is abstract and fluid.
> - They support all successful open projects 
> - Theories about what they represent:
>   - Arise when participants organise themselves, to get things done.
>   - Subgroups of individuals interacting regularly in shared project activity
>   - An important part of the success of any open project.
> - Benefits for community contributors:
>   - Opportunities for personal and professional development.  
>   - Extending skills through supported problem solving.
>   - Satisfaction in helping make open science better.
> - Participants interact more often, if they share an important distinguishing interest.
> - New communities appear when interactions are frequent and distinct from existing communities
>     - Functional communities organise publicised activities for interested participants
>     - Leadership and participation are responsive to the real needs of the community.
>     - The project can help with publicity, but the community must sustain itself.
>     - Anyone can participate in, or initiate new community activities.
>     - Developers contributing to the code repository are an example of an important community.
>     - *Communities of practice* maintain "best practice" tools, workflows and training for specific kinds of open science analysis.
>     - Regional communities support local training and other distributed activities.
> - Open source shared values help make communities safe, productive and enjoyable for everyone
>     - inclusive, participatory, professional behaviours are encouraged
>     - Project strives to welcome and engage contributors and users
>     - Community members are encouraged to help improve the project
> - Communities are a core resource, adding substantial value to project grant resources.
> - The project succeeds because it is a global, self governing, collaborating community.
{:  .tip}


### Communities of Practice

- Participants sharing an interest in a specific field form field-focussed communities
- Galaxy began with a focus on tools for analysing commodity sequencing data.
- Many new fields have adopted it for their own data and analyses
- Some of these are described in the companion stories tutorial.
- Best practice tool kits are maintained by contributors working in many fields:
    - [Proteomics](../../../proteomics)
    - [Climate science](../../../climate)
    - [Microbiome](../../../metagenomics)
    - public health 
    - small server administrators
    - add your own here...

### Regions

- Participants also tend to create geographically local communities and activities
- North America
- Europe
- Australia
- India
- add your own here...

### Developers

- Communities of developers are recognisable in the various source code repositories
- Interact through github issues, pull requests, and matrix chat channels

### Collaborating investigators

- Independent principal investigators
- Multiple independent grants that support interacting project deliverables
- Support source code, hub, usegalaxy.* analysis services, Galaxy Training Network, communities...
- Provide resources at scale for "free" Galaxy analysis services


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
 
- Galaxy is an open project
    - contributors are welcomed and supported
    - Project success and sustainability depends on contributions from community members.
- Project complexity offers many opportunities for contribution
- Suggestions are always welcomed to make it easier to navigate.
- In the past, community members turn into highly valued contributors
- Finding their own particular way to join in according to their situation.
- More recently, effort is being devoted to documenting and clarifying the process, to make it far easier to navigate, for all kinds of contributors. 
- This training module is part of that effort
- See the contributor stories material to learn how useful contributions have been made.

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
    - NHGRI AnVIL


### Developers

- Source code. Team and many independent contributors
- Tool wrappers. IUC. Communities of Practice
- Analysis packages - external authors.

### Educators

- GTN!

### Governance and organisation

- See hub.
- Self governing, open source style
- Unusual challenge: Respect independent grant holder responsibilities
- Executive Board
- GGSC
- Working groups
- Road maps
- IUC/IWC/IDC and other project initiatives - open but initiated by the Team?
- ....
- Your ideas here please?

### Your ideas for people things to describe here please


## 2. Resource inputs
- Many different kinds of resources are needed to support all the "free" services
- Large scale computational infrastructure and professionals 
- Individual skilled contributors
- Engaged communities
- Other interacting open communities

### Source code
- Open source community provides
    - Galaxy source code dependencies
    - Command line analysis package code for tools

### Computational infrastructure, support and analysis services

- Australia
- Europe
- USA

### Collaborating grants
- Would PI's want to advertise their deliverables here?
- One way of letting potential collaborators know what you're doing?
- Does not need budgets :)
- AnVIL
- Cancer....
- Your ideas here please?

### Community contributors
- Very important resource
- Bug fixes, ideas...
- 

## 3. Project outputs

- Most visible parts of the project to researchers using Galaxy for analyses
- Made possible by all the people and the resources
- Downstream products of Galaxy are important but hard to quantify
    - Increased productivity in open science analysis
    - Improved trustworthiness of sharable, replicable computation for analyses.
    - Scientific discovery

### Source code

- Generic analysis framework
- Utilities - planemo etc
- Your ideas here please?

### Tool wrappers to "flavour" framework instances

- +8k tools in public toolshed
- Your ideas here please?

### Free analysis services

- usegalaxy.*
- Large Australian, European and US research infrastructure allocations
- Tens of thousands of users
- Professional user support
- Stress testing framework code and tools
- Your ideas here please?

### Capacity building: training resources and services

- GTN offers training to enhance open science research capacity
    - Using Galaxy itself
    - Specific kinds of analyses
    - System administrators
    - Software developers 
    - Trainers
- Your ideas here please?

### Downstream: repeatable open science analyses

- Most important and probably largest deliverable
- 10k+ publications
- Tens of thousands of scientists trained
- Your ideas here please?
