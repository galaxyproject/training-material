---
layout: tutorial_hands_on
title: FAIR-by-design methodology
abbreviations:
  FAIR: Findable, Accessible, Interoperable, Reusable
zenodo_link: https://zenodo.org/records/11548062
questions:
  - How to develop FAIR learning materials?
objectives:
  - Define FAIR learning objects
  - Adapt and mix FAIR learning objects
  - Identify licenses and attribute correspondingly
  - Structure comprehensive learning materials
  - Manage file formats and tools
  - Define metadata using a schema
  - Create and publish FAIR-by-Design learning materials
  - Collaborate with other instructors
  - Assess FAIR-ness of existing learning objects
time_estimation: 60M
key_points:
  - FAIR data are data which meet principles of findability, accessibility, interoperability, and reusability (FAIR).
  - FAIR data are as open as possible, and as closed as necessary.
  - The main objective of FAIR data is to increase data reuse by researchers.
tags:
  - 'FAIR-by-Design Learning Materials'
  - "FAIR Learning Objects"
  - FAIR
  - Metadata
  - Learning objectives
  - Blooms taxonomy
priority: 1
contributions:
  authorship:
    - sonjafiliposka
    - amisev1971
    - korvoj
  funding:
    - skills4eosc

---


Welcome to the FAIR-by-Design GTN course.

Below are the main stages of the FAIR-by-Design Methodology that will help guide you on your journey of creating FAIR learning materials.


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Main stages

## Stage 1 - Prepare
> Before anything else, preparation is the key to success!
{: .quote author="Alexander Graham Bell"}

### First things first: What is FAIR?


> <question-title>Findable</question-title>
> The editable learning material has a unique and persistent identifier (PID) and is described with sufficiently detailed metadata.
{: .question}


> <question-title>Accessible</question-title>
>
> The human and machine readable metadata and object are stored in a trusted repository with clear authentication and authorization procedures.
>
{: .question}


> <question-title>Interoperable</question-title>
> The metadata describing the learning material follows a the RDA minimum metadata schema combined with agreed-upon controlled vocabularies.
> Formal, accessible, shared, and broadly applicable language(s) and format(s) are used to develop the material.
{: .question}


> <question-title>Reusable</question-title>
>The learning material has a clear usage license (CC-BY-4.0 recommended) and accurate information on provenance.
{: .question }

### Adopt a metadata schema
  <a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%201%20%E2%80%93%20Prepare/02-Preparing%20FAIR%20Learning%20Objects/02-Preparing%20FAIR%20Learning%20Objects/#rda-minimal-metadata-for-learning-resources" class="btn btn-primary btn-lg btn-block">The RDA Minimal Metadata Set for Learning Resources</a>


{% include _includes/tab-choices.html option1="Descriptive info fields" option2="Access info fields" option3="Educational info fields" default="Descriptive info fields" text="Metadata fields categories" %} 

<div class="Descriptive-info-fields" markdown="1">
- **Title** =	The human readable name of the resource
- **Abstract / Description** =	A brief synopsis about or description of the learning resource
- **Author(s)** =	Name of entity(ies) authoring the resource
- **Primary Language** =	Language in which the resource was originally published or made available
- **Keyword(s)** =	Keywords or tags used to describe the resource
- **Version Date** =	Version date for the most recently published or broadcast resource
</div>
<div class="Access-info-fields" markdown="1">
- **URL to Resource** =	URL that resolves to the learning resource or to a "landing page" for the resource that contains important contextual information including the direct resolvable link to the resource, if applicable.
- **Resource URL Type**	 = Designation of the identifier scheme used for the resource URL, e.g., DOI, ARK, Handle
- **License** =	A license document that applies to this content, typically indicated by URL
- **Access Cost**	= Choice stating whether or not there is a fee for use of the resource (yes, no, maybe)
</div>
<div class="Educational-info-fields" markdown="1">
- **Target Group (Audience)** =	Principal users(s) for which the resource was designed
- **Learning Resource Type** =	The predominant type or kind that characterizes the learning resource
- **Learning Outcome** =	Descriptions of what knowledge, skills or abilities a learner should acquire on completion of the resource
- **Expertise (Skill) Level** =	Target skill level in the topic being taught; example values include beginner, intermediate, advanced
</div>



More on the metadata here {% cite hoebelheinrich_2022 %}

More on FAIR learning objects definition: 
    - [FAIR instructional design skills](https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%201%20%E2%80%93%20Prepare/01-FAIR%20skills%20%26%20principles/01-FAIR%20skills%20%26%20principles/#fair-instructional-design-skills)
    - [FAIR guiding principles](https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%201%20%E2%80%93%20Prepare/01-FAIR%20skills%20%26%20principles/01-FAIR%20skills%20%26%20principles/#fair-guiding-principles)
    - [FAIR learning objects](https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%201%20%E2%80%93%20Prepare/01-FAIR%20skills%20%26%20principles/01-FAIR%20skills%20%26%20principles/#fair-learning-objects)

### Start Ideating...

<div class="row">
  <div class="col-sm-4">
   <!-- <div class="card text-white bg-secondary mb-3" > -->
   <div class="card bg-light mb-3" >
      <!-- <div class="card-header text-white"> -->
      <div class="card-header"><i class="fa fa-fast-backward" aria-hidden="true"></i> Think backward</div>
      <div class="card-body">
        <h5 class="card-title">Step 1</h5>
        <p class="card-text">What are your desired effects, i.e. learning outcomes...</p>
      </div>
    </div>
  </div>

  <div class="col-sm-4">
    <div class="card bg-light mb-3" >
      <div class="card-header"><i class="fa fa-fast-backward" aria-hidden="true"></i> Think backward</div>
      <div class="card-body">
        <h5 class="card-title">Step 2</h5>
        <p class="card-text">How are you going to assess their achievement....</p>
      </div>
    </div>
  </div>

  <div class="col-sm-4">
    <div class="card bg-light mb-3" >
      <div class="card-header"><i class="fa fa-fast-backward" aria-hidden="true"></i> Think backward</div>
      <div class="card-body">
        <h5 class="card-title">Step 3</h5>
        <p class="card-text">How should you structure the material to reach them...</p>
      </div>
    </div>
  </div>

</div>

More on backward learning:
    [Steps of the backward learning process](https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%201%20%E2%80%93%20Prepare/02-Preparing%20FAIR%20Learning%20Objects/02-Preparing%20FAIR%20Learning%20Objects_cont/#backward-instructional-design-process)


## Stage 2 - Discover

## Stage 3 - Design

### Time to brainstorm

### Structure is everything

### How to develop the learning content


# Conclusion

