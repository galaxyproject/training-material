---
layout: tutorial_hands_on
title: FAIR-by-Design methodology
abbreviations:
  FAIR: Findable, Accessible, Interoperable, Reusable
  OER: Open Educational Repositories
  EOSC: European Open Science Cloud
  CC: Creative Commons
  MVS: Minimum Viable Skills
  TeSS: Training eSupport System
zenodo_link: https://zenodo.org/records/11548062
questions:
  - How to develop FAIR learning materials?
  - How to incorporate the FAIR principles in the instructional design process?
objectives:
  - Design FAIR learning materials
  - Structure FAIR learning materials
  - Create and publish FAIR-by-Design learning materials
  - Assess FAIR-ness of learning objects
time_estimation: 20M
key_points:
  - The FAIR-by-Design Methodology ensures that the output of the instructional design process are FAIR learning materials.
  - Each stage of the FAIR-by-Design Methodology incorporates specific parts of the FAIR principles.
  - The ouput of the FAIR-by-Design Methodology are FAIR learning materials from both the perspective of instructors and of learners.
tags:
  - 'FAIR-by-Design Learning Materials'
  - "FAIR Learning Objects"
  - FAIR-by-Design Methodology
priority: 1
contributions:
  authorship:
    - sonjafiliposka
    - amisev1971
    - CarolinLeister
    - korvoj
    - domlgreen
    - isouyioul
    - achim-kit
    - agCleo
    - CRHAD26
  funding:
    - skills4eosc
follow_up_training:
 - type: external
   title: Skill4EOSC FAIR-by-design Methodology Training of Trainers
   link: "https://learning.skills4eosc.eu/course/view.php?id=19"
lang: "en"
level: "Introductory"
extra:
  description: "Short how to guide that guides you through the stages of the FAIR-by-Design methodology without any specific choice on tools and formats."
  target_audience: "Instructional designers"
  URL_resource_type: "Permanent URL"
  learning_resource_type: "self-paced microlearning unit"
  access_cost: "No"
  attribution: "This microlearning unit is based on Filiposka, Sonja, Mishev, Anastas. (2024). FAIR-by-Design Microlearning (1.0.0). https://doi.org/10.5281/zenodo.11548062"
---


Welcome to the FAIR-by-Design Methodology Microlearning GTN adapted tutorial.

> <agenda-title></agenda-title>
>
>Below are the main stages of the FAIR-by-Design Methodology {% cite filiposka_2023_8305540 %} that will help guide you on your journey of creating FAIR learning materials.
>
>![Visual representation of the FAIR-by-Design Methodology flow with its separate stages](../../images/methodology.png)
>
>The FAIR-by-Design Methodology is defined in [a number of stages](#main-stages) that help you incorporate the FAIR principles into the backward learning instructional design process so that the final output are FAIR learning materials from both learners and instructors perspective:
>
>1. [Prepare](#stage-1---prepare)
>	- Do you understand FAIR and its implications
>	- Define purpose, learning objectives, target audience
>2. [Discover](#stage-2---discover)
>	- Find existing resources
>	- Identify potential for reuse
>3. [Design](#stage-3---design)
>	- Define the syllabus and structure
>	- Identify granularity
>	- Define facilitation materials
>4. [Produce](#stage-4---produce)
>	- Develop content using common apps and file formats
>	- Define machine readable metadata
>	- Perform internal QA check
>5. [Publish](#stage-5---publish)
>	- Define license and other related info
>	- Release for instructors and learners
>	- Enable feedback gathering
>6. [Verify](#stage-6---verify)
>	- Final QA check
>	- Use gathered feedback for continuous improvement
>7. [Continuous Improvement](#stage-7---continuous-improvement)
>	- Create a list of potential improvements
>	- Choose a set of improvements to be implemented
>	- Start a new release cycle
>
>Review each stage and the essential steps that it includes.
>
>Essential information about the original microlearning unit that this adaptation is based on can be found in the [About](#about) section.
>
{: .agenda}

If you have any questions at any stage on your journey do not hesitate to [contact the FAIR-by-Design methodology team](mailto:sonja.filiposka@finki.ukim.mk).

May all your materials be FAIR!


# Main stages

## Stage 1 - Prepare

> Before anything else, preparation is the key to success!
{: .quote author="Alexander Graham Bell"}

### First things first: What is FAIR?

><question-title> What is FAIR? </question-title>
>
>> <solution-title>Findable</solution-title>
>> The editable learning material has a unique and persistent identifier (PID) and is described with sufficiently detailed metadata.
>{: .solution}
>
>
>> <solution-title>Accessible</solution-title>
>> The human and machine readable metadata and object are stored in a trusted repository with clear authentication and authorization procedures.
>{: .solution}
>
>
>> <solution-title>Interoperable</solution-title>
>> The metadata describing the learning material follows a the RDA minimum metadata schema combined with agreed-upon controlled vocabularies.
>> Formal, accessible, shared, and broadly applicable language(s) and format(s) are used to develop the material.
>{: .solution}
>
>
>> <solution-title>Reusable</solution-title>
>>The learning material has a clear usage license (CC-BY-4.0 recommended) and accurate information on provenance.
>{: .solution }
{: .question }

> <details-title>More details on the FAIR principles</details-title>
> [Follow the FAIR in a nutshell tutorial available on GTN]( {% link topics/fair/tutorials/fair-intro/tutorial.md %} )
> 
> [Follow the FAIR Galaxy Training Material tutorial available on GTN]( {% link topics/fair/tutorials/fair-gtn/tutorial.md %} )
{: .details}

### Adopt a metadata schema

If you are not using a discipline specific metadata schema, then, to ensure that your learning materials are appropriately described using a common approach, you should adopt the:

  <a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%201%20%E2%80%93%20Prepare/02-Preparing%20FAIR%20Learning%20Objects/02-Preparing%20FAIR%20Learning%20Objects/#rda-minimal-metadata-for-learning-resources" class="btn btn-primary stretched-link">Go to the RDA Minimal Metadata Set for Learning Resources description details</a>

The following are the fields described in the RDA Minimal Metadata Set divided into three categories:

> <details-title>Descriptive info metadata fields</details-title>
> - **Title** =	The human readable name of the resource
> - **Abstract / Description** =	A brief synopsis about or description of the learning resource
> - **Author(s)** =	Name of entity(ies) authoring the resource
> - **Primary Language** =	Language in which the resource was originally published or made available
> - **Keyword(s)** =	Keywords or tags used to describe the resource
> - **Version Date** =	Version date for the most recently published or broadcast resource
{: .details}

> <details-title>Access info metadata fields</details-title>
> - **URL to Resource** =	URL that resolves to the learning resource or to a "landing page" for the resource that contains important contextual information including the direct resolvable link to the resource, if applicable.
> - **Resource URL Type**	 = Designation of the identifier scheme used for the resource URL, e.g., DOI, ARK, Handle
> - **License** =	A license document that applies to this content, typically indicated by URL
> - **Access Cost**	= Choice stating whether or not there is a fee for use of the resource (yes, no, maybe)
{: .details}

> <details-title>Educational info metadata fields</details-title>
> - **Target Group (Audience)** =	Principal users(s) for which the resource was designed
> - **Learning Resource Type** =	The predominant type or kind that characterizes the learning resource
> - **Learning Outcome** =	Descriptions of what knowledge, skills or abilities a learner should acquire on completion of the resource
> - **Expertise (Skill) Level** =	Target skill level in the topic being taught; example values include beginner, intermediate, advanced
{: .details}

More on the RDA minimal metadata schema here {% cite hoebelheinrich_2022 %}

> <tip-title>Using metadata in GTN tutorials</tip-title>
>
>- When developing materials in GTN, this information should be included in the Tutorial.md file header metadata.
>- Most of the fields are already defined in the GTN Tutorial metadata schema. The ones that are missing can be added using the "extra" field.
>
> > <details-title>More details on the GTN Metadata</details-title>
> > [Follow the GTN Metadata tutorial]( {% link topics/contributing/tutorials/schemas/tutorial.md %} )
> {: .details}
>
{: .tip}


><details-title>Learn more about FAIR learning objects</details-title>
>- [FAIR instructional design skills](https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%201%20%E2%80%93%20Prepare/01-FAIR%20skills%20%26%20principles/01-FAIR%20skills%20%26%20principles/#fair-instructional-design-skills)
 >- [FAIR guiding principles](https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%201%20%E2%80%93%20Prepare/01-FAIR%20skills%20%26%20principles/01-FAIR%20skills%20%26%20principles/#fair-guiding-principles)
 >- [FAIR learning objects](https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%201%20%E2%80%93%20Prepare/01-FAIR%20skills%20%26%20principles/01-FAIR%20skills%20%26%20principles/#fair-learning-objects)
{: .details}

### Start Ideating

<div class="row">
  <div class="col-sm-4">
   <!-- <div class="card text-white bg-secondary mb-3" > -->
   <div class="card bg-light mb-3" >
      <!-- <div class="card-header text-white"> -->
      <div class="card-header"><i class="fa fa-fast-backward" aria-hidden="true"></i> Think backward</div>
      <div class="card-body">
        <h5 class="card-title">Step 1</h5>
        <p class="card-text">What are your desired effects, i.e. learning outcomes?</p>
      </div>
    </div>
  </div>

  <div class="col-sm-4">
    <div class="card bg-light mb-3" >
      <div class="card-header"><i class="fa fa-fast-backward" aria-hidden="true"></i> Think backward</div>
      <div class="card-body">
        <h5 class="card-title">Step 2</h5>
        <p class="card-text">How are you going to assess the learners achievement?</p>
      </div>
    </div>
  </div>

  <div class="col-sm-4">
    <div class="card bg-light mb-3" >
      <div class="card-header"><i class="fa fa-fast-backward" aria-hidden="true"></i> Think backward</div>
      <div class="card-body">
        <h5 class="card-title">Step 3</h5>
        <p class="card-text">How should you structure the material to provide effective learning experience?</p>
      </div>
    </div>
  </div>

</div>
><tip-title>Learn more about the backward learning process</tip-title>
> - [Steps of the backward learning process](https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%201%20%E2%80%93%20Prepare/02-Preparing%20FAIR%20Learning%20Objects/02-Preparing%20FAIR%20Learning%20Objects_cont/#backward-instructional-design-process)
{: .tip }

### Define

> <question-title>Purpose</question-title>
> When and how the learning materials can be used and for what purposes?
{: .question }

> <question-title>Target Audience</question-title>
> Is there anything specific that needs to be taken into account, such as cultural context?
{: .question }

> <question-title>Prerequisites</question-title>
> What does the target audience need to know or understand before starting the learning process?
{: .question }

> <question-title>Scope</question-title>
> Is it going to be a single learning unit, or a group such as a course?
{: .question }

> <question-title>Learning Objectives</question-title>
> What competences will be gained after successful completing of the learning process?
>
> ><tip-title>Be SMART</tip-title>
> >Objectives should be **s**pecific, **m**easurable, **a**ttainable, **r**elevant and **t**ime-bound.
> {: .tip}
>
> ><tip-title>Use Blooms Taxonomy</tip-title>
> >Formulate the objectives as actionable verb + observable knowledge, skill, attitude, behavior or ability.
> {: .tip}
> 
> [Read more](https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%201%20%E2%80%93%20Prepare/02-Preparing%20FAIR%20Learning%20Objects/02-Preparing%20FAIR%20Learning%20Objects_cont/#defining-learning-objectives)
{: .question }

> <details-title>More details on ideating learning outcomes and related information</details-title>
> [Read about defining intended Learning Outcomes in the Design and plan session, course, materials tutorial available on GTN]( {% link topics/contributing/tutorials/design/tutorial.md#define-intended-learning-outcomes-los %} )
{: .details}

<a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_Book/4%20-%20FAIR-by-design%20learning%20materials%20creation/4.1%20-%20Workflow%20stages%20description/411-prepare/" class="btn btn-dark text-white stretched-link">Go to the full description of FAIR-by-Design Methodology: Prepare stage ...</a>

## Stage 2 - Discover

> Greater even than the greatest discovery is to keep open the way to future discovery.
{: .quote author="John Jacob Abel"}
### Get inspired

Reusable materials can be found anywhere. These are just some examples:

> <code-in-title>GTN</code-in-title>
> [GTN learning resources](https://training.galaxyproject.org)
> 
> [TeSS Catalogue by Elixir](https://tess.elixir-europe.org)
{: .code-in}

> <code-in-title>OER</code-in-title>
> [DOAB](https://directory.doabooks.org/)
> 
> [MERLOT](https://www.merlot.org/merlot/index.htm)
> 
> [OASIS](https://oasis.geneseo.edu/index.php)
> 
> [OER Commons](https://www.oercommons.org/)
> 
> [OERTX CORA](https://www.projectcora.org/)
> 
> [GALILEO](https://oer.galileo.usg.edu/)
> 
> [FORRT](https://forrt.org/)
{: .code-in}

> <code-in-title>EOSC</code-in-title>
> [EOSC Training catalogue on the EOSC Marketplace](https://search.marketplace.eosc-portal.eu/search/training?q=*)
> 
> Most EOSC projects have their own training catalogues and/or platforms.
{: .code-in}

> <code-in-title>General</code-in-title>
> [Creative Commons Search](https://search.creativecommons.org/) - content provided under a CC license
> 
> [Zenodo](https://zenodo.org/) - a multi-disciplinary open repository
> 
> [OSF](https://osf.io/) - a free, open research platform
{: .code-in}

### Potential for reuse

> <warning-title> Respect the licenses, to be respected!</warning-title>
> - Materials with non-permissible licenses can be used for inspiration only.
> - Materials with permissible licenses should be reused based on the license rules.
{: .warning}

### Don't forget the multimedia search

Different learners have different learning modalities (read/write, auditory, visual, kinesthetic). To elevate the learning experience you should use all types of multimedia in your learning materials.

[Go to the full description of FAIR-by-Design Methodology: Discover stage ...](https://fair-by-design-methodology.github.io/FAIR-by-Design_Book/4%20-%20FAIR-by-design%20learning%20materials%20creation/4.1%20-%20Workflow%20stages%20description/412-discover/){: .btn.btn-dark.text-white.btn-lg.btn-block}

## Stage 3 - Design

> Design is intelligence made visible.
{: .quote author="Alina Wheeler"}

Now it's time to brainstorm.

### Concept map

Step 1: Build a concept map of your learning materials.

> <details-title>More details on concept maps</details-title>
> [Read about Concept Maps in the Design and plan session, course, materials tutorial available on GTN]( {% link topics/contributing/tutorials/design/tutorial.md#concept-maps %} )
{: .details}

Step 2: Make sure you align your concept map with the MVS profiles.

The aligned MVS profile can help you crystallize the learning objectives using the MVS taxonomy.

Each MVS profile defines a list of technical and soft skills required for the profile. Think on how to incorporate both aspects in your learning materials.

<a href="https://fair-by-design-methodology.github.io/MVS/latest/MVS%20Profiles/Civil%20Servant/civil_servant/" class="btn btn-primary stretched-link">Go to the MVS profiles catalogue</a>

### Structure is everything

><comment-title></comment-title>
>{%icon galaxy-gear %} Use an intuitive logical organisation of all learning materials.
{: .comment}

><comment-title></comment-title>
>{%icon tool-versions %} The goal is for other people to easily reuse a single item (plan, activity, unit, assessment, ...).
{: .comment}

><comment-title></comment-title>
>{%icon galaxy-dataset-map %} Take advantage of a hierarchical structure to combine learning units into larger compositions.
{: .comment}

><hands-on-title>How to organise the files</hands-on-title>
>
> GTN defines a specific hierarchical structure that needs to be followed when developing learning materials.
> 
> > <details-title>Check out the following tutorials</details-title>
> > [How to create a skeleton for a new learning topic]( {% link topics/contributing/tutorials/create-new-topic/tutorial.md#creating-a-new-topic-with-its-own-materials %}  )
> > 
> > [How to create the skeleton of a new tutorial]( {% link topics/contributing/tutorials/create-new-tutorial/tutorial.md %}  )
> {: .details}
> 
{: .hands-on}

><comment-title>Syllabus is ready</comment-title>
> You should by now have the first draft of your <a href='https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%203%20–%20Design/04-Conceptualisation/04-Conceptualisation/#syllabus/'>syllabus</a>. It contains all the fields from the RDA min metadata set plus the high level topics covered by the learning material.
{: .comment}

><tip-title>Available feedback form</tip-title>
> GTN provides a readily available feedback form that is automatically added at the end of each tutorial. The feedback form is used to gather quantitative and qualitative feedback.
> 
> > <details-title> Learn more about feedback in training</details-title>
> > [Follow the Assessment and feedback in training and teachings tutorial available on GTN]( {% link topics/teaching/tutorials/assessment/tutorial.md#dealing-with-feedback %}  )
> {: .details}
> 
{: .tip}

><tip-title>Facilitation guide kit</tip-title>
> A <a href='https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%203%20–%20Design/07-Facilitation/07-Facilitation/'> facilitation guide </a> should help prepare for the actual training.
> 
> The facilitation guide kit includes documentation of the process of organising and running a training. What different people need to have, know and do so that everything runs smoothly.
> 
> > <details-title> Learn more about facilitation</details-title>
> > [Follow the Organizing a workshop tutorial available on GTN]( {% link topics/teaching/tutorials/organize-workshop/tutorial.md %}  )
> > 
> > [Follow the Running a workshop as instructor tutorial available on GTN]( {% link topics/teaching/tutorials/running-workshop/tutorial.md %}  )
> {: .details}
> 
> Another option is to use something like the [TRIPLE project TRAINING TOOLKIT](https://project.gotriple.eu/project-deliverables/triple-training-toolkit/).
> 
{: .tip}

><question-title>What about instructor notes?</question-title>
>  They need to be detailed enough so that anyone can reuse the learning content, especially slides properly.
>  
>  Any specific information relevant for instructors that would like to organise a training based on your GTN tutorial content should be added in Details box.
>  
>  This practice significantly increases the reuse potential of the material.
{: .question}

### How to design the learning content

<a href="https://www.csun.edu/sites/default/files/Holle-Lesson-Planning.pdf" class="btn btn-primary stretched-link">Read more about the Hunter Model.</a>

<div class="row">
  <div class="col-sm-4">
   <div class="card bg-light mb-3" >
      <div class="card-body">
        <h5 class="card-title">{% icon curriculum %} 1. Set Learning Objectives</h5>
        <p class="card-text">
         What is the goal?
        </p>
      </div>
    </div>
  </div>
<div class="col-sm-4">
   <div class="card bg-light mb-3" >
      <div class="card-body">
        <h5 class="card-title">{% icon pref-identities %} 2. Identify Needs</h5>
        <p class="card-text">
        How to get there?
        </p>
      </div>
    </div>
  </div>
<div class="col-sm-4">
   <div class="card bg-light mb-3" >
      <div class="card-body">
        <h5 class="card-title">{% icon galaxy-panelview %} 3. Plan</h5>
        <p class="card-text">
        Share the agenda.
        </p>
      </div>
    </div>
  </div>
</div>

<div class="row">
  <div class="col-sm-4">
   <div class="card bg-light mb-3" >
      <div class="card-body">
        <h5 class="card-title">{% icon version %} 4. Hook</h5>
        <p class="card-text">
        Why is the content important?
        </p>
      </div>
    </div>
  </div>
<div class="col-sm-4">
   <div class="card bg-light mb-3" >
      <div class="card-body">
        <h5 class="card-title">{% icon galaxy-wf-edit %} 5. Instruct</h5>
        <p class="card-text">
        Watch how I do it.
        </p>
      </div>
    </div>
  </div>
<div class="col-sm-4">
   <div class="card bg-light mb-3" >
      <div class="card-body">
        <h5 class="card-title">{% icon galaxy-rulebuilder-history %} 6. Practise</h5>
        <p class="card-text">
        You help me do it, I'll watch you do it.
        </p>
      </div>
    </div>
  </div>
</div>

<div class="row">
  <div class="col-sm-4">
   <div class="card bg-light mb-3" >
      <div class="card-body">
        <h5 class="card-title">{% icon galaxy-barchart %} 7. Wrap-Up</h5>
        <p class="card-text">
         Foster retention and reinforcement.
        </p>
      </div>
    </div>
  </div>
<div class="col-sm-4">
   <div class="card bg-light mb-3" >
      <div class="card-body">
        <h5 class="card-title">{% icon license %} 8. Evaluate</h5>
        <p class="card-text">
        Monitor progress.
        </p>
      </div>
    </div>
  </div>
<div class="col-sm-4">
   <div class="card bg-light mb-3" >
      <div class="card-body">
        <h5 class="card-title">{% icon galaxy-history-answer %} 9. Reflect</h5>
        <p class="card-text">
        How did it go?
        </p>
      </div>
    </div>
  </div>
</div>

<a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_Book/4%20-%20FAIR-by-design%20learning%20materials%20creation/4.1%20-%20Workflow%20stages%20description/413-design/" class="btn btn-dark text-white stretched-link">Go to the full description of FAIR-by-Design Methodology: Design stage ...</a>

## Stage 4 - Produce

> To contrive is nothing! To construct is something! To produce is everything!
{: .quote author="Edward Rickenbacker"}

### Choose Tools & Formats

> <tip-title>Collaborative environment for team work</tip-title>
>  Choose an environment for producing the learning material that will enable multiple people to work on the same material at one.
>  
>  - GitHub is one of the most popular options at the moment (find out more <a href='https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%204%20–%20Produce/08-Development%20Tools/08-Introduction%20to%20Markdown%20and%20Git/'>here</a>)
>  - GTN is also based on GitHub and defines [specific procedures on how one can contribute and collaborate](https://github.com/galaxyproject/training-material/blob/main/CONTRIBUTING.md) using GitHub
>  
>  > <details-title> GTN: collaboration using GitHub</details-title>
> > [Follow the Contributing with GitHub via its interface tutorial available on GTN]( {% link topics/contributing/tutorials/github-interface-contribution/tutorial.md %}  )
> > 
> > [Follow the Contributing with GitHub via command-line tutorial available on GTN]( {% link topics/contributing/tutorials/github-command-line-contribution/tutorial.md %}  )
> {: .details}
{: .tip}

> <tip-title>Granular versioning for easy rollback</tip-title>
>    - Versioning helps you maintain control over your changes.
>    - GitHub natively provides versioning and history retention that help easy roll back to an earlier stable state.
{: .tip}

> <tip-title>Open file formats to foster reuse</tip-title>
> - For other people to reuse the materials they should be made available using open file formats.
> - GTN promotes the use of the MD open file format for the main learning content empowered with open scientific notebooks and workflows.
> - If you use close file formats then you MUST clearly state the tools that have been used in more details.
{: .tip}

> <tip-title>Multimodal content to reach all audience</tip-title>
>  - Don't forget to include different types of multimedia to provide support for different learning modalities: read/write, auditory, visual, kinesthetic.
>  - In addition to hands-on exercises, GTN also provides support for audio and video modalities.
{: .tip}

> <tip-title>Two file sets: editable + final</tip-title>
> - The MD and supporting files are the main files used for development of the content. These files are what matters for you and other instructors.
> - Based on the editable files, GTN automatically generates the final (non-editable version) in HTML (and PDF). These  are shared with the learners.
> - GTN takes care of the revision numbering for you.
{: .tip}

> <tip-title>Don't forget to take advantage of co-creation</tip-title>
> Truly FAIR learning materials enable co-creation with external parties.
>  GitHub is a collaborative environment that supports co-creation in every step of the learning materials development and revision process.
>  The GTN fork-and-pull process of contributing to the learning material is a clear example of co-creation implementation.
{: .tip}


### Plan to reuse existing material?

#### Check the license

The existing materials you are reusing are available under a [CC license](https://creativecommons.org/share-your-work/cclicenses/), but it is different than the one you plan to use for your materials. Depending on how you want to reuse the material, you will need to consider the following aspects.

><details-title>I want to reuse it as a whole</details-title>
>- You can't use something that is licensed with ND (no derivatives).
>- In this case you must follow the rules on combining and adapting CC material.
>
>[Read more](https://creativecommons.org/faq/#combining-and-adapting-cc-material)
{: .details}

><details-title>I want to reuse a small part of it</details-title>
> No problem, you can reuse any existing CC licensed material in your learning materials as long as the reused portion is used as a showcase or to make a specific point and it is not the core of your work.
> 
> Remember that if the work is licensed with ND, you can not modify it while reusing.
{: .details}

><tip-title>Learn more about IPR</tip-title>
>[How Intellectual Property Rights (IPR) protect the interests of the creators and owners by providing them with rights over their creation?](https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%201%20%E2%80%93%20Prepare/02-Preparing%20FAIR%20Learning%20Objects/02-Preparing%20FAIR%20Learning%20Objects/#intellectual-property-rights-ipr)
{: .tip}

Still need help? [{% icon point-right %} Go to CC licensing FAQ](https://creativecommons.org/faq/#before-using-cc-licensed-material)


#### Attribute

All CC licenses require that you attribute the author, and this rule is recommended even if the license is public domain CC-0.

If the work you are reusing has a copyright notice ('© some text') you need to reproduce it while you credit the work.

You should also be able to remove attribution upon request.

##### How to attribute?

> <solution-title>Use the authors recommended attribution</solution-title>
> If the original author has provided a cite-as information, use it to attribute the work.
{: .solution}

> <solution-title>Use TASL</solution-title>
> Provide the Title, Author, Source and License of the work that you are reusing.
>
>- Source is the URL to the original work.
>- If there is a URL to the author personal pages, provide it together with the name.
>- Provide the name of the license and a URL to the license.
{: .solution}

> <details-title>Examples</details-title>
> Examples are taken from [Best Practices for Creative Commons attributions - how to attribute works you reuse under a Creative Commons license](https://www.newmediarights.org/guide/how_to/creative_commons/best_practices_creative_commons_attributions) submitted by [New Media Rights](https://www.newmediarights.org/) available under a [CC BY-NC 3.0 US DEED](https://creativecommons.org/licenses/by-nc/3.0/us/).
> - Webpage/Blog - Title (with link to original work), author (or username) (with link to author's website), and license (with link).
> 	- [Undercover Vampire Policeman](https://chriszabriskie.bandcamp.com/album/undercover-vampire-policeman) by [Chris Zabriskie](https://chriszabriskie.bandcamp.com/), available under a [Creative Commons Attribution 4.0 License](http://creativecommons.org/licenses/by/4.0/)
> - Book – Title, author, license written somewhere near the title and author if it’s a hard copy or if it’s an online book you should include a link to the licensed terms.
> 	- [From Dust to Digital: Ten Years of the Endangered Archives Programme](https://books.google.com/books?id=ImO3BgAAQBAJ&pg=PR4&dq=creative+commons+4.0&hl=en&sa=X&ved=0CEoQ6AEwCGoVChMIspCXhPPxxgIVSF0eCh27NA5X#v=onepage&q=creative%20commons%204.0&f=false) by Maja Kominko under a Creative Commons Attribution Non-commercial Non-Derivative 4.0 International license (CC BY-NC-ND 4.0)
> - Online Video - Title, author, license written into credits at end of video.  Ideally make the text clickable to the original work.  Put links to the original work and the license terms in the information section for the particular work (i.e. on the right in YouTube).
> 	- [http://www.youtube.com/watch?v=fDbbdeIXO0w#t=3m0s](http://www.youtube.com/watch?v=fDbbdeIXO0w#t=3m0s)
> - Podcast/Audio - Title, author, license read at the end of the entire work.
> 	- [“Je Suis Rick Springfield”](http://www.jonathancoulton.com/wiki/Je_Suis_Rick_Springfield) from the album [Artificial Heart](http://www.jonathancoulton.com/wiki/Artificial_Heart_(album)), by Jonathan Coulton, used under a [Creative Commons Attribution-Noncommercial 3.0 Unported License](http://creativecommons.org/licenses/by-nc/3.0/)
>- Photo/Drawing/Illustration – Title, author, license (with link online) or in close proximity to the tangible work (either in the border or directly on the work, if applicable).
> 	- [Comcast protest](http://www.flickr.com/photos/ari/8503459/in/set-214952/)” by Flikr user [Steve Rhodes](http://www.flickr.com/photos/ari/) used under [Creative Commons Attribution 2.0 license](http://creativecommons.org/licenses/by/2.0/deed.en)
{: .details}

> <details-title>Read more about how to attribute</details-title>
> [{% icon point-right %}  Attribution and Citing](https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%201%20%E2%80%93%20Prepare/02-Preparing%20FAIR%20Learning%20Objects/02-Preparing%20FAIR%20Learning%20Objects/#attribution-and-citing)
> 
> [{% icon point-right %} How to handle attribution?](https://courses.lumenlearning.com/suny-oerguide/chapter/how-to-handle-attribution/)
{: .details}


### Accessibility

The developed learning materials should cover the widest range of learner variability including the ones that use or do not use assistive technology.

![Universal access logo](../../images/universal-access-6602642_640.png)  
<a href="https://pixabay.com/vectors/universal-access-human-icon-6602642/">"Universal Access Human"</a> by <a href="https://pixabay.com/users/inspire-studio-22128832/">J S</a> from <a href="https://pixabay.com/">Pixabay</a>


#### Standards

There are several standards that govern the rules on the level of accessibility.

Most commonly used is the <a href="https://www.w3.org/WAI/standards-guidelines/wcag/">W3C Web Content Accessibility Guidelines (WCAG) standard version 2.1</a>. Three conformance levels exist, you should aim for AA which is the middle one.

PDF document accessibility is measured with a separate technical specification <a href="https://pdfua.foundation/en/why-pdf-ua/">PDF/UA (Universal Accessibility)</a>.


#### Accessibility guidelines

<a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%204%20%E2%80%93%20Produce/11-Accessibility/11-Checking_accessibility/#general-guidelines-for-development-of-accessible-materials" class="btn btn-primary stretched-link">Review some general guidelines for development of accessible material</a>

><tip-title>GTN Accessibility</tip-title>
> GTN has developed an <a href="https://training.galaxyproject.org/training-material/accessibility.html"> accessibility mission </a> aiming for WCAG 2.0 AA compliance by implementing a number of accessibility features.
> 
{: .tip}


### Internal QA

<div class="row">
  <div class="col-sm-4">
   <!-- <div class="card text-white bg-secondary mb-3" > -->
   <div class="card bg-light mb-3" >
      <!-- <div class="card-header text-white"> -->
      <div class="card-body">
        <h5 class="card-title">{% icon congratulations %} QA Self-assessment</h5>
        <p class="card-text">
         Check if everything is as it should be.
        </p>
      </div>
    </div>
  </div>
<div class="col-sm-4">
   <!-- <div class="card text-white bg-secondary mb-3" > -->
   <div class="card bg-light mb-3" >
      <!-- <div class="card-header text-white"> -->
      <div class="card-body">
        <h5 class="card-title">{% icon pref-list %} Quantitative</h5>
        <p class="card-text">
        Are all required elements produced?
        </p>
      </div>
    </div>
  </div>
<div class="col-sm-4">
   <!-- <div class="card text-white bg-secondary mb-3" > -->
   <div class="card bg-light mb-3" >
      <!-- <div class="card-header text-white"> -->
      <div class="card-body">
        <h5 class="card-title">{% icon solution %} Qualitative</h5>
        <p class="card-text">
        Do all learning units provide materials to reach the learning objectives with different modalities?
        </p>
      </div>
    </div>
  </div>
</div>

<a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_Book/4%20-%20FAIR-by-design%20learning%20materials%20creation/4.1%20-%20Workflow%20stages%20description/414-produce/" class="btn btn-dark text-white stretched-link">Go to the full FAIR-by-Design Methodology: Produce stage ...</a>

## Stage 5 - Publish

> Publishing is the art of working on a creative idea and turning it into a masterpiece!
{: .quote author="Unknown"}

><warning-title>Publishing closed FAIR materials</warning-title>
> Having FAIR learning materials does not always mean that the materials are open to everyone and there are no costs or access rules attached. In this case the bundle that is going to be published in an open repository such as Zenodo should contain the following:
>  1. Syllabus, that contains all metadata that describe the materials. Metadata should always be open.
>  2. Accompanying information (optional) to augment the description of the materials and describe the details when it comes to accessing and using the materials from a trainer perspective.
>  
> The complete learning materials package itself should be published in a closed repository where the corresponding access rules (and costs) can be implemented.
{: .warning}

### Final preparations

><tip-title>Time to provide the accompanying information</tip-title>
> Usually provided at (or near) the root of the learning materials, the purpose of the accompanying information is to further describe the content and define the rules of reuse.
{: .tip}

<div class="row">
  <div class="col-sm-4">
   <!-- <div class="card text-white bg-secondary mb-3" > -->
   <div class="card bg-light mb-3" >
      <!-- <div class="card-header text-white"> -->
      <div class="card-body">
        <h5 class="card-title">{% icon text-document %} LICENSE</h5>
        <p class="card-text">
         The license defines the potential of reuse of your learning materials together with the rules regarding their attribution by others.<br>
         CC-BY-4.0 is the recommended license, which is default in GTN. <br>
         If you choose a different license, it must be supplied in the tutorial header as metadata.
         <br>
         <a href="https://creativecommons.org/licenses" > {% icon point-right %} Browse through the available CC license types
         </a>
        </p>
      </div>
    </div>
  </div>
<div class="col-sm-4">
   <!-- <div class="card text-white bg-secondary mb-3" > -->
   <div class="card bg-light mb-3" >
      <!-- <div class="card-header text-white"> -->
      <div class="card-body">
        <h5 class="card-title">{% icon sticky-note %} README</h5>
        <p class="card-text">
        A README is a text file that introduces and explains the contents of the learning topics and materials. It usually describes the context and written in a plain text format. <br> GTN provides a README file on the level of each topic, where context information about the topic and its tutorials can be provided.
        <br>
		<a href="https://www.makeareadme.com/">
        {% icon point-right %}  Make a README
        </a>
        </p>
      </div>
    </div>
  </div>
  <div class="col-sm-4">
   <!-- <div class="card text-white bg-secondary mb-3" > -->
   <div class="card bg-light mb-3" >
      <!-- <div class="card-header text-white"> -->
      <div class="card-body">
        <h5 class="card-title">{% icon cofest %} CONTRIBUTORS</h5>
        <p class="card-text">
        GTN CONTRIBUTORS is a YAML file that lists the information about all contributors of learning materials on the platform.
        <br>You must have a GitHub user to be listed as a contributor. Additional information such as ORCID can also be provided.
        <br>
        <a href="{% link topics/contributing/tutorials/create-new-tutorial-content/tutorial.md#listing-contributors %}">
		{% icon point-right %} Add to CONTRIBUTORS.yaml
		</a>
        </p>
      </div>
    </div>
  </div>
</div>

<div class="row">
<div class="col-sm-4">
   <!-- <div class="card text-white bg-secondary mb-3" > -->
   <div class="card bg-light mb-3" >
      <!-- <div class="card-header text-white"> -->
      <div class="card-body">
        <h5 class="card-title">{% icon help %} FAQs</h5>
        <p class="card-text">
	        In addition to the learning content, the Frequently Asked Questions can further help both learners and instructors with specific information about the learning context or practicalities.
	        <br> Consider adding a tutorial-specific FAQ to the GTN tutorial, where you can answer questions about the challenges of working with the hands on activities, or provide other hints and guidelines to instructors and learners.
	        <br>
		<a href="{% link topics/contributing/tutorials/create-new-tutorial-content/tutorial.md#creating-new-faqssnippets %}">
        {% icon point-right %}  Create snippets
        </a>
        </p>
      </div>
    </div>
  </div>
  <div class="col-sm-4">
   <!-- <div class="card text-white bg-secondary mb-3" > -->
   <div class="card bg-light mb-3" >
      <!-- <div class="card-header text-white"> -->
      <div class="card-body">
        <h5 class="card-title">{% icon event-cost %} FUNDERS</h5>
        <p class="card-text">
        If the creation of the learning materials is funded by an organisation or another body, they should be attributed accordingly. <br>
        In GTN this can be accomplished by adding the funding information into the special FUNDERS.yaml file.
         <br>
         <a href="https://github.com/galaxyproject/training-material/blob/main/FUNDERS.yaml">
         {% icon point-right %} Add to the Funders file
         </a>
        </p>
      </div>
    </div>
  </div>
<div class="col-sm-4">
   <!-- <div class="card text-white bg-secondary mb-3" > -->
   <div class="card bg-light mb-3" >
      <!-- <div class="card-header text-white"> -->
      <div class="card-body">
        <h5 class="card-title">{% icon new-history %} OTHER</h5>
        <p class="card-text">
        Additional miscellaneous information should be provided if possible.
        <br> For example, GTN also supports adding organisations information.
        <br>
        <a href="https://github.com/galaxyproject/training-material/blob/main/ORGANISATIONS.yaml">
        {% icon point-right %} Add org info
        </a>
		</p>
      </div>
    </div>
  </div>
</div>

<div class="row">
<div class="col-sm-4">
   <!-- <div class="card text-white bg-secondary mb-3" > -->
   <div class="card bg-light mb-3" >
      <!-- <div class="card-header text-white"> -->
      <div class="card-body">
        <h5 class="card-title">{% icon references %} CITATION</h5>
        <p class="card-text">
	        It is best practice to provide information on how you want others to cite your learning materials when they are referenced or reused. <br> GTN does this automatically by appending the "Citing this Tutorial" section at the end of each tutorial.
	        <br>
		<a href="https://training.galaxyproject.org/training-material/faqs/gtn/gtn_citing.html">
        {% icon point-right %}  See an example citation
        </a>
        </p>
      </div>
    </div>
  </div>
  <div class="col-sm-4">
   <!-- <div class="card text-white bg-secondary mb-3" > -->
   <div class="card bg-light mb-3" >
      <!-- <div class="card-header text-white"> -->
      <div class="card-body">
        <h5 class="card-title">{% icon license %} CODE OF CONDUCT</h5>
        <p class="card-text">
        A code of conduct defines the rules for how to engage in a co-creation community.
        <br>It is based on the premise of an inclusive environment that respects all contributions.
         <br>
         <a href="https://galaxyproject.org/community/coc/">
         {% icon point-right %} Read the Galaxy Project Code of Conduct
         </a>
        </p>
      </div>
    </div>
  </div>
<div class="col-sm-4">
   <!-- <div class="card text-white bg-secondary mb-3" > -->
   <div class="card bg-light mb-3" >
      <!-- <div class="card-header text-white"> -->
      <div class="card-body">
        <h5 class="card-title">{% icon comment %} REVISIONS</h5>
        <p class="card-text">
        Information about new or updated versions of your learning materials helps others keep track of the changes more easily. <br> GTN automatically updates the Revision number of each tutorial after each release.
        <br>
        <a href="https://github.com/galaxyproject/training-material/releases">
        {% icon point-right %} Compare changes between GTN releases
        </a>
		</p>
      </div>
    </div>
  </div>
</div>

### Storing and indexing

><tip-title>Automated publishing</tip-title>
>If you are working on GTN following the provided rules and procedures, GTN will automatically publish the new tutorial once the full review and release process is completed.
>
>Just follow the [GTN guide to contributing via GitHub]( {% link topics/contributing/#st-contribute %}).
{: .tip}

><tip-title>To training catalogue</tip-title>
>You are at the point when you should also consider making a record in a relevant training catalogue.
>
>GTN automatically creates a new record in the [Elixir TeSS training catalogue](https://tess.elixir-europe.org) that is the most relevant catalogue for its community.
{: .tip}

<a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_Book/4%20-%20FAIR-by-design%20learning%20materials%20creation/4.1%20-%20Workflow%20stages%20description/415-publish/" class="btn btn-dark text-white stretched-link">Go to the full FAIR-by-Design Methodology: Publish stage ...</a>

## Stage 6 - Verify

> You will find it a very good practice always to verify your references, sir!
{: .quote author="Martin Routh"}

### External QA

> <tip-title> A fresh set of eyes</tip-title>
> - Have someone who has not participated in the development of the learning materials review the final work. This will guarantee a review free of cognitive bias.
> - GTN implements this automatically via the fork-and-pull request process. Before publication your new contribution is being reviewed by GTN peers.
{: .tip}

> <tip-title> Go through the QA checklists</tip-title>
> In Skills4EOSC T2.4 has developed a number of QA checklists that you and your external reviewer need to go through so that you can ensure high-quality learning materials {% cite sanchez_2023_8305482 %}.
{: .tip}

### FAIR or not FAIR, that is the question

#### Measure FAIRness

Use the [FAIR-by-Design methodology QA checklist](https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%206%20–%20Verify/19-Final%20QA%20check/19-finalQA/#fair-by-design-methodology-qa-checklist) to check if you have followed the most important aspects of the methodology and managed to produce FAIR learning materials.

The questions marked as essential achieve bare minimum FAIRness.

><details-title>Essential requirements</details-title>
> - **Findable**
>   - Is the complete learning resource (including instructors info) registered or indexed in at least one searchable repository? Is it in a FAIR repository?
>   - Is metadata for the resource provided in both human- and machine-readable format (e.g JSON, XMLor YAML)?
> - **Accessible**
>   - Has an accessibility checker tool been utilised to improve the accessibility of all learning resource files (PDF, HTML, video, etc.)?
>   - Are access rules (authentication & authorisation) implemented for the learning resource?
> - **Interoperable**
>   - Is the RDA minimal (or domain specific) metadata schema used for the learning material description?
>   - Is the resource available in open file formats which are tool agnostic and compatible with a wide variety of existing software?
> - **Reusable**
>   - Is there clear attribution for all reused resources with compatible licenses?
>   - Has the learning resource been made available for use by defining a permissable license or policy information that allows derivations?
{: .details}

><details-title>Optional requirements</details-title>
> - Did you follow the stages of the backward instructional design process while developing the learning resource?
> - Are controlled vocabularies (CVs) used for describing the resource characteristics aligned with the chosen metadata schema?
> - Does the learning resource represent a complete learning object defined around minimum one learning objective?
> - Does the resource incorporate an instructor kit that aids in facilitating the process of others reusing learning material by offering helpful how-to guides?
>    - facilitator guide
>    - activities description
>    - assessment activities and strategy to assess
>    - general learning content or instructor notes
>    - lesson unit plan
>    - syllabus
> - Have you employed a versioning system to track and control changes in your materials?
> - Are the resource access rules (how to access, e.g. registration procedure) explicitly communicated to learners?
> - Is the learning resource searchable in at least one relevant catalogue?
>    - Is it FAIR (can be searched based on metadata)?
> - Does the course include the possibility to provide feedback or comments from users and-or trainers/designers?
>    - If so, do you regularly gather and analyse that feedback?
> - Does the resource adopt an open community approach regarding its quality and reachability?
> - Has the learning resource been checked by a third party regarding its learning experience quality?
{: .details}

### Feedback QA

> <comment-title>Regularly gather feedback from learners and instructors</comment-title>
> Ensure that you actively and regularly gather feedback from both perspectives: the learners and the instructors.
{: .comment}

<a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_Book/4%20-%20FAIR-by-design%20learning%20materials%20creation/4.1%20-%20Workflow%20stages%20description/416-verify/" class="btn btn-dark text-white stretched-link">Go to the full FAIR-by-Design Methodology: Verify stage ...</a>

## Stage 7 - Continuous Improvement
> To improve is to change; to be perfect is to change often
{: .quote author="Winston Churchill"}

### {% icon galaxy-download %} Gather

- Gather feedback from all available internal & external sources.
- Potential sources:
	- Feedback form
	- QA recommendations
	- Self-reflection after training
	- Git Issues
	- Matrix Chat
	- Direct mail contact
	- Other means of communication

### {% icon galaxy-barchart %} Analyse

- Analyse the gathered information in a structured way.
- Create a list of potential improvements with impact level (high, moderate, low).

### {% icon galaxy-wf-edit %} Improve

- Select items from the list that will be part of a new version.
- Choose items that make sense to be in the same new release.

### {% icon galaxy-history-refresh %} Repeat

- Start a new cycle of the FAIR-by-Design methodology that will implement the selected items.
- After the Verify stage, you will reenter continuous improvement with the newly gathered information.

<a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_Book/4%20-%20FAIR-by-design%20learning%20materials%20creation/4.2%20-%20Continuous%20Improvement/417-improvement/" class="btn btn-dark text-white stretched-link">Go to the full FAIR-by-Design Methodology: Continuous Improvement ...</a>

# About

## FAIR-by-Design Microlearning

### Self-paced how to FAIR-by-Design guide

### Original Location

- online: [https://fair-by-design-methodology.github.io/microlearning/latest/](https://fair-by-design-methodology.github.io/microlearning/latest/)
- this tutorial is an adapted version of the original microlearning guide

### Training Description

Short how-to that guides you through the stages of the FAIR-by-Design methodology without any specific choice on tools and formats.

### Target audience

All interested parties who need to develop learning materials for any type of project-related training.

### Expertise Level / Skill Level: Beginner

### Primary Language: English

### Access Cost: No

### Prerequisites

No prior knowledge is required.

### Learning Objectives

- Design FAIR learning materials
- Structure FAIR learning materials
- Create and publish FAIR-by-Design learning materials
- Assess FAIR-ness of learning objects

### Keywords

FAIR, learning objects, methodology, practical implementation

### Certification Information

This training has no certification or any recognition mechanism included.

### Original Author(s)

- Sonja Filiposka, Dominique Green, Anastas Mishev, Vojdan Kjorveziroski, Andrea Corleto, Eleonora Napolitano, Gabriella Paolini, Sara di Giorgio, Joanna Janik, Luca Schirru, Arnaud Gingold, Christine Chardosek, Irakleitos Souyioultzoglou, Carolin Leister, Emma Lazzeri


### Contact information

- For more information please contact the T2.3 FAIR-by-Design Methodology Task Leader Sonja Filiposka using [sonja.filiposka@finki.ukim.mk](mailto:sonja.filiposka@finki.ukim.mk).

### Original License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

### Original microlearning unit DOI

[https://doi.org/10.5281/zenodo.11548062](https://doi.org/10.5281/zenodo.11548062)

#### Acknowledgement

These learning materials have been developed by following the [FAIR-by-Design Methodology](https://zenodo.org/records/8305540).

