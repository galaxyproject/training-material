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


{% include _includes/tab-choices.html option1="Descriptive info fields" option2="Access info fields" option3="Educational info fields" default="Descriptive info fields" title="Metadata fields categories" %} 

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
> >Objectives should be specific, ​measurable, ​attainable, ​relevant and ​time-bound​.
> {: .tip}
>
> ><tip-title>Use Blooms Taxonomy</tip-title>
> >Formulate the objectives as actionable verb + observable knowledge, skill, attitude, behavior or ability.
> >
> >[See more](https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%201%20%E2%80%93%20Prepare/02-Preparing%20FAIR%20Learning%20Objects/02-Preparing%20FAIR%20Learning%20Objects_cont/#defining-learning-objectives)
> {: .tip}
{: .question }

<a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%201%20%E2%80%93%20Prepare/01-FAIR%20skills%20%26%20principles/01-FAIR%20skills%20%26%20principles/" class="btn btn-dark text-white btn-lg btn-block">Start an in-depth training on the Prepare stage....</a>
<a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_Book/4%20-%20FAIR-by-design%20learning%20materials%20creation/4.1%20-%20Workflow%20stages%20description/411-prepare/" class="btn btn-dark text-white btn-lg btn-block">FAIR-by-Design Methodology: Prepare stage....</a>

## Stage 2 - Discover
> Greater even than the greatest discovery is to keep open the way to future discovery.
{: .quote author="John Jacob Abel"}
### Get inspired
Reusable materials can be found anywhere. These are just some examples:
> <code-in-title>OER</code-in-title>
> [DOAB](https://directory.doabooks.org/)
> [MERLOT](https://www.merlot.org/merlot/index.htm)
> [OASIS](https://oasis.geneseo.edu/index.php)
> [OER Commons](https://www.oercommons.org/)
> [OERTX CORA](https://www.projectcora.org/)
> [GALILEO](https://oer.galileo.usg.edu/)
> [FORRT](https://forrt.org/)
{: .code-in}

> <code-in-title>EOSC</code-in-title>
> [EOSC Training catalogue on the EOSC Marketplace](https://search.marketplace.eosc-portal.eu/search/training?q=*)
> 
> Most EOSC projects have their own training catalogues and/or platforms...
{: .code-in}

> <code-in-title>General</code-in-title>
> [Creative Commons Search](https://search.creativecommons.org/) - content provided under a CC license
> [Zenodo](https://zenodo.org/) - a multi-disciplinary open repository
> [OSF](https://osf.io/) - a free, open research platform
{: .code-in}

### Potential for reuse
> <warning-title> Respect the licenses, to be respected!</warning-title>
> Materials with non-permissible licenses can be used for inspiration only. Materials with permissible licenses should be reused based on the license rules.
{: .warning}

### Don't forget the multimedia search
Different learners have different learning modalities (read/write, auditory, visual, kinesthetic). You should use all types of multimedia in your learning materials.
<a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%202%20%E2%80%93%20Discover/03-Existing%20learning%20materials/03-Existing%20learning%20materials/" class="btn btn-dark text-white btn-lg btn-block">Start an in-depth training on the Discover stage....</a>
<a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_Book/4%20-%20FAIR-by-design%20learning%20materials%20creation/4.1%20-%20Workflow%20stages%20description/412-discover/" class="btn btn-dark text-white btn-lg btn-block">FAIR-by-Design Methodology: Discover stage....</a>

## Stage 3 - Design
> Design is intelligence made visible.
{: .quote author="Alina Wheeler"}
### Time to brainstorm

### Structure is everything

### How to develop the learning content

## Stage 4 - Produce
> To contrive is nothing! To construct is something! To produce is everything! 
{: .quote author="Edward Rickenbacker​"}

### Choose Tools & Formats

### Plan to reuse existing material? 

### Accessibility

### Internal QA

<div class="row">
  <div class="col-sm-4">
   <!-- <div class="card text-white bg-secondary mb-3" > -->
   <div class="card bg-light mb-3" >
      <!-- <div class="card-header text-white"> -->
      <div class="card-body">
        <h5 class="card-title">{% icon congratulations %} QA Self-assessment</h5>
        <p class="card-text">
         to check if everything is as it should be. 
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
        are all required elements produced
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
        do all learning units provide materials to reach the learning objectives with different modalities. 
        </p>
      </div>
    </div>
  </div>
</div>


​<a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%204%20%E2%80%93%20Produce/08-Development%20Tools/08-Introduction%20to%20Markdown%20and%20Git/" class="btn btn-dark text-white btn-lg btn-block">Start an in-depth training on the Produce stage....</a>

<a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_Book/4%20-%20FAIR-by-design%20learning%20materials%20creation/4.1%20-%20Workflow%20stages%20description/414-produce/" class="btn btn-dark text-white btn-lg btn-block">FAIR-by-Design Methodology: Produce stage....</a>
## Stage 5 - Publish
> Publishing is the art of working on a creative idea and turning it into a masterpiece​!
{: .quote author="Unknown"}

><warning-title>Publishing closed FAIR materials</warning-title>
> Having FAIR learning materials does not always mean that the materials are open to everyone and there are no costs or access rules attached. In this case the bundle that is going to be published in an open repository such as Zenodo should contain the following:
>  1. Syllabus - that contains all metadata that describe the materials and metadata should always be open
>  2. Accompanying files - optional - augment the description of the materials and describe the details when it comes to accessing and using the materials from a trainer perspective 
> Another alternative is to publish the materials in a closed repository where the corresponding access rules can be implemented.
{: .warning}

### Final preparations
><tip-title>Time to create the accompanying files</tip-title>
>These are <a href='https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%205%20–%20Publish/16-Publishing%20Preparations/16-Publishing%20Preparations/'>files</a> that are provided in the root of your learning materials and describe the content and define the rules of reuse.
{: .tip}

<div class="row">
  <div class="col-sm-4">
   <!-- <div class="card text-white bg-secondary mb-3" > -->
   <div class="card bg-light mb-3" >
      <!-- <div class="card-header text-white"> -->
      <div class="card-body">
        <h5 class="card-title">{% icon text-document %} LICENSE</h5>
        <p class="card-text">
         Plaintext file that defines the license of your learning materials. Just copy paste it from the official CC website. CC-BY-4.0 is the recommended license. 
         <br>
         <a href="https://creativecommons.org/licenses/by/4.0/legalcode.txt" > {% icon point-right %} CC license
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
        A README is a text file that introduces and explains the contents of your learning materials. It usually describes the context and defines how the materials may be reused or co-created. It is usually written in a plain text format.
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
        <h5 class="card-title">{% icon references %} CITATION.cff</h5>
        <p class="card-text">
        Citation files are plain text files with human- and machine-readable citation information that tells others how to cite or attribute your work. 
        <br>
        <a href="https://citation-file-format.github.io/">
         {% icon point-right %} Create a citation file
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
        <h5 class="card-title">{% icon license %} CODE_OF_CONDUCT</h5>
        <p class="card-text">
        A code of conduct defines the rules for how to engage in a co-creation community. It is based on a premise of an inclusive environment that respects all contributions.
         <br>
         <a href="https://github.com/probot/template/blob/master/CODE_OF_CONDUCT.md">
         {% icon point-right %} Code of conduct template
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
        <h5 class="card-title">{% icon comment %} RELEASE NOTES</h5>
        <p class="card-text">
        A release note is a report published alongside new or updated version of your learning materials that details the changes in the new version. 
        <br>
        <a href="https://slite.com/templates/release-notes">
        {% icon point-right %} Creating Release Notes
        </a>
		</p>
      </div>
    </div>
  </div>
</div>

### Store in a repository

![Zenodo logo](../../images/Zenodo_logo.png)
<small>
<br>
<a href="https://upload.wikimedia.org/wikipedia/commons/5/58/Zenodo_logo.png">Zenodo logo</a> by <a href="https://twitter.com">a Twitter user</a> from <a href="https://commons.m.wikimedia.org/wiki/File:Zenodo_logo.png">Wikimedia</a> licensed under the <a href="https://creativecommons.org/licenses/by-sa/4.0/deed.en">Creative Commons Attribution-Share Alike 4.0 International license</a>
</small>
#### To Zenodo
Deposit your editable learning materials set to make available to other designers and instructors.
1. Create an archive of all the files in your logical hierarchical structure.
2. Create a new Zenodo record with the archive.
3. Provide a rich metadata description and link to any related resources.

<a href="https://help.zenodo.org/docs/deposit/create-new-upload/" class="btn btn-primary stretched-link">How to deposit in Zenodo</a>

><tip-title>Automated publishing to Zenodo</tip-title>
>If you are working on GitHub using the provided templates repository then the "publish to Zenodo" step is fully automated for you. Just follow the [guide to publishing](https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%205%20%E2%80%93%20Publish/17-Zenodo%20Publishing/17-Zenodo%20Publishing/).
{: .tip}

><tip-title>To training catalogue</tip-title>
>You are at the point when you should also consider making a record in a relevant training catalogue such as the EOSC training catalogue.
{: .tip}

### Provide to learners

![Moodle logo](../../images/Moodle-logo.svg.png)
<small>
<br>
<a href="https://upload.wikimedia.org/wikipedia/commons/thumb/c/c6/Moodle-logo.svg/320px-Moodle-logo.svg.png">Moodle logo</a> by <a href="https://moodle.org/">Moodle.org</a> from <a href="https://en.m.wikipedia.org/wiki/File:Moodle-logo.svg">Wikipedia</a> licensed under the <a href="https://en.wikipedia.org/wiki/en:GNU_General_Public_License">GNU General Public License</a>
</small>
#### To LMS
Generate the final versions from your editable content and add it to a course on the Skills4EOSC Learning Platform to make it available for learners.

0. Provide the course metadata.

1. Add the learning content in general non-editable file formats.

2. Add assessments such as quizzes or assignments.

3. Setup feedback gathering.

4. Define recognition mechanism such as open digital badges for successful completion.

<a href="https://docs.moodle.org/403/en/Table_of_Contents#Managing_a_Moodle_course" class="btn btn-primary stretched-link">Managing a Moodle course</a>


​​<a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%205%20%E2%80%93%20Publish/16-Publishing%20Preparations/16-Publishing%20Preparations/" class="btn btn-dark text-white btn-lg btn-block">Start an in-depth training on the Publish stage....</a>

<a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_Book/4%20-%20FAIR-by-design%20learning%20materials%20creation/4.1%20-%20Workflow%20stages%20description/415-publish/" class="btn btn-dark text-white btn-lg btn-block">FAIR-by-Design Methodology: Publish stage....</a>

## Stage 6 - Verify
> You will find it a very good practice always to verify your references, sir!
{: .quote author="Martin Routh"}
### External QA
> <tip-title> A fresh set of eyes</tip-title>
> Have someone who has not participated in the development of the learning materials review the final work. This will guarantee a review free of cognitive bias.
{: .tip}

> <tip-title> Don't forget to QA the LMS</tip-title>
>  The reviewer should play the role of a new learner in the LMS and check everything from the learner perspective.
{: .tip}

> <tip-title> Go through the QA checklists</tip-title>
> In Skills4EOSC T2.4 has developed a number of QA checklists that you and your external reviewer need to go through so that you can ensure high-quality learning materials.
{: .tip}

### FAIR or not FAIR, that is the question...

{% include _includes/tab-choices.html option1="Essential requirements" option2="Optional requirements"  default="Essential requirements" title="Measure FAIRness"  disambiguation="second" text="Use the FAIR-by-Design methodology QA checklist to check if you have followed the most important aspects of the methodology and managed to produce FAIR learning materials.
<br>
The questions marked as essential achieve bare minimum FAIRness." %} 

<div class="Essential-requirements" markdown="1">
- **Findable** =	Is the complete learning resource (including instructors info) registered or indexed in at least one searchable repository? Is it in a FAIR repository? Is metadata for the resource provided in both human- and machine-readable format (e.g JSON, XMLor YAML?
- **Accessible** =	Has an accessibility checker tool been utilised to improve the accessibility of all learning resource files (PDF, HTML, video, etc.)? Are access rules (authentication & authorisation) implemented for the learning resource?
- **Interoperable** = Is the RDA minimal (or domain specific) metadata schema used for the learning material description?</br> Is the resource available in open file formats which are tool agnostic and compatible with a wide variety of existing software?
- **Resuable** = Is there clear attribution for all reused resources with compatible licenses? Has the learning resource been made available for use by defining a permissable license or policy information that allows derivations?
</div>
<div class="Optional-requirements" markdown="1">
-  Did you follow the stages of the backward instructional design process while developing the learning resource?
-  Are controlled vocabularies (CVs) used for describing the resource characteristics aligned with the chosen metadata schema?
-  Does the learning resource represent a complete learning object defined around minimum one learning objective?
-  Does the resource incorporate an instructor kit that aids in facilitating the process of others reusing learning material by offering helpful how-to guides?
    -  facilitator guide
    -  activities description
    -  assessment activities and strategy to assess
    -  general learning content or instructor notes
    -  lesson unit plan
    -  syllabus
-  Have you employed a versioning system to track and control changes in your materials?
-  Are the resource access rules (how to access, e.g. registration procedure) explicitly communicated to learners?
-  Is the learning resource searchable in at least one relevant catalogue? 
    -  Is it FAIR (can be searched based on metadata)?
-  Does the course include the possibility to provide feedback or comments from users and-or trainers/designers? 
    -  If so, do you regularly gather and analyse that feedback?
-  Does the resource adopt an open community approach regarding its quality and reachability?
-  Has the learning resource been checked by a third party regarding its learning experience quality?
</div>

### Feedback QA 

> <comment-title>Regularly gather feedback from learners and instructors</comment-title>
> Ensure that you actively and regularly gather feedback from both perspectives: the learners and the instructors.
{: .comment}

<a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%206%20%E2%80%93%20Verify/19-Final%20QA%20check/19-finalQA/" class="btn btn-dark text-white btn-lg btn-block">Start an in-depth training on the Verify stage....</a>

<a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_Book/4%20-%20FAIR-by-design%20learning%20materials%20creation/4.1%20-%20Workflow%20stages%20description/416-verify/" class="btn btn-dark text-white btn-lg btn-block">FAIR-by-Design Methodology: Verify stage....</a>

## Stage 7 - Continuous Improvement
> To improve is to change; to be perfect is to change often
{: .quote author="Winston Churchill"}

#### {% icon galaxy-download %} Gather
Gather feedback from all available internal & external sources.

Potential sources:  
- Feedback form   
- QA recommendations   
- self-reflection after training   
- Survey   
- Direct mail contact   
- Other means of communication

<div class="row">
  <div class="col-sm-4">
   <!-- <div class="card text-white bg-secondary mb-3" > -->
   <div class="card bg-light mb-3" >
      <!-- <div class="card-header text-white"> -->
      <div class="card-body">
        <h5 class="card-title">{% icon galaxy-download %} Gather</h5>
        <p class="card-text">Potential sources:  </p>
        <ul>
		<li> Feedback form   </li>
		<li>QA recommendations   </li>
		<li> self-reflection after training </li>   
		<li> Survey </li>
		<li> Direct mail contact </li>
		<li>Other means of communication </li>
		</ul>
      </div>
    </div>
  </div>

<div class="col-sm-4">
   <!-- <div class="card text-white bg-secondary mb-3" > -->
   <div class="card bg-light mb-3" >
      <!-- <div class="card-header text-white"> -->
      <div class="card-body">
        <h5 class="card-title">{% icon galaxy-barchart %} Analyse</h5>
        <p class="card-text">Analyse the gathered information in a structured way. <br>
        Create a list of potential improvements with impact level (high, moderate, low).
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
        <h5 class="card-title">{% icon galaxy-wf-edit %} Improve</h5>
        <p class="card-text">Select items from the list that will be part of a new version. <br> Choose items that make sense to be in the same new release  </p>
      </div>
    </div>
  </div>

<div class="col-sm-4">
   <!-- <div class="card text-white bg-secondary mb-3" > -->
   <div class="card bg-light mb-3" >
      <!-- <div class="card-header text-white"> -->
      <div class="card-body">
        <h5 class="card-title">{% icon galaxy-history-refresh %} Repeat</h5>
        <p class="card-text">Start a new cycle of the FAIR-by-Design methodology that will implement the selected items. <br>After the Verify stage, you will reenter continuous improvement with newly gathered information....  </p>
      </div>
    </div>
  </div>
</div>

​<a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_ToT/latest/Stage%206%20%E2%80%93%20Verify/20-Continuous%20Improvement/20-CI/" class="btn btn-dark text-white btn-lg btn-block">Start an in-depth training on the Continuous Improvement stage....</a>

<a href="https://fair-by-design-methodology.github.io/FAIR-by-Design_Book/4%20-%20FAIR-by-design%20learning%20materials%20creation/4.2%20-%20Continuous%20Improvement/417-improvement/" class="btn btn-dark text-white btn-lg btn-block">FAIR-by-Design Methodology: Continuous Improvement stage....</a>

# About

