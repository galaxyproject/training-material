---
layout: tutorial_hands_on
title: REMBI: Recommended Metadata for Biological Images – metadata guidelines for bioimaging data

zenodo_link: ''

questions:

objectives:
- Organise bioimage metadata
- Find out what REMBI is and why it is useful
- Categorise what metadata belongs to each of the submodules of REMBI
- Gather the metadata for an example bioimage dataset
  
time_estimation: "15min"

key_points:

tags:
- fair
- data management
- bioimaging
  
priority: 5

contributions:
  authorship:
    - wee-snufkin
    - kkamieniecka
    - poterlowicz-lab

subtopic: fair-data

requirements:
  - type: "internal"
    topic_name: fair
    tutorials:
      - fair-intro
      - data-management
      - bioimage-metadata
---

# Metadata guidelines for bioimaging data

REMBI (Recommended Metadata for Biological Images) was proposed as a draft metadata guidelines to begin addressing the needs of diverse communities within light and electron microscopy. Currently, these guidelines are in draft form to encourage discussion within the community, but they provide a useful guide as to what metadata should be gathered to make your image data FAIR. They divide the metadata requirements into eight submodules, which breaks up what can seem to be a daunting task! To find out more, check the [publication](https://www.nature.com/articles/s41592-021-01166-8).


> <question-title></question-title>
>
> In the [REMBI paper](https://www.nature.com/articles/s41592-021-01166-8), the authors consider three potential user groups who require different metadata. Find out what these three groups are and what their metadata requirements are.
>
> > <solution-title></solution-title>
> > The identified three user groups are: Biologists, Imaging scientists, Computer-vision researchers. 
> > - A research biologist may be interested in the biological sample that has been imaged to compare it to similar samples that they are working with.
> > - An imaging scientist may be interested in how the image was acquired so they can improve upon current image acquisition techniques.
> > - A computer vision researcher may be interested in annotated ground-truth segmentations, that can be obtained from the image, so they can develop faster and more accurate  algorithms.
> {: .solution}
>
{: .question}

> <tip-title>Instructor Note</tip-title>
>
> If you're an instructor leading this training, you might ask people to work in small groups for this exercise and encourage discussion within them. Ask group members to share with their group which of the user groups they identify as and what metadata they would want.
> 
{: .tip}

# Categories of metadata 
REMBI covers different categories of metadata, such as: study, study component, biosample, specimen, image acquisition, image data, image correlation, analyzed data. Within each module, there are attributes that should be included to make the published data FAIR. We will explore all the modules and attributes suggested by REMBI and we'll show some examples as well. 

## Study 
This first module of REMBI metadata describes the Study and should include:
- Study type
- Study description
- General dataset information

### Study type
Ideally, the study type will be part of an ontology. You can look up the main subject of your study using a tool like [OLS](https://www.ebi.ac.uk/ols/index) to find a suitable ontology. This will help others to see where your study sits within the wider research area.

> <comment-title>Example</comment-title>

{: .comment}



