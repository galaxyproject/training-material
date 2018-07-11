---
layout: tutorial_hands_on
topic_name: contributing
tutorial_name: create-new-tutorial-metadata
---

The `metadata.yaml` file is a file located in the topic directory. It describes the metadata related to a topic and its material.

> ### Agenda
>
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Make topic aware of a tutorial

To make a tutorial appear within a topic, we need to add metadata for the tutorial in the topic's `metadata.yaml` file:

- `title`: title of the tutorial (it will appear on the tutorial page and the topic page)
- `type: "tutorial"`
- `name`: name of the tutorial (name of the subdirectory where the files related to the tutorial will be stored)
- `enable`: `false` to hide your tutorial from the topic page
- `hands_on`(`yes` or `no`): says if the tutorial includes hands-on sections in the `tutorial.md`
- `slides` (`yes` or `no`): says if slides are available for this material

> ### {% icon hands_on %} Hands-on: Fill the basic metadata
>
> 1. Open the `metadata.yaml` of `sequence-analysis` topic (`/topics/sequence-analysis/metadata.yaml`)
> 2. Update the tutorial information in the `material` section:
>
>     ```
>     -
>       title: "Similarity search with BLAST"
>       type: "tutorial"
>       name: "similarity-search"
>       hands_on: yes
>       slides: no
>     ```
>
{: .hands_on}

This information is used to automatically make the tutorial available on the website: [{{site.url}}{{ site.baseurl}} ]({{site.url}}{{ site.baseurl}})

# Technical support

After the definitining the tutorial in the topic metadata, we will also add some information about the technical support available for the tutorial:

- `zenodo_link`: link on Zenodo to the input data for the tutorial
- `workflows` (`yes` or `no`): says if a workflow is available for this material (in the `workflows` subdirectory)
- `galaxy_tour`(`yes` or `no`): says if an interactive tour is available for the tutorial (in the `tours` subdirectory)


> ### {% icon hands_on %} Hands-on: Fill the technical metadata
>
> 2. Add in the `material` section after `slides: no` the new tutorial:
>
>     ```
>       zenodo_link: ""
>       workflows: no
>       galaxy_tour: no
>     ```
>
{: .hands_on}

This information is used to display the data, workflow and tour links from the topic and tutorial page. They are also used to check which information are missing for the tutorials.

# Pedagogical support

We also define metadata related to the pedagogical content of the tutorial, which will appear in the top ("Overview" box) and bottom of the online tutorial:

- `requirements`: list of resources that the reader of the material should be familiar with before starting this training:
    - `title`
    - `link`: relative for internal (inside training material) requirement or full for external requirement)
    - `type`: the type of link (`internal` or `external`)
- `time_estimation`: an estimation of the time needed to complete the hands-on
- `questions`: list of questions that will be addressed in the tutorial
- `objectives`: list of learning objectives for the tutorial

    A learning objective is a single sentence describing what a learner will be able to do once they have done the tutorial

- `key_points`: list of take-home messages

    This information will appear at the end of the tutorial

For this category of metadata, we have taken inspiration from what Software Carpentry has done and particularly what they described in their [Instructor training](https://swcarpentry.github.io/instructor-training/).

> ### {% icon hands_on %} Hands-on: Fill out the pedagogical metadata
>
> 1. Define 2 questions that will be addressed during the tutorial and add them to the metadata
> 2. Define 2 learning objectives for the tutorial and add them to the metadata
{: .hands_on}

> ### {% icon comment %} When filling the pedagogical metadata
> We recommend that you fill out the *questions* and the *learning objectives* before starting writing the tutorial content. You can still refine them afterwards, but it will help to guide you in developing your tutorial, and gives you some time to think beforehand on what topics are worth being covered.
>
> For the take-home messages, it is easier to define them once the tutorial is written and you identified the issues.
{: .comment}

# Conclusion
{:.no_toc}
