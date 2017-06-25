---
layout: tutorial_hands_on
topic_name: training
tutorial_name: create-new-tutorial-metadata
---

# Introduction

Galaxy is a great solution to train the bioinformatics concepts:

- numerous bioinformatics tools are available (almost 5,000 in the ToolShed)
- it can be used by people without amy computer science skills
- it trains to use technology, outlining available resources and efforts that have made them accessible to researchers
- it is scalable

In 2016, the Galaxy Training Network decide to set up a new infrastructure for delivering easily Galaxy related training material. The idea was to develop something open and online based on a community effort, as always in Galaxy. 

We took inspiration from [Software Carpentry](https://software-carpentry.org) and collected everything on a GitHub repository: [https://github.com/galaxyproject/training-material](https://github.com/galaxyproject/training-material).
We decided on a structure based on tutorials with hands-on, fitting both for online self-training but also for workshops, grouped in topics. Each tutorial follows the same structure and comes with a virtualised isntance to run the training everywhere.

In this tutorial, you will learn how to annotate your training material with a lot of metadata, so that it can be reused and empower other services.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. [Setting up a new tutorial](#setting-up-a-new-tutorial)
> 2. [Filling the metadata](#filling-the-metadata)
> 3. [Checking the tutorial generation](#checking-the-tutorial-generation)
> 5. [Filling the tutorial content](#filling-the-tutorial-content)
> 6. [Adding slides](#adding-slides)
> 7. [Setting up the technical support](#setting-up-the-technical-support)
> 8. [Submitting the new tutorial to the GitHub repository](#submitting-the-new-tutorial-to-the-github-repository)
> {: .agenda}


# Filling the metadata

The first file we will fill is the `metadata.yaml` file. 

This file define the metadata related to a tutorial: 

- `title`: title of the tutorial
- `type: "tutorial"`
- `name`: name of the tutorial (name of the subdirectory where the files related to the tutorial will be stored)
- `zenodo_link`: link on Zenodo to the input data for the tutorial (not ideal but it can be empty)
- `galaxy_tour`: name of the galaxy tour
- `hands_on`(`"yes"` or `"no"`): tell if an hands on is available for this material
- `slides` (`"yes"` or `"no"`): tell if slides are available for this materialits title, its type, ...
- `requirements`: list of requirements specific to this tutorial (in addition to the one of the topic), with:
    - `title`
    - `link`: relative for internal (inside training material) requirement or full for external requirement)
    - `type`: the type of link (`internal` or `external`)

This information is used to automatically make the tutorial available on the online website: [http://galaxyproject.github.io/training-material/](http://galaxyproject.github.io/training-material/)

> ### :pencil2: Hands-on: Fill the basic metadata
>
> 1. Fill the basic metadata for our tutorial
>   - `title: Similarity search with BLAST`
>   - `type: "tutorial"`
>   - `name: "similarity-search"`
>   - `zenodo_link: ""` (we do not have data currently)
>   - `galaxy_tour: ""` (we do not have Galaxy Interactive Tour currently)
>   - `hands_on: "yes"`
>   - `slides: "no"`
>   - `requirements`: a requirement to "Galaxy introduction" with internal link to `introduction`
{: .hands_on}

In the second part of the metadata, we define metadata related to the content of the tutorial, which will appear in the top and bottom of the online tutorial:

- `time_estimation`: an estimation of the time needed to complete the hands-on
- `questions`: list of questions that will be addressed in the tutorial
- `objectives`: list of learning objectives of the tutorial

    A learning objective is a single sentence describing what a learner will be able to do once they have deone the tutorial

- `key_points`: list of take-home messages

    This information will appear at the end of the tutorial 

For this metadata, we take inspiration from what Software Carpentry is doing and particularly what they describe in their [Instructor training](http://swcarpentry.github.io/instructor-training/) and the section ["Lessons and Objectives"](http://swcarpentry.github.io/instructor-training/19-lessons/). 

> ### :pencil2: Hands-on: Fill the pedagogical metadata
>
> 1. Define 2 questions that will be addressed during the tutorial and add them to the metadata
> 2. Define 2 learning objectives for the tutorial and add them to the metadata
{: .hands_on}

We recommend you to fill the questions and the learning objectives before starting writing the tutorial content. You can still refine them afterwards, but it will help to design your tutorial and think beforehands what is worth training.

For the take-home messages, it is easier to define them once the tutorial is written and you identified the issues.

> ### :nut_and_bolt: Comment
>
> Temporarly, we need a duplication of the metadata (because of our temporary templating system) in the topic `metadata.yaml` in the section `material`. So you need to copy that.
> {: .comment}

