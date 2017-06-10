# Introduction

Galaxy is a great solution to train the bioinformatics concepts:

- numerous bioinformatics tools are available (almost 5,000 in the ToolShed)
- it can be used by people without amy computer science skills
- it trains to use technology, outlining available resources and efforts that have made them accessible to researchers
- it is scalable

In 2016, the Galaxy Training Network decide to set up a new infrastructure for delivering easily Galaxy related training material. The idea was to develop something open and online based on a community effort, as most of the time in Galaxy. 

We take inspiration from [Software Carpentry](https://software-carpentry.org). We collected everything on a GitHub repository: [https://github.com/galaxyproject/training-material](https://github.com/galaxyproject/training-material). We decided a structure based on tutorials with hands-on, fitting both for online self-training but also for workshops, grouped in topics. Each tutorial follows the same structure and comes with a technical support to be able to run. 

In this tutorial, you will understand how to design and develop a new tutorial fitting in this training material repository. As doing helps to understand, we will develop a small tutorial to explain BLAST with the full infrastructure to be able to run this tutorial anywhere. 

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. [Setting up a new tutorial](#setting-up-a-new-tutorial)
> 2. [Filling the metadata](#filling-the-metadata)
> 3. [Analysis of the differential expression](#analysis-of-the-differential-expression)
> {: .agenda}

# Setting up a new tutorial

Here, we want to develop a small tutorial to explain how to use BLAST. The first step we need to define is in which topic putting our tutorial. This first step can be tricky. 

When we structured the repository, we decided here to use as topic the names of the categories in the [ToolShed](https://toolshed.g2.bx.psu.edu/). So when decided where to put your tutorial, you can look in which ToolShed's category are the main tools used in the tutorial and use this category as topic. For example, this tutorial will rely on the NCBI Blast+ tool.

> ### :pencil2: Hands-on: Defining the topic for the tutorial
>
> 1. Search for NCBI Blast+ on the [ToolShed](https://toolshed.g2.bx.psu.edu/)
> 2. Check in which category it is
>
>    > ### :question: Questions
>    >
>    > In which topic will you put the tutorial?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    If we search for [NCBI Blast+ in the ToolShed](https://toolshed.g2.bx.psu.edu/view/devteam/ncbi_blast_plus/7538e2bfcd41), it is attributed to 2 categories (bottom): "Next Gen Mappers" and "Sequence Analysis".
>    >    We decided to put it in "Sequence analysis" because this is the most general one for this tutorial.
>    >    </details>
>    {: .question}
{: .hands_on}

Once the topic is chosen, serious things can start: creating the tutorial. It is meaning the tutorial content, the metadata related to the tutorial but also the technical support for the tutorial with the description of the needed tool and dataset, a workflow of the tutorial and also a Galaxy Interactive Tour.

To help you, we created a template for a tutorial with the different required files.

> ### :pencil2: Hands-on: Copy the needed file
>
> 1. Copy the `tutorial1` folder (you can find it in `templates/tutorials/`) in `topics/sequence-analysis/topics`
> 2. Rename the folder into `similarity-search`
{: .hands_on}

We will now start to fill the different files together.

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

> ### :pencil2: Hands-on: Fill the basic metadata
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

# Filling the tutorial

# Setting up the technical support


# Filling the tutorial

## Finding a good toy dataset

## Filling the tutorial content

Short introduction about this subpart.

> ### :pencil2: Hands-on: Data upload
>
> 1. Step1
> 2. Step2
>
>    > ### :question: Question
>    >
>    > Question?
>    >
>    > <details>
>    > <summary>Click to view answers</summary>
>    > Answer to question
>    > </details>
>    {: .question}
{: .hands_on}

Some blabla
> ### :pencil2: Hands-on: Data upload
>
> 1. Step1
> 2. Step2
>
>    > ### :question: Questions
>    >
>    > 1. Question1?
>    > 2. Question2?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>Answer for question1</li>
>    >    <li>Answer for question2</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 3. Step3
{: .hands_on}

## Testing the tutorial

## Extracting the workflow

## Creating a tour

# Adding slides

(If needed)

# Building the infrastructure to run the tutorial

## Filling the `tools.yaml`

## Filling the `data-library.yaml`

## Filling the `data-manager.yaml`

> ### :nut_and_bolt: Comment
>
> Do you want to learn more about the principles behind mapping? Follow our [training](../../NGS-mapping)
> {: .comment}

# Submitting the new tutorial

# Conclusion

Conclusion about the technical key points. And then relation between the technics and the biological question to end with a global view.
