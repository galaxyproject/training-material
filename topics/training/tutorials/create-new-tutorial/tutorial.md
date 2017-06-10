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
> 2. [Mapping](#mapping)
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


# Copy the template

Introduction about this part

# Filling the metadata 

## Defining the learning objectives, the questions, the take-home messages

Link to Software Carpentry?

> ### :pencil2: Hands-on: Data upload
>
> 1. Step1
> 2. Step2
>
>    > ### :nut_and_bolt: Comments
>    > A comment
>    {: .comment}
>
>    > ### :bulb: Tip: A tip
>    >
>    > * Step1
>    > * Step2
>    {: .tip}
{: .hands_on}

## Decide the time for the tutorial

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
