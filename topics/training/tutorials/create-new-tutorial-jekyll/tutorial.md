---
layout: tutorial_hands_on
topic_name: training
tutorial_name: create-new-tutorial-jekyll
---

# Introduction
{:.no_toc}

Galaxy is a great solution to train the bioinformatics concepts:

- numerous bioinformatics tools are available (almost 5,000 in the ToolShed)
- it can be used by people without amy computer science skills
- it trains to use technology, outlining available resources and efforts that have made them accessible to researchers
- it is scalable

In 2016, the Galaxy Training Network decide to set up a new infrastructure for delivering easily Galaxy related training material. The idea was to develop something open and online based on a community effort, as always in Galaxy.

We took inspiration from [Software Carpentry](https://software-carpentry.org) and collected everything on a GitHub repository: [https://github.com/galaxyproject/training-material ](https://github.com/galaxyproject/training-material).
We decided on a structure based on tutorials with hands-on, fitting both for online self-training but also for workshops, grouped in topics. Each tutorial follows the same structure and comes with a virtualised isntance to run the training everywhere.

In this tutorial, you will learn how to run a local instance of the GTN webiste with all materials to test and develop new training sessions.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

> ### Devloping GTN training material
>
> This tutorial is part of a series to develop GTN training material, feel free to also look at:
>
> 1. [Writing content in markdown](../create-new-tutorial-content/tutorial.html)
> 1. [Defining metadata](../create-new-tutorial-metadata/tutorial.html)
> 1. [Setting up the infrastructure](../create-new-tutorial-jekyll/tutorial.html)
> 1. [Creating Interactive Galaxy Tours](../create-new-tutorial-tours/tutorial.html)
> 1. [Building a Docker flavor](../create-new-tutorial-docker/tutorial.html)
> 1. [Submitting the new tutorial to the GitHub repository](../../../dev/tutorials/github-contribution/slides.html)
{: .agenda}

# Checking the website generation

If you want to run the entire GTN website locally or test your new training material we can do this! Currently, the website is generated from the metadata and the tutorials using Jekyll, a simple static site builder.
We can use Jekyll to run a server to check if the tutorial is correctly added and rendered.

> ### :pencil2: Hands-on: Checking the website generation locally
>
> 1. Install Jekyll using [RubyGems](https://rubygems.org/pages/download): `make install`
>
> 2. Run a local Jekyll server: `make serve`
> 3. Visualize at [http://localhost:4000/ ](http://localhost:4000/)
>
>    > ### :question: Questions
>    >
>    > How to check if the server was started and if all topics are included?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    Please check [http://localhost:4000/topics/ ](http://localhost:4000/topics/) to get a list of topics.
>    >    </details>
>    {: .question}
{: .hands_on}

# Conclusion
{:.no_toc}

> ### Developing GTN training material
>
> This tutorial is part of a series to develop GTN training material, feel free to also look at:
>
> 1. [Writing content in markdown](../create-new-tutorial-content/tutorial.html)
> 1. [Defining metadata](../create-new-tutorial-metadata/tutorial.html)
> 1. [Setting up the infrastructure](../create-new-tutorial-jekyll/tutorial.html)
> 1. [Creating Interactive Galaxy Tours](../create-new-tutorial-tours/tutorial.html)
> 1. [Building a Docker flavor](../create-new-tutorial-docker/tutorial.html)
> 1. [Submitting the new tutorial to the GitHub repository](../../../dev/tutorials/github-contribution/slides.html)
{: .agenda}