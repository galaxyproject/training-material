---
layout: tutorial_hands_on
topic_name: contributing
tutorial_name: create-new-topic
---

# Introduction
{:.no_toc}

Each training material is related to a topic. All training materials (slides, tutorials, ...) related to a topic are found in a dedicated directory (*e.g.* `transcriptomics` directory contains the material related to exome sequencing analysis). Each topic have the following structure:

```
├── README.md
├── metadata.yaml
├── images
├── docker
│   ├── Dockerfile
├── slides
│   ├── index.html
├── tutorials
│   ├── tutorial1
│   │   ├── tutorial.md
│   │   ├── slides.html
│   │   ├── tools.yaml
│   │   ├── data-library.yaml
│   │   ├── workflows
│   │   │   ├── workflow.ga
│   │   ├── tours
│   │   │   ├── tour.yaml
```

## `images` directory
{:.no_toc}

The `images` directory collects all images/pictures needed for the training materials related to the topic, *i.e* pictures for the slides or the tutorials.

Images shared between several topics are in the `shared/images` directory at the root.

All images for the slides must be in `images` directory. The images must be in good quality. The sources (`svg` or other) of the images must also be added to the `images` directory. We encourage you to use [yEd](https://www.yworks.com/products/yed) to easily generate diagrams and [Inkscape](https://inkscape.org/en/) for any other images.

## `slides` directory
{:.no_toc}

This directory contains introduction slide deck. There could be several slide decks, to cover different aspect. The slides are rendered using `remark.js` but written in Markdown to facilitate collaboration.

## `tutorials` directory
{:.no_toc}

This directory collects the tutorials related to the topic, one per subdirectory.

The tutorials are hands-on built for workshop and self-training, with description of the whole infrastructure needed to run the tutorial on any Galaxy instance (tools, data library, etc).

The templates for the tutorials are different from the other pages to help users to focus on the content of the tutorial. To improve the output of the tutorial, several metadata are mandatory for every tutorials, such as the requirements or the objectives of the tutorials. Boxes are also used to highlight some key points as the hands-on or the tips.

The content of each tutorial is generated with [Jekyll](https://jekyllrb.com/) from a Markdown file and some metadata (e.g. the requirements, the Zenodo link, the questions) defined inside the metadata of the related topic.

> Want to contribute to a tutorial? [Check out our training materials about that]({{ site.baseurl }}/topics/contributing/)

Sometimes, an hands-on tutorial is not the most appropriate format for a tutorial and slides are better. The content must be then added in the `slides` directory.

## `docker` directory
{:.no_toc}

For each topic, a flavored Docker image must integrate the tools needed for
the tutorials. The corresponding image must be based on official Galaxy Docker
images. We recommend to use the content of [`templates/docker`](templates/docker) as a template.

The `docker` image will also integrate the Galaxy tours available for each topics and the workflows.

> Want to learn more? [Check out our tutorial to build a Docker flavor for a tutorial]({{ site.baseurl }}/topics/training/tutorials/create-new-tutorial-docker/tutorial.html)

# Creating a new topic
{:.no_toc}

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Defining the topic

When we structured the repository, we decided to use as topics the categories that are used in the [ToolShed](https://toolshed.g2.bx.psu.edu/). The ToolShed assigns a category to each tool. Therefore, to decide where to put your tutorial, have a look at which ToolShed's category the main tools in your tutorial belong. For example, herethis tutorial will rely on the NCBI Blast+ tool.

> ### {% icon hands_on %} Hands-on: Defining the topic for the tutorial
>
> 1. Search for NCBI Blast+ on the [ToolShed](https://toolshed.g2.bx.psu.edu/)
> 2. Check in which category it has been placed
>
>    > ### {% icon question %} Questions
>    >
>    > In which topic will you put the new tutorial?
>    >
>    >    <details>
>    >    <summary>Click to view the answers</summary>
>    >    If we search for [NCBI Blast+ in the ToolShed](https://toolshed.g2.bx.psu.edu/view/devteam/ncbi_blast_plus/7538e2bfcd41), it is placed in 2 categories (bottom): "Next Gen Mappers", and "Sequence Analysis".
>    >    We decided to put it in "Sequence analysis" because this is the most general one for this tutorial.
>    >    </details>
>    {: .question}
{: .hands_on}

In this tutorial, we want to add a new topic called about "my-favorite-topic".

# Creating the directory for the topic

Once the topic name has been chosen, we can create it.

> ### {% icon hands_on %} Hands-on: Copy the required files
>
> 1. Copy the `templates` directory in `topics`
> 2. Rename the copied directory to `my-favorite-topic`
{: .hands_on}

# Make the templating system aware about the topic

We use Jekyll to generate the website out of the Markdown and YAML files. We need to tell Jekyll that there is a new topic by adding a symbolic link to the `metadata.yaml` inside the `metadata` folder.

> ### {% icon hands_on %} Hands-on: Add the new topic to the website
>
> 1. Go the `metadata` folder using the terminal: `cd metadata`
> 2. Add a symbolic link to the `metadata.yaml` file on our new topic: `ln -s ../topics/my-favorite-topic/metadata.yaml my-favorite-topic.yaml`
> 3. Move back to the root: `cd ..`
> 2. Make sure that Jekyll is running
>
>    > Want to learm how to start Jekyll? [Check out our tutorial to serve the website locally]({{ site.baseurl }}/topics/contributing/tutorials/running-jekyll/tutorial.html)
>
> 2. Check if the tutorial has been correctly added at [http://localhost:4000/training-material/](http://localhost:4000/training-material/)
>
{: .hands_on}

# Conclusion
{:.no_toc}

We just created a new topic. We can now fill it with new tutorials. To do that, you can follow our dedicated tutorials:

> ### Developing GTN training material
>
> This tutorial is part of a series to develop GTN training material, feel free to also look at:
>
> 1. [Setting up the tutorial infrastructure](../running-jekyll/tutorial.html)
> 1. [Writing content in markdown](../create-new-tutorial-content/tutorial.html)
> 1. [Defining metadata](../create-new-tutorial-metadata/tutorial.html)
> 1. [Creating a new topic](../create-new-topic/tutorial.html)
> 1. [Generating PDF handouts](../generating-pdf/tutorial.html)
> 1. [Creating Interactive Galaxy Tours](../create-new-tutorial-tours/tutorial.html)
> 1. [Defining technical requirements for a tutorial](../create-new-tutorial-technical/tutorial.html)
> 1. [Setting up Galaxy for training](../setup-galaxy-for-training/tutorial.html)
> 1. [Submitting the new tutorial to the GitHub repository](../github-command-line-contribution/slides.html)
> 1. [Good practices for teaching](../good-practices/slides.html)
{: .agenda}
