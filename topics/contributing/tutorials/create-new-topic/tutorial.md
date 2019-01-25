---
layout: tutorial_hands_on

title: "Including a new topic"
questions:
  - "How to include a new topic?"
objectives:
  - "Create a new topic"
  - "Set up the metadata for a topic"
time_estimation: "30m"
key_points:
  - "A new topic can be easily added for new tutorials"
contributors:
  - bebatut
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
images. We recommend to use the content of [`templates/docker`]({{ site.github_repository }}/tree/master/templates/docker) as a template.

The `docker` image will also integrate the Galaxy tours available for each topics and the workflows.

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

When we structured the repository, we decided to use as topics the categories that are used in the [ToolShed](https://toolshed.g2.bx.psu.edu/). The ToolShed assigns a category to each tool. Therefore, to decide where to put your tutorial, have a look at which ToolShed's category the main tools in your tutorial belong. For example, this tutorial will rely on the NCBI Blast+ tool.

> ### {% icon hands_on %} Hands-on: Defining the topic for the tutorial
>
> 1. Search for NCBI Blast+ on the [ToolShed](https://toolshed.g2.bx.psu.edu/)
> 2. Check in which category it has been placed
>
>    > ### {% icon question %} Questions
>    >
>    > In which topic will you put the new tutorial?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > If we search for [NCBI Blast+ in the ToolShed](https://toolshed.g2.bx.psu.edu/view/devteam/ncbi_blast_plus/7538e2bfcd41), it is placed in 2 categories (bottom): "Next Gen Mappers", and "Sequence Analysis".
>    > > We decided to put it in "Sequence analysis" because this is the most general one for this tutorial.
>    > {: .solution}
>    {: .question}
{: .hands_on}

In this tutorial, we want to add a new topic called about "my-favorite-topic".

# Creating the skeleton for the topic

Once the topic name has been chosen, we can create it.

> ### {% icon hands_on %} Hands-on: Create all the required files and folders structures automatically
>
> 1. Open a terminal
> 2. Run (by adapting the information between the quotes)
>
>    ```
>    $ planemo training_init \
>             --topic_name "my-favorite-topic" \
>             --topic_title "Test" \
>             --topic_target "use" \
>             --topic_summary "Summary of the topic"
>    ```
>
> 3. Check that a new directory has been generated in `topics`
> 4. Check that a YAML file with your topic name has been generated in `metadata` folder
> 5. Make sure that Jekyll is running
>
>    > Want to learn how to start Jekyll? [Check out our tutorial to serve the website locally]({{ site.baseurl }}/topics/contributing/tutorials/running-jekyll/tutorial.html)
>
> 6. Check if the topic has been correctly added at [http://localhost:4000/training-material/](http://localhost:4000/training-material/)
>
{: .hands_on}

# Adapt the metadata for your topic

Several metadata are defined in `metadata.yaml` file in your topic folder to :

- `name`: name of the topic (name of the folder)
- `title`: title of the topic (the one displayed on the webpage)
- `type`: target for the topic ('use', 'admin-dev', 'instructors')
- `summary`: summary of the focus of the topic
- `requirements`: list of resources that the reader of the material should be familiar with before starting any tutorial in this topic:
    - `type`: the type of link (`internal` or `external`)

    For internal, i.e. inside the Galaxy Training Material:
    - `topic_name`: name of the topic
    - `tutorials`: list of required tutorials inside of the topic

    For external:
    - `title`: title of the external resource
    - `link`: URL to the external resource

- `docker_image`: name of the Docker image for the topic

    If no Docker image exists for this topic, let this information empty

- `maintainers`: GitHub username of people maintaining the topic

> ### {% icon hands_on %} Hands-on: Update the new topic to the website
>
> 1. Open the `metadata.yaml` file in your topic folder
> 2. Fill the correct metadata of the topic
> 3. Make sure that Jekyll is running
>
>    > Want to learn how to start Jekyll? [Check out our tutorial to serve the website locally]({{ site.baseurl }}/topics/contributing/tutorials/running-jekyll/tutorial.html)
>
> 4. Check how it changes the local website
>
{: .hands_on}

# Conclusion
{:.no_toc}

We just created a new topic. We can now fill it by [creating new tutorials]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial/tutorial.html)
