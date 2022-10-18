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
contributions:
  authorship:
  - bebatut
  editing:
  - hexylena
---

# Introduction


Each training material is related to a topic. All training materials (slides, tutorials, ...) related to a topic are found in a dedicated directory (*e.g.* `transcriptomics` directory contains the material related to exome sequencing analysis).

## Directory structure

Each topic has the following structure:

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
│   │   ├── data-library.yaml
│   │   ├── workflows
│   │   │   ├── index.md
│   │   │   ├── workflow.ga
│   │   ├── tours
│   │   │   ├── tour.yaml
```


### `images` directory


The `images` directory collects all images/pictures needed for the training materials related to the topic, *i.e* pictures for the slides or the tutorials.

Images shared between several topics are in the `shared/images` directory at the root.

All images for the slides must be in `images` directory. The images must be in good quality. The sources (`svg` or other) of the images must also be added to the `images` directory. We encourage you to use [yEd](https://www.yworks.com/products/yed) to easily generate diagrams and [Inkscape](https://inkscape.org/en/) for any other images.

### `slides` directory


This directory contains introduction slide deck. There could be several slide decks, to cover different aspect. The slides are rendered using `remark.js` but written in Markdown to facilitate collaboration.

### `tutorials` directory


This directory collects the tutorials related to the topic, one per subdirectory.

The tutorials are hands-on built for workshop and self-training, with description of the whole infrastructure needed to run the tutorial on any Galaxy instance (tools, data library, etc).

The templates for the tutorials are different from the other pages to help users to focus on the content of the tutorial. To improve the output of the tutorial, several metadata are mandatory for every tutorials, such as the requirements or the objectives of the tutorials. Boxes are also used to highlight some key points as the hands-on or the tips.

The content of each tutorial is generated with [Jekyll](https://jekyllrb.com/) from a Markdown file and some metadata (e.g. the requirements, the Zenodo link, the questions) defined inside the metadata of the related topic.

> <comment-title>Contributing</comment-title>
> Want to contribute to a tutorial? Check out [our training materials about that]({% link topics/contributing/index.md %}).
{: .comment}

Sometimes, an hands-on tutorial is not the most appropriate format for a tutorial and slides are better. The content must be then added in the `slides` directory.

### `docker` directory


For each topic, a flavored Docker image must integrate the tools needed for
the tutorials. The corresponding image must be based on official Galaxy Docker
images.

The `docker` image will also integrate the Galaxy tours available for each topics and the workflows.

# Creating a new topic


> <agenda-title></agenda-title>
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Defining the topic

When we structured the repository, we decided to use as topics the categories that are used in the [ToolShed](https://toolshed.g2.bx.psu.edu/). The ToolShed assigns a category to each tool. Therefore, to decide where to put your tutorial, have a look at which ToolShed's category the main tools in your tutorial belong. For example, this tutorial will rely on the NCBI Blast+ tool.

> <hands-on-title>Defining the topic for the tutorial</hands-on-title>
>
> 1. Search for NCBI Blast+ on the [ToolShed](https://toolshed.g2.bx.psu.edu/)
> 2. Check in which category it has been placed
>
>    > <question-title></question-title>
>    >
>    > In which topic will you put the new tutorial?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > If we search for [NCBI Blast+ in the ToolShed](https://toolshed.g2.bx.psu.edu/view/devteam/ncbi_blast_plus/7538e2bfcd41), it is placed in 2 categories (bottom): "Next Gen Mappers", and "Sequence Analysis".
>    > > We decided to put it in "Sequence analysis" because this is the most general one for this tutorial.
>    > {: .solution}
>    {: .question}
{: .hands_on}

In this tutorial, we want to add a new topic called about "my-favorite-topic".

# Creating the skeleton for the topic

Once the topic name has been chosen, we can create it.

> <hands-on-title>Create all the required files and folders structures automatically</hands-on-title>
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
>    > <comment-title>Jekyll</comment-title>
>    > Want to learn how to start Jekyll? [Check out our tutorial to serve the website locally]({% link topics/contributing/tutorials/running-jekyll/tutorial.md %})
>    {: .comment}
>
> 6. Check if the topic has been correctly added at [http://localhost:4000/training-material/](http://localhost:4000/training-material/)
>
{: .hands_on}

# Adapt the metadata for your topic

Several metadata are defined in `metadata.yaml` file in your topic folder to :

{% assign kid_key = "Topic Schema" %}
{% assign kid_val = site.data['schema-topic'] %}
{% include _includes/schema-render.html key=kid_key value=kid_val %}

- `subtopics`: for large topics, we can define subtopics and create multiple tutorial lists:
  ```
  subtopics:
    singlecell:
      title: "Single Cell Analysis"
      description: "These tutorials cover single cell analysis"
    small:
      title: "Small RNA"
      description: "These Tutorial"
  ```
  tutorials can be assigned to subtopics by adding e.g. `subtopic: singlecell` to the tutorial metadata. An example of this subtopic division can be found in the [admin section]({{site.baseurl}}/topics/admin/ )

> <hands-on-title>Update the new topic to the website</hands-on-title>
>
> 1. Open the `metadata.yaml` file in your topic folder
> 2. Fill the correct metadata of the topic
> 3. Make sure that Jekyll is running
>
>    > <comment-title>Jekyll</comment-title>
>    > Want to learn how to start Jekyll? [Check out our tutorial to serve the website locally]({% link topics/contributing/tutorials/running-jekyll/tutorial.md %})
>    {: .comment}
>
> 4. Check how it changes the local website
>
{: .hands_on}

# Conclusion


We just created a new topic. We can now fill it by [creating new tutorials]({% link topics/contributing/tutorials/create-new-tutorial/tutorial.md %})
