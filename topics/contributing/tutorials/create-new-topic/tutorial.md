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

> Want to add a new topic? Please contact us before: open an issue on this GitHub repository to discuss and we will help you in this process

## `images` directory

The `images` directory collects all images/pictures needed for the training materials related to the topic, *i.e* pictures for the slides or the tutorials.

Images shared between several topics are in the `shared/images` directory at the root.

All images for the slides must be in `images` directory. The images must be in good quality. The sources (`svg` or other) of the images must also be added to the `images` directory. We encourage you to use [yEd](https://www.yworks.com/products/yed) to easily generate diagrams and [Inkscape](https://inkscape.org/en/) for any other images.

## `slides` directory

A slide deck is expected for every topic: the one with a general introduction of the topic. The slides are rendered using `remark.js` but written in Markdown to facilitate collaboration.

## `tutorials` directory

This directory collects the tutorials related to the topic, one per subdirectory. The tutorials are hands-on built for workshop and self-training, with description of the whole infrastructure needed to run the tutorial on any Galaxy instance (tools, data library, etc).

The templates for the tutorials are different from the other pages to help users to focus on the content of the tutorial. To improve the output of the tutorial, several metadata are mandatory for every tutorials, such as the requirements or the objectives of the tutorials. Boxes are also used to highlight some key points as the hands-on or the tips.

The content of each tutorial is generated with [Jekyll](https://jekyllrb.com/) from a Markdown file and some metadata (e.g. the requirements, the Zenodo link, the questions) defined inside the metadata of the related topic.

> Want to contribute to a tutorial? [Check out our training content about that](http://galaxyproject.github.io/training-material/topics/contributing/)

Sometimes, an hands-on tutorial is not the most appropriate format for a tutorial and slides are better. The content must be then added in the `slides` directory.

## `docker` directory

For each topic, a flavored Docker image must integrate the tools needed for
the tutorials. The corresponding image must be based on official Galaxy Docker
images. We recommend to use the content of [`templates/docker`](templates/docker) as a template.

The `docker` image must also integrate a Galaxy tour from the [`tours` repository](https://github.com/galaxyproject/galaxy-tours)

> Want to learn more? [Check out our tutorial to build a Docker flavor for a tutorial](https://galaxyproject.github.io/training-material//topics/training/tutorials/create-new-tutorial-docker/tutorial.html)

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Creating a new topic

symbolic links

# Conclusion
{:.no_toc}

