---
layout: tutorial_hands_on
topic_name: training
tutorial_name: create-new-tutorial-tours
---

# Introduction
{:.no_toc}

Galaxy is a great solution to train the bioinformatics concepts:

- numerous bioinformatics tools are available (almost 5,000 in the ToolShed)
- it can be used by people without amy computer science skills
- it trains to use technology, outlining available resources and efforts that have made them accessible to researchers
- it is scalable

In 2016, the Galaxy Training Network decide to set up a new infrastructure for delivering easily Galaxy related training material. The idea was to develop something open and online based on a community effort, as most of the time in Galaxy.

We take inspiration from [Software Carpentry](https://software-carpentry.org). We collected everything on a GitHub repository: [https://github.com/galaxyproject/training-material ](https://github.com/galaxyproject/training-material). We decided a structure based on tutorials with hands-on, fitting both for online self-training but also for workshops, grouped in topics. Each tutorial follows the same structure and comes with a technical support to be able to run.

In this tutorial, you will understand how to design and develop a new tutorial fitting in this training material repository. As doing helps to understand, we will develop a small tutorial to explain BLAST with the full infrastructure to be able to run this tutorial anywhere.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# A Galaxy Interactive Tour

A Galaxy Interactive Tour is a way to go through an entire analysis, step by step inside Galaxy in an interactive and explorative way. It is a great pedagogic way to run the tutorial directly inside Galaxy.

![Demonstration of an Interactive Tour](../../../dev/images/galaxy_tour_demo.gif "Demonstration of an Interactive Tour")

A tour is a YAML file like:

```
id: galaxy_ui
name: Galaxy UI
description: A gentle introduction to the Galaxy User Interface
title_default: "Welcome to Galaxy"

steps:
    - title: "Welcome to Galaxy"
      content: "This short tour will guide you through Galaxy's user interface.<br>
                You can navigate with your arrow keys and leave the tour at any time point
                with 'Escape' or the 'End tour' button."
      backdrop: true

    - title: "Upload your data"
      element: ".upload-button"
      intro: "Galaxy supports many ways to get in your data.<br>
              Use this button to upload your data."
      position: "right"
      postclick:
        - ".upload-button"
```

- at the top some metadata related to the Tour:
    - `id`: ID of the tour
    - `name`: name of the tour
    - `description`: a short description of the tour
    - `title_default`: a title
- several steps corresponding to the different boxes

    Each step is beginning with a dash `-` and with possible arguments

    Argument | Description
    ---  | ---
    `title`  | Header of each step-container
    `content` | Text that is shown to the user
    `element` | [JQuery Selector](https://api.jquery.com/category/selectors/) of the element you want to describe / click
    `placement` | Placement of the text box relative to the selected element
    `preclick` or `postclick` | Elements that recieve a click() event before (`preclick`) or after (`postclick`) the step is shown
    `textinsert` | Text to insert if element is a text box (e.g. tool search or upload)
    `backdrop` | `true/false`:  Show a dark backdrop behind the popover and its element, highlighting the current step

    [Full reference of the properties](https://bootstraptour.com/api/)

The YAML file of a tour can be integrated in a Galaxy instance by placing the YAML file in the `config/plugins/tours` directory of the Galaxy code and restarting the Galaxy instance

# Creating a Galaxy Interactive Tour

[A Web browser plugin](https://github.com/TailorDev/galaxy-tourbuilder) is available to help the creation and the test (on the fly) of an interactive tour.

<blockquote class="imgur-embed-pub" lang="en" data-id="a/0YVvz"><a href="//imgur.com/a/0YVvz">Galaxy Tour Builder by TailorDev</a></blockquote><script async src="//s.imgur.com/min/embed.js" charset="utf-8"></script>

> ### {% icon hands_on %} Hands-on: Install and start the plugin
>
> 1. Install the plugin using the app store of your web-browser:
>     - [Chrome Web Store](https://chrome.google.com/webstore/detail/galaxy-tour-builder/mdfbapknmcpnbmggahhaegehbbbmhmgg)
>     - [Mozilla Add-ons (Firefox)](https://addons.mozilla.org/en-US/firefox/addon/galaxy-tour-builder/)
>     - [Opera add-ons](https://addons.opera.com/en/extensions/details/galaxy-tour-builder/)
>
> 2. Load the webpage of any Galaxy instance
> 3. Start the plugin by clicking on the icon with Galaxy icon close to the address bar
{: .hands_on}

We can now create easily a Galaxy Interactive Tour and test it on the fly.

> ### {% icon hands_on %} Hands-on: Create a Galaxy Interactive Tour
>
> 1. Create a Galaxy Interactive Tour for "BLAST" tutorial
> 2. Test it with the plugin
> 3. Copy the YAML content and add it to a file
> 2. Add the file to the `tours` directory of the tutorial
> 3. Test it on a local Galaxy instance
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
