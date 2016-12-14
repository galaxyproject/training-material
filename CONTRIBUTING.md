gitContributing to Galaxy Training material
===

:+1::tada: First off, thanks for taking the time to contribute! :tada::+1:

The following is a set of guidelines for contributing to this training material on GitHub.

If you have any questions, you can reach us using the [Gitter chat](https://gitter.im/Galaxy-Training-Network/Lobby?utm_source=share-link&utm_medium=link&utm_campaign=share-link).

# Table of contents

- [What should I know before I get started?](#what-should-i-know-before-i-get-started)
- [How can I contribute?](#how-can-i-contribute)
    - [Reporting mistakes or errors](#reporting-mistakes-or-errors)
    - [Your first content contribution](#your-first-content-contribution)
    - [Pull Requests](#pull-requests)
- [How is the training material structured?](#how-the-training-material-is-structured)
    - [`images` directory](#images-directory)
    - [`tutorials` directory](#tutorials-directory)
    - [`slides` directory](#slides-directory)
    - [`docker` directory](#docker-directory)
    - [`README.md` file](#README.md-file)
- [How do I add new content?](#how-do-i-add-new-content)
    - [How do I add a new topic?](#how-do-i-add-a-new-topic)
    - [How do I add a new tutorial?](#how-do-i-add-a-new-tutorial)
    - [How do I fill a tutorial hands-on?](#how-do-i-fill-a-tutorial-hands-on)
    - [How do I fill introduction slides?](#how-do-i-fill-introduction-slides)
- [How is the training material maintained?](#how-the-training-material-is-maintained)
    - [Maintainers](#maintainers)
    - [Labels](#labels)

# What should I know before I get started?

This repository is a project of unification of the Galaxy training material. You can find more information about this project in this [small presentation](http://bgruening.github.io/training-material/shared/slides/project_presentation#/)

By contributing, you agree that we may redistribute your work under [our license](LICENSE.md).  

We will address your issues and/or assess your change proposal as promptly as we can, and help you become a member of our community.

# How can I contribute?

## Reporting mistakes or errors

The easiest way to start contributing is to file an issue to tell us about a spelling mistake or a factual error. You can then introduce yourself and meet some of our community members.

## Your first content contribution

Once you are feeling more comfortable, you can propose changes to this training material.

In [issues](https://github.com/bgruening/training-material/issues) and [project management system](https://github.com/bgruening/training-material/projects), you will find lists of issues to fix and features to change (with the "newcomer-friendly" label for example). Feel free to solve them.

We strongly recommend you read and follow Software Carpentry's recommendations on [lesson design](https://swcarpentry.github.io/lesson-example/01-design/) and [lesson writing](http://swcarpentry.github.io/instructor-training/19-lessons/) if you plan to add or change some training materials, and also to check the [structure of the training material](#how-the-training-material-is-structured).

## Pull Requests

To manage changes, we use [GitHub flow](https://guides.github.com/introduction/flow/) based on Pull Requests:

1. [Create a fork](https://help.github.com/articles/fork-a-repo/) of this repository on GitHub
2. Clone your fork of this repository to create a local copy on your computer and init needed submodules (`git submodule init` and `git submodule update`)
3. Create a new branch in your local copy for each significant change
4. Commit the changes in that branch
5. Push that branch to your fork on GitHub
6. Submit a pull request from that branch to the [master repository](https://github.com/bgruening/training-material)
7. If you receive feedback, make changes in your local clone and push them to your branch on GitHub: the pull request will update automatically

For beginners, the GitHub interface will help you in the process of editing a file. It will automatically create a fork of this repository where you can safely work and then submit the changes as a pull request without having to touch the command line.

# How is the training material structured?

Each training material is related to a topic. All training materials (slides, tutorials, ...) related to a topic are found in a dedicated directory (*e.g.* `Exome-seq` directory contains the material related to exome sequencing analysis). These repositories have to have the following structures (as in `Exome-seq` directory):

```
├── docker
│   ├── Dockerfile
│   ├── README.md
│   ├── tools.yaml
├── images
├── README.md
├── slides
├── tutorials
```

> Want to add a new topic? Check out [how to add new topic](#how-do-i-add-content).

## `images` directory

The `images` directory collects all images/pictures needed for the training materials related to the topic, *i.e* pictures for the slides or the tutorials.

Images shared between several topics are in the `shared/images` directory at the root.

All images for the slides must be in `images` directory. The images must be in good quality. The sources (`svg` or other) of the images must also be added to the `images` directory. We encourage you to use [yEd](http://www.yworks.com/products/yed) to easily generate diagrams and [Inkscape](https://inkscape.org/en/) for any other images.

## `tutorials` directory

This directory collects the tutorials related to the topic. The tutorials are hands-on built for workshop and self training. 

The template for the tutorials are different from the other pages to help users to focus on the content of the tutorial. To improve the output of the tutorial, several metadata are mandatory for every tutorials such as the requirements or the objectives of the tutorials. Boxes are also used to highlight some key points as the hands-on or the tips.  

The content of tutorial is generated with [Jekyll](http://jekyllrb.com/) from a markdown files and some metadata (*e.g.* the requirements, the Zenodo link, the questions) defined with the metadata of the related topic

> Want to contribute to a tutorial?
> - [Check out how to add a new tutorial?](#how-do-i-add-a-tutorial-hands-on)
> - [Check out how to fill a new tutorial?](#how-do-i-fill-a-tutorial-hands-on)

## `slides` directory

A slide deck is expected for every topic: the one with a general introduction of the topic. The slides are rendered using `remark.js` but written in markdown to facilate any contribution. [Check out how to fill introduction slides](#how-do-i-fill-introduction-slides).

## `docker` directory

For each topic, a flavored Docker image must integrate the tools needed for
the tutorials. The corresponding image must be based on official Galaxy Docker
images. We recommend to use the content of [`templates/docker`](templates/docker) as a template.

The `docker` image must also integrate a Galaxy tour from the [`galaxy-tours` repository](https://github.com/galaxyproject/galaxy-tours)

# How do I add new content?

Most of the content is written in markdown with some metadata (or variables) stored in YAML. To generate the website, we are using [Jekyll](http://jekyllrb.com/) and its templating system. 

So if you want to visualise locally how the website will look like, you need to run a local Jekyll server. So, Jekyll must be installed using [RubyGems](https://rubygems.org/pages/download):

```
$ make install
```

To run a local Jekyll server and visualize the changes, launch using the [Makefile](Makefile):

```
$ make serve
```

You can then visualize locally ([http://localhost:4000/](http://localhost:4000/)) the website before pushing your changes. 

## How do I add a new topic?

1. Add a `yml` file into the `metadata` directory similar to the one for [`RNA-Seq`](metadata/RNA-Seq.yml) and fill it with meta information on the topic
    - `name`: name of the topic (same name as the `yml` file and the directory)
    - `title`: title of the topic
    - `type`: targeted users (`"use"` or `""`) 
    - `summary`: summary of the content of the topic
    - `docker_image`: name of the [Docker image](#docker-directory) with the tools for this topic
    - `requirements`: list of requirements general for this topic, with a `title`, a `link` (relative for internal (inside training material) requirement or full for external requirement) and the type of link (`internal` or `external`)
    - `material`: list of material available for this topic

        For each material, you need to fill at least:

        - `title`
        - `type`: two possible types `introduction` or `tutorial`
            
            > For introduction material, [how to fill introduction slides](#how-do-i-fill-introduction-slides)

            > For tutorial material, check out [how to add a new tutorial](#how-do-i-add-a-new-tutorial)
        - `slides` (`"yes"` or `"no"`): tell if slides are available for this material

        
    - `maintainers`: the two maintainers of the topic with their `name`, `github_username`, `email`
    - `contributors`: list of people who contributed to the topic with `name`, `github_username`, `email`
    - `references`: list of references for this topic with `authors`, `title`, `link`, `summary`

    This information is used with [Jekyll](http://jekyllrb.com/) to generate the webpage related to the topic

2. Copy the template directory, rename it (with the same name as the `yml` file one) and fill it
    1. Change the `topic_name` in the `index.md` to fit the name of the directory and the name in `yml` file
    2. Add introduction slides

        > Check out [how to fill introduction slides](#how-do-i-fill-introduction-slides)

    3. Add tutorials

        > Check out [how to add a new tutorial](#how-do-i-add-a-new-tutorial)

## How do I add a new tutorial?

1. Add the metadata about the tutorial in `material` section in the `yml` file of the related topic that is in `metadata` directory
    - `title`: title of the tutorial
    - `type: "tutorial"`
    - `name`: name of the tutorial (name of the subdirectory where the files related to the tutorial will be stored)
    - `zenodo_link`: link on Zenodo to the input data for the tutorial (not ideal but it can be empty)
    - `galaxy_tour`: name of the galaxy tour
    - `hands_on`(`"yes"` or `"no"`): tell if an hands on is available for this material
    - `slides` (`"yes"` or `"no"`): tell if slides are available for this material
    - `questions`: list of questions that are adressed in the tutorial
    - `objectives`: list of objectives of the tutorial
    - `requirements`: list of requirements specific to this tutorial (in addition to the one of the topic), with a `title`, a `link` (relative for internal (inside training material) requirement or full for external requirement) and the type of link (`internal` or `external`)
    - `time_estimation`: estimation of the time needed to complete the hands-on
    - `key_points`: take home messages

    This information will appear in the top and bottom of the online hands-on generated using [Jekyll](http://jekyllrb.com/)

    ![](shared/images/tutorial_header.png)

2. Fill the hands-on file

    > Check out [how to fill it](#how-do-i-fill-a-tutorial-hands-on)

3. Add yourself as contributor for the topic in the `yml` file of the related topic that is in `metadata` directory

## How do I fill a tutorial hands-on?

1. Check that the metadata about the tutorial in the `yml` file in `metadata` directory are filled and correct

    They are used to generate the header and the footer of the tutorials. 

2. Fill the markdown file with the tutorial (after changing the `topic_name` and `tutorial_name`)

The content of a tutorial hands-on is written in `markdown`. They are rendered by [Jekyll](http://jekyllrb.com/) into the webpage for the tutorial. 

To improve the learning experience, we strongly recommend you to: 
- Add boxes to highlight:
    - Hands-on parts

        ```
        > ### :pencil2: Hands-on:
        >
        > 1. **Sort BAM dataset** :wrench:: Sort the paired-end BAM file by "Read names" with **Sort BAM dataset**
        {: .hands_on}
        ```

        ![](shared/images/tutorial_hand_on_box.png)

    - Questions (to make the learners think about what they are doing) and the collapsing and expanding answers

        ```
        > ### :question: Questions
        >
        > 1. Why are some tests filtered?
        > 2. Does it improve the *p*-value distribution?
        > 
        >    <details>
        >    <summary>Click to view answers</summary>
        >    Content goes here.
        >    </details>
        {: .question}
        ```

        ![](shared/images/tutorial_question_box.png)

    - Tips

        ```
        > ### :bulb: Tip: Importing data via links
        >
        > * Copy the link location
        > * Open the Galaxy Upload Manager
        > * Select **Paste/Fetch Data**
        > * Paste the link into the text field
        > * Press **Start**    
        {: .tip}
        ```

        ![](shared/images/tutorial_tip_box.png)

    - Comments

        ```
        > ### :nut_and_bolt: Comments
        > - Edit the "Database/Build" to select "dm3"
        > - Rename the datasets according to the samples
        {: .comment}
        ```

        ![](shared/images/tutorial_comment_box.png)

    To render the boxes correctly, the previous syntaxes have to be followed. The boxes can be nested, .e.g. for having tips inside hands-on.

- Add an agenda at the end of the introduction to indicate the plan of the tutorial
    
    ```
    > ### Agenda
    >
    > In this tutorial, we will analyze the data with:
    >
    > 1. [Pretreatments](#pretreatments)
    > 2. [Mapping](#mapping)
    > 3. [Analysis of the differential expression](#analysis-of-the-differential-expression)
    > {: .agenda}
    ```

    ![](shared/images/tutorial_agenda_box.png)

- Add pictures of the expected results
- Add at least one scheme or diagram to sum up the pipeline used at the end.

The input data required for the tutorials must be upload on [Zenodo](https://zenodo.org/) to obtain a dedicated DOI (in the [Galaxy training network community](https://zenodo.org/communities/galaxy-training/?page=1&size=20)).

You can also add yourself as contributor for the topic in the `yml` file of the related topic that is in `metadata` directory

## How do I fill introduction slides? 

Before starting filling the slides, you have to add the metadata about the tutorial in `material` section in the `yml` file of the related topic that is in `metadata` directory:

- `title`
- `type: "introduction"`
- `slides` (`"yes"` or `"no"`): tell if slides are available for this material

The introduction slides must be in the `index.html` file in `slides` directory for each topic. Even if the extension is `html`, slides are written in markdown. `---` is used to separate the slides.

```
---
layout: introduction_slides
topic_name: RNA-Seq
logo: "GTN"
---

# What is RNA sequencing?

---

### Where my data comes from?

![](../images/ecker_2012.jpg)

<small>[*Ecker et al, Nature, 2012*](http://www.nature.com/nature/journal/v489/n7414/full/489052a.html)</small>

---
```

The first slides (with the title, the requirements,...) are automatically generated using the metadata of the topic. Then the content to fill starts with the introduction. 

They are then rendered with [`Remark`](https://remarkjs.com/). Template for the `html` files can be found in
[`templates/slides/`](templates/slides/). Once the slides are on the `master` branch, they will be available at `http://bgruening.github.io/training-material/<topic>/slides/<slide_name>.html`

You can also add yourself as contributor for the topic in the `yml` file of the related topic that is in `metadata` directory

# How is the training material maintained?

## Maintainers

Each training topic has one or two maintainers who act as editors. They are responsible for making sure issues and change requests are looked at. They have the final say over what is included in the training material related to the topic. But they are not responsible for writing training material content or deciding what lessons ought to exist, both will be coming from the community.

## Labels

This repository is using the following labels for issues, pull requests and project management:

- Type
    - `bug`: errors to be fixed
    - `improvement`: enhancement to an existing functionality
    - `feature`: new functionality
    - `discussion`: discussion threads
    - `question`: often turn into discussion threads
- Status
    - `help-wanted`: requests for assistance
    - `newcomer-friendly`: suitable for people who want to start contributing
    - `work-in-progress`: someone is working on this
    - `review-needed`: requests for review
- Topic: each topic has its own label to easily relate the issue or pull request to the topic, but the label `template-and-tools` is used for questions/issues/pull requests related to the templates and tools rather than the lessons themselves.
