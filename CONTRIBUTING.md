Contributing to Galaxy Training material
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
- [How is the training material structured?](#how-is-the-training-material-structured)
    - [`images` directory](#images-directory)
    - [`tutorials` directory](#tutorials-directory)
    - [`slides` directory](#slides-directory)
    - [`docker` directory](#docker-directory)
- [How do I add new content?](#how-do-i-add-new-content)
- [How is the training material maintained?](#how-is-the-training-material-maintained)
    - [Maintainers](#maintainers)
    - [Labels](#labels)

# What should I know before I get started?

This repository is a project of unification of the Galaxy training material. You can find more information about this project in this [small presentation](https://galaxyproject.github.io/training-material/topics/training/slides/)

By contributing, you agree that we may redistribute your work under [this repository's license](LICENSE.md).

We will address your issues and/or assess your change proposal as promptly as we can, and help you become a member of our community.

# How can I contribute?

## Reporting mistakes or errors

The easiest way to start contributing is to file an issue to tell us about a spelling mistake or a factual error. You can then introduce yourself and meet some of our community members.

## Your first content contribution

Once you are feeling more comfortable, you can propose changes to this training material.

In [issues](https://github.com/galaxyproject/training-material/issues) and [project management system](https://github.com/galaxyproject/training-material/projects), you will find lists of issues to fix and features to change (with the "newcomer-friendly" label for example). Feel free to solve them.

We strongly recommend you read and follow Software Carpentry's recommendations on [lesson design](https://swcarpentry.github.io/lesson-example/01-design/) and [lesson writing](https://swcarpentry.github.io/instructor-training/19-lessons/) if you plan to add or change some training materials, and also to check the [structure of the training material](#how-the-training-material-is-structured).

## Pull Requests

To manage changes, we use [GitHub flow](https://guides.github.com/introduction/flow/) based on Pull Requests:

1. [Create a fork](https://help.github.com/articles/fork-a-repo/) of this repository on GitHub
2. Clone your fork of this repository to create a local copy on your computer and init needed submodules (`git submodule init` and `git submodule update`)
3. Create a new branch in your local copy for each significant change
4. Commit the changes in that branch
5. Push that branch to your fork on GitHub
6. Submit a pull request from that branch to the [master repository](https://github.com/galaxyproject/training-material)
7. If you receive feedback, make changes in your local clone and push them to your branch on GitHub: the pull request will update automatically

For beginners, the GitHub interface will help you in the process of editing a file. It will automatically create a fork of this repository where you can safely work and then submit the changes as a pull request without having to touch the command line.

Pull requests will be merged by the training team members after at least one other person has reviewed the Pull request and approved it.

# How is the training material structured?

Each training material is related to a topic. All training materials (slides, tutorials, ...) related to a topic are found in a dedicated directory (e.g. `Exome-seq` directory contains the material related to exome sequencing analysis). Each topic have the following structure:

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

> Want to contribute to a tutorial? [Check out our training content about that](http://galaxyproject.github.io/training-material/topics/training/)

Sometimes, an hands-on tutorial is not the most appropriate format for a tutorial and slides are better. The content must be then added in the `slides` directory.

## `docker` directory

For each topic, a flavored Docker image must integrate the tools needed for
the tutorials. The corresponding image must be based on official Galaxy Docker
images. We recommend to use the content of [`templates/docker`](templates/docker) as a template.

The `docker` image must also integrate a Galaxy tour from the [`tours` repository](https://github.com/galaxyproject/galaxy-tours)

> Want to learn more? [Check out our tutorial to build a Docker flavor for a tutorial](https://galaxyproject.github.io/training-material//topics/training/tutorials/create-new-tutorial-docker/tutorial.html)

# How do I add new content?

Most of the content is written in Markdown with some metadata (or variables) stored in YAML. To learn how to add new content, check out our [series of tutorials to create a new tutorial](http://galaxyproject.github.io/training-material/topics/training/):

- [Writing content in markdown to create a new tutorial](https://galaxyproject.github.io/training-material//topics/training/tutorials/create-new-tutorial-content/tutorial.html)
- [Defining metadata](https://galaxyproject.github.io/training-material//topics/training/tutorials/create-new-tutorial-metadata/tutorial.html)
- [Setting up the infrastructure](https://galaxyproject.github.io/training-material//topics/training/tutorials/create-new-tutorial-jekyll/tutorial.html) to run Jekyll and check the website generation
- [Creating Interactive Galaxy Tours](https://galaxyproject.github.io/training-material//topics/training/tutorials/create-new-tutorial-tours/tutorial.html)
- [Building a Docker flavor for a tutorial](https://galaxyproject.github.io/training-material//topics/training/tutorials/create-new-tutorial-docker/tutorial.html)

If you want to add a new topic, please contact us before: open an issue on this GitHub repository to discuss!

# How is the training material maintained?

## Maintainers

Each training topic and tutorial has one or two maintainers who act as editors.

They are responsible for making sure issues and change requests are looked at. They have the final say over what is included in the training material. But they are not responsible for writing/keeping up-to-date training material content or deciding what lessons ought to exist, both will be coming from the community.

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
