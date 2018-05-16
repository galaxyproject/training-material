---
layout: faq
---

1. TOC
{:toc}

# What is this website?

This website is a collection of tutorials developed and maintained by the worldwide Galaxy community. You can learn more about this effort with reading our [article](https://www.biorxiv.org/content/early/2018/04/05/225680).

# How can I get some help?

If you have any questions, you can reach us using the [Gitter chat](https://gitter.im/Galaxy-Training-Network/Lobby?utm_source=share-link&utm_medium=link&utm_campaign=share-link).

# How can I use these ressources to learn about Galaxy and data analyses?

Annotations

# How can I use these ressources as an instructor?

Good practices



## Docker

## Training Infrastructure as a Service

# How can I make sure a tutorial is running on a Galaxy instance?

Badge?
Installations


# How can I contribute?

First off, thanks for taking the time to contribute!

You can report mistakes or errors, create more contents, etc. Whatever is your background, there is probably a way to do it: via the GitHub website, via command-line or even without dealing with GitHub.

You can check our [tutorials](http://galaxyproject.github.io/training-material/topics/contributing) for more details.

We will address your issues and/or assess your change proposal as promptly as we can, and help you become a member of our community.

## How can I report mistakes or errors?

The easiest way to start contributing is to [file an issue]({{ site.github_repository }}/issues/new) to tell us about a spelling mistake or a factual error. You can then introduce yourself and meet some of our community members.

## How can I fix mistakes or expand an existing tutorial using the GitHub interface?

Check our tutorial to learn how to use the GitHub interface

## How can I create new content without dealing with git?

If you feel uncomfortable with using the git and the [GitHub flow](https://guides.github.com/introduction/flow/), you can write a new tutorial with any text editor.  and contact us: we will work together to integrate it.

## How can I contribute in "advanced" mode?

Most of the content is written in Markdown with some metadata (or variables) stored in YAML. To learn how to add new content, check out our [series of tutorials to create new content]({{ site.baseurl }}/topics/contributing/):

- [Writing content in markdown to create a new tutorial]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)
- [Defining metadata for a new tutorial]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-metadata/tutorial.html)
- [Setting up the infrastructure]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-jekyll/tutorial.html) to run Jekyll and check the website generation
- [Creating Interactive Galaxy Tours]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-tours/tutorial.html)
- [Building a Docker flavor for a tutorial]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-docker/tutorial.html)
- [Including a new topic]({{ site.baseurl }}/topics/contributing/tutorials/create-new-topic/tutorial.html)

To manage changes, we use [GitHub flow](https://guides.github.com/introduction/flow/) based on Pull Requests (check our [tutorial]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)):

1. [Create a fork](https://help.github.com/articles/fork-a-repo/) of this repository on GitHub
2. Clone your fork of this repository to create a local copy on your computer and init needed submodules (`git submodule init` and `git submodule update`)
3. Create a new branch in your local copy for each significant change
4. Commit the changes in that branch
5. Push that branch to your fork on GitHub
6. Submit a pull request from that branch to the [master repository](https://github.com/galaxyproject/training-material)
7. If you receive feedback, make changes in your local clone and push them to your branch on GitHub: the pull request will update automatically

Pull requests will be merged by the training team members after at least one other person has reviewed the Pull request and approved it.

## What can I do to help the project?

In [issues](https://github.com/galaxyproject/training-material/issues), you will find lists of issues to fix and features to change (with the "newcomer-friendly" label for example). Feel free to solve them.

We strongly recommend you read and follow Software Carpentry's recommendations on [lesson design](https://swcarpentry.github.io/lesson-example/01-design/) and [lesson writing](https://swcarpentry.github.io/instructor-training/19-lessons/) if you plan to add or change some training materials, and also to check the [structure of the training material](#how-the-training-material-is-structured).

# How is the training material structured?

Each training material is related to a topic. All training materials (slides, tutorials, ...) related to a topic are found in a dedicated directory (*e.g.* `transcriptomics` directory contains the material related to transcriptomis analysis). Each topic have the following structure:

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

You want to add a new topic? You can check our [dedicated tutorial]({{ site.baseurl }}/topics/contributing/tutorials/create-new-topic/tutorial.html) and also contact us to discuss about it.

# How is the content licensed?

The content of this website is Creative Commons Attribution 4.0 ([License]({{ site.github_repository }}/LICENSE.md)).

# How can I cite this effort?

We wrote an [article](https://www.biorxiv.org/content/early/2018/04/05/225680) about our effort there. You can cite it
