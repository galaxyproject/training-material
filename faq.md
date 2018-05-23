---
layout: faq
---

1. TOC
{:toc}

# Overview Questions

## What is this website?

This website is a collection of tutorials developed and maintained by the [worldwide Galaxy community](https://galaxyproject.org/). You can learn more about this effort by reading our [article](https://www.biorxiv.org/content/early/2018/04/05/225680).

## What is Galaxy?

[Galaxy](https://galaxyproject.org/) is an open data integration and analysis platform for the life sciences, and it is particularly well-suited for data analysis training in life science research.

## What are the tutorials for?

These tutorials can be used for learning and teaching how to use Galaxy for general data analysis, and for specific domains from assembly to differential gene expression analysis with RNA-Seq data.

## What audiences are the tutorials for?

There are two distinct audiences for these materials.

1. **Self-paced individual learners.** These tutorials provide everything you need to learn a topic, from explanations of basic concepts to detailed hands-on exercises.
2. **Instructors.** They are also designed to be used by instructors in teaching/training settings.  Slides, and detailed tutorials are provided.  Some tutorials also include computational support such as Docker images that can be used to scale the lessons up to many participants.

## How is the content licensed?

The content of this website is Creative Commons Attribution 4.0 ([License]({{ site.github_repository }}/LICENSE.md)).

## How can I cite this effort?

We wrote an [article](https://www.biorxiv.org/content/early/2018/04/05/225680) about our effort there. You can cite it and then giving us some credit.


# For Individual Learners

Learning Galaxy and data analysis on your own, at your own pace?  This material is for you.

## Where do I start?

If you are new to Galaxy then start with one of the introductory topics.  These introduce you to concepts that are useful in Galaxy, no matter what domain you are doing analysis in.

If you are already familiar with Galaxy basics and want to learn how to use it in a particular domain (for example, ChIP-Seq), then start with one of those topics.

## How do I use this material?

Many topics include slide decks and if the topic you are interested in has slides then start there.  These will introduce the topic and important concepts.

Most of your learning will happen in the next step - the hands-on tutorials.  This is where you'll become familiar with using the Galaxy interface and experiment with different ways to use Galaxy and the tools in Galaxy.

## Where can I run the hands-on tutorials?

To run the hands-on tutorials you need a Galaxy server to run them on.

Some topics have a Docker image that can be installed and run on all participants' laptops.  These images contain Galaxy instances that include all tools and datasets used in a tutorial, as well as saved analyses and repeatable workflows that are relevant.

You can also run many tutorials on public Galaxy servers.  These servers are available to anyone on the world wide web and have all the tools that are needed by a specific tutorial.

If your organization (or consortia) has its own Galaxy server, then you may be able to run tutorials on that. Before you start you should confirm that all needed tools and reference genomes are available on your server.

Finally, you can also run your tutorials on cloud-based infrastructures.  Galaxy is available on many national research infrastructures such as Jetstream (United States), GenAP (Canada), GVL (Australia), CLIMB (United Kingdom), and more.  These instances are typically easy to launch before you start, and easy to shut down when you are done.

If you are already familiar with, and have an account on Amazon Web Services then you can also launch a Galaxy server there.

## How can I get help?

**HELP FOR INDIVIDUAL LEARNERS**

If you have questions about this training material, you can reach us using the [Gitter chat](https://gitter.im/Galaxy-Training-Network/Lobby).  You'll need a [GitHub](https://github.com/) or [Twitter](https://twitter.com/) account to post questions.  If you have questions about Galaxy outside the context of training, see the [Galaxy Support page](https://galaxyproject.org/support/).

# For Instructors

This material can also be used to teach Galaxy and data analysis in a group setting to students and researchers.

## Where do I start?

Spend some time exploring the different tutorials and the different resources that are available. Become familiar with the structure of the tutorials and think about how you might use them in your teaching.

## What Galaxy instance should I use for my training?

To teach the hands-on tutorials you need a Galaxy server to run the examples on

Some training topics have a Docker image that can be installed and run on all participants' laptops.  These images contain Galaxy instances that include all tools and datasets used in a tutorial, as well as saved analyses and repeatable workflows that are relevant.

You can also run many tutorials on several public Galaxy servers.  These servers are available to anyone on the world wide web and will also have all the tools that are needed by a specific tutorial.  If you choose this option then you should work with that server's admins to confirm that think the server can handle the workload for a workshop.

If your organization (or consortia) has its own Galaxy server, then you may be able to run tutorials on that.  This can be ideal because then the instance you are teaching on is the same you your participants will be using after the training.  They'll also be able to revisit any analysis they did during the training.  If you pursue this option you'll need to work with your organization's Galaxy Admins to confirm that

- the server can support a room full of people all doing the same analysis at the same time.
- all tools and reference datasets needed in the tutorial are locally installed.  Some tutorials include scripts that install all needed tools and data on an existing server.
- all participants will be able to create/use accounts on the system.

TODO: Training Infrastructure as a Service on useGalaxy.eu

Finally, you can also run your training on cloud-based infrastructures.  Galaxy has been available on Amazon Web Services for years, and is also available on many national research infrastructures such as Jetstream (United States), GenAP (Canada), GVL (Australia), CLIMB (United Kingdom), and so on.

## What are the best practices for teaching with Galaxy?

- Have a look at our [Good practices slides]()


## How do I get help?

The support channel for instructors is the same as for individual learners.  We suggest you start by posting a question to the Galaxy Training Network Gitter channel.  Anyone can view the discussion, but you'll need to login (using your GitHub or Twitter account) to add to the discussion.

If you have questions about Galaxy in general (that are not training-centric) then there are numerous support options.

# Contributing

First off, thanks for your interest in contributing to the Galaxy training materials!

Individual learners and instructors can make these training more effective by contributing back to them. You can report mistakes and errors, create more content, etc. Whatever is your background, there is a way to contribute.  See the [contributing page]() for how.

TODO: @tnabtaf suggests moving rest of this section to CONTRIBUTING.md

: via the GitHub website, via command-line or even without dealing with GitHub.


You can check our [tutorials](http://galaxyproject.github.io/training-material/topics/contributing) for more details.

We will address your issues and/or assess your change proposal as promptly as we can, and help you become a member of our community.

## How can I give feedback?

- tutorial level
- contributions
- issues


## How can I report mistakes or errors?

The easiest way to start contributing is to [file an issue]({{ site.github_repository }}/issues/new) to tell us about a problem such as a typo, spelling mistake, or a factual error. You can then introduce yourself and meet some of our community members.

## How can I fix mistakes or expand an existing tutorial using the GitHub interface?

Check our tutorial to learn how to use the GitHub interface

## How can I create new content without dealing with git?

If you feel uncomfortable with using the git and the [GitHub flow](https://guides.github.com/introduction/flow/), you can write a new tutorial with any text editor and then contact us. We will work together to integrate the new content.

## How can I contribute in "advanced" mode?

Most of the content is written in [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/) with some metadata (or variables) stored in [YAML](http://yaml.org/). To learn how to add new content, check out our [series of tutorials on creating new content]({{ site.baseurl }}/topics/contributing/):

- [Writing content in markdown to create a new tutorial]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)
- [Defining metadata for a new tutorial]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-metadata/tutorial.html)
- [Setting up the infrastructure]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-jekyll/tutorial.html) to run [Jekyll](https://jekyllrb.com/) and check the website generation
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

In [issues](https://github.com/galaxyproject/training-material/issues), you will find lists of issues to fix and features to implement (with the "newcomer-friendly" label for example). Feel free to solve them.

We strongly recommend you read and follow [Software Carpentry's](https://software-carpentry.org/) recommendations on [lesson design](https://swcarpentry.github.io/lesson-example/01-design/) and [lesson writing](https://swcarpentry.github.io/instructor-training/19-lessons/) if you plan to add or change some training materials, and also to check the [structure of the training material](#how-the-training-material-is-structured).

# Other Questions

## Are there any upcoming events focused on Galaxy Training?

Yes.  As of May 2018, these events are on the horizon:

* [CarpentryCon 2018](http://www.carpentrycon.org/), 30 May - 1 June, Dublin, Ireland
* Not specifically about Galaxy Training, but an excellent opportunity to gather with other computational science educators. [Bérénice Batut](@bebatut) will present a poster and lightning talk on [Community-Driven Training for Biological Data Analysis with the Galaxy Training Network](https://github.com/carpentries/carpentrycon/blob/master/Sessions/2018-05-30/13-Poster-Session/abstract-berenice-batut.md)
* [GCCBOSC 2018](https://gccbosc2018.sched.com/), June 25-30, Portland, Oregon, United States
  * The annual gathering of the Galaxy Community is an opportunity to learn from experienced Galaxy trainers and to contribute to these efforts:
    * [Bioinformatics Training and Education with the Galaxy Training Network](http://sched.co/Drp9), training session on how to use and contribute to these materials, presented by [Bérénice Batut](@bebatut)
    * [A fruitful year for the Galaxy Training materials], conference talk by [Bérénice Batut](@bebatut)
    * The *Galaxy documentation, analysis, and training (Galaxy DAT)* track of  [CollaborationFest](https://galaxyproject.org/events/gccbosc2018/collaboration/), June 29 - July 2. Focus on expanding Galaxy community resources like training materials and documentation.
    * Quarterly online training material Contribution Fests: The training community will meet online on the 3rd Friday of every 3rd month to focus on enhancing particular areas of the training material.

Is the above list now out of data (it happens).  See the [Galaxy Community Events Calendar](https://galaxyproject.org/events/) for what coming up right now.

# How is the training material structured?

TODO: Suggest moving this to contributing.

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



# What's not yet listed above.

- Interactive tours (where to cover these?)


