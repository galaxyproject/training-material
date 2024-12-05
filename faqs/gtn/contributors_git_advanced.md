---
title: How can I contribute in "advanced" mode?
area: contributors
layout: faq
box_type: tip
contributors: [shiltemann,hexylena]
---


Most of the content is written in [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/) with some metadata (or variables) found in [YAML](http://yaml.org/) files. Everything is stored on our [GitHub repository]({{ site.github_repository }}). Each training material is related to a topic. All training materials (slides, tutorials, etc) related to a topic are found in a dedicated directory (*e.g.* `transcriptomics` directory contains the material related to transcriptomic analysis). Each topic has the following structure:

![Structure of the repository]({% link shared/images/repo_organization.png %}){: width="400px"}

- a metadata file in YAML format
- a directory with the topic introduction slide deck in Markdown with introductions to the topic
- a directory with the tutorials:

    Inside the tutorials directory, each tutorial related to the topic has its own subdirectory with several files:
    - a tutorial file written in Markdown with hands-on
    - an optional slides file in Markdown with slides to support the tutorial
    - a directory with Galaxy Interactive Tours to reproduce the tutorial
    - a directory with workflows extracted from the tutorial
    - a YAML file with the links to the input data needed for the tutorial
    - a YAML file with the description of needed tools to run the tutorial

- a directory with the Dockerfile describing the details to build a container for the topic (self-study environments).

To manage changes, we use [GitHub flow](https://guides.github.com/introduction/flow/) based on Pull Requests (check our [tutorial]({% link topics/contributing/tutorials/create-new-tutorial-content/tutorial.md %})):

1. [Create a fork](https://help.github.com/articles/fork-a-repo/) of this repository on GitHub
2. Clone your fork of this repository to create a local copy on your computer and initialize the required submodules (`git submodule init` and `git submodule update`)
3. Create a new branch in your local copy for each significant change
4. Commit the changes in that branch
5. Push that branch to your fork on GitHub
6. Submit a pull request from that branch to the [original repository]({{ site.github_repository }})
7. If you receive feedback, make changes in your local clone and push them to your branch on GitHub: the pull request will update automatically
8. Pull requests will be merged by the training team members after at least one other person has reviewed the Pull request and approved it.

Globally, the process of development of new content is open and transparent:

1. Creation of a branch derived from the main branch of the GitHub repository
2. Initialization of a new directory for the tutorial
3. Filling of the metadata with title, questions, learning objectives, etc
4. Generation of the input dataset for the tutorial
5. Filling of the tutorial content
6. Extraction of the workflows of the tutorial
7. Automatic extraction of the required tools to populate the tool file
8. Automatic annotation of the public Galaxy servers
9. Generation of an interactive tour for the tutorial with the [Tourbuilder](https://tailordev.fr/blog/2017/07/19/the-galaxy-tour-builder-extension/) web-browser extension
10. Upload of the datasets to Zenodo and addition of the links in the data library file.
11. Once ready, opening a Pull Request
12. Automatic checks of the changes are automatically checked for the right format and working links using continuous integration testing on Travis CI
13. Review of the content by several other instructors via discussions
14. After the review process, merge of the content into the main branch, starting a series of automatic steps triggered by Travis CI
15. Regeneration of the website and publication on [{{ site.url }}{{ site.baseurl }}/]({{ site.baseurl }}/)
16. Generation of PDF artifacts of the tutorials and slides and upload on the FTP server
18. Population of [TeSS](https://tess.elixir-europe.org/), the ELIXIR’s Training Portal, via the metadata

![Development process]({% link shared/images/development_process.png %})

To learn how to add new content, check out our [series of tutorials on creating new content]({% link topics/contributing/index.md %}):

{% assign topic = site.data["contributing"] %}
<ol>
{% assign topic_material = site | topic_filter:'contributing' %}
{% for material in topic_material %}
 {% if material.enable != "false" %}
  {% if material.type == "introduction" %}
<li><a href="{{ site.baseurl }}/topics/{{ topic.name }}/slides/{{ material.tutorial_name }}.html">{{ material.title }}</a></li>
 {% elsif material.type == "tutorial" %}
  {% if material.hands_on %}
<li><a href="{{ site.baseurl }}/topics/{{ topic.name  }}/tutorials/{{ material.tutorial_name }}/tutorial.html">{{ material.title }}</a></li>
  {% elsif material.slides %}
<li><a href="{{ site.baseurl }}/topics/{{ topic.name }}/tutorials/{{ material.tutorial_name }}/slides.html">{{ material.title }}</a></li>
   {% endif %}
  {% endif %}
 {% endif %}
{% endfor %}
 </ol>

We also strongly recommend you read and follow [The Carpentries](https://carpentries.org/) recommendations on [lesson design](https://carpentries.github.io/lesson-example/01-design/) and [lesson writing](https://carpentries.github.io/instructor-training/15-lesson-study/) if you plan to add or change some training materials, and also to check the structure of the training material above.


