---
layout: tutorial_hands_on
topic_name: contributing
tutorial_name: generating-pdf
---

# Introduction
{:.no_toc}

The website with the training material can be run locally. Sometimes, it is also interesting to freeze the tutorials or to get PDFs of the tutorials.

> ### Agenda
>
> In this tutorial, you will learn how to run a local instance of the GTN website:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Generate PDFs artifact

To generate the PDFs, a command `make pdf` is given. This command:

- Launches a detached Jekyll server to serve the website
- Generates the PDFs of the tutorials by calling Chrome via command line
- Generates the PDFs of the slide decks by calling decktage

> ### {% icon hands_on %} Hands-on: Checking the website generation locally
>
> 1. (If not done) Activate the conda environment: `source activate galaxy_training_material`
> 2. Install Chrome
>    - For OSX, install the [Chrome browser]()
>    - For Ubuntu, follow [these instructions]()
> 1. Generate the PDFs: `make pdf`
> 2. Check the generated PDFs in `_pdf` folder
{: .hands_on}

# Conclusion
{:.no_toc}

> ### Developing GTN training material
>
> This tutorial is part of a series to develop GTN training material, feel free to also look at:
>
> 1. [Setting up the tutorial infrastructure](../running-jekyll/tutorial.html)
> 1. [Writing content in markdown](../create-new-tutorial-content/tutorial.html)
> 1. [Defining metadata](../create-new-tutorial-metadata/tutorial.html)
> 1. [Creating a new topic](../create-new-topic/tutorial.html)
> 1. [Generating PDF handouts](../generate-pdf/tutorial.html)
> 1. [Creating Interactive Galaxy Tours](../create-new-tutorial-tours/tutorial.html)
> 1. [Defining technical requirements for a tutorial](../create-new-tutorial-technical/tutorial.html)
> 1. [Setting up Galaxy for training](../setup-galaxy-for-training/tutorial.html)
> 1. [Submitting the new tutorial to the GitHub repository](../github-command-line-contribution/slides.html)
{: .agenda}
