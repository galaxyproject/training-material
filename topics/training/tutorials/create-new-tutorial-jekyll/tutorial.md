---
layout: tutorial_hands_on
topic_name: training
tutorial_name: create-new-tutorial-jekyll
---

# Introduction
{:.no_toc}

Galaxy is a great solution to train the bioinformatics concepts:

- numerous bioinformatics tools are available (almost 5,000 in the ToolShed)
- it can be used by people without amy computer science skills
- it trains to use technology, outlining available resources and efforts that have made them accessible to researchers
- it is scalable

In 2016, the Galaxy Training Network decide to set up a new infrastructure for delivering easily Galaxy related training material. The idea was to develop something open and online based on a community effort, as always in Galaxy.

We took inspiration from [Software Carpentry](https://software-carpentry.org) and collected everything on a GitHub repository: [https://github.com/galaxyproject/training-material ](https://github.com/galaxyproject/training-material).
We decided on a structure based on tutorials with hands-on, fitting both for online self-training but also for workshops, grouped in topics. Each tutorial follows the same structure and comes with a virtualised isntance to run the training everywhere.

If you want to run the entire GTN website locally or test your new training material you can do this! 

Currently, the website is generated from the metadata and the tutorials using Jekyll, a simple static site builder.
We can use Jekyll to run a server to check if the tutorial is correctly added and rendered.


> ### Agenda
>
> In this tutorial, you will learn how to run a local instance of the GTN website:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Installation of the requirements

The first step is to install the requirements inside a conda environment. This step has to be done once.

> ### {% icon hands_on %} Hands-on: Install the requirements
>
> 1. Install [conda](https://conda.io/miniconda.html)
> 2. Navigate to `training-material/` folder
> 3. Create the conda environment: `conda env create -f environment.yml`
> 4. Activate the conda environment: `source activate galaxy_training_material`
> 5. Install Jekyll and related modules using [RubyGems](https://rubygems.org/pages/download): `make install`
{: .hands_on}

# Checking the website generation

Once Jekyll and its modules are installed in our conda environment, we can check the generation of the website locally:

> ### {% icon hands_on %} Hands-on: Checking the website generation locally
> 
> 1. (If not done) Activate the conda environment: `source activate galaxy_training_material`
> 1. Run a local Jekyll server with `make serve`
> 2. Visualize at [http://localhost:4000/ ](http://localhost:4000/)
>
>    > ### {% icon question %} Questions
>    >
>    > How to check if the server was started and if all topics are included?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    Please check [http://localhost:4000/topics/ ](http://localhost:4000/topics/) to get a list of topics.
>    >    </details>
>    {: .question}
{: .hands_on}

With `make serve`, a local Jekyll server will run in background. It will check the changes and regenerate the website accordingly. You may need to reload the page to see the changes (and sometimes to wait 1-2 minutes)

Once you are done, you can stop the server with `ctrl-c` and deactivate your conda environment with `source deactivate`.


# Generating PDF of the tutorials and slides

A PDF file of every tutorials and slide decks can be generated, using Chrome on command line and [Decktape](https://github.com/astefanutti/decktape).

> ### {% icon hands_on %} Hands-on: Checking the website generation locally
>
> 1. (If not done) Activate the conda environment: `source activate galaxy_training_material`
> 1. Install Chrome
>    - For OSX, install the [Chrome browser](https://www.google.com/chrome/browser/desktop/index.html)
>    - For Ubuntu, follow [these instructions](https://askubuntu.com/questions/510056/how-to-install-google-chrome#510186)
> 2. Run `make pdf`
> 3. Check the PDF in `_pdf` folder
{: .hands_on}


# (Not recommended) Using your fork for deployment and checking

We do not recommend it, but sometimes, you may want to have a version of the website on your fork (on the `gh-pages` branch) or to have the `gh-pages` in your fork is automatically updated once something is pushed on the master branch of your fork, as it is currently for the main GitHub repository (to check the changes, for example)

> ### {% icon hands_on %} Pushing the website on your fork
>
> 1. Activate `gh-pages` as branch on which the website is build:
>    1. Go on your GitHub repository
>    2. Go to "Settings"
>    3. Select "gh-pages branch" as source in the "GitHub Pages" section
> 2. Add your fork as a remote (if not already done): `git remote add fork "https://github.com/<your-github-username>/training-material"`
> 3. Build the website locally: `make build`
> 4. Prepare the `gh-pages` branch locally: `bin/publish`
> 5. Push on the `gh-pages` on your fork: `git push -f fork gh-pages`
> 6. Go back to your working branch: `git checkout --`
{: .hands_on}

> ### {% icon hands_on %} Checking the website on your fork
>
> 1. Activate `gh-pages` as branch on which the website is build
>    1. Go on your GitHub repository
>    2. Go to "Settings"
>    3. Select "gh-pages branch" as source in the "GitHub Pages" section
> 2. Generate a personal access token
>    1. Go on the "Settings" of your GitHub account (after clicking on your avatar on the top right of the GitHub interface)
>    2. Go the "Developer settings"
>    3. Go to "Personal access tokens"
>    4. Generate a new token by selecting "repo"
>    5. Copy the generated token
> 3. Encode your token
>    1. Install travis: `gem install travis`
>    2. Encrypt the token: `travis encrypt -r <your-github-username>/training-material`
>    3. Replace the keys in `secure` by the generated encrypted token into `.travis.yml`
> 4. Replace `bebatut` by your GitHub username in `bin/config-travis-git`
> 5. Push everything on the master branch of your fork
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
