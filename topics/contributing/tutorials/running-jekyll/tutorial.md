---
layout: tutorial_hands_on
topic_name: contributing
tutorial_name: running-jekyll
---

# Introduction
{:.no_toc}

If you want to run the entire GTN material website locally or test your new training material you can do this!

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
> 2. Visualize at [http://localhost:4000/training-material/ ](http://localhost:4000/training-material/)
>
{: .hands_on}

With `make serve`, a local Jekyll server will run in background. It will check the changes and regenerate the website accordingly. You may need to reload the page to see the changes (and sometimes to wait 1-2 minutes).

# Stoping the server

Once you are done, you can:
- stop the server with <kbd>CTRL</kbd>-<kbd>C</kbd>
- deactivate your conda environment with `source deactivate`
- clean the repository: `make clean`

# Conclusion
{:.no_toc}
