---
layout: tutorial_hands_on

title: "Running the Galaxy Training material website locally"
questions:
  - "How to setup the infrastructure to build training webpages?"
objectives:
  - "Installing packages needed for rendering the webpage"
  - "Running the GTN material website locally"
  - "Tracking changes to the content live in the webbrowser"
time_estimation: "15min"
key_points:
  - "Checking the generated website can be done locally"
contributors:
  - bebatut
  - bgruening
  - shiltemann
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

The first step is to install the needed tools inside a conda environment. A conda environment is a directory that contains a specific collection of packages. For example here to run the website, we need ruby, pandas, requests, etc. By creating a conda environment and installing the needed tools there, we do not affect your main installation.

This step has to be done once.

> ### {% icon hands_on %} Hands-on: Install the requirements
>
> 1. Open a Terminal
> 1. Navigate to the `training-material/` folder with `cd`
> 3. Set up the conda environment
>     
>     It will install some needed tools (ruby, nodejs, etc) in a protected environment, without interfering with the existing tools or versions.
>
>     1. Install conda (if not already installed): `make install-conda`
>     2. Add conda to the path: `export PATH="$HOME/miniconda3/bin:$PATH"`
>     3. Create conda environment: `make create-env`
> 
> 4. Activate the conda environment: `source activate galaxy_training_material`
> 5. Install Jekyll and related modules using [RubyGems](https://rubygems.org/pages/download): `make install`
{: .hands_on}

> ### {% icon tip %} Troubleshooting
> If you encounter an error about libxml2 on linux, please try installing `libxml2-dev`, e.g. by the command `sudo apt-get install libxml2-dev`, and re-run `make install`
{: .tip}


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

# Stopping the server

Once you are done, you can stop the server and cleaning your repository.

> ### {% icon hands_on %} Hands-on: Stoping the server
>
> 1. Stop the server with <kbd>CTRL</kbd>-<kbd>C</kbd>
> 2. Deactivate your conda environment with `source deactivate`
> 3. Clean the repository: `make clean`
{: .hands_on}

# Conclusion
{:.no_toc}
