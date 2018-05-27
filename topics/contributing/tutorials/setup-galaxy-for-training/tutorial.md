---
layout: tutorial_hands_on
topic_name: contributing
tutorial_name: setup-galaxy-for-training
---

# Introduction
{:.no_toc}

In this tutorial, you will learn how to provision your Galaxy instance to support training modules from the GTN repository.

Tutorials in this repo are all supplemented with files describing the technical requirements to run them. This makes it easy to automate installation of tutorial requirements.
  - `tools.yaml`: describes the Tool Shed tools used in the tutorial
  - `data-library.yaml`: describes the input datasets
  - `data-manager.yaml`: describes the reference data required by tools
  - `workflows` folder: contains one or more workflows with all steps in the tutorial
  - `tours` folder: contains one or more yaml files describing interactive tours

For more information about how to create these files, please see our module on [specifying the technical requirements for your tutorial]({{ site.baseurl }}/topics/contributing/create-new-tutorial-technical/).

For just the list of ephemeris commands needed for installation, skip to the [Quickstart section](#Quickstart) at the end of this tutorial.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Prepare for installation

If you have a Galaxy server already running somewhere and would like to support one or more training modules, [ephemeris](https://ephemeris.readthedocs.io) can be used to easily install all the required tools, reference data, data libraries, tours and workflows.

For the purposes of this training, we will create a Galaxy instance running on `localhost:8080`, but these instructions can be adapted for use with any Galaxy instance by replacing this value with the URL or IP address of your Galaxy server.

## Start a local Galaxy

To setup a Galaxy server locally, we will first clone the Galaxy github repository, make a few small edits to the `galaxy.yaml` configuration file, and then start the server.


> ### {% icon hands_on %} Hands-on: Setup a local Galaxy instance
>
> 1. Clone the github repository
>    ```bash
>    git clone https://github.com/galaxyproject/galaxy.git
>    cd galaxy
>    ```
>
> 2. Add yourself as admin user in `config/galaxy.yaml`
>    ```bash
>    cp config/galaxy.yml.sample config/galaxy.yml
>    ```
>    open the `galaxy.yml` file with your favorite editor and edit the following line with your email address:
>    ```yaml
>    admin_users: user@example.eu
>    ```
>
> 3. Start Galaxy
>    ```bash
>    sh run.sh
>    ```
>    Galaxy will now install all its requirements, which may take a few minutes,
>    when all is finished installing, you should see something like this in your screen:
>
>    ```
>    Starting server in PID 9560.
>    serving on http://localhost:8080
>    ```
>
> 4. Open Galaxy
>     - Open a web browser
>     - Navigate to `localhost:8080` to access Galaxy
>
{: .hands_on}


### Find your Galaxy API key

In order to install the tutorial requirements, we will need the API key of an admin user.

> ### {% icon hands_on %} Hands-on: Obtain Galaxy API key
>
> 1. Register an account on Galaxy using the email address you added to the `config/galaxy.yml` file
>    - Once logged in, verify that you have a menu item named `Admin` in your top menu bar.
> 2. Go to your `User -> preferences -> Manage API key` in the top menu bar
> 3. Click on `Create a new key` to generate an API key
>    - Copy your API key to somewhere convenient, you will need it throughout this tutorial
>
{: .hands_on}


## Install Ephemeris

To install to training requirements to our Galaxy, we will use ephemeris, let's install it now:


> ### {% icon hands_on %} Hands-on: Install Ephemeris
>
> ```bash
> # optional: create a virtual environment
> virtualenv .venv; . .venv/bin/activate
>
> # install ephemeris
> pip install ephemeris
> ```
{: .hands_on}


# Installing tutorial requirements

We have created a small bash script to automatically install all the tutorial :

[`bin/install_tutorial_requirements.sh`]({{ site.github_repository }}/tree/master/bin/install_tutorial_requirements.sh)

> ### {% icon hands_on %} Hands-on: Obtain Galaxy API key
>
> 1. If you have not done so yet, clone the training material github repo:
>    ```bash
>    git clone https://github.com/galaxyproject/training-material.git
>    cd training-material
>    ```
>



Example usage from root of the repository:

```bash
bin/install_tutorial_requirements.sh topics/transcriptomics/tutorials/ref-based -g <galaxy url> -a <api key>
```

### Installing an entire topic

If you would like to install all the requirements for every tutorial within an entire topic, you can use the script in [`bin/install_topic_requirements.sh`]({{ site.github_repository }}/tree/master/bin/install_topic_requirements.sh)

Example usage from root of the repository:

```bash
bin/install_topic_requirements.sh topics/metagenomics -g <galaxy url> -a <api key>
```



### Installing a subset of components

If you would like to pick and choose what to install for each tutorial, below are descriptions of the commands used to install each of the components (tools, workflows, reference data, data libraries, tours)

**Installing Tools**

The ephemeris command to install tools defined in a `tools.yaml` file to a running Galaxy instance is:

```bash
shed-tools install -g <Galaxy url> -a <API key> -t topics/<topic>/tutorials/<tutorial>/tools.yaml
```


**Installing Workflows**

The ephemeris command to install a workflow or directory of workflows:

```bash
workflow-install --publish-workflows -g <Galaxy url> -a <API key> -w topics/<topic>/tutorials/<tutorial>/workflows
```

This command will install all the workflows in the `workflows` directory, but you may also specify a single workflow file here.

The `--publish_workflow` parameter will make the workflows available to anybody on the Galaxy instance.


**Installing Data Libraries**

The ephemeris command to populate a data library with the input datasets from Zenodo is:

```
setup-data-libraries -g <Galaxy url> -a <API key> -i topic/<topic>/tutorial/<tutorial>/data-library.yaml
```

**Installing Reference Data**

To run the data manager that installs and configures the reference data needed for a tutorial, run the following ephemeris command:

```bash
run-data-managers -g <Galaxy url> -a <API key> --config topic/<topic>/tutorial/<tutorial>/data-manager.yaml
```

**Installing Tours**

This is currently not possible using ephemeris, however, these can be installed by copying the files to the `config/plugins/tours/` directory of your Galaxy instance.

# Quickstart

Below is the list of commands used in this tutorial.

Using the scripts in this repository:

```bash
# Make sure you are in the root of the training-material repo
cd <training-materials repo root>

# install single tutorial
bin/install_tutorial_requirements.sh <Galaxy url> <API key> topics/<yourtopic>/tutorial/<yourtutorial

# install entire topic
bin/install_topic_requirements.sh <Galaxy url> <API key> topics/<yourtopic>

```

Using ephemeris directly:

```
# install tools
shed-tools install -g <Galaxy url> -a <API key> -t topics/<topic>/tutorials/<tutorial>/tools.yaml

# create data library with input datasets
setup-data-libraries -g <Galaxy url> -a <API key> -i topic/<topic>/tutorial/<tutorial>/data-library.yaml

# install refefernce data
run-data-managers -g <Galaxy url> -a <API key> --config topic/<topic>/tutorial/<tutorial>/data-manager.yaml

# install workflows
workflow-install --publish-workflows -g <Galaxy url> -a <API key> -w topics/<topic>/tutorials/<tutorial>/workflows

# install tours
copy the contents of the "tours" directory for the tutorial to Galaxy's "config/plugins/tours/"
```

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
