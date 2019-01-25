---
layout: tutorial_hands_on

title: Set up a Galaxy for Training
time_estimation: 2h
questions:
  - How do I prepare my Galaxy instance to support a training module?
  - How can I generate a Docker Galaxy instance for my topic?
objectives:
  - Use ephemeris to install the training requirements to a Galaxy instance
  - Create a docker image for a training topic
key_points:
  - Technical requirements have been defined for all the training materials
  - Ephemeris can be used to automatically install these requirements to an existing Galaxy
  - Convenience scripts are provided in the training material repository allow for easy installation
  - Docker images can easily be created per topic
contributors:
  - shiltemann
  - bebatut
---

# Introduction
{:.no_toc}

In this tutorial, you will learn how to provision your Galaxy instance to support training modules from the GTN repository.

Tutorials in this repository are all supplemented with files describing the technical requirements to run them. This makes it easy to automate installation of tutorial requirements.
  - `tools.yaml`: describes the Tool Shed tools used in the tutorial
  - `data-library.yaml`: describes the input datasets
  - `data-manager.yaml`: describes the reference data required by tools
  - `workflows` folder: contains one or more workflows with all steps in the tutorial
  - `tours` folder: contains one or more yaml files describing interactive tours

> For more information about how to create these files, please see our module on [specifying the technical requirements for your tutorial]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-technical/tutorial.html).

For just the list of Ephemeris commands needed for installation, skip to the [Quickstart section](#quickstart) at the end of this tutorial.

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
>    admin_users: user@example.org
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


## Find your Galaxy API key

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

To install to training requirements to our Galaxy, we will use Ephemeris, let's install it now:


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

We have created a small bash script to automatically install all of a tutorial's requirements to an existing Galaxy. It's located in this repository under: [`bin/install_tutorial_requirements.sh`]({{ site.github_repository }}/tree/master/bin/install_tutorial_requirements.sh)

In this example we will install the requirements for the [*Reference-based RNASeq*]({{ site.baseurl }}/topics/transcriptomics/tutorials/ref-based/tutorial.html) tutorial to the Galaxy instance running on localhost.


> ### {% icon hands_on %} Hands-on: Install a tutorial
>
> 1. If you have not done so yet, clone the training material github repo:
>    ```bash
>    git clone https://github.com/galaxyproject/training-material.git
>    cd training-material
>    ```
>
> 2. Run the script to install the RNASeq tutorial
>
>    ```bash
>    bin/install_tutorial_requirements.sh topics/transcriptomics/tutorials/ref-based http://localhost:8080 <api key>
>    ```
>
{: .hands_on}

Installation may take some time. This script will automatically install the tools, create a data library and populate it with the input datasets from Zenodo, install and publish the workflows, and run any data managers that might be required.

The only thing the script cannot currently automate is the installation of the interactive tours. We will now do this manually by copying the contents of the `tours` folder to our Galaxy instance, in the folder `$GALAXY_ROOT/config/plugins/tours`

> ### {% icon hands_on %} Hands-on: Install the interactive Tours
>
> 1. Copy the `tour.yaml` file from the training materials repo to Galaxy
>    ```bash
>    cp -r topics/transcriptomics/tutorials/ref-based/tours/ $GALAXY_ROOT/config/plugins/tours
>    ```
>
{: .hands_on}




### Installing an entire topic

If you would like to install all the requirements for every tutorial within an entire topic, you can use the script in [`bin/install_topic_requirements.sh`]({{ site.github_repository }}/tree/master/bin/install_topic_requirements.sh)


### Installing a subset of components

If you would like to pick and choose what to install for each tutorial, below are descriptions of the commands used to install each of the components (tools, workflows, reference data, data libraries, tours) please see the [Quickstart section](#quickstart) for the individual commands used by the script

# Quickstart

Below is the list of commands used in this tutorial.

Using the scripts in this repository:

```bash
# Make sure you are in the root of the training-material repo
cd <training-materials repo root>

# install single tutorial
bin/install_tutorial_requirements.sh topics/<yourtopic>/tutorials/<yourtutorial> <Galaxy URL> <API key>

# install entire topic
bin/install_topic_requirements.sh topics/<yourtopic> <Galaxy URL> <API key>

```

Using ephemeris directly:

```
# install tools
shed-tools install -g <Galaxy URL> -a <API key> -t topics/<topic>/tutorials/<tutorial>/tools.yaml

# create data library with input datasets
setup-data-libraries -g <Galaxy URL> -a <API key> -i topics/<topic>/tutorials/<tutorial>/data-library.yaml

# install reference data
run-data-managers -g <Galaxy URL> -a <API key> --config topics/<topic>/tutorials/<tutorial>/data-manager.yaml

# install workflows
workflow-install --publish-workflows -g <Galaxy URL> -a <API key> -w topics/<topic>/tutorials/<tutorial>/workflows

# install tours
copy the contents of the "tours" directory for the tutorial to Galaxy's "config/plugins/tours/"
```

# Conclusion
{:.no_toc}
