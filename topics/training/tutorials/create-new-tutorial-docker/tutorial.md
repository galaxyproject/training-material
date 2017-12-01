---
layout: tutorial_hands_on
topic_name: training
tutorial_name: create-new-tutorial-docker
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

In this tutorial, you will learn how to create a virtualised Galaxy instance, based on Docker, to run your training - either on normal computers or cloud environments.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

> ### Devloping GTN training material
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


# Building a Galaxy instance specifically for your training

To able to run the tutorial, we need a Galaxy instance where the needed tools are installed and the data. We need then to describe the needed technical infrastructure.

This description will be used to automatically set up a Docker Galaxy flavour and also to test if a public Galaxy instance is able to run the tool.

## Filling the `tools.yaml`

The first file to fill is the `tools.yaml` file, containing the description of the required tools that could be installed from the ToolShed.

This file looks like:

```
---
api_key: admin
galaxy_instance: http://localhost:8080
tools:
- name: tool1
  owner: owner
  tool_panel_section_label: "Section1"
- name: tool2
  owner: owner
  tool_panel_section_label: "Section2"
```

with:

- `name`: the name of the wrapper of the tool in the ToolShed
- `owner`: the owner of the wrapper of the tool in the ToolShed
- `tool_panel_section_label`: section where to put the tool (in the left panel in the Galaxy instance)

> ### {% icon hands_on %} Hands-on: Fill the `tools.yaml`
>
> 1. Add the BLAST tool into the `tools.yaml` file
{: .hands_on}

## Filling the `data-library.yaml`

The data can also be integrated in the Galaxy instance inside a data libraries and then make the data shared between the users. It lets then avoid every trainees to redownload the input data.

Such data are described in the `data-library.yaml`:

```
libraries:
    - name: Name of the tutorial
      files:
        - url: "https://raw.githubusercontent.com/bgruening/galaxytools/master/tools/rna_tools/sortmerna/test-data/read_small.fasta"
          file_type: fasta
        - url: ""
          file_type: ""
```

with:

- `name`: name of the tutorial, where to put the data in the data libraries
- `files`: list of the files to download
    - `url`: URL to the input file
    - `file-type`: type of the input file

The URL must refer to the URL of the files in Zenodo.

> ### {% icon hands_on %} Hands-on: Fill the `data-library.yaml`
>
> 1. Add the input files into the `data-library.yaml` file
> 2. Add the link to Zenodo in the `metadata.yaml` file
{: .hands_on}

## Filling the `data-manager.yaml`

Some of the tools require specific databases, specifically prepared for the tool. Then some Galaxy tools come with data managers to manage these databases.

If you need such data managers for your tool, you can describe their running with the `data-manager.yaml` file:

```
data_managers:
    - id: url to data manager on ToolShed
      params:
        - 'param1': '{{ item }}'
        - 'param2': 'value'
      # Items refere to a list of variables you want to run this data manager. You can use them inside the param field with {{ item }}
      # In case of genome for example you can run this DM with multiple genomes, or you could give multiple URLs.
      items:
        - item1
        - item2
      # Name of the data-tables you want to reload after your DM are finished. This can be important for subsequent data managers
      data_table_reload:
        - all_fasta
        - __dbkeys__
```

## Extracting workflows

Once the tutorial is ready, we need to extract workflows with the different steps of the tutorial and add them to the `workflows` directory in the tutorial with some explanation about the tutorial in a `README.md` file

> ### {% icon hands_on %} Hands-on: Extract the workflow
>
> 1. Extract the workflow for the tutorial
> 2. Add some description about the tutorial in a `README.md` file with the workflow file
{: .hands_on}

## Adding a Galaxy Interactive Tour

A Galaxy Interactive Tour is a way to go through an entire analysis, step by step inside Galaxy in an interactive and explorative way.
It is a great way to run the tutorial directly inside Galaxy. To learn more about creating a Galaxy tour please have a look at our [dedicated tour training]({{site.baseurl}}/topics/training/tutorials/create-new-tutorial-tours/tutorial.html).

## Testing the technical infrastructure

Once we defined all the requirements for running the tutorial, we can test these requirements.

Every topic will come with a Docker image containing the tools, data, workflows and Galaxy Interactive Tours required by each tutorial of this topic. The Docker image is described in the Dockerfile found in the `docker` directory of each topic. This file uses scripts to automatically add the files for each tutorial. The only thing to change is the name of the topic in the Dockerfile copied from the templates.

> ### {% icon hands_on %} Hands-on: Testing the Docker
>
> 1. Check that the Dockerfile uses 'sequence-analysis' as topic name
> 2. Move to the root of the training material repository
> 3. Build the Docker image for the topic with: `docker build -f topic/sequence-analysis/docker/Dockerfile -t training-sequence-analysis .`
>     
>    This command needs to be launched a the root of training material repository because the Dockerfile uses some scripts available there to install the tools, import the data and the workflows
>
> 4. Launch the Docker container: `docker run -d -p 8080:80 training-sequence-analysis`
> 5. Check the Galaxy instance on [http://localhost:8080/](http://localhost:8080/):
>     1. Check the installed tools
>     2. Check the data libraries in "Shared data"
>     3. Check the workflows
>     4. Check the Galaxy Interactive Tours in "Help"
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
