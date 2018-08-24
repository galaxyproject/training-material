---
layout: tutorial_hands_on

title: "Creating a new tutorial - Defining the technical infrastructure"
questions:
  - "How can we define the technical infrastructure for a tutorial?"
  - "How to define the tools needed for a tutorial?"
  - "How to add the needed data directly in an instance?"
  - "How to add the workflows related to a tutorial?"
  - "How can we check the technical infrastructure is working?"
  - "How can we make an existing Galaxy instance able to run a tutorial?"
objectives:
  - "Extracting the technical description for a tutorial"
  - "Populating an existing instance with the needed tools, data and workflows for a tutorial"
  - "Creating a Galaxy Docker flavor with the needed tools, data and workflows for a tutorial"
  - "Testing the Galaxy Docker flavor of a tutorial"
time_estimation: "30min"
key_points:
  - "Tools, data and workflows can be easily integrated in a Docker flavor to have a useful technical support for a tutorial"
  - "A Galaxy Docker flavor is a great support for training"
  - "A Galaxy Docker flavor can be deployed 'anywhere' and is scalable"
contributors:
  - bebatut
  - bgruening
  - shiltemann
---

# Building a Galaxy instance specifically for your training
{:.no_toc}

To be able to run the tutorial, we need a Galaxy instance where all of the needed tools and data are available. Thus we need to describe the needed technical infrastructure.

This files we define in this tutorial will be used to automatically build a Docker Galaxy flavour, and also to test if a public Galaxy instance is able to run the tool.

In this tutorial, you will learn how to create a virtualised Galaxy instance, based on Docker, to run your training - either on normal computers or cloud environments.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Extracting workflows

Once the tutorial is ready, we need to develop a workflow that represents the steps taken in the tutorial, and then extract these workflow(s) and add them to the `workflows` directory in the tutorial. Additionally we will need to add some explanation about the workflow(s) in a `README.md` file

> ### {% icon hands_on %} Hands-on: Extract the workflow
>
> 1. Download the workflow for the tutorial
> 2. Save it in the `workflow` directory of the tutorial
> 3. Set `workflow` to `yes` in the appropriate section of the topic's `metadata.yaml`.
{: .hands_on}

# Creating the `tools.yaml` (recommended)

The first file to fill out is the `tools.yaml` file which contains the list of the required tools that could be installed from the ToolShed.

This file looks like:

```yaml
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

This list of tools can be automatically extracted from the workflow using [Ephemeris](https://ephemeris.readthedocs.io/en/latest/index.html) (which should be in the conda environment):

```console
$ workflow-to-tools -w path/to/workflow -o path/to/tools.yaml
```

After the extraction, some formatting is needed:

1. Add at the beginning:

    ```yaml
    ---
    api_key: admin
    galaxy_instance: http://localhost:8080
    ```

2. Change the `tool_panel_section_label` to something more informative

> ### {% icon hands_on %} Hands-on: Creating the `tools.yaml` from your workflow
>
> 1. Create the `tools.yaml` file using your workflow and Ephemeris
> 2. Correct the formatting of the `tools.yaml` file
{: .hands_on}


# Creating the `data-library.yaml` (recommended)

The datasets needed for a tutorial can also be integrated in the Galaxy instance inside of data libraries. These allow the datasets to be easily shared with all users of a Galaxy instance. Additionally it lets trainees avoid each re-downloading the input data.

These datasets are described in the `data-library.yaml` files:

```yaml
libraries:
    - name: Name of the tutorial
      files:
        - url: "https://raw.githubusercontent.com/bgruening/galaxytools/master/tools/rna_tools/sortmerna/test-data/read_small.fasta"
          file_type: fasta
        - url: ""
          file_type: ""
```

where:

- `name`: name of the tutorial, i.e. where to put the data in the data libraries
- `files`: list of the files to download
    - `url`: URL to the input file
    - `file-type`: type of the input file

The URL must refer to the URL of the files in [Zenodo](https://zenodo.org), do not just link to the overview page.

> ### {% icon hands_on %} Hands-on: Creating the `data-library.yaml`
>
> 1. Add the input files from your Zenodo dataset to the `data-library.yaml` file
> 2. Add the link to Zenodo in the `metadata.yaml` file
{: .hands_on}

# Creating the `data-manager.yaml` (optional)

Some of the tools may require specific databases, specifically prepared for the tool. In this case, some Galaxy tools come with "data managers" to simplify this process.

If you need such data managers for your training, then you should describe how to run them in the `data-manager.yaml` file:

```yaml
data_managers:
    - id: url to data manager on ToolShed
      params:
        - 'param1': '{{ item }}'
        - 'param2': 'value'
      # Items refer to a list of variables you want to run this data manager. You can use them inside the param field with {{ item }}
      # In case of genome for example you can run this DM with multiple genomes, or you could give multiple URLs.
      items:
        - item1
        - item2
      # Name of the data-tables you want to reload after your DM are finished. This can be important for subsequent data managers
      data_table_reload:
        - all_fasta
        - __dbkeys__
```

# Creating the Galaxy Interactive Tour (optional)

A Galaxy Interactive Tour is a way to go through an entire analysis, step by step inside Galaxy in an interactive and explorative way.
It is a great way to help users run the tutorial directly inside Galaxy. To learn more about creating a Galaxy tour please have a look at our [dedicated tour training]({{site.baseurl}}/topics/contributing/tutorials/create-new-tutorial-tours/tutorial.html).

# Testing the technical infrastructure

Once we have defined all the requirements for running the tutorial, we can test these requirements, either in a locally running Galaxy or in a Docker container. Please see our tutorial about [Setting up Galaxy for Training]({{site.baseurl}}/topics/contributing/tutorials/setup-galaxy-for-training/tutorial.html) about how to test your tutorial requirements.


# Conclusion
{:.no_toc}
