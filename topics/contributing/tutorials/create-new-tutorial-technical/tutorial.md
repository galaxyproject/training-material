---
layout: tutorial_hands_on
topic_name: contributing
tutorial_name: create-new-tutorial-technical
---

# Introduction
{:.no_toc}

In this tutorial, you will learn how to create a virtualised Galaxy instance, based on Docker, to run your training - either on normal computers or cloud environments.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
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

This list of tools can be automatically extracted from the workflow using [Ephemeris](https://ephemeris.readthedocs.io/en/latest/index.html) (which should be in the conda environment):

```
$ workflow-to-tools -w path/to/workflow -o path/to/tools.yaml
```

After the extraction, some formatting is needed:

1. Add at the beginning:

    ```
    ---
    api_key: admin
    galaxy_instance: http://localhost:8080
    ```

2. Change the `tool_panel_section_label` to something more informative

> ### {% icon hands_on %} Hands-on: Fill the `tools.yaml` from your workflow
>
> 1. Fill the `tools.yaml` file using your workflow and Ephemeris
> 2. Format the `tools.yaml` file correctly
{: .hands_on}


## Filling the `data-library.yaml`

The data can also be integrated in the Galaxy instance inside a data libraries and then make the data shared between the users. It lets then avoid every trainees to re-download the input data.

Such data are described in the `data-library.yaml`:

```
destination:
  description: Galaxy Training Network Material
  name: GTN - Material
  synopsis: Galaxy Training Network Material. See https://training.galaxyproject.org
  type: library
items:
- name: <Topic title>
  description: <Topic summary>
  items:
  - url: <URL to one of the file on Zenodo>
    ext: <Format of the file>
    info: <DOI on Zenodo>
    src: url
    
```

with:

- `name`: name of the topic
- `items`: list of the input files to download and put on the data library
    - `url`: URL to the input file on [Zenodo](https://zenodo.org)
    - `ext`: format of the input file
    - `info`: URL to the DOI on Zendo

This file can be automatically generated from the metadata information:

```
$ python bin/create-data-library.py --tutorial <tutorial-name> --topic <topic-name>
```

> ### {% icon hands_on %} Hands-on: Fill the `data-library.yaml`
>
> 1. Add the input files into the `data-library.yaml` file
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

## Extracting workflows

Once the tutorial is ready, we need to extract workflows with the different steps of the tutorial and add them to the `workflows` directory in the tutorial with some explanation about the tutorial in a `README.md` file

> ### {% icon hands_on %} Hands-on: Extract the workflow
>
> 1. Extract the workflow for the tutorial
> 2. Add some description about the tutorial in a `README.md` file with the workflow file
{: .hands_on}

## Adding a Galaxy Interactive Tour

A Galaxy Interactive Tour is a way to go through an entire analysis, step by step inside Galaxy in an interactive and explorative way.
It is a great way to run the tutorial directly inside Galaxy. To learn more about creating a Galaxy tour please have a look at our [dedicated tour training]({{site.baseurl}}/topics/contributing/tutorials/create-new-tutorial-tours/tutorial.html).

## Testing the technical infrastructure

Once we have defined all the requirements for running the tutorial, we can test these requirements, either in a locally running Galaxy or in a Docker container. Please see our tutorial about [Setting up Galaxy for Training](../setup-galaxy-for-training/tutorial.html) about how to test your tutorial requirements.


# Conclusion
{:.no_toc}
