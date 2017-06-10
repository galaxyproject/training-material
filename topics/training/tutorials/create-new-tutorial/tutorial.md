---
layout: tutorial_hands_on
topic_name: training
tutorial_name: create-new-tutorial
---

# Introduction

Galaxy is a great solution to train the bioinformatics concepts:

- numerous bioinformatics tools are available (almost 5,000 in the ToolShed)
- it can be used by people without amy computer science skills
- it trains to use technology, outlining available resources and efforts that have made them accessible to researchers
- it is scalable

In 2016, the Galaxy Training Network decide to set up a new infrastructure for delivering easily Galaxy related training material. The idea was to develop something open and online based on a community effort, as most of the time in Galaxy. 

We take inspiration from [Software Carpentry](https://software-carpentry.org). We collected everything on a GitHub repository: [https://github.com/galaxyproject/training-material](https://github.com/galaxyproject/training-material). We decided a structure based on tutorials with hands-on, fitting both for online self-training but also for workshops, grouped in topics. Each tutorial follows the same structure and comes with a technical support to be able to run. 

In this tutorial, you will understand how to design and develop a new tutorial fitting in this training material repository. As doing helps to understand, we will develop a small tutorial to explain BLAST with the full infrastructure to be able to run this tutorial anywhere. 

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. [Setting up a new tutorial](#setting-up-a-new-tutorial)
> 2. [Filling the metadata](#filling-the-metadata)
> 3. [Analysis of the differential expression](#analysis-of-the-differential-expression)
> {: .agenda}

# Setting up a new tutorial

Here, we want to develop a small tutorial to explain how to use BLAST. The first step we need to define is in which topic putting our tutorial. This first step can be tricky. 

When we structured the repository, we decided here to use as topic the names of the categories in the [ToolShed](https://toolshed.g2.bx.psu.edu/). So when decided where to put your tutorial, you can look in which ToolShed's category are the main tools used in the tutorial and use this category as topic. For example, this tutorial will rely on the NCBI Blast+ tool.

> ### :pencil2: Hands-on: Defining the topic for the tutorial
>
> 1. Search for NCBI Blast+ on the [ToolShed](https://toolshed.g2.bx.psu.edu/)
> 2. Check in which category it is
>
>    > ### :question: Questions
>    >
>    > In which topic will you put the tutorial?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    If we search for [NCBI Blast+ in the ToolShed](https://toolshed.g2.bx.psu.edu/view/devteam/ncbi_blast_plus/7538e2bfcd41), it is attributed to 2 categories (bottom): "Next Gen Mappers" and "Sequence Analysis".
>    >    We decided to put it in "Sequence analysis" because this is the most general one for this tutorial.
>    >    </details>
>    {: .question}
{: .hands_on}

Once the topic is chosen, serious things can start: creating the tutorial. It is meaning the tutorial content, the metadata related to the tutorial but also the technical support for the tutorial with the description of the needed tool and dataset, a workflow of the tutorial and also a Galaxy Interactive Tour.

To help you, we created a template for a tutorial with the different required files.

> ### :pencil2: Hands-on: Copy the needed file
>
> 1. Copy the `tutorial1` folder (you can find it in `templates/tutorials/`) in `topics/sequence-analysis/topics`
> 2. Rename the folder into `similarity-search`
{: .hands_on}

We will now start to fill the different files together.

# Filling the metadata

The first file we will fill is the `metadata.yaml` file. 

This file define the metadata related to a tutorial: 

- `title`: title of the tutorial
- `type: "tutorial"`
- `name`: name of the tutorial (name of the subdirectory where the files related to the tutorial will be stored)
- `zenodo_link`: link on Zenodo to the input data for the tutorial (not ideal but it can be empty)
- `galaxy_tour`: name of the galaxy tour
- `hands_on`(`"yes"` or `"no"`): tell if an hands on is available for this material
- `slides` (`"yes"` or `"no"`): tell if slides are available for this materialits title, its type, ...
- `requirements`: list of requirements specific to this tutorial (in addition to the one of the topic), with:
    - `title`
    - `link`: relative for internal (inside training material) requirement or full for external requirement)
    - `type`: the type of link (`internal` or `external`)

This information is used to automatically make the tutorial available on the online website: [http://galaxyproject.github.io/training-material/](http://galaxyproject.github.io/training-material/)

> ### :pencil2: Hands-on: Fill the basic metadata
>
> 1. Fill the basic metadata for our tutorial
>   - `title: Similarity search with BLAST`
>   - `type: "tutorial"`
>   - `name: "similarity-search"`
>   - `zenodo_link: ""` (we do not have data currently)
>   - `galaxy_tour: ""` (we do not have Galaxy Interactive Tour currently)
>   - `hands_on: "yes"`
>   - `slides: "no"`
>   - `requirements`: a requirement to "Galaxy introduction" with internal link to `introduction`
{: .hands_on}

In the second part of the metadata, we define metadata related to the content of the tutorial, which will appear in the top and bottom of the online tutorial:

- `time_estimation`: an estimation of the time needed to complete the hands-on
- `questions`: list of questions that will be addressed in the tutorial
- `objectives`: list of learning objectives of the tutorial

    A learning objective is a single sentence describing what a learner will be able to do once they have deone the tutorial

- `key_points`: list of take-home messages

    This information will appear at the end of the tutorial 

For this metadata, we take inspiration from what Software Carpentry is doing and particularly what they describe in their [Instructor training](http://swcarpentry.github.io/instructor-training/) and the section ["Lessons and Objectives"](http://swcarpentry.github.io/instructor-training/19-lessons/). 

> ### :pencil2: Hands-on: Fill the basic metadata
>
> 1. Define 2 questions that will be addressed during the tutorial and add them to the metadata
> 2. Define 2 learning objectives for the tutorial and add them to the metadata
{: .hands_on}

We recommend you to fill the questions and the learning objectives before starting writing the tutorial content. You can still refine them afterwards, but it will help to design your tutorial and think beforehands what is worth training.

For the take-home messages, it is easier to define them once the tutorial is written and you identified the issues.

> ### :nut_and_bolt: Comment
>
> Temporarly, we need a duplication of the metadata (because of our temporary templating system) in the topic `metadata.yaml` in the section `material`. So you need to copy that.
> {: .comment}

Once the metadata are filled, we can test if the metadata were correctly defined to serve the online website. Currently, the website is generated from the metadata and the tutorials using Jekyll, a templating system. We can run locally a Jekyll server to check if the tutorial is correctly added and rendered

> ### :pencil2: Hands-on: Checking the website generation
>
> 1. Install Jekyll using [RubyGems](https://rubygems.org/pages/download): `make install`
>
>    > ### :bulb: Tip: Installation on MacOS
>    >
>    > If you are installing it on Mac OSX, you need to install it this way as `/usr/bin/` is not writable. You need then to run:
>    >
>    > - `$ sudo gem update —system`
>    > - `$ sudo gem install -n /usr/local/bin/ gem name`
>    > - `$ sudo gem install -n /usr/local/bin/ jemoji`
>    > - `$ sudo gem install -n /usr/local/bin/ jekyll`
>    > - `$ sudo gem install -n /usr/local/bin/ jekyll-feed`
>    > - `$ sudo gem install -n /usr/local/bin/ bundler`
>    {: .tip}
>   
> 2. Run a local Jekyll server: `make serve`
> 3. Visualize at [http://localhost:4000/](http://localhost:4000/)
> 4. Check if the tutorial was correctly added
{: .hands_on}

# Filling the tutorial

## Finding a good toy dataset

## Filling the tutorial content

Short introduction about this subpart.

> ### :pencil2: Hands-on: Data upload
>
> 1. Step1
> 2. Step2
>
>    > ### :question: Question
>    >
>    > Question?
>    >
>    > <details>
>    > <summary>Click to view answers</summary>
>    > Answer to question
>    > </details>
>    {: .question}
{: .hands_on}

Some blabla
> ### :pencil2: Hands-on: Data upload
>
> 1. Step1
> 2. Step2
>
>    > ### :question: Questions
>    >
>    > 1. Question1?
>    > 2. Question2?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>Answer for question1</li>
>    >    <li>Answer for question2</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 3. Step3
{: .hands_on}

# Adding slides

(If needed)

# Setting up the technical support

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

> ### :pencil2: Hands-on: Fill the `tools.yaml`
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
        - url: "http://raw.githubusercontent.com/bgruening/galaxytools/master/tools/rna_tools/sortmerna/test-data/read_small.fasta"
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

> ### :pencil2: Hands-on: Fill the `data-library.yaml`
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

> ### :pencil2: Hands-on: Extract the workflow 
>
> 1. Extract the workflow for the tutorial
> 2. Add some description about the tutorial in a `README.md` file with the workflow file
{: .hands_on}

## Creating a Galaxy Interactive Tour

A Galaxy Interactive Tour is a way to go through an entire analysis, step by step inside Galaxy in an interactive and explorative way. It is a great pedogogic way to run the tutorial directly inside Galaxy

> ### :pencil2: Hands-on: Create a Galaxy Interactive Tour
>
> 1. Create a Galaxy Interactive Tour for the tutorial
{: .hands_on}

# Submitting the new tutorial

# Conclusion

Conclusion about the technical key points. And then relation between the technics and the biological question to end with a global view.
