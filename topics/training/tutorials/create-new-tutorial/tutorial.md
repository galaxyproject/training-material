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
> 3. [Checking the tutorial generation](#checking-the-tutorial-generation)
> 5. [Filling the tutorial content](#filling-the-tutorial-content)
> 6. [Adding slides](#adding-slides)
> 7. [Setting up the technical support](#setting-up-the-technical-support)
> 8. [Submitting the new tutorial to the GitHub repository](#submitting-the-new-tutorial-to-the-github-repository)
> {: .agenda}

# Setting up a new tutorial

Here, we want to develop a small tutorial to explain how to use BLAST.

## Clone the Galaxy Training material repository

Before anything, we need to get a local copy of the content of the GitHub repository by cloning it

> ### :pencil2: Hands-on: Clone the GitHub repository
>
> 1. Clone the repository locally with: `git clone https://github.com/galaxyproject/training-material.git`
> 2. Check that you have the same structure as the one in [GitHub](https://github.com/galaxyproject/training-material)
{: .hands_on}

## Defining the topic

The first step we need to define is in which topic putting our tutorial. This first step can be tricky. 

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

## Creating the directory for the tutorial

Once the topic is chosen, serious things can start: creating the tutorial. It is meaning the tutorial content, the metadata related to the tutorial but also the technical support for the tutorial with the description of the needed tool and dataset, a workflow of the tutorial and also a Galaxy Interactive Tour.

To help you, we created a template for a tutorial with the different required files.

> ### :pencil2: Hands-on: Copy the needed file
>
> 1. Copy the `tutorial1` folder (you can find it in `templates/tutorials/`) in `topics/sequence-analysis/topics`
> 2. Rename the folder into `similarity-search`
{: .hands_on}

## Keeping track of the changes

Once you started to change something, we need to keep track of these changes with a version control system (CVS). We are using Git as CVS and GitHub as hosting service.

This repository is developed collaboratively with more than 40 contributors. For the collaboration, we are using the [GitHub flow](https://guides.github.com/introduction/flow/) based on fork, branches and Pull Requests. We will explain the latter concept latter but now we will show how to start keeping track of the changes:

> ### :pencil2: Hands-on: Start keeping track of the changes
>
> 1. [Create a fork](https://help.github.com/articles/fork-a-repo/) of this repository on GitHub
> 2. Add your fork to the current local copy: `git remote add fork https://github.com/galaxyproject/training-material`
> 3. Create a new branch called "similarity-search" in your local copy: `git checkout -b similarity-search`
> 4. Commit the changes in that branch with
>     - `git add topics/sequence-analysis/tutorials/similarity-search`
>     - `git commit -m "Set up the similarity search tutorial"`
> 5. Push that branch to your fork on GitHub: `git push fork similarity-search`
{: .hands_on}

The GitHub interface can also help you in the process of editing a file. It will automatically create a fork of this repository where you can safely work.


We will now start to fill the different files together. We recommend you to commit regurlarly your changes. It help to follow them but also revert them if needed.

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

> ### :pencil2: Hands-on: Fill the pedagogical metadata
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

# Checking the tutorial generation

Once the metadata are filled, we can test if the metadata were correctly defined to serve the online website. Currently, the website is generated from the metadata and the tutorials using Jekyll, a templating system. We can run locally a Jekyll server to check if the tutorial is correctly added and rendered

> ### :pencil2: Hands-on: Checking the website generation locally
>
> 1. Install Jekyll using [RubyGems](https://rubygems.org/pages/download): `make install`
>
>    > ### :bulb: Tip: Installation on MacOS
>    >
>    > If you are installing it on Mac OSX, you need to install it this way as `/usr/bin/` is not writable. You need then to run:
>    >
>    >      $ sudo gem update —system
>    >      $ sudo gem install -n /usr/local/bin/ gem name
>    >      $ sudo gem install -n /usr/local/bin/ jemoji
>    >      $ sudo gem install -n /usr/local/bin/ jekyll
>    >      $ sudo gem install -n /usr/local/bin/ jekyll-feed
>    >      $ sudo gem install -n /usr/local/bin/ bundler
>    {: .tip}
>   
> 2. Run a local Jekyll server: `make serve`
> 3. Visualize at [http://localhost:4000/](http://localhost:4000/)
>
>    > ### :question: Questions
>    >
>    > Was the tutorial correctly added?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    If you check at [http://localhost:4000/topics/sequence-analysis/similarity-search](http://localhost:4000/topics/sequence-analysis/similarity-search) (the expected URL to the tutorial), we can not found the tutorial. The metadata at the top of the tutorial are not correctly defined. We will discuss about that latter.
>    >    </details>
>    {: .question}
{: .hands_on}

# Filling the tutorial content

## Finding a good toy dataset

## Filling the tutorial content

The content of the tutorial will go in the `tutorial.md`. The syntax of this file is really simple., as well as its structure:

```
---
layout: tutorial_hands_on
topic_name: training
tutorial_name: create-new-tutorial
---
# Introduction

blabla

# Section 1

blabla

## Subsection 1

blabla

# Section 2

blabla

## Subsection 2

blabla

# Conclusion

blabla
```

### Some metadata on the top

On the top, there is some metadata:

- `layout: tutorial_hands_on` (to let as it is)
- `topic_name: training` with the name of the topic
- `tutorial_name: create-new-tutorial` with the name of tutorial

These metadata are there to help the templating system to make the connection between the file and the metadata we defined before. If they are not correctly defined the tutorial can not be found on the website, as for our current tutorial.

> ### :pencil2: Hands-on: Fix the top metadata
>
> 1. Change the `tutorial-name` and the `topic_name` to fit to the ones defined in the metadata
> 2. Check if the tutorial has been correctly added at [http://localhost:4000/topics/sequence-analysis/similarity-search](http://localhost:4000/topics/sequence-analysis/similarity-search)
{: .hands_on}

### Content of the tutorial

After the metadata comes the content of the tutorial written in Markdow, a simple markup langage. 

> ### :bulb: Tip: Markdown
>
> Check [this cheatsheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet) to learn more how to use Markdown.
{: .tip}

The content in Markdown is transformed by our templating system into the nice webpage which add the metadata. Indeed in the `tutorial.md` file, no need to add the name of the tutorial: it is automatically added based on the title defined in the metadata.

We recommend to structure the tutorials like this

- An introdcution to introduce the tutorial with the use case, the data, the methods
- Several sections with the content of the tutorial and some hands-on parts (because we think that doing is an important part of the learning process)
- A conclusion to summarize what has been done in the tutorial (with a scheme)

> ### :pencil2: Hands-on: Structuring the tutorial
>
> 1. Add a small introduction about the dataset
> 2. Add one or two sections with ideas for the tutorial
> 3. Add a small conclusion
{: .hands_on}

### Improving the learning experience

To improve the learning experience in our tutorial, we defined some boxes to highlight.

They are defined always with the same structure:

```
> ### <an emoji> Type of boxe: Name of the box
> list
{: .type_of_box}
```

This structure needs to be respected otherwise it would not be interpreted correctly by the templating system. The different defined boxes are:

- Overview

    This box at the top of each tutorial is automatically generated using the metadata we defined

    > ### :pencil2: Hands-on: Checking the metadata
    >
    > 1. Check that the metadata added previously are correctly filling the overview box
    >
    >    > ### :question: Questions
    >    >
    >    > Which pedogical metadata are not added to this box?
    >    >
    >    >    <details>
    >    >    <summary>Click to view answers</summary>
    >    >    The take-home messages are not added to this box but into the last box of the tutorial
    >    >    </details>
    >    {: .question}
    {: .hands_on}

- Agenda

    The second box in most of the tutorial is the agenda box at the end of the introduction. It indicates the plan of the tutorial

        > ### Agenda
        >
        > In this tutorial, we will analyze the data with:
        >
        > 1. [Pretreatments](#pretreatments)
        > 2. [Mapping](#mapping)
        > 3. [Analysis of the differential expression](#analysis-of-the-differential-expression)
        > {: .agenda}

    ![](../../../../shared/images/tutorial_agenda_box.png)

    > ### :pencil2: Hands-on: Add an agenda box to the tutorial
    >
    > 1. Add an agenda box to the tutorial that fit the structure we defined previously
    {: .hands_on}

- Hands-on

    We think that doing is important in the learning process. So we emphasize it by adding regularly some hands-on sections where the trainees can do by themselves some analyses. We designed some special boxes to make these sections easy to find.

        > ### :pencil2: Hands-on: Sorting BAM dataset
        >
        > 1. **Sort BAM dataset** :wrench:: Sort the paired-end BAM file by "Read names" with **Sort BAM
        {: .hands_on}

    ![](../../../../shared/images/tutorial_hand_on_box.png)

    with the

    - `:pencil2:` emoji to define that is an hands-on
    - Short imperative sentence to make easy to identify the tasks
    - Name of the tool in bold to make it easy to identify and with `:wrench:` to insist
    - Parameters for the tool as a sublist<br/>

    > ### :pencil2: Hands-on: Add an hands-on box
    >
    > 1. Add an hands-on box to run a BLAST of the small sequence dataset against the chosen database 
    {: .hands_on}

-  Questions 
    
    The questions are then to force the trainees to think about what they are currently doing and to put things in perspective. They are also a way to help the instructors to expose and rectify in direct the misconceptions.

        > ### :question: Questions
        >
        > 1. Why are some tests filtered?
        > 2. Does it improve the *p*-value distribution?
        >
        >    <details>
        >    <summary>Click to view answers</summary>
        >    Content goes here.
        >    </details>
        {: .question}

    ![](../../../../shared/images/tutorial_question_box.png)

    They has to be quick to administer and evaluate. They can be small questions or also multiple choice questions (MCQs). If well designed with righlty chosen wrong asnwer, the latter solution can do much more than just measure how much someone knows by giving valuable insight.

    In the box, we add also the answer so the self-trainees can check the solution and its explanation. 

    > ### :pencil2: Hands-on: Add a question box
    >
    > 1. Add an hands-on box to construct the BLAST database
    {: .hands_on}

- Tips

    

        > ### :bulb: Tip: Importing data via links
        >
        > * Copy the link location
        > * Open the Galaxy Upload Manager
        > * Select **Paste/Fetch Data**
        > * Paste the link into the text field
        > * Press **Start**
        {: .tip}

    ![](../../../../shared/images/tutorial_tip_box.png)

- Comments

        > ### :nut_and_bolt: Comments
        > - Edit the "Database/Build" to select "dm3"
        > - Rename the datasets according to the samples
        {: .comment}

    ![](../../../../shared/images/tutorial_comment_box.png)

- Key points
    
    This last box of the tutorial is automatically filled with the take-home messages defined in the metadata


To render the boxes correctly, the previous syntaxes have to be followed. The boxes can be nested, *e.g.*for having tips inside hands-on. For example:

```
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
```

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
> 2. Add it to the `tours` directory
{: .hands_on}

## Testing the technical infrastructure

Once we defined all the requirements for running the tutorial, we can test these requirements. 

Every topic will come with a Docker image containing the tools, data, workflows and Galaxy Interactive Tours required by each tutorial of this topic. The Docker image is described in the Dockerfile found in the `docker` directory of each topic. This file uses scripts to automatically add the files for each tutorial. The only thing to change is the name of the topic in the Dockerfile copied from the templates.

> ### :pencil2: Hands-on: Testing the Docker
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

# Submitting the new tutorial to the GitHub repository

Once we are happy with your tutorial, we want to submit it to the GitHub repository. So we will open a Pull Request to propose our changes to the original GitHub repository

> ### :pencil2: Hands-on: Submitting a Pull Request with the new tutorial
>
> 1. Go on [GitHub](https://github.com/galaxyproject/training-material)
> 2. Create a new Pull Request and detail the proposed changes
> 3. If you receive feedback, make changes in your local clone and push them to your branch of your fork
>    The pull request will update automatically
{: .hands_on}

# Conclusion

Conclusion about the technical key points. And then relation between the technics and the biological question to end with a global view.
