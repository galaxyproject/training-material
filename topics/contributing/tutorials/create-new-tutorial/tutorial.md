---
layout: tutorial_hands_on

title: "Creating a new tutorial"
questions:
  - "How to integrate a new tutorial?"
  - "How to make a tutorial robust and reproducible?"
objectives:
  - "Create a tutorial from scratch"
  - "Link a tutorial to a topic"
  - "Create hands-on"
  - "Add technical support for a tutorial"
time_estimation: "15m"
key_points:
  - "Finding good training datasets is hard!"
  - "Creating a new tutorial involves several steps: some are mandatory, some can be skipped even if they are recommended"
subtopic: writing
contributors:
  - bebatut
  - hexylena
  - shiltemann
  - lldelisle
---

# Introduction


Galaxy is a great solution to train bioinformatics concepts:

- numerous bioinformatics tools are available (over 8,000 in the ToolShed)
- it can be used by people without any computer science skills
- it trains to use technology, outlining available resources and efforts that have made them accessible to researchers
- it is scalable

In 2016, the Galaxy Training Network decide to set up a new infrastructure for easily delivering Galaxy related training material. The idea was to develop something open, online, based on a community effort, and on top of the Galaxy platform.

We took inspiration from [Software Carpentry](https://software-carpentry.org) and collected everything on a GitHub repository: [https://github.com/galaxyproject/training-material](https://github.com/galaxyproject/training-material).
We decided on a structure focusing on tutorials with hands-on activities; fitting both for online self-training but also for workshops. Each tutorial follows the same structure and comes with a virtualised instance allowing you to run the training anywhere you have resources available.

Here you will learn how to create a new tutorial by developing a small tutorial explaining how to retrieve climate data from Copernicus (using **Copernicus Climate Data Store** tool).

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

> <comment-title></comment-title>
> This tutorial explains the different steps to create a tutorial for the Galaxy Training Material.
> It may require some knowledge that you may not have or do not have the time to learn. If this is the case, you can create a skeleton of a tutorial with whatever existing materials you have, using your prefered text editor, and then share it with us by opening [issue on GitHub]({{ site.github_repository }}/issues/new), writing us on [Gitter]({{ site.gitter_url }}), or sending us an [email](mailto:{{ site.email }}).
{: .comment}

# Define the topic

The first step we need to do is to identify in which topic to place our new tutorial. This can be tricky. When we structured the GTN, we decided that each training material should be related to a topic.

We decided to use the categories from the [ToolShed](https://toolshed.g2.bx.psu.edu/) as our initial list of topics. Since every tool uploaded to the ToolShed must be in at least one category, you can look at the main tools in your tutorial and see which categories they are placed in within the ToolShed. This can provide a guide for where you might put your new tutorial. For example, if your tutorial will rely on the **Copernicus** tool:

> <hands-on-title>Defining the topic for the tutorial</hands-on-title>
>
> 1. Search for **Copernicus** on the [ToolShed](https://toolshed.g2.bx.psu.edu/)
> 2. Check in which category it has been placed
>
> > <solution-title></solution-title>
> > There are a couple steps to reaching the answer:
> > 1. Search for `Copernicus`
> > 2. Press the <kbd>Enter</kbd> key to search
> > 3. Click on the result named `c3s`
> > 4. At the bottom of this page there is a box labelled "Categories"
> >
> > It is placed in one category: "Climate Analysis"
> {: .solution}
>
> 3. Look if the categories fit to existing topics in the Galaxy Training Material
{: .hands_on}

> <comment-title>No fitting topic for the tools in your tutorial?</comment-title>
> If the categories in ToolShed do not fit to any existing topics, we recommend to use your better judgment to identify in which topic your tutorial should go. You can also ask us on [Gitter]({{ site.gitter_url }}) or raise an issue on [GitHub]({{ site.github_repository }}) explaining the aim of the tutorial. We will be happy to help you there.
{: .comment}

> <comment-title>Creating a new topic</comment-title>
> Want to create a new topic? [Check out our tutorial to create a new topic]({% link topics/contributing/tutorials/create-new-topic/tutorial.md %})
{: .comment}

# Store a tutorial

All training materials (slides, tutorials, ...) related to a topic are found in a dedicated directory (*e.g.* `climate` directory contains the material related to Climate analyses). Each topic has the following structure:

```
├── README.md
├── metadata.yaml
├── images
├── docker
│   ├── Dockerfile
├── slides
│   ├── index.html
├── tutorials
│   ├── tutorial1
│   │   ├── tutorial.md
│   │   ├── slides.html
│   │   ├── data-library.yaml
│   │   ├── workflows
│   │   │   ├── workflow.ga
│   │   ├── tours
│   │   │   ├── tour.yaml
```

Once the topic has been chosen, we can create the tutorial. An ideal tutorial in the Galaxy Training Network contains:
- a tutorial file `tutorial.md` written in Markdown with hands-on
- an optional slides file `slides.html` in Markdown with slides to support the tutorial
- a directory `workflows` with workflows extracted from the tutorial
- a YAML file `data-library.yaml`  with the links to the input data needed for the tutorial

The most important file is the `tutorial.md` where the content of the tutorial is. The other files are not mandatory. They are there to support the tutorial and make it robust and usable across many environments. But they can help a lot in the development of new tutorial.

# Create the workflow

Our tutorials try to follow the "learn by doing" approach; they combine both theoretical and practical sections done on Galaxy.

Most tutorials explain step-by-step a data analysis by running the tools. The steps taken in the tutorial can be represented as a workflow.

Before writing the tutorial, it is a good practice to get a workflow with the different steps that will be run during the tutorial. The workflow does not have to be the final one but at least the major steps. It helps to get a direction for the tutorial but also to generate a skeleton of the tutorial as we will see later.

> <comment-title>Use tools that are available on the ToolShed</comment-title>
> We recommmend you to use in your workflows, specially for training, tools that are available on the Galaxy ToolShed.
{: .comment}


> <hands-on-title>Prepare the workflow for the training</hands-on-title>
>
> 1. Go to your favorite Galaxy server
> 2. Create a workflow with the different steps (tools) of your tutorial either from scratch or from an existing history
>
>    {% snippet faqs/galaxy/workflows_create_new.md %}
>
>    {% snippet faqs/galaxy/workflows_extract_from_history.md %}
>
>    For the Copernicus tutorial:
>    1. Add an **Input Dataset**
>    2. Add **Copernicus Climate Data Store**
>    3. Link the **Input Dataset** to **Copernicus Climate Data Store** on `API Request filename`
>    4. Rename **Input Dataset** to `API Request file`
>    5. Save the workflow
>
> 3. Add the topic name as a Tag and the tutorial title as Annotation/Notes to the workflow using the workflow editor
>
>    {% snippet faqs/galaxy/workflows_annotate.md %}
>
>    For the Copernicus tutorial:
>    - **Annotation**: `Retrieve climate data from Copernicus`
>    - **Tag**: `climate`
>
> 4. Make the workflow accessible (publishing is not necessary)
>
>    {% snippet faqs/galaxy/workflows_publish.md %}
>
{: .hands_on}


# Get a toy dataset

To run the different steps, the tutorial needs some data. The selected data must be informative enough to illustrate the meaning of using a tool or a given technique, but not too big to require long waiting times for processing during a workshop. Upload and download of files into and out of Galaxy is usually quick, but the time taken for a tool to run can be long. Tool run times of no more than 10-15 mins are recommended. Typically, the selected data should be the informative subset of a full real-life dataset.

Below we describe two examples of how toy datasets were generated for tutorials:

- **Example 1**: creating a toy dataset from scratch
  - Take one 16S sequence (for example found in the test case of a Galaxy tool):
  - Generate a reference database
      - Blast it on the NR database on [NCBI Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome)
      - Extracting one similar sequence found with Blast
      - Search and extract 2 other sequences of the same species using the [NCBI Nucleotide database](https://www.ncbi.nlm.nih.gov/nuccore)

- **Example 2**: creating a toy dataset from an existing larger one
  - When the experiment takes a FASTQ as input and a few reads are sufficient:
    - Use **seqtk_sample** {% icon tool %} to extract randomly reads from your input fastq.
  - However, when it requires a lot of reads to be meaningful, you can use the following strategy (used for the ATAC-seq tutorial using [this workflow](./workflows/Galaxy-Workflow-MakeAFakeInput.ga)):
    - Run the workflow until the mapping step on the full dataset (or big enough to have good results).
    - Select IDs of reads which map on the smallest chromosome (for example chr22 for human data).
    - In order to keep in the toy dataset enough diversity, you can also take randomly 1% of the reads IDs.
    - Concatenate the two lists and remove the duplicated IDs.
    - Use **seqtk_subseq** {% icon tool %} to sample your original FASTQ with the list of IDs.

We would then develop the tutorial and test it on this toy dataset. Once we were ready to share it, we would upload the datasets to [Zenodo](https://zenodo.org/) to store them on long-term and obtain a dedicated DOI in the [Galaxy training network community](https://zenodo.org/communities/galaxy-training/?page=1&size=20).

> <hands-on-title>Upload the dataset to Zenodo</hands-on-title>
>
> 1. Go to [Zenodo](https://zenodo.org/)
> 2. Log in using your GitHub credentials
>
>    You may need to authorize Zenodo to access your GitHub account (only to read your information)
>
> 3. Click on **Upload** (top panel)
> 4. Start a new upload
> 5. Upload the files corresponding to your datasets
>
>     > <comment-title>No possible changes in the files after publication</comment-title>
>     > File addition, removal or modification are not allowed after you have published your upload.
>     > So be careful when you start your upload that all your needed files are ready.
>     >
>     > The metadata can be changed after publication.
>     {: .comment}
>
> 6. Search for and Select *Galaxy Training Network* in **Communities**
> 7. Select *Dataset* in **Upload type**
> 8. Use the title of your tutorial and mention also Galaxy Training Material
> 9. Add all the persons who contributed to the tutorial as authors
> 10. Add a short description of the tutorial and a link to the training material website
> 11. Keep *Open Access* as **Access right** and *Creative Commons Attribution 4.0* as **License**
> 12. Fill out any remaining information
> 13. Click on **Save**
> 14. Click on **Publish**
> 15. Copy the DOI link in the new page
{: .hands_on}

# Create the skeleton of the tutorial

Once we have the workflows and toy dataset (optionally already on Zenodo), we can generate the tutorial folder, including the `tutorial.md`.

Writing the tutorial while adding the different tools and their parameters and following good formatting can be quite cumbersome. To help in that process, we developed a training development kit available with [Planemo](https://planemo.readthedocs.io/en/latest/index.html). The tool has several commands. But the main one can be used to create the skeleton of a tutorial (folder, well located, with good structure). It can also take a workflow as input, add it to the `workflows` folder. But more interestingly, the tool can extract from the workflow and add in the `tutorial.md` file the different steps with which tools to run and which parameters to select. If a Zenodo URL is provided, it also creates the `data-library.yaml` file.

This tool can be used via the command-line but also via a [webserver](https://ptdk.apps.galaxyproject.eu/). The webserver can take a public workflow that is avialable on [usegalaxy.org](https://usegalaxy.org), [usegalaxy.eu](https://usegalaxy.eu) or [usegalaxy.org.au](https://usegalaxy.org.au).

> <hands-on-title>Create the skeleton of a tutorial using the webserver</hands-on-title>
>
> 1. Make the workflow public
>
>    {% snippet faqs/galaxy/workflows_publish.md %}
>
> 2. Copy the workflow id that can be found in the URL of the current page (after `?id=`)
> 3. Open the [PTDK webserver](https://ptdk.apps.galaxyproject.eu/)
> 4. Fill in the information
>    - Tutorial name (the name will be the name of the folder of the tutorial)
>    - Tutorial title
>    - Galaxy instance with the public workflow
>    - Id of the workflow
>    - (Not mandatory) Zenodo URL with the input data
>
>    For the Copernicus tutorial:
>    - **Tutorial name** (the name will be the name of the folder of the tutorial): `climate-data-retrieval-copernicus`
>    - **Tutorial title**: `Retrieve climate data from Copernicus`
>    - **Galaxy instance with the public workflow**: `usegalaxy.eu`
>    - **Id of the workflow**: `ac5b66c42681e7a8`
>
> 5. Click on **Submit**
> 6. Download the generated archive
>
>    This archive contains the tutorial skeleton including:
>    - tutorial content, `tutorial.md` file, filled with all steps from the workflow and their parameters
>    - its workflow (`workflow` folder)
>    - a `data_library.yaml` file if Zenodo link was provided
>
> 7. Add the new material to Galaxy Training Material by unzip the downloaded archive in the tutorials folder of the topic for the new tutorial
>
>    > <comment-title>Using the GitHub interface</comment-title>
>    >
>    > Prefer to use the GitHub interface?
>    >
>    > 1. Unzip the downloaded archive
>    > 2. Edit the content of the `tutorial.md` (as explained below)
>    > 3. Go to the GitHub repository of the Training Material
>    > 4. Fork the GitHub repository
>    > 5. Click on `topics`
>    > 6. Select the topic for the new tutorial
>    > 7. Go to `tutorials`
>    > 8. Click on **Create new file**
>    > 9. Type `name/tutorial.md`, replacing "name" by the name of your tutorial (not the title)
>    > 10. Copy the content of downloaded and edited `tutorial.md` file there
>    > 11. Fill the **Commit new file** form
>    > 12. Create a new branch using the name of the tutorial
>    > 13. Click on **Propose new file**
>    > 14. Open a Pull Request, as explained in our tutorial
>    > 15. Add the workflow file (in `workflow` folder) by updating the Pull Request
>    {: .comment}
>
{: .hands_on}

If the workflow is not available on one of the previously listed Galaxy servers, we recommend you to run the tool via the command line:

> <hands-on-title>Create the skeleton of a tutorial via the command line</hands-on-title>
>
> 1. Get the workflow id
>
>    The id can be found on URL when running, editing or sharing the workflow (after `?id=`)
>
> 2. Get your API key on the Galaxy instance
>
>    {% snippet faqs/galaxy/preferences_admin_api_key.md %}
>
> 3. (If not done yet) Get the Galaxy Training Material repository locally and move in it
>
>     1. (If not done yet) Clone the training material GitHub repository: `git clone https://github.com/galaxyproject/training-material.git`
>     2. Navigate to the `training-material/` folder with `cd`
>
> 4. Get planemo
>
>     - Option 1: Using Conda:
>       1. Set up the conda environment
>
>          It will install some needed tools (ruby, nodejs, etc) in a protected environment, without interfering with the existing tools or versions.
>
>          1. Install conda (if not already installed): `make install-conda`
>          2. (You may need to exit the terminal and re-open for conda to be recognised. Navigate back to the same place.)
>          3. Create the `galaxy_training_material` conda environment: `make create-env`
>
>       2. Activate the conda environment with `conda activate galaxy_training_material`
>
>     -  Option 2: Using pip by running `pip install planemo`
>
> 6. Generate the skeleton of your tutorial (by adapting the information between the quotes)
>
>      ```
>      $ planemo training_init \
>             --topic_name "my-topic" \
>             --tutorial_name "my-new-tutorial" \
>             --tutorial_title "Title of the tutorial" \
>             --galaxy_url "URL to Galaxy instance in which you created the workflow" \
>             --galaxy_api_key "Your API key on the Galaxy instance" \
>             --workflow_id "ID of the workflow on the Galaxy instance" \
>             --zenodo_link "URL to the Zenodo record"
>      ```
>
> 4. Check that a new directory (with your tutorial name) has been generated in the topic folder
{: .hands_on}

# Write the tutorial

Now that you have the structure in place, you can then fill the tutorial per se.

> <hands-on-title>Write the tutorial</hands-on-title>
>
> 1. Open the `tutorial.md` file with your favorite text editor
> 2. Fill out the tutorial by following the [dedicated tutorial]({% link topics/contributing/tutorials/create-new-tutorial-content/tutorial.md %})
>
>    1. Add metadata on the top of the tutorial
>    2. Add a proper introduction
>    3. Organize the tutorial in different sections
>    4. Introduce the different concepts and the different steps of the analysis
>    5. Check the hands-on boxes
>    6. Comment on the outputs of the different steps
>    7. Add some question/solution boxes for assessment
>    8. Add a conclusion
>
> 3. (Optional) Build the website locally and check that the tutorial is there by following the [Jekyll tutorial]({% link topics/contributing/tutorials/running-jekyll/tutorial.md %})
{: .hands_on}

# Keep track of the changes

The material is stored in a [GitHub repository]({{ site.github_repository }}), a code hosting platform for version control and collaboration. So to develop training material, we are following the [GitHub flow](https://guides.github.com/introduction/flow/), which is based on fork, branches, and pull requests.

This can be done online via the GitHub interface or locally on your computer via command-line.

> <comment-title>Learning how to contribute</comment-title>
> Want to learn how to contribute? Check our tutorials:
> - [Contributing with GitHub via its interface]({% link topics/contributing/tutorials/github-interface-contribution/tutorial.md %})
> - [Contributing with GitHub via command-line]({% link topics/contributing/tutorials/github-command-line-contribution/tutorial.md %})
{: .comment}

# Add slides (optional)

Sometimes, you may want to have slides to support a tutorial and introduce it during a workshop. Sometimes, a set of slides is better than a tutorial to cover a specific topic.

> <hands-on-title>Add slides</hands-on-title>
>
> 1. Create a slide deck in `slides.html` following the [Slide tutorial]({% link topics/contributing/tutorials/create-new-tutorial-slides/slides.html %})
{: .hands_on}

# Conclusion


To develop a new tutorial:

1. Determine the topic
2. Create a workflow
3. Find a good toy dataset and upload it to Zenodo
4. Create the skeleton for the tutorial
5. Add the skeleton to the training material
6. Write the tutorial
7. Keep track of the changes
8. Add slides (optional)
9. Submit as a Pull Request to GitHub

For the next times, you can make it quicker.

> <hands-on-title>Generation of a tutorial</hands-on-title>
>
> 1. Determine the topic
> 2. Create your workflow on a running Galaxy instance
> 3. Add the topic name as a Tag and the tutorial title as Annotation/Notes to the workflow using the workflow editor
> 4. Create a Zenodo record with the input data
> 5. Generate the skeleton of your tutorial
>
>    - Option 1: from the [PTDK webserver](https://ptdk.apps.galaxyproject.eu/) and get the skeleton in the training material
>    - Option 2: from a workflow located on a Galaxy
>
>      ```
>      $ planemo training_init \
>             --topic_name "my-topic" \
>             --tutorial_name "my-new-tutorial" \
>             --tutorial_title "Title of the tutorial" \
>             --galaxy_url "URL to Galaxy instance in which you created the workflow" \
>             --galaxy_api_key "Your API key on the Galaxy instance" \
>             --workflow_id "ID of the workflow on the Galaxy instance" \
>             --zenodo_link "URL to the Zenodo record"
>      ```
>
> 6. Fill the remaining metadata in the `tutorial.md`
> 7. Fill the content of the `tutorial.md`
> 8. Check it by serving the website locally
>
>    > <comment-title>Serving the website locally</comment-title>
>    > Want to learn how to see the change on the website locally? [Check out our dedicated tutorial]({% link topics/contributing/tutorials/running-jekyll/tutorial.md %})
>    {: .comment}
>
> 9. Submit it to GitHub
{: .hands_on}
