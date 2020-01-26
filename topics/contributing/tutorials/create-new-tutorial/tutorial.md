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
contributors:
  - bebatut
  - hexylena
  - shiltemann
  - lldelisle
---

# Introduction
{:.no_toc}

Galaxy is a great solution to train bioinformatics concepts:

- numerous bioinformatics tools are available (almost 6,000 in the ToolShed)
- it can be used by people without any computer science skills
- it trains to use technology, outlining available resources and efforts that have made them accessible to researchers
- it is scalable

In 2016, the Galaxy Training Network decide to set up a new infrastructure for easily delivering Galaxy related training material. The idea was to develop something open, online, based on a community effort, and on top of the Galaxy platform.

We took inspiration from [Software Carpentry](https://software-carpentry.org) and collected everything on a GitHub repository: [https://github.com/galaxyproject/training-material](https://github.com/galaxyproject/training-material).
We decided on a structure focusing on tutorials with hands-on activities; fitting both for online self-training but also for workshops. Each tutorial follows the same structure and comes with a virtualised instance allowing you to run the training anywhere you have resources available.

Here you will learn how to create a new tutorial by developing a small tutorial that explains how to use BLAST.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

> ### {% icon comment %} Comment
> This tutorial explains the different steps to create a tutorial for the Galaxy Training Material.
> It may require some knowledge that you may not have or do not have the time to learn. If this is the case, you can create a skeleton of a tutorial with whatever existing materials you have, using your prefered text editor, and then share it with us by opening [issue on GitHub]({{ site.github_repository }}/issues/new), writing us on [Gitter]({{ site.gitter_url }}), or sending us an [email](mailto:{{ site.email }}).
{: .comment}

# Define the topic

The first question we need to answer is in which topic to place our new tutorial. This can be tricky. When we structured the repository, we decided to use the categories that are used in the [ToolShed](https://toolshed.g2.bx.psu.edu/) as our initial list of topics. Since every tool uploaded to the ToolShed must be in at least one category, you can look at the main tools in your tutorial and see which categories they are placed in within the ToolShed. This can provide a guide for where you might put your new tutorial. For example, this tutorial will rely on the NCBI Blast+ tool:

> ### {% icon hands_on %} Hands-on: Defining the topic for the tutorial
>
> 1. Search for NCBI Blast+ on the [ToolShed](https://toolshed.g2.bx.psu.edu/)
> 2. Check in which category it has been placed
>
> > ### {% icon solution %} Solution
> > There are a couple steps to reaching the answer:
> > 1. Search for `ncbi blast+`
> > 2. Press the <kbd>Enter</kbd> key to search
> > 3. Click on the result named `ncbi_blast_plus`
> > 4. At the bottom of this page there is a box labelled "Categories"
> >
> > It is placed in two categories, "Next Gen Mappers" and "Sequence Analysis"
> {: .solution}
>
{: .hands_on}

> ### {% icon comment %} Creating a new topic
> Want to create a new topic? [Check out our tutorial to create a new topic]({% link topics/contributing/tutorials/create-new-topic/tutorial.md %})
{: .comment}

# Keep track of the changes

The material is stored in a [GitHub repository]({{ site.github_repository }}), a code hosting platform for version control and collaboration. So to develop training material, we are following the [GitHub flow](https://guides.github.com/introduction/flow/), which is based on fork, branches, and pull requests.

This can be done online via the GitHub interface or locally on your computer via command-line.

> ### {% icon comment %} Learning how to contribute
> Want to learn how to contribute? Check our tutorials:
> - [Contributing with GitHub via its interface]({% link topics/contributing/tutorials/github-interface-contribution/tutorial.md %})
> - [Contributing with GitHub via command-line]({% link topics/contributing/tutorials/github-command-line-contribution/tutorial.md %})
{: .comment}

# Create the directory for the tutorial

Each training material is related to a topic. All training materials (slides, tutorials, ...) related to a topic are found in a dedicated directory (*e.g.* `transcriptomics` directory contains the material related to exome sequencing analysis). Each topic have the following structure:

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

Once the topic has been chosen and you set up your contribution environment, you can create the tutorial. An ideal tutorial in the Galaxy Training Network contains:
- a tutorial file `tutorial.md` written in Markdown with hands-on
- an optional slides file `slides.md` in Markdown with slides to support the tutorial
- a directory `tours` with Galaxy Interactive Tours to reproduce the tutorial
- a directory `workflows` with workflows extracted from the tutorial
- a YAML file `data-library.yaml`  with the links to the input data needed for the tutorial

The most important file is the `tutorial.md` where the content of the tutorial is. The other files are there to support the tutorial and make it robust and usable across many environments.

> ### {% icon hands_on %} Hands-on: Create all the required files and folders structures automatically
>
> 1. Run (by adapting the information between the quotes)
>
>    ```
>    $ planemo training_init \
>             --topic_name "my-topic" \
>             --tutorial_name "my-new-tutorial" \
>             --tutorial_title "Title of the tutorial" \
>             --hands_on
>    ```
>
> 2. Check that a new directory (with your tutorial name) has been generated in the topic folder
> 3. Make sure that Jekyll is running
>
>    > ### {% icon comment %} Jekyll
>    > Want to learn how to start Jekyll? [Check out our tutorial to serve the website locally]({% link topics/contributing/tutorials/running-jekyll/tutorial.md %})
>    {: .comment}
>
> 2. Check if the tutorial has been correctly added at [http://localhost:4000/training-material/](http://localhost:4000/training-material/)
{: .hands_on}

# A toy dataset

Our tutorials try to follow the "learn by doing" approach; they combine both theoretical and practical sections. The practical sections (or hands-on) are supposed to be done on Galaxy.

The first task is to select some data to use for the Hands-on sections. The selected data must be informative enough to illustrate the meaning of using a tool or a given technique, but not too big to require long waiting times for processing during a workshop. Upload and download of files into and out of Galaxy is usually quick, but the time taken for a tool to run can be long. Tool run times of no more than 10-15 mins are recommended. Typically, the selected data should be the informative subset of a full real-life dataset.

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

We would then develop the tutorial and test it on this toy dataset. Once we were ready to share it, we would upload the datasets on [Zenodo](https://zenodo.org/) to store them on long-term and obtain a dedicated DOI in the [Galaxy training network community](https://zenodo.org/communities/galaxy-training/?page=1&size=20).

> ### {% icon hands_on %} Hands-on: Upload the dataset to Zenodo
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
>     > ### {% icon comment %} No possible changes in the files after publication
>     > File addition, removal or modification are not allowed after you have published your upload.
>     > So be careful when you start your upload that all your needed files are ready.
>     >
>     > The metadata can be changed after publication.
>     {: .comment}
>
> 6. Search for and Select *Galaxy training network* in **Communities**
> 7. Select *Dataset* in **Upload type**
> 8. Use the title of your tutorial and mention also Galaxy Training Material
> 9. Add all the persons who contributed to the tutorial as authors
> 10. Add a short description of the tutorial and a link to the training material website
> 11. Keep *Open Access* as **Access right** and *Creative Commons Attribution 4.0* as **License**
> 12. Fill out any remaining information
> 13. Click on **Publish**
> 14. Copy the DOI link in the new page
> 15. Paste the link in `zenodo_link` in the tutorial header
{: .hands_on}

# Write the tutorial

Now that you have the structure in place, you can then fill the tutorial per se.

> ### {% icon hands_on %} Hands-on: Write the tutorial
>
> 1. Open the `tutorial.md` file with your favorite text editor
> 2. Fill out the tutorial by following the [dedicated tutorial]({% link topics/contributing/tutorials/create-new-tutorial-content/tutorial.md %})
> 2. (Optional) Build the website locally and check that the tutorial is there by following the [Jekyll tutorial]({% link topics/contributing/tutorials/running-jekyll/tutorial.md %})
{: .hands_on}

# Add some technical support (recommended)

To able to run the tutorial, we need a Galaxy instance where the needed tools and the data are available. We need then to describe the required technical infrastructure. Tools are installed based on the workflows in the `workflows` directory.

This description will be used to automatically set up a Docker Galaxy flavour, to set up an existing Galaxy instance and also to test if a public Galaxy instance is able to run the tool.

The technical support are different files:

- workflow file(s) in the `workflows` directory
- the `data-library.yaml` file with the links to the input data needed for the tutorial
- interactive tour file in the directory `tours` directory

> ### {% icon hands_on %} Hands-on: Add technical support for the tutorial
>
> 1. Add some technical support for the tutorial following the [tutorial]({% link topics/contributing/tutorials/create-new-tutorial-technical/tutorial.md %})
>    - Add the workflow
>    - (Recommended) Generate the `data-library.yaml`
>    - (Optional) Create an interactive tour
{: .hands_on}

# Add slides (optional)

Sometimes, you may want to have slides to support a tutorial and introduce it during a workshop. Sometimes, a set of slides is better than a tutorial to cover a specific topic.

> ### {% icon hands_on %} Hands-on: Add slides
>
> 1. Create a slide deck in `slides.html` following the [Slide tutorial]({% link topics/contributing/tutorials/create-new-tutorial-slides/slides.html %})
{: .hands_on}

# Conclusion
{:.no_toc}

To develop a new tutorial:

1. Determine the topic
2. Create the directory for the tutorial
3. Add some metadata
4. Find a good toy dataset and upload it on Zenodo
5. Write the tutorial
6. Add some technical support (recommended)
7. Add slides (optional)

For the next times, you can make it quicker.

> ### {% icon hands_on %} Hands-on: Generation of a tutorial
>
> 1. Determine the topic
> 2. Create your workflow on a running Galaxy instance
> 3. Create a Zenodo record with the input data
> 4. Generate the skeleton of your tutorial
>    - option 1: from a workflow located on a Galaxy
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
>    - option 2: from a local workflow file (`.ga`) (use only if your workflow is composed of tools from the main ToolShed)
>
>      ```
>      $ planemo training_init \
>             --topic_name "my-topic" \
>             --tutorial_name "my-new-tutorial" \
>             --tutorial_title "Title of the tutorial" \
>             --workflow "path/to/workflow" \
>             --zenodo_link "URL to the Zenodo record"
>      ```
>      You can use the example workflow file located in `topics/contributing/tutorials/create-new-tutorial/workflows/example-workflow.ga` if
>      you do not have a workflow of your own. This is the workflow belonging to the *Galaxy 101* introduction tutorial.
>
> 5. Fill the remaining metadata in the `tutorial.md`
> 6. Fill the content of the `tutorial.md`
> 7. Check it using Jekyll
{: .hands_on}
