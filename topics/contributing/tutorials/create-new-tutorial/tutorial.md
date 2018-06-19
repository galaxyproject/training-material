---
layout: tutorial_hands_on
topic_name: contributing
tutorial_name: create-new-tutorial
---

# Introduction
{:.no_toc}

Galaxy is a great solution to train the bioinformatics concepts:

- numerous bioinformatics tools are available (almost 5,000 in the ToolShed)
- it can be used by people without any computer science skills
- it trains to use technology, outlining available resources and efforts that have made them accessible to researchers
- it is scalable

In 2016, the Galaxy Training Network decide to set up a new infrastructure for delivering easily Galaxy related training material. The idea was to develop something open and online based on a community effort, as always in Galaxy.

We took inspiration from [Software Carpentry](https://software-carpentry.org) and collected everything on a GitHub repository: [https://github.com/galaxyproject/training-material ](https://github.com/galaxyproject/training-material).
We decided on a structure based on tutorials with hands-on, fitting both for online self-training but also for workshops, grouped in topics. Each tutorial follows the same structure and comes with a virtualised instance to run the training everywhere.

In this tutorial, you will learn how to create a new tutorial by developing a small tutorial to explain how to use BLAST.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

> ### {% icon comment %} Comment
> This tutorial explains the different steps to create a tutorial for the Galaxy Training Material.
> It may require some knowledge that you may not have or do not have the time to learn. If it is the case, you can create a skeleton of a tutorial using your prefered text editor and giving us a link to it by opening [issue on GitHub]({{ site.github_repository }}/issues/new), writing us on [Gitter]({{ site.gitter_url }}) or sending us an [email](mailto:{{ site.email }}).
{: .comment}

# Define the topic

The first step we need to define is in which topic to place our new tutorial. This can be tricky: when we structured the repository, we decided to use as topics the categories that are used in the [ToolShed](https://toolshed.g2.bx.psu.edu/). The ToolShed assigns a category to each tool. Therefore, to decide where to put your tutorial, have a look at which ToolShed's category the main tools in your tutorial belong. For example, this tutorial will rely on the NCBI Blast+ tool.

> ### {% icon hands_on %} Hands-on: Defining the topic for the tutorial
>
> 1. Search for NCBI Blast+ on the [ToolShed](https://toolshed.g2.bx.psu.edu/)
> 2. Check in which category it has been placed
{: .hands_on}

> ### {% icon comment %} Creating a new topic
> Want to create a new topic? [Check out our tutorial to create a new topic]({{ site.baseurl }}/topics/contributing/tutorials/create-new-topic/tutorial.html)
{: .comment}

# Keep track of the changes

The material is stored in a [GitHub repository]({{ site.github_repository }}), a code hosting platform for version control and collaboration. So to develop training material, we are following the [GitHub flow](https://guides.github.com/introduction/flow/), which is based on fork, branches, and pull requests.

It can be done online via the GitHub interface or locally on your computer via command-line.

> ### {% icon comment %} Learning how to contribute
> Want to learn how to contribute? Check our tutorials:
> - [Contributing with GitHub via its interface]({{ site.baseurl }}/topics/contributing/tutorials/github-interface-contribution/tutorial.html)
> - [Contributing with GitHub via command-line]({{ site.baseurl }}/topics/contributing/tutorials/github-command-line-contribution/tutorial.html)
{: .comment}

# Create the directory for the tutorial

Once the topic has been chosen and you set up your contribution environment, you can create the tutorial. An ideal tutorial in the Galaxy Training Network contains:
- a tutorial file `tutorial.md` written in Markdown with hands-on
- an optional slides file `slides.md` in Markdown with slides to support the tutorial
- a directory `tours` with Galaxy Interactive Tours to reproduce the tutorial
- a directory `workflows` with workflows extracted from the tutorial
- a YAML file `tools.yaml` with the description of needed tools to run the tutorial
- a YAML file `data-library.yaml`  with the links to the input data needed for the tutorial

The most important file is for sure the `tutorial.md` where the content of the tutorial is. The other files are there to support the tutorial and make it robust and usable in most environment.

To ease this process, we created a template for new tutorials, complete with all aforementioned requirements.

> ### {% icon hands_on %} Hands-on: Copy the required files
>
> 1. Copy the `tutorial1` directory (you can find it in `templates/tutorials/`) in `topics/sequence-analysis/tutorials`
> 2. Rename the copied directory to `similarity-search`
{: .hands_on}

# Add metadata

To link a tutorial to a topic and also define some technical and pedogical support, you need to add some metadata related to the tutorial in the `metadata.yaml` file of the topic. Once it is filled you can run the Galaxy Training material website locally to check that the new tutorial is accessible.

> ### {% icon hands_on %} Hands-on: Add metadata
>
> 1. Check out and run our [metadata tutorial]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-metadata/tutorial.html)
> 2. (Optional) Build the website locally by following the [Jekyll tutorial]({{ site.baseurl }}/topics/contributing/tutorials/running-jekyll/tutorial.html) and check that the tutorial is referenced in the topic page  
{: .hands_on}

# Find a good toy dataset and upload it on Zenodo

The tutorials are developed following the "learning by doing" approach. They combine both some theoretical and some practical sections. The practical sections (or hands-on) are supposed to be done on Galaxy.

The first question to come is what data to use for walking the tutorial through the hands-on sections. The selected data must be informative enough to illustrate the meaning of using a tool or a given technique, but not too big to require long waiting times for its processing during a workshop. Typically, the selected data should be the informative subset of a full real-life dataset.

For example for our tutorial, we generated a small dataset by

- Taking one 16S sequences (used for test of a Galaxy tool)
- Generating a reference database
    - Blasting it on the NR database on [NCBI Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome)
    - Extracting one similar sequence found with Blast
    - Searching and extracting 2 other sequences of the same species using the [NCBI Nucleotide database](https://www.ncbi.nlm.nih.gov/nuccore)

We then developed the tutorial and tested it on this toy dataset. Once we were ready to share it, we uploaded the datasets on [Zenodo](https://zenodo.org/) to store them on long-term and obtain a dedicated DOI (in the [Galaxy training network community](https://zenodo.org/communities/galaxy-training/?page=1&size=20)).

> ### {% icon hands_on %} Hands-on: Add the dataset on Zenodo
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
> 12. Fill any remaining information
> 13. Click on **Publish**
> 14. Copy the DOI link in the new page
> 15. Paste the link in `zenodo_link` in the tutorial section of the `metadata.yaml` file
{: .hands_on}

# Fill the tutorial

Now that you have the structure in place, you can the fill the tutorial per se. 

> ### {% icon hands_on %} Hands-on: Fill the tutorial
>
> 1. Open the `tutorial.md` file with your favorite text editor
> 2. Fill the tutorial by following the [dedicated tutorial]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html) 
> 2. (Optional) Build the website locally and check that the tutorial is there by following the [Jekyll tutorial]({{ site.baseurl }}/topics/contributing/tutorials/running-jekyll/tutorial.html)
{: .hands_on}

# Add some technical support (recommended)

To able to run the tutorial, we need a Galaxy instance where the needed tools and the data are available. We need then to describe the required technical infrastructure. 

This description will be used to automatically set up a Docker Galaxy flavour, to set un an existing Galaxy instance and also to test if a public Galaxy instance is able to run the tool.

The technical support are different files:

- workflow file(s) in the `workflows` directory
- the `tools.yaml` file with the description of needed tools to run the tutorial
- the `data-library.yaml` file with the links to the input data needed for the tutorial
- interactive tour file in the directory `tours` directory

> ### {% icon hands_on %} Hands-on: Add technical support for the tutorial
>
> 1. Add some technical support for the tutorial following the [tutorial]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)
>    - Add the workflow
>    - (Recommended) Generate the `tools.yaml`
>    - (Recommended) Generate the `data-library.yaml`
>    - (Optional) Create an interactive tour
> 2. Update the `metadata.yaml` file given the technical support added
{: .hands_on}

# Add slides (optional)

Sometimes, you may want to have slides to support a tutorial and introduce it during a workshop. Sometimes, a set of slides is better than a tutorial to cover a tutorial.

> ### {% icon hands_on %} Hands-on: Add slides
>
> 1. Create a slide deck in `slides.html` following the [Slide tutorial]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-slides/slides.html)
> 2. Update the `metadata.yaml` file by putting `yes` in `slides`
{: .hands_on}

# Conclusion
{:.no_toc}

To develop a new tutorial:

1. Define the topic
2. Create the directory for the tutorial
3. Add some metadata
4. Find a good toy dataset and upload it on Zenodo
5. Fill the tutorial
6. Add some technical support (recommended)
7. Add slides (optional)
