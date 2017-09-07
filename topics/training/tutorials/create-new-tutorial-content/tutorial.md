---
layout: tutorial_hands_on
topic_name: training
tutorial_name: create-new-tutorial-content
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

In this tutorial, you will learn how to write your first tutorial in markdown and contribute it to the Galaxy Training Network.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Setting up a new tutorial

Here, we want to develop a small tutorial to explain how to use BLAST.

## Clone the Galaxy Training material repository

Before anything, we need to get a local copy of the content of the GitHub repository by cloning it

> ### {% icon hands_on %} Hands-on: Clone the GitHub repository
>
> 1. Clone the repository locally with: `git clone https://github.com/galaxyproject/training-material.git`
> 2. Check that you have the same structure as the one in [GitHub](https://github.com/galaxyproject/training-material)
{: .hands_on}

## Defining the topic

The first step we need to define is in which topic putting our tutorial. This first step can be tricky.

When we structured the repository, we decided here to use as topic the names of the categories in the [ToolShed](https://toolshed.g2.bx.psu.edu/). So when decided where to put your tutorial, you can look in which ToolShed's category are the main tools used in the tutorial and use this category as topic. For example, this tutorial will rely on the NCBI Blast+ tool.

> ### {% icon hands_on %} Hands-on: Defining the topic for the tutorial
>
> 1. Search for NCBI Blast+ on the [ToolShed](https://toolshed.g2.bx.psu.edu/)
> 2. Check in which category it is
>
>    > ### {% icon question %} Questions
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

> ### {% icon hands_on %} Hands-on: Copy the needed file
>
> 1. Copy the `tutorial1` folder (you can find it in `templates/tutorials/`) in `topics/sequence-analysis/topics`
> 2. Rename the folder into `similarity-search`
{: .hands_on}

## Keeping track of the changes

Once you started to change something, we need to keep track of these changes with a version control system (VCS). We are using Git as VCS and GitHub as hosting service.

This repository is developed collaboratively with more than 40 contributors. For the collaboration, we are using the [GitHub flow](https://guides.github.com/introduction/flow/), which is based on forks, branches and pull requests.
We will explain GitHub flow later and show you now how to start keeping track of the changes:

> ### {% icon hands_on %} Hands-on: Start keeping track of the changes
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


# Filling the tutorial content

Once we set up the infrastructure, we are ready to write the tutorial

## Finding a good toy dataset

The first question to come is which data using for the tutorial and particularly for the hands-on parts. This data must be informative enough to illustrate the use of a tool or a technic, but it must not be too big to be able to be run during workshop or locally. Typically, this is a subset of a full dataset where the informative data has been extracted.

For example for our tutorial, we generated a small dataset by

- Taking one 16S sequences (used for test of a Galaxy tool)
- Generating a reference database
    - Blasting it on the NR database on [NCBI Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome)
    - Extracting one similar sequence found with Blast
    - Searching and extracting 2 other sequences of the species using the [NCBI Nucleotide database](https://www.ncbi.nlm.nih.gov/nuccore)

We then developed the tutorial and tested it on this dataset. Once we were ready to share it, we uploaded the datasets on [Zenodo](https://zenodo.org/) to obtain a dedicated DOI (in the [Galaxy training network community](https://zenodo.org/communities/galaxy-training/?page=1&size=20)).

## Filling the tutorial content

The content of the tutorial will go in the `tutorial.md`. The syntax of this file is really simple., as well as its structure:

```
---
layout: tutorial_hands_on
topic_name: training
tutorial_name: create-new-tutorial
---
# Introduction
{:.no_toc}

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
{:.no_toc}

blabla
```

### Some metadata on the top

On the top, there is some metadata:

- `layout: tutorial_hands_on` (keep the default)
- `topic_name: training` with the name of the topic
- `tutorial_name: create-new-tutorial` with the name of tutorial

These metadata are there to help the templating system to make the connection between the file and the global [metadata]({{site.url}}/topics/training/tutorials/create-new-tutorial-metadata/tutorial.html).
If they are not correctly defined the tutorial can not be found on the website.

> ### {% icon hands_on %} Hands-on: Fix the top metadata
>
> 1. Change the `tutorial-name` and the `topic_name` to fit to the ones defined in the metadata
> 2. Check if the tutorial has been correctly added at [http://localhost:4000/topics/sequence-analysis/similarity-search ](http://localhost:4000/topics/sequence-analysis/similarity-search)
{: .hands_on}

### Content of the tutorial

Directly after the short metadata section on top the content of your tutorial starts. It is written in Markdow - a simple markup langage.

> ### {% icon tip %} Tip: Markdown
>
> Check [this cheatsheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet) to learn more how to use Markdown.
{: .tip}

The content in Markdown is transformed by our templating system into the nice webpage which add the metadata. Indeed in the `tutorial.md` file,
no need to add the name of the tutorial: it is automatically added based on the title defined in the metadata.

We recommend to structure the tutorials like this

- An introduction to introduce the tutorial with the use case, the data, the methods
- Several sections with the content of the tutorial and some hands-on parts (practicing is an important part of the learning process)
- A conclusion to summarize what has been done in the tutorial (with a scheme)

> ### {% icon hands_on %} Hands-on: Structuring the tutorial
>
> 1. Add a small introduction about the dataset
> 2. Add one or two sections with ideas for the tutorial
> 3. Add a small conclusion
{: .hands_on}

> ### {% icon tip %} Tip: Adding images, with caption
>
> To add an image in markdown file, we need to use `![](../../images/image.png)`.
>
> On the top of that, we added a small plugin to add a caption for each image:
>
> ![This figure shows an example of a figure with a caption](../../images/image_caption_screenshot.png "Example of a figure with a caption")
>
> "Figure" and the number are automatically added and the caption is added by adding the information in the markdown call of the image:
>
>   ```
>   ![A textual description of the image](../images/image.png "This is my super caption")
>   ```
>
> We can also cross-reference the figure inside our markdown with an anchor. For example, we can link to [the previous figure](#figure-1) using `[the display text](#figure-nb)` (with changing `nb` to the figure number).
{: .tip}

> ### {% icon tip %} Tip: Writting mathematical expressions
>
> Mathematical expressions can be written in LaTeX: they will be rendered with [MathJax](https://www.mathjax.org/).
>
> It is easy: surround your math content with two `$` signs, like with a math block:
>
> - inline, *e.g.* `$$ 5 + 5 $$` will be rendered as $$ 5 + 5 $$
> - not inline, *e.g.* `$$ 5 + 5 $$` alone in new line will be rendered as
>
>   $$ 5 + 5 $$
>
>
> If you don't want to start an inline math statement, just escape the dollar signs and they will be treated as simple dollar signs.
>
>    > ### {% icon comment %} Comments
>    > LaTeX code that uses the pipe symbol `|` in inline math statements may lead to a line being recognized as a table line.
>    > This can be avoided by using the `\vert` command instead of `|`
>    {: .comment}
{: .tip}

### Improving the learning experience

To improve the learning experience in our tutorial, we defined some boxes to highlight.

They are defined always with the same structure:

```
> ### <an emoji> Type of box: Name of the box
> list
{: .type_of_box}
```

This structure needs to be respected otherwise it would not be interpreted correctly by the templating system. The different defined boxes are:

- Overview

    This box at the top of each tutorial is automatically generated using the metadata we defined

    > ### {% icon hands_on %} Hands-on: Checking the metadata
    >
    > 1. Check that the metadata added previously are correctly filling the overview box
    >
    >    > ### {% icon question %} Questions
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
        > In this tutorial, we will deal with:
        >
        > 1. TOC
        > {:toc}
        >
        {: .agenda}

    No need to fill the list: it will be done automatically based on the title of the sections.

    To avoid to add the "Introduction" and "Conclusion", we add `{:.no_toc}` below the section name.

    ![Example of agenda box](../../../../shared/images/tutorial_agenda_box.png "Example of agenda box")

    > ### {% icon hands_on %} Hands-on: Add an agenda box to the tutorial
    >
    > 1. Add an agenda box to the tutorial that fit the structure we defined previously
    {: .hands_on}

- Hands-on

    We think that doing is important in the learning process. So we emphasize it by adding regularly some hands-on sections where the trainees can do by themselves some analyses. We designed some special boxes to make these sections easy to find.

        > ### {% icon hands_on %} Hands-on: Sorting BAM dataset
        >
        > 1. **Sort BAM dataset** {% icon tool %}: Sort the paired-end BAM file by "Read names" with **Sort BAM
        {: .hands_on}

    ![Example of hands-on box](../../../../shared/images/tutorial_hand_on_box.png "Example of hands-on box")

    with the

    - `{% icon hands_on %}` emoji to define that is an hands-on
    - Short imperative sentence to make it easy to identify the tasks
    - Name of the tool in bold with the `{% icon tool %}` emoji to make it easy to identify a Galaxy tool
    - Parameters for the tool as a sublist<br/>
    <br/>

    > ### {% icon hands_on %} Hands-on: Add an hands-on box
    >
    > 1. Add an hands-on box to run a BLAST of the small sequence dataset against the chosen database
    {: .hands_on}

-  Questions

    The questions are then to force the trainees to think about what they are currently doing and to put things in perspective.
    They are also a way to help the instructors to expose and clearify misunderstanding earily on.

        > ### {% icon question %} Questions
        >
        > 1. Why are some tests filtered?
        > 2. Does it improve the *p*-value distribution?
        >
        >    <details>
        >    <summary>Click to view answers</summary>
        >    Content goes here.
        >    </details>
        {: .question}

    ![Example of question box](../../../../shared/images/tutorial_question_box.png "Example of question box")

    The questions has to be quick to answer. They can be small or also multiple choice (MCQs).
    With well choosen wrong answers MCQs can do much more than just measure how much someone knows.

    In the box below and hiffen we add also the correct answer, so that self-trainees can check the solution and its explanation.

    > ### {% icon hands_on %} Hands-on: Add a question box
    >
    > 1. Add an hands-on box to construct the BLAST database
    {: .hands_on}

- Tips

        > ### {% icon tip %} Tip: Importing data via links
        >
        > * Copy the link location
        > * Open the Galaxy Upload Manager
        > * Select **Paste/Fetch Data**
        > * Paste the link into the text field
        > * Press **Start**
        {: .tip}

    ![Example of tip box](../../../../shared/images/tutorial_tip_box.png "Example of tip box")

- Comments

        > ### {% icon comment %} Comments
        > - Edit the "Database/Build" to select "dm3"
        > - Rename the datasets according to the samples
        {: .comment}

    ![Example of comment box](../../../../shared/images/tutorial_comment_box.png "Example of comment box")

- Key points

    This last box of the tutorial is automatically filled with the take-home messages defined in the metadata


To render the boxes correctly, the syntax needs to be correct. If it does not work have a look at similar tutorials and get inspiration.
The boxes can be nested, *e.g.* for having tips inside a hands-on:

```
> ### {% icon hands_on %} Hands-on: Defining the topic for the tutorial
>
> 1. Search for NCBI Blast+ on the [ToolShed](https://toolshed.g2.bx.psu.edu/)
> 2. Check in which category it is
>
>    > ### {% icon question %} Questions
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

# Adding slides (optional)

Sometimes we also need slides to support the tutorial. With the current infrastructure, we also provide this possibility
to serve on the website slides related to the tutorial.

The slides are written in Markdown (only the file extension is .html), as the tutorial and are rendered as a webpage thanks to [`Remark`](https://remarkjs.com). However this is not done automatically. We first need to tell the templating system to search for the slides by changing `slides` in the metadata from `no` to `yes`.

Once it is done, the slides for our tutorial will be accessible at [http://localhost:4000/topics/sequence-analysis/tutorials/similarity-search/slides.html ](http://localhost:4000/topics/sequence-analysis/tutorials/similarity-search/slides.html)

We can now fill the `slides.html` file:


```
---
layout: tutorial_slides
topic_name: "sequence-analysis"
tutorial_name: search-similarity
logo: "GTN"
---

# What is similarity search?

---

### What is blast?

![](../../images/ecker_2012.jpg)

---

### Next slide
```

On the top, as for the `tutorial.md`, we need to define some metadata to link the slides to the correct tutorial. The metadata of the tutorial can then be used also to automatically fill the first slides (with the title, the requirements,...)  and the last slide (take-home message).

After, each new slide is introduced by `---`, and the content of each slide is written in Markdown.

> ### {% icon hands_on %} Hands-on: Add some slides for the tutorial
>
> 1. Add a few slides to explain the similarity search
> 2. Make sure they are accessible and correctly generated
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
