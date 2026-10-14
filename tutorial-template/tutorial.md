---
layout: tutorial_hands_on

title: Title of the tutorial
level: Introductory # or Intermediate or Advanced

zenodo_link: ''

questions:
- Which scientific questions are addressed by the tutorial?
- Which concepts are important to know for this type of data?

objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentences describing what a learner should be able to do once they
  have completed the tutorial
- You can use Bloom's Taxonomy to write effective learning objectives

time_estimation: 2H30M

key_points:
- The take-home messages
- They will appear at the end of the tutorial
- They must be in the form of statements.

contributions:
  authorship:
    - contributor1
    - contributor2
  editing:
    - contributor1

---


Start with a general introduction about the topic and then an introduction of the tutorial (the questions and the objectives).
It is nice also to have an overview to sum up the pipeline used during the tutorial.
The idea is to give to trainees insight into the content of the tutorial and the (theoretical and technical)
key concepts they will learn.

You may want to **cite some publications**; this can be done by adding citations to the bibliography file
(`tutorial.bib` file next to your `tutorial.md` file). These citations must be in bibtex format.
If you have the DOI for the paper you wish to cite, you can get the corresponding bibtex entry using
[doi2bib.org](https://doi2bib.org).

With the example you will find in the `tutorial.bib` file, you can add a citation to
this article here in your tutorial like this:

{% raw %} `{% cite Batut2018 %}`{% endraw %}.

This will be rendered like this: {% cite Batut2018 %}, and links to a
[bibliography section](#bibliography) which will automatically be created at the end of the
tutorial.

<!-- This is a comment. -->

**Please follow our
[tutorial to learn how to fill the Markdown]({% link topics/contributing/tutorials/create-new-tutorial-content/tutorial.md %})**

And have a look at the tutorial about how to [preview the GTN website as you edit your training material]({% link topics/contributing/tutorials/rendering_gtn/tutorial.md %}).

<!-- then we add the agenda, this will be automatically generated from your section headings -->

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}



# Title for your first section

Give some background about what the trainees will be doing in the section.
Remember that many people reading your materials will likely be novices,
so make sure to explain all the relevant concepts.

Tutorials are written in Markdown formats, so **this is bold**, and *this is italics*.
This is a hyperlinked URL: <https://usegalaxy.eu> or with link text: [Galaxy EU](https://usegalaxy.eu)


## Get data

Provide a bit of background about the data, what it is, where it came from, etc.

Then we usually start with a data upload hands-on box:

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    https://zenodo.org/record/XXXXXXX/files/myfile.ext
>    https://zenodo.org/record/XXXXXXX/files/myfile2.ext
>    [...]
>    ```
>
>    **TODO** update the zenodo links in the box above
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. **Rename** {% icon galaxy-pencil %} dataset `XX` to `YY`
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
> 4. Check that the datatype is `mytype`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="mytype" %}
>
> 5. **Add a tag** {% icon galaxy-tags %} named `#mytag` to dataset `XX`
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}


# Title of your second section

Below are a series of example boxes and other tutorial elements for you to use as a starting point.

Have fun!


## Title for a subsection
Section and subsection titles will be displayed in the tutorial index on the left side of
the page and in the agenda, so try to make them informative and concise!

### Title of a sub-subsection

You can nest sections as needed



# Title of the section usually corresponding to a conceptual step in the analysis

Start with a description of the step: some background and some theory.
Some images can be added there to support the theory explanation:

![Alternative text]({% link topics/contributing/images/website_screenshot.png "Legend of the image (optional)")

The idea is to keep the theory description before quite simple to focus more on the practical part.

*Consider adding a detail boxes to expand on the theory* these are collapsed by default
so that they can be easily skipped by learners, but expanded by anybody who really wants
to dive into the details.

> <details-title> More details about the theory </details-title>
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


## Subsection Title

Some intro text learners need to know about the next hands-on step

> <hands-on-title> Task description </hands-on-title>
>
> Here we describe the tools or other actions learners should perform
>
> 1. {% tool [Tool Name](toolshed.g2.bx.psu.edu/repos/UUU/VVV/WWW/v1.0+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `myfile.txt`
>    - *"Parameter name"*: `value`
>    - *"Another Parameter"*: `value`
>
> 2. {% tool [My Second Tool](toolshed.g2.bx.psu.edu/repos/UUU/VVV/WWW/v1.0+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: File
>    - *"Parameter name"*: `value`
>    - *"Another Parameter"*: `value`
>
>    ***TODO***: *Get tool's toolshed ID from Galaxy and replace above*
>
>    ***TODO***: *you can add a comment, tip box or FAQ with more info*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
{: .hands_on}


Discuss the outputs of the hands on steps, what it means scientifically, etc

*Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}



# Conclusion

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
