---
layout: tutorial_hands_on

title: "InterMine integration with Galaxy"
questions:
    - How to export your query results from your InterMine of choice to Galaxy?
    - How to export your data sets from Galaxy to your InterMine of choice?
objectives:
    - Learn how to import/export data from/to InterMine instances
    - Understand the InterMine Interchange Dataset
time_estimation: 1h
key_points:
    - TODO
contributors:
    - danielabutano
    - yochannah
---

# Introduction
{:.no_toc}

<!-- This is a comment. -->

General introduction about the topic and then an introduction of the
tutorial (the questions and the objectives). It is nice also to have a
scheme to sum up the pipeline used during the tutorial. The idea is to
give to trainees insight into the content of the tutorial and the (theoretical
and technical) key concepts they will learn.

You may want to cite some publications; this can be done by adding citations to the
bibliography file (`tutorial.bib` file next to your `tutorial.md` file). These citations
must be in bibtex format. If you have the DOI for the paper you wish to cite, you can
get the corresponding bibtex entry using [doi2bib.org](https://doi2bib.org).

With the example you will find in the `tutorial.bib` file, you can add a citation to
this article here in your tutorial like this:
{% raw %} `{% cite Batut2018 %}`{% endraw %}.
This will be rendered like this: {% cite Batut2018 %}, and links to a
[bibliography section](#bibliography) which will automatically be created at the end of the
tutorial.


**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Import data from InterMine

> ### {% icon hands_on %} Hands-on: Import
> 1. Search Galaxy for “InterMine” (not case sensitive; “intermine” is fine too), and click on “InterMine Server” under “Get Data”.
> 2. This will redirect you to the InterMine registry, which shows a full list of InterMines and the various organisms they support. Find an InterMine that has the organism type you’re working with, and > > click on it to redirect to that InterMine
> 3. Once you arrive at your InterMine of choice, run a query as normal - this could be a search, a list results page, a template, or a query in the query builder. Eventually you’ll be presented with an > > InterMine results table.
> 4. Click on Export (top right). This will bring up a modal window.
> 5. Select “Send to Galaxy” and double-check the Galaxy Location is correct.
> 6. Click on the “Send to Galaxy” button on the bottom right of the pop-up window.
>
>    > ### {% icon comment %} Allow popups
>    >
>    > If you get an error when you click on the “Send to Galaxy” button, please make sure to allow popups and try again.
>    {: .comment}
>
{: .hands_on}
You have now exported your query results from InterMine to Galaxy.


# Export data into InterMine
Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

## Create InterMine Interchange

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Create InterMine Interchange** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Tabular file"*: `output` (Input dataset)
>    - *"Feature Type"*: `Gene`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}


## Send data to InterMine

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
