---
layout: tutorial_hands_on

title: "Galaxy 101 for everyone"
zenodo_link: ""
level: Introductory
questions:
  - "What are the differences between the Iris species?"
objectives:
  - "Familiarize yourself with the basics of Galaxy"
  - "Learn how to obtain data from external sources"
  - "Learn how to tag datasets"
  - "Learn how to run tools"
  - "Learn how histories work"
  - "Learn how to create a workflow"
  - "Learn how to share your work"
time_estimation: "1H30M"
key_points:
  - "Galaxy provides an easy-to-use graphical user interface for often complex command-line tools"
  - "Galaxy keeps a full record of your analysis in a history"
  - "Workflows enable you to repeat your analysis on different data"
  - "Galaxy can connect to external sources for data import and visualization purposes"
  - "Galaxy provides ways to share your results and methods with others"
contributors:
  - chrisbarnettster
  - michelemaroni89
  - annefou
  - nagoue
  - olanag1
  - tnabtaf
---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

This practical aims to familiarize you with the Galaxy user interface. 
It will teach you how to perform basic tasks such as importing data, 
running tools, working with histories, 
creating workflows, and sharing your work.

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

> ### {% icon comment %} Background
> The Iris flower data set or Fisher’s Iris data set is a multivariate data set introduced by the British statistician and biologist Ronald Fisher in his 1936 paper ({% cite Fisher1936 %}). 
> Each row of the table represents an iris flower, including its species and dimensions of its botanical parts, sepal and petal, in centimeters.
> For more history of this dataset read here [Wikipedia](https://en.wikipedia.org/wiki/Iris_flower_data_set).
{: .comment}

# Create a new history

Galaxy allows you to create histories.
They are the pipeline of operations performed on a certain dataset in order to achieve the desired result.
After that your pipeline has been tested and leads to the correct result, 
it can be serialized as a workflow.
A workflow allows to obtain the same result with a different dataset, just by changing the input.

Overall a history represent and experimental lab book, or a recipe. 
A workflow is the concatenation of one or multiple histories as a series of building blocks 
for replicating an experimetal result or a recipe.

## Inspect Iris dataset

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial and give it a proper name such as `Galaxy 101 for everyone`
>
>    {% include snippets/create_new_history.md %}
>    {% include snippets/rename_history.md %}
>
> 2. Import `iris.csv` from [Zenodo](https://zenodo.org/record/1319069/files/iris.csv) or from the data library (ask your instructor)
>
>    ```
>    https://zenodo.org/record/1319069/files/iris.csv
>    ```
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
>    As default, Galaxy takes the link as name, so rename them.
>
> 3. Rename the dataset to `iris`
>
>    {% include snippets/rename_dataset.md %}
>
> 4. Check the datatype 
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 5. Add the tag `iris` to the dataset
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

# What does the dataset contain?

Now we are going to inspect the dataset and count how many different species are in the dataset and how many sample for each species.

## How many different species are in the dataset?

## How many sample for each species?

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


## Sub-step with **My Tool**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **My Tool** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: File
>    - *"Parameter"*: `a value`
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


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
