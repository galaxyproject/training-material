---
layout: tutorial_hands_on

title: Running Molecular Dynamics simulations
zenodo_link: ''
questions:
- Which biological questions are addressed by the tutorial?
- Which bioinformatics techniques is important to know for this type of data?
objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentence describing what a learner will be able to do once they
  have done the tutorial
- You can use the Bloom's Taxonomy to write effective learning objectives
requirements:
  -
    title: "Setting up molecular systems"
    type: "internal"
    link: "/computational-chemistry/tutorials/setting-up-molecular-systems/tutorial.html"
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- chrisbarnettster
- tsenapathi

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

General introduction about the topic and then an introduction of the
tutorial (the questions and the objectives). It is nice also to have a
scheme to sum up the pipeline used during the tutorial. The idea is to
give to trainees insight into the content of the tutorial and the (theoretical
and technical) key concepts they will learn.

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

# Title for your first section

Give some background about what the trainees will be doing in the section.

Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

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

# Analysis with NAMD (imported from CHARMM-GUI)

This workflow uses NAMD as a molecular dynamics engine.
An NVT simulation is followed by an NPT simulation.


Is this converged?

If you've imported from CHARMM-GUI read through the theory then skip to [namd_nvt](#namd_nvt) otherwise carry on and start with [setup](#setup)


It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Snapshot of CHARMM and NAMD analysis workflow](images/NAMD_CHARMMGUI_workflow.png"A simple simulation workflow starting with CHARMM for setup and NAMD to continue the production simulation")
![Snapshot of NAMD analysis workflow](images/NAMD_workflow.png "A simple NAMD simulation workflow")



The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> ### {% icon details %} NVT theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

> ### {% icon details %} NPT theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

## Step with **setup**
{: #setup}

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **setup** {% icon tool %} with the following parameters:
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

## Step with **minimizer**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **minimizer** {% icon tool %} with the following parameters:
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

## Step with **namd_nvt**
{: #namd_nvt}

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **namd_nvt** {% icon tool %} with the following parameters:
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

## Step with **namd_npt**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **namd_npt** {% icon tool %} with the following parameters:
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

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
