---
layout: tutorial_hands_on

title: Biodiversity data exploration
zenodo_link: https://zenodo.org/record/5930763/files/Reel_life_survey_fish_sample.tabular?download=1
questions:
- Which biological questions are addressed by the tutorial?
- Which bioinformatics techniques are important to know for this type of data?
objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentences describing what a learner should be able to do once they
  have completed the tutorial
- You can use Bloom's Taxonomy to write effective learning objectives
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- contributor1
- contributor2

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

# Title for your first section

Give some background about what the trainees will be doing in the section.
Remember that many people reading your materials will likely be novices,
so make sure to explain all the relevant concepts.

## Title for a subsection
Section and subsection titles will be displayed in the tutorial index on the left side of
the page, so try to make them informative and concise!

# Hands-on Sections
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
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

# Title of the section usually corresponding to a big step in the analysis

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


## Sub-step with **Spatial coordinates anonymization**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Spatial coordinates anonymization](toolshed.g2.bx.psu.edu/repos/ecology/tool_anonymization/tool_anonymization/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: `output` (Input dataset)
>    - *"Select column containing latitudes in decimal degrees"*: `c9`
>    - *"Select column containing longitudes in decimal degrees"*: `c10`
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

## Sub-step with **Homoscedasticity and normality**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Homoscedasticity and normality](toolshed.g2.bx.psu.edu/repos/ecology/ecology_homogeneity_normality/ecology_homogeneity_normality/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: `output` (Input dataset)
>    - *"Select column containing temporal data (year, date, ...)"*: `c11`
>    - *"Select column containing species"*: `c16`
>    - *"Select column containing numerical values (like abundances)"*: `c18`
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

## Sub-step with **Variables exploration**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Variables exploration](toolshed.g2.bx.psu.edu/repos/ecology/ecology_link_between_var/ecology_link_between_var/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: `output` (Input dataset)
>    - *"Variables links exploration"*: `Collinearity between selected numerical variables for each species`
>        - *"Select column containing species"*: `c16`
>        - *"Select columns containing numerical values"*: `c['12', '17', '18']`
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

## Sub-step with **Presence-absence and abundance**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Presence-absence and abundance](toolshed.g2.bx.psu.edu/repos/ecology/ecology_presence_abs_abund/ecology_presence_abs_abund/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: `output` (Input dataset)
>    - *"Variables presence, absence and abundance"*: `Abundance map in the environment `
>        - *"Select column containing latitudes "*: `c9`
>        - *"Select column containing longitudes"*: `c10`
>        - *"What do you study in this analysis ?"*: `fishes`
>        - *"Select column containing taxon "*: `c16`
>    - *"Select column containing abundances "*: `c18`
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

## Sub-step with **Statistics on presence-absence**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Statistics on presence-absence](toolshed.g2.bx.psu.edu/repos/ecology/ecology_stat_presence_abs/ecology_stat_presence_abs/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: `output` (Input dataset)
>    - *"Select a column containing numerical values (such as the abundance) "*: `c18`
>    - *"Select the column of the x-axis : most commonly species"*: `c16`
>    - *"Select column containing locations "*: `c8`
>    - *"Select column containing temporal data (year, date, ...) "*: `c11`
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

## Sub-step with **Local Contributions to Beta Diversity (LCBD)**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Local Contributions to Beta Diversity (LCBD)](toolshed.g2.bx.psu.edu/repos/ecology/ecology_beta_diversity/ecology_beta_diversity/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: `output` (Input dataset)
>    - *"Select column with abundances"*: `c18`
>    - *"Select column with locations"*: `c8`
>    - *"Select column containing taxon"*: `c16`
>    - *"Select column containing dates"*: `c11`
>    - *"Other LCBD : spatialized representation or xy plot."*: `Spatialized representation`
>        - *"Select column containing latitudes in decimal degrees"*: `c9`
>        - *"Select column containing longitudes in decimal degrees"*: `c10`
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
