---
layout: tutorial_hands_on

title: Development of statistical analysis and visualization workflows for metagenomic
  amplicon data using the Galaxy framework
level: Intermediate
zenodo_link: https://zenodo.org/records/11281381
questions:
- If we generated amplicon data, how can we analyse it with Galaxy ampvis2 tools?
- How can we visualise amplicon data by using heatmap, ordination plot, boxplot, rarefraction curve or timeseries?
objectives:
- use heatmap workflow to analyse and visualise amplicon data
- or use ordination plot, or boxplot, or rarefraction curve, or timeseries
- use grouped data or grouped data with facets
time_estimation: 3H
key_points:
- using different visualisation methods can present data from other points of view
- with enough metadata the data can be visualised w.r.t. groups and facets 
contributors:
- lenaarenot
- paulzierep

---



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

<!-- This is a comment. -->

**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

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

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    https://zenodo.org/api/records/11281381/files/Galaxy11-[MiDAS_otushort_table.tsv].mothur.axes/content
>    https://zenodo.org/api/records/11281381/files/Galaxy1-[MiDAS_metadata.tsv].tabular/content
>    https://zenodo.org/api/records/11281381/files/Galaxy3-[MiDAS_taxtable.tsv].tabular/content
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

> <details-title> More details about the theory </details-title>
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


## Sub-step with **ampvis2 load**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [ampvis2 load](toolshed.g2.bx.psu.edu/repos/iuc/ampvis2_load/ampvis2_load/2.8.6+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"OTU table"*: `output` (Input dataset)
>    - {% icon param-file %} *"Sample metadata"*: `output` (Input dataset)
>    - {% icon param-file %} *"Taxonomy table"*: `output` (Input dataset)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

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

## Sub-step with **ampvis2 subset samples**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [ampvis2 subset samples](toolshed.g2.bx.psu.edu/repos/iuc/ampvis2_subset_samples/ampvis2_subset_samples/2.8.6+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Ampvis2 RDS dataset"*: `ampvis` (output of **ampvis2 load** {% icon tool %})
>    - {% icon param-file %} *"Metadata list"*: `metadata_list_out` (output of **ampvis2 load** {% icon tool %})
>    - *"Metadata variable"*: ``
>    - *"Metadata value(s)"*: ``
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

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

## Sub-step with **ampvis2 heatmap**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [ampvis2 heatmap](toolshed.g2.bx.psu.edu/repos/iuc/ampvis2_heatmap/ampvis2_heatmap/2.8.6+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Ampvis2 RDS dataset"*: `ampvis` (output of **ampvis2 subset samples** {% icon tool %})
>    - {% icon param-file %} *"Metadata list"*: `metadata_list_out` (output of **ampvis2 subset samples** {% icon tool %})
>    - *"Group samples"*: ``
>    - *"Facet the samples"*: ``
>    - *"The taxonomic level to aggregate the OTUs"*: `Species`
>    - *"Additional taxonomic level(s) to display"*: ``
>    - *"Limit the number of shown taxa"*: `Select a number of taxa to show`
>    - *"Plot the values on the heatmap"*: `Yes`
>    - *"Sort heatmap by most abundant taxa"*: `No`
>    - *"Show functional information about the Genus-level OTUs"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

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


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.