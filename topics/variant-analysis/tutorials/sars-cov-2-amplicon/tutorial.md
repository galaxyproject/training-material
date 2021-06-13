---
layout: tutorial_hands_on

title: SARS-CoV-2 ARTIC sequence analysis
zenodo_link: 'https://zenodo.org/record/4899860'
questions:
- Which biological questions are addressed by the tutorial? Change this here
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
- pvanheus
- mudiboevans

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/Introduction.md %}

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
>https://zenodo.org/record/4899860/files/ERR4970105_1.fastq.gz
>https://zenodo.org/record/4899860/files/ERR4970105_2.fastq.gz
>https://zenodo.org/record/4899860/files/ERR4970106_1.fastq.gz
>https://zenodo.org/record/4899860/files/ERR4970106_2.fastq.gz
>https://zenodo.org/record/4899860/files/ERR4970107_1.fastq.gz
>https://zenodo.org/record/4899860/files/ERR4970107_2.fastq.gz
>https://zenodo.org/record/4899860/files/MN908947_3_Wuhan-Hu-1.fasta
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets (to get ride of the `http://zenodo.org...` prefix)
>
> 4. Check that the datatype
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to the sample ID
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
> 6. Convert the FASTQ datasets to a [List Paired collection](topics/galaxy-interface/tutorials/collections/tutorial.html).
>
{: .hands_on}

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/rbu.md %}

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/fastp.md %}

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/bwa_mem.md %}

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/samtools_stats.md %}

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/bam_filter.md %}

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/bamqc.md %}

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/ivar_trim.md %}

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/Variant_annotation.md %}

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/Consensus_genome_characterization.md %}

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/Lineage_designation.md %}

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


## Sub-step with **My Tool**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [My Tool]() %} with the following parameters:
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