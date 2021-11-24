---
layout: tutorial_hands_on

title: Masking repeats with RepeatMasker
zenodo_link: https://zenodo.org/record/5721490
tags:
  - eukaryote
questions:
- How to mask repeats in a genome?
- What is the difference between hard and soft masking?
objectives:
- Use RepeatMasker to soft mask a newly assmbled genome
time_estimation: 3H
level: Introductory
key_points:
- RepeatMasker can be used to soft-mask a genome
- It is an essential first step before running structural annotation pipelines
contributors:
- abretaud
- alexcorm
- lleroi
- r1corre
- stephanierobin

---


# Introduction
{:.no_toc}

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
>    https://zenodo.org/api/files/8458db9c-d517-4c6f-8a91-0ebadbb5a722/genome_raw.fasta
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
{: .hands_on}

# Soft-masking using RepeatMasker

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [RepeatMasker](toolshed.g2.bx.psu.edu/repos/bgruening/repeat_masker/repeatmasker_wrapper/4.1.2-p1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Genomic DNA"*: `genome_raw.fasta` (Input dataset)
>    - *"Repeat library source"*: `DFam (curated only, bundled with RepeatMasker)`
>        - *"Select species name from a list?"*: `Yes`
>            - *"Species"*: `Human (Homo sapiens)`
>    - *"Perform softmasking instead of hardmasking"*: `Yes`
>
{: .hands_on}


# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
