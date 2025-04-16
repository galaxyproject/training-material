---
layout: tutorial_hands_on

title: Phylogenetic analysis for bacterial comparative genomics
zenodo_link: 'https://zenodo.org/records/1'
draft: true
questions:
- to do
- to do
objectives:
- to do
- to do
time_estimation: 1H
key_points:
- to do
- to do
contributions:
  authorship:
  - hchiapello
  - bebatut
  - vloux
edam_ontology:
- topic_0622 # Genomics
- topic_3301 # Microbiology
- topic_0084 # Phylogeny
level: Introductory
tags:
- microgalaxy
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
bibliography section which will automatically be created at the end of the
tutorial.


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

# Galaxy and data preparation

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this analysis
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Rename the history
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
> 3. {% tool [Import](upload1) %} the contig file from [Zenodo]({{ page.zenodo_link }}) or from Galaxy shared data libraries:
>
>    ```
>    {{ page.zenodo_link }}/files/DRR187559_contigs.fasta
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
>
{: .hands_on}

# Tree building

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [FASTTREE](toolshed.g2.bx.psu.edu/repos/iuc/fasttree/fasttree/2.1.10+galaxy1) %} with the following parameters:
>    - *"Aligned sequences file (FASTA or Phylip format)"*: `fasta`
>        - {% icon param-file %} *"FASTA file"*: `output` (Input dataset)
>        - *"Set starting tree(s)"*: `No starting trees`
>    - *"Show advanced options"*: `No`
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


# Conclusion

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.