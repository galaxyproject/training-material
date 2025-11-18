---
layout: tutorial_hands_on

title: Bacterial genome quality control
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

# Quality control

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> <details-title> More details about the theory </details-title>
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [checkm2](toolshed.g2.bx.psu.edu/repos/iuc/checkm2/checkm2/1.0.2+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input MAG/SAG datasets"*: Input dataset collection
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

# Dereplication

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [dRep dereplicate](toolshed.g2.bx.psu.edu/repos/iuc/drep_dereplicate/drep_dereplicate/3.5.0+galaxy1) %} with the following parameters:
>    - {% icon param-collection %} *"Genomes"*: `output` (Input dataset collection)
>    - *"Genome quality filtering"*: `Run checkM`
>    - In *"Genome comparison and clustering"*:
>        - *"Steps in genome comparison"*: `Default: Run MASH clustering and a secondary clustering`
>            - *"Algorithm for secondary clustering comparisons"*: `ANImf: Align whole genomes with nucmer; filter alignment; compare aligned regions - RECOMMENDED`
>    - *"Select outputs"*: ``
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