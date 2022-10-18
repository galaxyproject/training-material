---
layout: tutorial_hands_on

title: Functional annotation of protein sequences
zenodo_link: https://zenodo.org/record/6861851
tags:
  - eukaryote
questions:
  - How to perform functional annotation on protein sequences?
objectives:
  - Perform functional annotation using EggNOG-mapper and InterProScan
time_estimation: 1h
level: Introductory
key_points:
  - EggNOG Mapper compares sequences to a database of annotated orthologous sequences
  - InterProScan detects known motifs in protein sequences
contributions:
  authorship:
    - abretaud
  funding:
    - erasmusplus
subtopic: eukaryote
priority: 6
---


# Introduction

When performing the structural annotation of a genome sequence, you get the position of each gene, but you don't have information about their name of their function. That's the goal of **functional annotation**.

In this short tutorial, we will run the most commonly used tools to perform functional annotation, starting from the predicted protein sequences of a few example genes.

For a more complete view of how this step integrates into a whole genome sequencing and annotation process, you can have a look at the [Funannotate tutorial]({% link topics/genome-annotation/tutorials/funannotate/tutorial.md %}).

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data upload

We will annotate a small set of **protein sequences**. These sequences were predicted from the gene structures obtained in the [Funannotate tutorial]({% link topics/genome-annotation/tutorials/funannotate/tutorial.md %})? Though these sequences from from a fungal species, you can run the same tools on proteins from any organisms, including prokaryotes.

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    https://zenodo.org/api/files/6628a5e4-d6be-47bd-bdaa-f2646112578e/proteins.fasta
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
{: .hands_on}

{% include {{ page.dir }}content.md short=true %}

# Conclusion

Congratulations for reaching the end of this tutorial! Now you know how to perform the functional annotation of a set of protein sequences, using EggNOG mapper and InterProScan.

If you want to collect more functional annotation, you can try to run the {% tool [NCBI BLAST+ blastp](toolshed.g2.bx.psu.edu/repos/devteam/ncbi_blast_plus/ncbi_blastp_wrapper/2.10.1+galaxy2) %} or {% tool [Diamond](toolshed.g2.bx.psu.edu/repos/bgruening/diamond/bg_diamond/2.0.15+galaxy0) %} tools against the UniProt or NR databases (Diamond runs much faster on big datasets). These tools will search for similarities between your protein sequences and the ones already described in big international databases.

Also note that many other more specialised tools exist to collect even more functional annotation, in particular for certain species (prokaryotes forexample), or enzyme/protein families.
