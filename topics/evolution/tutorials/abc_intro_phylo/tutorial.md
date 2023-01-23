---
layout: tutorial_hands_on
title: Phylogenetics - Back to basics
zenodo_link: 'https://zenodo.org/record/6010176'
enable: false
tags:
  - prokaryotes
  - phylogenetics
  - evolution
level: Introductory
questions:
- What information can I get from a phylogenetic tree?
- How do I estimate a phylogeny?
objectives:
- Understand the basic concepts behind phylogenetic trees
- Be able to read and interrogate a phylogeny encountered in the literature
time_estimation: 1H
contributors:
- mcharleston
- adamtaranto

---

# Introduction
Intro text


# Basic concepts: this is a major section
This is a link [introduction to phylogenetics](https://www.ebi.ac.uk/training/online/courses/introduction-to-phylogenetics/).

> <comment-title>This is a comment section</i></comment-title>
>
> Phylogenetics with MTB has some particularities rarely encountered with other organisms.
>
> 1. item 1
>
> 2. item 2
>
{: .comment}


## This is a subsection

Here is an embedded image

![Alignment](./images/MEGA_alignment.png){:width="500"}


## "Hands on task section" Get the data
> <hands-on-title>Obtain your data</hands-on-title>
>
> 1. Make sure you have an empty analysis history. Give it a name.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the following files from [Zenodo](https://zenodo.org/record/6010176) or from the shared data library
>
>    ```
>    https://zenodo.org/record/6010176/files/SNP_alignment.fasta
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
{: .hands_on}

# another section 
Internal link [Resources](#resources) 
An external link[RAxML](https://cme.h-its.org/exelixis/web/software/raxml/)


> <question-title>Exercise 1</question-title>
>
> 1. question?
>
> > <solution-title>1</solution-title>
> >
> > 1. Solution
> >
> {: .solution}
>
{: .question}



> <question-title>Exercise 2</question-title>
>
> example table
>
> | Sample       | Cluster_id | DR profile | Clustering  |
> |--------------|------------|------------|-------------|
> | ERR5987352   | 10         | Pre-MDR    | Clustered   |
> | ERR6362484   | 10         | Pre-MDR    | Clustered   |
> | ERR6362138   | 12         | MDR        | Clustered   |
> | ERR6362156   | 12         | Pre-XDR    | Clustered   |
> | ERR6362253   | 12         | MDR        | Clustered   |
>
> > <solution-title>4</solution-title>
> >
> > Solution
> >
> {: .solution}
>
{: .question}


# Resources
To develop a deeper understanding of phylogenetic trees, there is no better way than estimating phylogenies yourself --- and work through a book on the topic in your own mind's pace.

## Books
- *Phylogenetics in the genomics era*, 2020. An [open access book](https://hal.inria.fr/PGE) covering a variety of contemporary topics.
- *Tree Thinking*, 2013, by David A. Baum & Stacey D. Smith
- *Molecular Evolution*, 2014, by Ziheng Yang

## Useful links
- [MEGA Software](https://megasoftware.net/)
- [Tutorial on how to read a tree, with a virus example](https://artic.network/how-to-read-a-tree.html)
- [Tree Of Life web project](http://tolweb.org)
- [Phylogenetic Inference in the Stanford Encyclopedia](https://plato.stanford.edu/entries/phylogenetic-inference/)
