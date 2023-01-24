---
layout: tutorial_hands_on
title: Phylogenetics - Back to basics
zenodo_link: 'https://tinyurl.com/phylo-trees-1-data'
enable: false
tags:
  - phylogenetics
  - evolution
level: Introductory
questions:
- What information can I get from a phylogenetic tree?
- How do I estimate a phylogeny?
objectives:
- Understand the basic concepts behind phylogenetic trees
- Be able to read and interrogate a phylogeny encountered in the literature
time_estimation: 2H
contributors:
- mcharleston
- adamtaranto

---

# Introduction
Intro text

## Overview

This tutorial has the following structure:

- Introduction and motivation: why is phylogenetic inference important?
- A general overview of phylogenetic inference: from sequence data onward.
- Obtaining the data for this tutorial + exercise
- Sequence alignment (including manual methods, automatic methods, complexity issues / heuristics) + exercise
- Distances based on sequence alignment
- The Neighbor-Joining method and others (BioNJ and the other one that's better) **needs an update**
- Building the first tree (on Galaxy)
- Models of sequence evolution: from the sublime to the ridiculous
- Phylogenetic Networks (**SplitsTree needs install**), Neighbor-Net
- Assessing the quality of the tree(s): Bootstrapping, branch lengths; conflict in the networks
- Maximum Likelihood with IQTree


## Motivation

There are many ways in which we can use phylogenetic analyses: from the most fundamental understanding of the evolutionary relationships that exist between a set of species, as in Charles Darwin's famous sketch in Origin of Species:

![IThink](./images/Darwin_tree.png){:width="400"}

# Basic Methodology
This is a link [introduction to phylogenetics](https://www.ebi.ac.uk/training/online/courses/introduction-to-phylogenetics/).

> <comment-title>This is a comment section</i></comment-title>
>
> This is a comment section
>
> 1. item 1
>
> 2. item 2
>
{: .comment}

## Data upload

Background on the data used in this workshop.

## Get the data
> <hands-on-title>Obtain your data</hands-on-title>
>
> 1. Make sure you have an empty analysis history. Give it a name.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the following files from [Zenodo](https://tinyurl.com/phylo-trees-1-data) or from the shared data library
>     Note: Current link is to google drive. Update to Zenodo for final release.
>
>    ```
>    exon7-unaligned.fst
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
{: .hands_on}

**IMAGE HERE: View of unaligned fasta sequences**

# Sequence Alignment

Methods
Gaps
Scoring an alignment

## Aligning sequences with MAFFT

> <hands-on-title>Sequence alignment with MAFFT</hands-on-title>
>
> 1. Align seqs with settings: xxx
> 2. View alignment
>
{: .hands_on}

Here is an embedded image

![Alignment](./images/MEGA_alignment.png){:width="500"}


# Distance-based phylogenetic inference

Calculating distances from an alignment

## Building a Neighbor-Joining Tree 

> <hands-on-title>Build a Neighbour-Joining Tree with Splitstree</hands-on-title>
>
> 1. Step 1
> 2. Step 2
>
{: .hands_on}

**IMAGE HERE: NJ Tree image**


# Phylogenetic Networks

Intro to phylogenetic networks as an alternative to trees

## Building a Neighbor-Net phylogenetic network

> <hands-on-title>Build a Neighbor-Net with Splitstree</hands-on-title>
>
> 1. Step 1
> 2. Step 2
>
{: .hands_on}


**IMAGE HERE: Neighbour net image**

# Basics of Maximum Likelihood

Background on ML. How it works. Software. Bootstrap values.

Info on models


##  Estimating a Maximum Likelihood tree 

> <hands-on-title>Estimating a Maximum Likelihood tree with RaXML and IQTree</hands-on-title>
>
> 1. Step 1
> 2. Step 2
>
{: .hands_on}


**IMAGE HERE: ML tree with bootstrap values**


# Possible new section on formatting and exporting a pretty tree



# Summary 

- Phylogeny provides the statistical framework that is essential to comparing biological organisms
- There are so many trees, choosing the best one is very hard.
- There are many data options but one of the best is to use molecular
sequences.
- Optimality criteria (e.g., MP, ML) help us decide which trees are “good” – by how well they explain the data.
- We can search tree space for medium-sized problems with branch-and-bound, and bigger problems with heuristics.
- Bayesian analysis is a way of incorporating prior knowledge in ML analyses.
- Trees can be assessed for robustness by resampling. testing.


# Troubleshooting

Here are a few things that can still catch us out:

  - **Long Branch Attraction (LBA):** 
    Be wary of long branches that come out together in the estimated phylogeny: this can be the result pairs of sequences that are very different from the rest, so match each other “by chance” more than they match the rest.
    
    **Fix:** break up these long branches by adding in some taxa that are closely related to one or the other; remove one long branch at a time to see where the remaining one fits best; consider other methods that are more robust to LBA.

  - **Very “Gappy” Sequences:** 
    Sequences that are hard to align might contain many gaps and many equally “good” alignments.
  
    **Fix:** Try different multiple alignment programs; consider using “alignment-free” methods such as k-mer distances; remove very problematic regions using programs such as GBlocks.

  - **Low resolution:** 
    Low bootstrap support or lots of conflict in a network.
    
    **Fix:** Look at which sites support which splits (internal branches); consider sliding window approaches or check that your sequences don’t span regions with different selection pressures; consider using PartitionFinder or similar methods to work out which sets of sites have similar evolutionary dynamics.

  - **The gene trees are different!**
    
    **Fix:** Well they might not need fixing: it might just be that the genes’ evolutionary histories aren’t the same as those of the species that host them. Look at all the gene trees and see what other events might have led to the differences between them.

  - **I can’t find an outgroup!**

    **Fix:** Consider mid-point rooting: it is in most cases pretty good.


# Markdown hints
Internal link [Resources](#resources) 
An external link[RAxML](https://cme.h-its.org/exelixis/web/software/raxml/)

> example table
>
> | Sample       | Cluster_id | DR profile | Clustering  |
> |--------------|------------|------------|-------------|
> | ERR5987352   | 10         | Pre-MDR    | Clustered   |
> | ERR6362484   | 10         | Pre-MDR    | Clustered   |



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
> A question
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
