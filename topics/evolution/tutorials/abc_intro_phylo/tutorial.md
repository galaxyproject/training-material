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


> <details-title>Further reading</details-title>
> Here is a link to the PLoS article on Galaxy tutorials:
> - [Galaxy Training: A powerful framework for teaching!](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010752)
{: .details}

## Overview

This tutorial has the following structure:

- Introduction and motivation: why is phylogenetic inference important?
- A general overview of phylogenetic inference: from sequence data onward.
- Obtaining the data for this tutorial + exercise
- Sequence alignment (including manual methods, automatic methods, complexity issues / heuristics) + exercise
- Distances based on sequence alignment
- The Neighbor-Joining method & FastME2.0 (https://doi.org/10.1093/molbev/msv150) **needs an update** **not sure that there is a FastME implementation in Galaxy; there is an R package though**
- Building the first tree (on Galaxy)
- Models of sequence evolution: from the sublime to the ridiculous
- Phylogenetic Networks (**SplitsTree needs install**), Neighbor-Net
- Assessing the quality of the tree(s): Bootstrapping, branch lengths; conflict in the networks
- Maximum Likelihood with IQTree


## Motivation

There are many ways in which we can use phylogenetic analyses: from the most fundamental understanding of the evolutionary relationships that exist between a set of species, as in Charles Darwin's famous sketch in Origin of Species

![IThink](./images/Darwin_tree.png){:width="400"}

**needs a reference**

all birds
![AllBirds](./images/nature11631-f2.jpg){:width="400"}
(from from Jetz *et al.* 2012, Nature (491):444–448)

and much bigger projects across all of life:

![UnderstandinEvolTree](./images/nmicrobiol201648-f1.jpg){:width="500"}

(from Understanding Evolution. 2019. University of California Museum of Paleontology. 4th November 2019; http://evolution.berkeley.edu)

Aside from fundamental understanding, there are other strong motivators for inferring phylogenetic relationships:

- Designing vaccines, for example for SARS-CoV2 and influenza;
- Measuring phylogenetic diversity for guiding conservation efforts;
- Understanding coevolution: around 70% of emergent human diseases have come from other species;
- Dating major evolutionary events, to study the effects of environmental change on different species.

> <comment-title>Gene trees, species trees</comment-title>
> 
> It's worth noting that getting the phylogeny from a set of genes -- what we often call a *gene tree* -- might *not* > give us the phylogeny of the species that house those genes, *even if we get everything right!*
>
> This happens because there are other processes that can lead to the so-called gene tree not being the same as the species tree:
>	- lateral gene transfer events
>	- gene duplication
>	- gene loss and incomplete lineage sorting
>	- recombination
>
> This could send us off down a very deep and difficult rabbit-hole: that of the gene-tree / species tree problem: but today we will work under the assumption (which is reasonable in this case) that the gene tree will reflect the species relationships.
>
> **The situation where gene trees and species trees differ is often called the "gene tree / species tree reconciliation problem", and while it is very interesting and important, it is beyond the scope of this tutorial.  The interested reader is directed to...**
>
{: .comment}

# Basic Methodology

First and foremost, **phylogenetic inferences is a statistic estimation process.**

It is not generally possible to prove that any tree inferred is *correct* -- since we cannot go back in time and observe speciation events.
One obvious consequence of this is that different estimates of the phylogenetic tree relating a given set of species may differ, even if no errors were made.  
Finding an optimal tree is hard!

So, how do we do it?

There are several ways to estimate a tree, such as:

1. Go with what we think is the case already (not recommended!)
2. Attempt to build a tree based on similarity and dissimilarity, with such tools as Neighbor-Joining (NJ) or FastME;
3. Choose some kind of score function such as Parsimony, or Maximum Likelihood to potential trees and find the best one
4. Something else entirely (networks? inference based on their parasites)


This is a link [introduction to phylogenetics](https://www.ebi.ac.uk/training/online/courses/introduction-to-phylogenetics/).

> <comment-title>Common Evolutionary Assumptions used in Phylogenetic Estimation</comment-title>
> 
>  These may help your understanding of why things are done this way:
>
> 1. Evolution is “memoryless” (which means we can use the powerful mathematics of Markov processes);
> 
> 2. Phylogenetic relationships can be correctly represented by a tree! (This isn't always assumed, but it is very common.)
> 3. The *Molecular clock* assumption: sequences in a clade evolve at about the same rate as each other (this is easily tested);
> 4. Lineages don’t interact – once they have speciated, they are independent of each other.
>
> We do know that these assumptions do not always hold!  For instance, there is commonly variation in evolutionary rate between lineages, but if this variation is not significant, we can ignore it and use simpler models, to better leverage the phylogenetic information there is in the data.
> We also know that biological lineages *do* interact with each other
> 
{: .comment}

## Challenges


> <comment-title>This is a comment section</comment-title>
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
