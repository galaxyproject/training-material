---
title: Kraken2 and the k-mer approach for taxonomy classification
area: format
box_type: details
layout: faq
contributors: [bebatut]
---

In the $$k$$-mer approach for taxonomy classification, we use a database containing DNA sequences of genomes whose taxonomy we already know. On a computer, the genome sequences are broken into short pieces of length $$k$$ (called $$k$$-mers), usually 30bp.

**Kraken** examines the $$k$$-mers within the query sequence, searches for them in the database, looks for where these are placed within the taxonomy tree inside the database, makes the classification with the most probable position, then maps $$k$$-mers to the lowest common ancestor (LCA) of all genomes known to contain the given $$k$$-mer.

![Kraken2]({{site.baseurl}}/topics/metagenomics/faqs/images/kmers-kraken.jpg "Kraken sequence classification algorithm. To classify a sequence, each k-mer in the sequence is mapped to the lowest common ancestor (LCA, i.e. the lowest node) of the genomes that contain that k-mer in the database. The taxa associated with the sequence's k-mers, as well as the taxa's ancestors, form a pruned subtree of the general taxonomy tree, which is used for classification. In the classification tree, each node has a weight equal to the number of k-mers in the sequence associated with the node's taxon. Each root-to-leaf (RTL) path in the classification tree is scored by adding all weights in the path, and the maximal RTL path in the classification tree is the classification path (nodes highlighted in yellow). The leaf of this classification path (the orange, leftmost leaf in the classification tree) is the classification used for the query sequence. Source: {% cite Wood2014 %}")

__Kraken2__ uses a compact hash table, a probabilistic data structure that allows for faster queries and lower memory requirements. It applies a spaced seed mask of _s_ spaces to the minimizer and calculates a compact hash code, which is then used as a search query in its compact hash table; the lowest common ancestor (LCA) taxon associated with the compact hash code is then assigned to the k-mer.

You can find more information about the __Kraken2__ algorithm in the paper [_Improved metagenomic analysis with Kraken 2_](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0).