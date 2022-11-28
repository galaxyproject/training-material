---
layout: tutorial_hands_on

title: "Viewing Cancer Alignments in a Genome Browser"
questions:
objectives:
  - Very quickly navigate around the genome
  - Visualise HTS read alignments
  - Examine SNP calls and structural re-arrangements by eye
  - Tell the difference between germline and somatic variants
key_points:
time_estimation:
contributions:
  authorship: [rdeborjas]
  editing: [shiltemann]
  funding: [bioinformatics-ca,erasmusplus]


---

<!-- TODO add contributors. From bioinf.ca: "This lab is based on the HTS IGV lab originally by
Sorana Morrissy and was updated and modified by Heather Gibling for the Cancer Analysis workshop. " -->

# Introduction

This tutorial will introduce you to the Genome browsers; powerful tools for viewing many kinds of
genomic data, including data for DNA sequencing, RNA sequencing, microarrays, epigenetics, and
copy number alteration. For this tutorial we will use JBrowse, since it has a nice integration
with Galaxy, but many other Genome browsers exist. TODO: link some examples/comparisons

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Data upload

Before we view our data in the Genome browser, let's upload it to Galaxy

> <hands-on-title> Upload data </hands-on-title>
>
> 1. Make sure you have an empty analysis history. Give it a name.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import data from Zenodo
>
>    ```
>     TODO: tumor.bam, normal.bam
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
{: .hands_on}


These files are *alignment* files. This means the raw sequence reads have been mapped to the human reference genome. To learn more about this process of mapping, see our [dedicated tutorial]({% link topics/sequence-analysis/tutorials/mapping/tutorial.md %})


To view alignment files in a genome browser, you usually provide the BAM file, along with an *index* file, which allows the tools to quickly navigate the usually very large alignment files.
However, in our case, Galaxy automatically creates the index file whenever you upload a BAM file, and will supply it to the tools that need it behind the scenes.

TODO: note about test data we are using, scaled down for tutorial reasons



# The JBrowse interface

Let's start by loading the
