---
layout: tutorial_hands_on

title: Identification of yeasts in a beer
zenodo_link: https://doi.org/10.5281/zenodo.6620778
level: Introductory
questions:
- How can yeast strains in a beer sample be identified?
- How can we process metagenomic data sequenced using Nanopore?
objectives:
- Inspect metagenomics data
- Run metagenomics tools
- Identify yeast strains contained in a sequenced beer sample using DNA
- Visualize the microbiome community of a beer sample
time_estimation: 1H
key_points:
- Galaxy provides various tools for data analysis, dataset alteration, and visualization
- Inputting correct values for the parameters of proper tools is important
- Profiling...
- Microbiome
tags:
- nanopore
- beer
- citizen science
- metagenomics
contributors:
- plushz
- chensy96
- bebatut
---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

What is a microbiome? There are collections of small living creatures.
These small creatures are called bacteria and they are everywhere. In our gut,
in the soil, on vending machines, and even inside the beer. Most of these bacteria are
actually very good for us, but some can make us very ill.

Bacteria come in different shapes and sizes, but they have the same components.
One crucial component is the DNA, the blueprint of life. The DNA encodes the
shape and size and many other characteristics unique to a bacterial species. Because of
the encoding information the DNA can be used to identify what kind of bacteria
the DNA is from. Therefore, within a metagenomic sample, e.g. form soil, gut, or beer, one can
identify what kind of species are inside the sample.

In this tutorial, we will use data generated via the [BeerDEcoded project](https://streetscience.community/projects/beerdecoded/).

> ### {% icon comment %} The BeerDEcoded project
>
> The BeerDEcoded project are workshops organized with and for schools and general
> audience, to introduce biology and genomic science. People will learn about
> DNA, sequencing technologies, bioinformatics, open science, how these technologies and concepts are
> applied and how they are impacting their daily life.
>
> The 1-2 days continuous (or divided over several days) workshops include the following steps:
> 1. Extract yeasts and their DNA from beer bottle,
> 2. Sequence the extracted DNA using a MinION sequencer to obtain the sequence of bases/nucleotides (A, T, C and G) for each DNA fragment in the sample,
> 3. Analyze the sequenced data in order to know which organisms this DNA is from
>
> ![The image represents a BeerDEcoded workshop. On the left, there is a beer glass. An arrow goes from the bottle to DNA with "Extraction" written on the below. An arrow goes from DNA to DNA sequences with "Sequencing" written on the below. An arrow goes from the DNA sequences to Yeasts with "Data analysis" written on the below](./images/beerprocess.png)
>
{: .comment}

DNA of yeasts in a bottle of La Trappe beer has been extracted and sequenced using a
MinION to obtain sequences of DNA of the extracted yeasts. Now, for each obtained sequence, we would like to identify the yeast species to which it belongs, and thereby outline the diversity of organisms (the microbiome community) in the beer sample.

To get this information, we need to process the sequenced data in a few steps:
1. Check the quality of the data
2. Assign taxonomic label, i.e. assigh 'species' to the sequences
3. Visualize the species distribution

This type of data analysis requires running several bioinformatics tools and
usually requires a computer science background. [Galaxy](https://galaxyproject.org/) is
an open-source platform for data analysis that enables users to use bioinformatics
tools through its graphical web interface, accessible via any Web browser.

So, in this tutorial, we will use Galaxy to extract and visualize the community
of yeasts from a beer bottle.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Prepare Galaxy and data

First of all, this tutorial will get you hands on with some basic Galaxy tasks, including creating a history and importing data.

## Get familiar with Galaxy

> ### {% icon hands_on %} Hands-on: Open Galaxy
>
> 1. Open your favorite browser (Chrome, Safari or Firefox as your browser, not Internet Explorer!)
> 2. Create a Galaxy account if you do not have one
>
>    {% snippet faqs/galaxy/galaxy_creating_an_account.md %}
>
{: .hands_on}

The Galaxy homepage is divided into three panels:
* Tools on the left
* Viewing panel in the middle
* History of analysis and files on the right

![Galaxy interface screenshot showing history panel on the right, tools panel on the left, and main panel at the center](./images/galaxy_interface.png "The Galaxy interface")

The first time you use Galaxy, there will be no files in your history panel.

Any analysis should get its own Galaxy history. So let's start by creating a new one:

> ### {% icon hands_on %} Hands-on: Prepare the Galaxy history
>
> 1. Create a new history for this analysis
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Rename the history
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
{: .hands_on}

## Get data

Before we can begin any Galaxy analysis, we need to upload the input data: FASTQ files

> ### {% icon hands_on %} Hands-on: Upload your dataset
>
> 1. Import the sequenced data
>
>    - Option 1 [{% icon video %}](https://youtu.be/FFCDx1rMGAQ): Your own local data using **Upload Data** (recommended for 1-10 datasets).
>
>      {% snippet faqs/galaxy/datasets_upload.md %}
>
>    - Option 2: From Zenodo, an external server via URL
>
>      ```text
>      {{ page.zenodo_link }}
>      ```
>
>      {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    Your uploaded file is now in your current history. When the file has uploaded to Galaxy, it will turn green. But, what is this file?
>
> 2. Click on the {% icon galaxy-eye %} (eye) icon next to the dataset name to look at the file content
>
{: .hands_on}

The contents of the file will be displayed in the central Galaxy panel.

This file contains the sequences, also called **reads**, of DNA, *i.e.* succession of nucleotides, for all fragments from the yeasts in the beer, in FASTQ format.

{% snippet topics/sequence-analysis/faqs/fastq.md %}

# Data quality

## Assess data quality

Before starting to work on our data, it is necessary to assess its quality. This is an essential step if we aim to obtain a **meaningful downstream analysis**.

**FastQC** is one of the most widely used tools to **check the quality** of data generated by High Throughput Sequencing (HTS) technologies.

> ### {% icon hands_on %} Hands-on: Quality check
>
> 1. {% tool [FASTQC](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.73+galaxy0) %} with the following parameters
>    - {% icon param-file %} *"Raw read data from your current history"*: `Reads`
>
> 2. Inspect the generated HTML file
>
{: .hands_on}

> ### {% icon question %} Questions
>
> Given the Basic Statistics table on the top,
> 1. How many sequences are in the FASTQ file?
> 2. How long are the sequences?
>
> > ### {% icon solution %} Solution
> >
> > 1. There are 1876 sequences.
> > 2. The sequences range from 130 nucleotides to 2327 nucleotides. Not all sequences have then the same length.
> >
> {: .solution}
>
{: .question}

**FastQC** provides information on various parameters, such as the range of quality values across all bases at each position:

![FastQC Per base sequence quality with scores below 200](./images/fastqc_1.png "Per base sequence quality")

We can see that the quality of our sequencing data grows after the first few bases, stays around a score of 18, which is a relatively low value compared to other sequencing technologies, and then decreases again at the end of the sequences.

For more detailed information about the other plots in the FASTQC report, check out our [dedicated tutorial]({% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %}).

## Improve the dataset quality

In order to improve the quality of our data, we will use two tools:
- **porechop** ({% cite rrwick2017 %}) to remove adapters that were added for sequencing and chimera (contaminant)
- **fastp** ({% cite Chen_2018 %}) to filter sequences with low quality scores (below 9)

> ### {% icon hands_on %} Hands-on: Improve the dataset quality
>
> 1. {% tool [Porechop](toolshed.g2.bx.psu.edu/repos/iuc/porechop/porechop/0.2.3) %} with the following parameters:
>    - {% icon param-file %} *"Input FASTA/FASTQ"*:  `Reads`
>    - *"Output format for the reads"*: `fastq`
>
> 2. {% tool [fastp](toolshed.g2.bx.psu.edu/repos/iuc/fastp/fastp/0.23.2+galaxy0) %} with the following parameters:
>    - *"Single-end or paired reads"*: `Single-end`
>    - {% icon param-file %} *"Dataset collection"*: output of **Porechop**
>    - In *"Adapter Trimming Options"*:
>        - *"Disable adapter trimming"*: `Yes`
>    - In *"Filter Options"*:
>        - In *"Quality filtering options"*:
>            - *"Qualified quality phred"*: `9`
>    - In *"Read Modification Options"*:
>        - *"PolyG tail trimming"*: `Disable polyG tail trimming`
>
> 3. Inspect the HTML report of **fastp** to see how the quality has been improved
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. How many sequences are there before filtering? Is it the same number as in FASTQC report?
> 2. How many sequences are there after filtering? How many sequences have then been removed by filtering?
> 3. What is the mean length before filtering? And after filtering?
>
> > ### {% icon solution %} Solution
> >
> > 1. There are 1,869 reads before filtering. The number is lower than in the FASTQC report. Some reads may have been discarded via Porechop
> > 2. There are 1,350 reads after filtering. So the filtering step has removed $$1869-1350 = 519$$ sequences.
> > 3. The mean length is 314 nucleotide before filtering and 316bp after filtering.
> {: .solution}
>
{: .question}

# Assign taxonomic classifications

One of the main aim in microbiome data analysis is to identify the organisms sequenced. For that we try to identify the taxon to which each individual reads belong.

{% snippet topics/metagenomics/faqs/taxon.md %}

Taxonomic assignment or classification is the process of assigning an **Operational Taxonomic Unit** (OTUs, that is, groups of related individuals / taxon) to sequences. To assign an OTU to a sequence it is compared against a database, but this comparison can be done in different ways, with different bioinformatics tools. Here we will use **Kraken2** ({% cite wood2019improved %}).

> ### {% icon details %} Kraken2 and the k-mer approach for taxonomy classification
>
> In the $$k$$-mer approach for taxonomy classification, we use a database containing the DNA sequences of genomes we know the taxonomy. The genome database is broken into pieces of length $$k$$ (called $$k$$-mers), usually 30bp.
>
> **Kraken** examines the $$k$$-mers within the query sequence, searches for them in the database, looks for where these are placed within the taxonomy tree inside the database and make the classification with the most probable position.  the maps $$k$$-mers to the lowest common ancestor (LCA) of all genomes known to contain the given $$k$$-mer.
>
> ![Kraken2](../../images/metagenomics-nanopore/kmers-kraken.jpg "Kraken sequence classification algorithm. To classify a sequence, each k-mer in the sequence is mapped to the lowest common ancestor (LCA, i.e. lowest node) of the genomes that contain that k-mer in a database. The taxa associated with the sequence's k-mers, as well as the taxa's ancestors, form a pruned subtree of the general taxonomy tree, which is used for classification. In the classification tree, each node has a weight equal to the number of k-mers in the sequence associated with the node's taxon. Each root-to-leaf (RTL) path in the classification tree is scored by adding all weights in the path, and the maximal RTL path in the classification tree is the classification path (nodes highlighted in yellow). The leaf of this classification path (the orange, leftmost leaf in the classification tree) is the classification used for the query sequence. Source: {% cite Wood2014 %}")
>
{: .details}

<!-- Add something about database -->

> ### {% icon hands_on %} Hands-on: Kraken2
>
> 1. {% tool [Kraken2](toolshed.g2.bx.psu.edu/repos/iuc/kraken2/kraken2/2.0.8_beta+galaxy0) %} with the following parameters:
>    - *"Single or paired reads"*: `Single`
>        - {% icon param-file %} *"Input sequences"*: Output of **fastp**
>    - *"Print scientific names instead of just taxids"*: `Yes`
>    - In *"Create Report"*:
>        - *"Print a report with aggregrate counts/clade to file"*: `Yes`
>        - *"Format report output like Kraken 1's kraken-mpa-report"*: `Yes`
>    - *"Select a Kraken2 database"*: `Prebuilt Refseq indexes: PlusPF`
>
> 2. Inspect the report file
{: .hands_on}

The Kraken report in MPA format is a tabular file with 2 columns:

```
d__Eukaryota	756
d__Eukaryota|k__Metazoa	402
d__Eukaryota|k__Metazoa|p__Chordata	402
d__Eukaryota|k__Metazoa|p__Chordata|c__Mammalia	402
d__Eukaryota|k__Metazoa|p__Chordata|c__Mammalia|o__Primates	402
d__Eukaryota|k__Metazoa|p__Chordata|c__Mammalia|o__Primates|f__Hominidae	402
d__Eukaryota|k__Metazoa|p__Chordata|c__Mammalia|o__Primates|f__Hominidae|g__Homo	402
d__Eukaryota|k__Metazoa|p__Chordata|c__Mammalia|o__Primates|f__Hominidae|g__Homo|s__Homo sapiens	402
d__Eukaryota|k__Fungi	342
d__Eukaryota|k__Fungi|p__Ascomycota	341
d__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes	336
...
```

1. The first column lists clades, ranging from taxonomic domains (Bacteria, Archaea, etc.) through species.

    The taxonomic level of each clade is prefixed to indicate its level:
    - Domain: `d__`
    - Kingdom: `k__`
    - Phylum: `p__`
    - Class: `c__`
    - Order: `o__`
    - Family: `f__`
    - Genus: `g__`
    - Species: `s__`

2. The second column gives the number of reads assigned to the clade rooted at that taxon


> ### {% icon question %} Questions
>
> 1. How many reads are assigned Homo sapiens?
> 2. What is the most found taxon?
> 3. Which fungi species are found?
>
> > ### {% icon solution %} Solution
> >
> > 1. 402 reads are assigned to **Homo sapiens**
> > 2. Eukaryota is the most found taxon
> > 3. Saccharomyces cerevisiae, Saccharomyces paradoxus, Saccharomyces eubayanus, Kluyveromyces marxianus, Sugiyamaella lignohabitans
> >
> {: .solution}
>
{: .question}

# Visualize the community

Once we have assigned the corresponding taxa to the sequences, the next step is to properly visualize the data: visualize the diversity of taxons at different levels.

To do that, we will use the tool **Krona** ({% cite Ondov_2011 %}). But before that, we need to adjust the format of the data output from Kraken2.

> ### {% icon hands_on %} Hands-on: Prepare dataset for Krona
>
> 1. {% tool [Reverse](toolshed.g2.bx.psu.edu/repos/iuc/datamash_reverse/datamash_reverse/1.1.0) %} with the following parameters:
>    - {% icon param-file %} *"Input tabular dataset"*: output of **Kraken2** {% icon tool %})
> 2. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_line/1.1.2) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `out_file` (output of **Reverse** {% icon tool %})
>    - In *"Replacement"*:
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `\|`
>            - *"Replace with:"*: `\t`
{: .hands_on}

**Krona** creates an interactive report that allows hierarchical data (like taxonomy) to be explored with zooming, multi-layered pie charts. With this tool, we can easily visualize the composition of the microbiome communities and also compare how the populations of microorganisms are modified between different samples.

> ### {% icon hands_on %} Hands-on: Krona pie chart
>
> 1. {% tool [Krona pie chart](toolshed.g2.bx.psu.edu/repos/crs4/taxonomy_krona_chart/taxonomy_krona_chart/2.7.1) %} with the following parameters:
>    - *"What is the type of your input data"*: `Tabular`
>        - {% icon param-file %} *"Input file"*:  output of **Replace Text** {% icon tool %}
>
> 2. Inspect the generated file
{: .hands_on}

Let's take a look at the result.

<!--
<iframe id="krona" src="krona_all.html" frameBorder="0" width="100%" height="900px"> ![Krona](./images/krona.png) </iframe>

Add data interpretation and questions to understand the data

> ### {% icon question %} Questions
>
> 1.
>
> > ### {% icon solution %} Solution
> >
> > 1.
> >
> {: .solution}
>
{: .question}

--->

# Conclusion
{:.no_toc}


## (Optional) Sharing your history

One of the most important features of Galaxy comes at the end of an analysis: sharing your histories with others so they can review them.

> ### {% icon hands_on %} Hands-on: Share history and workflow
>
> 1. Click on the {% icon galaxy-gear %} (gear) symbol in the history panel
> 2. Select `Share or Publish`
> 3. Make your history accessible via a link
> 4. Copy the link
>
{: .hands_on}

> ### {% icon comment %} The 3 ways to share in Galaxy
>
> 1. **Make accessible via Link**
>
>     This generates a link that you can give out to others. Anybody with this link will be able to view your history.
>
> 2. **Publish History**
>
>     This will not only create a link, but will also publish your history. This means your history will be listed under `Shared Data → Published Histories` in the top menu.
>
> 3. **Share with Individual Users**
>
>     This will share the history only with specific users on the Galaxy instance.
>
{: .comment}




