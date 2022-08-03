---
layout: tutorial_hands_on

title: 'Genome Assembly Quality Control'
zenodo_link: 'https://zenodo.org/record/6947782'
tags:
  - assembly
  - quality control
questions:
- Is my genome assembly ready for annotation and/or scaffolding?
objectives:
- Assess assembly quality
time_estimation: 2h
level: Intermediate
key_points:
- Quast, BUSCO and Merqury make it easy to assess the quality of an assembly
- Post-processing can be necessary (purging, scaffolding) before annotation
contributors:
- abretaud
- alexcorm
- r1corre
- lleroi
- stephanierobin

follow_up_training:
 - type: internal
   topic_name: genome-annotation
   tutorials:
     - repeatmasker

---

# Introduction
{:.no_toc}

In this tutorial, we will assess the genome assembly quality of 2 assemblies generated with Hifiasm and Flye using PacBio HiFi reads of a species of fungi, *Saccharomyces cerevisiae* INSC1019 and compare the results with the actual reference genome [*Saccharomyces cerevisiae* S288C](https://www.ncbi.nlm.nih.gov/genome/?term=Saccharomyces%20cerevisiae).

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get data

We will use long quality reads sequencing data from PacBio HiFi sequencing of *Saccharomyces cerevisiae* INSC1019 genome. This data is a subset of data from ENA repository [SRR13577847](https://www.ebi.ac.uk/ena/browser/view/SRR13577847?show=reads). We will also use the reference genome assembly downloaded from the [NCBI website](https://www.ncbi.nlm.nih.gov/genome/?term=Saccharomyces%20cerevisiae) and we will use it as a comparison with our own assemblies.

## Get data from Zenodo

> ### {% icon hands_on %} Hands-on: Data upload from Zenodo
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo](https://zenodo.org/record/5702408)
>
>    ```
>    https://zenodo.org/record/6947782/files/GCA_000146045.2_genomic.fna
>    https://zenodo.org/record/6947782/files/Scerevisiae-INSC1019.flye.30x.fa
>    https://zenodo.org/record/6947782/files/Scerevisiae-INSC1019.hifiasm.30x.fa
>    https://zenodo.org/record/6947782/files/SRR13577847_subreads.30x.fastq.gz
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Check that the datatype is `fastqsanger.gz` for SRR13577847_subreads.30x.fastq.gz
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

# Evaluation of assembly contiguity with ***Quast*

> ### {% icon hands_on %} Hands-on: assembly evaluation with Quast
>
> 1. {% tool [Quast](toolshed.g2.bx.psu.edu/repos/iuc/quast/quast/5.0.2+galaxy1) %} with the following parameters:
>    - *"Use customized names for the input files?"*: `No, use datasets names`
>    - {% icon param-files %} *"Contigs/scaffolds file"*: `Scerevisiae-INSC1019.flye.30x.fa` and `Scerevisiae-INSC1019.hifiasm.30x.fa`
>    - *"Reads options"*: `Pacbio SMRT reads`
>        - {% icon param-file %} *"FASTQ file"*: `SRR13577847_subreads.30x.fastq.gz`
>    - *"Type of assembly"*: `Genome`
>        - *"Use a reference genome?"*: `Yes`
>            - *"Reference genome"*: `GCA_000146045.2_genomic.fna`
>    - *"Is genome large (>100Mpb)?"*: `No`
>
> 2. Rename the HTML report as `QUAST initial report`
>
{: .hands_on}

# Evaluation of assembly completness

## Core genes completness with **BUSCO**

> ### {% icon hands_on %} Hands-on: assessing assembly completness with BUSCO
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.0.0+galaxy0) %} with the following parameters:
>    - {% icon param-files %} *"Sequences to analyse"*:`Scerevisiae-INSC1019.flye.30x.fa`, `Scerevisiae-INSC1019.hifiasm.30x.fa` and `GCA_000146045.2_genomic.fna`
>    - *"Mode"*: `Genome assemblies (DNA)`
>        - *"Use Augustus instead of Metaeuk"*: `Use Metaeuk`
>    - *"Auto-detect or select lineage?"*: `Select lineage`
>       - *"Lineage"*: `Saccharomycetes`
>    - *"Which outputs should be generated"*: `short summary text` and `summary image`
>
> 2. Rename the summary as `BUSCO initial report`
>
{: .hands_on}

## k-mer based assembly evaluation with **Merqury**


> ### {% icon hands_on %} Hands-on: Generate k-mers count distribution
>
> 1. {% tool [Meryl](toolshed.g2.bx.psu.edu/repos/iuc/meryl/meryl/1.3+galaxy2) %} with the following parameters:
>    - *"Operation type selector"*: `Count operations`
>        - *"Count operations"*: `Count: count the ocurrences of canonical k-mers`
>        - {% icon param-file %} *"Input sequences"*: `SRR13577847_subreads.30x.fastq.gz`
>        - *"K-mer size selector"*: `Estimate the best k-mer size`
>            - "*Genome size*": `12000000`
>
> 2. Rename it as `Meryldb from HiFi`.
>
{: .hands_on}

> ### {% icon hands_on %} Hands-on: K-mer based evaluation with Merqury
>
> 1. {% tool [Merqury](toolshed.g2.bx.psu.edu/repos/iuc/merqury/merqury/1.3+galaxy2) %} with the following parameters:
>    - *"Evaluation mode"*: `Default mode`
>        - {% icon param-file %} *"K-mer counts database"*: `Meryldb from HiFi`
>        - *"Number of assemblies"*: `One assembly`
>            - {% icon param-files %} *"Genome assembly"*: `Scerevisiae-INSC1019.flye.30x.fa`, `Scerevisiae-INSC1019.hifiasm.30x.fa` and `GCA_000146045.2_genomic.fna`
>
{: .hands_on}

# Evaluation against reference genome (Chromeister)

> ### {% icon hands_on %} Hands-on: Synteny evaluation with Chromeister
>

https://toolshed.g2.bx.psu.edu/repos/iuc/chromeister

> 1. {% tool [Chromeister](toolshed.g2.bx.psu.edu/repos/iuc/chromeister/chromeister/1.5.a) %} with the following parameters:
>    - {% icon param-files %} *"Query sequence"*: `Scerevisiae-INSC1019.flye.30x.fa` and `Scerevisiae-INSC1019.hifiasm.30x.fa`
>    - {% icon param-file %} *"Reference sequence"*: `GCA_000146045.2_genomic.fna`
>
{: .hands_on}


# Conclusion
{:.no_toc}

This pipeline shows how to evaluate a genome assembly. Once you are satisfied with your genome sequence, you might want to purge it, make scaffolding and directly starting the annotation process!
