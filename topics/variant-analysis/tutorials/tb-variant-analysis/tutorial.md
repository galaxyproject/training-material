---
layout: tutorial_hands_on

title: "M. tuberculosis Variant Analysis"
zenodo_link: https://doi.org/10.5281/zenodo.3496437
tags:
  - prokaryote
questions:
  - "How do we detect differences between a set of reads from *M. tuberculosis* and a TB reference genome"
objectives:
  - "How should we filter those variants"
  - "How can we predict drug resistance from those variants"
  - "How do we annotate those variants"
time_estimation: "45m"
key_points:
  - "Something"
contributors:
  - pvanheus
---
# Introduction
{:.no_toc}

Tuberculosis (TB) is a infectious disease caused by the bacterium *Mycobacterium tuberculosis*. According to the [WHO](https://www.who.int/tb/publications/global_report/en/), in 2018 there were 10.0 million new cases of TB worldwide and 1.4 million deaths due to the disease, making TB the world's most deadly infectious disease. The [publication](https://www.ncbi.nlm.nih.gov/pubmed/9634230) of the genome of *M. tuberculosis H37Rv* in 1998 gave researchers a powerful new tool in understanding this pathogen. This genome has been revised since then, with the latest version being available
as RefSeq entry [NC_000962.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3/). The genome comprises a single circular chromosome of some 4.4 megabases. The H37Rv strain that the genome was sequenced from is a long-preserved laboratory strain, originally [isolated](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2132400) from a patient in 1905 and [named](https://journals.sagepub.com/doi/abs/10.3181/00379727-33-8330P) as H37Rv in 1935. It is notably different in some genomic [regions](https://www.sciencedirect.com/science/article/pii/S0888754317300617?via%3Dihub) from some modern clinical strains but remains the standard reference sequence for *M. tuberculosis* (Mtb). In a larger context *M. tuberculosis* is a prominent member of the Mycobacterium Tuberculosis Complex (MTBC).

This group of related species comprises of the 7 [lineages](https://www.ncbi.nlm.nih.gov/pubmed/29456241) of human-infecting *M. tuberculosis* as well as predominantly animal-infecting species such as *M. bovis* and *M. pinnipedii*. Two other close relatives of Mtb, *M. leprae* and *M. lepromatosis* circulate between humans, causing the disease leprosy. Finally amongst the Mycobacteria there are several other species that live in the environment and can cause human disease. These are the [Nontuberculous Mycobacteria](https://www.ncbi.nlm.nih.gov/pubmed/28345639).

Variation in the genome of *M. tuberculosis* (Mtb) is associated with changes in phenotype, for example [drug resistance](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0660-8) and virulence. It is also useful for [outbreak investigation](https://www.frontiersin.org/articles/10.3389/fpubh.2019.00087/full) as the single nucleotide polymorphisms (SNPs) in a sample can be used to build a phylogeny.

This tutorial will focus on identifying genomic variation in Mtb and using that do explore drug resistance and other aspects of the bacteria.

* Get data
* QC with FastQC
* MultiQC
* Trimmomatic
* QC and MultiQC
* Kraken2 report
* snippy
* TB variant filter
* TB profiler
* TB variant report

# Get your data

The data for today is a sample of *M. tuberculosis* [collected](https://www.ebi.ac.uk/ena/data/view/PRJEB6945) from a German [outbreak](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1001387) of tuberculosis. In addition to the bacterial sequence sample we will work with a Genbank format version of the genome of the inferred most recent common [ancestor](https://zenodo.org/record/3497110) of the M. tuberculosis complex which is combined with the annotation of the H37Rv reference sequence. This ancestral genome only differs from the H37Rv version 3 genome ([NC_000962.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3)) by the insertion of SNPs to try and model the ancestor of all lineages of Mtb.

This data is available at Zenodo using the following [link](http://doi.org/10.5281/zenodo.3531703).

> ### {% icon hands_on %} Hands-on: Get the data
>
> 1. Import the following three files into a new history
>   - [ERR550641_1.fastq.gz](https://zenodo.org/record/3531703/files/ERR550641_1.fastq.gz)
>   - [ERR550641_2.fastq.gz](https://zenodo.org/record/3531703/files/ERR550641_2.fastq.gz)
>   - [Mycobacterium_tuberculosis_ancestral_reference.gbk](https://zenodo.org/record/3531703/files/Mycobacterium_tuberculosis_ancestral_reference.gbk)
>
>    {% include snippets/import_via_link.md %}
>
{: .hands_on}

If you are using usegalaxy.eu for this tutorial, you can start by importing this [history](https://usegalaxy.eu/u/pvanheus/h/m-tuberculosis-variant-analysis-tutorial). Use the '+' button in the top right hand corner and select a meaningful name for the imported history before clicking 'Import'.




