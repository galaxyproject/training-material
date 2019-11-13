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
  - "Good example of methods used to perform variant calling in bacteria"
contributors:
  - pvanheus
  - slugger70
---
# Introduction
{:.no_toc}

Tuberculosis (TB) is a infectious disease caused by the bacterium *Mycobacterium tuberculosis*. According to the [WHO](https://www.who.int/tb/publications/global_report/en/), in 2018 there were 10.0 million new cases of TB worldwide and 1.4 million deaths due to the disease, making TB the world's most deadly infectious disease. The [publication](https://www.ncbi.nlm.nih.gov/pubmed/9634230) of the genome of *M. tuberculosis H37Rv* in 1998 gave researchers a powerful new tool in understanding this pathogen. This genome has been revised since then, with the latest version being available
as RefSeq entry [NC_000962.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3/). The genome comprises a single circular chromosome of some 4.4 megabases. The H37Rv strain that the genome was sequenced from is a long-preserved laboratory strain, originally [isolated](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2132400) from a patient in 1905 and [named](https://journals.sagepub.com/doi/abs/10.3181/00379727-33-8330P) as H37Rv in 1935. It is notably different in some genomic [regions](https://www.sciencedirect.com/science/article/pii/S0888754317300617?via%3Dihub) from some modern clinical strains but remains the standard reference sequence for *M. tuberculosis* (Mtb). In a larger context *M. tuberculosis* is a prominent member of the Mycobacterium Tuberculosis Complex (MTBC).

This group of related species comprises of the 7 [lineages](https://www.ncbi.nlm.nih.gov/pubmed/29456241) of human-infecting *M. tuberculosis* as well as predominantly animal-infecting species such as *M. bovis* and *M. pinnipedii*. Two other close relatives of Mtb, *M. leprae* and *M. lepromatosis* circulate between humans, causing the disease leprosy. Finally amongst the Mycobacteria there are several other species that live in the environment and can cause human disease. These are the [Nontuberculous Mycobacteria](https://www.ncbi.nlm.nih.gov/pubmed/28345639).

Variation in the genome of *M. tuberculosis* (Mtb) is associated with changes in phenotype, for example [drug resistance](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0660-8) and virulence. It is also useful for [outbreak investigation](https://www.frontiersin.org/articles/10.3389/fpubh.2019.00087/full) as the single nucleotide polymorphisms (SNPs) in a sample can be used to build a phylogeny.

This tutorial will focus on identifying genomic variation in Mtb and using that do explore drug resistance and other aspects of the bacteria.

TO DO:

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

# Quality control

This step serves the purpose of identifying possible issues with the raw
sequenced reads input data before embarking on any "real" analysis steps.

Some of the typical problems with NGS data can be mitigated by preprocessing
affected sequencing reads before trying to map them to the reference genome.
Detecting some other, more severe problems early on may at least save you a lot
of time spent on analyzing low-quality data that is not worth the effort.

Here, we will perform a standard quality check on our input data and only point
out a few interesting aspects about that data. For a more thorough explanation
of NGS data quality control, you may want to have a look at the dedicated
tutorial on [Quality control]({{ site.baseurl }}{% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %}).

> ### {% icon hands_on %} Hands-on: Quality control of the input datasets
> 1. Run **FastQC** {% icon tool %} on each of your six fastq datasets
>       - {% icon param-files %} *"Short read data from your current history"*: both FASTQ datasets selected with **Multiple datasets**
>
>    {% include snippets/select_multiple_datasets.md %}
>
>    When you start this job, four new datasets (one with the calculated raw
>    data, another one with an html report of the findings for each input
>    dataset) will get added to your history.
>
> 2. Use **MultiQC** {% icon tool %} to aggregate the raw **FastQC** data of all input datasets into one report
>      - In *"Results"*
>        - *"Which tool was used generate logs?"*: `FastQC`
>        - In *"FastQC output"*
>           - *"Type of FastQC output?"*: `Raw data`
>           - {% icon param-files %} *"FastQC output"*: both *RawData*
>             outputs of **FastQC** {% icon tool %})
>
> 3. Inspect the *Webpage* output produced by the tool
>
>    > ### {% icon question %} Questions
>    >
>    > 1. Based on the report, do you think preprocessing of the reads
>    >    (trimming and/or filtering) will be necessary before mapping?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > 1. Sequence quality is quite good overall. If anything you might
>    > >    consider trimming the 3' ends of reads (base qualities decline
>    > >    slightly towards the 3' ends) or to filter out the small fraction
>    > >    of reads with a mean base quality < 5.
>    > >    We will run **Trimmomatic** {% icon tool %} on the
>    > >    fastq datasets in the next step.
>    > >
>    > {: .solution}
>    {: .question}
>
> 4. Use **Trimmomatic** {% icon tool %} to clean up the reads and remove the poor quality sections.
>       - *"Single-end or paired-end reads?"*: `Paired End (two separate inut files)`
>       - {% icon param-files %} *"Input FASTQ file (R1/first of pair)"*: `ERR550641_1.fastq.gz`
>       - {% icon param-files %} *"Input FASTQ file (R2/second of pair)"*: `ERR550641_2.fastq.gz`
>       - *"+Insert Trimmomatic Operation"*
>           - *"Select Trimmomatic operation to perform"*: `Drop reads below a specified length (MINLEN)`
>           - *"Minimum length of reads to be kept"*: `20`
>
> 5. Inspect the output produced by Trimmomatic
>
>    > ### {% icon question %} Questions
>    >
>    > 1. Why are there 4 output read files instead of 2?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > 1. There are 4 output files: Forwards paired and single reads and reverse paired and single reads. The single reads come about when one read in a pair of reads has failed the quality checks and so is deleted. The other half of the pair may still be good and so it is put into the single reads file for the appropriate direction.
>    > >
>    > {: .solution}
>    {: .question}
{: .hands_on}

# Look for contamination with Kraken2

We should also look for contamination in our reads. Sometimes, other sources of DNA accidentally or inadvertantly get mixed in with our sample. Any reads from non-sample sources will confound our snp analysis. **Kraken 2** is an effective way of looking and which species is represented in our reads and so we can easily spot possible contamination of our sample.

> ### {% icon hands_on %} Hands-on: Run Kraken2
>
> 1. **Kraken2** {% icon tool %} with the following parameters
>   - *"Single or paired reads"*: `Paired`
>       - *"Forward Strand"*: `Trimmomatic on https://zenodo.org/record/3531703/files/ERR550641_1.fastq.gz (R1 paired)`
>       - *"Reverse Strand"*: `Trimmomatic on https://zenodo.org/record/3531703/files/ERR550641_2.fastq.gz (R2 paired)`
>   - *"Print scientific names instead of just taxids"*: `Yes`
>   - *"Enable quick operation"*: `Yes`
>   - Under *"Create reports"*:
>       - *"Print a report with aggregrate counts/clade to file"*: `Yes`
>   - *"Select a Kraken2 database"*: `Standard`
>
> 2. Inspect the report produced by Kraken
>
>    > ### {% icon question %} Questions
>    >
>    > 1. Was there any significant contamination of the sample?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > 1. 95.5% of the reads here have been positively identified as *Mycobacterium*. The others found were bacteria from the same kingdom. There were no contaminating human or viral sequences detected.
>    > >
>    > {: .solution}
>    {: .question}
{: .hands_on}

# Find variants with Snippy

We will now run the Snippy tool on our reads, comparing it to the reference.

Snippy is a tool for rapid bacterial SNP calling and core genome alignments. Snippy finds SNPs between a haploid reference genome and your NGS sequence reads. It will find both substitutions (snps) and insertions/deletions (indels).

If we give Snippy an annotated reference, it will silently run a tool called SnpEff which will figure out the effect of any changes on the genes and other features. If we just give Snippy the reference sequence alone without the annotations, it will not run SnpEff.

We have an annotated reference and so will use it in this case.

> ### {% icon hands_on %} Hands-on: Run Snippy
>
> 1. **Snippy** {% icon tool %} with the following parameters
>   - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>   - *"Use the following dataset as the reference sequence"*: `Mycobacterium_tuberculosis_ancestral_reference.gbk`
>   - *"Single or Paired-end reads"*: `Paired`
>       - *"Select first set of reads"*: `Trimmomatic on https://zenodo.org/record/3531703/files/ERR550641_1.fastq.gz (R1 paired)`
>       - *"Select second set of reads"*: `Trimmomatic on https://zenodo.org/record/3531703/files/ERR550641_2.fastq.gz (R2 paired)`
>   - Under *"Advanced parameters"*
>       - *"Minimum proportion for variant evidence"*: `0.1` (This is so we can see possible rare variants in our sample)
>   - Under *"Output selection"* select the following:
>       - *"The final annotated variants in VCF format"*
>       - *"The alignments in BAM format"*
>       - Deselect any others.
>
> 2. Inspect the Snippy VCF output
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What type of variant is the first one in the list?
>    >
>    > 2. What was the effect of this variant on the coding region it was found in?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > 1.
>    > >
>    > > 2.
>    > >
>    > {: .solution}
>    {: .question}
{: .hands_on}

RECAP: So far we have taken our sample reads, cleaned them up a bit, compared them with our reference sequence and then called variants (SNPs) between our sample and the reference genome.



#
