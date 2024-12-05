---
layout: tutorial_hands_on
title: Using the VGP workflows to assemble a vertebrate genome with HiFi and Hi-C
  data
zenodo_link: https://zenodo.org/record/5887339
level: Intermediate
tags:
- pacbio
- eukaryote
- VGP
questions:
- What combination of tools can produce the highest quality assembly of vertebrate
  genomes?
- How can we evaluate how good it is?
objectives:
- Learn the tools necessary to perform a de novo assembly of a vertebrate genome
- Evaluate the quality of the assembly
time_estimation: 2h
key_points:
- The VGP pipeline allows to generate error-free, near gapless reference-quality genome
  assemblies
- 'The assembly can be divided in four main stages: genome profile analysis, HiFi
  long read phased assembly with hifiasm, Bionano hybrid scaffolding and Hi-C hybrid
  scaffolding'
contributors:
- delphine-l
- astrovsky01
- gallardoalba
- pickettbd
- abueg
- nekrut
abbreviations:
  primary assembly: homozygous regions of the genome plus one set of alleles for the
    heterozygous loci
  alternate assembly: alternate loci not represented in the primary assembly
  QV: assembly consensus quality
  unitig: A uniquely assembleable subset of overlapping fragments. A unitig is an
    assembly of fragments for which there are no competing internal overlaps. A unitig
    is either a correctly assembled portion of a contig or a collapsed assembly of
    several high-fidelity copies of a repeat.
  contigs: contiguous sequences in an assembly
  collection: Galaxy's way to represent multiple datasets as a single interface entity
  collections: Galaxy's way to represent multiple datasets as a single interface entity
  scaffold: one or more contigs joined by gap sequence
  scaffolds: one or more contigs joined by gap sequence
  Hi-C: all-versus-all chromatin conformation capture
  HiFi: high fidelity reads
  GWS: Galaxy Workflow System
  VGP: Vertebrate Genome Project
  G10K: Genome 10K
recordings:
- youtube_id: _LO-migvwcM
  length: 1H56M
  galaxy_version: 24.1.2.dev0
  date: '2024-09-12'
  speakers:
  - abueg
  captioners:
  - abueg
  bot-timestamp: 1726177737
- youtube_id: TODO
  length: 1H30M
  galaxy_version: 24.1.3.dev0
  date: '2024-09-30'
  speakers:
  - mschatz
  - delphine-l
  captioners:
  - mschatz
  - delphine-l
  bot-timestamp: 1727666641


---




The {VGP}, a project of the {G10K} Consortium, aims to generate high-quality, near error-free, gap-free, chromosome-level, haplotype-phased, annotated reference genome assemblies for every vertebrate species ({% cite Rhie2021 %}). The VGP has developed a fully automated *de-novo* genome assembly pipeline, which uses a combination of three different technologies: Pacbio {HiFi}, {Hi-C} data, and (optionally) Bionano optical map data. The pipeline consists of nine distinct workflows. This tutorial provides a quick example of how to run these workflows for one particular scenario, which is, based on our experience, the most common: assembling genomes using {HiFi} Reads combined with {Hi-C} data (both generated from the same individual).

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Getting started on Galaxy

This tutorial assumes you are comfortable getting data into Galaxy, running jobs, managing history, etc. If you are unfamiliar with Galaxy, we recommend you visit the [Galaxy Training Network](https://training.galaxyproject.org). Consider starting with the following trainings:
- [Introduction to Galaxy]({% link topics/introduction/tutorials/introduction/slides.html %})
- [Galaxy 101]({% link topics/introduction/tutorials/galaxy-intro-101/tutorial.md %})
- [Getting Data into Galaxy]({% link topics/galaxy-interface/tutorials/get-data/slides.html %})
- [Using Dataset Collections]({% link topics/galaxy-interface/tutorials/collections/tutorial.md %})
- [Understanding the Galaxy History System]({% link topics/galaxy-interface/tutorials/history/tutorial.md %})
- [Introduction to Galaxy Analyses]({% link topics/introduction/index.md %})
- [Downloading and Deleting Data in Galaxy]({% link topics/galaxy-interface/tutorials/download-delete-data/tutorial.md %})

# The VGP-Galaxy pipeline

The {VGP} assembly pipeline has a modular organization, consisting in ten workflows (Fig. 1). It can used with the following types of input data:

|  Input data | Assembly quality  | Analysis trajectory <br>([Fig. 1)](#figure-1)|
|------|---------------|-----|
| HiFi | The minimum requirement | A |
| HiFi + HiC| Better continuity | B |
| HiFi + Bionano | Better continuity | C |
| HiFi + Hi-C + Bionano | Even better continuity | D |
| HiFi + parental data| Better haplotype resolution | E |
| HiFi + parental data + Hi-C| Better haplotype resolution and improved continuity | F |
| HiFi + parental + Bionano | Better haplotype resolution and improved continuity | G |
| HiFi + parental data + Hi-C + Bionano | Better haplotype resolution and ultimate continuity | H |

In this table, HiFi and Hi-C data are derived from the individual whose genome is being assembled. "Parental data" refers to high coverage whole genome resequencing data from the parents of the individual being assembled. Assemblies utilizing parental data for phasing are also called "*Trios*". Each combination of input datasets is supported by an *analysis trajectory*: a combination of workflows designed for generating assembly given a particular combination of inputs. These trajectories are listed in the table above and shown in the figure below. 

![The nine workflows of Galaxy assembly pipeline](../../images/vgp_assembly/VGP_workflow_modules.svg "Eight analysis trajectories are possible depending on the combination of input data. A decision on whether or not to invoke Workflow 6 is based on the analysis of QC output of workflows 3, 4, or 5. Thicker lines connecting Workflows 7, 8, and 9 represent the fact that these workflows are invoked separately for each phased assembly (once for maternal and once for paternal).")
<br>
The first stage of the pipeline is the generation of a *k*-mer profile of the raw reads to estimate genome size, heterozygosity, repetitiveness, and error rate necessary for parameterizing downstream workflows. The generation of *k*-mer counts can be done from HiFi data only (Workflow 1) or include data from parental reads for trio-based phasing (Workflow 2). The second stage is contig assembly. In addition to using only {HiFi} reads (Workflow 3), the contig-building (contiging) step can leverage {Hi-C} (Workflow 4) or parental read data (Workflow 5) to produce fully-phased haplotypes (hap1/hap2 or parental/maternal assigned haplotypes), using [`hifiasm`](https://github.com/chhylp123/hifiasm). The contiging workflows also produce a number of critical quality control (QC) metrics such as *k*-mer multiplicity profiles. Inspection of these profiles provides information to decide whether the third stage—purging of false duplication—is required. Purging (Workflow 6), using [`purge_dups`](https://github.com/dfguan/purge_dups) identifies and resolves haplotype-specific assembly segments incorrectly labeled as primary contigs, as well as heterozygous contig overlaps. This increases continuity and the quality of the final assembly. The purging stage is generally unnecessary for trio data, as thses assemblies achieve reliable haplotype resolution using *k*-mer set operations based on parental reads. Generally, HiC-phased contigs also do not need to undergo purging, but if so, then there exists Workflow 6b to purge a single haplotype at a time. The fourth stage, scaffolding, produces chromosome-level scaffolds using information provided by Bionano (Workflow 7), with [`Bionano Solve`](https://Bionano.com/software-downloads/) (optional) and Hi-C (Workflow 8) data and [`YaHS`](https://github.com/c-zhou/yahsscaffolding) algorithms. A final stage of decontamination (Workflow 9) removes exogenous sequences (*e.g.*, viral and bacterial sequences) from the scaffolded assembly. A separate workflow (WF0) is used to assemble and annotate a mitogenome.

> <comment-title>A note on data quality</comment-title>
> We suggest at least 30✕ PacBio HiFi coverage and 60✕ Hi-C coverage (both diploid coverage).
{: .comment}

# Getting the data

The following steps use PacBio {HiFi} and Illumina {Hi-C} data from baker's yeast ([*Saccharomyces cerevisiae*](https://en.wikipedia.org/wiki/Saccharomyces_cerevisiae)). The tutorial represents trajectory **B** from Fig. 1 above. For this tutorial, the first step is to get the datasets from Zenodo. Specifically, we will be uploading two datasets:

1. A set of PacBio {HiFi} reads in `fasta` format. Please note that your HiFi reads received from a sequencing center will usually be fastqsanger.gz format, but the dataset used in this tutorial has been converted to fasta for space..
2. A set of Illumina {Hi-C} reads in `fastqsanger.gz` format.

## Uploading `fasta` datasets from Zenodo

The following two steps demonstrate how to upload three PacBio {HiFi} datasets into your Galaxy history.


> <hands-on-title> Uploading <tt>FASTA</tt> datasets from Zenodo </hands-on-title>
>
> 1. Create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Copy the following URLs into clipboard.
>
>    (You can do this by clicking on {% icon copy %} button in the right upper corner of the box below. It will appear if you mouse over the box.)
>
>    ```
>    https://zenodo.org/record/6098306/files/HiFi_synthetic_50x_01.fasta
>    https://zenodo.org/record/6098306/files/HiFi_synthetic_50x_02.fasta
>    https://zenodo.org/record/6098306/files/HiFi_synthetic_50x_03.fasta
>    ```
>
> 3. Upload datasets into Galaxy.
>    - set the datatype to `fasta`
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md format="fasta" %}
>
>    {% snippet topics/assembly/tutorials/vgp_genome_assembly/faqs/dataset_upload_fasta_via_urls.md %}
>
{: .hands_on}


## Uploading `fastqsanger.gz` datasets from Zenodo

Illumina {Hi-C} data is uploaded in essentially the same way as shown in the following two steps.

> <warning-title> Make sure you the choose correct format!</warning-title>
> When selecting datatype in "**Type (set all)**" drop-down, make sure you select `fastqsanger` or `fastqsanger.gz` BUT NOT `fastqcssanger` or anything else!
{: .warning}

> <hands-on-title> Uploading <tt>fastqsanger.gz</tt> datasets from Zenodo </hands-on-title>
>
> 1. Copy the following URLs into clipboard.
>     - you can do this by clicking on {% icon copy %} button in the right upper corner of the box below. It will appear if you mouse over the box.
>
>    ```
>    https://zenodo.org/record/5550653/files/SRR7126301_1.fastq.gz
>    https://zenodo.org/record/5550653/files/SRR7126301_2.fastq.gz
>    ```
>
> 2. Upload datasets into Galaxy.
>    - set the datatype to `fastqsanger.gz`
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md format="fasta" %}
>
>    {% snippet topics/assembly/tutorials/vgp_genome_assembly/faqs/dataset_upload_fastqsanger_via_urls.md %}
>
> 3. Rename the datasets as follow: 
>    - Rename `SRR7126301_1.fastq.gz` as `Hi-C forward reads`
>    - Rename `SRR7126301_2.fastq.gz` as `Hi-C reverse reads`
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
{: .hands_on}



> <warning-title>These datasets are large!</warning-title>
> Hi-C datasets are large. It will take some time (~15 min) for them to be fully uploaded. Please, be patient.
{: .warning}

## Organizing the data

If everything goes smoothly, your history will look like shown in Fig. 4 below. The three {HiFi} fasta files are better represented as a collection: {collection}. Also, importantly, the workflow we will be using for the analysis of our data takes a collection as an input (it does not access individual datasets). So let's create a collection using steps outlines in the Tip {% icon tip %} "Creating a dataset collection" that you can find below Fig. 4.

![AfterUpload](../../images/vgp_assembly/making_list.svg "History after uploading HiFi and HiC data (left). Creation of a list (collection) combines all HiFi datasets into a single history item called 'HiFi data' (right). See below for instruction on how to make this collection.")

{% snippet faqs/galaxy/collections_build_list.md %}

> <details-title>Other ways to upload the data</details-title>
> You can obviously upload your own datasets via URLs as illustrated above or from your own computer. In addition, you can upload data from a major repository called [GenomeArk](https://genomeark.org). GenomeArk is integrated directly into Galaxy Upload. To use GenomeArk following the steps in the Tip {% icon tip %} below:
>
> {% snippet faqs/galaxy/datasets_upload_from_genomeark.md %}
{: .details}


Once we have imported the datasets, the next step is to import the workflows necessary for the analysis of our data from [DockStore](https://dockstore.org).

# Importing workflows

All analyses described in this tutorial are performed using *workflows*--chains of tools--shown in [Fig. 1](#figure-1). Specifically, we will use four workflows corresponding to analysis trajectory **B**: 1, 4, 6, and 8. To use these four workflows you need to import them into your Galaxy account following the steps below. Note: these are not necessarily the latest versions of the actual workflows, but versions that have been tested for this tutorial. To see the latest versions, see the [Galaxy Project VGP workflows page](https://galaxyproject.org/projects/vgp/workflows/) and click on the Dockstore links to import workflows. **Alternatively, for each section of the tutorial, there will be a "Launch [workflow] (View on Dockstore)" link at the beginning, which you can use to launch the workflow.**

# Performing the assembly

Workflows listed in [Fig. 1](#figure-1) support a variety of "analysis trajectories". The majority of species that were sequenced by the {VGP} usually contain {HiFi} reads for the individual being sequenced supplemented with {Hi-C} data. As a result most assemblies performed by us follow the trajectory **B**. This is why this tutorial was designed to follow this trajectory as well.

## Genome profile analysis (WF1)

Before the assembly can be run, we need to collect metrics on the properties of the genome under consideration, such as the expected genome size according to our data. The present pipeline uses **Meryl** for generating the *k*-mer database and **Genomescope2** for determining genome characteristics based on a *k*-mer analysis.

### Launching the workflow

{% snippet faqs/galaxy/workflows_run_ds.md title="Genome profile analysis (WF1)" dockstore_id="github.com/iwc-workflows/kmer-profiling-hifi-VGP1/main" version="v0.1.7" %}

> <hands-on-title> Running <i>K</i>-mer profile analysis workflow </hands-on-title>
>
> 1. **Identify inputs**
>
>    The profiling workflow takes the following inputs:
>   
>    1. {HiFi} reads as a collection
>    2. *K*-mer length
>    3. Ploidy
>
> 2. **Launch *k*-mer profiling workflow**
>
>    1. Click on the **Workflow** menu, located in the top bar
>    2. Click on the {% icon workflow-run %} **Run workflow** buttom corresponding to `kmer-profiling-hifi-VGP1`
>    3. In the **Workflow: kmer-profiling-hifi-VGP1** menu:
>        - {% icon param-collection %} "*Collection of Pacbio Data*": `HiFi data` collection
>        - "*K-mer length*": `31`
>        - "*Ploidy*": `2`
>    4. Click on the {% icon workflow-run %} **Run workflow** button
>   
>    This should look like this:
>
>    ![Parameters of *k*-mer profiling workflow](../../images/vgp_assembly/wf1_launch_ui.png  "Workflow main menu. The workflow menu lists all the workflows that have been imported. It provides useful information for organizing the workflows, such as last update and the tags. The worklows can be run by clicking in the play icon, marked in red in the image.")
>
>    > <comment-title>K-mer length</comment-title>
>    > In this tutorial, we are using a *k*-mer length of 31. This can vary, but the VGP pipeline tends to use a *k*-mer length of 21, which tends to work well for most mammalian-size genomes. There is more discussion about *k*-mer length trade-offs in the extended VGP pipeline tutorial.
>    {: .comment}
>
>
> 3. **Refill your coffee**
>
>    Assembly is not exactly an instantaneous type of analysis - this workflow will take approximately 15 minutes to complete. The same is true for all analyses in tutorial.
>
{: .hands_on}

### Interpreting the results

Once the workflow has finished, we can evaluate the transformed linear plot generated by [**Genomescope**](https://github.com/schatzlab/Genomescope), which includes valuable information such as the observed *k*-mer profile, fitted models, and estimated parameters. This file corresponds to the dataset `18` in this [history](https://usegalaxy.org/u/delphinel/h/vertebrate-genome-assembly-training). 

![Genomescope plot described further in caption](../../images/vgp_assembly/genomescope_plot.png "GenomeScope2 <i>k</i>-mer profile. The first peak located at about 25&times; corresponds to the heterozygous peak. The second peak at 50&times; corresponds to the homozygous peak. The plot also includes information about the inferred total genome length (len), genome unique length percent (uniq), overall heterozygosity rate (ab), mean <i>k</i>-mer coverage for heterozygous bases (kcov), read error rate (err), average rate of read duplications (dup) and <i>k</i>-mer size (k).")

This distribution is the result of the Poisson process underlying the generation of sequencing reads. As we can see, the *k*-mer profile follows a bimodal distribution, indicative of a diploid genome. The distribution is consistent with the theoretical diploid model (model fit > 93%). Low frequency *k*-mers are the result of sequencing errors, and are indicated by the red line. Genomescope2 estimated a haploid genome size of around 11.7 Mbp, a value reasonably close to the *Saccharomyces* genome size. We are using the transformed linear plot because it shows the ploidy better by reducing low ccoverage peaks (like *k*-mers containing errors) and increasing higher coverage peaks ({% cite Ranallo_Benavidez_2020 %}). 

It is worth noting that the genome characteristics such as length, error percentage, etc., are based on the GenomeScope2 model, which is the black line in the plot. If the model (black line) does not fit your observed data (blue bars), then these estimated characteristics might be very off. In the case of this tutorial, the model is a good fit to our data, so we can trust the estimates. 

## Assembly (contiging) with `hifiasm` (WF4)

To generate {contigs} we will use the [**hifiasm**](https://github.com/chhylp123/hifiasm) assembler. It is a part of the "Assembly with HiC (WF4)" workflow . This workflow uses **hifiasm** (HiC mode) to generate HiC-phased haplotypes (hap1 and hap2). This is in contrast to its default mode, which generates primary and alternate pseudohaplotype assemblies. This workflow includes three tools for evaluating assembly quality: [**gfastats**](https://github.com/vgl-hub/gfastats), [**BUSCO**](https://busco.ezlab.org/) and [**Merqury**](https://github.com/marbl/merqury).

### Launching the workflow

{% snippet faqs/galaxy/workflows_run_ds.md title="Assembly HiFi-HiC phasing (WF4)" dockstore_id="github.com/iwc-workflows/Assembly-Hifi-HiC-phasing-VGP4/main" version="v0.2.1" %}

> <hands-on-title> Launching assembly (contiging) workflow </hands-on-title>
>
> 1. **Identify inputs**
>
>    The assembly workflow takes the following inputs:
>
>    1. {HiFi} reads as a collection
>    2. Forward Hi-C reads
>    3. Reverse Hi-C reads
>    4. `Genomescope` Model Parameters generated by previous (*k*-mer profiling) workflow
>    5. `Genomescope` Summary generated by previous (*k*-mer profiling) workflow
>    6. Meryl *k*-mer database generated by previous (*k*-mer profiling) workflow
>    7. Busco Database
>    8. Busco lineage
>
> 2. **Launch the workflow**
>
>    1. Click on the **Workflow** menu, located in the top bar
>    2. Click on the {% icon workflow-run %} **Run workflow** button corresponding to `Assembly-Hifi-HiC-phasing-VGP4`
>    3. In the **Workflow: Assembly-Hifi-HiC-phasing-VGP4** menu fill the following parameters:
>        - {% icon param-collection %} "*Pacbio Reads Collection*": `HiFi data` collection
>        - {% icon param-file %} "*HiC forward reads*": `Hi-C forward reads`
>        - {% icon param-file %} "*HiC reverse reads*": `Hi-C reverse reads`
>        - {% icon param-file %} "*GenomeScope Summary*": `GenomeScope on data X Summary` (contains tag "`GenomeScopeSummary`")
>        - {% icon param-file %} "*Meryl database*": `Meryl on data X: read-db.mertyldb`: the Meryl *k*-mer database (contains tag "`MerylDatabase`")
>        - "*Database for Busco Lineage*": `Busco v5 Lineage Datasets` (or the latest version available in the instance)
>        - "*Provide lineage for BUSCO (e.g., Vertebrata)*": `Ascomycota`
>        - {% icon param-file %} "*GenomeScope Model Parameters*": `GenomeScope on data X Model parameters` (contains tag "`GenomeScopeParameters`")
>    4. Click on the {% icon workflow-run %} **Run workflow** button
> 
{: .hands_on}

> <comment-title>A note about "Homozygous Read Coverage"</comment-title>
>
>  Hifiasm tries to estimate the homozygous read coverage, but sometimes it can mistake a different peak as the homozygous coverage peak. To get around this, the VGP workflows calculate the homozygous coverage from the haploid coverage estimate from the GenomeScope outputs generated in the *k*-mer profiling workflow. However, sometimes GenomeScope can also misidentify the haploid peak, leading to it feeding an incorrect homozygous read coverage value into hifiasm. *If this is the case, you can run the hifiasm workflows but specify the homozygous read coverage yourself. Otherwise, this parameter is fine to leave blank to let the workflows try to get it right on their own.*
>
{: .comment}



### Interpreting the results

> <warning-title>There will be two assemblies!</warning-title>
> Because we are running hifiasm in HiC-phasing mode, it will produce two assemblies: hap1 and hap2!
{: .warning}

Let's have a look at the stats generated by **gfastats**. This output summarizes some main assembly statistics, such as contig number, N50, assembly length, etc. The workflow provide a joined table to display the statistics for both haplotype assemblies side-by-side. Below we provide a partial output of this file called `Assembly statistics for Hap1 and Hap2` in your history:

>| Statistic | Hap 1 | Hap 2 |
>|-----------|----------:|------:|
>| # contigs | 16 | 17 |
>| Total contig length | 11,304,507 | 12,160,985 |
>| Average contig length | 706,531.69  | 715,352.06 |
>| Contig N50 |  922,430 | 923,452 |
>| Contig N50 | 922,430 | 923,452 |
>| Contig auN | 895,018.82  | 904,515.40 |
>| Contig L50 | 5 | 6 |
>| Contig L50 | 5 | 6 |
>| Contig NG50 | 813,311 | 923,452 |
>| Contig NG50 | 813,311 | 923,452 |
>| Contig auNG | 861,278.90 | 936,364.06 |
>| Contig LG50 | 6 | 6 |
>| Contig LG50 | 6 | 6 |
>| Largest contig | 1,532,843 | 1,531,728 |
>| Smallest contig | 185,154 | 85,850 |
{: .matrix}

According to the report, both assemblies are quite similar; the hap1 assembly includes 16 {contigs}, whose cumulative length is around 11.3 Mbp. The hap2 assembly includes 17 contigs, whose total length is 12.1Mbp. Both assemblies come close to the estimated genome size, which is as expected since we used hifiasm-HiC mode to generate phased assemblies, and this lowers the chance of false duplications that can inflate assembly size.

> <comment-title>Are you working with pri/alt assemblies?</comment-title>
> This tutorial uses the hifiasm-HiC workflow, which generates phased hap1 and hap2 assemblies. The phasing helps lower the chance of false duplications, since the phasing information helps the assembler know which genomic variation is heterozygosity at the same locus versus being two different loci entirely. If you are working with primary/alternate assemblies (especially if there is no internal purging in the initial assembly), you can expect higher false duplication rates than we observe here with the yeast HiC hap1/hap2.
{: .comment}

> <question-title></question-title>
>
> 1. What is the longest contig in the hap1 assembly? And in the hap2 one?
> 2. What is the N50 of the hap2 assembly?
>
> > <solution-title></solution-title>
> >
> > 1. The longest contig in the hap1 assembly is 1,532,843 bp, and 1,531,728 bp in the hap2 assembly.
> > 2. The N50 of the hap2 assembly is 923,452 bp.
> >
> {: .solution}
>
{: .question}

Next, we are going to evaluate the outputs generated by **BUSCO**. This tool provides qualitative assessment of the completeness of a genome assembly in terms of expected gene content. It relies on the analysis of genes that should be present only once in a complete assembly or gene set, while allowing for rare gene duplications or losses ({% cite Simo2015 %}).

<br>

![BUSCO assessment](../../images/vgp_assembly/busco_after_contiging.svg "A composite of BUSCO completeness summaries for hap1 (left) and hap2 (right)")

<br>

As we can see in the report, the results are simplified into four categories: *complete and single-copy*, *complete and duplicated*, *fragmented* and *missing*.

> <question-title></question-title>
>
> 1. How many complete BUSCO genes have been identified in the hap1 assembly?
> 2. How many BUSCOs genes are absent in the hap1 assembly?
>
> > <solution-title></solution-title>
> >
> > 1. According to the report, our assembly contains the 1,436 complete BUSCO genes -- this includes ones that are single-copy and ones that are duplicated.
> > 2. 208 BUSCO genes are missing.
> >
> {: .solution}
>
{: .question}

Despite **BUSCO** being robust for species that have been widely studied, it can be inaccurate when the newly assembled genome belongs to a taxonomic group that is not well represented in [OrthoDB](https://www.orthodb.org/). `Merqury` provides a complementary approach for assessing assembly quality in a reference-free manner via *k*-mer copy number analysis. Specifically, it takes our hap1 as the first genome assembly, hap2 as the second genome assembly, and the merylDB generated previously from our sequencing reads for *k*-mer counts.

By default, `Merqury` generates three collections as output: stats (completeness stats), plots and {QV} stats. Let's first have a look at the copy number (CN) spectrum plot, known as the *spectra-cn* plot. The spectra-cn plot looks at both of your assemblies (here, your two haplotypes) taken *together* (fig. 6a). 

![Figure 6: Merqury spectra-cn plot for initial yeast contigs](../../images/vgp_assembly/yeast_c_merqury_cn.svg "Merqury CN plot for yeast assemblies. The plot tracks the multiplicity of each <i>k</i>-mer found in the read set and colors it by the number of times it is found in a given assembly. Merqury connects the midpoint of each histogram bin with a line, giving the illusion of a smooth curve. <b>a)</b>. <i>K</i>-mer distribution of both haplotypes. <b>b)</b>. <i>K</i>-mer distribution of an individual haplotype (hap2)."){:width="100%"}

Our haplotypes look clean! In the spectra-cn plot for both haplotypes, the peaks are where we should expect them. There is a peak of *k*-mers present at 1-copy (so, only in either hap1 or hap2) with *k*-mer multiplicity of ~25, corresponding to heterozygous coverage. These are our heterozygous alleles being represented only once in our diplolid assemblies, which matches with their haploid coverage. We also have a peak of 2-copy *k*-mers present at diploid coverage, around ~50. This looks good so far, but we should also look at *k*-mer multiplicity for each individual haplotype separately, not just combined as they were in the previous plot. Figure 6b shows that most of the *k*-mers in our assembly are 1-copy -- we have one copy of all the homozygous regions (the ones with 50x coverage) and of about half the heterozygous regions (the ones with 25x coverage). At the 25x coverage point, there is a "read-only" peak, which is expected, as these are the alternative alleles for those heterozygous loci, and those *k*-mers are likely in the other assembly, since they were not missing from the overall spectra-cn plot.


## Hi-C scaffolding (WF8)

In this final stage, we will run the **Scaffolding HiC YAHS (WF8)** workflow, which uses the fact that the contact frequency between a pair of loci strongly correlates with the one-dimensional distance between them (*i.e.*, linear distance on a chromosome). This information allows [**YAHS**](https://github.com/c-zhou/yahs) -- the main tool in this workflow -- to generate scaffolds that are often chromosome-sized.

### Launching Hi-C scaffolding workflow

> <warning-title>The scaffolding workflow is run on <b>ONE</b> haplotype at a time.</warning-title>
> Contiging (VGP4) works with both (hap1/hap2, primary/alternate) assemblies simultaneously. This is not the case for contiging -- it has to be run independently for each haplotype assembly. In this example (below) we run contiging on the hap1 assembly only. You can run the same analysis on the second haplotype by replacing hap1 with hap2. 
{: .warning}

{% snippet faqs/galaxy/workflows_run_ds.md title="Scaffolding HiC (WF8)" dockstore_id="github.com/iwc-workflows/Scaffolding-HiC-VGP8/main" version="v0.2.7" %}

> <hands-on-title> Launching Hi-C scaffolding workflow </hands-on-title>
>
> 1. **Identify inputs**
>
>    The scaffolding workflow takes the following inputs:
>
>    1. An assembly graph
>    2. Forward Hi-C reads
>    3. Reverse Hi-C reads
>    4. Estimated genome size parsed from GenomeScope summary by the previous run of assembly workflow (VGP4).
>    5. Restriction enzymes used in Hi-C library preparation procedure
>    6. Busco Database
>    7. Busco lineage
>
> 2. **Launch scaffolding workflow (WF8)**
>
>    1. Click on the **Workflow** menu, located in the top bar
>    2. Click on the {% icon workflow-run %} **Run workflow** button corresponding to `Scaffolding-HiC-VGP8`
>    3. In the **Workflow: Scaffolding-HiC-VGP8** menu:
>        - {% icon param-file %} "*input GFA*": `usable hap1 gfa` Output of the contiging workflow (VGP4) with a tag `hic_hap1_gfa` for hap1 assembly.
>        - "*Database for Busco Lineage*": `Busco v5 Lineage Datasets` (or the latest version available in the instance)
>        - "*Provide lineage for BUSCO (e.g., Vertebrata)*": `Ascomycota`
>        - {% icon param-file %} "*HiC forward reads*": `Hi-C forward reads`
>        - {% icon param-file %} "*HiC reverse reads*": `Hi-C reverse reads`
>        - "*Restriction enzymes*": `Dovetail Omni-C: enzyme-free prep` For this tutorial, we'll use the Omni-C option as it is the equivalent of not specifying any restriction enzyme cutsites, but for your own data you would want to select the appropriate option. 
>        - {% icon param-file %} "*Estimated genome size - Parameter File*": `Estimated genome size`: An output of the contiging workflow (VGP4) with a tag `estimated_genome_size`.
>    4. Click on the {% icon workflow-run %} **Run workflow** button
{: .hands_on}

### Interpreting the results

In order to evaluate the Hi-C hybrid scaffolding, we are going to compare the contact maps before and after running the HiC hybrid scaffolding workflow (Fig. below). They will have the following tags:

- Before scaffolding: `pretext_s1`
- After scaffolding: `pretext_s2`

Below is the comparison of the two maps obtained from the yeast tutorial data versus real example from a zebra finch (*Taeniopygia guttata*) genome assembly:

![Pretext final contact map](../../images/vgp_assembly/hi-c_pretext_final.svg "Hi-C maps generated by Pretext before and after scaffolding with Hi-C data. The red circles indicate the differences between the contact maps generated before and after Hi-C hybrid scaffolding. The bottom two panels show results of scaffolding on zebra finch where scaffolding dramatically decreases the number of segments by merging multiple contigs into scaffolds.")

The regions marked with red circles highlight the most notable difference between the two contact maps, where inversion has been fixed.

# Conclusion

To sum up, it is worthwhile to compare the final assembly with the [_S. cerevisiae_ S288C reference genome](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_assembly_stats.txt):

![Quast plot](../../images/vgp_assembly/quast_plot.png "Cumulative continuity plot comparing assembly generated here (red line) with existing yeast reference (black dotted line). Our assembly is slightly smaller (11,304,507 bp versus 12,071,326. Our assembly is lacking the mitochondrial genome (~86 kb) because the initial data does not include mitochondrial reads. This is partially responsible for this discrepancy. ")

With respect to the total sequence length, we can conclude that the size of our genome assembly is very similar to the reference genome. It is noteworthy that the reference genome consists of 17 sequences, while our assembly includes only 16 chromosomes. This is due to the fact that the reference genome also includes the sequence of the mitochondrial DNA, which consists of 85,779 bp. (The above comparison is performed using {% tool [Quast](toolshed.g2.bx.psu.edu/repos/iuc/quast/quast/5.2.0+galaxy1) %} using Primary assembly generated with scaffolding workflow (WF8) and yeast reference.)

![Comparison reference genome](../../images/vgp_assembly/hi-c_pretext_conclusion.svg "Comparison between contact maps generated using the final Primary assembly from this tutorial (left) and the reference genome (right).")

If we compare the contact map of our assembled genome with the reference assembly (Fig. above), we can see that the two are indistinguishable, suggesting that we have generated a chromosome level genome assembly.


