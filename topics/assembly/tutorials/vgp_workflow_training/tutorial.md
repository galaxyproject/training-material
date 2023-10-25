---
layout: tutorial_hands_on
title: VGP assembly pipeline - short version
zenodo_link: 'https://zenodo.org/record/5887339'
level: Intermediate
tags:
 - pacbio
 - eukaryote
 - VGP
questions:
- "what combination of tools can produce the highest quality assembly of vertebrate genomes?"
- "How can we evaluate how good it is?"
objectives:
- "Learn the tools necessary to perform a de novo assembly of a vertebrate genome"
- "Evaluate the quality of the assembly"
time_estimation: '1h'
key_points:
- "The VGP pipeline allows to generate error-free, near gapless reference-quality genome assemblies"
- "The assembly can be divided in four main stages: genome profile analysis, HiFi long read phased assembly with hifiasm, Bionano hybrid scaffolding and Hi-C hybrid scaffolding"
contributors:
- delphine-l
- astrovsky01
- gallardoalba
- pickettbd
- abueg
- nekrut
abbreviations:
  primary assembly: homozygous regions of the genome plus one set of alleles for the heterozygous loci
  alternate assembly: alternate loci not represented in the primary assembly
  QV: assembly consensus quality
  unitig: A uniquely assembleable subset of overlapping fragments. A unitig is an assembly of fragments for which there are no competing internal overlaps. A unitig is either a correctly assembled portion of a contig or a collapsed assembly of several high-fidelity copies of a repeat.
  contigs: contiguous sequences in an assembly
  scaffolds: one or more contigs joined by gap sequence
  Hi-C: all-versus-all chromatin conformation capture
  HiFi: high fidelity reads
  GWS: Galaxy Workflow System
  VGP: Vertebrate Genome Project
  G10K: Genome 10K
---

# Introduction

The {VGP}, a project of the {G10K} Consortium, aims to generate high-quality, near error-free, gap-free, chromosome-level, haplotype-phased, annotated reference genome assemblies for every vertebrate species ({% cite Rhie2021 %}). The VGP has developed a fully automated *de-novo* genome assembly pipeline, which uses a combination of three different technologies: Pacbio {HiFi}, {Hi-C} data, and (optionally) BioNano optical map data. The pipeline consists of nine distinct workflows. This tutorial provides a quick explanation on how these workflows can be used for producing high quality assemblies. It is intended for the "inpatient" types. 


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Getting started on Galaxy

This tutorial assumes you are comfortable getting data into Galaxy, running jobs, managing history, etc. If you are unfamiliar with Galaxy, we recommed you visit the [Galaxy Training Network](https://training.galaxyproject.org). Consider starting with the following trainings:
- [Introduction to Galaxy]({% link topics/introduction/tutorials/introduction/slides.html %})
- [Galaxy 101]({% link topics/introduction/tutorials/galaxy-intro-101/tutorial.md %})
- [Getting Data into Galaxy]({% link topics/galaxy-interface/tutorials/get-data/slides.html %})
- [Using Dataset Collections]({% link topics/galaxy-interface/tutorials/collections/tutorial.md %})
- [Introduction to Galaxy Analyses](https://training.galaxyproject.org/training-material/topics/introduction)
- [Understanding the Galaxy History System]({% link topics/galaxy-interface/tutorials/history/tutorial.md %})
- [Downloading and Deleting Data in Galaxy]({% link topics/galaxy-interface/tutorials/download-delete-data/tutorial.md %})


# VGP assembly workflow structure

The {VGP} assembly pipeline has a modular organization, consisting in five main subworkflows (fig. 1), each one integrated by a series of data manipulation steps. Firstly, it allows the evaluation of intermediate steps, which facilitates the modification of parameters if necessary, without the need to start from the initial stage. Secondly, it allows to adapt the workflow to the available data.

![Figure 1: The nine workflows of Galaxy assembly pipeline](../../images/vgp_assembly/VGP_workflow_modules.svg "Eight analysis trajectories are possible depending on the combination of input data. A decision on whether or not to invoke Workflow 6 is based on the analysis of QC output of workflows 3, 4, or 5. Thicker lines connecting Workflows 7, 8, and 9 represent the fact that these workflows are invoked separately for each phased assembly (once for maternal and once for paternal).")

----

The first stage of the pipeline is the generation of *k*-mer profiles of the raw reads to estimate genome size, heterozygosity, repetitiveness, and error rate necessary for parameterizing downstream workflows. The generation of *k*-mer counts can be done from HiFi data only (Workflow 1) or include data from parental reads for trio-based phasing (Workflow 2; trio is a combination of paternal sequencing data with that from an offspring that is being assembled). The second stage is the phased contig assembly. In addition to using only {HiFi} reads (Workflow 3), the contig building (contiging) step can leverage {Hi-C} (Workflow 4) or parental read data (Workflow 5) to produce fully-phased haplotypes (hap1/hap2 or parental/maternal assigned haplotypes), using hifiasm16. The contiging workflows also produce a number of critical quality control (QC) metrics such as *k*-mer multiplicity profiles7. Inspection of these profiles provides information to decide whether the third stage—purging of false duplication— is required. Purging (Workflow 6), using [`purge_dups`](https://github.com/dfguan/purge_dups) identifies and resolves haplotype-specific assembly segments incorrectly labeled as primary contigs, as well as heterozygous contig overlaps. This increases continuity and the quality of the final assembly. The purging stage is generally unnecessary for trio data for which reliable haplotype resolution is performed using *k*-mer profiles obtained from parental reads. The fourth stage, scaffolding, produces chromosome-level scaffolds using information provided by Bionano (Workflow 7), with [`Bionano Solve`](https://bionano.com/software-downloads/) (optional) and Hi-C (Workflow 8) data and [`YaHS`](https://github.com/c-zhou/yahsscaffolding) algorithms. A final stage of decontamination (Workflow 9) removes exogenous sequences (e.g., viral and bacterial sequences) from the scaffolded assembly. 

> <comment-title>A note on data quality</comment-title>
> For high quality genome assemblies, we recommends at least 30✕ of each HiFi and Hi-C per haplotype (i.e. 60✕ coverage for diploid genomes). The higher the quality and the coverage, the higher the continuity and the quality of the final assembly. 
{: .comment}

# Assembling the yeast genome

The following steps use PacBio {HiFi} and Illumina {Hi-C} data from baker's yeast ([*Saccharomyces cerevisiae*](https://en.wikipedia.org/wiki/Saccharomyces_cerevisiae)). The tutorial represents trajectory **B** from Fig. 1 above.

## Getting the data

For this tutorial, the first step is to get the datasets from Zenodo. Specifically, we will be uploading two datasets:

1. A set of PacBio {HiFi} reads in `fasta` format
2. A set of Illumina {Hi-C} reads in `fastqsanger.gz` format

### Uploading `fasta` datasets from Zenodo

The following two steps demonstrate how to upload three PacBio {HiFi} datasets into you Galaxy history.

#### **Step 1**: Copy the following URLs into clipboard
(you can do this by clicking on {% icon copy %} button in the right upper corner of the box below. It will appear if you mouse over the box.)

```
https://zenodo.org/record/6098306/files/HiFi_synthetic_50x_01.fasta
https://zenodo.org/record/6098306/files/HiFi_synthetic_50x_02.fasta
https://zenodo.org/record/6098306/files/HiFi_synthetic_50x_03.fasta
```

#### **Step 2**: Upload datasets into Galaxy

These datasets are in `fasta` format. Upload by following the steps shown in the figure below.

![Figure 2: Uploading fasta files in Galaxy]({{site.baseurl}}/faqs/galaxy/images/upload_fasta_via_url.png "Here we upload three fasta files. Compressed (.gz or .bz2) datasets are uploaded in exactly the same fashion by selecting an appropriate datatype (fasta.gz or fasta.bz2)")


### Uploading `fastqsanger.gz` datasets from Zenodo

Illumina {Hi-C} data is uploaded in essentially the same way as shown in the following two steps.

> <warning-title>DANGER: Make sure you choose correct format!</warning-title>
> When selecting datatype in "**Type (set all)**" drop-down, make sure you select `fastaqsanger` or `fastqsanger.gz` BUT NOT `fastqcssanger` or anything else!
{: .warning}

#### **Step 1**: Copy the following URLs into clipboard
(you can do this by clicking on {% icon copy %} button in the right upper corner of the box below. It will appear if you mouse over the box.)

```
https://zenodo.org/record/5550653/files/SRR7126301_1.fastq.gz
https://zenodo.org/record/5550653/files/SRR7126301_2.fastq.gz
```

#### **Step 2**: Upload datasets into Galaxy

These datasets are in `fastqsanger.gz` format. Upload by following the steps shown in the figure below.

![Figure 3: Uploading Fasta files in Galaxy]({{site.baseurl}}/faqs/galaxy/images/upload_fastqsanger_via_url.png "Here we upload two fastqsanger.gz files. Uncompressed or bz2 compressed (.bz2) detests are uploaded in exactly the same fashion by selecting an appropriation datatype (fastqsanger or fastasanger.bz2)")

### Other ways to upload data

You can obviously upload your own datasets via URLs as illustrated above or from your local system. In addition, you can upload data from a major repository called [GenomeArk](https://genomeark.org). GenomeArk is integrated directly into Galaxy Upload. To use follow these steps:

{% snippet faqs/galaxy/dataset_upload_from_genomeark.md %}


{% snippet faqs/galaxy/collections_build_list.md %}


Once we have imported the datasets, the next step is to import the VGP workflows from the WorkflowHub.

## Import workflows from WorkflowHub

{% snippet faqs/galaxy/workflows_import_from_workflowhub.md filter="name:vgp" %}

The workflows imported are marked with a red square in the following figure:

![Figure 2: Workflow menu](../../images/vgp_assembly/imported_workflows.png  "Workflow main menu. The workflow menu lists all the workflows that have been imported. It provides useful information for organizing the workflows, such as last update and the tags. The worklows can be run by clicking in the play icon, marked in red in the image.")

Once we have imported the datasets and the workflows, we can start with the genome assembly.

> <comment-title>Workflow-centric Research Objects</comment-title>
>
> In WorkflowHub, workflows are packaged, registered, downloaded and exchanged as Research Objects using the RO-Crate specification, with test and example data, managed metadata profiles, citations and more.
>
{: .comment}

# Genome profile analysis

[{% icon exchange %} Switch to step by step version]({% link topics/assembly/tutorials/vgp_genome_assembly/tutorial.md %}#genome-profile-analysis)

Now that our data and workflows are imported, we can run our first workflow. Before the assembly can be run, we need to collect metrics on the properties of the genome under consideration, such as the expected genome size according to our data. The present pipeline uses **Meryl** for generating the k-mer database and **Genomescope2** for determining genome characteristics based on a k-mer analysis.

> <hands-on-title>VGP genome profile analysis workflow</hands-on-title>
>
> 1. Click in the **Workflow** menu, located in the top bar
> 2. Click in the {% icon workflow-run %} **Run workflow** buttom corresponding to `VGP genome profile analysis`
> 3. In the **Workflow: VGP genome profile analysis** menu:
>   - {% icon param-collection %} "*Collection of Pacbio Data*": `7: HiFi_collection`
>   - "*K-mer length*": `31`
>   - "*Ploidy*": `2`
> 4. Click on the <kbd>Run workflow</kbd> buttom
>
> > <comment-title>K-mer length</comment-title>
> > In this tutorial, we are using a k-mer length of 31. This can vary, but the VGP pipeline tends to use a k-mer length of 21, which tends to work well for most mammalian-size genomes. There is more discussion about k-mer length trade-offs in the extended VGP pipeline tutorial.
> {: .comment}
>
{: .hands_on}

Once the workflow has finished, we can evaluate the linear plot generated by **Genomescope** (fig. 3), which includes valuable information such as the observed k-mer profile, fitted models and estimated parameters. This file corresponds to the dataset `26`.

![Figure 3: Genomescope plot](../../images/vgp_assembly/genomescope_plot.png "GenomeScope2 k-mer profile. The first peak located at about 25x corresponds to the heterozygous peak. The second peak at 50x, corresponds to the homozygous peak. The plot also includes information about the the inferred total genome length (len), genome unique length percent (uniq), overall heterozygosity rate (ab), mean k-mer coverage for heterozygous bases (kcov), read error rate (err), average rate of read duplications (dup) and k-mer size (k).")

This distribution is the result of the Poisson process underlying the generation of sequencing reads. As we can see, the k-mer profile follows a bimodal distribution, indicative of a diploid genome. The distribution is consistent with the theoretical diploid model (model fit > 93%). Low frequency *k*-mers are the result of sequencing errors, and are indicated by the red line. GenomeScope2 estimated a haploid genome size of around 11.7 Mbp, a value reasonably close to the *Saccharomyces* genome size.


# Assembly with hifiasm

[{% icon exchange %} Switch to step by step version]({% link topics/assembly/tutorials/vgp_genome_assembly/tutorial.md %}#assembly-with-hifiasm)

To generate {contigs}, the VGP pipeline uses **hifiasm**. After genome profiling, the next step is to run the **VGP HiFi phased assembly with hifiasm and HiC data workflow**. This workflow uses **hifiasm** (HiC mode) to generate HiC-phased haplotypes (hap1 and hap2). This is in contrast to its default mode, which generates primary and alternate pseudohaplotype assemblies. This workflow includes three tools for evaluating assembly quality: **gfastats**, **BUSCO** and **Merqury**.

> <hands-on-title>VGP HiFi phased assembly with hifiasm and HiC data workflow</hands-on-title>
> 1. Click in the **Workflow** menu, located in the top bar
> 2. Click in the {% icon workflow-run %} **Run workflow** buttom corresponding to `VGP HiFi phased assembly with hifiasm and HiC data`
> 3. In the **Workflow: VGP HiFi phased assembly with hifiasm and HiC data** menu:
>   - {% icon param-collection %} "*Pacbio Reads Collection*": `7. HiFi_collection`
>   - {% icon param-file %} "*Meryl database*": `12: Meryl on data 11, data 10, data 9: read-db.meryldb`
>   - {% icon param-file %} "*HiC forward reads*": `3: Hi-C_dataset_F`
>   - {% icon param-file %} "*HiC reverse reads*": `2: Hi-C_dataset_R`
>   - {% icon param-file %} "*Genomescope summary dataset*": `19: Genomescope on data 13 Summary`
> 4. Click on the <kbd>Run workflow</kbd> button
>
> > <comment-title>Input option order</comment-title>
> > Note that the order of the inputs may differ slightly.
> {: .comment}
>
{: .hands_on}

Let's have a look at the stats generated by **gfastats**. This output summarizes some main assembly statistics, such as contig number, N50, assembly length, etc.

According to the report, both assemblies are quite similar; the primary assembly includes 18 {contigs}, whose cumulative length is around 12.2Mbp. The alternate assembly includes 17 contigs, whose total length is 11.3Mbp. As we can see in figure 4a, both assemblies come close to the estimated genome size, which is as expected since we used hifiasm-HiC mode to generate phased assemblies which lowers the chance of false duplications that can inflate assembly size.

> <comment-title>Are you working with pri/alt assemblies?</comment-title>
> This tutorial uses the hifiasm-HiC workflow, which generates phased hap1 and hap2 assemblies. The phasing helps lower the chance of false duplications, since the phasing information helps the assembler know which genomic variation is heterozygosity at the same locus versus being two different loci entirely. If you are working with primary/alternate assemblies (especially if there is no internal purging in the initial assembly), you can expect higher false duplication rates than we observe here with the yeast HiC hap1/hap2.
{: .comment}

> <question-title></question-title>
>
> 1. What is the longest contig in the primary assembly? And in the alternate one?
> 2. What is the N50 of the primary assembly?
> 3. Which percentage of reads mapped to each assembly?
>
> > <solution-title></solution-title>
> >
> > 1. The longest contig in the primary assembly is 1.532.843 bp, and 1.532.843 bp in the alternate assembly.
> > 2. The N50 of the primary assembly is 922.430 bp.
> > 3. According to the report, 100% of reads mapped to both the primary assembly and the alternate assembly.
> >
> {: .solution}
>
{: .question}

Next, we are going to evaluate the outputs generated by **BUSCO**. This tool provides quantitative assessment of the completeness of a genome assembly in terms of expected gene content. It relies on the analysis of genes that should be present only once in a complete assembly or gene set, while allowing for rare gene duplications or losses ({% cite Simo2015 %}).

![Figure 5 : BUSCO](../../images/vgp_assembly/BUSCO_full_table.png "BUSCO full table. It contains the complete results in a tabular format with scores and lengths of BUSCO matches, and coordinates.")

As we can see in the report, the results are simplified into four categories: *complete and single-copy*, *complete and duplicated*, *fragmented* and *missing*.

> <question-title></question-title>
>
> 1. How many complete BUSCO genes have been identified in the primary assembly?
> 2. How many BUSCOs genes are absent?
>
> > <solution-title></solution-title>
> >
> > 1. According to the report, our assembly contains the complete sequence of 2080 complete BUSCO genes (97.3%).
> > 2. 19 BUSCO genes are missing.
> >
> {: .solution}
>
{: .question}

Despite **BUSCO** being robust for species that have been widely studied, it can be inaccurate when the newly assembled genome belongs to a taxonomic group that is not well represented in [OrthoDB](https://www.orthodb.org/). Merqury provides a complementary approach for assessing genome assembly quality metrics in a reference-free manner via *k*-mer copy number analysis. Specifically, it takes our hap1 as the first genome assembly, hap2 as the second genome assembly, and the merylDB generated previously for k-mer counts. Like the other QC metrics we have been looking at, the VGP Hifiasm-HiC workflow will automatically generate the Merqury analysis.

By default, **Merqury** generates three collections as output: stats, plots and {QV} stats. The "stats" collection contains the completeness statistics, while the "QV stats" collection contains the quality value statistics. Let's have a look at the copy number (CN) spectrum plot, known as the *spectra-cn* plot. The spectra-cn plot looks at both of your assemblies (here, your haplotypes) taken *together* (fig. 6a).  We can see a small amount of false duplications here: at the 50 mark on the x-axis, there is a small amount of k-mers present at 3-copy across the two assemblies (the green bump).

![Figure 6: Merqury spectra-cn plot for initial yeast contigs](../../images/vgp_assembly/yeast_c_merqury_cn.png "Merqury CN plot for these yeast haplotypes. The plot tracks the multiplicity of each k-mer found in the read set and colors it by the number of times it is found in a given assembly. Merqury connects the midpoint of each histogram bin with a line, giving the illusion of a smooth curve. K-mer distribution of both haplotypes (a). K-mer distribution of an individual haplotype (b)"){:width="100%"}

 Thus, we know there is some false duplication (the 3-copy green bump) present as 2-copy in one of our assemblies, but we don't know which one. We can look at the individual copy number spectrum for each haplotype in order to figure out which one contains the 2-copy k-mers (*i.e.*, the false duplications). In the Merqury spectra-CN plot for hap2 we can see the small bump of 2-copy k-mers at around the 50 mark on the x-axis (fig. 6b).

Now that we know which haplotype contains the false duplications, we can run the purging workflow to try to get rid of these duplicates.

# Purging with purge_dups

[{% icon exchange %} Switch to step by step version]({% link topics/assembly/tutorials/vgp_genome_assembly/tutorial.md %}#purging-with-purgedups)

An ideal haploid representation would consist of one allelic copy of all heterozygous regions in the two haplotypes, as well as all hemizygous regions from both haplotypes ({% cite Guan2019 %}). However, in highly heterozygous genomes, assembly algorithms are frequently not able to identify the highly divergent allelic sequences as belonging to the same region, resulting in the assembly of those regions as separate contigs. In order to prevent potential issues in downstream analysis, we are going to run the **VGP purge assembly with purge_dups workflow**, which will allow to identify and reassign heterozygous contigs. This step is only necessary if haplotypic duplications are observed, and the output should be carefully checked for overpurging.

> <hands-on-title>VGP purge assembly with purge_dups pipeline  workflow</hands-on-title>
>
> 1. Click in the **Workflow** menu, located in the top bar
> 2. Click in the {% icon workflow-run %} **Run workflow** buttom corresponding to `VGP purge assembly with purge_dups pipeline`
> 3. In the **Workflow: VGP purge assembly with purge_dups pipeline** menu:
>   - {% icon param-file %} "*Hifiasm Primary assembly*": `39: Hifiasm HiC hap1`
>   - {% icon param-file %} "*Hifiasm Alternate assembly*": `40: Hifiasm HiC hap2`
>   - {% icon param-collection %} "*Pacbio Reads Collection - Trimmed*": `22: Cutadapt`
>   - {% icon param-file %} "*Genomescope model parameters*": `20: Genomescope on data 13 Model parameters`
> 4. Click in the <kbd>Run workflow</kbd> buttom
>
{: .hands_on}

This workflow generates a large number of outputs, among which we should highlight the datasets `74` and `91`, which correspond to the purged primary and alternative assemblies respectively.

# Hybrid scaffolding with Bionano optical maps

[{% icon exchange %} Switch to step by step version]({% link topics/assembly/tutorials/vgp_genome_assembly/tutorial.md %}#hybrid-scaffolding-with-bionano-optical-maps)

Once the assemblies generated by **hifiasm** have been purged, the next step is to run the **VGP hybrid scaffolding with Bionano optical maps workflow**. It will integrate the information provided by optical maps with primary assembly sequences in order to detect and correct chimeric joins and misoriented contigs. In addition, this workflow includes some additonal steps for evaluating the outputs.

> <hands-on-title>VGP hybrid scaffolding with Bionano optical maps workflow</hands-on-title>
>
> 1. Click in the **Workflow** menu, located in the top bar
> 2. Click in the {% icon workflow-run %} **Run workflow** buttom corresponding to `VGP hybrid scaffolding with Bionano optical maps`
> 3. In the **Workflow: VGP hybrid scaffolding with Bionano optical maps** menu:
>   - {% icon param-file %} "*Bionano data*": `1: Bionano_dataset`
>   - {% icon param-file %} "*Hifiasm Purged Assembly*": `90: Purge overlaps on data 88 and data 33: get_seqs purged sequences` (note: the data numbers may differ in your workflow, but it should still be tagged `p1` from the purge_dups workflow)
>   - {% icon param-file %} "*Estimated genome size - Parameter File*": `60: Estimated Genome size`
>   - "*Is genome large (>100Mb)?*": `No`
> 4. Click on the <kbd>Run workflow</kbd> buttom
{: .hands_on}

Once the workfow has finished, let's have a look at the assembly reports.

As we can observe in the cumulative plot of the file `119` (fig. 7a), the total length of the assembly (12.160.926 bp) is slightly larger than the expected genome size. With respect to the NG50 statistic (fig. 7b), the value is 922.430 bp, which is significantly higher than the value obtained during the first evaluation stage (813.039 bp). This increase in NG50 means this scaffolded assembly is more contiguous compared to the non-scaffolded contigs.

It is also recommended to examine **BUSCO** outputs. In the summary image (fig. 7c), which can be found in the daset `117`, we can appreciate that most of the universal single-copy orthologs are present in our assembly at the expected single-copy.

> <question-title></question-title>
>
> 1. How many scaffolds are in the primary assembly after the hybrid scaffolding?
> 2. What is the size of the largest scaffold? Has this changed with respect to the previous evaluation?
> 3. What is the percentage of completeness on the core set genes in BUSCO? Has Bionano scaffolding increased the completeness?
>
> > <solution-title></solution-title>
> >
> > 1. The number of contigs is 17.
> > 2. The largest contig is 1.531.728 bp long. This value hasn't changed. This is expected, as the VGP pipeline implementation of Bionano scaffolding does not allow for breaking contigs.
> > 3. The percentage of complete BUSCOs is 95.7%. Yes, it has increased, since in the previous evaluation the completeness percentage was 88.7%.
> >
> {: .solution}
>
{: .question}


# Hi-C scaffolding

[{% icon exchange %} Switch to step by step version]({% link topics/assembly/tutorials/vgp_genome_assembly/tutorial.md %}#hi-c-scaffolding)

In this final stage, we will run the **VGP hybrid scaffolding with HiC data**, which exploits the fact that the contact frequency between a pair of loci strongly correlates with the one-dimensional distance between them. This information allows further scaffolding the Bionano scaffolds using **SALSA2**, usually generating chromosome-level scaffolds.


> <hands-on-title>VGP hybrid scaffolding with HiC data</hands-on-title>
>
> 1. Click in the **Workflow** menu, located in the top bar
> 2. Click in the {% icon workflow-run %} **Run workflow** buttom corresponding to `VGP hybrid scaffolding with HiC data`
> 3. In the **Workflow: VGP hybrid scaffolding with HiC data** menu:
>   - {% icon param-file %} "*Scaffolded Assembly*": `114: Concatenate datasets on data 110 and data 109`
>   - {% icon param-file %} "*HiC Forward reads*": `3: Hi-C_dataset_F (as fastqsanger)`
>   - {% icon param-file %} "*HiC Reverse reads*": `2: Hi-C_dataset_R (as fastqsanger)`
>   - {% icon param-file %} "*Estimated genome size - Parameter File*": `50: Estimated Genome size`
>   - "*Is genome large (>100Mb)?*": `No`
> 4. Click in the <kbd>Run workflow</kbd> buttom
{: .hands_on}

In order to evaluate the Hi-C hybrid scaffolding, we are going to compare the contact maps before and after running the HiC hybrid scaffolding workflow (fig. 8), corresponding to the datasets `130` and `141` respectively.

![Figure 8: Pretext final contact map](../../images/vgp_assembly/hi-c_pretext_final.png "Hi-C maps generated by Pretext using Hi-C data. The red circles indicate the  differences between the contact maps generated after (a) and before (b) Hi-C hybrid scaffolding.")


The regions marked with red circles highlight the most notable difference between the two contact maps, where inversion has been fixed.


# Conclusion

To sum up, it is worthwhile to compare the final assembly with the [_S. cerevisiae_ S288C reference genome](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_assembly_stats.txt).

![Figure 9: Final stats](../../images/vgp_assembly/stats_conclusion.png "Comparison between the final assembly generated in this training and the reference genome. Contiguity plot using the reference genome size (a). Assemby statistics (b).")

With respect to the total sequence length, we can conclude that the size of our genome assembly is almost identical to the reference genome (fig.9a,b). It is noteworthy that the reference genome consists of 17 sequences, while our assembly includes only 16 chromosomes. This is due to the fact that the reference genome also includes the sequence of the mitochondrial DNA, which consists of 85,779 bp. The remaining statistics exhibit very similar values (fig. 9b).

![Figure 10: Comparison reference genome](../../images/vgp_assembly/hi-c_pretext_conclusion.png "Comparison between contact maps generated using the final assembly (a) and the reference genome (b).")

If we compare the contact map of our assembled genome (fig. 10a) with the reference assembly (fig. 10b), we can see that the two are indistinguishable, suggesting that we have generated a chromosome level genome assembly.

# FAQs



