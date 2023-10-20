---
layout: tutorial_hands_on

title: 'Genome assembly using PacBio data'
zenodo_link: 'https://zenodo.org/record/5702408#.YZUb5uvjIiU'
tags:
  - assembly
  - pacbio
questions:
- How to perform a genome assembly with PacBio data ?
- How to check assembly quality ?
objectives:
- Assemble a Genome with PacBio data
- Assess assembly quality
time_estimation: 6h
level: Intermediate
key_points:
- PacBio data allows to perform good quality genome assembly
- Quast and BUSCO make it easy to compare the quality of assemblies
contributors:
- abretaud
- alexcorm
- r1corre
- lleroi
- stephanierobin
- erasmusplus

follow_up_training:
 - type: internal
   topic_name: genome-annotation
   tutorials:
     - repeatmasker

---


# Introduction


In this tutorial, we will assemble a genome of a species of fungi in the family Mucoraceae, *Mucor mucedo*, from PacBio sequencing data. These data were obtained from NCBI ([SRR8534473](https://www.ncbi.nlm.nih.gov/sra/?term=SRR8534473), [SRR8534474](https://www.ncbi.nlm.nih.gov/sra/?term=SRR8534474) and [SRR8534475](https://www.ncbi.nlm.nih.gov/sra/?term=SRR8534475)). The quality of the assembly obtained will be analyzed, in particular by comparing it to a reference assembly, obtained with Falcon assembler, and available on the [JGI website](https://mycocosm.jgi.doe.gov/Mucmuc1/Mucmuc1.info.html).


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get data

We will use long reads sequencing data: CLR (continuous long reads) from PacBio sequencing of *Mucor mucedo* genome. This data is a subset of data from NCBI. We will also use later a reference genome assembly downloaded from the [JGI website](https://mycocosm.jgi.doe.gov/Mucmuc1/Mucmuc1.info.html). This reference genome was assembled using the same PacBio data, we will use it as a comparison with our own assembly.

## Get data from Zenodo

> <hands-on-title>Data upload from Zenodo</hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo](https://zenodo.org/record/5702408)
>
>    ```
>    https://zenodo.org/api/files/d010d8f1-a1fd-4366-991f-916c2f0c55db/SRR8534473_subreads.fastq.gz
>    https://zenodo.org/api/files/d010d8f1-a1fd-4366-991f-916c2f0c55db/SRR8534474_subreads.fastq.gz
>    https://zenodo.org/api/files/d010d8f1-a1fd-4366-991f-916c2f0c55db/SRR8534475_subreads.fastq.gz
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>
> 3. Rename the datasets
> 4. Check that the datatype is `fastqsanger.gz`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
>
{: .hands_on}

## Get data from JGI website

> <hands-on-title>Data upload from JGI website</hands-on-title>
>
> 1. Create a JGI account in registration page of JGI: [JGI registration](https://contacts.jgi.doe.gov/registration/new)
> 2. Sign in JGI Genome Portal [JGI Genome Portal](https://genome.jgi.doe.gov/portal/)
> 3. Genome assembly is available here: [JGI Mucor mucedo](https://genome.jgi.doe.gov/portal/Mucmuc1/Mucmuc1.download.html)
> 4. Import fasta assembly file `Mucmuc1_AssemblyScaffolds.fasta` on your computer locally
> 5. Upload this file on Galaxy
> 6. Check that the datatype is `fasta`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

# Genome Assembly with **Flye**

We will use *Flye*, a de novo assembler for single molecule sequencing reads, such as those produced by PacBio and Oxford Nanopore Technologies. It is designed for a wide range of datasets, from small bacterial projects to large mammalian-scale assemblies. The package represents a complete pipeline: it takes raw PacBio / ONT reads as input and outputs polished contigs. Flye also has a special mode for metagenome assembly. All informations about Flye assembler are here: [Flye](https://github.com/fenderglass/Flye/).

> <hands-on-title>Assembly</hands-on-title>
>
> 1. {% tool [Flye](toolshed.g2.bx.psu.edu/repos/bgruening/flye/flye/2.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input reads"*: the three sequencing datasets
>    - *"Mode"*: `PacBio raw`
>    - *"Number of polishing iterations"*: `1`
>    - *"Reduced contig assembly coverage"*: `Disable reduced coverage for initial disjointing assembly`
>
>     The tool produces four datasets: consensus, assembly graph, graphical fragment assembly and assembly info
{: .hands_on}

> <question-title></question-title>
>
> What are the different output datasets?
>
> > <solution-title></solution-title>
> >
> > - The first dataset (consensus) is a fasta file containing the final assembly (1461 contigs). You may notice that the result (contigs number) you obtained is sligthy different from the one presented here. This is due to the Flye assembly algorithm which doesn't always give the eact same results.
> > - The second and third dataset are assembly graph files. These graphs are used to represent the final assembly of a genome, they are based on reads and their overlap information. Some tools such as [Bandage](http://rrwick.github.io/Bandage/) allow to visualize the assembly graph.
> > - The fourth dataset is a tabular file (assembly_info) containing extra information about contigs/scaffolds.
> >
> {: .solution}
>
{: .question}

# Quality assessment

## Genome assembly metrics with **Fasta Statistics**

***Fasta statistics*** displays the summary statistics for a fasta file. In the case of a genome assembly, we need to calculate different metrics such as assembly size, scaffolds number or N50 value. These metrics will allow us to evaluate the quality of this assembly.

> <hands-on-title>Fasta statistics on Flye assembly</hands-on-title>
>
> 1. {% tool [Fasta Statistics](toolshed.g2.bx.psu.edu/repos/iuc/fasta_stats/fasta-stats/2.0) %} with the following parameters:
>    - {% icon param-file %} *"fasta or multifasta file"*: `consensus` (output of **Flye** {% icon tool %})
>
{: .hands_on}

> <hands-on-title>Fasta statistics on the reference assembly</hands-on-title>
>
> 1. {% tool [Fasta Statistics](toolshed.g2.bx.psu.edu/repos/iuc/fasta_stats/fasta-stats/2.0) %} with the following parameters:
>    - {% icon param-file %} *"fasta or multifasta file"*: `Mucmuc1_AssemblyScaffolds.fasta`
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Compare the different metrics obtained for Flye assembly and reference genome.
> 2. What can you conclude about the quality of this new assembly ?
>
> > <solution-title></solution-title>
> >
> > 1. We compare the metrics of the two genome assembly:
> > - The Flye assembly: 1461 contigs/scaffolds, N50 = 222 kb, length max = 897 kb, size = 48.6 Mb, 36.6% GC
> > - The reference genome: 456 contigs/scaffolds, N50 = 202 kb, length max = 776 kb, size = 46.1 Mb, 36.7% GC
> >
> > 2. Metrics are very similar, Flye generated an assembly with a quality similar to that of the reference genome.
> >
> {: .solution}
>
{: .question}

## Genome assemblies comparison with **Quast**

Another way to calculate metrics assembly is to use ***QUAST = QUality ASsessment Tool***. Quast is a tool to evaluate genome assemblies by computing various metrics and to compare genome assembly with a reference genome. The manual of Quast is here: [Quast](http://quast.sourceforge.net/docs/manual.html#sec3)

> <hands-on-title>Task description</hands-on-title>
>
> 1. {% tool [Quast](toolshed.g2.bx.psu.edu/repos/iuc/quast/quast/5.0.2+galaxy3) %} with the following parameters:
>    - *"Use customized names for the input files?"*: `No, use dataset names`
>        - {% icon param-file %} *"Contigs/scaffolds file"*: `consensus` (output of **Flye** {% icon tool %})
>    - *"Type of assembly"*: `Genome`
>        - *"Use a reference genome?"*: `Yes`
>        - {% icon param-file %} *"Reference genome"*: `Mucmuc1_AssemblyScaffolds.fasta`
>        - *"Type of organism"*: `Fungus: use of GeneMark-ES for gene finding, ...`
>
{: .hands_on}

> <question-title></question-title>
>
> What additional informations are generated by Quast, compared to the **Fasta Statistics** outputs?
>
> > <solution-title></solution-title>
> >
> > Quast allows us to compare Flye assembly to the reference genome:
> > 1. Genome fraction (90.192 %) is the percentage of aligned bases in the reference genome.
> > 2. Duplication ratio (1.094) is the total number of aligned bases in the assembly divided by the total number of aligned bases in the reference genome.
> > 3. Largest alignment (698452) is the length of the largest continuous alignment in the assembly.
> > 4. Total aligned length (45.2 Mb) is the total number of aligned bases in the assembly.
> >
> > Quast also generates some plots:
> > 1. Cumulative length plot shows the growth of contig lengths. On the x-axis, contigs are ordered from the largest to smallest. The y-axis gives the size of the x largest contigs in the assembly.
> > 2. GC content plot shows the distribution of GC content in the contigs.
> >
> {: .solution}
>
{: .question}

## Genome assembly assessment with **BUSCO**

***BUSCO (Benchmarking Universal Single-Copy Orthologs)*** allows a measure for quantitative assessment of genome assembly based on evolutionarily informed expectations of gene content. Details for this tool are here: [Busco website](https://busco.ezlab.org/)

> <hands-on-title>BUSCO on Flye assembly</hands-on-title>
>
> **First on the Flye assembly:**
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.2.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `consensus` (output of **Flye** {% icon tool %})
>    - *"Auto-detect or select lineage"*: `Select lineage`
>        - *"Lineage"*: `Mucorales`
>
> **Then, on the reference assembly:**
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.2.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `Mucmuc1_AssemblyScaffolds.fasta`
>    - *"Auto-detect or select lineage"*: `Select lineage`
>        - *"Lineage"*: `Mucorales`
>
{: .hands_on}

> <question-title></question-title>
>
> Compare the number of BUSCO genes identified in the Flye assembly and the reference genome. What do you observe ?
>
> > <solution-title></solution-title>
> >
> > Short summary generated by BUSCO indicates that reference genome contains:
> > 1. 2327 Complete BUSCOs (of which 2302 are single-copy and 25 are duplicated),
> > 2. 13 fragmented BUSCOs,
> > 3. 109 missing BUSCOs.
> >
> > Short summary generated by BUSCO indicates that Flye assembly contains:
> > 1. 2348 complete BUSCOs (2310 single-copy and 38 duplicated),
> > 2. 8 fragmented BUSCOs
> > 3. 93 missing BUSCOs.
> >
> > BUSCO analysis confirms that these two assemblies are of similar quality, with similar number of complete, fragmented and missing BUSCOs genes.
> >
> {: .solution}
>
{: .question}



# Conclusion


This pipeline shows how to generate and evaluate a genome assembly from long reads PacBio data. Once you are satisfied with your genome sequence, you might want to annotate it: have a look at the [RepeatMasker]({% link topics/genome-annotation/tutorials/repeatmasker/tutorial.md %}) and [Funannoate]({% link topics/genome-annotation/tutorials/funannotate/tutorial.md %}) tutorials to learn how to do it!
