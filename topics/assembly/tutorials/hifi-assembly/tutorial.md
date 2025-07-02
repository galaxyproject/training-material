---
layout: tutorial_hands_on

title: 'Genome assembly using PacBio data'
zenodo_link: ''
tags:
  - assembly
  - pacbio
  - biodiversity
questions:
- How to perform a genome assembly with PacBio data ?
- How to check assembly quality ?
objectives:
- Assemble a Genome with PacBio data
- Assess assembly quality
time_estimation: 3h
level: Intermediate
key_points:
- PacBio data allows to perform good quality genome assembly
- Quast and BUSCO make it easy to compare the quality of assemblies
contributions:
  authorship:
  - abretaud
  - alexcorm
  - r1corre
  - lleroi
  - stephanierobin
  - scorreard
  funding:
  - gallantries
redirect_from:
  - topics/assembly/tutorials/flye-assembly/tutorial
follow_up_training:
 - type: internal
   topic_name: genome-annotation
   tutorials:
     - repeatmasker
edam_ontology:
- topic_0622 # Genomics
- topic_0196 # Assembly
- topic_3050 # Biodiversity

---



In this tutorial, we will assemble a genome of a species of fungi in the family Aspergillus, *Aspergillus niger*, from PacBio sequencing data. These data will be downloaded from ENA. The quality of the assembly obtained will be analyzed, in particular by comparing it to a reference assembly, available on ENA.


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get data

We will use long reads sequencing data: HiFi (High Fidelity long reads) from PacBio sequencing of *Aspergillus niger* genome. This data is available on ENA. We will also use later a reference genome assembly downloaded from ENA.

## Get data from ENA

> <hands-on-title>Data upload from ENA</hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [ENA](https://www.ebi.ac.uk/ena/browser/)
>
>    ```
>    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR317/012/SRR31719412/SRR31719412_subreads.fastq.gz
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>
> 3. Rename the datasets
> 4. Check that the datatype is `fastq.gz`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
>
{: .hands_on}

## Get reference genome from ENA

> <hands-on-title>Data upload from ENA</hands-on-title>
>
> 1. Reference genome is available here: [ASM4765177v1 assembly for Aspergillus niger](https://www.ebi.ac.uk/ena/browser/view/GCA_047651775.1)
> 2. Download the `WGS Set FASTA (JBKZXA01.fasta.gz)` on your computer
> Make sure you download `WGS Set FASTA (JBKZXA01.fasta.gz)` and *NOT*  `ALL Set FASTA`
> 4. Upload this file on Galaxy (Upload --> Choose local file --> Select the file --> Start --> Close)
> 5. Check that the datatype is `fasta.gz`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

# Genome Assembly

{% include _includes/cyoa-choices.html option1="Hifiasm" option2="Flye" default="Hifiasm"
       text="[Hifiasm](https://github.com/chhylp123/hifiasm) and [Flye](https://github.com/mikolmogorov/Flye) are two well known assembler." %}
<div class="Hifiasm" markdown="1">

Hifiasm is a fast haplotype-resolved de novo assembler initially designed for PacBio HiFi reads. In general, hifiasm generates the assembly graphs in the GFA format, so a step of conversion to fasta is necessary. The GFA 1.0 format is a tab-delimited text format for describing a set of sequences and their overlap.
<br>
Hifiasm produces arguably the best single-sample telomere-to-telomere assemblies combing HiFi, ultralong and Hi-C reads, and it is one of the best haplotype-resolved assemblers for the trio-binning assembly given parental short reads. For a human genome, hifiasm can produce the telomere-to-telomere assembly in one day.

> <hands-on-title>Assembly with Hifiasm</hands-on-title>
>
> 1. {% tool [Hifiasm](toolshed.g2.bx.psu.edu/repos/bgruening/hifiasm/hifiasm/0.25.0+galaxy0) %} with the following parameters:
>    - *"Mode"*: `Standard`
>    - {% icon param-file %} *"Input reads"*: the raw data (fastq.gz)
>    - *"Output log file"*: Set to yes
>
>     The tool produces five datasets: Haplotype-resolved raw unitig graph, Haplotype-resolved processed unitig graph without small bubbles, Primary assembly contig graph, Alternate assembly contig graph, [hap1]/[hap2] contig graph.
>
> 2. {% tool [GFA to FASTA](toolshed.g2.bx.psu.edu/repos/iuc/gfa_to_fa/gfa_to_fa/0.1.2) %} with the following parameters:
>    - {% icon param-file %} *"Input GFA file"*: primary assembly contig graph
{: .hands_on}

> <question-title></question-title>
>
> What are the different output datasets from Hifiasm?
>
> > <solution-title></solution-title>
> >
> > - Haplotype-resolved raw unitig graph: This graph keeps all haplotype information, including somatic mutations and recurrent sequencing errors.
> > - Haplotype-resolved processed unitig graph without small bubbles: This graph 'pops' small bubbles in the raw unitig graph; small bubbles might be caused by somatic mutations or noise in data, which are not the real haplotype information.
> > - Primary assembly contig graph: This graph includes a complete assembly with long stretches of phased blocks, though there may be some haplotype collapse.
> > - Alternate assembly contig graph: This graph consists of all contigs that are discarded from the primary contig graph.
> > - [hap1]/[hap2] contig graph: Each graph consists of phased contigs (output only with Hi-C phasing enabled).
> >
> {: .solution}
>
{: .question}

</div>

<div class="Flye" markdown="1">

We will use *Flye*, a de novo assembler for single molecule sequencing reads, such as those produced by PacBio and Oxford Nanopore Technologies. It is designed for a wide range of datasets, from small bacterial projects to large mammalian-scale assemblies. The package represents a complete pipeline: it takes raw PacBio / ONT reads as input and outputs polished contigs. Flye also has a special mode for metagenome assembly. All informations about Flye assembler are here: [Flye](https://github.com/fenderglass/Flye/).

> <hands-on-title>Assembly with Flye</hands-on-title>
>
> 1. {% tool [Flye](toolshed.g2.bx.psu.edu/repos/bgruening/flye/flye/2.9.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input reads"*: the raw data (fastq.gz)
>    - *"Mode"*: `PacBio HiFi`
>    - *"Number of polishing iterations"*: `1`
>    - *"Reduced contig assembly coverage"*: `Disable reduced coverage for initial disjointing assembly`
>    - *"Generate a log file"*: Set to yes
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
  
</div>

# Quality assessment

## Genome assemblies comparison with **Quast**

A way to calculate metrics assembly is to use ***QUAST = QUality ASsessment Tool***. Quast is a tool to evaluate genome assemblies by computing various metrics. The manual of Quast is here: [Quast](http://quast.sourceforge.net/docs/manual.html#sec3)
Quast allows to compare the assembly carried out with the reference genome and produces statistics such as the genome fraction.
However, when comparing to a reference genome, Quast results do not display the assembly N50, therefore, in this tutorial, you will generate two quast reports : One specifying a reference genome and one not specifying a reference genome.

> <hands-on-title>Quast</hands-on-title>
>
> 1. {% tool [Quast](toolshed.g2.bx.psu.edu/repos/iuc/quast/quast/5.0.2+galaxy3) %} without specifying a reference genome with the following parameters:
>    - *"Assembly mode?"*: `Individual assembly (1 contig file per sample)`
>    - *"Use customized names for the input files?"*: `No, use dataset names`
>        - {% icon param-collection %} *"Contigs/scaffolds file"*: `JBKZXA01.fasta.gz` (reference assembly), `fasta file` (output of **GFA to FASTA** {% icon tool %}) and/or `consensus` (output of **Flye** {% icon tool %})
>    - *"Type of assembly"*: `Genome`
>        - *"Use a reference genome?"*: `No`
>        - *"Type of organism"*: `Fungus: use of GeneMark-ES for gene finding, ...`
>
> 2. {% tool [Quast](toolshed.g2.bx.psu.edu/repos/iuc/quast/quast/5.0.2+galaxy3) %} specifying a reference genome with the following parameters:
>    - *"Assembly mode?"*: `Individual assembly (1 contig file per sample)`
>    - *"Use customized names for the input files?"*: `No, use dataset names`
>        - {% icon param-collection %} *"Contigs/scaffolds file"*: `fasta file` (output of **GFA to FASTA** {% icon tool %}) and/or `consensus` (output of **Flye** {% icon tool %})
>    - *"Type of assembly"*: `Genome`
>        - *"Use a reference genome?"*: `Yes`
>          - *"Reference genome"*: `JBKZXA01.fasta.gz` (reference assembly)
>        - *"Type of organism"*: `Fungus: use of GeneMark-ES for gene finding, ...`
>
> 
{: .hands_on}

> <question-title></question-title>
>
> Compare the different metrics obtained for the assembly you generated and the reference genome.
>
> > <solution-title></solution-title>
> >
> > We compare the metrics of the three assemblies:
> > - The reference genome (JBKZXA01.fasta.gz): 8 contigs, Contig N50 = 3.5Mb, Contig length max = 6.2 Mb, size = 35.4Mb, 49.51% GC
> > - The Hifiasm assembly: 105 contigs, N50 = 4.9 Mb, length max = 7.1Mb, assembly size = 42.1 Mb, 47.91% GC
> > - The Flye assembly: 13 contigs, N50 = 4.9Mb, length max = 6.8Mb, size = 38.7 Mb, 49.31% GC
> >
> {: .solution}
>
{: .question}

## Genome assembly assessment with **BUSCO**

***BUSCO (Benchmarking Universal Single-Copy Orthologs)*** allows a measure for quantitative assessment of genome assembly based on evolutionarily informed expectations of gene content. Details for this tool are here: [Busco website](https://busco.ezlab.org/)

> <hands-on-title>BUSCO on assembly</hands-on-title>
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.8.0+galaxy0) %} with the following parameters:
>    - *"Tool version"*: `Galaxy Version 5.8.0+galaxy0`
>    - {% icon param-file %} *"Sequences to analyse"*: Multiple datasets
>    - {% icon param-collection %} *"Sequences to analyse"*: `JBKZXA01.fasta.gz` (reference assembly), `fasta file` (output of **GFA to FASTA** {% icon tool %}) and/or `consensus` (output of **Flye** {% icon tool %})
>    - *"Lineage data source"*: `Use Cached Lineage Data`
>    - *"Auto-detect or select lineage"*: `Select lineage`
>        - *"Lineage"*: `Ascomycota`
>        - *"Which outputs should be generated"*: `short summary text; summary image`
{: .hands_on}

> <question-title></question-title>
>
> Compare the number of BUSCO genes identified in the generated assembly and the reference genome. What do you observe ?
>
> > <solution-title></solution-title>
> >
> > Short summary generated by BUSCO indicates that reference genome (JBKZXA01.fasta.gz) contains:
> > 1. 1678 Complete BUSCOs (of which 1675 are single-copy and 3 are duplicated),
> > 2. 7 fragmented BUSCOs,
> > 3. 21 missing BUSCOs.
> >
> > Short summary generated by BUSCO indicates that Hifiasm assembly contains:
> > 1. 1685 complete BUSCOs (1679 single-copy and 6 duplicated),
> > 2. 6 fragmented BUSCOs
> > 3. 15 missing BUSCOs.
> >
> > Short summary generated by BUSCO indicates that Flye assembly contains:
> > 1. 1685 complete BUSCOs (1679 single-copy and 6 duplicated),
> > 2. 6 fragmented BUSCOs
> > 3. 15 missing BUSCOs.
> >
> > BUSCO analysis confirms that these two assemblies are of similar quality, with similar number of complete, fragmented and missing BUSCOs genes.
> >
> {: .solution}
>
{: .question}



# Conclusion


This pipeline shows how to generate and evaluate a genome assembly from long reads PacBio data. Once you are satisfied with your genome sequence, you might want to annotate it: have a look at the [RepeatMasker]({% link topics/genome-annotation/tutorials/repeatmasker/tutorial.md %}) and [Funannoate]({% link topics/genome-annotation/tutorials/funannotate/tutorial.md %}) tutorials to learn how to do it!
