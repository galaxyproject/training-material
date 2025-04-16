---
layout: tutorial_hands_on

title: Genome annotation with Braker3
zenodo_link: https://zenodo.org/records/14770765
tags:
  - eukaryota
  - Braker3
  - jbrowse1

edam_ontology:
- topic_0196 # Sequence Assembly
- topic_0622 # Genomics
- topic_0362 # Genome annotation


questions:
    - How to annotate an eukaryotic genome with Braker3?
    - How to evaluate and visualize annotated genomic features?

objectives:
    - Load genome into Galaxy
    - Annotate genome with Braker3
    - Evaluate annotation quality with BUSCO
    - View annotations in JBrowse

time_estimation: 8h
level: Intermediate
key_points:
    - Braker3 allows to perform structural annotation of an eukaryotic genome
    - BUSCO and OMArk allow to inspect the quality of an annotation
contributions:
  authorship:
    - rlibouba
    - abretaud
  reviewing:
    - deeptivarshney
    - shiltemann

  funding:
    - eurosciencegateway

requirements:
 - type: internal
   topic_name: genome-annotation
   tutorials:
     - repeatmasker

subtopic: eukaryote
---

Annotating the eukaryotic genome represents a somewhat more complex challenge than that of prokaryotes, mainly due to the generally larger size of eukaryotic genomes and their greater number of genes. This annotation can be carried out at different levels of precision, ranging from simple identification of coding and non-coding parts to detailed structural labeling, including for example the precise location of exons, introns and other regulatory elements.

In this tutorial we will use a software tool called Braker3 to annotate the genome sequence of a small eukaryote: [*Mucor mucedo*](https://en.wikipedia.org/wiki/Mucor_mucedo) (a fungal plant pathogen).

[Braker3](https://github.com/Gaius-Augustus/BRAKER) is an automated bioinformatics tool that uses RNA-seq and protein data to annotate genomes. It integrates GeneMark-ETP and AUGUSTUS software to predict genes with a high degree of precision. By combining the results of these two tools, Braker3 generates a final file containing gene annotations with strong extrinsic support (i.e. based on external experimental data).

Braker3 facilitates genome annotation by leveraging transcriptomic and protein data to produce more reliable and robust gene predictions.

In this tutorial, you will learn how to perform structural annotation of the genome and assess its quality. To simplify learning and reduce analysis time, we will use a lightweight dataset, which may result in a less complete annotation than in the context of an analysis under real conditions.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data upload

To annotate our genome using Braker3, we will use the following files:

- The **genome sequence** in fasta format. For best results, the sequence should be soft-masked beforehand. You can learn how to do it by following the [RepeatMasker tutorial]({% link topics/genome-annotation/tutorials/repeatmasker/tutorial.md %}). For this tutorial, we will try to annotate the genome assembled in the [Flye assembly tutorial]({% link topics/assembly/tutorials/flye-assembly/tutorial.md %}). The size of the genome has been reduced to reduce Braker3 execution time.
- Some **RNAseq data** in bam format. These reads are already aligned on the genome, and Braker3 will use it as evidences to annotate genes.
- A set of **protein sequences**, like UniProt/SwissProt. It is important to have good quality, curated sequences here, and the UniProt/SwissProt databank fits very well. For this tutorial, we have prepared a subset of this databank to speed up computing, but you should use UniProt/SwissProt for real life analysis.

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    https://zenodo.org/records/14770765/files/genome.fasta
>    https://zenodo.org/records/14770765/files/RNASeq.bam
>    https://zenodo.org/records/14770765/files/protein_sequences.fasta
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
{: .hands_on}

# Preparing RNAseq data

In this tutorial, the alignments from RNA-seq are already prepared.
This data was obtained using the RNA STAR tool, as explained in detail in the [Funannotate]({% link topics/genome-annotation/tutorials/funannotate/tutorial.md %}) tutorial.

Note that when using Braker3 with STAR to align RNA-Seq data, the **--outSAMstrandField intronMotif** parameter must be added.
This parameter adds specific intron information to the alignment files (BAM).
This information is necessary for Braker3 to correctly understand and use the alignments to annotate genes.
Without this parameter, Braker3 may not function correctly or may produce incomplete results.

Here's how to run STAR with this specific parameter (**"Read alignement tags to include in the BAM output"** option), before annotating with Braker3:

> <hands-on-title></hands-on-title>
>
> 1. {% tool [RNA STAR](toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.8a+galaxy0) %} with the following parameters:
>    - *"Single-end or paired-end reads"*: `Paired-end (as individual datasets)`
>        - {% icon param-file %} *"RNA-Seq FASTQ/FASTA file, forward reads"*: `rnaseq_R1.fq.gz` (Input dataset)
>        - {% icon param-file %} *"RNA-Seq FASTQ/FASTA file, reverse reads"*: `rnaseq_R2.fq.gz` (Input dataset)
>    - *"Custom or built-in reference genome"*: `Use reference genome from history and create temporary index`
>        - {% icon param-file %} *"Select a reference genome"*: `genome_masked.fasta` (Input dataset)
>        - *"Length of the SA pre-indexing string"*: `11`
>    - In *"BAM output format specification"*:
>        - In *"Read alignement tags to include in the BAM output"*:
>             - Select `XS (strand flag, see parameter help below)`
>
{: .hands_on}

# Structural annotation

We now have the following required inputs:

- A masked genome in FASTA format
- RNAseq sequence alignments in BAM format
- Protein annotation in FASTA format

With these files, We can run [**Braker3**](https://github.com/Gaius-Augustus/BRAKER) to perform the structural annotation of the genome.


> <hands-on-title>Genome annotation with Braker3</hands-on-title>
>
> 1. {% tool [Braker3](toolshed.g2.bx.psu.edu/repos/iuc/braker3/braker3/3.0.8+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Assembly to annotate"*: `genome_masked.fasta` (Input dataset)
>    - *"Species name"*: `Mucor mucedo`
>    - {% icon param-file %} *"RNA-seq mapped to genome to train Augustus/GeneMark"*: `rnaseq.bam` (Input dataset)
>    - {% icon param-file %} *"Proteins to map to genome"*: `SwissProt_subset.fasta` (Input dataset)
>    - {% icon param-file %} *"Fungal genome"*: select `Yes`
>    - *"Output format"*: `GFF3`
>
{: .hands_on}

> <comment-title>On parameters</comment-title>
>
> - If you are working with fungi, we recommend that you activate the **--fungus** parameter. This enables Braker3 to run the GeneMark-ETP algorithm with a branch point model.
{: .comment}

> <comment-title>Don't wait</comment-title>
>
> - This step will take a bit of time to run: it could take around 2 hours for this tutorial. While it runs, we can already schedule the following functional annotation steps. Galaxy will run them automatically as soon as the structural annotation is ready.
> - There are several reasons for this: the size of the genome, the quality and quantity of the RNA-Seq data, the supervised learning model for prediction, and more.
{: .comment}

Braker3 generates an output file in GTF format by default. But it is possible to generate the output file in GFF3 output. This file contains informations on gene locations, exons, introns, etc.
GFF3 files are generally preferred for complex and varied annotations, while GTF is often used for transcript annotations in RNA-seq studies.

The GFF3 format is a standard bioinformatics format for storing genome annotations. Each row describes a genomic entity, with columns detailing its identifier, location, score and other attributes.

## Evaluation with **Busco**

[BUSCO](http://busco.ezlab.org/) (Benchmarking Universal Single-Copy Orthologs) is a widely used tool to evaluate the quality of a genome assembly and annotation. By comparing genomes from
various related and distantly related species, the authors determined sets of ortholog genes that are present in single copy in (almost) all the species of a clade (Bacteria, Fungi, Plants, Insects, Mammalians, …).
Most of these genes are essential for the organism to live, and are expected to be found in any newly sequenced and annotated genome from the corresponding clade. Using this data, BUSCO evaluates the "completeness" of genome annotation by assessing the proportion of these essential genes (also named BUSCOs) found in a set of (predicted) transcript or protein sequences.

We want to run BUSCO on the protein sequences predicted from gene sequences of the Braker3 annotation.

So first generate these sequences:

> <hands-on-title>Extract protein sequences</hands-on-title>
>
> 1. {% tool [GFFread](toolshed.g2.bx.psu.edu/repos/devteam/gffread/gffread/2.2.1.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input GFF3 or GTF feature file"*: output of {% tool [Braker](toolshed.g2.bx.psu.edu/repos/iuc/braker3/braker3/3.0.8+galaxy1) %}
>    - In *"Reference Genome"* select: `From your history` (Input dataset)
>    - *"Genome Reference Fasta"*: `masked genome` (Input dataset)
>    - In *"Select fasta outputs"* select: `fasta file with spliced exons for each GFF transcript (-y)`
>    - *"full GFF attribute preservation (all attributes are shown)"*: `Yes`
>    - *"decode url encoded characters within attributes"*: `Yes`
>    - *"warn about duplicate transcript IDs and other potential problems with the given GFF/GTF records"*: `Yes`
>
{: .hands_on}

The parameters for running BUSCO on these protein sequences:

> <hands-on-title>BUSCO in proteome mode</hands-on-title>
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `gffread: pep.fa`
>    - *"Mode"*: `annotated gene sets (protein)`
>    - *"Auto-detect or select lineage?"*: `Select lineage`
>        - *"Lineage"*: `Mucorales`
>    - In *"Advanced Options"*:
>        - *"Which outputs should be generated"*: select all outputs
{: .hands_on}

Several output files are generated:

- short summary: statistical summary of the quality of genomic assembly or annotation, including total number of genes evaluated, percentage of complete genes, percentage of partial genes, etc.
- full table: list of universal orthologs found in the assembled or annotated genome, with information on their completeness, location in the genome, quality score, etc.
- missing buscos: list of orthologs not found in the genome, which may indicate gaps in assembly or annotation.
- summary image: graphics and visualizations to visually represent the results of the evaluation, such as bar charts showing the proportion of complete, partial and missing genes.
- GFF: contain information on gene locations, exons, introns, etc.

This gives information about the completeness of the Braker3 annotation. A good idea is to compare this first result with the one you get on the initial genome sequence,and see if the annotation tool found all the genes that BUSCO finds in the raw genome sequence.

The parameters for running BUSCO on the masked genome:

> <hands-on-title>BUSCO in genome mode</hands-on-title>
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `masked genome` (Input dataset)
>    - *"Mode"*: `Genome assemblies (DNA)`
>    - *"Auto-detect or select lineage?"*: `Select lineage`
>        - *"Lineage"*: `Mucorales`
>    - In *"Advanced Options"*:
>        - *"Which outputs should be generated"*: select all outputs
{: .hands_on}


> <comment-title>What can we deduce from these results?</comment-title>
>
> - Completeness: The results of 50.3% complete BUSCOs for the genome and 43.9% for proteins are relatively low. These results can be explained by the reduction in file size in order to speed up analyses with Braker3.
> - Duplication and fragmentation: The duplication rate is low and the proportion of fragmented genes is small, suggesting that the annotation is correct overall, but incomplete.
{: .comment}

## Evaluation with **OMArk**

[OMArk](https://github.com/DessimozLab/OMArk) is proteome quality assessment software.
It provides measures of proteome completeness, characterises the consistency of all
protein-coding genes with their homologues and identifies the presence of contamination
by other species. OMArk is based on the OMA orthology database, from which it exploits
orthology relationships, and on the OMAmer software for rapid placement of all proteins
in gene families.

OMArk's analysis is based on HOGs (Hierarchical Orthologous Groups), which play a central
role in its assessment of the completeness and coherence of gene sets. HOGs make it possible
to compare the genes of a given species with groups of orthologous genes conserved across a
taxonomic clade.

> <hands-on-title>OMArk on extracted protein sequences</hands-on-title>
>
> 1. {% tool [OMArk](toolshed.g2.bx.psu.edu/repos/iuc/omark/omark/0.3.0+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Protein sequences"*: `gffread: pep.fa`
>    - *"OMAmer database*: select `LUCA-v2.0.0`
>    - In *"Which outputs should be generated"*: select `Detailed summary`
>
{: .hands_on}

The OMArk tool generated an output file in .txt format containing detailed information on
the assessment of the completeness, consistency and species composition of the proteome
analysed. This report includes statistics on conserved genes, the proportion of duplications,
missing genes and the identification of reference lineages.

> <comment-title>What can we deduce from these results?</comment-title>
>
> - Number of conserved HOGs: OMArk has identified a set of 5622 HOGs which are thought to be conserved in the majority of species in the Mucorineae clade.
> - However, 54.43% of these HOGs were absent from the analysed proteome, which may indicate a potential incompleteness of the annotation. This observation could be linked to the reduction in the number of RNA sequences aligned on the genome when Braker was run, a choice made to optimise analysis time.
> - 80.27% of genes are complete, so the annotation is of good quality in terms of genomic completeness.
> - Number of proteins in the whole proteome: 5642. Of which 80.27% are present and 15.37% of the proteome does not share sufficient similarities with known gene families.
> - No contamination detected.
> - The OMArk analysis is based on the Mucorineae lineage, a more recent and specific clade than that used in the BUSCO assessment, which selected the Mucorales as the reference group.
{: .comment}

# Visualisation with a genome browser

You can visualize the annotation generated using a genomic browser like [JBrowse](https://jbrowse.org/jbrowse1.html).
This browser enables you to navigate along the chromosomes of the genome and view the structure of each predicted gene.

> <hands-on-title>JBrowse visualisation</hands-on-title>
>
> 1. {% tool [JBrowse](toolshed.g2.bx.psu.edu/repos/iuc/jbrowse/jbrowse/1.16.11+galaxy1) %} with the following parameters:
>    - *"Reference genome to display"*: `Use a genome from history`
>        - {% icon param-file %} *"Select the reference genome"*: `genome_masked.fasta` (Input dataset)
>    - In *"Track Group"*:
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - *"Track Category"*: `Annotation`
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `GFF/GFF3/BED Features`
>                        - {% icon param-file %} *"GFF/GFF3/BED Track Data"*: `gff3` (output of **Braker3** {% icon tool %})
{: .hands_on}

Click on the newly created dataset's eye to display it. You will see a JBrowse genome browser. You can have a look at the [JBrowse tutorial]({% link topics/visualisation/tutorials/jbrowse/tutorial.md %}) for a more in-depth description of JBrowse.

We have embedded a copy of the resulting JBrowse here, if something went wrong during one of the steps you can always just check this output:


{% capture jbtracks %}DNA%2Caa6aca08d6fd82b94c7a86712a632760_0{% endcapture %}
{% snippet topics/visualisation/faqs/visualizations_jbrowse.html datadir="data" tracks=jbtracks %}

Check that **Braker3 annotation** is selected in the top left-hand corner. If not, select it to view the annotation.

# Conclusion

Congratulations on reaching the end of this tutorial! You now know how to perform a structural annotation of a new eukaryotic genome, using Braker3. And you've learned how to evaluate its quality and how to visualize it using JBrowse.

If you'd like to complete this annotation, we recommend you to follow the tutorial on [functional annotation]({% link topics/genome-annotation/tutorials/functional/tutorial.md %}) with EggNOG Mapper and InterProScan. You can follow it with the protein sequences we generated earlier with gffread.
