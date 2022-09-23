---
layout: tutorial_hands_on

title: Long non-coding RNAs (lncRNAs) annotation with FEELnc
zenodo_link: https://zenodo.org/record/7107050
tags:
  - eukaryote
questions:
  - How to annotate lncRNAs with FEELnc ?
  - How to classify lncRNAs according to their localisation and direction of transcription of proximal RNA transcripts ?
  - How to update genome annotation with these annotated lncRNAs ?

objectives:
  - Load data (genome assembly, annotation and mapped RNASeq) into Galaxy
  - Perform a transcriptome assembly with StringTie
  - Annotate lncRNAs with FEELnc
  - Classify lncRNAs according to their location
  - Update genome annotation with lncRNAs
  #- Evaluate annotation quality with BUSCO
  #- View annotations in JBrowse

time_estimation: 2h
level: Intermediate
key_points:
  - StringTie allows to perform a transcriptome assembly using mapped RNASeq data and provides an annotation file containing trancripts description.
  - FEELnc pipeline allows to perform annotation of long non-coding RNAs (lncRNAs).
  - Annotation is based on reconstructed transcripts from RNA-seq data (either with or without a reference genome)
  - Annotation can be performed without any training set of non-coding RNAs.
  - FEELnc provides the localisation and the direction of transcription of proximal RNA transcripts of lncRNAs.
contributions:
  authorship:
    - stephanierobin
  editing:
    - abretaud
  #funding:

abbreviations:
    lncRNAs: long non-coding RNAs
    mRNAs: messenger RNAs

requirements:
 - type: internal
   topic_name: genome-annotation
   tutorials:
     - funannotate

subtopic: eukaryote
priority: 2
---


# Introduction
{:.no_toc}

Messenger RNAs (mRNAs) are not the only type of RNAs present in different organisms (like mammals, insects or plants) and represent only a small fraction of the transcripts. A vaste repertoire of small (miRNAs, snRNAs) and long non-coding RNAs (lncRNAs) are also present. LncRNAs are generally defined as transcripts longer than 200 nucleotides that are not translated into functional proteins. LncRNAs are important because of their major roles in cellular machinery and their presence in large numbers. Indeed, they are notably involved in gene expression regulation, control of translation or imprinting. Statistics from the [GENCODE project](https://www.gencodegenes.org/human/stats_41.html) reveals that the human genome contains more than 19,095 lncRNA genes, almost as much as the 19,370 protein-coding genes.

Using RNASeq data, we can reconstruct assembled transcripts (with ou without any reference genome) which must then be annotated and identified in particular as mRNAs or lncRNAs.

In this tutorial, we will use a software tool called StringTie({% cite Pertea %}) to assemble the transcripts and a software tool called FEELnc ({% cite Wucher %}) to annotate the assembled transcripts of a small eukaryote: [*Mucor mucedo*](https://en.wikipedia.org/wiki/Mucor_mucedo) (a fungal plant pathogen).

[*StringTie*](https://github.com/gpertea/stringtie) is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts


[*FEELnc (FlExible Extraction of Long non-coding RNA)*](https://github.com/tderrien/FEELnc) is a pipeline to annotate lncRNAs from RNASeq assembled transcripts. FEELnc is composed of 3 modules:
* FEELnc_filter.pl	: Extract, filter candidate transcripts.
* FEELnc_codpot.pl	: Compute the coding potential of candidate transcripts.
* FEELnc_classifier.pl: Classify lncRNAs based on their genomic localization wrt others transcripts.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data upload

To assemble transcriptome with StringTie and annotate lncRNAs with FEELnc, we will use the following files :

- The **genome sequence** in fasta format. For this tutorial, we will use the genome assembled in the [Flye assembly tutorial]({% link topics/assembly/tutorials/flye-assembly/tutorial.md %}).
- The **genome annotation** in gff3 format. For this tutorial, we will use the genome annotation obtained in the [Funannotate tutorial]({% link topics/genome-annotation/tutorials/funannotate/tutorial.md %}).
- Some aligned **RNASeq data** in bam format. For this tutorial, we will use the mapped RNASeq data where mapping was done using STAR.

> ### {% icon hands_on %} Hands-on: Data upload
>
### Get data
{:.no_toc}

> ### {% icon hands_on %} Hands-on: Data upload
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
>    https://zenodo.org/api/files/0f8d27c5-8c8d-4379-90c4-c3cd950de391/genome_assembly.fasta
>    https://zenodo.org/api/files/0f8d27c5-8c8d-4379-90c4-c3cd950de391/genome_annotation.gff3
>    https://zenodo.org/api/files/0f8d27c5-8c8d-4379-90c4-c3cd950de391/all_RNA_mapped.bam
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
{: .hands_on}

# Transcripts assembly with StringTie

StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts.
StringTie takes as input a SAM, BAM or CRAM file sorted by coordinate (genomic location). This file should contain spliced RNA-seq read alignments such as the ones produced by TopHat, HISAT2 or STAR. The TopHat output is already sorted, but the SAM ouput from other aligners should be sorted using the samtools program.

A reference annotation file in GTF or GFF3 format can be provided to StringTie which can be used as 'guides' for the assembly process and help improve the transcript structure recovery for those transcripts.

> ### {% icon hands_on %} Transcripts assembly
>
> {% tool [StringTie](toolshed.g2.bx.psu.edu/repos/iuc/stringtie/stringtie/2.1.7+galaxy1) %} with the following parameters:
>    - *"Input options"*: `Short reads`
>    - {% icon param-file %} *"Input short mapped reads"*: `all_RNA_mapped.bam`
>    - *"Specify strand information"*: Unstranded
>    - *"Use a reference file to guide assembly?"*: Use reference GTF/GFF3
>    - *"Reference file"*: Use a file from history
>        - {% icon param-file %} *"GTF/GFF3 dataset to guide assembly"*: `genome_annotation.gff3`
>    - *"Use Reference transcripts only?"*: `No`
>    - *"Output files for differential expression?"*: `No additional output`
>    - *"Output coverage file"*: `No`
>
{: .hands_on}

We obtained an annotation file (gtf format) which contained all assembled transcripts present in the RNASeq data.

 After this step, the transcriptome is assembled and ready for lncRNAs annotation.

> ### {% icon question %} Question
>
> How many transcripts are assembled  ?
>
> > ### {% icon solution %} Solution
> >
> > Specific features can be extracted from the gtf file using for example `Extract features from GFF data (Galaxy Version 1.0.0)`. By selecting `transcript` From `column 3 / Feature`, we can select only transcripts described in this annotation file.
Assembly contains 14,877 transcripts (corresponding to the lines number of the filtered gtf file).
> >
> {: .solution}
>
{: .question}

# lncRNAs annotation with FEELnc

FEELnc is a pipeline which is composed of 3 steps. These 3 steps run directly when running FEELnc in Galaxy. The first step (FEELnc_filter) consists in filtering out unwanted/spurious transcripts and/or transcripts overlapping (in sense) exons of the reference annotation and especially protein_coding exons as they more probably correspond to new mRNA isoforms.


To use FEELnc, we need to have a reference annotation file in gtf format, which contains protein-coding genes annotation. Presently, we downloaded only the reference annotation file in gff3 format (annotation.gff3). To generate gtf format, we will use **gffread** which converts GFF3/GT2 records.


> ### {% icon hands_on %} Hands-on

> 1. {% tool [gffread](toolshed.g2.bx.psu.edu/repos/devteam/gffread/gffread/2.2.1.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input BED, GFF3 or GTF feature file"*: `genome_annotation.gff3`
>    - *"Feature File Output"*: `GTF`
> Rename output in
>
> 2. {% tool [FEELnc](toolshed.g2.bx.psu.edu/repos/iuc/feelnc/feelnc/0.2) %} with the following parameters:
>    - {% icon param-file %} *"Transcripts assembly"*: `Assembled transcript` (output of **StringTie** {% icon tool %})
>    - {% icon param-file %} *"Reference annotation"*: `genome_annotation.gtf` (Output of **gffread** {% icon tool %})
>    - {% icon param-file %} *"Genome sequence"*: `genome_assembly.fasta`
>
{: .hands_on}

FEELnc provides 3 output files
* lncRNA annotation file	: Annotation file in gtf format which contains the final set of lncRNAs
* mRNA annotation file	: Annotation file in gtf format which contains the final set of mRNAs
* Classifier output file : Table containing classification of lncRNAs based on their genomic localisation w.r.t others transcripts (direction : sense or antisense, type : genic -when the lncRNA gene overlaps an RNA gene from the reference annotation file- or intergenic (lincRNA) -otherwise-).

FEELnc provides also summary file in stdout.

> ### {% icon question %} Question
>
> How many RNAs does this annotation contain ?
> How many interactions between lncRNAs and mRNAs have been identified ?
> Can you describe the different types of lncRNAs ?
>
> > ### {% icon solution %} Solution
> >
> > The summary file indicates 104 lncRNAs and 0 new mRNAs were annotated by FEELnc. The initial annotation contains 13,795 mRNAs annotated. Therefore, a total of 13,898 RNAs are currently annotated.
> > The summary file indicates 652 interactions between lncRNAs and mRNAs. These interactions are described in the Classifier output file.
> > The different types of lncRNAs (intergenic (sense and antisense), intragenic (sense)) are described in the Classifier output file. We observe that the majority of the lncrnas are intergenic. These lncRNAs can each have interactions with several mRNAs. Only 7 lncRNAs are genic. These lncRNAs have only one interaction with the mRNA that contains it.
> >
> {: .solution}
>
{: .question}

# Conclusion
{:.no_toc}

Congratulations for reaching the end of this tutorial! Now you know how to perform an annotation of lncRNAs by using RNASeq data.
