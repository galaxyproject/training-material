---
layout: tutorial_hands_on

title: Long non-coding RNAs (lncRNAs) annotation with FEELnc
zenodo_link: https://zenodo.org/record/7107050
tags:
  - eukaryote
questions:
  - How to annotate lncRNAs with FEELnc?
  - How to classify lncRNAs according to their localisation and direction of transcription of proximal RNA transcripts?
  - How to update genome annotation with these annotated lncRNAs?

objectives:
  - Load data (genome assembly, annotation and mapped RNASeq) into Galaxy
  - Perform a transcriptome assembly with StringTie
  - Annotate lncRNAs with FEELnc
  - Classify lncRNAs according to their location
  - Update genome annotation with lncRNAs

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

abbreviations:
  lncRNAs: long non-coding RNAs
  LncRNAs: Long non-coding RNAs
  mRNAs: Messenger RNAs

requirements:
 - type: internal
   topic_name: genome-annotation
   tutorials:
     - funannotate

subtopic: eukaryote
priority: 7
---


# Introduction


{mRNAs} are not the only type of RNAs present in organisms (like mammals, insects or plants) and represent only a small fraction of the transcripts. A vast repertoire of small (miRNAs, snRNAs) and {lncRNAs} are also present. {LncRNAs} are generally defined as transcripts longer than 200 nucleotides that are not translated into functional proteins. They are important because of their major roles in cellular machinery and their presence in large number. Indeed, they are notably involved in gene expression regulation, control of translation or imprinting. Statistics from the [GENCODE project](https://www.gencodegenes.org/human/stats_41.html) reveals that the human genome contains more than 19,095 lncRNA genes, almost as much as the 19,370 protein-coding genes.

Using RNASeq data, we can reconstruct assembled transcripts (with ou without any reference genome) which can then be annotated and identified individually as {mRNAs} or {lncRNAs}.

In this tutorial, we will use a software tool called *StringTie* ({% cite Pertea %}) to assemble the transcripts and then *FEELnc* ({% cite Wucher %}) to annotate the assembled transcripts of a small eukaryote: [*Mucor mucedo*](https://en.wikipedia.org/wiki/Mucor_mucedo) (a fungal plant pathogen).

[*StringTie*](https://github.com/gpertea/stringtie) is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts.


[*FEELnc (FlExible Extraction of Long non-coding RNA)*](https://github.com/tderrien/FEELnc) is a pipeline to annotate {lncRNAs} from RNASeq assembled transcripts. It is composed of 3 modules:
* FEELnc_filter: Extract, filter candidate transcripts.
* FEELnc_codpot: Compute the coding potential of candidate transcripts.
* FEELnc_classifier: Classify {lncRNAs} based on their genomic localization wrt others transcripts.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data upload

To assemble transcriptome with StringTie and annotate {lncRNAs} with FEELnc, we will use the following files :

- The **genome sequence** in fasta format. For this tutorial, we will use the genome assembled in the [Flye assembly tutorial]({% link topics/assembly/tutorials/flye-assembly/tutorial.md %}).
- The **genome annotation** in GFF3 format. We will use the genome annotation obtained in the [Funannotate tutorial]({% link topics/genome-annotation/tutorials/funannotate/tutorial.md %}).
- Some aligned **RNASeq data** in bam format. Here, we will use some mapped RNASeq data where mapping was done using STAR.


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

> <hands-on-title>Transcripts assembly</hands-on-title>
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

We obtain an annotation file (GTF format) which contained all assembled transcripts present in the RNASeq data.

 After this step, the transcriptome is assembled and ready for {lncRNAs} annotation.

> <question-title></question-title>
>
> How many transcripts are assembled  ?
>
> > <solution-title></solution-title>
> >
> > Specific features can be extracted from the GTF file using for example {% tool [Extract features from GFF data](Extract_features1) %}. By selecting `transcript` From `column 3 / Feature`, we can select only the transcript elements present in this annotation file. Assembly contains 14,877 transcripts (corresponding to the number of lines in the filtered GTF file).
> >
> {: .solution}
>
{: .question}

# lncRNAs annotation with FEELnc

FEELnc is a pipeline which is composed of 3 steps. These 3 steps are run automatically when running FEELnc within Galaxy. The first step (FEELnc_filter) consists in filtering out unwanted/spurious transcripts and/or transcripts overlapping (in sense) exons of the reference annotation, and especially protein coding exons as they more probably correspond to new mRNA isoforms.

To use FEELnc, we need to have a reference annotation file in GTF format, which contains protein-coding genes annotation. Presently, we downloaded only the reference annotation file in GFF3 format (`annotation.gff3`). To convert from GFF3 to GTF format, we will use [*gffread*](https://github.com/gpertea/gffread).

> <hands-on-title>FEELnc</hands-on-title>
>
> 1. {% tool [gffread](toolshed.g2.bx.psu.edu/repos/devteam/gffread/gffread/2.2.1.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input BED, GFF3 or GTF feature file"*: `genome_annotation.gff3`
>    - *"Feature File Output"*: `GTF`
>
> 2. {% tool [FEELnc](toolshed.g2.bx.psu.edu/repos/iuc/feelnc/feelnc/0.2) %} with the following parameters:
>    - {% icon param-file %} *"Transcripts assembly"*: `Assembled transcript` (output of **StringTie** {% icon tool %})
>    - {% icon param-file %} *"Reference annotation"*: `genome_annotation.gtf` (Output of **gffread** {% icon tool %})
>    - {% icon param-file %} *"Genome sequence"*: `genome_assembly.fasta`
>
{: .hands_on}

FEELnc provides 3 output files
* lncRNA annotation file: annotation file in GTF format which contains the final set of {lncRNAs}
* mRNA annotation file: annotation file in GTF format which contains the final set of {mRNAs}
* Classifier output file: table containing classification of {lncRNAs} based on their genomic localisation w.r.t other transcripts (direction: `sense` or `antisense`, type: `genic`, if the lncRNA gene overlaps an RNA gene from the reference annotation file or `intergenic` (lincRNA) if not).

FEELnc provides also summary file in stdout.

> <question-title></question-title>
>
> How many RNAs does this annotation contain ?
> How many interactions between {lncRNAs} and {mRNAs} have been identified ?
> Can you describe the different types of {lncRNAs} ?
>
> > <solution-title></solution-title>
> >
> > The summary file indicates 104 {lncRNAs} and 0 new {mRNAs} were annotated by FEELnc. The initial annotation contains 13,795 {mRNAs} annotated. Therefore, a total of 13,898 RNAs are currently annotated.
> >
> > The summary file indicates 652 interactions between {lncRNAs} and {mRNAs}. These interactions are described in the Classifier output file.
> >
> > The different types of {lncRNAs} (intergenic (sense and antisense), intragenic (sense)) are described in the Classifier output file. We observe that the majority of the {lncRNAs} are intergenic. These {lncRNAs} can each have interactions with several {mRNAs}. Only 7 {lncRNAs} are genic. These {lncRNAs} have only one interaction with the mRNA that contains it.
> >
> {: .solution}
>
{: .question}

For future analyses, it would be interesting to use an updated annotation containing {mRNAs} and {lncRNAs} annotations. Thus, we will merge the reference annotation with those obtained with FEELnc.

> <hands-on-title>Merge the annotations</hands-on-title>
>
> {% tool [concatenate](https://toolshed.g2.bx.psu.edu/view/bgruening/text_processing/f46f0e4f75c4) %} with the following parameters:
>    - {% icon param-file %} *"Datasets to concatenate"*: `genome_annotation.gtf`
>    - Insert Dataset
>    - {% icon param-file %} *"Dataset"*: `lncRNA annotation with FEELnc`
>
{: .hands_on}


# Conclusion


Congratulations for reaching the end of this tutorial! Now you know how to perform an annotation of {lncRNAs} by using RNASeq data.
