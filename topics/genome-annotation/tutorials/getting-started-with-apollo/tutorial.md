---
layout: tutorial_hands_on
topic_name: genome_annotation
tutorial_name: getting_started_with_apollo
---

# Introduction
{:.no_toc}

<!-- This is a comment. -->

The Center for Phage Technology Galaxy program is an instance of the Galaxy Project, and provides a web interface for bioinformatics tools. Apollo, a genome browser, can be accessed through Galaxy in order to create, edit, and annotate genomes.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> * Background Information
>
>    > * History of Genome Browsers
>    >    > 1. Artemis
>    >    > 2. GBrowse
>    >    > 3. JBrowse
>    >    > 4. Apollo
>    > * Genome File Formats
>    >    > 1. FASTA
>    >    > 2. GFF3
>    >    > 3. GenBank
>
> * Annotation Within Apollo
>
{: .agenda}

# Background Information

Although it is not the first genome browser, Apollo is the first collaborative genomic annotation editor available solely on the Internet.

> ### {% icon details %} Important Definitions 
> * **Static:** Unmodifiable, specifically in the context of a computer resource that you are accessing. The website that you see cannot be modified by you, the user accessing them. This is opposed to “dynamic” where you can interact with the files or service, and your interactions can persist.
> * **Instance:** A specific copy of a web service made available over the internet. Given that the administrators can run 0-N copies of the same web service, we use the term “instance” to refer to a specific copy of a service.
> * **Tracks:** A set of analysis results that can be shown or hidden depending on the annotator’s needs.
> * **Feature:** Conceptually, a region of a genome with some annotations (such as a Name, Product, or Dbxref). Visually, a rectangular box in a track.
> * **Evidence:** Tracks contain evidence; these are results of specific computer methods (which are documented and citable), which we use to make annotations. Annotations should not be made without evidence. Evidence allows us to move the annotation process from an art to a science.
> * **Annotations:** Annotation is the addition of descriptive features to a DNA sequence, such as a protein’s function, or locating tRNAs, and terminators. The annotation process we do is 100% computer based, so keep in mind that until an annotation is experimentally tested in the lab, it is putative or assumed based on an educated hypothesis.
{: .details}

## History of Genome Browsers

This section will cover a bit of history about Genome Browsers. While not useful to the annotation process, it is important to know what the terms mean and how the parts all fit together, so that the developers and annotators can have a common language.

We use a lot of software under the umbrella term of GMOD, the [Generic Model Organism Database.](http://www.gmod.org/wiki/Main_Page “Generic Model Organism Database”)

GMOD is a collection open source software for maintaining Model Organism Databases (MODs). Having a common platform for MODs is important, as historically individual labs spent effort building their own, custom organism databases, and then faced challenges trying to interoperate with other databases. With GMOD and the associated tools, software that talks to one MOD can be re-used when talking to another MOD. We can use the same tools to work with the CPT’s Phage Database, as people use to access data in Yeast genome databases.

### 1. Artemis

Artemis, the first genome browser discussed here, is not actually a GMOD project. Artemis was an older, desktop-based genome browser. You had to install the software on your computer in order to run it. All of the other genome browsers observed here are web-based. Artemis allows for annotation, but those annotations were only stored on your local device. Artemis featured a three-pane view consisting of a high-level overview of the genome, a DNA-level view, and a list of all the features in the genome.

![](../../images/getting-started-with-apollo-screenshots/(2) Artemis.png)

### 2. GBrowse

One of the earlier genome browsers, Grows did *not* support annotation. Think of it like the old Yahoo-maps. Instead of just clicking and dragging the map, you had to click where you wanted to go, wait a few seconds, and the new map would be displayed. It makes the process tedious.

![](../../images/getting-started-with-apollo-screenshots/(3) GBrowse from GMOD Wiki.png)

### 3. JBrowse

JBrowse is used in Galaxy workflows for genome visualization. JBrowse is a more modern re-implementation of GBrowse. JBrowse is much more like Google Maps (or any other current web map service; you click and drag and can quickly browse around the genome, turning evidence tracks on and off where relevant. Many labs have deployed JBrowse instances to help showcase their annotation efforts to the community, and to make their data accessible. FlyBase has produced a demo in JBrowse, displaying *Drosophila melanogaster*. Note that JBrowse is a **_static_** visualization tool. You cannot make any changes to the data, and you cannot make annotations and save them. It is a “Read Only” view of genomes and annotations.

![](../../images/getting-started-with-apollo-screenshots/(4) JBrowse.png)

### 4. Apollo

Apollo takes JBrowse one step further and adds support for community annotation; it provides a “Read+Write” view of genomes. You can create new annotations on new gene features, and these are shared with everyone who has access to the Apollo server. From a computer perspective, Apollo embeds a copy of JBrowse. For the annotation workflow, we will use both Apollo and JBrowse.

![](../../images/getting-started-with-apollo-screenshots/(5) Compare and Contrast Browsers.png)

## Genome File Formats

There are three file formats to be aware of during genome annotation.
> * **FASTA:** stores genomic sequence information in the form of nucleotide sequences.
> * **GFF3:** stands for general feature format; stores genome annotations.
> * **GenBank:** an older format containing annotations and sequence information.

### 1. FASTA

The sequence contained in a FASTA file may be DNA, RNA, or protein sequences; they may contain unspecified bases (N/Y/X) or gaps (-). Within a FASTA file, each sequence begins with a “>,” which is immediately followed by the “FASTA ID.” Some sequences have a “description” after *any whitespace character*, such as the given example.

![](../../images/getting-started-with-apollo-screenshots/(6) FASTA Sequence Example.png)

### 2. GFF3

The eukaryotic gene model captures a lot of information about the biological process behind producing proteins from DNA, such as mRNAs, transcription, and alternative splicing. GFF3 files thus have to encode these complex, hierarchical, parent-child relationships. Characteristics of GFF3 files, such as the tab separated, key-value pairs, allot simplicity and ease of use.

![](../../images/getting-started-with-apollo-screenshots/(7) GFF3.png)

![](../../images/getting-started-with-apollo-screenshots/(8) GFF3 Visually.png "A visual representation of a GFF3 file.")

<!-- HOW TO PHYSICALL CAPTION SECOND PICTURE, NOT JUST HOVER TEXT --> 

At the top level is a “gene” (3rd column) spanning from 1000 to 9000, on the forward strand (7th column), with an ID of gene00001 and a Name of EDEN.

“Below” the gene is an mRNA feature. The mRNA has a Parent attribute (Parent=gene00001) set to the ID of the “parent” gene feature. This makes it a child of the gene feature.

Similarly all four exons and all four CDSs have a Parent of mRNA00001. ID, Name, and Parent are all known as feature attributes - metadata about a feature. Feature attributes also contain other fields; you will see sometimes see Notes, Products, and many others. Only a handful of these attributes have standards defining what information they contain, and the rest are free to be used as you like. From a computational perspective, we prefer the fields with standardized meanings. If they conform to a standard, we can apply automation in our processing. If they are a free-form “notes” field, we need a human to interpret and codify the evidence.

> ### {% icon tip %} Note that…
> All of this is a little bit excessive for phages (where exons are rare, and mRNAs not involved). Nevertheless, it is important to ensure that the data is accessible to other researchers so they can do experiments building on the work. Part of this requires that conformation to standard formats and conventions used by other groups. It is more important that to understand the format exists and that it encodes parent-child biological relationships, than the precise specifics of what each column means.
{: .tip}

### 3. GenBank

GenBank files are a fixed-width format which displays a “flat” gene model and lacks any way to represent the hierarchical relationships that are biologically relevant. There are a few major regions of a GenBank file:
> * The header (starting with LOCUS…); gives information on
>    > 1. Sequence ID (BK000583 in the example shown)
>    > 2. Genome or chromosome length (41724 bp)
>    > 3. Annotation set version (1, from VERSION BK000583.1)
>    > 4. References

![](../../images/getting-started-with-apollo-screenshots/(10) GenBank Header P22.png)

> * The feature table (starting with FEATURES) usually begins with a “source” type feature, which contains information about the chromosome/genome. Features consist of a feature type key on the left, and key-value pairs on the right formatted as */key=“Value…”

![](../../images/getting-started-with-apollo-screenshots/(11) GenBank Features P22.png)

> * The sequence data, which is displayed as six separated columns of ten characters each, with the sequence index annotated on the left.

![](../../images/getting-started-with-apollo-screenshots/(12) GenBank Sequence Data P22.png)

![](../../images/getting-started-with-apollo-screenshots/(13) Compare/Contrast File Formats.png)

# Annotation Within Apollo

Continuing on to actually using Apollo, this section will go through an example annotation. Additionally, characteristics of the program will be described to assist in the navigation of the program and allow greater ease of use. There are two primary components to annotation:

1. Structural annotation, which consists of locations of genomic features, such as genes and terminators. Several gene callers will identify possible genes in the phage genome. Putative genes in Apollo will be annotated based on these results.
2. Functional annotation, which entails identifying possible gene functions base on multiple sources of evidence.

<!-- LINK STRUCTURAL AND FUNCTIONAL TUTORIALS UPON COMPLETION. -->

## Apollo in Galaxy - General Use

The CPT developed a tool called JBrowse-in-Galaxy (JiG), which allows the building of JBrowse instances within Galaxy; this contrasts with how JBrowse instances are traditionally configured, through a complex and annual process at the command line. 

![](../../images/getting-started-with-apollo-screenshots/(14) JBrowse in Galaxy.png)

The CPT uses JBrowse as a tool for displaying the results of a bioinformatic analysis in a standardized way; instead of having to digest and understand 20+ different report formats, images, output files, tables, etc., all of our analysis is presented as easy-to-grasp features in evidence tracks. As its input, Apollo takes complete JBrowse instances. To view any data in Apollo, a JBrowse instance needs to be configured first. On the far left side of the Galaxy web page is a “Tools” column with a search bar. Search “JBrowse genome browser,” and click on the synonymous link underneath “CPT: Genomic Viz.”

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Step1
> 2. **My Tool** {% icon tool %} with the following parameters
>   - *"param1"*: the file `myfile`
>   - *"param2"*: `42`
>   - *"param3"*: `Yes`
>
> 3. **My Tool** {% icon tool %} with the following parameters
>   - {% icon param-text %} *"My text parameter"*: `my value`
>   - {% icon param-file %} *"My input file"*: `my file`
>   - {% icon param-files %} *"My multiple file input or collection"*: `my collection`
>   - {% icon param-select %} *"My select menu"*: `my choice`
>   - {% icon param-check %} *"My check box"*: `yes`
>
>    > ### {% icon question %} Questions
>    >
>    > 1. Question1?
>    > 2. Question2?
>    >
>    >    > ### {% icon solution %} Solution
>    >    >
>    >    > 1. Answer for question1
>    >    > 2. Answer for question2
>    >    >
>    >    {: .solution}
>    >
>    {: .question}
>
> 3. Step3
{: .hands_on}

# Conclusion
{:.no_toc}

Conclusion about the technical key points. And then relation between the techniques and the biological question to end with a global view.
