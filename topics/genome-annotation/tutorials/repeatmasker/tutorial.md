---
layout: tutorial_hands_on

title: Masking repeats with RepeatMasker
zenodo_link: https://zenodo.org/record/7085837
tags:
  - eukaryote
questions:
  - How to mask repeats in a genome?
  - What is the difference between hard and soft-masking?
objectives:
  - Use Red and RepeatMasker to soft-mask a newly assembled genome
time_estimation: 1H
level: Introductory
key_points:
  - RepeatMasker can be used to soft-mask a genome
  - It is an essential first step before running structural annotation pipelines
contributors:
  - abretaud
  - alexcorm
  - lleroi
  - r1corre
  - stephanierobin
  - erasmusplus
  - jkreplak

abbreviations:
  SINEs: Short Interspersed Nuclear Elements
  LINEs: Long Interspersed Nuclear Elements

follow_up_training:
- type: "internal"
  topic_name: genome-annotation
  tutorials:
    - funannotate
subtopic: eukaryote
priority: 1
---

# Introduction



When you assemble a new genome, you get its full sequence in FASTA format, in the form of contigs, scaffolds, or even whole chromosomes if you are lucky. However genomes, in particular for eukaryote organisms, contain a varying but significant proportion of repeated elements all along the sequence. These elements belong to different classes, including:

- Tandem repeats: small sequences (<60 base pairs) repeated next to each other, found in many places in the genome, in particular centromeres and telomeres
- Interspersed repeats: sequences repeated in distant positions, including transposons, {SINEs} or {LINEs}

These repeats are interesting on their own: they can originate from transposons or viral insertions, and they can have direct effects on the expression of genes. But they are also the source of a lot of trouble when you work on genomics data. First when sequencing a genome, assembly tools often have problems reconstructing the genome sequence in regions containing repeats (in particular when repeats are longer than the read size). Then, when you have a good assembly, you want to annotate it to find the location of genes. Unfortunately annotation tools have trouble identifying gene locations in regions rich in repeats.

The aim of repeat masking is to identify the location of all repeated elements along a genome sequence. Other tools (like annotation pipelines) can then take this information into account when producing their results.

The output of repeat masking tools is most often composed of a fasta file (with sometimes a GFF file containing the position of each repeat). There is two types of masking, producing slightly different fasta output:

- Soft-masking: repeat elements are written in lower case
- Hard-masking: repeat elements are replaced by stretches of the letter N

Normal (non-repeated) sequences are always kept in uppercase. Doing hard-masking is destructive because you lose large parts of the sequence which are replaced by stretches of N. If you want to perform an annotation, it is best to choose soft-masking.

We call this operation "masking" because, by making repeats lowercase, or replacing them with Ns, you make them "invisible" by annotation tools (they are written to mostly ignore the regions marked like this).

Multiple tools exist to perform the masking: [RepeatMasker](https://www.repeatmasker.org/), [RepeatModeler](https://www.repeatmasker.org/RepeatModeler/), [REPET](https://urgi.versailles.inra.fr/Tools/REPET), ... Each one have specificities: some can be trained on specific genomes, some rely on existing databases of repeated elements signatures ([Dfam](https://www.dfam.org/), [RepBase](https://www.girinst.org/repbase/)).

In this tutorial you will learn how to soft-mask the genome sequence of a small eukaryote: Mucor mucedo (a fungal plant pathogen). You can learn how this genome sequence was assembled by following the [Flye assembly tutorial]({% link topics/assembly/tutorials/flye-assembly/tutorial.md %}). We will use two different tools, Red and RepeatMasker, which are probably two of the simplest solutions giving an acceptable result before annotating the genome in the [Funannotate annotation tutorial]({% link topics/genome-annotation/tutorials/funannotate/tutorial.md %}).

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get data

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial
>
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>     ```
>     https://zenodo.org/api/files/debdbbfd-4739-4f2d-bb79-814ac032c8b5/genome_raw.fasta
>     https://zenodo.org/api/files/debdbbfd-4739-4f2d-bb79-814ac032c8b5/Muco_library_RM2.fasta
>     https://zenodo.org/api/files/debdbbfd-4739-4f2d-bb79-814ac032c8b5/Muco_library_EDTA.fasta
>     ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
{: .hands_on}

# Soft-masking using Red

First let's try Red, a tool than can mask repeats *de novo*. For that, select the input assembly in fasta format.

> <comment-title>*Ab initio* tool</comment-title>
>
> Red is an *ab initio* tool, it means that it will try to predict repeat elements using only the genomic sequence. It's perfect when you know nothing about the organism that you are working on.
>
{: .comment}

> <hands-on-title>Hands-on</hands-on-title>
>
> 1. {% tool [Red](toolshed.g2.bx.psu.edu/repos/iuc/red/red/2018.09.10+galaxy1) %} with the following parameters:
>   - {% icon param-file %} *"Genome sequence to mask"*: `genome_raw.fasta` (Input dataset)
>
{: .hands_on}

Red produces 2 output files :

- A *fasta* file: this is the soft-masked genome that you can use for future analysis. If you display it, you will notice that some portions of the sequence are in lowercase: these are the regions that were identified as repeats
- A *bed* file: this one contains the coordinates on the genome of each repeated loci

> <question-title></question-title>
>
> What proportion of the whole genome sequence is masked?
>
> > <solution-title></solution-title>
> >
> > You need to click on the {% icon galaxy-info %} on one of the output. You should find it at the end of the extended `Tool Standard Output` in **Job Information**. It should be ~ 30.62%
> >
> {: .solution}
>
{: .question}

> <question-title></question-title>
>
> How to hard-mask a genome with **Red** ?
>
> > <solution-title></solution-title>
> >
> > As you can see **Red** has no option to hard-mask your genome. However, one of the output is a `bed` file, so you can use {% tool bedtools MaskFastaBed %} to replace repeated regions with stretches of N:
> >
> > > <hands-on-title>Hands-on</hands-on-title>
> > >
> > > 1. {% tool [bedtools MaskFastaBed](toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_maskfastabed/2.30.0) %} with the following parameters:
> > >   - {% icon param-file %} *"BED/bedGraph/GFF/VCF/EncodePeak file"*: `Red on data` (`bed` file produced by red)
> > >   - {% icon param-file %} *"FASTA file"*: `genome_raw.fasta`
> > >
> > {: .hands_on}
> >
> {: .solution}
>
{: .question}

Red uses only the sequence of the genome to detect repeated regions, and does not provide a detailed classification of the detected repeats. Let's use another tool that works differently: RepeatMasker.

# Soft-masking using RepeatMasker

Let's run RepeatMasker, by selected the input assembly in fasta format. We select the soft-masking option, and we choose to use the Dfam database.

> <comment-title>Choosing the right species</comment-title>
>
> We select the `Human (Homo sapiens)` species here, even though we are masking a fungi genome. It means RepeatMasker will identify very common repeats found in many organisms. For more precise results, you can consider selecting a species closer to the one you analyse in the drop down list, or using other more advanced tools like RepeatModeler.
>
{: .comment}

> <hands-on-title>Hands-on</hands-on-title>
>
> 1. {% tool [RepeatMasker](toolshed.g2.bx.psu.edu/repos/bgruening/repeat_masker/repeatmasker_wrapper/4.1.2-p1+galaxy1) %} with the following parameters:
>   - {% icon param-file %} *"Genomic DNA"*: `genome_raw.fasta` (Input dataset)
>   - *"Repeat library source"*: `DFam (curated only, bundled with RepeatMasker)`
>     - *"Select species name from a list?"*: `Yes`
>       - *"Species"*: `Human (Homo sapiens)`
>   - *"Output annotation of repeats in GFF format"*: `Yes`
>   - *"Perform soft-masking instead of hard-masking"*: `Yes`
>
{: .hands_on}

RepeatMasker produces 4 output files:

- `masked sequence`: this is the fasta file that you will use for future analysis. If you display it, you will notice that some portions of the sequence are in lowercase: these are the regions that were identified as repeats.
- `repeat statistics`: this one contains some statistics on the number of repeats found in each category, and the total number of base pairs masked.
- `output log`: this is a tabular file listing all repeats.
- `repeat catalogue`: this one contains the list of all repeat sequences that were identified, with their position, and their similarity with known repeats from the Dfam database.
- `repeat annotation` : this one contains the coordinate of each repeat element in GFF2 format.

> <question-title></question-title>
>
> What proportion of the whole genome sequence is masked?
>
> > <solution-title></solution-title>
> >
> > You should find it in the `repeat statistics` output. It should be ~2.41%.
> {: .solution}
>
{: .question}

As we have used a generic species (Human), we only identified the most common and simple repeats, not very specific to this species. If you compare with Red results, your are missing at least ~28% of repeated content in the genome. However, RepeatMasker gives interesting information about repeat classification which could be interesting for future analysis.

To boost RepeatMasker performance, we need a tailored repeat library for *Mucor mucedo*. This step can take from a few hours to a few days and a large number of tools could be used. We pre-computed two librairies :

- `Muco_library_RM2.fasta` using [RepeatModeler](https://doi.org/10.1073/pnas.1921046117)
- `Muco_library_EDTA.fasta` using [EDTA](https://doi.org/10.1186/s13059-019-1905-y)


> <hands-on-title>Hands-on</hands-on-title>
>
> 1. {% tool [RepeatMasker](toolshed.g2.bx.psu.edu/repos/bgruening/repeat_masker/repeatmasker_wrapper/4.1.2-p1+galaxy1) %} with the following parameters:
>   - {% icon param-file %} *"Genomic DNA"*: `genome_raw.fasta` (Input dataset)
>   - *"Repeat library source"*: `Custom library of repeats`
>     - *"Custom library of repeats"*
>       - *"One of the two pre-computed libraires"*: `Muco_library_RM2.fasta` or `Muco_library_EDTA.fasta`
>   - *"Output annotation of repeats in GFF format"*: `Yes`
>   - *"Perform soft-masking instead of hard-masking"*: `Yes`
{: .hands_on}

> <question-title></question-title>
>
> Compare the different `repeat statistics` files produced, what is the highest library for RepeatMasker?
>
> > <solution-title></solution-title>
> >
> > The RepeatModeler library seems to have the highest percentage of repeats found with ~ 34.89%. It could be explained as RepeatModeler is specifically made to work with RepeatMasker.
> {: .solution}
>
{: .question}

Other tools might mask a greater proportion of the genome, at the cost of a more complex workflow with training steps. But masking more isn't always a positive results! In fact, large family of genes could be considered as a repeat by some tools or certain library. Only a manual curation can correct those mistakes, but this result is sufficient to perform an annotation by following the [Funannotate annotation tutorial]({% link topics/genome-annotation/tutorials/funannotate/tutorial.md %}).

# Conclusion



By following this tutorial you have learn how to mask an eukaryotic genome using Red and RepeatMasker, after assembling ([Flye assembly tutorial]({% link topics/assembly/tutorials/flye-assembly/tutorial.md %})) and before annotating it ([Funannotate annotation tutorial]({% link topics/genome-annotation/tutorials/funannotate/tutorial.md %})).

Often times, annotation tools prefer to use soft-masked genomes, as they primarily search for genes in non repeated regions, but tolerate that some genes overlap partially with these regions.
