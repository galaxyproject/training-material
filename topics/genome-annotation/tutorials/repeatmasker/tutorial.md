---
layout: tutorial_hands_on

title: Masking repeats with RepeatMasker
zenodo_link: https://zenodo.org/record/5726723
tags:
  - eukaryote
questions:
  - How to mask repeats in a genome?
  - What is the difference between hard and soft masking?
objectives:
  - Use RepeatMasker to soft mask a newly assmbled genome
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

abbreviations:
    SINEs: Short Interspersed Nuclear Elements
    LINEs: Long Interspersed Nuclear Elements

follow_up_training:
- type: "internal"
  topic_name: genome-annotation
  tutorials:
    - funannotate
---


# Introduction
{:.no_toc}

When you assemble a new genome, you get its full sequence in FASTA format, in the form of contigs, scaffolds, or even whole chromosomes if you are lucky. However genomes, in particular for eukaryote organisms, contain a varying but significant proportion of repeated elements all along the sequence. These elements belong to different classes, including:

- Tandem repeats: small sequences (<60 base pairs) repeated next to each other, found in many places in the genome, in particular centromeres and telomeres
- Interspersed repeats: sequences repeated in distant positions, including transposons, {SINEs} or {LINEs}

These repeats are interesting on their own: they can originate from transposons or viral insertions, and they can have direct effects on the expression of genes. But they are also the source of a lot of trouble when you work on genomics data. First when sequencing a genome, assembly tools often have problems reconstructing the genome sequence in regions containing repeats (in particular when repeats are longer than the read size). Then, when you have a good assembly, you want to annotate it to find the location of genes. Unfortunately annotation tools have trouble identifying gene locations in regions rich in repeats.

The aim of repeat masking is to identify the location of all repeated elements along a genome sequence. Other tools (like annotation pipelines) can then take this information into account when producing their results.

The output of repeat masking tools is most often composed of a fasta file (with sometimes a GFF file containing the position of each repeat). There is two types of masking, producing slightly different fasta output:

- Soft masking: repeat elements are written in lower case
- Hard masking: repeat elements are replaced by stretches of the letter N

Normal (non-repeated) sequences are always kept in uppercase. Doing hard masking is destructive because you lose large parts of the sequence which are replaced by stretches of N. If you want to perform an annotation, it is best to choose soft masking.

We call this operation "masking" because, by making repeats lowercase, or replacing them with Ns, you make them "invisible" by annotation tools (they are written to mostly ignore the regions marked like this).

Multiple tools exist to perform the masking: [RepeatMasker](https://www.repeatmasker.org/), [RepeatModeler](https://www.repeatmasker.org/RepeatModeler/), [REPET](https://urgi.versailles.inra.fr/Tools/REPET), ... Each one have specificities: some can be trained on specific genomes, some rely on existing databases of repeated elements signatures ([Dfam](https://www.dfam.org/), [RepBase](https://www.girinst.org/repbase/)).

In this tutorial you will learn how to soft mask the genome sequence of a small eukaryote: Mucor mucedo (a fungal plant pathogen). You can learn how this genome sequence was assembled by following the [Flye assembly tutorial]({% link topics/assembly/tutorials/flye-assembly/tutorial.md %}). We will use RepeatMasker, which is probably the simplest solution giving an acceptable result before annotating the genome in the [Funannotate annotation tutorial]({% link topics/genome-annotation/tutorials/funannotate/tutorial.md %}).

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    https://zenodo.org/api/files/71333591-99bd-4d99-bbdf-664cc18fd422/genome_raw.fasta
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
{: .hands_on}

# Soft-masking using RepeatMasker

Let's run RepeatMasker, by selected the input assembly in fasta format. We select the soft masking option, and we choose to use the Dfam database.

> ### {% icon comment %} Choosing the right species
>
> We select the `Human (Homo sapiens)` species here, even though we are masking a fungi genome. It means RepeatMasker will identify very common repeats found in many organisms. For more precise results, you can consider selecting a species closer to the one you analyse in the drop down list, or using other more advanced tools like RepeatModeler.
{: .comment}

> ### {% icon hands_on %} Hands-on
>
> 1. {% tool [RepeatMasker](toolshed.g2.bx.psu.edu/repos/bgruening/repeat_masker/repeatmasker_wrapper/4.1.2-p1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Genomic DNA"*: `genome_raw.fasta` (Input dataset)
>    - *"Repeat library source"*: `DFam (curated only, bundled with RepeatMasker)`
>        - *"Select species name from a list?"*: `Yes`
>            - *"Species"*: `Human (Homo sapiens)`
>    - *"Perform softmasking instead of hardmasking"*: `Yes`
>
{: .hands_on}

RepeatMasker produces 4 output files:

- `masked sequence`: this is the fasta file that you will use for future analysis. If you display it, you will notice that some portions of the sequence are in lowercase: these are the regions that were identified as repeats.
- `repeat statistics`: this one contains some statistics on the number of repeats found in each category, and the total number of base pairs masked.
- `output log`: this is a tabular file listing all repeats.
- `repeat catalogue`: this one contains the list of all repeat sequences that were identified, with their position, and their similarity with known repeats from the Dfam database.

> ### {% icon question %} Question
>
> What proportion of the whole genome sequence is masked?
>
> > ### {% icon solution %} Solution
> >
> > You should find it in the `repeat statistics` output. It should be ~2.41%.
> >
> {: .solution}
>
{: .question}

As we have used a generic species (Human), we only identified the most common repeats, not very specific to this species. Other tools might mask a greater proportion of the genome, at the cost of a more complex workflow with training steps. But this result is sufficient to perform an annotation by following the [Funannotate annotation tutorial]({% link topics/genome-annotation/tutorials/funannotate/tutorial.md %}).

# Conclusion
{:.no_toc}

By following this tutorial you have learn how to mask an eukaryotic genome using RepeatMasker, after assembling ([Flye assembly tutorial]({% link topics/assembly/tutorials/flye-assembly/tutorial.md %})) an before annotating it ([Funannotate annotation tutorial]({% link topics/genome-annotation/tutorials/funannotate/tutorial.md %})).

Often times, annotation tools prefer to use soft masked genomes, as they primarily search for genes in non repeated regions, but tolerate that some genes overlap partially with these regions.
