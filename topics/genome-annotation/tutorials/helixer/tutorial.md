---
layout: tutorial_hands_on

title: Genome annotation with Helixer
zenodo_link: https://zenodo.org/record/7867921
tags:
  - eukaryota
  - helixer
  - jbrowse1

edam_ontology:
- topic_0362 # Genome annotation

questions:
    - How to annotate an eukaryotic genome with Helixer?
    - How to evaluate and visualize annotated genomic features?

objectives:
    - Load genome into Galaxy
    - Annotate genome with Helixer
    - Evaluate annotation quality with BUSCO
    - View annotations in JBrowse

time_estimation: 4h
level: Intermediate
key_points:
    - Helixer allows to perform structural annotation of an eukaryotic genome
    - BUSCO allows to inspect the quality of an annotation
contributions:
  authorship:
    - rlibouba
    - abretaud
  funding:
    - eurosciencegateway
abbreviations:
    GPU: Graphics Processing Unit
    UTR: Untranslated region

requirements:
 - type: internal
   topic_name: genome-annotation
   tutorials:
     - repeatmasker

subtopic: eukaryote
---

Annotating the eukaryotic genome represents a somewhat more complex challenge than that of prokaryotes, mainly due to the generally larger size of eukaryotic genomes and their greater number of genes, but also to the complexity of eukaryotic genes structure (e.g. exons and {UTR}). This annotation can be carried out at different levels of precision, ranging from simple identification of coding and non-coding parts to detailed structural labeling, including for example the precise location of exons, introns and other regulatory elements.


In this tutorial we will use a software tool called Helixer to annotate the genome sequence of a small eukaryote: [*Mucor mucedo*](https://en.wikipedia.org/wiki/Mucor_mucedo) (a fungal plant pathogen).

[Helixer](https://github.com/weberlab-hhu/Helixer) is an annotation software with a new and different approach: it performs evidence-free predictions (no need for RNASeq data or sequence aligments), using {GPU}, with a much faster execution time. The annotation is based on the development and use of a cross-species deep learning model. The software is used to configure and train models for *ab initio* prediction of gene structure. In other words, it identifies the base pairs in a genome that belong to the UTR/CDS/Intron genes.

In this tutorial, you'll learn how to perform a structural annotation of the genome and how to assess its quality.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data upload

To annotate our genome using Helixer, we will use the following files:

- The **genome sequence** in fasta format. For best results, the sequence should be soft-masked beforehand. You can learn how to do it by following the [RepeatMasker tutorial]({% link topics/genome-annotation/tutorials/repeatmasker/tutorial.md %}). For this tutorial, we will try to annotate the genome assembled in the [Flye assembly tutorial]({% link topics/assembly/tutorials/flye-assembly/tutorial.md %}) and already masked for you using RepeatMasker.

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
>    https://zenodo.org/record/7867921/files/genome_masked.fasta
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
{: .hands_on}

# Structural annotation

We can run [**Helixer**](https://github.com/weberlab-hhu/Helixer) to perform the structural annotation of the genome.

We need to input the genome sequence we want to annotate.

We also need to choose between 4 different lineages: *invertebrate*, *vertebrate*, *land plant* or *fungi*. Select the one that fits the best to the species you're studying: *fungi* in our case. Helixer is shipped with these 4 models that were trained specifically to annotate genes from each of these lineages. Advanced users can upload their own lineage model in .h5 format with the *"Lineage model"* option.

As an option, we can also enter a species name.


> <hands-on-title></hands-on-title>
>
> 1. {% tool [Helixer](toolshed.g2.bx.psu.edu/repos/genouest/helixer/helixer/0.3.3+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Genomic sequence"*: `genome_masked.fasta` (Input dataset)
>    - In *"Available lineages"*: "*select*" `fungi`
>    - In *"Species name"*: `Mucor mucedo`
>
{: .hands_on}

> <comment-title>Advanced parameters</comment-title>
>
> Depending on the lineage, the parameters *"Subsequence length"*, *"Overlap offset"* and *"Overlap corelength"* are adjusted to corresponding default values (listed in the help of each option).
>
> This is due in particular to the size of the genomes. Indeed, it is recommended to increase the value of *"Subsequence length"* for genomes containing large genes. This is particularly important for vertebrates and invertebrates.
>
> The default values used by Galaxy are the ones recommended by Helixer authors. If you wish to modify these default values, you can do so by entering your values in the *"Subsequence length"*, *"Overlap offset"* and *"Overlap corelength"* parameters.
>
{: .comment}

> <comment-title>Don't wait</comment-title>
>
> This step can take a bit of time to run: although Helixer runs much faster than many other annotation tools (typically <20min for this tutorial), it requires specific hardware (GPU) that is often available in limited quantity on computing systems. It means your job can be placed in queue for a longer time than a more standard Galaxy job.
>
> While it runs, we can already schedule the following steps. Galaxy will run them automatically as soon as the Helixer annotation is ready.
{: .comment}

Helixer produces a single output dataset: a GFF3 file. The GFF3 format is a standard bioinformatics format for storing genome annotations. Each row describes a genomic entity, with columns detailing its identifier, location, score and other attributes.

# Quality evaluation

## General statistics

Genome Annotation Statistics is a program designed to analyze and provide statistics on genomic annotations. This software performs its analyses from a GFF3 file.

> <hands-on-title>Genome Annotation Statistics</hands-on-title>
>
> 1. {% tool [Genome Annotation Statistics](toolshed.g2.bx.psu.edu/view/iuc/jcvi_gff_stats/jcvi_gff_stats/0.8.4) %} with the following parameters:
>    - {% icon param-file %} *"Annotation to analyse "*: `GFF3 file` (Output of **Helixer**)
>    - In *"Reference genome"*: select `Use a genome from history`
>    - In *"Corresponding genome sequence"*: `genome_masked.fasta` (Input dataset)
{: .hands_on}

Two output files are generated:
- a file containing graphs in pdf format
- a summary in txt format

> <comment-title>What can we deduce from these results?</comment-title>
>
> - The summary file provides statistics on the genome annotation and gives a complete overview of the genomic structure and characteristics of the genes, exons and introns in the analysed genome.
> - We can see that there are 19,299 genes, 77% of which are multi-exons (i.e. 14,860) and 23% single-exons (i.e. 4,439).
> - We can obtain other information such as the average size of exons, the percentage in GC or the average size of transcripts.
>
{: .comment}

These statistics are interesting on their own: you often have a rough idea of the expected number of genes or mean length when annotating a new genome, by comparing with published similar species. You can also use them to compare the quality of annotations produced by different tools.

## Evaluation with **Busco**

[BUSCO](http://busco.ezlab.org/) (Benchmarking Universal Single-Copy Orthologs) is a tool allowing to evaluate the quality of a genome assembly or of a genome annotation. By comparing genomes from various more or less related species, the authors determined sets of ortholog genes that are present in single copy in (almost) all the species of a clade (Bacteria, Fungi, Plants, Insects, Mammalians, ...). Most of these genes are essential for the organism to live, and are expected to be found in any newly sequenced and annotated genome from the corresponding clade. Using this data, BUSCO is able to evaluate the proportion of these essential genes (also named BUSCOs) found in a set of (predicted) transcript or protein sequences. This is a good evaluation of the "completeness" of the annotation.

As an alternative for genomes only one can use [**compleasm**](https://github.com/huangnengCSU/compleasm) with the same BUSCO gene sets, as compleasm is a bit more sensitive and thus allows finding slightly more conserved genes. 

We want to run BUSCO on the protein sequences predicted from gene sequences of the Helixer annotation. So first generate these sequences:

> <hands-on-title>Extract protein sequences</hands-on-title>
>
> 1. {% tool [GFFread](toolshed.g2.bx.psu.edu/repos/devteam/gffread/gffread/2.2.1.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input GFF3 or GTF feature file"*: output of {% tool [Helixer](toolshed.g2.bx.psu.edu/repos/genouest/helixer/helixer/0.3.3+galaxy1)) %}
>    - In *"Reference Genome"* select: `From your history` (Input dataset)
>    - *"Genome Reference Fasta"*: `masked genome` (Input dataset)
>    - In *"Select fasta outputs"* select: `fasta file with spliced exons for each GFF transcript (-y)`
>    - *"full GFF attribute preservation (all attributes are shown)"*: `Yes`
>    - *"decode url encoded characters within attributes"*: `Yes`
>    - *"warn about duplicate transcript IDs and other potential problems with the given GFF/GTF records"*: `Yes`
>
{: .hands_on}

To run BUSCO on these protein sequences:

> <hands-on-title>BUSCO in proteome mode</hands-on-title>
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `gffread: pep.fa`
>    - *"Mode"*: `annotated gene sets (protein)`
>    - *"Auto-detect or select lineage?"*: `Select lineage`
>        - *"Lineage"*: `Mucorales`
>    - In *"Advanced Options"*:
>        - *"Which outputs should be generated"*: select all outputs
>
{: .hands_on}

Several output files are generated:

- short summary : statistical summary of the quality of genomic assembly or annotation, including total number of genes evaluated, percentage of complete genes, percentage of partial genes, etc.
- full table : list of universal orthologs found in the assembled or annotated genome, with information on their completeness, location in the genome, quality score, etc.
- missing buscos : list of orthologs not found in the genome, which may indicate gaps in assembly or annotation.
- summary image : graphics and visualizations to visually represent the results of the evaluation, such as bar charts showing the proportion of complete, partial and missing genes.
- GFF : contain information on gene locations, exons, introns, etc.

This gives information about the completeness of the Helixer annotation. A good idea is to compare this first result with the one you get on the initial genome sequence, and see if the annotation tool found all the genes that BUSCO finds in the raw genome sequence. So run BUSCO in genome mode:

> <hands-on-title>BUSCO in genome mode</hands-on-title>
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `masked genome` (Input dataset)
>    - *"Mode"*: `Genome assemblies (DNA)`
>    - *"Auto-detect or select lineage?"*: `Select lineage`
>        - *"Lineage"*: `Mucorales`
>    - In *"Advanced Options"*:
>        - *"Which outputs should be generated"*: select all outputs
>
{: .hands_on}

> <comment-title>What can we deduce from these results?</comment-title>
>
> - 94.6% of genes are complete, so the annotation is of high quality in terms of genomic completeness.
> - It is a little bit lower than what BUSCO is able to find in genome mode (95.7%), but the difference is quite small so Helixer seems to have generated a quite good result.
> - The duplication rate is low, with 1.3% and 1.5% of genes duplicated.
> - So the Helixer annotation looks like a good one, with high completeness and low duplication.
>
{: .comment}

# Visualisation with a genome browser

You can visualize the annotation generated using a genomic browser like [JBrowse](https://jbrowse.org/jbrowse1.html). This browser enables you to navigate along the chromosomes of the genome and view the structure of each predicted gene.

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
>                        - {% icon param-file %} *"GFF/GFF3/BED Track Data"*: `gff3` (output of **Helixer** {% icon tool %})
{: .hands_on}

Click on the newly created dataset's eye to display it. You will see a JBrowse genome browser. You can have a look at the [JBrowse tutorial]({% link topics/visualisation/tutorials/jbrowse/tutorial.md %}) for a more in-depth description of JBrowse.

# Conclusion

Congratulations on reaching the end of this tutorial! You now know how to perform a structural annotation of a new eukaryotic genome, using Helixer. And you've learned how to evaluate its quality and how to visualize it using JBrowse.

If you'd like to complete this annotation, we recommend you to follow the tutorial on [functional annotation]({% link topics/genome-annotation/tutorials/functional/tutorial.md %}) with EggNOG Mapper and InterProScan. You can follow it with the protein sequences we generated earlier with gffread.
