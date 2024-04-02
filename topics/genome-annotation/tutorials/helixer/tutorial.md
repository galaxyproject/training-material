---
layout: tutorial_hands_on

title: Genome annotation with Helixer
zenodo_link: https://zenodo.org/record/7867921
tags:
  - eukaryota
  - helixer
  - jbrowse

edam_ontology:
- topic_0196 # Sequence Assembly
- topic_0622 # Genomics

questions:
    - How to annotate an eukaryotic genome with Helixer?
    - How to evaluate and visualize annotated genomic features?

objectives:
    - Load genome into Galaxy
    - Annotate genome with Helixer
    - Evaluate annotation quality with BUSCO and Compleasm
    - View annotations in JBrowse

time_estimation: 4h
level: Intermediate
key_points:
    - Helixer allows to perform structural annotation of an eukaryotic genome
    - BUSCO ans Compleasm allow to inspect the quality of an annotation
contributions:
  authorship:
    - rlibouba
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


In this tutorial we will use a software tool called Helixer to annotate the genome sequence of a small eukaryote: [*Mucor mucedo*](https://en.wikipedia.org/wiki/Mucor_mucedo) (a fungal plant pathogen).

[Helixer](https://github.com/weberlab-hhu/Helixer) is a software package providing a framework for the development and use of a cross-species deep learning model that dramatically improves performance. The software enables models to be configured and trained for *ab initio* prediction of gene structure. In other words, to identify which base pairs in a genome belong to the UTR/CDS/Intron genes.

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

- The **genome sequence** in fasta format. For best results, the sequence should be soft-masked beforehand. You can learn how to do it by following the [RepeatMasker tutorial]({% link topics/genome-annotation/tutorials/repeatmasker/tutorial.md %}). For this tutorial, we will try to annotate the genome assembled in the [Flye assembly tutorial]({% link topics/assembly/tutorials/flye-assembly/tutorial.md %}).

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

We can run [**Helixer**](https://github.com/weberlab-hhu/Helixer) to perfom the structural annotation of the genome.

We need to input the genome sequence. We also specify the lineage. There are 4 differents lineages: invertebrate, vertebrate, land plant ans fungi. As an option, we can also enter the species name.

Depending on the lineage selected, default values for certain parameters change. These parameters are:
- 


> <hands-on-title></hands-on-title>
>
> 1. {% tool [Helixer](toolshed.g2.bx.psu.edu/repos/genouest/helixer/helixer/0.3.2+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Genomic sequence"*: `genome_masked.fasta` (Input dataset)
>    - In *"Available lineages"*: "*select*" `fungi`
>    - In *"Species name"*: `Mucor mucedo`

{: .hands_on}

> <comment-title>on parameters</comment-title>
>
> - Specie name is not a default setting. You don't have to enter every time.
> - Depending on the lineage, certain parameters have different default values. The parameters concerned are: subsequence-length, overlap-offset and overlap-core-length.
> - This is due in particular to the size of the genomes.
> - Indeed, it is recommended to increase the value of subsequence-length when genomes contain large genes. This is particularly important for vertebrates and invertebrates.
> - Recommendations have therefore been defined.
On Galaxy, if you wish to modify these default values, you can do so by entering your value in the parameters: subsequence length, overlep offset and overlap core legth.
{: .comment}


> <comment-title>Don't wait</comment-title>
>
> This step will take a bit of time to run. While it runs, we can already schedule the following functional annotation steps. Galaxy will run them automatically as soon as the structural annotation is ready.
{: .comment}

Helixer produces one output dataset: a .GFF3 file. The GFF3 format is a standard bioinformatics format for storing denomination annotations. Each row describes a genomic entity, with columns detailling ist identifier, location, score and other attributes. Thes fileexons and ohter genomic features.s are widely used yo exchange and store information on genes.

## Evaluate the quality of structural annotation with **Genome Annotation Statistics**

Genome Annotation Statistics is a program designed to analyze and provide statistics on genomic annotations. This software performs its analyses from a GFF3 file.

> <hands-on-title></hands-on-title>
>
> 1. {% tool [Genome Annotation Statistics](toolshed.g2.bx.psu.edu/view/iuc/jcvi_gff_stats/jcvi_gff_stats/0.8.4) %} with the following parameters:
>    - {% icon param-file %} *"Annotation to analyse "*: `GFF3 file` (Output of **Helixer**)
>    - In *"Reference genome"*: select `Use a genome from history`
>    - In *"Corresponding genome sequence"*: `genome_masked.fasta` (Input dataset)`
{: .hands_on}

Two output files are generated:
- a file containing graphs in pdf format
- a summary in txt format



## Evaluation with **Busco**

[BUSCO](http://busco.ezlab.org/) (Benchmarking Universal Single-Copy Orthologs) is a tool allowing to evaluate the quality of a genome assembly or of a genome annotation. By comparing genomes from various more or less related species, the authors determined sets of ortholog genes that are present in single copy in (almost) all the species of a clade (Bacteria, Fungi, Plants, Insects, Mammalians, ...). Most of these genes are essential for the organism to live, and are expected to be found in any newly sequenced and annotated genome from the corresponding clade. Using this data, BUSCO is able to evaluate the proportion of these essential genes (also named BUSCOs) found in a set of (predicted) transcript or protein sequences. This is a good evaluation of the "completeness" of the annotation.

> <hands-on-title></hands-on-title>
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.4.6+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `masked genome` (Input dataset)
>    - *"Mode"*: `Genome assemblies (DNA)`
>    - *"Auto-detect or select lineage?"*: `Select lineage`
>        - *"Lineage"*: `Mucorales`
>    - In *"Advanced Options"*:
>        - *"Which outputs should be generated"*: select all outputs
>
{: .hands_on}

Several output files will be generated:
- short summary : statistical summary of the quality of genomic assembly or annotation, including total number of genes evaluated, percentage of complete genes, percentage of partial genes, etc.
- full table : list of universal orthologs found in the assembled or annotated genome, with information on their completeness, location in the genome, quality score, etc.
- missing buscos : list of orthologs not found in the genome, which may indicate gaps in assembly or annotation.
- summary image : graphics and visualizations to visually represent the results of the evaluation, such as bar charts showing the proportion of complete, partial and missing genes.
- GFF : contain information on gene locations, exons, introns, etc.


## Evaluation with **Compleasm**

[Compleam](https://academic.oup.com/bioinformatics/article/39/10/btad595/7284108) is described as a faster and more accurate reimplementation of Busco. It serves as an efficient tool for assessing the completeness of genome assemblies. Compleam utilizes the protein-genome aligner from Miniprot and the conserved orthologous genes from Busco. On human assemblies, Compleam has been demonstrated to be 14 times faster than Busco, while maintaining precision, which increased from 95.7% to 99.6%. (mettre ref article)

> <hands-on-title></hands-on-title>
>
> 1. {% tool [Compleasm](toolshed.g2.bx.psu.edu/repos/iuc/compleasm/compleasm/0.2.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `masked genome` (Input dataset)
>    - In *"Which outputs should be generated"*: select all outputs
>    - *"The mode of evaluation"*: `BUSCO`
>    - *"Choose the BUSCO database to be used"*: select `Busco v5 Linage Datasets (all+2024+03+21)`
>    - In *"lineage?"*: `Mucorales`
>
{: .hands_on}

Several output files will be generated:
- full table busco
- full table
- miniprot (gff3 format)
- translated proteins (fasta format)

# Visualisation with a genome browser

You can visualize the annotation generated using a genomic browser like [JBrowse](https://jbrowse.org/jbrowse1.html). This browser enables you to navigate along the chromosomes of the genome and view the structure of each predicted gene.


> <hands-on-title></hands-on-title>
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

Congratulations on reaching the end of this tutorial! You now know how to perform a structural annotation of a new eukaryotic genome, using Helixer. You've learned how to visualize your new annotation using JBrowse.

If you're interested in using other genomes, take a look at the tutorial [Masking repeats with RepeatMasker](https://training.galaxyproject.org/training-material/topics/genome-annotation/tutorials/repeatmasker/tutorial.html). This tutorial explains how to mask a genome. 

If you'd like to complete this annotation, we recommend you follow the tutorial on annotation with [Funannotate](https://training.galaxyproject.org/training-material/topics/genome-annotation/tutorials/funannotate/tutorial.html) and take a look at the sections on functional annotation with EggNOG Mapper and InterProScan. The datasets are the same.


