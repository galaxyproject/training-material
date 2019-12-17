---
layout: tutorial_hands_on

title: Isoform-level ribosome occupancy estimation
zenodo_link: 'https://figshare.com/s/fcd05b448006437143f0'
questions:
- What are difficulties in the process to calculate gene isoforms expression at the translational level?
- What can we learn from the from gene expression at a more refine level?
objectives:
- Learning how to calculate gene expression at the isoform level
- Aware the fine tuning mechanism that occurred in the gene regulation at the isoform level
- You will know the complexity and precision of gene regualtion at the translational level
time_estimation: ''
key_points:
- You can learn about the complexity of gene expression due to various isoforms of them through this tutorial, hence the quatification of genes and their isoforms expression accurately is very necessary. 
contributors:
- ldyang14
- IceApink
---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

A gene usually contains many isoforms due to alternative splicing sites. Although these isoforms have overlapping regions, they can still perform diverse functions by translating to various polypeptide chains (As shown below). In general, the quantitative of gene expression is on account of the counts of reads on a gene, hence we neglected the question that gene may express a variety of isoforms. Thereby the more accurate and precise details of cell processes and molecular mechanisms were ignored.

![Gene isofroms](../../images/isoform-level/gene_isoforms.png "Gene isoforms (cited from {% cite aguiar2018bayesian %})")

We can know the precise location of reads through paired-end sequencing, thereby we should reconstruct the actual abundance of a gene and its isoforms in a more refined level with the help of some custom tools. For enhancing the reliability of analysis results, we generally remove ambiguously aligned reads such as multiple mapped reads according to the experience in the process of data analysis for transcriptome. If not the estimation of actual gene abundance would be affected. However, if you checked the mapping ratio of translatomics data, you will be amazed to find that multiple mapped reads occupy a large proportion compared to transcriptomics. The reason is that the length of reads from translatome is usually about 30nt owning to enclose by the ribosomes, thus a large number of reads produced by Ribo-Seq will align to repeat regions on the genome or overlapping areas of different gene isoforms. Therefore, we can not migrate the tools and strategies that take in the analysis of transcriptome to apply to the translatome. 

Therefore, a series of tools have been developed for quantification of genes and its isoforms abundance more accurately at the translatome level, for now, thus to explore gene function and molecular mechanisms more precisely. 




> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Import data 

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>    https://ndownloader.figshare.com/files/20042474?private_link=fcd05b448006437143f0
>    https://ndownloader.figshare.com/files/20042471?private_link=fcd05b448006437143f0
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
{: .hands_on}

# Calculate gene expression at the isoform level with **Ribomap**

> ### {% icon hands_on %} Hands-on: Calculate gene expression at the isoform level with **Ribomap**
>
> 1. **Ribomap** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input RNA-seq read fastq.gz file for transcript abundance estimation"*: `output` (Input dataset)
>    - {% icon param-file %} *"Input ribosome profiling (riboseq) read fastq.gz file"*: `output` (Input dataset)
>    - {% icon param-file %} *"Input trascriptome reference fasta file"*: `output` (Input dataset)
>    - *"Use contaminant sequence fasta file?"*: `no`
>    - *"Use coding sequence (CDS) range for all transcripts?"*: `no`
>    - *"The linker sequence attached to the 5' end of the Ribo-Seq reads."*:  `CTGTAGGCACCATCAAT`  
>    - *"Minimun read length to keep for downstream analysis."*: `25`
>    - *"Maximum riboseq read length to keep for downstream analysis."*: `34`
{: .hands_on}

# Conclusion

{:.no_toc}

In short,  we can know more details about mechanism of gene regulation through dectecting expression of gene isoforms. However, we have to face certain new challenges as higher resolution provided by ribosome profiling. For example, the determination of reads mapped on exon-exon junctions. Hence, a certain of new tools have been developing to resolve these issues for helping us to explore and understand the physiological process. 