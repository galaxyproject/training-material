---
layout: tutorial_hands_on

title: Genome annotation with Funannotate
zenodo_link: https://zenodo.org/record/xxxx
tags:
  - eukaryote
questions:
  - How to annotate an eukaryotic genome with Funannotate?
  - How to perform functional annotation?
  - How to evaluate and visualize annotated genomic features?
  - How to format the annotation for submission at NCBI?
objectives:
  - Load genome into Galaxy
  - Annotate genome with Funannotate
  - Perform functional annotation using EggNOG-mapper and InterProScan
  - Evaluate annotation quality with BUSCO
  - View annotations in JBrowse
time_estimation: 8H
key_points:
  - Funannotate allows to perform structural annotatation of a eukaryotic genome.
  - Functional annotation can be performed using EggNOG-mapper and InterProScan.
  - BUSCO and JBrowse allow to inspect the quality of an annotation.
  - Funannotate allows to format an annotation for sumission at NCBI.
contributors:
- abretaud
- alexcorm
- lleroi
- r1corre
- stephanierobin

---


# Introduction
{:.no_toc}


Genome annotation of eukaryotes is a little more complicated than for prokaryotes: eukaryotic genomes are usually larger than prokaryotes, with more genes. The sequences determining the beginning and the end of a gene are generally less conserved than the prokaryotic ones. Many genes also contain introns, and the limits of these introns (acceptor and donor sites) are not highly conserved.

In this tutorial we will use a software tool called Funannotate ({% cite jonathan_m_palmer_2020_4054262 %}) to annotate the genome sequence of a small eukaryote: [*Mucor mucedo*](https://en.wikipedia.org/wiki/Mucor_mucedo) (a fungal plant pathogen).

As explained on [Funannotate's website](https://funannotate.readthedocs.io/), "it was originally written to annotate fungal genomes (small eukaryotes ~ 30 Mb genomes), but has evolved over time to accomodate larger genomes". As other annotation tools like Maker or Braker, it works by aligning as many evidences as possible along the genome sequence, and then reconciliating all these signals to determine probable gene structures.

The evidences can be transcript or protein sequences from the same (or closely related) organism. These sequences can come from public databases (like NR or GenBank) or from your own experimental data (transcriptome assembly from an RNASeq experiment for example). Funannotate is also able to take into account repeated elements.

Funannotate uses ab-initio predictors ([Augustus](http://bioinf.uni-greifswald.de/augustus/), [SNAP](https://github.com/KorfLab/SNAP), [glimmerHMM](https://bio.tools/glimmer-hmm), [CodingQuarry](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1344-4) and [GeneMark-ES/ET](http://exon.gatech.edu/GeneMark/) (optional due to licensing)) to improve its predictions: these software tools are able to make gene structure predictions by analysing only the genome sequence with a statistical model.

While for [Maker]({% link topics/genome-annotation/tutorials/annotation-with-maker/tutorial.md %}) you need to perform training steps for the ab-initio predictors, Funannotate is able to take care of that for you, which makes it much easier to use.

In this tutorial you will learn how to perform a structural genome annotation, and how to evaluate its quality. Then you will learn how to run functional annotation, using EggNOG-mapper and InterProScan to automatically assign names and functions to the annotated genes. And you will also learn how Funannotate can prepare files ready for submission of your annotation to the NCBI.

Finally, you will learn how to use the [JBrowse](http://jbrowse.org/) genome browser to visualise your new annotation.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data upload

To annotate our genome using Funannotate, we will use the following files:

- The **genome sequence** in fasta format. For best results, the sequence should be soft-masked beforehand. You can learn how to do it by following the [RepeatMasker tutorial](todo). For this tutorial we will try to annotate the genome assemble in the [Flye assembly tutorial](todo).
- A set of pre-aligned **RNASeq data** in BAM format. Typically, you will run a mapper like [**STAR**](https://github.com/alexdobin/STAR) ({% cite dobin2013star %}) to align on the genome RNASeq reads from multiple conditions, and merge it into a single BAM file.
- A set of **protein sequences**, like UniProt/SwissProt. It is important to have good quality, curated sequences here, that's why, by default, Funannotate will use the UniProt/SwissProt databank. In this tutorial we have prepared a subset of this databank to speed up computing, but you should use UniProt/SwissProt for real life analysis.

 Funannotate will take iinto account the position of mapped RNASeq reads, and the alignment of protein sequences on the genome sequence to determine gene positions.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo](https://doi.org/10.5281/zenodo.5653163) or from the shared data library
>
>    ```
>    https://zenodo.org/api/files/xxxx/genome_masked.fasta
>    https://zenodo.org/api/files/xxxx/rnaseq_R1.fq.gz
>    https://zenodo.org/api/files/xxxx/rnaseq_R2.fq.gz
>    https://zenodo.org/api/files/xxxx/SwissProt_subset.fasta
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
{: .hands_on}

# Preparing the genome sequence

Before annotating the genome, we want to make sure that the fasta file is properly formatted. We do it now to make sure we will not encounter unexpected errors later in the annotation process.

Funannotate provides two little tools to help us. Let's run the two tools, one after the other.

The first one (**Funannotate assembly clean**) compares all the sequences between them, and removes the shorter ones that are already included in longer ones. This is to reduce unexpected redundancy in the genome. This step is recommended only for haploid genomes (we know our organism is haploid). This first tool will also remove any suspicious sequence (like sequences made only of 1 or 2 letters, instead of the 5 expected (ATGCN).

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Funannotate assembly clean](toolshed.g2.bx.psu.edu/repos/iuc/funannotate_clean/funannotate_clean/1.8.9+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Assembly to clean"*: `assembly_masked.fasta` (Input dataset)
>
{: .hands_on}

The second tool will ensure that our fasta file is sorted, based on the length of the contigs (the longest ones first). It will also rename contigs to make sure the name are standard (they will all begin with `scaffold_`, then a number).

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Sort assembly](toolshed.g2.bx.psu.edu/repos/iuc/funannotate_sort/funannotate_sort/1.8.9+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Assembly to sort"*: `output` (output of **Funannotate assembly clean** {% icon tool %})
>
{: .hands_on}

After this step, the genome is clean, sorted, and ready for the structural annotation.


> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [RNA STAR](toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.8a+galaxy0) %} with the following parameters:
>    - *"Single-end or paired-end reads"*: `Paired-end (as individual datasets)`
>        - {% icon param-file %} *"RNA-Seq FASTQ/FASTA file, forward reads"*: `output` (Input dataset)
>        - {% icon param-file %} *"RNA-Seq FASTQ/FASTA file, reverse reads"*: `output` (Input dataset)
>    - *"Custom or built-in reference genome"*: `Use reference genome from history and create temporary index`
>        - {% icon param-file %} *"Select a reference genome"*: `output` (Input dataset)
>        - *"Length of the SA pre-indexing string"*: `11`
>        - *"Build index with or without known splice junctions annotation"*: `build index without gene-model`
>    - *"Use 2-pass mapping for more sensitive novel splice junction discovery"*: `No`
>    - *"Per gene/transcript output"*: `No per gene or transcript output`
>    - In *"Output filter criteria"*:
>        - *"Exclude the following records from the BAM output"*: ``
>        - *"Would you like to set additional output filters?"*: `No`
>    - In *"Algorithmic settings"*:
>        - *"Configure seed, alignment and limits options"*: `Use Defaults`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

# Strutural annotation

We can now run **Funannotate predict annotation** to perform the structural annotation of the genome.

We need to input the genome sequence, the RNASeq data, and the proteins to align on the genome. We also specify the name of the species and strain (they can be used later for submission to NCBI). Finally, as

There are other parameters to finely tune how Funannotate will run ab-initio predictors to predict genes, and to filter the final results based on various criteria. As Funannotate uses [BUSCO](http://busco.ezlab.org/) (Benchmarking Universal Single-Copy Orthologs) for initial training of ab-initio predictors, we select a dataset close to the species we are annotating: `rhizopus_oryzae`.

Funannotate is also able to use GeneMark to predict new genes, but to due to licensing restrictions, this software is not available on every Galaxy instance. We will ignore this for this tutorial, it will not impact the results too much.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Funannotate predict annotation](toolshed.g2.bx.psu.edu/repos/iuc/funannotate_predict/funannotate_predict/1.8.9+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Assembly to annotate"*: `output` (output of **Sort assembly** {% icon tool %})
>    - *"Funannotate database"*: select the latest version available
>    - In *"Organism"*:
>        - *"Name of the species to annotate"*: `Mucor`
>        - *"Strain name"*: `muc1`
>        - *"Is it a fungus species?"*: `Yes`
>    - In *"Evidences"*:
>        - {% icon param-file %} *"RNA-seq mapped to genome to train Augustus/GeneMark-ET"*: `mapped_reads` (output of **RNA STAR** {% icon tool %})
>        - *"Select protein evidences"*: `Custom protein sequences`
>            - {% icon param-file %} *"Proteins to map to genome"*: `SwissProt_subset.fasta` (Input dataset)
>    - In *"Busco"*:
>        - *"Initial Augustus species training set for BUSCO alignment"*: `rhizopus_oryzae`
>    - *"Which outputs should be generated"*: Select all
>
>    > ### {% icon comment %} Comment
>    >
>    > For *"Select protein evidences"* we select `Custom protein sequences` to reduce the computing time, but for real data analyis, you should select the default value: `Use UniProtKb/SwissProt (from selected Funannotate database)`.
>    {: .comment}
>
{: .hands_on}

This tool produces several output dataset, in particular:

- the full structural annotation in Genbank, GFF3 or NCBI tbl formats: these files contain the position of all the genes that were found on the analysed genome.
- the  CDS, transcript and protein sequences of all the genes predicted by Funannotate
- some statistics and reports

TODO explain validation report + question on this?

This step will take a bit of time to run. While it runs, we can already schedule the following functional annotation steps. Galaxy will run them automatically as soon as the structural annotation is ready.

# Functional annotation

The aim of the previous step is to predict the position of the genes on the genome (structural annotation). Now we want to assign names and functions to the predicted genes. We can do this automatically using specialised tools: **EggNOG Mapper** and **InterProScan**.

## **EggNOG Mapper**

**EggNOG Mapper** compares each protein sequence of the annotation to a huge set of ortholog groups from the [EggNOG database](http://eggnog5.embl.de). In this database, each ortholog group is associated with functional annotation like [Gene Ontology (GO)](http://www.geneontology.org/) terms or [KEGG pathways](https://www.genome.jp/kegg/pathway.html). When the protein sequence of a new gene is found to be very similar to one of these ortholog groups, the corresponding functional annotation is transfered to this new gene.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [eggNOG Mapper](toolshed.g2.bx.psu.edu/repos/galaxyp/eggnog_mapper/eggnog_mapper/2.0.1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Fasta sequences to annotate"*: `fasta_proteins` (output of **Funannotate predict annotation** {% icon tool %})
>    - *"Version of eggNOG Database"*: select the latest version available
>    - In *"Output Options"*:
>        - *"Exclude header lines and stats from output files"*: `No`
>
{: .hands_on}

The output of this tool is a tabular file, where each line represents a gene from our annotation, with the functional annotation that was found by EggNOG-mapper. It includes a predicted protein name, GO terms, EC numbers, KEGG identifiers, ...

## **InterProScan**

[InterPro](https://www.ebi.ac.uk/interpro/) is a huge integrated database of protein families. Each family is characterized by one or muliple signatures (i.e. sequence motifs) that are specific to the protein family, and corresponding functional annotation like protein names or [Gene Ontology (GO)](http://www.geneontology.org/). A good proportion of the signatures are manually curated, which means they are of very good quality.

**InterProScan** is a tool that analyses each protein sequence from our annotation to determine if they contain one or several of the signatures from InterPro. When a protein contains a known signature, the corresponding functional annotation will be assigned to it by **InterProScan**.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [InterProScan](toolshed.g2.bx.psu.edu/repos/bgruening/interproscan/interproscan/5.52-86.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Protein FASTA File"*: `fasta_proteins` (output of **Funannotate predict annotation** {% icon tool %})
>    - *"InterProScan database"*: select the latest version available
>    - *"Output format"*: `XML`
>    - *"Use applications with restricted license, only for non-commercial use?"*: `Yes` (set it to `No` if you use InterProScan for commercial use)
>
{: .hands_on}

The output of this tool is both a tabular file and an XML file. Both contain the same information, but the tabular one is more readable for a Human: each line represents a gene from our annotation, with the different domains and motifs that ere found by InterProScan.

## Integrating the results

Now we have a structural annotation, and the results of both **EggNOG Mapper** and **InterProScan**. Each one is in a separate file, we will now combine all this data into a single file that will contain the structural *and* the functional annotation. This will be the final output of our annotation pipeline, ready to be submitted to the NCBI reference database.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Funannotate functional](toolshed.g2.bx.psu.edu/repos/iuc/funannotate_annotate/funannotate_annotate/1.8.9+galaxy1) %} with the following parameters:
>    - *"Input format"*: `GenBank (from 'Funannotate predict annotation' tool)`
>        - {% icon param-file %} *"Genome annotation in genbank format"*: `annot_gbk` (output of **Funannotate predict annotation** {% icon tool %})
>    - {% icon param-file %} *"NCBI submission template file"*: `output` (Input dataset)
>    - {% icon param-file %} *"Eggnog-mapper annotations file"*: `annotations` (output of **eggNOG Mapper** {% icon tool %})
>    - {% icon param-file %} *"InterProScan5 XML file"*: `outfile_xml` (output of **InterProScan** {% icon tool %})
>    - *"Strain name"*: `muc1`
>    - *"Which outputs should be generated"*: ``
>
{: .hands_on}

# Evaluation and visualisation

TODO look at the "stats" output of funannotate predict

## Sub-step with **Busco**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.2.2+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `fa_proteins` (output of **Funannotate functional** {% icon tool %})
>    - *"Mode"*: `annotated gene sets (protein)`
>    - In *"Advanced Options"*:
>        - *"Which outputs should be generated"*: ``
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

## Sub-step with **JBrowse**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [JBrowse](toolshed.g2.bx.psu.edu/repos/iuc/jbrowse/jbrowse/1.16.11+galaxy1) %} with the following parameters:
>    - *"Reference genome to display"*: `Use a genome from history`
>        - {% icon param-file %} *"Select the reference genome"*: `output` (Input dataset)
>    - *"JBrowse-in-Galaxy Action"*: `New JBrowse Instance`
>    - In *"Track Group"*:
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - *"Track Category"*: `Annotation`
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `GFF/GFF3/BED Features`
>                        - {% icon param-file %} *"GFF/GFF3/BED Track Data"*: `gff3` (output of **Funannotate functional** {% icon tool %})
>                        - *"This is match/match_part data"*: `Yes`
>                        - *"JBrowse Track Type [Advanced]"*: `Neat HTML Features`
>                        - In *"JBrowse Feature Score Scaling & Coloring Options [Advanced]"*:
>                            - *"Color Score Algorithm"*: `Ignore score`
>                                - *"Color Selection"*: `Automatically selected`
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - *"Track Category"*: `RNASeq`
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `BAM Pileups`
>                        - {% icon param-file %} *"BAM Track Data"*: `mapped_reads` (output of **RNA STAR** {% icon tool %})
>                        - *"Autogenerate SNP Track"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

# Comparing annotations

# Submission to NCBI


# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
