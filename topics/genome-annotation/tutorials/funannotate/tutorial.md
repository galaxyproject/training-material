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
  - Perform functional annotation using eggNOG-mapper and InterProScan
  - Evaluate annotation quality with BUSCO
  - View annotations in JBrowse
time_estimation: 8H
key_points:
  - Funannotate allows to perform structural annotatation of a eukaryotic genome.
  - Functional annotation can be performed using eggNOG-mapper and InterProScan.
  - BUSCO and JBrowse allow to inspect the quality of an annotation.
  - Funannotate allows to format an annotation for sumission at NCBI.
contributors:
- abretaud

---


# Introduction
{:.no_toc}


Genome annotation of eukaryotes is a little more complicated than for prokaryotes: eukaryotic genomes are usually larger than prokaryotes, with more genes. The sequences determining the beginning and the end of a gene are generally less conserved than the prokaryotic ones. Many genes also contain introns, and the limits of these introns (acceptor and donor sites) are not highly conserved.

In this tutorial we will use a software tool called Funannotate ({% cite jonathan_m_palmer_2020_4054262 %}) to annotate the genome sequence of a small eukaryote: [*Mucor mucedo*](https://en.wikipedia.org/wiki/Mucor_mucedo) (a fungal plant pathogen).

As explained on [Funannotate's website](https://funannotate.readthedocs.io/), "it was originally written to annotate fungal genomes (small eukaryotes ~ 30 Mb genomes), but has evolved over time to accomodate larger genomes". As other annotation tools like Maker or Braker, it works by aligning as many evidences as possible along the genome sequence, and then reconciliating all these signals to determine probable gene structures.

The evidences can be transcript or protein sequences from the same (or closely related) organism. These sequences can come from public databases (like NR or GenBank) or from your own experimental data (transcriptome assembly from an RNASeq experiment for example). Funannotate is also able to take into account repeated elements.

Funannotate uses ab-initio predictors ([Augustus](http://bioinf.uni-greifswald.de/augustus/), [SNAP](https://github.com/KorfLab/SNAP), [glimmerHMM](https://bio.tools/glimmer-hmm), [CodingQuarry](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1344-4) and [GeneMark-ES/ET](http://exon.gatech.edu/GeneMark/) (optional due to licensing)) to improve its predictions: these software tools are able to make gene structure predictions by analysing only the genome sequence with a statistical model.

While for [Maker]({% link topics/genome-annotation/tutorials/annotation-with-maker/tutorial.md %}) you need to perform training steps for the ab-initio predictors, Funannotate is able to take care of that for you, which makes it much easier to use.

In this tutorial you will learn how to perform a structural genome annotation, and how to evaluate its quality. Then you will learn how to run functional annotation, using eggNOG-mapper and InterProScan to automatically assign names and functions to the annotated genes. And you will also learn how Funannotate can prepare files ready for submission of your annotation to the NCBI.

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
>    https://zenodo.org/api/files/xxxx/assembly_masked.fasta
>    https://zenodo.org/api/files/xxxx/RNASeq_downsampled_0.1.bam
>    https://zenodo.org/api/files/xxxx/SwissProt_subset.fasta
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
{: .hands_on}

# Preparing the genome sequence

Before annotating the genome, we want to make sure that the fasta file is properly formatted. We do it now to make sure we will not encounter unexpected errors later in the annotation process.

Funannotate provides two little tools to help us, let's run them one after the other.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Funannotate assembly clean](toolshed.g2.bx.psu.edu/repos/iuc/funannotate_clean/funannotate_clean/1.8.9+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Assembly to clean"*: `assembly_masked.fasta` (Input dataset)
>
{: .hands_on}


> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Sort assembly](toolshed.g2.bx.psu.edu/repos/iuc/funannotate_sort/funannotate_sort/1.8.9+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Assembly to sort"*: `output` (output of **Funannotate assembly clean** {% icon tool %})
>
{: .hands_on}

## Sub-step with **Busco**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Busco** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `output` (Input dataset)
>    - *"Mode"*: `Genome assemblies (DNA)`
>        - *"Use Augustus instead of Metaeuk"*: `Use Metaeuk`
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

## Sub-step with **Funannotate predict annotation**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Funannotate predict annotation** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Assembly to annotate"*: `output` (output of **Sort assembly** {% icon tool %})
>    - In *"Organism"*:
>        - *"Name of the species to annotate"*: `Mucor`
>        - *"Strain name"*: `muc1`
>        - *"Is it a fungus species?"*: `Yes`
>    - In *"Evidences"*:
>        - {% icon param-file %} *"RNA-seq mapped to genome to train Augustus/GeneMark-ET"*: `output` (Input dataset)
>        - *"Select protein evidences"*: `Custom protein sequences`
>            - {% icon param-file %} *"Proteins to map to genome"*: `output` (Input dataset)
>    - In *"EVM settings (advanced)"*:
>        - *"Split contigs into partitions for EVM processing?"*: `Yes`
>    - *"Which outputs should be generated"*: ``
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

## Sub-step with **eggNOG Mapper**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **eggNOG Mapper** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Fasta sequences to annotate"*: `fasta_proteins` (output of **Funannotate predict annotation** {% icon tool %})
>    - In *"Diamond Options"*:
>        - *"Scoring matrix and gap costs"*: `BLOSUM62`
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

## Sub-step with **Interproscan functional predictions of ORFs**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Interproscan functional predictions of ORFs** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Protein Fasta File"*: `fasta_proteins` (output of **Funannotate predict annotation** {% icon tool %})
>    - *"Output format"*: `XML`
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

## Sub-step with **Funannotate functional**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Funannotate functional** {% icon tool %} with the following parameters:
>    - *"Input format"*: `GenBank (from 'Funannotate predict annotation' tool)`
>        - {% icon param-file %} *"Genome annotation in genbank format"*: `annot_gbk` (output of **Funannotate predict annotation** {% icon tool %})
>    - {% icon param-file %} *"NCBI submission template file"*: `output` (Input dataset)
>    - {% icon param-file %} *"Eggnog-mapper annotations file"*: `annotations` (output of **eggNOG Mapper** {% icon tool %})
>    - {% icon param-file %} *"InterProScan5 XML file"*: `outfile` (output of **Interproscan functional predictions of ORFs** {% icon tool %})
>    - *"Strain name"*: `muc1`
>    - *"Which outputs should be generated"*: ``
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

## Sub-step with **JBrowse**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **JBrowse** {% icon tool %} with the following parameters:
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
>                        - {% icon param-file %} *"BAM Track Data"*: `output` (Input dataset)
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

## Sub-step with **Busco**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Busco** {% icon tool %} with the following parameters:
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

## Sub-step with **Genome annotation statistics**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Genome annotation statistics** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotation to analyse"*: `gff3` (output of **Funannotate functional** {% icon tool %})
>    - *"Reference genome"*: `Use a genome from history`
>        - {% icon param-file %} *"Corresponding genome sequence"*: `output` (Input dataset)
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


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
