---
layout: tutorial_hands_on

title: 'Genome assembly using PacBio data'
zenodo_link: 'https://zenodo.org/record/5702408#.YZUb5uvjIiU'
questions:
- How to perform genome assembly with PacBio data ?
- How to check assembly quality ?
objectives:
- Assemble a Genome
- Assess assembly quality
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- Anthony Bretaudeau
- Alexandre Cormier
- Erwan Corre
- Laura Leroi
- Stéphanie Robin

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->
In this tutorial, we will assemble a Mucor mucedo genome from PacBio sequencing data. These data were obtained from NCBI (SRR8534473 SRR8534474 SRR8534475). The quality of the assembly obtained will be analyzed, in particular by comparing it to a reference assembly.


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}



# Get data
We will use long reads sequencing data (PacBio sequencing) which are a subset data from NCBI

## Get data from Zenodo

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo](https://zenodo.org/record/5702408)
>
>    ```
>   https://zenodo.org/record/5702408/files/SRR8534473_subreads.fastq.gz?download=1
>   https://zenodo.org/record/5702408/files/SRR8534474_subreads.fastq.gz?download=1
>   https://zenodo.org/record/5702408/files/SRR8534475_subreads.fastq.gz?download=1
>    ```
>    
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype are  `fastqsanger.gz`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
>
{: .hands_on}

# Genome Assembly




## Assembly with **Flye**

We will use *Flye*, a de novo assembler for single molecule sequencing reads, such as those produced by PacBio and Oxford Nanopore Technologies. It is designed for a wide range of datasets, from small bacterial projects to large mammalian-scale assemblies. The package represents a complete pipeline: it takes raw PacBio / ONT reads as input and outputs polished contigs. Flye also has a special mode for metagenome assembly. All informations about Flye assembler are here : [Flye](https://github.com/fenderglass/Flye/)

> ### {% icon hands_on %} Hands-on: Assembly
>
> 1. {% tool [Flye](toolshed.g2.bx.psu.edu/repos/bgruening/flye/flye/2.8.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input reads"*: the three sequencing datasets
>    - *"Mode"*: `PacBio raw`
>    - *"Number of polishing iterations"*: `1`
>    - *"Reduced contig assembly coverage"*: `Disable reduced coverage for initial disjointing assembly`
>    
>     The tool produces four dataset outputs : consensus, assembly graph, graphical fragment assembly and assembly info


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

## Assembly metrics with **Fasta Statistics**

*Fasta statistics* displays the summary statistics for a fasta file. In the case of a genome assembly, we need to calculate different metrics such as assembly size, scaffolds number or N50 value. These metrics will allow us to evaluate the quality of this assembly.

> ### {% icon hands_on %} Hands-on: Fasta statistics
>
> 1. {% tool [Fasta Statistics](toolshed.g2.bx.psu.edu/repos/iuc/fasta_stats/fasta-stats/1.0.3) %} with the following parameters:
>    - {% icon param-file %} *"fasta or multifasta file"*: `consensus` (output of **Flye** {% icon tool %})
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

## Assemblies comparison with **Quast**

QUAST = QUality ASsessment Tool is a tool to evaluate genome assemblies by computing various metrics. We will use to compare our genome assembly with a reference genome.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Quast](toolshed.g2.bx.psu.edu/repos/iuc/quast/quast/5.0.2+galaxy1) %} with the following parameters:
>    - *"Use customized names for the input files?"*: `No, use dataset names`
>        - {% icon param-file %} *"Contigs/scaffolds file"*: `consensus` (output of **Flye** {% icon tool %})
>    - *"Type of assembly"*: `Genome`
>        - *"Use a reference genome?"*: `Yes`
>        - {% icon param-file %} *"Reference genome"*: `Mucmuc1_AssemblyScaffolds.fasta"
>        - *"Type of organism"*: `Fungus: use of GeneMark-ES for gene finding, ...`
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

## Genome assembly assessment with **BUSCO**

> ### {% icon hands_on %} Hands-on: BUSCO
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.2.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `consensus` (output of **Flye** {% icon tool %})
>    - *"Mode"*: `Genome assemblies (DNA)`
>        - *"Use Augustus instead of Metaeuk"*: `Use Metaeuk`
>    - *"Auto-detect or select lineage"*: `Select lineage`
>        - *"Lineage"*: `Mucorales`
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
