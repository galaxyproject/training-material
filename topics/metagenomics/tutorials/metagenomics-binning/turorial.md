---
layout: tutorial_hands_on

title: "Binning of metagenomic sequencing data"
zenodo_link: "https://zenodo.org/record/7818827"
zenodo_link_results: "https://zenodo.org/record/7845138" 
level: Introductory
questions:
  - "What is metagenomic binning refers to?"
  - "Which tools should be used for metagenomic binning?"
  - "How to assess the quality of metagenomic data binning?"
# objectives:
#   - "Describe what an assembly is"
#   - "Describe what de-replication is"
#   - "Explain the difference between co-assembly and individual assembly"
#   - "Explain the difference between reads, contigs and scaffolds"
#   - "Explain how tools based on De Bruijn graph work"
#   - "Apply appropriate tools for analyzing the quality of metagenomic data"
#   - "Construct and apply simple assembly pipelines on short read data"
#   - "Apply appropriate tools for analyzing the quality of metagenomic assembly"
#   - "Evaluate the Quality of the Assembly with Quast, Bowtie2, and CoverM-Genome"
# time_estimation: "2H"
# key_points:
#   - "Assembly groups reads into contigs and scafolds."
#   - "De Brujin Graphs use k-mers to assembly reads"
#   - "MetaSPAdes and MEGAHIT are assemblers"
#   - "Quast is the tool to assess the assembly quality"
contributions:
  authorship:
    - npechl
    - fpsom
tags:
  - binning
  - metagenomics
---

# Introduction

Metagenomics is the study of genetic material recovered directly from environmental samples, such as soil, water, or gut contents, without the need for isolation or cultivation of individual organisms. Metagenomics binning is a process used to classify DNA sequences obtained from metagenomic sequencing into discrete groups, or bins, based on their similarity to each other.

The goal of metagenomics binning is to assign the DNA sequences to the organisms or taxonomic groups that they originate from, allowing for a better understanding of the diversity and functions of the microbial communities present in the sample. This is typically achieved through computational methods that use sequence similarity, composition, and other features to group the sequences into bins.

There are two main types of metagenomics binning: *reference-based* and *de novo*. *Reference-based binning* involves aligning the sequences to a database of known genomes or reference sequences, while *de novo binning* involves clustering the sequences based on similarity without prior knowledge of the organisms or reference sequences present in the sample.

Both methods have their strengths and limitations, and researchers often use a combination of approaches to improve the accuracy of their binning results. Metagenomics binning is an important tool for understanding the functional potential of microbial communities in various environments and has applications in fields such as biotechnology, environmental science, and human health.


**Metagenomics binning is a complex process that involves many steps and can be challenging due to several problems that can occur during the process**. Some of the most common problems encountered in metagenomics binning include:

- **High complexity**: Metagenomic samples contain DNA from multiple organisms, which can lead to high complexity in the data.
- **Fragmented sequences**: Metagenomic sequencing often generates fragmented sequences, which can make it difficult to assign reads to the correct bin.
- **Uneven coverage**: Some organisms in a metagenomic sample may be more abundant than others, leading to uneven coverage of different genomes.
- **Incomplete or partial genomes**: Metagenomic sequencing may not capture the entire genome of a given organism, which can make it difficult to accurately bin sequences from that organism.
- **Horizontal gene transfer**: Horizontal gene transfer (HGT) can complicate metagenomic binning, as it can introduce genetic material from one organism into another.
- **Chimeric sequences**: Sequences that are the result of sequencing errors or contamination can lead to chimeric sequences, which can make it difficult to accurately bin reads.
- **Strain variation**: Organisms within a species can exhibit significant genetic variation, which can make it difficult to distinguish between different strains in a metagenomic sample.


There are plenty of computational tools to perform metafenomics binning. Some of the most widely used include:

- **MaxBin**: A popular de novo binning algorithm that uses a combination of sequence features and marker genes to cluster contigs into genome bins.
- **MetaBAT**: Another widely used de novo binning algorithm that employs a hierarchical clustering approach based on tetranucleotide frequency and coverage information.
- **CONCOCT**: A de novo binning tool that uses a clustering algorithm based on sequence composition and coverage information to group contigs into genome bins.
- **MyCC**: A reference-based binning tool that uses sequence alignment to identify contigs belonging to the same genome or taxonomic group.
- **GroopM**: A hybrid binning tool that combines reference-based and de novo approaches to achieve high binning accuracy.
- **MetaWRAP**: A comprehensive metagenomic analysis pipeline that includes various modules for quality control, assembly, binning, and annotation.
- **Anvi'o**: A platform for visualizing and analyzing metagenomic data, including features for binning, annotation, and comparative genomics.

**In this tutorial, we will learn how to run metagenomic binning tools and evaluate the quality of the results**. In order to do that, we will use data from the study: [Temporal shotgun metagenomic dissection of the coffee fermentation ecosystem](https://www.ebi.ac.uk/metagenomics/studies/MGYS00005630#overview) and MetaBAT2 algorithm. For an in-depth analysis of the structure and functions of the coffee microbiome, a temporal shotgun metagenomic study (six time points) was performed. The six samples have been sequenced with Illumina MiSeq utilizing whole genome sequencing.

Based on the 6 original dataset of the coffee fermentation system, we generated mock datasets for this tutorial.

<!-- > <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda} -->

# Prepare analysis history and data

To run binning, we first need to get the data into Galaxy. Any analysis should get its own Galaxy history. So let's start by creating a new one:

> <hands-on-title>Prepare the Galaxy history</hands-on-title>
>
> 1. Create a new history for this analysis
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Rename the history
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
{: .hands_on}

We need to get the data into our history.

In case of a not very large dataset it's more convenient to upload data directly from your computer to Galaxy.

> <hands-on-title>Upload data into Galaxy</hands-on-title>
>
> 2. Import the sequence read data (\*.fasta) from [Zenodo]({{ page.zenodo_link }}) or a data library:
>
>    ```text
>    {{ page.zenodo_link }}/files/Assembly_with_MEGAHIT_on_ERR2231567.fasta
>    {{ page.zenodo_link }}/files/Assembly_with_MEGAHIT_on_ERR2231568.fasta
>    {{ page.zenodo_link }}/files/Assembly_with_MEGAHIT_on_ERR2231569.fasta
>    {{ page.zenodo_link }}/files/Assembly_with_MEGAHIT_on_ERR2231570.fasta
>    {{ page.zenodo_link }}/files/Assembly_with_MEGAHIT_on_ERR2231571.fasta
>    {{ page.zenodo_link }}/files/Assembly_with_MEGAHIT_on_ERR2231572.fasta
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
>    > <comment-title></comment-title>
>    > In case of large dataset, we can use FTP server or the [Galaxy Rule-based Uploader]({% link topics/galaxy-interface/tutorials/upload-rules/tutorial.md %}).
>    {: .comment}
>
> 3. Create a paired collection named `Raw reads`, rename your pairs with the sample name
>
>    {% snippet faqs/galaxy/collections_build_list_paired.md %}
>
{: .hands_on}

# Binning

As explained before, there are many challenges to metagenomics binning. The most common of them are listed below:

- High complexity.
- Fragmented sequences.
- Uneven coverage.
- Incomplete or partial genomes.
- Horizontal gene transfer.
- Chimeric sequences.
- Strain variation.

![Image show the binning process where sequences are grouped together based on genome signatures like the kmer profiles of each contig, contig coverage, or GC content](./images/binning.png "Binning"){:width="60%"}

In this tutorial we will learn how to use **MetaBAT2** tool through Galaxy:

- **MetaBAT2**: Metagenome Binning based on Abundance and Tetranucleotide frequency

  *Grouping large fragments assembled from shotgun metagenomic sequences to deconvolute complex microbial communities, or metagenome binning, enables the study of individual organisms and their interactions. Here we developed automated metagenome binning software, called MetaBAT, which integrates empirical probabilistic distances of genome abundance and tetranucleotide frequency. On synthetic datasets MetaBAT on average achieves 98percent precision and 90% recall at the strain level with 281 near complete unique genomes. Applying MetaBAT to a human gut microbiome data set we recovered 176 genome bins with 92% precision and 80% recall. Further analyses suggest MetaBAT is able to recover genome fragments missed in reference genomes up to 19%, while 53 genome bins are novel. In summary, we believe MetaBAT is a powerful tool to facilitate comprehensive understanding of complex microbial communities.*

We will use the uploaded assembled fasta files as input to the algorithm (For simplicity reasons all other parameters will be preserved with their default values).

> <hands-on-title>Individual assembly of short-reads with MEGAHIT</hands-on-title>
> 1.  {% tool [ MetaBAT2](toolshed.g2.bx.psu.edu/repos/iuc/megahit/megahit/1.2.9+galaxy0) %} with parameters:
>     - *"Fasta file containing contigs"*: `assembly fasta files`
>     - In *Advanced options*
>       - *"Percentage of good contigs considered for binning decided by connection among contigs"*: `default 95`
>       - *"Minimum score of an edge for binning"*: `default 60`
>       - *"Maximum number of edges per node"*: `default 200`
>       - *"TNF probability cutoff for building TNF graph"*: `default 0`
>       - *"Turn off additional binning for lost or small contigs?"*: `default "No"`
>       - *"Minimum mean coverage of a contig in each library for binning"*: `default 1`
>       - *"Minimum total effective mean coverage of a contig for binning "*: `default 1`
>       - *"For exact reproducibility"*: `default 0`
>     - In *Output options*
>       - *"Minimum size of a bin as the output"*: `default 200000`
>       - *"Output only sequence labels as a list in a column without sequences?"*: `default "No"`
>       - *"Save cluster memberships as a matrix format?"*: `"Yes"`
>
{: .hands_on}

The output files generated by MetaBAT2 include (some of the files below are optional and not produced unless it is required by the user):

1. The final set of genome bins in FASTA format (`.fa`)
2. A summary file with information on each genome bin, including its length, completeness, contamination, and taxonomy classification (`.txt`)
3. A file with the mapping results showing how each contig was assigned to a genome bin (`.bam`)
4. A file containing the abundance estimation of each genome bin (`.txt`)
5. A file with the coverage profile of each genome bin (`.txt`)
6. A file containing the nucleotide composition of each genome bin (`.txt`)
7. A file with the predicted gene sequences of each genome bin (`.faa`)

These output files can be further analyzed and used for downstream applications such as functional annotation, comparative genomics, and phylogenetic analysis.

> <comment-title></comment-title>
>
> Since the binning process would take some we are just going to import the results of the binning previously run.
>
> > <hands-on-title>Import generated assembly files</hands-on-title>
> >
> > 1. Import the six folders containg binning result files from [Zenodo]({{ page.zenodo_link_results }}) or the Shared Data library:
> >
> >    ```text
> >    {{ page.zenodo_link_results }}/files/26_%20MetaBAT2%20on%20data%20ERR2231567_%20Bins.zip
> >    {{ page.zenodo_link_results }}/files/38_%20MetaBAT2%20on%20data%20ERR2231568_%20Bins.zip
> >    {{ page.zenodo_link_results }}/files/47_%20MetaBAT2%20on%20data%20ERR2231569_%20Bins.zip
> >    {{ page.zenodo_link_results }}/files/57_%20MetaBAT2%20on%20data%20ERR2231570_%20Bins.zip
> >    {{ page.zenodo_link_results }}/files/65_%20MetaBAT2%20on%20data%20ERR2231571_%20Bins.zip
> >    {{ page.zenodo_link_results }}/files/74_%20MetaBAT2%20on%20data%20ERR2231572_%20Bins.zip
> >    ```
> >
> >
> > 2. Create a collection named `MEGAHIT Contig`, rename your pairs with the sample name
> >
> {: .hands_on}
{: .comment}

> <question-title></question-title>
>
> 1. How many bins has been for ERR2231567 sample?
> 2. How many sequences are contained in the second bin?
>
> > <solution-title></solution-title>
> >
> > 1. There are 6 bins identified
> > 2. 167 sequences are classified into the second bin.
> >
> {: .solution}
>
{: .question}

<!-- > <details-title>Co-assembly with MetaSPAdes</details-title>
>
> > <hands-on-title>Assembly with MetaSPAdes</hands-on-title>
> >
> > 1. {% tool [MetaSPAdes](toolshed.g2.bx.psu.edu/repos/nml/metaspades/metaspades/3.15.4+galaxy2) %} with following parameters
> >    - *"Pair-end reads input format"*: `Paired-end: list of dataset pairs`
> >        - {% icon param-collection %} *"FASTQ file(s): collection"*: `Raw reads`
> >     - *"Select k-mer detection option"*: `User specific`
> >        - *"K-mer size values"*: `21,33,55,77`
> {: .hands_on}
{: .details} -->

# Checking the quality of the bins

Once binning is done, it is important to check its quality.

Binning results can be evaluated with **CheckM**. CheckM is a software tool used in metagenomics binning to assess the completeness and contamination of genome bins. Metagenomics binning is the process of separating DNA fragments from a mixed community of microorganisms into individual bins, each representing a distinct genome.

CheckM compares the genome bins to a set of universal single-copy marker genes that are present in nearly all bacterial and archaeal genomes. By identifying the presence or absence of these marker genes in the bins, CheckM can estimate the completeness of each genome bin (i.e., the percentage of the total set of universal single-copy marker genes that are present in the bin) and the degree of contamination (i.e., the percentage of marker genes that are found in more than one bin).

This information can be used to evaluate the quality of genome bins and to select high-quality bins for further analysis, such as genome annotation and comparative genomics. CheckM is widely used in metagenomics research and has been shown to be an effective tool for assessing the quality of genome bins. Some of the key functionalities of CheckM are:

- *Estimation of genome completeness*: CheckM uses a set of universal single-copy marker genes to estimate the completeness of genome bins. The completeness score indicates the proportion of these marker genes that are present in the bin, providing an estimate of how much of the genome has been recovered.

- *Estimation of genome contamination*: CheckM uses the same set of marker genes to estimate the degree of contamination in genome bins. The contamination score indicates the proportion of marker genes that are present in multiple bins, suggesting that the genome bin may contain DNA from more than one organism.

- *Identification of potential misassemblies*: CheckM can identify potential misassemblies in genome bins based on the distribution of marker genes across the genome.

- *Visualization of results*: CheckM can generate various plots and tables to visualize the completeness, contamination, and other quality metrics for genome bins, making it easier to interpret the results.

- *Taxonomic classification*: CheckM can also be used to classify genome bins taxonomically based on the presence of specific marker genes associated with different taxonomic groups.

Based on the previous analysis we will use **CheckM lineage_wf**: *Assessing the completeness and contamination of genome bins using lineage-specific marker sets*

`CheckM lineage_wf` is a specific workflow within the CheckM software tool that is used for taxonomic classification of genome bins based on their marker gene content. This workflow uses a reference database of marker genes and taxonomic information to classify the genome bins at different taxonomic levels, from domain to species.

<!-- > <hands-on-title>Evaluation assembly quality with metaQUAST</hands-on-title>
>
> 1. {% tool [Quast](toolshed.g2.bx.psu.edu/repos/iuc/quast/quast/5.2.0+galaxy0) %} with parameters:
>    - *"Data structure for bins"*: `In collection`
>    - *"Bins "*: `Bins produced by MetaBAT2`
>
> 2. Inspect produced table
{: .hands_on} -->

> <comment-title></comment-title>
>
> Since the CheckM process would take some time we are just going to import the results:
>
> > <hands-on-title>Import generated `CheckM lineage_wf` results</hands-on-title>
> >
> > 1. Import the `CheckM lineage_wf` report files from [Zenodo]({{ page.zenodo_link }}) or the Shared Data library:
> >
> >    ```text
> >    {{ page.zenodo_link_results }}/files/CheckM_lineage_wf_on_data_ERR2231567__Bin_statistics.txt
> >    {{ page.zenodo_link_results }}/files/CheckM_lineage_wf_on_data_ERR2231568__Bin_statistics.txt
> >    {{ page.zenodo_link_results }}/files/CheckM_lineage_wf_on_data_ERR2231569__Bin_statistics.txt
> >    {{ page.zenodo_link_results }}/files/CheckM_lineage_wf_on_data_ERR2231570__Bin_statistics.txt
> >    {{ page.zenodo_link_results }}/files/CheckM_lineage_wf_on_data_ERR2231571__Bin_statistics.txt
> >    {{ page.zenodo_link_results }}/files/CheckM_lineage_wf_on_data_ERR2231572__Bin_statistics.txt
> >    ```
> >
> {: .hands_on}
{: .comment}

The output of "CheckM lineage_wf" includes several (optional) files and tables that provide information about the taxonomic classification and quality assessment of genome bins. Here are some of the key outputs:

- "checkm_taxonomy.tsv": This is a tab-separated file that lists the taxonomic assignments for each genome bin. The file includes columns for the bin ID, the domain, phylum, class, order, family, genus, and species (if available) of the closest reference genome(s) to the bin.

- "checkm_taxonomy.log": This file provides detailed information about the taxonomic classification process, including the marker genes used, the reference database, and the classification algorithm.

- "checkm_ms.txt": This file contains a summary of the completeness and contamination scores for each genome bin, along with other quality metrics.

- "checkm_qa.tsv": This file provides more detailed information about the completeness and contamination scores, including the number of marker genes used, the number of genome fragments included in the bin, and the size of the bin.

- "checkm_tree.nwk": This is a Newick format tree file that shows the phylogenetic relationship between the genome bins and the reference genomes used for classification.

<!-- # Visualization of the *de novo* assembly graph

Current metagenome assemblers like MEGAHIT and MetaSPAdes use **graphs**, most typically a de Bruijn graph to stich reads together. In an ideal case, the graph would contain one distinct path for each genome of each micro-organisms, but complexities such as repeated sequences usually prevent this.

Assembly graphs contain then **branching structures**: one node may lead into multiple others. **Contigs** correspond to the longest sequences in the graph that can be determined unambiguously. They are the final results of most assembler. But the assembly graph contains more information. It can be useful for finding sections of the graph, such as rRNA, or to try to find parts of a genome.

**Bandage** ({% cite wick2015bandage %}) is a tool creating interactive visualisations of assembly graphs.

> <hands-on-title>Visualization the assembly graph</hands-on-title>
>
> 1. {% tool [megahit contig2fastg](toolshed.g2.bx.psu.edu/repos/iuc/megahit_contig2fastg/megahit_contig2fastg/1.1.3+galaxy10) %} with parameters:
>    - {% icon param-collection %} *"Contig file"*: Output of **MEGAHIT**
>    - *"K-mer length"*: `91`
>
>      > <comment-title></comment-title>
>      > To get the value, you need to
>      > 1. Go into the **MEGAHIT** output collection
>      > 2. Expand one of the contig file by clicking on it in the history
>      > 3. Check in the dataset peek the name of the contig
>      > 4. Extract the value after the first `k` in the contig names
>      {: .comment}
>
> 2. {% tool [Bandage Image](toolshed.g2.bx.psu.edu/repos/iuc/bandage/bandage_image/0.8.1+galaxy4) %} with parameters:
>    - {% icon param-collection %} *"Graphical Fragment Assembly"*: Output of **megahit contig2fastg**
>
> 3. Inspect the generated image for ERR2231571
{: .hands_on}

![Image shows the assembly graphs, with longer stretch on the top and many small contigs on the bottom](./images/ERR2231571_graph.jpg "Assembly graph for ERR2231571 sample")

The graph is quite disconnected. On the top, we can see the longer stretches, that includes multiples contigs (each contig having a different color). On the bottom are the shortest stretches or single contigs.

But it is really hard to read or extract any information from the graph. Let's inspect the information about the assembly graph

> <hands-on-title>Visualization the *de novo* assembly graph</hands-on-title>
>
> 1. {% tool [Bandage Info](toolshed.g2.bx.psu.edu/repos/iuc/bandage/bandage_info/0.8.1+galaxy2) %} with parameters:
>    - {% icon param-collection %} *"Graphical Fragment Assembly"*: Output of **megahit contig2fastg**
>
> 2. {% tool [Column join](toolshed.g2.bx.psu.edu/repos/iuc/collection_column_join/collection_column_join/0.0.3) %} with parameters:
>    - {% icon param-collection %} *"Tabular files"*: Output of **Bandage Info**
>
> 3. Inspect the generated output
{: .hands_on}

> <question-title></question-title>
>
> 1. How many nodes are in the graph for ERR2231568? And for ERR2231572? What does they correspond to?
> 2. How many edges are in the graph for ERR2231568? And for ERR2231572? What is the impact of these numbers in relation to the number of nodes on the graph?
> 3. How many connected components are there for ERR2231568? And for ERR2231572? What does they correspond to?
> 4. What is the percentage of dead ends are there for ERR2231568? And for ERR2231572?
> 5. What are the smallest and larges edge overlaps?
> 6. What is the largest component? For which sample?
> 7. What is the shortest node? What does they correspond to?
>
> > <solution-title></solution-title>
> >
> > 1. There are 228,719 nodes for ERR2231568 and 122,526 for ERR2231572. They correspond to the number of contigs
> > 2. There are 16,580 edges for ERR2231568 and 13,993 for ERR2231572. There are less edges than nodes in the graph. It means that many nodes/contigs are disconnected
> > 3. There are 212,598 connected components, i.e. number of regions of the graph which are disconnected from each other, for ERR2231568 and 109,044 for ERR2231572
> > 4. There are 94.0702% dead ends, i.e. the end of a node not connected to any other nodes, for ERR2231568 and 90.7032% for ERR2231572. It confirms the previous observation
> > 5. The smallest and larges edge overlaps are 91bp, i.e. the k-mer length
> > 6. The largest component is 340,003 bp for ERR2231567
> > 7. The shortest node is 200 bp, i.e. the minimal size for a contig
> >
> {: .solution}
>
{: .question} -->

<!--# De-replication

De-replication is the process of identifying sets of genomes that are the "same" in a list of genomes, and removing all but the “best” genome from each redundant set. How similar genomes need to be to be considered “same”, how to determine which genome is “best”, and other important decisions are discussed in [Important Concepts](https://drep.readthedocs.io/en/latest/choosing_parameters.html).

A common use for genome de-replication is the case of individual assembly of metagenomic data. If metagenomic samples are collected in a series, a common way to assemble the short reads is with a “co-assembly”. That is, combining the reads from all samples and assembling them together. The problem with this is assembling similar strains together can severely fragment assemblies, precluding recovery of a good genome bin. An alternative option is to assemble each sample separately, and then “de-replicate” the bins from each assembly to make a final genome set.

![Image shows the process of individual assembly on two strains and five samples, after individual assembly of samples two samples are chosen for de-replication process. In parallel, co-assembly on all five samples is performed](./images/individual-assembly.png "Individual assembly followed by de-replication vs co-assembly"){:width="80%"}

If we have different samples, then we do an individual assembly for each sample. In the figure above we see that after individual assembly we have results for every individual assembly represented with pie charts. Different colours on these charts show different strains (organisms). Every chart has a different percentage of every strain which means that the assemblies contain different strains in different proportions in each sample.

Afterwards, we do the process of de-replication. We try to combine all the assemblies and try to identify which genomes are the most proper.

Individual assembly is a good practice as well as co-assembly. They both have pros and cons and that are just different techniques.

Co-assembly is a more common practice. But in case of co-assembly the genome might be more fragmented afterwards (like it is shown in the figure)  and sometimes it can be less proper. However, it should be decided in every single case which approach to use (co- or individual). More comprehensive information about de-replication you can learn from paper {%cite evans2020%} to get more intuition about how de-replication works.

> <hands-on-title>General list of actions for de-replication</hands-on-title>
> 1. Create new history
> 2. Assemble each sample separately using your favorite assembler
> 3. Perform a co-assembly to catch low-abundance microbes
> 4. Bin each assembly separately using your favorite binner
> 5. Bin co-assembly using your favorite binner
> 6. Pull the bins from all assemblies together
> 7. rRun **dRep** on them
> 8. Perform downstream analysis on the de-replicated genome list
>
{: .hands_on}

We will perform steps from 1 to 3 in this tutorial a bit later while steps 4 - 8 will be considered in the following tutorial - Binning tutorial.-->

# Conclusions

In summary, this tutorial shows a step-by-step on how to bin metagenomic contigs using MetaBAT2.

It is critical to select the appropriate binning tool for a specific metagenomics study, as different binning methods may have different strengths and limitations depending on the type of metagenomic data being analyzed. By comparing the outcomes of several binning techniques, researchers can increase the precision and accuracy of genome binning.

There are various binning methods available for metagenomic data, including reference-based, clustering-based, and hybrid approaches. Each method has its advantages and disadvantages, and the selection of the appropriate method depends on the research question and the characteristics of the data.

Comparing the outcomes of multiple binning methods can help to identify the most accurate and reliable method for a specific study. This can be done by evaluating the quality of the resulting bins in terms of completeness, contamination, and strain heterogeneity, as well as by comparing the composition and functional profiles of the identified genomes.

Overall, by carefully selecting and comparing binning methods, researchers can improve the quality and reliability of genome bins, which can ultimately lead to a better understanding of the functional and ecological roles of microbial communities in various environments.