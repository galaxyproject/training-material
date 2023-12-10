---
layout: tutorial_hands_on

title: "Binning of metagenomic sequencing data"
zenodo_link: "https://zenodo.org/record/7818827"
extra:
  zenodo_link_results: "https://zenodo.org/record/7845138" 
level: Introductory
questions:
  - "What is metagenomic binning refers to?"
  - "Which tools should be used for metagenomic binning?"
  - "How to assess the quality of metagenomic data binning?"
objectives:
  - "Describe what metagenomics binning is"
  - "Describe common problems in metagenomics binning"
  - "What software tools are available for metagenomics binning"
  - "Binning of contigs into metagenome-assembled genomes (MAGs) using MetaBAT 2 software"
  - "Evaluation of MAG quality and completeness using CheckM software"
time_estimation: "2H"
key_points:
  - "Metagenomics binning is a computational approach to grouping together DNA sequences from a mixed microbial sample into metagenome-assembled genomes (MAGs)"
  - "The metagenomics binning workflow involves several steps, including preprocessing of raw sequencing data, assembly of sequencing reads into contigs, binning of contigs into MAGs, quality assessment of MAGs, and annotation of functional genes and metabolic pathways in MAGs"
  - "The quality and completeness of MAGs can be evaluated using standard metrics, such as completeness, contamination, and genome size"
  - "Metagenomics binning can be used to gain insights into the composition, diversity, and functional potential of microbial communities, and can be applied to a range of research areas, such as human health, environmental microbiology, and biotechnology"
contributions:
  authorship:
    - npechl
    - fpsom
tags:
  - binning
  - metagenomics
  - microgalaxy
---

Metagenomics is the study of genetic material recovered directly from environmental samples, such as soil, water, or gut contents, without the need for isolation or cultivation of individual organisms. Metagenomics binning is a process used to classify DNA sequences obtained from metagenomic sequencing into discrete groups, or bins, based on their similarity to each other.

The goal of metagenomics binning is to assign the DNA sequences to the organisms or taxonomic groups that they originate from, allowing for a better understanding of the diversity and functions of the microbial communities present in the sample. This is typically achieved through computational methods that include sequence similarity, composition, and other features to group the sequences into bins.

There are several approaches to metagenomics binning, including:

- **Sequence composition-based binning**: This method is based on the observation that different genomes have distinct sequence composition patterns, such as GC content or codon usage bias. By analyzing these patterns in metagenomic data, sequence fragments can be assigned to individual genomes or groups of genomes.

- **Coverage-based binning**: This method uses the depth of coverage of sequencing reads to group them into bins. Sequencing reads that originate from the same genome are expected to have similar coverage, and this information can be used to identify groups of reads that represent individual genomes or genome clusters.

- **Hybrid binning**: This method combines sequence composition-based and coverage-based binning to increase the accuracy of binning results. By using multiple sources of information, hybrid binning can better distinguish closely related genomes that may have similar sequence composition patterns.

- **Clustering-based binning**: This method groups sequence fragments into clusters based on sequence similarity, and then assigns each cluster to a genome or genome cluster based on its sequence composition and coverage. This method is particularly useful for metagenomic data sets with high levels of genomic diversity.

- **Supervised machine learning-based binning**: This method uses machine learning algorithms trained on annotated reference genomes to classify metagenomic data into bins. This approach can achieve high accuracy but requires a large number of annotated genomes for training.

Each of these methods has its strengths and limitations, and the choice of binning method depends on the specific characteristics of the metagenomic data set and the research question being addressed.


**Metagenomics binning is a complex process that involves many steps and can be challenging due to several problems that can occur during the process**. Some of the most common problems encountered in metagenomics binning include:

- **High complexity**: Metagenomic samples contain DNA from multiple organisms, which can lead to high complexity in the data.
- **Fragmented sequences**: Metagenomic sequencing often generates fragmented sequences, which can make it difficult to assign reads to the correct bin.
- **Uneven coverage**: Some organisms in a metagenomic sample may be more abundant than others, leading to uneven coverage of different genomes.
- **Incomplete or partial genomes**: Metagenomic sequencing may not capture the entire genome of a given organism, which can make it difficult to accurately bin sequences from that organism.
- **Horizontal gene transfer**: Horizontal gene transfer (HGT) can complicate metagenomic binning, as it can introduce genetic material from one organism into another.
- **Chimeric sequences**: Sequences that are the result of sequencing errors or contamination can lead to chimeric sequences, which can make it difficult to accurately bin reads.
- **Strain variation**: Organisms within a species can exhibit significant genetic variation, which can make it difficult to distinguish between different strains in a metagenomic sample.

There are plenty of computational tools to perform metafenomics binning. Some of the most widely used include:

- **MaxBin**  ({%cite maxbin2015%}): A popular de novo binning algorithm that uses a combination of sequence features and marker genes to cluster contigs into genome bins.
- **MetaBAT** ({%cite Kang2019%}): Another widely used de novo binning algorithm that employs a hierarchical clustering approach based on tetranucleotide frequency and coverage information.
- **CONCOCT** ({%cite Alneberg2014%}): A de novo binning tool that uses a clustering algorithm based on sequence composition and coverage information to group contigs into genome bins.
- **MyCC** ({%cite Lin2016%}): A reference-based binning tool that uses sequence alignment to identify contigs belonging to the same genome or taxonomic group.
- **GroopM** ({%cite Imelfort2014%}): A hybrid binning tool that combines reference-based and de novo approaches to achieve high binning accuracy.
- **MetaWRAP** ({%cite Uritskiy2018%}): A comprehensive metagenomic analysis pipeline that includes various modules for quality control, assembly, binning, and annotation.
- **Anvi'o** ({%cite Eren2015%}): A platform for visualizing and analyzing metagenomic data, including features for binning, annotation, and comparative genomics.
- **SemiBin** ({%cite Pan2022%}): A command tool for metagenomic binning with deep learning, handles both short and long reads.

A benchmark study of metagenomics software can be found at {%cite Sczyrba2017%}. MetaBAT 2 outperforms previous MetaBAT and other alternatives in both accuracy and computational efficiency . All are based on default parameters ({%cite Sczyrba2017%}).

**In this tutorial, we will learn how to run metagenomic binning tools and evaluate the quality of the results**. In order to do that, we will use data from the study: [Temporal shotgun metagenomic dissection of the coffee fermentation ecosystem](https://www.ebi.ac.uk/metagenomics/studies/MGYS00005630#overview) and MetaBAT 2 algorithm. MetaBAT is a popular software tool for metagenomics binning, and there are several reasons why it is often used:
- *High accuracy*: MetaBAT uses a combination of tetranucleotide frequency, coverage depth, and read linkage information to bin contigs, which has been shown to be highly accurate and efficient.
- *Easy to use*: MetaBAT has a user-friendly interface and can be run on a standard desktop computer, making it accessible to a wide range of researchers with varying levels of computational expertise.
- *Flexibility*: MetaBAT can be used with a variety of sequencing technologies, including Illumina, PacBio, and Nanopore, and can be applied to both microbial and viral metagenomes.
- *Scalability*: MetaBAT can handle large-scale datasets, and its performance has been shown to improve with increasing sequencing depth.
- *Compatibility*: MetaBAT outputs MAGs in standard formats that can be easily integrated into downstream analyses and tools, such as taxonomic annotation and functional prediction.

For an in-depth analysis of the structure and functions of the coffee microbiome, a temporal shotgun metagenomic study (six time points) was performed. The six samples have been sequenced with Illumina MiSeq utilizing whole genome sequencing.

Based on the 6 original dataset of the coffee fermentation system, we generated mock datasets for this tutorial.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Prepare analysis history and data

MetaBAT 2 takes metagenomic sequencing data as input, typically in the form of assembled contigs in fasta format and coverage information in bam format. Specifically, MetaBAT 2 requires two input files:

- A fasta file containing the assembled contigs, which can be generated from raw metagenomic sequencing reads using an assembler such as MEGAHIT, SPAdes, or IDBA-UD.

- A bam file containing the read coverage information for each contig, which can be generated from the same sequencing reads using mapping software such as Bowtie2 or BWA.

MetaBAT 2 also requires a configuration file specifying various parameters and options for the binning process, such as the minimum contig length, the maximum number of clusters to generate, and the maximum expected contamination level.

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
>    {{ page.zenodo_link }}/files/contigs_ERR2231567.fasta
>    {{ page.zenodo_link }}/files/contigs_ERR2231568.fasta
>    {{ page.zenodo_link }}/files/contigs_ERR2231569.fasta
>    {{ page.zenodo_link }}/files/contigs_ERR2231570.fasta
>    {{ page.zenodo_link }}/files/contigs_ERR2231571.fasta
>    {{ page.zenodo_link }}/files/contigs_ERR2231572.fasta
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
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

![Image show the binning process where sequences are grouped together based on genome signatures like the kmer profiles of each contig, contig coverage, or GC content](./binning.png "Binning"){:width="60%"}

In this tutorial we will learn how to use **MetaBAT 2** tool through Galaxy. **MetaBAT** stands for "Metagenome Binning based on Abundance and Tetranucleotide frequency". It is:

> *Grouping large fragments assembled from shotgun metagenomic sequences to deconvolute complex microbial communities, or metagenome binning, enables the study of individual organisms and their interactions. Here we developed automated metagenome binning software, called MetaBAT, which integrates empirical probabilistic distances of genome abundance and tetranucleotide frequency. On synthetic datasets MetaBAT on average achieves 98percent precision and 90% recall at the strain level with 281 near complete unique genomes. Applying MetaBAT to a human gut microbiome data set we recovered 176 genome bins with 92% precision and 80% recall. Further analyses suggest MetaBAT is able to recover genome fragments missed in reference genomes up to 19%, while 53 genome bins are novel. In summary, we believe MetaBAT is a powerful tool to facilitate comprehensive understanding of complex microbial communities ({%cite Kang2019%}).*

We will use the uploaded assembled fasta files as input to the algorithm (For simplicity reasons all other parameters will be preserved with their default values).

> <hands-on-title>Individual binning of short-reads with MetaBAT 2</hands-on-title>
> 1.  {% tool [MetaBAT 2](toolshed.g2.bx.psu.edu/repos/iuc/megahit/megahit/1.2.9+galaxy0) %} with parameters:
>     - *"Fasta file containing contigs"*: `assembly fasta files`
<!-- >     - In *Advanced options*
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
>       - *"Save cluster memberships as a matrix format?"*: `"Yes"` -->
>
{: .hands_on}

The output files generated by MetaBAT 2 include (some of the files below are optional and not produced unless it is required by the user):

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
> > 1. Import the six folders containg binning result files from [Zenodo]({{ page.extra.zenodo_link_results }}) or the Shared Data library:
> >
> >    ```text
> >    {{ page.extra.zenodo_link_results }}/files/26_%20MetaBAT 2%20on%20data%20ERR2231567_%20Bins.zip
> >    {{ page.extra.zenodo_link_results }}/files/38_%20MetaBAT 2%20on%20data%20ERR2231568_%20Bins.zip
> >    {{ page.extra.zenodo_link_results }}/files/47_%20MetaBAT 2%20on%20data%20ERR2231569_%20Bins.zip
> >    {{ page.extra.zenodo_link_results }}/files/57_%20MetaBAT 2%20on%20data%20ERR2231570_%20Bins.zip
> >    {{ page.extra.zenodo_link_results }}/files/65_%20MetaBAT 2%20on%20data%20ERR2231571_%20Bins.zip
> >    {{ page.extra.zenodo_link_results }}/files/74_%20MetaBAT 2%20on%20data%20ERR2231572_%20Bins.zip
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

# De-replication

De-replication is the process of identifying sets of genomes that are the "same" in a list of genomes, and removing all but the “best” genome from each redundant set. How similar genomes need to be to be considered “same”, how to determine which genome is “best”, and other important decisions are discussed in [Important Concepts](https://drep.readthedocs.io/en/latest/choosing_parameters.html).

A common use for genome de-replication is the case of individual assembly of metagenomic data. If metagenomic samples are collected in a series, a common way to assemble the short reads is with a “co-assembly”. That is, combining the reads from all samples and assembling them together. The problem with this is assembling similar strains together can severely fragment assemblies, precluding recovery of a good genome bin. An alternative option is to assemble each sample separately, and then “de-replicate” the bins from each assembly to make a final genome set.

![Image shows the process of individual assembly on two strains and five samples, after individual assembly of samples two samples are chosen for de-replication process. In parallel, co-assembly on all five samples is performed](./individual-assembly.png "Individual assembly followed by de-replication vs co-assembly"){:width="80%"}

MetaBAT 2 does not explicitly perform dereplication in the sense of identifying groups of identical or highly similar genomes in a given dataset. Instead, MetaBAT 2 focuses on improving the accuracy of binning by leveraging various features such as read coverage, differential coverage across samples, and sequence composition. It aims to distinguish between different genomes present in the metagenomic dataset and assign contigs to the appropriate bins.

Several tools have been designed for the proccess of de-replication. **`dRep`** is a software tool designed for the dereplication of genomes in metagenomic datasets. The goal is to retain a representative set of genomes to improve downstream analyses, such as taxonomic profiling and functional annotation.

An typical workflow of how `dRep` works for dereplication in metagenomics includes:

- *Genome Comparison*: `dRep` uses a pairwise genome comparison approach to assess the similarity between genomes in a given metagenomic dataset.

- *Clustering*: Based on the genome similarities, `dRep` performs clustering to group similar genomes into "genome clusters." Each cluster represents a group of closely related genomes.

- *Genome Quality Assessment*: `dRep` evaluates the quality of each genome within a cluster. It considers factors such as completeness, contamination, and strain heterogeneity.

- *Genome Selection*: Within each genome cluster, `dRep` selects a representative genome based on user-defined criteria. This representative genome is considered as the "dereplicated" version of the cluster.

- *Dereplication Output*: The output of `dRep` includes information about the dereplicated genomes, including their identity, completeness, and contamination. The user can choose a threshold for genome similarity to control the level of dereplication.

> <hands-on-title>General list of actions for de-replication</hands-on-title>
> 1. Create new history
> 2. Assemble each sample separately using your favorite assembler
> 3. Perform a co-assembly to catch low-abundance microbes
> 4. Bin each assembly separately using your favorite binner
> 5. Bin co-assembly using your favorite binner
> 6. Pull the bins from all assemblies together
> 7. rRun **`dRep`** on them
> 8. Perform downstream analysis on the de-replicated genome list
>
{: .hands_on}


# Checking the quality of the bins

Once binning is done, it is important to check its quality.

Binning results can be evaluated with **CheckM** ({%cite Parks2015%}). CheckM is a software tool used in metagenomics binning to assess the completeness and contamination of genome bins. Metagenomics binning is the process of separating DNA fragments from a mixed community of microorganisms into individual bins, each representing a distinct genome.

CheckM compares the genome bins to a set of universal single-copy marker genes that are present in nearly all bacterial and archaeal genomes. By identifying the presence or absence of these marker genes in the bins, CheckM can estimate the completeness of each genome bin (i.e., the percentage of the total set of universal single-copy marker genes that are present in the bin) and the degree of contamination (i.e., the percentage of marker genes that are found in more than one bin).

This information can be used to evaluate the quality of genome bins and to select high-quality bins for further analysis, such as genome annotation and comparative genomics. CheckM is widely used in metagenomics research and has been shown to be an effective tool for assessing the quality of genome bins. Some of the key functionalities of CheckM are:

- *Estimation of genome completeness*: CheckM uses a set of universal single-copy marker genes to estimate the completeness of genome bins. The completeness score indicates the proportion of these marker genes that are present in the bin, providing an estimate of how much of the genome has been recovered.

- *Estimation of genome contamination*: CheckM uses the same set of marker genes to estimate the degree of contamination in genome bins. The contamination score indicates the proportion of marker genes that are present in multiple bins, suggesting that the genome bin may contain DNA from more than one organism.

- *Identification of potential misassemblies*: CheckM can identify potential misassemblies in genome bins based on the distribution of marker genes across the genome.

- *Visualization of results*: CheckM can generate various plots and tables to visualize the completeness, contamination, and other quality metrics for genome bins, making it easier to interpret the results.

- *Taxonomic classification*: CheckM can also be used to classify genome bins taxonomically based on the presence of specific marker genes associated with different taxonomic groups.

Based on the previous analysis we will use **CheckM lineage_wf**: *Assessing the completeness and contamination of genome bins using lineage-specific marker sets*

`CheckM lineage_wf` is a specific workflow within the CheckM software tool that is used for taxonomic classification of genome bins based on their marker gene content. This workflow uses a reference database of marker genes and taxonomic information to classify the genome bins at different taxonomic levels, from domain to species.

> <hands-on-title>Assessing the completeness and contamination of genome bins using lineage-specific marker sets with `CheckM lineage_wf`</hands-on-title>
> 1.  {% tool [CheckM lineage_wf](toolshed.g2.bx.psu.edu/repos/iuc/checkm_lineage_wf/checkm_lineage_wf/1.2.0+galaxy0) %} with parameters:
>     - *"Bins"*: `Folder containing the produced bins`
<!-- >     - In *Advanced options*
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
>       - *"Save cluster memberships as a matrix format?"*: `"Yes"` -->
>
{: .hands_on}

> <comment-title></comment-title>
>
> Since the CheckM process would take some time we are just going to import the results:
>
> > <hands-on-title>Import generated `CheckM lineage_wf` results</hands-on-title>
> >
> > 1. Import the `CheckM lineage_wf` report files from [Zenodo]({{ page.extra.zenodo_link_results }}) or the Shared Data library:
> >
> >    ```text
> >    {{ page.extra.zenodo_link_results }}/files/CheckM_lineage_wf_on_data_ERR2231567__Bin_statistics.txt
> >    {{ page.extra.zenodo_link_results }}/files/CheckM_lineage_wf_on_data_ERR2231568__Bin_statistics.txt
> >    {{ page.extra.zenodo_link_results }}/files/CheckM_lineage_wf_on_data_ERR2231569__Bin_statistics.txt
> >    {{ page.extra.zenodo_link_results }}/files/CheckM_lineage_wf_on_data_ERR2231570__Bin_statistics.txt
> >    {{ page.extra.zenodo_link_results }}/files/CheckM_lineage_wf_on_data_ERR2231571__Bin_statistics.txt
> >    {{ page.extra.zenodo_link_results }}/files/CheckM_lineage_wf_on_data_ERR2231572__Bin_statistics.txt
> >    ```
> >
> {: .hands_on}
{: .comment}

The output of "CheckM lineage_wf" includes several files and tables that provide information about the taxonomic classification and quality assessment of genome bins. Here are some of the key outputs:

- **CheckM Lineage Workflow Output Report**: This report provides a summary of the quality assessment performed by CheckM. It includes statistics such as the number of genomes analyzed, their completeness, contamination, and other quality metrics.

- **Lineage-specific Quality Assessment**: CheckM generates lineage-specific quality assessment files for each analyzed genome. These files contain detailed information about the completeness and contamination of the genome based on its taxonomic lineage.

- **Marker Set Analysis**: CheckM uses a set of marker genes to estimate genome completeness and contamination. The tool produces marker-specific analysis files that provide details on the presence, absence, and copy number of each marker gene in the analyzed genomes.

- **Visualizations**: CheckM generates various visualizations to aid in the interpretation of the results. These include plots such as the lineage-specific completeness and contamination plots, scatter plots, and other visual representations of the data.

- **Tables and Data Files**: CheckM generates tabular data files that contain detailed information about the analyzed genomes, including their names, taxonomic assignments, completeness scores, contamination scores, and other relevant metrics. These files are useful for further downstream analysis or data manipulation.

It should be noted that "CheckM lineage_wf" offers a range of optional outputs that can be generated to provide additional information to the user.

<!-- If we have different samples, then we do an individual assembly for each sample. In the figure above we see that after individual assembly we have results for every individual assembly represented with pie charts. Different colours on these charts show different strains (organisms). Every chart has a different percentage of every strain which means that the assemblies contain different strains in different proportions in each sample.

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

We will perform steps from 1 to 3 in this tutorial a bit later while steps 4 - 8 will be considered in the following tutorial - Binning tutorial. -->

# Conclusions

In summary, this tutorial shows a step-by-step on how to bin metagenomic contigs using MetaBAT 2.

It is critical to select the appropriate binning tool for a specific metagenomics study, as different binning methods may have different strengths and limitations depending on the type of metagenomic data being analyzed. By comparing the outcomes of several binning techniques, researchers can increase the precision and accuracy of genome binning.

There are various binning methods available for metagenomic data, including reference-based, clustering-based, hybrid approaches, and machine learning. Each method has its advantages and disadvantages, and the selection of the appropriate method depends on the research question and the characteristics of the data.

Comparing the outcomes of multiple binning methods can help to identify the most accurate and reliable method for a specific study. This can be done by evaluating the quality of the resulting bins in terms of completeness, contamination, and strain heterogeneity, as well as by comparing the composition and functional profiles of the identified genomes.

Overall, by carefully selecting and comparing binning methods, researchers can improve the quality and reliability of genome bins, which can ultimately lead to a better understanding of the functional and ecological roles of microbial communities in various environments.
