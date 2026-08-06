---
layout: tutorial_hands_on

title: "Pangenome analysis with Roary"
zenodo_link: "https://doi.org/10.5281/zenodo.21339040"
tags:
  - pangenome
  - prokaryote
  - phylogeny
questions:
  - "How can we identify core and accessory genes across multiple bacterial strains?"
  - "How do we generate a pangenome using Prokka and Roary?"
  - "How can we infer the evolutionary relationship of these strains?"
objectives:
  - "Annotate genomes using Prokka."
  - "Use Roary to identify core and accessory genes."
  - "Construct a phylogenetic tree of core genes with IQ-TREE."
time_estimation: "1h"
level: Introductory
key_points:
  - "Prokka is a useful tool to annotate a bacterial genome."
  - "Roary compares annotated genomes to find core and accessory genes."
  - "Core genes can be used to build a phylogenetic tree with IQ-TREE to show relationships."
contributions:
  authorship:
    - Maed0x
    - SaimMomin12
edam_ontology:
- topic_0362 # Genome annotation
- topic_0622 # Genomics

---

# Pangenomics

When the sequencing era truly began in the 1970s with the establishment of Sanger sequencing, it laid the foundation for the Human Genome Project, which in turn led to the human reference genome used widely today.

In traditional bioinformatic workflows, a single, linear reference genome (like the human reference genome) is dominantly used as the standard approach for comparative genomics. While this has driven great advancements in genetics, it suffers from a fundamental limitation known as reference bias ({% cite Ballouz2019 %}). Because a single reference genome only represents the genetic assembly of one individual (or a mosaic of a few), it fails to capture the broader genetic profile of an entire species. Consequently, structural variations or novel sequences present in other strains are often misaligned, unaligned, or completely discarded during mapping ({% cite Matthews2024 %}). Therefore, linear reference genomes are insufficient to represent true genetic diversity.

To overcome this, pangenomes can be used. The field of Pangenomics is an approach that captures the natural genomic variation within a species far better than a single linear reference genome ({% cite Liao2023 %}). Ideally, it represents the entire genomic repertoire of a species. By capturing all of this information, the pangenomic approach aims to eliminate reference bias, ensuring that no genetic information is ignored during mapping simply because it does not exist in a single reference assembly (see Figure 1).

<figure style="text-align: center;">
  <img src="./images/pangenome_comparison.png" alt="Comparison of a pangenome approach with traditional reference approach" style="width:100%;">
  <figcaption>
    <strong>Figure 1: Comparison of a pangenome approach with a traditional reference approach</strong> ({% cite Matthews2024 %}). With a traditional reference approach, any sample reads that cannot be mapped to the reference are discarded. A pangenomic approach solves this limitation by capturing a wider range of genetic diversity. This allows a higher percentage of sample reads to map, meaning more of the actual sample genome can be used for downstream analysis.
  </figcaption>
</figure>

The concept of a pangenome is not novel. However, when the term "Pangenome" was born, it was not financially or computationally feasible to pursue it further. "Pangenome" was first used by {% cite Sigaux2000 %} to describe a database of genomic and transcriptomic alterations found in tumor cells versus normal cells. {% cite Tettelin2005 %} also brought the concept of pangenomics into the spotlight while sequencing multiple strains of Group B Streptococcus (GBS) to help develop a vaccine. During this process, they noticed that with every new strain they sequenced, they continued to identify novel genes absent from previously sequenced strains. They concluded that a single reference genome was insufficient to capture the full genetic diversity of a species. If they had used the first sequenced strain as their sole reference, all the novel genes from the new strains would have been unmapped and discarded during analysis.

With the advent of Next-Generation Sequencing technologies around 2010, sequencing became significantly faster and more affordable. Today, the reduced sequencing costs, combined with increased computational resources compared to the early 2000s, make it possible to properly realize true pangenomes. This technological leap led to the release of the first human pangenome draft in 2023 ({% cite Liao2023 %}).

According to {% cite Matthews2024 %}, there exist three major types of pangenomes. Here we focus on the **Presence-Absence Variation Pangenome** (see Figure 2) that was originally introduced by {% cite Tettelin2005 %} and whose description has been expanded by others over the years. This consists of:
- The **Core Genome**: Set of genes shared by all individuals of a species. This is sometimes strictly defined as 100% presence, or more loosely as a "**soft-core**" with 95% or 99% presence.
- And the **Accessory Genome**: Set of genes present in only a subset of individuals. This is further divided into
  - **Shell Genome**: Set of genes present in 1% to 99% of individuals 
  - and the **Cloud Genome**: Set of genes present in less than 1%.

<figure style="text-align: center;">
  <img src="./images/presence_absence_variation_pangenome.png" alt="Presence-Absence Variation Pangenome" style="width:100%;">
  <figcaption>
    <strong>Figure 2: A Presence-Absence Variation Pangenome</strong> (Modified after {% cite Matthews2024 %}). (A) A comparative genomic analysis of four related sequences is shown, utilizing one reference genome alongside three variants from the same population. Genetic variations from the reference are highlighted in color, and specific genes are denoted by black shapes. (B) The total pool of genes is divided into two distinct groups in a Presence-Absence Variation Pangenome: the core genes and accessory genes.
  </figcaption>
</figure>

> <agenda-title></agenda-title>
>
> In this tutorial, you will apply the pangenomic approach to your own data. You are going to analyze five *Listeria monocytogenes* strains by first annotating their genomes with Prokka. Afterward, you will build a Presence-Absence Variation Pangenome using Roary to identify the core genes and the accessory genes. To finish, you will use IQ-TREE on the core genes to construct a phylogenetic tree, showing how these five strains are evolutionarily related.
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data upload

To identify core and accessory genes across multiple bacterial strains, you will need the genome sequences in FASTA format of each strain.

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial. Give it a name like `Pangenome Analysis`.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
> 2. {% tool [Import](upload1) %} the following files from [Zenodo](https://doi.org/10.5281/zenodo.21339040).
>
>    ```
>    https://zenodo.org/records/21339041/files/10403S.fasta
>    https://zenodo.org/records/21339041/files/AUF.fasta
>    https://zenodo.org/records/21339041/files/EGD-e.fna
>    https://zenodo.org/records/21339041/files/EGD.fasta
>    https://zenodo.org/records/21339041/files/F6900.fna
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Create a dataset collection containing the uploaded FASTA files. Give it a name like `Assemblies`.
>
>    {% snippet faqs/galaxy/collections_build_list.md %}
>
{: .hands_on}

# Annotate the Strains with Prokka

You will use Prokka to annotate all bacterial strains simultaneously using the dataset collection. Prokka analyzes the DNA sequences to identify genomic features, such as genes, and assigns each a name and predicted biological function.

> <hands-on-title>Annotate Strains</hands-on-title>
>
> 1. Select `Tools` in the left sidebar and search for {% tool [Prokka](toolshed.g2.bx.psu.edu/repos/crs4/prokka/prokka/1.14.6+galaxy2) %} in the list that appears. Select it to open the tool. 
>
> 2. Within Prokka, select the following parameters (leave everything else unchanged):
>    - {% icon param-file %} *"Contigs to annotate"*: Click on the *"Dataset collection"* icon and select your `Assemblies` dataset collection
>    - *"Genus name"*: Enter `Listeria`
>    - *"Species name"*: Enter `monocytogenes` 
>    - *"Additional outputs"*: Deselect all options except `Annotation in GFF3 format, containing both sequences and annotations (.gff)`
>
> 3. Run the tool. 
>
{: .hands_on}

## Examine the Output

Once Prokka has finished, two new dataset collections will appear in your history:

- One containing the GFF3 output files for each strain. This is the dataset collection you will use further on. Examine each output file in this collection.
- One containing the logs of each run.

> <details-title>What is a GFF3 file?</details-title>
>
> **GFF3** stands for **General Feature Format (version 3)**. It is a plain text file used in bioinformatics to describe the exact locations of genomic features. The file contains 9 tab-separated columns that provide information for each feature ([official specification](http://www.sequenceontology.org/gff3.shtml)).
>
> Standard GFF3 files only contain the coordinates of genomic features. However, Prokka attaches the raw DNA sequence in FASTA format at the bottom of the GFF3 file (starting with a ##FASTA line in the file). This is exactly what the next tool requires to build the pangenome.
>
{: .details}

# Calculate the Pangenome with Roary

Now you will run Roary to compare the genes identified across all strains by Prokka. Roary clusters the genes according to sequence similarity, thereby identifying the core and accessory genes that constitute the pangenome.

> <hands-on-title>Identify Core and Accessory Genes</hands-on-title>
>
> 1. Select `Tools` in the left sidebar and search for {% tool [Roary](toolshed.g2.bx.psu.edu/repos/iuc/roary/roary/3.13.0+galaxy3) %} in the list that appears. Select it to open the tool. 
>
> 2. Within Roary, select the following parameters (leave everything else unchanged):
>    - *"Individual gff files or a dataset collection"*: Select `Collection` in the dropdown
>      - {% icon param-file %} *"Dataset collection to submit to Roary"*: Select the collection containing the GFF3 output files from the previous step.
>    - *"Additional outputs"*: Tick `Accessory binary genes in newick format`
>
> 3. Run the tool. 
>
{: .hands_on}

## Examine the Output

Once Roary has finished, four different output files will appear in your history. Examine each output file:

- {% icon param-file %} **Summary statistics** file: Gives you a quick overview of how many genes belong to the core genome and the accessory genome.
- {% icon param-file %} **Gene Presence Absence** file: Shows which genes are found in which strains. 
- {% icon param-file %} **Core Gene Alignment** file: Contains the aligned sequences of the core genes in FASTA format.
- {% icon param-file %} **Accessory Binary Genes (Newick)** file: Contains a phylogenetic tree based on the presence and absence of the accessory genes.

> <question-title></question-title>
>
> 1. How many genes are in the pangenome in total and how many of them are shell genes?
> 2. In which strains can the gene for the Putative Cyclic di-GMP Phosphodiesterase *PdeC* be found?
>
> > <solution-title></solution-title>
> >
> > 1. There are 3417 genes in total in the pangenome and, of these, 829 are shell genes. You can see this by viewing the {% icon param-file %} **Summary statistics** file. You will see that Roary defines the shell genome as a set of genes present in 15% to 95% of individuals and the cloud genome as a set of genes present in less than 15% of individuals.
> > 2. *PdeC* can be found in all five strains: *10403S*, *AUF*, *EGD*, *F6900*, and *EGD-e*. You can see this by viewing the {% icon param-file %} **Gene Presence Absence** file. Each row corresponds to an annotated gene and provides information on how many isolates contain it. At the very end of the row, you can identify the specific strains by the presence of their locus tags. For the full view of the file, you may need to download the file and view it on your computer.
> >
> {: .solution}
>
{: .question}

# Construct the Phylogenetic Tree with IQ-TREE

You can now use the {% icon param-file %} **Core Gene Alignment** file to build a phylogenetic tree from it. This will show you the evolutionary relationship between all strains based on their shared genetic backbone. For this, you will use **IQ-TREE**.

> <hands-on-title>Build a phylogenetic tree</hands-on-title>
>
> 1. Select `Tools` in the left sidebar and search for {% tool [IQ-TREE](toolshed.g2.bx.psu.edu/repos/iuc/iqtree/iqtree/2.4.0+galaxy2) %} in the list that appears. Select it to open the tool. 
>
> 2. Within IQ-TREE, select the following parameters (leave everything else unchanged):
>    - {% icon param-file %} *"Specify input alignment file in PHYLIP, FASTA, NEXUS, CLUSTAL or MSF format."*: Select the `Core Gene Alignment` file
>
> 3. Run the tool. 
>
{: .hands_on}

## Examine the Output

Once IQ-TREE has finished, five different output files will appear in your history:

- {% icon param-file %} **BIONJ Tree** file: This is an initial draft tree.
- {% icon param-file %} **MaxLikelihood Tree** file: Your main result. It contains the final phylogenetic tree in Newick format.
- {% icon param-file %} **MaxLikelihood Distance Matrix** file: Shows the calculated evolutionary distance between every possible pair of your strains.
- {% icon param-file %} **Occurrence Frequencies in Bootstrap Trees** file: Contains the statistical support for the tree's branches.
- {% icon param-file %} **Report and Final Tree** file: A summary of the entire run and a text-based visual of the tree.

## Visualize the Phylogenetic Tree

Now that you have a phylogenetic tree from your bacterial strains, you can view it visually. 

> <hands-on-title>Visualize the Phylogenetic Tree</hands-on-title>
>
> 1. Select the {% icon param-file %} **MaxLikelihood Tree** file in your history and click the {% icon galaxy-eye %} **View** icon.
>
> 2. By default, the file is displayed in {% icon galaxy-eye %} **Preview** mode. In the top tab bar, select {% icon galaxy-barchart %} **Visualize**. 
>
> 3. Select the `Phylogenetic Tree Browser`
>
{: .hands_on}

The resulting phylogenetic tree should resemble Figure 3. Zoom in to examine individual branches.

<figure style="text-align: center;">
  <img src="./images/phylo_tree.png" alt="Phylogenetic tree of the five analyzed bacterial strains based on core genome" style="width:45%;">
  <figcaption>
    <strong>Figure 3: Phylogenetic tree of the five analyzed bacterial strains based on the core genome.</strong> Visualization of the evolutionary relationship and genetic divergence among the selected strains. The topology reveals two distinct groupings: (A) Strains 10403S, AUF, and EGD cluster closely together (their short branch lengths confirm that they share a near-identical core genome), and (B) Strains F6900 and EGD-e, which branch off into a separate cluster. This demonstrates a higher genetic distance and evolutionary divergence, both from the first group and from one another.
  </figcaption>
</figure>

> <comment-title>What can we deduce from these results?</comment-title>
>
> - Strains that appear highly similar in this tree may nonetheless differ substantially in phenotype. This is because variations can exist within their accessory genes.
> - For example, one strain might carry plasmids with antibiotic resistance, while another might harbor toxin genes introduced by phages.
> - These differences are hidden in this specific phylogenetic tree because it was built on core genes.
> - A core genome tree is precise for tracing evolutionary descent and relatedness, but it does not reveal the complete individual genetic repertoire.
{: .comment}

# Conclusion
This tutorial provided a step-by-step guide on how to construct a bacterial Presence-Absence Variation Pangenome by annotating the genomes with Prokka and extracting the core and accessory genes with Roary. Additionally, it provided a guide for creating a phylogenetic tree with IQ-TREE based on the core genes. By following these steps, you should be able to identify core and accessory genes, and uncover the evolutionary relationships among your own bacterial strains.