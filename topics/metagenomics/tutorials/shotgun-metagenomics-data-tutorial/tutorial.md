---
layout: tutorial_hands_on

title: "Analysis of shotgun metagenomics data"
zenodo_link: "https://doi.org/10.5281/zenodo.815875"
questions:
  - "How to analyze metagenomics shotgun data?"
objectives:
  - "Selection of tools to analyze shotgun data"
  - "Extracting taxonomic and functional information"
  - "Visualisation of a community structure"
time_estimation: "30M"
key_points:
  - "With amplicon data, we can extract information about the studied community structure"
  - "With shotgun data, we can extract information about the studied community structure and also the functions realised by the community"
  - "The tools used to analyze amplicon and shotgun data are different, except for the visualisation"
  - "Metagenomics data analyses are complex and time-consuming"
contributors:
  - shiltemann
  - bebatut
  - EngyNasr
---

# Introduction
{:.no_toc}

In metagenomics, information about micro-organisms in an environment can be extracted with two main techniques:

- [Amplicon sequencing](https://training.galaxyproject.org/training-material/topics/metagenomics/tutorials/mothur-miseq-sop/tutorial.html), which sequences only the rRNA or ribosomal DNA of organisms
- __Shotgun sequencing__, which sequences full genomes of the micro-organisms in the environment

In this tutorial, we will introduce the second type of analysis, shotgun data, with its general principles. For a more in-depth look at these analyses, we recommend our detailed tutorials on each analysis, Amplicon and Shotgun.

We will use the dataset from [project on the Argentinean agricultural pampean soils](https://www.ebi.ac.uk/metagenomics/projects/SRP016633). In this project, three different geographic regions that are under different types of land uses and two soil types (bulk and rhizospheric) were analyzed using shotgun and amplicon sequencing. We will focus on data from the Argentina Anguil and Pampas Bulk Soil (the original study included one more geographical regions, [see](https://doi.org/10.1186/2049-2618-1-21)).

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Shotgun metagenomics data

In [16S Microbial Analysis with mothur tutorial](https://training.galaxyproject.org/training-material/topics/metagenomics/tutorials/mothur-miseq-sop/tutorial.html), we saw how to analyze amplicon data to extract the community structure. Such information can also be extracted from shotgun metagenomic data.

In shotgun data analysis, full genomes of the micro-organisms in the environment are sequenced (not only the 16S or 18S). We can then have access to the rRNA (only a small part of the genomes), but also to the other genes of the micro-organisms. Using this information, we can try to answer questions such as "What are the micro-organisms doing?" in addition to the question "What micro-organisms are present?".

We will use a metagenomic sample of the Pampas Soil ([SRR606451](https://www.ebi.ac.uk/metagenomics/projects/SRP016633/samples/SRS372043/runs/SRR606451/results/versions/2.0)).

## Data upload

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history
> 2. Import the `SRR606451_pampa` Fasta file from [Zenodo](http://zenodo.org/record/815875) or from the data library (in "Analyses of metagenomics data")
>
>    ```
>    https://zenodo.org/record/815875/files/SRR606451_pampa.fasta
>    ```
>
{: .hands_on}

## Extraction of taxonomic information

As for amplicon data, we can extract taxonomic and community structure information from shotgun data. Different approaches can be used:

- Same approach as for amplicon data with identification and classification of OTUs

    Such an approach requires a first step of sequence sorting to extract only the 16S and 18S sequences, then using the same tools as for amplicon data. However, rRNA sequences represent a low proportion (< 1%) of the shotgun sequences so such an approach is not the most statistically supported

- Assignation of taxonomy on the whole sequences using databases with marker genes

In this tutorial, we use the second approach with MetaPhlAn2. This tools is using a database of ~1M unique clade-specific marker genes (not only the rRNA genes) identified from ~17,000 reference (bacterial, archeal, viral and eukaryotic) genomes.

> ### {% icon hands_on %} Hands-on: Taxonomic assignation with MetaPhlAn2
>
> 1. **MetaPhlAN2** {% icon tool %} with
>    - "Input file" to the imported file
>    - "Database with clade-specific marker genes" to `locally cached`
>    - "Cached database with clade-specific marker genes" to `MetaPhlAn2 clade-specific marker genes`
>
> This step may take a couple of minutes.
{: .hands_on}

3 files are generated:

- A tabular file with the community structure

    ```
    #SampleID   Metaphlan2_Analysis
    k__Bacteria 100.0
    k__Bacteria|p__Proteobacteria   86.20712
    k__Bacteria|p__Actinobacteria   13.79288
    k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria    86.20712
    k__Bacteria|p__Actinobacteria|c__Actinobacteria 13.79288
    ```

    Each line contains a taxa and its relative abundance found for our sample. The file starts with high level taxa (kingdom: `k__`) and go to more precise taxa.


- A BIOM file with the same information as the previous file but in BIOM format

    It can be used then by mothur and other tools requiring community structure information in BIOM format

- A SAM file with the results of the mapping of the sequences on the reference database

> ### {% icon question %} Questions
>
> 1. What is the most precise level we have access to with MetaPhlAn2?
> 2. What are the two orders found in our sample?
> 3. What is the most abundant family in our sample?
>
> > ### {% icon solution %} Solution
> > 1. We have access to species level
> > 2. Pseudomonadales and Solirubrobacterales are found in our sample
> > 3. The most abundant family is Pseudomonadaceae with 86.21 % of the assigned sequences
> {: .solution }
{: .question}

Even if the output of MetaPhlAn2 is bit easier to parse than the BIOM file, we want to visualize and explore the community structure with KRONA

> ### {% icon hands_on %} Hands-on: Interactive visualization with KRONA
>
> 1. **Format MetaPhlAn2 output for Krona** {% icon tool %} with
>    - "Input file" to `Community profile` output of `MetaPhlAn2`
>
> 2. **KRONA pie chart** {% icon tool %} with
>    - "What is the type of your input data" as `MetaPhlan`
>    - "Input file" to the output of `Format MetaPhlAn2`
>
{: .hands_on}

## Extraction of functional information

We would now like to answer the question "What are the micro-organisms doing?" or "Which functions are performed by the micro-organisms in the environment?".

In the shotgun data, we have access to the gene sequences from the full genome. We use that to identify the genes, associate them with a function, build pathways, etc., to investigate the functional part of the community.

> ### {% icon hands_on %} Hands-on: Metabolism function identification
>
> 1. **HUMAnN2** {% icon tool %} with
>    - "Input sequence file" to the imported sequence file
>    - "Use of a custom taxonomic profile" to `Yes`
>    - "Taxonomic profile file" to `Community profile` output of `MetaPhlAn2`
>    - "Nucleotide database" to `Locally cached`
>    - "Nucleotide database" to `Full`
>    - "Protein database" to `Locally cached`
>    - "Protein database" to `Full UniRef50`
>    - "Search for uniref50 or uniref90 gene families?" to `uniref50`
>    - "Database to use for pathway computations" to `MetaCyc`
>    - "Advanced Options"
>    - "Remove stratification from output" to `Yes`
>
>    This step is long so we generated the output for you!
>
> 2. Import the 3 files whose the name is starting with "humann2"
>
>    ```
>    https://zenodo.org/record/815875/files/humann2_gene_families_abundance.tsv
>    https://zenodo.org/record/815875/files/humann2_pathways_abundance.tsv
>    https://zenodo.org/record/815875/files/humann2_pathways_coverage.tsv
>    ```
{: .hands_on}

HUMAnN2 generates 3 files

- A file with the abundance of gene families

    Gene family abundance is reported in RPK (reads per kilobase) units to normalize for gene length. It reflects the relative gene (or transcript) copy number in the community.

    The "UNMAPPED" value is the total number of reads which remain unmapped after both alignment steps (nucleotide and translated search). Since other gene features in the table are quantified in RPK units, "UNMAPPED" can be interpreted as a single unknown gene of length 1 kilobase recruiting all reads that failed to map to known sequences.

- A file with the coverage of pathways

    Pathway coverage provides an alternative description of the presence (1) and absence (0) of pathways in a community, independent of their quantitative abundance.

- A file with the abundance of pathways

> ### {% icon question %} Questions
>
> How many gene families and pathways have been identified?
>
> > ### {% icon solution %} Solution
> > 44 gene families but no pathways are identified
> {: .solution }
{: .question}

The RPK for the gene families are quite difficult to interpret in term of relative abundance. We decide then to normalize the values

> ### {% icon hands_on %} Hands-on: Normalize the gene family abundances
>
> 1. **Renormalize a HUMAnN2 generated table** {% icon tool %} with
>    - "Gene/pathway table" to the gene family table generated with `HUMAnN2`
>    - "Normalization scheme" to `Relative abundance`
>    - "Normalization level" to `Normalization of all levels by community total`
>
>  > ### {% icon question %} Questions
>  >
>  > 1. What percentage of sequences has not be assigned to a gene family?
>  > 2. What is the most abundant gene family?
>  >
>  > > ### {% icon solution %} Solution
>  > > 1. 55% of the sequences have not be assigned to a gene family
>  > > 2. The most abundant gene family with 25% of sequences is a putative secreted protein
>  > {: .solution }
>  {: .question}
{: .hands_on}

With the previous analyses, we investigate "Which micro-organims are present in my sample?" and "What function are performed by the micro-organisms in my sample?". We can go further in these analyses (for example, with a combination of functional and taxonomic results). We did not detail that in this tutorial but you can find more analyses in our tutorials on shotgun metagenomic data analyses.

# Conclusion
{:.no_toc}

We can summarize the analyses with amplicon and shotgun metagenomic data:

![Scheme to sum up the analysis](../../images/general-tutorial-scheme.png)

Both analyses are quite complex! However, in this tutorial, we only showed simple cases of metagenomics data analysis with subset of real data.

Check our other tutorials to learn more in detail of how to analyze metagenomics data.
