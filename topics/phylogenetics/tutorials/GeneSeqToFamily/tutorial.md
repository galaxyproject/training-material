---
layout: tutorial_hands_on
topic_name: sequence-analysis
tutorial_name: GeneSeqToFamily
---

# Introduction
{:.no_toc}

GeneSeqToFamily: Discovery and Visualisation of Homologous Genes and Gene Families Using Galaxy

A gene family is a set of several similar genes formed by duplication of a single original gene.

The study of homologous genes helps to understand the evolution of gene families
Various tools available but they provide overview of homology at family level information such as, MSOAR, OrthoMCL, HomoloGene

They can not provide information about structural changes within gene


# Ensembl GeneTree pipeline    

The Ensembl GeneTrees computational pipeline generates gene families based on coding sequences. 
It uses various tools: 
BLAST
hcluster_sg
T-Coffee
TreeBeST


# Importing dataset

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history
>
> 2. Import the FASTA file: 
> 3. Import the GFF file: 
> 4. Import the species tree in Newick file format: 
>
>	> ### {% icon tip %} Tip: Importing data via links
>	>
>	> * Copy the link locations
>	> * Open the Galaxy Upload Manager
>	> * Select Paste/Fetch Data
>	> * Paste the link into the text field
>	> * Press Start
>	{: .tip}
>
>	> ### {% icon tip %} 
>	>
>	> Tip: Change the file type text to nhx once the data file is in your history
>	>
>	> Click on the pencil button displayed in your data file in the history
>	> Choose Datatype on the top
>	> Select nhx
>	> Press save
>	{: .tip}
>
>	> ### {% icon comment %} 
> 	> Rename the dataset to “First dataset”
>	{: .comment}
> By default, when data is imported via its link, Galaxy names it with its URL.
{: .hands_on}


# Data Preparation

To convert uploaded data into the format acceptable by GeneSeqToFamily workflow:

GeneSeqToFamily preparation is an open-source tool that converts genomic information in GFF/JSON format to SQLite format for easy access during the workflow, it also adds species information in Fasta header.

## GeneSeqToFamily preparation 
FastQC : Run GeneSeqToFamily preparation on the imported GFF/JSON and Fasta files

    Tip: Inspecting the content of a file in Galaxy
        
	    Open the GeneSeqToFamily preparation tool
	    Select JSON (if available)
	    Select GFFs one by one
	    Add specific species name
	    Select All Fasta Files 
	    Add Chromosome/Reference names to filter out sequences such as, MT, Chloroplast etc
	    Press Run

	Filter out sequences which has special codon translation.



## GeneSeqToFamily workflow 
Find GeneSeqToFamily workflow
Load appropriate data in the fields
Set Max target per seq = 4
Run the workflow
Wait for it…… 

# Visualise

## Aequatus.js
New JavaScript library for visualisation of homologous genes, 
extracted from the standalone Aequatus software package

Detailed view of gene structure across gene families

Shared exons use the same colour in each representation



# Conclusion
{:.no_toc}

blabla
