---
layout: tutorial_hands_on
topic_name: phylogenetics
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
{:.no_toc}

The Ensembl GeneTrees computational pipeline generates gene families based on coding sequences. 
It uses various tools: 
BLAST
hcluster_sg
T-Coffee
TreeBeST

> ## Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Importing dataset

> ## {% icon hands_on %} Hands-on: Data upload
>
> Make sure you have an empty analysis history. Give it a name.
>	> ### {% icon tip %} Tip: Starting a new history
>	>
>	> * Click the gear icon at the top of the history panel
>	> * Select the option Create New from the menu
>	{: .tip}
> 
> 2. Import Sample Data
>	> * FASTA file
>	> * GFF file: 
>	> * Species tree
>
>	> ### {% icon tip %} Tip: Importing data via links
>	>
>	> * Copy the link locations
>	> * Open the Galaxy Upload Manager
>	> * Select Paste/Fetch Data
>	> * Paste the link into the text field
>	> * Press Start
>	>
>	{: .tip}
>
>	> ### {% icon tip %} Tip: Change the file type text to nhx once the data file is in your history
>	>
>	> Click on the pencil button {% icon hands_on %} displayed in your data file in the history
>	> * Choose Datatype on the top
>	> * Select nhx
>	> * Press save
>	>
>	{: .tip}
>
>	> ### {% icon comment %} 
> 	> Rename the dataset to “First dataset”
>	{: .comment}
> By default, when data is imported via its link, Galaxy names it with its URL.
{: .hands_on}


# Data Preparation

To convert uploaded data into the format acceptable by GeneSeqToFamily workflow:

## GeneSeqToFamily preparation 
GeneSeqToFamily preparation is an open-source tool that converts genomic information in GFF/JSON format to SQLite format for easy access during the workflow, it also adds species information in Fasta header.

> ### {% icon hands_on %} Hands-on: GeneSeqToFamily preparation : Run GeneSeqToFamily preparation on the imported GFF/JSON and Fasta files
>
> 1. GeneSeqToFamily preparation {% icon tool %}
>	> * Select JSON and/or GFFs files
>	> * Add specific species name (in-case of GFFs)
>	> * Select All Fasta Files 
>	> * Keep Longest CDS per gene: Yes
>	> * Change the header line of the FASTA sequences to the >TranscriptId_species format: Yes
>	> * Comma-separated list of region IDs (e.g. chromosomes or scaffolds) for which FASTA sequences should be filtered:
>	> * Run tool
>
>	> ### {% icon comment %} Comments
>	>
>	> Filter out sequences which has special codon translation.
>	{: .comment}
{: .hands_on}

# CDS Translation

We use Transeq to convert a CDS to protein sequences in order to run BLASTP and find protein clusters. However, since downstream tools in the pipeline, such as TreeBeST, require nucleotide sequences to generate a gene tree, the protein sequences cannot be directly used as workflow input and are instead generated with Transeq.

> ### {% icon hands_on %} Hands-on: Transeq
>
> 1. Transeq {% icon tool %}
>	> * Frame(s) to translate: 1
>	> * Code to use: Standard
>	> * Remove all 'X' and '*' characters from the right end of the translation
>	> * Change all STOP codon positions from the '*' character to 'X'
>	> * Output sequence file format
>
>	> ### {% icon tip %} Tip
>	>
>	> Change Code to use from standard to other based on sequence type.
>	{: .tip}
{: .hands_on}

# Preclustering alignment

We are using BLASTP to run over the set of sequences against the database of the same input, as is the case with BLAST-all, in order to form clusters of related sequences.

## BLAST Database

> ### {% icon hands_on %} Hands-on: makeblastdb
>
> 1. NCBI BLAST+ makeblastdb {% icon tool %}
>	> * Molecule type of input: protein
>	> * Input FASTA files(s):
>	> * Run tool
{: .hands_on}

## BLASTP

> ### {% icon hands_on %} Hands-on: BLASTP : Run BLASTP
>
> 1. NCBI BLAST+ blastp {% icon tool %}
>	> * Protein query sequence(s): translated sequences
>	> * Subject database / sequence(s): BLAST database from your history
>	> * Protein BLAST database
>	> * Type of BLAST: blastp - Traditional BLASTP to compare a protein query to a protein database
>	> * Set expectation value cutoff: 1e-10
>	> * Output format: Tabular (extended 25 columns)
>	> * Advanced Options: Show Advanced Options
>	> * Maximum number of HSPs (alignmnets) to keep for any single query-subject pair
>	> * Run tool
>
>	> ### {% icon tip %} Tip
>	>
>	> Change BLASTP parameters to achieve expected result.
>	{: .tip}
{: .hands_on}

## BLAST parser

BLAST parser is a small Galaxy tool to convert the BLAST output into the input format required by hcluster_sg. 
It takes the BLAST 12 or 25-column output as input and generates a 3-column tabular file, comprising the BLAST query, the hit result, and the edge weight. The weight value is simply calculated as minus log10 of the BLAST e-value divided by 2, replacing this with 100 if this value is greater than 100. It also removes the self-matching BLAST results and lets the user filter out non-Reciprocal Best Hits (if selected).

> ### {% icon hands_on %} Hands-on: BLAST Parser 
>
> 1. BLAST parser tool {% icon tool %}
>	> * Reciprocal results: Yes
>	> * Run tool
{: .hands_on}

# Cluster generation

## hcluster_sg

hcluster_sg performs clustering for sparse graphs. It reads an input file that describes the similarity between 2 sequences, and iterates through the process of grouping 2 nearest nodes at each iteration. hcluster_sg outputs a single list of gene clusters.

> ### {% icon hands_on %} Hands-on: hcluster_sg 
>
> 1. hcluster_sg {% icon tool %}
>	> * Only find single-linkage clusters: No
>	> * Minimum edge density between a join: 0.34
>	> * Maximum size: 500
>	> * Run tool
{: .hands_on}

## hcluster_sg parser

hcluster_sg parser tool creates collection of files each containing sequence IDs for cluster.

> ### {% icon hands_on %} Hands-on: hcluster_sg parser
>
> 2. hcluster_sg parser {% icon tool %}
>	> * Minimum number of cluster elements: 3
>	> * Maximum number of cluster elements: 200
>	> * Run tool
>
>	> ### {% icon comment %} Comment
>	>
>	> Minimum number of cluster elements set to 3, because TreeBeST will not run if sequences are less than 3.
>	{: .comment}
{: .hands_on}

## Filter by FASTA IDs
Filter by FASTA IDs is used to create separate FASTA files using the sequence IDs listed in each gene cluster.

> ### {% icon hands_on %} Hands-on: hcluster_sg parser
>
> 1. filter_by_fasta_id {% icon tool %}
>	> * FASTA sequences
>	> * List of IDs to extract sequences for
>	> * Remove duplicate sequences: yes
>	> * Run tool
{: .hands_on}

# Cluster Alignment

## T-Coffee

T-Coffee is a MSA package, it can align both nucleotide and protein sequences. We use it to align the protein sequences in each cluster generated by hcluster_sg.

> ### {% icon hands_on %} Hands-on: T-Coffee 
>
> 1. T-Coffee {% icon tool %}
>	> * Filter FASTA input?: Yes
>	> * Multiple Sequence Alignment Methods: clustalw_msa
>	> * Output formats: fasta_aln
>	> * Run tool
{: .hands_on}

# Gene tree construction

## Tranalign

Tranalign is a tool that reads a set of nucleotide sequences and a corresponding aligned set of protein sequences and returns a set of aligned nucleotide sequences. Here, we use it to generate CDS alignments of gene sequences using the protein alignments produced by T-Coffee.

> ### {% icon hands_on %} Hands-on: Tranalign
>
> 1. Tranalign {% icon tool %}
>	> * Nucleic sequences: FASTA sequences generated from GeneSeqToFamily preparation
>	> * Protein sequences: Alignment generated from T-Coffee
>	> * Code to use: standard
>	> * Output sequence file format: FASTA (m)
>	> * Run tool
>
>	> ### {% icon tip %} Tip
>	>
>	> Change Code to use from standard to other based on sequence type.
>	{: .tip}
{: .hands_on}

## TreeBeST "best"

TreeBeST (Tree Building guided by Species Tree) is a tool to generate, manipulate, and display phylogenetic trees and can be used to build gene trees based on a known species tree.

> ### {% icon hands_on %} Hands-on: TreeBeST best
>
> 1. TreeBeST best {% icon tool %}
>	> * Species file in Newick format: Select input species file
>	> * CDS alignment: FASTA alignment generated from Tranalign 
>	> * Run tool
{: .hands_on}

# Gene alignment and family aggregation

## Gene Align and Family Aggregator 

Gene alignment and family aggregator (GAFA) is a Galaxy tool that generates a single SQLite database containing the gene trees and MSAs, along with gene features, in order to provide a reusable, persistent data store for visualization of synteny information with Aequatus.

> ### {% icon hands_on %} Hands-on: Gene Align and Family Aggregator 
>
> 1. Gene Align and Family Aggregator {% icon tool %}
>	> * Gene tree: GeneTrees generated from TreeBeST best
>	> * Protein alignments: Alignments generated from T-Coffee
>	> * Gene features: Gene features SQLite generated from GeneSeqToFamily preparation
>	> * Run tool 
{: .hands_on}

# Visualisation

## Aequatus visualisation Plugin 

The SQLite database generated by the GAFA tool can be rendered using a new visualization plugin, Aequatus.js. The Aequatus.js library, developed as part of the Aequatus project, has been configured to be used within Galaxy to visualize homologous gene structure and gene family relationships. 

> ### {% icon hands_on %} Hands-on: Aequatus visualization plugin 
>
> 1. Aequatus visualisation Plugin {% icon tool %}
>	> * Click on Result generated by previous step
>	> * Choose GeneTree from side panel
>	> * Visualise different GeneTrees
{: .hands_on}


# Conclusion
{:.no_toc}

Here, we convered all steps which makes GeneSeqToFamily workflow. In this tutorial we used default configuration of the workflow. They might need to be changed for various source of data.
