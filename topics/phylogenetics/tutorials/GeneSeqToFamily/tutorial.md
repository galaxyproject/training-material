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

The Ensembl GeneTrees computational pipeline generates gene families based on coding sequences. 
It uses various tools: 
BLAST
hcluster_sg
T-Coffee
TreeBeST

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

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

## GeneSeqToFamily preparation 
GeneSeqToFamily preparation is an open-source tool that converts genomic information in GFF/JSON format to SQLite format for easy access during the workflow, it also adds species information in Fasta header.


> ### {% icon hands_on %} Hands-on: GeneSeqToFamily preparation : Run GeneSeqToFamily preparation on the imported GFF/JSON and Fasta files
>
> 1. Open the GeneSeqToFamily preparation tool
> 2. Select JSON and/or GFFs files
> 3. Add specific species name (in-case of GFFs)
> 4. Select All Fasta Files 
> 5. Add Chromosome/Reference names to filter out sequences such as, MT, Chloroplast etc
> 6. Select Keep Longest CDS per gene
> 7. Run tool
>
>	> ### {% icon tip %} Tip: Importing data via links
>	>
>	> * Open the GeneSeqToFamily preparation tool
>	> * Select JSON (if available)
>	> * Select GFFs one by one
>	> * Add specific species name
>	> * Select All Fasta Files 
>	> * Add Chromosome/Reference names to filter out sequences such as, MT, Chloroplast etc
>	> * Press Run
>	{: .tip}
>
>
>	> ### {% icon comment %} Comments
>	>
>	> Filter out sequences which has special codon translation.
>	{: .comment}
{: .hands_on}

# CDS Translation
## Transeq
> ### {% icon hands_on %} Hands-on: Transeq
>
> 1. Frame(s) to translate: 1
> 2. Code to use: Standard
> 3. Remove all 'X' and '*' characters from the right end of the translation
> 4. Change all STOP codon positions from the '*' character to 'X'
> 5. Output sequence file format
>
> 
{: .hands_on}
# Preclustering alignment
## BLAST DB
> ### {% icon hands_on %} Hands-on: Create BLAST DB: Create BLAST database
>
> 1. Open the createblastdb tool
> 2. Select translated sequences
> 3. Run tool
{: .hands_on}
## BLASTP
> ### {% icon hands_on %} Hands-on: BLASTP : Run BLASTP
>
> 1. Open the blastp tool
> 2. Select protein query sequence(s): translated sequences
> 3. Subject database / sequence(s): BLAST database from your history
> 4. Protein BLAST database
> 5. Type of BLAST: blastp - Traditional BLASTP to compare a protein query to a protein database
> 6. Set expectation value cutoff: 1e-10
> 7. Output format: Tabular (extended 25 columns)
> 8. Advanced Options: Show Advanced Options
> 9. Maximum number of HSPs (alignmnets) to keep for any single query-subject pair
> 7. Run tool
{: .hands_on}
## BLAST parser
> ### {% icon hands_on %} Hands-on: BLAST Parser 
>
> 1. Open the BLAST parser tool
> 2. Reciprocal results: Yes
> 3. Run tool
{: .hands_on}
# Cluster generation
## hcluster_sg
> ### {% icon hands_on %} Hands-on: hcluster_sg 
>
> 1. Open the hcluster_sg
> 2. Only find single-linkage clusters: No
> 3. Minimum edge density between a join: 0.34
> 4. Maximum size: 500
> 5. Run tool
{: .hands_on}
## hcluster_sg parser
> ### {% icon hands_on %} Hands-on: hcluster_sg parser
>
> 1. Open the hcluster_sg parser
> 2. Minimum number of cluster elements: 3
> 3. Maximum number of cluster elements: 200
> 5. Run tool
{: .hands_on}
## Filter by FASTA IDs
# Cluster Alignment
## T-Coffee
> ### {% icon hands_on %} Hands-on: T-Coffee 
>
> 1. Open the T-Coffee
> 2. Filter FASTA input?: Yes
> 3. Multiple Sequence Alignment Methods: clustalw_msa
> 4. Output formats: fasta_aln
> 5. Run tool
{: .hands_on}
# Gene tree construction
## Tranalign
> ### {% icon hands_on %} Hands-on: Tranalign
>
> 1. Open the Tranalign
> 2. Nucleic sequences: FASTA sequences generated from GeneSeqToFamily preparation
> 3. Protein sequences: Alignment generated from T-Coffee
> 4. Code to use: standard
> 5. Output sequence file format: FASTA (m)
> 6. Run tool
{: .hands_on}
## TreeBeST "best"
> ### {% icon hands_on %} Hands-on: TreeBeST best
>
> 1. Open the TreeBeST best
> 2. Species file in Newick format: Select input species file
> 3. CDS alignment: FASTA alignment generated from Tranalign 
> 4. Run tool
{: .hands_on}

# Gene alignment and family aggregation
## Gene Align and Family Aggregator 
> ### {% icon hands_on %} Hands-on: Gene Align and Family Aggregator 
>
> 1. Open the Gene Align and Family Aggregator 
> 2. Gene tree: GeneTrees generated from TreeBeST best
> 3. Protein alignments: Alignments generated from T-Coffee
> 2. Gene features: Gene features SQLite generated from GeneSeqToFamily preparation
> 6. Run tool 
{: .hands_on}
# Visualization
## Aequatus visualization plugin
> ### {% icon hands_on %} Hands-on: equatus visualization plugin
>
> 1. Click on Result generated by previous step.
> 2. Choose GeneTree from side panel
> 4. Visualise different GeneTrees
{: .hands_on}


# Conclusion
{:.no_toc}

blabla
