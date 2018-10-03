---
layout: tutorial_hands_on
topic_name: genome_annotation
tutorial_name: phage-comparative-genomics
---
> ### Agenda
>
> In this tutorial, we will deal with:
>
> * Background
>    > 1. DNA Sequence Comparisons
>    > 2. Protein Sequence Comparisons
> * Workflow
>
{: .agenda}

# Background

The functional annotation workflow has been run, and a grip has been developed around annotating the novel phage genome, looking at relationships between the novel phage and other known phages can commence. The best way to do this is to compare the novel phage genome to sequences that are already deposited in the sequence databases. This comparison can be executed in two ways:

> * Compare the entire **DNA sequence** of the novel phage to the DNA sequences of other organisms
> * Compare the **protein sequences** of all the genes in the novel phage genome to the protein sequences of other organisms

As part of the functional workflow, BLASTp (for proteins) and BLASTn (for nucleotide sequences) have already been run; thus, both of these kinds of comparisons have already been made. The BLAST results will be compared against the NCBI NR and NT databases for proteins and DNA, respectively. These are the most current and comprehensive databases as they reflect the internationally shared INSDC database.