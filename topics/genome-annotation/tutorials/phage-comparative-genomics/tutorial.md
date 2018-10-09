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

As part of the functional workflow, BLASTp (for proteins) and BLASTn (for nucleotide sequences) have already been run; thus, both of these kinds of comparisons have already been made. The BLAST results will be compared against the NCBI NR and NT databases for proteins and DNA, respectively. These are the most current and comprehensive databases as they reflect the internationally shared [INSDC database.](http://www.insdc.org/)

> ### {% icon comment %} A Brief Center for Phage Technology Precedent
> An overview of how the [CPT](https://cpt.tamu.edu/) currently organizes and classifies genomes is provided in [this publication.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5408676/pdf/viruses-09-00070.pdf)
{: .tip}

## 1. DNA Sequence Comparisons

Because the triplet code that encodes protein sequences is *degenerate* (as in, there are multiple possible triplet codons for most amino acids), a DNA sequence can drift and still encode the same protein sequence by accumulation of silent mutations. This means that DNA sequences encoding similar proteins can diverge at a relatively high rate, thus DNA sequence comparisons are only particularly useful for organisms that are **closely related.** Between related organisms, DNA sequence comparisons can provide good overall [parallels, as this single analysis can demonstrate both conservation of DNA sequences and genome *synteny* (the order of genes in the genome). Additionally, high conservation of DNA sequences automatically means that proteins encoded by that sequence must also be similar. Once DNA sequence similarity drops below ≈ 30%, it is no longer very useful for comparisons.

## 2. Protein Sequence Comparisons

Comparison of organisms by protein sequence is much more sensitive, as there is stronger pressure to conserve a protein sequence for the protein to retain its function. Phage genome organization is considered to be *modular*, meaning that individual genes or groups of genes can be shared across otherwise very different phages. It is not unusual for two phages to be very similar across the genome with only particular genes - such as phage tail fibers - being different.

> ### {% icon comment %} Modular Phage Gene Examples
> Comparing phages lambda and P22, one will see that they share similar integration and lysogen genes, *but* different morphogenesis genes (siphophage versus podophage, respectively).
> Comparing phages lambda and T1, one will see that they have related morphogenesis genes, *but* different modules for control of replication and lysis.
{: .tip}

# Workflow

> * Open CPT Galaxy ([CPT Public Galaxy](https://cpt.tamu.edu/galaxy-pub), [CPT TAMU Galaxy](https://cpt.tamu.edu/galaxy)), and find the history that contains the results of the [functional workflow](LINK TUTORIAL!) for the desired phage genome.
> * Before retrieving the workflow, there are workflow outputs that need to be unhidden for the BLASTn analysis. At the top of the history, click on the “hidden” hyperlink to see all workflow outputs that have been hidden. Scroll down and find the “NT” dataset; this is the output from BLASTn against the NCBI nt database. Click “Unhide it” on this dataset and it will become available for analysis. If preferred, click “hide hidden” at the top of the history to re-hide the other datasets to prevent clutter.
> * At the top of the web page, click on the Shared Data drop-down menu and select workflows.

![](../../images/phage-comparative-genomics-screenshots/1_go_to_workflows.png)

> * Find the most recent version of the “Phage comparative genomics (v1.#) workflow,” where # is the highest number indicating the most recent version of this workflow. Click on the drop-down menu for that workflow and select “Import.” After this, a green box containing a message will appear to inform the user of a successful import.

![](../../images/phage-comparative-genomics-screenshots/2_import_workflow.png)

![](../../images/phage-comparative-genomics-screenshots/3_successful_import.png)

> * From there, one can click on “start using this workflow” within the message box to be brought to the page containing all of the user’s imported workflows. Find the Phage comparative genomics workflow, click the drop-down menu, and select “Run.”

![](../../images/phage-comparative-genomics-screenshots/4_run_workflow.png)

