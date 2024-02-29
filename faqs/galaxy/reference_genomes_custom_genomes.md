---
title: How to use Custom Reference Genomes?
area: reference genomes
box_type: tip
layout: faq
contributors: [jennaj, Nurzhamalyrys]
---


A **reference genome** contains the nucleotide sequence of the chromosomes, scaffolds, transcripts, or contigs for single species. It is representative of a specific genome assembly build or release.

 
There are two options for reference genomes in Galaxy.
* **Native**
   * Index provided by the server administrators.
   * Found on tool forms in a drop down menu.
   * A database key is automatically assigned. See tip 1.
   * The database is what links your data to a FASTA index. Example: used with BAM data
* **Custom** 
   * FASTA file uploaded by users. 
   * Input on tool forms then indexed at runtime by the tool.
   * An optional custom database key can be created and [assigned by the user]({% link faqs/galaxy/datasets_change_dbkey.md %}).

There are five basic steps to use a **Custom Reference Genome**, plus one optional.
1. Obtain a FASTA copy of the target genome. See tip 2.
2. Upload the genome to Galaxy and to add it as a dataset in your history.
3. [Clean up the format]({% link faqs/galaxy/datasets_working_with_fasta.md %}) with the tool **NormalizeFasta** using the options to wrap sequence lines at 80 bases and to trim the title line at the first whitespace.
4. Make sure the [chromosome identifiers]({% link faqs/galaxy/datasets_chromosome_identifiers.md %}) are a match for other inputs.
5. Set a tool form's options to use a custom reference genome from the history and select the loaded genome FASTA.
6. (Optional) Create a [custom genome build's database]({% link faqs/galaxy/analysis_add_custom_build.md %}) that you can [assign to datasets]({% link faqs/galaxy/datasets_change_dbkey.md %}).

{% icon tip %} TIP 1: Avoid [assigning a native database]({% link faqs/galaxy/datasets_change_dbkey.md %}) to uploaded data unless you confirmed the data are based on the [same exact genome assembly]({% link faqs/galaxy/datasets_chromosome_identifiers.md %}) or you [adjusted the data to be a match]({% link topics/introduction/tutorials/data-manipulation-olympics/tutorial.html %}) **first**!

{% icon tip %} TIP 2: When choosing your reference genome, consider [choosing your reference annotation]({% link faqs/galaxy/analysis_differential_expression_help.md %}) at the same time. Standardize the format of both as a preparation step. Put the files in a dedicated "reference data" history for easy reuse.
