---
title: How to use Custom Reference Genomes?
area: reference genomes
box_type: tip
layout: faq
contributors: [jennaj, Nurzhamalyrys]
---


A reference genome contains the nucleotide sequence of the chromosomes, scaffolds, transcripts, or contigs for single species. It is representative of a specific genome build or release. There are two options to use reference genomes in Galaxy: _native_ (provided by the server administrators and used by most of the tools) and _custom_ (uploaded by users in FASTA format). 

There are five basic steps to use a Custom Reference Genome:

1. Obtain a FASTA copy of the target genome.
2. Use FTP to upload the genome to Galaxy and load into a history as a dataset.
3. Clean up the format with the tool **NormalizeFasta** using the options to wrap sequence lines at 80 bases and to trim the title line at the first whitespace.
4. Make sure the chromosome identifiers are a match for other inputs.
5. Set a tool form's options to use a custom reference genome from the history and select the loaded genome.
