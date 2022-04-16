---
title: How to use Custom Reference Genome?
area: datasets
box_type: tip
layout: faq
contributors: [jennaj, Nurzhamalyrys]
---

- A reference genome contains the nucleotide sequence of the chromosomes, scaffolds, transcripts, or contigs for single species. It is representative of a specific genome build or release.

- In Galaxy, a custom reference genome is a FASTA formatted dataset that can be used in place of a native reference genome with most tools.

 - custom: a dataset from the history loaded by users
 - native: local or cached by administrators 


There are five basic steps to using a Custom Reference Genome:

1. Obtain a FASTA copy of the target genome.
2. FTP the genome to Galaxy and load into a history as a dataset.
3. Clean up the format with the tool NormalizeFasta using the options to wrap sequence lines at 80 bases and to trim the title line at the first whitespace.
4. Make sure the chromosome identifiers are a match for other inputs.
5. Set a tool form's options to use a custom reference genome from the history and select the loaded genome.
