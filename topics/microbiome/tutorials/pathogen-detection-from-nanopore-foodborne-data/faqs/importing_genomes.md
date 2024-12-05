---
title: From where can I import other genomes?
box_type: question
layout: faq
contributors: [EngyNasr]
redirect_from:
- /topics/metagenomics/tutorials/pathogen-detection-from-nanopore-foodborne-data/faqs/importing_genomes
---

In this tutorial, we used kalamari DB with the [full list of possible host sequences](https://github.com/lskatz/Kalamari/blob/master/src/Kalamari_v3.9.1.tsv) that can be removed. Reads are either tagged to map one of those species or are left unassigned. If the task at hand in the real world cannot be covered by those, you can also try another DB for Kraken2 that includes your species (or maybe retain unmapped reads from a read aligner such as Bowtie2, Minimap2…).

