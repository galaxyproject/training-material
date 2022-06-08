---
title: Mismatched Chromosome identifiers and how to avoid them
area: datasets
box_type: tip
layout: faq
contributors: [jennaj, Melkeb]
---


- The methods listed here help to identify and correct errors or unexpected results linked to inputs having non-identical chromosome identifiers and/or different chromosome sequence content.

- **If using a Custom Reference genome**, the methods below also apply, but the first step is to make certain that the Custom Genome is formatted correctly. Improper formatting is the most common root cause of CG related errors.

Method 1: [Finding BAM dataset identifiers](https://training.galaxyproject.org/training-material/faqs/galaxy/datasets_BAM_dataset_identifiers.html)

Method 2: [Directly obtaining UCSC sourced *genome* identifiers](https://training.galaxyproject.org/training-material/faqs/galaxy/datasets_UCSC_sourced_genome.html)

Method 3: [Adjusting identifiers for UCSC sourced data used with other sourced data](https://galaxyproject.org/support/chrom-identifiers/#adjusting-identifiers-or-input-source)

Method 4: [Adjusting identifiers or input source for any mixed sourced data](https://galaxyproject.org/support/chrom-identifiers/#any-mixed-sourced-data)

**A Note on Built-in Reference Genomes**

- The default variant for all genomes is "Full", defined as all primary chromosomes (or scaffolds/contigs) including mitochondrial plus associated unmapped, plasmid, and other segments.
- When only one version of a genome is available for a tool, it represents the default "Full" variant.
- Some genomes will have more than one variant available.

  - The "Canonical Male" or sometimes simply "Canonical" variant contains the primary chromosomes for a genome. For example a human "Canonical" variant contains chr1-chr22, chrX, chrY, and chrM.
  - The "Canonical Female" variant contains the primary chromosomes excluding chrY.
