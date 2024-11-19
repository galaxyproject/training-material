---
title: Mismatched Chromosome identifiers and how to avoid them
area: datasets
box_type: tip
layout: faq
contributors: [jennaj, Melkeb]
---

Reference data mismatches are similiar to bad reagents in a wet lab experiment: all sorts of odd problems can come up!

You inputs must be all based on an identical genome assembly build to achieve correct scientific results.

There are two areas to review for data to be considered **identical**.
1. The data are based on the same exact genome **assembly** (or "assembly release").
  * The "assembly" refers to the nucleotide sequence of the genome.
  * If the base order and length of the chromosomes are not the same, then your coordinates will have scientific problems.
  * Converting coordinates between assemblies may be possible. Search tool panel with `CrossMap`.
2. The data are based on the same exact genome assembly **build**.
  * The "build" refers to the labels used inside the file. In this context, pay attention to the chromosome identifiers.
  * These all may mean the same thing to a person but not to a computer or tool: chr1, Chr1, 1, chr1.1
  * Converting identifiers between builds may be possible. Search tool panel with `Replace`.
     
The methods listed below help to identify and correct errors or unexpected results when the underlying genome assembly build for all inputs are **not identical**.

**Method 1**: [Finding BAM dataset identifiers]({% link faqs/galaxy/datasets_BAM_dataset_identifiers.md %})

**Method 2**: [Directly obtaining UCSC sourced *genome* identifiers]({% link faqs/galaxy/datasets_UCSC_sourced_genome.md %})

**Method 3**: [Adjusting identifiers for UCSC sourced data used with other sourced data](https://galaxyproject.org/support/chrom-identifiers/#adjusting-identifiers-or-input-source)

**Method 4**: [Adjusting identifiers or input source for any mixed sourced data](https://galaxyproject.org/support/chrom-identifiers/#any-mixed-sourced-data)


{% icon tip %} Reference data is self referential. [More help for your genome, transcriptome, and annotation]({% link faqs/galaxy/analysis_differential_expression_help.md %})

{% icon tip %} Genome not available as a native index? [Use a custom genome fasta]({% link faqs/galaxy/reference_genomes_custom_genomes.md %}) and [create a custom build database]({% link faqs/galaxy/analysis_add_custom_build.md %}) instead.

{% icon tip %} More notes on Native Reference Genomes

* Native **reference genomes** (FASTA) are built as pre-computed indexes on the Galaxy server where you are working.
* Different servers host both common *and* different reference genome data.
* Most **reference annotation** (tabular, GTF, GFF3) is supplied from the history by the user, even when the genome is indexed.
* Public Galaxy servers source reference genomes preferentially from [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html).
* A **reference transcriptome** (FASTA) is supplied from the history by the user.
* Many experiements use a combination of all three types of reference data. Consider pre-preparing your files at the start!
* The default variant for a native genome index is "Full". Defined as: all primary chromosomes (or scaffolds/contigs) including mitochondrial plus associated unmapped, plasmid, and other segments.
* When only one version of a genome is available for a tool, it represents the default "Full" variant.
* Some genomes will have more than one variant available.
* The "Canonical Male" or sometimes simply "Canonical" variant contains the primary chromosomes for a genome. For example a human "Canonical" variant contains chr1-chr22, chrX, chrY, and chrM.
* The "Canonical Female" variant contains the primary chromosomes excluding chrY.
