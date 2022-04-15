---
title: Find BAM dataset identifiers
area: datasets
box_type: tip
layout: faq
contributors: [jennaj, Melkeb]
---


Quickly learn what the identifiers are in any **BAM** dataset that is the result from mapping.

Note: This will *not* work for "sequence-only" `bam` datasets, as these usually have no header.

1. Run **Samtools: IdxStats** on the aligned data (`bam` dataset).
2. The "index header" chromosome names and lengths will be listed in the output (along with read counts).
3. Compare the chromosome identifiers to the chromosome (aka "chrom") field in all other inputs: VCF, GTF, GFF(3), BED, Interval, etc.
4. Note: The orignal mapping target may have been a built-in genome index, custom genome (transcriptome, exome, other) -- the same `bam` data will still be summarized.
