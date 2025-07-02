---
title: Finding BAM dataset identifiers
area: datasets
description: How to find the reference sequence identifiers inside of a BAM file
box_type: tip
layout: faq
contributors: [jennaj, Melkeb]
---

Explore the content of your BAM.
1. Run **Samtools: IdxStats** on your `bam` dataset.
2. The reference sequence identifiers inside the "BAM header" will be listed in the result report. 
3. The report is a summary of the BAM content that includes: reference sequence identifiers (chromosome names), their lengths, and a count of the reads mapping to that reference sequence within the BAM file.
4. Compare the sequence identifiers in your BAM file to the the sequence identifiers (aka "chrom" field) field in all other inputs: VCF, GTF, GFF3, BED, Interval, Tabular.
5. It is usually important to use the same reference assembly for all steps within the same analysis. If you discover differences, you may need to choose different reference data.

{% icon tip %} Notes
* This method will *not* work for "sequence-only" `bam` datasets, as these usually have no header and are not associated with a reference assembly yet.
