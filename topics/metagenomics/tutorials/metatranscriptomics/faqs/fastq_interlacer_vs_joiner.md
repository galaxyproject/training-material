---
title: Why do we use FASTQ interlacer and not the FASTQ joiner?
box_type: question
layout: faq
contributors: [subinamehta]
---

The reason ASaiM-MT uses FASTQ-interlacer than FASTQ-joiner for combining forward and reverse reads is because the joiner tool combines the forward and reverse read sequence together while the interlacer puts the forward and reverse read sequences in the same file while retaining the entity of each read along with an additional file with unpaired sequences and it maintains the integrity of the reads while helping us distinguish between the forward and reverse reads.


