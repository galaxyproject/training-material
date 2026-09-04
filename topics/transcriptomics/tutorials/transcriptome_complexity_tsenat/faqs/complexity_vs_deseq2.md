---
title: Why is transcriptome complexity important if DESeq2 already detects changes?
area: conceptual
layout: faq
box_type: tip
contributors: [gallardoalba]
---

DESeq2 and edgeR measure **total gene abundance** changes. They miss situations where:
- A gene maintains constant total expression but **reorganizes its isoforms**
- Rare isoforms become more common (or vice versa) without affecting overall counts
- The **diversity** of the isoform landscape changes

These are valid biological signals (isoform switching via splicing regulation) that complement abundance-based methods.
