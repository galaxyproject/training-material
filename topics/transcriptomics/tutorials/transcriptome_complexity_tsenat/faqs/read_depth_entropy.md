---
title: Does read depth affect entropy estimates?
area: technical
layout: faq
box_type: tip
contributors: [gallardoalba]
---

Yes. TSENAT uses **automatic pseudocount selection** to handle this:
- Pseudocounts scale with library size
- Ensures entropy estimates are comparable across samples with different sequencing depths
- Breaks ties in rank-based validation methods
- Standard practice in RNA-seq analysis
