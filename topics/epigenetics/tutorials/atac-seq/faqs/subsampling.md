---
title: Is it possible to subsample some samples if you have more reads?
box_type: question
layout: faq
contributors: [lldelisle]
---

Yes, we would recommend to process all reads and just before the peak calling. You can use {% icon tool %} **Samtools view** to sample the BAM file.

