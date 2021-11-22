---
title:  I have RNAseq data for pilot experiment for differential expression in host associated bacterium. One dataset is obtained from bacterial culture, but the other comes from bacteria obtained from the host (plant). I expect strong contamination of the second sample with host RNA reads. Should I filter out reads from the host before performing the analysis (if so, what tools I could use for that), or could I just ignore the contamination (since I will use the bacterial genome to map the reads, it will disregard any host associated reads)?
area: RNA-seq 
box_type: tip
layout: faq
---

You could map both sets of reads to the reference genome. You are right - no pre-filtering required as the host reads shouldn't map to the ref and will be excluded.