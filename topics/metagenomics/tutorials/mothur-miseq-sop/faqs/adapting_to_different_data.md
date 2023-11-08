---
title: How can I adapt this tutorial to my own data?
area: analysis
layout: faq
box_type: question
contributors: [shiltemann]
---

If you would like to run this analysis on your own data, make sure to check which V-region was sequenced. In this tutorial, we sequenced the V4 region, and used a corresponding reference for just this region. If you sequenced another V-region, please use an appropriate reference (either the full SILVA reference, or the SILVA reference specific for your region). Similarly, the `Screen.seqs` step after the alignment filtered on start and end coordinates of the alignments. These will have to be adjusted to your V-region.
