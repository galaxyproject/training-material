---
title: I am using the RNA STAR tool to map my reads on the reference genome. while filling in the details for , I have to fill in "Length of the genomic sequence around annotated junctions", which apparently has to be 36. I'm lost for a moment why this is 36, what is meant by it, and why is it relevant? Does anyone have any ideas? Why should it be the length of the reads -1?
area: mapping
box_type: tip
layout: faq
---

RNA STAR is using the gene model to create the database of splice junctions, and that these don't "need" to have a length longer than the reads (37bp).
