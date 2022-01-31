---
title:  Regarding DESeq2, in the tutorial you used the normalised count table. Some people use VST normalised counts or rlog normalised counts for visualisation (heatmaps),  would you recommend it ? And second question, regarding the heatmap2, I think this depends on the data you analyse but do you have any advise on how to select the clustering method and the distance method?
area: DESeq2
box_type: tip
layout: faq
---

this depends on what you would like to do with the table. The DESeq2 wrapper in Galaxy can output all of these, and there is a nice discussion in the [DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-data-transformations) about this topic.