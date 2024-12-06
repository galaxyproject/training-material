---
title: "MSStats: what does ‘compare groups = yes’ mean? And the comparison matrix to define the contrast between the 2 groups?"
area: analysis
layout: faq
box_type: question
contributors: [foellmelanie]
---

MSstats consists of three parts:
- Reading the input files and converting them into an MSstats compatible format, doing some processing of the data at the same time
- Data processing: such as protein inference (summary), log2 transformation, normalization and missing value imputation
- `compare groups = yes`, means that the third step is performed, which is statistical analysis: Statistical modelling to find differentially abundant protein between different groups. The groups should be specified as “condition” in the annotation file and the group comparison matrix file specifies which groups to compare against each other. In the example this is quite simple because there are only 2 groups, with 3 or more groups the comparison matrix could become more complex.



