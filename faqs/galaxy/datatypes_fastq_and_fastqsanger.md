---
title: How do `fastq.gz` datasets relate to the `.fastqsanger` datatype metadata assignment?
area: datatypes
layout: faq
box_type: tip
contributors: [jennaj, Melkeb]
---

Before assigning `fastqsanger` or `fastqsanger.gz`, be sure to confirm the format.

**TIP:** 
- Using *non-fastqsanger* scaled quality values will cause scientific problems with tools that expected `fastqsanger` formatted input.
- Even if the tool does not fail, get the format right from the start to avoid problems. Incorrect format is still one of the most common reasons for tool errors or unexpected results (within Galaxy or not).
- For more information on [How to format fastq data for tools that require .fastqsanger format?](https://training.galaxyproject.org/training-material/faqs/galaxy/datatypes_fastq_format.html)
