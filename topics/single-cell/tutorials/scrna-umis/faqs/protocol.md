---
title: How do I know what protocol my data was sequenced with?
layout: faq
box_type: question
contributors: [nomadscientist,mtekman]
redirect_from:
  - /topics/transcriptomics/tutorials/scrna-umis/faqs/protocol
---

If you have 10x data, then you just need to count the length of the R1 reads to guess the Chromium version (see [this tutorial]({% link topics/single-cell/tutorials/scrna-preprocessing-tenx/tutorial.md %})). For other types of data, you must know the protocol in advance, and even then you must also know the multiplexing strategy and the list of expected (whitelisted) barcodes. The whitelist may vary from sequencing lab to sequencing lab, so always ask the wetlab people how the FASTQ data was generated.

