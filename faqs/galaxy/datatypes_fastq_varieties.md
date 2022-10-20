---
title: "FASTQ files: `fastq` vs `fastqsanger` vs .."
area: datatypes
layout: faq
box_type: tip
contributors: [shiltemann]
---

FASTQ files come in various flavours. They differ in the encoding scheme they use. See our [QC tutorial]({% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %}#inspect-a-raw-sequence-file) for a more detailed explanation of encoding schemes.

Nowadays, the most commonly used encoding scheme is sanger. In Galaxy, this is the `fastqsanger` datatype. If you are using older datasets, make sure to verify the FASTQ encoding scheme used in your data.

**Be Careful: choosing the wrong encoding scheme can lead to incorrect results!**

**Tip:** There are 2 Galaxy datatypes that have similar names, but are not the same, please make sure you `fastqsanger` and `fastqcssanger` (not the additional `cs`).

**Tip:** When in doubt, choose `fastqsanger`


