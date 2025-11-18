---
title: Merge Reads for Co/Grouped Assembly
description: How to merge reads from multiple samples for co- or group-assembly
area: analysis
box_type: tip
layout: faq
contributors: [paulzierep]
---

**Co-assembly** is the process of assembling reads from multiple metagenomic samples together into contigs or genomes (MAGs), instead of assembling each sample individually. If only subsets of samples are assembled together, this can also be referred to as **grouped assembly**.

In principle, co- or grouped assembly can be performed with any assembler by merging reads from the desired samples before running the assembly.

To simplify this process, we provide the {% tool [Fastq GroupMerge](toolshed.g2.bx.psu.edu/repos/iuc/fastq_groupmerge/fastq_groupmerge/1.0.1+galaxy0) %} tool. This tool allows flexible grouping (including multiple groupings) using a simple tabular file, as shown below:

| sample  | group |
| ------- | ----- |
| reads_A | 1     |
| reads_B | 1     |
| reads_B | 2     |
| reads_C | 3     |

**Explanation of the example:**

* Reads from **A** and **B** will be merged into a new synthetic sample called **1**.
* Reads from **B** will also be kept individually as **sample 2**.
* Reads from **C** will remain as **sample 3**.

The tool supports merging of **both single-end and paired-end reads**, making it flexible for various sequencing datasets.

