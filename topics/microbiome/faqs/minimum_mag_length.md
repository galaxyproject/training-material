---
title: "Minimum MAG length"
box_type: details
layout: faq
contributors: [bebatut]
---

In **dRep**, the set of **Metagenome-Assembled Genomes (MAGs)** undergoes **quality filtering** before any comparisons are performed. This critical step ensures that only **high-quality genomes** are retained for downstream analysis, improving both **accuracy and computational efficiency**.

One of the first steps in quality filtering is **length-based filtering**. Genomes that do not meet the **minimum length threshold** are **filtered out upfront**. This avoids unnecessary **CheckM computations** and ensures that only sufficiently large genomes are processed.

**Important Note:** All genomes must contain **at least one predicted Open Reading Frame (ORF)**. Without an ORF, **CheckM may stall** during processing. To prevent this, a **minimum genome length of 10,000 base pairs (bp)** is recommended.

The **default minimum length threshold in dRep is 50,000 bp**. This threshold strikes a balance between **computational efficiency** and **genome quality**, ensuring that only meaningful and sufficiently complete genomes are included in the de-replication process.
