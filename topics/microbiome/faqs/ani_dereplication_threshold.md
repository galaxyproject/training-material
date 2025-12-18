---
title: "ANI threshold for dereplication"
box_type: details
layout: faq
contributors: [bebatut]
---

Selecting an **appropriate ANI threshold** for de-replication depends on the specific goals of your analysis. ANI measures the similarity between genomes, but the **threshold for considering genomes "the same"** varies by application. Two key parameters influence this decision:
- **Minimum secondary ANI:** The lowest ANI value at which genomes are considered identical.
- **Minimum aligned fraction:** The minimum genome overlap required to trust the ANI calculation.

These parameters are determined by the **secondary clustering algorithm** used in tools like {% tool [dRep](toolshed.g2.bx.psu.edu/repos/iuc/drep_dereplicate/drep_dereplicate/3.6.2+galaxy1) %}

**Common Use Cases and Thresholds**

1. **Species-Level De-Replication**

    For generating **Species-Level Representative Genomes (SRGs)**, a **95% ANI threshold** is widely accepted. This threshold was used in studies like {% cite almeida2019new %} to create a comprehensive set of high-quality reference genomes for human gut microbes. While there is debate about whether bacterial species exist as discrete units or a continuum, **95% ANI** remains a practical standard for species-level comparisons. For context, see {% cite olm2017drep %}.

2. **Avoiding Read Mis-Mapping**

    If the goal is to **prevent mis-mapping of short metagenomic reads** (typically \(( \sim \)) 150 bp), a stricter **98% ANI threshold** is recommended. At this level, reads are less likely to map equally well to multiple genomes, ensuring clearer distinctions between closely related strains.

**Default and Practical Considerations**
- **Default ANI threshold in dRep:** **95%** (suitable for most species-level applications).
- **Upper limit for detection:** Thresholds up to **99.9%** are generally reliable, but higher values (e.g., 100%) are impractical due to algorithmic limitations (e.g., genomic repeats). Self-comparisons typically yield \(( \sim \)) 99.99% ANI due to these constraints.
- For **strain-level comparisons**, consider using {% tool [InStrain](toolshed.g2.bx.psu.edu/repos/iuc/instrain_profile/instrain_profile/1.5.3+galaxy0) %} ({% cite olm2021instrain %}), which provides detailed strain-resolution analyses.


**Important Notes**
- De-replication **collapses closely related but non-identical strains/genomes**. This is inherent to the process.
- To explore strain-level diversity after de-replication, map original reads back to the dereplicated genomes and visualize the **strain cloud** using tools like **InStrain** ([Bendall et al. 2016](#)).

For further discussion on ANI thresholds, see [Are these microbes the “same”? (Blog Post)](https://microbe.net/2017/02/15/are-these-microbes-the-same/)