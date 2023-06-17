---
title: Can we also use this workflow on Illumina raw reads?
box_type: question
layout: faq
contributors: [EngyNasr]
---

Yes, some tools would need to be changed or removed:

- For the Preprocessing workflow, plotting with Nanoplot shall be removed and keep only FastQC, MultiQC and Fastp.
- For the mapping in the SNP based pathogen detection workflow, instead of Minimap2, Bowtie can be used.
