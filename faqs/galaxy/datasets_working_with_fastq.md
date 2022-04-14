---
title: Working with very large fastq datasets
area: datasets
layout: faq
box_type: tip
contributors: [jennaj, garimavs]
---

- Run **FastQC** on your data to make sure the format/content is what you expect. Run more QA as needed. 
    - Search [GTN](https://training.galaxyproject.org/) tutorials with the keyword “qa-qc” for examples.
    - Search [Galaxy Help](https://help.galaxyproject.org/) with the keywords “qa-qc” and “fastq” for more help.
- How to create a single smaller input. Search the tool panel with the keyword “subsample” for tool choices.
- How to create multiple smaller inputs. Start with **Split file to dataset collection**, then merge the results back together using a tool specific for the datatype. Example: BAM results? Use **MergeSamFiles**.
