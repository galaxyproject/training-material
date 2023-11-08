---
title: Data retrieval with “NCBI SRA Tools” (fastq-dump)
area: data upload
box_type: tip
layout: faq
contributors: [jennaj, Melkeb]
---

This section will guide you through downloading experimental metadata, organizing the metadata to short lists corresponding to conditions and replicates, and finally importing the data from NCBI SRA in collections reflecting the experimental design.

**Downloading metadata**

- It is critical to understand the condition/replicate structure of an experiment before working with the data so that it can be imported as collections ready for analysis. Direct your browser to [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/) and in the search box enter GEO data set identifier (for example: GSE72018). Once the study appears, click the box to download the “RunInfo Table”.

**Organizing metadata**

- The “RunInfo Table” provides the experimental condition and replicate structure of all of the samples. Prior to importing the data, we need to parse this file into individual files that contain the sample IDs of the replicates in each condition. This can be achieved by using a combination of the ‘group’, ‘compare two datasets’, ‘filter’, and ‘cut’ tools **to end up with single column lists of sample IDs (SRRxxxxx) corresponding to each condition.**

**Importing data**

- Provide the files with SRR IDs to NCBI SRA Tools (fastq-dump) to import the data from SRA to Galaxy. By organizing the replicates of each condition in separate lists, **the data will be imported as “collections” that can be directly loaded to a workflow or analysis pipeline.**
