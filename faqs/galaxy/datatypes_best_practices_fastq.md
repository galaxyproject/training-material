---
title: Best practices for loading fastq data into Galaxy
area: datatypes
layout: faq
box_type: tip
contributors: [jennaj, Melkeb]
---

- As of release `17.09`, `fastq` data will have the datatype `fastqsanger` auto-detected when that quality score scaling is detected and "autodetect" is used within the Upload tool. Compressed `fastq` data will be converted to uncompressed in the history.
- To preserve `fastq` compression, directly assign the appropriate datatype (eg: `fastqsanger.gz`).
- If the data is close to or over 2 GB in size, be sure to use FTP.
- If the data was already loaded as `fastq.gz`, don't worry! Just test the data for correct format (as needed) and assign the [metadata type]({% link faqs/galaxy/datatypes_fastq_and_fastqsanger.md %}).
