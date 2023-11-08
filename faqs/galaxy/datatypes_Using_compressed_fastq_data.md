---
title: Using compressed fastq data as tool inputs
area: datatypes
box_type: tip
layout: faq
contributors: [jennaj, Melkeb]
---
 

- If the tool accepts `fastq` input, then `.gz` compressed data assigned to the datatype `fastq.gz` is appropriate.
- If the tool accepts `fastqsanger` input, then `.gz` compressed data assigned to the datatype `fastqsanger.gz` is appropriate.
- Using uncompressed `fastq` data is still an option with tools. The choice is yours.

**TIP:** Avoid labeling compressed data with an uncompressed datatype, and the reverse. Jobs using mismatched datatype versus actual format will fail with an error.
