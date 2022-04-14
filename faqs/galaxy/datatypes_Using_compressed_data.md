---
title: Using compressed data as tool inputs
area: 

- If the tool accepts fastq input, then gz compressed data assigned the datatype `fastq.gz` is appropriate.
- If the tool accepts fastqsanger input, then gz compressed data assigned the datatype `fastqsanger.gz` is appropriate.
- Using uncompress fastq data is still an option with tools. The choice is yours.

**TIP** Avoid labeling compressed data with an uncompressed datatype, and the reverse. Jobs using mismatched datatype versus actual format will fail with an error.
