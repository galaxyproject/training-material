---
title: Working with GFF GFT GTF2 GFF3 reference annotation
area: datasets
layout: faq
box_type: tip
contributors: [jennaj, garimavs]
---

- All annotation datatypes have a distinct format and content specification.
    - Data providers may release variations of any, and tools may produce variations. 
    - GFF3 data may be labeled as GFF.
    - Content can overlap but is generally not understood by tools that are expecting just one of these specific formats. 
- Best practices
    - The sequence identifiers must exactly match between reference annotation and reference genomes transcriptomes exomes.
    - Most tools expect GFT format unless the tool form specifically notes otherwise.
        - Get the GTF version from the data providers if it is available.
        - If only GFF3 is available, you can attempt to transform it with the tool **gffread**.
    - Was GTF data detected as GFF during Upload? It probably has headers.
        -Remove the headers (lines that start with a "#") with the **Select** tool using the option "NOT Matching" with the regular expression: ^#
        - [Redetect the datatype](https://training.galaxyproject.org/training-material/faqs/galaxy/#detecting-the-datatype-file-format). It should be GTF once corrected.
    - UCSC annotation
        - Find annotation under their Downloads area. The path will be similar to: `https://hgdownload.soe.ucsc.edu/goldenPath/<database>/bigZips/genes/`
        - Copy the URL from UCSC and paste it into the Upload tool, allowing Galaxy to detect the datatype.
