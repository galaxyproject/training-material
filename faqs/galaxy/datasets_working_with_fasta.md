---
title: Working with very large fasta datasets
area: datasets
layout: faq
box_type: tip
contributors: [jennaj, garimavs]
---

- Run **FastQC** on your data to make sure the format/content is what you expect. Run more QA as needed. 
    - Search [GTN](https://training.galaxyproject.org/) tutorials with the keyword “qa-qc” for examples.
    - Search [Galaxy Help](https://help.galaxyproject.org/) with the keywords “qa-qc” and “fasta” for more help.
- Assembly result?
    - Consider filtering by length to remove reads that did not assemble.
    - Formatting criteria:
        - All sequence identifiers must be unique.
        - Some tools will require that there is no description line content, only identifiers, in the fasta title line (“>” line). Use **NormalizeFasta** to remove the description (all content after the first whitespace) and wrap the sequences to 80 bases.
- [Custom genome](https://training.galaxyproject.org/training-material/faqs/galaxy/analysis_add_custom_build.html), transcriptome exome?
    - Only appropriate for smaller genomes (bacterial, viral, most insects).
    - Not appropriate for any mammalian genomes, or some plants/fungi.
    - Sequence identifiers must be an exact match with all other inputs or expect problems. See **GFF GFT GFF3**.  
    - Formatting criteria:
        - All sequence identifiers must be unique.
        - ALL tools will require that there is no description content, only identifiers, in the fasta title line (“>” line). Use **NormalizeFasta** to remove the description (all content after the first whitespace) and wrap the sequences to 80 bases. 
        - The only exception is when executing the **MakeBLASTdb** tool and when the input fasta is in NCBI BLAST format (see the tool form).