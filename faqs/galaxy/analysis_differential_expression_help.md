---
redirect_from: [/faqs/galaxy/analysis_extended_Extended_help_differential_expression_analysis_tools]
title: Extended Help for Differential Expression Analysis Tools
area: analysis
box_type: tip
layout: faq
contributors: [jennaj, Melkeb]
---

The error and usage help in this FAQ applies to:

- Deseq2
- Limma
- edgeR
- goseq
- DEXSeq
- Diffbind
- StringTie
- Featurecounts
- HTSeq
- Kalisto
- Salmon
- Sailfish
- DexSeq-count

Expect odd errors or content problems if any of the usage requirements below are not met:

- Differential expression tools all require count dataset replicates when used in Galaxy. At least two per factor level and the same number per factor level. These must all contain unique content.
- Factor/Factor level names should only contain alphanumeric characters and optionally underscores. Avoid starting these with a number and do not include spaces.
- If the tool uses `Conditions`, the same naming requirements apply. `DEXSeq` additionally requires that the first Condition is labeled as `Condition`.
- Reference annotation should be in GTF format for most of these tools, with no header/comment lines. Remove all GTF header lines with the tool **Remove beginning of a file**. If any are comment lines are internal to the file, those should be removed. The tool **Select** can be used.
- Make sure that if a GTF dataset is used, and tool form settings are expecting particular attributes, those are actually in your annotation file (example: gene_id).
- GFF3 data (when accepted by a tool) should have single `#` comment line and any others (at the start or internal) that usually start with a `##` should be removed. The tool **Select** can be used.
- If a GTF dataset is not available for your genome, a two-column tabular dataset containing `transcript <tab> gene` can be used instead with most of these tools. Some reformatting of a different annotation file type might be needed. Tools in the groups under **GENERAL TEXT TOOLS** can be used.
- Make sure that if your count inputs have a header, the option **Files have header?** is set to **Yes**. If no header, set to **No**.
- Custom genomes/transcriptomes/exomes must be formatted correctly before mapping.
- Any reference annotation should be an exact match for any genome/transcriptome/exome used for mapping. Build and version matter.
- Avoid using [UCSC's](https://genome.ucsc.edu/) annotation extracted from their Table Browser. All GTF datasets from the UCSC Table Browser have the same content populated for the transcript_id and gene_id values. Both are the "transcript_id", which creates scientific content problems, effectively meaning that the counts will be summarized "by transcript" and not "by gene", even if labeled in a tool's output as being "by gene". It is usually possible to extract gene/transcript in tabular format from other related tables. Review the Table Browser usage at [UCSC](https://genome.ucsc.edu/) for how to link/extract data or ask them for guidance if you need extra help to get this information for a specific data track.

Note: Selected genomes at UCSC do have a reference anotatation GTF pre-computed and available with a Gene Symbol populated into the "gene_id" value. Find these in the UCSC "Downloads" area. When available, the link can be directly copy/pasted into the Upload tool in Galaxy. Allow Galaxy to *autodetect the datatype* to produce an uncompressed GTF dataset in your history ready to use with tools.
