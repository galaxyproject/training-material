---
redirect_from: [/faqs/galaxy/analysis_extended_Extended_help_differential_expression_analysis_tools]
title: Extended Help for Differential Expression Analysis Tools
area: analysis
box_type: tip
layout: faq
contributors: [jennaj, Melkeb]
---

The error and usage help in this FAQ applies to most if not all Bioconductor tools.

- DEseq2
- Limma
- edgeR
- goseq
- Diffbind
- StringTie
- Featurecounts
- HTSeq-count
- HTseq-clip
- Kalisto
- Salmon
- Sailfish
- DEXSeq
- DEXSeq-count

{% icon galaxy-info %}  Review your error messages and you'll find some clues about what may be going wrong and what needs to be adjusted in your rerun. If you are getting a message from `R`, that usually means the underlying tool could not read in or understand your inputs. This can be a labeling problem (what was typed on the form) or a content problem (data within the files).

Expect odd errors or content problems if any of the usage requirements below are not met.

General

- Are your reference genome, reference transcriptome, and reference annotation all based on the same genome assembly?
    * Check the identifiers in all inputs and adjust as needed.
    * These all may mean the same thing to a person but not to a computer or tool: chr1, Chr1, 1, chr1.1
- Differential expression tools all require sample count replicates. [Rationale from two of the DEseq tool authors](https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/26388-deseq2-without-biol-replicates).
    * At least two factor levels/groups/conditions with two samples each.
    * All must all contain unique content for valid scientific results.
- Factor/Factor level names should only contain alphanumeric characters and optionally underscores.
    * Avoid starting these with a number and do not include spaces.
    * Galaxy *may* be able to normalize these values for you, but if you are getting an error: standardize the format yourself.
- **DEXSeq** additionally requires that the first Condition is labeled as `Condition`.
- If your count inputs have a header, the option **Files have header?** is set to **Yes**. If no headers, set to **No**. 
    * If your files have more than one header line: keep the sample header line, remove all extra line(s).
- Make sure that tool form settings match your annotation content or the tool cannot match up the inputs!
    * If you are counting by **gene_id**, your annotation should contain gene_id attributes (9th column)
    * If you are summarizing by **exon**, your annotation should contain exon features (3rd column)

Reference genome (fasta)

- Can be a server reference genome (hosted index in the pull down menu) or a custom reference genome (fasta from the history).
- Custom reference genomes must be [formatted correctly]({% link faqs/galaxy/reference_genomes_custom_genomes.md %}).
- If you are using **Salmon**, you probably don't need a reference genome but a reference transcriptome instead!
- More about understanding and [working with large fasta datasets]({% link faqs/galaxy/datasets_working_with_fasta.md %}).

Reference transcriptome (fasta)

- Fasta file containing assembled transcripts.
- Unassembled short or long reads will not work as a substitute.
- The transcript identifiers on the `>seq` fasta lines must exactly match the `transcript_id` values in your annotation or tabular mapping file.
- Sometimes **Salmon** or **DESeq2** (when comparing TMP values) does not understand `transcript_id.N` (where N is a version number). Try removing `.N` from all inputs. 
  
Reference annotation (tabular, GTF, GFF3)

- Reference annotation [in GTF format]({% link faqs/galaxy/datasets_working_with_reference_annotation.md %}) works best.
- If a GTF dataset is not available for your genome, a two-column tabular dataset containing `transcript <tab> gene` can be used instead with most of these tools. 
- **HTseq-count** requires GTF attributes. Featurecounts is an alternative tool choice.
- Sometimes the tool **gffread** is used to transform GFF3 data to GTF. 
- DO use UCSC's reference annotation (GTF) and reference transcriptome (fasta) data from their [Downloads](https://hgdownload.soe.ucsc.edu/downloads.html) area.
    * These are a match for the UCSC genomes indexed at public Galaxy servers.
    * Links can be directly copy/pasted into the Upload tool.
    * Allow Galaxy to *autodetect the datatype* to produce an uncompressed dataset in your history ready to use with tools.
- Avoid GTF data from the UCSC Table Browser: this leads to scientific problems. GTFs will have the same content populated for both the transcript_id and gene_id values. See the note at UCSC for more about why.
- Still have problems? Try removing all GTF header lines with the tool **Remove beginning of a file**.
- More about understanding and [working with GTF/GFF/GFF3 reference annotation]({% link faqs/galaxy/datasets_working_with_reference_annotation.md %})
