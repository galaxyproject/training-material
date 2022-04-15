---
title: Adjusting identifiers for UCSC sourced data used with other sourced data
area: datasets
box_type: tip
layout: faq
contributors: [jennaj, Melkeb]
---


- A GTF formatted dataset (potentially a "reference annotation dataset"), with Ensembl/UCSC/Other based chromosome identifiers, is to be used with a tool that also makes use of a different sourced reference genome or the reverse may be true, Ensembl/UCSC/Other sourced reference genome and a differnt sourced reference annotation.

- The underlying genome sequence content is otherwise identical. If not, see [Adjusting Identifiers or Input source for any mixed sourced data]({% link faqs/galaxy/datasets_adjusting_identifiers_mixed_data.md %})

- To adjust the Ensembl/Other reference annotation to match a UCSC-sourced reference genome (or another source that uses UCSC-style chromosome names), add a "chr" to the chromosome name, so that "N" becomes "chrN". Using tools from the group "Text Manipulation".

For **bed** data:

1. Tool **Add column**: add "chr" to the original dataset as a new column.
2. Tool **Merge Columns**: merge "c7" with "c1".
3. Tool **Cut**: cut "c8,c2,c3,c4,c5,c6" (replace c1 & c7 - with merged c8 - the new chrom identifier).
4. Click on the pencil icon for the result dataset, then the tab for "Datatype". Assign "bed" and save. Allow the metadata to complete assignment (the "yellow" dataset state).
5. Now click on the tab for "Attributes" and assign the remaining columns. Strand = 6, name = 4, and score = 5. Save. For best results with certain downstream tools, allow the metadata to complete assignment.

For **wig/wiggle** data (NOT compressed bigWig):

1. Tool **Replace parts of text**.
2. File to process: Use Multi-select select wig datasets to fix (one or more).
3. Find pattern: chrom=.
4. Replace with: chrom=chr.
5. Remainder of options left at default.

