---
title: Mismatched Chromosome identifiers and how to avoid them
area: datasets
box_type: tip
layout: faq
contributors: [jennaj, Melkeb]
---


- The methods listed help to identify and correct errors or unexpected results linked to inputs having non-identical chromosome identifiers and/or different chromosome sequence content.

- **If using a Custom Reference genome**, the methods below also apply, but the first step is to make certain that the Custom Genome is formatted correctly. Improper formating is the most common root cause of CG related errors.

Method 1: [Find BAM dataset identifiers]({% link faqs/galaxy/datasets_BAM_dataset_identifiers.md %})

Method 2: [Directly obtaining UCSC sourced *genome* identifiers]({% link faqs/galaxy/datasets_UCSC_sourced_genome.md %})

Method 3: [UCSC sourced data used with other sourced data to adjust identifiers]({% link faqs/galaxy/datasets_adjust_identifiers_UCSC.md %})

Method 4: [Adjusting Identifiers or Input source for any mixed sourced data]({% link faqs/galaxy/datasets_adjusting_identifiers_mixed_data.md %})

Method 5: [A Note on Built-in Reference Genomes]({% link faqs/galaxy/datasets_builtin_genomes.md %})
