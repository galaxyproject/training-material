---
title: Sorting Reference Genome
area: reference genome 
box_type: tip        
layout: faq        
contributors: [jennaj, nurzhamalyrys] 
---

Certain tools expect that reference genomes are sorted in [lexicographical order](https://en.wikipedia.org/wiki/Lexicographic_order). These tools are often downstream of the initial mapping tools, which means that a large investment in a project has already been made (i.e. a long mapping process), before a problem with sorting pops up in conclusion layer tools. No one likes to start over!

How to avoid? Always sort your FASTA reference genome dataset at the beginning of a project. Many sources only provide sorted genomes, but double checking is your own responsibility, and super easy in Galaxy. So easy that there isn't even a shared workflow, just a recipe (but feel free to make your own):

Quick lexicographical sort recipe:
1. Convert **Formats** -> **FASTA-to-Tabular**
2. Filter and Sort -> Sort
       on column: c1 
       with flavor: Alphabetical
       everything in: Ascending order
3. Convert Formats -> Tabular-to-FASTA

Note: The above sorting method is for most tools, but not all. In particular, GATK tools have a tool-specific sort order requirement. [The Broad Institute](https://www.broadinstitute.org/) has FAQ with input instructions. 
