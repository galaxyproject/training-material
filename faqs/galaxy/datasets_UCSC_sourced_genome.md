---
title: Directly obtaining UCSC sourced *genome* identifiers
area: data upload
box_type: tip
layout: faq
contributors: [jennaj, Melkeb]
---


**Option 1**

1. Go to [UCSC Genome Browser](http://genome.ucsc.edu/), navigate to "genomes", then the species of interest.
2. On the home page for the genome build, immediately under the top navigation box, in the blue bar next to the full genome build name, you will find **View sequences** button.
3. Click on the **View sequences** button and it will take you to a detail page with a table listing out the contents.

**Option 2**

1. Use the tool **Get Data -> UCSC Main**.
2. In the Table Browser, choose the target genome and build.
3. For "group" choose the last option "All Tables".
4. For "table" choose "chromInfo".
5. Leave all other options at default and send the output to Galaxy.
6. This new dataset will load as a tabular dataset into your history.
7. It will list out the contents of the genome build, including the chromosome identifiers (in the first column).
