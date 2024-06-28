---
title: Name tags make histories clear
description: Name tags allow you to track analyses in large histories
area: histories
layout: faq
box_type: tip
contributors: [nekrut]
---

**Name tags** (or "hash tags") are also called "propagating tags": once a dataset is tagged with a name tag all datasets that are derived from it will also carry that tag. Figure below explains why this is so useful. Consider the following analysis (numbers in parenthesis correspond to dataset numbers in the figure below): 

1. a set of forward and reverse reads (datasets 1 and 2) is mapped against a reference using {% tool [Bowtie2](Bowtie2) %} generating dataset 3;
1. dataset 3 is used to calculate read coverage using {% tool [BedTools Genome Coverage](Bedtools Genomecoverage) %} *separately* for `+` and `-` strands. This generates two datasets (4 and 5 for plus and minus, respectively);
1. datasets 4 and 5 are used as inputs to {% tool [Macs2 broadCall](Macs2) %} datasets generating datasets 6 and 8;
1. datasets 6 and 8 are intersected with coordinates of genes (dataset 9) using {% tool [BedTools Intersect](Bedtools Intersect) %} generating datasets 10 and 11.

![A history without name tags versus history with name tags]({% link shared/images/histories_why_nametags.svg %})

Now consider that this analysis is done without name tags. This is shown on the left side of the figure. It is hard to trace which datasets are containing "plus" data versus "minus" data. For example, does dataset 10 contain "plus" data or "minus" data? Probably "minus" but are you sure? In the case of a small history like the one shown here it is possible to trace this manually but as size of a history grows it will become very challenging.

Right side of the figure shows exactly the same analysis, but using name tags. When the analysis was conducted datasets 4 and 5 were tagged with `#plus` and `#minus`, respectively. When they were used as inputs to Macs2 resulting datasets 6 and 8 automatically inherited them and so on... As a result it is straightforward to trace both branches (plus and minus) of this analysis. 

More informatin is in a [dedicated nametag tutorial]({% link topics/galaxy-interface/tutorials/name-tags/tutorial.md %}).

<!-- Image is here = https://docs.google.com/drawings/d/1iiNsau6ddiE2MV9qMyekUq2mrpDHHcc02bXtcFEAnhY/edit?usp=sharing -->