---
title: Adding a tag
description: Tags can help you to better organize your history and track datasets.
area: datasets
layout: faq
box_type: tip
contributors: [bebatut,wm75,hexylena,shiltemann,nekrut]
---

Datasets can be tagged. This simplifies the tracking of datasets across the Galaxy interface. Tags can contain any combination of letters or numbers but cannot contain spaces. 

**To tag a dataset**:
 
1. Click on the dataset to expand it
2. Click on **Add Tags** {% icon galaxy-tags %}
3. Add {% if include.tag %} a tag named `{{include.tag}}`{% else %} tag text{% endif %}. Tags starting with `#` will be automatically propagated to the outputs of tools using this dataset (see below).
4. Press <kbd>Enter</kbd>
5. Check that the tag appears below the dataset name

**Tags beginning with `#` are special!**

They are called **Name tags**. The unique feature of these tags is that they *propagate*: if a dataset is labelled with a name tag, all derivatives (children) of this dataset will automatically inherit this tag (see below).
The figure below explains why this is so useful. Consider the following analysis (numbers in parenthesis correspond to dataset numbers in the figure below): 

1. a set of forward and reverse reads (datasets 1 and 2) is mapped against a reference using {% tool Bowtie2 %} generating dataset 3;
1. dataset 3 is used to calculate read coverage using {% tool BedTools Genome Coverage %} *separately* for `+` and `-` strands. This generates two datasets (4 and 5 for plus and minus, respectively);
1. datasets 4 and 5 are used as inputs to {% tool Macs2 broadCall %} datasets generating datasets 6 and 8;
1. datasets 6 and 8 are intersected with coordinates of genes (dataset 9) using {% tool BedTools Intersect %} generating datasets 10 and 11.

![A history without name tags versus history with name tags]({% link shared/images/histories_why_nametags.svg %})

Now consider that this analysis is done without name tags. This is shown on the left side of the figure. It is hard to trace which datasets contain "plus" data versus "minus" data. For example, does dataset 10 contain "plus" data or "minus" data? Probably "minus" but are you sure? In the case of a small history like the one shown here, it is possible to trace this manually but as the size of a history grows it will become very challenging.

The right side of the figure shows exactly the same analysis, but using name tags. When the analysis was conducted datasets 4 and 5 were tagged with `#plus` and `#minus`, respectively. When they were used as inputs to Macs2 resulting datasets 6 and 8 automatically inherited them and so on... As a result it is straightforward to trace both branches (plus and minus) of this analysis. 

More information is in a [dedicated #nametag tutorial]({% link topics/galaxy-interface/tutorials/name-tags/tutorial.md %}).


<!-- Image is here = https://docs.google.com/drawings/d/1iiNsau6ddiE2MV9qMyekUq2mrpDHHcc02bXtcFEAnhY/edit?usp=sharing -->

