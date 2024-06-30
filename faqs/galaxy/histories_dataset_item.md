---
title: Dataset snippet
description: Describes features of a single dataset element in the history
area: histories
box_type: tip
layout: faq
contributors: [nekrut]
---

A single Galaxy dataset can either be "collapsed" or "expanded".

**Collapsed dataset view**

Datasets in the panel are initially shown in a "collapsed" view:

![Collapsed view of a single Galaxy dataset]({% link shared/images/history_item_collapsed.png %})

It contains the following elements:

- **Dataset number**: ("1") order of dataset in the history;
- **Dataset name**: ("M117-bl_1.fq.gz") its name;
- {% icon galaxy-eye %}: click this to view the dataset contents;
- {% icon galaxy-pencil %}: click this to edit dataset properties;
- {% icon galaxy-delete %}: click this to delete the dataset from the history (*don't worry*, you can undo this action!).

Clicking on a collapsed dataset will expand it.

> <details-title>Some buttons can be disabled.</details-title>
> Some of the buttons above may be disabled if the dataset is in a state that doesn't allow the
> action. For example, the 'edit' button is disabled for datasets that are still queued or running
>
{: .details}

**Expanded dataset view**

Expanded dataset view adds a preview element and many additional controls. 

![Expanded view of a single Galaxy dataset]({% link shared/images/history_item_expanded.png %})

In addition to the elements described above for the collapsed dataset, its expanded view contains:

- **Add tags** {% icon galaxy-tags %}: click on this to tag this dateset;
- **Dataset size**: ("2 variants, 18 comments") lists the size of the dataset. When datasets are small (like in this example) the exact size is shown. For large datasets, Galaxy gives an approximate estimate.
- **format**: ("VCF") lists the datatype;
- **database**: ("?") lists which genome built this dataset corresponds to. This usually lists "?" unless the genome build is set explicitly or the dataset is derived from another dataset with defined genome build information;
- **info field**: ("INFO [2024-03-26 12:08:53,435]...") displays information provided by the tool that generated this dataset. This varies widely and depends on the type of job that generated this dataset.
- {% icon dataset-save %}: Saves dataset to disk;
- {% icon dataset-link %}: Copies dataset link into clipboard;
- {% icon dataset-info %}: Displays additional details about the dataset in the center pane;
- {% icon dataset-rerun %}: Reruns job that generated this dataset. This button is unavailable for datasets uploaded into history because they were not produced by a Galaxy tool;
- {% icon dataset-visualize %}: Displays visualization options for this dataset. The list of options is dependent on the datatype;
- {% icon dataset-related-datasets %}: Shows datasets related to this dataset. This is useful for tracking down parental datasets - those that were used as inputs into a job that produced this particular dataset.
