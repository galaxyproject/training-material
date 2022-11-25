---
layout: tutorial_hands_on
subtopic: deconvo
priority: 5
title: Benchmarking MuSiC Deconvolution
zenodo_link: https://zenodo.org/record/NUMBERHERE
tags:
  - single-cell
  - human
  - deconvolution
  - bulk
  - transcriptomics
questions:
- How well does my MuSiC deconvolution work? How accurate is it?
objectives:
- Generate pseudobulk expression matrices, wherein the transcript counts from known single cells are summed up
- Apply deconvolution to these pseudobulks to infer cell proportions
- Compare these inferred proportions to the actual cell proportions used to make the pseudobulks
- Repeat this at scale
time_estimation: 1H
key_points:
- Deconvolution might be better on some tissues than others, so it's a good idea to test your reference accuracy
contributors:
- nomadscientist
- mtekman
requirements:
  -
    type: "internal"
    topic_name: single-cell
    tutorials:
      - bulk-music
      - bulk-music-2-preparescref
      - bulk-music-3-preparebulk

gitter: Galaxy-Training-Network/galaxy-single-cell
---

# Introduction

So you want to do some deconvolution, but how accurate is it likely to be? Is your reference good enough? Is the MuSiC {% cite wang2019bulk %} algorithm itself good enough? Well, we can't answer these questions specifically, but we can give you a good indication of how well the deconvolution is working in a given situation.

To do this, we will generate pseudobulk matrices from our single cell reference, wherein we sum all the transcripts in a random sample of genes to create a 'bulk' sample. This strategy, and some of our downstream statistics (specifically the use of RSME, the root square mean error), parallel current benchmarking strategies in deconvolution {% cite Cobos2020 %}

![3 Panels: Panel 1 shows Input and images of people and tissues and the label "Human pancreas; Panel 2 shows construction of two matrices, one for single cells and one combining cell expression into a column; Panel 3 shows a chart of MAPE and RMSE statistical outputs as well as a heatmap comparing actual and inferred samples"](../../single-cell/images/bulk-music-benchmark/Overall_benchmark.png "Overall strategy for tutorial benchmarking (Adapted from Cobos et al) ")

> <comment-title></comment-title>
> - We will not be running tool by tool in this tutorial. The reason is, none of the tools are new for this tutorial, and mostly involve a lot of data manipulation (you can do a tutorial on that [here](https://training.galaxyproject.org/training-material/topics/introduction/tutorials/data-manipulation-olympics/tutorial.html) if you like!)
> - Therefore, instead, we've gone through the hard work of generating reusable workflows that you can simply run. We'll focus on data inputs, outputs, and interpretation, rather than running tool by tool.
{: .comment}

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data

We previously retrieved this pancreas single cell dataset {% cite Segerstolpe2016 %} from the [single cell expression atlas](https://www.ebi.ac.uk/gxa/sc/experiments) and manipulated it into an Expression Set object. As a reminder, these cells come from 6 healthy individuals and 4 individuals with Type II diabetes, and we combined all 10 samples into a single Expression Set object. To do this, we made a metadata file with details about each cell (including cell type), and a cell x gene matrix with cells as rows and genes as columns. We'll need all these files to generate our pseudobulks and test them against their reference.

# 1: Generating pseudobulks

Our first goal is to make pseudobulk matrices, with a list of gene expressions for a single 'pseudobulk' sample that we generate from the expression of known single cells. We also need to generate a list of those known single cells, which we'll use later to see how good our deconvolution was. To do this, we need two files:

 - Metadata file (which label each cell barcode with an author-declared cell type)
 - Matrix file (which has cells as rows and genes as columns)

## Get data FIXME

Input: https://usegalaxy.eu/u/wendi.bacon.training/h/deconvolution-benchmarking--1---scrna-matrix-to-pseudobulk-input

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial *"Deconvolution: Benchmarking"*
> 2. Import the files from [Zenodo]({{ page.zenodo_link }})
>
>    * Human single cell RNA ESet objects (tag: `#singlecell`)
>
>      ```
>    {{ page.zenodo_link }}/files/sc_metadata.FIXME
>    {{ page.zenodo_link }}/files/sc_matrix.FIXME
>    {{ page.zenodo_link }}/files/
>      ```
>
>    * Human bulk RNA ESet objects (tag: `#bulk`)
>      ```
>    {{ page.zenodo_link }}/files/ESet_object_bulk_healthy.rdata
>    {{ page.zenodo_link }}/files/ESet_object_bulk_T2D.rdata
>      ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Rename the datasets
>
> 5. Add to each file a tag corresponding to `#bulk` and `#scrna`
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

Let's look at the datasets.

Inspect {% icon galaxy-eye %} the `Metadata` file in your Galaxy history.

> <question-title></question-title>
>
> 1. How many cells are in your sample?
>
> > <div id="solution-1" class="box-title"><button type="button" aria-controls="solution-1-contents" aria-expanded="true" aria-label="Toggle solution box: "><i class="far fa-eye" aria-hidden="true"></i><span class="visually-hidden"></span> Solution<span role="button" class="fold-unfold fa fa-minus-square"></span></button></div>
> >
> > 1. If you look at the dataset in the Galaxy window, you'll find that it has 2915 lines. One of those is a header, so you have 2914 cells. Each row has a cell barcode identifying it as a cell.
> >
> {: .solution}
{: .question}

Now check the `Matrix` file in your Galaxy history - do not inspect {% icon galaxy-eye %} it as it is a massive file!

> <question-title></question-title>
>
> 1. How are the genes labelled in your matrix, by gene symbol or Ensembl ID?
>
> > <div id="solution-1" class="box-title"><button type="button" aria-controls="solution-1-contents" aria-expanded="true" aria-label="Toggle solution box: "><i class="far fa-eye" aria-hidden="true"></i><span class="visually-hidden"></span> Solution<span role="button" class="fold-unfold fa fa-minus-square"></span></button></div>
> >
> > 1. If you look at the dataset in the Galaxy window, you'll find that the rows correspond to the metadata file as cell barcodes, while the columns contain long ENSG##### names. These are Ensembl IDs and will have to be changed in the workflow to gene symbols.
> >
> {: .solution}
{: .question}

Ultimately, we'll use both these files to generate our pseudbulks.

![2 Panels: Panel 1 shows a metadata table with each cell as a row and information such as individual, sex, and cell type as columns. Panel 2 shows a matrix of cell barcodes as rows and genes as columns. Arrows connect cell barcodes between both panels with either cell types or the expression matrix."](../../single-cell/images/bulk-music-benchmark/retrieving_actuals.png "Part 1: Using cell metadata to generate known cell proportions within a pseudobulk sample)

## Get the workflow - Part 1

> <hands-on-title>Generate pseudbulk populations using workflow</hands-on-title>
>
> 1. Import this [workflow](https://usegalaxy.eu/u/wendi.bacon.training/w/deconvolution-benchmarking--1---scrna-matrix-to-pseudobulks).
>
>    {% snippet faqs/galaxy/workflows_import.md %}
>
> 2. Run the workflow on your sample with the following parameters:
>    - {% icon param-file %}*"Formatted metadata"*:  `sc_metadata`
>    - {% icon param-file %}*"Formatted metadata"*:  `sc_matrix`
>    - *"Cell Sample Expr"*: `table.sample(100)` (This makes pseudobulks of 100 cells, although you could change this to suit your needs)
>    - *"Label the second sample"*: `BEGIN { print "column\tB_infer" { print $0 }` (This labels the second sample B}
>    - *"Organism"*: `Hs` (This stipulates your organism for the Ensembl ID collapsing.)
>    - *"Label the first sample"*: `BEGIN { print "column\tA_infer" { print $0 }` (This labels the first sample A}
>
>    {% snippet faqs/galaxy/workflows_run.md %}
>
{: .hands_on}

Workflow:


## Run the workflow - Part 1

FAQ here
Inputs set-up

## Inspect outputs - Part 1

> <warning-title>Danger: Your samples will look different!</warning-title>
> Because you're randomly generating pseudobulk samples, your data will look different from what is in this tutorial. We'll catch this up in Part 2 by giving you our inputs, but you're also welcome to have a go with your own datasets with slightly altered numbers.
{: .warning}

Answer key:
https://usegalaxy.eu/u/wendi.bacon.training/h/deconvolution-benchmarking--1---scrna-matrix-to-pseudobulk-answer-key



# 2: Inferring cellular proportions

## Get data - Part 2

# 3: Statistical analysis

## Get data - Part 3

# 4: Doing this... at scale!

## Get data - Part 4

# Conclusion

## Final histories for each part

Workflow 1
 -



BELOW IS ALL OLD

# Introduction

The goal of this tutorial is to apply bulk RNA deconvolution techniques to a problem with multiple variables - in this case, a model of diabetes is compared with its healthy counterparts. All you need to compare inferred cell compositions are well-annotated, high quality reference scRNA-seq datasets, transformed into MuSiC-friendly Expression Set objects, and your bulk RNA-samples of choice (also transformed into MuSiC-friendly Expression Set objects). For more information on how MuSiC works, you can check out their github site [MuSiC](https://xuranw.github.io/MuSiC/articles/MuSiC.html) or published article {% cite wang2019bulk %}





# Data

In the standard MuSiC tutorial, we used human pancreas data. We will now use the same single cell reference dataset {%cite segerstolpe2016single %} withits 10 samples of 6 healthy subjects and 4 with Type-II diabetes (T2D), as well as the bulk RNA-samples from the same lab (3 healthy, 4 diseased). Both of these datasets were accessed from the public EMBL-EBI repositories and transformed into Expression Set objects in the previous two tutorials. For both the single cell reference and the bulk samples of interest, you have generated Expression Set objects with only T2D samples, only healthy samples, and a final everything-combined sample for the scRNA reference. The plan is to analyse this data in three ways: using a combined reference (altogether); using only the healthy single cell reference (healthyscref); or using a healthy and combined reference separately (like4like), all to identify differences in cellular composition.

![Three colours of arrows connect bulk healthy & diseased data sets to a combined single cell (altogether); bulk healthy and single cell healthy & bulk diseased with single cell diseased (like4like); and bulk diseased and healthy with the single cell healthy reference (healthyscref).](../../images/bulk-music/comparison.png "Plan of analysis")

## Get data

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial *"Deconvolution: Compare"*
> 2. Import the files from [Zenodo]({{ page.zenodo_link }})
>
>    * Human single cell RNA ESet objects (tag: `#singlecell`)
>
>      ```
>    {{ page.zenodo_link }}/files/ESet_object_sc_combined.rdata
>    {{ page.zenodo_link }}/files/ESet_object_sc_T2D.rdata
>    {{ page.zenodo_link }}/files/ESet_object_sc_healthy.rdata
>      ```
>
>    * Human bulk RNA ESet objects (tag: `#bulk`)
>      ```
>    {{ page.zenodo_link }}/files/ESet_object_bulk_healthy.rdata
>    {{ page.zenodo_link }}/files/ESet_object_bulk_T2D.rdata
>      ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Rename the datasets
>
> 5. Add to each file a tag corresponding to `#bulk` and `#scrna`
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

# Infer cellular composition & compare

It's finally time!

## Altogether: Deconvolution with a combined sc reference

> <hands-on-title> Comparing: altogether </hands-on-title>
>
> 1. {% tool [MuSiC Compare](toolshed.g2.bx.psu.edu/repos/bgruening/music_compare/music_compare/0.1.1+galaxy4) %} with the following parameters:
>    - In *"New scRNA Group"*:
>        - {% icon param-repeat %} *"Insert New scRNA Group"*
>            - *"Name of scRNA Dataset"*: `scRNA_set`
>            - {% icon param-file %} *"scRNA Dataset"*: `ESet_object_sc_combined.rdata` (Input dataset)
>            - In *"Advanced scRNA Parameters"*:
>                - *"Cell Types Label from scRNA dataset"*: `Inferred cell type - author labels`
>                - *"Samples Identifier from scRNA dataset"*: `Individual`
>                - *"Comma list of cell types to use from scRNA dataset"*: `alpha cell,beta cell,delta cell,gamma cell,acinar cell,ductal cell`
>            - In *"Bulk Datasets in scRNA Group"*:
>                - {% icon param-repeat %} *"Insert Bulk Datasets in scRNA Group"*
>                    - *"Name of Bulk Dataset"*: `Bulk_set:Normal`
>                    - {% icon param-file %} *"Bulk RNA Dataset"*: `ESet_object_bulk_healthy.rdata` (Input dataset)
>                    - *"Factor Name"*: `Disease`
>                - {% icon param-repeat %} *"Insert Bulk Datasets in scRNA Group"*
>                    - *"Name of Bulk Dataset"*: `Bulk_set:T2D`
>                    - {% icon param-file %} *"Bulk RNA Dataset"*: `ESet_object_bulk_T2D.rdata` (Input dataset)
>                    - *"Factor Name"*: `Disease`
>
> 2. To each of the outputs, add the `#altogether` tag.
>
{: .hands_on}

There are four sets of output files.
1. *Summarised Plots* <- This is the most interesting output, because it has the pretty pictures!
2. *Individual Heatmaps* <- This kind of does what standard (non-Comparing) MuSiC does for each sample, rather than combining them.
3. *Stats* <- This will be very handy if you want to make any statistical calculations, as it contains medians and quartiles
4. *Tables* <- This contains the cell proportions found within each sample as well as the number of reads.

### Summarised Plots

Examine {% icon galaxy-eye %} the output file *Summarised Plots (MuSiC)*. Now the first few pages are similar to the standard deconvolution tool, but now comparing across the factor of interest (disease). Among the myriad of visualisations available, our favourite is on page 5 - a comparison of inferred cell proportions across disease.

![Graph showing two rows for each cell type (gamma, ductal, delta, beta, alpha, and acinar cells) comparison normal or T2D proportions by either read or by sample.](../../images/bulk-music/plot_altogether.png "Altogether: Resulting Graph")

Here we can see that the bulk-RNA seq samples from the T2D patients contain markedly fewer beta cells as compared with their healthy counterparts. This makes sense, so that's good!

### Individual Heatmaps

Examine {% icon galaxy-eye %} the output file *Individual heatmaps (MuSiC)*. This shows the cell distribution across each of the individual samples, separated out by disease factor into two separate plots, but ultimately isn't particularly informative.

![Heatmap showing vertically four samples (Bulk) and horizontally each cell type (from alpha to ductal), with gradations of colour referring to the proportion of that cell type within each sample.](../../images/bulk-music/individual_heatmaps.png "Altogether: Individual heatmaps")

### Stats

If you select the *Stats* dataset, you'll find it contains four sets of data, `Bulk_disease: Read Props`, `Bulk_disease: Sample Props`, `Bulk_healthy: Read Props` and `Bulk_healthy: Sample Props`. Examine {% icon galaxy-eye %} the file `Bulk_disease: Sample Props`. This contains summary statistics (Min, quartiles, median, mean, etc.) for each phenotype. This could be quite helpful if you're trying to statistically identify differences across samples.

![Table with rows as cell types and columns as statistics: Min, 1st Qu., Median, Mean.](../../images/bulk-music/summary_statistics.png "Altogether: Stats")

### Tables

Finally, if you select the *Tables* dataset, you'll find it contains three sets of data, `Data Table`, `Matrix of Cell Type Read Counts`, and `Matrix of Cell Type Sample Proportions`.

Examine {% icon galaxy-eye %} the file `Data Table`. This contains the inferred proportions and reads associated with each sample and cell type, alongwith its important factor of interest (Disease). In this tutorial, we tend to use sample proportions rather than read count, but either works. The two other matrix files are just portions of this data table.

![Table with individual rows for each cell type within each sample, with columns of Cell, Factor, CT Prop in Sample and Number of Reads.](../../images/bulk-music/data_table.png "Altogether: Data table")

> <question-title></question-title>
>
> 1. Why does the data table contain 42 rows?
>
> > <solution-title></solution-title>
> >
> > 1. The data table contains a row for each cell type within each sample. Since there are 6 cell types and 7 samples, `6*7 = 42` rows.
> >
> {: .solution}
>
{: .question}

Hopefully, this has been illumating! Now let's try two other ways of inferring from a reference and see if it makes a difference.

## Like4like: Deconvolution of healthy samples with a healthy reference and diseased samples with a diseased reference

> <hands-on-title> Like4like Inferrence </hands-on-title>
>
> 1. {% tool [MuSiC Compare](toolshed.g2.bx.psu.edu/repos/bgruening/music_compare/music_compare/0.1.1+galaxy4) %} with the following parameters:
>    - In *"New scRNA Group"*:
>        - {% icon param-repeat %} *"Insert New scRNA Group"*
>            - *"Name of scRNA Dataset"*: `scRNA_set:Normal`
>            - {% icon param-file %} *"scRNA Dataset"*: `ESet_object_sc_healthy.rdata` (Input dataset)
>            - In *"Advanced scRNA Parameters"*:
>                - *"Cell Types Label from scRNA dataset"*: `Inferred cell type - author labels`
>                - *"Samples Identifier from scRNA dataset"*: `Individual`
>                - *"Comma list of cell types to use from scRNA dataset"*: `alpha cell,beta cell,delta cell,gamma cell,acinar cell,ductal cell`
>            - In *"Bulk Datasets in scRNA Group"*:
>                - {% icon param-repeat %} *"Insert Bulk Datasets in scRNA Group"*
>                    - *"Name of Bulk Dataset"*: `Bulk_set:Normal`
>                    - {% icon param-file %} *"Bulk RNA Dataset"*: `ESet_object_bulk_healthy.rdata` (Input dataset)
>                    - *"Factor Name"*: `Disease`
>        - {% icon param-repeat %} *"Insert New scRNA Group"*
>            - *"Name of scRNA Dataset"*: `scRNA_set:T2D`
>            - {% icon param-file %} *"scRNA Dataset"*: `ESet_object_sc_T2D.rdata` (Input dataset)
>            - In *"Advanced scRNA Parameters"*:
>                - *"Cell Types Label from scRNA dataset"*: `Inferred cell type - author labels`
>                - *"Samples Identifier from scRNA dataset"*: `Individual`
>                - *"Comma list of cell types to use from scRNA dataset"*: `alpha cell,beta cell,delta cell,gamma cell,acinar cell,ductal cell`
>            - In *"Bulk Datasets in scRNA Group"*:
>                - {% icon param-repeat %} *"Insert Bulk Datasets in scRNA Group"*
>                    - *"Name of Bulk Dataset"*: `bulk_set:T2D`
>                    - {% icon param-file %} *"Bulk RNA Dataset"*: `ESet_object_bulk_T2D.rdata` (Input dataset)
>                    - *"Factor Name"*: `Disease`
>
> 2. Add the `#like4like` tag to each of the outputs.
>
{: .hands_on}

> <question-title></question-title>
>
> 1. How have the cell inferrences changed, now that we have changed the scRNA references used?
>
> > <solution-title></solution-title>
> >
> > ![Two graphs showing two rows for each cell type (gamma, ductal, delta, beta, alpha, and acinar cells) comparison normal or T2D proportions by either read or by sample, with the top graph labelled #altogether and the bottom labelled #like4like. Differences are more pronounced in the top #altogether graph.](../../images/bulk-music/plot_altogether.png "Altogether vs like4like")
> > 1. Overall, our interpretation here is that the differences are less pronounced. It's interesting to conjecture whether this is an artefact of analysis, or whether - possibly - the beta cells in the diseased samples are not only fewer, but also contain fewer beta-cell specific transcripts (and thereby inhibited beta cell function), thereby lowering the bar for the inferrence of a beta cell and leading to a higher proportion of interred B-cells.
> >
> {: .solution}
>
{: .question}

Let's try one more inferrence - this time, we'll use only healthy cells as a reference, to (theoretically) make a more consistent analysis across the two phenotypes.

## healthyscref: Deconvolution using only healthy cells as a reference

> <hands-on-title> Healthy sc reference only inferrence </hands-on-title>
>
> 1. {% tool [MuSiC Compare](toolshed.g2.bx.psu.edu/repos/bgruening/music_compare/music_compare/0.1.1+galaxy4) %} with the following parameters:
>    - In *"New scRNA Group"*:
>        - {% icon param-repeat %} *"Insert New scRNA Group"*
>            - *"Name of scRNA Dataset"*: `scRNA_set:Normal`
>            - {% icon param-file %} *"scRNA Dataset"*: `ESet_object_sc_healthy.rdata` (Input dataset)
>            - In *"Advanced scRNA Parameters"*:
>                - *"Cell Types Label from scRNA dataset"*: `Inferred cell type - author labels`
>                - *"Samples Identifier from scRNA dataset"*: `Individual`
>                - *"Comma list of cell types to use from scRNA dataset"*: `alpha cell,beta cell,delta cell,gamma cell,acinar cell,ductal cell`
>            - In *"Bulk Datasets in scRNA Group"*:
>                - {% icon param-repeat %} *"Insert Bulk Datasets in scRNA Group"*
>                    - *"Name of Bulk Dataset"*: `Bulk_set:Normal`
>                    - {% icon param-file %} *"Bulk RNA Dataset"*: `ESet_object_bulk_healthy.rdata` (Input dataset)
>                    - *"Factor Name"*: `Disease`
>                - {% icon param-repeat %} *"Insert Bulk Datasets in scRNA Group"*
>                    - *"Name of Bulk Dataset"*: `Bulk_set:T2D`
>                    - {% icon param-file %} *"Bulk RNA Dataset"*: `ESet_object_bulk_T2D.rdata` (Input dataset)
>                    - *"Factor Name"*: `Disease`
>
> 2. Add the `#healthyscref` tag to each of the outputs.
>
{: .hands_on}

> <question-title></question-title>
>
> 1. How have the cell inferrences changed this time?
>
> > <solution-title></solution-title>
> >
> > > > ![Three graphs showing two rows for each cell type (gamma, ductal, delta, beta, alpha, and acinar cells) comparison normal or T2D proportions by either read or by sample, with the top graph labelled #altogether; the middle labelled #like4like; and the bottom labelled #healthyscref. Differences are most pronounced in the bottom #healthyscref graph.](../../images/bulk-music/compare_3.png "The impact of the single cell reference")
> >
> > 1. If using a like4like inferrence reduced the difference between the phenotype, aligning both phenotypes to the same (healthy) reference exacerbated them - there are even fewer beta cells in the output of this analysis.
> >
> {: .solution}
>
{: .question}

Overall, it's important to remember how the inferrence changes depending on the reference used - for example, a combined reference might have majority healthy samples or diseased samples, so that would impact the inferred cellular compositions.

# Conclusion

Congrats! You've made it to the end of this suite of deconvolution tutorials! You've learned how to find quality data for reference and for analysis, how to reformat it for deconvolution using MuSiC, and how to compare cellular inferrences using multiple kinds of reference datasets. You can find the workflow for this tutorial [here](https://usegalaxy.eu/u/wendi.bacon.training/w/music-deconvolution-compare-1) and an example history [here](https://usegalaxy.eu/u/wendi.bacon.training/h/music-deconvolution-comparisons--final).

 We hope this helps you in your research!

![Workflow editor showing 5 inputs and 3 runs of the MuSiC Compare tool](../../images/bulk-music/compare_workflow.png "MuSiC Compare Tutorial Workflow")
