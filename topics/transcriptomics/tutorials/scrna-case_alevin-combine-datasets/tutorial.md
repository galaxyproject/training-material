---
layout: tutorial_hands_on

title: "Combining datasets after pre-processing"
subtopic: single-cell-CS
priority: 2

zenodo_link: 'https://zenodo.org/record/4574153'

questions:
  - I have some single cell FASTQ files I want to analyse. Where do I start?

objectives:
  - Repeat matrix generation for any droplet-based single cell sequencing data
  - Apply data combination and metadata editing for particular experimental designs
  - Interpret quality control (QC) plots to make informed decisions on cell thresholds
  - Find relevant information in GTF files for the particulars of their study, and include this in data matrix metadata

time_estimation: 3H

key_points:
  - Create a scanpy-accessible AnnData object from FASTQ files, including relevant cell and gene metadata
  - Combine multiple samples and label according to study design

tags:
  - single-cell
  - 10x
  - paper-replication
  - español

contributors:
  - nomadscientist
  - pinin4fjords

requirements:
  - type: "internal"
    topic_name: transcriptomics
    tutorials:
        - scrna-intro
        - scrna-umis

translations:
  - es

gitter: Galaxy-Training-Network/galaxy-single-cell

---

# Introduction
{:.no_toc}

<!-- This is a comment. -->

This tutorial will take you from raw FASTQ files to a cell x gene data matrix in AnnData format. What's a data matrix, and what's AnnData format? Well you'll find out! Importantly, this is the first step in processing single cell data in order to start analysing it. Currently you have a bunch of strings of `ATGGGCTT` etc. in your sequencing files, and what you need to know is how many cells you have and what genes appear in those cells. In the second part of this tutorial, we will also look at combining FASTQ files and adding in metadata (for instance, SEX or GENOTYPE) for analysis later on. These steps are the most computationally heavy in the single cell world, as you're starting with 100s of millions of reads, each 4 lines of text. Later on in analysis this data becomes simple gene counts such as 'Cell A has 4 GAPDHs', which is a lot easier to store! Because of this data overload, we have downsampled the FASTQ files to speed up the analysis a bit. Saying that, you're still having to map loads of reads to the massive murine genome, so get yourself a cup of coffee and prepare to analyse!

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Generating a matrix

In this section, we will show you the principles of the initial phase of single-cell RNA-seq analysis: generating expression measures in a matrix. We'll concentrate on droplet-based (rather than plate-based) methodology, since this is the process with most differences with respect to conventional approaches developed for bulk RNA-seq.

Droplet-based data consists of three components: cell barcodes, unique molecular identifiers (UMIs) and cDNA reads. To generate cell-wise quantifications we need to:

 * Process cell barcodes, working out which ones correspond to 'real' cells, which to sequencing artefacts, and possibly correct any barcodes likely to be the product of sequencing errors by comparison to more frequent sequences.
 * Map biological sequences to the reference genome or transcriptome.
 * 'De-duplicate' using the UMIs.

This used to be a complex process involving multiple algorithms, or was performed with technology-specific methods (such as 10X's 'Cellranger' tool)  but is now much simpler thanks to the advent of a few new methods. When selecting methodology for your own work you should consider:

# Combining FASTQ files

This sample was originally one of seven. So to run the other [12 downsampled FASTQ files](https://humancellatlas.usegalaxy.eu/u/wendi.bacon.training/h/alevin-tutorial---all-samples---400k), you can use a [workflow](https://humancellatlas.usegalaxy.eu/u/wendi.bacon.training/w/pre-processing-with-alevin---part-1-imported-from-uploaded-file)! Note - the N705 subsample is unluckily largely junk reads, so emptyDrops doesn't work. Instead, I processed it with Alevin. The total sample runs fine on emptyDrops of course. All these samples are going to take a while, so go and have several cups of tea... Or, better yet, I have [run them myself](https://humancellatlas.usegalaxy.eu/u/wendi.bacon.training/h/pre-processing-with-alevin---part-2---input-generation), and plopped them in a [new clean history](https://humancellatlas.usegalaxy.eu/u/wendi.bacon.training/h/pre-processing-with-alevin---part-2---input) for you to import as a fresh history. Alternatively, you can get data with zenodo.

## Data

> ### {% icon hands_on %} Hands-on: Data upload - Combining files
>
> 1. Create a new history for this tutorial (if you're not importing the history above)
> 2. Import the different AnnData files and the experimental design table from [Zenodo](https://zenodo.org/record/4574153#.YD56YS-l2uU)
>
>    ```
>    {{ page.zenodo_link }}/files/Experimental_Design.tabular
>    {{ page.zenodo_link }}/files/N701-400k-AnnData.h5ad
>    {{ page.zenodo_link }}/files/N702-400k-AnnData.h5ad
>    {{ page.zenodo_link }}/files/N703-400k-AnnData.h5ad
>    {{ page.zenodo_link }}/files/N704-400k-AnnData.h5ad
>    {{ page.zenodo_link }}/files/N705-400k-AnnData.h5ad
>    {{ page.zenodo_link }}/files/N706-400k-AnnData.h5ad
>    {{ page.zenodo_link }}/files/N707-400k-AnnData.h5ad
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype is `h5ad`, otherwise you will need to change each file to `h5ad`!
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

Inspect the {% icon galaxy-eye %} `Experimental Design` text file. This shows you how each `N70X` corresponds to a sample, and whether that sample was from a male or female. This will be important metadata to add to our sample, which we will add very similarly to how you added the `gene_name` and `mito` metadata above!

## Concatenating objects
> ### {% icon hands_on %} Hands-on: Concatenating AnnData objects
>
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.7.5+galaxy0){% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `N701-400k-AnnData`
>    - *"Function to manipulate the object"*: 'Concatenate along the observations axis'
>    - {% icon param-file %} *"Annotated data matrix to add"*: 'Select all the other matrix files from bottom to top'
>    - *"Join method"*: `Intersection of variables`
>    - *"Key to add the batch annotation to obs"*: `batch`
>    - *"Separator to join the existing index names with the batch category"*: `-`
{: .hands_on}

Now let's look at what we've done! Unfortunately, AnnData objects are quite complicated, so the {% icon galaxy-eye %} won't help us too much here. Instead, we're going to use a tool to look into our object from now on.

> ### {% icon hands_on %} Hands-on: Inspecting AnnData Objects
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: output of **Manipulate AnnData** {% icon tool %}
>    - *"What to inspect?"*: `General information about the object`
> 2. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: output of **Manipulate AnnData** {% icon tool %}
>    - *"What to inspect?"*: `Key-indexed observations annotation (obs)`
> 3. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: output of **Manipulate AnnData** {% icon tool %}
>    - *"What to inspect?"*: `Key-indexed annotation of variables/features (var)`
{: .hands_on}

Now have a look at the three {% icon tool %} **Inspect AnnData** outputs.

> ### {% icon question %} Question
>
> 1. How many cells do you have now?
> 2. Where is `batch` information stored?
>
> > ### {% icon solution %} Solution
> >
> > 1. If you look at the **General information** {% icon tool %} output, you can see there are now `4079 cells`, as the matrix is now 4079 cells x 35734 genes. You can see this as well in the **obs** {% icon tool %} (cells) and **var** {% icon tool %} (genes) file sizes.
> > 2. Under **Key-indexed observations annotation (obs)**. Different version of the Manipulate tool will put the `batch` columns in different locations. The tool version in this course has the `9th` column at the farthest right is `batch`. Batch refers to the order in which the matrices were added. The files are added from the bottom of the history upwards, so be careful how you set up your histories when running this!
> {: .solution}
>
{: .question}

# Adding batch metadata

I set up the example history with the earliest indices at the bottom.

![Ordered history](../../images/wab-history-files-ascending.png "Note how N701 is lowest, ordered ascending to N707")

Therefore, when it is all concatenated together, the `batch` appears as follows:

| Index | Batch | Genotype | Sex |
|------ |--------------------|
| N701 | 0    | wildtype    | male    |
| N702 | 1    | knockout   | male    |
| N703 | 2    | knockout   | female    |
| N704 | 3    | wildtype    | male    |
| N705 | 4    | wildtype    | male    |
| N706 | 5    | wildtype    | male    |
| N707 | 6    | knockout    | male    |

If you used Zenodo to import files, they may not have imported in order (i.e. N701 to N707, ascending). In that case, you will need to tweak the parameters of the next tools appropriately to label your batches correctly!

The two critical pieces of metadata in this experiment are **sex** and **genotype**. I will later want to color my cell plots by these parameters, so I want to add them in now!

> ### {% icon hands_on %} Hands-on: Labelling sex
>
> 1. {% tool [Replace Text in a specific column](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: output of **Inspect AnnData: Key-indexed observations annotation (obs)** {% icon tool %})
>    - *"1. Replacement"*
>
>         - *"in column"*: `Column: 9` - or whichever column `batch` is in
>         - *"Find pattern"*: `0|1|3|4|5|6`
>         - *"Replace with"*: `male`
>    - **+ Insert Replacement**
>    - *"2. Replacement"*
>
>         - *"in column"*: `Column: 9`
>         - *"Find pattern"*: `2`
>         - *"Replace with"*: `female`
>    - **+ Insert Replacement**
>    - *"3. Replacement"*
>
>         - *"in column"*: `Column: 9`
>         - *"Find pattern"*: `batch`
>         - *"Replace with"*: `sex`
>
>    Now we want only the column containing the sex information - we will ultimately add this into the cell annotation in the AnnData object.
>
> 2. {% tool [Cut columns from a table](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c9`
>    - *"Delimited by"*: `Tab`
>    - % icon param-file %} *"From"*: output of **Replace text** {% icon tool %}
>
> 3. Rename {% icon galaxy-pencil %} output `Sex metadata`
{: .hands_on}

That was so fun, let's do it all again but for genotype!

> ### {% icon hands_on %} Hands-on: Labelling genotype
>
> 1. {% tool [Replace Text in a specific column](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: output of **Inspect AnnData: Key-indexed observations annotation (obs)** {% icon tool %}
>    - *"1. Replacement"*
>
>         - *"in column"*: `Column: 9`
>         - *"Find pattern"*: `0|3|4|5`
>         - *"Replace with"*: `wildtype`
>    - **+ Insert Replacement**
>    - *"2. Replacement"*
>
>         - *"in column"*: `Column: 9`
>         - *"Find pattern"*: `1|2|6`
>         - *"Replace with"*: `knockout`
>    - **+ Insert Replacement**
>    - *"3. Replacement"*
>
>         - *"in column"*: `Column: 9`
>         - *"Find pattern"*: `batch`
>         - *"Replace with"*: `genotype`
>
>    Now we want only the column containing the genotype information - we will ultimately add this into the cell annotation in the AnnData object.
>
> 2. {% tool [Cut columns from a table](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c9`
>    - *"Delimited by"*: `Tab`
>    - {% icon param-file %} *"From"*: output of **Replace text** {% icon tool %}
>
> 3. Rename {% icon galaxy-pencil %} output `Genotype metadata`
{: .hands_on}

You might want to do this with all sorts of different metadata - which labs handled the samples, which days they were run, etc. Once you've added all your metadata columns, we can add them together before plugging them into the AnnData object itself.

> ### {% icon hands_on %} Hands-on: Combining metadata columns
>
> 1. {% tool [Paste two files side by side](Paste1) %} with the following parameters:
>    - {% icon param-file %} *"Paste"*: `Genotype metadata`
>    - {% icon param-file %} *"and"*: `Sex metadata`
>    - *"Delimit by"*: `Tab`
> 2. Rename {% icon galaxy-pencil %} output `Cell Metadata`
{: .hands_on}

Let's add it to the AnnData object!

> ### {% icon hands_on %} Hands-on: Adding metadata to AnnData object
>
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.7.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: output of previous **Manipulate AnnData** {% icon tool %}
>    - *"Function to manipulate the object"*: `Add new annotation(s) for observations or variables`
>    - *"What to annotate?"*: `Observations (obs)``
>    - {% icon param-file %} *"Table with new annotations"*: `Cell Metadata`
{: .hands_on}

Woohoo! We're there! You can run an **Inspect AnnData** to check now, but I want to clean up this AnnData object just a bit more first. It would be a lot nicer if 'batch' meant something, rather than 'the order in which the Manipulate AnnData tool added my datasets'.

> ### {% icon hands_on %} Hands-on: Labelling batches
>
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.7.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: output of **Manipulate AnnData - Add new annotations** {% icon tool %}
>    - *"Function to manipulate the object"*: `Rename categories of annotation`
>    - *"Key for observations or variables annotation"*: `batch`
>    - *"Comma-separated list of new categories"*: `N701,N702,N703,N704,N705,N706,N707`
{: .hands_on}


Huzzah! We are JUST about there. However, while we've been focussing on our cell metadata (sample, batch, genotype, etc.) to relabel the 'observations' in our object...

# Mitochondrial reads

Do you remember when we mentioned mitochondria early on in this tutorial? And how often in single cell samples, mitochondrial RNA is often an indicator of stress during dissociation? We should probably do something with our column of true/false in the gene annotation that tells us information about the cells. You will need to do this whether you have combined FASTQ files or are analysing just one (and thus skipping sections 4 & 5).

> ### {% icon hands_on %} Hands-on: Calculating mitochondrial RNA in cells
>
> 1. {% tool [AnnData Operations](toolshed.g2.bx.psu.edu/repos/ebi-gxa/anndata_ops/anndata_ops/0.0.3+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input object in hdf5 AnnData format"*: output of **Manipulate AnnData - Rename categories** {% icon tool %}
>    - *"Format of output object"*: `AnnData format`
>    - *"Copy AnnData to .raw"*: `No`
>    - *"Gene symbols field in AnnData"*: `NA.`
>    - *"Flag genes that start with these names"*: `Insert Flag genes that start with these names`
>    - *"Starts with"*: `True`
>    - *"Var name"*: `mito`
>    - *"Number of top genes"*: `50`
{: .hands_on}

{% icon congratulations %}Well done!  I strongly suggest have a play with the **Inspect AnnData** {% icon tool %} on your final `Pre-processed object` to see the wealth of information that has been added. You are now ready to move along to further filtering! There is a cheat that may save you time in the future though...

# Pulling single cell data from public resources

If you happen to be interested in analysing publicly available data, particularly from the [Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/home), you may be interested in the following tool {% cite Moreno2020.04.08.032698 %} which rather skips forward all these steps in one! For this tutorial, the dataset can be seen [here](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-6945/downloads) with experiment id of `E-MTAB-6945`.

> ### {% icon hands_on %} Hands-on: Retrieving data from Single Cell Expression Atlas
>
> 1. {% tool [EBI SCXA Data Retrieval](toolshed.g2.bx.psu.edu/repos/ebi-gxa/retrieve_scxa/retrieve_scxa/v0.0.2+galaxy2) %} with the following parameters:
>      - *"SC-Atlas experiment accession"*: `E-MTAB-6945`
>      - *"Choose the type of matrix to download"*: `Raw filtered counts`
>
>    Now we need to transform this into an AnnData objects
>
> 2. {% tool [Scanpy Read10x](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_read_10x/scanpy_read_10x/1.6.0+galaxy0) %} with the following parameters:
>    - *"Expression matrix in sparse matrix format (.mtx)"*: `EBI SCXA Data Retrieval on E-MTAB-6945 matrix.mtx (Raw filtered counts)`
>    - *"Gene table"*:  `EBI SCXA Data Retrieval on E-MTAB-6945 genes.tsv (Raw filtered counts)`
>    - *"Barcode/cell table"*: `EBI SCXA Data Retrieval on E-MTAB-6945 barcodes.tsv (Raw filtered counts)`
>    - *"Cell metadata table"*: `EBI SCXA Data Retrieval on E-MTAB-6945 exp_design.tsv`
{: .hands_on}

It's important to note that this matrix is processed somewhat through the SCXA pipeline, which is quite similar to this tutorial, and it contains any and all metadata provided by their pipeline as well as the authors (for instance, more cell or gene annotations).

# Conclusion
{:.no_toc}

![Workflow Part 1](../../images/wab-alevin-part1workflow.png "Workflow  - Steps 1-3")

![Workflow Part 2](../../images/wab-alevin-part2workflow.png "Workflow  - Steps 4-6")

You've reached the end of this session!
You may be interested in seeing an [example history](https://humancellatlas.usegalaxy.eu/u/wendi.bacon.training/h/pre-processing-with-alevin---part-2---answer-key-1) and [Part 2 workflow](https://humancellatlas.usegalaxy.eu/u/wendi.bacon.training/w/pre-processing-with-alevin---part-2). Note that the workflow will require changing of the `column` containing the batch metadata depending on how you are running it. The final object containing all the reads can be found in [here](https://humancellatlas.usegalaxy.eu/u/wendi.bacon.training/h/pre-processing-with-alevin---part-2---total-anndata-example).

We have:

 * Taken raw read data and annotations and necessary input files for quantification.
 * Run Alevin in two different parameterisations, both allowing Alevin to make its own calls on what constitutes empty droplets, and applying emptyDrops instead.
 * Deployed barcode rank plots as a way of quickly assessing the signal present in droplet datasets.
 * Applied the necessary conversion to pass these data to downstream processes.
 * Retrieved partially analysed data from the Single Cell Expression Atlas

 To discuss with like-minded scientists, join our Gitter channel for all things Galaxy-single cell!
 [![Gitter](https://badges.gitter.im/Galaxy-Training-Network/galaxy-single-cell.svg)](https://gitter.im/Galaxy-Training-Network/galaxy-single-cell?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
