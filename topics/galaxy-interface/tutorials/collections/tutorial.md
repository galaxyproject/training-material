---
layout: tutorial_hands_on
redirect_from:
  - /topics/galaxy-data-manipulation/tutorials/collections/tutorial
  - /topics/galaxy-data-manipulation/tutorials/processing-many-samples-at-once/tutorial
  - /topics/galaxy-interface/tutorials/processing-many-samples-at-once/tutorial

title: "Using dataset collections"
zenodo_link: 'https://zenodo.org/record/5119008'
tags:
- collections
level: Intermediate
questions:
  - "How to manipulate large numbers of datasets at once?"
objectives:
  - "Understand and master dataset collections"
time_estimation: "30m"
key_points:
  - "Multiple datasets can be combined in a collection."
  - "This significantly simplifies the analysis."
  - "This tutorial showed how to (1) create collection, (2) run tools on a collection, (3) combine collection elements into a final analysis results."
  - "There are a variety of **Collection operation** tools allowing to perform a variety of transformations."
contributions:
  authorship:
    - nekrut
    - tflowers15
  funding:
    - unimelb
    - melbournebioinformatics
    - AustralianBioCommons

subtopic: manage

recordings:
- youtube_id: uN6nP3I7QLE
  captioners:
  - delphine-l
  date: '2023-05-09'
  length: 13M
  speakers:
  - delphine-l
- captioners:
  - annefou
  - mariipia10
  date: '2021-08-09'
  length: 11M
  youtube_id: uZUt9XIHUQo
  speakers:
  - nekrut

---

{% snippet faqs/gtn/galaxy_tested_with.md version="26.1" %}

Here we will show Galaxy features designed to help with the analysis of large numbers of samples. When you have just a few samples - clicking through them is easy. But once you've got hundreds - it becomes very annoying. In Galaxy we have introduced **Dataset collections** that allow you to combine numerous datasets in a single entity that can be easily manipulated.


# Getting data

First, we need to upload datasets. Cut and paste the following URLs to Galaxy upload tool (see a {% icon tip %} **Tip** on how to do this [below](#tip-upload-fastqsanger-datasets-via-links)).

```
https://zenodo.org/record/5119008/files/M117-bl_1.fq.gz
https://zenodo.org/record/5119008/files/M117-bl_2.fq.gz
https://zenodo.org/record/5119008/files/M117-ch_1.fq.gz
https://zenodo.org/record/5119008/files/M117-ch_2.fq.gz
https://zenodo.org/record/5119008/files/M117C1-bl_1.fq.gz
https://zenodo.org/record/5119008/files/M117C1-bl_2.fq.gz
https://zenodo.org/record/5119008/files/M117C1-ch_1.fq.gz
https://zenodo.org/record/5119008/files/M117C1-ch_2.fq.gz
```

> <hands-on-title>Set format to `fastqsanger.gz`</hands-on-title>
> The above datasets are in `fastqsanger.gz` format. Galaxy is fairly good at correctly identifying the datatype/format for imported/uploaded datasets. However, if Galaxy does not default to the correct datatype, this can be set manually. Another option available is to specify the datatype when uploading/importing the datasets. The {% icon tip %} **Tip** sections below explains how to upload these data and set the correct format. There is a variety of [fastq format flavors](https://en.wikipedia.org/wiki/FASTQ_format) and it can be difficult to guess them automatically.
>
> {% snippet faqs/galaxy/datasets_import_via_link.md format="fastqsanger.gz" %}
>
> {% snippet faqs/galaxy/datasets_change_datatype.md format="fastqsanger.gz" %}
>
> {% snippet topics/assembly/tutorials/vgp_genome_assembly/faqs/dataset_upload_fastqsanger_via_urls.md %}
>
{: .hands_on}

## About these datasets

These datasets represent genomic DNA (enriched for mitochondria via a long range PCR) isolated from blood (`bl`) and cheek (buccal swab, `ch`) of mother (`M117`) and her child (`M117C1`) that was sequenced on an Illumina miSeq machine as paired-read library (250-bp reads; see our [2014](http://www.pnas.org/content/111/43/15474.abstract) manuscript for **Methods**):

- `M117-bl_1` - family 117, mother, *forward* (**F**) reads from **blood**
- `M117-bl_2` - family 117, mother, *reverse* (**R**) reads from **blood**
- `M117-ch_1` - family 117, mother, *forward* (**F**) reads from **cheek**
- `M117-ch_2` - family 117, mother, *reverse* (**R**) reads from **cheek**
- `M117C1-bl_1`- family 117, child, *forward* (**F**) reads from **blood**
- `M117C1-bl_2`- family 117, child, *reverse* (**R**) reads from **blood**
- `M117C1-ch_1`- family 117, child, *forward* (**F**) reads from **cheek**
- `M117C1-ch_2`- family 117, child, *reverse* (**R**) reads from **cheek**

# Creating a paired dataset collection

You can see that there are eight datasets forming four pairs. Obviously, we can manipulate them one-by-one (e.g. start four mapping jobs, call variants four times and so on), but this will unnecessarily tedious. Moreover, imagine if you have 100s or 1,000s of pairs: it will be impossible to process them individually.

This is exactly why we developed collections. Dataset collections allow combining multiple datasets into a single entity. Thus, instead of dealing with four, a hundred, or a thousand of individual datasets you have only one item in Galaxy history to deal with.

Because our data is *paired* we need to create a hierarchical collection called a **Paired Dataset Collection** or **Paired Collection**. In such collections there are two layers. The first layer corresponds to individual samples (e.g. `M117-bl`). The second layer represent `forward` and `reverse` reads corresponding to each sample:


![paired collection](../../images/collections/paired_collection.svg "The logic of Paired Collection. Here <em>N</em> datasets are bundled into a paired collection with two layers. The first layer corresponds to samples and the second to forward and reverse reads within each sample.")


## Selecting Datasets and Building a Collection

To begin creating a collection we need to select datasets we would like to bundle. This is done using checkbox button {% icon galaxy-selector %} of Galaxy's history menu. Figure 2 below shows this process.

![selecting multiple datasets for list builder](../../images/collections/collections_select_datasets_build.gif  "Selecting multiple datasets and creating a paired collection using the list builder.")

After selecting the datasets to be included in the collection, there are two options for building the collection: (1) **Auto Build List** and (2) **Advanced Build List**. In the above figure, the **Advanced Build List** option was selected, which opens the Galaxy collection builder at the first step which shows all the collection builder options. In this case Galaxy has automatically identified that the selected datasets are paired-end as the dataset names contain the text `_1` or `_2` and has selected the **List of Paired Datasets** option. The selected option is indicated by the box outlined in dark blue. If the selected datasets are not recognised as paired-end by Galaxy, then the **Flat List** option will be selected by default.

> <comment-title>Using Auto Build List</comment-title>
>
> If using the **Auto Build List** option after selecting datasets, the collection builder will determine if the selected datasets are best suited for the **Flat List** or **List of Paired Datasets** option based on the dataset names. If all datasets have a paired-end indicators (e.g. `_F` and `_R`, `_1` and `_2`, or `.1` and `.2`) and Galaxy is able to match up pairs of datasets with identical names only differing by the paired-end indicator then Galaxy will select **List of Paired Datasets**. Alternatively, if there are no paired-end indicators identified then Galaxy will select **Flat List**. Figure 3 below shows the **Auto Build List** functionality.
>
> ![selecting multiple datasets for list autobuilder](../../images/collections/collections_select_datasets_autobuild.gif  "Selecting multiple datasets and creating a paired collection using the list auto builder.")
> 
{: .comment}

> <comment-title>A Note on Galaxy Terminology</comment-title>
>
> Within Galaxy the terms "Collection" and "List" are used interchangeably.
> 
{: .comment}

> <comment-title>Collection/List Builder</comment-title>
>
> There are six options for creating collections available through the Galaxy **List Builder**:
> - **Flat List**: This option creates a simple flat list of files. If your data isn't nested and does not contained paired datasets, this is the the option to choose.
> - **List of Paired Datasets**: This option creates a list of paired datasets (e.g., paired-end data from Illumina and ElementBio), use this option if all of your data should be paired off.
> - **Rule Builder**: The rule builder allows building arbitrarily nested collections of data with reproducible metadata parsing and other advanced features.
> - **Mixed List of Paired and Unpaired Datasets**: For studies that have a mix of paired and unpaired data, this may be the right choice. Existing tools and workflows may have to be updated to allow this modality though.
> - **Nested List of Datasets**: A nested list of datasets. This type of collection is a list where each element is in turn a list of datasets. The outer list identifier might group samples together by condition, sample type, etc..
> - **Nested List of Dataset Pairs**: A nested list of datasets. This type of collection is a list where each element is in turn a list of pairs of datasets.
>
> The most commonly used options are **Flat List** and **List of Paired Datasets**.
> 
{: .comment}

After selecting the type of **List Builder**, in this instance the **List of Paired Datasets**, click **Next** at the bottom of the page to progress to the **Auto Pairing** page. In this case Galaxy automatically assigned pairs using the `_1` and `_2` endings of dataset names. The paired datasets are shown under **Auto-matched Pairs**. Let's however pretend that this did not happen. Delete the text in the two **Filter text** text boxes (free text boxes pre-filled with `_1` and `_2`). The interface will change into its unpaired state (all datasets will be listed under **Un-matched Datasets**).

Here datasets containing the first (forward) and the second (reverse) read are differentiated by having `_1` and `_2` in the filename. We can use this feature in the dataset collection builder to pair our datasets.  Type `_1` in the left **Filter text** text box and `_2` in the right. You will see that the dataset collection builder will automatically pair the datasets once both filters are set. Once all datasets have been assigned to a pair, click on **Next** at the bottom of the page to progress to the **Builder** page.

The **Builder** page again shows the datasets to be assigned to each pair and the corresponding *List Identifier* that will be used in the created collection (this will generally be the base of the dataset name, i.e. without the paired-end indicator). There are also options to remove the extensions from the dataset names (e.g. `.fq.gz`) and to hide the original dataset elements in the history, which are both toggled on by default. The new collection must be given a name before it can be created. Type `M117-collection` in the **Name** text box and create the collection by clicking **Build**. A new collection item will appear in the history and the original `.fq.gz` datasets will be hidden. This process is demonstrated below in Figure 4.

![using collection builder](../../images/collections/collections_create_collection.gif "Working with the collection builder. The text above this figure explains each step.")

## Viewing Datasets in a Paired-end Collection

Clicking on collection will expand it to show four pairs it contains. Clicking individual pairs will expand them further to reveal **forward** and **reverse** datasets. Expanding these further will enable one to see individual datasets, as shown below in Figure 5.

![expanding collection](../../images/collections/collections_viewing_datasets.gif "To look at what is inside a collection, just click on it.")

# Processing data organized as a collection

By now we see that a collection can be used to bundle a large number of items into a single history item. Galaxy tools take collection as input. Let's map reads contained in collection `M117-collection` against human mitochondrial genome. Before we can do this we need to upload mitochondrial genome using the following URL (see a {% icon tip %} **Tip** on how to do this [below](#tip-importing-via-links)):

```
https://zenodo.org/record/5119008/files/chrM.fa.gz
```

> <hands-on-title>Set format to `fasta.gz`</hands-on-title>
> The above dataset is in `fasta.gz` format. The {% icon tip %} **Tip** section below explains how to upload these data and set the correct format.
>
> {% snippet faqs/galaxy/datasets_import_via_link.md reset_form="True" link="https://zenodo.org/record/5119008/files/chrM.fa.gz" format="fasta.gz" %}
{: .hands_on}

## Mapping reads

{% tool Map with BWA-MEM %} is a widely used sequence aligner for short-read sequencing datasets such as those we are analysing in this tutorial. (You can find the tool by typing `BWA MEM` in the search box at the top left corner of Galaxy interface).

> <hands-on-title>Map sequencing reads to reference genome</hands-on-title>
>
> Run {% tool [Map with BWA-MEM](toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.17.1) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `chrM.fa.gz` (The mitochondrial genome we just uploaded)
>    - *"Single or Paired-end reads"*: `Paired Collection`
>        - {% icon param-file %} *"Select a paired collection"*: `M117-collection` (the collection we built at the beginning of this tutorial.)
>    - *"Set read groups information?"*: `Do not set`
>    - *"Select analysis mode"*: `1.Simple Illumina mode`
>
> The interface should look like this:
>
> ![bwa_mem_interface](../../images/collections/bwa_mem_interface_coll_tut.png "Tool interface")
>
>    - Click **Run Tool** button
>
{: .hands_on}

You will see jobs being submitted and new datasets appearing in the history. Because our collection contains four paired datasets Galaxy will actually generate four separate {% tool Map with BWA-MEM %} jobs. In the end this {% tool Map with BWA-MEM %} run will produce a new collection containing four (4) BAM datasets. Let's look at this collection by clicking on it (panel **A** in the figure below). You can see that now this collection is no longer paired (compared to the collection we created in the beginning of this tutorial). This is because {% tool Map with BWA-MEM %} takes forward and reverse data as input, but produces only a single BAM dataset as the output. So what we have in the result is a *list* of four dataset (BAM files; panel **B**). If you click on any of the datasets you will see that it is indeed a BAM dataset (panel **C**).

![bwa_memCollection_ABC](../../images/collections/collection_expansion.svg)

## Calling variants

After we mapped reads against the mitochondrial genome, we can now call variants. In this step a variant calling tool {% tool Call variants with LoFreq %} will take a collection of BAM datasets (the one produced by {% tool Map with BWA-MEM %}), identify differences between reads and the reference, and output these differences as a collection of [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) datasets.

> <hands-on-title>Call variants  </hands-on-title>
>
> Run {% tool [Call variants](toolshed.g2.bx.psu.edu/repos/iuc/lofreq_call/lofreq_call/2.1.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input reads in BAM format"*: `Map with BWA-MEM...` (output of {% tool Map with BWA-MEM %})
>    - *"Choose the source for the reference genome"*: `History`
>        - {% icon param-file %} *"Reference"*: `chrM.fa.gz (as fasta)` (Input dataset)
>    - *"Call variants across"*: `Whole reference`
>    - *"Types of variants to call"*: `SNVs and indels`
>
> The interface should look like this:
>
> ------
>
> ![lofreq_interface](../../images/collections/lofreq_interface.png)
>
> ------
>
>    - Click **Run Tool** button
>
{: .hands_on}

## Create table of variants using **SnpSift Extract Fields**

We will now convert VCF datasets into tab delimited format as it will be easier to work with. This will be done with {% tool SnpSift %}: a tool specifically designed for manipulation of tab-delimited data.


> <hands-on-title>Create table of variants</hands-on-title>
>
> Run {% tool [SnpSift Extract Fields](toolshed.g2.bx.psu.edu/repos/iuc/snpsift/snpSift_extractFields/4.3+t.galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Variant input file in VCF format"*: `Call variants on collection...` (output of {% tool Call variants with LoFreq %})
>    - *"Fields to extract"*: `CHROM POS REF ALT QUAL DP AF SB DP4`
>    - *"One effect per line"*: `Yes`
>
> The interface should look like this:
>
> ------
>
> ![snpsift_interface](../../images/collections/snpsift_interface.png)
>
> ------
>
>    - Click **Run Tool** button
>
{: .hands_on}


As a result of this operation we now have a collection of four tab delimited files. Yet, ultimately we want to summaize these data as one final table. The next step does just that.

## Collapse data into a single dataset

We now extracted meaningful fields from VCF datasets. But they still exist as a collection. To move towards secondary analysis we need to **collapse** this collection into a single dataset with {% tool Collapse Collection %}. For more information about collapsing collections see this video:

{% include _includes/youtube.html id="ypuFZ1RKMIY" title="Collapsing collections" %}

> <hands-on-title>Collapse a collection</hands-on-title>
>
> Run {% tool [Collapse Collection](toolshed.g2.bx.psu.edu/repos/nml/collapse_collections/collapse_dataset/4.0) %} with the following parameters:
>    - {% icon param-collection %} *"Collection of files to collapse into single dataset"*: `SnpSift Extract Fields ...` (output of {% tool SnpSift Extract Fields %})
>    - "*Keep one header line*": `Yes`
>    - "*Prepend File name*": `Yes`
>    - "*Where to add dataset name*": `Same line and each line in dataset`
>
> The interface should look like this:
>
> ------
>
> ![snpsift_interface](../../images/collections/collapse_collection.png)
>
> ------
>
>    - Click **Run Tool** button
>
{: .hands_on}

You can see that this tool takes lines from all collection elements (in our case we have two), add element name as the first column, and pastes everything together. So if we have a collection as an input:

> <code-in-title>A collection with two items</code-in-title>
> A collection element named `M117-bl.fq`
>
>```
>chrM   152 T C  3707.0 1242 0.99 0 2,2,540,697
>chrM 16519 T C 35149.0 1033 0.99 0 1,1,611,420
>```
>
>A collection element named `M117-ch.fq`:
>
>```
chrM   152 T C  4098.0 1440 0.99 0 0,1,575,863
chrM 16519 T C 36574.0 1039 0.99 2 3,0,713,321
>```
>
>A collection element named `M117C1-bl.fq`:
>
>```
>chrM   152 T C  4888.0 1235 1.00 0 0,0,548,687
>chrM 16519 T C 35220.0 1042 0.99 0 0,0,598,443
>```
>
>A collection element named `M117C1-ch.fq`:
>
>```
>chrM   152 T C  2757.0 1413 0.99 0 2,2,576,833
>chrM 16455 G A    54.0  100 0.04 0 89,7,4,0
>chrM 16519 T C 36363.0 1061 0.99 6 3,4,691,362
>```
>
{: .code-in}

We will have a single dataset as the output:

> <code-out-title>A single dataset</code-out-title>
>
>then the {% tool Collapse Collection %} will produce this:
>
>```
>M117-bl.fq   chrM   152 T C  3707.0 1242 0.99 0 2,2,540,697
>M117-bl.fq   chrM 16519 T C 35149.0 1033 0.99 0 1,1,611,420
>M117-ch.fq   chrM   152 T C  4098.0 1440 0.99 0 0,1,575,863
>M117-ch.fq   chrM 16519 T C 36574.0 1039 0.99 2 3,0,713,321
>M117C1-bl.fq chrM   152 T C  4888.0 1235 1.00 0 0,0,548,687
>M117C1-bl.fq chrM 16519 T C 35220.0 1042 0.99 0 0,0,598,443
>M117C1-ch.fq chrM   152 T C  2757.0 1413 0.99 0 2,2,576,833
>M117C1-ch.fq chrM 16455 G A    54.0  100 0.04 0 89,7,4,0
>M117C1-ch.fq chrM 16519 T C 36363.0 1061 0.99 6 3,4,691,362
>```
{: .code-out}

you can see that added a column with dataset ID taken from collection element name.

## We did not fake this:
The history described in this page is accessible directly from here:

* History [**Collections (full analysis)**](https://usegalaxy.org/u/cartman/h/collection-tutorial)

From there you can import histories to make them your own.


# Collection operations

In this brief analysis we took four paired datasets, created a collection, analyzed this collection and finally created a single report. Such "lifecycle" is shown in the figure below. Here we started with eight fastq datasets representing four paired end samples. A paired collection was reduced to a list of BAM datasets by {% tool Map with BWA-MEM %}. Varinat calling by {% tool Call variants with LoFreq %} and field extraction with {% tool SnpSift %} maintained collection structure: these tools processed four individual datasets changing their formats from BAM to VCF, and from VCF to Tab-delimited. Finally, we collapsed collection by merging its content into a single dataset.

![Collection lifecycle](../../images/collections/collection_lifecycle.svg "Collection lifecycle. Arrows = individual fastq datasets; Four shades of yellow = four samples analyzed in this example. ")

The last step of our analysis, collapsing a collection, is an example of a *collection operation*. Galaxy contains an entire section of tools designed for handling of collection data. These can be classified as:

 - Tools that manipulate elements within a collection
 - Tools that change collection structure
 - Tools that combine elements of a collection

Let's look at these categories in more detail:

## Tools that manipulate elements within a collection

### Extract dataset

{% tool Extract dataset %} extracts datasets from a collection based on either position or identifier.

The tool allow extracting datasets based on position (**The first dataset** and **Select by index** options) or name (**Select belement identifier** option). This tool effectively collapses the inner-most collection into a dataset. For nested collections (e.g a list of lists of lists: outer:middle:inner, extracting the inner dataset element) a new list is created where the selected element takes the position of the inner-most collection (so outer:middle, where middle is not a collection but the inner dataset element).

### Filter empty

{% tool Filter empty %} removes empty elements from a collection.

This tool takes a dataset collection and filters out (removes) empty datasets. This is useful for continuing a multi-sample analysis when downstream tools require datasets to have content.

![Filtering empty datasets](../../images/collections/filter_empty.svg)

### Filter failed datasets

{% tool Filter failed datasets %} removes datasets in error (red) from a collection.

This tool takes a dataset collection and filters out (removes) datasets in the failed (red) state. This is useful for continuing a multi-sample analysis when one or more of the samples fails at some point.

![filter_failed](../../images/collections/filter_error.svg)

### Build list

{% tool Build list %} creates a new list collection from individual datasets or collections.

This tool combines individual datasets or collections into a new collection. The simplest scenario is building a new colection from individual datasets (case **A** in the image below). You can merge a collection with individual dataset(s). In this case (see **B** in the image below) the individual dataset(s) will be merged with each element of the input collection to create a nested collection. Finally, two or more collection can be merged together creating a nested collection (case **C** in the image below).

![Build list](../../images/collections/build_list.svg)

### Filter collection

{% tool Filter collections %} removes elements from a collection using a list supplied in a file.

This tools allow filtering elements from a data collection.  It takes an input collection and a text file with names (i.e. identifiers). The tool behaviour is controlled by **How should the elements to remove be determined?** drop-down. It has the following options:

#### Remove if identifiers are ABSENT from file

Given a collection:

```
 Collection: [Dataset A]
             [Dataset B]
             [Dataset X]
```

and a text file:

```
             A
             B
             Z
```

the tool will return two collections:


```

 (filtered):  [Dataset A]
              [Dataset B]

 (discarded): [Dataset X]
```

#### Remove if identifiers are PRESENT in file

Given a collection:

```
 Collection: [Dataset A]
             [Dataset B]
             [Dataset X]
```
and a text file:

```
             A
             B
             Z
```

the tool will return two collections:

```
 (filtered):  [Dataset X]

 (discarded): [Dataset A]
              [Dataset B]

```

### Relabel identifiers

{% tool Relabel identifiers %} changes identifiers of datasets within a collection using identifiers from a supplied file.

New identifiers can be supplied as either a simple list or a tab-delimited file mapping old identifier to the new ones. This is controlled using **How should the new identifiers be specified?** drop-down:

#### Using lines in a simple text file

Given a collection:

```
 Collection: [Dataset A]
             [Dataset B]
             [Dataset X]
```

and a simple text file:


```
             Alpha
             Beta
             Gamma
```

the tool will return:

```
 Collection: [Dataset Alpha]
             [Dataset Beta]
             [Dataset Gamma]
```

#### Map original identifiers to new ones using a two column table

Given a collection:

```
 Collection: [Dataset A]
             [Dataset B]
             [Dataset X]
```

and a simple text file (you can see that entries do not have to be in order here):


```
             B Beta
             X Gamma
             A Alpha
```

the tool will return:

```
 Collection: [Dataset Alpha]
             [Dataset Beta]
             [Dataset Gamma]
```

### Sort collection

{% tool Sort collection %} sorts dataset collection alphabetically, numerically, or using predetermined order from a supplied file.

#### Numeric sort

The tool sort in ascending order. When *numeric* sort is chosen, the tool ignores non-numeric characters. For example, if a collection contains the following elements:

```
 Collection: [Horse123]
             [Donkey543]
             [Mule176]
```

The tool will output:

```
 Collection: [Horse123]
             [Mule176]
             [Donkey543]
```

#### Sorting from file

Alternative, one can supply a single column text file containing elements identifiers in the desired sort order. For example, suppose there a collection:

```
 Collection: [Horse123]
             [Donkey543]
             [Mule176]
```

and a file specifying sort order:

```
 Donkey543
 Horse123
 Mule176
```

the output will predictably look like this:

```
 Collection: [Donkey543]
             [Horse123]
             [Mule176]
```

### Tag elements

{% tool Tag elements %} adds tags (including name: and group: tags) to collection elements.

The relationship between element names and tags is specified in a two column tab-delimited file. This file may contain less entries than elements in the collection. In that case only matching list identifiers will be tagged.

To create name: or group: tags prepend them with `#` (you can also use `name:`) or `group:`, respectively.

More about tags

> <tip-title>More about tags</tip-title>
> Galaxy allows tagging datasets to facilitate analyses. There are several types of tags including simple tags, name tags, and group tags. **Simple** tags allow you to attach an alternative label to a dataset, which will make it easier to find it later. **Name** tags allow you to track propagation of a dataset through the analyses: all datasets derived from the initial dataset labeled with a name tag will inherit it. Finally, **group** tags allow you to label group of datasets. This is useful. for example, for differential expression analysis where you can have two groups of datasets labeled as "treatment" and "control".
>
>To learn mote about tags go to [training site]({% link search.md %}?query=tags).
{: .tip}

## Tools that change collection structure

### Flatten collection

{% tool Flatten collection %} collapses nested collection into a simple list.

This tool takes nested collections such as a list of lists or a list of dataset pairs and produces a flat list from the inputs. It effectively "flattens" the hierarchy. The collection identifiers are merged together (using `_` as default) to create new collection identifiers in the flattened result:

![Flatten collection](../../images/collections/flatten.svg)

### Merge collections

{% tool Merge collections %} takes two or more collections and creates a single collection from them.

By default the tool assumes that collections that are being merged have unique dataset names. If it not the case only one (the first) of the datasets with a repeated name will be included in the merged collection. For example, suppose you have two collections. Each has two datasets named "A" and "B":
```
 Collection 1: [Dataset A]
               [Dataset B]
               [Dataset X]
 Collection 2: [Dataset A]
               [Dataset B]
               [Dataset Y]
```

Merging them will produce a single collection with only two datasets:

```
 Merged Collection: [Dataset A]
                    [Dataset B]
                    [Dataset X]
                    [Dataset Y]
```

This behavior can be changed by clicking on "*Advanced Options*" link. The following options are available:

#### Keep first instance (Default behavior)

Input:

```
 Collection 1: [Dataset A]
               [Dataset B]
               [Dataset X]
 Collection 2: [Dataset A]
               [Dataset B]
               [Dataset Y]
```

Output:

```
 Merged Collection: [Dataset A]
                    [Dataset B]
                    [Dataset X]
                    [Dataset Y]
```

Here if two collection have identical dataset names, a dataset is chosen from the *first* collection.


#### Keep last instance

Input:

```
 Collection 1: [Dataset A]
               [Dataset B]
               [Dataset X]
 Collection 2: [Dataset A]
               [Dataset B]
               [Dataset Y]
```

Output:

```
 Merged Collection: [Dataset A]
                    [Dataset B]
                    [Dataset X]
                    [Dataset Y]
```

Here if two collection have identical dataset names, a dataset is chosen from the *last* collection.


#### Append suffix to conflicted element identifiers

Input:

```
 Collection 1: [Dataset A]
               [Dataset B]
               [Dataset X]
 Collection 2: [Dataset A]
               [Dataset B]
               [Dataset Y]
```

Output:

```
 Merged Collection: [Dataset A_1]
                    [Dataset B_1]
                    [Dataset A_2]
                    [Dataset B_2]
                    [Dataset X]
                    [Dataset Y]
```


#### Append suffix to conflicted element identifiers after first on encountered

Input:

```

 Collection 1: [Dataset A]
               [Dataset B]
               [Dataset X]
 Collection 2: [Dataset A]
               [Dataset B]
               [Dataset Y]
```

Output:

```
 Merged Collection: [Dataset A]
                    [Dataset B]
                    [Dataset A_2]
                    [Dataset B_2]
                    [Dataset X]
                    [Dataset Y]
```

#### Append suffix to every element identifier

Input:

```
 Collection 1: [Dataset A]
               [Dataset B]
               [Dataset X]
 Collection 2: [Dataset A]
               [Dataset B]
               [Dataset Y]
```

Output:

```
 Merged Collection: [Dataset A_1]
                    [Dataset B_2]
                    [Dataset A_2]
                    [Dataset B_2]
                    [Dataset X_1]
                    [Dataset Y_2]
```

#### Fail collection creation

This option will simply trigger an error.

### Zip collection

{% tool Zip collection %} takes two collections and creates a paired collection from them.

If you have one collection containing only forward reads and one containing only reverse, this tools will "zip" them together into a simple paired collection. For example, given two collections with `forward` and `reverse` reads they can be "zipped" into a single paired collection:

![Zip collection](../../images/collections/zip.svg)

### Unzip collection

{% tool Unzip collection %} takes a paired collection and "unzips" it into two simple dataset collections (lists of datasets).

Given a paired collection of forward and reverse reads this tool will "unzip" it into two collections containing forward and reverse reads, respectively:

![Unzip collection](../../images/collections/unzip.svg)

## Tools that combine elements of a collection

### Column join

{% tool Column join %} merges elements of a collection on a given column.

If you have a collection with three elements (image below), merging it on the first column will first produce a union on values found in the first column of each elements and then paste elements having the same value side-by-side:

![Column join](../../images/collections/join_on_column.svg)

### Collapse collection

{% tool Collapse collection %} merges elements together (head-to-tail) in the order of the collection. Its power comes from the ability to add identifiers when it performs the merge. Identifiers can be added in variety of ways specified by the **Prepend File name** option as shown in the figure below (we used option **A** in the last step of this tutorial). **A** = `Same line and each line in dataset`; **B** = `Same line and only once per dataset`; **C** =  `Line above`

![Collapse collection](../../images/collections/collapse_collection.svg)

<!-- GTN:IGNORE:002 -->
