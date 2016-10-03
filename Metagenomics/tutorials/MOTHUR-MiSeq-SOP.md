Metagenomics with Mothur: MiSeq SOP
==================================

:grey_question: ***Questions***

- *What is the effect of normal variation in the gut microbiome on host health?*

:dart: ***Objectives***

- *Learn how to analyse 16S rRNA sequences in Galaxy using the Mothur toolsuite*
- *Learn how to use dataset collections to process a large number of samples at once*

:heavy_check_mark: ***Requirements***

- *Galaxy introduction*

:hourglass: ***Time estimation*** *1d/3h/6h*

# Introduction

This tutorial will demonstrate how to perform the *standard operating procedure (SOP)* for the analysis of 16S rRNA gene sequences generated using Illumina's MiSeq platform, with the [Mothur toolsuite](http://www.mothur.org/wiki) within Galaxy. This SOP was developed by the Schloss lab and described [here](http://www.mothur.org/wiki/MiSeq_SOP) on the mother wiki.

<!-- TODO: add citation
Kozich JJ, Westcott SL, Baxter NT, Highlander SK, Schloss PD. (2013): Development of a dual-index sequencing strategy and curation pipeline for analyzing amplicon sequence data on the MiSeq Illumina sequencing platform. Applied and Environmental Microbiology. 79(17):5112-20. -->

## Overview
In this tutorial we will perform the following steps:

1. Obtaining and preparing input data
2. Quality Control
3. OTU-based analysis
4. Phylotype-based analysis
5. Phylogeny-based analysis


# Part 1: Obtaining and preparing data

## Understanding our input data
In this tutorial we are interested in understanding the effect of normal variation in the gut microbiome on host health. To that end fresh feces from mice were collected on a daily basis for 365 days post weaning. During the first 150 days post weaning (dpw), nothing was done to our mice except allow them to eat, get fat, and be merry. We were curious whether the rapid change in weight observed during the first 10 dpw affected the stability microbiome compared to the microbiome observed between days 140 and 150. We will address this question in this tutorial using a combination of OTU, phylotype, and phylogenetic methods. To make this tutorial easier to execute, we are providing only part of the data - you are given the flow files for one animal at 10 time points (5 early and 5 late). In addition, to sequencing samples from mice fecal material, we resequenced a mock community composed of genomic DNA from 21 bacterial strains. We will use the 10 fecal samples to look at how to analyze microbial communities and the mock community to measure the error rate and its effect on other analyses.

:pencil2: ***Hands on!***

So now that we know what our input data is, let's get it into our history:

1. Create a **new history** and name it "Mothur MiSeq SOP"

2. **Import** the training data **to your history**. There are two ways to do this. The easiest is using the data available from a *shared data library*, if this is not possible you can download the data yourself and upload it to your Galaxy instance.
  - From data libray:
      - Navigate to the shared data library named *Galaxy training: Mothur MiSeq SOP* and import all fastq files you encounter there
  - From your computer:
      - obtain data directly from [here](http://www.mothur.org/w/images/d/d6/MiSeqSOPData.zip)  <!-- TODO: zenodo link-->
      - unzip it
      - upload all fastq files to your history

3. Create a **paired collection**. Since we have paired end data, each sample consist of two seperate fastq files, one containing the forward reads, and one containing the reverse reads. We can recognize the pairing from the file names, they will differ only by `_R1` or `_R2` in the filename. We can tell Galaxy about this pairing naming convention, so that our tools will know which files belong together.

 - Click on the **checkmark icon** at top of your history
 - Select all fastq files, then click on **for all selected..** and select **Build list of dataset pairs** from the dropdown menu.
 - In the next dialog window you can create the list of pairs. By default Galaxy will look for pairs of files that differ only by a `_1` and `_2` part in their names. In our case however, these should be `_R1` and `_R2`. Please change these values accordingly. You should now see a list of pairs detected by galaxy, examine the file names, if it looks good, you can click on **auto-pair** to create the suggested pairs.
 - Enter a name for your new collection at the bottom right of the screen.
 - Click the **Create List** button. A new dataset collection item should now appear in your history.

TODO: screenshots for collection creation

## Reference Data

Apart from our input data we will also need some additional reference data.

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

Some blabla

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

# Part 2

Short introduction about this subpart.

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

## Subpart 2

Short introduction about this subpart.

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

Some blabla

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

# Conclusion

Conclusion about the technical key points. And then relation between the technics and the biological question to end with a global view.

:grey_exclamation: ***Key Points***

- *Simple sentence to sum up the first key point of the tutorial (Take home message)*
- *Second key point*
- *Third key point*
- *...*

# :clap: Thank you
