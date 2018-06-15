---
layout: tutorial_hands_on
topic_name: variant-analysis
tutorial_name: snpeff-variant-analysis
---

# Introduction
{:.no_toc}

In this tutorial we are going to use data from a study on antibiotic resistance ([Barbosa et al. 2017](http://academic.oup.com/mbe/article/34/9/2229/3829862/Alternative-Evolutionary-Paths-to-Bacterial)).
The study grew a strain of Pseudomonas aeruginosa grown in media containing different antibiotics to study apparition of resistance, its impact on growth and phenomena of cross-resistance. The data are composed of sequencing of 83 samples for 10 conditions:
 - Control samples (con): 10 samples
 - Wild Type (wt): 1 sample
 - Ciprofloxacin (cip): 10 samples
 - Gentamicin (gen): 10 samples
 - Streptomycin (str): 10 samples
 - Piperacillin (pip): 10 samples
 - Carbenicillin (car): 10 samples
 - Doripenem (dor): 10 samples
 - Imipenem (imi): 2 samples
 - Cefsulodin (cef): 10 samples




> ### Agenda
>
> In this tutorial, we will deal with:
>
>
> {:toc}
>
{: .agenda}

# Get the data

## Create reads collection

First we are going to use the SRA uploading tool to upload a collection of paired end reads in Galaxy.

> ### {% icon hands_on %} Hands-on: Get the read data
>
> 1. Open [bioproject of the study](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA355367)
> 2. In *Related information* click on SRA, you have now access to each sample of the experiment.
> 3. For a more interactive access, click on _Send Result to RunSelector_ on the top of the page
> 4. In the download section, click on _Runinfo Table_ and _Accession List_
> 5. Upload the two files in Galaxy
> 6. Run the tool **Download and Extract Reads in FASTA/Q format from NCBI SRA** with parameters as follows :
>   - *select input type* as *List of SRA accession, one per line*
>   - *sra accession list* as the Accession List file we donloaded from NCBI
{: .hands_on}

The downlad of the data can take a little time, and will create two collections, one of paired end data, one of single end data. When you click on the Sungle end collection, you can notice it is empty. To keep the history as clean as possible you can delete it by clicking on the cross.

When you click on the paired end collection you can notice it contains a list of pairs of files (forward and reverse), labeled with samples IDs. The names are not very informative, and in order to facilitate the analysis we want to change the names so it includes the condition of growth of the sample.

> ### {% icon hands_on %} Hands-on: Change the collections names
>We are goin go use a tool to replace the old names by new ones including the condition information. This starts by creating a file containing the old and new names.
>
> 1. Run the tool **Cut columns from a table (cut)** with parameters as follows:
>   - *File to cut* set as the *RunInfo Table* we downloaded from NCBI
>   - *List of Fields* set as *column:1* and *column:3* , these columns contain the sample ID and the condition
> 2.  Run the tool **Add column to an existing dataset** with parameters as follows:
>   - Set  *Add this value* with `_` (This will serve as a separator in the new names)
>   - Set *to Dataset* to the file we generated at the previous step
> 3. Run the tool **Merge Columns together** with parameters as follows:
{: .hands_on}

## Get reference genome data

Once we got a clean collection with our sequencing data, we need to upload the genome and annotations file of the reference genome.
You can find these files on the [Pseudomonas Genome Database Website](http://pseudomonas.com/strain/download). Download the GEnomic sequence file and the gbk annotation file for *Pseudomonas aeruginosa UCBPP-PA14* (second line) and upload them in Galaxy.

> ### {% icon tip %} Tip: Upload file in Galaxy through their url
>
>![Tool search](../../images/tool_search.png "Use search box to find tools!")
{: .tip}


# Align the reads to the reference genome

# Search variations

## Plot coverage to identify structural variations

## Identify and annotation small variations

# Identify relevant variants

## Hierarchical clustering

## Functional analysis


# Conclusion
{:.no_toc}
