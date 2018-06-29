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
> In this tutorial, you will learn perform the annotation and analysis of variants:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get the data

## Create reads collection

First we are going to use the SRA uploading tool to upload a collection of paired end reads in Galaxy.

> ### {% icon hands_on %} Hands-on: Get the read data
>
> 1. Open the [bioproject of the study](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA355367)
> 2. In *Related information* click on SRA, you have now access to each sample of the experiment.
> 3. For a more interactive access, click on _Send Result to RunSelector_ on the top of the page
> 4. In the download section, click on _Runinfo Table_ and _Accession List_
> 5. Upload the two files in Galaxy
> 6. **Download and Extract Reads in FASTA/Q format from NCBI SRA** {% icon tool %} : download SRA files with
>   - {% icon param-text %} *"select input type"* : `List of SRA accession, one per line`
>   - {% icon param-text %} *"sra accession list"* : the `Accession List` file we downloaded from NCBI
{: .hands_on}

The download of the data can take a little time, and will create two collections, one of paired end data, one of single end data. When you click on the single end collection, you can notice it is empty. To keep the history as clean as possible you can delete it by clicking on the cross.

When you click on the paired end collection you can notice it contains a list of pairs of files (forward and reverse), labeled with samples IDs. The names are not very informative, and in order to facilitate the analysis we want to change the names so it includes the condition of growth of the sample.

> ### {% icon hands_on %} Hands-on: Change the collection names
>We are goin go use a tool to replace the old names by new ones including the condition information. This starts by creating a file containing the old and new names.
>
> 1. **Cut columns from a table (cut)** {% icon tool %} : Remove columns with
>   - {% icon param-text %} *"File to cut"* : the `RunInfo Table` we downloaded from NCBI
>   - {% icon param-text %} *"List of Fields"* : `column:9` and `column:10` , these columns contain the sample ID and the condition
> 2. **Add column to an existing dataset** {% icon tool %} : Add column with
>   - {% icon param-text %}  *"Add this value"* : `_` (This will serve as a separator in the new names)
>   - {% icon param-text %}  *"to Dataset"* : the file we generated at the previous step
> 3. **Replace Text in a specific column** {% icon tool %} : Replace text with
>   - {% icon param-text %}  *"File to Process"* : the file we generated at the previous steps
>   - {% icon param-text %}  *"in column"* : `Column:2`
>   - {% icon param-text %}  *"Find Pattern"* : `(.{3}).*` to select the first three letters corresponding to the condition of growth
>   - {% icon param-text %}  *"Replace with"* : `\\1` to replace the whole string with only the condition name
> 4. **Merge Columns together**  {% icon tool %} : Merge columns with
>   - {% icon param-text %} *"Select data"* : the file we just generated
>   - {% icon param-text %} *"Merge column"* : `Column:1`
>   - {% icon param-text %} *"with column"*  to `Column:3`
>   - {% icon param-text %} *"Insert column"* :
>   >    *  {% icon param-text %}  *"Add column"* : `Column:2`
> 5. **Cut columns from a table (cut)**  {% icon tool %} : Remove column with
>   - {% icon param-text %} *"File to cut"* : the file we generated at the previous step
>   - {% icon param-text %} *"List of Fields"* as `column:1` and `column:4` , these columns contain old and new names of our files
> 6. **Relabel List Identifiers from contents of a file**  {% icon tool %} : Rename the files in your collectioni with
>   - {% icon param-text %} *"Input Collection"* : the collection of paired-end data downloaded from SRA accessions
>   - {% icon param-text %} *"How should the new labels be specified?"* as `Maps original identifiers to new ones using a two column table`
>   - {% icon param-text %} *New identifiers* : the two column datasets we generated at the previous step
> 7. Rename your collection with a meaningful name, for example `Pa14 Experiment paired Reads`
{: .hands_on}

## Get reference genome data

Once we got a clean collection with our sequencing data, we need to upload the genome and annotations file of the reference genome.
You can find these files on the [Pseudomonas Genome Database Website](http://pseudomonas.com/strain/download). Download the Genomic sequence file and the gbk annotation file for *Pseudomonas aeruginosa UCBPP-PA14* (second line) and upload them in Galaxy.

> ### {% icon tip %} Tip: Importing data via links
>
> 1. Copy the link location
> 2. Open the Galaxy Upload Manager
> 3. Select **Paste/Fetch Data**
> 4. Paste the link into the text field
> 5. Press **Start**
> 6. Click on **Edit attribute** to give more meaningful names to the datasets.
{: .tip}

Now that you have the reference files in your history we can start the analysis by aligning the reads to the reference genome.
This will allow us to assign a position for each reads and therefore detect where variations occurs.

# Align the reads to the reference genome

## Quality Control

Before we can align the reads on the reference genome, we need to check the quality and the average length of the reads.

> ### {% icon hands_on %} Hands-on: Quality Control
>We are going to use two tools to perform quality control on our reads : FastQC that is going to generate one result per file and MultiQC to aggregate all the results
>
> 1. **FastQC Read Quality reports**  {% icon tool %} : Evaluate the quality of your reads with
>   - {% icon param-text %}  *"Short Read data from your history"* : the collection {% icon param-files %} of paired-end data
> 2. **Flatten Collection into a flat list of datasets** {% icon tool %} : MultiQC need a simple collection as an input, we therefore need to flatten our list of pairs of files with
>   - {% icon param-text %}  *"Input Collection"* : the raw data output of Fastqc `FastQC on collection [...]:RawData`
> 3.  **MultiQC** {% icon tool %} : Aggregate quality control outputs with :
>   - {% icon param-text %}  *Which tool was used generate logs?* : `FastQC`
>   - {% icon param-text %}  *FastQC output* : the flatten collection we generated at the previous step
{: .hands_on}

MultiQC provides two outputs : a Webpage aggregating all the results from FastQC and a collection containing statistics in the text format.
As we said earlier, we are looking for two informations here: the quality and the average length of reads.

You can find the read quality on the web page output.

![MultiQC Quality](../../images/multiqc_quality.png "MultiQC Output shows read quality")

We can see there that the quality of the reads are good over the whole length of the reads. We therefore don't need to perform quality treatment and we can use them as they are.

We can find the average read length in the file *general_stats* of the stats file collection output of MultiQC.


> ### {% icon question %} Questions
>
> 1. What is the average read length(rounded to the closest integer)?
> 2. MultiQC grouped our 80 samples in two samples, based on which criteria?
> 3. Can that grouping be a problem?  
>    > ### {% icon solution %} Solution
>    >
>    > 1. The average read length is 101 bases.
>    > 2. MultiQC grouped the samples in *forward* and  *reverse* samples.
>    > 3. We are applying the same treatment to all our samples, so we are only interested in the average quality accross all samples. Therefore grouping them is not a problem.
>    >
>    {: .solution}
{: .question}


## Mapping

Now that we have good quality reads we are going to map them against the reference genome.

> ### {% icon hands_on %} Hands-on: Mapping
>Now that we have good quality reads we are going to map them against the reference genome.
> 1.   **BWA-MEM**  {% icon tool %} : Map your reads against the reference with
 
>   -  
{: .hands_on}

# Search variations

## Plot coverage to identify structural variations

## Identify and annotation small variations

# Identify relevant variants

## Hierarchical clustering

## Functional analysis


# Conclusion
{:.no_toc}
