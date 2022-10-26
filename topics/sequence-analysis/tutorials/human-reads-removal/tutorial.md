---
layout: tutorial_hands_on

title: "Removal of human reads from SARS-CoV-2 sequencing data"
questions:
  - How can you remove traces of human sequence data from a sequenced viral sample?
objectives:
  - Obtain viral (SARS-CoV-2) sequencing data with contaminating human reads from public sources
  - Organize the data into a collection
  - Preprocess and map the data against the human reference genome
  - Eliminate reads/read pairs that map to the human genome from the original data
time_estimation: "1h"
level: Intermediate
key_points:
  - Before submitting raw viral sequencing data to public databases you will want to remove human sequence traces
  - Human reads can be identified by mapping the data to the human reference genome
  - After mapping, you can extract the read identifiers of the non-human reads and use these to extract just the desired original sequenced reads from the input data
requirements:
  -
    type: "internal"
    topic_name: sequence-analysis
    tutorials:
      - quality-control
      - mapping
  -
    type: "internal"
    topic_name: galaxy-interface
    tutorials:
      - collections
tags:
  - covid19
contributors:
  - wm75

---

# Introduction


Patient samples for pathogen detection are usually "contaminated" with human host DNA.
Such contamination, if not removed from sequencing data, may pose an issue with certain types of
data analyses.
Another issue is that submitting the raw sequenced reads with the human sequence "contamination" included to a public database may simply violate national or institutional regulations for handling patient data.
For this second reason in particular, researchers will regularly have to remove traces of human sequencing data even from samples that, by their nature, should be highly enriched for just pathogen reads, as is the case, for example, for ampliconic viral samples, in which viral sequences have been PCR-amplified with virus-specific primers before sequencing.

This tutorial will guide you through a typical workflow for clearing human sequences from any kind (ampliconic or not) of viral sequenced sample, which retains non-human reads in a format ready to be submitted to public databases.

> <comment-title>Nature of the input data</comment-title>
> We will use sequencing data of bronchoalveolar lavage fluid (BALF) samples obtained from early COVID-19 patients in China as our input data.
> Since such samples are expected to be contaminated signficantly with human sequenced reads, the effect of the cleaning steps will be much more apparent than for ampliconic sequencing data from diagnostic patient swabs.
> Nevertheless the cleaning processs does not rely on any particular sample isolation or preprcocessing method and could be used unaltered on, for example, ARTIC amplified SARS-CoV-2 samples.
>
{: .comment}


> <agenda-title></agenda-title>
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get data

As always, it is best to give each analysis you are performing with Galaxy its own history, so let's start with preparing this:

> <hands-on-title>Prepare new history</hands-on-title>
>
> 1. Create a new history for this tutorial and give it a proper name
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>    {% snippet faqs/galaxy/histories_rename.md %}
>
{: .hands_on}

## Getting the sequencing data

The sequenced reads datasets used in this tutorial (4 files representing 2 Illumina paired-end sequenced samples) have been deposited on [Zenodo](https://zenodo.org/record/3732359/) and can be uploaded to Galaxy via their URLs.
After that, we will arrange the uploaded data into a collection in our Galaxy history to facilitate its handling in the analysis.

> <hands-on-title>Data upload and rearrangement</hands-on-title>
>
> 1. Upload the data from Zenodo via URLs
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    The URLs for our example data are these:
>
>    ```
>    https://zenodo.org/record/3732359/files/SRR10903401_r1.fq.gz
>    https://zenodo.org/record/3732359/files/SRR10903401_r2.fq.gz
>    https://zenodo.org/record/3732359/files/SRR10903402_r1.fq.gz
>    https://zenodo.org/record/3732359/files/SRR10903402_r2.fq.gz
>    ```
>
> 2. Arrange the data into a paired dataset collection
>
>    {% snippet faqs/galaxy/collections_build_list_paired.md %}
>
>    For the example datasets this means:
>    - You need to tell Galaxy about the suffix for your forward and reverse reads, respectively:
>      - set the text of *unpaired forward* to: `_r1.fq.gz`
>      - set the text of *unpaired reverse* to: `_r2.fq.gz`
>      - click: `Auto-pair`
>
>      All datasets should now be moved to the *paired section* of the dialog, and the middle column there should show that only the sample accession numbers, *i.e.* `SRR10903401` and `SRR10903402`, will be used as the pair names.
>
>    - Make sure *Hide original elements* is checked to obtain a cleaned-up history after building the collection.
>    - Click *Create Collection*
>
{: .hands_on}

# Read trimming and mapping

In order to learn which of the reads of each of the input samples are of likely human origin, we need to map the reads to the human reference genome. To improve the quality of this step we will also trim low-quality bases from the ends of all reads before passing them to the read mapping software.

> <comment-title></comment-title>
>
> Do you want to learn more about the principles behind quality control (including trimming) and mapping? Follow our dedicated tutorials:
> - [Quality Control]({% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %}).
> - [Mapping]({% link topics/sequence-analysis/tutorials/mapping/tutorial.md %}).
{: .comment}

> <hands-on-title>Trim reads and map them to the human reference genome</hands-on-title>
>
> 1. {% tool [Trimmomatic](toolshed.g2.bx.psu.edu/repos/pjbriggs/trimmomatic/trimmomatic/0.38.0) %} with default settings and:
>    - *"Single-end or paired-end reads?"*: `Paired-end (as collection)`
>      - *"Select FASTQ dataset collection with R1/R2 pair"*: the collection of sequenced reads input datasets
> 2. {% tool [Map with BWA-MEM](toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.17.2) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a built-in genome index`
>      - *"Using reference genome"*: `Human Dec. 2013 (GRCh38/hg38)(hg38)`
>
>        The exact name of this version of the human reference genome may vary between Galaxy servers. Start typing `hg38` in the select box to reveal available choices on your server.
>
>    - *"Single or Paired-end reads"*: `Paired Collection`
>      - *"Select a paired collection"*: the `paired` output of **Trimmomatic**
>
{: .hands_on}

# Obtaining the identifiers of non-human read pairs

The mapping step above has produced a collection of BAM datasets that store, for each sample, which of its reads map to the human genome (and if so to which position on the genome though we are not interested in this information), and which ones could not be mapped to it. It is these *unmapped* reads that we are interested in, but how do we get at them?

The logic of the following steps is:

1. Filter the BAM dataset for each sample for only those reads that could not be mapped to the human reference genome, and for which also the read mate in the pair could not be mapped to this genome.

   In other words, if only one read of a read pair can be mapped to the human genome we play it safe and discard the pair.

2. Emit the remaining reads in a format (we choose FASTA here) from which it is easy to extract the read identifiers

3. Extract the read identifiers

> <hands-on-title>Filter for read pairs not mapping to the human genome and extract their identifiers</hands-on-title>
>
> 1. {% tool [Samtools fastx](toolshed.g2.bx.psu.edu/repos/iuc/samtools_fastx/samtools_fastx/1.9+galaxy1) %} with the following parameter settings:
>    - {% icon param-collection %} *"BAM or SAM file to convert"*: the output of **Map with BWA-MEM**
>    - *"Output format"*: `FASTA`
>    - *"outputs"*: select `READ1`
>
>      With this setting we are choosing to retrieve only the forward reads of each pair as a single result dataset.
>      Since the forward and reverse reads of a pair share their read identifier, we do not need the reverse reads.
>
>    - *"Omit or append read numbers"*: `Do not append /1 and /2 to read names`
>
>      Adding this first or second read of the pair indicator to the read name would make the result dependent on which half of the reads we retrieve. We do not want this.
>
>    - *"Require that these flags are set"*: 
>      - `read is unmapped`
>      - `mate is unmapped`
>
>      This setting makes sure we retain only read pairs for which none of the two reads could be mapped to the human genome.
>
> 2. {% tool [Select lines that match an expression](Grep1) %} with the following parameters:
>    - {% icon param-collection %} *"Select lines from"*: the output of **Samtools fastx**
>    - *"that"*: `Matching`
>    - *"the pattern"*: `^>.+`
>
>    This will select only those lines from the FASTA inputs that start with a `>` character, *i.e.* will retain the identifier lines, but not the sequence lines.
>
> 3. {% tool [Replace Text in entire line](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_line/1.1.2) %} with the following parameters:
>    - {% icon param-collection %} *"File to process"*: the output of **Select**
>    - {% icon param-repeat %} *"Replacement"*:
>      - *"Find Pattern"*: `^>(.+)`
>      - *"Replace with"*: `\1`
>
>      This find and replacement pattern will work on lines starting with a `>` character and replace each such line with the content following the `>` character (the `(.+)` part that gets referenced as `\1` in the replacement pattern).
>
{: .hands_on}

# Use the non-human read identifiers to extract the reads of interest from the original inputs

We are almost there! We already know the identifiers of all reads we would like to retain, we only need a way to filter the original input datasets for these identifiers.

> <comment-title>Why not use the unmapped reads directly?</comment-title>
>
> When you followed the previous steps you may have noticed that **Samtools fastx** offers *fastq* as an output format.
> Since we used this tool to report the unmapped reads, you may be wondering why we did not simply ask the tool to write the reads in the format we want and be done with it.
>
> The answer is: we could have done that (and some people follow this approach), but remember that we trimmed the reads before mapping them, which means the extracted unmapped reads would not be identical to the input reads.
>
> It is good practice to submit sequencing data to public databases in as raw a state as possible to leave it up to the people downloading that data to apply processing steps as they see fit.
>
{: .comment}

The **seqtk_subseq** tool offfers the identifier-based filtering functionality we are looking for.
As a slight complication, however, the tool is not prepared to handle paired read collections like the one we organized our input data into.
For this reason, we need to first unzip our collection into two simple list collections - one with the forward reads, the other one with the reverse reads of all samples.
Then, after processing both with **seqtk_subseq**, we zip the filtered collections back into one.

> <hands-on-title>Identifier-based extraction of non-human reads from the input data</hands-on-title>
>
> 1. {% tool [Unzip Collection](__UNZIP_COLLECTION__) %} with the following setting:
>    - *"Input Paired Dataset"*: the collection of original sequenced reads input datasets
>
> 2. {% tool [seqtk_subseq](toolshed.g2.bx.psu.edu/repos/iuc/seqtk/seqtk_subseq/1.3.1) %} with the following parameters:
>    - {% icon param-collection %} *"Input FASTA/Q file"*: the first (forward reads) output of **Unzip Collection**
>    - *"Select source of sequence choices"*: `FASTA/Q ID list`
>      - {% icon param-collection %} the collection of read identifiers; output of **Replace Text**
>
> 3. {% tool [seqtk_subseq](toolshed.g2.bx.psu.edu/repos/iuc/seqtk/seqtk_subseq/1.3.1) %} with the following parameters:
>    - {% icon param-collection %} *"Input FASTA/Q file"*: the second (reverse reads) output of **Unzip Collection**
>    - *"Select source of sequence choices"*: `FASTA/Q ID list`
>      - {% icon param-collection %} the collection of read identifiers; output of **Replace Text**
>
> 4. {% tool [Zip Collection](__ZIP_COLLECTION__) %} with the following settings:
>    - {% icon param-collection %} *"Input Dataset (Forward)"*: the collection of cleaned forward reads; output of the first run of **seqtk_subseq**
>    - {% icon param-collection %} *"Input Dataset (Reverse)"*: the collection of cleaned reverse reads; output of the second run of **seqtk_subseq**
>
{: .hands_on}

# Conclusion


Data cleaning is a standard procedure for clinical sequencing datasets and is a rather straightforward process in Galaxy.
You just learnt how to remove human reads from any number of paired-end sequenced samples in parallel using collections and just a handful of tools.
Cleaning of single-end data would use the same tools (with adjusted settings), but would work on data arranged in a simple list collection instead of in a paired collection.
