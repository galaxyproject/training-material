---
layout: tutorial_hands_on
topic_name: transcriptomics
tutorial_name: clipseq
---

# Introduction
{:.no_toc}

<!-- This is a comment. -->

General introduction about the topic and then an introduction of the
tutorial (the questions and the objectives). It is nice also to have a
scheme to sum up the pipeline used during the tutorial. The idea is to
give to trainees insight into the content of the tutorial and the (theoretical
and technical) key concepts they will learn.

**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Title for your first section

Give some background about what the trainees will be doing in the section.

Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Import the following files from [Zenodo](https://zenodo.org/record/1327423) or from a data
>    library named `TODO` if available (ask your instructor)
>
>    ```
>    https://zenodo.org/api/files/102d29d5-2180-490b-be7c-bb0e4ca7b109/hg19_chr_sizes.txt
>    https://zenodo.org/api/files/102d29d5-2180-490b-be7c-bb0e4ca7b109/Homo_sapiens.GRCh37.74.gtf
>    https://zenodo.org/api/files/102d29d5-2180-490b-be7c-bb0e4ca7b109/RBFOX2-204-CLIP_S1_R1_RBFOX2.fastq
>    https://zenodo.org/api/files/102d29d5-2180-490b-be7c-bb0e4ca7b109/RBFOX2-204-CLIP_S1_R2_RBFOX2.fastq
>    https://zenodo.org/api/files/102d29d5-2180-490b-be7c-bb0e4ca7b109/RBFOX2-204-INPUT_S2_R1.fastq
>    https://zenodo.org/api/files/102d29d5-2180-490b-be7c-bb0e4ca7b109/RBFOX2-204-INPUT_S2_R2.fastq
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    > ### {% icon tip %} Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * Press **Start**
>    >
>    > By default, Galaxy uses the url as the name, so please rename them to something more pleasing.
>    {: .tip}
>
>    > ### {% icon tip %} Tip: Importing data from a data library
>    >
>    > * Go into "Shared data" (top panel) then "Data libraries"
>    > * Click on "Training data" and then "Transcriptomics"
>    > * Select interesting file
>    > * Click on "Import selected datasets into history"
>    > * Import in a new history
>    {: .tip}
>
{: .hands_on}

# Title of the section usually corresponding to a big step in the analysis

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


## Sub-step with **Unzip Collection**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Unzip Collection** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"Input Paired Dataset"*: `output` (Input dataset collection)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **FastQC**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **FastQC** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"Short read data from your current history"*: `output` (Input dataset collection)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Unzip Collection**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Unzip Collection** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"Input Paired Dataset"*: `output` (Input dataset collection)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **FastQC**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **FastQC** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"Short read data from your current history"*: `output` (Input dataset collection)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Cutadapt**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Cutadapt** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Fastq file to trim"*: `forward` (output of **Unzip Collection** {% icon tool %})
>    - *"Track Paired Reads"*: `Yes`
>    - In *"3' Adapters"*:
>        - Click on *"Insert 3' Adapters"*:
>        - In *"1: 3' Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 3' adapter sequence"*: `AACTTGTAGATCGGA`
>        - Click on *"Insert 3' Adapters"*:
>        - In *"2: 3' Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 3' adapter sequence"*: `AGGACCAAGATCGGA`
>    - In *"5' (Front) Adapters"*:
>        - Click on *"Insert 5' (Front) Adapters"*:
>        - In *"1: 5' (Front) Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 5' adapter sequence"*: `CTTCCGATCTACAAGTT`
>        - Click on *"Insert 5' (Front) Adapters"*:
>        - In *"2: 5' (Front) Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 5' adapter sequence"*: `CTTCCGATCTTGGTCCT`
>    - *"Maximum error rate"*: `"0.1"`
>    - *"Do not allow indels (Use ONLY with anchored 5' (front) adapters)."*: `Yes`
>    - *"Match times"*: `"1"`
>    - *"Minimum overlap length"*: `"5"`
>    - *"Match Read Wildcards"*: `Yes`
>    - *"Output filtering options"*: `Set Filters`
>        - *"Minimum length"*: `10`
>    - *"Additional output options"*: `Default`
>    - *"Additional modifications to reads"*: `Set Modification Options`
>        - *"Cut bases from reads before adapter trimming"*: `-5`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Cutadapt**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Cutadapt** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Fastq file to trim"*: `forward` (output of **Unzip Collection** {% icon tool %})
>    - *"Track Paired Reads"*: `Yes`
>    - In *"3' Adapters"*:
>        - Click on *"Insert 3' Adapters"*:
>        - In *"1: 3' Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 3' adapter sequence"*: `AACTTGTAGATCGGA`
>        - Click on *"Insert 3' Adapters"*:
>        - In *"2: 3' Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 3' adapter sequence"*: `AGGACCAAGATCGGA`
>    - In *"5' (Front) Adapters"*:
>        - Click on *"Insert 5' (Front) Adapters"*:
>        - In *"1: 5' (Front) Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 5' adapter sequence"*: `CTTCCGATCTACAAGTT`
>        - Click on *"Insert 5' (Front) Adapters"*:
>        - In *"2: 5' (Front) Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 5' adapter sequence"*: `CTTCCGATCTTGGTCCT`
>    - *"Maximum error rate"*: `"0.1"`
>    - *"Do not allow indels (Use ONLY with anchored 5' (front) adapters)."*: `Yes`
>    - *"Match times"*: `"1"`
>    - *"Minimum overlap length"*: `"5"`
>    - *"Match Read Wildcards"*: `Yes`
>    - *"Output filtering options"*: `Set Filters`
>        - *"Minimum length"*: `10`
>    - *"Additional output options"*: `Default`
>    - *"Additional modifications to reads"*: `Set Modification Options`
>        - *"Cut bases from reads before adapter trimming"*: `-5`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Cutadapt**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Cutadapt** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Fastq file to trim"*: `paired_output` (output of **Cutadapt** {% icon tool %})
>    - *"Track Paired Reads"*: `Yes`
>    - In *"3' Adapters"*:
>        - Click on *"Insert 3' Adapters"*:
>        - In *"1: 3' Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 3' adapter sequence"*: `AACTTGTAGATCGGA`
>        - Click on *"Insert 3' Adapters"*:
>        - In *"2: 3' Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 3' adapter sequence"*: `AGGACCAAGATCGGA`
>    - In *"5' (Front) Adapters"*:
>        - Click on *"Insert 5' (Front) Adapters"*:
>        - In *"1: 5' (Front) Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 5' adapter sequence"*: `CTTCCGATCTACAAGTT`
>        - Click on *"Insert 5' (Front) Adapters"*:
>        - In *"2: 5' (Front) Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 5' adapter sequence"*: `CTTCCGATCTTGGTCCT`
>    - *"Maximum error rate"*: `"0.1"`
>    - *"Do not allow indels (Use ONLY with anchored 5' (front) adapters)."*: `Yes`
>    - *"Match times"*: `"1"`
>    - *"Minimum overlap length"*: `"5"`
>    - *"Match Read Wildcards"*: `Yes`
>    - *"Output filtering options"*: `Set Filters`
>        - *"Minimum length"*: `10`
>    - *"Additional output options"*: `Default`
>    - *"Additional modifications to reads"*: `No Read Modifications`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Cutadapt**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Cutadapt** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Fastq file to trim"*: `paired_output` (output of **Cutadapt** {% icon tool %})
>    - *"Track Paired Reads"*: `Yes`
>    - In *"3' Adapters"*:
>        - Click on *"Insert 3' Adapters"*:
>        - In *"1: 3' Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 3' adapter sequence"*: `AACTTGTAGATCGGA`
>        - Click on *"Insert 3' Adapters"*:
>        - In *"2: 3' Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 3' adapter sequence"*: `AGGACCAAGATCGGA`
>    - In *"5' (Front) Adapters"*:
>        - Click on *"Insert 5' (Front) Adapters"*:
>        - In *"1: 5' (Front) Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 5' adapter sequence"*: `CTTCCGATCTACAAGTT`
>        - Click on *"Insert 5' (Front) Adapters"*:
>        - In *"2: 5' (Front) Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 5' adapter sequence"*: `CTTCCGATCTTGGTCCT`
>    - *"Maximum error rate"*: `"0.1"`
>    - *"Do not allow indels (Use ONLY with anchored 5' (front) adapters)."*: `Yes`
>    - *"Match times"*: `"1"`
>    - *"Minimum overlap length"*: `"5"`
>    - *"Match Read Wildcards"*: `Yes`
>    - *"Output filtering options"*: `Set Filters`
>        - *"Minimum length"*: `10`
>    - *"Additional output options"*: `Default`
>    - *"Additional modifications to reads"*: `No Read Modifications`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **UMI-tools extract**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **UMI-tools extract** {% icon tool %} with the following parameters:
>    - *"Library type"*: `Paired-end`
>        - *"Barcode on both reads?"*: `Barcode on first read only`
>    - *"Use Known Barcodes?"*: `No`
>    - *"Method to extract barcodes"*: ``
>    - *"Barcode pattern for first read"*: `"NNNNN"`
>    - *"Is the barcode at the 5' end?"*: `Yes`
>    - *"Output log?"*: `Yes`
>    - *"Enable quality filter?"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **UMI-tools extract**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **UMI-tools extract** {% icon tool %} with the following parameters:
>    - *"Library type"*: `Paired-end`
>        - *"Barcode on both reads?"*: `Barcode on first read only`
>    - *"Use Known Barcodes?"*: `No`
>    - *"Method to extract barcodes"*: ``
>    - *"Barcode pattern for first read"*: `"NNNNN"`
>    - *"Is the barcode at the 5' end?"*: `Yes`
>    - *"Output log?"*: `Yes`
>    - *"Enable quality filter?"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **RNA STAR**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **RNA STAR** {% icon tool %} with the following parameters:
>    - *"Single-end or paired-end reads"*: `Paired-end (as individual datasets)`
>    - *"Custom or built-in reference genome"*: `Use a built-in index`
>        - *"Reference genome with or without an annotation"*: `use genome reference with builtin gene-model`
>            - *"Select reference genome"*: `Homo sapiens (hg19+GRCh37.75)`
>    - *"Count number of reads per gene"*: `Yes`
>    - *"Would you like to set output parameters (formatting and filtering)?"*: `Yes`
>        - *"Extra SAM attributes to include"*: `All`
>        - *"Include strand field flag XS"*: `Yes -- and reads with inconsistent and/or non-canonical introns are filtered out`
>        - *"Would you like to set additional output parameters (formatting and filtering)?"*: `Yes`
>    - *"Other parameters (seed, alignment, limits and chimeric alignment)"*: `Extended parameter list`
>        - In *"Alignment parameters"*:
>            - *"Use end-to-end read alignments, with no soft-clipping?"*: `Yes`
>        - *"Would you like to set chimeric alignment parameters?"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **RNA STAR**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **RNA STAR** {% icon tool %} with the following parameters:
>    - *"Single-end or paired-end reads"*: `Paired-end (as individual datasets)`
>    - *"Custom or built-in reference genome"*: `Use a built-in index`
>        - *"Reference genome with or without an annotation"*: `use genome reference with builtin gene-model`
>            - *"Select reference genome"*: `Homo sapiens (hg19+GRCh37.75)`
>    - *"Count number of reads per gene"*: `Yes`
>    - *"Would you like to set output parameters (formatting and filtering)?"*: `Yes`
>        - *"Extra SAM attributes to include"*: `All`
>        - *"Include strand field flag XS"*: `Yes -- and reads with inconsistent and/or non-canonical introns are filtered out`
>        - *"Would you like to set additional output parameters (formatting and filtering)?"*: `Yes`
>    - *"Other parameters (seed, alignment, limits and chimeric alignment)"*: `Extended parameter list`
>        - In *"Alignment parameters"*:
>            - *"Use end-to-end read alignments, with no soft-clipping?"*: `Yes`
>        - *"Would you like to set chimeric alignment parameters?"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **UMI-tools deduplicate**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **UMI-tools deduplicate** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Reads to deduplicate in SAM or BAM format"*: `mapped_reads` (output of **RNA STAR** {% icon tool %})
>    - *""*: ``
>    - *"Separator between read id and UMI."*: `"_"`
>    - *"Tag which contains UMI."*: `""`
>    - *"Method used to identify PCR duplicates within reads."*: ``
>    - *"Edit distance threshold"*: `"1"`
>    - *"BAM is paired end"*: `Yes`
>    - *"Spliced reads are unique"*: `Yes`
>    - *"Soft clip threshold"*: `"4"`
>    - *"Use the read length as as a criterion when deduping"*: `Yes`
>    - *"Consider all alignments to a single contig together"*: `Yes`
>    - *"Only consider a random selection of the reads"*: `"1.0"`
>    - *"Only consider a single chromosome"*: `Yes`
>    - *"Deduplicate per contig"*: `Yes`
>    - *"Deduplicate per gene"*: `Yes`
>    - *"Deduplicate by this gene tag"*: `""`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **multiBamSummary**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **multiBamSummary** {% icon tool %} with the following parameters:
>    - *"Sample order matters"*: `No`
>    - *"Choose computation mode"*: `Bins`
>    - *"Region of the genome to limit the operation to"*: `""`
>    - *"Show advanced options"*: `no`
>    - *"Save raw counts (coverages) to file"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **plotFingerprint**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **plotFingerprint** {% icon tool %} with the following parameters:
>    - *"Sample order matters"*: `No`
>    - *"Region of the genome to limit the operation to"*: `""`
>    - *"Show advanced options"*: `no`
>    - *"Show advanced output settings"*: `no`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **UMI-tools deduplicate**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **UMI-tools deduplicate** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Reads to deduplicate in SAM or BAM format"*: `mapped_reads` (output of **RNA STAR** {% icon tool %})
>    - *""*: ``
>    - *"Separator between read id and UMI."*: `"_"`
>    - *"Tag which contains UMI."*: `""`
>    - *"Method used to identify PCR duplicates within reads."*: ``
>    - *"Edit distance threshold"*: `"1"`
>    - *"BAM is paired end"*: `Yes`
>    - *"Spliced reads are unique"*: `Yes`
>    - *"Soft clip threshold"*: `"4"`
>    - *"Use the read length as as a criterion when deduping"*: `Yes`
>    - *"Consider all alignments to a single contig together"*: `Yes`
>    - *"Only consider a random selection of the reads"*: `"1.0"`
>    - *"Only consider a single chromosome"*: `Yes`
>    - *"Deduplicate per contig"*: `Yes`
>    - *"Deduplicate per gene"*: `Yes`
>    - *"Deduplicate by this gene tag"*: `""`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **FastQC**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **FastQC** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Short read data from your current history"*: `output` (output of **UMI-tools deduplicate** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **plotCorrelation**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **plotCorrelation** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Matrix file from the multiBamSummary tool"*: `outFile` (output of **multiBamSummary** {% icon tool %})
>    - *"Correlation method"*: ``
>    - *"Plotting type"*: `Heatmap`
>    - *"Skip zeros"*: `Yes`
>    - *"Image file format"*: ``
>    - *"Remove regions with very large counts"*: `Yes`
>    - *"Save the matrix of values underlying the heatmap"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **FastQC**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **FastQC** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Short read data from your current history"*: `output` (output of **UMI-tools deduplicate** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **PEAKachu**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PEAKachu** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Experiment Libraries"*: `output` (output of **UMI-tools deduplicate** {% icon tool %})
>    - {% icon param-file %} *"Control Libraries"*: `output` (output of **UMI-tools deduplicate** {% icon tool %})
>    - *"Pairwise Replicates"*: `Yes`
>    - *"Paired End"*: `Yes`
>    - *"Maximum Insert Size"*: `"200"`
>    - *"Features"*: `""`
>    - *"Sub-Features"*: `""`
>    - *"Select Mode"*: `Adaptive`
>        - *"Normalisation Method."*: `DESeq2`
>    - *"Mad Multiplier"*: `"0.0"`
>    - *"Fold Change Threshold"*: `"2.0"`
>    - *"Adjusted p-value Threshold"*: `"0.05"`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Extract alignment ends**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Extract alignment ends** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Alignments in SAM or BAM format"*: `output` (output of **UMI-tools deduplicate** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Text reformatting**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Text reformatting** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `peak_tables` (output of **PEAKachu** {% icon tool %})
>    - *"AWK Program"*: `"NR>1{\nif ($3 < $4) {\n   print $1,$3,$4,\"clip_peak_\"NR-1,$9,$5;\n}\nelse {\n   print $1,$4,$3,\"clip_peak_\"NR-1,$9,$5;\n}\n}"`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SortBED**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **SortBED** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sort the following BED file"*: `alignment_ends` (output of **Extract alignment ends** {% icon tool %})
>    - *"Sort by"*: ``
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Get crosslinked nucleotides**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Get crosslinked nucleotides** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Alignments in BED format"*: `alignment_ends` (output of **Extract alignment ends** {% icon tool %})
>    - *"Set position one nt downstream of 3'-end as crosslinked nucleotide"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SlopBed**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **SlopBed** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"BED/VCF/GFF file"*: `outfile` (output of **Text reformatting** {% icon tool %})
>    - *"Genome file"*: `Genome file from your history`
>    - *"Define -l and -r as a fraction of the feature’s length"*: `Yes`
>    - *"Define -l and -r based on strand"*: `Yes`
>    - *"Choose what you want to do"*: `Increase the BED/GFF/VCF entry by the same number base pairs in each direction.`
>        - *"Number of base pairs"*: `20`
>    - *"Print the header from the A file prior to results"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Create a BedGraph of genome coverage**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Create a BedGraph of genome coverage** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"The BAM or BED file from which coverage should be computed"*: `output` (output of **SortBED** {% icon tool %})
>    - *"Report regions with zero coverage"*: `Yes`
>    - *"Treat split/spliced BAM or BED12 entries as distinct BED intervals when computing coverage."*: `Yes`
>    - *"Calculate coverage based on"*: ``
>    - *"Scale the coverage by a constant factor"*: `""`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Text reformatting**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Text reformatting** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `crosslinking_coordinates` (output of **Get crosslinked nucleotides** {% icon tool %})
>    - *"AWK Program"*: `"$2!=0"`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Text reformatting**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Text reformatting** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `output` (output of **SlopBed** {% icon tool %})
>    - *"AWK Program"*: `"{print $0}"`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Extract Genomic DNA**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Extract Genomic DNA** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Fetch sequences for intervals in"*: `output` (output of **SlopBed** {% icon tool %})
>    - *"Interpret features when possible"*: ``
>    - *"Choose the source for the reference genome"*: `locally cached`
>        - *"Using reference genome"*: ``
>    - *"Select output format"*: ``
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Wig/BedGraph-to-bigWig**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Wig/BedGraph-to-bigWig** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Convert"*: `output` (output of **Create a BedGraph of genome coverage** {% icon tool %})
>    - *"Converter settings to use"*: `Default`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SortBED**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **SortBED** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sort the following BED file"*: `outfile` (output of **Text reformatting** {% icon tool %})
>    - *"Sort by"*: ``
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **RNA Centric Annotation System**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **RNA Centric Annotation System** {% icon tool %} with the following parameters:
>    - *"Genome Version"*: ``
>    - {% icon param-file %} *"Target regions in BED format"*: `outfile` (output of **Text reformatting** {% icon tool %})
>    - {% icon param-file %} *"Reference annotation in ENSEMBL GTF format"*: `output` (Input dataset)
>    - *"Run annotation."*: `Yes`
>    - *"Run GO term enrichment"*: `Yes`
>    - *"Run gene set enrichment"*: `No`
>    - *"Run motif search"*: `Yes`
>    - *"Downsampling (N)"*: `"0"`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **MEME-ChIP**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **MEME-ChIP** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Primary sequences"*: `output` (output of **Extract Genomic DNA** {% icon tool %})
>    - *"Sequence alphabet"*: ``
>    - *"Options Configuration"*: `Advanced`
>        - *"Limit of sequences to pass to MEME"*: `100`
>        - *"Should subsampling be random?"*: `Yes`
>            - *"Seed for the randomized selection of sequences"*: `123`
>        - *"maximum size of a sequence before it is cut down to a centered section"*: `0`
>        - *"Search given strand only"*: `Yes`
>        - *"What is the expected motif site distribution?"*: `Zero or one occurances per sequence`
>        - *"Minimum motif width"*: `5`
>        - *"Maximum motif width"*: `20`
>        - *"Maximum number of motifs to find"*: `20`
>        - *"Stop DREME searching after finding this many motifs"*: `5`
>    - *"I certify that I am not using this tool for commercial purposes."*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Create a BedGraph of genome coverage**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Create a BedGraph of genome coverage** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"The BAM or BED file from which coverage should be computed"*: `output` (output of **SortBED** {% icon tool %})
>    - *"Report regions with zero coverage"*: `Yes`
>    - *"Treat split/spliced BAM or BED12 entries as distinct BED intervals when computing coverage."*: `Yes`
>    - *"Calculate coverage based on"*: ``
>    - *"Scale the coverage by a constant factor"*: `""`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Wig/BedGraph-to-bigWig**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Wig/BedGraph-to-bigWig** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Convert"*: `output` (output of **Create a BedGraph of genome coverage** {% icon tool %})
>    - *"Converter settings to use"*: `Default`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}


## Re-arrange

To create the template, each step of the workflow had its own subsection. 

***TODO***: *Re-arrange the generated subsections into sections or other subsections. 
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.