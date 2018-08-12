---
layout: tutorial_hands_on
topic_name: transcriptomics
tutorial_name: scrna_preprocessing
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
> 1. Import the following files from [Zenodo]() or from a data
>    library named `TODO` if available (ask your instructor)
>
>    ```
>    
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


## Sub-step with **UMI-tools extract**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **UMI-tools extract** {% icon tool %} with the following parameters:
>    - *"Library type"*: `Paired-end Dataset Collection`
>        - {% icon param-collection %} *"Reads in FASTQ format"*: `output` (Input dataset collection)
>        - *"Barcode on both reads?"*: `Barcode on first read only`
>    - *"Use Known Barcodes?"*: `Yes`
>        - {% icon param-file %} *"Barcode File"*: `output` (Input dataset)
>    - *"Barcode pattern for first read"*: `NNNNNNCCCCCC`
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
>    - *"Single-end or paired-end reads"*: `Single-end`
>        - {% icon param-file %} *"RNA-Seq FASTQ/FASTA file"*: `out2` (output of **UMI-tools extract** {% icon tool %})
>    - *"Custom or built-in reference genome"*: `Use a built-in index`
>        - *"Reference genome with or without an annotation"*: `use genome reference without builtin gene-model`
>            - {% icon param-file %} *"Gene model (gff3,gtf) file for splice junctions"*: `output` (Input dataset)
>    - *"Would you like to set output parameters (formatting and filtering)?"*: `No`
>    - *"Other parameters (seed, alignment, limits and chimeric alignment)"*: `Use Defaults`
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

## Sub-step with **featureCounts**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **featureCounts** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Alignment file"*: `mapped_reads` (output of **RNA STAR** {% icon tool %})
>    - *"Gene annotation file"*: `in your history`
>        - {% icon param-file %} *"Gene annotation file"*: `output` (Input dataset)
>    - In *"Options for paired-end reads"*:
>        - *"Count fragments instead of reads"*: `Disabled; all reads/mates will be counted individually`
>    - In *"Advanced options"*:
>        - *"Count multi-mapping reads/fragments"*: `Disabled; multi-mapping reads are excluded (default)`
>        - *"Exon-exon junctions"*: `Yes`
>        - *"Annotates the alignment file with 'XS:Z:'-tags to described per read or read-pair the corresponding assigned feature(s)."*: `Yes`
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

## Sub-step with **Filter**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Filter** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"BAM dataset(s) to filter"*: `output_bam` (output of **featureCounts** {% icon tool %})
>    - In *"Condition"*:
>        - Click on *"Insert Condition"*:
>        - In *"1: Condition"*:
>            - In *"Filter"*:
>                - Click on *"Insert Filter"*:
>                - In *"1: Filter"*:
>                    - *"Select BAM property to filter on"*: `alignmentFlag`
>                        - *"Filter on this alignment flag"*: `0`
>        - Click on *"Insert Condition"*:
>        - In *"2: Condition"*:
>            - In *"Filter"*:
>                - Click on *"Insert Filter"*:
>                - In *"1: Filter"*:
>                    - *"Select BAM property to filter on"*: `alignmentFlag`
>                        - *"Filter on this alignment flag"*: `16`
>        - Click on *"Insert Condition"*:
>        - In *"3: Condition"*:
>            - In *"Filter"*:
>                - Click on *"Insert Filter"*:
>                - In *"1: Filter"*:
>                    - *"Select BAM property to filter on"*: `tag`
>                        - *"Filter on a particular tag"*: `nM:<3`
>        - Click on *"Insert Condition"*:
>        - In *"4: Condition"*:
>            - In *"Filter"*:
>                - Click on *"Insert Filter"*:
>                - In *"1: Filter"*:
>                    - *"Select BAM property to filter on"*: `tag`
>                        - *"Filter on a particular tag"*: `NH:<2`
>    - *"Would you like to set rules?"*: `Yes`
>        - *"Enter rules here"*: `(1 | 2) & 3 & 4`
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

## Sub-step with **UMI-tools count**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **UMI-tools count** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sorted BAM file"*: `out_file1` (output of **Filter** {% icon tool %})
>    - *"Umi Extract Method"*: `Barcodes are contained at the end of the read seperated by a delimiter`
>    - *"Method to identify group of reads"*: `Unique`
>    - *"Prepend a label to all column headers"*: `Dataset Name`
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