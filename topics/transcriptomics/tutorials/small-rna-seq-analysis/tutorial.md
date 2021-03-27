---
layout: tutorial_hands_on

title: Small RNA-seq analysis
zenodo_link: https://zenodo.org/record/2650182#.XMB3UKbgpp8
questions:
- Which biological questions are addressed by the tutorial?
- Which bioinformatics techniques are important to know for this type of data?
objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentences describing what a learner should be able to do once they
  have completed the tutorial
- You can use Bloom's Taxonomy to write effective learning objectives
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- stephanierobin
- abretaud

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

# Analysis strategy and upload data

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
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/A12.fastq.gz
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/A16.fastq.gz
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/A18.fastq.gz
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/A19.fastq.gz
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/A22.fastq.gz
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/A24.fastq.gz
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/A27.fastq.gz
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/A29.fastq.gz
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/api-mir-mature.fa
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/api-mir-stem-loop.fa
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/genome.fa
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/K12.fastq.gz
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/K16.fastq.gz
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/K18.fastq.gz
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/K19.fastq.gz
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/K22.fastq.gz
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/K24.fastq.gz
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/K27.fastq.gz
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/K29.fastq.gz
>    https://zenodo.org/api/files/d2858738-eb3d-4953-8173-9ea539451452/utr.fa
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Rename each dataset according to the sample id (e.g. A12.fastq.gz)
>    {% include snippets/rename_dataset.md %}
> 4. Check that the datatype is fastq.gz for sequencing files and fasta for other files
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 5. Add to each dataset 2 tags corresponding to the 2 conditions  (A, K, T1, T2)
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

# Identification and quantification of miRNA

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


## Mapping with **MiRDeep2 Mapper**

The miRDeep2 Mapper module process and map small RNA sequencing data to the reference genome. First step is adapter removal, length filtering and collapsing of identical read sequences. Processed reads are mapped against the reference genome.
The input files are the reads files and the genome sequence file. If clip 3' adapter sequence is selected, adapter Sequence should be added. The output files are a file with the processed sequencing reads (fasta format) and a file with the reads mapped against the genome (arf format : proprietary file format generated and processed by miRDeep2).

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **MiRDeep2 Mapper** {% icon tool %} with the following parameters:
>    - *"Pool multiple read sets"*: `Yes`
>        - In *"Reads"*:
>            - {% icon param-repeat %} *"Insert Reads"*
>                - *"Sample name"*: `A12`
>            - {% icon param-repeat %} *"Insert Reads"*
>                - *"Sample name"*: `A16`
>            - {% icon param-repeat %} *"Insert Reads"*
>                - *"Sample name"*: `A18`
>            - {% icon param-repeat %} *"Insert Reads"*
>                - *"Sample name"*: `A19`
>            - {% icon param-repeat %} *"Insert Reads"*
>                - *"Sample name"*: `A22`
>            - {% icon param-repeat %} *"Insert Reads"*
>                - *"Sample name"*: `A24`
>            - {% icon param-repeat %} *"Insert Reads"*
>                - *"Sample name"*: `A27`
>            - {% icon param-repeat %} *"Insert Reads"*
>                - *"Sample name"*: `A29`
>            - {% icon param-repeat %} *"Insert Reads"*
>                - *"Sample name"*: `K12`
>            - {% icon param-repeat %} *"Insert Reads"*
>                - *"Sample name"*: `K16`
>            - {% icon param-repeat %} *"Insert Reads"*
>                - *"Sample name"*: `K18`
>            - {% icon param-repeat %} *"Insert Reads"*
>                - *"Sample name"*: `K19`
>            - {% icon param-repeat %} *"Insert Reads"*
>                - *"Sample name"*: `K22`
>            - {% icon param-repeat %} *"Insert Reads"*
>                - *"Sample name"*: `K24`
>            - {% icon param-repeat %} *"Insert Reads"*
>                - *"Sample name"*: `K27`
>            - {% icon param-repeat %} *"Insert Reads"*
>                - *"Sample name"*: `K29`
>    - *"Remove reads with non-standard nucleotides"*: `Yes`
>    - *"Convert RNA to DNA alphabet (to map against genome)"*: `Yes`
>    - *"Clip 3' Adapter Sequence"*: `Clip Sequence`
>        - *"Sequence to clip"*: `TGGAATTCTCGGGTGCCAAGG`
>    - *"Collapse reads and/or Map"*: `Collapse reads and Map`
>        - *"Will you select a reference genome from your history or use a built-in index?"*: `Use one from the history`
>        - *"Select the reference genome"*: `genome.fa`
>
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

## Identification of new miRNAs with **MiRDeep2**

The miRDeep2 module identifies known and novel miRNAs in small RNA-seq data.
The input files are the genome sequence file, the mapping file, and optionally, the sequence files of known miRNAs : mature miRNAs for the species analyzed, and/or for related species, precursor miRNAs and star miRNAs.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **MiRDeep2** {% icon tool %} with the following parameters:
>    - *"Collapsed deep sequencing reads"*: ``
>    - *"Detailed fasta output"*: `Yes`
>    - *"Detailed fasta output"*: `Yes`
>    - *"Detailed fasta output"*: `Yes`
>    - *"Detailed fasta output"*: `Yes`
>    - *"Detailed fasta output"*: `Yes`
>    - *"Detailed fasta output"*: `Yes`
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

## Concatenate known and novel mature miRNA with **Concatenate multiple datasets**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Concatenate multiple datasets** {% icon tool %} with the following parameters:
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

## Concatenate known and novel precursor miRNA with **Concatenate multiple datasets**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Concatenate multiple datasets** {% icon tool %} with the following parameters:
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


## Quantification of miRNA with **MiRDeep2 Quantifier**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **MiRDeep2 Quantifier** {% icon tool %} with the following parameters:
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

## Remove columns with **Remove_columns**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Remove_columns** {% icon tool %} with the following parameters:
>    - *"Remove column by"*: `column number`
>        - In *"Column number"*:
>            - {% icon param-repeat %} *"Insert Column number"*
>                - *"Column number"*: `1`
>            - {% icon param-repeat %} *"Insert Column number"*
>                - *"Column number"*: `2`
>            - {% icon param-repeat %} *"Insert Column number"*
>                - *"Column number"*: `4`
>            - {% icon param-repeat %} *"Insert Column number"*
>                - *"Column number"*: `21-36`
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

## Analysis of differential gene expression with **edgeR**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **edgeR** {% icon tool %} with the following parameters:
>    - *"Count Files or Matrix?"*: `Single Count Matrix`
>        - *"Input factor information from file?"*: `No`
>            - In *"Factor"*:
>                - {% icon param-repeat %} *"Insert Factor"*
>                    - *"Factor Name"*: `treatment`
>                    - *"Groups"*: `A,A,A,A,A,A,A,A,K,K,K,K,K,K,K,K`
>    - *"Use Gene Annotations?"*: `No`
>    - In *"Contrast"*:
>        - {% icon param-repeat %} *"Insert Contrast"*
>            - *"Contrast of Interest"*: `A-K`
>    - In *"Filter Low Counts"*:
>        - *"Filter lowly expressed genes?"*: `Yes`
>            - *"Filter on CPM or Count values?"*: `CPM`
>                - *"Minimum CPM"*: `10.0`
>                - *"Minimum Samples"*: `3`
>    - In *"Output Options"*:
>        - *"Output Normalised Counts Table?"*: `Yes`
>        - *"Output Rscript?"*: `Yes`
>        - *"Output RData file?"*: `Yes`
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
# Identification of miRNA targets with **miRanda**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **miRanda** {% icon tool %} with the following parameters:
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
