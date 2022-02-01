---
layout: tutorial_hands_on

title: Metabarcoding/eDNA through Obitools
zenodo_link: https://zenodo.org/record/5932108/files/wolf_tutorial.zip?download=1
questions:
- how to analyze DNA metabarcoding / eDNA data produced on Illumina sequencers using the OBITools?
objectives:
- Deal with paired-end data to create consensus sequences
- Clean, filter and anlayse data to obtain strong results
time_estimation: 1H
key_points:
- From raw reads you can process, clean and filter data to obtain a list of species from environmental DNA (eDNA) samples.
contributors:
- colineroyaux
- onorvez
- obitools_team
- yvanlebras

---


# Introduction
{:.no_toc}

The data used in this tutorial correspond to the [official OBITools tutorial](https://pythonhosted.org/OBITools/wolves.html) and show how to analyse four wolf scats, using the protocol published in Shehzad et al. (2012) for assessing carnivore diet. After extracting DNA from the faeces, the DNA amplifications were carried out using the primers TTAGATACCCCACTATGC and TAGAACAGGCTCCTCTAG amplifiying the 12S-V5 region (Riaz et al. 2011), together with a wolf blocking oligonucleotide.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# OBITools overview
It is always a good idea to have a look at the intermediate results or to evaluate the best parameter for each step. Some commands are designed for that purpose, for example you can use :

- obicount to count the number of sequence records in a file
- obihead and obitail to view the first or last sequence records of a file
- obistat to get some basic statistics (count, mean, standard deviation) on the attributes (key=value combinations) in the header of each sequence record (see The extended OBITools fasta format in the fasta format description)
- any Galaxy tools corresponding to classical unix command such as less, awk, sort, wc to check your files.

# Manage input data
The data needed to run the tutorial are the following:

- fastq files resulting of a GA IIx (Illumina) paired-end (2 x 108 bp) sequencing assay of DNA extracted and amplified from four wolf faeces:
    - wolf_F.fastq
    - wolf_R.fastq
- the file describing the primers and tags used for all samples sequenced:
    - wolf_diet_ngsfilter.txt The tags correspond to short and specific sequences added on the 5’ end of each primer to distinguish the different samples
- the file containing the reference database in a fasta format:
    -db_v05_r117.fasta This reference database has been extracted from the release 117 of EMBL using ecoPCR

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the zip archive containing input files from [Zenodo](https://zenodo.org/record/5932108/files/wolf_tutorial.zip?download=1) 
>
>    ```
>    https://zenodo.org/record/5932108/files/wolf_tutorial.zip?download=1
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the dataset, here a zip archive, if needed
> 4. Check that the datatype is `zip`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>

## **Unzip** the downloaded archive

> ### {% icon hands_on %} Hands-on: Unzip the downladed .zip archive and prepare unzipped files to be used by OBITools Galaxy tools
>
> 1. {% tool [Unzip](toolshed.g2.bx.psu.edu/repos/imgteam/unzip/unzip/0.2) %} with the following parameters:
>    - *"Extract single file"*: `All files`
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
> 2. Add to each datafile a tag and/or modify names (*optional*)
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
> 3. Unhide all dataset from the resulting data collection so you can use these files independently.
>
> 4. Modify datatype from txt to tabular for the `wolf_diet_ngsfilter` dataset
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




# Use OBITools

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


## Sub-step with **illuminapairedend**

> ### {% icon hands_on %} Hands-on: Recover consensus sequences from overlapping forward and reverse reads.
> 
> When using the result of a paired-end sequencing assay with supposedly overlapping forward and reverse reads, the first step is to recover the assembled sequence.
>
> The forward and reverse reads of the same fragment are at the same line position in the two fastq files obtained after sequencing. Based on these two files, the assembly of the forward and reverse reads is done with the illuminapairedend utility that aligns the two reads and returns the reconstructed sequence.
>
> 1. {% tool [illuminapairedend](toolshed.g2.bx.psu.edu/repos/iuc/obi_grep/obi_illuminapairedend/1.2.13) %} with the following parameters:
>    - *"Read from file"*: `wolf_F` for the 3p file
>    - *"Read from file"*: `wolfRF` for the 5p file
>    - *"minimum score for keeping aligment"*: `40.0`
>
>    > ### {% icon comment %} Comment
>    >
>    > Sequence records corresponding to the same read pair must be in the same order in the two files !
>    > 
>    > If the alignment score is below the defined score, here 40, the forward and reverse reads are not aligned but concatenated, and the value of the mode attribute in the sequence header is set to joined instead of alignment
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


## Sub-step with **obigrep**

> ### {% icon hands_on %} Hands-on: Remove unaligned sequence records
> We here use the value of the mode attribute in the sequence header to discard sequences not "joined" (see explanation about this mode on the previous step)
>
> 1. {% tool [obigrep](toolshed.g2.bx.psu.edu/repos/iuc/obi_grep/obi_grep/1.2.13) %} with the following parameters:
>    - *"Choose the sequence record selection option"*: `predicat`
>        - *"Python boolean expression to be evaluated for each sequence record."*: `mode!=joined`
>
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}


> ### {% icon question %} Questions
>
> 1. How do you verify the operation is successfull?
> 2. How many sequences are kept?
>
> > ### {% icon solution %} Solution
> >
> > 1. you can search in the input file content the presence of `mode=joined` and same on the output file (just clicking the eye to visualize the content of each file and typing CTRL+C for example to search `mode=joined` in the file, or using a regex Galaxy tool for example). You can also at least look at the size of the output file, if smaller than input file, this is a first good indication.
> > 2. You can use a Galaxy tool like `Line/Word/Character count of a dataset` to count the number of lines of each dataset (input and output of obigrep) and divided by 4 (as in a FastQ file, each sequence is represented by a block of 4 lines). 45 276 sequences for input file. 44 717 for output file. Thus  559 sequences.
> >
> {: .solution}
>
{: .question}

## Sub-step with **NGSfilter**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [NGSfilter](toolshed.g2.bx.psu.edu/repos/iuc/obi_ngsfilter/obi_ngsfilter/1.2.13) %} with the following parameters:
>    - *"Output data type"*: `fastq`
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

## Sub-step with **obiuniq**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [obiuniq](toolshed.g2.bx.psu.edu/repos/iuc/obi_uniq/obi_uniq/1.2.13) %} with the following parameters:
>    - *"Attribute to merge"*: `sample`
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
> 1. {% tool [FastQC](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.73+galaxy0) %} with the following parameters:
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

## Sub-step with **obiannotate**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [obiannotate](toolshed.g2.bx.psu.edu/repos/iuc/obi_annotate/obi_annotate/1.2.13) %} with the following parameters:
>    - In *"Keep only attribute with key"*:
>        - *"key"*: `count`
>        - *"if you want to specify a second key"*: `merged_sample`
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

## Sub-step with **obistat**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [obistat](toolshed.g2.bx.psu.edu/repos/iuc/obi_stat/obi_stat/1.2.13) %} with the following parameters:
>    - In *"Category attribute"*:
>        - {% icon param-repeat %} *"Insert Category attribute"*
>            - *"How would you specify the category attribute key?"*: `simply by a key of an attribute`
>                - *"Attribute used to categorize the sequence records"*: `count`
>    - *"Use a specific option"*: `no`
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

## Sub-step with **obigrep**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [obigrep](toolshed.g2.bx.psu.edu/repos/iuc/obi_grep/obi_grep/1.2.13) %} with the following parameters:
>    - *"Choose the sequence record selection option"*: `predicat`
>        - *"Python boolean expression to be evaluated for each sequence record."*: `count>=10`
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

## Sub-step with **obigrep**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [obigrep](toolshed.g2.bx.psu.edu/repos/iuc/obi_grep/obi_grep/1.2.13) %} with the following parameters:
>    - *"Choose the sequence record selection option"*: `lmin`
>        - *"lmin"*: `80`
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

## Sub-step with **obiclean**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [obiclean](toolshed.g2.bx.psu.edu/repos/iuc/obi_clean/obi_clean/1.2.13) %} with the following parameters:
>    - *"Threshold ratio between counts (rare/abundant counts) of two sequence records so that the less abundant one is a variant of the more abundant (default: 1, i.e. all less abundant sequences are variants)"*: `0.05`
>    - *"Do you want to select only sequences with the head status in a least one sample?"*: `Yes`
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

## Sub-step with **NCBI BLAST+ blastn**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [NCBI BLAST+ blastn](toolshed.g2.bx.psu.edu/repos/devteam/ncbi_blast_plus/ncbi_blastn_wrapper/2.10.1+galaxy0) %} with the following parameters:
>    - *"Subject database/sequences"*: `FASTA file from your history (see warning note below)`
>    - *"Set expectation value cutoff"*: `0.0001`
>    - *"Output format"*: `Tabular (extended 25 columns)`
>    - *"Advanced Options"*: `Hide Advanced Options`
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
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - *"With following condition"*: `c3>99.99`
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
