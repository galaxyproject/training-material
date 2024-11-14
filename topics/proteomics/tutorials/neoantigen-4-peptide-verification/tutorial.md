---
layout: tutorial_hands_on

title: Neoantigen_PepQuery2_Verification
zenodo_link: ''
questions:
- How can neoantigens be verified using bioinformatics tools?
- What is the role of mass spectrometry and peptide sequence databases in neoantigen discovery?
objectives:
- Understand the workflow for neoantigen validation.
- Apply bioinformatics tools to validate peptides and proteins.
- Interpret the results from various analytical steps.
time_estimation: 3H
key_points:
- Understand the workflow for neoantigen validation.
- Apply bioinformatics tools to validate peptides and proteins.
- Interpret the results from various analytical steps.
contributions:
  authorship:
    - subinamehta
    - katherine-d21
    - jj-umn
  editing:
    - pratikdjagtap
    - timothygriffin
requirements:
  -
    type: "internal"
    topic_name: proteomics
subtopic: neoantigen
follow_up_training:

    -
        type: "internal"
        topic_name: proteomics
        tutorials:
            - neoantigen-non-normal-database
tags: [label-free]
redirect_from:
- proteomics/tutorials/neoantigen-4-peptide-verification/tutorial

---


# Introduction

<!-- This is a comment. -->

General introduction about the topic and then an introduction of the
tutorial (the questions and the objectives). It is nice also to have a
scheme to sum up the pipeline used during the tutorial. The idea is to
give to trainees insight into the content of the tutorial and the (theoretical
and technical) key concepts they will learn.

You may want to cite some publications; this can be done by adding citations to the
bibliography file (`tutorial.bib` file next to your `tutorial.md` file). These citations
must be in bibtex format. If you have the DOI for the paper you wish to cite, you can
get the corresponding bibtex entry using [doi2bib.org](https://doi2bib.org).

With the example you will find in the `tutorial.bib` file, you can add a citation to
this article here in your tutorial like this:
{% raw %} `{% cite Batut2018 %}`{% endraw %}.
This will be rendered like this: {% cite Batut2018 %}, and links to a
[bibliography section](#bibliography) which will automatically be created at the end of the
tutorial.


**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Title for your first section

Give some background about what the trainees will be doing in the section.
Remember that many people reading your materials will likely be novices,
so make sure to explain all the relevant concepts.

## Title for a subsection
Section and subsection titles will be displayed in the tutorial index on the left side of
the page, so try to make them informative and concise!

# Hands-on Sections
Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

## Get data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

# Title of the section usually corresponding to a big step in the analysis

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> <details-title> More details about the theory </details-title>
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


## Sub-step with **msconvert**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [msconvert](toolshed.g2.bx.psu.edu/repos/galaxyp/msconvert/msconvert/3.0.20287.2) %} with the following parameters:
>    - {% icon param-file %} *"Input unrefined MS data"*: `output` (Input dataset)
>    - *"Do you agree to the vendor licenses?"*: `Yes`
>    - *"Output Type"*: `mgf`
>    - In *"Data Processing Filters"*:
>        - *"Apply peak picking?"*: `Yes`
>        - *"(Re-)calculate charge states?"*: `no`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **PepQuery2**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [PepQuery2](toolshed.g2.bx.psu.edu/repos/galaxyp/pepquery2/pepquery2/2.0.2+galaxy2) %} with the following parameters:
>    - *"Validation Task Type"*: `novel peptide/protein validation`
>    - In *"Input Data"*:
>        - *"Input Type"*: `peptide`
>            - *"Peptides?"*: `Peptide list from your history`
>                - {% icon param-file %} *"Peptide Sequences (.txt)"*: `output` (Input dataset)
>        - *"Protein Reference Database from"*: `history`
>            - {% icon param-file %} *"Protein Reference Database File"*: `output` (Input dataset)
>        - *"MS/MS dataset to search"*: ` Spectrum Datasets from history`
>            - {% icon param-file %} *"Spectrum File"*: `output` (output of **msconvert** {% icon tool %})
>        - *"Report Spectrum Scan as"*: `spectrum title in MGF`
>    - In *"Modifications"*:
>        - *"Fixed modification(s)"*: ``
>        - *"Variable modification(s)"*: ``
>    - In *"Digestion"*:
>        - *"Enzyme"*: `Non enzyme`
>    - In *"Mass spectrometer"*:
>        - In *"Tolerance"*:
>            - *"Precursor Unit"*: `ppm`
>        - In *"PSM"*:
>            - *"Fragmentation Method"*: `CID/HCD`
>            - *"Minimum Charge"*: `2`
>            - *"Maximum Charge"*: `3`
>            - *"Minimum length of peptide"*: `8`
>            - *"Maximum length of peptide"*: `9`
>    - *"Select outputs"*: ``
>    - *"Use fast mode for searching"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Query Tabular**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.3.2) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `psm_rank_txt` (output of **PepQuery2** {% icon tool %})
>    - *"SQL Query to generate tabular output"*: `SELECT c1,c4
FROM t1
WHERE (c20 = 'Yes')
`
>    - *"include query result column headers"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Tabular-to-FASTA**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Tabular-to-FASTA](toolshed.g2.bx.psu.edu/repos/devteam/tabular_to_fasta/tab2fasta/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Tab-delimited file"*: `output` (output of **Query Tabular** {% icon tool %})
>    - *"Title column(s)"*: `c['2']`
>    - *"Sequence column"*: `c1`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **NCBI BLAST+ blastp**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NCBI BLAST+ blastp](toolshed.g2.bx.psu.edu/repos/devteam/ncbi_blast_plus/ncbi_blastp_wrapper/2.14.1+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Protein query sequence(s)"*: `output` (output of **Tabular-to-FASTA** {% icon tool %})
>    - *"Subject database/sequences"*: `Locally installed BLAST database`
>        - *"Protein BLAST database"*: ``
>    - *"Type of BLAST"*: `blastp-short - BLASTP optimized for queries shorter than 30 residues`
>    - *"Set expectation value cutoff"*: `200000.0`
>    - *"Output format"*: `Tabular (extended 25 columns)`
>    - *"Advanced Options"*: `Show Advanced Options`
>        - *"Scoring matrix and gap costs"*: `PAM30`
>            - *"Gap Costs"*: `Existence: 9  Extension: 1`
>        - *"Maximum number of HSPs (alignments) to keep for any single query-subject pair"*: `1`
>        - *"Word size for wordfinder algorithm"*: `2`
>        - *"Multiple hits window size: use 0 to specify 1-hit algorithm, leave blank for default"*: `15`
>        - *"Minimum score to add a word to the BLAST lookup table."*: `16`
>        - *"Composition-based statistics"*: `0: No composition-based statistics`
>        - *"Restrict search of database to a given set of ID's"*: `Taxonomy identifiers (TaxId's)`
>            - {% icon param-file %} *"Restrict search of database to list of TaxId's"*: `output` (Input dataset)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Query Tabular**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.3.2) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `output1` (output of **NCBI BLAST+ blastp** {% icon tool %})
>            - In *"Table Options"*:
>                - *"Specify Name for Table"*: `blast`
>                - *"Specify Column Names (comma-separated list)"*: `qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,sallseqid,score,nident,positive,gaps,ppos,qframe,sframe,qseq,sseq,qlen,slen,salltitles`
>                - In *"Table Index"*:
>                    - {% icon param-repeat %} *"Insert Table Index"*
>                        - *"Index on Columns"*: `qseqid`
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `output` (output of **Query Tabular** {% icon tool %})
>            - In *"Table Options"*:
>                - *"Specify Name for Table"*: `pep`
>                - *"Specify Column Names (comma-separated list)"*: `pep,seq`
>                - In *"Table Index"*:
>                    - {% icon param-repeat %} *"Insert Table Index"*
>                        - *"Index on Columns"*: `pep`
>    - *"SQL Query to generate tabular output"*: `SELECT DISTINCT pep.*
FROM pep 
JOIN blast ON pep.pep = blast.qseq
WHERE pep.pep NOT IN (
    SELECT qseq 
    FROM blast 
    WHERE pident = 100
)
AND (blast.pident < 100 
     OR blast.gapopen >= 1 
     OR blast.length < blast.qlen)
ORDER BY pep.pep`
>    - *"include query result column headers"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
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

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.

# Disclaimer 

Please note that all the software tools used in this workflow are subject to version updates and changes. As a result, the parameters, functionalities, and outcomes may differ with each new version. Additionally, if the protein sequences are downloaded at different times, the number of sequences may also vary due to updates in the reference databases or tool modifications. We recommend the users to verify the specific versions of software tools used to ensure the reproducibility and accuracy of results.
