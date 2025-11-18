---
layout: tutorial_hands_on

title: Dataset construction for bacterial comparative genomics
zenodo_link: 'https://zenodo.org/records/1'
draft: true
questions:
- to do
- to do
objectives:
- to do
- to do
time_estimation: 1H
key_points:
- to do
- to do
contributions:
  authorship:
  - hchiapello
  - bebatut
  - vloux
edam_ontology:
- topic_0622 # Genomics
- topic_3301 # Microbiology
level: Introductory
tags:
- microgalaxy
---

Introduction


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get genomes from GTDB

> <hands-on-title> Task description </hands-on-title>
>
> 1. 
>
{: .hands_on}



> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this analysis
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Rename the history
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
> 3. {% tool [Import](upload1) %} the contig file from [Zenodo]({{ page.zenodo_link }}) or from Galaxy shared data libraries:
>
>    ```
>    {{ page.zenodo_link }}/files/DRR187559_contigs.fasta
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
>
{: .hands_on}



# Extract IDs for NCBI

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1`
>    - *"Delimited by"*: `Comma`
>    - {% icon param-file %} *"From"*: `output` (Input dataset)
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


> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c6`
>    - {% icon param-file %} *"From"*: `output` (Input dataset)
>
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


> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Concatenate datasets](cat1) %} with the following parameters:
>    - {% icon param-file %} *"Concatenate Dataset"*: `out_file1` (output of **Cut** {% icon tool %})
>    - In *"Dataset"*:
>        - {% icon param-repeat %} *"Insert Dataset"*
>            - {% icon param-file %} *"Select"*: `out_file1` (output of **Cut** {% icon tool %})
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

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Sort](sort1) %} with the following parameters:
>    - {% icon param-file %} *"Sort Dataset"*: `out_file1` (output of **Concatenate datasets** {% icon tool %})
>    - *"on column"*: `c1`
> 1. {% tool [Unique](toolshed.g2.bx.psu.edu/repos/bgruening/unique/bg_uniq/0.3) %} with the following parameters:
>    - {% icon param-file %} *"from query"*: `out_file1` (output of **Sort** {% icon tool %})
>    - *"Advanced Options"*: `Hide Advanced Options`
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

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Remove beginning](Remove beginning1) %} with the following parameters:
>    - *"Remove first"*: `2`
>    - {% icon param-file %} *"from"*: `outfile` (output of **Unique** {% icon tool %})
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

# Download genomes from NCBI

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NCBI Datasets Genomes](toolshed.g2.bx.psu.edu/repos/iuc/ncbi_datasets/datasets_download_genome/16.42.0+galaxy0) %} with the following parameters:
>    - In *"Query"*:
>        - *"Choose how to find genomes to download"*: `By NCBI assembly or BioProject accession`
>            - *"Enter accession or read from file ?"*: `Read a list of NCBI Assembly accessions from a dataset`
>                - {% icon param-file %} *"Select dataset with list of NCBI Assembly accessions"*: `out_file1` (output of **Remove beginning** {% icon tool %})
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


# Conclusion

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.