---
layout: tutorial_hands_on

title: Neoantigen_PepPointer_Annotation
zenodo_link: ''
questions:
- How can neoantigens be identified in cancer genomes?
- What role do neoantigens play in personalized immunotherapy?
- How can mutations in cancer cells be used to predict immune system responses?
objectives:
- Identify potential neoantigens from sequencing data.
- Annotate somatic mutations and predict peptide sequences.
- Predict MHC binding affinities for neoantigens.
- Interpret data using bioinformatics tools for cancer immunotherapy applications.
time_estimation: 3H
key_points:
- Neoantigen discovery involves the identification of somatic mutations that generate novel peptide sequences recognized by the immune system.
- Accurate prediction of MHC binding affinity is critical for evaluating the immunogenic potential of neoantigens.
- Immunoinformatics tools and pipelines are integral for analyzing neoantigens from genomic data and visualizing results.
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
- proteomics/tutorials/neoantigen-5-variant-annotation/tutorial

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


## Sub-step with **Query Tabular**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.3.2) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `output` (Input dataset)
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `output` (Input dataset)
>    - *"SQL Query to generate tabular output"*: `SELECT t1.c1,t2.c13,t2.c5,t2.c6
FROM t1
INNER JOIN t2
ON t1.c1 = t2.c1`
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

## Sub-step with **Convert**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Convert](Convert characters1) %} with the following parameters:
>    - *"Convert all"*: `Pipes`
>    - {% icon param-file %} *"in Dataset"*: `output` (output of **Query Tabular** {% icon tool %})
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

## Sub-step with **Column Regex Find And Replace**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Column Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regexColumn1/1.0.3) %} with the following parameters:
>    - {% icon param-file %} *"Select cells from"*: `out_file1` (output of **Convert** {% icon tool %})
>    - *"using column"*: `c3`
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `u_`
>            - *"Replacement"*: `u:`
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

## Sub-step with **Convert**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Convert](Convert characters1) %} with the following parameters:
>    - *"Convert all"*: `Colons`
>    - {% icon param-file %} *"in Dataset"*: `out_file1` (output of **Column Regex Find And Replace** {% icon tool %})
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
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `out_file1` (output of **Convert** {% icon tool %})
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `output` (Input dataset)
>    - *"SQL Query to generate tabular output"*: `SELECT t1.c1,t1.c8,t1.c9,t1.c11,t1.c12,t2.c5,t2.c6,t1.c7
FROM t1
INNER JOIN t2
ON t1.c3 = t2.c4`
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

## Sub-step with **Query Tabular**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.3.2) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `output` (output of **Query Tabular** {% icon tool %})
>    - *"SQL Query to generate tabular output"*: `SELECT t1.*, 
       t1.c4 * 3 AS c2_multiplied, 
       t1.c5 * 3 AS c3_multiplied
FROM t1`
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

## Sub-step with **Query Tabular**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.3.2) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `output` (output of **Query Tabular** {% icon tool %})
>    - *"SQL Query to generate tabular output"*: `SELECT t1.*, 
       CASE 
           WHEN t1.c7 = '-' THEN t1.c3 - t1.c9
           WHEN t1.c7 = '+' THEN t1.c2 + t1.c9
       END AS start,
       CASE 
           WHEN t1.c7 = '-' THEN t1.c3 - t1.c10
           WHEN t1.c7 = '+' THEN t1.c2 + t1.c10
       END AS stop
FROM t1`
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

## Sub-step with **Query Tabular**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.3.2) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `output` (output of **Query Tabular** {% icon tool %})
>    - *"SQL Query to generate tabular output"*: `SELECT 
    c8 AS `chromosome`,      
   c11  AS `chromStart`, 
   c12 AS `chromEnd`,         
   c1 AS `name`,                  
   c6 AS `score`,                    
    c7 AS `strand`           
FROM 
   t1`
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

## Sub-step with **PepPointer**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [PepPointer](toolshed.g2.bx.psu.edu/repos/galaxyp/pep_pointer/pep_pointer/0.1.3+galaxy1) %} with the following parameters:
>    - *"Choose the source of the GTF file"*: `From history`
>        - {% icon param-file %} *"GTF file with the genome of interest"*: `output` (Input dataset)
>    - {% icon param-file %} *"BED file with chromosomal coordinates of peptide"*: `output` (output of **Query Tabular** {% icon tool %})
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
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `classified` (output of **PepPointer** {% icon tool %})
>    - *"SQL Query to generate tabular output"*: `SELECT 
c4 AS Peptide,
c1 AS Chromosome,
c2 AS Start,
c3 AS End,
c6 AS Strand,
c7 AS Annotation,
c1||':'||c2||'-'||c3 AS IGV_Genome_Coordinate,
'https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position='||c1||'%3A'||c2||'-'||c3 AS UCSC_Genome_Browser
FROM  t1`
>    - *"include query result column headers"*: `Yes`
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

In this tutorial, we have covered the process of annotating neoantigen peptides using the PepPointer tool in a bioinformatics pipeline. The main focus was to guide you through the various steps involved in the analysis, including querying tabular datasets, converting data formats, and using the PepPointer tool to classify peptides based on genomic information.

By the end of this tutorial, you should now be able to:
- Understand the steps involved in preparing and processing peptide data for neoantigen prediction.
- Use tools like Query Tabular, Convert, and PepPointer to manipulate and classify peptide sequences.
- Recognize the importance of annotating neoantigens for immunotherapy research and how bioinformatics tools can facilitate this process.

The key takeaway from this tutorial is the ability to manipulate complex genomic data and generate insights into potential therapeutic targets using bioinformatics techniques. You should now be better equipped to handle similar datasets and apply these techniques to other types of genomic and proteomic research.
