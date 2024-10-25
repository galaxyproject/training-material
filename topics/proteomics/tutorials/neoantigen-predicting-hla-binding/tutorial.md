---
layout: tutorial_hands_on

title: Neoantigen_Predicting_HLA_Binding
zenodo_link: ''
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
- contributor1
- contributor2

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


## Sub-step with **OptiType**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [OptiType](toolshed.g2.bx.psu.edu/repos/iuc/optitype/optitype/1.3.5+galaxy0) %} with the following parameters:
>    - *"Single or Paired-end reads"*: `Paired`
>        - {% icon param-file %} *"Select first set of reads"*: `output` (Input dataset)
>        - {% icon param-file %} *"Select second set of reads"*: `output` (Input dataset)
>    - *"Enumerations"*: `3`
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

## Sub-step with **seq2HLA**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [seq2HLA](toolshed.g2.bx.psu.edu/repos/iuc/seq2hla/seq2hla/2.3+galaxy0) %} with the following parameters:
>    - *"Name prefix for this analysis"*: `STS26TGen`
>    - *"Paired-end reads"*: `Paired`
>        - {% icon param-file %} *"Select first set of reads"*: `output` (Input dataset)
>        - {% icon param-file %} *"Select second set of reads"*: `output` (Input dataset)
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

## Sub-step with **Text reformatting**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Text reformatting](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_awk_tool/9.3+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `result` (output of **OptiType** {% icon tool %})
>    - *"AWK Program"*: `$1 ~ /[0-9]/{ 
    for (i = 2; i <=7; i++) { allele[$i]++}
}
END {
    for (i in allele) {
        print i
    }
}`
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
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `outfile` (output of **Text reformatting** {% icon tool %})
>            - In *"Filter Dataset Input"*:
>                - In *"Filter Tabular Input Lines"*:
>                    - {% icon param-repeat %} *"Insert Filter Tabular Input Lines"*
>                        - *"Filter By"*: `by regex expression matching`
>                            - *"regex pattern"*: `^(\w+[*]\d\d:\d\d\d?).*$`
>                            - *"action for regex match"*: `include line on pattern match`
>                    - {% icon param-repeat %} *"Insert Filter Tabular Input Lines"*
>                        - *"Filter By"*: `regex replace value in column`
>                            - *"enter column number to replace"*: `c1`
>                            - *"regex pattern"*: `^(\w+[*]\d\d:\d\d\d?).*$`
>                            - *"replacement expression"*: `HLA-\1`
>            - In *"Table Options"*:
>                - *"Specify Name for Table"*: `optitype`
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `c1_genotype4digits` (output of **seq2HLA** {% icon tool %})
>            - In *"Filter Dataset Input"*:
>                - In *"Filter Tabular Input Lines"*:
>                    - {% icon param-repeat %} *"Insert Filter Tabular Input Lines"*
>                        - *"Filter By"*: `skip leading lines`
>                            - *"Skip lines"*: `1`
>                    - {% icon param-repeat %} *"Insert Filter Tabular Input Lines"*
>                        - *"Filter By"*: `select columns`
>                            - *"enter column numbers to keep"*: `2,4`
>                    - {% icon param-repeat %} *"Insert Filter Tabular Input Lines"*
>                        - *"Filter By"*: `regex replace value in column`
>                            - *"enter column number to replace"*: `1`
>                            - *"regex pattern"*: `^(\w+[*]\d\d:\d\d\d?).*$`
>                            - *"replacement expression"*: `HLA-\1`
>                    - {% icon param-repeat %} *"Insert Filter Tabular Input Lines"*
>                        - *"Filter By"*: `regex replace value in column`
>                            - *"enter column number to replace"*: `2`
>                            - *"regex pattern"*: `^(\w+[*]\d\d:\d\d\d?).*$`
>                            - *"replacement expression"*: `HLA-\1`
>            - In *"Table Options"*:
>                - *"Specify Name for Table"*: `seq2hla`
>    - *"SQL Query to generate tabular output"*: `SELECT hla
FROM
(SELECT c1 AS hla FROM optitype
UNION
SELECT c1 AS hla FROM seq2hla WHERE c1 LIKE '%*%:%'
UNION 
SELECT c2 AS hla FROM seq2hla WHERE c2 LIKE '%*%:%') 
ORDER BY hla`
>    - *"include query result column headers"*: `No`
>    - In *"Additional Queries"*:
>        - In *"SQL Query"*:
>            - {% icon param-repeat %} *"Insert SQL Query"*
>                - *"SQL Query to generate tabular output"*: `SELECT c1 FROM optitype
ORDER BY c1`
>                - *"include query result column headers"*: `No`
>            - {% icon param-repeat %} *"Insert SQL Query"*
>                - *"SQL Query to generate tabular output"*: `SELECT c1 FROM seq2hla WHERE c1 LIKE '%*%:%' 
UNION 
SELECT c2 FROM seq2hla WHERE c2 LIKE '%*%:%'`
>                - *"include query result column headers"*: `No`
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