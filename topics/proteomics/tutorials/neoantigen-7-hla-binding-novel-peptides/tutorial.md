---
layout: tutorial_hands_on

title: Neoantigen_IEDB_binding_PepQuery_novel_peptides
zenodo_link: ''
questions:
- What are neoantigens, and why are they significant in cancer immunotherapy?
- How can binding predictions and validation help in distinguishing strong and weak binders?
- What tools and techniques are commonly used for neoantigen identification and validation?
objectives:
- Understand the process of neoantigen identification and the role of peptide binding predictions.
- Learn how to use IEDB to predict the binding affinity of peptides to MHC molecules.
- Gain practical experience using PepQuery to validate novel peptides from proteomics data.
- Distinguish between strong and weak binders based on predicted binding affinity.
time_estimation: 3H
key_points:
- Neoantigens are crucial targets for cancer immunotherapy, with their binding strength to MHC molecules playing a critical role in immune response.
- IEDB allows for the prediction of peptide-MHC binding affinities, helping to identify candidate neoantigens.
- PepQuery enables validation of these peptides in proteomics data, adding confidence to neoantigen identification.
- Separating strong and weak binders allows for prioritization in therapeutic applications.
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
tags: [label-free]
redirect_from:
- proteomics/tutorials/neoantigen-7-hla-binding-novel-peptides/tutorial

---


# Introduction

Neoantigens are peptides derived from tumor-specific mutations, which are recognized by the immune system as foreign and can stimulate an immune response against cancer cells. Identifying these neoantigens is a crucial step in the development of personalized cancer immunotherapies, as they serve as targets for T-cell mediated immune responses. However, predicting which peptides from the tumor genome will bind effectively to major histocompatibility complex (MHC) molecules—key proteins that present antigens to immune cells—remains a significant challenge.

This tutorial outlines a comprehensive workflow for the identification, prediction, and validation of potential neoantigens. We begin by using the Immune Epitope Database (IEDB) to predict the binding affinity of peptide sequences to MHC molecules. IEDB provides powerful tools to model how peptides interact with different MHC alleles, helping to prioritize peptides that are most likely to be presented by the immune system. Next, we validate these peptides using PepQuery, a tool that allows for the comparison of predicted neoantigens with experimental proteomics data, providing an additional layer of confidence in their relevance. Finally, we categorize the peptides into strong and weak binders, based on their predicted affinity, which helps in identifying the most promising candidates for cancer immunotherapy.


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
>
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


## Sub-step with **IEDB**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [IEDB](toolshed.g2.bx.psu.edu/repos/iuc/iedb_api/iedb_api/2.15.2) %} with the following parameters:
>    - *"Prediction"*: `MHC-I Binding`
>        - *"prediction method"*: `netmhcpan_el`
>        - *"Alleles"*: `From history`
>            - {% icon param-file %} *"Alleles file"*: `output` (Input dataset)
>        - *"peptide lengths for prediction"*: ``
>    - *"Peptide sequences"*: `Fasta file`
>        - {% icon param-file %} *"Peptide Sequence Fasta"*: `output` (Input dataset)
>
>
>
{: .hands_on}



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

## Sub-step with **Filter**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `output` (output of **IEDB** {% icon tool %})
>    - *"With following condition"*: `c11>0.5 and c11<=2`
>    - *"Number of header lines to skip"*: `1`
>
>
>
{: .hands_on}



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

## Sub-step with **Filter**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `output` (output of **IEDB** {% icon tool %})
>    - *"With following condition"*: `c11<=0.5`
>    - *"Number of header lines to skip"*: `1`
>
>
>
{: .hands_on}



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

## Sub-step with **Table Compute**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Table Compute](toolshed.g2.bx.psu.edu/repos/iuc/table_compute/table_compute/1.2.4+galaxy0) %} with the following parameters:
>    - *"Input Single or Multiple Tables"*: `Single Table`
>        - {% icon param-file %} *"Table"*: `out_file1` (output of **Filter** {% icon tool %})
>        - *"Type of table operation"*: `Perform a full table operation`
>            - *"Operation"*: `Pivot`
>                - *"Index"*: `icore`
>                - *"Column"*: `allele`
>                - *"Values"*: `percentile_rank`
>                - *"Aggregator Function"*: `Maximum`
>
>
>
{: .hands_on}



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

## Sub-step with **Table Compute**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Table Compute](toolshed.g2.bx.psu.edu/repos/iuc/table_compute/table_compute/1.2.4+galaxy0) %} with the following parameters:
>    - *"Input Single or Multiple Tables"*: `Single Table`
>        - {% icon param-file %} *"Table"*: `out_file1` (output of **Filter** {% icon tool %})
>        - *"Type of table operation"*: `Perform a full table operation`
>            - *"Operation"*: `Pivot`
>                - *"Index"*: `icore`
>                - *"Column"*: `allele`
>                - *"Values"*: `percentile_rank`
>                - *"Aggregator Function"*: `Maximum`
>
>
>
{: .hands_on}



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

## Sub-step with **Remove beginning**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Remove beginning](Remove beginning1) %} with the following parameters:
>    - {% icon param-file %} *"from"*: `table` (output of **Table Compute** {% icon tool %})
>
>
>
{: .hands_on}



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

## Sub-step with **Remove beginning**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Remove beginning](Remove beginning1) %} with the following parameters:
>    - {% icon param-file %} *"from"*: `table` (output of **Table Compute** {% icon tool %})
>
>
>
{: .hands_on}


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

## Sub-step with **Cut**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1`
>    - {% icon param-file %} *"From"*: `out_file1` (output of **Remove beginning** {% icon tool %})
>
>
>
{: .hands_on}



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

## Sub-step with **Cut**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1`
>    - {% icon param-file %} *"From"*: `out_file1` (output of **Remove beginning** {% icon tool %})
>
>
{: .hands_on}


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
>                - {% icon param-file %} *"Peptide Sequences (.txt)"*: `out_file1` (output of **Cut** {% icon tool %})
>        - *"Protein Reference Database from"*: `download`
>            - *"Public protein sequence database"*: `gencode:human`
>        - *"MS/MS dataset to search"*: ` Spectrum Datasets from history`
>            - {% icon param-file %} *"Spectrum File"*: `output` (Input dataset)
>        - *"Report Spectrum Scan as"*: `spectrum title in MGF`
>
>
>
{: .hands_on}



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
>                - {% icon param-file %} *"Peptide Sequences (.txt)"*: `out_file1` (output of **Cut** {% icon tool %})
>        - *"Protein Reference Database from"*: `download`
>            - *"Public protein sequence database"*: `gencode:human`
>        - *"MS/MS dataset to search"*: ` Spectrum Datasets from history`
>            - {% icon param-file %} *"Spectrum File"*: `output` (Input dataset)
>        - *"Report Spectrum Scan as"*: `spectrum title in MGF`
>
>
>
{: .hands_on}



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

## Sub-step with **Filter**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `psm_rank_txt` (output of **PepQuery2** {% icon tool %})
>    - *"With following condition"*: `c20=='Yes'`
>
>
>
{: .hands_on}



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

## Sub-step with **Filter**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `psm_rank_txt` (output of **PepQuery2** {% icon tool %})
>    - *"With following condition"*: `c20=='Yes'`
>
>
>
{: .hands_on}



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
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `table` (output of **Table Compute** {% icon tool %})
>            - In *"Table Options"*:
>                - *"Use first line as column names"*: `Yes`
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `out_file1` (output of **Filter** {% icon tool %})
>    - *"SQL Query to generate tabular output"*: `SELECT *
FROM t1 
WHERE t1.icore IN (SELECT c1 FROM t2)
`
>    - *"include query result column headers"*: `Yes`
>    - In *"Additional Queries"*:
>        - In *"SQL Query"*:
>            - {% icon param-repeat %} *"Insert SQL Query"*
>                - *"SQL Query to generate tabular output"*: `SELECT icore
FROM t1 
WHERE t1.icore IN (SELECT c1 FROM t2)
ORDER BY icore`
>                - *"include query result column headers"*: `No`
>
>
{: .hands_on}


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
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `table` (output of **Table Compute** {% icon tool %})
>            - In *"Table Options"*:
>                - *"Use first line as column names"*: `Yes`
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `out_file1` (output of **Filter** {% icon tool %})
>    - *"SQL Query to generate tabular output"*: `SELECT *
FROM t1 
WHERE t1.icore IN (SELECT c1 FROM t2)
`
>    - *"include query result column headers"*: `Yes`
>    - In *"Additional Queries"*:
>        - In *"SQL Query"*:
>            - {% icon param-repeat %} *"Insert SQL Query"*
>                - *"SQL Query to generate tabular output"*: `SELECT icore
FROM t1 
WHERE t1.icore IN (SELECT c1 FROM t2)
ORDER BY icore`
>                - *"include query result column headers"*: `No`
>
>
{: .hands_on}


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

The ability to predict and validate neoantigens is a cornerstone of modern cancer immunotherapy, particularly in personalized immunotherapies such as cancer vaccines and adoptive T-cell therapies. The identification of high-affinity neoantigens, which can effectively stimulate the immune system, is key to developing therapeutic strategies that selectively target tumor cells while sparing healthy tissue.

This workflow is particularly relevant in the neoantigen discovery process, as it combines computational prediction with experimental validation, ensuring that only the most promising peptides are considered for therapeutic development. By integrating IEDB’s predictive capabilities with PepQuery’s validation approach, this workflow bridges the gap between in silico predictions and in vitro experimental data, allowing for a more robust and reliable neoantigen identification pipeline. Separating peptides into strong and weak binders based on their predicted MHC affinity further refines the selection process, ensuring that high-priority candidates are prioritized for downstream therapeutic application.

Given the increasing demand for personalized cancer treatments, this workflow represents a vital approach for accelerating the identification of clinically relevant neoantigens, thus advancing the field of cancer immunotherapy and personalized medicine.

# Disclaimer 

Please note that all the software tools used in this workflow are subject to version updates and changes. As a result, the parameters, functionalities, and outcomes may differ with each new version. Additionally, if the protein sequences are downloaded at different times, the number of sequences may also vary due to updates in the reference databases or tool modifications. We recommend the users to verify the specific versions of software tools used to ensure the reproducibility and accuracy of results.
