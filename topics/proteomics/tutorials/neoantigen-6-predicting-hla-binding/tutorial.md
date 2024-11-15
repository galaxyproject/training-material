---
layout: tutorial_hands_on

title: "Neoantigen 6: Predicting HLA Binding"
zenodo_link: ''
questions:
- How can we predict the neoantigens presented by tumor cells?
- How does the HLA genotype affect the immune response to cancer?
objectives:
The following tools will be used in this tutorial:
  - Predict potential neoantigens based on HLA binding affinity.
  - Understand the role of HLA genotyping in predicting personalized immune responses.
  - Use specific tools for processing sequence data to predict HLA-binding peptides.
time_estimation: 3H
key_points:
- Neoantigen prediction is a key step in personalized medicine.
- HLA binding affinity helps identify potential neoantigens.
- Tools like OptiType, seq2HLA, and query tabular tools are essential for the workflow.
contributors:
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
            - neoantigen-7-hla-binding-novel-peptides
tags: [label-free]
redirect_from:
- proteomics/tutorials/neoantigen-6-predicting-hla-binding/tutorial

---


# Introduction

Neoantigen prediction for HLA binding is a critical component of personalized cancer immunotherapy. Neoantigens, which are tumor-specific antigens resulting from mutations in cancer cells, can be recognized by the immune system, making them promising targets for tailored immunotherapies. Human leukocyte antigen (HLA) molecules play a key role in presenting these neoantigens on the surface of cells, where they can be detected by T-cells, triggering an immune response.

This tutorial focuses on predicting HLA binding affinities for potential neoantigens, an essential step in identifying effective targets for immunotherapy. The workflow includes the use of tools such as OptiType and seq2HLA to determine HLA genotypes, as well as data reformatting and querying methods to manage complex data outputs.

By completing this tutorial, learners will gain an understanding of:

- The bioinformatics approaches required to predict HLA binding,
- Practical skills for neoantigen identification, including managing and analyzing genomic data,
- How to apply specific bioinformatics tools to analyze HLA binding predictions effectively.

This hands-on tutorial is designed to offer practical experience for learners aiming to integrate bioinformatics techniques into their immunotherapy research or clinical practice, emphasizing the theoretical and technical concepts essential for neoantigen prediction.

![HLAtyping-overview]({% link topics/proteomics/images/neoantigen/Prediction_of_HLA_Binding.PNG %})

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Predicting HLA binding
This tutorial provides a step-by-step guide for predicting HLA binding of neoantigens, a crucial part of personalized immunotherapy research. Using OptiType and seq2HLA, we will perform HLA typing and analyze which neoantigens are likely to bind to a specific individual’s HLA molecules, potentially driving an immune response. This process is essential for identifying candidate peptides that could serve as effective targets in immunotherapy.

### 1. Get Data
The tutorial begins with the acquisition and upload of necessary sequencing files. These files typically contain raw paired-end reads or other sequencing outputs, which will be processed to identify HLA alleles. Proper organization and tagging of data are crucial for maintaining an efficient workflow.

### 2. HLA Typing with OptiType
OptiType is used first in this tutorial for HLA typing. This tool processes the uploaded sequencing data to identify HLA class I alleles, which are essential for the subsequent neoantigen binding predictions. OptiType’s output will include the predicted HLA alleles that are specific to the individual’s immune profile.

### 3. HLA Typing with seq2HLA
As a complementary approach, we will use seq2HLA for HLA typing to cross-validate results. By comparing predictions from both OptiType and seq2HLA, we can ensure greater accuracy in identifying relevant HLA alleles, allowing us to confidently proceed with neoantigen predictions.

### 4. Reformatting and Filtering HLA Alleles
After obtaining results from both HLA typing tools, the next step is to reformat and filter the data. This involves removing redundant entries and preparing the data in a format that’s optimized for binding prediction. Cleaning and filtering ensure that only the most relevant HLA alleles are used in the neoantigen analysis.

### 5. Querying and Validating HLA Typing Results
The final step is to use SQL-like queries to further refine and validate the HLA typing results. This allows us to filter and focus on the HLA alleles most relevant to our research question, ensuring the quality and relevance of the dataset for predicting neoantigen binding.

This structured workflow enables a streamlined approach for accurate HLA typing and binding prediction, setting the foundation for personalized immunotherapy and advancing the development of effective cancer treatments.

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


## HLA typing with **OptiType**

OptiType is a bioinformatics tool designed specifically for HLA class I typing using paired-end sequencing reads. In this workflow, OptiType identifies the HLA alleles that a person possesses, which are crucial for understanding which neoantigens may effectively bind to an individual’s immune receptors. This tool enhances the personalization of immunotherapy by providing highly accurate HLA allele predictions. In turn, this allows for more targeted analysis of potential neoantigen binding candidates.

> <hands-on-title> OptiType </hands-on-title>
>
> 1. {% tool [OptiType](toolshed.g2.bx.psu.edu/repos/iuc/optitype/optitype/1.3.5+galaxy0) %} with the following parameters:
>    - *"Single or Paired-end reads"*: `Paired`
>        - {% icon param-file %} *"Select first set of reads"*: `output` (Input dataset)
>        - {% icon param-file %} *"Select second set of reads"*: `output` (Input dataset)
>    - *"Enumerations"*: `3`
>
>
{: .hands_on}


> <question-title></question-title>
>
> 1. Why do we use paired-end reads in OptiType?
> 2. What does the "Enumerations" parameter mean in OptiType, and why set it to 3?
>
> > <solution-title></solution-title>
> >
> > 1. Paired-end reads are used because they provide more comprehensive information about the DNA sequence. This improves the accuracy of HLA typing, as OptiType can analyze data from both directions of a DNA segment, reducing ambiguity and increasing the chance of correctly identifying HLA alleles.
> > 2. The "Enumerations" parameter controls the number of best-matching HLA allele combinations that OptiType will report. Setting it to 3 provides multiple possible allele combinations for verification and ensures that we capture the most likely HLA types while still keeping the output manageable for downstream analysis.
> >
> {: .solution}
>
{: .question}

## HLA typing with **seq2HLA**

seq2HLA is a computational tool for identifying HLA types from RNA-Seq or DNA-Seq data. It is especially useful in immunogenomics research, as it allows for the prediction of an individual’s HLA class I and II alleles, which play a key role in immune response. In this workflow, seq2HLA provides essential data for identifying potential neoantigens by aligning sequencing reads to reference HLA alleles, allowing researchers to understand the specific immune profile of each sample.

> <hands-on-title> seq2HLA </hands-on-title>
>
> 1. {% tool [seq2HLA](toolshed.g2.bx.psu.edu/repos/iuc/seq2hla/seq2hla/2.3+galaxy0) %} with the following parameters:
>    - *"Name prefix for this analysis"*: `STS26TGen`
>    - *"Paired-end reads"*: `Paired`
>        - {% icon param-file %} *"Select first set of reads"*: `output` (Input dataset)
>        - {% icon param-file %} *"Select second set of reads"*: `output` (Input dataset)
>
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Why is it important to specify a "Name prefix for this analysis"?
> 2. What is the significance of using paired-end reads in seq2HLA?
>
> > <solution-title></solution-title>
> >
> > 1. The "Name prefix for this analysis" parameter allows users to easily identify the output files related to a specific sample or experiment. By assigning a unique prefix, you can better organize and retrieve results, especially in workflows with multiple samples.
> > 2. Paired-end reads enhance seq2HLA's accuracy by providing complementary information from both ends of the DNA or RNA fragment. This improves the resolution and reliability of the HLA typing results, as more context is available to match against the reference HLA alleles.
> >
> {: .solution}
>
{: .question}

## Reformatting and Filtering HLA Alleles 

This step involves reformatting and filtering the HLA alleles output from the OptiType tool to retain only the relevant and unique alleles. The process is essential for streamlining the data and ensuring that only the significant allele information is passed on for further analysis. Using a simple AWK program, we can filter out redundant data and focus on distinct HLA alleles, which are crucial for downstream neoantigen discovery and immunogenomics applications.

> <hands-on-title> Text reformatting </hands-on-title>
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
>
{: .hands_on}

> <question-title></question-title>
>
> 1. What is the purpose of using the AWK program in this step?
> 2. Why is it important to filter out redundant HLA alleles?
>
> > <solution-title></solution-title>
> >
> > 1. The AWK program here helps identify and filter unique HLA alleles from the output of OptiType. This step removes duplicates and ensures that only relevant alleles are retained for further analysis, improving data quality and accuracy.
> > 2. Filtering out redundant alleles simplifies the dataset and enhances computational efficiency, making downstream analyses quicker and more focused on significant variants. This is essential for accurately interpreting immune responses in personalized medicine and neoantigen discovery.
> >
> {: .solution}
>
{: .question}

## Querying and Validating HLA Typing Results

In this step, we use the Query Tabular tool to validate and organize the HLA typing results from OptiType and seq2HLA outputs. By querying and filtering the data, we create a consolidated list of HLA alleles, formatted consistently for further analysis. This process ensures that we have a clean dataset of HLA alleles that will improve the accuracy of downstream analyses, such as neoantigen prediction or immune response studies.

> <hands-on-title> Query tabular </hands-on-title>
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
>
{: .hands_on}


> <question-title></question-title>
>
> 1. Why is it necessary to use SQL queries on HLA typing results?
> 2. How do we know the query results are accurate for the next steps?
>
> > <solution-title></solution-title>
> >
> > 1. SQL queries allow us to filter, organize, and validate the HLA typing results efficiently. This approach ensures consistent formatting and removes unnecessary or redundant data, making the dataset more manageable for further analysis.
> > 2. By carefully setting regex patterns and filtering options, we can ensure that only valid HLA allele formats are included. Reviewing the query results for consistency and completeness further ensures accuracy before moving to the next steps in the pipeline.
> >
> {: .solution}
>
{: .question}

# Conclusion

This tutorial covers the workflow for identifying and validating HLA alleles using OptiType and seq2HLA. Each step in the process, from reformatting and filtering to querying, is essential for generating accurate and consistent allele data. This validated data can then be used confidently in downstream immunological analyses, including neoantigen prediction and personalized medicine applications.


# Disclaimer 

Please note that all the software tools used in this workflow are subject to version updates and changes. As a result, the parameters, functionalities, and outcomes may differ with each new version. Additionally, if the protein sequences are downloaded at different times, the number of sequences may also vary due to updates in the reference databases or tool modifications. We recommend the users to verify the specific versions of software tools used to ensure the reproducibility and accuracy of results.
