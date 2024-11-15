---
layout: tutorial_hands_on

title: Neoantigen 7:IEDB binding PepQuery Validated Neopeptides
zenodo_link: ''
questions:
- What are neoantigens, and why are they significant in cancer immunotherapy?
- How can binding predictions and validation help distinguish strong and weak binders?
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

![HLA-binding-peptides-overview]({% link topics/proteomics/images/neoantigen/Predict_MHC_Binding_for_Novel_Neoantigens.PNG %})


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Neoantigen HLA binding Novel peptides

This workflow outlines a structured approach to predicting and analyzing peptide-HLA (MHC-I) binding affinities, focusing on neoantigens for immunotherapy. Here's an overview of the workflow, broken into key steps:

### 1. Get Data
The first step involves gathering and uploading the necessary proteomics data files into the analysis environment. These files typically contain protein sequences or raw spectrum data that will be processed throughout the tutorial. Proper data organization and tagging are essential to ensure smooth workflow execution.

### 2. Predicting MHC Binding with IEDB
In this step, the goal is to predict how peptide sequences will bind to MHC-I molecules, an essential part of identifying potential neoantigens. The IEDB tool is launched with parameters such as prediction method (netmhcpan_el), selecting relevant alleles, and uploading peptide sequences in FASTA format.

### 3. Filtering Results Based on Affinity
Once the MHC-I binding affinities are predicted, it's time to filter the results based on specific criteria to narrow down the most promising peptides.The peptides are categorized into strong and weak binders: an affinity value of less than 0.5 is a strong binder, while a value between 0.5 and 2 is a weak binder.

### 4. Refining Results Using Table Operations
To further refine the results, the data is aggregated, focusing on the highest-ranked peptides based on various parameters like affinity and percent rank. This tool is used to aggregate data and perform operations such as selecting the highest-ranking peptides based on specified criteria (e.g., using a pivot operation).

### 5. Novel peptide verification with PepQuery2 
The PepQuery2 tool is used to validate the identified peptides by checking them against a reference protein database to determine whether they are novel or already known. In this step, peptides are compared to the gencode protein reference database, and their sequences are verified using MS/MS spectral data from previous steps. The tool helps ensure the authenticity of the identified peptides, which is critical for confirming their relevance in the biological context. By determining whether peptides are novel or known, this step increases confidence in the results, supporting the identification of unique biomarkers or therapeutic targets and contributing to the overall quality of the research.

### 6. Filtering confident PepQuery2 validated peptides and annotating their binding affinity to HLA.
In this step, the peptides that have been confidently validated by PepQuery2 are filtered to focus on those that exhibit strong evidence of binding affinity to Human Leukocyte Antigen (HLA) molecules. This process involves selecting peptides with high confidence in their validation and annotating them with their corresponding HLA binding affinity data. By filtering and annotating these peptides, the analysis narrows down the most promising candidates for further investigation, ensuring that only those with reliable peptide-HLA interactions are considered for potential clinical or research applications. This step is critical for identifying peptides with the potential to serve as targets for immunotherapy or vaccine development.

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



## Predicting MHC Binding with IEDB
In this step of the workflow, the IEDB (Immune Epitope Database) tool is used to predict peptide binding to MHC-I molecules, which is crucial for identifying potential neoantigens. This is accomplished by using the netmhcpan_el prediction method to analyze peptide sequences and their binding affinity to specific MHC alleles. The netmhcpan_el method is a state-of-the-art tool that predicts how likely a peptide sequence is to bind to a given MHC-I allele. The alleles to be tested against are selected from a history file, ensuring that the analysis is tailored to specific genetic variations. Peptides that are predicted to bind MHC-I molecules are provided in a FASTA format, which is the standard for representing protein sequences.
The IEDB tool is integral in neoantigen discovery because it allows researchers to identify which peptides (out of many potential candidates) have a high likelihood of binding to MHC-I molecules. These peptides are of interest as they could trigger an immune response against cancer cells or pathogens, making them candidates for immunotherapy.


> <hands-on-title> IEDB </hands-on-title>
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
> 1. Why do I need to specify the netmhcpan_el method for peptide-MHC binding prediction?
> 2. What should I do if my peptide sequence file is not in the correct FASTA format?
>
> > <solution-title></solution-title>
> >
> > 1. The netmhcpan_el method is a highly accurate and widely used algorithm for predicting peptide binding to MHC-I molecules. By specifying this method, you're ensuring that the prediction is based on a well-established and reliable tool, increasing the likelihood of accurate results for identifying potential neoantigens.
> > 2. If your peptide sequence file is not in the FASTA format, you'll need to reformat the file before uploading it to the IEDB tool. FASTA format is a plain text format used to represent nucleotide or peptide sequences, where each sequence begins with a > symbol followed by a description line and then the sequence itself. Many bioinformatics tools, including IEDB, require this format to interpret the sequences correctly.
> >
> {: .solution}
>
{: .question}

## Filtering Results Based on Affinity

In this step, the Filter tool is used to refine the results from the previous IEDB step by selecting peptides based on their binding affinity to MHC-I molecules. Peptides are categorized into separate groups based on their binding affinities, with those that bind with lower or higher affinity being filtered accordingly. We utilize the affinity values to identify strong and weak binders. Specifically, peptides with an affinity value greater than 0.5 and less than or equal to 2 are classified as weak binders, while those with an affinity value less than 0.5 are classified as strong binders. This distinction is important, as peptides with very low or extremely high affinities may be as relevant for further analysis.

> <hands-on-title> Filter weak binders </hands-on-title>
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `output` (output of **IEDB** {% icon tool %})
>    - *"With following condition"*: `c11>0.5 and c11<=2`
>    - *"Number of header lines to skip"*: `1`
>
>
>
{: .hands_on}


> <hands-on-title> Filter strong binders </hands-on-title>
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
> 1. Why do we need to filter strong and weak binder
>
> > <solution-title></solution-title>
> >
> > 1. Filtering weak and strong binders is crucial for optimizing the identification of peptides that bind to MHC-I molecules with an optimal affinity, which is key for selecting potential neoantigens. Weak binders (with affinities greater than 0.5 but less than or equal to 2) are less likely to trigger a robust immune response, while strong binders (with affinities less than or equal to 0.5) might bind too tightly, potentially causing inefficient or excessive immune reactions. By filtering both extremes, the analysis focuses on peptides with moderate affinity that are more likely to stimulate an effective immune response without causing adverse effects. This targeted filtering not only enhances the relevance of the peptides for therapeutic applications but also improves the quality of predictions by reducing the inclusion of ineffective or irrelevant peptides, ultimately increasing the chances of identifying candidates suitable for vaccine development or immunotherapy.
> >
> {: .solution}
>
{: .question}

## Refining Results Using Table Operations

### Pivoting the table to aggregate affinity scores
In this step, the Table Compute tool is used to refine and aggregate the data from the filtered MHC-I peptide binding results. Specifically, the tool is performing a Pivot operation, which transforms the data for easier analysis and interpretation. The tool organizes the data by a specific index, such as icore, which could represent a unique identifier or a feature of interest in the dataset. The column to be analyzed is allele, which groups the data by different alleles (i.e., the genetic variants of the MHC molecules). The values of interest, which represent the binding affinity scores or ranks (e.g., percentile_rank), are then aggregated using the Maximum function, ensuring that for each group of alleles, the highest binding affinity is retained.

This operation is important because it helps prioritize the peptides that exhibit the strongest binding potential for each allele. By focusing on the best binders per allele, the data becomes more manageable and relevant for further analysis, such as identifying the most promising neoantigen candidates.

This Table Compute tool is used to perform a data manipulation operation called Pivot. This operation transforms the table format and aggregates data based on specific parameters to provide clearer insights. 
Specifically:
- The Index is set to icore, which could refer to a specific identifier or feature within the dataset. This helps group the data by this particular attribute.
- The Column parameter is set to allele, so the data will be organized based on different alleles of the MHC molecules, grouping them accordingly.
- The Values being analyzed are the percentile_rank, which represent the binding affinities of the peptides.
- The Aggregator Function is set to Maximum, meaning that for each allele, only the highest binding affinity (best percentile rank) is retained in the output.

> <hands-on-title> Weak-Table Compute </hands-on-title>
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


> <hands-on-title> Strong-Table Compute </hands-on-title>
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
> 1. Why is the Pivot operation used in this step?
> 2. Why do we use the 'Maximum' aggregator function?
>
> > <solution-title></solution-title>
> >
> > 1. The Pivot operation is used to reorganize the data so that each allele is associated with the highest percentile rank (maximum value) for the peptides. This helps in focusing on the top-ranked peptides for each allele, ensuring that the analysis highlights the most promising candidates for immunotherapy or vaccine development.
> > 2. The 'Maximum' function ensures that, for each allele, we keep the peptide with the highest binding affinity (best percentile rank). This step prioritizes the most relevant peptides for further analysis, eliminating lower-affinity binders that are less likely to have therapeutic potential.
> >
> {: .solution}
>
{: .question}

### Removing headers

In this sub-step, the Remove Beginning tool is used to clean the data by removing unnecessary rows or headers from the start of the table. After performing the Pivot operation in the previous step, the dataset may include extra header rows or metadata that aren't needed for analysis. This tool helps streamline the data by removing these initial rows, ensuring that only relevant information remains for further processing.

By applying the Remove Beginning tool, the user ensures that any unwanted starting rows—such as those containing column names, labels, or metadata that might have been carried over from previous operations—are removed, leaving the dataset clean and ready for the next analysis step.

> <hands-on-title> Weak- Remove beginning </hands-on-title>
>
> 1. {% tool [Remove beginning](Remove beginning1) %} with the following parameters:
>    - {% icon param-file %} *"from"*: `table` (output of **Table Compute** {% icon tool %})
>
>
>
{: .hands_on}

> <hands-on-title> Strong- Remove beginning </hands-on-title>
>
> 1. {% tool [Remove beginning](Remove beginning1) %} with the following parameters:
>    - {% icon param-file %} *"from"*: `table` (output of **Table Compute** {% icon tool %})
>
>
>
{: .hands_on}


## Extract Peptide column from the tabular 

> <hands-on-title> Weak Peptide extraction  </hands-on-title>
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1`
>    - {% icon param-file %} *"From"*: `out_file1` (output of **Remove beginning** {% icon tool %})
>
>
>
{: .hands_on}


> <hands-on-title> Strong Peptide extraction </hands-on-title>
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1`
>    - {% icon param-file %} *"From"*: `out_file1` (output of **Remove beginning** {% icon tool %})
>
>
{: .hands_on}


## Novel peptide verification with PepQuery2

In this step, the PepQuery2 tool is used to validate the identified peptides by checking them against a protein reference database to verify whether they are novel or already known. This validation is a critical step in proteomics workflows, as it ensures the authenticity and novelty of the peptides identified in earlier steps. PepQuery2 is utilized to perform novel peptide/protein validation, which involves checking peptide sequences against a reference protein database to confirm if the peptide is novel (not previously identified in the database) or known (matches a protein sequence in the database). The PepQuery2 tool is used in the final step to validate whether the peptides identified in earlier stages are novel or known. This is an essential step to confirm the relevance and accuracy of the findings in a broader biological context. If a peptide is novel, it could represent a potential new target for further research. If it is already known, it helps contextualize the results within existing knowledge. By comparing the identified peptides to a comprehensive database of known human proteins, this step ensures that only relevant and novel peptides are prioritized for further study, enhancing the overall quality and focus of the research.

> <hands-on-title> Weak Peptide verification - PepQuery 2 </hands-on-title>
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


> <hands-on-title> Strong Peptide verification - PepQuery 2 </hands-on-title>
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


## Filtering confident PepQuery2 validated peptides and annotating their binding affinity to HLA.

The first step is to filter the peptides based on the confidence column. The confident peptides are then annotated with the corresponding HLA they bind to.

### Filtering confident peptides 
> <hands-on-title> Filter weak confident peptides </hands-on-title>
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `psm_rank_txt` (output of **PepQuery2** {% icon tool %})
>    - *"With following condition"*: `c20=='Yes'`
>
>
>
{: .hands_on}


> <hands-on-title> Filter strong confident peptides </hands-on-title>
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `psm_rank_txt` (output of **PepQuery2** {% icon tool %})
>    - *"With following condition"*: `c20=='Yes'`
>
>
>
{: .hands_on}


### Annotating peptides based on their binding affinity

The peptides, including both strong and weak binders, are annotated with their respective HLA types.

 > <hands-on-title> Annotating weak binder peptides </hands-on-title>
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


> <hands-on-title> Annotating strong binder peptides </hands-on-title>
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

Please note that all the software tools used in this workflow are subject to version updates and changes. As a result, the parameters, functionalities, and outcomes may differ with each new version. Additionally, if the protein sequences are downloaded at different times, the number of sequences may also vary due to updates in the reference databases or tool modifications. We recommend the users to verify the specific versions of software tools used to ensure the reproducibility and accuracy of results.
