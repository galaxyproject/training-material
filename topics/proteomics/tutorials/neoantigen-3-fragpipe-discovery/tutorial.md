---
layout: tutorial_hands_on

title: Neoantigen_Database_merge_FragPipe_discovery
zenodo_link: ''
questions:
- What are the key features and unique sequences in protein datasets that contribute to neoantigen discovery?
- How can we identify neoantigens from proteomic data?
objectives:
- Understand the process of merging neoantigen databases with human protein sequences.
- Learn to use FragPipe for proteomics data analysis.
- Gain hands-on experience with bioinformatics tools such as FASTA file processing, database validation, and peptide identification.
time_estimation: 3H
key_points:
- Merging of non-normal databases with known human protein sequences.
- Processing and validation of FASTA files for proteomics analysis.
- Utilization of FragPipe to identify neoantigens and perform downstream analysis.
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
- proteomics/tutorials/neoantigen-3-fragpipe-discovery/tutorial

---


# Introduction

In this tutorial, we will guide you through a bioinformatics workflow aimed at merging neoantigen databases with known human protein sequences, preparing the data for proteomics analysis using FragPipe. This workflow involves processing FASTA files, filtering for unique sequences, and validating the FASTA databases before using FragPipe to perform peptide identification and validation of neoantigens.

Throughout the tutorial, you will learn how to integrate multiple datasets, ensuring that you can perform analyses such as the identification of potential neoantigens, which are critical for cancer immunotherapy and vaccine development. The tools and steps covered here are important for any bioinformatics pipeline dealing with proteomics and neoantigen discovery.


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Neoantigen Fragpipe Discovery

This tutorial guides users through the process of performing database searching or neoantigen protein/peptide discovery. It encompasses essential bioinformatics steps to identify and variant-specific peptides for immunological studies. Below is an overview of each major stage:

### 1. Get Data
The first step involves gathering and uploading the necessary proteomics data files into the analysis environment. These files typically contain protein sequences or raw spectrum data that will be processed throughout the tutorial. Proper data organization and tagging are essential to ensure smooth workflow execution.

### 2. Merging FASTA Files and Filtering for Unique Sequences
In this step, multiple FASTA files containing protein sequences are merged into a single file. After merging, sequences are filtered to retain only the unique ones, ensuring that redundancy is removed and only relevant protein data is used for downstream analysis.

### 3. Validating FASTA Databases
Once the FASTA files are merged and filtered, it's important to validate the database to ensure that the protein sequences are correctly formatted and usable for analysis. This step helps identify and correct any issues in the dataset before performing more complex analysis tasks.

### 4. Running FragPipe for Neoantigen Discovery
FragPipe, a proteomics analysis tool, is then employed to process the data further. This involves peptide identification, protein quantification, and running specialized workflows such as nonspecific-HLA searches to identify potential neoantigens that may be targeted for immunotherapy.

### 5. Collapsing Datasets and Selecting Relevant Data
After the analysis, the resulting datasets are collapsed into a single dataset to simplify further processing. This step helps streamline the data, making it easier to select and focus on the relevant sequences that match the biological question being addressed.

### 6. Querying Tabular Results for Further Analysis
In the final step, tabular results from the analysis are queried using SQL-like commands to filter and extract the most relevant data. This allows for focused analysis on specific protein sequences or neoantigens that were identified, enabling further downstream analysis and interpretation.

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

# Merging FASTA Files and Filtering for Unique Sequences

Next, we will merge the FASTA files, ensuring that any redundant sequences are removed. This step ensures that we only work with unique sequences, improving the quality and accuracy of the subsequent analysis. In this step, we combine the fusion database generated from the Arriba Pipeline (first neoantigen workflow) with the non-normal database created from HISAT, Freebayes, CustomPRODB, and the Stringtie Pipeline (second neoantigen workflow). Once merging is done, we validate the database to ensure that the sequences are in the right format.


## Merging fusion and non-normal databases with **FASTA Merge Files and Filter Unique Sequences**

> <hands-on-title> FASTA Merge Files and Filter Unique Sequences</hands-on-title>
>
> 1. {% tool [FASTA Merge Files and Filter Unique Sequences](toolshed.g2.bx.psu.edu/repos/galaxyp/fasta_merge_files_and_filter_unique_sequences/fasta_merge_files_and_filter_unique_sequences/1.2.0) %} with the following parameters:
>    - *"Run in batch mode?"*: `Merge individual FASTAs (output collection if input is collection)`
>        - In *"Input FASTA File(s)"*:
>            - {% icon param-repeat %} *"Insert Input FASTA File(s)"*
>                - {% icon param-file %} *"FASTA File"*: `output` (Input dataset)
>
>
{: .hands_on}


## Sequence database parsing with **Validate FASTA Database**

> <hands-on-title> Validate FASTA Database </hands-on-title>
>
> 1. {% tool [Validate FASTA Database](toolshed.g2.bx.psu.edu/repos/galaxyp/validate_fasta_database/validate_fasta_database/0.1.5) %} with the following parameters:
>    - {% icon param-file %} *"Select input FASTA dataset"*: `output` (output of **FASTA Merge Files and Filter Unique Sequences** {% icon tool %})
>
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Why is it important to validate a FASTA database before using it for further analysis in a proteomics pipeline?
> 2. What potential issues might be identified by the Validate FASTA Database tool?
>
> > <solution-title></solution-title>
> >
> > 1. Validating a FASTA database ensures that the sequences are correctly formatted and that no errors are present. It helps ensure the accuracy of the data used in downstream analysis, which is critical for correct protein identification and neoantigen discovery.
> > 2. The Validate FASTA Database tool might identify issues such as incorrect sequence formatting, missing information (e.g., headers or sequence data), or invalid characters in the sequences, which could cause errors during further analysis.
> >
> {: .solution}
>
{: .question}

##  Running FragPipe for Neoantigen Discovery
The FragPipe tool is used for processing and analyzing proteomics data in mass spectrometry experiments. This tool integrates multiple components such as MSFragger, IonQuant, and Percolator, allowing users to perform tasks such as peptide-spectrum match (PSM) identification, protein quantification, and post-translational modification (PTM) analysis. In this task, FragPipe is being used to identify potential neoantigens by analyzing mass spectrometry (MS) data and correlating it with a validated FASTA sequence database.
The FragPipe tool is executed with the MS data and a validated FASTA sequence database to identify peptides, proteins, and potential neoantigens from the raw mass spectrometry data. Users must provide MS spectrum files, a manifest file, and the validated FASTA database (created in the previous step using the Validate FASTA Database tool). The tool runs a series of analyses, including non-specific enzyme digestion, quantification of protein abundance, and validation of peptide-spectrum matches (PSMs) with Percolator and Protein Prophet. Output includes identified peptides, proteins, and quantification results, which are essential for downstream neoantigen discovery.

In this workflow, FragPipe is used after FASTA database validation to ensure that the sequence database is correctly formatted and ready for mass spectrometry analysis. It integrates several tools in a single workflow to process the raw MS data, identify peptides and proteins, and provide validation for the results. It also supports isobaric and label-free quantification methods for protein and peptide quantification, important for understanding relative protein abundance and identifying potential neoantigens.

> <hands-on-title> Fragpipe </hands-on-title>
>
> 1. {% tool [FragPipe -  Academic Research and Education User License (Non-Commercial)](toolshed.g2.bx.psu.edu/repos/galaxyp/fragpipe/fragpipe/20.0+galaxy2) %} with the following parameters:
>    - *"I understand that these tools, including MSFragger, IonQuant, Bruker, and Thermo Raw File Reader, are available freely for academic research and educational purposes only, and agree to the following terms."*: `Yes`
>    - {% icon param-file %} *"Proteomics Spectrum files"*: `output` (Input dataset)
>    - {% icon param-file %} *"Manifest file"*: `output` (Input dataset)
>    - {% icon param-file %} *"Proteomics Search Database in FASTA format"*: `goodFastaOut` (output of **Validate FASTA Database** {% icon tool %})
>    - *"Split database"*: `200`
>    - *"Workflow"*: `Nonspecific-HLA`
>        - In *"MSFragger Options"*:
>            - In *"Search Tolerances"*:
>                - *"Set Precursor Mass tolerances"*: `Use default`
>                - *"Set Precursor True tolerance"*: `Use default`
>                - *"Set Fragment Mass tolerances"*: `ppm`
>            - In *"In-silico Digestion Parameters"*:
>                - *"Protein Digestion Enzyme"*: `nonspecific`
>                - *"Second Enzyme Digest"*: `No`
>                - *"Maximum length of peptides to be generated during in-silico digestion"*: `20`
>            - In *"Variable Modifications"*:
>                - *"Select Variable Modifications"*: ``
>            - In *"Spectrum Processing"*:
>                - *"Precursor mass mode"*: `corrected`
>                - *"Precursor Charge Override"*: `Use default`
>        - In *"Validation"*:
>            - *"Run Validation"*: `Yes`
>                - *"PSM Validation"*: `Run Percolator`
>                - *"Run Protein Prophet"*: `Yes`
>                - *"Generate Philosopher Reports"*: `Yes`
>        - In *"Quant (MS1)"*:
>            - *"Perform Label-Free Quantification"*: `Use workflow default`
>        - In *"PTMs"*:
>            - *"Run PTM Shepherd"*: `no`
>        - In *"Quant (Isobaric)"*:
>            - *"Perform Isobaric Quantification"*: `no`
>    - *"Additional outputs"*: ``
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

## Sub-step with **Collapse Collection**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Collapse Collection](toolshed.g2.bx.psu.edu/repos/nml/collapse_collections/collapse_dataset/5.1.0) %} with the following parameters:
>    - {% icon param-file %} *"Collection of files to collapse into single dataset"*: `output_peptide` (output of **FragPipe -  Academic Research and Education User License (Non-Commercial)** {% icon tool %})
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

## Sub-step with **Select**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Select](Grep1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `output` (output of **Collapse Collection** {% icon tool %})
>    - *"that"*: `NOT Matching`
>    - *"the pattern"*: `(HUMAN)|(contam_)|(con_)|(ENSP)`
>    - *"Keep header line"*: `Yes`
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
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `out_file1` (output of **Select** {% icon tool %})
>            - In *"Table Options"*:
>                - *"Specify Name for Table"*: `fp`
>    - *"SQL Query to generate tabular output"*: `SELECT c1
FROM fp
WHERE (c16 IS NULL OR c16 = '')
AND (c18 IS NULL OR c18 = '')`
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
