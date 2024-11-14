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

Verification is an essential step in the process of analyzing neoantigens to ensure that the identified peptides or proteins are accurate and reliable. Neoantigens, which are novel antigens arising from tumor-specific mutations, are crucial for cancer immunotherapy. However, the prediction and identification of these neoantigens must be rigorously validated to ensure they are truly present and capable of eliciting an immune response. Without proper verification, the therapeutic potential of neoantigens cannot be realized, leading to wasted resources and failed treatments.

In this tutorial, we will guide you through a series of bioinformatics tools used to verify novel peptides and proteins predicted from mass spectrometry data. These tools help confirm the presence and specificity of the identified neoantigens by comparing them to protein databases, checking for validation against spectral data, and performing sequence analysis.


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## PepQuery2 Verification Workflow

This workflow involves several key steps to validate neoantigens identified through mass spectrometry data. Each step leverages bioinformatics tools to ensure that the predicted peptides or proteins are accurate and can be used for further therapeutic research. Below are the detailed steps of this workflow:

### 1. Get Data
The first step in the workflow involves gathering and uploading the necessary mass spectrometry (MS) data files. These files, which include raw spectra or peptide data, are essential for identifying potential neoantigens. Proper data organization and renaming ensure that the workflow executes smoothly.

### 2. Preprocessing Raw MS Data with msconvert
Next, the raw MS data files are preprocessed using msconvert to convert them into a suitable format for analysis, specifically the Mascot Generic Format (MGF). This step is critical for transforming the data into a structured format that can be analyzed using peptide validation tools.

### 3. Peptide Validation with PepQuery2
Once the data is preprocessed, PepQuery2 is employed to validate the peptides. This tool compares the identified peptides against a protein reference database, using MS/MS data to confirm their authenticity. It ensures that the identified peptides are truly valid and not false positives, which is essential for identifying reliable neoantigens.

### 4. Filtering Tabular Data
After peptide validation, Query Tabular is used to filter the results. This step narrows down the data by focusing on specific criteria such as peptide match status and spectrum quality, ensuring only the most relevant peptides are retained for further analysis.

### 5. Converting Data to FASTA Format
The validated and filtered peptides are then converted into a FASTA format using the Tabular-to-FASTA tool. This format is commonly used for sequence analysis and provides a standardized structure for further processing of the data.

### 6. BLAST-P Sequence Alignment
In the next step, the peptides are subjected to sequence alignment using NCBI BLAST+ blastp against a protein database. This step helps confirm whether the peptides are novel or match known proteins, adding an extra layer of validation to ensure their potential as therapeutic neoantigens.

### 7. Refining Data with Query Tabular
Following the BLAST alignment, Query Tabular is used again to filter out peptides that do not meet the criteria or are not novel. This step helps ensure that only the most promising neoantigen candidates remain for further investigation.

### 8. Summary of Results
The final step involves summarizing the results, which include the validated and refined peptides identified as potential neoantigens. These findings are organized and presented for interpretation, providing valuable insights into the therapeutic potential of the identified neoantigens for immunotherapy.

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

# Preprocessing Raw MS Data with msconvert
msconvert is a tool used to preprocess raw mass spectrometry (MS) data files by converting them into a compatible format for downstream analysis. In this workflow, msconvert is employed to convert the input unrefined MS data into Mascot Generic Format (MGF), a standardized format for storing peptide and protein data. Additionally, during this step, peak picking is applied to identify the most relevant peaks from the raw spectra, which improves the quality of the data. However, the charge state recalculation is disabled as it's not necessary for the subsequent analysis. By converting the raw data into MGF format, msconvert prepares the dataset for further processing in tools like PepQuery2 and BLAST, ensuring that the data is in the right format for peptide identification and validation.

## Converting RAW to MGF using **msconvert**

> <hands-on-title> msconvert</hands-on-title>
>
> 1. {% tool [msconvert](toolshed.g2.bx.psu.edu/repos/galaxyp/msconvert/msconvert/3.0.20287.2) %} with the following parameters:
>    - {% icon param-file %} *"Input unrefined MS data"*: `output` (Input dataset)
>    - *"Do you agree to the vendor licenses?"*: `Yes`
>    - *"Output Type"*: `mgf`
>    - In *"Data Processing Filters"*:
>        - *"Apply peak picking?"*: `Yes`
>        - *"(Re-)calculate charge states?"*: `no`
>
>
{: .hands_on}


> <question-title></question-title>
>
> 1. What does "peak picking" mean in the context of msconvert?
> 2. Why is the charge state recalculation turned off in msconvert?
>
> > <solution-title></solution-title>
> >
> > 1. Peak picking refers to the process of identifying the most significant peaks in the mass spectrometry data, which correspond to the ions of interest. This step is important because it helps filter out noise and retain only the relevant peaks that will be used for peptide identification.
> > 2. The charge state recalculation is turned off because the data being processed is already assumed to have accurate charge assignments from the initial raw spectra. Recalculating the charge states could lead to unnecessary complications or errors in the context of this workflow, which focuses on peptide identification rather than charge state assignment.
> >
> {: .solution}
>
{: .question}

## Peptide Validation with PepQuery2 with **PepQuery2**

PepQuery2 is a tool used to validate novel peptides and proteins by searching mass spectrometry data against a reference protein database. In this workflow, PepQuery2 is used for peptide validation, specifically focusing on the detection of novel peptides or proteins. The input data includes a list of peptides (obtained from a previous step) and a protein reference database, which is also provided from the workflow history. The spectrum file, generated by msconvert, is used to search for peptide-spectrum matches (PSMs). The tool is configured with specific parameters, such as peptide charge and length, and the fragmentation method (CID/HCD). The output will help identify novel peptides that may not be present in the reference database, aiding in the discovery of potential new biomarkers or proteins.

> <hands-on-title> PepQuery2 </hands-on-title>
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
>
{: .hands_on}



> <question-title></question-title>
>
> 1. Why is "Non enzyme" selected as the digestion method in PepQuery2?
> 2. What is the significance of selecting "CID/HCD" as the fragmentation method?
>
> > <solution-title></solution-title>
> >
> > 1. "Non enzyme" is selected because, in this workflow, the peptides are already prepared and may not require digestion by a specific enzyme. This option allows for the validation of the peptides without applying any enzyme-specific cleavage rules.
> > 2. CID (Collision-Induced Dissociation) and HCD (Higher Energy Collisional Dissociation) are fragmentation techniques used to break peptides into smaller ions for analysis. These methods are commonly employed to improve the identification of peptides in mass spectrometry by generating clearer and more interpretable spectra, which is essential for accurate peptide validation.
> >
> {: .solution}
>
{: .question}

## Filtering Tabular Data with **Query Tabular**

Query Tabular is a tool used to query tabular datasets using SQL-like commands. In this step, the tool is used to filter and extract specific data from the results generated by PepQuery2 (specifically the psm_rank_txt dataset). The query provided selects specific columns (e.g., c1 and c4) from the table where a condition(e.g., c20 = 'Yes') is met. This allows the user to filter the peptide-spectrum matches (PSMs) based on certain criteria, such as identifying PSMs that were validated with a "Yes" condition. The tool helps organize the data into a more manageable and relevant format for further analysis.

> <hands-on-title> Query tabular </hands-on-title>
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
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Why is the option "include query result column headers" set to 'No'?
> 2. What does the SQL query SELECT c1, c4 FROM t1 WHERE (c20 = 'Yes') do?
>
> > <solution-title></solution-title>
> >
> > 1. The option is set to 'No' because, in this case, the user may not need column headers in the output file. This can be useful when the output is being further processed or integrated into another system where the headers are not required.
> > 2. This SQL query selects specific columns (c1 and c4) from the table t1 where the value in column c20 is equal to 'Yes'. Essentially, it filters the dataset to retrieve only the rows where a certain condition (e.g., peptide validation) has been met.
> >
> {: .solution}
>
{: .question}

## Converting Data to FASTA Format with **Tabular-to-FASTA**
Tabular-to-FASTA is a tool used to convert tabular data into the FASTA format, commonly used for storing sequence data. In this step, the tool takes the output of the Query Tabular step, which contains the filtered peptide data, and converts it into a FASTA format. The parameters specify which columns in the tabular data correspond to the sequence (c1) and the title (c['2']). The output will be a FASTA file where each peptide sequence is accompanied by its associated title or identifier. This conversion is crucial for further analysis, such as alignment or database searches, using sequence data in FASTA format. This database is being generated for BLAST-P searching.

> <hands-on-title> Tabular to FASTA </hands-on-title>
>
> 1. {% tool [Tabular-to-FASTA](toolshed.g2.bx.psu.edu/repos/devteam/tabular_to_fasta/tab2fasta/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Tab-delimited file"*: `output` (output of **Query Tabular** {% icon tool %})
>    - *"Title column(s)"*: `c['2']`
>    - *"Sequence column"*: `c1`
>
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Why is there a need to convert the files to FASTA for BLASTP searches?
>
> > <solution-title></solution-title>
> >
> > 1. We need to convert data to FASTA format for BLASTP because it is the required input format for the tool. FASTA provides a standardized way to represent protein sequences with a header and sequence lines, making it compatible with BLASTP’s alignment algorithms. This format allows BLASTP to properly parse, compare, and search the query sequence against protein databases.
> >
> {: .solution}
>
{: .question}

## BLAST-P Sequence Alignment
In this step, the NCBI BLAST+ blastp tool is used for performing protein sequence alignment. It compares the query protein sequence (generated from the Tabular-to-FASTA tool) against a locally installed BLAST database of protein sequences. The alignment is optimized for shorter protein queries (less than 30 residues) using the blastp-short option. The tool uses a PAM30 scoring matrix, along with specific gap costs, and generates tabular output with extended 25 columns, allowing users to examine the results in detail.

**Parameters:** 
- Protein query sequence(s): The input protein sequence in FASTA format.
- Subject database/sequences: The BLAST database against which the query sequence will be compared.
- Type of BLAST: Optimized for short protein queries.
- Expectation value cutoff: Specifies the cutoff for statistical significance.
- Advanced options: Includes scoring matrix (PAM30), gap costs, and word size for the BLAST algorithm.

> <hands-on-title> NCBI BLAST+ blastp </hands-on-title>
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
>
{: .hands_on}

> <question-title></question-title>
>
> 1. What is the significance of using the blastp-short option in this tool?
> 2. Why is the PAM30 matrix used in this step?
>
> > <solution-title></solution-title>
> >
> > 1. The blastp-short option is used because the query protein sequence is shorter than 30 residues, and this setting optimizes the alignment for such sequences.
> > 2. The PAM30 matrix is a substitution matrix suitable for shorter protein sequences and is specifically designed for those with fewer than 50 amino acids.
> >
> {: .solution}
>
{: .question}

## Refining Data with **Query Tabular**
The Query Tabular tool is used in this workflow to filter and organize data from previous steps. In this case, it processes the output from the NCBI BLAST+ blastp tool, which contains protein sequence alignment results. The tool allows the user to execute SQL-like queries on the tabular data, enabling filtering and sorting based on specific criteria. For example, the query selects distinct peptides from the pep table that are not perfectly aligned (i.e., with a sequence identity of 100%) and have certain mismatches or gaps in the alignment. The filtered results are then used for further analysis.

The tool helps to refine the data by removing sequences that meet specific criteria, such as perfect alignment, and ensuring that the remaining sequences meet the necessary quality standards for further exploration.

> <hands-on-title> Query Tabular</hands-on-title>
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
>
{: .hands_on}


> <question-title></question-title>
>
> 1. Why is the query filtering out sequences with 100% identity from the BLAST results?
> 2. What is the significance of using SQL-like queries in Query Tabular?
>
> > <solution-title></solution-title>
> >
> > 1. Sequences with 100% identity are filtered out to focus on those that show variability in the alignment, which might indicate novel or significant biological differences worth exploring.
> > 2. SQL-like queries allow for efficient data manipulation and filtering, enabling users to extract only the most relevant results based on complex conditions (such as mismatches, gaps, or alignment length).
> >
> {: .solution}
>
{: .question}


# Conclusion

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.

# Disclaimer 

Please note that all the software tools used in this workflow are subject to version updates and changes. As a result, the parameters, functionalities, and outcomes may differ with each new version. Additionally, if the protein sequences are downloaded at different times, the number of sequences may also vary due to updates in the reference databases or tool modifications. We recommend the users to verify the specific versions of software tools used to ensure the reproducibility and accuracy of results.
