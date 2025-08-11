---
layout: tutorial_hands_on

title: "Clinical Metaproteomics 3: Verification"
zenodo_link: "https://doi.org/10.5281/zenodo.10105821"
questions:
- Why do we need to verify our identified peptides
- What is the importance of making a new database for quantification
objectives:
- Verification of peptides helps in confirming the presence of the peptides in our samplle
time_estimation: 3H
key_points:
- Perform verification
- Extraction of accession numbers for getting protein sequences
contributions:
  authorship:
    - subinamehta
    - katherine-d21
    - dechendb
  editing:
    - pratikdjagtap
    - timothygriffin
requirements:
  -
    type: "internal"
    topic_name: proteomics
subtopic: clinical-metaproteomics
follow_up_training:

    -
        type: "internal"
        topic_name: proteomics
        tutorials:
            - clinical-mp-4-quantitation
tags: [label-TMT11]
redirect_from:
- /topics/proteomics/tutorials/clinical-mp-data-verification/tutorial

recordings:
- captioners:
  - katherine-d21
  date: '2024-06-21'
  galaxy_version: '23.1'
  length: 16M
  youtube_id: wCVzWMQANGE
  speakers:
  - katherine-d21
---


In proteomic research, the primary goal is to obtain accurate and meaningful insights into the proteome of a biological system. Verifying the presence of peptides or proteins is a critical step in achieving this goal, ensuring the quality and reliability of the data and the biological relevance of the findings. This tutorial is a sequel to the clinical metaproteomics discovery workflow. Once you have identified microbial peptides, the next step is to verify these peptides, for which we use PepQuery.

The PepQuery tool is used to validate the identified microbial peptides from SearchGUI/PeptideShaker and MaxQuant, to ensure that they are indeed of microbial origin and that human peptides were not misassigned. To do this, all confident microbial peptides from the two database search algorithms were merged and searched against the Human UniProt Reference proteome (with Isoforms) and cRAP databases.

Interestingly, the PepQuery tool does not rely on searching peptides against a reference protein sequence database as “traditional” shotgun proteomics does, which enables it to identify novel, disease-specific sequences with sensitivity and specificity in its protein validation (Figure A). Then we extract microbial protein sequences that are assigned to the PepQuery verified peptides. To this, we again add the Human UniProt Reference proteome (with Isoforms) and cRAP databases for creating a database for quantitation purposes (Figure B).

![Peptide Verification]({% link topics/proteomics/images/clinical-mp/clinical-mp-verification-1.JPG %})

![Database generation from verified peptides]({% link topics/proteomics/images/clinical-mp/clinical-mp-verification-2.JPG %})


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Get data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    https://zenodo.org/records/10105821/files/PTRC_Skubitz_Plex2_F10_9Aug19_Rage_Rep-19-06-08.mgf
>    https://zenodo.org/records/10105821/files/PTRC_Skubitz_Plex2_F11_9Aug19_Rage_Rep-19-06-08.mgf
>    https://zenodo.org/records/10105821/files/PTRC_Skubitz_Plex2_F13_9Aug19_Rage_Rep-19-06-08.mgf
>    https://zenodo.org/records/10105821/files/PTRC_Skubitz_Plex2_F15_9Aug19_Rage_Rep-19-06-08.mgf
>    https://zenodo.org/records/10105821/files/SGPS_Peptide_Report.tabular
>    https://zenodo.org/records/10105821/files/MaxQuant_Peptide_Report.tabular
>    https://zenodo.org/records/10105821/files/Distinct_Peptides_for_PepQuery.tabular
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
> 5. Add to each database a tag corresponding to input files.
> 6. Users can create a database collection of the MGF files.
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}


# Import Workflow

> <hands-on-title>Running the Workflow</hands-on-title>
>
> 1. **Import the workflow** into Galaxy:
>
>    {% snippet faqs/galaxy/workflows_run_trs.md path="topics/proteomics/tutorials/clinical-mp-3-verification/workflows/WF3_Verification_Workflow.ga" title="Verification Workflow" %}
>
> 2. Run **Workflow** {% icon workflow %} using the following parameters:
>    - *"Send results to a new history"*: `No`
>    - {% icon param-file %} *" Input MGFs Dataset Collection "*: `MGF dataset collection`
>    - {% icon param-file %} *" SGPS_peptide-report"*: `SGPS_Peptide_Report.tabular`
>    - {% icon param-file %} *" Distinct Peptides for PepQuery"*: `Distinct_Peptides_for_PepQuery.tabular`
>    - {% icon param-file %} *" MaxQuant-peptide-report "*: `MaxQuant_Peptide_Report.tabular`
>
>    {% snippet faqs/galaxy/workflows_run.md %}
>
{: .hands_on}


# Extraction of Microbial Peptides from SearchGUI/PeptideShaker and MaxQuant
Now that we have identified microbial peptides from SearchGUI/PeptideShaker and MaxQuant, we need to extract the microbial peptide sequences and group them to obtain a list of distinct microbial peptides. This list of distinct peptides will be used as input for PepQuery2 to verify confident microbial peptides.

First, we will use the Cut tool to select the peptide and protein columns from the SearchGUI/PeptideShaker and MaxQuant Peptide Reports. Then we use Remove header lines from SGPS and MaxQuant to prepare for concatenation with Remove beginning.

> <hands-on-title> Extracting peptides </hands-on-title>
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c6,c2`
>    - {% icon param-file %} *"From"*: `output` (Input dataset)
>
>
> 2. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1,c35`
>    - {% icon param-file %} *"From"*: `output` (Input dataset)
>
>
> 3. {% tool [Remove beginning](Remove beginning1) %} with the following parameters:
>    - {% icon param-file %} *"from"*: `out_file1` (output of **Cut** {% icon tool %})
>
>
> 4. {% tool [Remove beginning](Remove beginning1) %} with the following parameters:
>    - {% icon param-file %} *"from"*: `out_file1` (output of **Cut** {% icon tool %})
>
>
{: .hands_on}


## Concatenate peptides from MaxQuant and SGPS for PepQuery2

We will now concatenate the peptide and protein datasets from SearchGUI/PeptideShaker and MaxQuant. Later, we will generate a list of confident peptides using PepQuery2. The list of confident peptides will be searched against the concatenated peptide-protein datasets from SearchGUI/PeptideShaker and MaxQuant to generate a list of verified peptides.

> <hands-on-title> Concatenate SGPS and MaxQuant peptides </hands-on-title>
>
> 1. {% tool [Concatenate datasets](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cat/0.1.1) %} with the following parameters:
>    - {% icon param-files %} *"Datasets to concatenate"*: `out_file1` (output of **Remove beginning** {% icon tool %}), `out_file1` (output of **Remove beginning** {% icon tool %})
>
>
{: .hands_on}


## Creating input database for PepQuery2

We generate and merge Human UniProt (with Isoforms) and contaminants (cRAP) to make an input database for PepQuery2.

> <hands-on-title> FASTA Merge Files and Filter Unique Sequences</hands-on-title>
>
> 1. {% tool [FASTA Merge Files and Filter Unique Sequences](toolshed.g2.bx.psu.edu/repos/galaxyp/fasta_merge_files_and_filter_unique_sequences/fasta_merge_files_and_filter_unique_sequences/1.2.0) %} with the following parameters:
>    - *"Run in batch mode?"*: `Merge individual FASTAs (output collection if input is collection)`
>        - In *"Input FASTA File(s)"*:
>            - {% icon param-repeat %} *"Insert Input FASTA File(s)"*
>                - {% icon param-file %} *"FASTA File"*: `Human UniProt+Isoforms FASTA` (output of **Protein Database Downloader** {% icon tool %})
>                - - {% icon param-file %} *"FASTA File"*: `cRAP database` (output of **Protein Database Downloader** {% icon tool %})
>
>
{: .hands_on}



## Peptide verification
The PepQuery2 tool will be used to validate the identified microbial peptides from SearchGUI/PeptideShaker and MaxQuant to ensure that they are indeed of microbial origin and that human peptides were not misassigned. We will use the list of `Distinct Peptides` (from the [Discovery Module](https://github.com/subinamehta/training-material/blob/main/topics/proteomics/tutorials/clinical-mp-discovery/tutorial.md)), `Human UniProt+Isoforms+cRAP database`, and our `MGF file collection` as inputs for PepQuery2. The outputs we are interested in are the four PSM Rank (txt) files (one for each MGF file).

Interestingly, the PepQuery2 tool does not rely on searching peptides against a reference protein sequence database as “traditional” shotgun proteomics does, which enables it to identify novel, disease-specific sequences with sensitivity and specificity in its protein validation. More information about PepQuery is available, including the first {% cite Wen2019 %} and second iterations {% cite Wen2023 %}.


> <hands-on-title> Peptide verification </hands-on-title>
>
> 1. {% tool [PepQuery2](toolshed.g2.bx.psu.edu/repos/galaxyp/pepquery2/pepquery2/2.0.2+galaxy0) %} with the following parameters:
>    - *"Validation Task Type"*: `novel peptide/protein validation`
>    - In *"Input Data"*:
>        - *"Input Type"*: `peptide`
>            - *"Peptides?"*: `Peptide list from your history`
>                - {% icon param-file %} *"Peptide Sequences (.txt)"*: `output` (Input dataset)
>        - *"Protein Reference Database from"*: `history`
>            - {% icon param-file %} *"Protein Reference Database File"*: `output` (output of **FASTA Merge Files and Filter Unique Sequences** {% icon tool %})
>        - *"MS/MS dataset to search"*: ` Spectrum Datasets from history`
>            - {% icon param-collection %} *"Spectrum File"*: `output` (Input dataset collection)
>        - *"Report Spectrum Scan as"*: `spectrum title in MGF`
>    - In *"Modifications"*:
>        - *"Fixed modification(s)"*: `1: Carbamidomethylation of C [57.02146372057]` `13: TMT 11-plex of K [229.16293213472]` `14: TMT 11-plex of peptide N-term [229.16293213472]`
>        - *"Variable modification(s)"*: `2: Oxidation of M [15.99491461956]`
>        - *"Use more stringent criterion for unrestricted modification searching"*: `Yes`
>        - *"Consider amino acid substitution modifications?"*: `Yes`
>    - In *"Digestion"*:
>        - *"Enzyme"*: `Trypsin`
>        - *"Max Missed Cleavages"*: `2`
>    - In *"Mass spectrometer"*:
>        - In *"Tolerance"*:
>            - *"Precursor Tolerance"*: `10`
>            - *"Precursor Unit"*: `ppm`
>            - *"Tolerance"*: `0.6`
>        - In *"PSM"*:
>            - *"Fragmentation Method"*: `CID/HCD`
>            - *"Scoring Method"*: `HyperScore`
>            - *"Minimum Charge"*: `2`
>            - *"Maximum Charge"*: `6`
>
>
{: .hands_on}

# Collapsing all the data

Remember that PepQuery2 generates a PSM Rank file for each input MGF file, so we will have four PSM Rank files. To make the analysis more efficient, we will collapse these four PSM Rank files into one dataset.


> <hands-on-title> Collasping PSM rank files into a singular dataset using Collapse Collection </hands-on-title>
>
> 1. {% tool [Collapse Collection](toolshed.g2.bx.psu.edu/repos/nml/collapse_collections/collapse_dataset/5.1.0) %} with the following parameters:
>    - {% icon param-file %} *"Collection of files to collapse into single dataset"*: `psm_rank_txt` (output of **PepQuery2** {% icon tool %})
>    - *"Keep one header line"*: `Yes`
>
>
{: .hands_on}



## Filtering out confident peptides

Now, we want to filter for confident peptides from PepQuery2 and prepare them for the Query Tabular tool.

> <hands-on-title> Filter </hands-on-title>
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `output` (output of **Collapse Collection** {% icon tool %})
>    - *"With following condition"*: `c20=='Yes'`
>    - *"Number of header lines to skip"*: `1`
>
>
{: .hands_on}


> <hands-on-title> Remove header line from filtered PepQuery peptides with Remove beginning </hands-on-title>
>
> 1. {% tool [Remove beginning](Remove beginning1) %} with the following parameters:
>    - {% icon param-file %} *"from"*: `out_file1` (output of **Filter** {% icon tool %})
>
>
{: .hands_on}


> <hands-on-title> Cut (select out) peptide sequences from PepQuery output with Cut </hands-on-title>
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1`
>    - {% icon param-file %} *"From"*: `out_file1` (output of **Remove beginning** {% icon tool %})
>
>
{: .hands_on}


##  **Querying verified peptides**

We will use the Query Tabular tool {% cite Johnson2019 %} to search the PepQuery-verified peptides against the concatenated dataset that contains peptides and proteins from SearchGUI/Peptide and MaxQuant. This step ensures all the PepQuery-verified peptides are assigned to their protein/protein groups.

> <hands-on-title> Querying verified peptides </hands-on-title>
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.3.0) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `out_file1` (output of **Cut** {% icon tool %})
>            - In *"Table Options"*:
>                - *"Specify Name for Table"*: `pep`
>                - *"Specify Column Names (comma-separated list)"*: `mpep`
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `out_file1` (output of **Concatenate datasets** {% icon tool %})
>            - In *"Table Options"*:
>                - *"Specify Name for Table"*: `prot`
>                - *"Specify Column Names (comma-separated list)"*: `pep,prot`
>    - *"SQL Query to generate tabular output"*: `select pep.mpep, prot.prot
>    `FROM pep
>    INNER JOIN  prot on pep.mpep=prot.pep `
>    - *"include query result column headers"*: `Yes`
>    `
>
>
> > <comment-title>SQL Query information</comment-title>
> >  The query input files are the list of peptides and the peptide report we obtained from MaxQuant and SGPS. The query is matching each peptide (m.pep) from the PepQuery results to the peptide reports so that each verified peptide has its protein/protein group assigned to it.
> {: .comment}
{: .hands_on}

> <hands-on-title> Remove Header with Remove beginning </hands-on-title>
>
> 1. {% tool [Remove beginning](Remove beginning1) %} with the following parameters:
>    - {% icon param-file %} *"from"*: `output` (output of **Query Tabular** {% icon tool %})
>
>
{: .hands_on}


Using the Group tool, we can select distinct (unique) peptides and proteins from the Query Tabular tool.

> <hands-on-title> Extract distinct peptides with Group </hands-on-title>
>
> 1. {% tool [Group](Grouping1) %} with the following parameters:
>    - {% icon param-file %} *"Select data"*: `out_file1` (output of **Remove beginning** {% icon tool %})
>    - *"Group by column"*: `c1`
>    - In *"Operation"*:
>        - {% icon param-repeat %} *"Insert Operation"*
>            - *"Type"*: `Concatenate Distinct`
>            - *"On column"*: `c2`
>
>
{: .hands_on}


## Retrieve UniProt IDs for distinct peptides

Again, we will use the Query Tabular tool to retrieve UniProt IDs (accession numbers) for the distinct (grouped) peptides.

> <hands-on-title> Query Tabular </hands-on-title>
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.3.0) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `out_file1` (output of **Group** {% icon tool %})
>            - In *"Filter Dataset Input"*:
>                - In *"Filter Tabular Input Lines"*:
>                    - {% icon param-repeat %} *"Insert Filter Tabular Input Lines"*
>                        - *"Filter By"*: `normalize list columns, replicates row for each item in list`
>                            - *"enter column numbers to normalize"*: `2`
>                            - *"List item delimiter in column"*: `;`
>                    - {% icon param-repeat %} *"Insert Filter Tabular Input Lines"*
>                        - *"Filter By"*: `regex replace value in column`
>                            - *"enter column number to replace"*: `2`
>                            - *"regex pattern"*: `(tr|sp)[|]`
>                    - {% icon param-repeat %} *"Insert Filter Tabular Input Lines"*
>                        - *"Filter By"*: `regex replace value in column`
>                            - *"enter column number to replace"*: `2`
>                            - *"regex pattern"*: `[ ]+`
>                    - {% icon param-repeat %} *"Insert Filter Tabular Input Lines"*
>                        - *"Filter By"*: `regex replace value in column`
>                            - *"enter column number to replace"*: `2`
>                            - *"regex pattern"*: `[|].*$`
>            - In *"Table Options"*:
>                - *"Specify Name for Table"*: `t1`
>                - *"Use first line as column names"*: `Yes`
>                - *"Specify Column Names (comma-separated list)"*: `pep,prot`
>    ` *"SQL Query to generate tabular output"*: `SELECT distinct(prot) AS Accession
>       from t1`
>        - *"include query result column headers"*: `No`
>
>
{: .hands_on}


> <question-title></question-title>
>
> 1. What is the accession number of a protein?
> 2. Can there be multiple accession numbers for one peptide or protein?
>
> > <solution-title></solution-title>
> >
> > 1. An accession number of a protein, also called a protein accession number, is a unique identifier assigned to a specific protein sequence in a protein sequence database. These accession numbers are used to reference and catalog proteins in a standardized and systematic manner
> >
> > 2. Yes, a single peptide or protein can have multiple accession numbers, particularly when dealing with different protein sequence databases, databases for specific species, or different versions of the same database. That's the reason in our workflow we merge both accession and sequences.
> {: .solution}
>
{: .question}


## Generate FASTA database from UniProt IDs
Using the UniProt IDs from Query Tabular, we will be able to generate a FASTA database for our PepQuery-verified peptides.

> <hands-on-title> UniprotXML-downloader </hands-on-title>
>
> 1. {% tool [UniProt](toolshed.g2.bx.psu.edu/repos/galaxyp/uniprotxml_downloader/uniprotxml_downloader/2.4.0) %} with the following parameters:
>    - *"Select"*: `A history dataset with a column containing Uniprot IDs`
>        - {% icon param-file %} *"Dataset (tab separated) with ID column"*: `output` (Input dataset)
>        - *"Column with ID"*: `c1`
>        - *"Field"*: `Accession`
>    - *"uniprot output format"*: `fasta`
>
>
{: .hands_on}


## Generating compact database

Lastly, we will merge the Human UniProt (with isoforms), contaminants (cRAP) and the PepQuery-verified FASTA databases into one Quantitation Database that will be used as input for the [Quantification Module](https://github.com/subinamehta/training-material/blob/main/topics/proteomics/tutorials/clinical-mp-quantitation/tutorial.md).

> <hands-on-title> Generation of Compact Verified Database with UniProt </hands-on-title>
>
> 1. {% tool [FASTA Merge Files and Filter Unique Sequences](toolshed.g2.bx.psu.edu/repos/galaxyp/fasta_merge_files_and_filter_unique_sequences/fasta_merge_files_and_filter_unique_sequences/1.2.0) %} with the following parameters:
>    - *"Run in batch mode?"*: `Merge individual FASTAs (output collection if input is collection)`
>        - In *"Input FASTA File(s)"*:
>            - {% icon param-repeat %} *"Insert Input FASTA File(s)"*
>                - {% icon param-file %} *"FASTA File"*: `proteome` (output of **UniProt** {% icon tool %})
>
>
{: .hands_on}


# Conclusion

A peptide verification workflow is a critical step in proteomic research that enhances data reliability, quantitative accuracy, and biological understanding by confirming the presence and validity of selected peptides. It is a pivotal quality control process that ensures the trustworthiness of proteomic findings and supports downstream investigations. By completing this tutorial, you have not only verified the microbial peptides but also created a database consisting of protein sequences from the PepQuery-verified peptides. The need of such a database is to ensure that when we quantify our proteins and peptides we are reducing the introduction of false positives. This database will be now used for quantitation purposes.
