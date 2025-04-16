---
layout: tutorial_hands_on

title: Data submission using ENA upload Tool
questions:
- How to prepare sequences for submission to ENA?
- How to upload raw sequences to ENA?
objectives:
- Manage sequencing files (ab1, FASTQ, FASTA, FASTQ.GZ)
- Clean sequences in an automated and reproducible manner
- Perform alignments for each sequence
- Have the necessary sequence format to submit to ENA
- Submit raw reads to ENA using the ENA upload Tool
time_estimation: 2h
key_points:
- Clean raw ab1 sequences and compare filtered sequences to NCBI nucleotidic database
- Submit cleaned and unique sequences to European Nucleotide Archive (ENA) resource
contributions:
    authorship:
        - Najatamk
        - yvanlebras
    testing:
        - PaulineSGN
        - yvanlebras
    editing:
        - VerenaMoo
    funding:
        - pndb

subtopic: ecologymetadatamgt
---


This tutorial will guide you through the necessary steps to manage and prepare sequencing files (ab1, FASTQ, FASTA) for submission to the genomic database ENA.
This workflow will take you from raw sequences in AB1 format through all the necessary steps to integrate these sequences into the ENA genomic database. We will convert the files into FASTQ and FASTA formats after performing quality control.
Additionally, we will perform alignments with the NCBI database to ensure the accuracy of your sequences.You will then need to fill a metadata Excel template to use the ENA upload Tool.
The worklow is made of 17 Galaxy tools, we will  present them and explain what they do.
The goal is to present an accessible and reproductible workflow for data submission.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Prepare raw data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. **Create a new history** for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. **Import** the raw sequences files.
>
>    ```
>    https://data.indores.fr/api/access/datafile/3673
>    https://data.indores.fr/api/access/datafile/3609
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. **Rename** {% icon galaxy-pencil %} your datafiles
>    - `3673` becomes `A2_RC_8F2_B.pl_HCOI.ab1`
>    - `3609` becomes `A12_RC_9G4_B.md_HCOI.ab1`
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
> 4. **Check the datatype**
>    - Make sure it is `ab1`, and change it if not.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md %}
>
> 5. **Build a Collection** containing these two files, you can ame i "ab1" for example
>
>    {% snippet faqs/galaxy/collections_build_list.md %}
>
{: .hands_on}


## Tools used in the "Prepare Data submission" Workflow

Following steps take as input ab1 sequences files and produce filtered FastQ and Fasta files so sequences passing the quality checks are compared to NCBI nucleotidic database using Blastn operation.

### Converting Ab1 files to FASTQ

> <hands-on-title> ab1 to FASTQ converter </hands-on-title>
>
> 1. {% tool [ab1 to FASTQ converter](toolshed.g2.bx.psu.edu/repos/ecology/ab1_fastq_converter/ab1_fastq_converter/1.20.0) %} with the following parameters:
>    - {% icon param-collection %} *"Input ab1 file"*: `ab1` data collection created at the previous step
>
{: .hands_on}


### Quality Control

We are doing a first Quality control on the raw files using [Falco](https://falco.readthedocs.io/en/latest/) and MultiQC. Falco is an efficiency optimized rewrite of [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

> <hands-on-title> Falco </hands-on-title>
> 1. {% tool [Falco](toolshed.g2.bx.psu.edu/repos/iuc/falco/falco/1.2.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Raw read data from your current history"*: `ab1.fastq` data collection created at the previous step
>
> 2. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.11+galaxy1) %} with the following parameters:
>    - In *"Results"*:
>        - {% icon param-repeat %} *"Insert Results"*
>            - *"Which tool was used generate logs?"*: `FastQC` (no matter whether FastQC or Falco was used)
>    - In *"FastQC output"*:
>         - {% icon param-file %} *"RawData FastQC output"*: `Falco on collection X:` data collection created at the previous step
>
> 3. *Check on the HTML files the general quality statistics of your sequences*
>
{: .hands_on}


> <question-title> Question </question-title>
>
> 1. What is the quality of your sequences?
> 2. Do you have adapters?
>
> > <solution-title></solution-title>
> >
> > 1. Quality is quite good looking at the "status checks" section of MultiQC. As expected (because here we only have one sequence by file) "Per base sequence Content" and "Overrepresented sequences" Sections are "bad" for both sequences files. "adapter content" section also show a "bad" result for A2_RC_8F2_B.pl_HCOI.ab1 file.
> >
> > 2. A2_RC_8F2_B.pl_HCOI.ab1 file seems to have adapters in it.
> >
> {: .solution}
{: .question}

# Cleaning the Data

## Cutadapt

Cutadapt enables the removal of adapters, polyA tails, and other artifacts from sequences. The tool also filters reads based on quality.

> <hands-on-title> Cutadapt </hands-on-title>
>
> 1. {% tool [Cutadapt](toolshed.g2.bx.psu.edu/repos/lparsons/cutadapt/cutadapt/4.8+galaxy0) %} with the following parameters:
>     - {% icon param-collection %} *"FASTQ/A file"*: the collection with your data (output of {% icon tool %} **ab1 to FastQ converter**)
>     - **"Single-end or Paired-end reads?"**: `Single-end`
>     - In **"Other Read Trimming Options"**:
>         - **"Quality cutoff(s) (R1)"**: `30`
>         - **"Shortening reads to a fixed length"**: `Disabled`
>
> > <comment-title> Suggestions </comment-title>
> >
> > You may consider changing these parameters depending on the quality of your dataset.
> {: .comment}
>
{: .hands_on}

> <comment-title> Quality Control </comment-title>
>
> We do a second quality control similar to the first one to check the quality of the sequences after cleaning them.
{: .comment}


## Quality Control with Falco and MultiQC

> <hands-on-title> Falco </hands-on-title>
>
> 1. {% tool [Falco](toolshed.g2.bx.psu.edu/repos/iuc/falco/falco/1.2.4+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Raw read data from your current history"*: output from {% icon tool%} **Cutadapt**
>
> 2. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.11+galaxy1) %} with the following parameters:
>    - In *"Results"*:
>        - *"Which tool was used generate logs?"*: `FastQC` (no matter whether FastQC or Falco was used)
>        - {% icon param-repeat %} *"Insert FastQC output"*
>          - {%  icon param-collection %} *"FastQC output"*: the `raw` output from {% icon tool %} **Falco**
>
>    > <comment-title> Comment </comment-title>
>    >
>    > You should notice an improvement on the quality of your sequences.
>    {: .comment}
>
{: .hands_on}

## Filtering the collection


> <hands-on-title> Filter empty datasets </hands-on-title>
>
> 1. {% tool [Filter empty datasets](__FILTER_EMPTY_DATASETS__) %} with the following parameters
>    - {% icon param-collection %} *"Input Collection"*: output collection from Cutadapt step
>
> 2. {% tool [FASTQ Groomer](toolshed.g2.bx.psu.edu/repos/devteam/fastq_groomer/fastq_groomer/1.1.5+galaxy2) %} with the following parameters:
>    - {% icon param-collection %} *"File to groom"* : output collection from the {% icon tool %} **Filter empty datasets**
>
>    This step is notably there to produce "standardized" fastqsanger sequences files so we can then use other tools accepting only such data format.
>
> 3. {% tool [Filter FASTQ](toolshed.g2.bx.psu.edu/repos/devteam/fastq_filter/fastq_filter/1.1.5) %} with the following parameters:
>    - *"FASTQ File"*: output collecton from  {% icon tool %} **FastQ Groomer**
>    - *"Minimum size"*: `300`
>
>    > <comment-title> Comment </comment-title>
>    >
>    > Here we descide to only keep sequences of 300bp or above, you may change this parameter depending on your dataset
>    {: .comment}
>
{: .hands_on}


### Changing files names

> <hands-on-title>Extract element identifiers and remove extensions</hands-on-title>
>
> 1. {% tool [Extract element identifiers](toolshed.g2.bx.psu.edu/repos/iuc/collection_element_identifiers/collection_element_identifiers/0.0.2) %}
>    - {% icon param-collection %} *"Dataset collection"*: output from the previous step
>
> 2. {% tool [Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.3) %} with the following parameters:
>    - *"Select lines from"*: output of the previous step
>    - In *"Check"*:
>      - {% icon param-repeat %} *"Insert Check"*
>        - *"Find Regex"*: `.ab1`
>        - *"Replacement"*: ``
>
>    > <comment-title> Comment </comment-title>
>    >
>    > This is to ensure that all your files names end with .fastq.gz
>    {: .comment}
>
> 3. {% tool [Paste](Paste1) %} with the following parameters:
>    - {% icon param-file %} *"Paste"*: the file from {% icon tool %} **Extract element identifiers**
>    - {% icon param-file %} *"and"*: the file from {% icon tool %} **Regex Find And Replace**
>    - {% icon param-select %} *"Delimited by"*: Tab
>
> 4. **Check the datatype**
>    - should be 'tabular'. If not, change it now.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md %}
>
{: .hands_on}



> <hands-on-title> Relabel identifiers </hands-on-title>
>
> 1. {% tool [Relabel identifiers](__RELABEL_FROM_FILE__) %} with the following parameters:
>    - {% icon param-collection %} *"Input Collection"*: output from {% icon tool %} **Filter FastQ**
>    - *"How should the new labels be specified?"*: `Map original identifiers to new ones using a two column table.`
>
{: .hands_on}


## Alignments on NCBI database

> <hands-on-title> NCBI BLAST alignment </hands-on-title>
>
> 1. {% tool [FASTQ to FASTA](toolshed.g2.bx.psu.edu/repos/devteam/fastqtofasta/fastq_to_fasta_python/1.1.5) %} with the following parameters:
>    - {% icon param-collection %} *"Input FASTQ File"*: output collection from {% icon tool %} **Relabel Identifiers**
>
> 2. {% tool [NCBI BLAST+ blastn](toolshed.g2.bx.psu.edu/repos/devteam/ncbi_blast_plus/ncbi_blastn_wrapper/2.14.1+galaxy2) %} with the following parameters:
>    - {% icon param-collection %} Nucleotide query sequence(s): output from the previous step
>    - *"Subject database/sequences"*: `Locally installed BLAST database`
>        - *"Nucleotide BLAST database"*: `NCBI NT (01 Sep 2023)`
>    - *"Output format"*: `Tabular (extended 25 columns)`
>    - *"Advanced Options"*: `Hide Advanced Options`
{: .hands_on}


> <hands-on-title> Extracting best hits </hands-on-title>
>
> 1. {% tool [Unique](toolshed.g2.bx.psu.edu/repos/bgruening/unique/bg_uniq/0.3) %} with the following parameters:
>    - {% icon param-collection %} *"File to scan for unique values"*: output from the previous step
>    - *"Advanced Options"*: `Show Advanced Options`
>        - *"Column start"*: `c1`
>        - *"Column end"*: `c1`
>
{: .hands_on}

## Workflow Outputs

1. **Collection of raw FASTQ files:** Input AB1 files converted into FASTQ files.

2. **Collection of FASTQ files (after quality control)**: Renamed Fastq files ready for submission after quality control and filtering.

3. **Collection of FASTA files**: FASTQ files converted into FASTA format. Used for conducting BLAST alignments.

4. **Falco Quality Control Results** before and after cleaning: Both raw Falco results and HTML reports are created

5. **MultiQC Quality Control Results** before and after cleaning: Both raw MultiQC statistics and HTML report are created

6. **Raw Blast Results**: Results of BLAST alignments conducted on our sequences. Columns names are:

   | Column | NCBI name    | Description                                    |
   |--------|--------------|------------------------------------------------|
   | 1      | qaccver      | Query accession dot version                    |
   | 2      | saccver      | Subject accession dot version (database hit)   |
   | 3      | pident       | Percentage of identical matches                |
   | 4      | length       | Alignment length                               |
   | 5      | mismatch     | Number of mismatches                           |
   | 6      | gapopen      | Number of gap openings                         |
   | 7      | qstart       | Start of alignment in query                    |
   | 8      | qend         | End of alignment in query                      |
   | 9      | sstart       | Start of alignment in subject (database hit)   |
   | 10     | send         | End of alignment in subject (database hit)     |
   | 11     | evalue       | Expectation value (E-value)                    |
   | 12     | bitscore     | Bit score                                      |
   | 13     | sallseqid    | All subject Seq-id(s), separated by a ';'      |
   | 14     | score        | Raw score                                      |
   | 15     | nident       | Number of identical matches                    |
   | 16     | positive     | Number of positive-scoring matches             |
   | 17     | gaps         | Total number of gaps                           |
   | 18     | ppos         | Percentage of positive-scoring matches         |
   | 19     | qframe       | Query frame                                    |
   | 20     | sframe       | Subject frame                                  |
   | 21     | qseq         | Aligned part of query sequence                 |
   | 22     | sseq         | Aligned part of subject sequence               |
   | 23     | qlen         | Query sequence length                          |
   | 24     | slen         | Subject sequence length                        |
   | 25     | salltitles   | All subject title(s), separated by a '<>'      |

7. **Filtered Blast Results**
Files containing the closest homologous sequences.

8. **Collection of Fastq files**
Contains filtered sequences.

# How to use ENA upload Tool

## Adding ENA "Webin" credentials to your Galaxy user information

> <comment-title> Having an ENA Submission Account </comment-title>
>
>  Make sure you have a submission account with the European Nucleotide Archive (ENA). You will need the identifier and the password, available through https://www.ebi.ac.uk/ena/submit/webin/login.
>
{: .comment}


>  <hands-on-title> Add your "WEBIN" credentials to your Galaxy account </hands-on-title>
>  **Instructions:**
>	- From the Menu, click on "User" > "Preferences". Click on "Manage Information". Scroll down to "Your ENA Webin account details" and enter your ENA "Webin" identifier and password.
>  ![Adding ENA Webin credentials](./images/imageI.png)
{: .hands_on}

## Submitting using a metadata template file

For this tutorial we will use the ENA default sample checklist.

![Excel Metadata template](./images/image2.png)

**Note:** It is crucial to fill in all the fields marked "Mandatory" and ensure that the sequence names match exactly those indicated in the Excel file.

> <comment-title> ENA Metadata Templates </comment-title>
>
> You can find metadata templates for each checklist in the [ELIXIR-Belgium GitHub repository](https://github.com/ELIXIR-Belgium/ENA-metadata-templates)
>
> 1. Direct download link of the [ENA default sample checklist]( https://github.com/ELIXIR-Belgium/ENA-metadata-templates/raw/main/templates/ERC000011/metadata_template_ERC000011.xlsx)
>
> 2. Direct download link of the [ENA default sample checklist filled with elements for the training](https://github.com/galaxyproject/training-material/raw/24776cf161e38ac0449755749d23e851400020aa/topics/ecology/tutorials/ENA_Biodiv_submission/metadata_GdBqCOI_ERC000011_Test.xlsx)
>
> You will need to import this file into your Galaxy history. Then, use the ENA Upload Tool to proceed with the submission.
>
{: .comment}



> <hands-on-title> Excel Metadata Template </hands-on-title>
>
> 1. Import the ENA default sample checklist file.
>
>    ```
>    https://github.com/galaxyproject/training-material/raw/24776cf161e38ac0449755749d23e851400020aa/topics/ecology/tutorials/ENA_Biodiv_submission/metadata_GdBqCOI_ERC000011_Test.xlsx
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 2. {% tool [ENA Upload tool](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.11+galaxy1) %} with the following parameters:
>    - *"Action to execute"*: `Add new (meta)data`
>    - *"Select the metadata input method"*: `Excel file`
>    - *"Select the ENA sample checklist"*: `ENA default sample checklist (ERC000011)`
>    - *"Select Excel file based on template"*: `metadata_GdBqCOI_ERC000011_Test.xlsx`
>    - *"Select input data"*: `Dataset or dataset collection`
>    - *"Add .fastq (.gz, .bz2) extension to the Galaxy dataset names to match the ones described in the input tables?"*: `Yes`
>
>    > <comment-title> Datatype </comment-title>
>    >
>    > The ENA upload tool will then automatically compress fastq sequences files into .fastq.gz format before submission
>    {: .comment}
>
>    > <warning-title> Danger: Submit to ENA test server! </warning-title>
>    > We suggest you first submit to the [ENA test server](https://wwwdev.ebi.ac.uk/ena/submit/webin/) before making a public submission! Submission can be seen in `Dashboard/Study Report`
>    {: .warning}
>
> ![ENA Upload tool](./images/image3.png)
>
{: .hands_on}



# Conclusion

This tutorial guides you through quality check and preparing raw data files for ENA submission. You can then verify that your sequences have been successfully sent by logging into the Test ENA portal (https://wwwdev.ebi.ac.uk/ena/submit/webin/login) and navigating to the Study Report section.
