---
layout: tutorial_hands_on

title: "Neoantigen 2: Non-normal-Database-Generation"
zenodo_link: ""
questions:
- Why must we generate a customized fusion database for Proteogenomics research?
objectives:
- Downloading databases related to 16SrRNA data
- For better neoantigen identification results.
time_estimation: 3H
key_points:
- Create a customized Variant proteomics database from 16SrRNA results.
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
- /topics/proteomics/tutorials/neoantigen-1-non-normal-database-generation/tutorial

---


# Introduction

Proteogenomics leverages mass spectrometry (MS)-based proteomics data alongside genomics and transcriptomics data to identify neoantigens—unique peptide sequences arising from tumor-specific mutations. In the initial section of this tutorial, we will construct a customized protein database (FASTA) using RNA-sequencing files (FASTQ) derived from tumor samples. Following this, we will conduct sequence database searches using the resultant FASTA file and MS data to identify peptides corresponding to novel proteoforms, specifically focusing on potential neoantigens. We will then assign genomic coordinates and annotations to these identified peptides and visualize the data, assessing both spectral quality and genomic localization.

In this framework, Proteogenomics incorporates RNA-Seq data to generate tailored protein sequence databases, enabling the identification of protein sequence variants, including neoantigens, through mass spectrometry analysis ({% cite Chambers_2017 %}).


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Neoantigen-Non-normal-database generation

This tutorial guides users through the process of generating a non-normal variant database. It encompasses essential bioinformatics steps to identify and prepare variant-specific peptides for immunological studies. Below is an overview of each major stage:

### 1. Get Data
The workflow begins with uploading raw sequencing data, followed by a quality assessment to ensure data integrity. This step establishes a solid foundation for subsequent analyses by addressing any issues in the initial dataset.

### 2. Variant Detection and Mapping
Next, the RNA sequencing data is aligned to a reference genome using tools like HISAT2 and StringTie. Alignment events are detected with specialized tools like Freebayes, CustomProDB, and GFFcompare, which identify non-normal gene transcripts. These tools analyze the resulting alignments to characterize the gene segments in CDS, single nucleotide variants, indels, UTRs, or frameshifts.

### 3. Text reformatting and Database generation
Once variants are identified, we generate a customized database and apply various reformatting techniques to tag it, ensuring optimal usability for downstream processing.

### 4. Addition of known protein sequences
Known proteomics databases are added to the variant database to create a comprehensive database.

### 5. Final Database Construction
The workflow concludes with applying regex adjustments and other formatting functions to standardize the output. This process culminates in creating a comprehensive database of potential non-normal protein sequences, making them ready for experimental validation and clinical exploration.


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


## **Convert compressed file to uncompressed.**

Uncompressing data is a crucial first step in many bioinformatics workflows because raw sequencing data files, especially from high-throughput sequencing, are often stored in compressed formats (such as `.gz` or `.zip`) to save storage space and facilitate faster data transfer. Compressed files need to be uncompressed to make the data readable and accessible for analysis tools, which generally require the data to be in plain text or other compatible formats. By uncompressing these files, we ensure that downstream applications can efficiently process and analyze the raw sequencing data without compatibility issues related to compression. In this workflow, we do that for both forward and reverse files.

> <hands-on-title> Converting compressed to uncompressed </hands-on-title>
>
> 1. {% tool [Convert compressed file to uncompressed.](CONVERTER_gz_to_uncompressed) %} with the following parameters:
>    - {% icon param-file %} *"Choose compressed file"*: `output` (Input dataset)
>
>
{: .hands_on}

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Convert compressed file to uncompressed.](CONVERTER_gz_to_uncompressed) %} with the following parameters:
>    - {% icon param-file %} *"Choose compressed file"*: `output` (Input dataset)
>
>
{: .hands_on}

# Extracting Single amino acid variants with HISAT and Freebayes

## Aligning to the reference genome with **HISAT2**
HISAT2 is a fast and efficient tool used in bioinformatics workflows for aligning sequence reads to a reference genome. In this task, HISAT2 is being utilized to align paired-end reads against the human genome version GRCh38 (hg38). This alignment is essential for downstream analyses such as variant calling or transcript quantification. HISAT2 is configured to use default alignment and scoring options to ensure simplicity and speed, which is often suitable for general-purpose analyses.

In this workflow, HISAT2 serves the critical role of mapping raw sequencing data (reads) to a reference genome. This step is a foundation for understanding genetic variation and gene expression in the sample. By aligning the reads to a reference, HISAT2 provides a structured output that can be further analyzed in various bioinformatics applications.


> <hands-on-title> HISAT2 </hands-on-title>
>
> 1. {% tool [HISAT2](toolshed.g2.bx.psu.edu/repos/iuc/hisat2/hisat2/2.2.1+galaxy1) %} with the following parameters:
>    - *"Source for the reference genome"*: `Use a built-in genome`
>        - *"Select a reference genome"*: `Human Dec. 2013 (GRCh38/hg38) (hg38)`
>    - *"Is this a single or paired library"*: `Paired-end`
>        - {% icon param-file %} *"FASTA/Q file #1"*: `output1` (output of **Convert compressed file to uncompressed.** {% icon tool %})
>        - {% icon param-file %} *"FASTA/Q file #2"*: `output1` (output of **Convert compressed file to uncompressed.** {% icon tool %})
>        - *"Paired-end options"*: `Use default values`
>    - In *"Advanced Options"*:
>        - *"Input options"*: `Use default values`
>        - *"Alignment options"*: `Use default values`
>        - *"Scoring options"*: `Use default values`
>        - *"Spliced alignment options"*: `Use default values`
>        - *"Reporting options"*: `Use default values`
>        - *"Output options"*: `Use default values`
>        - *"SAM options"*: `Use default values`
>        - *"Other options"*: `Use default values`
>
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Why do we need to align sequencing reads to a reference genome in this workflow?
> 2. What do paired-end reads represent, and why are they important in this alignment?
>
> > <solution-title></solution-title>
> >
> > 1. Aligning sequencing reads to a reference genome is essential for identifying the location of each read within the genome, which enables analysis of gene expression levels, detection of mutations, and understanding of overall genome structure in the sample.
> > 2. Paired-end reads are sequences generated from both ends of a DNA fragment, providing more information about the fragment's size and alignment accuracy. They improve alignment quality and facilitate more accurate mapping, especially for repetitive regions, by providing context from both ends of each read.
> >
> {: .solution}
>
{: .question}

## Variant Calling with **FreeBayes**
FreeBayes is a variant calling tool used in bioinformatics to identify genetic variants (such as SNPs, insertions, and deletions) from aligned sequence data. In this step, FreeBayes takes the output from the HISAT2 alignment (BAM or CRAM file) and identifies variations in the aligned reads compared to the human reference genome GRCh38. By running in "Simple diploid calling" mode, FreeBayes is configured to detect variants in diploid samples, assuming two sets of chromosomes as typically found in humans.


In this workflow, FreeBayes performs the essential function of variant calling, which is critical for identifying genetic differences that could be associated with diseases, traits, or other biological characteristics. The output from FreeBayes can then be used for downstream analyses such as functional annotation or association studies.

> <hands-on-title> FreeBayes </hands-on-title>
>
> 1. {% tool [FreeBayes](toolshed.g2.bx.psu.edu/repos/devteam/freebayes/freebayes/1.3.6+galaxy0) %} with the following parameters:
>    - *"Choose the source for the reference genome"*: `Locally cached`
>        - *"Run in batch mode?"*: `Run individually`
>            - {% icon param-file %} *"BAM or CRAM dataset"*: `output_alignments` (output of **HISAT2** {% icon tool %})
>        - *"Using reference genome"*: `Human Dec. 2013 (GRCh38/hg38) (hg38)`
>    - *"Limit variant calling to a set of regions?"*: `Do not limit`
>    - *"Read coverage"*: `Use defaults`
>    - *"Choose parameter selection level"*: `1. Simple diploid calling`
>
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Why is variant calling important after aligning reads to a reference genome?
> 2. What does "Simple diploid calling" mean in the context of FreeBayes?
>
> > <solution-title></solution-title>
> >
> > 1. Variant calling is important because it identifies specific genetic differences between the sample and the reference genome. These differences can provide insights into genetic variation, which may be associated with certain phenotypes, diseases, or treatment responses.
> > 2. "Simple diploid calling" refers to the assumption that each sample has two copies of each chromosome (diploid). FreeBayes uses this assumption to identify variants present in either one or both copies, which is appropriate for most human genetic analyses.
> >
> {: .solution}
>
{: .question}

## Customized Database generation using  **CustomProDB**
CustomProDB is a bioinformatics tool used to generate custom protein databases that incorporate specific genetic variants found in a sample. In this task, CustomProDB utilizes genome annotation data (GRCh38/hg38), a BAM file from HISAT2 alignment, and a VCF file from FreeBayes variant calling to create a variant-aware protein database. This database includes the sample’s unique variant proteins, which is essential for personalized proteomics research.

In this workflow, CustomProDB plays a critical role in translating genetic variants identified by FreeBayes into custom protein sequences. This variant-specific database is valuable for applications in proteomics, as it allows for the detection of variant-specific peptides in mass spectrometry data. The generated outputs, including a variant FASTA file and mapping files, support downstream analyses, such as studying how genetic variations may affect protein function or abundance.

> <hands-on-title> CustomProDB </hands-on-title>
> 
>
> 1. {% tool [CustomProDB](toolshed.g2.bx.psu.edu/repos/galaxyp/custom_pro_db/custom_pro_db/1.22.0) %} with the following parameters:
>    - *"Will you select a genome annotation from your history or use a built-in annotation?"*: `Use a built-in genome annotation`
>        - *"Select genome annotation"*: `Human (Ensembl 78 hsapiens) (hg38/GRCh38)`
>        - {% icon param-file %} *"BAM file"*: `output_alignments` (output of **HISAT2** {% icon tool %})
>        - {% icon param-file %} *"VCF file"*: `output_vcf` (output of **FreeBayes** {% icon tool %})
>    - *"Create a variant FASTA for short insertions and deletions"*: `Yes`
>    - *"Create SQLite files for mapping proteins to genome and summarizing variant proteins"*: `Yes`
>    - *"Create RData file of variant protein coding sequences"*: `Yes`
>
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Why is it necessary to create a custom protein database with genetic variants?
> 2. What is the purpose of creating SQLite and RData files in CustomProDB?
>
> > <solution-title></solution-title>
> >
> > 1. A custom protein database allows for the detection of sample-specific protein variants, which may affect protein structure and function. This is especially useful in personalized medicine, where identifying unique variants helps understand individual differences in disease mechanisms or treatment responses.
> > 2. The SQLite files map proteins to genome regions and summarize variant proteins, which supports efficient data storage and retrieval for variant analysis. The RData file contains coding sequences for variant proteins, which is useful for further statistical and bioinformatic analyses in R. Although we are not using this in our current workflow, the users can definitely use it to manipulate in their workflows.
> >
> {: .solution}
>
{: .question}

## Converting FASTA database from CustomPRODB to tabular for text processing with **FASTA-to-Tabular**
FASTA-to-Tabular is a tool that converts FASTA-formatted sequence files into tabular format, where each sequence entry is represented as a row with separate columns for identifiers and sequences. In this task, FASTA-to-Tabular processes the variant FASTA file generated by CustomProDB to create a tabular representation of the protein variants. This format is often easier to analyze or integrate into downstream data workflows.

In this workflow, FASTA-to-Tabular enables the conversion of variant protein sequences into a structured tabular format, which is helpful for subsequent data processing and analysis. This format allows researchers to efficiently filter, sort, or query specific sequence information and simplifies integration with other data analysis tools or databases. We do this for the indels, single nucleotide variants and rpkm databases.

> <hands-on-title> INDEL - FASTA-to-Tabular </hands-on-title>
>
> 1. {% tool [FASTA-to-Tabular](toolshed.g2.bx.psu.edu/repos/devteam/fasta_to_tabular/fasta2tab/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Convert these sequences"*: `output_indel` (output of **CustomProDB** {% icon tool %})
>
>
{: .hands_on}

> <hands-on-title> SNV - FASTA-to-Tabular </hands-on-title>
>
> 1. {% tool [FASTA-to-Tabular](toolshed.g2.bx.psu.edu/repos/devteam/fasta_to_tabular/fasta2tab/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Convert these sequences"*: `output_snv` (output of **CustomProDB** {% icon tool %})
>
>
{: .hands_on}

> <hands-on-title> RPKM - FASTA-to-Tabular </hands-on-title>
>
> 1. {% tool [FASTA-to-Tabular](toolshed.g2.bx.psu.edu/repos/devteam/fasta_to_tabular/fasta2tab/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Convert these sequences"*: `output_rpkm` (output of **CustomProDB** {% icon tool %})
>
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Why is it helpful to convert FASTA sequences to a tabular format?
> 2. What information is typically retained when converting a FASTA file to tabular format?
>
> > <solution-title></solution-title>
> >
> > 1. Converting FASTA sequences to tabular format makes it easier to visualize, filter, and manage the data. The tabular structure allows researchers to perform sorting, searching, and data manipulation more efficiently, which is beneficial for downstream analysis.
> > 2. The conversion usually retains the sequence identifiers and the sequences themselves. Each sequence entry in the FASTA file is represented in two columns: one for the identifier and another for the actual sequence, enabling easy access to both key elements of the data.
> >
> {: .solution}
>
{: .question}

## Manipulating the headers with **Column Regex Find And Replace**

Column Regex Find And Replace is a tool that applies regular expression (regex) patterns to specified columns in a tabular dataset to find and replace text patterns. In this task, the tool is used on the output from the FASTA-to-Tabular step to standardize or format specific patterns in the data. By applying regex patterns to column c1, the tool identifies and modifies specific sequence identifiers or annotations, making them easier to interpret or use in further analyses.

In this workflow, Column Regex Find And Replace cleans and formats the data in a way that makes identifiers or variant descriptions consistent. This is important for data compatibility, especially when the data needs to be used across different tools or integrated into larger datasets. It ensures that all sequence labels or variant annotations follow a uniform format, which reduces errors in downstream analyses.


> <hands-on-title> INDEL-Column Regex Find And Replace </hands-on-title>
>
> 1. {% tool [Column Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regexColumn1/1.0.3) %} with the following parameters:
>    - {% icon param-file %} *"Select cells from"*: `output` (output of **FASTA-to-Tabular** {% icon tool %})
>    - *"using column"*: `c1`
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `^(ENS[^_]+_\d+:)([ACGTacgt]+)>([ACGTacgt]+)\s*`
>            - *"Replacement"*: `\1\2_\3`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `,([A-Y]\d+[A-Y]?)\s*`
>            - *"Replacement"*: `.\1`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `^(ENS[^ |]*)\s*`
>            - *"Replacement"*: `\1`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `^(.*)$`
>            - *"Replacement"*: `generic|INDEL_\1`
>
>
{: .hands_on}

> <hands-on-title> SNV-Column Regex Find And Replace </hands-on-title>
>
> 1. {% tool [Column Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regexColumn1/1.0.3) %} with the following parameters:
>    - {% icon param-file %} *"Select cells from"*: `output` (output of **FASTA-to-Tabular** {% icon tool %})
>    - *"using column"*: `c1`
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `^(ENS[^_]+_\d+:)([ACGTacgt]+)>([ACGTacgt]+)\s*`
>            - *"Replacement"*: `\1\2_\3`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `,([A-Y]\d+[A-Y]?)\s*`
>            - *"Replacement"*: `.\1`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `^(ENS[^ |]*)\s*`
>            - *"Replacement"*: `\1`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `^(.*)$`
>            - *"Replacement"*: `generic|SNV_\1`
>
>
{: .hands_on}

> <hands-on-title> RPKM -Column Regex Find And Replace </hands-on-title>
>
> 1. {% tool [Column Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regexColumn1/1.0.3) %} with the following parameters:
>    - {% icon param-file %} *"Select cells from"*: `output` (output of **FASTA-to-Tabular** {% icon tool %})
>    - *"using column"*: `c1`
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `^(ENS[^_]+_\d+:)([ACGTacgt]+)>([ACGTacgt]+)\s*`
>            - *"Replacement"*: `\1\2_\3`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `,([A-Y]\d+[A-Y]?)\s*`
>            - *"Replacement"*: `.\1`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `^(ENS[^ |]*)\s*`
>            - *"Replacement"*: `\1`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `^(.*)$`
>            - *"Replacement"*: `generic|RPKM_\1`
>
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Why use regular expressions to modify text in a column?
> 2. What is the purpose of the generic|INDEL_\1 replacement pattern ( or any other pattern) in this tool?
>
> > <solution-title></solution-title>
> >
> > 1. Regular expressions allow for flexible, precise pattern matching, enabling the user to find and replace complex text patterns within a column. This is useful for standardizing data formats, extracting relevant portions, or appending additional information to entries.
> > 2. The generic|INDEL_\1 replacement pattern prepends generic|INDEL_ to each entry in the column, standardizing the label format and clarifying that each entry is related to an insertion or deletion (INDEL) variant. This makes it clear in downstream analysis that the sequence is variant-specific.
> >
> {: .solution}
>
{: .question}

## Converting the manipulated tabular files to FASTA with **Tabular-to-FASTA**
Tabular-to-FASTA is a tool that converts tabular data back into FASTA format, where each entry has a title and a sequence. In this step, the tool uses the processed output from the Column Regex Find And Replace step to create a FASTA file. Column c1 is set as the title (identifier) for each sequence, and column c2 contains the sequence data. This conversion is useful when the data needs to be returned to a standard sequence format for further bioinformatics analyses.

In this workflow, Tabular-to-FASTA converts the formatted tabular data back into a FASTA file, making it compatible with tools that require FASTA input for further analysis. This step enables the standardized, cleaned sequences from previous steps to be utilized in additional bioinformatics workflows or databases, maintaining the variant-specific information in a commonly used format. We do this for all the tabular files (SNV, INDEL, and RPKM).

> <hands-on-title> INDEL-Tabular-to-FASTA </hands-on-title>
>
> 1. {% tool [Tabular-to-FASTA](toolshed.g2.bx.psu.edu/repos/devteam/tabular_to_fasta/tab2fasta/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Tab-delimited file"*: `out_file1` (output of **Column Regex Find And Replace** {% icon tool %})
>    - *"Title column(s)"*: `c['1']`
>    - *"Sequence column"*: `c2`
>
>
{: .hands_on}

> <hands-on-title> SNV-Tabular-to-FASTA  </hands-on-title>
>
> 1. {% tool [Tabular-to-FASTA](toolshed.g2.bx.psu.edu/repos/devteam/tabular_to_fasta/tab2fasta/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Tab-delimited file"*: `out_file1` (output of **Column Regex Find And Replace** {% icon tool %})
>    - *"Title column(s)"*: `c['1']`
>    - *"Sequence column"*: `c2`
>
>
{: .hands_on}

> <hands-on-title> RPKM-Tabular-to-FASTA  </hands-on-title>
>
> 1. {% tool [Tabular-to-FASTA](toolshed.g2.bx.psu.edu/repos/devteam/tabular_to_fasta/tab2fasta/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Tab-delimited file"*: `out_file1` (output of **Column Regex Find And Replace** {% icon tool %})
>    - *"Title column(s)"*: `c['1']`
>    - *"Sequence column"*: `c2`
>
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Why is it important to convert tabular data back to FASTA format?
> 2. What information is typically included in the title and sequence columns in a FASTA file?
>
> > <solution-title></solution-title>
> >
> > 1. Converting tabular data back to FASTA format makes it compatible with a wide range of bioinformatics tools that require FASTA input. It allows the sequence data, now standardized and labeled, to be used in additional analyses, such as protein structure prediction or functional annotation.
> > 2. The title column usually contains a unique identifier or label for each sequence, often including relevant annotations. The sequence column contains the actual DNA, RNA, or protein sequence. Together, these provide essential details for recognizing and analyzing each sequence entry in downstream applications.
> >
> {: .solution}
>
{: .question}

## Merging Single amino acid variant databases with **FASTA Merge Files and Filter Unique Sequences**
FASTA Merge Files and Filter Unique Sequences is a tool that combines multiple FASTA files into a single file and removes any duplicate sequences, keeping only unique entries. In this task, the tool takes the FASTA file generated from the Tabular-to-FASTA step and merges it with any other FASTA files in the input list. The tool then filters the sequences to ensure that only unique sequences are retained in the final output, which is important for reducing redundancy in the dataset.

In this workflow, FASTA Merge Files and Filter Unique Sequences consolidate all sequence data into a single, non-redundant FASTA file. This step is essential for removing duplicate sequences, which helps streamline the dataset for further analysis. A unique sequence file reduces computational load and minimizes potential biases in downstream applications that could be affected by redundant data. We are merging the indel, snv, and rpkm databases in this step.


> <hands-on-title> FASTA Merge Files and Filter Unique Sequences </hands-on-title>
>
> 1. {% tool [FASTA Merge Files and Filter Unique Sequences](toolshed.g2.bx.psu.edu/repos/galaxyp/fasta_merge_files_and_filter_unique_sequences/fasta_merge_files_and_filter_unique_sequences/1.2.0) %} with the following parameters:
>    - *"Run in batch mode?"*: `Merge individual FASTAs (output collection if the input is a collection)`
>        - In *"Input FASTA File(s)"*:
>            - {% icon param-repeat %} *"Insert Input FASTA File(s)"*
>                - {% icon param-file %} *"FASTA File"*: `output` (output of **Tabular-to-FASTA** {% icon tool %})
>
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Why is it useful to filter out duplicate sequences from a FASTA file?
> 2. What does “Run in batch mode” mean in the context of merging FASTA files?
>
> > <solution-title></solution-title>
> >
> > 1. Filtering out duplicate sequences ensures that unique sequences are included, which reduces redundancy in the dataset. This can enhance computational efficiency, prevent overrepresentation of certain sequences, and ensure that analyses are based on a unique set of entries.
> > 2. In this context, “Run in batch mode” means merging multiple FASTA files into a single collection, even if the input consists of multiple individual FASTA files. This mode consolidates the files to make them easier to manage and analyze as a single dataset.
> >
> {: .solution}
>
{: .question}

# Extracting Non-Normal with Stringtie

## Sub-step with **StringTie**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [StringTie](toolshed.g2.bx.psu.edu/repos/iuc/stringtie/stringtie/2.2.3+galaxy0) %} with the following parameters:
>    - *"Input options"*: `Short reads`
>        - {% icon param-file %} *"Input short mapped reads"*: `output_alignments` (output of **HISAT2** {% icon tool %})
>    - *"Use a reference file to guide assembly?"*: `Use reference GTF/GFF3`
>        - *"Reference file"*: `Use a file from history`
>            - {% icon param-file %} *"GTF/GFF3 dataset to guide assembly"*: `output` (Input dataset)
>        - *"Output files for differential expression?"*: `No additional output`
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




## Sub-step with **GffCompare**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [GffCompare](toolshed.g2.bx.psu.edu/repos/iuc/gffcompare/gffcompare/0.12.6+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"GTF inputs for comparison"*: `output_gtf` (output of **StringTie** {% icon tool %})
>    - *"Use reference annotation"*: `Yes`
>        - *"Choose the source for the reference annotation"*: `History`
>            - {% icon param-file %} *"Reference annotation"*: `output` (Input dataset)
>    - *"Use sequence data"*: `No`
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



## Sub-step with **Convert gffCompare annotated GTF to BED**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Convert gffCompare annotated GTF to BED](toolshed.g2.bx.psu.edu/repos/galaxyp/gffcompare_to_bed/gffcompare_to_bed/0.2.1) %} with the following parameters:
>    - {% icon param-file %} *"GTF annotated by gffCompare"*: `transcripts_annotated` (output of **GffCompare** {% icon tool %})
>    - *"filter gffCompare class_codes to convert"*: ``
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





## Sub-step with **Translate BED transcripts**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Translate BED transcripts](toolshed.g2.bx.psu.edu/repos/galaxyp/translate_bed/translate_bed/0.1.0) %} with the following parameters:
>    - {% icon param-file %} *"A BED file with 12 columns"*: `output` (output of **Convert gffCompare annotated GTF to BED** {% icon tool %})
>    - *"Source for Genomic Sequence Data"*: `Locally cached twobit`
>        - *"Select reference 2bit file"*: `hg38`
>    - In *"Fasta ID Options"*:
>        - *"fasta ID source, e.g. generic"*: `generic`
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

## Sub-step with **bed to protein map**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [bed to protein map](toolshed.g2.bx.psu.edu/repos/galaxyp/bed_to_protein_map/bed_to_protein_map/0.2.0) %} with the following parameters:
>    - {% icon param-file %} *"A BED file with 12 columns, thickStart and thickEnd define protein coding region"*: `translation_bed` (output of **Translate BED transcripts** {% icon tool %})
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

# Extracting Protein Accession IDs


## Sub-step with **FASTA-to-Tabular**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [FASTA-to-Tabular](toolshed.g2.bx.psu.edu/repos/devteam/fasta_to_tabular/fasta2tab/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Convert these sequences"*: `output` (Input dataset)
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

## Sub-step with **Filter Tabular**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Filter Tabular](toolshed.g2.bx.psu.edu/repos/iuc/filter_tabular/filter_tabular/3.3.1) %} with the following parameters:
>    - {% icon param-file %} *"Tabular Dataset to filter"*: `output` (output of **FASTA-to-Tabular** {% icon tool %})
>    - In *"Filter Tabular Input Lines"*:
>        - {% icon param-repeat %} *"Insert Filter Tabular Input Lines"*
>            - *"Filter By"*: `select columns`
>                - *"enter column numbers to keep"*: `1`
>        - {% icon param-repeat %} *"Insert Filter Tabular Input Lines"*
>            - *"Filter By"*: `regex replace value in column`
>                - *"enter column number to replace"*: `1`
>                - *"regex pattern"*: `^[^|]+[|]([^| ]+).*$`
>                - *"replacement expression"*: `\1`
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

# Filtering RPKM accessions

## Sub-step with **Filter Tabular**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Filter Tabular](toolshed.g2.bx.psu.edu/repos/iuc/filter_tabular/filter_tabular/3.3.1) %} with the following parameters:
>    - {% icon param-file %} *"Tabular Dataset to filter"*: `output` (output of **FASTA-to-Tabular** {% icon tool %})
>    - In *"Filter Tabular Input Lines"*:
>        - {% icon param-repeat %} *"Insert Filter Tabular Input Lines"*
>            - *"Filter By"*: `select columns`
>                - *"enter column numbers to keep"*: `1`
>        - {% icon param-repeat %} *"Insert Filter Tabular Input Lines"*
>            - *"Filter By"*: `regex replace value in column`
>                - *"enter column number to replace"*: `1`
>                - *"regex pattern"*: `^([^ |]+).*$`
>                - *"replacement expression"*: `\1`
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

## Sub-step with **Concatenate datasets**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Concatenate datasets](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cat/9.3+galaxy1) %} with the following parameters:
>    - {% icon param-files %} *"Datasets to concatenate"*: `output` (output of **Filter Tabular** {% icon tool %}), `output` (output of **Filter Tabular** {% icon tool %})
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



## Sub-step with **FASTA Merge Files and Filter Unique Sequences**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [FASTA Merge Files and Filter Unique Sequences](toolshed.g2.bx.psu.edu/repos/galaxyp/fasta_merge_files_and_filter_unique_sequences/fasta_merge_files_and_filter_unique_sequences/1.2.0) %} with the following parameters:
>    - *"Run in batch mode?"*: `Merge individual FASTAs (output collection if input is collection)`
>        - In *"Input FASTA File(s)"*:
>            - {% icon param-repeat %} *"Insert Input FASTA File(s)"*
>                - {% icon param-file %} *"FASTA File"*: `translation_fasta` (output of **Translate BED transcripts** {% icon tool %})
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
