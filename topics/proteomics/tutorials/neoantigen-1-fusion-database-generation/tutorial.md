---
layout: tutorial_hands_on

title: "Neoantigen 1: Fusion-Database-Generation"
zenodo_link: ""
questions:
- Why do we need to generate a customized fusion database for proteogenomics research?
objectives:
- Downloading databases related to 16SrRNA data
- For better neoantigen identification results.
time_estimation: 2H
key_points:
- Create a customized fusion proteomics database from 16SrRNA results.
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
            - neoantigen-2-non-normal-database-generation
tags: [label-free]
redirect_from:
- /topics/proteomics/tutorials/neoantigen-1-fusion-database-generation/tutorial

---


A neoantigen is a novel peptide (protein fragment) that is produced by cancer cells due to mutations, including gene fusions, that alter the DNA sequence in a way that generates unique proteins not found in normal cells. Because these mutated proteins are unique to the tumor, they are recognized as "foreign" by the immune system. Neoantigens are valuable in immunotherapy because they can serve as specific targets for the immune system, allowing treatments to selectively attack cancer cells while sparing normal tissue. By stimulating an immune response specifically against these neoantigens, therapies like cancer vaccines or T-cell-based treatments can be developed to enhance the body’s natural defense mechanisms, making neoantigens a promising avenue for personalized cancer treatment.

Creating a fusion database is essential in cancer genomics and personalized medicine, as it enables the identification of crucial biomarkers, enhances diagnostic accuracy, and supports therapeutic development. Gene fusions, where parts of two previously separate genes merge, can produce abnormal proteins that drive cancer. Cataloging these fusion events in a database helps researchers identify specific biomarkers linked to cancer types and design more targeted treatments. Additionally, fusion events may lead to unique peptide sequences, known as neoantigens, which are found only in cancer cells. These neoantigens can be targeted by the immune system, making fusion databases valuable in designing personalized immunotherapies like cancer vaccines or T-cell therapies. Some gene fusions also create oncogenic proteins that promote tumor growth, such as the BCR-ABL fusion in chronic myeloid leukemia. Including such information in a database aids in identifying potential therapeutic targets and predicting treatment efficacy. On the diagnostic side, known gene fusions serve as reliable markers, helping clinicians better classify cancer types and choose the most effective treatments. Finally, fusion databases provide a critical reference for researchers studying fusion mechanisms, their impact on disease progression, and their prevalence across cancers, ultimately fueling the discovery of novel treatments and therapies.

To generate the fusion database, RNA star and Arriba tools are used in this workflow.


![Fusions_Protein_Database workflow]({% link topics/proteomics/images/neoantigen/Fusions_Protein_Database.PNG %})

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Neoantigen Database Generation

## Overview of Fusion Neoantigen Database Workflow

The workflow in this tutorial guides users through the generation of a fusion neoantigen database, covering key steps in bioinformatics to identify, filter, and prepare fusion-specific peptides for further immunological study. Below is an overview of each major stage:

### 1. Get Data
The process begins with the upload and quality assessment of raw sequencing data, which is then uncompressed. This stage sets the groundwork for all subsequent analyses.

### 2. Fusion Detection and Alignment
RNA sequencing data undergoes alignment to a reference genome using tools like **RNA STAR**, followed by **Arriba** to detect fusion events. These tools identify gene fusions and help characterize the gene segments that combine to form new fusion genes.

### 3. Filtering and Refinement
After identifying fusions, various filters are applied to remove non-specific or common fusion events using blacklist data and other criteria. This step ensures that only relevant, unique fusion events are retained for neoantigen prediction.

### 4. Peptide Sequence Extraction and Formatting
Potential neoantigen peptides are extracted from the fusion gene sequences. Using tools such as **Text Reformatting** and **Tabular-to-FASTA**, the data is transformed into formats suitable for further immunological analysis.

### 5. Final Database Formatting
The workflow concludes by applying regex adjustments and formatting functions to standardize the output, creating a database of potential fusion neoantigens.

### Summary
This workflow provides a structured approach to preparing fusion neoantigen data for downstream applications, such as immunotherapy research, by making fusion-derived peptides accessible in a database for experimental or clinical exploration.


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

# Data preparation


## Convert compressed file to uncompressed

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


## Alignment with **RNA STAR**

**RNA STAR** (Spliced Transcripts Alignment to a Reference) is a high-performance tool used to align RNA sequencing (RNA-seq) reads to a reference genome. It identifies the best matches between RNA reads and genome sequences by detecting exon-exon junctions, which are critical for accurately mapping reads from spliced transcripts. RNA STAR uses a "two-pass" mapping approach that first identifies splice junctions across all reads and then uses these junctions to guide a more accurate alignment on the second pass. This capability is especially valuable for studying gene expression, discovering novel splice variants, and identifying fusion genes in cancer and other disease research. The output includes aligned sequences that can be used in subsequent steps of bioinformatics pipelines, such as fusion detection and differential expression analysis.

> <hands-on-title> Spliced transcripts Alignment to a human reference </hands-on-title>
>
> 1. {% tool [RNA STAR](toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.10b+galaxy4) %} with the following parameters:
>    - *"Single-end or paired-end reads"*: `Paired-end (as individual datasets)`
>        - {% icon param-file %} *"RNA-Seq FASTQ/FASTA file, forward reads"*: `output1` (output of **Convert compressed file to uncompressed.** {% icon tool %})
>        - {% icon param-file %} *"RNA-Seq FASTQ/FASTA file, reverse reads"*: `output1` (output of **Convert compressed file to uncompressed.** {% icon tool %})
>    - *"Custom or built-in reference genome"*: `Use a built-in index`
>        - *"Reference genome with or without an annotation"*: `use genome reference without builtin gene-model but provide a gtf`
>            - *"Select reference genome"*: `Human Dec. 2013 (GRCh38/hg38) (hg38)`
>            - {% icon param-file %} *"Gene model (gff3,gtf) file for splice junctions"*: `output` (Input dataset)
>            - *"Per gene/transcript output"*: `No per gene or transcript output`
>    - *"Use 2-pass mapping for more sensitive novel splice junction discovery"*: `Yes, perform single-sample 2-pass mapping of all reads`
>    - *"Report chimeric alignments?"*: `Within the BAM output (together with regular alignments; WithinBAM SoftClip) soft-clipping in the CIGAR for supplemental chimeric alignments`
>    - In *"Output filter criteria"*:
>        - *"Would you like to set additional output filters?"*: `No`
>    - In *"Algorithmic settings"*:
>        - *"Configure seed, alignment and limits options"*: `Use parameters suggested for STAR-Fusion`
>    - *"Compute coverage"*: `No coverage`
>
>
{: .hands_on}


> <question-title></question-title>
>
> 1. What is RNA STAR, and what does it do?
> 2. How do I interpret the alignment statistics in STAR’s output?
>
> > <solution-title></solution-title>
> >
> > 1. STAR is a tool for aligning RNA-Seq reads to a reference genome, helping researchers understand gene expression and identify splice junctions. STAR requires RNA-Seq reads, usually in FASTQ format. It also needs a reference genome file in FASTA format and annotation files in GTF/GFF format to build an index. STAR outputs alignments in BAM/SAM format, as well as splice junction files. It can also provide additional alignment stats in log files.
> > 2. STAR provides logs with mapping statistics, such as the percentage of uniquely mapped reads, which can be useful for quality control. Aligned BAM files from STAR can be visualized in genome browsers like IGV (Integrative Genomics Viewer) to examine coverage and splicing.
> >
> {: .solution}
>
{: .question}

## Fusion detection with **Arriba**

**Arriba** is a specialized tool used for detecting gene fusions from RNA sequencing (RNA-seq) data. It is particularly focused on identifying fusion events in cancer, where gene fusions can drive oncogenic processes. Arriba uses the output from **RNA STAR** alignments, specifically looking at chimeric alignments that result from fusion transcripts, and applies a series of filtering steps to reduce false positives. 

Arriba’s pipeline includes features for:
- Filtering out common artifacts and false-positive fusions based on blacklisted regions.
- Annotating fusion breakpoints.
- Generating a visualization of detected fusion events.

The output includes a list of fusion candidates with key information like fusion partners, breakpoint locations, reading frames, and peptide sequences. Arriba’s results can provide insight into potential neoantigens, helping guide research into therapeutic targets or immune-based therapies for cancer.

> <hands-on-title> Fusion detection </hands-on-title>
>
> 1. {% tool [Arriba](toolshed.g2.bx.psu.edu/repos/iuc/arriba/arriba/2.4.0+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"STAR Aligned.out.sam"*: `mapped_reads` (output of **RNA STAR** {% icon tool %})
>    - *"Genome assembly fasta (that was used for STAR alignment)"*: `From your history`
>        - {% icon param-file %} *"Genome assembly fasta"*: `output` (Input dataset)
>    - *"Genome GTF annotation source"*: `From your history`
>        - {% icon param-file %} *"Gene annotation in GTF format"*: `output` (Input dataset)
>    - {% icon param-file %} *"File containing blacklisted ranges."*: `blacklist` (output of **Arriba Get Filters** {% icon tool %})
>    - {% icon param-file %} *"File containing protein domains"*: `protein_domains` (output of **Arriba Get Filters** {% icon tool %})
>    - {% icon param-file %} *"File containing known fusions"*: `known_fusions` (output of **Arriba Get Filters** {% icon tool %})
>    - *"Use whole-genome sequencing data"*: `no`
>    - *"Generate visualization"*: `Yes`
>        - {% icon param-file %} *"Cytobands"*: `cytobands` (output of **Arriba Get Filters** {% icon tool %})
>
>
{: .hands_on}


> <question-title></question-title>
>
> 1. What is ARRIBA, and what does it do?
> 2. How can I ensure ARRIBA finds specific known fusions?
>
> > <solution-title></solution-title>
> >
> > 1. ARRIBA is a tool for detecting gene fusions in RNA-Seq data, especially helpful for identifying cancer-associated fusions and other structural variations. ARRIBA needs:A sorted BAM file with RNA-Seq reads aligned by STAR; STAR’s chimeric output (Chimeric.out.junction) to identify candidate fusion junctions; Reference annotation files, like a gene annotation GTF file and a blacklist file to filter false positives.
> > 2. Ensure that the STAR alignment and ARRIBA parameters are optimized for sensitivity. Adjusting settings for segment length and alignment quality in STAR can improve detection of specific known fusions.
> >
> {: .solution}
>
{: .question}

## Clean up data using **Text reformatting**

**Text Reformatting** is a step used in bioinformatics workflows to manipulate and clean up data for easier downstream processing. In fusion detection workflows, text reformatting is often used to parse and restructure output files, making the data consistent and accessible for subsequent analysis steps.

In this workflow, text reformatting involves:
- Extracting specific columns or fields from tabular outputs, such as gene names, breakpoint coordinates, or fusion peptide sequences.
- Formatting peptide sequences and related information into specific columns or concatenating fields for unique identifiers.
- Converting the data into a consistent format that downstream tools can interpret, such as converting tab-separated values into a structured layout for database input or analysis.

The reformatting step ensures that the processed data adheres to the requirements of other tools, enabling seamless integration across the workflow and supporting reliable, interpretable final results.

> <hands-on-title> Formating Arriba output</hands-on-title>
>
> 1. {% tool [Text reformatting](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_awk_tool/1.1.2) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `fusions_tsv` (output of **Arriba** {% icon tool %})
>    - *"AWK Program"*:
> ```
> (NR==1){
>    for (i=1;i<=NF;i++) {
>        if ($i ~ gene1) { 
>            gene1 = i;
>        }
>        if ($i == gene2) { 
>            gene2 = i;
>        }
>        if ($i == breakpoint1) { 
>            breakpoint1 = i;
>        }
>        if ($i == breakpoint2) { 
>            breakpoint2 = i;
>        }
>        if ($i == reading_frame) { 
>            reading_frame = i;
>        }
>        if ($i == peptide_sequence) { 
>            pscol = i;
>        }
>    }
> }
> (NR>1){
>    pseq = $pscol
>    if (pseq != .) {
>        bp = index(pseq,|);
>        pos = bp - 8; 
>        n=split(pseq,array,|);
>        pep = toupper(array[1] array[2])
>        sub([*],,pep)
>        g1 = $gene1;
>        g2 = $gene2;
>        sub([(,].*,,g1);
>        sub([(,].*,,g2);
>        id = g1 _ g2
>        brkpnts = $breakpoint1 _ $breakpoint2 
>        neopep = substr(pep,pos)
>        if ($reading_frame == in-frame) {
>            neopep = substr(pep,pos,16)
>        }
>        print(id \t (NR-1) \t brkpnts \t neopep);  
>    }
> }
> ```
>
>
>
{: .hands_on}



## Data refinement with **Query Tabular**

**Query Tabular** is a bioinformatics tool used to extract and manipulate specific data from tabular datasets in workflows. This tool allows users to perform SQL-like queries on tabular data, enabling them to filter, aggregate, and transform datasets based on user-defined criteria.

In this workflow, the **Query Tabular** tool is employed for several purposes:

- **Data Filtering:** Users can select specific rows based on certain conditions (e.g., filtering fusions that meet particular criteria).
- **Column Manipulation:** Users can specify which columns to retain or create new columns by combining or transforming existing data.
- **Aggregation:** The tool allows for summarizing data, such as counting occurrences of specific fusion events or summarizing results based on particular categories.
- **Output Customization:** Users can format the output to suit downstream processing needs, making it easier to pass data to subsequent analysis tools.

By leveraging **Query Tabular**, researchers can efficiently refine and structure their data, ensuring that only relevant information is carried forward in the workflow, ultimately aiding in the identification and analysis of significant biological insights.

> <hands-on-title> Manipulating the data to extract fusions </hands-on-title>
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.3.1) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `outfile` (output of **Text reformatting** {% icon tool %})
>            - In *"Table Options"*:
>                - *"Specify Column Names (comma-separated list)"*: `c1,c2,c3,c4`
>    - *"SQL Query to generate tabular output"*:
> ``` sql
> SELECT t1.c1 || '__' || t1.c2  || '__' || t1.c3, t1.c4
>FROM t1
> ```
>    - *"include query result column headers"*: `No`
>
>
{: .hands_on}


##  **Tabular-to-FASTA**

Tabular to FASTA conversion is a common task in bioinformatics that transforms data structured in a tabular format (such as CSV or TSV) into FASTA format, widely used for representing nucleotide or protein sequences. This conversion is essential when sequence data needs to be input into various bioinformatics tools or databases that require FASTA-formatted files.

> <hands-on-title> Converting tabular to fasta </hands-on-title>
>
> 1. {% tool [Tabular-to-FASTA](toolshed.g2.bx.psu.edu/repos/devteam/tabular_to_fasta/tab2fasta/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Tab-delimited file"*: `output` (output of **Query Tabular** {% icon tool %})
>    - *"Title column(s)"*: `c['1']`
>    - *"Sequence column"*: `c2`
>
>
{: .hands_on}


## Using **Regex Find And Replace**

Using regex (regular expressions) for find and replace is a powerful technique for text manipulation, allowing you to search for patterns and replace them with desired text. Below is a guide on how to use regex for find and replace, including examples in different programming languages. In this context, we are adding "fusion" to the database header.

> <hands-on-title> Adding fusion tag in the fasta header </hands-on-title>
>
> 1. {% tool [Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.3) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `output` (output of **Tabular-to-FASTA** {% icon tool %})
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `>(\b\w+\S+)(.*$)`
>            - *"Replacement"*: `>generic|fusion_\1|\2`
>
>
{: .hands_on}


# Conclusion

The workflow outlined above demonstrates a systematic approach to processing biological data, emphasizing the importance of each step in ensuring accurate and reliable results. By integrating tools like RNA-STAR for alignment and Arriba for structural variant detection, researchers can effectively analyze complex genomic information. The transition from tabular data to FASTA format and the application of regex for find-and-replace operations further streamline data management, enhancing efficiency and clarity. Ultimately, this workflow not only facilitates the identification of neoantigens but also contributes to the broader goals of personalized medicine and targeted therapies. By leveraging these methodologies, researchers can gain deeper insights into the genetic underpinnings of diseases and advance the development of innovative treatments.

# Disclaimer

Please note that all the software tools used in this workflow are subject to version updates and changes. As a result, the parameters, functionalities, and outcomes may differ with each new version. Additionally, if the protein sequences are downloaded at different times, the number of sequences may also vary due to updates in the reference databases or tool modifications. We recommend the users to verify the specific versions of software tools used to ensure the reproducibility and accuracy of results.
