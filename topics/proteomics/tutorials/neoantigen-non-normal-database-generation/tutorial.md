---
layout: tutorial_hands_on

title: "Neoantigen 2: Non-normal-Database-Generation"
zenodo_link: ""
questions:
- Why do we need to generate a customized fusion database for Proteogenomics research?
objectives:
- Downloading databases related to 16SrRNA data
- For better neoantigen identification results.
time_estimation: 3H
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
            - neoantigen-non-normal-database
tags: [label-free]
redirect_from:
- /topics/proteomics/tutorials/neoantigen-1-fusion-database-generation/tutorial

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


## Sub-step with **Convert compressed file to uncompressed.**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Convert compressed file to uncompressed.](CONVERTER_gz_to_uncompressed) %} with the following parameters:
>    - {% icon param-file %} *"Choose compressed file"*: `output` (Input dataset)
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

## Sub-step with **Convert compressed file to uncompressed.**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Convert compressed file to uncompressed.](CONVERTER_gz_to_uncompressed) %} with the following parameters:
>    - {% icon param-file %} *"Choose compressed file"*: `output` (Input dataset)
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

## Sub-step with **HISAT2**

> <hands-on-title> Task description </hands-on-title>
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

## Sub-step with **FreeBayes**

> <hands-on-title> Task description </hands-on-title>
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

## Sub-step with **CustomProDB**

> <hands-on-title> Task description </hands-on-title>
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

## Sub-step with **FASTA-to-Tabular**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [FASTA-to-Tabular](toolshed.g2.bx.psu.edu/repos/devteam/fasta_to_tabular/fasta2tab/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Convert these sequences"*: `output_indel` (output of **CustomProDB** {% icon tool %})
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

## Sub-step with **FASTA-to-Tabular**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [FASTA-to-Tabular](toolshed.g2.bx.psu.edu/repos/devteam/fasta_to_tabular/fasta2tab/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Convert these sequences"*: `output_snv` (output of **CustomProDB** {% icon tool %})
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

## Sub-step with **FASTA-to-Tabular**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [FASTA-to-Tabular](toolshed.g2.bx.psu.edu/repos/devteam/fasta_to_tabular/fasta2tab/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Convert these sequences"*: `output_rpkm` (output of **CustomProDB** {% icon tool %})
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

## Sub-step with **SQLite to tabular**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SQLite to tabular](toolshed.g2.bx.psu.edu/repos/iuc/sqlite_to_tabular/sqlite_to_tabular/2.0.0) %} with the following parameters:
>    - {% icon param-file %} *"SQLite Database"*: `output_genomic_mapping_sqlite` (output of **CustomProDB** {% icon tool %})
>    - *"SQL query"*: `SELECT pro_name, chr_name, cds_chr_start - 1, cds_chr_end,strand,cds_start - 1, cds_end
  FROM genomic_mapping
  ORDER BY pro_name, cds_start, cds_end`
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

## Sub-step with **FASTA-to-Tabular**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [FASTA-to-Tabular](toolshed.g2.bx.psu.edu/repos/devteam/fasta_to_tabular/fasta2tab/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Convert these sequences"*: `output_rpkm` (output of **CustomProDB** {% icon tool %})
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

## Sub-step with **SQLite to tabular**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SQLite to tabular](toolshed.g2.bx.psu.edu/repos/iuc/sqlite_to_tabular/sqlite_to_tabular/2.0.0) %} with the following parameters:
>    - {% icon param-file %} *"SQLite Database"*: `output_variant_annotation_sqlite` (output of **CustomProDB** {% icon tool %})
>    - *"SQL query"*: `SELECT var_pro_name,pro_name,cigar,annotation
FROM variant_annotation`
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

## Sub-step with **Column Regex Find And Replace**

> <hands-on-title> Task description </hands-on-title>
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

## Sub-step with **Column Regex Find And Replace**

> <hands-on-title> Task description </hands-on-title>
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

## Sub-step with **Column Regex Find And Replace**

> <hands-on-title> Task description </hands-on-title>
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

## Sub-step with **Column Regex Find And Replace**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Column Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regexColumn1/1.0.3) %} with the following parameters:
>    - {% icon param-file %} *"Select cells from"*: `query_results` (output of **SQLite to tabular** {% icon tool %})
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

## Sub-step with **Column Regex Find And Replace**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Column Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regexColumn1/1.0.3) %} with the following parameters:
>    - {% icon param-file %} *"Select cells from"*: `query_results` (output of **SQLite to tabular** {% icon tool %})
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

## Sub-step with **Tabular-to-FASTA**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Tabular-to-FASTA](toolshed.g2.bx.psu.edu/repos/devteam/tabular_to_fasta/tab2fasta/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Tab-delimited file"*: `out_file1` (output of **Column Regex Find And Replace** {% icon tool %})
>    - *"Title column(s)"*: `c['1']`
>    - *"Sequence column"*: `c2`
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

## Sub-step with **Tabular-to-FASTA**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Tabular-to-FASTA](toolshed.g2.bx.psu.edu/repos/devteam/tabular_to_fasta/tab2fasta/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Tab-delimited file"*: `out_file1` (output of **Column Regex Find And Replace** {% icon tool %})
>    - *"Title column(s)"*: `c['1']`
>    - *"Sequence column"*: `c2`
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

## Sub-step with **Tabular-to-FASTA**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Tabular-to-FASTA](toolshed.g2.bx.psu.edu/repos/devteam/tabular_to_fasta/tab2fasta/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Tab-delimited file"*: `out_file1` (output of **Column Regex Find And Replace** {% icon tool %})
>    - *"Title column(s)"*: `c['1']`
>    - *"Sequence column"*: `c2`
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

## Sub-step with **Query Tabular**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.3.2) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `out_file1` (output of **Column Regex Find And Replace** {% icon tool %})
>            - In *"Table Options"*:
>                - *"Specify Name for Table"*: `variant_annotation`
>                - *"Specify Column Names (comma-separated list)"*: `name,reference,cigar,annotation`
>                - In *"Table Index"*:
>                    - {% icon param-repeat %} *"Insert Table Index"*
>                        - *"Index on Columns"*: `name,cigar`
>    - *"Save the sqlite database in your history"*: `Yes`
>    - *"include query result column headers"*: `Yes`
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
>    - {% icon param-files %} *"Datasets to concatenate"*: `out_file1` (output of **Column Regex Find And Replace** {% icon tool %}), `output` (output of **bed to protein map** {% icon tool %})
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
>                - {% icon param-file %} *"FASTA File"*: `output` (output of **Tabular-to-FASTA** {% icon tool %})
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
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `out_file1` (output of **Concatenate datasets** {% icon tool %})
>            - In *"Table Options"*:
>                - *"Specify Name for Table"*: `feature_cds_map`
>                - *"Specify Column Names (comma-separated list)"*: `name,chrom,start,end,strand,cds_start,cds_end`
>                - In *"Table Index"*:
>                    - {% icon param-repeat %} *"Insert Table Index"*
>                        - *"Index on Columns"*: `name,cds_start,cds_end`
>    - *"Save the sqlite database in your history"*: `Yes`
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
