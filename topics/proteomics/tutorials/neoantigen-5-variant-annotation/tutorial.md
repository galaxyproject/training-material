---
layout: tutorial_hands_on

title: "Neoantigen 5: Variant Annotation"
zenodo_link: ''
questions:
- How can neoantigens be identified in cancer genomes?
- What role do neoantigens play in personalized immunotherapy?
- How can mutations in cancer cells be used to predict immune system responses?
objectives:
- Identify potential neoantigens from sequencing data.
- Annotate somatic mutations and predict peptide sequences.
- Predict MHC binding affinities for neoantigens.
- Interpret data using bioinformatics tools for cancer immunotherapy applications.
time_estimation: 3H
key_points:
- Neoantigen discovery involves the identification of somatic mutations that generate novel peptide sequences recognized by the immune system.
- Accurate prediction of MHC binding affinity is critical for evaluating the immunogenic potential of neoantigens.
- Immunoinformatics tools and pipelines are integral for analyzing neoantigens from genomic data and visualizing results.
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
            - neoantigen-6-predicting-hla-binding
tags: [label-free]
redirect_from:
- proteomics/tutorials/neoantigen-5-variant-annotation/tutorial

---


Neoantigens are tumor-specific antigens that arise from somatic mutations in cancer cells. These mutations result in the generation of abnormal peptides that can be presented by the immune system’s Major Histocompatibility Complex (MHC) molecules. Neoantigens are gaining significant attention in cancer immunotherapy, as they hold the potential to be used in personalized vaccines and therapies aimed at stimulating the immune system to target and destroy tumor cells.

The process of identifying neoantigens begins with the analysis of genomic data to detect somatic mutations, followed by the prediction of the corresponding peptide sequences that can bind to MHC molecules. This process involves a combination of bioinformatics tools and techniques, including variant calling, peptide prediction, and MHC binding affinity analysis. Successful identification of neoantigens can pave the way for personalized treatment strategies that improve patient outcomes in cancer therapy.

This tutorial focuses on the Neoantigen Annotation pipeline, which is designed to assist researchers and clinicians in annotating and predicting neoantigens from high-throughput genomic data. By the end of this tutorial, you will be equipped with the skills to extract meaningful insights from genomic datasets, predict immunogenic peptides, and understand how these peptides interact with the immune system. These insights are essential for advancing the field of personalized immunotherapy and improving cancer treatment outcomes.

![Variant-annotation-overview-workflow]({% link topics/proteomics/images/neoantigen/PepPointer_Characterization_1.PNG %})

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Neoantigen Variant Annotation

This tutorial introduces the process of annotating neoantigens using the PepPointer tool. The aim is to identify potential neoantigens by predicting peptide sequences derived from somatic mutations in cancer genomes, followed by their binding affinity to MHC molecules. The workflow integrates several key bioinformatics steps, including variant analysis, peptide prediction, and immune system interaction analysis. Below is an overview of each major step involved in the tutorial:It is divided into two parts: A, which is the variant annotation, and B, which focuses on the database generation for IEDB.
# A
### 1. Get Data
The first step in the process is to gather the necessary data for analysis, which typically includes genomic data in VCF (Variant Call Format) or other formats containing information about somatic mutations. The data is then uploaded into the analysis environment. A well-organized dataset is crucial for smooth and efficient analysis, so it’s important to ensure proper file structure and metadata tagging before beginning.

### 2. Mutation to Peptide Mapping
The next step is to map these mutations to peptide sequences. This involves creating a list of peptides that contain the mutated residue, representing potential neoantigens. Peptide prediction tools are used to generate the corresponding peptide sequences from the mutated positions.

### 3. Annotation and Filtering
At this stage, the predicted peptides are annotated with relevant biological and immunological information, such as their predicted MHC class, binding affinity, and potential for being recognized by T-cells. Filtering is performed to retain only the most promising candidates, based on binding affinity thresholds and relevance to the tumor type being studied.

### 4. Mapping Peptide sequences with **PepPointer**
PepPointer is used to map the peptide sequences to their corresponding genomic coordinates. This tool helps align peptide sequences (often derived from proteomic data) to the genomic context, providing useful insights into where these peptides are located in the genome. It allows researchers to determine which genomic regions are associated with the peptides of interest, facilitating the study of their potential functional roles.

### 5. Visualization and Interpretation
The final step involves visualizing the results of the annotation and filtering steps. Various bioinformatics tools can be used to present the data in a way that is easy to interpret, such as visualizing peptide binding affinity scores or generating summary plots that highlight the most immunogenic neoantigens. This step helps in drawing meaningful conclusions about the potential of the identified peptides for cancer immunotherapy.

# B 
### 6. Generating FASTA for MHC binding tool
In this step, we prepare the peptide sequences in FASTA format to be used with an MHC binding prediction tool. MHC (Major Histocompatibility Complex) binding tools are often used in immunology research to predict which peptides can bind to specific MHC molecules and present them to T-cells. 

# A: Variant annotation
![Variant-annotation-overview-workflow]({% link topics/proteomics/images/neoantigen/PepPointer_Characterization_2.PNG %})
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


## Mutation to Peptide Mapping with **Query Tabular**

In this step, we will use the Query Tabular tool to perform a SQL query on a tabular dataset. This tool allows you to filter and extract data from a dataset by applying SQL commands, making it ideal for working with structured data such as CSV or TSV files. You will use it to join tables, select relevant columns, and generate a new output based on the criteria specified in the query. This step helps in refining and organizing the data for further analysis.
This step extracts information about novel peptides from Frapipe, which primarily includes protein details and the location of the peptide on the protein.

> <hands-on-title> Query Tabular </hands-on-title>
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.3.2) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `output` (Input dataset)
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `output` (Input dataset)
>    - *"SQL Query to generate tabular output"*:
> ``` sql
> SELECT t1.c1,t2.c13,t2.c5,t2.c6
> FROM t1
> INNER JOIN t2
> ON t1.c1 = t2.c1
> ```
>    - *"include query result column headers"*: `No`
>
{: .hands_on}



## Annotation and Filtering

### Converting delimiters 
In this step, we will use the Convert tool to modify characters in a dataset. Specifically, the tool will be used to convert all instances of pipe characters (|) to another format. The tool helps clean and standardize the data for subsequent processing or analysis. This is often a necessary step when preparing data for tools that require specific formats or when performing tasks such as parsing or file importation.

> <hands-on-title> Convert </hands-on-title>
>
> 1. {% tool [Convert](Convert characters1) %} with the following parameters:
>    - *"Convert all"*: `Pipes`
>    - {% icon param-file %} *"in Dataset"*: `output` (output of **Query Tabular** {% icon tool %})
>
>
{: .hands_on}


> <question-title></question-title>
>
> 1. What is the purpose of using the "Convert" tool in this step?
>
> > <solution-title></solution-title>
> >
> > 1. The "Convert" tool is used to replace or remove specific characters in a dataset, such as pipe characters (|). This ensures that the data conforms to the required format for subsequent analysis or tools. Pipe characters are often used as delimiters in tabular data. Replacing them might be necessary if the data needs to be formatted for a different tool or if the pipe character interferes with parsing in downstream steps.
> >
> {: .solution}
> 
{: .question}

### Editing certain "u" characters 
In this step, we will use the Column Regex Find And Replace tool to find and replace specific patterns in the data columns. Regular expressions (regex) allow for pattern matching, making it easier to identify and modify data based on specific patterns. In this case, the tool will search for occurrences of the pattern u_ and replace them with u: in column c3 of the dataset. This is commonly done to standardize data formatting or clean data for further analysis. In this step, any protein sequences with a "u_" will be replaced with "u".

> <hands-on-title> Column Regex Find And Replace </hands-on-title>
>
> 1. {% tool [Column Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regexColumn1/1.0.3) %} with the following parameters:
>    - {% icon param-file %} *"Select cells from"*: `out_file1` (output of **Convert** {% icon tool %})
>    - *"using column"*: `c3`
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `u_`
>            - *"Replacement"*: `u:`
>
>    
>
{: .hands_on}


### Removing Colons and converting into a different tabular 
In this step, we will use the Convert tool to remove colons from the dataset. This conversion is useful when you need to standardize the data or prepare it for downstream processing, where colons may interfere with the analysis or formatting. By selecting the "Convert all" option for colons, all occurrences of colons in the dataset will be removed, creating a clean tabular format for further analysis.

> <hands-on-title> Convert </hands-on-title>
>
> 1. {% tool [Convert](Convert characters1) %} with the following parameters:
>    - *"Convert all"*: `Colons`
>    - {% icon param-file %} *"in Dataset"*: `out_file1` (output of **Column Regex Find And Replace** {% icon tool %})
>
>
{: .hands_on}


### Extracting bed file information **Query Tabular**
In this step, we will use the Query Tabular tool to extract specific information from a dataset, such as a BED file containing genomic regions, and match it with novel peptides. This allows for identifying the relevant genomic and peptide information by querying data from two sources and combining them through an SQL query. By using an INNER JOIN operation, we can merge data from two tables based on shared columns, and retrieve the necessary information. This query extracts specific columns from both the BED file (such as genomic coordinates) and the novel peptide dataset (such as peptide sequences or identifiers), enabling the identification of peptides that correspond to specific genomic regions. These are the columns that will be extracted - 
- Chrom: Chromosome name (e.g., chr1).
- Start: Starting position of the feature (zero-based index).
- End: Ending position of the feature (one-based index).
- Name: Identifier for the feature (e.g., gene name).
- Score: Numeric score representing feature strength or relevance.
- Strand: Strand orientation (+ or -).
- ThickStart and ThickEnd: Define the start and end of the transcribed or relevant part of the feature, often the coding region.

> <hands-on-title> Query Tabular </hands-on-title>
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.3.2) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `out_file1` (output of **Convert** {% icon tool %})
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `output` (Input dataset)
>    - *"SQL Query to generate tabular output"*:
> ``` sql
> SELECT t1.c1,t1.c8,t1.c9,t1.c11,t1.c12,t2.c5,t2.c6,t1.c7
> FROM t1
> INNER JOIN t2
> ON t1.c3 = t2.c4
> ```
>    - *"include query result column headers"*: `No`
>
>
{: .hands_on}


### Performing calculations to convert proteomic coordinates to genomic coordinates.  
To convert proteomic coordinates to genomic coordinates, it is essential to account for the relationship between the protein sequence and its corresponding gene or genomic region. In this workflow, the proteomic coordinates have already been extracted at the amino acid level. Since each amino acid in the protein sequence corresponds to a triplet of nucleotides (a codon) in the mRNA, we need to multiply the proteomic coordinate by 3 to obtain the genomic coordinate. This conversion will give us the position of each amino acid within the genomic sequence. The resulting genomic coordinates are stored in a separate column for easy reference. Once this step is completed, we can extract and organize the information in the correct order for further analysis or mapping to the genomic reference.


> <hands-on-title> Query Tabular </hands-on-title>
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.3.2) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `output` (output of **Query Tabular** {% icon tool %})
>    - *"SQL Query to generate tabular output"*:
> ``` sql
> SELECT t1.*, t1.c4 * 3 AS c2_multiplied, t1.c5 * 3 AS c3_multiplied FROM t1
> ```
>    - *"include query result column headers"*: `No`
>
>
{: .hands_on}


### Annotating the genomic coordinate

The Query Tabular step in this workflow is used to extract and calculate genomic coordinates based on the proteomic data. The SQL query within the tool defines two calculations for genomic coordinates, start and stop, based on the strand information of the data. For each row in the input dataset (t1), if the strand (t1.c7) is "-" (negative), the genomic coordinates are calculated by subtracting the position from the given end (t1.c3 - t1.c9 for start, and t1.c3 - t1.c10 for stop). If the strand is "+" (positive), the genomic coordinates are calculated by adding the respective positions (t1.c2 + t1.c9 for start, and t1.c2 + t1.c10 for stop). These calculated coordinates are then returned in the query results, where they will be included as new columns (start and stop). This step is essential for transforming the proteomic information into genomic positions for further analysis.

> <hands-on-title> Query Tabular </hands-on-title>
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.3.2) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `output` (output of **Query Tabular** {% icon tool %})
>    - *"SQL Query to generate tabular output"*:
> ``` sql
> SELECT t1.*,
> CASE
> WHEN t1.c7 = '-' THEN t1.c3 - t1.c9
> WHEN t1.c7 = '+' THEN t1.c2 + t1.c9
> END AS start,
> CASE
> WHEN t1.c7 = '-' THEN t1.c3 - t1.c10
> WHEN t1.c7 = '+' THEN t1.c2 + t1.c10
> END AS stop FROM t1
> ```
>    - *"include query result column headers"*: `No`
>
>
{: .hands_on}


> <question-title></question-title>
>
> 1. Why do I need to differentiate between the positive and negative strands when calculating genomic coordinates?
> 2. What is the significance of the t1.c9 and t1.c10 columns in the SQL query?
>
> > <solution-title></solution-title>
> >
> > 1. The positive and negative strands represent the two directions in which DNA is read, and they affect how genomic coordinates are calculated. For the negative strand, the coordinates are subtracted, whereas for the positive strand, they are added. This ensures the correct mapping of the protein sequence to the genome.
> > 2. The t1.c9 and t1.c10 columns represent specific offsets or lengths within the dataset, which are used to adjust the calculated genomic start and stop coordinates. These offsets could correspond to peptide lengths or specific features of the protein that need to be accounted for in the conversion to genomic coordinates.
> >
> {: .solution}
>
{: .question}

### Generating BED file for Peppointer 
This step is necessary to extract and reorganize relevant genomic information from the dataset. By querying specific columns such as chromosome (chromosome), start (chromStart), end (chromEnd), and strand (strand), we are preparing the data for further analysis. These values are essential for mapping proteomic or peptide data to the genomic coordinates, ensuring accurate alignment and interpretation of the sequence in the context of its genomic location. Additionally, renaming columns enhances clarity and standardizes the format, making it easier to work with the data in subsequent steps.

> <hands-on-title> Query Tabular </hands-on-title>
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.3.2) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `output` (output of **Query Tabular** {% icon tool %})
>    - *"SQL Query to generate tabular output"*:
> ``` sql
> SELECT
> c8 AS `chromosome`,
> c11  AS `chromStart`,
> c12 AS `chromEnd`,
> c1 AS `name`,
> c6 AS `score`,
> c7 AS `strand`
> FROM t1
> ```
>    - *"include query result column headers"*: `No`
>
>
{: .hands_on}


> <question-title></question-title>
>
> 1. What is the purpose of renaming columns like chromosome, chromStart, and chromEnd in the SQL query?
> 2. Why is the strand column important in this query, and how does it affect genomic coordinate interpretation?
>
> > <solution-title></solution-title>
> >
> > 1. Renaming the columns helps in providing more meaningful names for the data, making it easier to understand the genomic features such as the chromosome and start/end coordinates. These renamed columns align with common genomic nomenclature, improving readability and analysis.
> > 2. The strand column indicates whether the sequence is on the positive or negative strand of the DNA. This is important for correctly interpreting the orientation of the sequence and calculating its genomic coordinates. A positive strand means the coordinates will be calculated from the start, while a negative strand requires reverse calculations to adjust the coordinates properly.
> >
> {: .solution}
>
{: .question}

## Mapping Peptide sequences with **PepPointer**

PepPointer is a tool designed to map peptide sequences to their respective genomic locations using data such as GTF and BED files. In this workflow, PepPointer takes the GTF file, which contains gene annotations and genomic coordinates, and combines it with a BED file, which provides chromosomal coordinates of the peptides. By doing this, PepPointer can identify the exact genomic locations of the peptides, linking them to genes and exons. This step is crucial for accurately correlating proteomic data to genomic sequences, enabling a better understanding of the genomic context in which the peptides are found.

> <hands-on-title> PepPointer </hands-on-title>
>
> 1. {% tool [PepPointer](toolshed.g2.bx.psu.edu/repos/galaxyp/pep_pointer/pep_pointer/0.1.3+galaxy1) %} with the following parameters:
>    - *"Choose the source of the GTF file"*: `From history`
>        - {% icon param-file %} *"GTF file with the genome of interest"*: `output` (Input dataset)
>    - {% icon param-file %} *"BED file with chromosomal coordinates of peptide"*: `output` (output of **Query Tabular** {% icon tool %})
>
>
{: .hands_on}


> <question-title></question-title>
>
> 1. How does the GTF file help in identifying the genomic location of peptides?
> 2. What happens if the coordinates in the BED file don't align with the genomic features in the GTF file?
>
> > <solution-title></solution-title>
> >
> > 1. The GTF file contains detailed annotations of the genome, including information about genes, exons, and other genomic features. By providing the chromosomal coordinates of these features, the GTF file allows PepPointer to match peptide sequences to the corresponding genes and exons. This ensures that the peptides are mapped accurately within the genome, linking them to the correct genomic context.
> > 2. If the coordinates in the BED file don't align with the genomic features in the GTF file, PepPointer may not be able to accurately map the peptides to their corresponding genes or exons. This misalignment could result in incorrect annotations or a failure to identify the proper genomic location of the peptides. It is important to ensure that both files are from the same genome assembly to avoid such discrepancies.
> >
> {: .solution}
>
{: .question}

## Visualization and Interpretation

In this step, we are using Query Tabular to extract and format relevant information from the results produced by PepPointer. The SQL query in this tool is designed to structure the data into a more readable format, providing key details such as the peptide ID, chromosome location, start and end positions, strand orientation, and any annotations. Additionally, it generates genome coordinates in a format suitable for viewing in genome browsers like IGV and UCSC Genome Browser. This formatting step helps to visualize and interpret the data more effectively by linking the peptides to their genomic context.

> <hands-on-title> Query Tabular </hands-on-title>
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.3.2) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `classified` (output of **PepPointer** {% icon tool %})
>    - *"SQL Query to generate tabular output"*:
> ``` sql
> SELECT 
> c4 AS Peptide,
> c1 AS Chromosome,
> c2 AS Start,
> c3 AS End,
> c6 AS Strand,
> c7 AS Annotation,
> c1||':'||c2||'-'||c3 AS IGV_Genome_Coordinate,
> 'https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position='||c1||'%3A'||c2||'-'||c3 AS UCSC_Genome_Browser
> FROM  t1
> ```
>    - *"include query result column headers"*: `Yes`
>
>
{: .hands_on}


> <question-title></question-title>
>
> 1. Why is it important to include both the "IGV_Genome_Coordinate" and "UCSC_Genome_Browser" columns?
> 2. What is the significance of the "Annotation" column, and what type of data might it contain?
>
> > <solution-title></solution-title>
> >
> > 1. The IGV_Genome_Coordinate column provides a formatted coordinate that is directly usable in genome viewers like IGV, helping users visualize the peptide's genomic location. The UCSC_Genome_Browser column provides a clickable URL to launch UCSC Genome Browser with the coordinates, offering another convenient method for visualizing the peptide's location in a widely used genome browser. Including both makes it easier for users to explore the data in different genomic tools.
> > 2. The Annotation column provides additional context about the genomic feature associated with the peptide, such as whether the peptide is part of a gene, exon, or other regulatory element. This column can help to interpret the biological relevance of the peptide in the genomic region and is important for understanding the functional role of the peptide within the genome. The data might include annotations like "gene", "exon", "intron", or other functional genomic elements.
> >
> {: .solution}
>
{: .question}


# B: Database for IEDB 
![Variant-annotation-overview-workflow]({% link topics/proteomics/images/neoantigen/PepPointer_Characterization_3.PNG %})

## Generating FASTA for MHC binding tool
This output is an input for the next workflow.

> <hands-on-title> Tabular-to-FASTA </hands-on-title>
>
> 1. {% tool [Tabular-to-FASTA](toolshed.g2.bx.psu.edu/repos/devteam/tabular_to_fasta/tab2fasta/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Tab-delimited file"*: `output` (output of **Query Tabular** {% icon tool %})
>    - *"Title column(s)"*: `c['2']`
>    - *"Sequence column"*: `c1`
>
>
{: .hands_on}


# Conclusion

In this tutorial, we have covered the process of annotating neoantigen peptides using the PepPointer tool in a bioinformatics pipeline. The main focus was to guide you through the various steps involved in the analysis, including querying tabular datasets, converting data formats, and using the PepPointer tool to classify peptides based on genomic information.

By the end of this tutorial, you should now be able to:
- Understand the steps involved in preparing and processing peptide data for neoantigen prediction.
- Use tools like Query Tabular, Convert, and PepPointer to manipulate and classify peptide sequences.
- Recognize the importance of annotating neoantigens for immunotherapy research and how bioinformatics tools can facilitate this process.

The key takeaway from this tutorial is the ability to manipulate complex genomic data and generate insights into potential therapeutic targets using bioinformatics techniques. You should now be better equipped to handle similar datasets and apply these techniques to other types of genomic and proteomic research.

# Disclaimer 

Please note that all the software tools used in this workflow are subject to version updates and changes. As a result, the parameters, functionalities, and outcomes may differ with each new version. We recommend that users verify the specific versions of software tools used to ensure the reproducibility and accuracy of results.

