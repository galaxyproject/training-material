---
layout: tutorial_hands_on

title: "Clinical Metaproteomics 2: Discovery"
zenodo_link: "https://doi.org/10.5281/zenodo.10105821"
questions:
- How to perform database searching?
- How to extract microbial and Human protein and peptide sequences from the results
objectives:
- Perform Database searching using two algorithms
- Extract confident peptides and proteins
- Generate a microbial peptide panel for verification
time_estimation: 3H
key_points:
- Employ SearchGUI/PeptideShaker and MaxQuant for database searching
- Extraction of confident microbial peptides for verification
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
            - clinical-mp-3-verification
tags: [label-TMT11]
redirect_from:
- /topics/proteomics/tutorials/clinical-mp-discovery/tutorial

recordings:
- captioners:
  - katherine-d21
  date: '2024-06-21'
  galaxy_version: '23.1'
  length: 16M
  youtube_id: 5Pg9LLfFGX4
  speakers:
  - katherine-d21
---



This tutorial can be followed with any user-defined database but would work best if the clinical metaproteomics database generation module was used (see [Database Generation tutorial](https://github.com/subinamehta/training-material/blob/main/topics/proteomics/tutorials/clinical-mp-database-generation/tutorial.md)). The MetaNovo tool generates a more manageable database that contains identified proteins. The MetaNovo-generated database merged with Human SwissProt (reviewed only) and contaminants (cRAP) databases to generate a compact database (~21.2k protein sequences) that will be used for peptide identification.


# Peptide identification
The MSMS data will be searched against the compact database `Human UniProt Microbial Proteins (from MetaNovo) and cRAP` to identify peptide and protein sequences via sequence database searching. For this tutorial, two peptide identification programs will be used: SearchGUI/PeptideShaker and MaxQuant. However, you could use other software too, such as Fragpipe or Scribe. For the purpose of this tutorial, a dataset of the 4 RAW/MGF files will be used as the MS/MS input.


![Discovery Workflow]({% link topics/proteomics/images/clinical-mp/clinical-mp-discovery.JPG %})


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Database Searching
This step is to identify proteins based on mass spectrometry data. The algorithms identify peptides in the spectra and search a protein sequence database to match observed peptide data with theoretical peptide masses and spectra. Scoring and false discovery rate control help assess the reliability of matches, followed by protein inference to determine the proteins present in the sample. These algorithms are essential for interpreting mass spectrometry data, aiding in protein identification, quantification, and insights into biological processes and disease mechanisms in proteomics research.

## Get data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    https://zenodo.org/records/10105821/files/Human_UniProt_Microbial_Proteins_(from_MetaNovo)_and_cRAP.fasta
>    https://zenodo.org/records/10105821/files/PTRC_Skubitz_Plex2_F10_9Aug19_Rage_Rep-19-06-08.raw
>    https://zenodo.org/records/10105821/files/PTRC_Skubitz_Plex2_F11_9Aug19_Rage_Rep-19-06-08.raw
>    https://zenodo.org/records/10105821/files/PTRC_Skubitz_Plex2_F13_9Aug19_Rage_Rep-19-06-08.raw
>    https://zenodo.org/records/10105821/files/PTRC_Skubitz_Plex2_F15_9Aug19_Rage_Rep-19-06-08.raw
>    https://zenodo.org/records/10105821/files/Experimental-Design_Discovery_MaxQuant.tabular
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
> 5. Add to each database a tag corresponding to user.
> 6. Create a dataset collection of all the raw files and MGF files.
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

# Import Workflow


> <hands-on-title>Running the Workflow</hands-on-title>
>
> 1. **Import the workflow** into Galaxy:
>
>    {% snippet faqs/galaxy/workflows_run_trs.md path="topics/proteomics/tutorials/clinical-mp-2-discovery/workflows/WF2_Discovery-Workflow.ga" title="Discovery Workflow" %}
>
> 2. Run **Workflow** {% icon workflow %} using the following parameters:
>    - *"Send results to a new history"*: `No`
>    - {% icon param-file %} *" RAW files"*: `RAW dataset collection`
>    - {% icon param-file %} *" Human UniProt Microbial Proteins (from MetaNovo) and cRAP"*: `Human_UniProt_Microbial_Proteins_(from_MetaNovo)_and_cRAP.fasta`
>    - {% icon param-file %} *" Experimental Design Discovery MaxQuant"*: `Experimental-Design_Discovery_MaxQuant.tabular`
>
>    {% snippet faqs/galaxy/workflows_run.md %}
>
{: .hands_on}


# Peptide identification
Using the compact database generated by MetaNovo as the input database, we will match MS/MS data to peptide sequences via sequence database searching.

For this tutorial, two peptide identification programs will be used: SearchGUI/PeptideShaker and MaxQuant. For both programs, the created dataset of the four MS datasets  in the history will be used as the MS/MS input. The RAW MS/MS data files will be converted into mascot generic format (MGF) files as that is the standard format in which MS/MS searches are performed.

Peptides identified from each program will be verified with the PepQuery tool to generate a master list of confident verified microbial peptides.

## Appending decoy sequenced to FASTA database with **FastaCLI**

Using the FastaCLI tool, decoy sequences will be appended to the FASTA database. Decoy sequences are protein sequences are not expected to be present in samples. For more information on how to generate and append decoy sequences, see [GTN Protein FASTA Database Handling](https://training.galaxyproject.org/archive/2019-12-01/topics/proteomics/tutorials/database-handling/tutorial.html#creating-a-decoy-database).

> <hands-on-title> FastaCLI </hands-on-title>
>
> 1. {% tool [FastaCLI](toolshed.g2.bx.psu.edu/repos/galaxyp/peptideshaker/fasta_cli/4.0.41+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Protein Database"*: `output` (Input dataset)
>
>
{: .hands_on}

## Converting RAW files to MGF files with **msconvert**

The msconvert tool allows for the conversion of mass spectrometry data files between different formats, such as thermo.raw, mgf, or mzml.

> <hands-on-title> msconvert: RAW to MGF </hands-on-title>
>
> 1. {% tool [msconvert](toolshed.g2.bx.psu.edu/repos/galaxyp/msconvert/msconvert/3.0.20287.2) %} with the following parameters:
>    - {% icon param-collection %} *"Input unrefined MS data"*: `output` (Input dataset collection)
>    - *"Do you agree to the vendor licenses?"*: `Yes`
>    - *"Output Type"*: `mgf`
>    - In *"Data Processing Filters"*:
>        - *"Apply peak picking?"*: `Yes`
>        - *"(Re-)calculate charge states?"*: `no`
>
>
>
{: .hands_on}


> <question-title></question-title>
>
> 1. Why do we need to use MGF instead of RAW files for Search GUI?
>
> > <solution-title></solution-title>
> >
> > 1. SearchGUI is compatible only with MGF files, hence you have to use msconvert or Thermofile converter tools to convert the RAW format to MGF fomat.
> >
> {: .solution}
>
{: .question}

## Perform Database searching with **Search GUI**
SearchGUI is a database-searching tool that comprises different search engines to match sample MS/MS spectra to known peptide sequences. In our analysis, we will use X!Tandem and MS-GF+ as search algorithms within SearchGUI for matching spectra from mass spectrometry data against peptides from the protein sequence database.

The SearchGUI tool will perform a database search based on the parameters we've set and will generate a file (called a SearchGUI archive file) that will serve as the input for the PeptideShaker tool. The SearchGUI archive file contains Peptide-Spectral Matches (PSMs), and PeptideShaker is a post-processing software that will assess the confidence of the data. PeptideShaker also infers the identities of proteins based on the matched peptide sequences, and users are able to visualize these outputs to interpret results.  More information about database searching using SearchGUI and PeptideShaker is accessible at [Metaproteomics tutorial](https://gxy.io/GTN:T00221).

> <hands-on-title> Peptide discovery using SearchGUI </hands-on-title>
>
> 1. {% tool [Search GUI](toolshed.g2.bx.psu.edu/repos/galaxyp/peptideshaker/search_gui/4.0.41+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Identification Parameters file"*: `Identification_Parameters_File` (output of **Identification Parameters** {% icon tool %})
>    - {% icon param-file %} *"Fasta file"*: `input_database_concatenated_target_decoy` (output of **FastaCLI** {% icon tool %})
>    - {% icon param-file %} *"Input Peak Lists"*: `output` (output of **msconvert** {% icon tool %})
>    - *"SearchGUI Options"*: `Default`
>
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Why do we need to add decoy sequences to our FASTA database for Search GUI? And how many do we need to add?
>
> > <solution-title></solution-title>
> >
> > 1. Adding decoy sequences helps in FDR estimation, discriminating true positives from false positives, and quality control of the data. The number of decoy sequences you need to add to your database depends on the desired FDR level you want to achieve. A common practice is to use a 1:1 ratio of target sequences to decoy sequences. In other words, for every real protein sequence in your database, you would add a decoy sequence. This allows you to estimate the FDR at 1%, 5%, or any other chosen threshold.
> >
> {: .solution}
>
{: .question}
> <question-title></question-title>
>
> 1. What is the Identification Parameters tool?
>
> > <solution-title></solution-title>
> >
> > 1. Identification Parameters tool is an input required by the search GUI tool, it contains all the parameters required to run the search algorithms.
> >
> {: .solution}
>
{: .question}

## Post-processing of SearchGUI output using with **Peptide Shaker**

> <hands-on-title> Peptide Shaker </hands-on-title>
>
> 1. {% tool [Peptide Shaker](toolshed.g2.bx.psu.edu/repos/galaxyp/peptideshaker/peptide_shaker/2.0.33+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Compressed SearchGUI results"*: `searchgui_results` (output of **Search GUI** {% icon tool %})
>    - In *"Exporting options"*:
>        - *"Follow-up analysis export options"*: `Do not export`
>        - *"Identification features reports to be generated"*: `PSM Report` `Peptide Report` `Protein Report` `Certificate of Analysis`
>
>
{: .hands_on}


> <question-title></question-title>
>
> 1. What are the differences between the following reports from PeptideShaker: PSM report, Peptide report, and Protein report?
>
> > <solution-title></solution-title>
> >
> > 1. PSM reports focus on individual peptide-spectrum matches, providing detailed information about each spectrum and its assigned peptide sequence. Peptide reports summarize information about unique peptides and their properties. Protein reports, on the other hand, focus on proteins, including protein inference, grouping, and quantification, making them more suitable for understanding the overall protein composition in a sample. These reports serve different purposes in proteomic data analysis and are used to extract various levels of information from mass spectrometry results.
> >
> {: .solution}
>
{: .question}

##  Using Text Manipulation Tools to Manage Microbial Outputs from SearchGUI/PeptideShaker

> <hands-on-title> Selecting microbial peptides from SearchGUI/PeptideShaker with Select tool </hands-on-title>
>
> 1. {% tool [Select](Grep1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `output_peptides` (output of **Peptide Shaker** {% icon tool %})
>    - *"that"*: `NOT Matching`
>    - *"the pattern"*: `(_HUMAN)|(_REVERSED)|(REV_)|(CON)|(con)`
>    - *"Keep header line"*: `Yes`
>
>
{: .hands_on}


> <question-title></question-title>
>
> 1. What is the purpose of this step?
>
> > <solution-title></solution-title>
> >
> > 1. This step is to extract microbial peptides or to remove any peptides that match humans, reverse, contaminants, etc.
> >
> {: .solution}
>
{: .question}

> <hands-on-title> Selecting microbial PSMs from SearchGUI/PeptideShaker with Select </hands-on-title>
>
> 1. {% tool [Select](Grep1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `output_psm` (output of **Peptide Shaker** {% icon tool %})
>    - *"that"*: `NOT Matching`
>    - *"the pattern"*: `(_HUMAN)|(_REVERSED)|(REV_)|(CON)|(con)`
>    - *"Keep header line"*: `Yes`
>
>
{: .hands_on}


> <hands-on-title> Filtering confident microbial peptides from SGPS  with Filter </hands-on-title>
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `out_file1` (output of **Select** {% icon tool %})
>    - *"With following condition"*: `c17=='Confident'`
>    - *"Number of header lines to skip"*: `1`
>
>
>
{: .hands_on}


> <question-title></question-title>
>
> 1. In the Filtering steps, what does “Confidence” mean quantitatively, i.e. what is the percentage cutoff?
>
> > <solution-title></solution-title>
> >
> > 1. The term "Confidence" in the context of proteomic data analysis often refers to a measure of how reliable or trustworthy a particular protein or peptide identification is. However, the specific numerical value or percentage cutoff for confidence can vary depending on the software or approach you are using and the goals of your analysis. In many proteomics studies, researchers use a false discovery rate (FDR) to set a quantitative confidence threshold. Here we have set it as 1%FDR, which means that you're accepting only 1% or less of your reported identifications as likely to be false positives.
> >
> {: .solution}
>
{: .question}


> <hands-on-title> Filtering confident microbial PSMs from SGPS with Filter </hands-on-title>
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `out_file1` (output of **Select** {% icon tool %})
>    - *"With following condition"*: `c24=='Confident'`
>    - *"Number of header lines to skip"*: `1`
>
>
{: .hands_on}


We will generate and merge the Human SwissProt Protein Database and contaminants (cRAP) and convert the resulting FASTA file to a tabular file that will be used in the Query Tabular tool to generate distinct microbial peptides from SearchGUI/PeptideShaker.

> <hands-on-title> Merging Human SwissProt and cRAP databases for Query Tabular with FASTA Merge Files and Filter Unique Sequences </hands-on-title>
>
> 1. {% tool [FASTA Merge Files and Filter Unique Sequences](toolshed.g2.bx.psu.edu/repos/galaxyp/fasta_merge_files_and_filter_unique_sequences/fasta_merge_files_and_filter_unique_sequences/1.2.0) %} with the following parameters:
>    - *"Run in batch mode?"*: `Merge individual FASTAs (output collection if input is collection)`
>        - In *"Input FASTA File(s)"*:
>            - {% icon param-repeat %} *"Insert Input FASTA File(s)"*
>                - {% icon param-file %} *"FASTA File"*: `Human Swissprot Protein Database` (output of **Protein Database Downloader** {% icon tool %})
>                - - {% icon param-file %} *"FASTA File"*: `Contaminants cRAP database` (output of **Protein Database Downloader** {% icon tool %})
>
>
{: .hands_on}


> <hands-on-title> Converting FASTA sequences to TAB-delimited file with FASTA-to-Tabular </hands-on-title>
>
> 1. {% tool [FASTA-to-Tabular](toolshed.g2.bx.psu.edu/repos/devteam/fasta_to_tabular/fasta2tab/1.1.0) %} with the following parameters:
>    - {% icon param-file %} *"Convert these sequences"*: `output` (output of **FASTA Merge Files and Filter Unique Sequences** {% icon tool %})
>
>
{: .hands_on}


> <hands-on-title> Filtering out accession numbers from TAB-delimited file  with Filter Tabular </hands-on-title>
>
> 1. {% tool [Filter Tabular](toolshed.g2.bx.psu.edu/repos/iuc/filter_tabular/filter_tabular/3.3.0) %} with the following parameters:
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
>
{: .hands_on}


> <question-title></question-title>
>
> 1. What’s the difference between a FASTA and Tabular output?
>
> > <solution-title></solution-title>
> >
> > 1. FASTA Output: Typically used to report identified peptide or protein sequences, which are useful for building or updating sequence databases, for downstream sequence analysis, or for re-searching against the sequences.
> > Tabular Output: Used for presenting various information related to identified peptides or proteins, such as accession numbers, scores, abundance values, and other attributes. Tabular output facilitates data analysis, comparisons, and custom data processing.
> >
> >
> {: .solution}
>
{: .question}


> <hands-on-title> Querying protein accession numbers and peptide sequences of confident microbial PSMs (from SGPS) with Query Tabular </hands-on-title>
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.3.0) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `out_file1` (output of **Filter** {% icon tool %})
>            - In *"Filter Dataset Input"*:
>                - In *"Filter Tabular Input Lines"*:
>                    - {% icon param-repeat %} *"Insert Filter Tabular Input Lines"*
>                        - *"Filter By"*: `skip leading lines`
>                            - *"Skip lines"*: `1`
>                    - {% icon param-repeat %} *"Insert Filter Tabular Input Lines"*
>                        - *"Filter By"*: `prepend a line number column`
>            - In *"Table Options"*:
>                - *"Specify Name for Table"*: `psms`
>                - *"Specify Column Names (comma-separated list)"*: `ln,id,Proteins,Sequence`
>                - *"Only load the columns you have named into database"*: `Yes`
>                - In *"Table Index"*:
>                    - {% icon param-repeat %} *"Insert Table Index"*
>                        - *"Index on Columns"*: `ln`
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `out_file1` (output of **Filter** {% icon tool %})
>            - In *"Filter Dataset Input"*:
>                - In *"Filter Tabular Input Lines"*:
>                    - {% icon param-repeat %} *"Insert Filter Tabular Input Lines"*
>                        - *"Filter By"*: `skip leading lines`
>                            - *"Skip lines"*: `1`
>                    - {% icon param-repeat %} *"Insert Filter Tabular Input Lines"*
>                        - *"Filter By"*: `prepend a line number column`
>                    - {% icon param-repeat %} *"Insert Filter Tabular Input Lines"*
>                        - *"Filter By"*: `normalize list columns, replicates row for each item in list`
>                            - *"enter column numbers to normalize"*: `3`
>            - In *"Table Options"*:
>                - *"Specify Name for Table"*: `prots`
>                - *"Specify Column Names (comma-separated list)"*: `ln,id,prot`
>                - *"Only load the columns you have named into database"*: `Yes`
>                - In *"Table Index"*:
>                    - {% icon param-repeat %} *"Insert Table Index"*
>                        - *"This is a unique index"*: `Yes`
>                        - *"Index on Columns"*: `prot,ln`
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `output` (output of **Filter Tabular** {% icon tool %})
>            - In *"Table Options"*:
>                - *"Specify Name for Table"*: `Uniprot`
>                - *"Specify Column Names (comma-separated list)"*: `prot`
>                - In *"Table Index"*:
>                    - {% icon param-repeat %} *"Insert Table Index"*
>                        - *"Index on Columns"*: `prot`
>    - *"SQL Query to generate tabular output"*: `SELECT id,Proteins,Sequence FROM psms WHERE psms.ln NOT IN (SELECT distinct prots.ln FROM prots JOIN Uniprot ON prots.prot = Uniprot.prot) ORDER BY psms.ln`
>    - *"include query result column headers"*: `Yes`
>
>
{: .hands_on}


> <hands-on-title> Cutting out peptide sequences from Query Tabular with Cut </hands-on-title>
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c3`
>    - {% icon param-file %} *"From"*: `output` (output of **Query Tabular** {% icon tool %})
>
>
{: .hands_on}


> <hands-on-title> Grouping distinct (unique) peptides from SGPS with Group </hands-on-title>
>
> 1. {% tool [Group](Grouping1) %} with the following parameters:
>    - {% icon param-file %} *"Select data"*: `out_file1` (output of **Cut** {% icon tool %})
>    - *"Group by column"*: `c1`
>
>
{: .hands_on}


## Perform peptide discovery with **MaxQuant**
MaxQuant is an MS-based proteomics platform that is capable of processing raw data and provides improved mass precision and high precursor mass accuracy (HPMA), which resulted in increased protein identification and more in-depth proteomic analysis. Raw MS/MS spectra will be searched against the reduced MetaNovo-generated database (~21.2k sequences). More information about analysis using MaxQuant is available, including [Label-free data analysis](https://gxy.io/GTN:T00218) and [MaxQuant and MSstats for the analysis of TMT data](https://gxy.io/GTN:T00220).

> <hands-on-title> Peptide discovery using MaxQuant </hands-on-title>
>
> 1. {% tool [MaxQuant](toolshed.g2.bx.psu.edu/repos/galaxyp/maxquant/maxquant/2.0.3.0+galaxy0) %} with the following parameters:
>    - In *"Input Options"*:
>        - {% icon param-file %} *"FASTA files"*: `output` (Input dataset)
>    - In *"Search Options"*:
>        - {% icon param-file %} *"Specify an experimental design template (if needed). For detailed instructions see the help text."*: `output` (Input dataset)
>        - *"minimum peptide length"*: `8`
>        - *"Match between runs"*: `Yes`
>        - *"Maximum peptide length for unspecific searches"*: `50`
>    - In *"Protein quantification"*:
>        - *"Use only unmodified peptides"*: `Yes`
>            - *"Modifications used in protein quantification"*: `Oxidation (M)`
>        - In *"LFQ Options"*:
>            - *"iBAQ (calculates absolute protein abundances by normalizing to copy number and not protein mass)"*: `No`
>    - In *"Parameter Group"*:
>        - {% icon param-repeat %} *"Insert Parameter Group"*
>            - {% icon param-collection %} *"Infiles"*: `output` (Input dataset collection)
>            - *"fixed modifications"*: `Carbamidomethyl (C)`
>            - *"variable modifications"*: `Oxidation (M)`
>            - *"enzyme"*: `Trypsin/P`
>            - *"Quantitation Methods"*: `reporter ion MS2`
>                - *"isobaric labeling"*: `TMT11plex`
>                - *"Filter by PIF"*: `True`
>    - *"Generate PTXQC (proteomics quality control pipeline) report? (experimental setting)"*: `False`
>    - In *"Output Options"*:
>        - *"Select the desired outputs."*: `Protein Groups` `mqpar.xml` `Peptides` `MSMS` `msms scans` `summary` `MaxQuant and PTXQC log` `yaml config file`
>
>
>
{: .hands_on}


> <question-title></question-title>
>
> 1. What is the Experimental Design file for MaxQuant?
> 
> > <solution-title></solution-title>
> >
> > 1. In MaxQuant, the **Experimental Design** file is used to specify the experimental conditions, sample groups, and the relationships between different samples in a proteomics experiment. This file is a crucial component of the MaxQuant analysis process because it helps the software correctly organize and analyze the mass spectrometry data. The **Experimental Design** file typically has a ".txt" extension and is a tab-delimited text file. Here's what you might include in an Experimental Design file for MaxQuant: **Sample Names** (You specify the names of each sample in your experiment. These names should be consistent with the naming conventions used in your raw data files.), **Experimental Conditions** (You define the experimental conditions or treatment groups associated with each sample. For example, you might have control and treated groups, and you would assign the appropriate condition to each sample.), **Replicates** (You indicate the replicates for each sample, which is important for assessing the statistical significance of your results. Replicates are typically denoted by numeric values (e.g., "1," "2," "3") or by unique identifiers (e.g., "Replicate A," "Replicate B")), **Labels** (If you're using isobaric labeling methods like TMT (Tandem Mass Tag) or iTRAQ (Isobaric Tags for Relative and Absolute Quantitation), you specify the labels associated with each sample. This is important for quantification.), **Other Metadata** (You can include additional metadata relevant to your experiment, such as the biological source, time points, or any other information that helps describe the samples and experimental conditions.)
> >
> {: .solution}
>
{: .question}


## Using Text Manipulation Tools to Manage Microbial Outputs from MaxQuant

> <hands-on-title> Selecting microbial peptides from MaxQuant with Select </hands-on-title>
>
> 1. {% tool [Select](Grep1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `peptides` (output of **MaxQuant** {% icon tool %})
>    - *"that"*: `NOT Matching`
>    - *"the pattern"*: `(_HUMAN)|(_REVERSED)|(REV_)|(CON)|(con)`
>    - *"Keep header line"*: `Yes`
>
>
{: .hands_on}


> <hands-on-title> Cutting out microbial peptide sequences with Cut </hands-on-title>
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1`
>    - {% icon param-file %} *"From"*: `out_file1` (output of **Select** {% icon tool %})
>
>
{: .hands_on}


> <hands-on-title> Remove header line from MaxQuant peptide output with Remove beginning </hands-on-title>
>
> 1. {% tool [Remove beginning](Remove beginning1) %} with the following parameters:
>    - {% icon param-file %} *"from"*: `out_file1` (output of **Cut** {% icon tool %})
>
>
{: .hands_on}


> <hands-on-title> Grouping distinct (unique) peptide sequences from MaxQuant with Group </hands-on-title>
>
> 1. {% tool [Group](Grouping1) %} with the following parameters:
>    - {% icon param-file %} *"Select data"*: `out_file1` (output of **Remove beginning** {% icon tool %})
>    - *"Group by column"*: `c1`
>
>
{: .hands_on}


> <question-title></question-title>
>
> 1. How case-sensitive is the Group tool? Can I only group by column values, and not row values?
>
> > <solution-title></solution-title>
> >
> > 1. You can make it case sensitive, by default it is not. The tool here does column grouping only.
> >
> {: .solution}
>
{: .question}


## Process SGPS and MaxQuant peptides to compile one list of unique microbial peptides


> <hands-on-title> Concatenate SGPS and MaxQuant peptides into a singular database with Concatenate datasets </hands-on-title>
>
> 1. {% tool [Concatenate datasets](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cat/0.1.1) %} with the following parameters:
>    - {% icon param-files %} *"Datasets to concatenate"*: `out_file1` (output of **Group** {% icon tool %}), `out_file1` (output of **Group** {% icon tool %})
>
>
{: .hands_on}


> <hands-on-title> Group the peptides from SGPS and MaxQuant to remove duplicates with Group </hands-on-title>
>
> 1. {% tool [Group](Grouping1) %} with the following parameters:
>    - {% icon param-file %} *"Select data"*: `out_file1` (output of **Concatenate datasets** {% icon tool %})
>    - *"Group by column"*: `c1`
>
>
{: .hands_on}


# Conclusion
By following this tutorial, you have effectively conducted a search of your MS/MS data against the compact database and successfully retrieved reliable microbial peptides. After identifying these microbial peptides with the assistance of MaxQuant and SearchGUI, the next step is to verify the presence of these peptides. This compiled list of unique peptides will serve as the input for PepQuery to validate the confident identification of microbial peptides with the help of the verification workflow.
