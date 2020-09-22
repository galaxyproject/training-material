---
layout: tutorial_hands_on

title: metaQuantome Data creation tutorial
zenodo_link: '(https://doi.org/10.5281/zenodo.4037137)'
questions:
 - "How do I perform functional and taxonomy analysis on metaproteomics data?"
 - "How can I perform quantitation on metaproteomics data?"
 - "How do I create inputs that can be used in metaquantome to examine differentially expressed proteins?"
objectives:
  - "A taxonomy, functional and quantitational analysis of metaproteomic mass spectrometry data."
time_estimation: "1 hr"
key_points:
  - "Use dataset collections"
  - "With SearchGUI and PeptideShaker you can gain access to multiple search engines"
  - "Learning the basics of SQL queries can pay off"
contributors:
  - subinamehta
  - timothygriffin
  - pratikdjagtap 
  - emmaleith
  - mariecrane

---


# Introduction
{:.no_toc}

Metaproteomics involves the identification and analysis of microbial proteins at community level. The data
is used to determine taxonomic and functional state of the microbiome.Currently, metaproteomics related
research and bioinformatics softwares have several limitations, such as, supporting only spectral counts or 
its inability to co-relate functional and taxonomy interaction. The Galaxy-P team published a [metaQuantome tool](https://www.mcponline.org/content/18/8_suppl_1/S82), 
a multifarious package suite that leverages the taxonomic, functional and peptide level quantitative information 
to analyze the microbial community in different conditions. ![Microbiome](../../images/microbiome.png)

 
Across multiple experimental conditions, metaQuantome offers differential abundance analysis, principal 
components analysis, and clustered heat map visualizations. metaQuantome is an open source tool and 
available on the command line and in Galaxy making it accessible, flexible and reproducible. However, 
creating the data that is compatible with the metaQuantome suite is also not trivial. Hence, we developed 
a metaQuantome data creation workflow, wherein we create the inputs that are compatible with the metaquantome workflow
> ![metaQuantome workflow](../../images/metaquantomeworkflow.png)


To demonstrate the use of the data creation workflow, the metaproteomics data set came from a thermophilic 
biogas reactor which digests municipal food waste and manure (Fig1). After one round in the reactor, the 
microbial community was simplified and enriched via serial dilution to extinction. This inoculum was then 
transferred to a solution of cellulose from Norwegian Spruce and incubated at 65°C. Triplicate mRNA samples 
were taken in a time series from 0 to 43 hours after inoculation. For this training, we chose two time points-13 hour and 38 hour. 


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Pretreatments

The first step in a tutorial is to get the data from the zenodo link provided and making sure that it is in the correct format.

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files ( 6 MZML files, a Protein FASTA file and Experimental Design file) from [Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.4037137.svg)](https://doi.org/10.5281/zenodo.4037137) or from
>    the shared data library (`GTN - Material` -> `Proteomics` ->`MetaQuantome Datacreation`
>     -> `{{ page.title }}`):
>
>    ```
>    [Experimental Design](https://zenodo.org/record/4037137/files/ExperimentalDesign.tsv?download=1)
>    [Protein Database](https://zenodo.org/record/4037137/files/ProteinDB_cRAP.fasta?download=1)
>    [T2_A1 MZML File](https://zenodo.org/record/4037137/files/T2_A1.mzml?download=1)
>    [T2_B1 MZML File](https://zenodo.org/record/4037137/files/T2_B1.mzml?download=1)
>    [T7A_1 MZML File](https://zenodo.org/record/4037137/files/T7A_1.mzml?download=1)
>    [T7B_1 MZML File](https://zenodo.org/record/4037137/files/T7B_1.mzml?download=1)
>    ```
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
>
> 3. Build a **Dataset list** for the four mzml files.
>    - Click the **Operations on multiple datasets** check box at the top of the history panel
>       ![Operations on multiple datasets button](../../images/operations_icon.png)
>    - Check the four boxes next to the mzml files.
>    - Click **For all selected...** and choose **Build dataset list**
>
> 4. Rename the datasets (IF needed)
> 5. Check that the datatype ( Make sure they are in the correct formats).
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 6. Add to each database a tag corresponding to the name of the input data (optional).
{: .hands_on}

# Analysis

## Match peptide sequences

For this, the sequence database-searching program called [SearchGUI](https://compomics.github.io/projects/searchgui.html) will be used.
The created dataset collection of the four *MZML files* in the history has to be first converted to MGF to be used as the MS/MS input. 


## Convert mzml to MGF with **msconvert**
msconvert is used in order to convert the input file type, a mzml data collection, to a mgf file type. 
The mgf file type can then be used as the Input Peak Lists when running SearchGUI. 

> ### {% icon hands_on %} Hands-on: mzml to MGF
>
> 1. {% tool [msconvert](toolshed.g2.bx.psu.edu/repos/galaxyp/msconvert/msconvert/3.0.19052.0) %} with the following parameters:
>    - {% icon param-collection %} *"Input unrefined MS data"*: `output` (Input dataset collection)
>    - *"Do you agree to the vendor licenses?"*: `Yes`
>    - *"Output Type"*: `mgf`
>    - In *"Data Processing Filters"*:
>        - *"Apply peak picking?"*: `Yes`
>        - *"Apply m/z refinement with identification data?"*: `Yes`
>        - *"(Re-)calculate charge states?"*: `no`
>        - *"Filter m/z Window"*: `Yes`
>        - *"Filter out ETD precursor peaks?"*: `Yes`
>        - *"De-noise MS2 with moving window filter"*: `Yes`
>    - In *"Scan Inclusion/Exclusion Filters"*:
>        - *"Filter MS Levels"*: `Yes`
>    - In *"General Options"*:
>        - *"Sum adjacent scans"*: `Yes`
>        - *"Output multiple runs per file"*: `Yes`
>
>
>    > ### {% icon comment %} Comment
>    >This is a critical step for running this workflow. 
>    {: .comment}
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. Why do we need to convert the files to MGF?
> 2. Can we use any other input format?
>
> > ### {% icon solution %} Solution
> >
> > 1. The files have to be converted to MGF for this workflow because we use SearchGUI as the searching tool and it can only read MGF files
> > 2. Yes, we can also use RAW files as input and just convert RAW files to MGF.
> >
> {: .solution}
>
{: .question}

## **Search GUI**

> ### {% icon hands_on %} Hands-on: 

SearchGUI is a tool that searches sequence databases on any number of MGF files. In this case, the previously made collection of three MGF files (entitles MGF files) will be used as the MS/MS input. This tool will produce an output file, called a SearchGUI archive file. This file will serve as in input for the next tool used, PeptideShaker.
>
> 1. {% tool [Search GUI](toolshed.g2.bx.psu.edu/repos/galaxyp/peptideshaker/search_gui/3.3.10.1) %} with the following parameters:
>    - {% icon param-file %} *"Protein Database"*: `output` (Input dataset collection)
>    - {% icon param-file %} *"Input Peak Lists (mgf)"*: `output` (output of **msconvert** {% icon tool %})
>    - In *"Search Engine Options"*:
>        - *"DB-Search Engines"*: ``
>    - In *"Protein Digestion Options"*:
>        - *"Digestion"*: `Select Enzymes`
>    - In *"Precursor Options"*:
>        - *"Fragment Tolerance"*: `0.2`
>        - *"Maximum Charge"*: `6`
>    - In *"Protein Modification Options"*:
>        - *"Fixed Modifications"*: ``
>        - *"Variable Modifications"*: ``
>    - In *"Andvanced Options"*:
>        - *"SearchGUI Options"*: `Default`
>        - *"X!Tandem Options"*: `Advanced`
>            - *"X!Tandem: Quick Acetyl"*: `Yes`
>            - *"X!Tandem: Quick Pyrolidone"*: `Yes`
>            - *"X!Tandem: Maximum Valid Expectation Value"*: `100.0`
>            - *"X!Tandem peptide model refinement"*: `Don't refine`
>        - *"OMSSA Options"*: `Default`
>        - *"MSGF Options"*: `Default`
>        - *"MS Amanda Options"*: `Default`
>        - *"TIDE Options"*: `Default`
>        - *"MyriMatch Options"*: `Default`
>        - *"Comet Options"*: `Default`
>        - *"DirectTag Options"*: `Default`
>        - *"Novor Options"*: `Default`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > > Note that sequence databases used for metaproteomics are usually much larger than the excerpt used in this tutorial. When using large databases, the peptide identification step can take much more time for computation. In metaproteomics, choosing the optimal database is a crucial step of your workflow, for further reading see [Timmins-Schiffman et al (2017)](https://www.ncbi.nlm.nih.gov/pubmed/27824341).
>
> To learn more about database construction in general, like integrating contaminant databases or using a decoy strategy for FDR searching, please consult our tutorial on [Database Handling]({{site.baseurl}}/topics/proteomics/tutorials/database-handling/tutorial.html).
>
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. How many Search Engines can be used?
> 2. Can the parameters be manipulated?
>
> > ### {% icon solution %} Solution
> >
> > 1. There are 8 database search algorithms, you can use as many as you want. Ideally, 4 database algorithms gives the best results.
> > 2. Yes, The parameters can be manipulated according to the experimental design of the datasets.
> >
> {: .solution}
>
{: .question}

## **Peptide Shaker**

> ### {% icon hands_on %} Hands-on: 
[PeptideShaker](https://compomics.github.io/projects/peptide-shaker.html) is a post-processing software tool that processes data from the SearchGUI software tool. PeptideShaker is a search engine for interpretation of proteomics identification results from multiple search engines, currently supporting X!Tandem, MS-GF+, MS Amanda, OMSSA, MyriMatch, Comet, Tide, Mascot, Andromeda and mzIdentML. More specifically, PeptideShaker processes data from  the SearchGUI tool through the organization of Peptide-Spectral Matches (PSMs) generated. In addition to organization, it provides an assessment of confidence of the data and generates outputs that can be visualized by users to interpret the results. 
>
> 1. {% tool [Peptide Shaker](toolshed.g2.bx.psu.edu/repos/galaxyp/peptideshaker/peptide_shaker/1.16.36.3) %} with the following parameters:
>    - {% icon param-file %} *"Compressed SearchGUI results"*: `searchgui_results` (output of **Search GUI** {% icon tool %})
>    - *"Specify Advanced PeptideShaker Processing Options"*: `Advanced Processing Options`
>        - *"The PTM probabilistic score to use for PTM localization"*: `A-score`
>    - *"Specify Advanced Filtering Options"*: `Advanced Filtering Options`
>        - *"Maximum Peptide Length"*: `60`
>    - *"Specify Contact Information for mzIdendML"*: `GalaxyP Project contact (Not suitable for PRIDE submission)`
>    - In *"Exporting options"*:
>        - *"Creates a mzIdentML file"*: `Yes`
>        - *"Compress results into a single zip file"*: `Yes`
>        - *"Reports to be generated"*: ``
>
>
>    > ### {% icon comment %} Comment
>    >
>    >  There are a number of choices for different data files that can be generated using
> PeptideShaker. A compressed file can be made containing all information needed to view the
> results in the standalone PeptideShaker viewer. A `mzidentML` file can be created that contains
> all peptide sequence matching information and can be utilized by compatible downstream
> software. Other outputs are focused on the inferred proteins identified from the PSMs, as well
> as phosphorylation reports, relevant if a phosphoproteomics experiment has been undertaken.
> More detailed information on peptide inference using SearchGUI and PeptideShaker can be found in our tutorial on [Peptide and Protein ID]({{site.baseurl}}/topics/proteomics/tutorials/protein-id-sg-ps/tutorial.html).
>    {: .comment}
>
{: .hands_on}


##  **Select**

> ### {% icon hands_on %} Hands-on: 
> This Select tool is used to remove all the contaminants from the Peptide Spectral Match (PSM) search results.
>
> 1. {% tool [Select](Grep1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `output_psm` (output of **Peptide Shaker** {% icon tool %})
>    - *"that"*: `NOT Matching`
>    - *"the pattern"*: `con_`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. Why is removing contaminants important?
>
> > ### {% icon solution %} Solution
> >
> > 1. Ideally, we would like to remove known contaminants from our samples just to maintain discovering novel proteoforms in our sample.
> >
> {: .solution}
>
{: .question}

## **Select**

> ### {% icon hands_on %} Hands-on: 
> This Select tool is used to remove all the contaminants from the Peptide report(PSM) obtained from Peptide Shaker.
>
> 1. {% tool [Select](Grep1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `output_peptides` (output of **Peptide Shaker** {% icon tool %})
>    - *"that"*: `NOT Matching`
>    - *"the pattern"*: `con_`
>
{: .hands_on}


## **Replace Text**

> ### {% icon hands_on %} Hands-on: 
This is a data manipulation step to make the data compatible with other downstream processing tools. The Replace text tool replaces the .mgf extention from the PSM report so that it can be used as an input for FlashLFQ.
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `out_file1` (output of **Select** {% icon tool %})
>    - In *"Replacement"*:
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"in column"*: `c10`
>            - *"Find pattern"*: `.mgf`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > Replace Text searches given columns and finds and replaces patterns provided by the user. 
This tool is removing the extensions (.raw,.mzml,.mgf) provided by the PeptideShaker tool. This step is critical for FlashLFQ to work.
>    {: .comment}
>
{: .hands_on}

## **Cut**

> ### {% icon hands_on %} Hands-on: 
This step cuts the peptide lest from the PSM report with no contaminants.
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c6`
>    - {% icon param-file %} *"From"*: `out_file1` (output of **Select** {% icon tool %})
>
>
{: .hands_on}

## **FlashLFQ**

> ### {% icon hands_on %} Hands-on: 
FlashLFQ can quantify MS peaks in order to find the abundances of peptides. Additionally, the abundances of peptides within the sample can be compared between samples as further analysis beyond this workflow.  
>
> 1. {% tool [FlashLFQ](toolshed.g2.bx.psu.edu/repos/galaxyp/flashlfq/flashlfq/1.0.3.0) %} with the following parameters:
>    - {% icon param-file %} *"identification file"*: `outfile` (output of **Replace Text** {% icon tool %})
>    - {% icon param-file %} *"spectrum files"*: `output` (output of **msconvert** {% icon tool %})
>    - *"match between runs"*: `Yes`
>    - *"Use experimnetal design for normalization or protein fold-change analysis"*: `Yes`
>        - {% icon param-file %} *"ExperimentalDesign.tsv"*: `output` (Input dataset)
>        - *"Perform Bayesian protein fold-change analysis"*: `Yes`
>            - *"control condition for Bayesian protein fold-change analysis"*: ``
>
>
>    > ### {% icon comment %} Comment
>    >
>    > [FlashLFQ](https://github.com/smith-chem-wisc/FlashLFQ) is a label-free quantification tool for mass-spectrometry proteomics. It supports both .mzML and Thermo .raw file formats. 
>    {: .comment}
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. Can it be used with fractionated data?
> 2. Does it perform peptide and protein level quantification?
>
> > ### {% icon solution %} Solution
> >
> > 1. Yes, this tool can be used with fractionated datasets and multiple conditions
> > 2. It performed both peptide level and protein level quantification, For protein level it used Bayesian Fold change analysis.
> >
> {: .solution}
>
{: .question}

## **Filter**

> ### {% icon hands_on %} Hands-on: 
> This is a data manipulation tool. Here, we select those peptides with less than 50 amino acids in length.
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `out_file1` (output of **Cut** {% icon tool %})
>    - *"With following condition"*: `len(c1)<=50`
>    - *"Number of header lines to skip"*: `1`
>
>    > ### {% icon comment %} Comment
>    > Unipept fails with peptides more than 50 amino acids in length, thus we decided to work with peptides that are less than 50 amino acids.
>    > 
>    {: .comment}
>

## Sub-step with **Regex Find And Replace**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `quantifiedPeptides` (output of **FlashLFQ** {% icon tool %})
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Base Sequence`
>            - *"Replacement"*: `peptide`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Intensity_`
>
{: .hands_on}


## Sub-step with **Unipept**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Unipept](toolshed.g2.bx.psu.edu/repos/galaxyp/unipept/unipept/4.3.0) %} with the following parameters:
>    - *"Unipept application"*: `peptinfo: Tryptic peptides and associated EC and GO terms and lowest common ancestor taxonomy`
>        - *"Equate isoleucine and leucine"*: `Yes`
>        - *"retrieve extra information"*: `Yes`
>        - *"group responses by GO namespace (biological process, molecular function, cellular component)"*: `Yes`
>        - *"allfields"*: `Yes`
>    - *"Peptides input format"*: `tabular`
>        - {% icon param-file %} *"Tabular Input Containing Peptide column"*: `out_file1` (output of **Filter** {% icon tool %})
>        - *"Select column with peptides"*: `c1`
>    - *"Choose outputs"*: ``
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Unipept**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Unipept](toolshed.g2.bx.psu.edu/repos/galaxyp/unipept/unipept/4.0.0) %} with the following parameters:
>    - *"Unipept application"*: `pept2lca: lowest common ancestor`
>        - *"Equate isoleucine and leucine"*: `Yes`
>        - *"allfields"*: `Yes`
>    - *"Peptides input format"*: `tabular`
>        - {% icon param-file %} *"Tabular Input Containing Peptide column"*: `out_file1` (output of **Filter** {% icon tool %})
>        - *"Select column with peptides"*: `c1`
>    - *"Choose outputs"*: ``
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}


## Sub-step with **Cut**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1,c3`
>    - {% icon param-file %} *"From"*: `output_ec_tsv` (output of **Unipept** {% icon tool %})
>
{: .hands_on}


## Sub-step with **Query Tabular**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.0.0) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `output_go_tsv` (output of **Unipept** {% icon tool %})
>            - In *"Table Options"*:
>                - *"Specify Name for Table"*: `Goterm`
>                - *"Specify Column Names (comma-separated list)"*: `peptide,total_protein_count,go_term,protein_count,go_name,go_funct`
>    - *"SQL Query to generate tabular output"*: `SELECT Goterm.*
```
FROM Goterm 
WHERE ((1.0*Goterm.protein_count)/(1.0*Goterm.total_protein_count)) >= 0.05
```

>    - *"include query result column headers"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}


## Sub-step with **Replace Text**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_line/1.1.2) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `output_tsv` (output of **Unipept** {% icon tool %})
>    - In *"Replacement"*:
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `#peptide`
>            - *"Replace with:"*: `peptide`
>
>
>    > ### {% icon comment %} Comment
>    >
>    > This step is to remove the hashtag (#) from the peptide info tsv output from Unipept.
>    {: .comment}
>
{: .hands_on}


## Sub-step with **Query Tabular**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.0.0) %} with the following parameters:
>    - In *"Database Table"*:
>        - {% icon param-repeat %} *"Insert Database Table"*
>            - {% icon param-file %} *"Tabular Dataset for Table"*: `out_file1` (output of **Cut** {% icon tool %})
>            - In *"Table Options"*:
>                - *"Specify Name for Table"*: `ec`
>                - *"Specify Column Names (comma-separated list)"*: `peptide,go_ec`
>    - *"SQL Query to generate tabular output"*: `SELECT *
FROM ec`
>    - *"include query result column headers"*: `Yes`
>
{: .hands_on}


## Sub-step with **Filter**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `output` (output of **Query Tabular** {% icon tool %})
>    - *"With following condition"*: `c5=='biological process'`
>    - *"Number of header lines to skip"*: `1`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > This tool helps extract peptides assigned to their  biological process.
>    {: .comment}
>
{: .hands_on}


## Sub-step with **Filter**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `output` (output of **Query Tabular** {% icon tool %})
>    - *"With following condition"*: `c5=='cellular component'`
>    - *"Number of header lines to skip"*: `1`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > This tool helps extract peptides assigned to their Cellular component.
>    {: .comment}
>
{: .hands_on}


## Sub-step with **Filter**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `output` (output of **Query Tabular** {% icon tool %})
>    - *"With following condition"*: `c5=='molecular function'`
>    - *"Number of header lines to skip"*: `1`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}



# Conclusion
{:.no_toc}

This metaQuantome data create tutorial is a guide to have compatible dataset that are metaQuantome ready. We have incorporated two conditions in this workflow but users can use as many as they want. 
