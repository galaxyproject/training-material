---
layout: tutorial_hands_on

title: metaQuantome Data creation tutorial
zenodo_link: "https://doi.org/10.5281/zenodo.4037137"
questions:
 - "How do I perform functional and taxonomy analysis on metaproteomics data?"
 - "How can I perform quantitation on metaproteomics data?"
 - "How do I create inputs that can be used in metaquantome to examine differentially expressed proteins?"
objectives:
  - "A taxonomy, functional and quantitational analysis of metaproteomic mass spectrometry data."
time_estimation: "1h"
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

Metaproteomics involves characterization of community level expression of microbial proteins from an environmental 
or clinical sample. Metaproteomics data is primarily used to determine the functional status of the microbiome under 
study along with its taxonomic composition. The Galaxy-P team published a software suite named [metaQuantome](https://www.mcponline.org/content/18/8_suppl_1/S82) to enable quantitative and statistical analysis and visualization of functional, taxonomic expression as well as functional and 
taxonomy interaction. metaQuantome leverages peptide level quantitative information to analyze the taxonomic, functional 
expression within the microbial community in different conditions.

<p align="center">
  <img width="850" height="500" src="../../images/microbiome.png" alt="Microbiome" title="Microbiome" />
</p>

 
metaQuantome offers differential abundance analysis, principal components analysis, and clustered heat map visualizations, 
across multiple experimental conditions. metaQuantome, an open source tool, is available via command line and also 
accessible via Galaxy platform for reproducible analysis. As a first step for metaQuantome analysis, metaproteomics 
data needs to be made compatible for subsequent analysis. With this in mind, we have developed a metaQuantome data 
generation workflow tutorial that ill help users generate inputs for metaQuantome analysis.

<p align="center">
  <img width="850" height="500" src="../../images/metaquantomeworkflow.png" alt="Workflow" title="Workflow">
</p>

To demonstrate the use of the data creation workflow, we have used a thermophilic biogas reactor dataset wherein municipal 
food waste and manure is digested to generate methane gas. After one round in the reactor, the microbial community was 
simplified and enriched via serial dilution. This inoculum was then transferred to a solution of cellulose from Norwegian 
Spruce and incubated at 65°C. Triplicate samples were taken in a time series from 0 to 43 hours after inoculation and mass 
spectrometry data was acquired on a Q-Exactive (Thermo) mass spectrometer. For this training, we have chosen two time points-8 hour and 33 hour.

<p align="center">
  <img width="800" height="300" src="../../images/biogasdataset.png" alt="Dataset" title="Dataset">
</p>



> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# **Pretreatments**

The first step in a tutorial is to get the data from the zenodo link provided and making sure that it is in the correct format.

## *Get data*

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


# **Match peptide sequences**

For this, the sequence database-searching program called [SearchGUI](https://compomics.github.io/projects/searchgui.html) will be used.
The created dataset collection of the four *MZML files* in the history has to be first converted to MGF to be used as the MS/MS input. 


### *Convert mzml to MGF with msconvert*

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
> > 1. The files have to be converted to MGF for this workflow because we use SearchGUI as the searching tool and it can only read MGF files.
> > 2. Yes, we can also use RAW files as input and just convert RAW files to MGF.
> >
> {: .solution}
>
{: .question}

##  *Search GUI*
SearchGUI is a tool that searches sequence databases on any number of MGF files. In this case, the previously made collection of three MGF files (entitles MGF files) will be used as the MS/MS input. This tool will produce an output file, called a SearchGUI archive file. This file will serve as in input for the next tool used, PeptideShaker.
>
> ### {% icon hands_on %} Hands-on: 
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
>
>
>    > ### {% icon comment %} Comment
>    >
>    >  Note that sequence databases used for metaproteomics are usually much larger than the excerpt used in this tutorial. When using large databases, the peptide identification step can take much more time for computation. In metaproteomics, choosing the optimal database is a crucial step of your workflow, for further reading see [Timmins-Schiffman et al (2017)](https://www.ncbi.nlm.nih.gov/pubmed/27824341). To learn more about database construction in general, like integrating contaminant databases or using a decoy strategy for FDR searching, please consult our tutorial on [Database Handling]({{site.baseurl}}/topics/proteomics/tutorials/database-handling/tutorial.html).
>    {: .comment}
>
{: .hands_on}


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

##  *Peptide Shaker*

[PeptideShaker](https://compomics.github.io/projects/peptide-shaker.html) is a post-processing software tool that processes data from the SearchGUI software tool. PeptideShaker is a search engine for interpretation of proteomics identification results from multiple search engines, currently supporting X!Tandem, MS-GF+, MS Amanda, OMSSA, MyriMatch, Comet, Tide, Mascot, Andromeda and mzIdentML. More specifically, PeptideShaker processes data from  the SearchGUI tool through the organization of Peptide-Spectral Matches (PSMs) generated. In addition to organization, it provides an assessment of confidence of the data and generates outputs that can be visualized by users to interpret the results. 
>
> ### {% icon hands_on %} Hands-on: 
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
>
>    > ### {% icon comment %} Comment
>    >
>    >  There are a number of choices for different data files that can be generated using
 PeptideShaker. A compressed file can be made containing all information needed to view the
results in the standalone PeptideShaker viewer. A `mzidentML` file can be created that contains
all peptide sequence matching information and can be utilized by compatible downstream
software. Other outputs are focused on the inferred proteins identified from the PSMs, as well
as phosphorylation reports, relevant if a phosphoproteomics experiment has been undertaken.
More detailed information on peptide inference using SearchGUI and PeptideShaker can be found in 
our tutorial on [Peptide and Protein ID]({{site.baseurl}}/topics/proteomics/tutorials/protein-id-sg-ps/tutorial.html).
>    {: .comment}
>
{: .hands_on}


##  *Select*

> ### {% icon hands_on %} Hands-on: 
> This Select tool is used to remove all the contaminants from the Peptide Spectral Match (PSM) search results.
>
> 1. {% tool [Select](Grep1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `output_psm` (output of **Peptide Shaker** {% icon tool %})
>    - *"that"*: `NOT Matching`
>    - *"the pattern"*: `con_`
>
>
>
>    > ### {% icon comment %} Comment
>    >
>    > In Proteomics, contamination is generally detected as peaks in spectra that did not originate 
from the samples and can be introduced in the sample from a variety of environmental sources or human error. Identification of these 
contaminants is critical to enable their removal before data analysis, mainly, to maintain the validity of conclusions 
drawn from statistical analyses. Thus, this selection tool helps us remove the contaminants that were identified in the spectral data. 
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

## *Select*

This Select tool is used to remove all the contaminants from the Peptide report obtained from Peptide Shaker.

> ### {% icon hands_on %} Hands-on: 
> 
>
> 1. {% tool [Select](Grep1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `output_peptides` (output of **Peptide Shaker** {% icon tool %})
>    - *"that"*: `NOT Matching`
>    - *"the pattern"*: `con_`
>
{: .hands_on}


## *Replace Text*
This is a data manipulation step to make the data compatible with other downstream processing tools. The Replace text tool replaces the .mgf extention from the PSM report so that it can be used as an input for FlashLFQ.
>
> ### {% icon hands_on %} Hands-on: 
>
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `out_file1` (output of **Select** {% icon tool %})
>    - In *"Replacement"*:
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"in column"*: `c10`
>            - *"Find pattern"*: `.mgf`
>
>
>
>    > ### {% icon comment %} Comment
>    >
>    > Replace Text searches given columns and finds and replaces patterns provided by the user. 
This tool is removing the extensions (.raw,.mzml,.mgf) in the spectral file column provided by the PeptideShaker tool. This step is critical for FlashLFQ to work.
>    {: .comment}
>
{: .hands_on}

## *Cut*

> ### {% icon hands_on %} Hands-on: 
This step selects the peptide column from the Select output ( where we have removed the contaminants)
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c6`
>    - {% icon param-file %} *"From"*: `out_file1` (output of **Select** {% icon tool %})
>
>
{: .hands_on}

# **Peptide Quantification**

In this tutorial, we are using FlashLFQ as the quantitation tool. The user can choose to work with other quantitation tools, For eg: moFF and MaxQuant are available in Galaxy. 

### *FlashLFQ*
FlashLFQ can quantify MS peaks in order to find the abundances of peptides. Additionally, the abundances of peptides within the sample can be compared between samples as further analysis beyond this workflow.  
>
> ### {% icon hands_on %} Hands-on: 
>
> 1. {% tool [FlashLFQ](toolshed.g2.bx.psu.edu/repos/galaxyp/flashlfq/flashlfq/1.0.3.0) %} with the following parameters:
>    - {% icon param-file %} *"identification file"*: `outfile` (output of **Replace Text** {% icon tool %})
>    - {% icon param-file %} *"spectrum files"*: `output` (output of **msconvert** {% icon tool %})
>    - *"match between runs"*: `Yes`
>    - *"Use experimental design for normalization or protein fold-change analysis"*: `Yes`
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

## *Filter*

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
{: .hands_on}


## *Regex Find And Replace*

> ### {% icon hands_on %} Hands-on: 
> Regex Find And Replace goes line by line through the input file and will remove any patterns specified by the user and replace them with expressions also specified by the user. In this case, Regex Find And Replace is being used on a FlashLFQ output file and manipulating the header to make it compatible with metaQuantome alongwith completely removing the N-terminus and C-terminus tag in the peptide sequences.
>
> 1. {% tool [Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `quantifiedPeptides` (output of **FlashLFQ** {% icon tool %})
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Base Sequence`
>            - *"Replacement"*: `peptide`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Intensity_`
>        - {% icon param-repeat %} *"Insert Check"*  
>            - *”Find Regex”*: `NH2-`
>            - *”Replacement”*: ``
>        - {% icon param-repeat %} *"Insert Check"*  
>            - *”Find Regex”*: `-COOH`
>            - *”Replacement”*: ``
>
{: .hands_on}

# **Functional and Taxonomy annotation**


## *Unipept* for taxonomy annotation

Unipept is used again to match tryptic peptides and find the taxonomy and lowest common ancestor of each peptide. 

> ### {% icon hands_on %} Hands-on: 
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
>
>
>    > ### {% icon comment %} Comment
>    >
>    > There are two Unipept in this workflow, One for taxonomy and other for function.
>    {: .comment}
>
{: .hands_on}

<p align="center">
  <img width="800" height="400" src="../../images/taxa.png" alt="Taxa" title="Taxonomy">
</p>

> ### {% icon question %} Questions
>
> 1. Can any other taxonomy and functional tool be used apart from Unipept?
>
> > ### {% icon solution %} Solution
> >
> > 1. Yes, any tool can be used for taxonomy and functional output. Please make sure the output has the information that incluldes peptide,taxon_name, taxon_id, genus, species etc.
> >
> {: .solution}
>
{: .question}

The JSON output from the Taxonomy can be visualized using the visualize option and Select the Unipept Taxonomyviewer. 
<p align="center">
  <img width="250" height="300" src="../../images/UnipeptJSON.png" alt="Unipept-JSON" title="Unipept-JSON">
</p>
<p align="center">
  <img width="250" height="125" src="../../images/unipept_taxonomy_viewer.png" alt="Taxa-viewer" title="Taxa-viewer">
</p>
<p align="center">
  <img width="700" height="800" src="../../images/UnipeptJSONoutput.png" alt="Output" title="Taxa-output">
</p>

## *Unipept* for Functional annotation

Unipept is used to match tryptic peptides and find the taxonomy and Functional annotation of the peptides. Unipept is used to match sample tryptic peptides to proteins using a fast-matching algorithm. Although Unipept can be accessed and used through the web page, the use of Unipept on galaxy allows the production of output datasets including the peptide information to be used in sequential steps. Unipept requires a list containing the peptide sequences which was generated by Query Tabular.

> ### {% icon hands_on %} Hands-on: 
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
>
>    > ### {% icon comment %} Comment
>    >
>    > There are two Unipept in this workflow, One for taxonomy and other for function. Please select all the output options from Unipept.
>    {: .comment}
>
The JSON output from the Taxonomy can be visualized using the visualize option and Select the Unipept Taxonomyviewer. 
>
>
{: .hands_on}

## *Cut*

> ### {% icon hands_on %} Hands-on: 
The cut tool cuts out specific columns from the dataset. In this case, the cut tool is being used to extract columns 1 (peptide) and 3 (EC number) from the dataset peptinfo EC.tsv output. This is a manipulation tool for metaQuantome's convinience.
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1,c3`
>    - {% icon param-file %} *"From"*: `output_ec_tsv` (output of **Unipept** {% icon tool %})
>
{: .hands_on}


## *Query Tabular*

Query Tabular is a tool that can load tabular data into a SQLite database. This step precedes UniPept, as a list containing the peptide sequences must be generated. In this step a list of gene ontology (GO) terms is being generated.

> ### {% icon hands_on %} Hands-on: 
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
>
>    > ### {% icon comment %} Comment
>    >
>    > In the Unipept API output, the threshold is set to 0.5% of the overall number of peptides unambiguously assigned to a taxon at a particular taxonomic rank level. Here in the Galaxy platform, we are using Query tabular to perform this filtering. 
>    {: .comment}
>
{: .hands_on}


## *Replace Text*

This step is to remove the hashtag from the Peptide header in the Unipept output.

> ### {% icon hands_on %} Hands-on: >
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_line/1.1.2) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `output_tsv` (output of **Unipept** {% icon tool %})
>    - In *"Replacement"*:
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `#peptide`
>            - *"Replace with:"*: `peptide`
>
>
{: .hands_on}


## *Query Tabular*

We are using this Query tabular ot rename the output that we obtained from the Cut column tool.

> ### {% icon hands_on %} Hands-on: 
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


## *Filter* - Biological Functions

The filter tool allows restriction of the dataset using simple conditional statements. This step is used to filter out the GO terms with biological processes and the corresponding number of peptides associated with these terms.


> ### {% icon hands_on %} Hands-on: 
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `output` (output of **Query Tabular** {% icon tool %})
>    - *"With following condition"*: `c5=='biological process'`
>    - *"Number of header lines to skip"*: `1`
>
{: .hands_on}

<p align="center">
  <img width="800" height="600" src="../../images/biologicalprocess.png" alt="Biological-Processes" title="Biological-Processes">
</p>



## *Filter* - Cellular components

This step is used to filter out the GO terms with cellular components and the corresponding number of peptides associated with these terms.

> ### {% icon hands_on %} Hands-on: 
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `output` (output of **Query Tabular** {% icon tool %})
>    - *"With following condition"*: `c5=='cellular component'`
>    - *"Number of header lines to skip"*: `1`
>
{: .hands_on}

<p align="center">
  <img width="800" height="500" src="../../images/cellularcomponent.png" alt="Cellular-Component" title="Cellular-Component">
</p>

## *Filter* - Molecular Function

This step is used to filter out the GO terms with molecular function and the corresponding number of peptides associated with these terms.

> ### {% icon hands_on %} Hands-on: 
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `output` (output of **Query Tabular** {% icon tool %})
>    - *"With following condition"*: `c5=='molecular function'`
>    - *"Number of header lines to skip"*: `1`
>
{: .hands_on}

<p align="center">
  <img width="800" height="500" src="../../images/molecularfunction.png" alt="Molecular-Function" title="Molecular-Function">
</p>


# **Conclusion**
{:.no_toc}

This completes the walkthrough of the metaQuantome data creation workflow .This tutorial is a guide to have datasets that are metaQuantome ready/compatible and can be used for metaproteomics research. We have incorporated only two conditions in this workflow but users can use as many as they want. Researchers can use this workflow with their data also, please make sure the tool parameters and the workflow will be needed to be modified accordingly.

Please look at the following tutorials in this Metaproteomics series:

[Metaproteomics]({% link topics/proteomics/tutorials/metaproteomics/tutorial.md %})
This workflow is also available at [proteomics.usegalaxy.eu](https://proteomics.usegalaxy.eu/)

This workflow was developed by the Galaxy-P team at the University of Minnesota. For more information about Galaxy-P or our ongoing work, please visit us at galaxyp.org

> ### {% icon comment %} References
>
>
> - [Galaxy workflows for metaproteomics](https://www.ncbi.nlm.nih.gov/pubmed/26058579)
>
> - [Metaproteomics community effort](https://z.umn.edu/gcc2017mporal)
>
> - [Unipept](https://www.ncbi.nlm.nih.gov/pubmed/28552653)
>
> - [Galaxy-P Metaproteomics instance](https://proteomics.usegalaxy.eu/)
>
> - [Metaproteomics video](http://z.umn.edu/mpvideo2018)
{: .comment}
