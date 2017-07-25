---
layout: tutorial_hands_on
topic_name: proteomics
tutorial_name: metaproteomics
---

# Introduction
{:.no_toc}

In this metaproteomics tutorial we will identify expressed proteins from a complex bacterial community sample.
For this MS/MS data will be matched to peptide sequences provided through a FASTA file.

Metaproteomics is the large-scale characterization of the entire protein complement of environmental microbiota
at a given point in time. It has the potential to unravel the mechanistic details of microbial interactions with
the host / environment by analyzing the functional dynamics of the microbiome.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Pretreatments

## Data upload

There are a many ways how you can upload your data. Three among these are:

*   Upload the files from your computer
*   Using a direct weblink
*   Import from the data library if your instance provides the files

In this tutorial, we will get the data from Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.820308.svg)](https://doi.org/10.5281/zenodo.820308). 

> ### :pencil2: Hands-on: Data upload and organization
>
> 1. Create a new history and name it something meaningful (e.g. *Metaproteomics tutorial*)
> 2. Import the three MGF MS/MS files and the FASTA sequence file from Zenodo.
>
>    > ### :bulb: Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field. You can add multiple links, each on a seperate line.
>    > * Press **Start**    
>    {: .tip}
>
>    As default, Galaxy takes the link as name.
>
>    > ### :nut_and_bolt: Comments
>    > - Rename the datasets to a more descriptive name
>    {: .comment}
>
> 3. Build a **Dataset list** for the three MGF files
>    - Click the **Operations on multiple datasets** check box at the top of the history panel
>       ![](../../images/operations_icon.png)
>    - Check the three boxes next to the MGF files
>    - Click **For all selected...** and choose **Build dataset list**
>    - Ensure the three control samples are the only ones selected, and enter a name for the new collection (e.g. *MGF files*)
>    - Click **Create list** and exit by clicking again the dataset operations icon
>
{: .hands_on}

# Analysis

## Match peptide sequences

The search database labelled `FASTA_Bering_Strait_Trimmed_metapeptides_cRAP.FASTA` is the input database that
will be used to match MS/MS to peptide sequences via a sequence database search.
For this, the sequence database-searching program called **SearchGUI** will be used.
The created dataset collection of the three *MGF files* in the history is used as the MS/MS input.

#### SearchGUI

> ### :pencil2: Hands-on: SearchGUI
>
> 1. **SearchGUI** :wrench:: Run **SearchGUI** with:
>    - **Protein Database**: `FASTA_Bering_Strait_Trimmed_metapeptides_cRAP.FASTA`(or however you named the `FASTA` file)
>    - **Input Peak lists (mgf)**: `MGF files` dataset collection. 
>
>    > ### :bulb: Tip: Select dataset collections as input
>    >
>    > * Click the **Dataset collection** icon on the left of the input field:
>    >  
>    >      ![](../../images/dataset_button.png)
>    > * Select the appropriate dataset collection from the list
>    {: .tip}    
>    
>    Section **Search Engine Options**:
>
>    - **B-Search Engines**: `X!Tandem`
>
>    > ### :nut_and_bolt: Comment
>    >
>    > The section **Search Engine Options** contains a selection of sequence database searching
>    > programs that are available in SearchGUI. Any combination of these programs can be used for
>    > generating PSMs from MS/MS data. For the purpose of this tutorial, **X!Tandem** we will be used.
>    {: .comment}
>
>    Section **Protein Digestion Options**:
>
>    - **Digestion**: `Trypsin`
>    - **Maximum Missed Cleavages**: `2`
>
>    Section **Precursor Options**:
>
>    - **Precursor Ion Tolerance Units**: `Parts per million (ppm)`
>    - **Precursor Ion Tolerance:**: `10`
>    - **Fragment Tolerance (Daltons)**: `0.02`- this is high resolution MS/MS data
>    - **Minimum charge**: `2`
>    - **Maximum charge**: `6`
>    - **Forward Ion**: `b`
>    - **Reverse Ion**: `y`
>    - **Minimum Precursor Isotope**: `0`
>    - **Maximum Precursor Isotope**: `1`
>
>    Section **Protein Modification Options**:
>
>    - **Fixed Modifications**: `Carbamidomethylation of C`
>    - **Variable modifications**: `Oxidation of M`
>
>    > ### :bulb: Tip: Search for options
>    >
>    > * For selection lists, typing the first few letters in the window will filter the available options.
>    {: .tip}
>
>    Section **Advanced Options**:
>    - **X!Tandem Options**: `Advanced`
>    - **X!Tandem: Total Peaks**: `50`
>    - **X!Tandem: Min Peaks**: `15`
>    - **X!Tandem: Min Frag m/z**: `200`
>    - **X!Tandem: Min Precursor Mass**: `200`
>    - **X!Tandem: Noise Suppression**: `Yes`
>    - **X!Tandem: Dynamic Range**: `100`
>    - **X!Tandem: Quick Acetyl**: `No`
>    - **X!Tandem: Quick Pyrolidone**: `No`
>    - **X!Tandem: Protein stP Bias**: `No`
>    - **X!Tandem: Maximum Valid Expectation Value**: `100`
>    - **X!Tandem: Output Proteins**: `No`
>    - **X!Tandem: Output Sequences**: `No`
>    - **X!Tandem: Output Spectra**: `Yes`
>    - **X!Tandem peptide model refinement**: `Don't refine`
>    - leave everything else as default
>
> 2. Click **Execute**.
>
{: .hands_on}

Once the database search is completed, the SearchGUI tool will output a file (called a
SearchGUI archive file) that will serve as an input for the next section, PeptideShaker.

> ### :nut_and_bolt: Comment
> In order to maintain accuracy and effectiveness in spectral / peptide / protein identification, a
> target-decoy search strategy can be used to discern how correct and incorrect a spectral or
> peptide or protein match is. The most popular approach for generating decoy databases is the
> *reverse database* approach. Essentially, protein sequences are reversed to generate a *decoy*
> database. Any matches and their associated scores against a target and decoy database are
> noted (with the premise that matches against decoy matches are incorrect). Later the matches
> are ranked according to descending scores and *decoy matches* are used to calculate false
> discovery rate (FDR) to set a threshold for valid identifications. The FDR approach allows for a
> fairer comparison of datasets across labs, machines and proteomic workflows. Please read
> the manuscript by Elias and Gygi (2010) for more information.
>
{: .comment}

#### PeptideShaker

**PeptideShaker** is a post-processing software tool that
processes data from the SearchGUI software tool. It serves to organize the Peptide-Spectral
Matches (PSMs) generated from SearchGUI processing and is contained in the SearchGUI archive.
It provides an assessment of confidence of the data, inferring proteins identified from the
matched peptide sequences and generates outputs that can be visualized by users to interpret
results. PeptideShaker has been wrapped in Galaxy to work in combination with SearchGUI
outputs.

> ### :nut_and_bolt: Comment
> There are a number of choices for different data files that can be generated using
> PeptideShaker. A compressed file can be made containing all information needed to view the
> results in the standalone PeptideShaker viewer. A `mzidentML` file can be created that contains
> all peptide sequence matching information and can be utilized by compatible downstream
> software. Other outputs are focused on the inferred proteins identified from the PSMs, as well
> as phosphorylation reports, relevant if a phosphoproteomics experiment has been undertaken.
{: .comment}

> ### :pencil2: Hands-on: PeptideShaker
>
> 1. **PeptideShaker** :wrench:: Run **PeptideShaker** with:
>   - **Compressed SearchGUI results**: The SearchGUI archive file
>   - **Specify Advanced PeptideShaker Processing Options**: `Default Processing Options`
>   - **Specify Advanced Filtering Options**: `Default Filtering Options`
>   - **Specify Contact Information for mzIdendML**: You can leave the default dummy options for now, but feel free to enter custom contact information.
>   - **Include the protein sequences in mzIdentML**: `No`
>   - **Output options**: Select the `PSM Report` (Peptide-Spectral Match) and the `Certificate of Analysis`
>
>       > ### :nut_and_bolt: Comment
>       > 
>       > The **Certificate of Analysis** provides details on all the parameters
>       > used by both SearchGUI and PeptideShaker in the analysis. This can be downloaded from the
>       > Galaxy instance to your local computer in a text file if desired.
>       {: .comment}
>
> 2. Click **Execute** and inspect the resulting files after they turned green with the **View data** icon:
>     ![](../../images/view_data_icon.png)
>
{: .hands_on}


A number of new items will appear in your history, each corresponding to the outputs selected
in the PeptideShaker parameters. Most relevant for this tutorial is the PSM report:

![](../../images/psm_report.png)

Scrolling at the bottom to the left will show the sequence for the PSM that matched to these
metapeptide entries. Column 3 is the sequence matched for each PSM entry. Every PSM is a
new row in the tabular output.

In the following steps of this tutorial, selected portions of this output will be extracted and used for
analysis of the taxonomic make-up of the sample as well as the biochemical functions
represented by the proteins identified.

## Unipept metaproteomic analysis

In the previous section, the genome sequencing and mass spectrometry data from
processing of biological samples was used to identify peptides present in those samples.
Now those peptides are used as evidence to infer which organisms are represented in the sample,
and what biological functions those peptides and associated proteins suggest are occurring.

The UniProt organization collects and annotates all known proteins for organisms. A UniProt
entry includes the protein amino acid sequence, the NCBI taxonomy, and any annotations
about structure and function of the protein. The UniPept web resource developed
by Ghent University will be used to match the sample peptides to proteins. UniPept indexes all Uniprot
proteins and provides a fast matching algorithm for peptides.

> ### :bulb: Tip: Unipept
>
> Users can access UniPept via a [web page](https://unipept.ugent.be) and paste peptide
> sequences into the search form to retrieve protein information. But we`ll use a Galaxy
> *Unipept* tool to automate the process. The *Unipept* tool sends the peptide list to the
> UniPept REST API service, then transforms the results into datasets that can be further analyzed
> or operated on within Galaxy.
{: .tip}

#### Recieving the list of peptides: Query Tabular

In order to use *Unipept*, a list containing the peptide sequences has to be generated.
The tool **Query Tabular** can load tabular data (the PSM report in this case) into a SQLite data base.
As a tabular file is being read, line filters may be applied and an SQL query can be performed.

> ### :pencil2: Hands-on: Query Tabular
>
> 1. **Query Tabular** :wrench:: Run **Query Tabular** with:
>
>    - **Database Table**: Click on `+ Insert Database Table`:
>    - **Tabular Dataset for Table**: The PSM report
>
>    Section **Filter Dataset Input**:
>   
>    - **Filter Tabular Input Lines**: Click on `+ Insert Filter Tabular Input Lines`:
>    - **Filter By**: Select `by regex expression matching`
>        - **regex pattern**: `^\d`
>        - **action for regex match**: `include line on pattern match`
>    
>    Section **Table Options**:
>
>    - **Specify Name for Table**: `psm`
>    - **Specify Column Names (comma-separated list)**: `id,,sequence,,,,,,,,,,,,,,,,,,,,confidence,validation`
>
>        > ### :nut_and_bolt: Comment
>        > 
>        > By default, table columns will be named: c1,c2,c3,...,cn (column names for a table must be unique).
>        > You can override the default names by entering a comma separated list of names, e.g. `,name1,,,name2`
>        > would rename the second and fifth columns.
>        >
>        > Check your input file to find the settings which best fits your needs.
>        {: .comment}
>
>    - **Only load the columns you have named into database**: `Yes`
>
>    - **Save the sqlite database in your history**: `Yes`
>
>        > ### :bulb: Tip
>        >
>        > * **Query Tabular** can also use an existing SQLite database. Activating `Save the sqlite database in your history`
>        > will store the created database in the history, allowing to reuse it directly. 
>        >
>        {: .tip}
>
>    - **SQL Query to generate tabular output**:
>         
>          SELECT distinct sequence
>         
>          FROM psm
>
>          WHERE validation IS NOT 'Confident' AND confidence >= 95
>          
>          ORDER BY sequence
>
>    > ### :question: Questions
>    >
>    > The SQL query might look confusing at first, but having a closer look should clarify a lot.
>    >
>    > 1. What does `FROM psm` mean?
>    > 2. What need to be changed if we only want peptides with a confidence higher then 98%?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>We want to read from table "psm". We defined the name before in the "Specify Name for Table" option.</li>
>    >    <li>We need to change the value in line 3: "WHERE validation IS NOT 'Confident' AND confidence >= 98"</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
>    - **Omit column headers from tabular output**: `Yes`
>
> 2. Click **Execute** and inspect the query results file after it turned green. If everything went well, it should look similiar:
>
>     ![](../../images/query_tabular_1.png)
>
{: .hands_on}

While we can proceed with this list of peptides, let's practice using the created SQLite database for further queries.
We might not only be interested in all the distinct peptides, but also on how many PSMs a single peptide had.
Therefore we can search the database for the peptides and count the occurrence without configuring the tables and columns again:

> ### :pencil2: Hands-on: SQLite to tabular
>
> 1. **SQLite to tabular** :wrench:: Run **SQLite to tabular** with:
>
>    - **SQL Query**:
>         
>          SELECT sequence as "peptide", count(id) as "PSMs" 
>         
>          FROM psm
>
>          WHERE validation IS NOT 'Confident' AND confidence >= 95
>          
>          GROUP BY sequence
>
>          ORDER BY sequence
>
> 2. Click **Execute**. The resulting file should have two columns, one with the distinct peptides, the other with the count number of PSMs.
>
{: .hands_on}


#### Retrieve taxonomy for peptides: Unipept

The generated list of peptides can now be used to search via *Unipept*.
We do a taxonomy analysis using the UniPept pept2lca function to return the taxonomic lowest common ancestor for each peptide:

> ### :pencil2: Hands-on: Unipept
>
> 1. **Unipept** :wrench:: Run **Unipept** with:
>
>    - **Unipept application**: `pept2lca: lowest common ancestor`
>    - **Peptides input format**: `tabular`
>    - **Tabular Input Containing Peptide column**: The query results file.
>    - **Select column with peptides**: `Column 1`
>    - **Choose outputs**: Select `tabular` and `JSON taxonomy tree`
>
> 2. Click **Execute**. The history should grow by two files. View each to see the difference.
>
>       > ### :nut_and_bolt: Comment
>       > 
>       > The JSON (JavaScript Object Notation) file contains the same information as the tabular file but is not comfortably human readable. 
>       > Instead, we can use it to use JavaScript libraries to visualize this data.
>       {: .comment}
>
> 3. Visualize the data:
>
>    - Click on the JSON output file from the *Unipept* tool to expand it. Click on the **Visualize** button and select **Unipept Tree viewer**:
>
>       ![](../../images/visualize_button.png)
>
>    - A new window should appear with a visualization of the taxonomy tree of your data. Use the mouse wheel to scroll in and out and click on nodes to expand or collapse them:
>
>       ![](../../images/unipept_tree_viewer.png)
>
{: .hands_on}

## Genus taxonomy level summary

The tabular *Unipept* output lists the taxonomy assignments for each peptide. To create a meaningful summary, the **Query Tabular** tool is
once again used, aggregating the number of peptides and PSMs for each genus level taxonomy assignment:

> ### :pencil2: Hands-on: Query Tabular
>
> 1. **Query Tabular** :wrench:: Run **Query Tabular** with:
>
>    - **Database Table**: Click on `+ Insert Database Table`
>    - **Tabular Dataset for Table**: The PSM report
>
>    Section **Filter Dataset Input**:
>   
>    - **Filter Tabular Input Lines**: Click on `+ Insert Filter Tabular Input Lines`:
>    - **Filter By**: Select `by regex expression matching`
>        - **regex pattern**: `^\d`
>        - **action for regex match**: `include line on pattern match`
>    
>    Section **Table Options**:
>
>    - **Specify Name for Table**: `psm`
>    - **Specify Column Names (comma-separated list)**: `,,sequence,,,,,,,,,,,,,,,,,,,,confidence,validation`
>
>    - **Only load the columns you have named into database**: `Yes`
>
> 2. Repeat this step to have a second **Database Table**:
>
>    - **Database Table**: Click on `+ Insert Database Table`
>    - **Tabular Dataset for Table**: The *Unipept* tabular/tsv output
>
>    Section **Filter Dataset Input**:
>   
>    - **Filter Tabular Input Lines**: Click on `+ Insert Filter Tabular Input Lines`:
>    - **Filter By**: Select `by regex expression matching`
>        - **regex pattern**: `#peptide`
>        - **action for regex match**: `exclude line on pattern match`
>    
>    Section **Table Options**:
>
>    - **Specify Name for Table**: `lca`
>    - **Specify Column Names (comma-separated list)**: `peptide,,,,,,,,,,,,,,,,,,,,,genus`
>
>    - **Only load the columns you have named into database**: `Yes`
>
>    - **Save the sqlite database in your history**: `No`
>
>    - **SQL Query to generate tabular output**:
>         
>          SELECT lca.genus,count(psm.sequence) as "PSMs",count(distinct psm.sequence) as "DISTINCT PEPTIDES" 
>         
>          FROM psm LEFT JOIN lca ON psm.sequence = lca.peptide 
>
>          WHERE validation IS NOT 'Confident' AND confidence >= 95
>          
>          GROUP BY lca.genus
>
>          ORDER BY PSMs desc, 'DISTINCT PEPTIDES' desc
>
>    - **Omit column headers from tabular output**: `No`
>
> 2. Click **Execute** and inspect the query results file after it turned green:
>
>     ![](../../images/metaproteomics_summary.png)
>
{: .hands_on}