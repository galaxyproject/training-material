---
layout: tutorial_hands_on

title: "Metaproteomics tutorial"
edam_ontology: "topic_0121"
zenodo_link: "https://doi.org/10.5281/zenodo.839701"
questions:
  - "How can I match metaproteomic mass spectrometry data to peptide sequences derived from shotgun metagenomic data?"
  - "How can I perform taxonomy analysis and visualize metaproteomics data?"
  - "How can I perform functional analysis on this metaproteomics data?"
objectives:
  - "A taxonomy and functional analysis of metaproteomic mass spectrometry data."
time_estimation: "2h"
key_points:
  - "Use dataset collections"
  - "With SearchGUI and PeptideShaker you can gain access to multiple search engines"
  - "Learning the basics of SQL queries can pay off"
contributors:
  - subinamehta
  - timothygriffin
  - pratikdjagtap 
  
  
---

# Introduction
{:.no_toc}

In this metaproteomics tutorial we will identify expressed proteins from a complex bacterial community sample.
For this MS/MS data will be matched to peptide sequences provided through a FASTA file.

Metaproteomics is the large-scale characterization of the entire protein complement of environmental microbiota
at a given point in time. It has the potential to unravel the mechanistic details of microbial interactions with
the host / environment by analyzing the functional dynamics of the microbiome.

In this tutorial, we will analyze a sample of sea water that was collected in August of 2013 from the Bering
Strait chlorophyll maximum layer (7m depth, 65° 43.44″ N, 168° 57.42″ W). The data were originally published in [May et al., 2016](https://www.ncbi.nlm.nih.gov/pubmed/27396978).

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

There are three ways to upload your data.

*   Upload/Import the files from your computer
*   Using a direct link
*   Import from the data library if your instance provides the files

In this tutorial, we will get the data from Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.839701.svg)](https://doi.org/10.5281/zenodo.839701).

> ### {% icon hands_on %} Hands-on: Data upload and organization
>
> 1. Create a new history and name it something meaningful (e.g. *Metaproteomics tutorial*)
>
>    {% include snippets/create_new_history.md %}
>    {% include snippets/rename_history.md %}
>
> 2. Import the three MGF MS/MS files and the FASTA sequence file from Zenodo.
>
>    {% include snippets/import_via_link.md %}
>
>    As default, Galaxy takes the link as name.
>
>    > ### {% icon comment %} Comments
>    > - Rename the datasets to a more descriptive name
>    > - There is a GO term file in the zenodo folder is for reference purposes, please do not load it on your account.
>    {: .comment}
>
> 3. Build a **Dataset list** for the three MGF files
>    - Click the **Operations on multiple datasets** check box at the top of the history panel
>       ![Operations on multiple datasets button](../../images/operations_icon.png)
>    - Check the three boxes next to the MGF files
>    - Click **For all selected...** and choose **Build dataset list**
>    - Ensure the three control samples are the only ones selected, and enter a name for the new collection (e.g. *MGF files*)
>    - Click **Create list** and exit by clicking again the dataset operations icon
>
{: .hands_on}

# Analysis

## Match peptide sequences

The search database labelled `FASTA_Bering_Strait_Trimmed_metapeptides_cRAP.FASTA` is the input database that
will be used to match MS/MS to peptide sequences via a sequence database search. It is a small excerpt of the original database, which was constructed based on a metagenomic screening of the sea water samples (see [May et al. (2016)](https://www.ncbi.nlm.nih.gov/pubmed/27396978)). The full original database can be accessed from [here](https://noble.gs.washington.edu/proj/metapeptide/data/metapeptides_BSt.fasta). The contaminant database (cRAP) was merged with the original database.

For this, the sequence database-searching program called [SearchGUI](https://compomics.github.io/projects/searchgui.html) will be used.
The created dataset collection of the three *MGF files* in the history is used as the MS/MS input.

#### SearchGUI

> ### {% icon hands_on %} Hands-on: SearchGUI
>
> 1. **SearchGUI** {% icon tool %}: Run **SearchGUI** with:
>    - **Protein Database**: `FASTA_Bering_Strait_Trimmed_metapeptides_cRAP.FASTA`(or however you named the `FASTA` file)
>    - **Input Peak lists (mgf)**: `MGF files` dataset collection.
>
>    > ### {% icon tip %} Tip: Select dataset collections as input
>    >
>    > * Click the **Dataset collection** icon on the left of the input field:
>    >
>    >      ![Dataset collection button](../../images/dataset_button.png)
>    > * Select the appropriate dataset collection from the list
>    {: .tip}
>
>    Section **Search Engine Options**:
>
>    - **Search Engines**: `X!Tandem`
>
>   - {% icon param-file %} *”Protein database”*: `Output dataset ‘output’ from step 1`
>   - {% icon param-check %} *”Create a concatenated target/decoy before running PeptideShaker”*: `Yes`
>   - {% icon param-check %} *”Gene mappings will be used and salves along with the project (UniProt databases only)”* `No`
>   - {% icon param-check %} *”Update gene mappings automatically from Esembl (UniProt databases only)”*: `No`
>   - {% icon param-file %} *”Input Peak Lists (mgf)”*: `Output dataset ‘output’ from step 2`
>   - {% icon param-check %} *”DB-Search engines”*:`X!Tandem`
>   - {% icon param-select %} *”Digestion”*: `Select enzymes`
>   - {% icon param-select %} *”Enzyme”* `Trypsin`
>   - {% icon param-text %} *”Maximum Missed Cleavages”*: `2`
>   - {% icon param-select %} *”Precursor Ion Tolerance Units”*: `Parts per million (ppm)`
>   - {% icon param-text %} *”Precursor Ion Tolerance”*: `10.0`
>   - {% icon param-select %} *”Fragment Tolerance Units”*: `Daltons`
>   - {% icon param-text %} *” Fragment Tolerance”*: `0.02`
>   - {% icon param-text %} *”Minimum Charge”*: `2`
>   - {% icon param-text %} *”Maximum Charge”*: `6`
>   - {% icon param-select %} *”Forward Ion”*: `b`
>   - {% icon param-text %} *”Reverse Ion”*: `y`
>   - {% icon param-text %} *”Minimum precursor isotope”*: `0`
>   - {% icon param-text %} *”Maximum precursor isotope”*: `1`
>   - {% icon param-select %} *”Fixed Modifications”*: `Carbamidomethylation of C`
>   - {% icon param-select %} *”Variable Modifications”*: `Oxidation of M`
>
>    > ### {% icon tip %} Tip: Search for options
>    >
>    > * For selection lists, typing the first few letters in the window will filter the available options.
>    {: .tip}
>
>   - {% icon param-select %} *”Search GUI Options”*: `Default`
>   - {% icon param-select %} *”X!Tandem Options”*: `Advanced`
>   - {% icon param-text %} *”X!Tandem: Total Peaks”*: `50`
>   - {% icon param-text %} *”X!Tandem: Min Peaks”*: `15`
>   - {% icon param-text %} *”X!Tandem: Min Frag m/z”*: `200`
>   - {% icon param-text %} *”X!Tandem: Min Precursor Mass”*: `200`
>   - {% icon param-check %} *”X!Tandem: Noise Suppression”*: `Yes`
>   - {% icon param-text %} *”X!Tandem: Dynamic Range”*: `100`
>   - {% icon param-check %} *”X!Tandem: Quick Acetyl”*: `No`
>   - {% icon param-check %} *”X!Tandem: Quick Pyrolidone”*: `No`
>   - {% icon param-check %} *”X!Tandem: Protein stP Bias”*: `No`
>   - {% icon param-text %} *”X!Tandem: Maximum Valid Expectation Value”*: `100.0`
>   - {% icon param-check %} *” X!Tandem: Output proteins”*: `No`
>   - {% icon param-check %} *”X!Tandem: Output sequences”*: `No`
>   - {% icon param-check %} *”X!Tandem: Output Spectra”*: `Yes`
>   - {% icon param-select %} *”X!Tandem peptide model refinement”*: `Don’t refine`
>   - leave everything else as default
>
> 2. Click **Execute**.
>
{: .hands_on}

Once the database search is completed, the SearchGUI tool will output a file (called a
SearchGUI archive file) that will serve as an input for the next section, PeptideShaker.

> ### {% icon comment %} Comment
> Note that sequence databases used for metaproteomics are usually much larger than the excerpt used in this tutorial. When using large databases, the peptide identification step can take much more time for computation. In metaproteomics, choosing the optimal database is a crucial step of your workflow, for further reading see [Timmins-Schiffman et al (2017)](https://www.ncbi.nlm.nih.gov/pubmed/27824341).
>
> To learn more about database construction in general, like integrating contaminant databases or using a decoy strategy for FDR searching, please consult our tutorial on [Database Handling]({{site.baseurl}}/topics/proteomics/tutorials/database-handling/tutorial.html).
>
{: .comment}

#### PeptideShaker

[PeptideShaker](https://compomics.github.io/projects/peptide-shaker.html) is a post-processing software tool that
processes data from the SearchGUI software tool. It serves to organize the Peptide-Spectral
Matches (PSMs) generated from SearchGUI processing and is contained in the SearchGUI archive.
It provides an assessment of confidence of the data, inferring proteins identified from the
matched peptide sequences and generates outputs that can be visualized by users to interpret
results. PeptideShaker has been wrapped in Galaxy to work in combination with SearchGUI
outputs.

> ### {% icon comment %} Comment
> There are a number of choices for different data files that can be generated using
> PeptideShaker. A compressed file can be made containing all information needed to view the
> results in the standalone PeptideShaker viewer. A `mzidentML` file can be created that contains
> all peptide sequence matching information and can be utilized by compatible downstream
> software. Other outputs are focused on the inferred proteins identified from the PSMs, as well
> as phosphorylation reports, relevant if a phosphoproteomics experiment has been undertaken.
> More detailed information on peptide inference using SearchGUI and PeptideShaker can be found in our tutorial on [Peptide and Protein ID]({{site.baseurl}}/topics/proteomics/tutorials/protein-id-sg-ps/tutorial.html).
{: .comment}

> ### {% icon hands_on %} Hands-on: PeptideShaker
>
> 1. **PeptideShaker** {% icon tool %}: Run **PeptideShaker** with:
>   - {% icon param-file %} *”Compressed SearchGUI results”*: `Output dataset ‘searchgui_results’ from step 3`
>   - {% icon param-select %} *”Specify Advanced Peptide Shaker Processing Options”*: `Advanced Processing Options`
>   - {% icon param-text %} *”FDR at the protein level”*: `1.0`
>   - {% icon param-text %} *”FDR at the peptide level”*: `1.0`
>   - {% icon param-text %} *”FDR at the PSM level”*: `1.0`
>   - {% icon param-text %} *”Minimum confidence required for a protein in the fraction MW plot”*: `95.0`
>   - {% icon param-select %} *”The PRM probabilistic score to use for PTM localization”*: `A-score`
>   - {% icon param-select %} *”The PTM to peptide sequence matching type”*: `Amino Acids`
>   - {% icon param-check %} *”Align peptide ambiguously localizes PRMs on confident sites”*: `Yes`
>   - {% icon param-select %} *”Specify Advanced FIltering Options”*: `Advanced Filtering Options`
>   - {% icon param-text %} *”Minimum Peptide Length”*: `6`
>   - {% icon param-text %} *”Maximum Peptide Length”*: `65`
>   - {% icon param-text %} *”Maximum Precursor Error”*: `10.0`
>   - {% icon param-select %} *”Maximum Precursor Error Type”*: `ppm`
>   - {% icon param-check %} *”Exclude Unknown PTMs”: `Yes`
>   - {% icon param-check %} *”Creates a mzldentML file*”: `No`
>   - {% icon param-check %} *”Compress results into single zip file”*: `No`
>   - {% icon param-check %} *”Exports the CPS file”*: `No`
>   - {% icon param-check %} *”Reports to be generated”*: `Select the ‘PSM Report’`
>
>
> 2. Click **Execute** and inspect the resulting files after they turned green with the **View data** icon:
>     ![View data button](../../images/view_data_icon.png)
>
{: .hands_on}


A number of new items will appear in your history, each corresponding to the outputs selected
in the PeptideShaker parameters. Most relevant for this tutorial is the PSM report:

![Display of the PSM report tabular file](../../images/psm_report.png "The PSM report")

Scrolling towards left will show the sequence for the PSM that matched to these
metapeptide entries. Column 3 is the sequence matched for each PSM entry. Every identified PSM is a
new row in the tabular output.

In the following steps of this tutorial, selected portions of this output will be extracted and used for
analysis of the taxonomic make-up of the sample as well as the biochemical functions
represented by the proteins identified.

## Taxonomy analysis

In the previous section, the genome sequencing and mass spectrometry data from
processing of biological samples was used to identify peptides present in those samples.
Now those peptides are used as evidence to infer which organisms are represented in the sample,
and what biological functions those peptides and associated proteins suggest are occurring.

The UniProt organization collects and annotates all known proteins for organisms. A UniProt
entry includes the protein amino acid sequence, the NCBI taxonomy, and any annotations
about structure and function of the protein. The UniPept web resource developed
by Ghent University will be used to match the sample peptides to proteins. UniPept indexes all Uniprot
proteins and provides a fast matching algorithm for peptides.

> ### {% icon comment %} Unipept
>
> Users can access UniPept via a [web page](https://unipept.ugent.be) and paste peptide
> sequences into the search form to retrieve protein information. But we'll use the Galaxy
> *Unipept* tool to automate the process. The *Unipept* tool sends the peptide list to the
> UniPept REST API service, then transforms the results into datasets that can be further analyzed
> or operated on within Galaxy.
{: .comment}

#### Recieving the list of peptides: Query Tabular

In order to use *Unipept*, a list containing the peptide sequences has to be generated.
The tool **Query Tabular** can load tabular data (the PSM report in this case) into a SQLite data base.
As a tabular file is being read, line filters may be applied and an SQL query can be performed.

> ### {% icon hands_on %} Hands-on: Query Tabular
>
> 1. **Query Tabular** {% icon tool %}: Run **Query Tabular** with:
>
>   - {% icon param-file %} *”Tabular Dataset for Table”*: `Output dataset ‘output_psm’ from step 4`
>   - {% icon param-text %} *”Specify Name for Table”*: `psm`
>   - {% icon param-check %} *”Use first line as column names”*: `No`
>   - {% icon param-text %} *”Specify Column Names (comma-separated list)”*: `id,,sequence,,,,,,,,,,,,,,,,,,,,confidence,validation`
>   - {% icon param-check %} *”Only load the columns you have named into database”*: `Yes`
>   - {% icon param-check %} *”Save the sqlite database in your history”*: `No`
>
>        > ### {% icon comment %} Comment
>        >
>        > By default, table columns will be named: c1,c2,c3,...,cn (column names for a table must be unique).
>        > You can override the default names by entering a comma separated list of names, e.g. `,name1,,,name2`
>        > would rename the second and fifth columns.
>        >
>        > Check your input file to find the settings which best fits your needs.
>        {: .comment}
>
>        > ### {% icon comment %} Querying SQLite Databases
>        >
>        > * **Query Tabular** can also use an existing SQLite database. Activating `Save the sqlite database in your history`
>        > will store the created database in the history, allowing to reuse it directly.
>        >
>        {: .comment}
>
>    - **SQL Query to generate tabular output**:
>
>          `SELECT distinct sequence
>
>          FROM psm
>
>          WHERE confidence >= 95
>
>          ORDER BY sequence`
>
>    > ### {% icon question %} Questions
>    >
>    > The SQL query might look confusing at first, but having a closer look should clarify a lot.
>    >
>    > 1. What does `FROM psm` mean?
>    > 2. What need to be changed if we only want peptides with a confidence higher then 98%?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. We want to read from table "psm". We defined the name before in the "Specify Name for Table" option.
>    > > 2. We need to change the value in line 3: "WHERE validation IS NOT 'Confident' AND confidence >= 98"
>    > {: .solution }
>    {: .question}
>
>    - **include query result column headers**: `No`
>
> 2. Click **Execute** and inspect the query results file after it turned green. If everything went well, it should look similiar:
>
>     ![Query Tabular output showing the peptides](../../images/query_tabular_1.png "Query Tabular output")
>
{: .hands_on}


#### Retrieve taxonomy for peptides: Unipept

The generated list of peptides can now be used to search via *Unipept*.
We do a taxonomy analysis using the UniPept peptinfo function to return the taxonomic lowest common ancestor for each peptide:

> ### {% icon hands_on %} Hands-on: Unipept
>
> 1. **Unipept** {% icon tool %}: Run **Unipept** with:
>   - {% icon param-select %} *”Unipept application”*: `peptinfo: Tryptic peptides and associated EC and GO terms and lowest common ancestor taxonomy`
>   - {% icon param-check %} *”Equate isoleucine and leucine”*: `No`
>   - {% icon param-check %} *”Retrieve extra information*”: `No`
>   - {% icon param-check %} *”Group responses by GO namespace (biological process, molecular function, cellular component)”*: `Yes`
>   - {% icon param-check %} *”Names”*: `No`
>   - {% icon param-check %} *”Allfields”*: `No`
>   - {% icon param-select %} *”Peptides input format”*: `Tabular`
>   - {% icon param-file %} *”Tabular Input Containing Peptide column”*: `Output dataset ‘output’ from step 5`
>   - {% icon param-text %} *”Select column with peptides”*: `1`
>   - {% icon param-check %} *”Choose outputs”*: `Select the ‘JSON Taxonomy Tree (for pepet2lca, pe2taxa, and peptinfo),’ ‘Peptide GO terms in normalized tabular (for pept2go, pept2funct, and peptinfo),’ ‘Peptide EC terms in normalized tabular (for pept2ec, pept2funct, and peptinfo),’ ‘JSON EC Coverage Tree (for pept2ec, pep2funct, and peptinfo)’`
>   - {% icon param-check %} *”Exit with error on invalid peptides, otherwise ignore them”*: `No`
>
> 2. Click **Execute**. The history should grow by two files. View each to see the difference.
>
>       > ### {% icon comment %} Comment
>       >
>       > The JSON (JavaScript Object Notation) file contains the same information as the tabular file but is not comfortably human readable.
>       > Instead, we can use it to use JavaScript libraries to visualize this data.
>       {: .comment}
>
> 3. Visualize the data:
>
>    - Click on the JSON output file from the *Unipept* tool to expand it. Click on the **Visualize** button and select **Unipept Tree viewer**:
>
>       ![Visualize button](../../images/Visualize_output.png)
>
>
>       ![Viewer](../../images/Unipept_viewer.png)
>
>
>    - A new window should appear with a visualization of the taxonomy tree of your data. Use the mouse wheel to scroll in and out and click on nodes to expand or collapse them:
>
>       ![Unipept Tree viewer visual output](../../images/Treeview_Unipet.png "Interactive Tree viewer visualization from the Unipept Tree viewer plugin")
>
>       ![Unipept Sunburst visual output](../../images/sunburst_unipept.png "Interactive Sunburst visualization from the Unipept Tree viewer plugin")
>
>       ![Unipept Tree map visual output](../../images/Treemap_unipept.png "Interactive Tree map visualization from the Unipept Tree viewer plugin")

> The user can perform the same functions for viewing the EC output.
> ![Visualize button](../../images/EC-viewer.png)
>
> The Unipept viewer also provides interactive viewer for EC number (Tree view, Sunburst and Treemap)
>
> ![Unipept Tree viewer visual output](../../images/EC_treemap.png "Interactive Tree viewer visualization from the Unipept Tree viever plugin")
{: .hands_on}

## Genus taxonomy level summary

The tabular *Unipept* output lists the taxonomy assignments for each peptide. To create a meaningful summary, the **Query Tabular** tool is
once again used, aggregating the number of peptides and PSMs for each genus level taxonomy assignment:

> ### {% icon hands_on %} Hands-on: Query Tabular
>
> 1. **Query Tabular** {% icon tool %}: Run **Query Tabular** with:
>   - {% icon param-file %} *”Add tables to this Database”*: ‘‘
>   - {% icon param-file %} *”Tabular Dataset for Table”*: `Output dataset ‘output_psm’ from step 4`
>   - {% icon param-select %} *”Filter By”*: `by regex expression matching`
>   - {% icon param-text %} *”Regex pattern”*: `^\d`
>   - {% icon param-select %} *”Action for regex match”*: `include line on pattern match`
>   - {% icon param-text %} *”Specify Name for Table”*: `psm`
>   - {% icon param-check %} *”Use first line as column names”*: `No`
>   - {% icon param-text %} *”Specify Column Names (comma-separated list)”*: `,,sequence,,,,,,,,,,,,,,,,,,,,confidence,validation`
>   - {% icon param-check %} *”Only load the columns you have named into database”*: `Yes`
>   - {% icon param-text %} *”Add an auto increment primary key column with this name”*: ``
>   - {% icon param-file %} *”Tabular Dataset for Table”*: `Output dataset ‘output_tsv’ from step 7`
>   - {% icon param-select %} *”Filter By”*: `by regex expression matching`
>   - {% icon param-text %} *”Regex pattern”*: `#peptide`
>   - {% icon param-select %} *”Action for regex match”*: `exclude line on pattern match`
>   - {% icon param-text %} *”Specify Name for Table”*: `lca`
>   - {% icon param-check %} *”Use first line as column names”*: `No`
>   - {% icon param-text %} *”Specify Column Names (comma-separated list)”*: `peptide,,,,,,,,,,,,,,,,,,,,,genus`
>   - {% icon param-check %} *”Only load the columns you have named into database”*: `Yes`
>   - {% icon param-text %} *”Add an auto increment primary key column with this name”*: ``
>   - {% icon param-check %} *”Save the sqlite database in your history”*: `Yes`
>   - {% icon param-text %} *”SQL Query to generate tabular output”*:
 `SELECT lca.genus,count(psm.sequence) as "PSMs",count(distinct psm.sequence) as "DISTINCT PEPTIDES"
FROM psm LEFT JOIN lca ON psm.sequence = lca.peptide
WHERE confidence >= 95
GROUP BY lca.genus
ORDER BY PSMs desc, 'DISTINCT PEPTIDES' desc`
>
>   - {% icon param-check %} *”Include query result column headers”*: `Yes`
>   - {% icon param-select %} *”Prefix character for column_header line”*: ‘#’
>
>
> 2. Click **Execute** and inspect the query results file after it turned green:
>
>     ![Query Tabular output showing gene, PSMs and distinct peptides](../../images/Genera_PSM.png "Query Tabular output")
>
{: .hands_on}

## Functional Analysis

Recent advances in microbiome research indicate that functional characterization via metaproteomics analysis has the potential to accurately measure the microbial response to perturbations. In particular, metaproteomics enables the estimation of the function of the microbial community based on expressed microbial proteome.

In the following chapter, a functional analysis will be performed using the **UniPept** application `peptinfo` in order to match the list of peptides with the correlated Gene Ontology terms.
This allows to get an insight of the **biological process**, the **molecular function** and the **cellular component** related to the sample data. We also performed EC number analysis on this dataset

> ### {% icon comment %} Gene Ontology (GO) Consortium
>
> The [Gene Ontology Consortium](http://www.geneontology.org/) provides with its Ontology a framework for the model of biology.
> The GO defines concepts/classes used to describe gene function, and relationships between these concepts. It classifies functions along three aspects:
>
>
> - **molecular function**
>
>   - molecular activities of gene products
>
> - **cellular component**
>
>   - where gene products are active
>
> - **biological process**
>
>   - pathways and larger processes made up of the activities of multiple gene products.
>
> [more information](http://geneontology.org/page/ontology-documentation)
>
{: .comment}


#### Retrieve GO and EC IDs for peptides: Unipept

The **UniPept** application `peptinfo` can be used to return the list of proteins containing each peptide.The option `retrieve extra information` option is set to `yes` so that we retrieve Gene Ontology assignments for each protein. Unipept 4.0 has inbuilt GO terms.

> ### {% icon hands_on %} Hands-on: Unipept
>
> 1. **Unipept** {% icon tool %}: Run **Unipept** with:
>   - {% icon param-select %} *”Unipept application”*: ‘peptinfo: Tryptic peptides and associated EC and GO terms and lowest common ancestor taxonomy’
>   - {% icon param-check %} *”Equate isoleucine and leucine”*: `No`
>   - {% icon param-check %} *”Retrieve extra information*”: `No`
>   - {% icon param-check %} *”Group responses by GO namespace (biological process, molecular function, cellular component)”*: `Yes`
>   - {% icon param-check %} *”Names”*: `No`
>   - {% icon param-check %} *”Allfields”*: `No`
>   - {% icon param-select %} *”Peptides input format”*: `Tabular`
>   - {% icon param-file %} *”Tabular Input Containing Peptide column”*: `Output dataset ‘output’ from step 5`
>   - {% icon param-text %} *”Select column with peptides”*: `1`
>   - {% icon param-check %} *”Choose outputs”*: `Select the ‘Tabular with one line per peptide,’ ‘JSON Taxonomy Tree (for pepet2lca, pe2taxa, and peptinfo),’ ‘Peptide GO terms in normalized tabular (for pept2go, pept2funct, and peptinfo),’ ‘JSON EC Coverage Tree (for pept2ec, pep2funct, and peptinfo)`
>   - {% icon param-check %} *”Exit with error on invalid peptides, otherwise ignore them”*: `No`
> 2. Click **Execute**.
>
> 3. inspect the result:
>
>    - The output should be a tabular file containing a column labeled `go_references`. This is what we're looking for.
>
{: .hands_on}


#### Combine all information to quantify the EC and GO results

As a final step we will use **Query Tabular** in a more sophisticated way to combine all information to quantify the EC analysis. 
> ### {% icon hands_on %} Hands-on: Query Tabular
>
>1.**Query Tabular** {% icon tool %} with the following parameters:
>   - {% icon param-file %} *”Add tables to this Database”*: ``
>   - {% icon param-file %} *”Tabular Dataset for Table”*: `Output dataset ‘output_psm’ from step 4`
>   - {% icon param-select %} *”Filter By”*: `by regex expression matching`
>   - {% icon param-text %} *”Regex pattern”*: `^\d`
>   - {% icon param-select %} *”Action for regex match”*: `include line on pattern match`
>   - {% icon param-text %} *”Specify Name for Table”*: `psm`
>   - {% icon param-check %} *”Use first line as column names”*: `No`
>   - {% icon param-text %} *”Specify Column Names (comma-separated list)”*: `,,sequence,,,,,,,,,,,,,,,,,,,,confidence,validation`
>   - {% icon param-check %} *”Only load the columns you have named into database”*: `Yes`
>   - {% icon param-text %} *”Add an auto increment primary key column with this name”*: ``
>   - {% icon param-file %} *”Tabular Dataset for Table”*: `Output dataset ‘output_ec_tsv’ from step 6`
>   - {% icon param-select %} *”Filter By”*: `by regex expression matching`
>   - {% icon param-text %} *”Regex pattern”*: `#peptide`
>   - {% icon param-select %} *”Action for regex match”*: `exclude line on pattern match`
>   - {% icon param-text %} *”Specify Name for Table”*: `goec`
>   - {% icon param-check %} *”Use first line as column names”*: `No`
>   - {% icon param-text %} *”Specify Column Names (comma-separated list)”*: `peptide,total_protein_count,ec_number,protein_count`
>   - {% icon param-check %} *”Save the sqlite database in your history”*: `No`
>   - {% icon param-text %} *”SQL Query to generate tabular output”*: 
`SELECT goec.ec_number,count(psm.sequence) as "PSMs",count(distinct psm.sequence) as "DISTINCT PEPTIDES" 
FROM psm LEFT JOIN goec ON psm.sequence = goec.peptide 
WHERE confidence >= 95 
GROUP BY goec.ec_number 
ORDER BY PSMs desc, 'DISTINCT PEPTIDES' desc`
>   - {% icon param-check %} *”Include query result column headers”*: `Yes`
>   - {% icon param-select %} *”Prefix character for column_header line”*: `#`
>
> 2. Click **Execute**.
>
{: .hands_on}


Here is another query tabular on extracting all the GO terms from the Unipept results

> ### {% icon hands_on %} Hands-on: Query Tabular
>
> 1. **Query Tabular** {% icon tool %}: Run **Query Tabular** with:
>
>   - {% icon param-file %} *”Add tables to this Database”*: ``
>   - {% icon param-file %} *”Tabular Dataset for Table”*: ‘Output dataset `output_psm’ from step 4`
>   - {% icon param-select %} *”Filter By”*: `by regex expression matching`
>   - {% icon param-text %} *”Regex pattern”*: `^\d`
>   - {% icon param-select %} *”Action for regex match”*: `include line on pattern match`
>   - {% icon param-text %} *”Specify Name for Table”*: `psm`
>   - {% icon param-check %} *”Use first line as column names”*: `No`
>   - {% icon param-text %} *”Specify Column Names (comma-separated list)”*: `,,sequence,,,,,,,,,,,,,,,,,,,,confidence,validation`
>   - {% icon param-check %} *”Only load the columns you have named into database”*: `Yes`
>   - {% icon param-text %} *”Add an auto increment primary key column with this name”*: ``
>   - {% icon param-file %} *”Tabular Dataset for Table”*: `Output dataset ‘output_go_tsv’ from step 6`
>   - {% icon param-select %} *”Filter By”*: `by regex expression matching`
>   - {% icon param-text %} *”Regex pattern”*: `#peptide`
>   - {% icon param-select %} *”Action for regex match”*: `exclude line on pattern match`
>   - {% icon param-text %} *”Specify Name for Table”*: `goterm`
>   - {% icon param-check %} *”Use first line as column names”*: `No`
>   - {% icon param-text %} *”Specify Column Names (comma-separated list)”*: `peptide,total_protein_count,go_term,protein_count,go_name`
>   - {% icon param-check %} *”Only load the columns you have named into database”*: `Yes`
>   - {% icon param-text %} *”Add an auto increment primary key column with this name”*: ``
>   - {% icon param-check %} *”Save the sqlite database in your history”*: `Yes`
>   - {% icon param-text %} *”SQL Query to generate tabular output”*: 
`SELECT goterm.go_term,count(psm.sequence) as "PSMs",count(distinct psm.sequence) as "DISTINCT PEPTIDES",goterm.go_name 
FROM psm LEFT JOIN goterm ON psm.sequence = goterm.peptide 
WHERE confidence >= 95 
GROUP BY goterm.go_term 
ORDER BY PSMs desc, 'DISTINCT PEPTIDES' desc`
>   - {% icon param-check %} *”Include query result column headers”*: `Yes`
>   - {% icon param-select %} *”Prefix character for column_header line”*: `#`
>
>
> 2. Click **Execute**.
>
{: .hands_on}

The next three steps are to filter out the three different Go terms. For that we use the Filter data on any column using simple expressions tool and extract the GO terms and the corresponding number of peptides associated with these terms.

> ### {% icon hands_on %} Hands-on: Filter data on any column using simple expressions
>
> **Filter** {% icon tool %} with the following parameters: 
>   - {% icon param-file %} *”Filter”*: `Output dataset ‘output’ from step 9`
>   - {% icon param-text %} *”With following condition”*: `c4=='biological process'`
>   - {% icon param-text %} *”Number of header lines to skip”*: `1`
>
> 2. Click **Execute**.
>
{: .hands_on}
> ![Filter output](../../images/Biological_process.png)

> ### {% icon hands_on %} Hands-on: Filter data on any column using simple expressions
>
> **Filter** {% icon tool %} with the following parameters: 
>   - {% icon param-file %} *”Filter”*: `Output dataset ‘output’ from step 9`
>   - {% icon param-text %} *”With following condition”*: `c4=='cellular component'`
>   - {% icon param-text %} *”Number of header lines to skip”*: `1`
>
> 2. Click **Execute**.
>
{: .hands_on}
> ![Filter output](../../images/Cellular_component.png)
> ### {% icon hands_on %} Hands-on: Filter data on any column using simple expressions
>
> **Filter** {% icon tool %} with the following parameters: 
>   - {% icon param-file %} *”Filter”*: `Output dataset ‘output’ from step 9`
>   - {% icon param-text %} *”With following condition”*: `c4==’molecular function’`
>   - {% icon param-text %} *”Number of header lines to skip”*: `1`
>
> 2. Click **Execute**.
>
{: .hands_on}
> ![Filter output](../../images/molecular_function.png)

With these three output files the functional analysis of this tutorial is finished. Each record contains the name of a GO term, the amount of peptides related to it and the amount of PSMs for these peptides.

This marks the end to the metaproteomics workflow! 

> ### {% icon comment %} References
>
> - [Dataset](https://www.ncbi.nlm.nih.gov/pubmed/27824341) and [SixGill software](https://www.ncbi.nlm.nih.gov/pubmed/27396978)
>
> - [Galaxy workflows for metaproteomics](https://www.ncbi.nlm.nih.gov/pubmed/26058579)
>
> - [Metaproteomics community effort](https://z.umn.edu/gcc2017mporal)
>
> - [Unipept](https://www.ncbi.nlm.nih.gov/pubmed/28552653)
>
> - [Galaxy-P Metaproteomics instance](http://z.umn.edu/metaproteomicsgateway)
>
> - [Metaproteomics video](http://z.umn.edu/mpvideo2018)
{: .comment}

