---
layout: tutorial_hands_on

title: Annotating a protein list identified by LC-MS/MS experiments

zenodo_link: 
- 'https://zenodo.org/record/2650868'
- 'https://zenodo.org/record/2650874'
- 'https://zenodo.org/record/2650872'

questions:
- How to filter out technical contaminants?
- How to check for tissue-specificity?
- How to perform enrichment analysis?
- How to map your protein list to pathways (Reactome)?
- How to compare your proteome with other studies?
objectives:
- Execute a complete annotation pipeline of a protein list identified by LC-MS/MS experiments
time_estimation: 60 minutes
key_points:

contributors:
- vloux
- combesf
- davidchristiany
- yvandenb

---

# Introduction
{:.no_toc}

[ProteoRE Galaxy instance](http://www.proteore.org) provides necessary tools to execute a whole annotation pipeline of a protein list identified by LC-MS/MS experiments. This activity introduces these tools and guides you through a simple pipeline using some example datasets based on the following study: [Proteomic characterization of human exhaled breath condensate](https://www.ncbi.nlm.nih.gov/pubmed/29189203) by Lacombe *et al., European Journal of Breath, 2018*. 


Once identified and/or quantified using a MS-based approach, interpreting the proteome in a sample is an important step to characterize its content in terms of functional properties in order to extend the biological knowledge related to this sample. In this activity, we illustrate the annotation and the exploration of the human exhaled breath condensate (EBC) proteome by performing the following steps:
> ### Agenda
>
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get Input Datasets

For this tutorial, we will use 3 datasets: the list of proteins identified by LC-MS/MS in the exhaled breath condensate (EBC) from Lacombe *et al.* and two others EBC proteomes previously published.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial and give it a name
>   {% include snippets/create_new_history.md %}
>
> 2. Import the files from [Zenodo](https://zenodo.org) or from the shared data library (ask your instructors). 
>
>   The datasets are available on Zenodo under the references: [2650868](https://zenodo.org/record/2650868) for Lacombe et
>   al., [2650874](https://zenodo.org/record/2650874) for Mucilli and [2650872](https://zenodo.org/record/2650872) for
>   Bredberg. 
>
> You can import the 3 data files directly with these links in the galaxy upload menu :
> 
> https://zenodo.org/record/2650868/files/Lacombe_et_al_2017_OK.txt
> 
> https://zenodo.org/record/2650874/files/Mucilli.txt
> 
> https://zenodo.org/record/2650872/files/Bredberg.txt
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
>
{: .hands_on}


# Filtering out technical contaminants

A group of 10 proteins were identified in both “technical” control samples with an enrichment in EBC samples below a fixed threshold. These proteins were thus considered to be technical contaminants (see list of proteins in Table 4 in [_Lacombe et al. 2018_](https://www.ncbi.nlm.nih.gov/pubmed/29189203)) and have to be removed from the initial dataset.

> ### {% icon hands_on %} Hands-on: Remove the contaminants
>
> 1. **Filter by keywords and/or numerical value** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `Lacombe_et_al_2017.txt` (Input dataset)
>    - Keep default option `Yes` for *"Does file contain header?"*
>    - *"Operation"*: `Discard`. We want to remove technical contaminants.
>    - Keep the `OR` option (by default) for the **operator** parameter. We don't need that parameter for a single filter.
>    - First add a  {% icon param-repeat %} *"Insert Filter by keywords"* box with a list of keywords to be filtered out.   In this case, keywords are list of Uniprot accession numbers.
>    - The *"Column number on which to apply the filter"*, in this case is the column that contains Uniprot accession numbers (`c1` as by default).
>    - *"Search for exact match ?"* You can perform exact or partial match with the keywords entered. Partial match is set by default. We keep default option `No` in this tutorial.
>    - To *"Enter keywords"*, you can either copy and paste list of keywords to the text area or choose a file that contains keywords in text format, in which each lines contains a keyword. Here we choose to `copy/paste` the following list of Uniprot accession number `P04264 P35908 P13645 Q5D862 Q5T749 Q8IW75 P81605 P22531 P59666 P78386`
>
>    > ### {% icon comment %} Outputs
>    > - **Filtered_Lacombe_et_al_2017.txt - Discarded_lines**: output list with the ten proteins (contaminants) removed from the original dataset (10 proteins)
>   > - **Filtered_Lacombe_et_al_2017.txt**: output contains the remaining proteins that will be considered for further analysis (151 proteins)
>    > 
>    {: .comment}
>
{: .hands_on}


# Check for the presence of biological contaminants

As EBC samples are obtained from air exhaled through the oral cavity, and even though the RTube collection device contained a saliva trap to separate saliva from the exhaled breath, contamination with salivary proteins had to be assessed. We decided to check the expression pattern for each protein of the "core" EBC proteome using the Human Protein Atlas (HPA). As HPA is indexed by Ensembl gene identifier (ENSG) we first need to convert Uniprot ID to Ensembl gene (ENSG). Secondly, check for proteins which are highly expressed in the salivary glands as reported by HPA, then in a third step, we filter out these proteins.

> ### {% icon hands_on %} Hands-on: Convert Uniprot ID to Ensembl gene ID
>
> 1. **ID Converter** {% icon tool %} with the following parameters:
>    - *"Enter IDs"*: `Input file containing IDs`
>        - {% icon param-file %} *"Select your file"*: `Filtered_Lacombe_et_al_2017.txt` (output of **Build tissue-specific expression dataset** {% icon tool %})
>    - *"Does file contain header"*: `Yes`
>    - *"Column number of IDs to map"*: `c1`
>    - *"Species"*: `Human (Homo sapiens)`
>        - *"Type/source of IDs"*: `Uniprot accession number (e.g. P31946)`
>        - In *"Target type"*:
>            - *"Target type of IDs you would like to map to"*: `Ensembl gene ID (e.g. ENSG00000166913)`
>
>    > ### {% icon comment %} Output
>    >
>    > - In the output file, a new column which contains Ensembl IDs was added (at the end)
>    {: .comment}
{: .hands_on}


> ### {% icon hands_on %} Hands-on: Check for proteins highly expressed in salivary glands
>
> 1. **Add expression data** {% icon tool %} with the following parameters:
>    - *"Enter your IDs (Ensembl gene IDs only, e.g. ENSG00000064787)"*: `Input file containing your IDs`
>        - {% icon param-file %} *"Select your file"*: `output of ID Converter` (output of **ID Converter** {% icon tool %})
>        - *"Column IDs: `c4`
>        - *"Does file contain header"*: `Yes`
>    - Numerous information can be extracted from the HPA source files, you can read user documentation at the end of the submission form of the tool for more detailed description. In this activity, in *"RNAseq/Ab-based expression data"*, we *"Select information to add to your list"*:  
>       - `Gene name`
>       - `Gene description`
>       - `RNA tissue category (according to HPA)`
>       - `RNA tissue specificity abundance in "Transcript Per Million`
>
>
>    > ### {% icon comment %} Outputs
>    >
>    > - In the outpu file, four columns were added (n°5, 6, 7 and 8) corresponding to the retrieved information from HPA.
>    {: .comment}
>
>   > ### {% icon comment %} Comments
>   >
>   > Scroll down the table, note at the end of the list (column n°8), that **AMY1B**, **CALML5**, **PIP**, **ZG16B**, **CST4**, **MUC7**, **CST1** and **CST2** have been reported as highly enriched in salivary gland with elevated RNA transcript specific TPM value for each, suggesting that these proteins may come from the saliva and not from the exhaled breath condensate. We thus will remove these biological contaminants from our initial protein set.
>   {: .comment}
>
{: .hands_on}


> ### {% icon hands_on %} Hands-on: Filter the data to remove the biological contaminants (i.e. proteins highly expressed in salivary glands)
>
>  **Filter by keywords and/or numerical value** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `Add expression data on data 6` (output of **Add expression data** {% icon tool %})
>    - Keep default option `Yes` for *"Does file contain header?"*
>    - *"Operation"*: `Discard`. What we want is to remove biological contaminants.
>    - Keep the `OR` option (by default) for the **operator** parameter. We don't need that parameter for a single filter.
>    - First add a  {% icon param-repeat %} *"Insert Filter by keywords"* box with a list of keywords to be filtered out.
>   In this step, we will filter out the lines that contain "salivary" in the column of RNA transcript specific TPM.
>    - *"Column number on which to apply the filter"* : `c8`, the column of RNA transcript specific TPM
>    - *"Search for exact match ?" : `No`
>    - To *"Enter keywords"*, you can either copy and paste list of keywords to the text area or choose a file that contains keywords in text format, in which each lines contains a keyword. Here we choose to `copy/paste` the keyword to be filtered out :  `salivary`
>
>    > ### {% icon comment %} Outputs
>   > Two output files are created:
>   > - **Filtered Add expression data on data 6 - Discarded lines** (12 proteins)
>   > - **FilteredAdd expression data on data 6** (151 proteins)
>{: .comment}
>
> > ### {% icon tip %} Tip: Using genes instead of keywords
> > Note also that a list of “gene” may have been entered (selected on the basis of their TPM value) applied to column n°5 instead of the keywords "salivary" to column n°8, as it has been done in [_Lacombe et al, 2018_](https://www.ncbi.nlm.nih.gov/pubmed/29189203).
> {: .tip}
>
{: .hands_on}


# Functional annotation of the EBC proteome (enrichment analysis)

The resulting list of 151 proteins identified in the two pooled EBC samples (excluding the 10 contaminants proteins) is now submitted to Gene Ontology (GO)-term enrichment analysis to determine functions that were significantly enriched in our EBC proteomic dataset compared to the lung proteome (corresponding to tissue-specific genes extracted from the Human Protein Atlas). To do so, we first build a lung reference proteome (that should be more representative of the studied sample rather than a full human proteome) that will be used for enrichment analysis performed with the ClusterProfiler tool (based on the R package clusterProfiler). 

> ### {% icon hands_on %} Hands-on: Build a lung reference proteome as a background for GO terms enrichment analysis
>
> 1. **Build tissue-specific expression dataset** {% icon tool %} with the following parameters:
>    - *"Experimental data source (antibody- or RNAseq-based)"*: `Expression profiles based on immunohistochemistry`
>    - *"Select tissue"*: `Lung` and `Bronchus`
>    - *"Expression level"*: `High`, `Medium` and `Low`
>    - *"Reliability score"*: `Enhanced` and `Supported`
>
>   > ### Output
>   > - **Tissue-specific expression from IHC**: List of the selected proteins. 
>   > 6 columns: 'Gene', 'Gene name' and the retrieved info from HPA. 
>   {: .comment}
{: .hands_on}


> ### {% icon tip %} Tip
> Note that expression information about respiratory cell types is retrieved (column 4; e.g. macrophages, pneumocytes, respiratory epithelial cells) 
> that could be used for further refinement of your reference background.
{: .tip}


As the ClusterProfiler tool (which we will use for the enrichment analysis) does not consider ENSG (Ensembl gene) identifiers as input, we need to convert IDs into either entrez Gene ID or Uniprot accession number. 

> ### {% icon hands_on %} Hands-on: Convert Ensembl ID to Uniprot and Entrez Gene ID
>
> 1. **ID Converter** {% icon tool %} with the following parameters:
>    - *"Enter IDs"*: `Input file containing IDs`
>        - {% icon param-file %} *"Select your file"*: `Tissue-specific expression from IHC` (output of **Build tissue-specific expression dataset** {% icon tool %})
>    - *"Does file contain header"*: `Yes`
>    - *"Column number of IDs to map"*: `c1`
>    - *"Species"*: `Human (Homo sapiens)`
>        - *"Type/source of IDs"*: `Ensembl gene ID (e.g. ENSG00000166913)`
>        - In *"Target type"*:
>            - *"Target type of IDs you would like to map to"*: `UniProt accession number (e.g. P31946)` and `Entrez gene ID (e.g. 7529)`
>
>    > ### {% icon comment %} Output
>    >
>    > - In the output file, 2 new columns have been added with the ID retrieved thanks to the conversion. 
>    {: .comment}
{: .hands_on}

Now we can perform the GO terms analysis. Input list is the EBC proteome to be analyzed after technical and biological contaminants 
removal, which is the output of biological contaminants filter step. 

> ### {% icon hands_on %} Hands-on: GO terms analysis
>
> 1. **GO terms classification and enrichment analysis** {% icon tool %} with the following parameters:
>    - *"Enter your IDs (UniProt Accession numer or Gene ID)"*: `Input file containing your IDs`
>        - {% icon param-file %} *"Choose a file that contains your list of IDs"*: `FilteredAdd expression data on data 6` (output of **Filter by keywords and/or numerical value** {% icon tool %})
>    - *"Select type/source of IDs"*: `UniProt accession number (e.g.:P31946)`
>    - - *"Species"*: `Homo sapiens`
>    - *"Select GO terms category"*: select all three options `Cellular Component`, `Biological process`, 
>    and `Molecular Function`
>    - *"Perform GO categories representation analysis?"*: `Yes`
>        - *"Ontology level (the higher this number, the deeper the GO level)"*: `3`
>    - *"Perform GO categories enrichment analysis?"*: `Yes`
>        - *"Define your own background IDs?"*: `Yes`
>            - *"Enter your background IDs (UniProt Accession number or Entrez Gene ID)"*: `Input file containing your background IDs`
>                - {% icon param-file %} *"Select file that contains your background IDs list"*: `ID converter on data 10` (output of **ID Converter** {% icon tool %})
>                - *"Column number of IDs"*: `c7`
>            - *"Select type of background IDs"*: `UniProt Accession number`
>        - *"Graphical display"*: `dotplot`
>
>   > ### {% icon comment %} Output
>   >
>   > Results created in History panel are the following:    
>   >   - Cluster profiler
>   >   - ClusterProfiler diagram outputs (collection dataset of all graphical outputs)
>   >   - ClusterProfiler text files (collection dataset of all text files) 
>
>   >   The suffix “GGO” (GroupGO) corresponds to the results “GO categories representation analysis” option
>   >   (performs a gene/protein classification based on GO distribution at a specific level). The suffix
>   >   “EGO” (EnrichGO) corresponds to the results from the enrichment analysis (based on an
>   >   over-representation
>   test of GO terms against the lung reference background). Two types of graphical output are provided either
>   >   in the form of bar-plot or dot-plot. 
>   >   According to this analysis, the main biological processes over-represented in EBC compared to
>   >   lung were some processes related to the immune system and exocytosis (see EGO.BP.dot.png, for Enriched
>   >   Biological Process GO terms dot-plot representation in png format).
>    {: .comment}
{: .hands_on}


# Visualize EBC proteome on biological pathways (using Reactome)

The 151 proteins identified in EBC samples are now mapped to biological pathways and visualized
via the web service of Reactome, an open access, manually curated and peer-reviewed human pathway
database that aims ti provide intuitive bioinformatics tools for the visualization,
interpretation and analysis of pathway knowledge.

> ### {% icon hands_on %} Hands-on: Protein list mapping on Reactome database
>
> 1. **Query Reactome pathway database** {% icon tool %} with the following parameters:
>    - *"Input IDs (UniProt Accession number, Entrez Gene ID or Gene Name"*: `Input file containing your IDs`
>    - {% icon param-file %} *"Input file containing your IDs"*: `Filtered_Add expression data on data 6` (output of **Filter by keywords and/or numerical value** {% icon tool %})
>    - *"Column number of IDs"*: `c1`
>    - *"Species"*: `Human (Homo sapiens)`
>   > ### {% icon comment %} Output
>   > 
>   > You can click on a link that opens the connection on Reactome and shows the image below: 
>   > ![the mapping of your IDs on the database](../../images/reactome.png).  
>   > 
>   {: .comment}
{: .hands_on}


Examine {% icon galaxy-eye %} the Reactome map of your IDs to see the context of your biological pathways. 


# Comparison with other proteomic datasets from previous studies

> ### {% icon hands_on %} Hands-on: Lists comparison with Venn diagramm tool
>
> 1. **Venn diagram** {% icon tool %} with the following parameters:
>    - In *"List to compare"*:
>        - {% icon param-repeat %} *"Insert List to compare"*
>            - *"Enter your list"*: `Input file containing your list`
>                - {% icon param-file %} *"Select your file"*: `kept_lines` (output of **Filter by keywords and/or numerical value** {% icon tool %})
>            - *"Enter the name of this list"*: `Lacombe et al`
>        - {% icon param-repeat %} *"Insert List to compare"*
>            - *"Enter your list"*: `Input file containing your list`
>                - {% icon param-file %} *"Select your file"*: `output` (Input dataset)
>            - *"Enter the name of this list"*: `Bredberg et al`
>        - {% icon param-repeat %} *"Insert List to compare"*
>            - *"Enter your list"*: `Input file containing your list`
>                - {% icon param-file %} *"Select your file"*: `output` (Input dataset)
>            - *"Enter the name of this list"*: `Mucilli et al`
>
>   > ### {% icon comment %} Output
>   > 
>   > The Venn diagram shows the number of proteins specific and in common between the 3 lists. 
>   > ![Graphical output of the Venn diagram](../../images/Venn-proteome-annot.png).  
>   > 
>   {: .comment}
{: .hands_on}

# Conclusion
{:.no_toc}

ProteoRE offers a panel of tools to annotate a protein list. 
We showed that it is possible to make ID conversion, perform tissu-expression annotation, but also 
Gene Ontology analysis as well as Reactome interogation. 
This allows the user to go deeper in the analysis of proteomics analyses results. 

