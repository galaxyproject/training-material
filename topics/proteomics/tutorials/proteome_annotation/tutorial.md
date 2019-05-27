---
layout: tutorial_hands_on

title: Annotating a protein list identified by LC-MS/MS experiments
zenodo_link: ''
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

---

# Introduction
{:.no_toc}

<!-- This is a comment. -->
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

For this tutorial, we will use three datasets, the list of proteins identified by LC-MS/MS in the exhaled breath condensate (EBC) from Lacombe *et al.* and two others EBC proteomes previously published.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial and give it a name
>   {% include snippets/create_new_history.md %}

> 2. Import the files from [Zenodo]() or from the shared data library (ask your instructors)
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>

>
{: .hands_on}


# Filter out technical contaminants

A group of 10 proteins were identified in both “technical” control samples with an enrichment in EBC samples below a fixed threshold. These proteins were thus considered to be technical contaminants (see list of proteins in Table 4 in [_Lacombe et al. 2018_](https://www.ncbi.nlm.nih.gov/pubmed/29189203)) and have to be removed from the initial dataset.

> ### {% icon hands_on %} Hands-on: Filter by keyword and  numerical values
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

<!-- ***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

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
{: .question} -->

# Check for the presence of biological contaminants

As EBC samples are obtained from air exhaled through the oral cavity, and even though the RTube collection device contained a saliva trap to separate saliva from the exhaled breath, contamination with salivary proteins had to be assessed. We decided to check the expression pattern for each protein of the "core" EBC proteome using the Human Protein Atlas (HPA). As HPA is indexed by Ensembl gene identifier (ENSG) we first need to convert Uniprot ID to Ensembl gene (ENSG). Secondly, check for proteins which are highly expressed in the salivary glands as reported by HPA, then in a third step, we filter out these proteins.

> ### {% icon hands_on %} Hands-on: Convert Uniprot ID to Ensembl gene
>
> 1. **ID Converter** {% icon tool %} with the following parameters:
>    - *"Enter IDs"*: `Input file containing IDs`
>        - {% icon param-file %} *"Select your file"*: `Filter_by_keywords_or_numerical_value_on_Lacombe_et_al_2017.txt` (output of **Build tissue-specific expression dataset** {% icon tool %})
>    - *"Does file contain header"*: `Yes`
>     - *"Column number of IDs to map"*: `c1`
>    - *"Species"*: `Human (Homo sapiens)`
>        - *"Type/source of IDs"*: `Uniprot accession number (e.g. P31946)`
>        - In *"Target type"*:
>            - *"Target type of IDs you would like to map to"*: `Ensembl gene ID (e.g. ENSG00000166913)`
>
>
>    > ### {% icon comment %} Output
>    >
>    > - **ID Converter on data 4**: In this file, a new column which contains Ensembl IDs was added.
>    {: .comment}
{: .hands_on}

<!-- ***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

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
{: .question} -->


> ### {% icon hands_on %} Hands-on: Check for proteins highly expressed in salivary glands
>
> 1. **Add expression data** {% icon tool %} with the following parameters:
>    - *"Enter your IDs (Ensembl gene IDs only, e.g. ENSG00000064787)"*: `Input file containing your IDs`
>        - {% icon param-file %} *"Select your file"*: `output of ID Converter step ID Converter on data 4` (output of **ID Converter** {% icon tool %})
>        - *"Column IDs: `c4`
>        - *"Does file contain header"*: `Yes`
>        - 
>    - Numerous information can be extracted from the HPA source files, you can read user documentation at the end of the submission form of the tool for more detailed description. In this activity, in *"RNAseq/Ab-based expression data"*, we *"Select information to add to your list"*:  
>       - `Gene name`
>       - `Gene description`
>       - `RNA tissue category (according to HPA)`
>       - `RNA tissue specificity abundance in "Transcript Per Million`
>
>
>    > ### {% icon comment %} Outputs
>    >
>    > - **Add expression data on data 6**: Four columns were added (n°5, 6, 7 and 8) corresponding to the HPA information previously selected.
>    {: .comment}
>
>   > ### {% icon comment %} Comments
>   >
>   > Scroll down the table, note at the end of the list (column n°8), that **AMY1B**, **CALML5**, **PIP**, **ZG16B**, **CST4**, **MUC7**, **CST1** and **CST2** have been reported as highly enriched in salivary gland with elevated RNA transcript specific TPM value for each, suggesting that these proteins may come from the saliva and not from the exhaled breath condensate. We thus will remove these biological contaminants from our initial protein set.
>   {: .comment}
>
{: .hands_on}

<!-- ***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

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
{: .question} -->


> ### {% icon hands_on %} Hands-on: Filter out the contaminant
>
>  **Filter by keywords and/or numerical value** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `Add expression data on data 6` (output of **Add expression data** {% icon tool %})
>    - Keep default option `Yes` for *"Does file contain header?"*
>    - *"Operation"*: `Discard`. We want to remove technical contaminants.
>    - Keep the `OR` option (by default) for the **operator** parameter. We don't need that parameter for a single filter.
>    - First add a  {% icon param-repeat %} *"Insert Filter by keywords"* box with a list of keywords to be filtered out.   In this step, we will filter out the lines that contain "salivary" in the column of RNA transcript specific TPM..
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
> > ### {% icon tip %} Tip
> > Note also that a list of “gene” may have been entered (selected on the basis of their TPM value) applied to column n°5 instead of the keywords "salivary" to column n°8, as it has been done in [_Lacombe et al, 2018_](https://www.ncbi.nlm.nih.gov/pubmed/29189203).
> {: .tip}
>
{: .hands_on}




## Sub-step with **Venn diagram**

> ### {% icon hands_on %} Hands-on: Task description
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



## Sub-step with **Filter by keywords and/or numerical value**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Filter by keywords and/or numerical value** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `output` (output of **Add expression data** {% icon tool %})
>    - *"Operation"*: `Discard`
>    - In *"Filter by keywords"*:
>        - {% icon param-repeat %} *"Insert Filter by keywords"*
>            - *"Column number on which to apply the filter"*: `c8`
>            - *"Enter keywords"*: `copy/paste`
>                - *"Copy/paste keywords to find (keep or discard)"*: `salivary`
>    - In *"Filter by numerical value"*:
>        - {% icon param-repeat %} *"Insert Filter by numerical value"*
>            - *"Column number on which to apply the filter"*: `[`
>            - *"Select operator"*: ``
>            - *"Value"*: `[`
>        - {% icon param-repeat %} *"Insert Filter by numerical value"*
>            - *"Column number on which to apply the filter"*: `]`
>            - *"Select operator"*: ``
>            - *"Value"*: `]`
>    - In *"Filter by range of numerical values"*:
>        - {% icon param-repeat %} *"Insert Filter by range of numerical values"*
>            - *"Column number on which to apply the filter"*: `[`
>            - *"Enter the bottom value"*: `[`
>            - *"Enter the top value"*: `[`
>            - *"inclusive range ?"*: `Yes`
>        - {% icon param-repeat %} *"Insert Filter by range of numerical values"*
>            - *"Column number on which to apply the filter"*: `]`
>            - *"Enter the bottom value"*: `]`
>            - *"Enter the top value"*: `]`
>            - *"inclusive range ?"*: `Yes`
>    - *"Sort by column ?"*: `Yes`
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

## Sub-step with **Query Reactome pathway database**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Query Reactome pathway database** {% icon tool %} with the following parameters:
>    - *"Input IDs (UniProt Accession number, Entrez Gene ID or Gene Name"*: `Input file containing your IDs`
>        - {% icon param-file %} *"Input file containing your IDs"*: `kept_lines` (output of **Filter by keywords and/or numerical value** {% icon tool %})
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

## Sub-step with **GO terms classification and enrichment analysis**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **GO terms classification and enrichment analysis** {% icon tool %} with the following parameters:
>    - *"Enter your IDs (UniProt Accession numer or Gene ID)"*: `Input file containing your IDs`
>        - {% icon param-file %} *"Choose a file that contains your list of IDs"*: `kept_lines` (output of **Filter by keywords and/or numerical value** {% icon tool %})
>    - *"Select type/source of IDs"*: `UniProt accession number (e.g.:P31946)`
>    - *"Select GO terms category"*: ``
>    - *"Perform GO categories representation analysis?"*: `Yes`
>        - *"Ontology level (the higher this number, the deeper the GO level)"*: `3`
>    - *"Perform GO categories enrichment analysis?"*: `Yes`
>        - *"Define your own background IDs?"*: `Yes`
>            - *"Enter your background IDs (UniProt Accession number or Entrez Gene ID)"*: `Input file containing your background IDs`
>                - {% icon param-file %} *"Select file that contains your background IDs list"*: `output` (output of **ID Converter** {% icon tool %})
>                - *"Column number of IDs"*: `c7`
>            - *"Select type of background IDs"*: `UniProt Accession number`
>        - *"Graphical display"*: ``
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


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.