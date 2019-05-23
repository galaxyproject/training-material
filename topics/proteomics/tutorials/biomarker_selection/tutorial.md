---
layout: tutorial_hands_on

title:  Biomarker candidate identification
zenodo_link: ''
questions:
- How to mine public databases to retrieve info?
- How to build a selection strategy by applying
  successive biochemical/cellular criteria to a list of gene/protein?
- How to select biomarkers candidates using experimental information (transcriptomics & proteomics)
  and annotation from public databases?

objectives:
- Build a workflow implementing a strategy for the selection of tissue-leakage biomarkers using ProteoRE

time_estimation: 3H

key_points:
- biomarker candidates selection workflow, public proteomics data retrieval and annotation

contributors:
- combesf
- davidchristiany
- vloux
- yvandenb

---


# Introduction
{:.no_toc}

A biomarker is a measurable biological component that can be routinely detected in clinical practice and reflects a disease state,
response to therapeutic treatment, or other relevant biological state.

[ProteoRE Galaxy instance](http://www.proteore.org) provides necessary tools to execute a complete biomarkers selection pipeline.
In this tutorial we introduce successively the tools of this pipeline, and guide you to execute them in order to complete the entire
pipeline on a concrete example. This strategy is described by {% cite Nguyen2019 %}.


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}



# Global view of the strategy
{:.no_toc}

For this tutorial, no input data are required as the first steps will be to select data from
public databases with ProteoRE tools.

The strategy consists in selecting, one step after another, the most interesting candidates biomarkers.
Our use-case here is to identify candidate biomarkers for myocardial infarction tissue-leakage.

Criteria candidate biomarkers have to fulfill through this pipeline are:
- heart-specificity
- cytoplasmic localization
- detection in LC-MS/MS experiments already done

![Overview of the pipeline described in this tutorial](../../images/pipeline.png "Pipeline of the tutorial, from Nguyen et al., as of in 2019.")


# Selection of tissue-specific proteins based on experimental data available in HPA


> ### {% icon hands_on %} Hands-on: Build tissue-specific expression dataset based on ImmunoHistoChemistry
>
> 1. **Create a new history** and give it a name.
>    {% include snippets/create_new_history.md %}
>
> 2. **Build tissue-specific expression dataset** {% icon tool %} with the following parameters:
>    - *"Experimental data source (antibody- or RNAseq-based)"*: `Expression profiles based on immunohistochemistry`
>    - *"Select tissue"*: `Heart muscle`
>    - *"Expression level"*: `High` and `Medium`
>    - *"Reliability score"*: `Enhanced` and `Supported`
>
>   > ### Output
>   > - **Tissue-specific expression from IHC** (1596 lines): List of the selected proteins.
>   > 6 columns: 'Gene', 'Gene name' and the retrieved info from HPA.
>   {: .comment}
{: .hands_on}

We will now rerun the same tool but to select transcripts according to their expression profile.


> ### {% icon hands_on %} Hands-on: Build tissue-specific expression dataset based on RNAseq
>
> 1. **Build tissue-specific expression dataset** {% icon tool %} with the following parameters:
>    - *"Experimental data source (antibody- or RNAseq-based)"*: `RNA levels based on RNA-seq experiments`
>    - *"Select tissue"*: `Heart muscle`
>
>   > ### Output
>   > - **Tissue-specific expression from RNAseq** (19613 lines): List of the selected transcripts.
>   > 4 columns: 'Gene', 'Gene name' and the retrieved info from HPA.
>   {: .comment}
{: .hands_on}


This second list must be reduced by removing transcripts that are not highly enriched in heart muscle.
To do so, a filter is applied on the expression value provided by HPA and measured in TPM (last column of the output file).
In ProteoRE we'll use the "Filter by keywords and/or numerical value" tool.


> ### {% icon hands_on %} Hands-on: Filter on expression value criterium
>
> 1. **Filter by keywords and/or numerical value** {% icon tool %} with the following parameters:
>    - *"Operation"*: `Discard`
>    - In *"Filter by numerical value"*:
>        - {% icon param-repeat %} *"Insert Filter by numerical value"*
>            - *"Column number on which to apply the filter"*: `c4`
>            - *"Select operator"*: `<=`
>            - *"Value"*: `10.0`
>    - *"Sort by column ?"*: `Yes`
>
>
> > ### {% icon question %} Questions
> >
> > How many lines are there in the file of heart transcripts
> > with a TPM value >10 ?
> >
> > > ### {% icon solution %} Solution
> > > 5257 lines.
> > {: .solution}
> {: .question}
>
>
{: .hands_on}

We have now 2 datasets of heart-muscle proteins/transcripts, based on IHC data or TPM value.

We want now to select candidate biomarkers that are expressed in the heart muscle according to **both** IHC and RNA-seq data, using the Jvenn tool.


> ### {% icon hands_on %} Hands-on: Venn diagram
>
> 1. **Venn diagram** {% icon tool %} with the following parameters:
>    - In *"List to compare"*:
>        - {% icon param-repeat %} *"Insert List to compare"*
>            - *"Enter your list"*: `Input file containing your list`
>                - *"Enter the name of this list"*: `heart IHC`
>        - {% icon param-repeat %} *"Insert List to compare"*
>            - *"Enter your list"*: `Input file containing your list`
>                - *"Enter the name of this list"*: `heart RNAseq`
>
>   > ### {% icon question %} Questions
>   > How many IDs are in common to both IHC and RNA-seq lists ?
>   >   > ### {% icon solution %} Solution
>   >   > 931 (you can visualize it in [the output of the Venn](#figure-1))
>   > {: .solution}
>   {: .question}
{: .hands_on}

You can see the the graphical output of the Venn:

![Venn diagram output](../../images/jVenn_chart-tuto2.png)

For greater clarity we'll keep only the column with those 931 IDs to continue our pipeline.


> ### {% icon hands_on %} Hands-on: Cut
>
> 1. **Cut** {% icon tool %} with the following parameters:
>    - *"Cut by"*: `fields`
>        - *"List of Fields"*: `c['3']`
>
{: .hands_on}

Now we'll filter this dataset not to keep the 'NA' lines.

> ### {% icon hands_on %} Hands-on: Filter by keywords and/or numerical value
>
> 1. **Filter by keywords and/or numerical value** {% icon tool %} with the following parameters:
>    - *"Operation"*: `Discard`
>    - In *"Filter by keywords"*:
>        - {% icon param-repeat %} *"Insert Filter by keywords"*
>            - *"Search for exact match?"*: `Yes`
>            - *"Enter keywords"*: `copy/paste`
>                - *"Copy/paste keywords to find (keep or discard)"*: `NA`
>    - *"Sort by column ?"*: `Yes`
>
>
{: .hands_on}

Let's rename the 931 IDs dataset in "heart931" for simplification.

Pipeline will then continue based on those 931, from which we have to select biomarkers that are
highly specific to the heart using additional expression data (still from HPA).


> ### {% icon hands_on %} Hands-on: Add expression data
>
> 1. **Add expression data** {% icon tool %} with the following parameters:
>    - *"Enter your IDs (Ensembl gene IDs only, e.g. ENSG00000064787)"*: `Input file containing your IDs`
>    - In *"RNAseq/Ab-based expression data"*:
>        - *"Select information to add to your list"*: `Gene name, Gene description, RNA tissue category, RNA tissue specificity abundance in ‘Transcript Per Million`
>
>
{: .hands_on}

We wish to focus on transcripts that have been classed as (according to the HPA definition):
* "tissue enriched" (expression in one tissue at least fivefold higher than all other tissues),
* "group enriched" (fivefold higher average TPM in a group of two to to seven tissues compared to all other tissues) and
* "tissue enhanced" (fivefold higher average TPM in one or more tissues/cell lines compared to the mean TPM for all tissues)

This information is listed in the column 4 : "RNA tissue category" of the result dataset.

Let's use the "Filter by keywords and/or numerical value" tool to select the candidate biomarkers based on this
"RNA tissue category" criterium.


> ### {% icon hands_on %} Hands-on: Filter by keywords and/or numerical value
>
> 1. **Filter by keywords and/or numerical value** {% icon tool %} with the following parameters:
>    - In *"Filter by keywords"*:
>        - {% icon param-repeat %} *"Insert Filter by keywords"*
>            - *"Column number on which to apply the filter"*: `c4`
>            - *"Enter keywords"*: `copy/paste`
>                - *"Copy/paste keywords to find (keep or discard)"*: `enriched enhanced`
>    - *"Sort by column ?"*: `Yes`
>
>   > ### Output
>   > - **Filtered Add_expression_data_on_data_8**: output list of the heart biomarkers with RNA tissue category
>   > containing "enriched" or "enhanced" (115 lines = what we are interested in)
>   > - **Filtered Add_expression_data_on_data_8 - discarded lines**: output list of the heart biomarkers with
>   > RNA tissue category NOT containing "enriched" or "enhanced" (not what we are interested in)
>   {: .comment}
{: .hands_on}

We now have identified 115 candidates considered to have significantly higher expression in heart muscle according to HPA criteria.
Let's call the dataset where are those 115 candidates '**heart115**'.


# Annotation of this protein list with biochemical and cellular features


Candidate biomarkers we want to identify have to be cytoplasmic and without transmembrane domains (TMD).
Thus we will retrieve protein features from neXtProt to retrieve those informations.

Since HPA only considers ENSG identifiers (related to the gene), although neXtProt uses UniProt
identifiers (related to proteins), first thing to do is to map the Ensembl identifiers contained our list of (115)
candidates to their corresponding UniProt accession number. The tool **ID Converter** is what we need to do so.


> ### {% icon hands_on %} Hands-on: ID Converter
>
> 1. **ID Converter** {% icon tool %} with the following parameters:
>    - *"Enter IDs"*: `Input file containing IDs`
>    - *"Species"*: `Human (Homo sapiens)`
>        - *"Type/source of IDs"*: `Ensembl gene ID (e.g. ENSG00000166913)`
>        - In *"Target type"*:
>            - *"Target type of IDs you would like to map to"*: `Uniprot accession number, Uniprot ID`
>
>   > ### Output
>   > **ID converter on data 11**: In this dataset, 2 columns (columns 6 and 7, at the end) which contain
>   > UniProt accession number and ID are added.
>   {: .comment}
{: .hands_on}

We have now UniProt IDs for the 115 candidate biomarkers: we are able to collect protein features from neXtProt. For this purpose,
we use the **Add protein features** ProteoRE tool.


> ### {% icon hands_on %} Hands-on: Add protein features
>
> 1. **Add protein features** {% icon tool %} with the following parameters:
>    - *"Enter your IDs (neXtProt or UniProt)"*: `Input file containing your IDs `
>        - *"Column IDs (e.g : Enter c1 for column n°1)"*: `c6`
>    - In *"Select features"*:
>        - *"Physico-Chemical Features"*: `Number of transmembrane domains`
>        - *"Localization"*: `Subcellular Location`
>        - *"Disease information"*: `Yes`
>
>   > ### Output
>   > **Add information from NextProt**: In this file (431 lines), 3 columns (columns 8, 9 and 10)
>   > were added (at the end). These columns present TMDomains, Subcell Location and Diseases info.
>   {: .comment}
{: .hands_on}

With this dataset, we can select proteins reported as localized in the cytoplasm and having
no transmembrane domains by running the Filter by keywords and/or numerical value tool.


> ### {% icon hands_on %} Hands-on: Filter by keywords and/or numerical value
>
> 1. **Filter by keywords and/or numerical value** {% icon tool %} with the following parameters:
>    - *"Select an operator to combine your filters (if more than one)"*: `AND`
>    - In *"Filter by keywords"*:
>        - {% icon param-repeat %} *"Insert Filter by keywords"*
>            - *"Column number on which to apply the filter"*: `c9`
>            - *"Enter keywords"*: `copy/paste`
>                - *"Copy/paste keywords to find (keep or discard)"*: `cytoplasm cytosol`
>    - In *"Filter by numerical value"*:
>        - {% icon param-repeat %} *"Insert Filter by numerical value"*
>            - *"Column number on which to apply the filter"*: `c8`
>            - *"Value"*: `0`
>    - *"Sort by column ?"*: `Yes`
>        - *"Sort result files by:"*: `c5`
>        - *"Sort in descending order ?"*: `Yes`
>
>
>   > ### Output
>   > - **Filtered Add_information_from_neXtProt**: output list of the proteins having a cytoplasmic
>   > location and no TMD (48 proteins)
>   > - **Filtered Add_information_from_neXtProt - discarded lines**: output list of the proteins NOT
>   > cytoplasmic and having at least 1 TMD.
>   {: .comment}
{: .hands_on}

We have now 48 proteins.

Next step : to identify proteins already seen in LS MS/MS experiments.


# Check whether these proteins have already been detected by LC-MS/MS experiments

> ### {% icon hands_on %} Hands-on: Get MS/MS observations in tissue/fluid
>
> 1. **Get MS/MS observations in tissue/fluid** {% icon tool %} with the following parameters:
>    - *"Enter your IDs (UniProt Accession number only)"*: `Input file containing your IDs `
>        - *"Column of IDs"*: `c6`
>    - *"Proteomics dataset (biological sample)"*:
>        - `Human Heart`
>        - `Human Plasma non glyco`
>
>   > ### Output
>   > **Get MS/MS observations in tissue/fluid on data 15**: In this file, 2 columns (11 and 12, at the end)
>   > were added with the info of number of times peptides were seen by MS/MS.
>   {: .comment}
{: .hands_on}

Let's now keep only proteins that have already been seen by MS/MS in the plasma (last column of the file).


> ### {% icon hands_on %} Hands-on: Filter by keywords and/or numerical value
>
> 1. **Filter by keywords and/or numerical value** {% icon tool %} with the following parameters:
>    - *"Operation"*: `Discard`
>    - In *"Filter by keywords"*:
>        - {% icon param-repeat %} *"Insert Filter by keywords"*
>            - *"Column number on which to apply the filter"*: `c12`
>            - *"Search for exact match?"*: `Yes`
>            - *"Enter keywords"*: `copy/paste`
>                - *"Copy/paste keywords to find (keep or discard)"*: `NA`
>    - *"Sort by column ?"*: `Yes`
>
>   > ### Output
>   > - **Filtered Get MS/MS observations in tissue/fluid on data 15**: output list of the proteins
>   > whose some peptides have been seen in plasma (21 proteins)
>   > - **Filtered Get MS/MS observations in tissue/fluid on data 15 - discarded lines**:
>   > output list of proteins with no peptides seen in the plassma
>   {: .comment}
{: .hands_on}



# Conclusion
{:.no_toc}
At the end of the process we end up with a list of 21 biomarkers that are **highly enriched in heart muscle**, localized
**in the cytosol** and **detectable by MS in the plasma**.

Briefly and from a biological point of view, 3 of these proteins exhibit a relative low detection level in the plasma
compared to heart muscle tissue, and are reported with a very high heart-muscle-specific RNA abundance.
These potential mechanistic biomarkers of myocardial infarction include (i) cardiac Troponin I type 3 (TNNI3 (P19429))
that is routinely used as the most specific marker of myocardial injury and (ii) the heart-type fatty acid-binding protein
(FABP3 (P05413)) that has been proposed as a diagnostic and prognostic marker for acute and chronic cardiac injury.

Extraction of workflow is first valuable for analyses reproductibility and tracability.

Moreover, a workflow can also be made reusable with modifiable parameters, and as a result this strategy can
be applied to other types of tissue injury (e.g. brain, liver, kidney).
