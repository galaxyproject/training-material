---
layout: tutorial_hands_on

title: Using MultiGSEA
zenodo_link: 'https://zenodo.org/records/14216972'
questions:
- How to use MultiGSEA for GSEA-based pathway enrichment for multiple omics layers?
objectives:
- Perform GSEA-based pathway enrichment for transcriptomics, proteomics, and metabolomics data.
- Understand how to combine p-values across multiple omics layers.
time_estimation: 1H
key_points:
- MultiGSEA provides an integrated workflow for pathway enrichment analysis across multi-omics data.
- Supports pathway definitions from several databases and robust ID mapping.
contributors:
- stehling


---


The multiGSEA package was designed to run a robust GSEA-based pathway enrichment for multiple omics layers. The enrichment is calculated for each omics layer separately and aggregated p-values are calculated afterwards to derive a composite multi-omics pathway enrichment.

Pathway definitions can be downloaded from up to eight different pathway databases by means of the graphite Bioconductor package (Sales, Calura, and Romualdi 2018). Feature mapping for transcripts and proteins is supported towards Entrez Gene IDs, Uniprot, Gene Symbol, RefSeq, and Ensembl IDs. The mapping is accomplished through the AnnotationDbi package (Pagès et al. 2019) and currently supported for 11 different model organisms including human, mouse, and rat. ID conversion of metabolite features to Comptox Dashboard IDs (DTXCID, DTXSID), CAS-numbers, Pubchem IDs (CID), HMDB, KEGG, ChEBI, Drugbank IDs, or common metabolite names is accomplished through the AnnotationHub package metabliteIDmapping. This package provides a comprehensive ID mapping for more than 1.1 million entries.

This tutorial covers a simple example workflow illustrating how the multiGSEA package works. The omics data sets that will be used throughout the example were originally provided by Quiros et al. (Quirós et al. 2017). In their publication the authors analyzed the mitochondrial response to four different toxicants, including Actinonin, Diclofenac, FCCB, and Mito-Block (MB), within the transcriptome, proteome, and metabolome layer.


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Preparing the Data

To perform pathway enrichment with MultiGSEA, you'll need omics datasets in the file type TSV . These datasets contain columns for feature Symbol, logFC pValue and adj.p-values. We'll use example data provided on Zenodo.

## Get data

> ### Data Upload
>
> 1. Create a new history for this tutorial.
    {% snippet faqs/galaxy/histories_create_new.md %}
> 2. Import the datasets from [Zenodo]({{ page.zenodo_link }})  into your Galaxy instance:
>   - **transcriptomics.tsv**
>   - **proteomics.tsv**
>   - **metabolomics.tsv**
>

# Running MultiGSEA

In this step, you'll use the MultiGSEA tool to perform GSEA-based pathway enrichment on the uploaded datasets.

>## Selecting parameters
><hands-on-title> Task description </hands-on-title>
>
> 1. Select the tool {% tool [multiGSEA](toolshed.g2.bx.psu.edu/repos/iuc/multigsea/multigsea/1.12.0+galaxy0) %} in Galaxy.
> 2. Configure the input parameters as follows:
>    - *"Select transcriptomics data"*: `Enabled`
>        - {% icon param-file %} *"Transcriptomics data"*: `Transcriptomics`
>    - *"Select proteomics data"*: `Enabled`
>        - {% icon param-file %} *"Proteomics data"*: `Proteomics`
>    - *"Select metabolomics data"*: `Enabled`
>        - {% icon param-file %} *"Metabolomics data"*: `Metabolomics`
> 3. You can also choose the Gene ID format for every data set. In this tutorial we will use the preset "SYMBOL" for transcriptomics and proteomics. For metabolomics we use HMDB.
> 4. Select in **Supported organisms** the organism of which the data is about. In our case we select `Homo sapiens (Human)`.
> 5. **Pathway databases**: Select relevant databases. For the tutorial we choose `KEGG`
> 6. **Combine p-values method**: Choose a method (here `Stouffer` for balanced weighting).
> 7. **P-value correction method** (for controlling false discovery rate): Choose `Holm`.
> 8. Click on `Run Tool`
>
>    {% snippet faqs/galaxy/tools_run.md %}
{: .hands_on}

---


> <question-title></question-title>
>
> 1. What file format is required for the input data in MultiGSEA?
> 2. What is the purpose of the “Combine p-values method” parameter, and which method was selected in this tutorial?
> 3. Why is it important to select pathway databases (e.g., KEGG) when using MultiGSEA?
>
> > <solution-title></solution-title>
> >
> > 1. The required file format is TSV.
> > 2. The “Combine p-values method” parameter is used to aggregate p-values across omics layers. In this tutorial, the method Stouffer was selected to apply balanced weighting.
> > 3. Selecting pathway databases ensures that the analysis uses appropriate and relevant pathway definitions for enrichment.
> >
> {: .solution}
>
{: .question}

> # Conclusion
> 
> In this tutorial, you explored the capabilities of MultiGSEA for performing pathway enrichment analysis across multiple omics layers, including transcriptomics, proteomics, and metabolomics data. By following the steps, you learned how to:
>
> - Prepare and upload the required omics datasets.
> - Configure and execute the MultiGSEA tool within Galaxy.
> -  Combine p-values from different omics layers to derive a unified perspective on pathway enrichment.