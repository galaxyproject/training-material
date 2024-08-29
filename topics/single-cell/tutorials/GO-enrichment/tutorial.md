---
layout: tutorial_hands_on

title: GO Enrichment Analysis on Single-Cell RNA-Seq Data
zenodo_link: 'https://zenodo.org/records/13343340'

questions:
- What is Gene Ontology (GO) enrichment analysis, and why should I perform it on my marker genes?
- How can I use GO enrichment analysis to better understand the biological functions of the genes in my clusters?
- Can GO enrichment analysis help me confirm that my clusters represent distinct cell types or states?
- How can I visualize my GO enrichment results to make them easier to understand and interpret?
  
objectives:
- Understand the Role of GO Enrichment in Single-Cell Analysis.
- Use marker genes from different cell clusters or conditions for GO enrichment analysis. 
- Compare enrichment across experimental conditions (e.g., wild type vs. knockout) to uncover functional changes associated with genetic or environmental perturbations.
- link GO enrichment results with a previously annotated cell clusters, providing a clearer picture of the functional roles of different cell populations.
  
requirements:
-
    type: "internal"
    topic_name: single-cell
    tutorials:
        - scrna-case_basic-pipeline
tags:
- 10x
- paper-replication
  
time_estimation: 3H

key_points:
- GO enrichment helps make sense of your data and understand what makes each cell cluster/condition unique.
- GO enrichment analysis is used to discover new insights about how cells work, which can lead to better understanding of biological processes and diseases.
  
contributors:
- nomadscientist
- MennaGamal
---


# Introduction

In the previous tutorial [Filter, plot and explore single-cell RNA-seq data with Scanpy], we took an important step in our single-cell RNA sequencing analysis by identifying marker genes for each of the clusters in our dataset. These marker genes are crucial, as they help us distinguish between different cell types and states, giving us a clearer picture of the cellular diversity within our samples.
However, simply identifying marker genes is just the beginning. To truly understand what makes each cluster unique, we need to dive deeper into the biological functions these genes are involved in. This is where Gene Ontology (GO) enrichment analysis comes into play.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Data description

In this tutorial will use the following datasets:
**Marker Genes:**
We'll start with two input datasets of marker genes (Study sets):
  *Marker genes per cell cluster: This dataset lists the genes that are significantly different in each cell cluster.
  *Marker genes per condition (wt and ko): This dataset lists the genes that are significantly different between the wild-type (wt) and knockout (ko) conditions.
**GO Enrichment Files:**
We'll also use three additional files for GO enrichment analysis. 
  *Gene Ontology: This file contains information about Gene Ontology terms.
  *GO Annotations: This file maps genes to their corresponding GO terms.
  *Population set: This file provides a list of genes used as a background gene set for the analysis.
*Note:* There are several online databases available for downloading GO and GO Annotations files, including the Gene Ontology website, Ensembl, and the UCSC Genome Browser.
  
> <comment-title>Concept behind GO Enrichment Analysis</comment-title>
> The goal of GO enrichment analysis is to interpret the biological significance of long lists of marker genes. By summarizing these genes into a shorter list of enriched GO terms. The analysis works by comparing each GO term between your list of marker genes and a background gene set. Statistical tests are then used to calculate a p-value that indicates whether a particular GO term is significantly enriched in the marker gene list compared to the background.
The concept of enrichment analysis is well explained in this [youtube video] (https://www.youtube.com/watch?v=H1cUs6pql9s&t=157s) about Pathway Enrichment Analysis, here we use GO terms instead of pathways.
{: .comment}

## Get data

You can access the data for this tutorial in multiple ways:

1. **Importing from a history** - You can import [this history](https://usegalaxy.eu/u/mennagamal/h/go-enrichment-analysis-on-scrna-seq-data)

   {% snippet faqs/galaxy/histories_import.md %}

2. **Uploading from Zenodo** (see below)

> <hands-on-title> Data Upload from Zenodo </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) 
>
>    ```
>   {{ page.zenodo_link }}/files/Markers_clusters.tabular
>   {{ page.zenodo_link }}/files/Markers_genotype.tabular
    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype is `tabular`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular" %}
>
{: .hands_on}

# Important tips for easier analysis

{% snippet faqs/galaxy/tutorial_mode.md %}

{% snippet topics/single-cell/faqs/single_cell_omics.md %}

{% snippet faqs/galaxy/analysis_troubleshooting.md sc=true %}


# Data processing

To perform GO enrichment analysis on each cell cluster individually, we need to separate our "Markers_clusters Dataset" into seven files, one for each cluster and the "Markers_genotype Dataset" into 2 files, one for each condition. We'll use the "Split file" tool for this step. 

> <hands-on-title>File splitting</hands-on-title>
>
> 1. {% tool [Split file](toolshed.g2.bx.psu.edu/repos/bgruening/split_file_on_column/tp_split_on_column/0.4) %} with the following parameters:
>    - {% icon param-file %} *"File to select"*: `output` (Input dataset)
>    - *"on column"*: `c1`
>    - *"Include the header in all splitted files?"*: `Yes`
>
>
>    > <comment-title> Input Dataset </comment-title>
>    >
>    > As we have two datasets, one with the marker genes for all the seven clusters and one with the marker genes for the knockout (KO) and wild-type (WT) conditions. Make sure to repeat the analysis twice for the two different datasets. Otherwise you can run this [workflow] (https://usegalaxy.eu/u/mennagamal/w/unnamed-workflow) for parallel analysis of the datasets 
>    {: .comment}
>
{: .hands_on}

Next, we need to isolate the Ensembl gene IDs column from each file. We'll use the "Cut Columns" tool to achieve this.

> <hands-on-title> Extract Ensembl IDs </hands-on-title>
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c4`
>    - {% icon param-file %} *"From"*: `split_output` (output of **Split file** {% icon tool %})
>
>
>    > <comment-title> The gene format to use </comment-title>
>    >
>    In this example we extract column 4 because it contains the Ensembl gene IDs on which the subsequent steps are ideally working. While there are other gene formats like gene symbols, Entrez gene IDs, and more, make sure to check the specific format accepted by the tool you are using. There are also tools available to convert between different gene formats if needed. 
> 
> {: .comment}
{: .hands_on}


# GO Analysis using **GOEnrichment** tool

Now we will perform the GO Enrichment analysis on the list of ensembl gene IDs.

> <hands-on-title> GOEnrichment </hands-on-title>
>
> 1. {% tool [GOEnrichment](toolshed.g2.bx.psu.edu/repos/iuc/goenrichment/goenrichment/2.0.1) %} with the following parameters:
>    - {% icon param-file %} *"Gene Ontology File"*: `output` (Input dataset)
>    - {% icon param-file %} *"Gene Product Annotation File"*: `output` (Input dataset)
>    - {% icon param-file %} *"Study Set File"*: `out_file1` (output of **Cut** {% icon tool %})
>    - {% icon param-file %} *"Population Set File (Optional)"*: `output` (Input dataset)
>
>
>    > <comment-title> Population Set File Selection for GO Enrichment </comment-title>
>    >
>    > When choosing a background for GO enrichment analysis (Population Set File), it's important to consider the context of your data. While using a broad background (like all genes in the organism) is common, it might be more informative to limit the background to genes expressed in the specific tissue or cell type being profiled. In this tutorial we used only genes involved in the experiment before selecting the marker genes.
>    {: .comment}
>
{: .hands_on}


> <question-title></question-title>
>
> 1. Take a look at the enriched terms for cluster 7, Can you find any GO terms that are specific to this cluster?
> 2. Can we perform manual anntation of the cell cluster based on GO enrichment results?
>
> > <solution-title></solution-title>
> >
> > 1. Cluster 7 is enriched for terms like "regulation of cell death","T cell mediated cytotoxicity", and "peptidase activator activity involved in apoptotic process".
> > 2. By looking at the most enriched functions and with biological knowledge, we are able to infer the cell types for most of the clusters. We know that data come from thymus tissue, so we have a prior idea of the cell types we would expect. For example, the enrichment for those terms in cluster 7 can confirm our results that the cell type is Macrophages, as they help create a supportive environment for thymocyte maturation by clearing apoptotic cells and debris.
> >
> {: .solution}
>
{: .question}

# GO Analysis **gProfiler GOSt** tool

The gProfiler GOSt (Gene Ontology Sequential Testing) is another popular tool used to perform gene ontology (GO) enrichment analysis. In addition to providing enrichment results for the standard GO categories of Biological Process (BP), Cellular Component (CC), and Molecular Function (MF), the tool also analyzes enrichment across several other functional annotation databases, including: KEGG Pathways, Reactome Pathways, WikiPathways and TF Targets. It also gives a plot to better visualize the results.

> <hands-on-title> gProfiler GOSt </hands-on-title>
>
> 1. {% tool [gProfiler GOSt](toolshed.g2.bx.psu.edu/repos/iuc/gprofiler_gost/gprofiler_gost/0.1.7+galaxy11) %} with the following parameters:
>    - {% icon param-file %} *"Input is whitespace-separated list of genes, proteins, probes, term IDs or chromosomal regions."*: `out_file1` (output of **Cut** {% icon tool %})
>    - *"Organism"*: `Common organisms`
>        - *"Common organisms"*: `Mus musculus (Mouse)`
>    - In *"Tool settings"*:
>        - *"Export plot"*: `Yes`
>
>
>    > <comment-title> Picking the right species matter </comment-title>
>    >
>    > The species you select should match the species your genes come from. If you choose the wrong species, the tool might use incorrect information, leading to inaccurate results. For example, human genes behave differently from mouse genes, so selecting the correct species ensures the analysis is relevant to your data.
>    {: .comment}
>
{: .hands_on}


> <question-title></question-title>
>
> 1. Can you find enriched GO terms that are inline with the [published study] (https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2018.02523/full) findings in KO vs WT results files?
> 2. What might be happening to the stem cells in the KO mice compared to the WT mice?
>
> > <solution-title></solution-title>
> >
> > 1. In the KO g:GOSt result file, enrichment for the GO term "Negative regulation of stem cell differentiation" is found.
> > 2. This suggests that the KO condition is causing a delay in the differentiation of stem cells into mature T cells in the thymus which is inline with the study findings.
> >
> {: .solution}
>
{: .question}



# Conclusion

In this tutorial we have performed GO enrichment analysis on the differentially expressed genes between 2 different conditions and between different cell types.
