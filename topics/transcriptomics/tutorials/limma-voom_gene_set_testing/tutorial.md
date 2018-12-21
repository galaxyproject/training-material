---
layout: tutorial_hands_on
title: RNA-seq genes to pathways
zenodo_link: "https://figshare.com/s/f5d63d8c265a05618137"
enable: "false"
questions:
  - "What are the differentially expressed pathways in the mammary gland of pregnant versus lactating mice?"
objectives:
  - "Identification of differentially expressed pathways"
time_estimation: "3h"
key_points:
  - "Multiple methods can be used to help identify differentially expressed pathways"
contributors:
  - mblue9
  - bphipson
  - annatrigos
  - mritchie
  - hdashnow
  - charitylaw
---


# Introduction
{:.no_toc}

Measuring gene expression on a genome-wide scale has become common practice over the last two decades or so, with microarrays predominantly used pre-2008. With the advent of next generation sequencing technology in 2008, an increasing number of scientists use this technology to measure and understand changes in gene expression in often complex systems. As sequencing costs have decreased, using RNA-Seq to simultaneously measure the expression of tens of thousands of genes for multiple samples has never been easier. The cost of these experiments has now moved from generating the data to storing and analysing it.

There are many steps involved in analysing an RNA-Seq experiment. Analysing an RNA-seq experiment begins with sequencing reads. These are aligned to a reference genome, then the number of reads mapped to each gene can be counted. This results in a table of counts, which is what we perform statistical analyses on to determine differentially expressed genes and pathways. The purpose of this tutorial is to demonstrate how to perform gene set testing using tools in Galaxy. How to generate differentially expressed genes from reads (FASTQs) is covered in the accompanying tutorials [RNA-seq reads to counts]({{ site.baseurl }}/topics/transcriptomics/tutorials/limma-voom_fastqs_to_counts/tutorial.html) and [RNA-seq counts to genes]({{ site.baseurl }}/topics/transcriptomics/tutorials/limma-voom/tutorial.html).

## Mouse mammary gland dataset

The data for this tutorial comes from a Nature Cell Biology paper, [EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival](https://www.ncbi.nlm.nih.gov/pubmed/25730472)), Fu et al. 2015. Both the raw data (sequence reads) and processed data (counts) can be downloaded from Gene Expression Omnibus database (GEO) under accession number [GSE60450](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450).

This study examined the expression profiles of basal stem-cell enriched cells (B) and committed luminal cells (L) in the mammary gland of virgin, pregnant and lactating mice. Six groups are present, with one for each combination of cell type and mouse status. Note that two biological replicates are used here, two independent sorts of cells from the mammary glands of virgin, pregnant or lactating mice, however three replicates is usually recommended as a minimum requirement for RNA-seq. In this tutorial we will use the GEO counts file as a starting point for our analysis. Alternatively, you could create a count matrix from the raw sequence reads, as demonstrated in the [RNA-seq reads to counts tutorial]({{ site.baseurl }}/topics/transcriptomics/tutorials/limma-voom_fastqs_to_counts/tutorial.html). The GEO count file was generated from aligning the reads to the mouse `mm10` genome with the [Rsubread](https://www.biorxiv.org/content/early/2018/08/15/377762) aligner, followed by counting reads mapped to RefSeq genes with [featureCounts](https://academic.oup.com/bioinformatics/article/30/7/923/232889) (Liao, Smyth, and Shi 2014), see the [Fu paper](https://www.nature.com/articles/ncb3117) for details.

We used **limma-voom** for identifying differentially expressed genes here. Other popular alternatives are edgeR and DESeq2. Limma-voom has been shown to be perform well in terms of precision, accuracy and sensitivity [Costa-Silva, Domingues and Lopes 2017](https://www.ncbi.nlm.nih.gov/pubmed/29267363) and, due to its speed, it's particularly recommended for large-scale datasets with 100s of samples [Chen, Lun, Smyth 2016](https://f1000research.com/articles/5-1438/v2).

This is a Galaxy tutorial based on material from the [COMBINE R RNAseq workshop](http://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html), first taught [here](http://combine-australia.github.io/2016-05-11-RNAseq/). Some of the gene set testing material is inspired by the Cancer Research UK workshop [here](https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html)

![Tutorial Dataset](../../images/limma-voom_f2c/mouse_exp.png "Tutorial Dataset")


> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Preparing the inputs

We will use three files for this analysis:

 * **Count matrix** (genes in rows, samples in columns)
 * **Sample information** file (sample id, group)
 * **Gene annotation** file (gene id, symbol, description)

## Import data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this RNA-seq exercise e.g. `RNA-seq with limma-voom`
> 2. Import the mammary gland counts table and the associated sample information file.
>
>     To import the files, there are two options:
>     - Option 1: From a shared data library if available (ask your instructor)
>     - Option 2: From [Figshare](https://figshare.com/s/1d788fd384d33e913a2a)
>
>         > ### {% icon tip %} Tip: Importing data via links
>         >
>         > * Copy the link location
>         > * Open the Galaxy Upload Manager
>         > * Select **Paste/Fetch Data**
>         > * Paste the link into the text field
>         > * Press **Start**
>         {: .tip}
>
>         - You can paste both links below into the **Paste/Fetch** box:
>
>           ```
>       https://ndownloader.figshare.com/files/5057929?private_link=1d788fd384d33e913a2a
>       https://ndownloader.figshare.com/files/5999829?private_link=1d788fd384d33e913a2a
>           ```
>
>         - Select *"Genome"*: `mm10`
>
> 2. Rename the counts dataset as `seqdata` and the sample information dataset as `sampleinfo` using the {% icon galaxy-pencil %} (pencil) icon.
> 3. Check that the datatype is `tabular`.
>    If the datatype is not `tabular`, please change the file type to `tabular`.
>
>    > ### {% icon tip %} Tip: Changing the datatype
>    > * Click on the {% icon galaxy-pencil %} (pencil) icon displayed in your dataset in the history
>    > * Choose **Datatype** on the top
>    > * Select `tabular`
>    > * Press **Save**
>    {: .tip}
{: .hands_on}


Let’s take a look at the data. The `seqdata` file contains information about genes (one gene per row), the first column has the Entrez gene id, the second has the gene length and the remaining columns contain information about the number of reads aligning to the gene in each experimental sample. There are two replicates for each cell type and time point (detailed sample info can be found in file “GSE60450_series_matrix.txt” from the GEO website). The first few rows and columns of the seqdata file are shown below.

![seqdata file](../../images/limma-voom/seqdata.png "Count file (before formatting)"){: width="50%"}

The `sampleinfo` file contains basic information about the samples that we will need for the analysis. See below.

![sampleinfo file](../../images/limma-voom/sampleinfo.png "Sample information file (before formatting)"){: width="50%"}


# Gene Set Testing

Sometimes there is quite a long list of differentially expressed genes to interpret after a differential expression analysis, and it is usually infeasible to go through the list one gene at a time trying to understand it’s biological function. A common downstream procedure is gene set testing, which aims to understand which pathways/gene networks the differentially expressed genes are implicated in. There are many different gene set testing methods that can be applied and it can be useful to try several.

## Gene Ontology testing with **goseq**

We would like to know if there are biological categories that are enriched among the differentially expressed genes. To do this we will perform a Gene Ontology analysis, similar to the [RNA-seq ref-based tutorial]({{ site.baseurl }}/topics/transcriptomics/tutorials/ref-based/tutorial.html).

[Gene Ontology (GO)](http://www.geneontology.org/) analysis is widely used to reduce complexity and highlight biological processes in genome-wide expression studies. However, standard methods give biased results on RNA-seq data due to over-detection of differential expression for long and highly-expressed transcripts.

[goseq tool](https://bioconductor.org/packages/release/bioc/vignettes/goseq/inst/doc/goseq.pdf) provides methods for performing GO analysis of RNA-seq data, taking length bias into account. The methods and software used by goseq are equally applicable to other category based tests of RNA-seq data, such as KEGG pathway analysis.

goseq needs 2 files as inputs:
- a **differentially expressed genes** file. Information for all genes tested for differential expression (all genes after filtering lowly expressed). This file should have 2 columns:
    - the Gene IDs (unique within the file)
    - True (differentially expressed) or False (not differentially expressed)
- a **gene lengths** file. Information to correct for potential length bias in differentially expressed genes. This file should have 2 columns:
    - the Gene IDs (unique within the file)
    - the gene lengths

We will use the tables of differentially expressed results output from the limma-voom tool, for both the basal and luminal constrasts, and call genes differentially expressed if they have an adjusted P value below 0.01 and a fold change of 1.5 (equivalent to a $$log_{2} FC$$ of 0.58), as in the Fu paper. We can use the gene lengths from the original table we imported from GEO (`seqdata`). But if we didn't have that we could use a tool like **featureCounts** {% icon tool %} to output a gene lengths file. The original file with gene lengths contains all >20k genes, but we only want the ~15k we have in our differentially expressed genes file after filtering low counts, and in the same order. So we will join the lengths file with the differentially expressed genes file, keeping only the lengths information for genes present in the differentially expressed genes file. We can then cut out the columns we need for the two inputs (gene id, length) (gene id, DE status) and as a bonus they will both be sorted in the same order, what we need for goseq.

To generate the two input files we will use:
* **Compute** to add a column to the limma-voom table that gives genes meeting our adj.P and lfc thresholds the value "True" and all other genes the value "False". We want genes that have a lfc < -0.58 (downregulated) and lfc > 0.58 (upregulated). We could use separate filters (e.g. `bool(c4<-0.58) or bool(c4>0.58)`) or more simply, we can use the absolute (abs) value 0.58, where the minus sign is ignored.
* **Join two Datasets** to add the gene lengths information to the differentially expressed genes, matching on gene ids
* **Cut** to extract the two columns for the differentially expressed genes information
* **Cut** to extract the two columns for the gene lengths information


> ### {% icon hands_on %} Hands-on: Prepare the two inputs for GOSeq
>
> 1. **Compute** {% icon tool %} with
>    - *"Add expression"*: `bool(c8<0.01) and bool(abs(c4)>0.58)` (adj.P < 0.01 and lfc of 0.58)
>    - {% icon param-collection %} *"as a new column to"*: the `DE tables` output of **limma** {% icon tool %} (containing both the basal and luminal contrasts)
> 2. **Join two Datasets** {% icon tool %} with
>    - *"Join"*: output of **Compute** {% icon tool %}
>    - *"using column"*: `Column: 1`
>    - {% icon param-file %} *"with"* the original GEO counts file `seqdata`
>    - *"and column"*: `Column: 1`
>    - *"Keep lines of first input that do not join with second input"*: `No`
>    - *"Keep lines of first input that are incomplete"*: `No`
>    - *"Fill empty columns"*: `No`
>    - *"Keep the header lines"*: `No`
> 3. **Cut columns from a table** {% icon tool %} with
>    - *"Cut columns"*: `c1,c9` (the gene ids and DE status)
>    - *"Delimited by"*: `Tab`
>    - {% icon param-file %} *"From"*: the output of **Join** {% icon tool %}
>    - Rename to `goseq DE status`
> 4. **Cut columns from a table** {% icon tool %} with
>    - *"Cut columns"*: `c1,c11` (the gene ids and lengths)
>    - *"Delimited by"*: `Tab`
>    - {% icon param-file %} *"From"*: the output of **Join** {% icon tool %}
>    - Rename to `goseq gene lengths`
{: .hands_on}

We now have the two required input files for goseq for both our basal and luminal contrasts.

> ### {% icon hands_on %} Hands-on: Perform GO analysis
>
> 1. **goseq** {% icon tool %} with
>    - {% icon param-collection %} *"Differentially expressed genes file"*: `goseq DE status`
>    - {% icon param-file %} *"Gene lengths file"*: `goseq gene lengths`
>    - *"Gene categories"*:  `Get categories`
>       - *"Select a genome to use"*:  `Mouse(mm10)`
>       - *"Select Gene ID format"*:  `Entrez Gene ID`
>       - *"Select one or more categories"*: `GO: Biological Process`
>    - *"Output Options"*
>        - *"Output Top GO terms plot?"* `Yes`
{: .hands_on}

goseq generates a big table with the following columns for each GO term:
1. `category`: GO category
2. `over_rep_pval`: *p*-value for over representation of the term in the differentially expressed genes
3. `under_rep_pval`: *p*-value for under representation of the term in the differentially expressed genes
4. `numDEInCat`: number of differentially expressed genes in this category
5. `numInCat`: number of genes in this category
6. `term`: detail of the term
7. `ontology`: MF (Molecular Function - molecular activities of gene products), CC (Cellular Component - where gene products are active), BP (Biological Process - pathways and larger processes made up of the activities of multiple gene products)
8. `p.adjust.over_represented`: *p*-value for over representation of the term in the differentially expressed genes, adjusted for multiple testing with the Benjamini-Hochberg procedure
9. `p.adjust.under_represented`: *p*-value for over representation of the term in the differentially expressed genes, adjusted for multiple testing with the Benjamini-Hochberg procedure

To identify categories significantly enriched/unenriched below some p-value cutoff, it is necessary to use the adjusted *p*-value.

A plot of the top 10 over-represented GO terms (by adjusted *p*-value) can be output from the goseq tool to help visualise results. Note that the top 10 are selected by adjusted p-value so if there are multiple terms with the same value there will be more than 10 terms in the plot. Click on the `Top over-represented GO terms plot` in the history. There should be 2 PDFs, one for each contrast, that look similar to below.

![Basal Plot](../../images/limma-voom/basal_top_GO.png "Basal pregnant vs lactating top 10 GO terms")

![Luminal Plot](../../images/limma-voom/luminal_top_GO.png "Luminal pregnant vs lactating top 10 GO terms")

The Fu paper also used goseq and found enrichment for cell contractility genes in the basal cells and enrichment in the luminal cells for general metabolic processes, lipid biosynthesis and transport proteins, and .

> ### {% icon question %} Questions
>
> Take a look at the top 10 GO plots for the luminal and basal contrast. How do you think they compare to what the authors found?
>
> > ### {% icon solution %} Solution
> >
> > The top 10 GO terms seem to describe similar processes to what the authors found.
> >
> {: .solution}
{: .question}

## Gene Set Enrichment Analysis with **fgsea**

Gene Set Enrichment Analysis (GSEA) [(Subramanian et al., 2005)](https://www.ncbi.nlm.nih.gov/pubmed/16199517) is a widely used method that determines whether a set of genes is enriched in a list of differentially expressed genes. Unlike the previous method with goseq, no threshold is applied for what is considered "differentially expressed", all genes are used. If a gene set falls at either the top (over-expressed) or bottom (under-expressed) of the list it is said to be enriched. [fgsea](https://www.biorxiv.org/content/early/2016/06/20/060012) is a faster implementation of the GSEA method. fgsea requires a ranked list of genes and some gene sets to test.

The [Molecular Signatures Database (MSigDB)](http://software.broadinstitute.org/gsea/msigdb/index.jsp) contains curated collections of gene sets that are commonly used in a GSEA analysis. They can be downloaded from the [Broad website](http://software.broadinstitute.org/gsea/downloads.jsp). But these collections are only of human gene sets. If working with another species you would need to first map the genes to their human orthologues. However, MSigDB versions for mouse are provided [by the Smyth lab here](http://bioinf.wehi.edu.au/software/MSigDB/index.html) so we'll use those. There are several MSigDB collections, we'll use the [Hallmark collection](https://www.cell.com/cell-systems/fulltext/S2405-4712(15)00218-5), which contains 50 gene sets. According to MSigDB, "each gene set in the hallmark collection consists of a “refined” gene set, derived from multiple “founder” sets, that conveys a specific biological state or process and displays coherent expression. The hallmarks effectively summarize most of the relevant information of the original founder sets and, by reducing both variation and redundancy, provide more refined and concise inputs for gene set enrichment analysis".

There are several ways we could choose to rank our genes, we could rank by log-fold change (most upregulated to most downregulated) but that doesn't take into account any error in the log fold change value. Another way is to use the "signed fold change" which is to rank by the sign of the fold change multiplied by the P value (as described [here](http://genomespot.blogspot.com/2014/09/data-analysis-step-8-pathway-analysis.html). We could also use the t statistic that's output from limma, as that takes into account the log-fold change and it's standard error, see [here](https://support.bioconductor.org/p/6124/) for more explanation on the t statistic. We'll use the t statistic to rank here.

> ### {% icon hands_on %} Hands-on: Perform gene set enrichment with fgsea
>
> 1. Import the mouse Hallmark collection of gene sets from `http://bioinf.wehi.edu.au/software/MSigDB/mouse_H_v5p2.rdata` using the Paste/Fetch upload box
>    - Set the file **Type** to `rdata`
>    - Rename file as `mouse_hallmark_sets`
> 2. **Cut columns from a table** {% icon tool %} with
>    - *"Cut columns"*: `c1,c6` (the Entrez gene ids and t-statistic)
>    - *"Delimited by"*: `Tab`
>    - {% icon param-collection %} *"From"*: the `DE tables` output of **limma** {% icon tool %}
> 3. **Sort data in ascending or descending order** {% icon tool %} with
>    - {% icon param-collection %} *"Sort Query"*: the output of **Cut** {% icon tool %}
>    - *"Number of header lines"*: `1`
>    - *"Column selections"*:
>        - *"on column"*: `Column: 2`
>        - *"in"*: `Descending order`
>        - *"Flavor"*: `Fast numeric sort (-n)`
> 4. **fgsea** {% icon tool %} with
>    - {% icon param-collection %} *"Ranked Genes"*: the output of **Sort** {% icon tool %}
>    - *"File has header?"*: `Yes`
>    - {% icon param-file %} *"Gene Sets"*: `mouse_hallmark_sets`
>    - *"Minimum Size of Gene Set"*: `15`
>    - *"Output plots"*: `Yes`
{: .hands_on}

fgsea outputs a table of results containing a list of pathways with P values and enrichment scores. It can also output a summary table plot of the top pathways like the one shown below for the `basallpregnant-basallactate` contrast.

![fgsea Table](../../images/limma-voom/fgsea_table.png "fgsea Summary table"){: width="50%"}

An enrichment plot of the each of the top pathways can also be produced, one is shown below. The barcode pattern shows where the genes in the set are found in the list of ranked genes. Most of the bars to the left indicate enrichment of the set at the top of the ranked list of genes (upregulated) and most bars towards the right indicate enrichment at the bottom of the list (downregulated). The enrichment score reflects the degree to which the genes are enriched at the top or bottom of the list.

![fgsea Enrichment](../../images/limma-voom/fgsea_enrichplot.png "fgsea Enrichment plot")

## Ensemble gene set enrichment analyses with **EGSEA**

The ensemble of genes set enrichment analyses (EGSEA) [(Alhamdoosh et al, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/27694195) is a method developed for RNA-sequencing data that combines results from multiple algorithms and calculates collective gene set scores, to try to improve the biological relevance of the highest ranked gene sets. EGSEA has built-in gene sets from MSigDB and KEGG for human and mouse. We'll show here how it can be used with the MSigDB Hallmark collection and KEGG pathways. For input we need a count matrix and EGSEA will perform a limma-voom analysis before gene set testing. We can use the filtered counts output from limma, where the low count genes have been filtered out, we just need to remove the gene symbol and description columns. We also need a symbols mapping file containing just the Entrez ids and symbols, which we can generate from the filtered counts file. The third input we need is a factors information file, containing what groups the samples belong to, we can use the one we used with limma.

> ### {% icon hands_on %} Hands-on: Perform ensemble gene set testing with EGSEA
>
> 1. Rerun **limma** selecting *"Output Filtered Counts Table?"*: `Yes`
> 2. **Cut** {% icon tool %}: Run **Cut columns from a table (cut)** with the following settings:
>      - {% icon param-file %} *"File to cut"*: `Filtered Counts` output from **limma**
>      - *"Operation"*: `Discard`
>      - *"List of fields"*: Select `Column:2`, `Column:3`
>      - Rename to `EGSEA counts`
> 3. **Cut** {% icon tool %}: Run **Cut columns from a table (cut)** with the following settings:
>      - {% icon param-file %} *"File to cut"*: `Filtered Counts` output from **limma**
>      - *"Operation"*: `Keep`
>      - *"List of fields"*: Select `Column:1`, `Column:2`
>      - Rename to `EGSEA anno`
> 4. **EGSEA** {% icon tool %} with
>      - *"Count Files or Matrix?*": `Single Count Matrix`
>          - {% icon param-file %} *"Count Matrix"*: Select `EGSEA counts`
>      - *"Input factor information from file?"*: `Yes`
>          - {% icon param-file %} *"Factor File"*: Select `factordata`
>      - {% icon param-file %} *"Symbols Mapping file"*: `EGSEA anno`
>      - *"Contrast of Interest"*: `basalpregnant-basallactate`
>      - *"Contrast of Interest"*: `luminalpregnant-luminallactate`
>      - *"Species"*: `mouse`
>      - *"Gene Set Testing Methods"*: Tick `camera`, `safe`, `gage`, `zscore`, `gsva`, `globaltest`, `ora`, `ssgsea`, `padog`, `plage`, `fry`
>      - *"MSigDB Gene Set Collections"*: `H: hallmark gene sets`
>      - *"KEGG Pathways"*: `Metabolism` and `Signalling`
>      - *"I certify that I am not using this tool for commercial purposes"*: `Yes`
{: .hands_on}

This generates a report like below.

![EGSEA report](../../images/limma-voom/EGSEA_report.png "EGSEA report"){: width="50%"}

In addition to a table of results, plots are generated like the heatmaps of the top ranked pathways, shown below. Note that we see some similar pathways in the results here as with the fgsea analysis.

![EGSEA heatmaps](../../images/limma-voom/EGSEA_heatmaps.png "EGSEA heatmaps"){: width="50%"}

KEGG pathway diagrams are generated if KEGG pathways are selected, as shown below.  These show the expression values of the genes, genes upregulated in the contrast are shown in red, downregulated in blue.

![EGSEA KEGG](../../images/limma-voom/EGSEA_KEGG.png "EGSEA KEGG pathways")

# Conclusion
{:.no_toc}

In this tutorial we have seen some gene set testing methods that can be used to help interpret lists of differentially expressed genes. This follows on from the accompanying tutorials, RNA-seq reads to counts and RNA-seq counts to genes, that showed how to turn reads (FASTQs) into differentially expressed genes for this dataset. For further reading on analysis of RNA-seq count data and the methods used here, see the articles; RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR [(Law et al. 2016)](https://f1000research.com/articles/5-1408/v2) and From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline [(Chen, Lun, Smyth 2016)](https://f1000research.com/articles/5-1438/v2).
