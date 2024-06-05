---
layout: tutorial_hands_on

title: Single-cell ATAC-seq standard processing with SnapATAC2
subtopic: scmultiomics
priority: 2
redirect_from:
  - /topics/transcriptomics/tutorials/scatac-standard-processing-snapatac2/tutorial
level: Intermediate
zenodo_link: https://zenodo.org/records/11369811
questions:
- What does ATAC-seq data tell us about the cell?
- Which steps are necessary for clustering of single-cell ATAC-seq data?
- Why is dimension reduction important for analysis of single-cell data?
objectives:
- Learn how ATAC-seq works
- Create a count-matrix from a 10X fragment file
- Perform filtering, dimension reduction and clustering on AnnData matrices
- Generate and filter a cell-by-gene matrix
- Identify marker genes for the clusters
time_estimation: 2H
key_points:
- Single-cell ATAC-seq can identify open chromatin-sites  
- Clustering groups similar cells together which can then be annotated to cell-types
requirements:
-
    type: "internal"
    topic_name: single-cell
    tutorials:
        - scatac-preprocessing-tenx
tags:
- 10x
- epigenetics
abbreviations:
    scATAC-seq: Single-cell Assay for Transposase-Accessible Chromatin using sequencing
    PBMC's: peripheral blood mononuclear cells

contributors:
- timonschlegel
gitter: Galaxy-Training-Network/galaxy-single-cell


---

{scATAC-seq} analysis is a method to decipher the chromatin states of the analyzed cells. In general, genes are only expressed in accessible (i.e. "open") chromatin and not in closed chromatin. 
By analyzing which genomic sites have an _open_ chromatin state, cell-type specific patterns of gene accessibility can be determined. 
{scATAC-seq} is particularly usefull for analyzing tissue containing different cell populations, such as {PBMC's}. 

In this tutorial we will analyze {scATAC-seq} data using the tool suites [SnapATAC2](https://kzhang.org/SnapATAC2/version/2.5/index.html) ({% cite Zhang2024 %}) and [Scanpy](https://scanpy.readthedocs.io/en/stable/index.html) ({%cite Wolf2018%}). 
With both of these tool suites we will perform preprocessing, clustering and identification of {scATAC-seq} datasets from [10x Genomics](https://www.10xgenomics.com/products/single-cell-atac). The analysis will be performed using a dataset of {PBMC's} containing ~4,620 single nuclei. 


<!-- This is a comment. -->


> <comment-title></comment-title>
>
> This tutorial is significantly based on ["Standard pipeline" tutorial from SnapATAC2](https://kzhang.org/SnapATAC2/version/2.5/tutorials/pbmc.html), and can be seen as the {scATAC-seq} counterpart to the scRNA-seq tutorial [Clustering 3K PBMCs with Scanpy]( {% link topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.md %} ).
>
{: .comment}


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# {scATAC-seq} with 10X Genomics
ATAC-seq utilizes a hyperactive Tn5 transposase ({% cite Kia2017 %}) to ligate adaptors to genome fragments, created by the transposase. Performing ATAC-seq on individual cells used to be an expensive and time consuming labour. The 10X Chromium NextGEM system made {scATAC-seq} a cost-effective method for gaining high-resolution data with a simple protocol. 
After transposition of nuclei in bulk, the nuclei are put into Gel beads in Emulsion (GEM), containing unique 10x cell barcodes and sequencing adaptors for Illumina sequencing. 
![Library Preparation]({% link topics/single-cell/images/scatac-pre-processing/tenx_libprep_scatac.png %} "An overview of the 10X single-nuclei ATAC-seq library preparation")


# Data

The 5k {PBMC's} dataset for this tutorial is available for free from [10X Genomics](https://www.10xgenomics.com/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-next-gem-v-1-1-1-1-standard-2-0-0). The blood samples were collected from a healthy donor and were prepared following the Chromium Next GEM scATAC-seq protocol. After sequencing on Illumina NovaSeq, the reads were processed by the **Cell Ranger ATAC 2.0.0** pipeline from 10X to generate a [*Fragments File*](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments). 

> <details-title>Fragments File </details-title>
>
>   The Fragments File is a tabular file in a [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)-like format, containing information about the position of the read on the chromosome. A Fragments File is created as part of the tutorial ["Pre-processing of 10X Single-Cell ATAC-seq Datasets"]( {% link topics/single-cell/tutorials/scatac-preprocessing-tenx/tutorial.md %} )
>
{: .details}
SnapATAC2 requires 3 input files for the standard pathway of processing:
- `fragments_file.tsv`: A tabular file containing the chromosome positions of the reads and their corresponding 10X cell barcodes. 
- `chrom_sizes.txt`: A tabular file of the number of bases of a chromosome
- `gene_annotation.gtf.gz`: A tabular file listing genomic features of the human genome (GENCODE GRCh38)

## Data upload
Section and subsection titles will be displayed in the tutorial index on the left side of
the page, so try to make them informative and concise!
> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the `fragments_file.tsv`, `chrom_sizes.txt` and `gene_annotation.gtf.gz` from [Zenodo]({{ page.zenodo_link }}) or from the shared data library
>
>
>    ```
>    {{ page.zenodo_link }}/files/atac_pbmc_5k_nextgem_fragments.tsv
>    {{ page.zenodo_link }}/files/chrom_sizes.txt
>    {{ page.zenodo_link }}/files/gencode.v46.annotation.gtf.gz
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets
>
>    > <details-title>Renaming the input datasets </details-title>
>    >- {% icon galaxy-pencil %} **Rename** the file `atac_pbmc_5k_nextgem_fragments.tsv` to `fragments_file.tsv`
>    >- {% icon galaxy-pencil %} **Rename** the file `gencode.v46.annotation.gtf.gz` to `gene_annotation.gtf.gz`
>    {: .details}
> 
> 4. Inspect `chrom_sizes` and `fragments_file` 
{: .hands_on}
> <question-title></question-title>
>
> 1. How many chromosomes are in `chrom_sizes`?
> 2. In which column are the cell barcodes stored in the `fragments_file`?
>
> > <solution-title></solution-title>
> >
> > 1. There are 25 chromosomes. The 22 autosomes (Chr. 1-22), both sex chromosomes (Chr. X and Y) and the small circular mitochondrial chromosome (Chr. M). 
> > 2. The cell barcodes are unique 16 bp oligos, located in the column `Name`.  
> >
> {: .solution}
>
{: .question}

# Preprocessing

Preprocessing of the scATAC-seq data contained in the `fragments_file` with SnapATAC2 begins with importing the files and computing basic quality control (QC) metrics. 

SnapATAC2 compresses and stores the fragments into an `AnnData` object. 

## AnnData
The [`AnnData`](https://anndata.readthedocs.io/en/latest/) format was initially developed for the [`Scanpy`](https://scanpy.readthedocs.io/en/stable/index.html) package and is now a widely accepted data format to store annotated data matrices in a space efficient manner. 

![Anndata format]({% link topics/single-cell/images/scatac-standard-snapatac2/anndata_schema.svg %} "<code>AnnData</code> format stores a count matrix <code>X</code> together with annotations of observations (i.e. cells) <code>obs</code>, variables (i.e. genes) <code>var</code> and unstructured annotations <code>uns</code>.")

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> <details-title> More details about the theory </details-title>
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


## Sub-step with **SnapATAC2 Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Import data fragment files and compute basic QC metrics, using 'pp.import_data'`
>        - {% icon param-file %} *"Fragment file, optionally compressed with gzip or zstd"*: `output` (Input dataset)
>        - {% icon param-file %} *"A tabular file containing chromosome names and sizes"*: `output` (Input dataset)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SnapATAC2 Plotting**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Plotting](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_plotting/snapatac2_plotting/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for plotting"*: `Plot fragment size distribution, using 'pl.frag_size_distr'`
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Preprocessing** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SnapATAC2 Plotting**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Plotting](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_plotting/snapatac2_plotting/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for plotting"*: `Plot fragment size distribution, using 'pl.frag_size_distr'`
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Preprocessing** {% icon tool %})
>        - *"Change the y-axis (fragment counts) to log scale"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SnapATAC2 Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Compute the TSS enrichment score (TSSe) for each cell, using 'metrics.tsse'`
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Preprocessing** {% icon tool %})
>        - {% icon param-file %} *"GTF/GFF file containing the gene annotation"*: `output` (Input dataset)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SnapATAC2 Plotting**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Plotting](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_plotting/snapatac2_plotting/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for plotting"*: `Plot the TSS enrichment vs. number of fragments density figure, using 'pl.tsse'`
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Preprocessing** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SnapATAC2 Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Filter cell outliers based on counts and numbers of genes expressed, using 'pp.filter_cells'`
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Preprocessing** {% icon tool %})
>        - *"Minimum number of counts required for a cell to pass filtering"*: `5000`
>        - *"Minimum TSS enrichemnt score required for a cell to pass filtering"*: `10.0`
>        - *"Maximum number of counts required for a cell to pass filtering"*: `100000`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SnapATAC2 Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Generate cell by bin count matrix, using 'pp.add_tile_matrix'`
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Preprocessing** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SnapATAC2 Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Perform feature selection, using 'pp.select_features'`
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Preprocessing** {% icon tool %})
>        - *"Number of features to keep"*: `250000`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SnapATAC2 Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Compute probability of being a doublet using the scrublet algorithm, using 'pp.scrublet'`
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Preprocessing** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SnapATAC2 Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Remove doublets according to the doublet probability or doublet score, using 'pp.filter_doublets'`
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Preprocessing** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SnapATAC2 Clustering**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Clustering](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_clustering/snapatac2_clustering/2.5.3+galaxy1) %} with the following parameters:
>    - *"Dimension reduction and Clustering"*: `Perform dimension reduction using Laplacian Eigenmap, using 'tl.spectral'`
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Preprocessing** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SnapATAC2 Clustering**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Clustering](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_clustering/snapatac2_clustering/2.5.3+galaxy1) %} with the following parameters:
>    - *"Dimension reduction and Clustering"*: `Compute Umap, using 'tl.umap'`
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Clustering** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SnapATAC2 Clustering**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Clustering](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_clustering/snapatac2_clustering/2.5.3+galaxy1) %} with the following parameters:
>    - *"Dimension reduction and Clustering"*: `Compute a neighborhood graph of observations, using 'pp.knn'`
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Clustering** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SnapATAC2 Clustering**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Clustering](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_clustering/snapatac2_clustering/2.5.3+galaxy1) %} with the following parameters:
>    - *"Dimension reduction and Clustering"*: `Cluster cells into subgroups, using 'tl.leiden'`
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Clustering** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SnapATAC2 Plotting**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Plotting](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_plotting/snapatac2_plotting/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for plotting"*: `Plot the UMAP embedding, using 'pl.umap'`
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Clustering** {% icon tool %})
>        - *"Color"*: `leiden`
>        - *"Height of the plot"*: `500`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SnapATAC2 Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Generate cell by gene activity matrix, using 'pp.make_gene_matrix'`
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Clustering** {% icon tool %})
>        - {% icon param-file %} *"GTF/GFF file containing the gene annotation"*: `output` (Input dataset)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Filter**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.9.6+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Preprocessing** {% icon tool %})
>    - *"Method used for filtering"*: `Filter genes based on number of cells or counts, using 'pp.filter_genes'`
>        - *"Filter"*: `Minimum number of cells expressed`
>            - *"Minimum number of cells expressed required for a gene to pass filtering"*: `5`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Normalize**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Normalize](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_normalize/scanpy_normalize/1.9.6+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Filter** {% icon tool %})
>    - *"Method used for normalization"*: `Normalize counts per cell, using 'pp.normalize_total'`
>        - *"Exclude (very) highly expressed genes for the computation of the normalization factor (size factor) for each cell"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Inspect and manipulate**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.9.6+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Normalize** {% icon tool %})
>    - *"Method used for inspecting"*: `Logarithmize the data matrix, using 'pp.log1p'`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Normalize**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Normalize](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_normalize/scanpy_normalize/1.9.6+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Inspect and manipulate** {% icon tool %})
>    - *"Method used for normalization"*: `Denoising using Markov Affinity-based Graph Imputation of Cells (MAGIC) API 'external.pp.magic'`
>        - *"Denoised genes to return"*: `PCA only`
>        - *"Which solver to use"*: `"approximate", is faster that performs imputation in the PCA space and then projects back to the gene space`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **AnnData Operations**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [AnnData Operations](toolshed.g2.bx.psu.edu/repos/ebi-gxa/anndata_ops/anndata_ops/1.9.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input object in hdf5 AnnData format"*: `anndata_out` (output of **Normalize** {% icon tool %})
>    - *"Copy embeddings (such as UMAP, tSNE)"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Plot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.9.6+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output_h5ad` (output of **AnnData Operations** {% icon tool %})
>    - *"Method used for plotting"*: `Embeddings: Scatter plot in UMAP basis, using 'pl.umap'`
>        - *"Keys for annotations of observations/cells or variables/genes"*: `{'id': 3, 'output_name': 'output'}`
>        - *"Show edges?"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
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

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
