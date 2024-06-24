---
layout: tutorial_hands_on

title: Multi-sample batch correction with Harmony and SnapATAC2
subtopic: scmultiomics
priority: 3
level: Intermediate
zenodo_link: ''
questions:
- Why is batch correction important when analyzing data from multiple samples?
- How is batch correction performed on single cell ATAC-seq data?
objectives:
- Perform batch correction on a collection of single cell ATAC-seq data
- Learn how Harmony integrates different samples
time_estimation: 3H
key_points:
- Batch correction is important for integration of data from multiple experiments
- How harmony works
requirements:
  -
    type: "internal"
    topic_name: single-cell
    tutorials:
      - scatac-preprocessing-tenx
      - scatac-standard-processing-snapatac2
tags:
- 10x
- epigenetics
abbreviations:
    scATAC-seq: Single-cell Assay for Transposase-Accessible Chromatin using sequencing
    QC: quality control
    TSSe: transcription start site enrichment
    TSS: transcription start sites
    UMAP: Uniform Manifold Approximation and Projection
contributors:
- timonschlegel
gitter: Galaxy-Training-Network/galaxy-single-cell


---


# Introduction

<!-- This is a comment. -->
Performing biological experiments in replicates is one of the cornerstones of modern science. However, when integrating data from multiple single-cell sequencing experiments, technical confounders might impact the results. 
To reduce technical confounders, such as different experimenters, experimental protocols, sequencing lanes or sequencing technologies, a batch correction might be beneficial. 

In this tutorial, we will perform batch correction on five datasets of {scATAC-seq} data with the algorithm *Harmony* ({% cite Korsunsky2019 %}) and the tool suite [**SnapATAC2**] (https://kzhang.org/SnapATAC2/version/2.5/index.html) ({% cite Zhang2024 %}). 

{% snippet topics/single-cell/faqs/single_cell_omics.md %}

{% snippet faqs/galaxy/tutorial_mode.md %}

> <comment-title></comment-title>
>
> This tutorial is significantly based on ["Multi-sample Pipeline" tutorial from SnapATAC2](https://kzhang.org/SnapATAC2/version/2.5/tutorials/integration.html). 
> The data analysis is performed with the same tools shown in the tutorial [Single-cell ATAC-seq standard processing with SnapATAC2]( {% link topics/single-cell/tutorials/scatac-standard-processing-snapatac2/tutorial.md %} ). 
> - That tutorial also explains the steps of the ATAC-seq analysis with SnapATAC2 in more detail. 
> - We recommend completing that tutorial before continuing with this one. 
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

# Data

The datasets for this tutorial are colon samples from multiple donors, provided by the [SnapATAC2 documentation](https://kzhang.org/SnapATAC2/version/2.5/tutorials/integration.html). 

## Get data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    {{ page.zenodo_link }}/files/colon_multisample.tar
>    {{ page.zenodo_link }}/files/chrom_sizes.txt
>    {{ page.zenodo_link }}/files/gencode.v46.annotation.gtf.gz
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets if necessary
>   
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
> 4. Check that the datatype of the `colon_multisample` files is set to `bed`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
> 5. Create a dataset collection with all `colon_multisample` datasets.
> 
>    {% snippet faqs/galaxy/collections_build_list.mdu name="Colon Multisample" %}
>
{: .hands_on}

# SnapATAC2 preprocessing and filtering

With our data imported and the collection built, we can now begin the {scATAC-seq} data preprocessing with SnapATAC2. 

The first step is importing the datasets into an AnnData object with the tool *pp.import_data*. Next, the {TSSe} will be calculated. The  {TSS} serves as a {QC} measurement to only filter droplets containing high-quality cells. 

> <hands-on-title> Preprocessing and Filtering </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Import data fragment files and compute basic QC metrics, using 'pp.import_data'`
>        - {% icon param-collection %} *"Fragment file, optionally compressed with gzip or zstd"*: `output` (Input dataset collection)
>        - {% icon param-file %} *"A tabular file containing chromosome names and sizes"*: `output` (Input dataset)
>        - *"Number of unique fragments threshold used to filter cells"*: `1000`
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

## Sub-step with **SnapATAC2 Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Filter cell outliers based on counts and numbers of genes expressed, using 'pp.filter_cells'`
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Preprocessing** {% icon tool %})
>        - *"Minimum TSS enrichemnt score required for a cell to pass filtering"*: `7.0`
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
>        - *"The size of consecutive genomic regions used to record the counts"*: `5000`
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
>        - *"Number of features to keep"*: `50000`
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

## Sub-step with **Extract element identifiers**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Extract element identifiers](toolshed.g2.bx.psu.edu/repos/iuc/collection_element_identifiers/collection_element_identifiers/0.0.2) %} with the following parameters:
>    - {% icon param-file %} *"Dataset collection"*: `anndata_out` (output of **SnapATAC2 Preprocessing** {% icon tool %})
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

## Sub-step with **Extract dataset**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Extract dataset](__EXTRACT_DATASET__) %} with the following parameters:
>    - {% icon param-file %} *"Input List"*: `anndata_out` (output of **SnapATAC2 Preprocessing** {% icon tool %})
>    - *"How should a dataset be selected?"*: `The first dataset`
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

## Sub-step with **Select first**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Select first](Show beginning1) %} with the following parameters:
>    - *"Select first"*: `1`
>    - {% icon param-file %} *"from"*: `output` (output of **Extract element identifiers** {% icon tool %})
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

## Sub-step with **Filter collection**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Filter collection](__FILTER_FROM_FILE__) %} with the following parameters:
>    - {% icon param-file %} *"Input Collection"*: `anndata_out` (output of **SnapATAC2 Preprocessing** {% icon tool %})
>    - *"How should the elements to remove be determined?"*: `Remove if identifiers are PRESENT in file`
>        - {% icon param-file %} *"Filter out identifiers present in"*: `out_file1` (output of **Select first** {% icon tool %})
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

## Sub-step with **Manipulate AnnData**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.10.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output` (output of **Extract dataset** {% icon tool %})
>    - *"Function to manipulate the object"*: `Concatenate along the observations axis`
>        - {% icon param-file %} *"Annotated data matrix to add"*: `output_filtered` (output of **Filter collection** {% icon tool %})
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
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata` (output of **Manipulate AnnData** {% icon tool %})
>        - *"Number of features to keep"*: `50000`
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

## Sub-step with **SnapATAC2 Plotting**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Plotting](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_plotting/snapatac2_plotting/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for plotting"*: `Plot the UMAP embedding, using 'pl.umap'`
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Clustering** {% icon tool %})
>        - *"Color"*: `batch`
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
>    - *"Method used for preprocessing"*: `Use harmonypy to integrate different experiments,using 'pp.harmony'`
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
>    - *"Dimension reduction and Clustering"*: `Compute Umap, using 'tl.umap'`
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Preprocessing** {% icon tool %})
>        - *"Use the indicated representation in `.obsm`"*: `X_spectral_harmony`
>        - *"`adata.obs` key under which t add cluster labels"*: `umap_harmony`
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
>        - *"Color"*: `batch`
>        - *"Use the indicated representation in .obsm"*: `X_umap_harmony`
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

## Sub-step with **SnapATAC2 Clustering**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Clustering](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_clustering/snapatac2_clustering/2.5.3+galaxy1) %} with the following parameters:
>    - *"Dimension reduction and Clustering"*: `Compute a neighborhood graph of observations, using 'pp.knn'`
>        - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **SnapATAC2 Clustering** {% icon tool %})
>        - *"The key for the matrix"*: `X_spectral_harmony`
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
>        - *"Use the indicated representation in .obsm"*: `X_umap_harmony`
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


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.