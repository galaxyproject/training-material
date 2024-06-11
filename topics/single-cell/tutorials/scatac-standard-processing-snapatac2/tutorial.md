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
    QC: quality control
    TSSe: transcription start site enrichment
    TSS: transcription start sites
    UMAP: Uniform Manifold Approximation and Projection

contributors:
- timonschlegel
gitter: Galaxy-Training-Network/galaxy-single-cell


---

{scATAC-seq} analysis is a method to decipher the chromatin states of the analyzed cells. In general, genes are only expressed in accessible (i.e. "open") chromatin and not in closed chromatin. 
By analyzing which genomic sites have an _open_ chromatin state, cell-type specific patterns of gene accessibility can be determined. 
{scATAC-seq} is particularly usefull for analyzing tissue containing different cell populations, such as {PBMC's}. 

In this tutorial we will analyze {scATAC-seq} data using the tool suites [SnapATAC2](https://kzhang.org/SnapATAC2/version/2.5/index.html) ({% cite Zhang2024 %}) and [Scanpy](https://scanpy.readthedocs.io/en/stable/index.html) ({%cite Wolf2018%}). 
With both of these tool suites we will perform preprocessing, clustering and identification of {scATAC-seq} datasets from [10x Genomics](https://www.10xgenomics.com/products/single-cell-atac). The analysis will be performed using a dataset of {PBMC's} containing ~4,620 single nuclei. 

{% snippet topics/single-cell/faqs/single_cell_omics.md %}

{% snippet faqs/galaxy/tutorial_mode.md %}


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
After transposition of nuclei in bulk, individual nuclei are put into Gel beads in Emulsion (GEM), containing unique 10x cell barcodes and sequencing adaptors for Illumina sequencing. 
![Library Preparation]({% link topics/single-cell/images/scatac-standard-snapatac2/tenx_libprep_scatac.png %} "An overview of the 10X single-nuclei ATAC-seq library preparation")


# Data

The 5k {PBMC's} dataset for this tutorial is available for free from [10X Genomics](https://www.10xgenomics.com/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-next-gem-v-1-1-1-1-standard-2-0-0). The blood samples were collected from a healthy donor and were prepared following the Chromium Next GEM scATAC-seq protocol. After sequencing on Illumina NovaSeq, the reads were processed by the **Cell Ranger ATAC 2.0.0** pipeline from 10X to generate a [*Fragments File*](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments). 

> <details-title>Fragments File </details-title>
>
>   The Fragments File is a tabular file in a [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)-like format, containing information about the position of the read on the chromosome. 
>
{: .details}
SnapATAC2 requires 3 input files for the standard pathway of processing:
- `fragments_file.tsv`: A tabular file containing the chromosome positions of the reads and their corresponding 10X cell barcodes. 
- `chrom_sizes.txt`: A tabular file of the number of bases of the human chromosomes
- `gene_annotation.gtf.gz`: A tabular file listing genomic features of the human genome (GENCODE GRCh38)

> <comment-title></comment-title>
> - This tutorial starts with a `fragments_file`. 
> - SnapATAC2 also accepts mapped reads in a `.bam` file.
> - To learn how to get a `fragments_file` or `.bam` file from raw `.FASTQ`-reads, please check out the tutorial ["Pre-processing of 10X Single-Cell ATAC-seq Datasets"]( {% link topics/single-cell/tutorials/scatac-preprocessing-tenx/tutorial.md %} )
{: .comment}

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

Preprocessing of the scATAC-seq data contained in the `fragments_file` with SnapATAC2 begins with importing the files and computing basic {QC} metrics. 

SnapATAC2 compresses and stores the fragments into an `AnnData` object. 

## AnnData
The [`AnnData`](https://anndata.readthedocs.io/en/latest/) format was initially developed for the [`Scanpy`](https://scanpy.readthedocs.io/en/stable/index.html) package and is now a widely accepted data format to store annotated data matrices in a space efficient manner. 

![Anndata format]({% link topics/single-cell/images/scatac-standard-snapatac2/anndata_schema.svg %} "<code>AnnData</code> format stores a count matrix <code>X</code> together with annotations of observations (i.e. cells) <code>obs</code>, variables (i.e. genes) <code>var</code> and unstructured annotations <code>uns</code>.")


## Import files to SnapATAC2

> <hands-on-title> Create an AnnData object </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Import data fragment files and compute basic QC metrics, using 'pp.import_data'`
>        - {% icon param-file %} *"Fragment file, optionally compressed with gzip or zstd"*: `fragments_file.tsv` (Input dataset)
>        - {% icon param-file %} *"A tabular file containing chromosome names and sizes"*: `chrom_sizes.txt` (Input dataset)
>        - {% icon param-toggle %} *"Whether the fragment file has been sorted by cell barcodes"*: `No` 
> 
> 2. Rename the generated file to `Anndata 5k PBMC`
>
> 3. Check that the format is `h5ad`
{: .hands_on}

Because the `AnnData` format is an extension of the HDF5 format, i.e. a binary format, an `AnnData` object can not be inspected directly in Galaxy by clicking on the {% icon galaxy-eye %} (**View data**) icon. Instead we need to use a dedicated tool from the **AnnData** suite.

> <hands-on-title>Inspect an AnnData object</hands-on-title>
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.10.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC`
>    - *"What to inspect?"*: `General information about the object`
>
> 2. {% icon galaxy-eye %} Inspect the generated file
>
>    > <question-title></question-title>
>    >
>    > ```
>    > AnnData object with n_obs × n_vars = 14232 × 0
>    >  obs: 'n_fragment', 'frac_dup', 'frac_mito'
>    >  uns: 'reference_sequences'
>    >  obsm: 'fragment_paired'
>    > ```
>    >
>    > 1. How many observations are there? What do they represent?
>    > 2. How many variables are there? 
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > 1. There are 14232 observations, representing the cells.
>    > > 2. There are 0 variables, representing genomic regions. This is because genome-wide 500-bp bins are only added after initial filtering. 
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
>    > <comment-title>Faster Method for General Information</comment-title>
>    > * Many toolsets producing outputs in *AnnData* formats in Galaxy, provide the general information by default:
>    >    * Click on the name of the dataset in the history to expand it.
>    >    * General Anndata information would be given in the expanded box:
>    >
>    >      e.g.
>    >
>    >      ```
>    >      [n_obs x n_vars]
>    >      -    14232 × 0
>    >      ```
>    >    * {% icon details %} This feature isn't the most reliable and might display:
>    >      
>    >      ```
>    >      [n_obs x n_vars]
>    >      -    1 × 1
>    >      ```
>    > * In such cases and for more specific queries, {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.10.3+galaxy0) %} is required.
>    {: .comment}
>
{: .hands_on}


## Calculate and visualize {QC} metrics

> <hands-on-title> Fragment-size distribution </hands-on-title>
>
> 1. {% tool [SnapATAC2 Plotting](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_plotting/snapatac2_plotting/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for plotting"*: `Plot fragment size distribution, using 'pl.frag_size_distr'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC` (output of **pp.import_data** {% icon tool %})
> 2. {% icon galaxy-eye %} Inspect the `.png` output
> 3. {% tool [SnapATAC2 Plotting](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_plotting/snapatac2_plotting/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for plotting"*: `Plot fragment size distribution, using 'pl.frag_size_distr'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC` (output of **pp.import_data** {% icon tool %})
>        - *"Change the y-axis (fragment counts) to log scale"*: `Yes`
> 4. {% icon galaxy-eye %} Inspect the `.png` output
>
>
>  ![fragment_size_distribution]({% link topics/single-cell/images/scatac-standard-snapatac2/pl.frag_size.png %})
>
>  ![log_fragment_size_distribution]({% link topics/single-cell/images/scatac-standard-snapatac2/log_pl.frag_size.png %})
> > <question-title></question-title>
> >
> > 1. What distinct features do the plots have? And what do they represent?
> > 2. Which fragments are generally from open chromatin?
> >
> > > <solution-title></solution-title>
> > >
> > > 1. 3 peaks are clearly visible (at <100-bp, ~200-bp and ~400-bp). The smallest fragments are from nucleosome-free regions, while the larger peaks (200- and 400-bp) contain mono- and di-nucleosom fragments, respectively. 
> > > 2. The small fragments (<100-bp) are open chromatin reads, since the Tn5 transposase could easily access the loosely packed DNA ({% cite Yan2020 %}). 
> > > 
> > {: .solution}
> >
> {: .question}
{: .hands_on}


The {TSSe} is another important {QC} metric. Nucleosome-free fragments are expected to be enriched at {TSS}. TSSe shows increased fragmentation of chromatin around the TSS. This suggests open and accessible nucleosome-free chromatin. 

{TSSe} is used as a QC metric, since an increased enrichment around TSS regions suggest that the experiment has captured biological meaningful genomic features. 
TSSe scores of individual cells can be calculated using SnapATAC2's *metrics.tsse* function. 

> <hands-on-title> Calculate and Plot TSSe </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Compute the TSS enrichment score (TSSe) for each cell, using 'metrics.tsse'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC` (output of **pp.import_data** {% icon tool %})
>        - {% icon param-file %} *"GTF/GFF file containing the gene annotation"*: `gene_annotation.gtf.gz` (Input dataset)
>
> 2. Rename the generated file to `Anndata 5k PBMC TSSe` or add tag `TSSe` to the dataset:
> 
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
>
> 3. {% tool [SnapATAC2 Plotting](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_plotting/snapatac2_plotting/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for plotting"*: `Plot the TSS enrichment vs. number of fragments density figure, using 'pl.tsse'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC TSSe` (output of **metrics.tsse** {% icon tool %})
> 4. {% icon galaxy-eye %} Inspect the `.png` output
> 
>
> ![TSSe plot against number of unique fragments]({% link topics/single-cell/images/scatac-standard-snapatac2/pl.tsse.png %})
High-quality cells can be identified in the plot of {TSSe} scores against number of unique fragments for each cell. 
>
> > <question-title></question-title>
> >
> > 1. Where are high-quality cells located in the plot?
> > 2. Based on this plot, how should the filter be set?
> >
> > > <solution-title></solution-title>
> > >
> > > 1. The cells in the upper right are high-quality cells, enriched for {TSS}. Fragments in the lower left represent low-quality cells or empty droplets and should be filtered out. 
> > > 2. Setting the minimum number of counts at 5,000 and the minimum TSS enrichment to 10.0 is an adequate filter. 
> > >
> > {: .solution}
> >
> {: .question}
{: .hands_on}

## Filtering
Based on the {TSSe} plot the cells can be filtered by TSSe and fragment counts.

> <hands-on-title> Filter cells </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Filter cell outliers based on counts and numbers of genes expressed, using 'pp.filter_cells'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC TSSe` (output of **metrics.tsse** {% icon tool %})
>        - *"Minimum number of counts required for a cell to pass filtering"*: `5000`
>        - *"Minimum TSS enrichemnt score required for a cell to pass filtering"*: `10.0`
>        - *"Maximum number of counts required for a cell to pass filtering"*: `100000`
>
> 2. Rename the generated file to `Anndata 5k PBMC TSSe filtered` or add the tag `filtered` to the dataset
> 3. {% icon galaxy-eye %} Inspect the general information of the `.h5ad` output
>
>    > <question-title></question-title>
>    >
>    > ```
>    > AnnData object with n_obs × n_vars = 4564 × 0
>    >  obs: 'n_fragment', 'frac_dup', 'frac_mito', 'tsse'
>    >  uns: 'reference_sequences'
>    >  obsm: 'fragment_paired'
>    > ```
>    >
>    > 1. How have the observations changed, compared to the first `Anndata 5k PBMC` AnnData file?
>    > 2. What does this tell us about the quality of the data?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > 1. There are only 4564 observations, compared to the initial 14232 observations. 
>    > >    
>    > >    And the `obs: 'tsse'` has been added (but already during **metrics.tsse**{% icon tool %})
>    > > 2. The empty droplets and low-quality cells have been filtered out, leaving us with 4564 high-quality cells. 
>    > >
>    > {: .solution}
>    >
>    {: .question}
> 
{: .hands_on}

## Feature selection
Currently, our AnnData matrix does not contain any variables. The variables will be added in the following step with the function *pp.add_tile_matrix*. This creates a cell by bin matrix containing insertion counts across genome-wide 500-bp bins. 

After creating the variables, the most accessible features are selected. 
> <hands-on-title> Select features </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Generate cell by bin count matrix, using 'pp.add_tile_matrix'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC TSSe filtered` (output of **pp.filter_cells** {% icon tool %})
>
> 2. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Perform feature selection, using 'pp.select_features'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata tile_matrix` (output of **pp.add_tile_matrix** {% icon tool %})
>        - *"Number of features to keep"*: `250000`
>
>    > <comment-title> Select features </comment-title>
>    >
>    > - Including more features improves resolution and can reveal finer details, but it may also introduce noise. 
>    >    - To optimize results, experiment with the `n_features` parameter to find the most appropriate value for your dataset. 
>    > - At this step you can provide a blacklist or whitelist to specifically select relevant features. 
>    >    - For example the [**ENCODE Blacklist**](https://github.com/Boyle-Lab/Blacklist) ({% cite Amemiya2019 %}) can be applied here.  
>    {: .comment}
>
> 3. Rename the generated file to `Anndata 5k PBMC select_features` or add the tag `select_features` to the dataset
> 4. {% icon galaxy-eye %} Inspect the general information of the `.h5ad` output
>
>
>    > <question-title></question-title>
>    >
>    > ```
>    > AnnData object with n_obs × n_vars = 4564 × 6062095
>    >  obs: 'n_fragment', 'frac_dup', 'frac_mito', 'tsse'
>    >  var: 'count', 'selected'
>    >  uns: 'reference_sequences'
>    >  obsm: 'fragment_paired'
>    > ```
>    >
>    > 1. How did `n_vars` change compared to `Anndata 5k PBMC TSSe filtered`
>    > 2. Where are the selected features stored in the count matrix?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > 1. There are 6,062,095 variables, compared to 0 in `TSSe filtered`. 
>    > > 2. Selected features are stored in `var: 'selected'`. 
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
{: .hands_on}


## Doublet removal

Doublets are removed by calling a customized [**scrublet**](https://github.com/swolock/scrublet) algorithm. *pp.scrublet* will identify potential doublets and the function *pp.filter_doublets* removes them. 
> <hands-on-title> Scrublet </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Compute probability of being a doublet using the scrublet algorithm, using 'pp.scrublet'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC select_features` (output of **pp.select_features** {% icon tool %})
>
> 2. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Remove doublets according to the doublet probability or doublet score, using 'pp.filter_doublets'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata scrublet` (output of **pp.scrublet** {% icon tool %})
> 3. Rename the generated file to `Anndata 5k PBMC filter_doublets` or add the tag `filter_doublets` to the dataset
> 4. {% icon galaxy-eye %} Inspect the general information of the `.h5ad` output
>
>    > <question-title></question-title>
>    >
>    > ```
>    > AnnData object with n_obs × n_vars = 4430 × 6062095
>    >  obs: 'n_fragment', 'frac_dup', 'frac_mito', 'tsse', 'doublet_probability', 'doublet_score'
>    >  var: 'count', 'selected'
>    >  uns: 'doublet_rate', 'reference_sequences', 'scrublet_sim_doublet_score'
>    >  obsm: 'fragment_paired'
>    > ```
>    >
>    > 1. What was removed by **pp.filter_doublets**?
>    > 2. Where are the new annotations stored?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > 1. Cell doublets were removed. `n_obs` was reduced from 4564 to 4430 cells. 
>    > > 2. The outputs of **pp.scrublet** are stored in observations `obs: 'doublet_probability', 'doublet_score'` and in unstructured annotations `uns: 'scrublet_sim_doublet_score'`. The output of **pp.filter_doublets** is stored in `uns: 'doublet_rate'`. 
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
>
{: .hands_on}


# Dimension reduction

Dimension reduction (also known as embedding) is a very important step during the analysis of single cell data. During this, the complex multi-dimensional data is projected into lower-dimensional space, while retaining as much information as possible. Dimension reduction enables quicker downstream analysis, since the data is more simplified and thus the memory usage is reduced. 

> <details-title>Dimension reduction with SnapATAC2</details-title>
>
> - Dimension reduction algorithms can be either linear or non-linear. 
> - Linear methods are computationally efficient and well scalable. 
>   A popular linear dimension reduction algorithm is: 
>     - **PCA** (Principle Component Analysis), implemented in **Scanpy** (please check out our [Scanpy]({% link topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.md %}) tutorial for an explanation). 
> - Nonlinear methods however are well suited for multimodal and complex datasets. 
>     - As such, they are implemented in many algorithms to visualize the data in 2 dimensions (f.ex. **UMAP** embedding).
> - The nonlinear dimension reduction algorithm, through *spectral embedding*, used in SnapATAC2 is currently the fastest and most memory efficient non-linear algorithm available ({% cite Zhang2024%}). 
>     - **Spectral embedding** utilizes an iterative algorithm to calculate the **spectrum** (*eigenvalues* and *eigenvectors*) of a matrix without computing the matrix itself. 
{: .details}

## Spectral embedding
The dimension reduction, produced by the algorithm *tl.spectral*, is required for later steps, such as plotting and clustering. 

> <hands-on-title> Spectral embedding </hands-on-title>
>
> 1. {% tool [SnapATAC2 Clustering](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_clustering/snapatac2_clustering/2.5.3+galaxy1) %} with the following parameters:
>    - *"Dimension reduction and Clustering"*: `Perform dimension reduction using Laplacian Eigenmap, using 'tl.spectral'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC filter_doublets` (output of **pp.filter_doublets** {% icon tool %})
>        - *"Distance metric"*: `cosine` 
>
> 2. Rename the generated file to `Anndata 5k PBMC spectral` or add the tag `spectral` to the dataset
> 3. {% icon galaxy-eye %} Inspect the general information of the `.h5ad` output
>
>    > <question-title></question-title>
>    >
>    > ```
>    > AnnData object with n_obs × n_vars = 4430 × 6062095
>    >  obs: 'n_fragment', 'frac_dup', 'frac_mito', 'tsse', 'doublet_probability', 'doublet_score'
>    >  var: 'count', 'selected'
>    >  uns: 'doublet_rate', 'reference_sequences', 'scrublet_sim_doublet_score', 'spectral_eigenvalue'
>    >  obsm: 'fragment_paired', 'X_spectral'
>    > ```
>    >
>    > Where are the new annotations stored?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > The outputs of **tl.spectral** are stored in unstructured annotations `uns: 'spectral_eigenvalue'` and as multidimensional observations `obsm: 'X_spectral'`. 
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
{: .hands_on}
## UMAP embedding
With the already reduced dimensionality of the data stored in `X_spectral`, the cells can be further embedded (i.e. transformed into lower dimensions) with {UMAP}. UMAP projects the cells and their relationship to each other into 2-dimensional space, which can be easily visualized. 

> <hands-on-title> UMAP embedding </hands-on-title>
>
> 1. {% tool [SnapATAC2 Clustering](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_clustering/snapatac2_clustering/2.5.3+galaxy1) %} with the following parameters:
>    - *"Dimension reduction and Clustering"*: `Compute Umap, using 'tl.umap'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC spectral` (output of **tl.spectral** {% icon tool %})
>
> 2. Rename the generated file to `Anndata 5k PBMC UMAP` or add the tag `UMAP` to the dataset
{: .hands_on}

# Clustering
During clustering, cells that share similar accessibility profiles are organized into clusters. **SnapATAC2** utilizes graph-based community clustering with the *Leiden* method. This method takes the results of k-nearest neighbor (KNN) method as input data and produces well-connected communities. 

## Community clustering

> <hands-on-title> Clustering analysis </hands-on-title>
>
> 1. {% tool [SnapATAC2 Clustering](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_clustering/snapatac2_clustering/2.5.3+galaxy1) %} with the following parameters:
>    - *"Dimension reduction and Clustering"*: `Compute a neighborhood graph of observations, using 'pp.knn'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC umap` (output of **tl.umap** {% icon tool %})
>
> 2. {% tool [SnapATAC2 Clustering](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_clustering/snapatac2_clustering/2.5.3+galaxy1) %} with the following parameters:
>    - *"Dimension reduction and Clustering"*: `Cluster cells into subgroups, using 'tl.leiden'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata knn` (output of **pp.knn** {% icon tool %})
>        - *"Whether to use the Constant Potts Model (CPM) or modularity"*: `modularity` 
> 
>    > <comment-title> CPM or modularity </comment-title>
>    > - make sure you selected `modularity`
>    > - the clusters produced by `CPM` are not represented well in the UMAP projections
>    > 
>    {: .comment}
> 2. Rename the generated file to `Anndata 5k PBMC leiden` or add the tag `leiden` to the dataset
> 3. {% icon galaxy-eye %} Inspect the general information of the `.h5ad` output
>
>    > <question-title></question-title>
>    >
>    > ```
>    > AnnData object with n_obs × n_vars = 4564 × 6062095
>    >  obs: 'n_fragment', 'frac_dup', 'frac_mito', 'tsse', 'doublet_probability', 'doublet_score', 'leiden'
>    >  var: 'count', 'selected'
>    >  uns: 'doublet_rate', 'reference_sequences', 'scrublet_sim_doublet_score', 'spectral_eigenvalue'
>    >  obsm: 'X_spectral', 'X_umap', 'fragment_paired'
>    >  obsp: 'distances'
>    > ```
>    >
>    > Where are the **leiden** clusters stored in the AnnData?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > The clusters are stored in `obs: 'leiden'`
>    > >
>    > {: .solution}
>    >
>    {: .question}
> 
{: .hands_on}


## Plotting of clusters

> <hands-on-title> Plotting the clusters </hands-on-title>
>
> 1. {% tool [SnapATAC2 Plotting](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_plotting/snapatac2_plotting/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for plotting"*: `Plot the UMAP embedding, using 'pl.umap'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC leiden` (output of **tl.leiden** {% icon tool %})
>        - *"Color"*: `leiden`
>        - *"Height of the plot"*: `500`
> 4. {% icon galaxy-eye %} Inspect the `.png` output
>
>  ![umap_leiden_clustering]({% link topics/single-cell/images/scatac-standard-snapatac2/pl.umap.png %})
>
> > <question-title></question-title>
> >
> > 1. How many leiden clusters were discovered?
> > 2. What does the distance of clusters to each other tell us about their chromatin states?
> >
> > > <solution-title></solution-title>
> > >
> > > 1. There are 12 leiden clusters. 
> > > 2. Clusters in close proximity (f.ex. clusters 0 and 5) share a similar chromatin accessibility profile, compared to a cluster further away (f.ex. cluster 9). 
> > > 
> > {: .solution}
> >
> {: .question}
{: .hands_on}

# Cell cluster annotation

After clustering the cells, they must be annotated. This categorizes the clusters into known cell types. **Manual Cell Annotation** requires known marker genes and varying expression profiles of the marker genes among clusters. As marker genes for {PBMC's} are [known](https://panglaodb.se/markers.html) ({% cite Franzn2019 %}), we can annotate our clusters manually. 

## Gene matrix
Since our data currently doesn't contain gene information, we have to create a cell by gene activity matrix using the function *pp.make_gene_matrix*. 


> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.5.3+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Generate cell by gene activity matrix, using 'pp.make_gene_matrix'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC leiden` (output of **tl.leiden** {% icon tool %})
>        - {% icon param-file %} *"GTF/GFF file containing the gene annotation"*: `gene_annotation.gtf.gz` (Input dataset)
> 2. Rename the generated file to `Anndata 5k PBMC gene_matrix` or add the tag `gene_matrix` to the dataset
>    > <tip-title> Gene matrix </tip-title>
>    >
>    > - Please note that *pp.make_gene_matrix* removes all annotations except those stored in `n_obs`. 
>    > - Therefore it might be necessary to remove propagating tags (tags starting with `#`) from `Anndata 5k PBMC gene_matrix`. 
>    >    - Tags can be removed by expanding the dataset with a tag and clicking the `x` next to the tag.
>    {: .tip}
>
> 3. {% icon galaxy-eye %} Inspect the general information of the `.h5ad` output
>
>    > <question-title></question-title>
>    >
>    > ```
>    > AnnData object with n_obs × n_vars = 4430 × 60606
>    >  obs: obs: 'n_fragment', 'frac_dup', 'frac_mito', 'tsse', 'doublet_probability', 'doublet_score', 'leiden'
>    > ```
>    >
>    > What does `n_vars` represent in `Anndata 5k PBMC gene_matrix` and what did it represent in `Anndata 5k PBMC leiden`?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > - The variables now represent accessible genes. There are `60606` accessible genes in our samples. In `Anndata 5k PBMC leiden` and all earlier AnnData the variables represented fixed-sized genomic bins. 
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
{: .hands_on}


## Imputation with scanpy
Similar to scRNA-seq data, the cell by gene activity matrix is very sparse. Additionally, high gene variance between cells, due to technical confounders, could impact the downstream analysis. In scRNA-seq, filtering and normalization are therefore required to produce a high-quality gene matrix. 

Since the cell by gene activity matrix resembles the cell by gene expression matrix of scRNA-seq, we can use the tools of the [Scanpy](https://scanpy.readthedocs.io/en/stable/index.html) ({%cite Wolf2018%}) toolsuite to continue with our data. 

> <hands-on-title> Filter and normalize </hands-on-title>
>
> 1. {% tool [Filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.9.6+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC gene_matrix` (output of **pp.make_gene_matrix** {% icon tool %})
>    - *"Method used for filtering"*: `Filter genes based on number of cells or counts, using 'pp.filter_genes'`
>        - {% icon param-select %} *"Filter"*: `Minimum number of cells expressed`
>            - *"Minimum number of cells expressed required for a gene to pass filtering"*: `5`
> 
> 2. {% tool [Normalize](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_normalize/scanpy_normalize/1.9.6+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Anndata filter_genes` (output of **pp.filter_genes** {% icon tool %})
>    - *"Method used for normalization"*: `Normalize counts per cell, using 'pp.normalize_total'`
>    - {% icon param-toggle %} *"Exclude (very) highly expressed genes for the computation of the normalization factor (size factor) for each cell"*: `No`
>
> 3. {% tool [Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.9.6+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Anndata normalize` (output of **Normalize** {% icon tool %})
>    - *"Method used for inspecting"*: `Logarithmize the data matrix, using 'pp.log1p'`
> 
> 4. {% tool [Normalize](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_normalize/scanpy_normalize/1.9.6+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Anndata log1p` (output of **log1p** {% icon tool %})
>    - *"Method used for normalization"*: `Denoising using Markov Affinity-based Graph Imputation of Cells (MAGIC) API 'external.pp.magic'`
>        - *"Denoised genes to return"*: `PCA only`
>        - *"Which solver to use"*: `"approximate", is faster that performs imputation in the PCA space and then projects back to the gene space`
>
>    > <warning-title>Large output files!</warning-title>
>     >   - The settings are important for this step!
>     >   - The setting `Denoised genes to return: 'All genes'` produces a very large output file
>     >      - **57.4 GB** compared to **1.5 GB** with `PCA only`
>     >   - The compute time for `Which solver to use: 'exact'` is very long and not necessary for our purposes. 
>     >
>     {: .warning}
>
> 5. Rename the generated file to `Anndata 5k PBMC gene_matrix magic` or add the tag `magic` to the dataset
>
> 6. {% icon galaxy-eye %} Inspect the general information of the `.h5ad` output
>
>
>    > <question-title></question-title>
>    >
>    > ```
>    > AnnData object with n_obs × n_vars = 4430 × 55106
>    >  obs: obs: 'n_fragment', 'frac_dup', 'frac_mito', 'tsse', 'doublet_probability', 'doublet_score', 'leiden'
>    >  var: 'n_cells'
>    >  uns: 'reference_sequences'
>    >  obsm: 'log1p'
>    > ```
>    >
>    > 1. How did `n_vars` change, compared to `Anndata 5k PBMC gene_matrix`?
>    > 2. Which data was 'lost', compared to `Anndata 5k PBMC leiden`, and must be added to the file in order to produce **UMAP** plots?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > 1. The number of accessible genes was reduced from 60606 to 55106 by the filtering. And additional annotations were added, such as: `var: 'n_cells'`, `uns: 'reference_sequences'` and `obsm: 'log1p'`
>    > > 2. The UMAP embeddings `obsm: 'X_umap'` are missing and should be added to the Anndata in the next step. 
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
{: .hands_on}

## Copy-over embeddings

> <hands-on-title> Copy UMAP embedding </hands-on-title>
>
> 1. {% tool [AnnData Operations](toolshed.g2.bx.psu.edu/repos/ebi-gxa/anndata_ops/anndata_ops/1.9.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input object in hdf5 AnnData format"*: `Anndata 5k PBMC gene_matrix magic` (output of **external.pp.magic** {% icon tool %})
>    - *"Copy embeddings (such as UMAP, tSNE)"*: `Yes`
>       - *"Keys from embeddings to copy"*: `X_umap`
>       - {% icon param-file %} *"IAnnData objects with embeddings to copy"*: `Anndata 5k PBMC leiden`
>
>    > <comment-title> Annotations to copy </comment-title>
>    >
>    > - This tutorial only focuses on producing an **UMAP** plot with marker-genes. 
>    > - If further analysis, with tools requiring more annotations, is intended, these can be added in a similar way as shown above.  
>    >     - f.ex. *Peak and Motif Analysis* with {% tool [Snapatac2 peaks and motif](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_peaks_and_motif/snapatac2_peaks_and_motif/2.5.3+galaxy1) %} requires annotations from `uns`. 
>    > - It is also possible to leave the input *"Keys from embeddings to copy"* empty, to copy all `obsm`. 
>    {: .comment}
>
> 5. Rename the generated file to `Anndata 5k PBMC gene_matrix magic UMAP` or add the tag `UMAP` to the dataset
>
> 6. {% icon galaxy-eye %} Inspect the general information of the `.h5ad` output, to check if `obsm` contains `X_umap`
>  
>   ```
>   AnnData object with n_obs × n_vars = 4430 × 55106
>    obs: obs: 'n_fragment', 'frac_dup', 'frac_mito', 'tsse', 'doublet_probability', 'doublet_score', 'leiden'
>    var: 'n_cells'
>    uns: 'reference_sequences'
>    obsm: 'log1p'
>   ```
>
{: .hands_on}

## Visualize gene activity of marker genes
The gene activity of selected marker genes can now be visualized with Scanpy. 

> <hands-on-title> Plot marker genes </hands-on-title>
>
> 1. {% tool [Plot with Scanpy](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.9.6+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output_h5ad` (output of **AnnData Operations** {% icon tool %})
>    - *"Method used for plotting"*: `Embeddings: Scatter plot in UMAP basis, using 'pl.umap'`
>        - *"Keys for annotations of observations/cells or variables/genes"*: `leiden, MS4A1, CD3D, LEF1, NKG7, TREM1, LYZ, PPBP`
>        - {% icon param-toggle %} *"Show edges?"*: `No`
>
>  ![umap_leiden_marker_gene_clustering]({% link topics/single-cell/images/scatac-standard-snapatac2/pl.umap.png %})
>
> > <question-title></question-title>
> >
> > 1. Are the marker genes selectively expressed in the clusters?
> > 2. Which marker genes have overlapping expression profiles? And what could that imply?
> >
> > > <solution-title></solution-title>
> > >
> > > 1. Answer 1
> > > 2. Answer 2
> > > 
> > {: .solution}
> >
> {: .question}
{: .hands_on}

## Manual cluster annotation
Comparison of marker gene expression in our clusters with a table of canonical marker genes, enables us to annotate the clusters manually. 

Cell type | Marker genes
--- | ---
memory T cells | CD3D
naive T cells | LEF1
Monocytes | LYZ
B cells | MS4A1
Natural killer (NK) cells | NKG7
Dendritic Cells | TREM1
Megakaryocytes | PPBP

These canonical marker genes can easily match the clusters to known cell types:

Cluster | Cell type
--- | ---
0 | Dendritic cells
1 | memory T cells
2 | memory T cells
3 | Dendritic cells
4 | naive T cells
5 | Monocytes
6 | B cells
7 | naive T cells
8 | naive T cells
9 | NK cells
10 | Dendritic cells
11 | B cells
12 | Megakaryocytes

> <comment-title></comment-title>
> Note that some clusters contain subtypes (f.ex. the annotated T cell clusters contain both CD4+ and CD8+ T cells). The cell-type annotation can be refined by choosing more specific marker genes. 
Hands-on: manually annotate the clusters
{: .comment}

To manually annotate the *Leiden* clusters, we will need to perform multiple steps: 
    1. **Inspect** the key-indexed observations of `Anndata 5k PBMC gene_matrix magic UMAP` 
    2. **Cut** the *Leiden* annotations out of the table
    3. Make a *replace file* containing the new cell type annotations for the *Leiden* clusters
    4. **Replace** the values of the cluster annotation with cell type annotation
    5. **Add** the cell type annotation to the AnnData
    6. **Plot** the annotated cell types

> <hands-on-title> Manual annotation </hands-on-title>
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.10.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC gene_matrix magic UMAP`
>    - *"What to inspect?"*: `Key-indexed observations annotation`
> 2. {% icon galaxy-eye %} Inspect the generated file
>
>    > <question-title></question-title>
>    > In which column is the `Leiden` annotation located?
>    > > <solution-title></solution-title>
>    > > The `Leiden` annotation is in column 8. 
>    > > ```
>    > > Column 1 | Column 2 | Column 3 | Column 4 | Column 5 | Column 6 | Column 7 | Column 8 | Column 9 | Column 10
>    > > --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
>    > > "" | n_fragment | frac_dup | frac_mito | tsse | doublet_probability | doublet_score | leiden | n_genes | n_counts
>    > > AAACGAAAGACGTCAG-1 | 22070 | 0.5219425551271499 | 0.0 | 30.43315066436454 | 0.004634433324822066 | 0.009276437847866418 | 8 | 52303 | 16521.599844068267
>    > > AAACGAAAGATTGACA-1 | 10500 | 0.5345125681606597 | 0.0 | 29.10551296093465 | 0.004668403569267374 | 0.001088139281828074 | 1 | 54501 | 15020.42495602328
>    > > AAACGAAAGGGTCCCT-1 | 19201 | 0.5101785714285714 | 0.0 | 19.90011850347046 | 0.004634433324822066 | 0.009276437847866418 | 5 | 54212 | 16294.751533305309
>    > > AAACGAACAATTGTGC-1 | 13242 | 0.487399837417257 | 0.0 | 29.060913705583758 | 0.004660125753854076 | 0.0022172949002217295 | 7 | 53530 | 15456.629863655084
>    > > ``` 
>    > {: .solution}
>    >
>    {: .question}
>
> 3. Rename the generated file to `5k PBMC observations` or add the tag `obs` to the dataset
>
> {% snippet  faqs/galaxy/analysis_cut.md %}
> 4. {% tool [Cut columns](toolshed.g2.bx.psu.edu/repos/devteam/cut_columns/Cut1/1.0.2) %} with the following parameters:
>    - {% icon param-select %} *"Cut columns"*: `c8`
>    - {% icon param-file %} *"From"*: `5k PBMC observations` (output of **Inspect AnnData** {% icon tool %})
>
> 5. Create a new **.csv** file from the following
>    ```
>    leiden, cell_type
>    0, Dendritic_cells
>    1, memory_Tcells
>    2, memory_Tcells
>    3, Dendritic_cells
>    4, naive_Tcells
>    5, Monocytes
>    6, Bcells
>    7, naive_Tcells
>    8, naive_Tcells
>    9, NKcells
>    10, Dendritic_cells
>    11, Bcells
>    12, Megakaryocytes
>    ```
>    {% snippet faqs/galaxy/datasets_create_new_file.md format="csv" name="replace_file"%}
> 
>    > <details-title>Replace file</details-title>
>    >
>    > - The first column of the replace file contains the "old" annotations and the second column contains the "new" annotation. 
>    > - {% icon warning %} Spaces between entries can lead to errors. Please use underscores (`_`) instead. 
>    >
>    {: .details}
> 6. {% tool [Replace column](toolshed.g2.bx.psu.edu/repos/bgruening/replace_column_by_key_value_file/replace_column_with_key_value_file/0.2) %} with the following parameters:
>    - {% icon param-file %} *"File in which you want to replace some values"*: `Cut columns leiden` (output of **Cut columns** {% icon tool %})
>    - {% icon param-file %} *"Replace information file"*: `replace_file` 
>    - *"Which column should be replaced?*: `Column: 1`
> 7. {% icon galaxy-eye %} Inspect the generated file to check if the replacement was successful
> 8. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.10.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC gene_matrix magic UMAP`
>    - *"Function to manipulate the object"*: `Add new annotation(s) for observations or variables`
>        - *"What to annotate"*: `Observations (obs)`
>        - {% icon param-file %} *"Table with new annotations"*: `Replace column`(output of **Replace column** {% icon tool %})
> 9. Rename the generated file to `Anndata 5k PBMC gene_matrix magic cell_type` or add the tag `cell_type` to the dataset
>
> 10. {% tool [Plot with Scanpy](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.9.6+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC gene_matrix magic cell_type` (output of **Replace column** {% icon tool %})
>    - *"Method used for plotting"*: `Embeddings: Scatter plot in UMAP basis, using 'pl.umap'`
>        - *"Keys for annotations of observations/cells or variables/genes"*: `cell_type`
>        - {% icon param-toggle %} *"Show edges?"*: `No`
>        - In *"Plot attributes"*
>           - *"Location of legend"*: `on data`
>           - {% icon param-toggle %} *"Draw a frame around the scatter plot?"*: `No`
> 6. {% icon galaxy-eye %} Inspect the `.png` output
>
>
>    > <question-title></question-title>
>    > 1. Question
>    > > <solution-title></solution-title>
>    > >
>    > > 1. solution
>    > {: .solution}
>    >
>    {: .question}
>
{: .hands_on}


# Conclusion
{% icon congratulations %} Well done, you’ve made it to the end! You might want to consult your results with this [control history](https://singlecell.usegalaxy.eu/u/timn/h/test-of-5k-pbmc-tutorial-workflow), or check out the [full workflow](https://singlecell.usegalaxy.eu/u/timn/w/2combined-snapatac2) for this tutorial.

In this tutorial we produced a count matrix of {scATAC-seq} reads in the `AnnData` format and performed: 
1. Preprocessing: 
   1. Plotting the fragment-size distributions
   2. Calculating and plotting {TSSe} scores
   3. Filtering cells and selecting features (fixed-size genomic bins)
2. Dimension reduction through **Spectral embedding** and **{UMAP}**
3. Clustering of cells via the *Leiden* method
4. Cluster annotation
   1. Producing and filtering a cell by gene activity matrix
   2. Data normalization and imputation with **Scanpy**
   3. Visualizing marker genes in the clusters
   4. Manually annotating the cell types

