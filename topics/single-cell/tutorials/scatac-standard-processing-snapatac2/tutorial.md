---
layout: tutorial_hands_on

title: Single-cell ATAC-seq standard processing with SnapATAC2
subtopic: scmultiomics
priority: 2
redirect_from:
  - /topics/transcriptomics/tutorials/scatac-standard-processing-snapatac2/tutorial
level: Intermediate
zenodo_link: https://zenodo.org/records/12707159
questions:
- What does ATAC-seq data tell us about the cell?
- Which steps are necessary to cluster the cells of single-cell ATAC-seq data?
- Why is dimension reduction important for the analysis of single-cell data?
objectives:
- Learn how single-cell ATAC-seq data is processed
- Create a count-matrix from a 10X fragment file
- Perform filtering, dimension reduction and clustering on AnnData matrices
- Generate and filter a cell-by-gene matrix
- Identify marker genes for the clusters and annotate the cell types
time_estimation: 4H
key_points:
- Single-cell ATAC-seq can identify open chromatin sites
- Dimension reduction is required to simplify the data while preserving important information about the relationships of cells to each other.
- Clusters of similar cells can be annotated to their respective cell-types
requirements:
  -
    type: "internal"
    topic_name: single-cell
    tutorials:
      - scatac-preprocessing-tenx
tags:
- 10x
- epigenetics
- single-cell
abbreviations:
    scATAC-seq: Single-cell Assay for Transposase-Accessible Chromatin using sequencing
    PBMCs: peripheral blood mononuclear cells
    QC: quality control
    TSSe: transcription start site enrichment
    TSS: transcription start sites
    UMAP: Uniform Manifold Approximation and Projection
contributions:
  authorship:
    - timonschlegel
  editing:
    - pavanvidem
    - bgruening
  testing:
    - pavanvidem
gitter: Galaxy-Training-Network/galaxy-single-cell


---

{scATAC-seq} analysis is a method to decipher the chromatin states of the analyzed cells. In general, genes are only expressed in accessible (i.e. "open") chromatin and not in closed chromatin.
By analyzing which genomic sites have an _open_ chromatin state, cell-type specific patterns of gene accessibility can be determined.
{scATAC-seq} is particularly useful for analyzing tissue containing different cell populations, such as {PBMCs}.

In this tutorial we will analyze {scATAC-seq} data using the tool suites [SnapATAC2](https://kzhang.org/SnapATAC2/version/2.5/index.html) ({% cite Zhang2024 %}) and [Scanpy](https://scanpy.readthedocs.io/en/stable/index.html) ({%cite Wolf2018%}).
With both of these tool suites we will perform preprocessing, clustering and identification of {scATAC-seq} datasets from [10x Genomics](https://www.10xgenomics.com/products/single-cell-atac).
The analysis will be performed using a dataset of {PBMCs} containing ~4,620 single nuclei.

{% snippet topics/single-cell/faqs/single_cell_omics.md %}

{% snippet faqs/galaxy/tutorial_mode.md %}


<!-- This is a comment. -->

> <comment-title></comment-title>
>
> This tutorial is significantly based on ["Standard pipeline" tutorial from SnapATAC2](https://kzhang.org/SnapATAC2/version/2.5/tutorials/pbmc.html). It can be seen as the {scATAC-seq} counterpart to the scRNA-seq tutorial [Clustering 3K PBMCs with Scanpy]( {% link topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.md %} ).
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
ATAC-seq utilizes a hyperactive Tn5 transposase ({% cite Kia2017 %}) to ligate adaptors to genome fragments, created by the transposase. Performing ATAC-seq on individual cells used to be an expensive and time-consuming labour.
The 10X Chromium NextGEM system made {scATAC-seq} a cost-effective method for gaining high-resolution data with a simple protocol.
After the transposition of nuclei in bulk, individual nuclei are put into Gel beads in Emulsion (GEM), containing unique 10x cell barcodes and sequencing adaptors for Illumina sequencing.
![Library Preparation]({% link topics/single-cell/images/scatac-standard-snapatac2/tenx_libprep_scatac.png %} "An overview of the 10X single-nuclei ATAC-seq library preparation")


# Data

The 5k {PBMCs} dataset for this tutorial is available for free from [10X Genomics](https://www.10xgenomics.com/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-next-gem-v-1-1-1-1-standard-2-0-0). The blood samples were collected from a healthy donor and were prepared following the Chromium Next GEM scATAC-seq protocol. After sequencing on Illumina NovaSeq, the reads were processed by the **Cell Ranger ATAC 2.0.0** pipeline from 10X to generate a [*Fragments File*](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments).

> <details-title>Fragments File </details-title>
>
>   The Fragments File is a tabular file in a [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)-like format, containing information about the position of the fragments on the chromosome and their corresponding 10x cell barcodes.
>
{: .details}
SnapATAC2 requires 3 input files for the standard pathway of processing:
- `fragments_file.tsv`: A tabular file containing the chromosome positions of the fragments and their corresponding 10x cell barcodes.
- `chrom_sizes.txt`: A tabular file of the number of bases of each chromosome of the reference genome
- `gene_annotation.gtf.gz`: A tabular file listing genomic features of the reference genome (GENCODE GRCh38)

> <details-title>Chromosome sizes </details-title>
>
> - A chromosome sizes file can be generated using the tool {% tool [Compute sequence length](toolshed.g2.bx.psu.edu/repos/devteam/fasta_compute_length/fasta_compute_length/1.0.3) %}.
> - The reference genome can either be selected from cached genomes or uploaded to the galaxy history.
>
{: .details}

> <comment-title></comment-title>
> - This tutorial starts with a `fragment` file.
> - SnapATAC2 also accepts mapped reads in a `BAM` file.
> - To learn how to get a `fragment` file or `BAM` file from raw `.FASTQ`-reads, please check out the tutorial ["Pre-processing of 10X Single-Cell ATAC-seq Datasets"]( {% link topics/single-cell/tutorials/scatac-preprocessing-tenx/tutorial.md %} )
> - If you would like to start the analysis with a `BAM` file, you can expand the details section ["Details: Creating a fragment file"]( {% link topics/single-cell/tutorials/scatac-standard-processing-snapatac2/tutorial.md %}#creating-a-fragment-file).
{: .comment}



## Get Data
> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the `fragments_file.tsv`, `chrom_sizes.txt` and `gene_annotation.gtf.gz` from [Zenodo]({{ page.zenodo_link }}) or from the shared data library
>
>
>    ```
>    {{ page.zenodo_link }}/files/atac_pbmc_5k_nextgem_fragments.tsv.gz
>    {{ page.zenodo_link }}/files/chrom_sizes.txt
>    {{ page.zenodo_link }}/files/gencode.v46.annotation.gtf.gz
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets
>   - {% icon galaxy-pencil %} **Rename** the file `atac_pbmc_5k_nextgem_fragments.tsv` to `fragments_file.tsv`
>   - {% icon galaxy-pencil %} **Rename** the file `gencode.v46.annotation.gtf.gz` to `gene_annotation.gtf.gz`
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
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
> > 2. The cell barcodes are unique 16 bp oligos, located in column 4.
> >
> {: .solution}
>
{: .question}

## Creating a fragment file
> <details-title>Creating a fragment file</details-title>
>  > <hands-on-title>fragment file</hands-on-title>
>  > 1. Import the file  `BAM_500-PBMC` from [Zenodo]({{ page.zenodo_link }}) or from the shared data library
>  > ```
>  > {{ page.zenodo_link }}/files/BAM_500-PBMC.bam
>  > ```
>  >   - This dataset contains mapped reads in the `BAM` format.
>  >   - It was generated by following the tutorial ["Pre-processing of 10X Single-Cell ATAC-seq Datasets"]( {% link topics/single-cell/tutorials/scatac-preprocessing-tenx/tutorial.md %} ) until the output of {% tool [Map with BWA-MEM](toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.18) %}
>  > 2. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.6.4+galaxy1) %} with the following parameters:
>  >    - *"Method used for preprocessing"*: `Convert a BAM file to a fragment file, using 'pp.make_fragment_file'`
>  >        - {% icon param-file %} *"File name of the BAM file"*: `BAM_500-PBMC` (Input dataset)
>  >        - {% icon param-toggle %} *"Indicate whether the BAM file contain paired-end reads"*: `Yes`
>  >        - *"How to extract barcodes from BAM records?"*: `From read names using regular expressions`
>  >          - *"Extract barcodes from read names of BAM records using regular expressions"*: `(................):`
>  >
>  >    > <comment-title></comment-title>
>  >    > - Not every regular expression type is supported.
>  >    > - This expression selects 16 characters if they are followed by a colon. Only the cell barcodes of the `BAM` file will match.
>  >    {: .comment}
>  >
>  > 3. Rename the generated file to `Fragments 500 PBMC`
>  > 4. Now you can continue with either the `fragments_file` from earlier or the new file `Fragments 500 PBMC`.
>  >    - {% icon galaxy-info %} The tool `pp.make_fragment_file` {% icon tool %} has implemented additional {QC} measures. These filter out larger fragments, which will be noticeable in the log-scale fragment size distribution.
>  >    - {% icon galaxy-info %} Please note that `Fragments 500 PBMC` only contains 500 {PBMCs} and thus the clustering will produce different outputs compared to the outputs generated by `fragments_file` (with 5k PBMC).
>  {: .hands_on}
{: .details}

# Preprocessing

Preprocessing of the scATAC-seq data contained in the `fragment` file with SnapATAC2 begins with importing the files and computing basic {QC} metrics.

SnapATAC2 compresses and stores the fragments into an `AnnData` object.

## AnnData
The [`AnnData`](https://anndata.readthedocs.io/en/latest/) format was initially developed for the [`Scanpy`](https://scanpy.readthedocs.io/en/stable/index.html) package and is now a widely accepted data format to
store annotated data matrices in a space-efficient manner.

![Anndata format]({% link topics/single-cell/images/scatac-standard-snapatac2/anndata_schema.svg %} "<code>AnnData</code> format stores a count matrix <code>X</code> together with annotations of
observations (i.e. cells) <code>obs</code>, variables (i.e. genes) <code>var</code> and unstructured annotations <code>uns</code>.")


## Import files to SnapATAC2

> <hands-on-title> Create an AnnData object </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.6.4+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Import data fragment files and compute basic QC metrics, using 'pp.import_data'`
>        - {% icon param-file %} *"Fragment file, optionally compressed with gzip or zstd"*: `fragments_file.tsv` (Input dataset)
>        - {% icon param-file %} *"A tabular file containing chromosome names and sizes"*: `chrom_sizes.txt` (Input dataset)
>        - {% icon param-toggle %} *"Whether the fragment file has been sorted by cell barcodes"*: `No`
>
>    > <details-title>Sorted by barcodes</details-title>
>    > - This tool requires the fragment file to be sorted according to cell barcodes.
>    > - If **pp.make_fragment_file** {% icon tool %} was used to generate the fragment file, this has automatically been done.
>    >   - Otherwise, the setting *"sorted by cell barcodes"* should remain `No`.
>    {: .details}
>
> 2. Rename the generated file to `Anndata 5k PBMC`
>
> 3. Check that the format is `h5ad`
{: .hands_on}

Because the `AnnData` format is an extension of the HDF5 format, i.e. a binary format, an `AnnData` object can not be inspected directly in Galaxy by clicking on the {% icon galaxy-eye %} (**View data**) icon.
Instead, we need to use a dedicated tool from the **AnnData** suite.

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
>    > <tip-title>Faster Method for General Information</tip-title>
>    > * Many toolsets producing outputs in *AnnData* formats in Galaxy, provide the general information by default:
>    >    * Click on the name of the dataset in the history to expand it.
>    >    * General Anndata information will be given in the expanded box:
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
>    {: .tip}
>
{: .hands_on}


## Calculate and visualize {QC} metrics

> <hands-on-title> Fragment-size distribution </hands-on-title>
>
> 1. {% tool [SnapATAC2 Plotting](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_plotting/snapatac2_plotting/2.6.4+galaxy1) %} with the following parameters:
>    - *"Method used for plotting"*: `Plot fragment size distribution, using 'pl.frag_size_distr'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC` (output of **pp.import_data** {% icon tool %})
> 2. {% icon galaxy-eye %} Inspect the `.png` output
> 3. {% tool [SnapATAC2 Plotting](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_plotting/snapatac2_plotting/2.6.4+galaxy1) %} with the following parameters:
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
> > > 1. 3 peaks are clearly visible (at <100-bp, ~200-bp and ~400-bp). The smallest fragments are from nucleosome-free regions, while the larger peaks of 200- and 400-bp contain mono- and di-nucleosome fragments, respectively.
> > > 2. The small fragments (<100-bp) are from open chromatin reads, since the Tn5 transposase could easily access the loosely packed DNA ({% cite Yan2020 %}).
> > >
> > {: .solution}
> >
> {: .question}
{: .hands_on}


The {TSSe} is another important {QC} metric. Nucleosome-free fragments are expected to be enriched at {TSS}. TSSe shows increased fragmentation of chromatin around the TSS. This suggests open and accessible nucleosome-free chromatin.

{TSSe} is used as a QC metric, since an increased enrichment around TSS regions suggests that the experiment has captured biological meaningful genomic features.
TSSe scores of individual cells can be calculated using SnapATAC2's *metrics.tsse* function.

> <hands-on-title> Calculate and Plot TSSe </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.6.4+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Compute the TSS enrichment score (TSSe) for each cell, using 'metrics.tsse'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC` (output of **pp.import_data** {% icon tool %})
>        - {% icon param-file %} *"GTF/GFF file containing the gene annotation"*: `gene_annotation.gtf.gz` (Input dataset)
>
> 2. Rename the generated file to `Anndata 5k PBMC TSSe` or add the tag {% icon galaxy-tags %} `TSSe` to the dataset:
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
>
> 3. {% tool [SnapATAC2 Plotting](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_plotting/snapatac2_plotting/2.6.4+galaxy1) %} with the following parameters:
>    - *"Method used for plotting"*: `Plot the TSS enrichment vs. number of fragments density figure, using 'pl.tsse'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC TSSe` (output of **metrics.tsse** {% icon tool %})
> 4. {% icon galaxy-eye %} Inspect the `.png` output
>
>
> ![TSSe plot against number of unique fragments]({% link topics/single-cell/images/scatac-standard-snapatac2/pl.tsse.png %})
> High-quality cells can be identified in the plot of {TSSe} scores against a number of unique fragments for each cell.
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
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.6.4+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Filter cell outliers based on counts and numbers of genes expressed, using 'pp.filter_cells'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC TSSe` (output of **metrics.tsse** {% icon tool %})
>        - *"Minimum number of counts required for a cell to pass filtering"*: `5000`
>        - *"Minimum TSS enrichemnt score required for a cell to pass filtering"*: `10.0`
>        - *"Maximum number of counts required for a cell to pass filtering"*: `100000`
>
> 2. Rename the generated file to `Anndata 5k PBMC TSSe filtered` or add the tag {% icon galaxy-tags %} `filtered` to the dataset
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
Currently, our AnnData matrix does not contain any variables. The variables will be added in the following step with the function *pp.add_tile_matrix*.
This creates a cell-by-bin matrix containing insertion counts across genome-wide 500-bp bins.

After creating the variables, the most accessible features are selected.
> <hands-on-title> Select features </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.6.4+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Generate cell by bin count matrix, using 'pp.add_tile_matrix'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC TSSe filtered` (output of **pp.filter_cells** {% icon tool %})
>        - *"The size of consecutive genomic regions used to record the counts"*: `500` 
>        - *"The strategy to compute feature counts"*: `"insertion": based on the number of insertions that overlap with a region of interest` 
>
> 2. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.6.4+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Perform feature selection, using 'pp.select_features'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata tile_matrix` (output of **pp.add_tile_matrix** {% icon tool %})
>        - *"Number of features to keep"*: `250000`
>
>    > <comment-title> Select features </comment-title>
>    >
>    > - Including more features improves resolution and can reveal finer details, but it may also introduce noise.
>    >    - To optimize results, experiment with the `n_features` parameter to find the most appropriate value for your dataset.
>    >
>    >  > <details-title> Different number of features </details-title>
>    >  > - To demonstrate the differences when selecting features, the following UMAP plots are the outputs from processing with a number of features between 1,000 and 500,000.
>    >  > - Fewer features result in fewer, but larger clusters. Selecting a lot of features will output more granular clusters and the compute time will increase.
>    >  > ![Different number of features UMAP]({% link topics/single-cell/images/scatac-standard-snapatac2/number_features.png %}"UMAP plots with different selected features")
>    >  >
>    >  {: .details}
>    > - At this step you can provide a blacklist or whitelist to specifically select relevant features.
>    >    - For example the [**ENCODE Blacklist**](https://github.com/Boyle-Lab/Blacklist) ({% cite Amemiya2019 %}) can be applied here.
>    {: .comment}
>
> 3. Rename the generated file to `Anndata 5k PBMC select_features` or add the tag {% icon galaxy-tags %} `select_features` to the dataset
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

Doublets are removed by calling a customized [**scrublet**](https://github.com/AllonKleinLab/scrublet) ({% cite Wolock2019 %}) algorithm. *pp.scrublet* will identify potential doublets and the function *pp.filter_doublets* removes them.

> <details-title>Doublets in Single-cell datasets</details-title>
>
> - During single-cell sequencing, multiple nuclei can be encapsulated into the same 10x gel bead. The resulting multiplets (>97% of which are doublets) produce sequences from both cells.
> - Doublets can confound the results by appearing as "new" clusters or artifactual intermediary cell states.
>    - These problematic doublets are called **neotypic** doublets, since they appear as "new" cell types.
> - **Scrublet** (Single-cell Remover of Doublets) is an algorithm which can detect neotypic doublets that produce false results.
>    - The algorithm first simulates doublets by combining random pairs of observed cell features.
>    - The observed features of the "cells" are then compared to the simulated doublets and scored on their doublet probability.
>    - SnapATAC2's *pp.filter_doublets* then removes all cells with a doublet probability >50%.
>
> ![Doublet removal with scrublet]({% link topics/single-cell/images/scatac-standard-snapatac2/doublets-and-scrublet.png %} "Scrublet simulates expected doublets and produces doublet scores for each cell. ({% cite Wolock2019 %})")
>
{: .details}

> <hands-on-title> Scrublet </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.6.4+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Compute probability of being a doublet using the scrublet algorithm, using 'pp.scrublet'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC select_features` (output of **pp.select_features** {% icon tool %})
>
> 2. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.6.4+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Remove doublets according to the doublet probability or doublet score, using 'pp.filter_doublets'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata scrublet` (output of **pp.scrublet** {% icon tool %})
>        - *"Threshold for doublet probability"*: `0.5`
> 3. Rename the generated file to `Anndata 5k PBMC filter_doublets` or add the tag {% icon galaxy-tags %} `filter_doublets` to the dataset
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

Dimension reduction is a very important step during the analysis of single cell data. During this, the complex multi-dimensional data is projected into lower-dimensional space, while the lower-dimensional embedding of the complex data retains as much information as possible. Dimension reduction enables batch correction, data visualization and quicker downstream analysis since the data is more simplified and the memory usage is reduced ({% cite Zhang2024%}).


> <details-title>Dimension reduction with SnapATAC2</details-title>
>
> - Dimension reduction algorithms can be either linear or non-linear.
> - Linear methods are generally computationally efficient and well scalable.
>
>   A popular linear dimension reduction algorithm is:
>     - **PCA** (Principle Component Analysis), implemented in **Scanpy** (please check out our [Scanpy]({% link topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.md %}) tutorial for an explanation).
> - Nonlinear methods however are well suited for multimodal and complex datasets.
>     - in contrast to linear methods, which often preserve global structures, non-linear methods have a locality-preserving character.
>     - This makes non-linear methods relatively insensitive to outliers and noise while emphasizing natural clusters in the data ({% cite Belkin2003%})
>     - As such, they are implemented in many algorithms to visualize the data in 2 dimensions (f.ex. **UMAP** embedding).
> - The nonlinear dimension reduction algorithm, through *matrix-free spectral embedding*, used in **SnapATAC2** is a very fast and memory efficient non-linear algorithm ({% cite Zhang2024%}).
>     - **Spectral embedding** utilizes an iterative algorithm to calculate the **spectrum** (*eigenvalues* and *eigenvectors*) of a matrix without computing the matrix itself.
> - For a simple introduction into *spectral embedding* and how it compares to *PCA*, please check out the blog post ["On Laplacian Eigenmaps for Dimensionality Reduction"](https://juanitorduz.github.io/laplacian_eigenmaps_dim_red/) by Juan Orduz.
>
{: .details}

## Spectral embedding
The dimension reduction, produced by the algorithm *tl.spectral*, is required for later steps, such as plotting and clustering.

> <hands-on-title> Spectral embedding </hands-on-title>
>
> 1. {% tool [SnapATAC2 Clustering](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_clustering/snapatac2_clustering/2.6.4+galaxy1) %} with the following parameters:
>    - *"Dimension reduction and Clustering"*: `Perform dimension reduction using Laplacian Eigenmap, using 'tl.spectral'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC filter_doublets` (output of **pp.filter_doublets** {% icon tool %})
>        - *"Distance metric"*: `cosine`
>
>    > <comment-title> Distance metric </comment-title>
>    >
>    > - The fast and well scalable *"matrix-free spectral embedding"* algorithm depends on the distance metric: `cosine`
>    {: .comment}
>
> 2. Rename the generated file to `Anndata 5k PBMC spectral` or add the tag {% icon galaxy-tags %} `spectral` to the dataset
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
With the already reduced dimensionality of the data stored in `X_spectral`, the cells can be further embedded (i.e. transformed into lower dimensions) with {UMAP}.
**UMAP** projects the cells and their relationship to each other into 2-dimensional space, which can be easily visualized ({% cite McInnes2018%}).

> <hands-on-title> UMAP embedding </hands-on-title>
>
> 1. {% tool [SnapATAC2 Clustering](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_clustering/snapatac2_clustering/2.6.4+galaxy1) %} with the following parameters:
>    - *"Dimension reduction and Clustering"*: `Compute Umap, using 'tl.umap'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC spectral` (output of **tl.spectral** {% icon tool %})
>        - *"Use the indicated representation in '.obsm'"*: `X_spectral`
>
> 2. Rename the generated file to `Anndata 5k PBMC UMAP` or add the tag  {% icon galaxy-tags %} `UMAP` to the dataset
{: .hands_on}

# Clustering
During clustering, cells that share similar accessibility profiles are organized into clusters. **SnapATAC2** utilizes graph-based community clustering with the *Leiden* algorithm ({% cite Traag2019%}).
This method takes the k-nearest neighbor (KNN) graph as input data and produces well-connected communities.


## Community clustering

> <hands-on-title> Clustering analysis </hands-on-title>
>
> 1. {% tool [SnapATAC2 Clustering](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_clustering/snapatac2_clustering/2.6.4+galaxy1) %} with the following parameters:
>    - *"Dimension reduction and Clustering"*: `Compute a neighborhood graph of observations, using 'pp.knn'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC UMAP` (output of **tl.umap** {% icon tool %})
>
> 2. {% tool [SnapATAC2 Clustering](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_clustering/snapatac2_clustering/2.6.4+galaxy1) %} with the following parameters:
>    - *"Dimension reduction and Clustering"*: `Cluster cells into subgroups, using 'tl.leiden'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata knn` (output of **pp.knn** {% icon tool %})
>        - *"Whether to use the Constant Potts Model (CPM) or modularity"*: `modularity`
>
>    > <comment-title> CPM or modularity </comment-title>
>    > - make sure you selected `modularity`
>    > - the clusters produced by `CPM` are not represented well in the UMAP projections
>    >
>    {: .comment}
> 2. Rename the generated file to `Anndata 5k PBMC leiden` or add the tag {% icon galaxy-tags %} `leiden` to the dataset
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
Now that we have produced **UMAP** embeddings of our cells and have organized the cells into **leiden** clusters, we can now visualize this information with *pl.umap*.

> <hands-on-title> Plotting the clusters </hands-on-title>
>
> 1. {% tool [SnapATAC2 Plotting](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_plotting/snapatac2_plotting/2.6.4+galaxy1) %} with the following parameters:
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
> > > 1. There are 13 leiden clusters.
> > > 2. Clusters in close proximity (f.ex. clusters 0 and 5) share a similar chromatin accessibility profile (and very likely also a similar cell type), compared to a cluster further away (f.ex. cluster 9).
> > >
> > {: .solution}
> >
> {: .question}
{: .hands_on}

# Cell cluster annotation

After clustering the cells, they must be annotated. This categorizes the clusters into known cell types. **Manual Cell Annotation** requires known marker genes and varying expression profiles of the marker genes among clusters.

Luckily, the marker genes for {PBMCs} are known ({% cite Sun2019 %}).
Marker genes for other single cell datasets can also be found in databases such as [PanglaoDB](https://panglaodb.se/markers.html) ({% cite Franzn2019 %}). Using the known marker genes we can now annotate our clusters manually.

## Gene matrix
Since our data currently doesn't contain gene information, we have to create a cell-by-gene activity matrix using the function *pp.make_gene_matrix*.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SnapATAC2 Preprocessing](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_preprocessing/snapatac2_preprocessing/2.6.4+galaxy1) %} with the following parameters:
>    - *"Method used for preprocessing"*: `Generate cell by gene activity matrix, using 'pp.make_gene_matrix'`
>        - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC leiden` (output of **tl.leiden** {% icon tool %})
>        - {% icon param-file %} *"GTF/GFF file containing the gene annotation"*: `gene_annotation.gtf.gz` (Input dataset)
> 2. Rename the generated file to `Anndata 5k PBMC gene_matrix` or add the tag {% icon galaxy-tags %} `gene_matrix` to the dataset
>    > <tip-title> Gene matrix </tip-title>
>    >
>    > - Please note that *pp.make_gene_matrix* removes all annotations except those stored in `obs`.
>    > - Therefore it might be necessary to remove propagating tags {% icon galaxy-tags %} (tags starting with `#`) from `Anndata 5k PBMC gene_matrix`.
>    >    - Tags can be removed by expanding the dataset with a tag and clicking {% icon galaxy-cross %} next to the tag.
>    {: .tip}
>
> 3. {% icon galaxy-eye %} Inspect the general information of the `.h5ad` output
>
>    > <question-title></question-title>
>    >
>    > ```
>    > AnnData object with n_obs × n_vars = 4430 × 60606
>    >  obs: 'n_fragment', 'frac_dup', 'frac_mito', 'tsse', 'doublet_probability', 'doublet_score', 'leiden'
>    > ```
>    >
>    > What does `n_vars` represent in `Anndata 5k PBMC gene_matrix` and what did it represent in `Anndata 5k PBMC leiden`?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > - The variables now represent accessible genes. There are 60606 accessible genes in our samples. In `Anndata 5k PBMC leiden` and all earlier AnnData the variables represented the 6062095 fixed-sized genomic bins.
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
{: .hands_on}


## Imputation with Scanpy and MAGIC
Similar to scRNA-seq data, the cell-by-gene-activity matrix is very sparse. Additionally, high gene variance between cells, due to technical confounders, could impact the downstream analysis.
In scRNA-seq, filtering and normalization are therefore required to produce a high-quality gene matrix.

Since the *cell-by-gene-activity* matrix resembles the *cell-by-gene-expression* matrix of scRNA-seq, we can use the tools of the [Scanpy](https://scanpy.readthedocs.io/en/stable/index.html) ({%cite Wolf2018%}) tool suite to continue with our data.

> <details-title>Imputation with MAGIC</details-title>
>
> - The count matrices of single-cell data are sparse and noisy.
> - Confounding issues, such as "dropout" effects, where some mRNA or DNA-segments are not detected although they are present in the cell, also result in some cells missing important cell-type defining features.
>    - These problems can obscure the data, as only the strongest gene-gene relationships are still detectable.
> - The *Markov Affinity-based Graph Imputation of Cells* (MAGIC) algorithm ({%cite vanDijk2018%}) tries to solve these issues by filling in missing data from some cells with transcript information from similar cells.
>    - The algorithm calculates the likely gene expression of a single cell, based on similar cells and fills in the missing data to produce the expected expression.
>      - *MAGIC* achieves this by building a graph from the data and using data diffusion to smooth out the noise.
>
> ![Imputation with the MAGIC algorithm]({% link topics/single-cell/images/scatac-standard-snapatac2/magic_method.png %} "MAGIC restores noisy and sparse single-cell data using diffusion geometry ({%cite vanDijk2018%})")
>
{: .details}

> <hands-on-title> Filter and normalize </hands-on-title>
>
> 1. {% tool [Filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.9.6+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC gene_matrix` (output of **pp.make_gene_matrix** {% icon tool %})
>    - *"Method used for filtering"*: `Filter genes based on number of cells or counts, using 'pp.filter_genes'`
>        - {% icon param-select %} *"Filter"*: `Minimum number of cells expressed`
>            - *"Minimum number of cells expressed required for a gene to pass filtering"*: `5`
>
> 2. {% tool [Normalize](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_normalize/scanpy_normalize/1.9.6+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Anndata filter_genes` (output of **Filter** {% icon tool %})
>    - *"Method used for normalization"*: `Normalize counts per cell, using 'pp.normalize_total'`
>    - {% icon param-toggle %} *"Exclude (very) highly expressed genes for the computation of the normalization factor (size factor) for each cell"*: `No`
>
> 3. {% tool [Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.9.6+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Anndata normalize` (output of **Normalize** {% icon tool %})
>    - *"Method used for inspecting"*: `Logarithmize the data matrix, using 'pp.log1p'`
>
> 4. {% tool [Normalize](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_normalize/scanpy_normalize/1.9.6+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Anndata log1p` (output of **log1p** {% icon tool %})
>    - *"Method used for normalization"*: `Denoising using Markov Affinity-based Graph Imputation of Cells (MAGIC) API 'external.pp.magic'`
>        - *"Denoised genes to return"*: `All genes`
>        - *"Which solver to use"*: `"approximate", is faster that performs imputation in the PCA space and then projects back to the gene space`
>
>    > <comment-title> </comment-title>
>     >   - Choosing the setting `Which solver to use: 'exact'` will result in a output file with better resolution.
>     >   - This is not necessary for our purposes, since the compute time also increases with this setting.
>     >
>     {: .comment}
>
> 5. Rename the generated file to `Anndata 5k PBMC magic` or add the tag {% icon galaxy-tags %} `magic` to the dataset
>
> 6. {% icon galaxy-eye %} Inspect the general information of the `.h5ad` output
>
>
>    > <question-title></question-title>
>    >
>    > ```
>    > AnnData object with n_obs × n_vars = 4430 × 55106
>    >  obs: 'n_fragment', 'frac_dup', 'frac_mito', 'tsse', 'doublet_probability', 'doublet_score', 'leiden', 'n_genes', 'n_counts'
>    >  var: 'n_cells', 'n_counts'
>    >  uns: 'log1p'
>    > ```
>    >
>    > 1. How did `n_vars` change, compared to `Anndata 5k PBMC gene_matrix`?
>    > 2. Which data was 'lost', compared to `Anndata 5k PBMC leiden`, and must be added to the file in order to produce **UMAP** plots?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > 1. The number of genes was reduced from 60606 to 55106 by the filtering. Additional annotations were added, such as: `obs: 'n_genes', 'n_counts'`, `var: 'n_cells', 'n_counts'` and `uns: 'log1p'`.
>    > > 2. The UMAP embeddings `obsm: 'X_umap'` are missing and should be added to the Anndata in the next step. Without `X_umap` it won't be possible to visualize the plots.
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
>    - {% icon param-file %} *"Input object in hdf5 AnnData format"*: `Anndata 5k PBMC magic` (output of **external.pp.magic** {% icon tool %})
>    - *"Copy embeddings (such as UMAP, tSNE)"*: `Yes`
>       - *"Keys from embeddings to copy"*: `X_umap`
>       - {% icon param-file %} *"IAnnData objects with embeddings to copy"*: `Anndata 5k PBMC leiden`
>
>    > <comment-title> Annotations to copy </comment-title>
>    >
>    > - This tutorial only focuses on producing an **UMAP** plot with marker-genes.
>    > - If further analysis, with tools requiring more annotations, is intended, these can be added in a similar way as shown above.
>    >     - f.ex. *Peak and Motif Analysis* with {% tool [Snapatac2 peaks and motif](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_peaks_and_motif/snapatac2_peaks_and_motif/2.6.4+galaxy1) %} requires annotations from `uns`.
>    > - It is also possible to leave the input *"Keys from embeddings to copy"* empty, to copy all annotations of a given category such as `obsm`.
>    {: .comment}
>
> 5. Rename the generated file to `Anndata 5k PBMC magic UMAP` or add the tag {% icon galaxy-tags %} `UMAP` to the dataset
>
> 6. {% icon galaxy-eye %} Inspect the general information of the `.h5ad` output, to check if `obsm` contains `X_umap`
>
>   ```
>   AnnData object with n_obs × n_vars = 4430 × 55106
>    obs: obs: 'n_fragment', 'frac_dup', 'frac_mito', 'tsse', 'doublet_probability', 'doublet_score', 'leiden', 'n_genes', 'n_counts'
>    var: 'n_cells', 'n_counts'
>    uns: 'log1p'
>    obsm: 'X_umap'
>   ```
>
{: .hands_on}

## Visualize gene activity of marker genes
The gene activity of selected marker genes can now be visualized with Scanpy.

> <hands-on-title> Plot marker genes </hands-on-title>
>
> 1. {% tool [Plot with Scanpy](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.9.6+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output_h5ad` (output of **AnnData Operations** {% icon tool %})
>    - *"Method used for plotting"*: `Embeddings: Scatter plot in UMAP basis, using 'pl.umap'`
>        - *"Keys for annotations of observations/cells or variables/genes"*: `leiden, CD3D, CD8A, CD4, MS4A1, NKG7, CD14, FCER1A`
>        - {% icon param-toggle %} *"Show edges?"*: `No`
>        - In *"Plot attributes"*
>           - *"Number of panels per row"*: `2`
>
>  ![umap_leiden_marker_gene_clustering]({% link topics/single-cell/images/scatac-standard-snapatac2/umap_leiden_marker-genes.png %})
>
> > <question-title></question-title>
> >
> > 1. Are the marker genes selectively expressed in the clusters?
> > 2. Which marker genes have overlapping expression profiles? And what could that imply?
> >
> > > <solution-title></solution-title>
> > >
> > > 1. Some marker genes, such as `MS4A1` or `CD8A`, are only expressed in a few clusters (clusters 6+11 and clusters 1+8, respectively).
> > > 2. The marker gene `CD3D` is expressed in multiple clusters (1, 2, 4, 7 and 8). Overlapping expression profiles imply similar cell types since similar cell types have similar marker genes upregulated. In this case, `CD3D` expression classifies the cells in these clusters as T-cells.
> > >
> > {: .solution}
> >
> {: .question}
{: .hands_on}

## Manual cluster annotation
Comparison of marker gene expression in our clusters with a table of canonical marker genes ({% cite Sun2019 %}), enables us to annotate the clusters manually.

Cell type | Marker genes
--- | ---
CD8+ T cells | CD3D+, CD8A+, CD4-
CD4+ T cells | CD3D+, CD8A-, CD4+
B cells | CD3D-, MS4A1+
Natural killer (NK) cells | CD3D-, NKG7+
Monocytes | CD3D-, CD14+
Dendritic cells | CD3D-, FCER1A+

These canonical marker genes can match the clusters to known cell types:

Cluster | Cell type
--- | ---
0 | Monocytes
1 | CD8+ T cells
2 | CD4+ T cells
3 | Monocytes
4 | CD4+ T cells
5 | Monocytes
6 | B cells
7 | CD4+ T cells
8 | CD8+ T cells
9 | NK cells
10 | Monocytes
11 | B cells
12 | Dendritic cells

> <comment-title></comment-title>
> Note that some clusters contain subtypes (f.ex. the annotated B cell clusters contain both naive and memory B cells). The cell-type annotation can be refined by choosing more specific marker genes.
{: .comment}

To manually annotate the *Leiden* clusters, we will need to perform multiple steps:

1. **Inspect** the key-indexed observations of `Anndata 5k PBMC gene_matrix magic UMAP`
2. **Cut** the *Leiden* annotations out of the table
3. **Upload** a *replace file* containing the new cell type annotations for the *Leiden* clusters
4. **Replace** the values of the cluster annotation with cell type annotation
5. **Add** the cell type annotation to the AnnData
6. **Plot** the annotated cell types

> <hands-on-title> Manual annotation </hands-on-title>
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.10.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC gene_matrix magic UMAP`
>    - *"What to inspect?"*: `Key-indexed observations annotation`
> 2. Rename the generated file to `5k PBMC observations`
> 3. {% icon galaxy-eye %} Inspect the generated file
>
>    > <question-title></question-title>
>    > In which column is the **Leiden** annotation located?
>    > > <solution-title></solution-title>
>    > > The **Leiden** annotation is in column 8.
>    > >
>    > > Column 1 | Column 2 | Column 3 | Column 4 | Column 5 | Column 6 | Column 7 | Column 8 | Column 9 | Column 10
>    > > --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
>    > > "" | n_fragment | frac_dup | frac_mito | tsse | doublet_probability | doublet_score | leiden | n_genes | n_counts
>    > > AAACGAAAGACGTCAG-1 | 22070 | 0.5219425551271499 | 0.0 | 30.43315066436454 | 0.004634433324822066 | 0.009276437847866418 | 8 | 52303 | 16521.599844068267
>    > > AAACGAAAGATTGACA-1 | 10500 | 0.5345125681606597 | 0.0 | 29.10551296093465 | 0.004668403569267374 | 0.001088139281828074 | 1 | 54501 | 15020.42495602328
>    > > AAACGAAAGGGTCCCT-1 | 19201 | 0.5101785714285714 | 0.0 | 19.90011850347046 | 0.004634433324822066 | 0.009276437847866418 | 5 | 54212 | 16294.751533305309
>    > > AAACGAACAATTGTGC-1 | 13242 | 0.487399837417257 | 0.0 | 29.060913705583758 | 0.004660125753854076 | 0.0022172949002217295 | 7 | 53530 | 15456.629863655084
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
>
> 4. {% tool [Cut columns](toolshed.g2.bx.psu.edu/repos/devteam/cut_columns/Cut1/1.0.2) %} with the following parameters:
>    - {% icon param-select %} *"Cut columns"*: `c8`
>    - {% icon param-file %} *"From"*: `5k PBMC observations` (output of **Inspect AnnData** {% icon tool %})
>
> 5. Create a new **tabular** file from the following
>    ```
>    leiden cell_type
>    0 Monocytes
>    1 CD8_Tcells
>    2 CD4_Tcells
>    3 Monocytes
>    4 CD4_Tcells
>    5 Monocytes
>    6 Bcells
>    7 CD4_Tcells
>    8 CD8_Tcells
>    9 NKcells
>    10 Monocytes
>    11 Bcells
>    12 Dendritic_cells
>    ```
>    {% snippet faqs/galaxy/datasets_create_new_file.md name='replace_file' format='tabular' convertspaces='True' %}
>
>    > <details-title>Replace file</details-title>
>    >
>    > - The first column of the replace file contains the "old" annotations and the second column contains the "new" annotation.
>    > - {% icon warning %} Spaces in the new annotations can lead to errors. Please use underscores (`_`) instead.
>    >
>    {: .details}
> 6. {% tool [Replace column](toolshed.g2.bx.psu.edu/repos/bgruening/replace_column_by_key_value_file/replace_column_with_key_value_file/0.2) %} with the following parameters:
>    - {% icon param-file %} *"File in which you want to replace some values"*: `Cut columns leiden` (output of **Cut columns** {% icon tool %})
>    - {% icon param-file %} *"Replace information file"*: `replace_file`
>    - *"Which column should be replaced?*: `Column: 1`
> 7. Rename the generated file to `Cell type annotation` and {% icon galaxy-eye %} inspect the file to check if the replacement was successful
>    - check if the datatype is set to `tabular`
>    - you may need to change the datatype manually
>
>      {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular"%}
>
> 8. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.10.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC gene_matrix magic UMAP`
>    - *"Function to manipulate the object"*: `Add new annotation(s) for observations or variables`
>        - *"What to annotate"*: `Observations (obs)`
>        - {% icon param-file %} *"Table with new annotations"*: `Cell type annotation`(output of **Replace column** {% icon tool %})
> 9. Rename the generated file to `Anndata 5k PBMC gene_matrix magic cell_type` or add the tag {% icon galaxy-tags %} `cell_type` to the dataset
>
> 10. {% tool [Plot with Scanpy](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.9.6+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Anndata 5k PBMC gene_matrix magic cell_type` (output of **Manipulate AnnData** {% icon tool %})
>    - *"Method used for plotting"*: `Embeddings: Scatter plot in UMAP basis, using 'pl.umap'`
>        - *"Keys for annotations of observations/cells or variables/genes"*: `cell_type`
>        - {% icon param-toggle %} *"Show edges?"*: `No`
>        - In *"Plot attributes"*
>           - *"Location of legend"*: `on data`
>           - {% icon param-toggle %} *"Draw a frame around the scatter plot?"*: `Yes`
> 11. {% icon galaxy-eye %} Inspect the `.png` output
>
>  ![UMAP annotated cell types]({% link topics/single-cell/images/scatac-standard-snapatac2/umap_cell-types.png %})
>
>    > <question-title></question-title>
>    > 1. Are clusters with the same assigned cell type located close to each other?
>    > 2. Are the [expected cell type percentages of {PBMCs}](https://www.akadeum.com/blog/what-are-pbmc-human-pbmc-cells/) visible in the annotated plot?
>    >
>    >    Cell Type | Expected Percentage
>    >    --- | ---
>    >    CD3+ T cells | 50-70%
>    >    CD4+ T cells | 25-45%
>    >    CD8+ T cells | 10-25%
>    >    B cells | 5-15%
>    >    NK cells | <5%
>    >    Monocytes | 10-30%
>    >    Dendritic cells | <2%
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > 1. Yes. B-cells are far away from NK and T cells. Only the myeloid lineage of monocytes and dendritic cells are located close to each other. This might be due to the common progenitor cell lineage and thus a similar chromatin profile.
>    > > 2. Yes. The most common cell types are T cells, followed by Monocytes and B cells. NK cells and Dendritic cells only make up a small percentage of PBMCs.
>    > {: .solution}
>    >
>    {: .question}
>
{: .hands_on}


# Conclusion
{% icon congratulations %} Well done, you’ve made it to the end! You might want to consult your results with this [control history](https://usegalaxy.eu/u/videmp/h/standard-processing-of-10x-single-cell-atac-seq-data-with-snapatac2), or check out the [full workflow](https://usegalaxy.eu/u/videmp/w/standard-processing-of-10x-single-cell-atac-seq-data-with-snapatac2) for this tutorial.

In this tutorial, we produced a count matrix of {scATAC-seq} reads in the `AnnData` format and performed:
1. Preprocessing:
   1. Plotting the fragment-size distributions
   2. Calculating and plotting {TSSe} scores
   3. Filtering cells and selecting features (fixed-size genomic bins)
2. Dimension reduction through **Spectral embedding** and **{UMAP}**
3. Clustering of cells via the *Leiden* method
4. Cluster annotation
   1. Producing and filtering a cell by gene activity matrix
   2. Data normalization and imputation with **Scanpy** and the *MAGIC* algorithm
   3. Visualizing marker genes in the clusters
   4. Manually annotating the cell types with selected marker genes

![SnapATAC2 processing pipeline]({% link topics/single-cell/images/scatac-standard-snapatac2/snapatac2-pipeline.png %})

The {scATAC-seq} analysis can now continue with downstream analysis, for example *differential peak analysis*. 

> <details-title>Differential peak analysis</details-title>
>
> The **SnapATAC2** tools for differential peak analysis are already accessible on Galaxy. However, there are no GTN trainings available yet. Until such a tutorial is uploaded, you can visit the **SnapATAC2** documentation for a [tutorial on differential peak analysis](https://kzhang.org/SnapATAC2/version/2.6/tutorials/diff.html). And check out our [example history](https://usegalaxy.eu/u/timonschlegel/h/differential-peak-analysis-with-snapatac2) and the exemplary [workflow](https://usegalaxy.eu/u/timonschlegel/w/copy-of-differential-peak-analysissnapatac2). 
>
> The tools are available in Galaxy under {% tool [SnapATAC2 Peaks and Motif](toolshed.g2.bx.psu.edu/repos/iuc/snapatac2_peaks_and_motif/snapatac2_peaks_and_motif/2.6.4+galaxy1) %}. 
> - If you want to continue with differential peak analysis, please make sure that the AnnData object with the annotated cell types contains unspecified annotations for the reference sequences (`uns: 'reference_sequences'`). 
>   - The section [Copy-over embeddings]( {% link topics/single-cell/tutorials/scatac-standard-processing-snapatac2/tutorial.md %}#copy-over-embeddings ) explains how to copy annotations from one AnnData object to another. 
>
{: .details}