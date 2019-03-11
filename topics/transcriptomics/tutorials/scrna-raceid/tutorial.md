---
layout: tutorial_hands_on
title: Downstream single-cell RNA analysis with RaceID
zenodo_link: 'https://zenodo.org/record/1511582'
tags:
  - single-cell
questions:
  - What are
  - What is the difference between PCA and tSNE?
objectives:
  - Clustering cells in a matrix
  - Assessing the quality of individual clusters
  - Inferring a lineage between cell types
  - Examining gene expression
  - Determining the top most expressive genes per cluster
  - Correcting for unwanted variation
#requirements:
#  -
#    type: "internal"
#    topic_name: transcriptomics
#    tutorials:
#        - scrna-introduction
#        - scrna-plates-batches-barcodes
#        - scrna-umis
#        - scrna_preprocessing

time_estimation: 2H
key_points:
  - The take-home messages
  - They will appear at the end of the tutorial
contributors:
  - mtekman

---


# Introduction
{:.no_toc}

The data provided here as part of this tutorial analyses single-cell RNA-seq data from a study published by [Grün et.al](https://doi.org/10.1016/j.stem.2016.05.010) in 2016. The data was used to cluster cells from *Lgr5*-positive intestinal stem cells of C57BL6/J mice, with the aim of discovering distinct cell sub-populations and deriving a lineage tree between them to find out how these sub-populations relate (or are derived from) one another.

The input data consists of a single count matrix consisting of ~21,000 genes (rows) and ~400 cells (columns), following the [tidy data](https://cran.r-project.org/web/packages/tidyr/vignettes/tidy-data.html) convention prevalent amongst the R data analysis community which assigns every value to a variable and an observation.

Here, the values are the number of reads which are assigned to a particular gene (a variable) that was observed within a specific cell (an observation).

Normally a count matrix consists of integers, but this matrix has undergone an UMI-to-transcript [count alteration](https://www.nature.com/articles/nmeth.2930#methods) to correct against UMI errors, yielding decimal values instead.

This tutorial will perform cell clustering and lineage construction, as well as exploring some genes of interest.


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Getting the Count Matrix

Every single-cell RNA analysis begins with a count matrix, which contains the raw data required for the downstream analysis. Annotation data such as cell phenotype or gene annotations are sometimes also used to enrich the analysis, but these are typically only useful later, and can be generated from external databases if required.

For now we will work with the count matrix alone.

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial and name it "RaceID on scRNA"
> 2. Import the file from [Zenodo](https://zenodo.org/record/1511582) or from the shared data library
>
>    ```
>    https://zenodo.org/record/1511582/files/intestinalData.tsv
>    ```
> 
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Rename the dataset to *"intestinal"*
> 4. Check that the datatype is a tab-seperated file
>
>    {% include snippets/change_datatype.md datatype="tabular" %}
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. How many genes are in the count matrix?
> 2. How many cells?
>
> > ### {% icon solution %} Solution
> >
> > A summary of our dataset is given by expanding the file preview window by clicking on the file name.
> > 
> > 1. 20,269 lines, where lines/rows denote our genes.
> > 2. Scroll the mini-preview window to the right to see that the number of cells are 431
> >
> {: .solution}
>
{: .question}



## Inspecting the Cell Labelling

We can inspect the dataset by clicking on the {% icon galaxy-eye %}  symbol. Immediately we can see that the header of this file contains cell names, following a naming convention that separates the cells into 5 phenotypes: *I5d, II5d, III5d, IV5d*, and *V5d*

We can see this for ourselves by extracting the headers, and reformatting them to see how many unique types we can detect:


> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Select first** {% icon tool %} with the following parameters:
>    - *"Number of lines"*: `1`
>
> 1. **Transpose** {% icon tool %} with the following parameters:
>
> 1. **Text transformation** {% icon tool %} with the following parameters:
>    - *"SED Program"*: `s/_[0-9]+//`
>
>    > ### {% icon comment %} Comment
>    >
>    > The above text is a regular expression used to match on anything that contains a `_` followed by a number, and removing it.
>    {: .comment}
>
> 1. **Sort** {% icon tool %} with the following parameters:
>    - In *"Column selections"*:
>        - {% icon param-repeat %} *"Insert Column selections"*
>            - *"on column"*: `c1`
>
> 1. **Unique lines** {% icon tool %} with the following parameters:
>    - *"Do you want to group each unique group"*: `No`
>       - *"Counting number of occurrences"*: `Yes`
>       - *"Only print duplicate lines"*: `No`
>       - *"Only print unique lines"*: `No`
>    - *"Ignore differences in case when comparing"*: `No`
>
{: .hands_on}


> ### {% icon question %} Questions
>
> 1. What did the transpose step do? Why is it neccesary?
> 2. How many rows remain after the **Text transformation** step?
> 3. How many unique cell phenotypes were identified in the cell headers?
> 4. Which cell phenotype is least represented in the count matrix?
>
> > ### {% icon solution %} Solution
> > 
> > 1. The sole purpose of the Transpose tool is to switch columns with rows (and vice versa), which will make it easier to inspect and sort data.
> > 2. The number of rows has not changed since the last step, but the cell names have lost their numbering and are identified purely by their phenotype.
> > 3. There are 5 types of cells in our count matrix: *I5d*,*II5d*,*III5d*,*IV5d*, and *V5d*.
> > 4. There are only 48 *IV5d* cells compared to the other types which have 95 or 96.
> {: .solution}
>
{: .question}

With these types already labelled in the header of our data, we can validate the clustering that we will perform later. Ideally, the cells described by these 5 different labels should cluster into 5 seperate clusters, with varying degrees of proximity to one another.

## Inspecting the Quality of the Count Matrix

Low quality cells and genes are sometimes caught at the preprocessing stage and removed, but sometimes more filtering is required to remove unwanted noise from the data.

A gene that has a low number of total counts across multiple cells might be differentially expressed across those cells (e.g. 1 count in CellA, and 10 counts in CellB, yielding a fold change of 10) but would not neccasarily be significant compared to a gene that has more total counts across cells (e.g. 100 counts in cellA, and 1000 counts in CellB, also yielding a fold change 10).

We can refine filtering thresholds by examining how much a histogram of our plots change before and after filtering using standard parameters

> ### {% icon hands_on %} Hands-on: Unique Cell Types
>
> 1. **Filtering, Normalisation, and Confounder Removal using RaceID** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Count Matrix"*: `intestinal` (Input dataset)
>    - In *"Filtering"*:
>        - *"Min Transcripts"*: `3000`
>        - *"Min Expression"*: `5`
>        - *"Min Cells"*: `5`
>        - *"Use Defaults?"*: `Yes`
>
>
{: .hands_on}

This tool generates four histograms with the top line giving the raw expression data fed into the tool, and the bottom line giving the filtered data. 

![Histograms of Raw and Filtered Data]({{site.baseurl}}{% link topics/transcriptomics/images/raceid_filter_plots.png %} "Histograms of Raw and Filtered Data")

<!-- Todo: update image with vertical red bar -->
<!-- Todo: add feature that looks for number of features >=1, instead of >0, because RaceID adds 0.1 to 0 counts -->

The top row shows the count distributions of the Library Size and Number of Features of the raw data:

* (Top-Left) Library Size (total number of transcripts per cell)
  * A lower-tail heavy distribution centred around $$10^{3-4}$$ counts per cell, with a few cells having library sizes containing a handful of counts (0-1000).
* (Top-Right)
  * Another lower-tail heavy distribution with a peak centred around $$10^3.5$$. The lower of the two profiles is defined somewhat sporadically, so we cannot assume that cells with less than 10 features (< 1.0) indicate anything meaningful except just being of low quality.

The bottom row shows the count distributions of the Library Size and Number of Features of the filtered data:

* (Bottom-Left) The lower-tail of our previous distribution has been trimmed off (note that we required cells with no less than 3000 total transcripts), which gives an even normal-looking distribution centred around $$10^3.4$$ transcripts per cell.
* (Bottom-Right) Instead of a distribution we have a single bar that indicates that all of our cells have the exact number of features.

> ### {% icon details %} Details: Why the Same Number of Features
> * RaceID normalises the data so that all cells are compared using the same features. If the features compared between cells are different, then it is hard to make a meaningful assessment of how much one cell differs from another.
> * For cells that are not filtered out during this stage which have less than the "required" number of features, a value of 0.1 is added to the count data so that these features are not lost during the analysis. This makes the assumption that the feature *is* detectable for that cell (i.e. no errors during sequencing) but that the transcript was very lowly expressed.
> * For a more 'realistic' distribution of features, re-run the tool with *"Count filtered features greater than or equal to 1"* enabled.
{: .details}

> ### {% icon question %} Questions
>
> 1. How many cells remain after filtering?
> 2. How many genes remain after filtering?
>
> > ### {% icon solution %} Solution
> >
> > The answer to both questions can be seen in Metrics file
> >
> > 1. 287 cells remain (66%)
> > 2. 2089 genes remain (10%)
> >
> {: .solution}
>
{: .question}


The filtered distributions are what we expect a well filtered and normalised dataset to resemble: a count matrix with all observations having roughly the same number of transcripts, but distributed differently across different sets of common features. With this we can now perform initial clustering to see whether we can cluster any of our cells into distinct cell types.



# Initial Clustering

There are two approaches to the analysis:

* *Naive Approach*, which assumes that the data contains only the biologically relevant information that we wish to cluster, and that rare cell subtypes have not been filtered out.

* *Refined Approach*, which takes into account that there may be some technical noise and unwanted sources of technological or biological variability in the data that needs to be accounted for.

### Biological Variation

![Sources of variation]({{site.baseurl}}{% link topics/transcriptomics/images/raceid_cellcycle.svg %} "Sources of unwanted biological variation: (Left) Transcriptional Bursting, and (Right) Cell-cycle Variation")

Transcriptional bursting is a stochastic model for the transcription process in a cell, where transcription does not occur as a smooth or continuous process but occurs in spontaneous and discrete 'bursts' thought to be only loosely assosciated with chromatin conformation/availability. It is an effect that is not seen in bulk RNA-seq due to the smoothing effect of measuring average gene expression across a tissue. However, the effect is more pronounced in single-cell and it is hard to model against.

On the other hand, cell-cycle variation is well defined and can be modelled against. As the cell grows from the G1 to the M phase, the amount of mRNA transcribed grows with it, meaning that cells in the later stages of their cycle are more likely to produce more transcripts of a given gene than a cell of the same type in the earlier stages of its cycle. Such differences can give false variation that would cluster two cells of the same type but at different time-points seperately. Fortunately, there are a well-defined set of genes whose expression is known to covary with the cell-cycle, and thus this effect can be modelled out.


### Technical Variation

![Sources of variation]({{site.baseurl}}{% link topics/transcriptomics/images/raceid_technical_variation.svg %} "Sources of unwanted technical variation")

Technical variation appears in three main forms: *Library size variation*, *Amplification bias*, and *Dropout events*.

1. **Library size variation** occurs where two cells of the same cell type may produce a different amount of total transcripts than one another (e.g. due to cell-cycle effects, or differing capture efficiencies), but have a similar *proportion* of transcripts for specific genes. For example, a neural cell with a library size of 100 may express 10 counts of SOX2 (10%), and another neural cell with a library size of 200 may express 20 counts of SOX2 (also 10%). The two cells have different library sizes, but harbor the same expression for that gene because they are of the same cell type.

1. **Amplification bias** stems from an uneven amplification of certain transcripts of a cell over others, giving a false number of reads for the number of mRNA molecules actually observed in the cell. Unique Moleculer Identifiers can significantly reduce this bias, and it is a topic that is covered more extensively in the [*Understanding Barcodes*]({{site.baseurl}}{% link topics/transcriptomics/tutorials/scrna-umis/tutorial.md %}) hands-on.

1. **Dropout events** are the zero counts that are prevalent in the data due to the reduced sequencing sensitivity in detecting reads, which yields many false negatives in the detection of genes, often resulting in over 80% of the count values in the count matrix being zero. A major point to take into account is that some of these zeroes are *real* (i.e. no transcripts of that gene were detected in that cell) and some of these are *false* (i.e. the transcripts were never captured due to the low sequencing depth). Modelling this duality in the data and mitigating against it is one of the biggest challenges of normalising single-cell data.


## Naive Approach

The naive clustering approach assumes that there is no unwanted technical or biological variability in the data and that the cells will cluster purely based on their phenotypes. This assumption is not completely without merit, since often the biological signal is strong enough to counter the lesser unwanted variation.

Here we will attempt to perform some filtering, normalisation, and clustering using the recommended default settings to see if we can detect different cell types. 


### Filtering, Normalisation, and Clustering

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Filtering, Normalisation, and Confounder Removal using RaceID** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Count Matrix"*: `output` (Input dataset)
>    - In *"Filtering"*:
>        - *"Min Transcripts"*: `3000`
>        - *"Min Expression"*: `5`
>        - *"Min Cells"*: `5`
>        - *"Use Defaults?"*: `Yes`
>
> 1. **Clustering using RaceID** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input RaceID RDS"*: `outrdat` (output of **Filtering, Normalisation, and Confounder Removal using RaceID** {% icon tool %})
>    - In *"Clustering"*:
>        - *"Use Defaults?"*: `Yes`
>    - In *"Outliers"*:
>        - *"Use Defaults?"*: `Yes`
>    - In *"tSNE and FR"*:
>        - *"Use Defaults?"*: `Yes`
> 
{: .hands_on}

The first three plots tell us about the stability/reliability of our clusters, and are more important indicators for the quality of our clustering than any of the resultant graph projections, such as PCA or tSNE. 

![Stability Plots]({{site.baseurl}}{% link topics/transcriptomics/images/raceid_sat_jacc.png %} "Stability Plots")

The first plot measures the levels of dispersion within each cluster and produces the mean over all clusters (as defined by the k parameter). As *k* increases, the reduction in this dispersion is measured for each increase of *k* until the change in the mean within-cluster dispersion no longer changes. Here we can see that reduction saturates at k=12, which is chosen to the be the number of clusters detected in our data for all further analysis. The second plot is the same as the first but with the actual dispersion plotted instead of the relative change of dispersion.

The third plot measures the direct stability of each of the derived (in this case, 12) clusters using the [Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index), which is a fraction that measures how many elements in two sets overlap divided by the union of both sets. Here the (top N) genes expressed by the cells in each cluster are intersected with cells in all other clusters to measure how unique the expression profile is to that cluster. The scales given by the plot are actually measuring the dissimalirity between sets, which is one minus the index.

Ideally the Jaccard should have above 0.6 in most clusters, but it is acceptable to have one or two more poorly defined clusters

----

### Outlier Detection

Outlier detection attempts to refine the initially detected clusters to find smaller (sub-)clusters that could be used to define rarer cell types.

The next three plots attempts to do this by describing the variation of the gene expression, given by; a Background plot, a Sensitivity plot, and an Outlier probability plot.


![Gene Expression Plots]({{site.baseurl}}{% link topics/transcriptomics/images/raceid_gexpr.png %} "Stability Plots")


1. A background model is calibrated and outliers are identified based on the distribution of transcript counts within a cluster. The counts for each gene are assumed to follow a negative binomial distribution determined by a mean (average expression of a gene across all cells in a cluster), and a dispersion parameter. The dispersion is dervied from the average variance-mean dependence, modelled as a logarithmic second order polynomial under the assumption that *most* genes are not differentially expressed between clusters, and that true biological variability should exceed this assumption.
   As we can see from the Background plot, the (red) regression of the variance on the mean (as approximated by a second-order polynomial in logarithmic space) is higher than the variance or most genes (all grey dots below the red curve) as expected, since they are not differentially expressed. The genes above this regression are therefore significant for the detection of outlier cells. The orange line is the local regression (moving average variance per mean) and is used purely for illustrative purposes.

1. Outlier cells are detected if the probability for that cell $$c$$, a minimum number of genes $$G_{min}$$ of observing total counts $$T_G_{min}$$ is less than a specific threshold $$P_{\textnormal{thr}}$$. To summarize formally: $$P(\sum_{g\epsilon \textnormal{ outlg}}T_{(c,g)}) < \textnormal{probthr}$$
  This is shown in the chart below as the number of outliers as a function of the probablity threshold, which is set to $$1e-3$$ by default. Ideally, this threshold should be chosen to that the tail of the distribution is captured as outliers to ensure a maximum sensitivity of this method. If the sensitivity of the sequencing was low, then only a few highly expressed genes would be reliably quantified, so the outlier probability threshold would need to be higher (e.g. up to 1).

1. A barplot of the outlier probabilities of all cells across all clusters. All outlier cells are merged into their own clusters if their similarity exceeds a quantile threshold of the similarity distribution for all pairs of cells within one of the original clusters. After the outlier cells are merged, then new cluster centers are defined for the original clusters after removing the outliers. Then, each cell is assigned to the nearest cluster center using k-partitioning.

The most differentially expressed genes (below a maximum p-value cutoff) in each of the clusters can be seen in the output file *Clustering using RaceID on data: Cluster - Genes per Cluster*, with 7 columns describing:

 | Column | Description |
 |--------|-------------|
 | **.** | Gene Name |
 | **n** | Cluster number |
 | **mean.ncl** | Mean expression of gene across cells *outside* the cluster |
 | **mean.cl** | Mean expression of gene across cells *inside* the cluster |
 | **fc** | Fold-change of mean expression in the cluster, against all remaining cells |
 | **pv** | Inferred p-value for differential expression |
 | **padj** | Adjusted p-value using the Benjami-Hochberg correction for the false discovery rate |


---
### Heatmaps

The remainder of the plots are heatmaps derived from k-medoids clustering, showing the similarity between clusters for both the initial clustering and the final (post outlier detection) clustering.

![Heatmaps]({{site.baseurl}}{% link topics/transcriptomics/images/raceid_heatmaps.png %} "Heatmaps for initial and final clusters")


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


### Clustering

> ### {% icon hands_on %} Hands-on: Task description
>

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

### Assessing the Quality of the Clusters

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Cluster Inspection using RaceID** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input RaceID RDS"*: `outrdat` (output of **Clustering using RaceID** {% icon tool %})
>    - *"Plot All Clusters?"*: `Yes`
>    - *"Perform Subset Analysis?"*: `No`
>    - *"Examine Genes of Interest"*: `No`
>    - *"Differential Gene Testing"*: `No`
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

### Observations to take into account

- Mention cell-cycle variation and transcriptional bursting (in slides?)


## Refined Approach


### Filtering, Normalisation, and Confounder Removal

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Filtering, Normalisation, and Confounder Removal using RaceID** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Count Matrix"*: `output` (Input dataset)
>    - In *"Filtering"*:
>        - *"Use Defaults?"*: `Yes`
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


### Clustering

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Clustering using RaceID** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input RaceID RDS"*: `outrdat` (output of **Filtering, Normalisation, and Confounder Removal using RaceID** {% icon tool %})
>    - In *"Clustering"*:
>        - *"Use Defaults?"*: `Yes`
>    - In *"Outliers"*:
>        - *"Use Defaults?"*: `Yes`
>    - In *"tSNE and FR"*:
>        - *"Use Defaults?"*: `Yes`
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

### Assessing the Quality of the Clusters

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Cluster Inspection using RaceID** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input RaceID RDS"*: `outrdat` (output of **Clustering using RaceID** {% icon tool %})
>    - *"Plot All Clusters?"*: `Yes`
>    - *"Perform Subset Analysis?"*: `No`
>    - *"Examine Genes of Interest"*: `No`
>    - *"Differential Gene Testing"*: `No`
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



## Learning to Play with the Data

blahblahblah

# Cluster Inspection

## Differential Gene Analysis Between Two Clusters

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Cluster Inspection using RaceID** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input RaceID RDS"*: `outrdat` (output of **Clustering using RaceID** {% icon tool %})
>    - *"Plot All Clusters?"*: `No`
>    - *"Perform Subset Analysis?"*: `No`
>    - *"Examine Genes of Interest"*: `No`
>    - *"Differential Gene Testing"*: `Yes`
>        - In *"Cells in Set A"*:
>            - *"Name of Set"*: `Cells in 1`
>            - *"Selection method"*: `Cluster Numbers`
>                - *"List of clusters"*: `1`
>        - In *"Cells in Set B"*:
>            - *"Name of Set"*: `Cells in 3`
>            - *"Selection method"*: `Cluster Numbers`
>                - *"List of clusters"*: `3`
>        - *"Use Defaults?"*: `Yes`
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

## Differential Gene Expression Across All Clusters

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Cluster Inspection using RaceID** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input RaceID RDS"*: `outrdat` (output of **Clustering using RaceID** {% icon tool %})
>    - *"Plot All Clusters?"*: `No`
>    - *"Perform Subset Analysis?"*: `No`
>    - *"Examine Genes of Interest"*: `Yes`
>        - *"Genes to Examine"*: `Ptma,Rps2`
>        - *"Use Defaults?"*: `Yes`
>    - *"Differential Gene Testing"*: `No`
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

# Lineage Tree Construction

## Lineage computation

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Lineage computation using StemID** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input RDS"*: `outrdat` (output of **Clustering using RaceID** {% icon tool %})
>    - In *"Compute transcriptome entropy of each cell"*:
>        - *"Use Defaults?"*: `Yes`
>    - In *"Compute Cell Projections for Randomized Background Distribution"*:
>        - *"Use Defaults?"*: `Yes`
>    - In *"StemID2 Lineage Graph"*:
>        - *"Use Defaults?"*: `Yes`
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

## Examining a Trajectory with StemID

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Lineage Branch Analysis using StemID** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input RDS"*: `outrdat` (output of **Lineage computation using StemID** {% icon tool %})
>    - In *"StemID Branch Link Examine"*:
>        - *"Perform StemID?"*: `Yes`
>            - *"Trajectory Path i, j, k"*: `1,3,5`
>            - *"Use Defaults?"*: `Yes`
>    - In *"FateID Branch Link Examine"*:
>        - *"Perform FateID?"*: `No`
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

## Examining a Trajectory with FateID

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Lineage Branch Analysis using StemID** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input RDS"*: `outrdat` (output of **Lineage computation using StemID** {% icon tool %})
>    - In *"StemID Branch Link Examine"*:
>        - *"Perform StemID?"*: `No`
>    - In *"FateID Branch Link Examine"*:
>        - *"Perform FateID?"*: `Yes`
>            - *"Cells from Clusters"*: `1,3,5`
>            - *"Use Defaults?"*: `Yes`
>            - *"Perform Additional FateID Analysis with Self-Organised Map?"*: `No`
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