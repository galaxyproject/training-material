---
layout: tutorial_hands_on
title: RNA-seq counts to genes
zenodo_link: "https://figshare.com/s/f5d63d8c265a05618137"
enable: "false"
tags:
  - limma-voom
  - mouse
  - QC
questions:
  - "What are the differentially expressed genes in the mammary gland of pregnant versus lactating mice?"
  - "How to analyze RNA count data using limma-voom?"
  - "How to perform quality control of RNA-seq count data?"
objectives:
  - "Analysis of RNA-seq count data using limma-voom"
  - "Quality control of count data"
  - "Visualisation and interactive exploration of count data"
  - "Identification of differentially expressed genes"
time_estimation: "3h"
key_points:
  - "The limma-voom tool can be used to perform differential expression and output useful plots"
  - "Multiple comparisons can be input and compared"
  - "Results can be interactively explored with limma-voom via Glimma"
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

There are many steps involved in analysing an RNA-Seq experiment. The analysis begins with sequencing reads (FASTQ files). These are usually aligned to a reference genome, if available. Then the number of reads mapped to each gene can be counted. This results in a table of counts, which is what we perform statistical analyses on to determine differentially expressed genes and pathways. The purpose of this tutorial is to demonstrate how to perform differential expression on count data with **limma-voom**. How to generate counts from reads (FASTQs) is covered in the accompanying tutorial [RNA-seq reads to counts]({{ site.baseurl }}/topics/transcriptomics/tutorials/rna-seq-reads-to-counts/tutorial.html).

**Mouse mammary gland dataset**

The data for this tutorial comes from a Nature Cell Biology paper, [EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival](https://www.ncbi.nlm.nih.gov/pubmed/25730472)), Fu et al. 2015. Both the raw data (sequence reads) and processed data (counts) can be downloaded from Gene Expression Omnibus database (GEO) under accession number [GSE60450](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450).

This study examined the expression profiles of basal and luminal cells in the mammary gland of virgin, pregnant and lactating mice. Six groups are present, with one for each combination of cell type and mouse status. Note that two biological replicates are used here, two independent sorts of cells from the mammary glands of virgin, pregnant or lactating mice, however three replicates is usually recommended as a minimum requirement for RNA-seq. In this tutorial we will use the GEO counts file as a starting point for our analysis. Alternatively, you could create a count matrix from the raw sequence reads, as demonstrated in the [RNA-seq reads to counts tutorial]({{ site.baseurl }}/topics/transcriptomics/tutorials/rna-seq-reads-to-counts/tutorial.html). The GEO count file was generated from aligning the reads to the mouse `mm10` genome with the [Rsubread](https://www.biorxiv.org/content/early/2018/08/15/377762) aligner, followed by counting reads mapped to RefSeq genes with [featureCounts](https://academic.oup.com/bioinformatics/article/30/7/923/232889) (Liao, Smyth, and Shi 2014), see the [Fu paper](https://www.nature.com/articles/ncb3117) for details.

We will use **limma-voom** for identifying differentially expressed genes here. Other popular alternatives are edgeR and DESeq2. Limma-voom has been shown to be perform well in terms of precision, accuracy and sensitivity ([Costa-Silva, Domingues and Lopes 2017](https://www.ncbi.nlm.nih.gov/pubmed/29267363)) and, due to its speed, it's particularly recommended for large-scale datasets with 100s of samples ([Chen, Lun, Smyth 2016](https://f1000research.com/articles/5-1438/v2)).

This is a Galaxy tutorial based on material from the [COMBINE R RNAseq workshop](http://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html), first taught [here](http://combine-australia.github.io/2016-05-11-RNAseq/).

![Tutorial Dataset](../../images/rna-seq-reads-to-counts/mouse_exp.png "Tutorial Dataset")


> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}


{% include snippets/warning_results_may_vary.md %}

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

![seqdata file](../../images/rna-seq-counts-to-genes/seqdata.png "Count file (before formatting)"){: width="50%"}

The `sampleinfo` file contains basic information about the samples that we will need for the analysis. See below.

![sampleinfo file](../../images/rna-seq-counts-to-genes/sampleinfo.png "Sample information file (before formatting)"){: width="50%"}

## Format the data

Let’s create a new file, `countdata`, that contains only the counts for the 12 samples i.e. we'll remove the gene length column with the **Cut columns from a table (cut)** tool. The sample names are also pretty long so we'll use the **Replace Text in entire line** tool to shorten these to contain only the relevant information about each sample.

> ### {% icon hands_on %} Hands-on: Format the counts data
>
> 1. **Cut columns from a table (cut)** {% icon tool %} with the following parameters:
>      - {% icon param-file %} *"File to cut"*: `seqdata`
>      - {% icon param-select %} *"Operation"*: `Discard`
>      - {% icon param-select %} *"List of fields"*: Select `Column:2`
> 2. **Replace Text in entire line** {% icon tool %} with the following parameters:
>      - {% icon param-file %} *"File to process"*: output of **Cut** {% icon tool %}
>      - {% icon param-text %} *"Find pattern"*: `_B[A-Z0-9_]+`
> 3. Rename file as `countdata` using the {% icon galaxy-pencil %} (pencil) icon. The file should look like below.
>    ![countdata file](../../images/rna-seq-counts-to-genes/countdata.png "Count file (after formatting)")
{: .hands_on}

Next, let's create a new file, `factordata`, that contains the groups information that we need for the limma-voom tool. We'll combine the cell type and mouse status to make 6 groups e.g. we'll combine the CellType `basal` with the Status `pregnant` for the group `basalpregnant`. We'll use the **Merge Columns** tool to combine the cell type and mouse status columns in the sample information file, making a column with the 6 group names.

> ### {% icon hands_on %} Hands-on: Format the sample information file
>
> 1. **Merge Columns together** {% icon tool %} with the following parameters:
>      - {% icon param-file %} *"Select data"*: `sampleinfo`
>      - {% icon param-select %} *"Merge column"*: `Column: 3`
>      - {% icon param-select %} *"with column"*: `Column: 4`
> 2. **Cut columns from a table (cut)** {% icon tool %} with the following parameters:
>      - {% icon param-file %} *"File to cut"*: output of **Merge Columns** {% icon tool %}
>      - {% icon param-select %} *"Operation"*: `Keep`
>      - {% icon param-select %} *"List of fields"*: Select `Column:2` and `Column:5`
> 3. Rename file as `factordata` using the {% icon galaxy-pencil %} (pencil) icon. The file should look like below.
>    ![factordata file](../../images/rna-seq-counts-to-genes/factordata.png "Sample information file (after formatting)")
{: .hands_on}

## Get gene annotations

Optionally, gene annotations can be provided to the limma-voom tool and if provided the annotation will be available in the output files. We'll get gene symbols and descriptions for these genes using the Galaxy **annotateMyIDs** tool, which provides annotations for human, mouse, fruitfly and zebrafish.

> ### {% icon hands_on %} Hands-on: Get gene annotations
>
> 1. **annotateMyIDs** {% icon tool %} with the following parameters:
>      - {% icon param-file %} *"File with IDs"*: `countdata`
>      - {% icon param-check %} *"File has header"*: `Yes`
>      - {% icon param-select %} *"Organism"*: `Mouse`
>      - {% icon param-select %} *"ID Type"*: `Entrez`
>      - {% icon param-check %} "*Output columns"*: tick
>          - `ENTREZID`
>          - `SYMBOL`
>          - `GENENAME`
> 2. Rename file as `annodata` using the {% icon galaxy-pencil %} (pencil) icon. The file should look like below.
>    ![annodata file](../../images/rna-seq-counts-to-genes/annodata.png "Gene annotation file"){: width="50%"}
{: .hands_on}

# Differential expression with limma-voom

## Filtering to remove lowly expressed genes

It is recommended to filter for lowly expressed genes when running the limma-voom tool. Genes with very low counts across all samples provide little evidence for differential expression and they interfere with some of the statistical approximations that are used later in the pipeline. They also add to the multiple testing burden when estimating false discovery rates, reducing power to detect differentially expressed genes. These genes should be filtered out prior to further analysis.

There are a few ways to filter out lowly expressed genes. When there are biological replicates in each group, in this case we have a sample size of 2 in each group, we favour filtering on a minimum counts-per-million (CPM) threshold present in at least 2 samples. Two represents the smallest sample size for each group in our experiment. In this dataset, we choose to retain genes if they are expressed at a CPM above 0.5 in at least two samples. The CPM threshold selected can be compared to the raw count with the CpmPlots (see below).

> ### {% icon details %} More details on filtering
>
> The limma tool uses the `cpm` function from the edgeR package [Robinson, McCarthy, and Smyth 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/) to generate the CPM values which can then be filtered. Note that by converting to CPMs we are normalizing for the different sequencing depths for each sample. A CPM of 0.5 is used as it corresponds to a count of 10-15 for the library sizes in this data set. If the count is any smaller, it is considered to be very low, indicating that the associated gene is not expressed in that sample. A requirement for expression in two or more libraries is used as each group contains two replicates. This ensures that a gene will be retained if it is only expressed in one group. Smaller CPM thresholds are usually appropriate for larger libraries. As a general rule, a good threshold can be chosen by identifying the CPM that corresponds to a count of 10, which in this case is about 0.5. You should filter with CPMs rather than filtering on the counts directly, as the latter does not account for differences in library sizes between samples.
{: .details}

## Normalization for composition bias

In an RNA-seq analysis, the counts are normalized for different sequencing depths between samples. Normalizing to eliminate composition biases between samples is also typically performed. Composition biases can occur, for example, if there are a few highly expressed genes dominating in some samples, leading to less reads from other genes. By default, TMM normalization [(Robinson and Oshlack 2010)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2864565/) is performed by the limma tool using the edgeR `calcNormFactors` function (this can be changed under **Advanced Options**). TMM stands for Trimmed Mean of M values, where a weighted trimmed mean of the log expression ratios is used to scale the counts for the samples. See the figure from the TMM paper below. Note the plot (Figure 1c) that shows how a few highly expressed genes in the liver sample (where the arrow is) results in the majority of other genes in the sample having the appearance of being expressed lower in liver. The mid-line through the points is offset from the expected zero and the TMM normalization factor (red line) scales the counts to adjust for this.

![TMM normalization](../../images/rna-seq-counts-to-genes/TMM.png "TMM normalization (Robinson and Oshlack 2010)"){: width="50%"}

## Specify Contrast(s) of interest

Since we are interested in differences between groups, we need to specify which comparisons we want to test. For example, if we are interested in knowing which genes are differentially expressed between the pregnant and lactating group in the basal cells we specify `basalpregnant-basallactate` for the *Contrast of Interest*. Note that the group names in the contrast must exactly match the names of the groups in the `factordata` file. More than one contrast can be specified using the `Insert Contrast` button, so we could look at more comparisons of the groups here, but first we'll take a look at `basalpregnant-basallactate`.

> ### {% icon hands_on %} Hands-on: Differential expression with limma-voom
>
> 1. **limma** {% icon tool %} with the following parameters:
>      - {% icon param-file %} *"Differential Expression Method"*: `limma-voom`
>      - {% icon param-select %} *"Count Files or Matrix?*": `Single Count Matrix`
>          - {% icon param-file %} *"Count Matrix"*: Select `countdata`
>      - {% icon param-select %} *"Input factor information from file?"*: `Yes`
>          - {% icon param-file %} *"Factor File"*: Select `factordata`
>      - {% icon param-select %} *"Use Gene Annotations?"*: `Yes`
>          - {% icon param-file %} *"Factor File"*: Select `annodata`
>      - {% icon param-text %} *"Contrast of Interest"*: `basalpregnant-basallactate`
>      - {% icon param-select %} *"Filter lowly expressed genes?"*: `Yes`
>          - {% icon param-select %} *"Filter on CPM or Count values?"*: `CPM`
>          - {% icon param-text %} *"Minimum CPM"*: `0.5`
>          - {% icon param-text %} *"Minimum Samples"*: `2`
> 2. Inspect the `Report` produced by clicking on the {% icon galaxy-eye %} (eye) icon
{: .hands_on}

# Quality Control of count data

Before we check out the differentially expressed genes, we can look at the `Report` information to check that the data is good quality and that the samples are as we would expect.

## Multidimensional scaling plot

By far, one of the most important plots we make when we analyse RNA-Seq data are MDS plots. An MDS plot is a visualisation of a principal components analysis, which determines the greatest sources of variation in the data. A principal components analysis is an example of an unsupervised analysis, where we don’t need to specify the groups. If your experiment is well controlled and has worked well, what we hope to see is that the greatest sources of variation in the data are the treatments/groups we are interested in. It is also an incredibly useful tool for quality control and checking for outliers. This Galaxy limma tool outputs an MDS plot by default in the `Report` and a link is also provided to a PDF version (`MDSPlot_CellTypeStatus.pdf`). A scree plot is also produced that shows how much variation is attributed to each dimension. If there was a batch effect for example, you may see high values for additional dimensions. The limma tool plots the first two dimensions by default (1 vs 2), however you can also plot additional dimensions 2 vs 3 and 3 vs 4 using under **Output Options** Additional Plots `MDS Extra` These are displayed in the `Report` along with a link to a PDF version (`MDSPlot_extra.pdf`).

![MDS Plot](../../images/rna-seq-counts-to-genes/mdsscree.png "MDS Plot")

Take a look at the MDS plot coloured by group.

> ### {% icon question %} Question
>
> Do you notice anything about the samples in this plot?
>
>    > ### {% icon solution %} Solution
>    >
>    > Two samples don't appear to be in the right place.
>    >
>    {: .solution}
{: .question}

It turns out that there has been a mix-up with two samples, they have been mislabelled in the sample information file. This shows how the MDS plot can also be useful to help identify if sample mix-ups may have occurred. We need to redo the limma-voom analysis with the correct sample information.

> ### {% icon hands_on %} Hands-on: Use the Rerun button to redo steps
>
> 1. Import the correct sample information file from `https://ndownloader.figshare.com/files/5999832?private_link=1d788fd384d33e913a2a`
> 2. Use the Rerun button in the History to redo the **Merge Columns** and **Cut** steps on the correct sample information file.
> 3. Delete the incorrect sample information datasets to avoid any confusion.
> 4. Rerun **limma** as before with the correct `sampleinfo` file and adding the following parameters:
>      - **Output Options**
>          - {% icon param-check %} *"Additional Plots"* tick:
>              - `Density Plots (if filtering)`
>              - `CpmsVsCounts Plots (if filtering on cpms)`
>              - `Box Plots (if normalising)`
>              - `MDS Extra (Dims 2vs3 and 3vs4)`
>              - `MD Plots for individual samples`
>              - `Heatmaps (top DE genes)`
>              - `Stripcharts (top DE genes)`
>          - {% icon param-check %} *"Output Library information file?"*: `Yes`
{: .hands_on}

In the `Report` you should then see the correct MDS plot as below.

![MDS Plot](../../images/rna-seq-counts-to-genes/mdsscree_corr.png "MDS Plot (correct samples)")

> ### {% icon details %} More details on MDS plots
>
> The distance between each pair of samples in the MDS plot is calculated as the leading fold change, defined as the root-mean-square of the largest 500 log2-fold changes between that pair of samples. Replicate samples from the same group cluster together in the plot, while samples from different groups form separate clusters. This indicates that the differences between groups are larger than those within groups, i.e., differential expression is greater than the variance and can be detected. In the MDS plot, the distance between basal samples on the left and luminal cells on the right is about 6 units, corresponding to a leading fold change of about 64-fold (2^6 = 64) between basal and luminal. The expression differences between virgin, pregnant and lactating are greater for luminal cells than for basal.
>
> Clustering in the MDS plot can be used to motivate changes to the analysis in light of potential batch effects. For example, imagine that the first replicate of each group was prepared at a separate time from the second replicate. If the MDS plot showed separation of samples by time, it might be worthwhile including time as an additional factor in the differential expression analysis, to account for the time-based effect.
{: .details}

> ### {% icon question %} Question
>
> What is the greatest source of variation in the data (i.e. what does dimension 1 represent)?
> What is the second greatest source of variation in the data?
>
>    > ### {% icon solution %} Solution
>    >
>    > Dimension 1 represents the variation due to cell type, basal vs luminal. Dimension 2 represents the variation due to the stages, virgin, pregnant or lactating.
>    >
>    {: .solution}
{: .question}

Next, scroll down the `Report` to take a look at the **Additional information** and **Summary of experimental data** sections near the bottom. It should look similar to below. Here you can check that the correct samples have been assigned to the correct groups, what settings were used (e.g. filters, normalization method) and also how many genes were filtered out due to low expression.

![Report summary](../../images/rna-seq-counts-to-genes/report_summary.png){: width="600px"}

> ### {% icon question %} Question
>
> How many genes have been filtered out for low expression?
>
>    > ### {% icon solution %} Solution
>    >
>    > 11375 genes were filtered out as insignificant as they were without more than 0.5 CPM in at least 2 samples.
>    >
>    {: .solution}
{: .question}


## Density plots

Density plots can be output in the `Report` if *Filter lowly expressed genes* is selected. A link is also provided in the `Report` to a PDF version (`DensityPlots.pdf`). These plots allow comparison of the counts distributions before and after filtering. The samples are coloured by the groups. Count data is not normally distributed, so if we want to examine the distributions of the raw counts we need to log the counts. We typically check the distribution of the read counts on the log2 scale. A CPM value of 1 is equivalent to a log-CPM value of 0 and the CPM we used of 0.5 is equivalent to a log-CPM of -1. It can be seen in the Raw counts (before filtering) plot below, that a large proportion of genes within each sample are not expressed or lowly-expressed, and the Filtered counts plot shows our filter of CPM of 0.5 (in at least 2 samples) removes a lot of these uninformative genes.

![Density Plots](../../images/rna-seq-counts-to-genes/densityplots.png "Density Plots")

We can also have a look more closely to see whether our threshold of 0.5 CPM does indeed correspond to a count of about 10-15 reads in each sample with the plots of CPM versus raw counts.

The `Report` provides links to PDFs of all plots shown in the `Report` and also to the rest of the additional plots selected to be output.

![Report Outputs](../../images/rna-seq-counts-to-genes/report_plots.png "Report outputs")

Click on the `CpmPlots.pdf` link in the `Report`. You should see 12 plots, one for each sample. Two of the plots are shown below. From these plots we can see that 0.5 CPM is equivalent to ~10 counts in each of the 12 samples, so 0.5 seems to be an appropriate threshold for this dataset (these samples all have sequencing depth of 20-30 million, see the `Library information` file below, so a CPM value of 0.5 would be ~10 counts).

![CPM threshold Plots](../../images/rna-seq-counts-to-genes/cpmsvscounts.png "CPM vs Raw Counts Plots")

> ### {% icon tip %} Tip
>
> * A threshold of 1 CPM in at least minimum group sample size is a good rule of thumb for samples with about 10 million reads. For larger library sizes increase the CPM theshold and for smaller library sizes decrease it. Check the CpmPlots to see if the selected threshold looks appropriate for the samples (equivalent to ~10 reads).
>
{: .tip}


## Box plots

We can also use box plots to check the distributions of counts in the samples. Box plots can be selected to be output by the Galaxy limma-voom tool if normalization is applied (TMM is applied by default). The plots are output in the `Report` and a link is also provided to a PDF version (`BoxPlots.pdf`). The samples are coloured by the groups. With the box plots for these samples we can see that overall the distributions are not identical but still not very different. If a sample is really far above or below the blue horizontal line we may need to investigate that sample further.

![Box Plots](../../images/rna-seq-counts-to-genes/boxplots.png "Box Plots")

> ### {% icon question %} Question
>
> Compare the box plots before and after TMM normalisation. Can you see any differences?
>
>    > ### {% icon solution %} Solution
>    >
>    > After the normalization more of the samples are closer to the median horizontal line.
>    >
>    {: .solution}
{: .question}


The TMM normalization generates normalization factors, where the product of these factors and the library sizes defines the effective library size. TMM normalization (and most scaling normalization methods) scale relative to one sample. The normalization factors multiply to unity across all libraries. A normalization factor below one indicates that the library size will be scaled down, as there is more suppression (i.e., composition bias) in that library relative to the other libraries. This is also equivalent to scaling the counts upwards in that sample. Conversely, a factor above one scales up the library size and is equivalent to downscaling the counts. We can see the normalization factors for these samples in the `Library information` file that we selected to output. Click on the {% icon galaxy-eye %} (eye) icon to view.

![Library Info file](../../images/rna-seq-counts-to-genes/libinfo.png "Library information file")

> ### {% icon question %} Question
>
> Which sample has the largest normalization factor? Which sample has the smallest?
>
>    > ### {% icon solution %} Solution
>    >
>    > MCL1.LA has the largest normalization factor and MCL1.LE the smallest.
>    >
>    {: .solution}
{: .question}

## MD plots for samples

It is considered good practice to make mean-difference (MD) plots for all the samples as a quality check, as described in this [edgeR workflow article](https://f1000research.com/articles/5-1438/v2). These plots allow expression profiles of individual samples to be explored more closely. An MD plot shows the log-fold change between a sample against the average expression across all the other samples. This visualisation can help you see if there are genes highly upregulated or downregulated in a sample. If we look at mean difference plots for these samples, we should be able to see the composition bias problem. The mean-difference plots show average expression (mean: x-axis) against log-fold-changes (difference: y-axis).

Click on the `MDPlots_Samples.pdf` link in the `Report`. You should see 12 MD plots, one for each sample. Let's take a look at the plots for the two samples MCL1.LA and MCL1.LE that had the largest and smallest normalization factors. The MD plots on the left below show the counts normalized for library size and the plots on the right show the counts after the TMM normalization has been applied. MCL1.LA had the largest normalization factor and was above the median line in the unnormalized by TMM box plots. MCL1.LE had the smallest normalization factor and was below the median line in the box plots. These MD plots help show the composition bias problem has been addressed.
![MD Plot LA](../../images/rna-seq-counts-to-genes/mdsampleLA.png "MD Plots for MCL1.LA before and after TMM normalization")
![MD Plot LE](../../images/rna-seq-counts-to-genes/mdsampleLE.png "MD Plots for MCL1.LE before and after TMM normalization")


## Voom variance plot

This plot is generated by the voom method and displayed in the `Report` along with a link to a PDF version (`VoomPlot.pdf`). It shows the mean-variance relationship of the genes in the dataset. It can help show if low counts have been filtered adequately and if there is a lot of variation in the data.

![Voom Plot](../../images/rna-seq-counts-to-genes/voomplot.png "Voom Plot")

> ### {% icon details %} More details on Voom variance plots
>
> If we didn't filter this dataset for the lowly expressed genes the variance plot would look like below.
>
>   ![Voom Nofilter Plot](../../images/rna-seq-counts-to-genes/voomplot_nofilt.png "Voom Plot unfiltered counts")
>
> If we look at the plot generated with the two samples mixed up we can see there's more variation.
>
>   ![Voom Mixedup Plot](../../images/rna-seq-counts-to-genes/voomplot_mixed.png "Voom Plot mixed-up samples")
>
> More examples of the variation this plot can show can be seen in Figure 1 from the [limma-voom](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29) paper, shown below.
>
>   ![Voom Plot Examples](../../images/rna-seq-counts-to-genes/voom_variance_examples.jpg "Voom Plot more examples")
>
>   *Figure 1: Mean-variance relationships. Gene-wise means and variances of RNA-seq data are represented by black points with a LOWESS trend. Plots are ordered by increasing levels of biological variation in datasets. (a) voom trend for HBRR and UHRR genes for Samples A, B, C and D of the SEQC project; technical variation only. (b) C57BL/6J and DBA mouse experiment; low-level biological variation. (c) Simulation study in the presence of 100 upregulating genes and 100 downregulating genes; moderate-level biological variation. (d) Nigerian lymphoblastoid cell lines; high-level biological variation. (e) Drosophila melanogaster embryonic developmental stages; very high biological variation due to systematic differences between samples. (f) LOWESS voom trends for datasets (a)–(e). HBRR, Ambion’s Human Brain Reference RNA; LOWESS, locally weighted regression; UHRR, Stratagene’s Universal Human Reference RNA*.
{: .details}

## MD and Volcano plots for DE results

Genome-wide plots that are useful for checking differentially expressed (DE) results are MD plots (or MA plots) and Volcano plots. There are functions in limma for generating these plots and they are used by this tool. These plots are output by default and shown in the `Report` along with a link to PDF versions (`MDPlot_basalpregnant-basallactate.pdf` and `VolcanoPlot_basalpregnant-basallactate.pdf`). In the volcano plot the top genes (by adjusted p-value) are highlighted. The number of top genes is 10 by default and the user can specify the number of top genes to view (up to 100) under **Advanced Options**.

![MDVol Plot](../../images/rna-seq-counts-to-genes/mdvolplot_basalpregnant-basallactate.png "MD Plot and Volcano Plot")

The MD Plot highlighted genes are significant at an adjusted p-value (adj.P) threshold of 0.05 and exhibit log2-fold-change (lfc) of at least 0. These thresholds can be changed under **Advanced Options**.


> ### {% icon question %} Question
>
> How many genes are differentially expressed at the default thresholds of adj.P=0.05 and lfc=0?
>
>    > ### {% icon solution %} Solution
>    >
>    > The number of DE genes at these adj.P and lfc thresholds is shown in the table in the `Report` as below.
>    >
>    > ![DE counts](../../images/rna-seq-counts-to-genes/decounts.png){: width="500px"}
>    >
>    {: .solution}
{: .question}


> ### {% icon comment %} Comment on adjusted P value
>
> **A note about deciding how many genes are significant**: In order to decide which genes are differentially expressed, we usually take a cut-off (e.g. 0.05 or 0.01) on the adjusted p-value, NOT the raw p-value. This is because we are testing many genes (more than 15000 genes here), and the chances of finding differentially expressed genes is very high when you do that many tests. Hence we need to control the false discovery rate, which is the adjusted p-value column in the results table. What this means is that, if we choose an adjusted p-value cut-off of 0.05, and if 100 genes are significant at a 5% false discovery rate, we are willing to accept that 5 will be false positives.
{: .comment}

# Testing relative to a threshold (TREAT)

When there is a lot of differential expression, sometimes we may want to cut-off on a fold change threshold, as well as a p-value threshold, so that we follow up on the most biologically significant genes. However, it is not recommended to simply rank by p-value and then discard genes with small logFC’s, as this has been shown to increase the false discovery rate. In other words, you are not controlling the false discovery rate at 5% any more. There is a function called `treat` in limma that performs this style of analysis correctly [(McCarthy and Smyth 2009)](https://www.ncbi.nlm.nih.gov/pubmed/19176553). TREAT will simply take a user-specified log fold change cut-off and recalculate the moderated t-statistics and p-values with the new information about logFC. There are thousands of genes differentially expressed in this `basalpregnant-basallactate` comparison, so let's rerun the analysis applying TREAT and similar thresholds to what was used in the Fu paper: an adjusted P value of 0.01 (1% false discovery rate) and a log-fold-change cutoff of 0.58 (equivalent to a fold change of 1.5).

> ### {% icon hands_on %} Hands-on: Testing relative to a threshold (TREAT)
>
> 1. Rerun **limma** {% icon tool %} with the following parameters:
>      - *"**Output Options**"*
>          - {% icon param-check %} "Output Library information file?": `No`
>      - *"**Advanced Options**"*
>          - {% icon param-text %} *"Minimum Log2 Fold Change"*: `0.58`
>          - {% icon param-text %} *"P-Value Adjusted Threshold"*: `0.01`
>          - {% icon param-check %} *"Test significance relative to a fold-change threshold (TREAT)"*: `Yes`
> 2. Inspect the `Report`
{: .hands_on}

We can see that much fewer genes are now highlighted in the MD plot and identified as differentially expressed.

![TREAT MDVol Plot](../../images/rna-seq-counts-to-genes/TREAT_mdvolplot_basalpregnant-basallactate.png "TREAT MD and Volcano plots")

![TREAT DE counts](../../images/rna-seq-counts-to-genes/TREATdecounts.png "TREAT differentially expressed genes"){: width="400px"}


# Visualising results

Before following up on the DE genes with further lab work, it is recommended to have a look at the expression levels of the individual samples for the genes of interest. The Galaxy limma tool can auto-generate heatmaps of the top genes to show the expression levels across the *samples*. This enables a quick view of the expression of the top differentially expressed genes. This can help show if expression is consistent between replicates in the groups.

## Heatmap of top genes

Click on the `Heatmap_basalpregnant-basallactate.pdf` link in the `Report`. You should see a plot like below.

![Heatmap](../../images/rna-seq-counts-to-genes/heatmap.png "Heatmap of top genes")

## Stripcharts of top genes

The limma-voom tool can also auto-generate stripcharts to view the expression of the top genes across the *groups*. Click on the `Stripcharts_basalpregnant-basallactate.pdf` link in the `Report`. You should see 10 plots, one for each top gene. Four are shown below. Note that here you can see if replicates tend to group together and how the expression compares to the other groups.

![Stripchart Plot](../../images/rna-seq-counts-to-genes/stripcharts.png "Stripcharts of top genes"){: width="950px"}


## Interactive plots (Glimma)

An interactive version of the mean-difference plots is output by the limma-voom tool via the [Glimma](https://github.com/Shians/Glimma) package, if a gene annotation file is provided. A link to a html page is generated in the `Report` that allows the user to search for their favourite gene.

Click on the `Glimma_MDPlot_basalpregnant-basallactate.html` link in the `Report`. You should see a two-panel interactive MD plot like below. The left plot shows the log-fold-change vs average expression. The right plot shows the expression levels of a particular gene of each sample by groups (similar to the stripcharts). Hovering over points on left plot will plot expression level for corresponding gene, clicking on points will fix the expression plot to gene. Clicking on rows on the table has the same effect as clicking on the corresponding gene in the plot.

<iframe src="../../images/rna-seq-counts-to-genes/glimma/MD-Plot.html" width="100%" height="1000"></iframe>

> ### {% icon hands_on %} Hands-on: Search for a gene of interest
>
> `Egf` was a gene identifed as very highly expressed in the Fu paper and confirmed with qRT-PCR, see Fig. 6c from the paper below.
> ![Fu EGF gene](../../images/rna-seq-counts-to-genes/fu_egf.png "Fu et al, Nat Cell Biol 2015")
>
> Search for `Egf` in the Glimma interactive table. You should see something similar to below.
>
> ![Glimma EGF gene](../../images/rna-seq-counts-to-genes/glimma_egf.png "Glimma EGF gene")
>
>  Notice that in the plot on the right above, showing `Egf` expression for all samples, we can see it is highly expressed in the luminal lactate group but not the other samples.
{: .hands_on}

Multiple contrasts can be run with the limma tool. For example, we can compare the pregnant and lactating conditions for both the basal and luminal cells. So let's rerun the limma-voom TREAT analysis (adj.P <0.01 and lfc=0.58) and this time use the `Insert Contrast` button to include the additional contrast `luminalpregnant - luminallactate`. We can then see how the number of differentially expressed genes in the luminal cells compares to the basal cells.

> ### {% icon hands_on %} Hands-on: Run multiple contrasts
>
> 1. Rerun **limma** {% icon tool %} adding the following parameters *(i.e. run with 2 contrasts)*:
>      - {% icon param-text %} *"Contrast of Interest"*: `basalpregnant-basallactate`
>      - {% icon param-text %} *"Contrast of Interest"*: `luminalpregnant-luminallactate`
>      - *"**Advanced Options**"*
>          - {% icon param-text %} *"Minimum Log2 Fold Change"*: `0.58`
>          - {% icon param-text %} *"P-Value Adjusted Threshold"*: `0.01`
>          - {% icon param-check %} *"Test significance relative to a fold-change threshold (TREAT)"*: `Yes`
> 2. Inspect the `Report`
>
> You should find that there are more genes differentially expressed for the luminal cells than the basal. There are ~274 genes DE in basal cells versus ~ 1610 in the luminal cells.
> ![Basal Luminal DE counts](../../images/rna-seq-counts-to-genes/basal_luminal_decounts.png "Basal vs Luminal DE counts"){: width="50%"}
>
> This is similar to what Fu et al found, many more genes differentially expressed in the luminal cells on lactation, compared to the basal cells.
> ![Fu DE genes](../../images/rna-seq-counts-to-genes/fu_degenes.png "Fu et al, Nat Cell Biol 2015"){: width="25%"}
{: .hands_on}

The tables of differentially expressed genes are output as links in the `Report` (`limma-voom_basalpregnant-basallactate.tsv` and `limma-voom_luminalpregnant-luminallactate.tsv`), see below, and also as datasets in the history (`DE tables`). With multiple contrasts, a plot for each contrast is generated for relevant plots, as shown below. This enables a quick and easy visual comparison of the contrasts.

![Multiple contrasts](../../images/rna-seq-counts-to-genes/multiple_contrasts.png "Multiple contrasts output")

> ### {% icon tip %} Tip
>
> The `Report` with all the plots and tables can be downloaded by clicking on the floppy disk icon on the dataset in the history as shown below.
> ![Report download](../../images/rna-seq-counts-to-genes/download_report.png "Download limma report"){: width="20%"}
>
{: .tip}

To see some methods for identifying differentially expressed pathways in this dataset, see the follow-on tutorial [RNA-seq genes to pathways]({{ site.baseurl }}/topics/transcriptomics/tutorials/rna-seq-genes-to-pathways/tutorial.html). To see how to create a heatmap of custom genes using this dataset, see the tutorial [Visualization of RNA-Seq results with heatmap2]({{ site.baseurl }}/topics/transcriptomics/tutorials/rna-seq-viz-with-heatmap2/tutorial.html)

# Conclusion
{:.no_toc}

In this tutorial we have seen how counts files can be converted into differentially expressed genes with limma-voom. This follows on from the accompanying tutorial, [RNA-seq reads to counts]({{ site.baseurl }}/topics/transcriptomics/tutorials/rna-seq-reads-to-counts/tutorial.html), that showed how to generate counts from the raw reads (FASTQs) for this dataset. In this part we have seen ways to visualise the count data, and QC checks that can be performed to help assess the quality and results. We have also reproduced results similar to what the authors found in the original paper with this dataset. For further reading on analysis of RNA-seq count data and the methods used here, see the articles; RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR [(Law et al. 2016)](https://f1000research.com/articles/5-1408/v2) and From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline [(Chen, Lun, Smyth 2016)](https://f1000research.com/articles/5-1438/v2).
