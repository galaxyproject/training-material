---
layout: tutorial_hands_on

title: "De novo transcriptome reconstruction with RNA-Seq"
subtopic: introduction
priority: 3

zenodo_link: "https://zenodo.org/record/254485#.WKODmRIrKRu"
questions:
  - "What genes are differentially expressed between G1E cells and megakaryocytes?"
  - "How can we generate a transcriptome de novo from RNA sequencing data?"
objectives:
  - "Analysis of RNA sequencing data using a reference genome"
  - "Reconstruction of transcripts without reference transcriptome (de novo)"
  - "Analysis of differentially expressed genes"
time_estimation: "6h"
key_points:
  - "De novo transcriptome reconstruction is the ideal approach for identifying differentially expressed known and novel transcripts."
  - "Differential gene expression testing is improved with the use of replicate experiments and deep sequence coverage."
  - "Visualizing data on a genome browser is a great way to display interesting patterns of differential expression."
contributions:
  authorship:
  - malloryfreeberg
  - moheydarian
  editing:
  - jxtx
---

# Introduction


The data provided here are part of a Galaxy tutorial that analyzes RNA-seq data from a study published by *Wu et al.* in 2014 [DOI:10.1101/gr.164830.113](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4248312/). The goal of this study was to investigate "the dynamics of occupancy and the role in gene regulation of the transcription factor Tal1, a critical regulator of hematopoiesis, at multiple stages of hematopoietic differentiation." To this end, RNA-seq libraries were constructed from multiple mouse cell types including G1E - a GATA-null immortalized cell line derived from targeted disruption of GATA-1 in mouse embryonic stem cells - and megakaryocytes. This RNA-seq data was used to determine differential gene expression between G1E and megakaryocytes and later correlated with Tal1 occupancy. This dataset (GEO Accession: GSE51338) consists of biological replicate, paired-end, poly(A) selected RNA-seq libraries. Because of the long processing time for the large original files, we have downsampled the original raw data files to include only reads that align to chromosome 19 and a subset of interesting genomic loci identified by Wu *et al*.

# Analysis strategy

The goal of this exercise is to identify what transcripts are present in the G1E and megakaryocyte cellular states and which transcripts are differentially expressed between the two states. We will use a *de novo* transcript reconstruction strategy to infer transcript structures from the mapped reads in the absence of the actual annotated transcript structures. This will allow us to identify novel transcripts and novel isoforms of known transcripts, as well as identify differentially expressed transcripts.

> <agenda-title></agenda-title>
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Data upload

Due to the large size of this dataset, we have downsampled it to only include reads mapping to chromosome 19 and certain loci with relevance to hematopoeisis. This data is available at [`Zenodo`](https://zenodo.org/record/583140#.WSW3NhPyub8), where you can find the forward and reverse reads corresponding to replicate RNA-seq libraries from G1E and megakaryocyte cells and an annotation file of RefSeq transcripts we will use to generate our transcriptome database.

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this RNA-seq exercise
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
> 2. Open the data upload manager (Get Data -> Upload file)
> 3. Copy and paste the links for the reads and annotation file
> 4. Select **Paste/Fetch Data**
> 5. Paste the link(s) into the text field
> 6. Change the datatype of the read files to **fastqsanger**
> 7. Change the datatype of the annotation file to **gtf** and assign the Genome as **mm10**
> 8. Press **Start**
> 9. Rename the files in your history to retain just the necessary information (*e.g.* "G1E R1 forward reads")
>
>
> 10. Import the data from [Zenodo](https://zenodo.org/record/583140) [![DOI](https://zenodo.org/badge/DOI/10.123/GTNdenovoRNAseq.svg)](https://zenodo.org/record/583140)
>
>     ```
>     https://zenodo.org/record/583140/files/G1E_rep1_forward_read_%28SRR549355_1%29
>     https://zenodo.org/record/583140/files/G1E_rep1_reverse_read_%28SRR549355_2%29
>     https://zenodo.org/record/583140/files/G1E_rep2_forward_read_%28SRR549356_1%29
>     https://zenodo.org/record/583140/files/G1E_rep2_reverse_read_%28SRR549356_2%29
>     https://zenodo.org/record/583140/files/Megakaryocyte_rep1_forward_read_%28SRR549357_1%29
>     https://zenodo.org/record/583140/files/Megakaryocyte_rep1_reverse_read_%28SRR549357_2%29
>     https://zenodo.org/record/583140/files/Megakaryocyte_rep2_forward_read_%28SRR549358_1%29
>     https://zenodo.org/record/583140/files/Megakaryocyte_rep2_reverse_read_%28SRR549358_2%29
>     https://zenodo.org/record/583140/files/RefSeq_reference_GTF_%28DSv2%29
>     ```
>
>     You will need to fetch the link to the annotation file yourself ;)
>
>     {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
{: .hands_on}

## Quality control

For quality control, we use similar tools as described in [NGS-QC tutorial]({{site.baseurl}}/topics/sequence-analysis/): [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).

> <hands-on-title>Quality control</hands-on-title>
>
> 1. **FastQC** {% icon tool %}: Run `FastQC` on the forward and reverse read files to assess the quality of the reads.
>
>    > <question-title></question-title>
>    >
>    > 1. What is the read length?
>    > 2. Is there anything interesting about the quality of the base calls based on the position in the reads?
>    >
>    > > <solution-title></solution-title>
>    > > 1. The read length is 99 bp
>    > > 2. The quality of base calls declines throughout a sequencing run.
>    > {: .solution }
>    {: .question}
>
> 2. **Trimmomatic** {% icon tool %}: Trim off the low quality bases from the ends of the reads to increase mapping efficiency. Run `Trimmomatic` on each pair of forward and reverse reads with the following settings:
>
>    - *"Single-end or paired-end reads?"*: `Paired-end (two separate input files)`
>       - {% icon param-file %} *"Input FASTQ file (R1/first of pair)"*: `G1E_rep1 forward read`
>       - {% icon param-file %} *"Input FASTQ file (R2/second of pair)"*: `G1E_rep1 reverse read`
>    - *"Perform initial ILLUMINACLIP step?"*: `No`
>
> 3. **FastQC** {% icon tool %}: Re-run `FastQC` on trimmed reads and inspect the differences.
>
>    > <question-title></question-title>
>    >
>    > 1. What is the read length?
>    > 2. Is there anything interesting about the quality of the base calls based on the position in the reads?
>    >
>    > > <solution-title></solution-title>
>    > > 1. The read lengths range from 1 to 99 bp after trimming
>    > > 2. The average quality of base calls does not drop off as sharply at the 3' ends of reads.
>    > {: .solution }
>    {: .question}
> ![Before and after trimming comparison](../../images/transcriptomics_images/BeforeAndAfterTrimming.png)
>
> 4. **Trimmomatic** {% icon tool %}: Run `Trimmomatic` on the remaining forward/reverse read pairs with the same parameters.
>
{: .hands_on}

Now that we have trimmed our reads and are fortunate that there is a reference genome assembly for mouse, we will align our trimmed reads to the genome.

> <comment-title></comment-title>
>
> Instead of running a single tool multiple times on all your data, would you rather run a single tool on multiple datasets at once? Check out the [dataset collections](https://galaxyproject.org/tutorials/collections/) feature of Galaxy!
{: .comment}

# Mapping

To make sense of the reads, their positions within mouse genome must be determined. This process is known as aligning or 'mapping' the reads to the reference genome.

> <comment-title></comment-title>
>
> Do you want to learn more about the principles behind mapping? Follow our [training]({{site.baseurl}}/topics/sequence-analysis/)
{: .comment}

In the case of a eukaryotic transcriptome, most reads originate from processed mRNAs lacking introns. Therefore, they cannot be simply mapped back to the genome as we normally do for reads derived from DNA sequences. Instead, the reads must be separated into two categories:

- Reads contained within mature exons - these align perfectly to the reference genome
- Reads that span splice junctions in the mature mRNA - these align with gaps to the reference genome

Spliced mappers have been developed to efficiently map transcript-derived reads against genomes. [`HISAT`](https://ccb.jhu.edu/software/hisat2/index.shtml) is an accurate and fast tool for mapping spliced reads to a genome. Another popular spliced aligner is [`TopHat`](https://ccb.jhu.edu/software/tophat/index.shtml), but we will be using `HISAT` in this tutorial.

> <comment-title></comment-title>
> As it is sometimes quite difficult to determine which settings correspond to those of other programs, the following table might be helpful to identify the library type:
>
> Library type | **Infer Experiment** | **TopHat** | **HISAT** | **htseq-count** | **featureCounts**
> --- | --- | --- | --- | --- | ---
> PE | 1++,1--,2+-,2-+ | FR Second Strand | FR | yes | 1
> PE | 1+-,1-+,2++,2-- | FR First Strand | RF | reverse | 2
> SE | ++,-- | FR Second Strand | F | yes | 1
> SE | +-,-+ | FR First Strand | R | reverse | 2
> SE,PE | undecided | FR Unstranded | default | no | 0
>
{: .comment}

> <hands-on-title>Spliced mapping</hands-on-title>
>
> 1. **HISAT2** {% icon tool %}: Run `HISAT2` on one forward/reverse read pair and modify the following settings:
>    - *"Source for the reference genome"*: `Use a built-in genome`
>       - {% icon param-file %} *"Select a reference genome"*: `Mouse (Mus Musculus): mm10`
>    - *"Single-end or paired-end reads?"*: `Paired-end`
>       - {% icon param-file %} *"FASTA/Q file #1"*: Trimmomatic on G1E_rep1 forward read (R1 paired)
>       - {% icon param-file %} *"FASTA/Q file #2"*: Trimmomatic on G1E_rep1 reverse read (R2 paired)
>       - *"Specify strand information"*: `Forward(FR)`
>    - *"Advanced options"*
>       - *"Spliced alignment options"*
>           - *"Penalty for non-canonical splice sites"*: '3'
>           - *"Penalty function for long introns with canonical splice sites"*: 'Constant [f(x) = B]''
>           - *"Constant term (B)"*: '0.0'
>           - *"Penalty function for long introns with non-canonical splice sites"*: 'Constant [f(x) = B]''
>           - *"Transcriptome assembly reporting"*: 'Report alignments tailored for transcript assemblers including StringTie'
>
> 2. **HISAT2** {% icon tool %}: Run `HISAT2` on the remaining forward/reverse read pairs with the same parameters.
>
{: .hands_on}

# De novo transcript reconstruction
Now that we have mapped our reads to the mouse genome with `HISAT`, we want to determine transcript structures that are represented by the aligned reads. This is called *de novo* transcriptome reconstruction. This unbiased approach permits the comprehensive identification of all transcripts present in a sample, including annotated genes, novel isoforms of annotated genes, and novel genes. While common gene/transcript databases are quite large, they are not comprehensive, and the *de novo* transcriptome reconstruction approach ensures complete transcriptome(s) identification from the experimental samples. The leading tool for transcript reconstruction is `Stringtie`. Here, we will use `Stringtie` to predict transcript structures based on the reads aligned by `HISAT`.

> <hands-on-title>Transcriptome reconstruction</hands-on-title>
>
> 1. **Stringtie** {% icon tool %}: Run `Stringtie` on the `HISAT2` alignments using the default parameters.
>    - Use batch mode to run all four samples from one tool form.
>    - {% icon param-file %} *"Specify strand information"*: `Forward (FR)`
>
{: .hands_on}

# Transcriptome assembly

We just generated four transcriptomes with `Stringtie` representing each of the four RNA-seq libraries we are analyzing. Since these were generated in the absence of a reference transcriptome, and we ultimately would like to know what transcript structure corresponds to which annotated transcript (if any), we have to make a **transcriptome database**. We will use the tool `Stringtie - Merge` to combine redundant transcript structures across the four samples and the RefSeq reference. Once we have merged our transcript structures, we will use `GFFcompare` to annotate the transcripts of our newly created transcriptome so we know the relationship of each transcript to the RefSeq reference.

> <hands-on-title>Transcriptome assembly</hands-on-title>
>
> 1. **Stringtie-merge** {% icon tool %}: Run `Stringtie-merge` on the `Stringtie` assembled transcripts along with the RefSeq annotation file we imported earlier.
>    - {% icon param-file %} *"Transcripts"*: `all four `Stringtie` assemblies`
>    - {% icon param-file %} *"Reference annotation to include in the merging"*: `RefSeq_reference_GTF`
>
> 2. **GFFCompare** {% icon tool %}: Run `GFFCompare` on the `Stringtie-merge` generated transcriptome along with the RefSeq annotation file.
>    - {% icon param-file %} *"GTF inputs for comparison"*: `output of Stringtie-merge`
>    - *"Use Reference Annotation"*: `Yes`
>       - *"Choose the source for the reference annotation"*: `History`
>           - {% icon param-file %} *"Reference Annotation"*: `RefSeq_reference_GTF`
>    - *"Use Sequence Data"*: `Yes`
>       - *"Choose the source for the reference list"*: `Locally cached`
>           - *"Using reference genome"*: 'Mouse (Mus Musculus): mm10'
>
{: .hands_on}

> <comment-title>Note: Transcript categorization used by `GFFcompare`</comment-title>
>
> |**Class code** | **Transcript category**|
> |:---:|:---|
> |= | Annotated in reference|
> |j | Novel isoform of reference|
> |u | Intergenic|
> |x | Anti-sense|
> |r | Repetitive|
> |c | Contained in exon of reference|
> |s | Anti-sense spliced intronic|
> |e | Single exon transfrag overlapping a reference exon and at least 10 bp of a reference intron, indicating a possible pre-mRNA fragment.|
> |i | A transfrag falling entirely within a reference intron|
> |o | Generic exonic overlap with a reference transcript|
> |p | Possible polymerase run-on fragment (within 2Kbases of a reference transcript)|
{: .comment}


# Analysis of the differential gene expression

We just generated a transriptome database that represents the transcripts present in the G1E and megakaryocytes samples. This database provides the location of our transcripts with non-redundant identifiers, as well as information regarding the origin of the transcript.

We now want to identify which transcripts are differentially expressed between the G1E and megakaryocyte cellular states. To do this we will implement a counting approach using `FeatureCounts` to count reads per transcript. Then we will provide this information to `DESeq2` to generate normalized transcript counts (abundance estimates) and significance testing for differential expression.

## Count the number of reads per transcript

To compare the abundance of transcripts between different cellular states, the first essential step is to quantify the number of reads per transcript. [`FeatureCounts`](http://bioinf.wehi.edu.au/featureCounts/ ) is one of the most popular tools for counting reads in genomic features. In our case, we'll be using `FeatureCounts` to count reads aligning in exons of our `GFFCompare` generated transcriptome database.

The recommended mode is "union", which counts overlaps even if a read only shares parts of its sequence with a genomic feature and disregards reads that overlap more than one feature.

> <hands-on-title>Counting the number of reads per transcript</hands-on-title>
>
> 1. **FeatureCounts** {% icon tool %}: Run `FeatureCounts` on the aligned reads (`HISAT2` output) using the `GFFCompare` transcriptome database as the annotation file.
>
>    - {% icon param-file %} *"Alignment file"*: `four HISAT2 aligned read files`
>    - *"Specify strand information"*: `Stranded (Forward)`
>    - {% icon param-file %} *"Gene annotation file"*: `in your history`
>       - *"Gene annotation file"*: `the annotated transcripts GTF file output by GFFCompare`
>    - *"Options for paired-end reads"*
>       - *"Count fragments instead of reads"*: `Enabled; fragments (or templates) wil be counted instead of reads`
>    - *"Advanced options"*
>       - *"GFF gene identifier"*: `transcript_id`
>
{: .hands_on}

## Perform differential gene expression testing

Transcript expression is estimated from read counts, and attempts are made to correct for variability in measurements using replicates. This is absolutely essential to obtaining accurate results. We recommend having at least two biological replicates.

[`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) is a great tool for differential gene expression analysis. It accepts read counts produced by `FeatureCounts` and applies size factor normalization:

- Computation for each gene of the geometric mean of read counts across all samples
- Division of every gene count by the geometric mean
- Use of the median of these ratios as sample's size factor for normalization

> <hands-on-title></hands-on-title>
>
> 1. **DESeq2** {% icon tool %}: Run `DESeq2` with the following parameters:
>    - *"1: Factor"*
>       - *"1: Factor level"*: `G1E`
>          - {% icon param-file %} *"Counts file(s)"*: `featureCount files corresponding to the two G1E replicates`
>       - *"2: Factor level"*: `Mega`
>          - {% icon param-file %} *"Counts file(s)"*: `featureCount files corresponding to the two Mega replicates`
>
>       > <comment-title></comment-title>
>       >
>       > You can select several files by holding down the CTRL (or COMMAND) key and clicking on the desired files
>       {: .comment}
>    - *"Visualising the analysis results"*: Yes
>    - *"Output normalized counts table"*: Yes
>
{: .hands_on}

The first output of `DESeq2` is a tabular file. The columns are:

1.	Gene identifiers
2.	Mean normalized counts, averaged over all samples from both conditions
3.	Logarithm (base 2) of the fold change (the values correspond to up- or downregulation relative to the condition listed as Factor level 1)
4.	Standard error estimate for the log2 fold change estimate
5.	[Wald](http://data.princeton.edu/wws509/notes/c2s3.html) statistic
6.	*p*-value for the statistical significance of this change
7.	*p*-value adjusted for multiple testing with the Benjamini-Hochberg procedure which controls false discovery rate ([FDR](http://www.biostathandbook.com/multiplecomparisons.html))


> <hands-on-title></hands-on-title>
>
>1. **Filter** {% icon tool %}: Run `Filter` to extract genes with a significant change in gene expression (adjusted *p*-value less than 0.05) between treated and untreated samples
>
>    > <question-title></question-title>
>    >
>    > How many transcripts have a significant change in expression between these conditions?
>    >
>    > > <solution-title></solution-title>
>    > > To filter, use "c7<0.05". And we get 249 transcripts with a significant change in gene expression between the G1E and megakaryocyte cellular states.
>    > {: .solution }
>    {: .question}
>
> 2. **Filter** {% icon tool %}: Determine how many transcripts are up or down regulated in the G1E state.
>
>    > <comment-title></comment-title>
>    > Rename your datasets for the downstream analyses
>    {: .comment}
>
>    > <question-title></question-title>
>    >
>    > Are there more upregulated or downregulated genes in the treated samples?
>    >
>    > > <solution-title></solution-title>
>    > > To obtain the up-regulated genes in the G1E state, we filter the previously generated file (with the significant change in transcript expression) with the expression "c3>0" (the log2 fold changes must be greater than 0). We obtain 102  genes (40.9% of the genes with a significant change in gene expression). For the down-regulated genes in the G1E state, we did the inverse and we find 149 transcripts (59% of the genes with a significant change in transcript expression).
>    > {: .solution }
>    {: .question}
{: .hands_on}

In addition to the list of genes, `DESeq2` outputs a graphical summary of the results, useful to evaluate the quality of the experiment:

1. Histogram of *p*-values for all tests

    ![DESeq2 histogram]({% link topics/transcriptomics/images/transcriptomics_images/Deseq2_histogram3.png %})

2. [MA plot](https://en.wikipedia.org/wiki/MA_plot): global view of the relationship between the expression change of conditions (log ratios, M), the average expression strength of the genes (average mean, A), and the ability of the algorithm to detect differential gene expression. The genes that passed the significance threshold (adjusted p-value < 0.1) are colored in red.

    ![MA plot]({% link topics/transcriptomics/images/transcriptomics_images/Deseq2_MAplot3.png %})

3. Principal Component Analysis ([PCA](https://en.wikipedia.org/wiki/Principal_component_analysis)) and the first two axes

    ![Principal component analysis output]({% link topics/transcriptomics/images/transcriptomics_images/Deseq2_PCA3.png %})

    Each replicate is plotted as an individual data point. This type of plot is useful for visualizing the overall effect of experimental covariates and batch effects.

4. Heatmap of sample-to-sample distance matrix: overview over similarities and dissimilarities between samples

    ![Heatmap of sample-to-sample distance matrix]({% link topics/transcriptomics/images/transcriptomics_images/Deseq2_heatmap3.png %})

5. Dispersion estimates: gene-wise estimates (black), the fitted values (red), and the final maximum a posteriori estimates used in testing (blue)

    ![Dispersion estimates]({% link topics/transcriptomics/images/transcriptomics_images/Deseq2_dispersion3.png %})

    This dispersion plot is typical, with the final estimates shrunk from the gene-wise estimates towards the fitted estimates. Some gene-wise estimates are flagged as outliers and not shrunk towards the fitted value. The amount of shrinkage can be more or less than seen here, depending on the sample size, the number of coefficients, the row mean and the variability of the gene-wise estimates.


For more information about `DESeq2` and its outputs, you can have a look at [`DESeq2` documentation](https://www.bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf).

# Visualization
Now that we have a list of transcript expression levels and their differential expression levels, it is time to visually inspect our transcript structures and the reads they were predicted from. It is a good practice to visually inspect (and present) loci with transcripts of interest. Fortunately, there is a built-in genome browser in Galaxy, **Trackster**, that make this task simple (and even fun!).

In this last section, we will convert our aligned read data from BAM format to bigWig format to simplify observing where our stranded RNA-seq data aligned to. We'll then initiate a session on Trackster, load it with our data, and visually inspect our interesting loci.

> <hands-on-title>Converting aligned read files to bigWig format</hands-on-title>
>
> 1. **bamCoverage** {% icon tool %}: Run `bamCoverage` on all four aligned read files (`HISAT2` output) with the following parameters:
>    - *"Bin size in bases"*: '1'
>    - *"Effective genome size"*: 'mm9 (2150570000)''
>    - *"Advanced options"*
>       - *"Only include reads originating from fragments from the forward or reverse strand"*: 'forward'
>
> 2. **Rename** {% icon tool %}: Rename the outputs to reflect the origin of the reads and that they represent the reads mapping to the PLUS strand.
>
> 3. **bamCoverage** {% icon tool %}: Repeat Step 1 except changing the following parameter:
>    - *"Only include reads originating from fragments from the forward or reverse strand"*: 'reverse'
>
> 4. **Rename** {% icon tool %}: Rename the outputs to reflect the origin of the reads and that they represent the reads mapping to the MINUS strand.
{: .hands_on}

> <hands-on-title>Trackster based visualization</hands-on-title>
>
> 1. **Viz** {% icon tool %}: On the center console at the top of the Galaxy interface, choose " Visualization" -> "New track browser"
>    - Name your visualization someting descriptive under "Browser name:"
>    - Choose "Mouse Dec. 2011 (GRCm38/mm10) (mm10)" as the "Reference genome build (dbkey)
>    - Click "Create" to initiate your Trackster session
> ![Trackster opening window](../../images/transcriptomics_images/Trackster_opening_window.png)
>
> 2. **Viz** {% icon tool %}: Click "Add datasets to visualization"
>    - Select the "RefSeq GTF mm10" file
>    - Select the output files from `Stringtie`
>    - Select the output file from `GFFCompare`
>    - Select the output files from `bamCoverage`
>
> 3. {% icon tool %}: Using the grey labels on the left side of each track, drag and arrange the track order to your preference.
>
> 4. {% icon tool %}: Hover over the grey label on the left side of the "RefSeq GTF mm10" track and click the "Edit settings" icon.
>    - Adjust the block color to blue (#0000ff) and antisense strand color to red (#ff0000)
>
> 5. {% icon tool %}: Repeat the previous step on the output files from `StringTie` and `GFFCompare`.
>
> 6. {% icon tool %}: Hover over the grey label on the left side of the "G1E R1 plus" track and click the "Edit settings" icon.
>    - Adjust the color to blue (#0000ff)
>
> 7. {% icon tool %}: Repeat the previous step on the other three bigWig files representing the plus strand.
>
> 8. {% icon tool %}: Hover over the grey label on the left side of the "G1E R1 minus" track and click the "Edit settings" icon.
>    - Adjust the color to red (#ff0000)
>
> 9. {% icon tool %}: Repeat the previous step on the other three bigWig files representing the minus strand.
>
> 10. {% icon tool %}: Adjust the track height of the bigWig files to be consistent for each set of plus strand and minus strand tracks.
> ![Adjusting the track height](../../images/transcriptomics_images/Hoxb13_locus_screenshot.png)
> 11. {% icon tool %}: Direct Trackster to the coordinates: chr11:96191452-96206029, what do you see?
>
>    > <question-title></question-title>
>    > what do you see?
>    >
>    > > <solution-title></solution-title>
>    > > 1. There are two clusters of transcripts that are exclusively expressed in the G1E background
>    > > 2. The left-most transcript is the Hoxb13 transcript
>    > > 3. The center cluster of transcripts are not present in the RefSeq annotation and are determined by `GFFCompare` to be "u" and "x"
>    > {: .solution }
>    {: .question}
>
{: .hands_on}

# Conclusion


In this tutorial, we have analyzed RNA sequencing data to extract useful information, such as which genes are expressed in the G1E and megakaryocyte cellular states and which of these genes are differentially expressed between the two cellular states. In addition, we identified unannotated genes that are expressed in a cell-state dependent manner and at a locus with relevance to differentiation and development. To identify these transcripts, we analyzed RNA sequence datasets using a de novo transcriptome reconstruction RNA-seq data analysis approach. This approach can be summed up with the following scheme:

![RNAseq de novo tutorial workflow](../../images/transcriptomics_images/schematic_for_RNAseq_de_novo_tutorial.png)
