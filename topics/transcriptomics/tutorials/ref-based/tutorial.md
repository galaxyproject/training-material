---
layout: tutorial_hands_on
topic_name: transcriptomics
tutorial_name: ref-based
---

# Introduction
{:.no_toc}

In the study of [Brooks *et al.* 2011](http://genome.cshlp.org/content/21/2/193.long), the Pasilla (PS) gene, *Drosophila* homologue of the Human splicing regulators Nova-1 and Nova-2 Proteins, was depleted in *Drosophila melanogaster* by RNAi. The authors wanted to identify exons that are regulated by Pasilla gene using RNA sequencing data.

Total RNA was isolated and used for preparing either single-end or paired-end RNA-seq libraries for treated (PS depleted) samples and untreated samples. These libraries were sequenced to obtain a collection of RNA sequencing reads for each sample. The effects of Pasilla gene depletion on splicing events can then be analyzed by comparison of RNA sequencing data of the treated (PS depleted) and the untreated samples.

The genome of *Drosophila melanogaster* is known and assembled. It can be used as reference genome to ease this analysis.  In a reference based RNA-seq data analysis, the reads are aligned (or mapped) against a reference genome, *Drosophila melanogaster* here, to significantly improve the ability to reconstruct transcripts and then identify differences of expression between several conditions.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Pretreatments

## Data upload

The original data is available at NCBI Gene Expression Omnibus (GEO) under accession number [GSE18508](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18508).

We will look at the 7 first samples:

- 3 treated samples with Pasilla (PS) gene depletion: [GSM461179](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461179), [GSM461180](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461180), [GSM461181](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4611810)
- 4 untreated samples: [GSM461176](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461176), [GSM461177](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461177), [GSM461178](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461178), [GSM461182](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461182)

Each sample constitutes a separate biological replicate of the corresponding condition (treated or untreated). Moreover, two of the treated and two of the untreated samples are from a paired-end sequencing assay, while the remaining samples are from a single-end sequencing experiment.

We have extracted sequences from the Sequence Read Archive (SRA) files to build FASTQ files.

> ### :pencil2: Hands-on: Data upload
>
> 1. Create a new history for this RNA-seq exercise
> 2. Import a FASTQ file pair (*e.g.*  `GSM461177_untreat_paired_chr4_R1.fastq` and `GSM461177_untreat_paired_chr4_R2.fastq`) from [Zenodo](https://dx.doi.org/10.5281/zenodo.290221)
>
>    > ### :bulb: Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * Press **Start**    
>    {: .tip}
>
>    > ### :bulb: Tip: Changing the file type `fastq` to `fastqsanger` once the data file is in your history
>    >
>    > * Click on the pencil button displayed in your dataset in the history
>    > * Choose **Datatype** on the top
>    > * Select `fastqsanger`
>    > * Press **Save**
>    {: .tip}
>
>    As default, Galaxy takes the link as name. It also do not link the dataset to a database or a reference genome.
>
>    > ### :nut_and_bolt: Comments
>    > - Edit the "Database/Build" to select "dm3"
>    > - Rename the datasets according to the samples
>    {: .comment}
>
{: .hands_on}

Both files contain the reads that belong to chromosome 4 of a paired-end sample. The sequences are raw sequences from the sequencing machine, without any pretreatments. They need to be controlled for their quality.

## Quality control

For quality control, we use similar tools as described in [NGS-QC tutorial]({{site.url}}topics/sequence-analysis): [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).

> ### :pencil2: Hands-on: Quality control
>
> 1. **FastQC** :wrench:: Run FastQC on both FastQ files to control the quality of the reads
>
>    > ### :question: Questions
>    >
>    > 1. What is the read length?
>    > 2. Is there anything what you find striking when you compare both reports?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>The read length is 37 bp</li>
>    >    <li>Both reports for GSM461177_untreat_paired_chr4_R1 and for GSM461177_untreat_paired_chr4_R2 are quite ok. For GSM461177_untreat_paired_chr4_R1, there is several warnings and an issue on the Kmer content. For GSM461177_untreat_paired_chr4_R2, the quality in the 2nd tile is bad (maybe because of some event during sequencing). We need to be careful for the quality treatment and to do it with paired-end information</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 2. **Trim Galore** :wrench:: Treat for the quality of sequences by running Trim Galore on the paired-end datasets
>
>    > ### :question: Questions
>    >
>    > Why is Trim Galore run once on the paired-end dataset and not twice on each dataset?
>    >
>    > <details>
>    > <summary>Click to view answers</summary>
>    > Trim Galore can remove sequences if they become too short during the trimming process. For paired-end files Trim Galore! removes entire sequence pairs if one (or both) of the two reads became shorter than the set length cutoff. Reads of a read-pair that are longer than a given threshold but for which the partner read has become too short can optionally be written out to single-end files. This ensures that the information of a read pair is not lost entirely if only one read is of good quality.
>    > </details>
>    {: .question}
>
> 3. **FastQC** :wrench:: Re-run FastQC on Trim Galore's outputs and inspect the differences
>
>    > ### :question: Questions
>    >
>    > 1. How are the changes in the read length?
>    > 2. Is there any characteristics impacted by Trim Galore?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>The read length is then now from 20 to 37 bp</li>
>    >    <li>For GSM461177_untreat_paired_chr4_R1, the per base sequence content is now red. For GSM461177_untreat_paired_chr4_R2, the per tile sequence quality is still bad but now also the per base sequence content and the Kmer Content</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
{: .hands_on}

As the genome of *Drosophila melanogaster* is known and assembled, we can use this information and map the sequences on this genome to identify the effects of Pasilla gene depletion on splicing events.

# Mapping

To make sense of the reads, their positions within *Drosophila melanogaster* genome must be determined. This process is known as aligning or 'mapping' the reads to the reference genome.

> ### :nut_and_bolt: Comment
>
> Do you want to learn more about the principles behind mapping? Follow our [training]({{site.url}}topics/sequence-analysis/)
> {: .comment}

Because in the case of a eukaryotic transcriptome, most reads originate from processed mRNAs lacking exons, they cannot be simply mapped back to the genome as we normally do for DNA data. Instead the reads must be separated into two categories:

- Reads that map entirely within exons
- Reads that cannot be mapped within an exon across their entire length because they span two or more exons

Spliced mappers have been developed to efficiently map transcript-derived reads against genomes. [TopHat](https://ccb.jhu.edu/software/tophat/index.shtml) was one of the first tools designed specifically to address this problem:

1. Identification of potential exons using reads that do map to the genome
2. Generation of possible splices between neighboring exons
3. Comparison of reads that did not initially map to the genome against these *in silico* created junctions

![Kim et al.](../../images/tophat2.png "Kim et al., Genome Biology, 2013")

Here, we will use HISAT2, a successor to TopHat2 that is faster with low memory requirements.

To run efficiently the mapping, HISAT2 needs to know on important parameters about the sequencing library: the library type.

This information should usually come with your FASTQ files, ask your sequencing facility! If not, try to find them on the site where you downloaded the data or in the corresponding publication. Another option is to estimate these parameters with a *preliminary mapping* of a *downsampled* file and some analysis programs. Afterward, the actual mapping can be redone on the original files with the optimized parameters.

## Preliminary mapping

In a preliminary mapping, we will estimate the library type to run HISAT2 efficiently afterwards.

> ### :nut_and_bolt: Comment
> This step is not necessary if you don't need to estimate the library type of your data.
{: .comment}

The library type corresponds to a protocol used to generate the data: which strand the RNA fragment is synthesized from.

![Credit: Zhao Zhang](../../images/strand_library_type.png "Credit: Zhao Zhang(http://onetipperday.sterding.com/2012/07/how-to-tell-which-library-type-to-use.html)")

In the previous illustration, you could see that for example dUTP method is to only sequence the strand from the first strand synthesis (the original RNA strand is degradated due to the dUTP incorporated).

Examples of protocol | Description | Library Type (HISAT2)
--- | --- | ---
Standard Illumina | Reads from the left-most end of the fragment (in transcript coordinates) map to the transcript strand, and the right-most end maps to the opposite strand | Unstranded (default)
dUTP, NSR, NNSR | Same as above except we enforce the rule that the right-most end of the fragment (in transcript coordinates) is the first sequenced (or only sequenced for single-end reads). Equivalently, it is assumed that only the strand generated during first strand synthesis is sequenced. | First strand (FR/F)
Ligation, Standard SOLiD | Same as above except we enforce the rule that the left-most end of the fragment (in transcript coordinates) is the first sequenced (or only sequenced for single-end reads). Equivalently, it is assumed that only the strand generated during second strand synthesis is sequenced. | Second strand (RF/R)

If you do not know the library type, you can find it by yourself by mapping the reads on the reference genome and infer the library type from the mapping results by comparing reads mapping information to the annotation of the reference genome.

![Type of library, depending also of the type of sequencing](../../images/library_type_mapping.png "Type of library, depending also of the type of sequencing")

The sequencer always read from 5' to 3'. So, in First Strand case, all reads from the left-most end of RNA fragment (always from 5' to 3') are mapped to transcript-strand, and (for pair-end sequencing) reads from the right-most end are always mapped to the opposite strand.

We can now try to determine the library type of our data.

> ### :pencil2: Hands-on: Determining the library type
>
> 1. Load the Ensembl gene annotation for *Drosophila melanogaster* ([`Drosophila_melanogaster.BDGP5.78.gtf`](https://zenodo.org/record/290221/files/Drosophila_melanogaster.BDGP5.78.gtf)) from [Zenodo](https://dx.doi.org/10.5281/zenodo.290221) into your current Galaxy history and rename it
>
> 2. **Select first** :wrench:: Downsample the both FastQ files generated by Trim Galore to 200k or 1M reads
>
>    > ### :question: Questions
>    >
>    > How many rows must be selected to conserve 200k reads?
>    >
>    > <details>
>    > <summary>Click to view answers</summary>
>    > In a FastQ file, a read corresponds to 4 lines. So to conserve 200,000 reads, 800,000 must be selected.
>    > </details>
>    {: .question}
>
> 3. **HISAT2** :wrench:: Run **HISAT2** with:
>    - "FASTQ" as "Input data format"
>    - "Individual paired reads"
>    - Downsampled "Trimmed reads pair 1" (Trim Galore output) as "Forward reads"
>    - Downsampled "Trimmed reads pair 2" (Trim Galore output) as "Reverse reads"
>    - "dm3" as reference genome
>    - Default values for other parameters
>
> 4. **Infer Experiment** :wrench:: Run **Infer Experiment** to determine the library type:
>    - HISAT2 output as "Input BAM/SAM file"
>    - Drosophila annotation as "Reference gene model"
>
> 5. Check the results and search the tool's documentation for help on the meaning
>    
>    > ### :nut_and_bolt: Comment
>    > As it is sometimes quite difficult to find out which settings correspond to those of other programs, the following table might be helpful to identify the library type:
>    >
>    > Sequencing type | **Infer Experiment** | **TopHat** | **HISAT2** | **htseq-count** | **featureCounts**
>    > --- | --- | --- | --- | --- | ---
>    > Paired-End (PE) | "1++,1--,2+-,2-+" | "FR Second Strand" | "Second Strand F/FR" | "yes" | "1"
>    > PE | "1+-,1-+,2++,2--" | "FR First Strand" | "First Strand R/RF" | "reverse" | "2"
>    > Single-End (SE) | "++,--" | "FR Second Strand" | "Second Strand F/FR" | "yes" | "1"
>    > SE | "+-,-+" | "FR First Strand" | "First Strand R/RF" | "reverse" | "2"
>    > SE,PE | undecided | "FR Unstranded" | default | "no" | "0"
>    >
>    {: .comment}
>    
>    > ### :question: Question
>    >
>    > 1. Which fraction of the different type of library type?
>    > 2. Which library type do you chose? What is the correspondance in **HISAT2**?
>    >
>    >    <details>
>    >    <summary>Click to view answer</summary>
>    >    <ol type="1">
>    >    <li>Fraction of reads explained by "1++,1--,2+-,2-+": 0.0151 and Fraction of reads explained by "1+-,1-+,2++,2--": 0.9843</li>
>    >    <li>The library type seems then to be First Strand.</li>
>    >    </ol>
>    >    </details>
>    {: .question}
{: .hands_on}

## Actual mapping

We can now map all the RNA sequences on the *Drosophila melanogaster* genome using HISAT2.

> ### :pencil2: Hands-on: Spliced mapping
>
> 1. **HISAT2** :wrench:: Run **HISAT2** with:
>    - "FASTQ" as "Input data format"
>    - "Individual paired reads"
>    - "Trimmed reads pair 1" (Trim Galore output) as "Forward reads"
>    - "Trimmed reads pair 2" (Trim Galore output) as "Reverse reads"
>    - "dm3" as reference genome
>    - Default values for other parameters except "Spliced alignment parameters"
>    - "Specify strand-specific information" to the previously determined value
>    - `Drosophila_melanogaster.BDGP5.78.gtf` as "GTF file with known splice sites"
>
> 2. Inspect the mapping statistics
>    - Click on "View details" (Warning icon)
>    - Click on "stderr" (Tool Standard Error)
>
>    > ### :question: Question
>    >
>    > 1. How many paired reads were mapped 1 time? And how many paired reads were mapped more than 1 time?
>    > 2. How many reads were mapped but without their mate pair?
>    > 3. What is the overall alignment rate?
>    >
>    >    <details>
>    >    <summary>Click to view answer</summary>
>    >    <ol type="1">
>    >    <li>29.06% and 6.10%</li>
>    >    <li>11,963 reads (the reads aligned discordantly 1 time)</li>
>    >    <li>The overall alignment rate is 75.77%. It counts proportion of mapped reads: (12017+2486+11963+7390/2+1413/2)/40736</li>
>    >    </ol>
>    >    </details>
>    {: .question}
{: .hands_on}

**HISAT2** generates a BAM file with the mapped reads.

> ### :question: Question
>
> 1. What is a BAM file?
> 2. What does such file contain?
>
>    <details>
>    <summary>Click to view answer</summary>
>    <ol type="1">
>    <li>BAM file is the binary version of SAM file</li>
>    <li>It contains information about the mapping: for each mapped read, the position on the reference genome, the mapping quality, ...</li>
>    </ol>
>    </details>
{: .question}

The mapping exercise worked for you? Great! :tada:

> ### :pencil2: (Optional) Hands-on: Map other datasets
>
> You can do the same process on the other sequence files available on [Zenodo](https://dx.doi.org/10.5281/zenodo.290221)
>
> - Paired-end data
>     - `GSM461178_untreat_paired_chr4_R1` and `GSM461178_untreat_paired_chr4_R2`
>     - `GSM461180_treat_paired_chr4_R1` and `GSM461180_treat_paired_chr4_R2`
>     - `GSM461181_treat_paired_chr4_R1` and `GSM461181_treat_paired_chr4_R2`
> - Single-end data
>     - `GSM461176_untreat_single_chr4`
>     - `GSM461179_treat_single_chr4`
>     - `GSM461182_untreat_single_chr4`
>
> This is really interesting to redo on the other datasets, specially to check how the parameters are inferred given the different type of data.
{: .hands_on}

## Inspection of HISAT2 results

The BAM file contains information about where the reads are mapped on the reference genome. But it is binary file and with the information for more than 3 millions of reads, it makes it difficult to visualize it, to inspect how the reads are mapped on the reference genome.

There is

> ### :pencil2: Hands-on: Inspection of HISAT2 results
>
> 1. **IGV** :wrench:: Visualize the HISAT2 output BAM file, particularly the region on chromosome 4 between 560 kb to 600 kb (`chr4:560,000-600,000`)
>
>    > ### :nut_and_bolt: Comment
>    >
>    > Check [IGV documentation](https://software.broadinstitute.org/software/igv/AlignmentData)
>    {: .comment}
>
>    > ### :question: Question
>    >
>    > 1. Which information does appear on the top in grey?
>    > 2. What does the horizontal bar in the following screenshot represent?
>    >
>    >    ![Screenshot of IGV on Chromosome 4](../../images/junction_igv_screenshot.png "Screenshot of IGV on Chromosome 4")
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>The coverage plot: the sum of mapped reads at each position</li>
>    >    <li>They represent the junctions events (or splice sites): when reads are mapped over an intron</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 3. **IGV** :wrench:: Unzoom to `chr4:560,000-600,000` and inspect the splice junctions using a **Sashimi plot**
>
>    > ### :bulb: Tip: Creation of a Sashimi plot
>    >
>    > * Right click on the BAM file
>    > * Select **Sashimi Plot** from the context menu
>    {: .tip}    
>
>    > ### :question: Question
>    >
>    > ![Screenshot of a Sashimi plot of Chromosome 4](../../images/hisat_igv_sashimi.png "Screenshot of a Sashimi plot of Chromosome 4")
>    >
>    > 1. What does the vertical bar graph represent? And the numbered arcs?
>    > 2. What does the number means?
>    > 3. Why do we observe 4 different blue group of linked boxes on the bottom?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>The coverage for each alignment track is plotted as a bar graph. Arcs representing splice junctions connecting exons</li>
>    >    <li>Arcs display the number of reads split across the junction (junction depth). </li>
>    >    <li>The 4 group of linked boxes on the bottom represent the 4 possible RNAs, changing because of different splicing events of the first intron.</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
>    > ### :nut_and_bolt: Comment
>    >
>    > Check [IGV documentation on Sashimi plots](https://software.broadinstitute.org/software/igv/Sashimi) to find some clues
>    {: .comment}
>
{: .hands_on}

After the mapping, we have in the generated mapping file the information about where the reads are mapped on the reference genome. So for each mapped read, we know where it is mapped and how good it was mapped.

The next step in the RNA-Seq data analysis is quantification of expression level of the genomic features (gene, transcript, exons, ...) to be able then to compare several samples for the different expression analysis. The quantification consist into taking each known genomic feature (*e.g.* gene) of the reference genome and then counting how many reads are mapped on this genomic feature. So, in this step, we start with an information per mapped reads to end with an information per genomic feature.

> ### :nut_and_bolt: Comment
>
> The quantification depends on the definition of the genomic features of the reference genome, and then on the annotations. We strongly recommend you to use an annotation corresponding to the same version of the reference genome you used for the mapping.
{: .comment}

To identify exons that are regulated by the Pasilla gene, we need to identify genes and exons which are differentially expressed between samples with PS gene depletion and control samples.
In this tutorial, we will then analyze the differential gene expression, but also the differential exon usage.

# Analysis of the differential gene expression

We will first investigate the differential gene expression to identify which genes are impacted by the Pasilla gene depletion

## Count the number of reads per annotated gene

To compare the expression of single genes between different conditions (*e.g.* with or without PS depletion), an first essential step is to quantify the number of reads per gene. [**HTSeq-count**](https://www-huber.embl.de/users/anders/HTSeq/doc/count.html) is one of the most popular tool for gene quantification.

To quantify the number of reads mapped to a gene, an annotation of the gene position is needed. We already upload on Galaxy the `Drosophila_melanogaster.BDGP5.78.gtf` with the Ensembl gene annotation for *Drosophila melanogaster*.

In principle, the counting of reads overlapping with genomic features is a fairly simple task. But there are some details that need to be decided, such how to handle multi-mapping reads. **HTSeq-count** offers 3 choices ("union", "intersection_strict" and "intersection_nonempty") to handle read mapping to multiple locations, reads overlapping introns, or reads that overlap more than one genomic feature:

![HT-Seq methods to handle overlapping reads](../../images/htseq_count.png "HT-Seq methods to handle overlapping reads (HTSeq documentation: https://www-huber.embl.de/users/anders/HTSeq/doc/count.html)")

The recommended mode is "union", which counts overlaps even if a read only shares parts of its sequence with a genomic feature and disregards reads that overlap more than one feature.

> ### :pencil2: Hands-on: Counting the number of reads per annotated gene
>
> 1. **HTSeq-count** :wrench:: Run **HTSeq-count** on the BAM file with
>    - `Drosophila_melanogaster.BDGP5.78.gtf` as "GFF file"
>    - The "Union" mode
>    - A "Minimum alignment quality" of 10
>    - Appropriate value for "Stranded" option
>       
> 2. Inspect the result files
>
>    > ### :question: Question
>    >
>    > 1. Which information does the result file contain?
>    > 2. Which feature has the most reads mapped on it?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>The useful result file is a tabular file with two columns: the gene id and the number of reads mapped on the corresponding gene</li>
>    >    <li>To display the most found feature, we need first to sort the output file with the feature and the number of reads found for these feature. We do that using Sort tool, on the second column and in descending order. And we found that FBgn0017545 is the feature with the most reads mapped on it with 4,030 reads.</li>
>    >    </ol>
>    >    </details>
>    {: .question}
{: .hands_on}

## Analysis of the differential gene expression

In the previous section, we counted only reads that mapped to genes of chromosome 4 and for only one sample. To be able to identify differential gene expression induced by PS depletion, all datasets (3 treated and 4 untreated) must be analyzed with the similar procedure and for whole genome.

For time and computer saving, we run the previous steps for you and obtain 7 count files, available on [Zenodo](https://dx.doi.org/10.5281/zenodo.290221).

These files contain for each gene of Drosophila the number of reads mapped to it. We could compare directly the files and then having the differential gene expression. But the number of sequenced reads mapped to a gene depends on:

- Its own expression level
- Its length
- The sequencing depth of the sample
- The expression of all other genes within the sample

Either for within or for inter-sample comparison, the gene counts need to be normalized. We can then use the Differential Gene Expression (DGE) analysis, whose two basic tasks are:

- Estimate the biological variance using the replicates for each condition
- Estimate the significance of expression differences between any two conditions

This expression analysis is estimated from read counts and attempts are made to correct for variability in measurements using replicates that are absolutely essential accurate results. For your own analysis, we advice you to use at least 3, better 5 biological replicates per condition. You can have different number of replicates per condition.

[**DESeq2**](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) is a great tool for DGE analysis. It takes read counts produced by **HTseq-count**, combine them into a big table (with gene in the rows and samples in the columns) and applies size factor normalization:

- Computation for each gene of the geometric mean of read counts across all samples
- Division of every gene count by the geometric mean
- Use of the median of these ratios as sample's size factor for normalization

Multiple factors with several levels can then be incorporated in the analysis. With the normalization, we can then compare with statistical liability the gene expression of different levels of a factor.

In our example, we have samples with two varying factors that can explain differences in gene expression:

- Treatment (either treated or untreated)
- Sequencing type (paired-end or single-end)

Here treatment is the primary factor which we are interested in. The sequencing type is some further information that we know about the data that might effect the analysis. This particular multi-factor analysis allows us to assess the effect of the treatment taking also the sequencing type into account.

> ### :nut_and_bolt: Comment
>
> We recommend you to add as many factors as you think that may affect the gene expression. It can be the sequencing type like here, but it can also be the manipulation (if different persons are involved in the library preparation), ...
{: .comment}

> ### :pencil2: Hands-on: Analysis of the differential gene expression (1)
>
> 1. Create a new history
> 2. Import the seven count files from [Zenodo](https://dx.doi.org/10.5281/zenodo.290221)
>    - `GSM461176_untreat_single.deseq.counts`
>    - `GSM461177_untreat_paired.deseq.counts`
>    - `GSM461178_untreat_paired.deseq.counts`
>    - `GSM461179_treat_single.deseq.counts`
>    - `GSM461180_treat_paired.deseq.counts`
>    - `GSM461181_treat_paired.deseq.counts`
>    - `GSM461182_untreat_single.deseq.counts`
>
> 3. **DESeq2** :wrench:: Run **DESeq2** with:
>    - "Treatment" as first factor with "treated" and "untreated" as levels and selection of count files corresponding to both levels
>
>       > ### :bulb: Tip
>       >
>       > You can select several files by keeping the CTRL (or COMMAND) key pressed and clicking on the interesting files
>       {: .tip}
>
>    - "Sequencing" as second factor with "PE" and "SE" as levels and selection of count files corresponding to both levels
>
>    > ### :nut_and_bolt: Comment
>    >
>    > File names have all information needed
>    {: .comment}
{: .hands_on}

The first output of **DESeq2** is a tabular file. The columns are:

1.	Gene identifiers
2.	Mean normalized counts, averaged over all samples from both conditions
3.	Logarithm (to basis 2) of the fold change


    The log2 fold changes are based on primary factor level 1 vs. factor level 2. The order of factor levels is then important. For example, for the factor 'Treatment', DESeq2 computes fold changes of 'treated' samples against 'untreated', *i.e.* the values correspond to up- or downregulations of genes in treated samples.

4.	Standard error estimate for the log2 fold change estimate
5.	[Wald](https://en.wikipedia.org/wiki/Wald_test) statistic
6.	*p*-value for the statistical significance of this change
7.	*p*-value adjusted for multiple testing with the Benjamini-Hochberg procedure which controls false discovery rate ([FDR](https://en.wikipedia.org/wiki/False_discovery_rate))

> ### :pencil2: Hands-on: Analysis of the differential gene expression (2)
>
> 1. **Filter** :wrench:: Run **Filter** to extract genes with a significant change in gene expression (adjusted *p*-value equal or below 0.05) between treated and untreated samples
>
>    > ### :question: Question
>    >
>    > How many genes have a significant change in gene expression between these conditions?
>    >
>    > <details>
>    > <summary>Click to view answers</summary>
>    > To filter, you need to add the expression "c7&lt;0.05". And we get 751 genes (5.05%) with a significant change in gene expression between treated and untreated samples.
>    > </details>
>    {: .question}
>
>    > ### :nut_and_bolt: Comment
>    >
>    > The file with the independent filtered results can be used for further downstream analysis as it excludes genes with only few read counts as these genes will not be considered as significantly differentially expressed.
>    {: .comment}
>
> 2. **Filter** :wrench:: Extract genes that are significantly up and downregulated in treated samples
>
>    > ### :nut_and_bolt: Comments
>    > Rename your datasets for the downstream analyses
>    {: .comment}
>
>    > ### :question: Question
>    >
>    > Are there more upregulated or downregulated genes in the treated samples?
>    >
>    > <details>
>    > <summary>Click to view answers</summary>
>    > To obtain the up-regulated genes, we filter the previously generated file (with the significant change in gene expression) with the expression "c3>0" (the log2 fold changes must be greater than 0). We obtain 331 genes (44.07% of the genes with a significant change in gene expression). For the down-regulated genes, we did the inverse and we 420 genes (55.93% of the genes with a significant change in gene expression)
>    > </details>
>    {: .question}
{: .hands_on}

In addition to the list of genes, **DESeq2** outputs a graphical summary of the results, useful to evaluate the quality of the experiment:

1. Histogram of *p*-values for all tests
2. [MA plot](https://en.wikipedia.org/wiki/MA_plot): global view of the relationship between the expression change of conditions (log ratios, M), the average expression strength of the genes (average mean, A), and the ability of the algorithm to detect differential gene expression. The genes that passed the significance threshold (adjusted p-value < 0.1) are colored in red.
3. Principal Component Analysis ([PCA](https://en.wikipedia.org/wiki/Principal_component_analysis)) and the first two axes

    Each replicate is plotted as an individual data point. This type of plot is useful for visualizing the overall effect of experimental covariates and batch effects.

    > ### :question: Questions
    >
    > 1. What is the first axis separating?
    > 2. And the second axis?    
    >
    >    <details>
    >    <summary>Click to view answers</summary>
    >    <ol type="1">
    >    <li>The first axis is seperating the treated samples from the untreated samples, as defined when DESeq2 was launched</li>
    >    <li>The second axis is separating the single-end datasets from the paired-end datasets</li>
    >    </ol>
    >    </details>
    {: .question}


4. Heatmap of sample-to-sample distance matrix: overview over similarities and dissimilarities between samples

    > ### :question: Questions
    >
    > How are the samples grouped?
    >
    >    <details>
    >    <summary>Click to view answers</summary>
    >    <ol type="1">    
    >    <li>They are first grouped depending on the treatment (the first factor) and after on the library type (the second factor), as defined when DESeq2 was launched</li>
    >    </ol>
    >    </details>
    {: .question}

5. Dispersion estimates: gene-wise estimates (black), the fitted values (red), and the final maximum a posteriori estimates used in testing (blue)

    This dispersion plot is typical, with the final estimates shrunk from the gene-wise estimates towards the fitted estimates. Some gene-wise estimates are flagged as outliers and not shrunk towards the fitted value. The amount of shrinkage can be more or less than seen here, depending on the sample size, the number of coefficients, the row mean and the variability of the gene-wise estimates.


For more information about **DESeq2** and its outputs, you can have a look at [**DESeq2** documentation](https://www.bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf).

## Analysis of the functional enrichment among differentially expressed genes

We have extracted genes that are differentially expressed in treated (with PS gene depletion) samples compared to untreated samples. We would like to know the functional enrichment among the differentially expressed genes.

The Database for Annotation, Visualization and Integrated Discovery ([DAVID](https://david.ncifcrf.gov/)) provides a comprehensive set of functional annotation tools for investigators to understand the biological meaning behind large lists of genes.

The query to DAVID can be done only on 100 genes. So, we will need to select the ones where the most interested in.

> ### :pencil2: Hands-on:
>
> 1. **Sort** :wrench:: Sort the 2 datasets generated previously (upregulated genes and downregulated genes) given the log2 fold change, in descending or ascending order (to obtain the higher absolute log2 fold changes on the top)
> 1. **Select first lines from a dataset** :wrench:: Extract the first 100 lines of sorted files
> 2. **DAVID** :wrench:: Run **DAVID** on these files with
>     - First column as "Column with identifiers"
>     - "ENSEMBL_GENE_ID" as "Identifier type"
>
>    The output of the **DAVID** tool is a HTML file with a link to the DAVID website.
>
> 2. Inspect the Functional Annotation Chart
>
>    > ### :question: Questions
>    >
>    > What functional categories are the most represented?
>    >  
>    > <details>
>    > <summary>Click to view answers</summary>
>    > The up-regulated genes are mostly related to membrane (in the number of genes). The most represented functional categories are linked to signal and pathways for the down-regulated genes.
>    > </details>
>    {: .question}
>
> 3. Inspect the Functional Annotation Clusterings
>
>    > ### :question: Questions
>    >
>    > What functional annotations are the first clusters related to?
>    >  
>    > <details>
>    > <summary>Click to view answers</summary>
>    > For the up-regulated genes, the first cluster is more composed of functions related to chaperone and stress response. The down-regulated genes are more linked to ligase activity.
>    > </details>
>    {: .question}
{: .hands_on}

# Inference of the differential exon usage

Now, we would like to know the differential exon usage between treated (PS depleted) and untreated samples using RNA-seq exon counts. We will rework on the mapping results we generated previously.

We will use [DEXSeq](https://www.bioconductor.org/packages/release/bioc/html/DEXSeq.html). DEXSeq detects high sensitivity genes, and in many cases exons, that are subject to differential exon usage. But first, as for the differential gene expression, we need to count the number of reads mapping the exons.

## Count the number of reads per exon

This step is similar to the step of [counting the number of reads per annotated gene](#count-the-number-of-reads-per-annotated-gene). Here instead of HTSeq-count, we are using DEXSeq-Count.

> ### :pencil2: Hands-on: Counting the number of reads per exon
>
> 1. **DEXSeq-Count** :wrench:: Use the **DEXSeq-Count** to prepare the *Drosophila* annotations (`Drosophila_melanogaster.BDGP5.78.gtf`) to extract only exons with corresponding gene ids
>     - "Prepare annotation" of "Mode of operation"
>
>    The output is again a GTF file that is ready to be used for counting
>
> 4. **DEXSeq-Count** :wrench:: Count reads using **DEXSeq-Count** with
>     - HISAT2 output as "Input bam file"
>     - The formatted GTF file
> 5. Inspect the result files
>
>    > ### :question: Question
>    >
>    > Which exon has the most read mapped on it? From which gene has this exon beed extracted? Is it similar to the previous result with HTSeq-count?
>    >
>    > <details>
>    > <summary>Click to view answers</summary>
>    > FBgn0017545:004 is the exon with the most read mapped on it. It is part of FBgn0017545, the feature with the most reads mapped with HTSeq-count
>    > </details>
>    {: .question}
{: .hands_on}

## Differential exon usage

DEXSeq usage is similar to DESeq2. It uses similar statistics to find differentially used exons.

As for DESeq2, in the previous step, we counted only reads that mapped to exons on chromosome 4 and for only one sample. To be able to identify differential exon usage induced by PS depletion, all datasets (3 treated and 4 untreated) must be analyzed with the similar procedure. For time saving, we did that for you. The results are available on [Zenodo](https://dx.doi.org/10.5281/zenodo.290221):

- [dexseq.gtf](https://zenodo.org/record/290221/files/dexseq.gtf): the results of the running DEXSeq-count in 'Prepare annotation' mode
- Seven files, counts files generated in 'Count reads' mode

> ### :pencil2: Hands-on:
>
> 1. Create a new history
> 2. Import the seven count files and the dexseq.gtf from [Zenodo](https://dx.doi.org/10.5281/zenodo.290221)
>    - `dexseq.gtf`
>    - `GSM461176_untreat_single.dexseq.counts`
>    - `GSM461177_untreat_paired.dexseq.counts`
>    - `GSM461178_untreat_paired.dexseq.counts`
>    - `GSM461179_treat_single.dexseq.counts`
>    - `GSM461180_treat_paired.dexseq.counts`
>    - `GSM461181_treat_paired.dexseq.counts`
>    - `GSM461182_untreat_single.dexseq.counts`
>
> 3. **DEXSeq** :wrench:: Run **DEXSeq** with
>    - "Condition" as first factor with "treated" and "untreated" as levels and selection of count files corresponding to both levels
>    - "Sequencing" as second factor with "PE" and "SE" as levels and selection of count files corresponding to both levels
>
>    > ### :nut_and_bolt: Comment
>    >
>    > Unlike DESeq2, DEXSeq does not allow flexible primary factor names. Always use your primary factor name as "condition"
>    {: .comment}
{: .hands_on}

Similarly to DESeq2, DEXSeq generates a table with:

1.  Exon identifiers
2.  Gene identifiers
3.  Exon identifiers in the Gene
4.  Mean normalized counts, averaged over all samples from both conditions
5.  Logarithm (to basis 2) of the fold change

    The log2 fold changes are based on primary factor level 1 vs. factor level 2. The order of factor levels is then important. For example, for the factor 'Condition', DESeq2 computes fold changes of 'treated' samples against 'untreated', *i.e.* the values correspond to up- or downregulations of genes in treated samples.

6.  Standard error estimate for the log2 fold change estimate
7.  *p*-value for the statistical significance of this change
8.  *p*-value adjusted for multiple testing with the Benjamini-Hochberg procedure which controls false discovery rate ([FDR](https://en.wikipedia.org/wiki/False_discovery_rate))

> ### :pencil2: Hands-on:
>
> 1. **Filter** :wrench:: Run **Filter** to extract exons with a significant usage (adjusted *p*-value equal or below 0.05) between treated and untreated samples
>
>    > ### :question: Question
>    >
>    > How many exons have a significant change in usage between these conditions?
>    >
>    > <details>
>    > <summary>Click to view answers</summary>
>    > We get 38 exons (12.38%) with a significant usage change between treated and untreated samples
>    > </details>
>    {: .question}
{: .hands_on}

# Annotation of the result tables with gene information

Unfortunately, in the process of counting, we loose all the information of the gene except its identifiant. In order to get the information back to our final counting tables, we can use a tool to make the correspondance between identifiant and annotation.

> ### :pencil2: Hands-on:
>
> 1. **Annotate DE(X)Seq result** :wrench:: Run **Annotate DE(X)Seq result** on a counting table (from DESeq or DEXSeq) using the `Drosophila_melanogaster.BDGP5.78.gtf` as annotation file
{: .hands_on}

# Conclusion
{:.no_toc}

In this tutorial, we have analyzed real RNA sequencing data to extract useful information, such as which genes are up- or downregulated by depletion of the Pasilla gene and which genes are regulated by the Pasilla gene. To answer these questions, we analyzed RNA sequence datasets using a reference-based RNA-seq data analysis approach. This approach can be sum up with the following scheme:

![Summary of the used pipeline](../../images/rna_quantification.png "Summary of the used pipeline")
