# GALAXY WORKSHOP on RNA-seq DATA ANALYSIS  

For a more detailed version have a look at: https://github.com/nekrut/galaxy/wiki/Reference-based-RNA-seq

##Hands on!

**Slides used in this workshop can be found [here](https://github.com/bgruening/training_data/raw/master/RNAseq/RNAseq_introduction.pdf).**

**Datasets used in this exercise**

This exercise uses RNA-seq data from the study by [Brooks et al. 2011](http://genome.cshlp.org/content/21/2/193.long), in which the pasilla gene in *Drosophila melanogaster* was depleted by RNAi and the effects on splicing events were analysed by RNA-seq. The data is available at NCBI Gene Expression Omnibus (GEO) under accession number [GSE18508](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18508). Here, we analyse three treated (pasilla-depleted) and four untreated (wt) samples. It should be highlighted that each sample constitutes a separate biological replicate of the corresponding condition (treated or untreated). Also note that two of the treated and two of the untreated samples are from a paired-end sequencing assay, while the remaining samples are from a single-end sequencing experiment. 


## Spliced Read Mapping to the Reference Genome
**Step 1: Inspecting the FASTQ files**

- Create a new history for this RNA-seq exercise.

- Import a FASTQ file pair (e.g. GSM461177_untreat_paired_subset_1 and 2) with sample id from [Zenodo](http://dx.doi.org/10.5281/zenodo.61771) (select data file with "fastq" ending). Load them into Galaxy by right-clicking →  copy link location and paste the link in Galaxy → Get Data →  Upload File from your computer →  paste/fetch data → Start. 

 (Recommended: Select the correct file type ("fastqsanger") and genome ("dm3") directly in the upload dialogue. A lot of downstream programs will require these information. With the upload you can assign the correct settings for all uploaded files at once!)

 Both files contain the first 100.000 paired-end reads of one untreated sample. Rename the datasets according to the samples (recommended). As default, Galaxy takes the link as name.

- Run the tool **FastQC** on one of the two FASTQ files to control the quality of the reads. What is the read length? Is there anything what you find striking?

- Trim low quality bases from the 3' end using **Trim Galore** on both paired-end datasets. In order to use Trim Galore make sure that the file type is set to *fastqsanger* (not *fastq*), change it if necessary: click on the pencil button displayed in your dataset in the history, choose Datatype → select fastqsanger → Save.

- Re-run **FastQC** and inspect the differences.


**Step 2: Mapping of the reads with TopHat (version 2)**

Preparation:

- Tophat needs information about the type of **quality scores** in the FASTQ files. The most common type nowadays is `fastqsanger`, signalling Sanger-scaled quality scores, which are also used by the current generation of Illumina high-throughput sequencers. Make sure that the  type is set correctly.

- If you want TopHat to take advantage from already known reference gene annotations, load a reference annotation file into your current Galaxy history. For this exercise please import the Ensembl gene annotation for Drosophila melanogaster (`Drosophila_melanogaster.BDGP5.78.gtf`). From [Zenodo](http://dx.doi.org/10.5281/zenodo.61771): right-click → copy link location and paste the link in Galaxy → Upload File from your computer → paste/fetch data → Start.

- TopHat also needs to know two important parameters about the sequencing library: 1) the *strandedness* being *unstranded* or *stranded* (if *stranded* there are many types) and 2) the *inner distance* between the two reads for paired-end data. **These information should usually come with your FASTQ files!!!** If not, try to find them on the site where you downloaded the data or in the corresponding publication.

 Another option is to estimate these parameters from a *preliminary mapping* of a *downsampled* file, then run some analysis programs, and eventually redo the actual mapping on the original files with the optimized parameters. 
 

*Preliminary mapping* (**NOT necessary if you don't need to estimate parameters!!!**):

- Downsample the FASTQ file. For the provided files downsampling is not necessary as they only contain 100k reads. For real data use the Galaxy tool `Select first` to downsample to 200k to 1M reads.

- Run **TopHat** with some minimal settings. Just stick with the defaults for the settings you don't know, *strandedness* and *insert size*. As this data comes from paired-end sequencing, switch from single-end to `paired-end (as individual datasets)` and specify the FASTQ files. Align the reads to the built in reference Drosophila melanogaster `dm3` genome → Execute.

- Let's estimate the *inner distance*: Run RSeQC `Inner Distance` on the BAM file using the `Drosophila_melanogaster.BDGP5.78.gtf` reference gene model. Inspect the resulting PDF. What is the mean value for the inner distance? (If you already have read the corresponding paper carefully you might know that the fragment size is ~200bp. With read lengths of 2x37bp an educated guess could also be `125` for the inner distance. It's up to you decision, which value you prefer...)

- Check the strandedness of the library: Run RSeQC `Infer Experiment` with the same files. Check the results and search the tool's documentation for help on the meaning. As it is sometimes quite difficult to find out which settings correspond to those of other programs a table might be helpful:

 |       | RSeQC Infer Experiment |      TopHat2      |    HISAT2    | htseq-count | featureCounts |
 |:-----:|:----------------------:|:-----------------:|:------------:|:-----------:|:-------------:|
 |   PE  |    "1++,1--,2+-,2-+"   | 'fr-secondstrand' |     'FR'     |    'yes'    |      '1'      |
 |   PE  |    "1+-,1-+,2++,2--"   |  'fr-firststrand' |     'RF'     |  'reverse'  |      '2'      |
 |   SE  |         "++,--"        | 'fr-secondstrand' |      'F'     |    'yes'    |      '1'      |
 |   SE  |         "+-,-+"        |  'fr-firststrand' |      'R'     |  'reverse'  |      '2'      |
 | SE,PE |        undecided       |  'fr-unstranded'  |    default   |     'no'    |      '0'      |

*Actual Mapping*:

- Run **TopHat** with the full parameter set to get the best mapping results. Use `paired-end (as individual datasets)` and specify the FASTQ files. Set `Mean Inner Distance` to `112` (or `125`?). Select the built in reference *Drosophila melanogaster* `dm3` genome. Allow `Tophat settings to use` → `Full parameter list`: Set the correct `library type` → `FR First Strand`. Supply `own junction data` → `Yes`, `Use Gene Annotation Model` → `Yes` and select the appropriate `Gene Model Annotations` → `Drosophila_melanogaster.BDGP5.78.gtf` for your organism to enable transcriptome alignment. Enable `coverage-based search for junctions` → `Yes (--coverage-search)` for novel splice junctions to increase sensitivity.  The TopHat algorithm splits reads into segments to map the reads across splice junctions. The default `minimum length of read segments` is 25 by default, but a value of `18` seems to be a more appropriate value for this input data. Why?

 
**Step 3: Inspecting TopHat results**

- **TopHat** returns a BAM file with the mapped reads and three bed files containing splice junctions, insertions and deletions. 

- The mapping exercise worked for you? Great! However, the datasets are too small to get you a good impression of how real data looks like.

- Please therefore import the following 4 files from Tophat outputs into your history (`GSM461177_untreat_paired_chr4.bam`, `GSM461177_untreat_paired_deletions_chr4.bed`,`GSM461177_untreat_paired_insertions_chr4.bed`, `GSM461177_untreat_paired_junctions_chr4.bed`) from [Zenodo](http://dx.doi.org/10.5281/zenodo.61771) by right-clicking →  copy link location and paste the link in Galaxy → Get data → Upload File from your computer → paste/fetch data → Start.
 These files are restricted to reads that map to chr4 of Drosophila dm3.

- Visualise this bam file and the three bed files in **IGV**. You might for example inspect the region on chr4 between 560 kb to 600 kb (`chr4:560,000-600,000`). Which information does each of the bed files contain? You may have to change the data type from "tabular" to "bed" (use the pencil button). 
- Inspect the results using a **Sashimi plot** (right-click on the bam file → select `Sashimi Plot` from the context menu). Look around to find other regions with in interesting junctions, e.g. `chr4:870,000-940,000`.


## Transcript Assembly
**Step 4: Predict novel transcripts with Cufflinks**

- Import the bam files (GSM461177_untreat_paired_chr4.bam,GSM461178_untreat_paired_chr4.bam) of all samples from [Zenodo](http://dx.doi.org/10.5281/zenodo.61771) into your history by right-clicking → copy link location and paste the link in Galaxy →  Upload File from your computer →  paste/fetch data →  start.

- The tool **Cufflinks** can be used to annotate novel transcripts and splice isoforms. We will assemble the reads from the two untreated paired-end samples into transcripts. Please run **Cufflinks** separately on the bam files GSM461177_untreat_paired_chr4.bam and GSM461178_untreat_paired_chr4.bam. For real-world applications, please carefully read the [manual](http://cole-trapnell-lab.github.io/cufflinks/cufflinks/index.html) and select the appropriate options (consider especially *Use Reference Annotation*, *Perform Bias Correction* and *Use multi-read correct*).

- Merge the two resulting Cufflinks assemblies with the tool **Cuffmerge**.

- Visualise the merged predicted transcripts in **IGV** and compare them with the annotated genes. You might either download the Cuffmerge GTF file and import it to **IGV** or convert the GTF file to BED using **GFF-to-BED** tool.

## Analysing Differential Gene Expression
**Step 5: Count the number of reads per annotated gene with htseq-count**

Methods that estimate the differential expression of genes across samples require an annotation of genomic features such as genes or exons. Some methods rely on unscaled read counts per gene to evaluate the differential gene expression.
- The tool **htseq-count** can be used to count reads per features in different samples. The tools **htseq-count** expects a BAM file as input. In case of paired-end reads, the alignments in BAM should be sorted by read name. First use the tool **sort** of NGS: **SAM Tools** to sort all the paired-end BAM files. Sort by *read names*.

- Then, we need a file in GFF/GTF format with feature, i.e. gene, annotations. If not done yet, please import the Ensembl genes for dm3 →                                                                           (Drosophila_melanogaster.BDGP5.78.gtf) from [Zenodo](http://dx.doi.org/10.5281/zenodo.61771) by right-clicking →  copy link location and paste the link in Galaxy →  Upload File from your computer →  paste/fetch data →  start.

- Apply the tool **htseq-count** to the all samples. Select *dm3.ensGene.gtf* file as feature file. Use the *Union* mode for reads overlapping more than one feature. Set the *Minimum Alignment Quality* to 10.

- Inspect the result files.


**Step 6: Analyse differential gene expression with DESeq2**

- In Step 5, we counted only reads that mapped to chr4. To get more meaningful results in the following analysis, please import the 3 treated and 4 untreated count files from [Zenodo](http://dx.doi.org/10.5281/zenodo.61771) by right-clicking →  copy link location and paste the link in Galaxy →  Upload File from your computer →  paste/fetch data →  start. These files contains the read counts for all Drosophila genes and not only for reads mapped to chr4.

In our example, we have samples with two varying factors: (1) condition (either treated or untreated) and (2) sequencing type (paired-end or single-end). A multi-factor analysis allows us to assess  the effect of the treatment taking also the sequencing type into account.

- Run the tool **DESeq2** using the count files as input. In addition to the first factor *condition* with the levels *treated* and *untreated*, please add a second factor *sequencing* with the levels *PE* and *SE*. Choose the corresponding count files for each factor and level. File names have all information needed.

- Inspect the result files. The **DESeq2** tool page describes the content of the columns in the two tabular output files. The file with the independent filtering results should be used for further downstream analysis as it excludes genes with only few read counts as these genes will not be called as significantly differentially expressed.

- **Filter** for all genes from the DESeq2 result file that have a significant adjusted p-value of 0.05 or below (**Filter** tool: condition *c7<=0.05*). Please note that the output was already sorted by adjusted p-value. 

- Similarly, separate the up and down regulated genes (3rd column contains fold changes).

- Select first 100 lines of the data set.


**Step 7: Functional enrichment among differentially expressed genes**

- Use the adjusted p-value filtered data from Step 6 as input data set for **DAVID**. The identifiers in the first column are Flybase gene ids. The output of the **DAVID** tool is a HTML file with a link to the DAVID website. There, you can for example analyse cluster of functional enrichment.


-------------------------------------------

<a name="literature"/></a>
##Useful literature 

<a name="chipseq"/></a>
###ChIP-seq in general:

**Nice graphical overview of the RNA-seq processing and analysis steps:** http://figshare.com/articles/RNA_seq_Workflows_and_Tools/662782

**Auer and Doerge (2010):** [Statistical Design and Analysis of RNA Sequencing Data](http://www.genetics.org/content/185/2/405), (doi:10.1534/genetics.110.114983) - Insights into proper planning of your RNA-seq run! Read **BEFORE** your RNA-seq experiment! 

**Ian Korf (2013):**[Genomics: the state of the art in RNA-seq analysis](http://www.ncbi.nlm.nih.gov/pubmed/24296473), (doi:10.1038/nmeth.2735) - A refreshingly honest view on the non-trivial aspects of RNA-seq analysis

**Dillies et al. (2012):** [A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis](http://bib.oxfordjournals.org/content/early/2012/09/15/bib.bbs046.long), (doi:10.1093/bib/bbs046) - Systematic comparison of seven representative normalization methods for the differential analysis of RNA-seq data (Total Count, Upper Quartile, Median (Med), DESeq, edgeR, Quantile and Reads Per Kilobase per Million mapped reads (RPKM) normalization)

**Rapaport et al. (2013):** [Comprehensive evaluation of differential gene expression analysis methods for RNA-seq data](http://www.genomebiology.com/2013/14/9/R95/abstract), (doi:10.1186/gb-2013-14-9-r95) - Evaluation of methods for differential gene expression analysis

**Soneson et al. (2013):** [A comparison of methods for differential expression analysis of RNA-seq data](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-91), (doi:10.1186/1471-2105-14-91)

**Roberts et al. (2011):** [Improving RNA-Seq expression estimates by correcting for fragment bias Fragment bias correction](http://www.ncbi.nlm.nih.gov/pubmed/21410973), (doi: 10.1186/gb-2011-12-3-r22)

**Garber et al. (2011):** [Computational methods for transcriptome annotation and quantification using RNA-seq](http://www.ncbi.nlm.nih.gov/pubmed/21623353), (doi: 10.1038/nmeth.1613) - Classical paper about the computational aspects of RNA-seq data analysis
