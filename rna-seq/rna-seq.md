# GALAXY WORKSHOP on RNA-seq DATA ANALYSIS  

##Literature

**Nice graphical overview of the RNA-seq processing and analysis steps:** http://figshare.com/articles/RNA_seq_Workflows_and_Tools/662782

**Auer and Doerge (2010):** [Statistical Design and Analysis of RNA Sequencing Data](http://www.genetics.org/content/185/2/405), (doi:10.1534/genetics.110.114983) - Insights into proper planning of your RNA-seq run! Read **BEFORE** your RNA-seq experiment! 

**Ian Korf (2013):**[Genomics: the state of the art in RNA-seq analysis](http://www.ncbi.nlm.nih.gov/pubmed/24296473), (doi:10.1038/nmeth.2735) - A refreshingly honest view on the non-trivial aspects of RNA-seq analysis

**Dillies et al. (2012):** [A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis](http://bib.oxfordjournals.org/content/early/2012/09/15/bib.bbs046.long), (doi:10.1093/bib/bbs046) - Systematic comparison of seven representative normalization methods for the differential analysis of RNA-seq data (Total Count, Upper Quartile, Median (Med), DESeq, edgeR, Quantile and Reads Per Kilobase per Million mapped reads (RPKM) normalization)

**Rapaport et al. (2013):** [Comprehensive evaluation of differential gene expression analysis methods for RNA-seq data](http://www.genomebiology.com/2013/14/9/R95/abstract), (doi:10.1186/gb-2013-14-9-r95) - Evaluation of methods for differential gene expression analysis

**Soneson et al. (2013):** [A comparison of methods for differential expression analysis of RNA-seq data](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-91), (doi:10.1186/1471-2105-14-91)

**Roberts et al. (2011):** [Improving RNA-Seq expression estimates by correcting for fragment bias Fragment bias correction](http://www.ncbi.nlm.nih.gov/pubmed/21410973), (doi: 10.1186/gb-2011-12-3-r22)

**Garber et al. (2011):** [Computational methods for transcriptome annotation and quantification using RNA-seq](http://www.ncbi.nlm.nih.gov/pubmed/21623353), (doi: 10.1038/nmeth.1613) - Classical paper about the computational aspects of RNA-seq data analysis


##Hands on!

**Slides used in this workshop can be found [here](https://github.com/bgruening/training_data/raw/master/RNAseq/RNAseq_introduction.pdf).**

**Datasets used in this exercise**

This exercise uses RNA-seq data from the study by [Brooks et al. 2011](http://genome.cshlp.org/content/21/2/193.long), in which the pasilla gene in *Drosophila melanogaster* was depleted by RNAi and the effects on splicing events were analysed by RNA-seq. The data is available at NCBI Gene Expression Omnibus under accession number [GSE18508](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18508). Here, we analyse three treated (pasilla-depleted) and four untreated (wt) samples. It should be highlighted that each sample constitutes a separate biological replicate of the corresponding condition (treated or untreated). Also note that two of the treated and two of the untreated samples are from a paired-end sequencing assay, while the remaining samples are from a single-end sequencing experiment. 


## Spliced Read Mapping to the Reference Genome
**Step 1: Inspecting the FASTQ files**

- Create a new history for this RNA-seq exercise.

- Import any two FASTQ files with sample id by right-clicking →  copy link location and paste the link in Galaxy →  Upload File from your computer →  paste/fetch data →  start:

      - [GSM461176_untreat_single_subset](https://www.googledrive.com/host/0B9urRnOAUUI8TEFrWTFHcy0tQWs)
      - [GSM461177_untreat_paired_subset_1](https://www.googledrive.com/host/0B9urRnOAUUI8NXB1dTVOdWRScDA)
      - [GSM461177_untreat_paired_subset_2](https://www.googledrive.com/host/0B9urRnOAUUI8MXlQYUk2LWhkUlk)
      - [GSM461178_untreat_paired_subset_1](https://www.googledrive.com/host/0B9urRnOAUUI8cnpsem9mWml4OVk)
      - [GSM461178_untreat_paired_subset_2](https://www.googledrive.com/host/0B9urRnOAUUI8Y1BuYVdQNEhYRE0)
      - [GSM461179_treat_single_subset](https://www.googledrive.com/host/0B9urRnOAUUI8cXd2Qmgya3NhdDA)
      - [GSM461180_treat_paired_subset_1](https://www.googledrive.com/host/0B9urRnOAUUI8SU1XdjJDX016MVk)
      - [GSM461180_treat_paired_subset_2](https://www.googledrive.com/host/0B9urRnOAUUI8a3NIbWR6dS1Vb0k)
      - [GSM461181_treat_paired_subset_1](https://www.googledrive.com/host/0B9urRnOAUUI8TGhhR01PemYtbzQ)
      - [GSM461181_treat_paired_subset_2](https://www.googledrive.com/host/0B9urRnOAUUI8VHN4RURLcVY3X28)
      - [GSM461182_untreat_single_subset](https://www.googledrive.com/host/0B9urRnOAUUI8YlBIcXhEclhINlk)


 These two files contain the first 100.000 paired-end reads of one untreated sample. Rename the datasets   according to the samples if you like. As default, Galaxy takes the link as name.

- Run the tool **FastQC** on one of the two FASTQ files to control the quality of the reads. What is the read length? Is there anything what you find striking?

- Trim low quality bases from the 3' end using **Trim Galore** on both paired-end datasets.

- Re-run **FastQC** and inspect the differences.


**Step 2: Mapping of the reads with TopHat**

- Annoying but necessary: Before continuing with the mapping step, change the file type of the FASTQ files to *fastqsanger*. This signals that the FASTQ contains Sanger-scaled quality values. **TopHat2** will otherwise not be able to accept the files for input! This can be done by clicking on the pencil item of the dataset (edit attributes) →  Datatype →  select fastqsanger

- If you want **TopHat** to take advantage from already known reference gene annotations, make sure that the data is available in your current Galaxy history. For this exercise, please import the Ensembl gene annotation for [dm3.ensGene.gtf](https://www.googledrive.com/host/0B9urRnOAUUI8Tml0MzBabEE5MUU) by right-clicking →  copy link location and paste the link in Galaxy →  Upload File from your computer →  paste/fetch data →  start.

- Run **TopHat2** with the two trimmed FASTQ files as input for forward and reverse reads as your reads are from a paired-end sequencing run. Align the reads to the Drosophila dm3 genome. The fragment size is 200 and the read length is 37, so why is 125 a good value for the mean inner distance? How do you know the fragment size (Hint: Go to the GEO entry of your data set, browse to the SRA entry of your sample and explore the library information)? 

  Enable the *Coverage Search* for novel splice junctions to increase sensitivity. Also enable *Use Own Junctions*, *Use Gene Annotation Model* and select the appropriate *Gene Model Annotations* (dm3.ensGene.gtf) for your organism to enable the transcriptome alignment of **TopHat**.
  
  The TopHat algorithm splits reads into segments to map the reads across splice junctions. The default minimum length of read segments is 25, but a value of 18 seems to be a more appropriate value for this input data. Why?


**Step 3: Inspecting the TopHat results**

- **TopHat2** returns a bam file with the mapped reads and three bed files containing splice junctions, insertions and deletions. 

- The exercise of step 2 worked for you? Fine! However the dataset might be a little too small to get you a good impression of how real data looks like. Please therefore import all the files from tophat2 outputs into your history by right-clicking →  copy link location and paste the link in Galaxy →  Upload File from your computer →  paste/fetch data →  start.

      - [GSM461177_untreat_paired_chr4.bam](https://www.googledrive.com/host/0B9urRnOAUUI8ZW5VWnRhUzNqUGc)
      - [GSM461177_untreat_paired_deletions_chr4.bed](https://www.googledrive.com/host/0B9urRnOAUUI8bHpSRW5ZUlR6S0k)
      - [GSM461177_untreat_paired_insertions_chr4.bed](https://www.googledrive.com/host/0B9urRnOAUUI8MTBQbVZYSGhIU3M)
      - [GSM461177_untreat_paired_junctions_chr4.bed](https://www.googledrive.com/host/0B9urRnOAUUI8ZmZDZkIxcU4zSFU)


These files contain the TopHat results for the sample GSM461177_untreat_paired, but are restricted to reads that map to chr4 of Drosophila dm3.

- Visualise this bam file and the three bed files in **IGV**. You might for example inspect the region between 560 kb to 600 kb on chr4. Which information does each of the bed files contain? 

- Also inspect the results using a **Sashimi plot** (activate by right-clicking on the reads).


## Transcript Assembly
**Step 4: Predict novel transcripts with Cufflinks**

- Import the bam files of all samples into your history by right-clicking →  copy link location and paste the link in Galaxy →  Upload File from your computer →  paste/fetch data →  start.

      - [GSM461177_untreat_paired_chr4.bam](https://www.googledrive.com/host/0B9urRnOAUUI8RlNiSkd0TG9TeE0)
      - [GSM461178_untreat_paired_chr4.bam](https://www.googledrive.com/host/0B9urRnOAUUI8b09JSF9rM3REeEU)


- The tool **Cufflinks** can be used to annotate novel transcripts and splice isoforms. We will assemble the reads from the two untreated paired-end samples into transcripts. Please run **Cufflinks** separately on the bam files GSM461177_untreat_paired_chr4.bam and GSM461178_untreat_paired_chr4.bam. For real-world applications, please carefully read the [manual](http://cole-trapnell-lab.github.io/cufflinks/cufflinks/index.html) and select the appropriate options (consider especially *Use Reference Annotation*, *Perform Bias Correction* and *Use multi-read correct*).

- Merge the two resulting Cufflinks assemblies with the tool **Cuffmerge**.

- Visualise the merged predicted transcripts in **IGV** and compare them with the annotated genes. You might either download the Cuffmerge GTF file and import it to **IGV** or convert the GTF file to BED using **GFF-to-BED** tool.

## Analysing Differential Gene Expression
**Step 5: Count the number of reads per annotated gene with htseq-count**

Methods that estimate the differential expression of genes across samples require an annotation of genomic features such as genes or exons. Some methods rely on unscaled read counts per gene to evaluate the differential gene expression.
- The tool **htseq-count** can be used to count reads per features in different samples. The tools **htseq-count** expects a BAM file as input. In case of paired-end reads, the alignments in BAM should be sorted by read name. First use the tool **sort** of NGS: **SAM Tools** to sort all the paired-end BAM files. Sort by *read names*.

- Then, we need a file in GFF/GTF format with feature, i.e. gene, annotations. If not done yet, please import the Ensembl genes for dm3 →                                                                           [dm3.ensGene.gtf](https://www.googledrive.com/host/0B9urRnOAUUI8Tml0MzBabEE5MUU) by right-clicking →  copy link location and paste the link in Galaxy →  Upload File from your computer →  paste/fetch data →  start.

- Apply the tool **htseq-count** to the all samples. Select *dm3.ensGene.gtf* file as feature file. Use the *Union* mode for reads overlapping more than one feature. Set the *Minimum Alignment Quality* to 10.

- Inspect the result files.


**Step 6: Analyse differential gene expression with DESeq2**

- In Step 5, we counted only reads that mapped to chr4. To get more meaningful results in the following analysis, please import the 3 treated and 4 untreated count files from the data library Galaxy course →  RNA-seq →  count_matrix. These files contains the read counts for all Drosophila genes and not only for reads mapped to chr4.

In our example, we have samples with two varying factors: (1) condition (either treated or untreated) and (2) sequencing type (paired-end or single-end). A multi-factor analysis allows us to assess  the effect of the treatment taking also the sequencing type into account.

- Run the tool **DESeq2** using the count files as input. In addition to the first factor *condition* with the levels *treated* and *untreated*, please add a second factor *sequencing* with the levels *PE* and *SE*. Choose the corresponding count files for each factor and level. File names have all information needed.

- Inspect the result files. The **DESeq2** tool page describes the content of the columns in the two tabular output files. The file with the independent filtering results should be used for further downstream analysis as it excludes genes with only few read counts as these genes will not be called as significantly differentially expressed.

- **Filter** for all genes from the DESeq2 result file that have a significant adjusted p-value of 0.05 or below (**Filter** tool: condition *c7<=0.05*). Please note that the output was already sorted by adjusted p-value. 

- Similarly, separate the up and down regulated genes (3rd column contains fold changes).

- Select first 100 lines of the data set.


**Step 7: Functional enrichment among differentially expressed genes**

- Use the adjusted p-value filtered data from Step 6 as input data set for **DAVID**. The identifiers in the first column are Flybase gene ids. The output of the **DAVID** tool is a HTML file with a link to the DAVID website. There, you can for example analyse cluster of functional enrichment.
