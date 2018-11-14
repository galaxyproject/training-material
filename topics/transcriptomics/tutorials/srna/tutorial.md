---
layout: tutorial_hands_on

title: "Differential abundance testing of small RNAs"
zenodo_link: "https://zenodo.org/record/826906"
questions:
  - "What small RNAs are expressed?"
  - "What RNA features have significantly different numbers of small RNAs targeting them between two conditions?"
objectives:
  - "Process small RNA-seq datasets to determine quality and reproducibility."
  - "Filter out contaminants (e.g. rRNA reads) in small RNA-seq datasets."
  - "Differentiate between subclasses of small RNAs based on their characteristics."
  - "Identify differently abundant small RNAs and their targets."
time_estimation: "3h"
key_points:
  - "Analysis of small RNAs is complex due to the diversity of small RNA subclasses."
  - "Both alignment to references and transcript quantification approaches are useful for small RNA-seq analyses."
contributors:
  - malloryfreeberg
---

# Introduction
{:.no_toc}

Small, noncoding RNA (sRNA) molecules, typically 18-40nt in length, are key features of post-transcriptional regulatory mechanisms governing gene expression. Through interactions with protein cofactors, these tiny sRNAs typically function by perfectly or imperfectly basepairing with substrate RNA molecules, and then eliciting downstream effects such as translation inhibition or RNA degradation. Different subclasses of sRNAs - *e.g.* microRNAs (miRNAs), Piwi-interaction RNAs (piRNAs), and endogenous short interferring RNAs (siRNAs) - exhibit unique characteristics, and their relative abundances in biological contexts can indicate whether they are active or not. In this tutorial, we will examine expression of the piRNA subclass of sRNAs and their targets in *Drosophila melanogaster*.

The data used in this tutorial are from polyphosphatase-treated sRNA sequencing (sRNA-seq) experiments in *Drosophila*. The goal of this study was to determine how siRNA expression changes in flies treated with RNA interference (RNAi) to knock down Symplekin, which is a component of the core cleaveage completx. To that end, in addition to sRNA-seq, mRNA-seq experiments were performed to determine whether targets of differentially expressed siRNAs were also differentially expressed. The original published study can be found [here](https://www.ncbi.nlm.nih.gov/pubmed/28415970). Because of the long processing time for the large original files - which contained 7-22 million reads each - we have downsampled the original input data to include only a subset of usable reads.

# Analysis strategy

In this exercise we will identify what small RNAs are present in *Drosophila* Dmel-2 tissue culture cells treated with either RNAi against  core cleavage complex component Symplekin or blank (control) RNAi (published data available in GEO at [GSE82128](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82128)). In this study, biological triplicate sRNA- and mRNA-seq libraries were sequenced for both RNAi conditions. After removing contaminant ribosomal RNA (rRNA) and miRNA reads, we will identify and quantify endogenous siRNAs from sequenced reads and test for differential abundance. We will follow a popular small RNA analysis pipeline developed for piRNAs by the Phillip Zamore Lab and the ZLab (Zhiping Weng) at UMass Med School called [PiPipes](https://github.com/bowhan/piPipes). Although PiPipes was developed for analysis of piRNAs, many of the basical principles can be applied to other classes of small RNAs.

It is of note that this tutorial uses datasets that have been de-multiplexed so that the input data files are a single FASTQ formatted file for each sample. Because sRNAs are typically much smaller than fragments generated for RNA-seq or other types of deep sequencing experiments, single-end sequencing strategies are almost always used to sequence sRNAs. This tutorial uses the *Collections* feature of Galaxy to orgainze each set of replicates into a single group, making tool form submission easier and ensuring reproducibility.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Data upload and organization

Due to the large size of the original sRNA-seq datasets, we have downsampled them to only inlcude a subset of reads. These datasets are avaialble at [`Zenodo`](https://doi.org/10.5281/zenodo.826906) where you can find the FASTQ files corresponding to replicate sRNA-seq experiments and additional annotation files for the *Drosophila melanogaster* genome version dm3.

> ### {% icon hands_on %} Hands-on: Data upload and organization
>
> 1. Create a new history and name it something meaningful (*e.g.* sRNA-seq tutorial)
>
>    {% include snippets/create_new_history.md %}
>    {% include snippets/rename_history.md %}
>
> 2. Import the 3 `Blank_RNAi_sRNA-seq` and 3 `Symp_RNAi_sRNA-seq` FASTQ files from [Zenodo](https://zenodo.org/record/1324070) or from the data library (ask your instructor)
>
>    ```
https://zenodo.org/record/826906/files/Blank_RNAi_sRNA-seq_rep1_downsampled.fastqsanger.gz
https://zenodo.org/record/826906/files/Blank_RNAi_sRNA-seq_rep2_downsampled.fastqsanger.gz
https://zenodo.org/record/826906/files/Blank_RNAi_sRNA-seq_rep3_downsampled.fastqsanger.gz
https://zenodo.org/record/826906/files/Symp_RNAi_sRNA-seq_rep1_downsampled.fastqsanger.gz
https://zenodo.org/record/826906/files/Symp_RNAi_sRNA-seq_rep2_downsampled.fastqsanger.gz
https://zenodo.org/record/826906/files/Symp_RNAi_sRNA-seq_rep3_downsampled.fastqsanger.gz
>    ```
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
>    Set the datatype of the read (.fastqsanger) files to **fastq**
>    As default, Galaxy takes the link as name, so rename them.
>
> 2. Build a dataset list for both sets (`Blank` and `Symp`) of replicate FASTQ files
>
>    {% include snippets/build_dataset_list.md %}
>
>    ![History section with all the loaded datasets](../../images/sRNA/Fig3_dataset_lists_made.png)
> 
> 3. Import the remaining files from [Zenodo](https://zenodo.org/record/1324070) or from the data library (ask your instructor)
>
>    ```
https://zenodo.org/record/826906/files/dm3_transcriptome_Tx2Gene_downsampled.tab.gz
https://zenodo.org/record/826906/files/dm3_transcriptome_sequences_downsampled.fa.gz
https://zenodo.org/record/826906/files/dm3_miRNA_hairpin_sequences.fa.gz
https://zenodo.org/record/826906/files/dm3_rRNA_sequences.fa.gz
>    ```
>
>    - Set the datatype of the annotation file to **tab** and assign the Genome as **dm3**
>    - Set the datatype of the reference files to **fasta** and assign the Genome as **dm3**
>
{: .hands_on}

## Read quality checking

Read quality scores (phred scores) in FASTQ-formatted data can be encoded by one of a few different encoding schemes. Most Galaxy tools assume that input FASTQ files are using the Sanger/Illumina 1.9 encoding scheme. If the input FASTQ files are using an alternate encoding scheme, then some tools will not interpret the quality score encodings correctly. It is good practice to confirm the quality encoding scheme of your data and then convert to Sanger/Illumina 1.9, if necessary. We can check the quality encoding scheme using the [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) tool (further described in the [NGS-QC tutorial]({{site.baseurl}}/topics/sequence-analysis)).


> ### {% icon hands_on %} Hands-on: Quality checking
>
> 1. Run **FastQC** {% icon tool %} on the `Blank` collection
>    - {% icon param-collection %} *"Short read data from your current history"*: `Blank` collection
>
>    {% include snippets/select_collection.md %}
>
> 2. Run **FastQC** {% icon tool %} for the *Symplekin* RNAi dataset collection
> 3. Inspect the **FastQC** {% icon tool %} reports
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What quality score encoding scheme is being used for each sample?
>    > 2. What is the read length for each sample?
>    > 3. What does the base/read quality look like for each sample?
>    >
>    >    ![The output of FastQC as a boxplot](../../images/sRNA/Fig5a_fastqc_result_stats_boxplot.png)
>    >
>    > 4. Are there any adaptors present in these reads? Which one(s)?
>    >
>    >    ![The output of FastQC stats adapted](../../images/sRNA/Fig5b_fastqc_result_stats_adapter.png)
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. All samples use the Illumina 1.9 quality encoding scheme, so we do **not** need to convert. 
>    > > 2. All samples have a read length of 50 nt. 
>    > > 3. The base quality across the entire length of the reads is good (phred score > 28 for the most part). 
>    > > 4. Yes, the "Illumina Universal Adapter" is present. 
>    > {: .solution }
>    {: .question}
{: .hands_on}

> ### {% icon hands_on %} Hands-on: Formatting 
>
> If your data are not in **Sanger / Illumina 1.9 format**:
> 1. **FASTQ Groomer** {% icon tool %} on each collection
>    - {% icon param-collection %} *"File to groom"*: input collection
>    - *"Input FASTQ quality scores type"*: `Illumina 1.3-1.7`
> 2. Check that the format is **fastqsanger** instead of **fastq**
>
> If your data are already in **Sanger / Illumina 1.9 format**:
> 1. Change the datatype to `fastqsanger`
>
>    {% include snippets/change_datatype.md datatype="fastqsanger" %}
>
{: .hands_on}

Go back to the **FastQC** {% icon tool %} and scroll down to the "Adapter Content" section. You can see that Illumina Universal Adapters are present in **>80%** of the reads. The next step is to remove these artificial adaptors because they are not part of the biological sRNAs. If a different adapter is present, you can update the **Adapter sequence to be trimmed off** in the **Trim Galore!** {% icon tool %} step described below.

## Adaptor trimming

sRNA-seq library preparation involves adding an artificial adaptor sequence to both the 5' and 3' ends of the small RNAs. While the 5' adaptor anchors reads to the sequencing surface and thus are not sequenced, the 3' adaptor is typically sequenced immediately following the sRNA sequence. In the example datasets here, the 3' adaptor sequence is identified as the Illumina Universal Adapter, and needs to be removed from each read before aligning to a reference. We will be using the Galaxy tool **Trim Galore!** {% icon tool %} which implements the [`cutadapt`](https://cutadapt.readthedocs.io/en/stable/) tool for adapter trimming.

> ### {% icon hands_on %} Hands-on: Adaptor trimming
>
> 1. **Trim Galore!** {% icon tool %} on each collection of FASTQ read files to remove Illumina adapters from the 3' ends of reads with the following parameters:
>    - *"Is this library paired- or single-end"*: `Single-end`
>       - {% icon param-collection %} *"Reads in FASTQ format?"*: `Blank` collection
>       - *"Adapter sequence to be trimmed"*: `Illumina universal`
>    - *"Trim Galore! advanced settings"*: `Full parameter list`
>       - *"Trim low-quality ends from reads in addition to adapter removal"*: `0`
>       - *"Overlap with adapter sequence required to trim a sequence"*: `6`
>       - *"Discard reads that became shorter than length N"*: `12`
>       - *"Generate a report file"*: `Yes`
>
> 2. **Trim Galore!** {% icon tool %}: Repeat for the *Symplekin* RNAi dataset collection
>
{: .hands_on}

We don't want to trim for quality because the adapter-trimmed sequences represent a full small RNA molecule, and we want to maintain the integrity of the entire molecule. We increase the minimum read length required to keep a read because small RNAs can potentially be shorter than 20 nt (the default value). We can check out the generated report file for any sample and see the command for the tool, a summary of the total reads processed and number of reads with an adapter identified, and histogram data of the length of adaptor trimmed. We also see that a very small percentage of low-quality bases have been trimmed

> ### {% icon hands_on %} Hands-on: Adaptor trimming check
>
> 1. **FastQC** {% icon tool %} on each collection of trimmed FASTQ read files to assess whether adapters were successfully removed.
> 2. Inspect the **FastQC** {% icon tool %} reports
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What is the read length?
>    >
>    >   ![FastQC output after trimming](../../images/sRNA/Fig8a_fastqc_post_trimming_result.png)
>    >
>    > 2. Are there any adaptors present in these reads? Which one(s)?
>    >
>    >   ![FastQC output after trimming](../../images/sRNA/Fig8b_fastqc_post_trimming_lengths.png)
>    > > ### {% icon solution %} Solution
>    > > 1. The read lengths range from 12 to 50 nt after trimming.
>    > > 2. No, Illumina Universal Adaptors are no longer present. No other adapters are present.
>    > {: .solution }
>    {: .question}
{: .hands_on}

An interesting thing to note from our **FastQC** {% icon tool %} results is the *Sequence Length Distribution* results. While many RNA-seq experiments have normal distribution of read lengths, an unusual spike at 22nt is observed in our data. This spike represents the large set of endogenous siRNAs that occur in the cell line used in this study, and is actually confirmation that our dataset captures the biological molecule we are interested in.

Now that we have converted to *fastqsanger* format and trimmed our reads of the Illumina Universal Adaptors, we will align our trimmed reads to reference *Drosophila* rRNA and miRNA sequences (dm3) to remove these artifacts. Interestingly, *Drosophila* have a short 2S rRNA sequence that is 30nt long and typically co-migrates with sRNA populations during gel electrophoresis. rRNAs make up a very large proportion of all non-coding RNAs, and thus need to be removed. Oftentimes, experimental approaches can be utilized to deplete or avoid capture of rRNAs, but these methods are not always 100% efficient. We also want to remove any miRNA sequences, as these are not relevant to our analysis. After removing rRNA and miRNA reads, we will analyze the remaining reads, which may be siRNA or piRNA sequences.

## Hierarchical read alignment to remove rRNA/miRNA reads

To first identify rRNA-originating reads (which we are not interested in in this case), we will align the reads to a reference set of rRNA sequences using [**HISAT2** {% icon tool %}](https://ccb.jhu.edu/software/hisat2/index.shtml), an accurate and fast tool for aligning reads to a reference.

> ### {% icon hands_on %} Hands-on: Hierarchical alignment to rRNA and miRNA reference sequences
>
> 1. **HISAT2** {% icon tool %} on each collection of trimmed reads to align to reference rRNA sequences with the following parameters
>    - *"Source for the reference genome to align against"*: `Use a genome from history`
>       - {% icon param-file %} *"Select the reference genome"*: `Dmel_rRNA_sequences.fa`
>    - *"Is this a single or paired library?"*: `Single-end`
>       - {% icon param-collection %} *"FASTA/Q file"*: output of **Trim Galore!** {% icon tool %} for `Blank` collection
>       - *"Specify strand-specific information"*: `Forward (F)`
>
> 2. **HISAT2** {% icon tool %}: Repeat for the *Symplekin* RNAi dataset collection
{: .hands_on}

We now need to extract the *unaligned* reads from the output BAM file for aligning to reference miRNA sequences. We can do this by using the `Filter SAM or BAM, output SAM or BAM` tool to obtain reads with a bit flag = 4 (meaning the read is unaligned) and then converting the filtered BAM file to FASTQ format with the `Convert from BAM to FastQ` tool.

> ### {% icon hands_on %} Hands-on: Hierarchical alignment to rRNA and miRNA reference sequences (2)
> 1. **Filter SAM or BAM, output SAM or BAM** {% icon tool %} on collection of HISAT2 output BAM files with the following parameters
>    - {% icon param-collection %} *"SAM or BAM file to filter"*: output of **HISAT2** {% icon tool %} for `Blank` collection
>    - *"Filter on bitwise flag"*: `Yes`
>       - *"Only output alignments with all of these flag bits set"*: `The read in unmapped`
>
> 2. **Filter SAM or BAM, output SAM or BAM** {% icon tool %}: Repeat for the *Symplekin* RNAi dataset collection
> 3. **Convert from BAM to FastQ** {% icon tool %} on each collection of filtered HISAT2 output BAM files with the following parameters:
>    - {% icon param-collection %} *"Convert the following BAM file to FASTQ"*: output of **Filter SAM or BAM** {% icon tool %} for `Blank` collection
>
> 4. **Convert from BAM to FastQ** {% icon tool %}: Repeat for the *Symplekin* RNAi dataset collection
{: .hands_on}

Next we will align the non-rRNA reads to a known set of miRNA hairpin sequences to identify miRNA reads.

> ### {% icon hands_on %} Hands-on: Hierarchical alignment to rRNA and miRNA reference sequences (3)
> 1. **HISAT2** {% icon tool %} on each collection of filtered HISAT2 output FASTQ files to align non-rRNA reads to reference miRNA hairpin sequences using the following parameters:
>    - *"Source for the reference genome to align against"*: `Use a genome from history`
>       - {% icon param-file %} *"Select the reference genome"*: `Dmel_miRNA_sequences.fa`
>    - *"Is this a single or paired library?"*: `Single-end`
>       - {% icon param-collection %} *"FASTA/Q file"*: output of **Convert from BAM to FastQ** {% icon tool %} for `Blank` collection
>       - *"Specify strand-specific information"*: `Forward (F)`
> 2. **HISAT2** {% icon tool %}: Repeat for the *Symplekin* RNAi dataset collection
{: .hands_on}

For this tutorial we are not interested in miRNA reads, so we need to extract *unaligned* reads from the output BAM files. To do this, repeat the **Filter SAM or BAM, output SAM or BAM** {% icon tool %} steps for each dataset collection. Finally, rename the converted FASTQ files something meaningful (*e.g.* "non-r/miRNA control RNAi sRNA-seq").

In some instances, miRNA-aligned reads are the desired output for downstream analyses, for example, if we are investigating how loss of key miRNA pathway components affect levels of pre-miRNA and mature miRNAs. In this case, the BAM output of HISAT2 can directly be used in subsequent tools as it contains the subset of reads aligned to miRNA sequences.

## Small RNA subclass distinction

In *Drosophila*, non-miRNA small RNAs are typically divided into two major groups: endogenous siRNAs which are 20-22nt long and piRNAs which are 23-29nt long. We want to analyze these sRNA subclasses independently, so next we are going to filter the non-r/miRNA reads based on length using the **Manipulate FASTQ** {% icon tool %} tool.

> ### {% icon hands_on %} Hands-on: Extract subclasses
>
> 1. **Manipulate FASTQ** {% icon tool %} on each collection of non-r/miRNA reads to identify siRNAs (20-22nt) using the following parameters.
>    - {% icon param-collection %} *"FASTQ File"*: output of last **HISAT** {% icon tool %} for `Blank` collection
>    - Click on *"Insert Match Reads"*
>       - In "1: Match Reads"
>         - *"Match Reads by"*: `Sequence Content`
>           - *"Match by"*: `^.{12,19}$|^.{23,50}$`
>    - Click on *"Manipulate Reads"*
>       - In "1: Manipulate Reads"
>         - *"Manipulate Reads on"*: `Miscellaneous actions`
>           - *"Miscellaneous Manipulation Type"*: `Remove Read`
> 2. **Manipulate FASTQ** {% icon tool %}: Repeat for the *Symplekin* RNAi dataset collection
> 3. Rename each resulting dataset collection something meaningful (*i.e.* "blank RNAi - siRNA reads (20-22nt)")
{: .hands_on}

The regular expression in the **Match by** parameter tells the tool to identify sequences that are length 12-19 or 23-50 (inclusive), and the **Miscellaneous Manipulation Type** parameter tells the tool to remove these sequences. What remains are sequences of length 20-22nt. We will now repeat these steps to identify sequences in the size-range of piRNAs (23-29nt).

> ### {% icon hands_on %} Hands-on: Extract subclasses
> 1. **Manipulate FASTQ** {% icon tool %}: Run `Manipulate FASTQ` on each collection of non-r/miRNA reads to identify sequences in the size-range of piRNAs (23-29nt)
>    - {% icon param-collection %} *"FASTQ File"*: output of last **HISAT** {% icon tool %} for `Blank` collection
>    - Click on *"Insert Match Reads"*
>       - In "1: Match Reads"
>         - *"Match Reads by"*: `Sequence Content`
>           - *"Match by"*: `^.{12,22}$|^.{30,50}$`
>    - Click on *"Manipulate Reads"*
>       - In "1: Manipulate Reads"
>         - *"Manipulate Reads on"*: `Miscellaneous actions`
>           - *"Miscellaneous Manipulation Type"*: `Remove Read`
> 2. **Manipulate FASTQ** {% icon tool %}: Repeat for the *Symplekin* RNAi dataset collection
> 3. Rename each resulting dataset collection something meaningful (*i.e.* "blank RNAi - piRNA reads (23-29nt)")
{: .hands_on}

> ### {% icon hands_on %} Hands-on: Check the read length
> 1. **FastQC** {% icon tool %}: Run `FastQC` on each collection of siRNA and piRNA read files to confirm the correct read lengths.
{: .hands_on}

![FastQC result on both collections](../../images/sRNA/Fig13_FastQC_piRNA_siRNA_result_Blank_rep1.png)

We see in the above image (for Blank RNAi, replicate 3) that we have ~32k reads that are 20-22nt long in our Reads_20-22_nt_length_siRNAs file and ~16k reads that are 23-29nt long in our Reads_23-29_nt_length_piRNAs file, as expected. Given that we had ~268k trimmed reads in the Blank RNAi replicate 3 file, siRNAs represent ~12.1% of the reads and piRNAs represent ~6.2% of the reads in this experiment. There is no distinct peak in the piRNA range (you can find the example of the real piRNA profile in Figure 1 from [Colin D. Malone et al, 2009](http://dx.doi.org/10.1016/j.cell.2009.03.040).

The next step in our analysis pipeline is to identify which RNA features - e.g. protein-coding mRNAs, transposable elements - the siRNAs and piRNAs align to. Most fly piRNAs originate from and thus align to transposable elements (TEs) but some also originate from genome piRNA clusters or protein-coding genes. siRNAs are thought to inhibit TE mobility and potentially target mRNAs for degradation via an RNA-interference (RNAi) mechanism. Ultimately, we want to know whether an RNA feature (TE, piRNA cluster, mRNA, etc.) has significantly different numbers of siRNAs or piRNAs targeting it. To determine this, we will use **Salmon** {% icon tool %} to simultaneously align and quantify siRNA reads against known target sequences to estimate siRNA abundance per target. For the remainder of the tutorial we will be focusing on siRNAs, but similar approaches can be done for piRNAs, miRNAs, or any other subclass of small, non-coding RNAs.

## siRNA abundance estimation

We want to identify which siRNAs are differentially abundance between the blank and *Symplekin* RNAi conditions. To do this we will implement a counting approach using **Salmon** {% icon tool %} to quantify siRNAs per RNA feature. Specifically, we will analyzing mRNA and TE elements by counting siRNAs against a FASTA file of transcript sequences, not the reference genome. This approach is especially useful in this case because we are interested in TEs, which occur at many copies (100s - 1000s) in the genome, and thus will result in multiple alignments if we align to the genome with a tool like **HISAT2** {% icon tool %}. Further, we will be counting abundance of siRNAs that align *sense* or *antisense* to an RNA feature independently, as opposite sense alignments are correlated with unique downstream silencing effects. Then, we will provide this information to **DESeq2** {% icon tool %} to generate normalized counts and significance testing for differential abundance of siRNAs per feature.

> ### {% icon hands_on %} Hands-on: siRNA abundance estimation
>
> 1. **Salmon** {% icon tool %} on each collection of siRNA reads (20-22nt) to quantify the abundance of *antisense* siRNAs at relevant targets
>    - *"Select a reference transcriptome from your history or use a built-in index?"*: `Use one from the history`
>      - {% icon param-file %} *"Select the reference transcriptome"*: `the reference mRNA and TE fasta file`
>      - *"The size should be odd number"*: `19`
>    - *"Is this library mate-paired?"*: `Single-end`
>      - {% icon param-collection %} *"FASTQ/FASTA file"*: collection with the Blank RNAi siRNA (20-22nt) reads
>      - *"Specify the strandedness of the reads"*: `read 1 (or single-end read) comes from the reverse strand (SR)`
>    - In *"Additional Options"*
>      - *"Incompatible Prior"*: `0`
>
> 1. **Salmon** {% icon tool %}: Repeat for the *Symplekin* RNAi siRNAs (20-22nt) reads dataset collection
> 1. **Salmon** {% icon tool %}: Repeat on each collection of siRNA reads (20-22nt) to quantify the abundance of *sense* siRNAs at relevant targets by changing the following parameters:
>    - *"Specify the strandedness of the reads"*: `read 1 (or single-end read) comes from the forward strand (SF)`
{: .hands_on}

The output of **Salmon** {% icon tool %} includes a table of RNA features, estimated counts, transcripts per million, and feature length. We will be using the TPM quantification column as input to **DESeq2** {% icon tool %} in the next section.

## siRNA differential abundance testing

[**DESeq2** {% icon tool %}](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) is a great tool for differential expression analysis, but we also employ it here for estimation of abundance of reads targeting each of our RNA features. As input, **DESeq2** {% icon tool %} can take *t*ranscripts *p*er *m*illion (TPM) counts produced by **Salmon** {% icon tool %} for each feature. TPMs are estimates of the relative abundance of a given transcript in units, and analogously can be used here to quantify siRNAs that align either sense or antisense to transcripts.

> ### {% icon hands_on %} Hands-on: siRNA differential abundance testing
>
> 1. **DESeq2** {% icon tool %} to test for differential abundance of *antisense* siRNAs at mRNA and TE features using the following parameters:
>    - In *"Factor"*:
>       - In "1: Factor"
>           - *"Specify a factor name"*: `RNAi`
>           - In *"Factor level"*:
>               - In *"1: Factor level"*:
>                   - *"Specify a factor level"*: `Symplekin`
>                   - {% icon param-collection %} *"Counts file(s)"*: the Symplekin RNAi siRNA counts from reverse strand (SR) collection
>               - In *"2: Factor level"*:
>                   - *"Specify a factor level"*: `RNAi`
>                   - {% icon param-collection %} *"Counts file(s)"*: the Symplekin RNAi siRNA counts from reverse strand (SR)
>    - *"Choice of Input data"*: `TPM values (e.g. from sailfish and salmon)`
>       - *"Program used to generate TPMs"*: `Salmon`
>       - *"Gene mapping format"*: `Transcript-ID and Gene-ID mapping file`
>         - {% icon param-file %} *"Tabular file with Transcript - Gene mapping"*: `reference tx2g (transcript to gene) file` (available through zenodo link)
>    - *"Output normalized counts table"*: `Yes`
>
> 1. **DESeq2** {% icon tool %} for the *sense* siRNAs at mRNA and TE features:
>    - In *"Factor"*:
>       - In "1: Factor"
>           - *"Specify a factor name"*: `RNAi`
>           - In *"Factor level"*:
>               - In *"1: Factor level"*:
>                   - *"Specify a factor level"*: `Symplekin`
>                   - {% icon param-collection %} *"Counts file(s)"*: the Symplekin RNAi siRNA counts from forward strand (SF)
>               - In *"2: Factor level"*:
>                   - *"Specify a factor level"*: `RNAi`
>                   - {% icon param-collection %} *"Counts file(s)"*: the Symplekin RNAi siRNA counts from forward strand (SF)
>    - *"Choice of Input data"*: `TPM values (e.g. from sailfish and salmon)`
>       - *"Program used to generate TPMs"*: `Salmon`
>       - *"Gene mapping format"*: `Transcript-ID and Gene-ID mapping file`
>         - {% icon param-file %} *"Tabular file with Transcript - Gene mapping"*: `reference tx2g (transcript to gene) file` (available through zenodo link)
>    - *"Output normalized counts table"*: `Yes`
>
> 1. **Filter** {% icon tool %} to extract features with a significantly different *antisense* siRNA abundance (adjusted *p*-value less than 0.05).
>    - {% icon param-file %} *"Filter"*: Select the `DESeq2` result file from testing *antisense* (SR) siRNA abundances
>    - *"With following condition"*: `c7<0.05`
>
> 1. **Filter** {% icon tool %} for the *sense* siRNAs at mRNA and TE features:
>    - {% icon param-file %} *"Filter"*: Select the `DESeq2` result file from testing *sense* (SF) siRNA abundances
>    - *"With following condition"*: `c7<0.05`
{: .hands_on}

> ### {% icon question %} Question
>
> How many features have a significant difference in *antisense* and *sense* siRNA abundance in the *Symplekin* RNAi condition?
>
> > ### {% icon solution %} Solution
> > There are 87 *antisense* and 15 *sense* features that have significanlty different siRNA abundances.
> {: .solution }
{: .question}

For more information about `DESeq2` and its outputs, have a look at [`DESeq2` documentation](https://www.bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf).

## Visualization

Coming soon!

## Conclusion
{:.no_toc}

Analysis of small RNAs is a complicated and intricate process due to the diversity in characteristics, functionality, and nuances of small RNA subclasses. The goal of this tutorial is to introduce you to a common small RNA workflow that specifically identifies changes in endogenous siRNA abundances that target protein-coding mRNAs and transposable elements. The steps presented here can be rearranged and modified based on small RNA features of specific systems and the needs of the user.

## Small RNA analysis (general) pipeline

![sRNA tutorial workflow](../../images/sRNA/schematic_for_sRNAseq_tutorial.png)
