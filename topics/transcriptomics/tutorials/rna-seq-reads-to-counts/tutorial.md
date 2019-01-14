---
layout: tutorial_hands_on
title: RNA-Seq reads to counts
zenodo_link: "https://figshare.com/s/f5d63d8c265a05618137"
enable: "false"
tags:
  - collections
  - mouse
  - QC
questions:
- How to convert RNA-seq reads into counts?
- How to perform quality control of RNA-seq reads?
- How to do this analysis efficiently in Galaxy?
objectives:
- Learn how RNA-seq reads are converted into counts
- Understand QC steps that can be performed on RNA-seq reads
- Generate interactive reports to summarise QC information with MultiQC
- Use the Galaxy Rule-based Uploader to import FASTQs from URLs
- Make use of Galaxy Collections for a tidy analysis
- Create a Galaxy Workflow that converts RNA-seq reads into counts
time_estimation: '3h'
key_points:
- In RNA-seq, reads (FASTQs) are mapped to a reference genome with a spliced aligner (e.g HISAT2, STAR)
- The aligned reads (BAMs) can then be converted to counts (e.g featureCounts, HTSeq)
- Many QC steps can be performed to help check the quality of the data
- MultiQC can be used to create a nice summary report of QC information
- The Galaxy Rule-based Uploader, Collections and Workflows can help make analysis more efficient and easier
contributors:
- mblue9
- bphipson
- hdashnow

---


# Introduction
{:.no_toc}

Measuring gene expression on a genome-wide scale has become common practice over the last two decades or so, with microarrays predominantly used pre-2008. With the advent of next generation sequencing technology in 2008, an increasing number of scientists use this technology to measure and understand changes in gene expression in often complex systems. As sequencing costs have decreased, using RNA-Seq to simultaneously measure the expression of tens of thousands of genes for multiple samples has never been easier. The cost of these experiments has now moved from generating the data to storing and analysing it.

There are many steps involved in analysing an RNA-Seq experiment. he analysis begins with sequencing reads (FASTQ files). These are usually aligned to a reference genome, if available. Then the number of reads mapped to each gene can be counted. This results in a table of counts, which is what we perform statistical analyses on to determine differentially expressed genes and pathways. The purpose of this tutorial is to demonstrate how to do read alignment and counting, prior to performing differential expression. Differential expression analysis with limma-voom is covered in an accompanying tutorial [RNA-seq counts to genes]({{ site.baseurl }}/topics/transcriptomics/tutorials/rna-seq-counts-to-genes/tutorial.html). The tutorial here shows how to start from FASTQ data and perform the mapping and counting steps, along with associated Quality Control.

**Mouse mammary gland dataset**

The data for this tutorial comes from a Nature Cell Biology paper, [EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival](https://www.ncbi.nlm.nih.gov/pubmed/25730472)), Fu et al. 2015. Both the raw data (sequence reads) and processed data (counts) can be downloaded from Gene Expression Omnibus database (GEO) under accession number [GSE60450](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450).

This study examined the expression profiles of basal and luminal cells in the mammary gland of virgin, pregnant and lactating mice. Six groups are present, with one for each combination of cell type and mouse status. Note that two biological replicates are used here, two independent sorts of cells from the mammary glands of virgin, pregnant or lactating mice, however three replicates is usually recommended as a minimum requirement for RNA-seq.

This is a Galaxy tutorial based on material from the [COMBINE R RNAseq workshop](http://combine-australia.github.io/RNAseq-R/07-rnaseq-day2.html), first taught [here](http://combine-australia.github.io/2016-05-11-RNAseq/).

![Tutorial Dataset](../../images/rna-seq-reads-to-counts/mouse_exp.png "Tutorial Dataset")


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


{% include snippets/warning_results_may_vary.md %}

# Preparing the reads

## Import data from URLs

Read sequences are usually stored in compressed (gzipped) FASTQ files. Before the differential expression analysis can proceed, these reads must be aligned to the reference genome and counted into annotated genes. Mapping reads to the genome is a very important task, and many different aligners are available, such as HISAT2 ([Kim, Langmead, and Salzberg, 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4655817/)), STAR ([Dobin et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/)) and Subread ([Liao, Smyth, and Shi 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3664803/)]). Most mapping tasks require larger computers than an average laptop, so usually read mapping is done on a server in a linux-like environment, requiring some programming knowledge. However, Galaxy enables you to do this mapping without needing to know programming and if you don't have access to a server you can try to use one of the publically available Galaxies e.g. [usegalaxy.org](https://usegalaxy.org), [usegalaxy.eu](https://usegalaxy.eu), [usegalaxy.org.au](https://usegalaxy.org.au/).

If you are sequencing your own data, the sequencing facility will almost always provide compressed FASTQ files which you can upload into Galaxy. If your FASTQs are provided through Galaxy's Shared Data, you can easily import them into a history. For publicly available sequence data, such as from GEO/SRA/ENA, Galaxy's Rule-based Uploader can be used to import the files from URLs, saving on the need to download to your computer and upload into Galaxy.

The raw reads used in this tutorial were obtained from SRA from the link given in GEO for the the mouse mammary gland dataset (Fu et al. 2015) (e.g `ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP045%2FSRP045534`). For the purpose of this tutorial we are going to be working with a small part of the FASTQ files. We are only going to be mapping 1000 reads from each sample to enable running through all the steps quickly. If working with your own data you would use the full data and some results for the full mouse dataset will be shown for comparison. The small FASTQ files are available in [Figshare](https://figshare.com/s/f5d63d8c265a05618137) and the links to the FASTQ files are provided below. We are going to import the files into a Collection. Using Galaxy Collections helps keep the datasets organised and tidy in the history. Collections also make it easier to maintain the sample names through tools and workflows. If you are not familiar with collections, see the [Galaxy Collections tutorial]({{ site.baseurl }}/topics/galaxy-data-manipulation/tutorials/collections/tutorial.html).

The sample information (sample ID, Group) and link to the FASTQ file (URL) are in the grey box below. To generate a file like this to import data from SRA/ENA see the [Galaxy Rule-based Uploader tutorial]({{ site.baseurl }}/topics/galaxy-data-manipulation/tutorials/upload-rules/tutorial.html).

```
SampleID	Group	URL
MCL1-DL	basallactate	https://ndownloader.figshare.com/files/5053573?private_link=f5d63d8c265a05618137
MCL1-DK	basallactate	https://ndownloader.figshare.com/files/5053570?private_link=f5d63d8c265a05618137
MCL1-DJ	basalpregnant	https://ndownloader.figshare.com/files/5053567?private_link=f5d63d8c265a05618137
MCL1-DI	basalpregnant	https://ndownloader.figshare.com/files/5053564?private_link=f5d63d8c265a05618137
MCL1-DH	basalvirgin	https://ndownloader.figshare.com/files/5053561?private_link=f5d63d8c265a05618137
MCL1-DG	basalvirgin	https://ndownloader.figshare.com/files/5053558?private_link=f5d63d8c265a05618137
MCL1-LF	luminalvirgin	https://ndownloader.figshare.com/files/5053555?private_link=f5d63d8c265a05618137
MCL1-LE	luminalvirgin	https://ndownloader.figshare.com/files/5053588?private_link=f5d63d8c265a05618137
MCL1-LD	luminalpregnant	https://ndownloader.figshare.com/files/5053579?private_link=f5d63d8c265a05618137
MCL1-LC	luminalpregnant	https://ndownloader.figshare.com/files/5053582?private_link=f5d63d8c265a05618137
MCL1-LB	luminalvirgin	https://ndownloader.figshare.com/files/5053585?private_link=f5d63d8c265a05618137
MCL1-LA	luminalvirgin	https://ndownloader.figshare.com/files/5053552?private_link=f5d63d8c265a05618137
```

In order to get these files into Galaxy, we will want to do a few things:

* Strip the *header* out of the sample information (it doesn’t contain a URL Galaxy can download).
* Define the file **Identifier** column (`SampleID`).
* Define the **URL** column (`URL`) (this is the location Galaxy can download the data from).
* Tell Galaxy to treat these files as `fastqsanger.gz` files
* Tell Galaxy to treat these files as `mm10` files

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial e.g. `RNA-seq reads to counts`
> 2. Import the files from Figshare using Galaxy's Rule-based Uploader.
>    - Open the Galaxy Upload Manager
>    - Click the tab **Rule-based**
>        - *"Upload data as"*: `Collection(s)`
>        - *"Load tabular data from"*: `Pasted Table`
>    - Paste the table from the grey box above. *(You should now see below)*
>    - Click **Build**
>
>       ![Rule-based Uploader](../../images/rna-seq-reads-to-counts/rule_uploader.png "Rule-based Uploader")
>
>    - In the `rules editor` that pops up:
>
>        - **Remove the header**. From the **Filter** menu select `First or Last N Rows`
>            - *"Filter which rows?"*: `first`
>            - *"Filter how many rows?"*: `1`
>            - Click `Apply`
>
>        - **Define the Identifier and URL columns**. From the **Rules** menu select `Add / Modify Column Definitions`
>            - Click `Add Definition` button and select List Identifier(s)
>                - *"List Identifier(s)"*: `A`
>            - Click `Add Definition` button again and select URL instead
>                - *"URL"*: `C`
>            - Click `Apply`, and you should see your new column definitions listed
>
>        - **Specify the file type**. Select *"Type"*: `fastqsanger.gz`
>        - **Specify the genome**. Select *"Genome"*: `mm10`
>        - **Name the collection**. For *"Name"* enter: `fastqs` *(You should now see below)*
>        - Click `Upload`
>
>        ![Rules Editor](../../images/rna-seq-reads-to-counts/rule_editor.png "Rules Editor")
>
>        You should see a collection (list) called `fastqs` in your history containing all 12 FASTQ files, like below.
>
>        ![Collection](../../images/rna-seq-reads-to-counts/collection.png "Collection of FASTQs")
>
{: .hands_on}

If your data is not accessible by URL, for example, if your FASTQ files are located on your laptop, you can upload into a collection as below.

> ### {% icon tip %} Tip: Upload local files into a collection
>
> - Open the Galaxy Upload Manager
> - Click the tab **Collection**
> - Click **Choose Local Files** and locate the files you want to upload
>     - *"Collection Type"*: `List`
>     - *"File Type"*: `fastqsanger.gz`
>     - *"Genome"*: `mm10`
> - In the pop up that appears:
>     - *"Name"*: `fastqs`
>     - Click `Create list`
{: .tip}

If your FASTQ files are located in Shared Data, you can import them into your history as a collection as below.

> ### {% icon tip %} Tip: Import files from Shared Data into a collection
>
> - In the Menu at the top go to Shared Data > Data Libraries
> - Locate your FASTQ files
> - Tick the checkboxes to select the files
> - From the **To History** menu select `as a Collection`
> - In the pop up that appears:
>     - *"Which datasets?"*: `current selection`
>     - *"Collection type"*: `List`
>     - *"Select history"*: `select your History`
>     - Click `Continue`
> - In the pop up that appears:
>     - *"Name"*: `fastqs`
>     - Click `Create list`
{: .tip}

Take a look at one of the FASTQ files to see what it contains.

> ### {% icon tip %} Tip: FASTQ format
> If you are not familiar with FASTQ format, see the [Quality Control tutorial]({{ site.baseurl }}/topics/sequence-analysis/tutorials/quality-control/tutorial.html)
{: .tip}

> ### {% icon hands_on %} Hands-on: Take a look at FASTQ format
>
> 1. Click on the collection name (`fastqs`)
> 2. Click on the name of one of the collection elements (e.g `MCL1-DL`). Check that the data type is `fastqsanger.gz` and the database `mm10`.
>
>    ![Collection Element](../../images/rna-seq-reads-to-counts/collection_element.png "Collection Element"){: width="50px" height="50px"}
>
> 3. Click on the {% icon galaxy-eye %} (eye) icon of one of the FASTQ files to have a look at what it contains
{: .hands_on}


## Check raw reads

During sequencing, errors are introduced, such as incorrect nucleotides being called. These are due to the technical limitations of each sequencing platform. Sequencing errors might bias the analysis and can lead to a misinterpretation of the data. Every base sequence gets a quality score from the sequencer and this information is present in the FASTQ file. A quality score of 30 corresponds to a 1 in 1000 chance of an incorrect base call (a quality score of 10 is a 1 in 10 chance of an incorrect base call). To look at the overall distribution of quality scores across the reads, we can use FastQC.

Sequence quality control is therefore an essential first step in your analysis. We will use similar tools as described in the ["Quality control" tutorial]({{site.baseurl}}/topics/sequence-analysis): [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html).

> ### {% icon hands_on %} Hands-on: Check raw reads with **FastQC**
>
> 1. **FastQC** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"Short read data from your current history"*: `fastqs` (Input dataset collection)
>
>    > ![FastQC Collection](../../images/rna-seq-reads-to-counts/run_fastqc.png "Run FastQC on collection"){: width="50px"}
>
> 2. Inspect the `Webpage` output of **FastQC** {% icon tool %} for the `MCL1-DL` sample by clicking on the {% icon galaxy-eye %} (eye) icon
{: .hands_on}

> ### {% icon question %} Questions
> 1. What is the read length?
> 2. What base quality score encoding is used?
>
> > ### {% icon solution %} Solution
> >
> > 1. The read length is 100 bp.
> > 2. Sanger quality score encoding is used.
> > This information can be seen at the top of the FastQC Webpage as below.
> >
> > ![FastQC Webpage](../../images/rna-seq-reads-to-counts/fastqc_webpage.png "FastQC Webpage"){: width="50px"}
> >
> {: .solution}
>
{: .question}

The FastQC report contains a lot of information and we can look at the report for each sample. However, that is quite a few reports, 12 for this dataset. If you had more samples it could be a lot more. Luckily, there is a very useful tool called MultiQC that can summarise QC information for multiple samples into a single report.

> ### {% icon hands_on %} Hands-on: Aggregate FastQC reports with **MultiQC**
>
> 1. **MultiQC** {% icon tool %} with the following parameters to aggregate the FastQC reports
>      - In *"Results"*
>        - {% icon param-select %}*"Which tool was used generate logs?"*: `FastQC`
>        - In *"FastQC output"*
>           - {% icon param-select %} *"Type of FastQC output?"*: `Raw data`
>           - {% icon param-collection %} *"FastQC output"*: `RawData` files (output of **FastQC** {% icon tool %} on trimmed reads)
>           ![MultiQC](../../images/rna-seq-reads-to-counts/multiqc.png "Run MultiQC on collection"){: width="50px"}
> 2. Inspect the `Webpage` output from MultiQC
{: .hands_on}

Note that these are the results for just 1000 reads. The FastQC results for the full dataset are shown below. These 1000 reads are the first reads from the FASTQ files, and the first reads usually originate from the flowcell edges, so we can expect that they may have lower quality and the patterns may be a bit different from the distribution in the full dataset.

Here, most of the plots in the small FASTQs look similar to the full dataset. However, in the small FASTQs, we see less duplication, some Ns in the reads and some overrepresented sequences.

![General Statistics](../../images/rna-seq-reads-to-counts/fastqc_table.png "General Statistics")
![Sequence Counts](../../images/rna-seq-reads-to-counts/fastqc_sequence_counts_plot.png "Sequence Counts")
![Sequence Quality](../../images/rna-seq-reads-to-counts/fastqc_per_base_sequence_quality_plot.png "Sequence Quality")
![Per Sequence Quality Scores](../../images/rna-seq-reads-to-counts/fastqc_per_sequence_quality_scores_plot.png "Per Sequence Quality Scores")
![Per Sequence GC Content](../../images/rna-seq-reads-to-counts/fastqc_per_sequence_gc_content_plot.png "Per Sequence GC Content")
![Per base N content](../../images/rna-seq-reads-to-counts/fastqc_per_base_n_content_plot.png "Per base N content")
![Sequence Duplication Levels](../../images/rna-seq-reads-to-counts/fastqc_sequence_duplication_levels_plot.png "Sequence Duplication Levels")
![Adapter Content](../../images/rna-seq-reads-to-counts/fastqc_adapter_content_plot.png "Adapter Content")

See the [Quality Control tutorial]({{ site.baseurl }}/topics/sequence-analysis/tutorials/quality-control/tutorial.html) for more information on FastQC plots.

> ### {% icon question %} Questions
>
> What do you think of the overall quality of the sequences?
>
> > ### {% icon solution %} Solution
> >
> > Overall, the samples look pretty good. The main things to note here are:
> > * The **base quality** is high in all samples.
> > * Some Illumina **adapter** has been detected.
> > * Some duplication in RNA-seq can be normal due to the presence of highly expressed genes. However, for some reason `MCL1-LE` and `MCL1-LF` have higher numbers of duplicates detected than the other samples.
> >
> {: .solution}
{: .question}

We will use Cutadapt to trim the reads to remove the Illumina adapter and any low quality bases at the ends (quality score < 20). We will also output the Cutadapt report for summarising with MultiQC.

The Cutadapt tool Help section provides the sequence we can use to trim this standard Illumina adapter `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC`, as given on the [Cutadapt website](https://cutadapt.readthedocs.io/en/stable/guide.html#illumina-truseq). Other Illumina adapter sequences (e.g. Nextera) can be found at the [Illumina website](http://sapac.support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html). We will discard any sequences that are too short (< 20bp) after trimming.

## Trim reads

> ### {% icon hands_on %} Hands-on: Trim reads with **Cutadapt**
>
> 1. **Cutadapt** {% icon tool %} with the following parameters:
>    - {% icon param-select %} *"Single-end or Paired-end reads?"*: `Single-end`
>        - {% icon param-collection %} *"FASTQ/A file"*: `fastqs` (Input dataset collection)
>        - In *"Read 1 Options"*:
>            - In *"3' (End) Adapters"*:
>                - Click on *"Insert 3' (End) Adapters"*:
>                - In *"1: 3' (End) Adapters"*:
>                    - {% icon param-select %} *"Source"*: `Enter custom sequence`
>                        - {% icon param-text %} *"Enter custom 3' adapter name (Optional)"*: `Illumina`
>                        - {% icon param-text %} *"Enter custom 3' adapter sequence"*: `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC`
>    - In *"Filter Options"*:
>        - {% icon param-text %} *"Minimum length"*: `20`
>    - In *"Read Modification Options"*:
>        - {% icon param-text %} *"Quality cutoff"*: `20`
>    - In *"Output Options"*:
>        - {% icon param-check %} *"Report"*: `Yes`
>
{: .hands_on}

We can take a look at the reads again now that they've been trimmed.

## Check trimmed reads

> ### {% icon hands_on %} Hands-on: Check trimmed reads with **FastQC**
>
> 1. **FastQC** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"Short read data from your current history"*: `RawData` (output of **Cutadapt** {% icon tool %})
> 2. **MultiQC** {% icon tool %} with the following parameters to aggregate the FastQC reports
>    - In *"Results"*
>        - {% icon param-select %} *"Which tool was used generate logs?"*: `FastQC`
>        - In *"FastQC output"*
>           - {% icon param-select %} *"Type of FastQC output?"*: `Raw data`
>           - {% icon param-collection %} *"FastQC output"*: `RawData` files (output of **FastQC** {% icon tool %})
> 3. Inspect the `Webpage` output from MultiQC
{: .hands_on}

The MultiQC plot below shows the result from the full dataset for comparison.

![Adapter Content post-trimming](../../images/rna-seq-reads-to-counts/post_cutadapt_adapter_content.png "Adapter Content post-trimming")
![Sequence Length post-trimming](../../images/rna-seq-reads-to-counts/post_cutadapt_fastqc_sequence_length_distribution_plot.png "Sequence Length post-trimming")

After trimming we can see that:

* No adapter is detected now.
* The reads are no longer all the same length, we now have sequences of different lengths detected.

# Mapping

Now that we have prepared our reads, we can align the reads for our 12 samples. There is an existing reference genome for mouse and we will map the reads to that. The current version of the mouse reference genome is `mm10/GRCm38`. Here we will use [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml), the descendent of widely-used TopHat, to align the reads, but alternative mappers could be used, such as STAR. See the [RNA-seq ref-based tutorial]({{ site.baseurl }}/topics/transcriptomics/tutorials/ref-based/tutorial.html#mapping) for more information on RNA-seq mappers. There are often numerous mapping parameters that we can specify, but usually the default mapping parameters are fine. However, library type (paired-end vs single-end) and library strandness (stranded vs unstranded) require some different settings when mapping and counting, so they are two important pieces of information to know about samples. The mouse data comprises unstranded, single-end reads so we will specify that where necessary. If we had paired-end data or stranded data, we would specify that. HISAT2 can output a mapping summary file that tells what proportion of reads mapped to the reference genome. Summary files for multiple samples can be summarised with MultiQC. As we’re only using a subset of 1000 reads per sample, aligning should just take a minute or so. To run the full samples from this dataset would take longer.

## Map reads to reference genome

> ### {% icon hands_on %} Hands-on: Map reads to reference with **HISAT2**
>
> **HISAT2** {% icon tool %} with the following parameters:
>    - {% icon param-select %} *"Source for the reference genome"*: `Use a built-in genome`
>        - {% icon param-select %} *"Select a reference genome"*: `mm10`
>    - {% icon param-select %} *"Is this a single or paired library?"*: `Single-end`
>        - {% icon param-collection %} *"FASTA/Q file"*: `Read 1 Output` (output of **Cutadapt** {% icon tool %})
>    - In *"Summary Options"*:
>        - {% icon param-check %} *"Output alignment summary in a more machine-friendly style."*: `Yes`
>        - {% icon param-check %} *"Print alignment summary to a file."*: `Yes`
{: .hands_on}

> ### {% icon tip %} Tip: Settings for Paired-end or Stranded reads
>
> - If you have **paired-end** reads
>     - Select *"Is this a single or paired library"* `Paired-end` or `Paired-end Dataset Collection` or `Paired-end data from single interleaved dataset`
> - If you have **stranded** reads
>     - Select *"Specify strand information"*: `Forward (FR)` or `Reverse (RF)`
{: .tip}

**HISAT2** generates a BAM file with mapped reads.

> ### {% icon hands_on %} Hands-on: Take a look at BAM format
>
> 1. Click on the collection name (`HISAT2 on collection N: aligned reads (BAM)`)
> 2. Click on the {% icon galaxy-eye %} (eye) icon of one of the BAM files to have a look at what it contains.
>
>    > ### {% icon question %} Questions
>    >
>    > What is BAM format?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > See the answer in the [RNA-seq ref-based tutorial]({{ site.baseurl }}/topics/transcriptomics/tutorials/ref-based/tutorial.html#pretreatments)
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
{: .hands_on}

## Check mapped reads

An important metric to check is the percentage of reads mapped to the reference genome. A low percentage can indicate issues with the data or analysis. We can use MultiQC again to summarise the QC information.

> ### {% icon hands_on %} Hands-on: Aggregate the HISAT2 summary files with **MultiQC**
>
> 1. **MultiQC** {% icon tool %} with the following parameters to aggregate the HISAT2 summary files
>    - In *"Results"*
>        - {% icon param-select %} *"Which tool was used generate logs?"*: `HISAT2`
>        - {% icon param-collection %} *"Output of HISAT2"*: `Mapping summary` (output of **HISAT2** {% icon tool %})
> 2. Inspect the `Webpage` output from MultiQC
{: .hands_on}

The MultiQC plot below shows the result from the full dataset for comparison.

![HISAT2 mapping](../../images/rna-seq-reads-to-counts/hisat2_se_plot.png "HISAT2 mapping")

Over 90% of reads have mapped in all samples, a good mapping rate. And the vast majority of reads have mapped uniquely, they haven't mapped to multiple locations in the reference genome.

It is also good practice to visualise the read alignments in the BAM file, for example using IGV, see the [RNA-seq ref-based tutorial]({{ site.baseurl }}/topics/transcriptomics/tutorials/ref-based/tutorial.html#pretreatments).

> ### {% icon tip %} Tip: Downloading a collection
>
> To download a collection of datasets (e.g. the collection of BAM files) click on the floppy icon as shown below. This will download a tar file containing all the datasets in the collection. Note that for BAM files the .bai indexes (required for IGV) will be included automatically in the download.
>
> ![Downloading collection](../../images/rna-seq-reads-to-counts/collection_download.png "Downloading a collection")
{: .tip}

## Check strandness

As far as we know this data is unstranded, but as a sanity check you can check the strandness. You can use RSeQC Infer Experiment tool to "guess" the strandness, as explained in the [RNA-seq ref-based tutorial]({{ site.baseurl }}/topics/transcriptomics/tutorials/ref-based/tutorial.html). This is done through comparing the "strandness of reads" with the "strandness of transcripts". For this tool, and many of the other RSeQC tools, a reference bed file of genes (`reference genes`) is required. RSeQC provides some reference BED files for model organisms. You can import the RSeQC mm10 RefSeq BED file from the link `https://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/mm10_RefSeq.bed.gz/download` (and rename to `reference genes`) or import a file from Shared data if provided. Alternatively, you can provide your own BED file of reference genes, for example from UCSC (see the [Peaks to Genes tutorial]({{ site.baseurl }}/topics/introduction/tutorials/galaxy-intro-peaks2genes/tutorial.html). Or the **Convert GTF to BED12** tool can be used to convert a GTF into a BED file.

{% include snippets/import_via_link.md %}

> ### {% icon hands_on %} Hands-on: Check strandness with **Infer Experiment**
>
> 1. **Infer Experiment** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"Input .bam file"*: `aligned reads (BAM)` (output of **HISAT2** {% icon tool %})
>    - {% icon param-file %} *"Reference gene model"*: `reference genes` (Reference BED file)
> 2. **MultiQC** {% icon tool %} with the following parameters:
>       - In *"1: Results"*:
>           - {% icon param-select %} *"Which tool was used generate logs?"*: `RSeQC`
>               - {% icon param-select %} *"Type of RSeQC output?"*: `infer_experiment`
>                   - {% icon param-collection %} *"RSeQC infer_experiment output"*: `Infer Experiment output` (output of **Infer Experiment** {% icon tool %})
> 3. Inspect the `Webpage` output from MultiQC
{: .hands_on}

The MultiQC plot below shows the result from the full dataset for comparison.

![Infer Experiment](../../images/rna-seq-reads-to-counts/rseqc_infer_experiment_plot.png "Infer Experiment")

> ### {% icon question %} Questions
>
> Do you think the data is stranded or unstranded?
>
> > ### {% icon solution %} Solution
> >
> > It is unstranded as approximately equal numbers of reads have aligned to the sense and antisense strands.
> >
> {: .solution}
>
{: .question}

## Check duplicate reads

Duplicate reads are usually kept in RNA-seq differential expression analysis as they can come from highly-expressed genes but it is still a good metric to check. A high percentage of duplicates can indicate a problem with the sample, for example, PCR amplification of a low complexity library (not many transcripts) due to not enough RNA used as input. FastQC gives us an idea of duplicates in the reads before mapping (note that it just takes a sample of the data). We can assess the numbers of duplicates in all mapped reads using the **Picard MarkDuplicates** tool. Picard considers duplicates to be reads that map to the same location, based on the start position of where the read maps.

> ### {% icon hands_on %} Hands-on: Check duplicate reads with **MarkDuplicates**
>
> 1. **MarkDuplicates** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"Select SAM/BAM dataset or dataset collection"*: `aligned reads (BAM)` (output of **HISAT2** {% icon tool %})
> 2. **MultiQC* {% icon tool %} with the following parameters:
>       - In *"1: Results"*:
>           - {% icon param-select %} *"Which tool was used generate logs?"*: `Picard`
>               - {% icon param-select %} *"Type of Picard output?"*: `Markdups`
>                   - {% icon param-collection %} *"Picard output"*: `MarkDuplicate metrics` (output of **MarkDuplicates** {% icon tool %})
> 3. Inspect the `Webpage` output from MultiQC
{: .hands_on}

The MultiQC plot below shows the result from the full dataset for comparison.

![MarkDups metrics](../../images/rna-seq-reads-to-counts/picard_deduplication.png "MarkDups metrics")

> ### {% icon question %} Questions
>
> Which two samples have the most duplicates detected?
>
> > ### {% icon solution %} Solution
> >
> > `MCL1-LE` and `MCL1-LF` have the highest number of duplicates in mapped reads compared to the other samples, similar to what we saw in the raw reads with FastQC.
> >
> {: .solution}
>
{: .question}

## Check number of reads for each chromosome

You can check the numbers of reads mapped to each chromosome with the **Samtools IdxStats** tool. This can help assess the sample quality, for example, if there is an excess of mitochondrial contamination. It could also help to check the sex of the sample through the numbers of reads mapping to X/Y or to see if any chromosomes have highly expressed genes.

> ### {% icon hands_on %} Hands-on: Count reads mapping to each chromosome with **IdxStats**
>
> 1. **IdxStats** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"BAM file"*: `aligned reads (BAM)` (output of **HISAT2** {% icon tool %})
> 2. **MultiQC* {% icon tool %} with the following parameters:
>       - In *"1: Results"*:
>           - {% icon param-select %} *"Which tool was used generate logs?"*: `Samtools`
>               - {% icon param-select %} *"Type of Samtools output?"*: `idxstats`
>                   - {% icon param-collection %} *"Samtools idxstats output"*: `IdxStats output` (output of **IdxStats** {% icon tool %})
> 3. Inspect the `Webpage` output from MultiQC
{: .hands_on}

The MultiQC plot below shows the result from the full dataset for comparison.

![IdxStats chromsome mappings](../../images/rna-seq-reads-to-counts/samtools-idxstats-mapped-reads-plot.png "IdxStats Chromosome Mappings")
![IdxStats XY ](../../images/rna-seq-reads-to-counts/samtools-idxstats-xy-plot.png "IdxStats X/Y Mappings")

> ### {% icon question %} Questions
>
> 1. What do you think of the chromosome mappings?
> 2. Are the samples male or female? *(If a sample is not in the XY plot it means no reads mapped to Y)*
>
> > ### {% icon solution %} Solution
> >
> > 1. Some of the samples have very high mapping on chromosome 5. What is going on there?
> > 2. The samples appear to be all female as there are few reads mapping to the Y chromosome. As this is a experiment studying virgin, pregnant and lactating mice if we saw large numbers of reads mapping to the Y chromosome in a sample it would be unexpected and a probable cause for concern.
> >
> {: .solution}
>
{: .question}

## Check coverage across gene bodies (5'-3')

The coverage of reads along gene bodies can be assessed to check if there is any bias in coverage. For example, a bias towards the 3' end of genes could indicate degradation of the RNA. Alternatively, a 3' bias could indicate that the data is from a 3' assay (e.g. oligodT-primed, 3'RNA-seq). You can use the RSeQC **Gene Body Coverage (BAM)** tool to assess gene body coverage in the BAM files.

> ### {% icon hands_on %} Hands-on: Check coverage of genes with **Gene Body Coverage (BAM)**
>
> 1. **Gene Body Coverage (BAM)** {% icon tool %} with the following parameters:
>    - *"Run each sample separately, or combine mutiple samples into one plot"*: `Run each sample separately`
>        - {% icon param-collection %} *"Input .bam file"*: `aligned reads (BAM)` (output of **HISAT2** {% icon tool %})
>    - {% icon param-file %} *"Reference gene model"*: `reference genes` (Input dataset)
> 2. **MultiQC* {% icon tool %} with the following parameters:
>       - In *"1: Results"*:
>           - {% icon param-select %} *"Which tool was used generate logs?"*: `RSeQC`
>               - {% icon param-select %} *"Type of RSeQC output?"*: `gene_body_coverage`
>                   - {% icon param-collection %} *"RSeQC gene_body_coverage output"*: `Gene Body Coverage (BAM) (text)` (output of **Gene Body Coverage (
BAM)** {% icon tool %})
> 3. Inspect the `Webpage` output from MultiQC
>
{: .hands_on}

The MultiQC plot below shows the result from the full dataset for comparison.

![Gene Body Coverage](../../images/rna-seq-reads-to-counts/rseqc_gene_body_coverage_plot.png "Gene Body Coverage")

The plot below from the RSeQC website shows what samples with 3'biased coverage would look like.
![Gene Body Coverage comparison](../../images/rna-seq-reads-to-counts/genebodycoverage.png "Gene Body Coverage comparison")

> ### {% icon question %} Questions
>
> What do you think of the coverage across gene bodies in these samples?
>
> > ### {% icon solution %} Solution
> >
> > It looks good. This plot looks a bit noisy in the small FASTQs but it still shows there's pretty even coverage from 5' to 3' ends with no obvious bias in all the samples.
> >
> {: .solution}
>
{: .question}


## Check distribution of reads across features (exons, introns, intergenic..)

We can also check the distribution of reads across known gene features, such as exons (CDS, 5'UTR, 3'UTR), introns and intergenic regions. In RNA-seq we expect most reads to map to exons rather than introns or intergenic regions. It is also the reads mapped to exons that will be counted so it is good to check what proportions of reads have mapped to those. High numbers of reads mapping to intergenic regions could indicate the presence of DNA contamination.

> ### {% icon hands_on %} Hands-on: Check distribution of reads with **Read Distribution**
>
> 1. **Read Distribution** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"Input .bam/.sam file"*: `aligned reads (BAM)` (output of **HISAT2** {% icon tool %})
>    - {% icon param-file %} *"Reference gene model"*: `reference genes` (Input dataset)
> 2. **MultiQC* {% icon tool %} with the following parameters:
>       - In *"1: Results"*:
>           - {% icon param-select %} *"Which tool was used generate logs?"*: `RSeQC`
>               - {% icon param-select %} *"Type of RSeQC output?"*: `read_distribution`
>                   - {% icon param-collection %} *"RSeQC read_distribution output"*: `Read Distribution output` (output of **Read Distribution** {% icon tool %})
> 3. Inspect the `Webpage` output from MultiQC
>
{: .hands_on}

The MultiQC plot below shows the result from the full dataset for comparison.

![Read Distribution](../../images/rna-seq-reads-to-counts/rseqc_read_distribution_plot.png "Read Distribution")

> ### {% icon question %} Questions
>
> What do you think of the read distribution?
>
> > ### {% icon solution %} Solution
> >
> > It looks good, most of the reads have mapped to exons and not many to introns or intergenic regions. The samples have pretty consistent read distribution, albeit with slightly higher numbers of reads mapping to CDS exons for `MCL1-LC` and `MCL1-LD`, and `MCL1-LE` and `MCL1-LF` have more reads mapping to CDS exons than the other samples.
> >
> {: .solution}
>
{: .question}

Now that we have checked the data and are happy with how it looks, we can proceed to counting.

# Counting

The alignment produces a set of BAM files, where each file contains the read alignments for each sample. In the BAM file, there is a chromosomal location for every read that mapped. Now that we have figured out where each read comes from in the genome, we need to summarise the information across genes or exons. The mapped reads can be counted across mouse genes by using a tool called featureCounts. featureCounts requires gene annotation specifying the genomic start and end position of each exon of each gene. For convenience, featureCounts contains built-in annotation for mouse (`mm10`, `mm9`) and human (`hg38`, `hg19`) genome assemblies, where exon intervals are defined from the NCBI RefSeq annotation of the reference genome. Reads that map to exons of genes are added together to obtain the count for each gene, with some care taken with reads that span exon-exon boundaries. For other species, users will need to read in a data frame in GTF format to define the genes and exons. Users can also specify a custom annotation file in SAF format. See the tool help in Galaxy, which has an example of what an SAF file should like like, or the Rsubread users guide for more information.

> ### {% icon comment %} Comment
>
> In this example we have kept many of the default settings, which are typically optimised to work well under a variety of situations. For example, the default setting for featureCounts is that it only keeps reads that uniquely map to the reference genome. For testing differential expression of genes, this is preferred, as the reads are unambigously assigned to one place in the genome, allowing for easier interpretation of the results. Understanding all the different parameters you can change involves doing a lot of reading about the tool that you are using, and can take a lot of time to understand! We won’t be going into the details of the parameters you can change here, but you can get more information from looking at the tool help.
{: .comment}

## Count reads mapped to genes

> ### {% icon hands_on %} Hands-on: Count reads mapped to genes with **featureCounts**
>
> **featureCounts** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"Alignment file"*: `aligned reads (BAM)` (output of **HISAT2** {% icon tool %})
>    - {% icon param-select %} *"Gene annotation file"*: `featureCounts built-in`
>        - {% icon param-select %} *"Select built-in genome"*: `mm10`
{: .hands_on}

> ### {% icon tip %} Tip: Settings for Paired-end or Stranded reads
>
> - If you have **paired-end** reads
>     - Click *"Options for paired-end reads"*
>         - {% icon param-select %} *"Count fragments instead of reads"*: `Enabled; fragments (or templates) will be counted instead of reads`
> - If you have **stranded** reads
>     - {% icon param-select %} Select *"Specify strand information"*: `Stranded (Forward)` or `Stranded (Reverse)`
{: .tip}

## Check assignments of reads to genes

featureCounts reports the numbers of unassigned reads and the reasons why they are not assigned (eg. ambiguity, multi-mapping, secondary alignment, mapping quality, fragment length, chimera, read duplicate, non-junction and so on), in addition to the number of successfully assigned reads for each library. The statistics of the read mapping are output as a file from featureCounts. MultiQC can be used to assess the numbers of reads assigned to features, genes in this case.

> ### {% icon hands_on %} Hands-on: Assess reads mapped to genes
>
> 1. **MultiQC** {% icon tool %} with the following parameters:
>    - {% icon param-select %} *"Which tool was used generate logs?"*: `featureCounts`
>        - {% icon param-collection %} *"Output of FeatureCounts"*: `featureCounts summary` (output of **featureCounts** {% icon tool %})
> 2. Inspect the `Webpage` output from MultiQC
{: .hands_on}

The MultiQC plot below shows the result from the full dataset for comparison.

![featureCounts assignments](../../images/rna-seq-reads-to-counts/featureCounts_assignment_plot.png "featureCounts assignments")

> ### {% icon question %} Questions
>
> What % reads are assigned to exons?
>
> > ### {% icon solution %} Solution
> >
> > ~60-70% of reads are assigned to exons. This is a fairly typical number for RNA-seq.
> >
> {: .solution}
>
{: .question}

The counts for the samples are output as tabular files. Take a look at one. The numbers in the first column of the counts file represent the Entrez gene identifiers for each gene, while the second column contains the counts for each gene for the sample.

## Create count matrix

The counts files are currently in the format of one file per sample. However, it is often convenient to have a count matrix. A count matrix is a single table containing the counts for all samples, with the genes in rows and the samples in columns. The counts files are all within a collection so we can use the Galaxy **Column Join on Collection** tool to easily create a count matrix from the single counts files.

> ### {% icon hands_on %} Hands-on: Create count matrix with **Column Join on Collection**
>
> **Column Join on Collection** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"Tabular files"*: `featureCounts output` (output of **featureCounts** {% icon tool %})
>    - {% icon param-text %} *"Identifier column"*: `1`
>    - {% icon param-text %} *"Number of Header lines in each item"*: `1`
>    - {% icon param-check %} *"Keep original column header"*: `No`
>
{: .hands_on}

Take a look at the output (the output for the full dataset is shown below).

![Count matrix](../../images/rna-seq-reads-to-counts/count_matrix.png "Count matrix")

Now it is easier to see the counts for a gene across all samples. The accompanying limma-voom tutorial shows how gene information (symbols etc) can be added to a count matrix.

# Generating an interactive QC summary report

So far we have run MultiQC on one step at a time which generates multiple reports for the analysis. However, it would be more convenient to have a single report containing all the QC information for the analysis. This can be easily achieved in Galaxy with the MultiQC tool.

> ### {% icon question %} Questions
>
> Can you create a MultiQC report that summarises all the QC information into a single report?
>
> * Summarise the QC info for **Cutadapt**, **FastQC**, **HISAT2**, **Infer Experiment**, **MarkDups**, **IdxStats**, **Gene Body Coverage**, **Read Distribution**, **featureCounts**.
> * Use just one of the two FastQC outputs e.g. the trimmed FastQC RawData and not the untrimmed aswell, as you currently can't have more than one run of the same tool in the same MultiQC report.
>
> > ### {% icon solution %} Solution
> >
> > MultiQC summary of all programs used (screenshot from the full dataset).
> >
> > ![MultiQC Summary](../../images/rna-seq-reads-to-counts/multiqc_summary.png "MultiQC Summary")
> >
> > **MultiQC** {% icon tool %} with the following parameters:
> >    - In *"Results"*:
> >        - In *"1: Results"*:
> >            - {% icon param-select %} *"Which tool was used generate logs?"*: `Cutadapt/Trim Galore!`
> >                - {% icon param-collection %} *"Output of Cutadapt"*: `Cutadapt Report` (output of **Cutadapt** {% icon tool %})
> >        - {% icon param-repeat %} Click on *"Insert Results"*:
> >        - In *"2: Results"*:
> >            - {% icon param-select %} *"Which tool was used generate logs?"*: `FastQC`
> >                - {% icon param-collection %} *"FastQC output"*: `FastQC RawData` (output of **FastQC** {% icon tool %})
> >        - {% icon param-repeat %} Click on *"Insert Results"*:
> >        - In *"3: Results"*:
> >            - {% icon param-select %} *"Which tool was used generate logs?"*: `HISAT2`
> >                - {% icon param-collection %} *"Output of HISAT2"*: `Mapping summary` (output of **HISAT2** {% icon tool %})
> >        - {% icon param-repeat %} Click on *"Insert Results"*:
> >        - In *"4: Results"*:
> >            - {% icon param-select %} *"Which tool was used generate logs?"*: `RSeQC`
> >                - {% icon param-select %} *"Type of RSeQC output?"*: `infer_experiment`
> >                    - {% icon param-collection %} *"RSeQC infer_experiment output"*: `Infer Experiment output` (output of **Infer Experiment** {% icon tool %})
> >        - {% icon param-repeat %} Click on *"Insert Results"*:  
> >        - In *"5: Results"*:
> >            - {% icon param-select %} *"Which tool was used generate logs?"*: `Picard`
> >                - {% icon param-select %} *"Type of Picard output?"*: `Markdups`
> >                    - {% icon param-collection %} *"Picard output"*: `MarkDuplicate metrics` (output of **MarkDuplicates** {% icon tool %})
> >        - {% icon param-repeat %} Click on *"Insert Results"*:
> >        - In *"6: Results"*:
> >            - {% icon param-select %} *"Which tool was used generate logs?"*: `Samtools`
> >               - {% icon param-select %} *"Type of Samtools output?"*: `idxstats`
> >                    - {% icon param-collection %} *"Samtools idxstats output"*: `IdxStats output` (output of **IdxStats** {% icon tool %})
> >        - {% icon param-repeat %} Click on *"Insert Results"*:
> >        - In *"7: Results"*:
> >            - {% icon param-select %} *"Which tool was used generate logs?"*: `RSeQC`
> >                - {% icon param-select %} *"Type of RSeQC output?"*: `gene_body_coverage`
> >                    - {% icon param-collection %} *"RSeQC gene_body_coverage output"*: `Gene Body Coverage (BAM) (text)` (output of **Gene Body Coverage (
BAM)** {% icon tool %})
> >        - {% icon param-repeat %} Click on *"Insert Results"*:
> >        - In *"8: Results"*:
> >            - {% icon param-select %} *"Which tool was used generate logs?"*: `RSeQC`
> >                - {% icon param-select %} *"Type of RSeQC output?"*: `read_distribution`
> >                    - {% icon param-collection %} *"RSeQC read_distribution output"*: `Read Distribution output` (output of **Read Distribution** {% icon tool %})
> >        - {% icon param-repeat %} Click on *"Insert Results"*:
> >        - In *"9: Results"*:
> >            - {% icon param-select %} *"Which tool was used generate logs?"*: `featureCounts`
> >                - {% icon param-collection %} *"Output of FeatureCounts"*: `featureCounts summary` (output of **featureCounts** {% icon tool %})
> >
> {: .solution}
>
{: .question}


The MultiQC report can be downloaded by clicking on the floppy disk icon on the dataset in the history as shown below.

![MultiQC download](../../images/rna-seq-reads-to-counts/multiqc_download.png "Download MultiQC report")

> ### {% icon question %} Questions
>
> Can you think of any other QCs that could be performed on RNA-seq reads?
>
> > ### {% icon solution %} Solution
> >
> > The reads could be checked for:
> > * Ribosomal contamination
> > * Contamination with other species e.g. bacteria
> > * GC bias of the mapped reads
> > * This is single-end data but paired-end mapped reads could be checked for fragment size (distance between the read pairs).
> >
> {: .solution}
>
{: .question}

# Creating a workflow of the analysis

You can now extract a workflow from the steps that have been performed, as shown in the [Peaks to Genes tutorial]({{ site.baseurl }}/topics/introduction/tutorials/galaxy-intro-peaks2genes/tutorial.html#extracting-workflow). This is a way to help keep a record of the steps performed. It also enables efficient rerunning of a multi-step analysis, such as RNA-seq.

![fastq to counts Workflow](../../images/rna-seq-reads-to-counts/workflow.png "FASTQ to counts Workflow")

# Conclusion
{:.no_toc}

In this tutorial we have seen how reads (FASTQ files) can be converted into counts. We have also seen QC steps that can be performed to help assess the quality of the data. A follow-on tutorial, [RNA-seq counts to genes]({{ site.baseurl }}/topics/transcriptomics/tutorials/rna-seq-counts-to-genes/tutorial.html), shows how to perform differential expression and QC on the counts for this dataset.
