---
layout: tutorial_hands_on

title: Limma-voom (fastqs to counts)
zenodo_link: "https://figshare.com/s/f5d63d8c265a05618137"
questions:
- How to convert RNA-seq reads into counts before performing differential expression?
- How to perform quality control of RNA-seq data?
- How to do this analysis efficiently in Galaxy?
objectives:
- Know how RNA-seq reads are converted into counts
- Understand QC steps that can be performed on RNA-seq data
- Generate interactive reports to summarise QC information with MultiQC
- Use the Galaxy Rule Uploader to import fastqs from URLs
- Make use of Galaxy Collections for a tidy analysis
- Create a Galaxy Workflow that converts RNA-seq reads into counts
time_estimation: '2h'
key_points:
- In RNA-seq, reads (fastqs) are mapped to a reference genome with a spliced aligner (e.g HISAT2, STAR)
- The aligned reads (BAMs) can then be converted to counts (e.g featureCounts, HTSeq)
- Many QC steps can be performed to help check the quality of the data
- MultiQC can be used to create a nice summary report of QC information
- The Galaxy Rule Uploader, Collections and Workflows can help make analysis more efficient and easier
contributors:
- mblue9
- Belinda Phipson
- Harriet Dashnow

---


# Introduction
{:.no_toc}

Measuring gene expression on a genome-wide scale has become common practice over the last two decades or so, with microarrays predominantly used pre-2008. With the advent of next generation sequencing technology in 2008, an increasing number of scientists use this technology to measure and understand changes in gene expression in often complex systems. As sequencing costs have decreased, using RNA-Seq to simultaneously measure the expression of tens of thousands of genes for multiple samples has never been easier. The cost of these experiments has now moved from generating the data to storing and analysing it.

There are many steps involved in analysing an RNA-Seq experiment. Analysing an RNAseq experiment begins with sequencing reads. These are aligned to a reference genome, then the number of reads mapped to each gene can be counted. This results in a table of counts, which is what we perform statistical analyses on to determine differentially expressed genes and pathways. The purpose of this tutorial is to demonstrate how to do read alignment and counting, prior to performing differential expression with limma-voom. Differential expression analysis with limma-voom is covered in an accompanying tutorial. This tutorial shows how to start from fastq data and perform the mapping and counting steps, along with associated Quality Control. The associated tutorial, limma-voom, shows how to perform differential expression and QC on the counts.

Mouse mammary gland dataset
The data for this tutorial comes from a Nature Cell Biology paper, EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival ([Fu et al. 2015](https://www.ncbi.nlm.nih.gov/pubmed/25730472)). Both the raw data (sequence reads) and processed data (counts) can be downloaded from Gene Expression Omnibus database (GEO) under accession number [GSE60450](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450).

This study examines the expression profiles of basal stem-cell enriched cells (B) and committed luminal cells (L) in the mammary gland of virgin, pregnant and lactating mice. Six groups are present, with one for each combination of cell type and mouse status. Three replicates is usually recommended as a *minimum* requirement for RNA-seq, however, here each group contains two biological replicates, two independent sorts of cells from the mammary glands of virgin, pregnant or lactating mice.

This is a Galaxy tutorial based on material from the [COMBINE R RNAseq workshop](http://combine-australia.github.io/RNAseq-R/07-rnaseq-day2.html), first taught [here](http://combine-australia.github.io/2016-05-11-RNAseq/).


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Preparing the reads

## Import data from URLs

Read sequences are usually stored in compressed (gzipped) FASTQ files. Before the differential expression analysis can proceed, these reads must be aligned to the reference genome and counted into annotated genes. Mapping reads to the genome is a very important task, and many different aligners are available, such as HISAT2 ([Kim, Langmead, and Salzberg, 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4655817/)), STAR ([Dobin et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/)) and Subread ([Liao, Smyth, and Shi 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3664803/)]). Most mapping tasks require larger computers than an average laptop, so usually read mapping is done on a server in a linux-like environment, requiring some programming knowledge. However, Galaxy enables you to do this mapping without needing to know programming and if you don't have access to a server you can try to use one of the publically available Galaxies e.g. [usegalaxy.org](https://usegalaxy.org), [usegalaxy.eu](https://usegalaxy.eu), [usegalaxy.org.au](https://usegalaxy.org.au/). 

If you are sequencing your own data, the sequencing facility will almost always provide FASTQ.gz files which you can upload into Galaxy. If your FASTQs are provided through Galaxy's Shared Data, you can easily import them into a history. For publicly available sequence data, such as from GEO/SRA/ENA, Galaxy's Rule-Based Uploader can be used to import the files from URLs, saving on the need to download to your computer and upload into Galaxy. For more information on the Rule-Based Uploader see the tutorial [here](http://galaxyproject.github.io/training-material/topics/galaxy-data-manipulation/tutorials/upload-rules/tutorial.html).

The raw reads used in this tutorial were obtained from SRA from the link given in GEO for the the mouse mammary gland dataset (Fu et al. 2015) (e.g `ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP045%2FSRP045534`). For the purpose of this tutorial we are going to be working with a small part of the fastq files. We are only going to be mapping 1000 reads from each sample to enable running through all the steps quickly. The small fastqs files are available in [Figshare](https://figshare.com/s/f5d63d8c265a05618137) and the links to the fastqs are provided below. We are going to import the files into a Collection. Using Galaxy Collections helps keep the datasets organised and tidy in the history. Collections also make it easier to maintain the sample names through tools and workflows. If you are not familiar with collections, see the [Galaxy Collections tutorial](http://galaxyproject.github.io/training-material/topics/galaxy-data-manipulation/tutorials/collections/tutorial.html).

The sample information (sample ID, Group) and link to the FASTQ file (URL) are in the grey box below.

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

* Strip the header out of the sample information (it doesn’t contain a URL Galaxy can download).
* Define the SampleID column as the dataset identifier.
* Define the URL column as the dataset URL (this is the location Galaxy can download the data from).
* Tell Galaxy to treat these files as `fastqsanger.gz` files.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from Figshare using Galaxy's Rule-Based Uploader.
>    - Open the Galaxy Upload Manager
>    - Click the tab **Rule-based**
>        - *"Upload data as"*: `Collection(s)`
>        - *"Load tabular data from"*: `Pasted Table`
>    - Paste the table below
>    - Click **Build**
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
>            - Repeat this again and select URL instead
>                - *"URL"*: `C`
>            - Click `Apply`, and you should see your new column definitions listed
>
>        - **Specify the file type**. Select *"Type"*: `fastqsanger.gz`
>        - **Specify the genome**. Select *"Genome"*: `mm10`
>        - **Name the collection**. For *"Name"* enter: `fastqs`
>        - Click `Upload`
>        - You should see a collection (list) called `fastqs` in your history containing all 12 FASTQ files, like below.
> ![Collection](../../images/limma-voom_f2c/collection.png "Collection of FASTQs")
>
{: .hands_on}

If your data is not accessible by URL, for example, if your fastq files are located on your laptop, you can upload into a collection as below.

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

If your fastqs are located in Shared Data, you can import them into your history as a collection as below.

> ### {% icon tip %} Tip: Import files from Shared Data into a collection
>
> - In the Menu at the top go to Shared Data > Data Libraries
> - Locate your fastqs
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
> If you are not familiar with FASTQ format, see the tutorial [here](https://galaxyproject.github.io/training-material/topics/introduction/tutorials/galaxy-intro-ngs-data-managment/tutorial.html)
{: .tip}

> ### {% icon hands_on %} Hands-on: Take a look at FASTQ format
>
> 1. Click on the collection name (`fastqs`)
> 2. Click on the eye icon of one of the FASTQ files to have a look at what it contains
>
>
{: .hands_on}

## Check raw reads

During sequencing, errors are introduced, such as incorrect nucleotides being called. These are due to the technical limitations of each sequencing platform. Sequencing errors might bias the analysis and can lead to a misinterpretation of the data. Every base sequence gets a quality score from the sequencer and this information is present in the FASTQ file. A quality score of 30 corresponds to a 1 in 1000 chance of an incorrect base call (a quality score of 10 is a 1 in 10 chance of an incorrect base call). To look at the overall distribution of quality scores across the reads, we can use FastQC.

Sequence quality control is therefore an essential first step in your analysis. We will use similar tools as described in the ["Quality control" tutorial]({{site.baseurl}}/topics/sequence-analysis): [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html).

> ### {% icon hands_on %} Hands-on: Check raw reads with **FastQC**
>
> 1. **FastQC** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"Short read data from your current history"*: `fastqs` (Input dataset collection)
> 2. Inspect the WebPage output of **FastQC** {% icon tool %} for the `MCL1-DL` sample by clicking on the eye icon
>
>    > ### {% icon question %} Questions
>    >
>    > What is the read length?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > The read length is 100 bp
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
{: .hands_on}

The FastQC report contains a lot of information and we can look at the report for each sample. However, that is quite a few reports, 12 for this dataset. If you had more samples it could be a lot more. Luckily, there is a very useful tool called MultiQC that can summarise QC information for multiple samples into a single report.

> ### {% icon hands_on %} Hands-on: Aggregate FastQC data with **MultiQC**
>
> 1. **MultiQC** {% icon tool %} with the following parameters to aggregate the FastQC reports
>      - In *"Results"*
>        - *"Which tool was used generate logs?"*: `FastQC`
>        - In *"FastQC output"*
>           - *"Type of FastQC output?"*: `Raw data`
>           - {% icon param-collection %} *"FastQC output"*: `Raw data` files (output of **FastQC** {% icon tool %})
>
> 2. Inspect the webpage output from MultiQC
{: .hands_on}

Note that this is the results for just 1000 reads. The results for the full dataset are shown below (most of the plots look similar to the small fastqs).

![General Statistics](../../images/limma-voom_f2c/fastqc_table.png "General Statistics")
![Sequence Counts](../../images/limma-voom_f2c/fastqc_sequence_counts_plot.png "Sequence Counts")
![Sequence Quality](../../images/limma-voom_f2c/fastqc_per_base_sequence_quality_plot.png "Sequence Quality")
![Per Sequence Quality Scores](../../images/limma-voom_f2c/fastqc_per_sequence_quality_scores_plot.png "Per Sequence Quality Scores")
![Per Sequence GC Content](../../images/limma-voom_f2c/fastqc_per_sequence_gc_content_plot.png "Per Sequence GC Content")
![Per base N content](../../images/limma-voom_f2c/fastqc_per_base_n_content_plot.png "Per base N content")
![Sequence Duplication Levels](../../images/limma-voom_f2c/fastqc_sequence_duplication_levels_plot.png "Sequence Duplication Levels")
![Adapter Content](../../images/limma-voom_f2c/fastqc_adapter_content_plot.png "Adapter Content")

See the [ref-based tutorial](http://galaxyproject.github.io/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html) for more information on FastQC plots.

> ### {% icon question %} Questions
>
> What do you think of the overall quality of the sequences?
>
> > ### {% icon solution %} Solution
> > Overall, the samples look pretty good. The main things to note are 
> > * The base quality is high in all samples. 
> > * Some Illumina adapter has been detected. 
> > * Some duplication in RNA-seq can be normal due to the presence of highly expressed genes. However, for some reason `MCL1-LE` and `MCL1-LF` have higher numbers of duplicates detected than the other samples.
> {: .solution}
{: .question}

We will use Cutadapt to trim the reads to remove the Illumina adapter and any low quality bases at the ends. We will discard any sequences that are too short (< 20bp) after trimming. The [Cutadapt website](https://cutadapt.readthedocs.io/en/stable/guide.html#illumina-truseq) provides the sequence we can use `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` to trim the Illumina adapter. We will also output the Cutadapt report for summarising with MultiQC.

## Trim reads

> ### {% icon hands_on %} Hands-on: Trim reads with **Cutadapt**
>
> 1. **Cutadapt** {% icon tool %} with the following parameters:
>    - *"Single-end or Paired-end reads?"*: `Single-end`
>        - {% icon param-collection %} *"FASTQ/A file"*: `output` (Input dataset collection)
>        - In *"Read 1 Options"*:
>            - In *"3' (End) Adapters"*:
>                - Click on *"Insert 3' (End) Adapters"*:
>                - In *"1: 3' (End) Adapters"*:
>                    - *"Source"*: `Enter custom sequence`
>                        - *"Enter custom 3' adapter name (Optional)"*: `Illumina`
>                        - *"Enter custom 3' adapter sequence"*: `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC`
>    - In *"Filter Options"*:
>        - *"Minimum length"*: `20`
>    - In *"Read Modification Options"*:
>        - *"Quality cutoff"*: `20`
>    - In *"Output Options"*:
>        - *"Report"*: `Yes`
>
{: .hands_on}

We can take a look at the FASTQs again now that they've been trimmed.

## Check trimmed reads

> ### {% icon hands_on %} Hands-on: Check trimmed reads with **FastQC**
>
> 1. **FastQC** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"Short read data from your current history"*: `out1` (output of **Cutadapt** {% icon tool %})
> 2. **MultiQC** {% icon tool %} with the following parameters to aggregate the FastQC reports
>    - In *"Results"*
>        - *"Which tool was used generate logs?"*: `FastQC`
>        - In *"FastQC output"*
>           - *"Type of FastQC output?"*: `Raw data`
>           - {% icon param-collection %} *"FastQC output"*: `Raw data` files (output of **FastQC** {% icon tool %})
{: .hands_on}

> ### {% icon question %} Questions
>
> What differences can you see after trimming?
>
> > ### {% icon solution %} Solution
> >
> > No adapter is detected now. 
> > The sequences are no longer all 100bp, we now have sequences of different lengths.
> >
> {: .solution}
>
{: .question} 

# Mapping

Now that we have prepared our reads, we can align the reads for our 12 samples. There is an existing reference genome for mouse and we will map the reads to that. The current version of the mouse reference genome is `mm10/GRCm38`. Here we will use [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml), the descendent of widely-used TopHat, to align the reads, but alternative mappers could be used, such as STAR. See the [RNA-seq ref-based tutorial](http://galaxyproject.github.io/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html#Mapping) for more information on RNA-seq mappers. There are often numerous mapping parameters that we can specify, but usually the default mapping parameters are fine. However, library type (paired-end vs single-end) and library strandedness (stranded vs unstranded) require some different settings when mapping and counting, so they are two important pieces of information to know about samples. The mouse data comprises unstranded, single-end reads so we will specify that where necessary. If we had paired-end data or stranded data, we would specify that. HISAT2 can output a mapping summary file that tells what proportion of reads mapped to the reference genome. Summary files for multiple samples can be summarised with MultiQC. As we’re only using a subset of 1000 reads per sample, aligning should just take a minute or so. To run the full samples from this dataset would take longer.

## Map reads to reference genome

> ### {% icon hands_on %} Hands-on: Map reads to reference with **HISAT2**
>
> **HISAT2** {% icon tool %} with the following parameters:
>    - *"Source for the reference genome"*: `Use a built-in genome`
>        - *"Select a reference genome"*: `mm10`
>    - *"Single-end or paired-end reads?"*: `Single-end`
>        - {% icon param-collection %} *"FASTA/Q file"*: `out1` (output of **Cutadapt** {% icon tool %})
>    - In *"Summary Options"*:
>        - *"Output alignment summary in a more machine-friendly style."*: `Yes`
>        - *"Print alignment summary to a file."*: `Yes`
>    - In *"Advanced Options"*:
>        - *"Input options"*: `Use default values`
>        - *"Alignment options"*: `Use default values`
>        - *"Scoring options"*: `Use default values`
>        - *"Spliced alignment options"*: `Use default values`
>        - *"Reporting options"*: `Use default values`
>        - *"Output options"*: `Use default values`
>        - *"Other options"*: `Use default values`
{: .hands_on}

**HISAT2** generates a BAM file with mapped reads.

> ### {% icon hands_on %} Hands-on: Take a look at BAM format
>
> 1. Click on the collection name (`HISAT2 on collection N: aligned reads (BAM)`)
> 2. Click on the eye icon of one of the BAM files to have a look at what it contains
>
>    > ### {% icon question %} Questions
>    >
>    > What is BAM format?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > See the answer in the [ref-based tutorial](http://galaxyproject.github.io/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html#%20Pretreatments)
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
> **MultiQC** {% icon tool %} with the following parameters to aggregate the HISAT2 summary files
>    - In *"Results"*
>        - *"Which tool was used generate logs?"*: `HISAT2`
>        - {% icon param-collection %} *"Output of HISAT2"*: `summary_file` (output of **HISAT2** {% icon tool %})
{: .hands_on}

> ### {% icon question %} Questions
>
> What do you think of the results? What % reads mapped? 
>
> > ### {% icon solution %} Solution
> >
> > The mapping rate looks good as over 90% of reads mapped in all samples.
> >
> {: .solution}
>
{: .question}

It is also good to visualise the read alignments in the BAM file, for example using IGV, see the [ref-based tutorial](http://galaxyproject.github.io/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html#%20Pretreatments).

## Check strandedness

As far as is known, this data is unstranded, but as a sanity check you can check strandness. You can use RSeQC Infer Experiment tool to do this, as explained in the [RNA-seq ref-based tutorial](http://galaxyproject.github.io/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html). For this tool, and many of the other RSeQC tools, we need a reference bed file of genes (`reference genes`). RSeQC provides some files for model organisms. You can import the RSeQC mm10 RefSeq `BED` file from the link `https://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/mm10_RefSeq.bed.gz/download` or a file from Shared data if provided.

{% include snippets/import_via_link.md %}

> ### {% icon hands_on %} Hands-on: Check strandedness with **Infer Experiment**
>
> 1. **Infer Experiment** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"Input .bam file"*: `output_alignments` (output of **HISAT2** {% icon tool %})
>    - {% icon param-file %} *"Reference gene model"*: `reference genes` (Reference BED file)
> 2. **MultiQC* {% icon tool %} with the following parameters:
>       - In *"1: Results"*:
>           - *"Which tool was used generate logs?"*: `RSeQC`
>               - In *"RSeQC output"*:
>                   - Click on *"Insert RSeQC output"*:
>                   - In *"1: RSeQC output"*:
>                       - *"Type of RSeQC output?"*: `infer_experiment`
>                           - {% icon param-collection %} *"RSeQC infer_experiment output"*: `output` (output of **Infer Experiment** {% icon 
tool %})
{: .hands_on}

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

Duplicate reads are usually kept in RNA-seq differential expression analysis but it is still a good metric to check. A high percentage of duplicates can indicate a problem with the sample, for example, a low complexity library with not many transcripts, due to not enough RNA used as input. FastQC gives us an idea of duplicate reads but it just takes a sample of the data. We can get a more accurate asessment of duplciate reads using the **Picard MarkDuplicates** tool. 

> ### {% icon hands_on %} Hands-on: Check duplicate reads with **MarkDuplicates**
>
> 1. **MarkDuplicates** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"Select SAM/BAM dataset or dataset collection"*: `output_alignments` (output of **HISAT2** {% icon tool %})
> 2. **MultiQC* {% icon tool %} with the following parameters:
>       - In *"1: Results"*:
>           - *"Which tool was used generate logs?"*: `Picard`
>               - In *"Picard output"*:
>                   - Click on *"Insert Picard output"*:
>                   - In *"1: Picard output"*:
>                       - *"Type of Picard output?"*: `Markdups`
>                       - {% icon param-collection %} *"Picard output"*: `metrics_file` (output of **MarkDuplicates** {% icon tool %})
{: .hands_on}

> ### {% icon question %} Questions
>
> What are the % of duplicate reads in the samples?
>
> > ### {% icon solution %} Solution
> >
> > The results are similar to what we saw with FastQC. Most of the samples have ~50-60% duplicates but 2 samples, `MCL1-LE` and `MCL1-LF` have a higher number of duplicates (~70%).
> >
> {: .solution}
>
{: .question}

## Check number of reads for each chromosome

You can check the numbers of reads mapped to each chromosome with the **Samtools IdxStats** tool. This can help assess the sample quality, for example, if there is high mitochondrial contamination. It could also help to check the gender of sample through the numbers of reads mapping to X/Y.

> ### {% icon hands_on %} Hands-on: Count reads mapping to each chromosome with **IdxStats**
>
> 1. **IdxStats** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"BAM file"*: `output_alignments` (output of **HISAT2** {% icon tool %})
> 2. **MultiQC* {% icon tool %} with the following parameters:
>       - In *"1: Results"*:
>           - *"Which tool was used generate logs?"*: `Samtools`
>               - In *"Samtools output"*:
>                   - Click on *"Insert Samtools output"*:
>                   - In *"1: Samtools output"*:
>                      - *"Type of Samtools output?"*: `idxstats`
>                           - {% icon param-collection %} *"Samtools idxstats output"*: `output` (output of **IdxStats** {% icon tool %})
>
{: .hands_on}

> ### {% icon question %} Questions
>
> What do you think of the chromosome mappings?
>
> > ### {% icon solution %} Solution
> >
> > Some of the samples have very high mapping on chromosome 5. What is going on there?
> >
> {: .solution}
>
{: .question}

## Check coverage across gene bodies (5'-3')

The coverage of reads along gene bodies can be assessed to check if there is any bias in coverage. For example, a bias towards the 3' end of genes could indicate degradation of the RNA. Alternatively, a 3' bias could indicate that the data is from a 3' assay (e.g. oligodT-primed, 3'RNA-seq). You can use the RSeQC Gene body coverage tool to assess gene body coverage.

> ### {% icon hands_on %} Hands-on: Check coverage of genes with **Gene Body Coverage (BAM)**
>
> 1. **Gene Body Coverage (BAM)** {% icon tool %} with the following parameters:
>    - *"Run each sample separately, or combine mutiple samples into one plot"*: `Run each sample separately`
>        - {% icon param-collection %} *"Input .bam file"*: `output_alignments` (output of **HISAT2** {% icon tool %})
>    - {% icon param-file %} *"Reference gene model"*: `output` (Input dataset)
> 2. **MultiQC* {% icon tool %} with the following parameters:
>       - In *"1: Results"*:
>           - *"Which tool was used generate logs?"*: `RSeQC`
>               - In *"RSeQC output"*:
>                   - Click on *"Insert RSeQC output"*:
>                   - In *"1: RSeQC output"*:
>                       - *"Type of RSeQC output?"*: `gene_body_coverage`
>                           - {% icon param-collection %} *"RSeQC gene_body_coverage output"*: `outputtxt` (output of **Gene Body Coverage (
BAM)** {% icon tool %})
>
{: .hands_on}

> ### {% icon question %} Questions
>
> What do you think of the gene body coverage in these samples?
>
> > ### {% icon solution %} Solution
> >
> > It looks good, no bias in coverage towards 5' or 3' ends.
> >
> {: .solution}
>
{: .question}

## Check distribution of reads across features (exons, introns, UTRs..)

> ### {% icon hands_on %} Hands-on: Check distribution of reads with **Read Distribution**
>
> 1. **Read Distribution** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"Input .bam/.sam file"*: `output_alignments` (output of **MultiQC2** {% icon tool %})
>    - {% icon param-file %} *"Reference gene model"*: `output` (Input dataset)
> 2. **MultiQC* {% icon tool %} with the following parameters:
>       - In *"1: Results"*:
>           - *"Which tool was used generate logs?"*: `RSeQC`
>               - In *"RSeQC output"*:
>                   - Click on *"Insert RSeQC output"*:
>                   - In *"1: RSeQC output"*:
>                       - *"Type of RSeQC output?"*: `read_distribution`
>                           - {% icon param-collection %} *"RSeQC read_distribution output"*: `output` (output of **Read Distribution** {% 
icon tool %})
>
{: .hands_on}

> ### {% icon question %} Questions
>
> What do you think of the read distribution?
>
> > ### {% icon solution %} Solution
> >
> > It looks good, most of the reads have mapped to exons (which is what we'll be counting) and not many to introns or intergenic regions. The samples have pretty consistent read distribution, albeit with slightly higher numbers of reads mapping to CDS exons for `MCL1-LC` and `MCL1-LD`, and `MCL1-LE` and `MCL1-LF` have more reads mapping to CDS exons than the other samples.
> >
> {: .solution}
>
{: .question}

Now that we have checked the data and are happy with how it looks, we can proceed to counting.

# Counting

The alignment produces a set of BAM files, where each file contains the read alignments for each sample. In the BAM file, there is a chromosomal location for every read that mapped. Now that we have figured out where each read comes from in the genome, we need to summarise the information across genes or exons. The mapped reads can be counted across mouse genes by using a tool called featureCounts. featureCounts requires gene annotation specifying the genomic start and end position of each exon of each gene. For convenience, featureCounts contains built-in annotation for mouse (mm9, mm10) and human (hg19) genome assemblies, where exon intervals are defined from the NCBI RefSeq annotation of the reference genome. Reads that map to exons of genes are added together to obtain the count for each gene, with some care taken with reads that span exon-exon boundaries. For other species, users will need to read in a data frame in GTF format to define the genes and exons. Users can also specify a custom annotation file in SAF format. See the tool help in Galaxy, which has an example of what an SAF file should like like, or the Rsubread users guide for more information.


## Count reads mapped to genes

> ### {% icon hands_on %} Hands-on: Count reads mapped to genes with **featureCounts**
>
> **featureCounts** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"Alignment file"*: `output_alignments` (output of **MultiQC2** {% icon tool %})
>    - *"Gene annotation file"*: `featureCounts built-in`
>        - *"Select built-in genome"*: `mm10`
{: .hands_on}

## Check assignments of reads to genes

The statistics of the read mapping are output as a file from featureCounts. MultiQC can be used to assess the numbers of reads assigned to features, genes in this case.

> ### {% icon hands_on %} Hands-on: Assess reads mapped to genes
>
> **MultiQC** {% icon tool %} with the following parameters:
>    - *"Which tool was used generate logs?"*: `featureCounts`
>        - {% icon param-collection %} *"Output of FeatureCounts"*: `output_summary` (output of **featureCounts** {% icon tool %})
{: .hands_on}

featureCounts reports the numbers of unassigned reads and the reasons why they are not assigned (eg. ambiguity, multi-mapping, secondary alignment, mapping quality, fragment length, chimera, read duplicate, non-junction and so on), in addition to the number of successfully assigned reads for each library.

> ### {% icon question %} Questions
>
> How many reads are assigned to exons?
>
> > ### {% icon solution %} Solution
> >
> > ~60-70% of reads are assigned to exons. This is a fairly typical number for RNA-seq.
> >
> {: .solution}
> 
{: .question}

In this example we have kept many of the default settings, which have been optimised to work well under a variety of situations. The default setting for featureCounts is that it only keeps reads that uniquely map to the reference genome. For testing differential expression of genes, this is what we want, as the reads are unambigously assigned to one place in the genome, allowing for easier interpretation of the results. Understanding all the different parameters you can change involves doing a lot of reading about the aligner that you are using, and can take a lot of time to understand! Today we won’t be going into the details of the parameters you can change, but you can get more information from looking at the tool help.

The counts for the samples are output as tabular files. Take a look at one. The numbers in the first column of the counts file represent the Entrez gene identifiers for each gene, while the second column contains the counts for each gene for the sample.

## Create count matrix

The counts files are currently in the format of one file per sample. However, it is often convenient to have a count matrix. A count matrix is a single table containing the counts for all samples, with the genes in rows and the samples in columns. The counts files are all within a collection so we can use the Galaxy **Column Join on Collection** tool to easily create a count matrix from the single counts files.

> ### {% icon hands_on %} Hands-on: Create count matrix with **Column Join on Collection**
>
> **Column Join on Collection** {% icon tool %} with the following parameters:
>    - {% icon param-collection %} *"Tabular files"*: `output_counts` (output of **featureCounts** {% icon tool %})
>    - *"Identifier column"*: `1`
>    - *"Number of Header lines in each item"*: `1`
>    - *"Keep original column header"*: `No`
>
{: .hands_on}

Take a look at the output. Now it is easier to see the counts for a gene across all samples.

# Generating an interactive QC summary report

So far we have run MultiQC on one step at a time, which generates multiple reports for the analysis. However, it would be more convenient to have a single report containing all the QC information for the analysis.

> ### {% icon question %} Questions
>
> 1. Can you create a MultiQC report that summarises all the QC information into a single report?
>    Choose the untrimmed FastQC data (As don't think you can view both - the pre-trimmed and trimmed)
>
> > ### {% icon solution %} Solution
> >
> > 1. **MultiQC** {% icon tool %} with the following parameters:
> >    - In *"Results"*:
> >        - In *"1: Results"*:
> >            - *"Which tool was used generate logs?"*: `FastQC`*:
> >                - {% icon param-collection %} *"FastQC output"*: `FastQC RawData` (output of **FastQC** {% icon tool %})
> >        - Click on *"Insert Results"*:
> >        - In *"2: Results"*:
> >            - *"Which tool was used generate logs?"*: `Cutadapt/Trim Galore!`
> >                - {% icon param-collection %} *"Output of Cutadapt"*: `summary_file` (output of **Cutadapt** {% icon tool %})
> >        - Click on *"Insert Results"*:
> >        - In *"3: Results"*:
> >            - *"Which tool was used generate logs?"*: `HISAT2`
> >                - {% icon param-collection %} *"Output of HISAT2"*: `summary_file` (output of **HISAT2** {% icon tool %})
> >        - Click on *"Insert Results"*:
> >        - In *"4: Results"*:
> >            - *"Which tool was used generate logs?"*: `RSeQC`
> >                - In *"RSeQC output"*:
> >                    - Click on *"Insert RSeQC output"*:
> >                    - In *"1: RSeQC output"*:
> >                        - *"Type of RSeQC output?"*: `infer_experiment`
> >                            - {% icon param-collection %} *"RSeQC infer_experiment output"*: `output` (output of **Infer Experiment** {% icon 
tool %})
> >        - In *"5: Results"*:
> >            - *"Which tool was used generate logs?"*: `Picard`
> >                - In *"Picard output"*:
> >                    - Click on *"Insert Picard output"*:
> >                    - In *"1: Picard output"*:
> >                        - *"Type of Picard output?"*: `Markdups`
> >                        - {% icon param-collection %} *"Picard output"*: `metrics_file` (output of **MarkDuplicates** {% icon tool %})
> >        - Click on *"Insert Results"*:
> >        - In *"6: Results"*:
> >            - *"Which tool was used generate logs?"*: `Samtools`
> >                - In *"Samtools output"*:
> >                    - Click on *"Insert Samtools output"*:
> >                    - In *"1: Samtools output"*:
> >                       - *"Type of Samtools output?"*: `idxstats`
> >                            - {% icon param-collection %} *"Samtools idxstats output"*: `output` (output of **IdxStats** {% icon tool %})
> >        - Click on *"Insert Results"*:
> >        - In *"7: Results"*:
> >            - *"Which tool was used generate logs?"*: `RSeQC`
> >                - In *"RSeQC output"*:
> >                    - Click on *"Insert RSeQC output"*:
> >                    - In *"1: RSeQC output"*:
> >                        - *"Type of RSeQC output?"*: `gene_body_coverage`
> >                            - {% icon param-collection %} *"RSeQC gene_body_coverage output"*: `outputtxt` (output of **Gene Body Coverage (
BAM)** {% icon tool %})
> >        - Click on *"Insert Results"*:
> >        - In *"8: Results"*:
> >            - *"Which tool was used generate logs?"*: `RSeQC`
> >                - In *"RSeQC output"*:
> >                    - Click on *"Insert RSeQC output"*:
> >                    - In *"1: RSeQC output"*:
> >                        - *"Type of RSeQC output?"*: `read_distribution`
> >                            - {% icon param-collection %} *"RSeQC read_distribution output"*: `output` (output of **Read Distribution** {% 
icon tool %})
> >        - Click on *"Insert Results"*:
> >        - In *"9: Results"*:
> >            - *"Which tool was used generate logs?"*: `featureCounts`
> >                - {% icon param-collection %} *"Output of FeatureCounts"*: `output_summary` (output of **featureCounts** {% icon tool %})
> >
> {: .solution}
> 
{: .question}

> ### {% icon question %} Questions
>
> Are there any other QCs that could be performed on RNA-seq data?
>
> > ### {% icon solution %} Solution
> >
> > The data could be checked for 
> > * Ribosomal contamination
> > * Contamination with other species e.g. bacteria
> > * This is single-end data but paired-end data could be checked for fragment size (distance between the read pairs).
> >
> {: .solution}
> 
{: .question}

# Creating a workflow of the analysis

You can now extract a workflow from the steps that have been performed, as shown in the [Galaxy Peaks to Genes tutorial](https://galaxyproject.github.io/training-material/topics/introduction/tutorials/galaxy-intro-peaks2genes/tutorial.html#extracting-workflow). This is a way to help keep a record of the steps performed. It also enables efficient rerunning of a multi-step analysis, such as RNA-seq.

![fastq to counts Workflow](../../images/limma-voom_f2c/workflow.png "FASTQ to counts Workflow")

# Conclusion
{:.no_toc}

In this tutorial we have seen how fastq files can be converted into counts. We have also seen QC steps that can be performed to help assess the quality of the data. The accompanying tutorial, limma-voom, shows how to perform differential expression and QC on the counts.