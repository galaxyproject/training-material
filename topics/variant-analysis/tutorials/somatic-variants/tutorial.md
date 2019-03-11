---
layout: tutorial_hands_on

title: Identification of somatic and germline variants from tumor and normal sample
  pairs
zenodo_link: ''
questions:
- What are the specific challenges in somatic variant calling that set it apart
  from regular diploid variant calling?
- How can you call variants and classify them according to their
  presence/absence in/from tumor and normal tissue of the same individual?
- How can you annotate variants and affected genes with prior knowledge from
  human genetic and cancer-specific databases to generate clinically relevant
  reports?
objectives:
- Call variants and their somatic status from whole-exome sequencing data
- Annotate variants with a wealth of human genetic and cancer-specific
  information extracted from public databases
- Add gene-level annotations and generate reports of annotated somatic and
  germline variants, loss-of-heterozygosity (LOH) events, and affected genes,
  ready for interpretation by clinicians
time_estimation: 7h
key_points:
- "Follow best practices for read mapping, quality control and mapped reads
  postprocessing to minimze false-positive variant calls."
- "Use a dedicated somatic variant caller to call variants and to classify them
  into somatic, germline and LOH event variants on statistical grounds."
- "Annotations and queries based on variant properties add relevance to variant
  and gene reports."
- "A framework like GEMINI is very helpful for managing, annotating and querying
  lists of variants in a flexible way."
- "Prefer public, free annotation sources to foster reproducibility and
  information sharing."
contributors:
- wm75

---


# Introduction
{:.no_toc}

When sequencing genomic material from a human tumor, the underlying clinical or
research question typically is **what spectrum of mutations distinguishes this
tumor from healthy tissue**.

This question cannot be answered adequately just by comparing the tumor tissue
to the human reference genome because even in healthy tissue there will be many
thousands of variants compared to the reference genome. This is because every
individual inherits a unique pattern of that many variants from her parents. A
fundamental difference between these variants and the tumor-specific mutations
is that the former are present in the carrier's germline, while the latter have
been acquired somatically and will, thus, not be transmitted to offspring.
Therefor, we talk of germline variants to refer to variants present in healthy
and tumor tissue alike, and of somatic variants to refer to tumor-specific
variants. To be able to distinguish between these two types of variants always
requires a direct comparison of data from tumor and normal tissue samples. 

In addition to acquiring new variants, tumors can also lose or gain chromosomal
copies of variants found heterozygously in an individual's germline. This
phenomenon is termed *loss of heterozygosity* (LOH) because only one of the two
original alleles persists in the tumor (either in a hemizygous state if the
other allele is simply dropped, or in a homozygous state in the case of a
duplication of one allelic copy accompanied by loss of the other).
The detection of LOH events, again, is dependent on a comparison of tumor and
normal tissue data.

In this tutorial we are going to identify somatic and germline variants, as
well as variants affected by LOH, from a tumor and a normal sample of the same
patient. Our goal is to report the variant sites, and the genes affected by
them, annotated with the content of general human genetic and cancer-specific
databases. Ideally, this may provide insight into the genetic events driving
tumor formation and growth in the patient, and might be of prognostic and even
therapeutic value by revealing variants known to affect drug
resistance/sensitivity, tumor aggressiveness, *etc*.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data Preparation

First we need to upload and prepare some input data to analyze. The sequencing
reads we are going to analyze are from real-world data from a cancer patient's
tumor and normal tissue samples.
For the sake of an acceptable speed of the analyis, the original data has been
downsampled though to include only the reads from human chromosomes 5, 7 and
12.

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial and give it a meaningful name
>
>    {% include snippets/create_new_history.md %}
>    {% include snippets/rename_history.md %}
>
> 2. Import the following four files from
>    [Zenodo](https://zenodo.org/record/2582555):
>
>    ```
>    https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
>    https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
>    https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
>    https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz
>    ```
>    where the first two files represent the forward and reverse reads sequence
>    data from a patient's **normal** tissue, and the last two represent the
>    data of the same patient's **tumor** tissue.
> 
>    Alternatively, the same files may be available on your Galaxy server
>    through a shared data library (your instructor may tell you so), in which
>    case you may prefer to import the data directly from there.
>
>    {% include snippets/import_via_link.md format="fastqsanger.gz" %}
>
>    {% include snippets/import_from_data_library.md %}
>
> 3. Import the reference genome
>
>    > ### {% icon comment %} Shortcut
>    > You can skip this step if the Galaxy server you are working on offers
>    > a `hg19` version of the human reference genome with prebuilt indexes for
>    > *bwa-mem* **and** *samtools* (ask your instructor, or check the tools
>    > **Map with BWA-MEM** {% icon tool %} and **VarScan Somatic**
>    > {% icon tool %} if they list a `hg19` version as an option under
>    > **(Using) reference genome**).
>    {: .comment}
>
>    Import the file:
>    `https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz`
>    from [Zenodo](https://zenodo.org/record/2582555), which contains the
>    `hg19` version of the sequences of human chromosomes 5, 12 and 17,
>    and make sure you set its datatype to `fasta`.
>
>    Alternatively, load the dataset from a shared data library.
>
> 4. At this point, check that all newly created datasets in your history
>    have their datatypes assigned correctly, and fix any missing or wrong
>    datatype assignment.
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 5. Rename the datasets and add appropriate tags to them
>
>    {% include snippets/rename_dataset.md %}
>    > 1. For datasets that you upload via a link, Galaxy will pick the link
>    >    address as the dataset name, which you will likely want to shorten to
>    >    just the file names.
>    > 2. The reference genome you have imported above came as a compressed
>    >    file, but got unpacked by Galaxy to plain `fasta` format according to
>    >    your datatype selection. You may now wish to remove the `.gz` suffix
>    >    from the dataset name.
>    > 3. Large parts of the analysis in this tutorial will consist of
>    >    identical steps performed on the normal and on the tumor tissue data
>    >    in parallel.
>    >
>    >    > ### {% icon tip %} Tip: Dataset tags
>    >    > To make it easier to keep track of which dataset represents which
>    >    > branch of an analysis in a linear history, Galaxy supports dataset
>    >    > tags. In particular, if you attach a tag starting with `#` to any
>    >    > dataset, that tag will automatically propagate to any new dataset
>    >    > derived from the tagged dataset.
>    >    >
>    >    > To add a tag, simply click on the dataset you want to tag, then
>    >    > select {% icon galaxy-tags %} **Edit dataset tags** and start
>    >    > typing in the text box that appears.
>    >    {: .tip}
>    >
>    >    Before starting our analysis it is, thus, a good idea to tag the two
>    >    fastq datasets representing the normal tissue with, *e.g.*, `#normal`
>    >    and the two datasets representing the tumor tissue with, *e.g.*,
>    >    `#tumor`.
>
{: .hands_on}

# Quality control and mapping of NGS reads

Before starting our analysis, we would like to make sure that the input data
is of good quality, *i.e.*, that there haven't been any major issues during
DNA preparation, exon capture, or during actual sequencing. To avoid spurious
variant calls due to low input quality, we can ensure that all sequencing reads
used in the analysis meet some minimal quality criteria by trimming low-quality
parts off of the ends of reads and/or discarding reads of poor quality
altogether. The resulting set of polished reads then needs to be mapped to the
human reference genome because knowing the genomic positions that the bases of
a read provide evidence for is, of course, a prerequisite for variant calling.

> ### {% icon comment %} More on quality control and mapping
> If you would like to explore the topics of quality control and read mapping
> in detail, you should take a look at the separate
> [Quality Control](../../../sequence-analysis/tutorials/quality-control/tutorial.html)
> and [Mapping](../../../sequence-analysis/tutorials/mapping/tutorial.html)
> tutorials.
> Here, we will only illustrate the concrete steps necessary for quality
> control and read mapping of our particular datasets.
{: .comment}

## Quality control
> ### {% icon hands_on %} Hands-on: Quality control of the input datasets
> 1. Run **FastQC** {% icon tool %} on each of your four fastq datasets
>
>    The easiest way to do this is by specifying multiple input files under
>    *“Short read data from your current history”*.
>    {% include snippets/select_multiple_datasets.md %}
>    When you start this job, eight new datasets (one with the calculated raw
>    data, another one with an html report of the findings for each input
>    dataset) will get added to your history.
> 2. Use **MultiQC** {% icon tool %} to aggregate the raw **FastQC** data of
>    all four input datasets into one report
>
>    by populating the *"Results"* block in the tool interface like this:
>    - *"Which tool was used generate logs?"*: `FastQC`
>    - Under *"FastQC output"*
>      - *"Type of FastQC output?"*: `Raw data`
>      - As {% icon param-files %} *"FastQC output"*, select all four *RawData*
>        datasets produced by **FastQC** {% icon tool %}
> 3. Inspect the *Webpage* output produced by the tool
>    > ### {% icon question %} Questions
>    >
>    > 1. What do you think of the base qualities of the sequences?
>    > 2. Which aspect of the quality report is most puzzling to you?
>    >    (Hint: Have a look at the GC content of the reads)
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > 1. The quality of the sequences looks promising. There are no
>    > >    discernible systematic problems with it.
>    > >
>    > >    Even the reverse reads, which are typically of somewhat poorer
>    > >    quality than the corresponding forward reads, look good on average. 
>    > > 2. The GC content plots of the forward and the reverse reads from both
>    > >    samples reveal a very peculiar bimodal distribution.
>    > >
>    > >    Typically, a non-normal distribution of the GC content of the reads
>    > >    from a sample are considered to hint at possible contamination.
>    > >    Here, however, we are dealing with sequencing data from captured
>    > >    exomes, *i.e*, the reads are not representing random sequences from
>    > >    a genome, but rather an arbitrary selection.
>    > >    In fact, the samples at hand were prepared using Agilent's
>    > >    SureSelect V5 technology for exome enrichment, and bimodal GC
>    > >    content distributions can be seen as a hallmark of that capture
>    > >    method in several publications.
>    > >
>    > {: .solution}
>    {: .question}
{: .hands_on}

## Read trimming and filtering
Although the raw reads used in this tutorial are of relatively good overall
quality already, we will apply read trimming and filtering to see if we can
improve things still a bit, but also to demonstrate the general concept.

> ### {% icon hands_on %} Hands-on: Read trimming and filtering
> We will use **Trimmomatic** {% icon tool %} for this task.
>
> 1. Trim and filter the **normal tissue** reads
>    - Select *"Single-end or paired-end reads?"*: `Paired-end (two separate
>      input files)`
>      
>      This makes the tool treat the forward and reverse reads simultaneously.
>
>      - As *"Input FASTQ file (R1/first of pair)"* choose the forward reads
>        (r1) dataset of the normal tissue sample
>      - As *"Input FASTQ file (R2/second of pair)"* choose the reverse reads
>        (r2) dataset of the normal tissue sample
>
>    As further settings use:
>    - *"Perform initial ILLUMINACLIP step?"*: `Yes`
>      - *"Select standard adapter sequences or provide custom?": `Standard`
>        - *"Adapter sequences to use"*: `TruSeq3 (paired-ended, for MiSeq and
>          HiSeq)`
>      - *"Maximum mismatch count which will still allow a full match to be
>        performed"*: `2`
>      - *"How accurate the match between the two 'adapter ligated' reads must
>        be for PE palindrome read alignment"*: `30`
>      - *"How accurate the match between any adapter etc. sequence must be
>        against a read"*: `10`
>      - *"Minimum length of adapter that needs to be detected (PE specific/
>        palindrome mode)"*: `8`
>      - *"Always keep both reads (PE specific/palindrome mode)?"*: `Yes`
>
>    to cut ILLUMINA-specific adapter sequences from the reads.
>
>    Then specify the following three trimming and filtering operations:
>    1. *"Select Trimmomatic operation to perform"*: `Cut the specified number
>       of bases from the start of the read (HEADCROP)`
>       - *"Number of bases to remove from the start of the read"*: `3`
>    2. *"Select Trimmomatic operation to perform"*: `Cut bases off the end of
>       a read, if below a threshold quality (TRAILING)`
>       - *"Minimum quality required to keep a base"*: `10`
>    3. *"Select Trimmomatic operation to perform"*: `Drop reads below a
>       specified length (MINLEN)`
>       - *"Minimum length of reads to be kept"*: `25`
>
>    that are to be applied to the reads in the given order after adapter
>    trimming.
>
>    Running this job will generate four output datasets:
>    - two for the trimmed forward and reverse reads that still have a proper
>      mate in the other dataset
>    - two more datasets of orphaned forward and reverse reads, for which the
>      corresponding mate got dropped because of insufficient length after
>      trimming; when you inspect these two files, however, you should find
>      that they are empty because none of our relatively high quality reads
>      got trimmed that excessively. You can delete the two datasets to keep
>      your history more compact.
>
>    > ### {% icon details %} More on handling unpaired reads in paired-end data
>    > Splitting out potential unpaired reads into separate datasets like this
>    > is important because read mappers, typically, expect reads in forward
>    > and reverse input datasets to be arranged in proper pairs, and reads
>    > in one of the datasets without a counterpart in the other would destroy
>    > that expected structure. Therefor, when your data is paired-end data
>    > always make sure you use Trimmomatic in paired-end mode, too!
>    {: .details}
>
> 2. Trim and filter the **tumor tissue** reads
>
>    - Follow the same steps as above - just change the two input datasets - to
>      treat the tumor tissue reads with identical settings.
>
>      As before, you should observe that the two unpaired reads datasets are
>      empty, and delete them.
>    
{: .hands_on}

> ### {% icon hands_on %} Exercise: Quality control of the polished datasets
> Use **FastQC** {% icon tool %} and **MultiQC** {% icon tool %} like before,
> but using the four trimmed datasets produced by Trimmomatic as input.
>
>    > ### {% icon question %} Questions
>    >
>    > How did read trimming affect the quality reports?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > As expected, trimming the relatively high-quality raw reads did not
>    > > have any substantial impact on average. A small fraction of successful
>    > > adapter removal is visible though.
>    > >
>    > {: .solution}
>    {: .question}
{: .hands_on}

## Read Mapping
> ### {% icon hands_on %} Hands-on: Read Mapping
> Use **Map with BWA-MEM** {% icon tool %} to map the reads from each sample
> to the reference genome.
>
> 1. Map the reads from the **normal tissue** sample
>
>    If you determined at the Data Preparation step that you will be able to
>    follow the tutorial with `hg19` reference genomes and indexes preinstalled
>    on your Galaxy server instance:
>
>    - Under *"Will you select a reference genome from your history or use a
>      built-in index?"* select `Use a built-in genome index`
>      - *"Using reference genome"*: `Human: hg19` (or a similarly named
>        option)
> 
>    If you have imported the `hg19` sequence as a fasta dataset into your
>    history instead:
>
>    - Under *"Will you select a reference genome from your history or use a
>      built-in index?"* select `Use a genome from history and build index`
>      - Under *"Use the following dataset as the reference sequence"* select
>        your imported `hg19` fasta dataset.
>
>    We have paired-end (forward and reverse reads) for our samples, so select:
>    - *"Single or Paired-end reads"*: `Paired`
>
>    Since we are going to map the reads of the normal tissue sample, under
>    - *"Select first set of reads"*: choose the forward reads (r1) dataset of
>      the normal tissue sample
>    - *"Select second set of reads"*: choose the reverse reads (r2) dataset of
>      the normal tissue sample
>    - Under *"Set read groups information?"* select `Set read groups (SAM/BAM
>      specification)`
>      - *"Auto-assign"*: `No`
>        - *"Read group identifier (ID)"*: `231335` (this value being taken
>          from the original file name of the normal tissue input files)
>        - *"Auto-assign"*: `No`
>          - *"Read group sample name (SM)"*: `Normal`
>
>    > ### {% icon details %} More on read group identifiers and sample names
>    > In general, you are free to choose ID and SM values to your liking,
>    > but the ID should unambiguously identify the sequencing run that
>    > produced the reads, while the SM value should identify the
>    > biological sample.
>    {: .details}
>
> 2. Map the reads from the **tumor tissue** sample
> 
>    - Follow the same steps as above, but when you get to *"Select first set
>      of reads"* and *"Select second set of reads"*, be sure to select the two
>      datasets of the tumor tissue sample.
>    - Additionaly, make these changes:
>      - *"Read group identifier (ID)"*: `231336` (this value again being taken
>          from the original file name of the tumor tissue input files)
>      - *"Read group sample name (SM)"*: `Tumor`
>
{: .hands_on}


# Mapped reads postprocessing

To ensure that we base our variant analysis only on unambiguous, high-quality
read mappings we will do some postprocessing now.

## Filtering on mapped reads properties

> ### {% icon hands_on %} Hands-on: Filtering for mapping status and quality
>
> To produce new filtered BAM datasets with only those reads retained that
> have a minimal mapping quality of 1 and are mapped in a proper pair, use
> **Filter SAM or BAM, output SAM or BAM** {% icon tool %} with the following
> parameters (leaving non-mentioned ones at their defaults):
>   - *"SAM or BAM file to filter"*: change to {% icon param-files %}
>     **Multiple datasets** mode, then select your mapped reads datasets from
>     the normal *and* the tumor tissue data
>   - *"Minimum MAPQ quality score"*: `1`
>   - *"Filter on bitwise flag"*: `yes`
>
>     Under *"Only output alignments with all of these flag bits set"*,
>     - Check {% icon param-check %} *"Read is mapped in a proper pair"*
>
> This will produce two new datasets, one for each of the normal and tumor
> data.
{: .hands_on}

## Removing duplicate reads

> ### {% icon hands_on %} Hands-on: Remove duplicates
>
> Use **RmDup** {% icon tool %} with the following parameters:
>   - *"BAM file"*: change to {% icon param-files %}
>     **Multiple datasets** mode, then select your filtered reads datasets from
>     the normal *and* the tumor tissue data
>   - *"Is this paired-end or single end data"*: `BAM is paired-end`
>   - *"Treat as single-end"*: `No`
>
> Again, this will produce two new datasets, one for each of the normal and
> tumor data.
{: .hands_on}

## Left-align reads around indels

> ### {% icon hands_on %} Hands-on: Left-align
>
> Use **BamLeftAlign** {% icon tool %} with the following parameters:
> - *"Choose the source for the reference genome"*:
>
>   - If you are following the tutorial using `hg19` reference genomes and
>     indexes preinstalled on your Galaxy server instance, select:
>     `Locally cached`
>
>     then select *"Using reference genome"*: `Human: hg19` (or a similarly
>     named choice).
>   - If you are using an imported `hg19` sequence in your history instead,
>     select: `History`
>
>     then under *"Using reference file"*: choose your reference genome
>     dataset.
>
>   - In either case, change *"BAM dataset to re-align"* to
>     {% icon param-files %} **Multiple datasets** mode, and select your
>     filtered and deduplicated reads datasets from the normal *and* the tumor
>     tissue data.
> - *"Maximum number of iterations"*: `5`
>
> As before, this will produce two new datasets, one for each of the normal and
> tumor data.
{: .hands_on}

## Recalibrate read base and mapping qualities

> ### {% icon hands_on %} Hands-on: Recalibrate read quality scores
>
> Use **CalMD** {% icon tool %} with the following parameters:
> - *"BAM file to recalculate"*: change to {% icon param-files %}
>   **Multiple datasets**, then select your left-aligned datasets from the
>   normal *and* the tumor tissue data
> - *"Choose the source for the reference genome"*:
>
>   - If you are following the tutorial using `hg19` reference genomes and
>     indexes preinstalled on your Galaxy server instance, select:
>     `Use a built-in genome`
>
>     then select *"Using reference genome"*: `Human: hg19` (or a similarly
>     named choice).
>   - If you are using an imported `hg19` sequence in your history instead,
>     select: `Use a genome from the history`
>
>     then under *"Using reference file"*: choose your reference genome
>     dataset.
>
> - *"Do you also want BAQ (Base Alignment Quality) scores to be calculated?"*:
>   `No`
>
> > ### {% icon details %} More on base quality adjustment
> > The *VarScan somatic* tool that we are going to use for calling variants at
> > the next step is typically used in combination with unadjusted base quality
> > scores because the general opinion is that the base quality downgrades
> > performed by *CalMD* and other tools from the *samtools* suite of tools
> > are too severe for *VarScan* to retain good sensitivity. We are sticking
> > with this practice in this tutorial.
> >
> > If, for your own data, you would like to experiment with adjusted base
> > quality scores, it is important to understand that *VarScan somatic* will
> > only make use of the adjusted scores if they are incorporated directly into
> > the read base qualities of a BAM input dataset, but not if they are written
> > to the dataset separately. Hence, should you ever decide to use:
> > - *"Do you also want BAQ (Base Alignment Quality) scores to be
> >   calculated?"*: `Yes, run BAQ calculation`
> >
> >   and you want this setting to affect downstream variant calling with
> >   *VarScan somatic* make sure you also set then:
> >   - *"Use BAQ to cap read base qualities"*: `Yes`
> >
> > Please also note that BAQ scores are quite expensive to calculate so be
> > prepared to see an enormous (up to 10x) increase in job run time when
> > enabling it.
> {: .details}
> - *"Additional options"*: `Advanced options`
>   - *"Change identical bases to '='"*: `No`
>   - *"Coefficient to cap mapping quality of poorly mapped reads"*: `50`
>
> This will, once more, produce two new datasets, one for each of the normal
> and tumor data.
{: .hands_on}


# Variant calling and classification

> ### {% icon hands_on %} Hands-on: Variant calling and classification
>
> Use **VarScan somatic** {% icon tool %} with the following parameters:
> - *"Will you select a reference genome from your history or use a built-in
>   genome?"*:
>
>   - If you are following the tutorial using `hg19` reference genomes and
>     indexes preinstalled on your Galaxy server instance, select:
>     `Use a built-in genome`
>
>     then select *"reference genome"*: `Human: hg19` (or a similarly
>     named choice).
>   - If you are using an imported `hg19` sequence in your history instead,
>     select: `Use a genome from my history`
>
>     then select the *"reference genome"* from your history.
>
> - As *"aligned reads from normal sample"* you want to select the mapped and
>   fully post-processed normal tissue dataset (*i.e.*, one of the two outputs
>   produced by CalMD)
> - As *"aligned reads from tumor sample"* select the tumor tissue dataset
>   produced by CalMD
>
> - *"Estimated purity (non-tumor content) of normal sample"*: `1`
> - *"Estimated purity (tumor content) of tumor sample"*: `0.5`
>
> - *"Generate separate output datasets for SNP and indel calls?"*: `No`
>
> Of the many, many additional options offered by the tool, we only want to
> change two:
> - Under *"Settings for Variant Calling"* select `Customize settings`
>
> then change the following two parameters from their defaults:
> - *"Minimum base quality"*: `28`
>
>   We have seen, at the quality control step, that our sequencing data is of
>   really good quality, and we have chosen not to downgrade base qualities at
>   the quality scores recalibration step above, so we can increase the base
>   quality required at any given position without throwing away too much of
>   our data. 
> - *"Minimum mapping quality"*: `1`
>
>   During postprocessing, we have filtered our reads for ones with a mapping
>   quality of at least one, so requiring this quality also here does not
>   actually change the results, but it makes the requirement more explicit.
>
> Leave all other settings at their default values and also go with:
> - *"Settings for Posterior Variant Filtering"*: `Use default values`
{: .hands_on}

# Variant annotation and reporting

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Import the following **variant annotation files** from
>    [Zenodo](https://zenodo.org/record/2581873):
>
>    ```
>    https://zenodo.org/record/2581873/files/hotspots.bed
>    https://zenodo.org/record/2581873/files/cgi_variant_positions.bed
>    https://zenodo.org/record/2581873/files/01-Feb-2019-CIVic.bed
>    https://zenodo.org/record/2582555/files/dbsnp.b147.chr5_12_17.vcf.gz
>    ```
>    Make sure you select `bed` as the datatype for the first three files and
>    `vcf` for the last file.
>
>    Alternatively, add the files from a shared data library on your Galaxy
>    server instance.
>
> 2. Import some **gene-level annotation files** from
>    [Zenodo](https://zenodo.org/record/2581881):
>
>    ```
>    https://zenodo.org/record/2581881/files/Uniprot_Cancer_Genes.13Feb2019.txt
>    https://zenodo.org/record/2581881/files/cgi_genes.txt
>    https://zenodo.org/record/2581881/files/01-Feb-2019-GeneSummaries.tsv
>    ```
>    and make sure you select `tabular` as their datatype, or add them from
>    the shared data library.
>
> 3. Download **SnpEff functional genomic annotations**
>
>    > ### {% icon comment %} Shortcut
>    > You can skip this step if the Galaxy server you are working on offers
>    > `Homo sapiens: hg19` as a locally installed snpEff database. You can
>    > check the **SnpEff eff** {% icon tool%} tool under **Genome source** to
>    > see if this is the case.
>    {: .tip}
>
>    Use **SnpEff Download** {% icon tool %} to download genome annotation
>    database `hg19`.
>
> 4. Check that the datatypes of all new datasets have been set correctly, and
>    change them if necessary. You may also want to shorten some of the
>    dataset names.
>
{: .hands_on}

## Adding annotations to the called variants

### Adding functional genomic annotations

Certainly, not all variants are equal. Many may just be silent mutations with
no effect at the amino acid level, while a few others may be disrupting the
coding sequence of a protein by introducing a premature stop codon or a
frameshift. Of course, it is also important to know *which* gene is affected
by a variant.
Such functional genomic annotations can be added to a VCF dataset of variants
with *SnpEff*.

> ### {% icon hands_on %} Hands-on: Adding annotations with SnpEff
> In **SnpEff eff** {% icon tool %}
> - set *"Sequence changes (SNPs, MNPs, InDels)"* to the output from VarScan
>   somatic
> - *"Input format"*: `VCF`
> - *"Output format"*: `VCF (only if input is VCF)`
>
> If you just downloaded a `hg19` SnpEff genome database to your history,
> select now:
> - *"Genome source"*: `Reference genome from your history` and
> - choose the *"SnpEff4.3 Genome Data"* from your history.
>
> To annotate your variants with a snpEff database pre-installed on the Galaxy
> server instead, select:
> - *"Genome source"*: `Locally installed reference genome` and
> - pick *"Genome"*: `Homo sapiens: hg19` (or a similarly named option).
>
> You can leave all other options at their defaults, but to keep your history
> a bit more compact, you may scroll down to near the end of the page and set:
> - *"Produce Summary Stats"*: `No`
{: .hands_on}

### Adding genetic and clinical evidence-based annotations

Other interesting pieces of information about a variant include aspects like
whether the variant has been observed before in the human population and,
if so, at which frequency. If a variant is known to be associated with specific
diseases, we would also very much like to know that.
To proceed with this kind of genetic and clinical evidence-based annotations,
we are going to convert our list of variants into a database that can be
handled more efficiently than a VCF dataset. We will use the *GEMINI* tool
suite for this task and for all further work with the variants.

> ### {% icon hands_on %} Hands-on: Creating a GEMINI database from a variants dataset
> Run **GEMINI load** {% icon tool %} with
> - *"VCF dataset to be loaded in the GEMINI database"* set to the output of
>   SnpEff eff.
>
> Choose the following further settings:
> - *"The variants in this input are"*: `annotated with snpEff`
>
>   because, obviously, that is the case.
> - *"This input comes with genotype calls for its samples"*: `Yes`
>
>   Calling sample genotypes was part of what we used VarScan somatic for.
> - *"Choose a gemini annotation source"*: select the latest available
>   annotations snapshot (most likely, there will be only one)
> - *"Sample and family information in PED format"*: leave at
>   `Nothing selected`
> - Under *"Load the following optional content into the database"* check the
>   following options:
>   - {% icon param-check %} *"GERP scores"*
>   - {% icon param-check %} *"CADD scores"*
>   - {% icon param-check %} *"Gene tables"*
>   - {% icon param-check %} *"Sample genotypes"*
>   - {% icon param-check %} *"variant INFO field"*
>
>   but leave **unchecked** the following:
>   - *"Genotype likelihoods (sample PLs)"*
>
>     because VarScan somatic does not generate these values
>   - *"only variants that passed all filters"*
>
>     because it is simple and more flexible to filter for this property later
{: .hands_on}

During the creation of the database *GEMINI* already (silently) adds an
impressive amount of annotations it knows about to our variants (including,
*e.g.*, the frequency at which every variant has been observed in large human
genome sequencing projects). We have gotten all of these for free, just by converting the variants to a GEMINI database!

*GEMINI* also extracts a lot of the information stored in the VCF input dataset for us (such as the functional genomic annotations that we added with
*SnpEff*).

However, there are typically additional annotations from other sources (not
incorporated into *GEMINI*) that one would like to add. In addition, *GEMINI*
is not prepared to extract some non-standard information from VCF datasets,
including some important bits added by *VarScan somatic*.

**GEMINI annotate** {% icon tool %} is the tool that is designed to help you
with these rather common issues. It lets you add further annotations to the
variants in an already loaded GEMINI database.

As a first step, we are going to use the tool to add some crucial information
generated by *VarScan somatic*, but not recognized by *GEMINI load*, to the
database. Specifically, we are interested in three values calculated by
*VarScan somatic* for each variant it detected:

- *Somatic status (SS)*

  This is a simple numeric code, in which a value of `1` indicates a germline
  variant, `2` a somatic variant and `3` an LOH event.
  
- *Germline p-value (GPV)*

  For variants with a somatic status of 1, this is the error probability
  associated with that status call.
  
- *Somatic p-value (SPV)*

  This is the error probability associated with status calls of 2 and 3
  (somatic and LOH calls).
  
These values are encoded in the *INFO* column of the VarScan-generated VCF
dataset and we are going to extract them from there and add them to the
*GEMINI* database.

> ### {% icon details %} Handling of the VCF INFO field by GEMINI load
> If you paid close attention to how we generated the GEMINI database, you
> might remember that, under *"Load the following optional content into the
> database"*, we checked the option {% icon param-check %} *"variant INFO
> field"*. Did the tool not do what it was supposed to do?
>
> The answer is that it did, but not in the way you may expect.
>
> In fact, *GEMINI load* always extracts INFO field elements it knows about and
> stores them into predefined columns of the variants table in the database.
> In addition, it can store the *whole* INFO field of each variant into a
> separate *info* column so that no information gets lost on the way from the
> VCF input to the GEMINI database, and that's the optional behavior we have
> been asking for. To keep things compact, however, this column stores the INFO
> content in compressed form, which is not readily accessible.
>
> That's why it is necessary to extract non-standard INFO field elements
> explicitly if they are supposed to be used for filtering and querying.
{: .details}

> ### {% icon hands_on %} Hands-on: Making variant call statistics accessible
> Configure these parameters of **GEMINI annotate** {% icon tool %}:
> - As *"GEMINI database"* select the dataset generated by GEMINI load in the
>   previous step.
> - *"Dataset to use as the annotation source"* is the VarScan-generated vcf
>   dataset from your history.
> - *"Strict variant-identity matching of database and annotation records
>   (VCF format only)"*: `Yes`
>
>   This setting does not really matter here since you have built the GEMINI
>   database from the exact same list of variants that we are now retrieving
>   annotations from and because VarScan somatic does not call multiple alleles
>   at single sites. Matching on variant-identity is the behavior we would
>   like to see though, so we may as well be explicit about it.
> - *"Type of information to add to the database variants"*: `Specific values
>   extracted from matching records in the annotation source (extract)`
>
>    Now fill in the following three {% icon param-repeat %}
>    *"Annotation extraction recipe[s]"*:
>    1. Recipe for extracting the VarScan-generated *SS* field and adding it as
>       a new column *somatic_status* to the GEMINI database
>       - *"Elements to extract from the annotation source"*: `SS`
>       - *"Database column name to use for recording annotations"*:
>         `somatic_status`
>       - *"What type of data are you trying to extract?"*: `Integer numbers`
>       - *"If multiple annotations are found for the same variant,
>         store ..."*: `the first annotation found`
>    2. Recipe for extracting the VarScan-generated *GPV* field and adding it
>       as a new column *germline_p* to the GEMINI database
>       - *"Elements to extract from the annotation source"*: `GPV`
>       - *"Database column name to use for recording annotations"*:
>         `germline_p`
>       - *"What type of data are you trying to extract?"*: `Numbers with
>         decimal precision`
>       - *"If multiple annotations are found for the same variant,
>         store ..."*: `the first annotation found`
>    3. Recipe for extracting the VarScan-generated *SPV* field and adding it
>       as a new column *somatic_p* to the GEMINI database
>       - *"Elements to extract from the annotation source"*: `SPV`
>       - *"Database column name to use for recording annotations"*:
>         `somatic_p`
>       - *"What type of data are you trying to extract?"*: `Numbers with
>         decimal precision`
>       - *"If multiple annotations are found for the same variant,
>         store ..."*: `the first annotation found`
>
>    *"Execute"* the job.
{: .hands_on}

Next, we are going to add additional annotations beyond the ones directly
obtainable through *GEMINI* or from the input VCF dataset.

> ### {% icon hands_on %} Hands-on: Adding further annotations
> Use **GEMINI annotate** {% icon tool %} to add further annotations from four
> different sources.
> 1. Add information from **dbSNP**
>
>    As part of the database creation process, GEMINI already checks all
>    variants whether they occur in dbSNP and, if so, stores their dbSNP IDs.
>    Here, we are only storing some additional information of interest.
>
>    Configure these parameters of **GEMINI annotate** {% icon tool %}:
>    - As *"GEMINI database"* select the dataset generated by GEMINI annotate
>      in the previous step.
>    - *"Dataset to use as the annotation source"* is the
>      `dbsnp.b147.chr5_12_17.vcf` dataset imported in the Get Data step of
>      this section.
>    - *"Strict variant-identity matching of database and annotation records
>      (VCF format only)"*: `Yes`
>
>      dbSNP stores information about specific SNPs observed in human
>      populations and we would like to record if any exact same SNPs are among
>      our variants.
>    - *"Type of information to add to the database variants"*: `Specific
>      values extracted from matching records in the annotation source
>      (extract)`
>
>    Now fill in the following three {% icon param-repeat %}
>    *"Annotation extraction recipe[s]"*:
>    1. Recipe for extracting the dbSNP *SAO* field and adding it as
>       *rs_ss* to the GEMINI database
>       - *"Elements to extract from the annotation source"*: `SAO`
>       - *"Database column name to use for recording annotations"*:
>         `rs_ss`
>       - *"What type of data are you trying to extract?"*: `Integer numbers`
>       - *"If multiple annotations are found for the same variant,
>         store ..."*: `the first annotation found`
>    2. Recipe for extracting the dbSNP *CFL* field and adding it as
>       *rs_cfl* to the GEMINI database
>       - *"Elements to extract from the annotation source"*: `CFL`
>       - *"Database column name to use for recording annotations"*:
>         `rs_cfl`
>       - *"What type of data are you trying to extract?"*: `Integer numbers`
>       - *"If multiple annotations are found for the same variant,
>         store ..."*: `the first annotation found`
>    3. Recipe for extracting the dbSNP *ASP* field and adding it as
>       *rs_asp* to the GEMINI database
>       - *"Elements to extract from the annotation source"*: `ASP`
>       - *"Database column name to use for recording annotations"*:
>         `rs_asp`
>       - *"What type of data are you trying to extract?"*: `Integer numbers`
>       - *"If multiple annotations are found for the same variant,
>         store ..."*: `the first annotation found`
>
>    *"Execute"* the job.
> 2. Add information from **cancerhotspots**
>
>    Configure these parameters of **GEMINI annotate** {% icon tool %}:
>    - As *"GEMINI database"* select the dataset generated by the previous run
>      of GEMINI annotate.
>    - *"Dataset to use as the annotation source"* is the
>      `cancerhotspots_v2.bed` dataset imported in the Get Data step of
>      this section.
>    - *"Strict variant-identity matching of database and annotation records
>      (VCF format only)"*: with input in BED format this setting will be
>      ignored
>
>      For the cancerhotspots data, we are simply going to record the best
>      q-value associated with any cancerhotspots variant overlapping one
>      of our variant sites.
>    - *"Type of information to add to the database variants"*: `Specific
>      values extracted from matching records in the annotation source
>      (extract)`
>
>    Now fill in the following *"Annotation extraction recipe"* for extracting
>    the *q-values* of overlapping cancerhotspots sites and adding them as
>    *hs_qvalue* to the GEMINI database:
>       - *"Elements to extract from the annotation source"*: `5`
>
>       The q-values are stored in the fifth column of the BED dataset.
>       - *"Database column name to use for recording annotations"*:
>         `hs_qvalue`
>       - *"What type of data are you trying to extract?"*: `Numbers with
>         decimal precision`
>       - *"If multiple annotations are found for the same variant,
>         store ..."*: `the smallest of the (numeric) values`
>
>    *"Execute"* the job.
> 3. Add links to **CIViC**
>
>    Configure these parameters of **GEMINI annotate** {% icon tool %}:
>    - As *"GEMINI database"* select the dataset generated by the previous run
>      of GEMINI annotate.
>    - *"Dataset to use as the annotation source"* is the
>      `CIViC.bed` dataset imported in the Get Data step of
>      this section.
>    - *"Strict variant-identity matching of database and annotation records
>      (VCF format only)"*: with input in BED format this setting will be
>      ignored
>
>      For the CIViC data, we are going to record the hyperlink to any variant
>      found in the CIViC database that overlaps one of our variant sites.
>    - *"Type of information to add to the database variants"*: `Specific
>      values extracted from matching records in the annotation source
>      (extract)`
>
>    Now fill in the following *"Annotation extraction recipe"* for extracting
>    the hyperlinks of overlapping CIViC sites and adding them as a list of
>    *overlapping_civic_urls* to the GEMINI database:
>       - *"Elements to extract from the annotation source"*: `4`
>
>       The hyperlinks are stored in the fourth column of the BED dataset.
>       - *"Database column name to use for recording annotations"*:
>         `overlapping_civic_url`
>       - *"What type of data are you trying to extract?"*: `Text (text)`
>       - *"If multiple annotations are found for the same variant,
>         store ..."*: `a comma-separated list of non-redundant (text) values`
>
>    *"Execute"* the job.
> 4. Add information from the **Cancer Genome Interpreter (CGI)**
>
>    Configure these parameters of **GEMINI annotate** {% icon tool %}:
>    - As *"GEMINI database"* select the dataset generated by the previous run
>      of GEMINI annotate.
>    - *"Dataset to use as the annotation source"* is the
>      `cgi_variant_positions.bed` dataset imported in the Get Data step of
>      this section.
>    - *"Strict variant-identity matching of database and annotation records
>      (VCF format only)"*: with input in BED format this setting will be
>      ignored
>
>      For the CGI data, we are going to record if any of our variant sites
>      is overlapped by a variant in the CGI biomarkers database.
>    - *"Type of information to add to the database variants"*: `Binary
>      indicator (1=found, 0=not found) of whether the variant had any match in
>      the annotation source (boolean)`
>      - *"Database column name to use for recording annotations"*: `in_cgidb`
>
>    *"Execute"* the job to obtain the final GEMINI database with all desired
>    annotations incorporated into it.
{: .hands_on}

## Reporting selected subsets of variants

Now that we have built our GEMINI database and enriched it with additional
annotations, it is time that we explore the wealth of information stored in it.
The goal is to produce filtered variant reports that list specific classes
(somatic, germline, LOH) of high-quality variants together with their most
relevant annotations. The way to achieve this is through **GEMINI queries**
that specify:

1. filters we want to apply to the variant list stored in the database
2. the pieces of information about the filtered variants that we would like to
   retrieve
   
> ### {% icon comment %} The GEMINI query language
> GEMINI queries are extremely flexible, enabling users to express many
> different ideas to explore the variant space, and as such, the complete
> query syntax can be a bit overwhelming for beginners.
>
> In this section, we will start off with a rather simple query, then build on
> it stepwise before trying a really complex query in the next section.
>
> For a more detailed explanation of the query syntax, you should consult the
> [query section of the GEMINI documentation](https://gemini.readthedocs.io/en/latest/content/querying.html).
> Since the GEMINI query syntax is built on the *SQLite* dialect of *SQL*, the
> *SQLite* documentation, in particular, its chapter on
> [SQLite core functions](https://sqlite.org/lang_corefunc.html), is another
> really helpful resource for understanding more sophisticated queries.
{: .comment}

Lets start by configuring a simple query to obtain a report of *bona fide*
somatic variants:

> ### {% icon hands_on %} Hands-on: Querying the GEMINI database for somatic variants
> Run **GEMINI query** {% icon tool %} with:
> - *"GEMINI database"* set to the fully annotated database created in the last
>   step
> - Select *"Build GEMINI query using"*: `Basic variant query constructor`
>   
> Next, we define a few criteria for variants to get included in the report.
>
> When looking for somatic variants we may want to disregard questionable
> variants, for which either a non-negligible amount of supporting sequencing
> reads is also found in the normal tissue data, or which are only supported by
> a very small fraction of the reads from the tumor sample.
>
> We can build corresponding filter criteria using a *"Genotype filter
> expression"*:
> - Click on {% icon param-repeat %} *"Insert Genotype filter expression"* and
> - Configure: {% icon param-text %} *"Restrictions to apply to genotype
    values"*: `gt_alt_freqs.NORMAL <= 0.05 AND gt_alt_freqs.TUMOR >= 0.10`
>
>   to retain only variants that are supported by less than 5% of the reads of
>   the normal sample, but by more than 10% of the reads of the tumor sample.
>
> Among the info stored in the GEMINI database is the somatic status that
> VarScan somatic has called for every variant. Remember how we used *GEMINI
> annotate* to add this info? Now it is time to use it:
> - Set {% icon param-text %} *"Additional constraints expressed in SQL
>   syntax"*: `somatic_status = 2`
>
>   to retain only those variants passing the genotype filter above, which are
>   **also** considered somatic variants by the variant caller.
>
> Finally, we need to specify what information about the variants that pass our
> filters we would like to have reported.
> - Go to the section *"Output format options"* and select:
>   - *"Type of report to generate"*: `tabular (GEMINI default)
>   - *"Add a header of column names to the output"*: `Yes`
>   - *"Set of columns to include in the variant report table"*: `Custom
>     (report user-specified columns)`
>
>   Then, from *"Choose columns to include in the report"* select:
>   - {% icon param-check %} *"chrom"*
>   - {% icon param-check %} *"start"*
>   - {% icon param-check %} *"ref"*
>   - {% icon param-check %} *"alt"*
>
>   and configure:
>   - *"Additional columns (comma-separated)"*: `gene, aa_change, rs_ids,
>     hs_qvalue, cosmic_ids`
>
> Run the job.
>
> > ### {% icon comment %} How am I supposed to know these column names?
> > Obviously, you need to know the names of the columns (in the tables of the
> > GEMINI database) to include them in the report, but how are you supposed to
> > know them?
> >
> > The standard ones (added by *GEMINI load* when building the database) are
> > listed in the
> > [GEMINI documentation](https://gemini.readthedocs.io/en/latest/content/database_schema.html).
> > The non-standard columns (the ones you added with *GEMINI annotate*) have
> > the names you gave them, when you added them.
> >
> > Alternatively, to get the tables and column names of a specific database
> > listed, you can use **GEMINI database info** {% icon tool %} like so:
> >
> > - *"GEMINI database"*: the database you want to explore
> > - *"Information to retrieve from the database"*: `Names of database tables
> > and their columns`
> {: .comment}
{: .hands_on}

> ### {% icon hands_on %} Hands-on: More complex filter criteria
> Run **GEMINI query** {% icon tool %} with the exact same settings as before,
> but:
> - Set {% icon param-text %} *"Additional constraints expressed in SQL
>   syntax"*: `somatic_status = 2 AND somatic_p <= 0.05 AND (filter IS NULL OR rs_ids IS NOT NULL) AND rs_cfl != 1 and rs_asp != 1`
>
>   This translates into "variants classified as somatic with a p-value <=
>   0.05, which haven't been flagged as likely false-positives or, if so, are
>   listed in dbSNP, but in there, are not flagged as being assembly-dependent
>   or -specific".
{: .hands_on}

If you've followed all steps up to here exactly, running this job should give
you a tabular dataset of 43 variants, and with the annotations in the report
it is relatively easy to pick out a few interesting ones.
Before we focus on the content of the report, however, we could enhance the
report format a bit more.

> ### {% icon hands_on %} Hands-on: SQL-based output formatting
> Run **GEMINI query** {% icon tool %} with the exact same settings as in the
> last example, but, in the *"Output format options"* section:
> - Set {% icon param-text %} *"Additional columns (comma-separated)"*:
>   `type, gt_alt_freqs.TUMOR, gt_alt_freqs.NORMAL,
>   ifnull(nullif(round(max_aaf_all,2),-1.0),0) as MAF, gene, impact_so,
>   aa_change, ifnull(round(cadd_scaled,2),'.') as cadd_scaled,
>   round(gerp_bp_score,2) as gerp_bp, ifnull(round(gerp_element_pval,2),'.')
>   as gerp_element_pval, ifnull(round(hs_qvalue,2), '.') as hs_qvalue,
>   in_omim, ifnull(clinvar_sig,'.') as clinvar_sig,
>   ifnull(clinvar_disease_name,'.') as clinvar_diesease_name,
>   ifnull(rs_ids,'.') as dbsnp_ids, rs_ss, ifnull(cosmic_ids,'.') as
>   cosmic_ids, ifnull(overlapping_civic_url,'.') as overlapping_civic_url,
>   in_cgidb`
>
> This adds a lot more annotations to the report and it also uses some
> [SQLite functions](https://sqlite.org/lang_corefunc.html) to clean up the
> output.
>
> Compare the new report to the previous one to see what has changed.
>
> > ### {% icon details %} Limitations of genotype column queries
> > In case you are wondering why the above query does not use rounding on
> > the alternate allele frequencies of the samples, *i.e.*, on
> > `gt_alt_freqs.TUMOR` and `gt_alt_freqs.NORMAL`, or why it does not rename
> > these columns, that is because it would break the query.
> >
> > As a general rule, note that all columns in the variants table starting
> > with `gt` (the genotype columns, calculated from the genotype fields of a
> > VCF dataset) cannot be used like regular SQLite columns, but are parsed by
> > GEMINI separately. That is why you cannot mix them with regular SQLite
> > elements like functions and alias specifications with `AS`.
> > You may also have noticed that `gt_alt_freqs.TUMOR` and
> > `gt_alt_freqs.NORMAL` do not obey the column order specification of the
> > query, but end up as the last columns of the tabular report. This is
> > another artefact of GEMINI's special treatment of `gt` columns.
> {: .details}
{: .hands_on}

## Generating reports of genes affected by variants

As a final step, let us now try to generate a gene-centered report based on the
same somatic variants we just selected above.

Such a gene- centered report would include annotations that apply to a whole
gene affected by a variant rather than to the variant itself. Examples of such
annotations include known synonyms of an affected gene, its NCBI entrez number,
the ClinVar phenotype, if any, associated with the gene, a hyperlink to the
gene's page at CIViC.org, *etc.*.#

Some of this information comes built-in into every GEMINI database, but it is
stored in a separate table called `gene_detailed`, while all information we
used and queried so far was from the `variants` table.

To access information from the `variants` and the `gene_detailed` table in the
same query we need to join the two tables. Such an operation is not possible
with the `Basic query constructor` we have used so far, but requires an
advanced mode for composing the query.

> ### {% icon hands_on %} Hands-on: Turning query results into gene-centered reports
> Run **GEMINI query** {% icon tool %} in advanced mode by choosing
> - *"Build GEMINI query using"*: `Advanced query constructor`
>
>   Then use the settings:
>   - *"The query to be issued to the database"*: `SELECT v.gene, v.chrom,
>     g.synonym, g.hgnc_id, g.entrez_id, g.rvis_pct, v.clinvar_gene_phenotype
>     FROM variants v, gene_detailed g WHERE v.chrom = g.chrom AND
>     v.gene = g.gene AND v.somatic_status = 2 AND v.somatic_p <= 0.05 AND
>     (v.filter IS NULL OR v.rs_ids IS NOT NULL) AND v.rs_cfl != 1 and
>     v.rs_asp != 1 GROUP BY g.gene`
>
> > ### {% icon comment %} Elements of the SQL query
> > In this full SQL query, the part between `SELECT` and `FROM` states which
> > columns from which database tables we wish to retrieve, while the part
> > between `FROM` and `WHERE` specifies the database tables that need to be
> > consulted and gives them simpler aliases (`v` becomes an alias for
> > the `variants` table, `g` for the `gene_detailed` table), which we can then
> > use everywhere else in the query.
> >
> > The part following `WHERE` are the filter criteria we want to apply. Note
> > that these criteria are almost the same as those we have used in our
> > earlier somatic variants query, but since we are working with two tables
> > instead of just one, we need to state which table the filter columns come
> > from through table prefixes. Thus, `somatic_status` becomes
> > `v.somatic_status`, *etc.*. In addition, we want to report, of course,
> > corresponding information from the two tables, which is ensured by the
> > additional criteria `v.chrom = g.chrom` and `v.gene = g.gene` (the SQL
> > terminology for this is: we want to join the `variants` and the
> > `gene_detailed` tables on their `chrom` and `gene` columns).
> >
> > The `GROUP BY` part, finally, specifies that we want to collapse records
> > affecting the same gene into one.
> {: .comment}
>
>  - *"Genotype filter expression"*: `gt_alt_freqs.NORMAL <= 0.05 AND
>    gt_alt_freqs.TUMOR >= 0.10`
>
>    This remains the same as in the previous somatic variants query.    
{: .hands_on}

## Adding additional annotations to the gene-centered report

Unfortunately, *GEMINI annotate* lets you add columns only to the variants
table of a GEMINI database, but there is no simple way to enhance the
`gene_detailed` table with additional annotations. That's why we are going to
add such extra annotations now to the tabular gene-centered report using more
general-purpose Galaxy tools.

> ### {% icon hands_on %} Hands-on: Join and rearrange to get a fully annotated gene report
> Annotating the tabular gene report produced with GEMINI is a two-step
> process, in which we first *join* the report and tabular annotation sources
> into a larger tabular dataset, from which we then eliminate redundant and
> unwanted columns, while *rearranging* the remaining ones.
>
> **Step 1** consists of three separate *join* operations that sequentially
> pull in the annotations found in the three gene-based tabular datasets
> obtained in the *Get Data* step of this section.
> 1. Adding **UniProt cancer genes** information
>
>    Use **Join two files** {% icon tool %} with these settings:
>    - *"1st file"*: the GEMINI-generated gene report from the previous step
>    - *"Column to use from 1st file"*: `Column: 1`
>    - *"2nd file"*: the `Uniprot_Cancer_Genes` dataset imported in the *Get
>      Data* step
>    - *"Column to use from 2nd file"*: `Column: 1`
>
>    Together these settings indicate that we want to join the two selected
>    datasets using the first column (the gene column) of both to determine
>    corresponding lines.
>
>    Next we need to indicate which lines from each dataset should be retained
>    in the output. Configure:
>    - *"Output lines appearing in"*: `Both 1st & 2nd file, plus unpairable
>      lines from 1st file. (-a 1)`
>
>      If a variant-affected gene is not listed as a Uniprot cancer gene, then,
>      obviously, we still want to have it included in the final report.
>    - *"First line is a header line"*: `Yes`
>    - *"Ignore case"*: `No`
>    - *"Value to put in unpaired (empty) fields"*: `0`
>
>      If you inspect the `Uniprot_Cancer_Genes` dataset, you will see that it
>      has two annotation columns - one indicating, using `1` and `0`, whether
>      a given gene is a *proto-oncogene* or not, the other one indicating
>      *tumor suppressor* genes the same way. For genes missing from this
>      annotation dataset, we want to fill the corresponding two columns of the
>      join result with `0` to indicate the common case that a gene affected by
>      a variant is neither a known proto-oncogene, nor a tumor suppressor
>      gene.
>
>    - *"Execute"* the job and inspect the result
>
> 2. Adding **CGI biomarkers** status column
>
>    Use **Join two files** {% icon tool %} with these settings:
>    - *"1st file"*: the partially annotated dataset from step 1
>    - *"Column to use from 1st file"*: `Column: 1`
>    - *"2nd file"*: the `cgi_genes` dataset imported in the *Get Data* step
>    - *"Column to use from 2nd file"*: `Column: 1`
>    - *"Output lines appearing in"*: `Both 1st & 2nd file, plus unpairable
>      lines from 1st file. (-a 1)`
>    - *"First line is a header line"*: `Yes`
>    - *"Ignore case"*: `No`
>    - *"Value to put in unpaired (empty) fields"*: `0`
>
>    - *"Execute"* the job and inspect the input and the result dataset to
>      make sure you understand what happened at this step.
>
> 3. Adding **CIViC gene summaries**
>
>    Use **Join two files** {% icon tool %} with these settings:
>    - *"1st file"*: the partially annotated dataset from step 2
>    - *"Column to use from 1st file"*: `Column: 1`
>    - *"2nd file"*: the `GeneSummaries` dataset imported in the *Get Data*
>      step
>    - *"Column to use from 2nd file"*: `Column: 3`
>
>      The gene column in the CIViC gene summaries annotation dataset is *not*
>      the first one!
>    - *"Output lines appearing in"*: `Both 1st & 2nd file, plus unpairable
>      lines from 1st file. (-a 1)`
>    - *"First line is a header line"*: `Yes`
>    - *"Ignore case"*: `No`
>    - *"Value to put in unpaired (empty) fields"*: `.`
>
>    - *"Execute"* the job and inspect the input and the result dataset to
>      make sure you understand what happened at this step.
>
> If you inspected all output datasets as suggested, you will have noticed
> that each of the *join* operations kept the gene columns from both input
> datasets. In addition, we had no control over the order, in which columns got
> added to the report, nor could we exclude columns.
>
> In **Step 2** of the report preparation we are going to address all of these
> issues.
>
> 4. Rearranging and filtering report columns
>
>    Use **Column arrange by header name** {% icon tool %} configured like
>    this:
>    - *"file to rearrange"*: the final *join* result dataset from step 3
>    - set up **12** {% icon param-repeat %} *"Specify the first few coulmns by
>      name"* elements and
>    - set them to `gene`, `chrom`, `synonym`, `hgnc_id`, `entrez_id`,
>      `rvis_pct`, `is_OG`, `is_TS`, `in_cgi_biomarkers`,
>      `clinvar_gene_phenotype`, `gene_civic_url`, `description` in this order
>    - set *"Discard unspecified columns"*: `Yes`
>
> > ### {% icon comment %} Alternative tool suggestion
> > If your Galaxy server does not offer the *Column arrange* tool, it will
> > almost certainly offer the **Cut columns from a table** {% icon tool %},
> > which can be used as a drop-in replacement.
> > Instead of column names, however, this tool expects a comma-separated list
> > of column indexes, like `c1,c2` for the first and second column, so you
> > will have to first figure out the column numbers of your columns of
> > interest.
> >
> > **Note**: Be sure not to confuse the suggested tool with *Cut
> > columns from a table (cut)*, which does *not* let you change the order of
> > columns!
> {: .comment}
{: .hands_on}

# Conclusion
{:.no_toc}

In addition to merely calling variants, *somatic variant calling* tries to
distinguish *somatic mutations*, which are private to tumor tissue, from
*germline* mutations, that are shared by tumor and healthy tissue, and
*loss-of-heterozygosity* events, which involve the loss, from tumor tissue, of
one of two alleles found at a biallelic site in healthy tissue.

Dedicated somatic variant callers can perform this classification on
statistical grounds, but the interpretation of any list of variants (somatic,
germline or LOH) also depends crucially on rich genetic and cancer-specific
variant and gene annotations.

