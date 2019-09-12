---
layout: tutorial_hands_on

title: Identification of somatic and germline variants from tumor and normal sample
  pairs
zenodo_link: https://doi.org/10.5281/zenodo.2582555
questions:
- What are the specific challenges in somatic variant calling that set it apart from
  regular diploid variant calling?
- How can you call variants and classify them according to their presence/absence
  in/from tumor and normal tissue of the same individual?
- How can you annotate variants and affected genes with prior knowledge from human
  genetic and cancer-specific databases to generate clinically relevant reports?
objectives:
- Call variants and their somatic status from whole-exome sequencing data
- Annotate variants with a wealth of human genetic and cancer-specific information
  extracted from public databases
- Add gene-level annotations and generate reports of annotated somatic and germline
  variants, loss-of-heterozygosity (LOH) events, and affected genes, ready for interpretation
  by clinicians
time_estimation: 7h
key_points:
- Follow best practices for read mapping, quality control and mapped reads postprocessing
  to minimze false-positive variant calls.
- Use a dedicated somatic variant caller to call variants and to classify them into
  somatic, germline and LOH event variants on statistical grounds.
- Annotations and queries based on variant properties add relevance to variant and
  gene reports.
- A framework like GEMINI is very helpful for managing, annotating and querying lists
  of variants in a flexible way.
- Prefer public, free annotation sources to foster reproducibility and information
  sharing.
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
For the sake of an acceptable speed of the analysis, the original data has been
downsampled though to include only the reads from human chromosomes 5, 12 and
17.

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
> 3. Check that all newly created datasets in your history have their datatypes assigned 
>    correctly, and fix any missing or wrong datatype assignment
>
>    {% include snippets/change_datatype.md datatype="fastqsanger.gz" %}
>
> 4. Rename the datasets and add appropriate tags to them
>
>    For datasets that you upload via a link, Galaxy will pick the link
>    address as the dataset name, which you will likely want to shorten to
>    just the file names.
>
>    {% include snippets/rename_dataset.md %}
>
>    Large parts of the analysis in this tutorial will consist of
>    identical steps performed on the normal and on the tumor tissue data
>    in parallel.
>
>    To make it easier to keep track of which dataset represents which
>    branch of an analysis in a linear history, Galaxy supports dataset
>    tags. In particular, if you attach a tag starting with `#` to any
>    dataset, that tag will automatically propagate to any new dataset
>    derived from the tagged dataset.
>
>    {% include snippets/add_tag.md %}
>
>    Before starting our analysis it is, thus, a good idea to tag the two
>    fastq datasets representing the normal tissue with, *e.g.*, `#normal`
>    and the two datasets representing the tumor tissue with, *e.g.*,
>    `#tumor`.
>
> 5. Import the reference genome with the `hg19` version of the sequences of
>    human chromosomes 5, 12 and 17:
>
>    ```
>    https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz
>    ```
>
>    Make sure you specify the datatype as `fasta` in the import dialog.
>
>    > ### {% icon comment %} Shortcut
>    > You can skip this step if the Galaxy server you are working on offers
>    > a `hg19` version of the human reference genome with prebuilt indexes for
>    > *bwa-mem* **and** *samtools* (ask your instructor, or check the tools
>    > **Map with BWA-MEM** {% icon tool %} and **VarScan Somatic**
>    > {% icon tool %} if they list a `hg19` version as an option under
>    > *"(Using) reference genome"*).
>    {: .comment}
>
>    Alternatively, load the dataset from a shared data library.
>
> 6. Rename the reference genome
>
>    The reference genome you have imported above came as a compressed
>    file, but got unpacked by Galaxy to plain `fasta` format according to
>    your datatype selection. You may now wish to remove the `.gz` suffix
>    from the dataset name.
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
> [Quality Control]({{ site.baseurl }}{% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %})
> and [Mapping]({{ site.baseurl }}{% link topics/sequence-analysis/tutorials/mapping/tutorial.md %})
> tutorials.
> Here, we will only illustrate the concrete steps necessary for quality
> control and read mapping of our particular datasets.
{: .comment}

## Quality control

> ### {% icon hands_on %} Hands-on: Quality control of the input datasets
> 1. Run **FastQC** {% icon tool %} on each of your four fastq datasets
>       - {% icon param-files %} *"Short read data from your current history"*: all 4 FASTQ  datasets selected with **Multiple datasets**
>
>    {% include snippets/select_multiple_datasets.md %}
>
>    When you start this job, eight new datasets (one with the calculated raw
>    data, another one with an html report of the findings for each input
>    dataset) will get added to your history.
>
> 2. Use **MultiQC** {% icon tool %} to aggregate the raw **FastQC** data of all four input datasets into one report
>      - In *"Results"*
>        - *"Which tool was used generate logs?"*: `FastQC`
>        - In *"FastQC output"*
>           - *"Type of FastQC output?"*: `Raw data`
>           - {% icon param-files %} *"FastQC output"*: all four *RawData*
>             outputs of **FastQC** {% icon tool %})
>
> 3. Inspect the *Webpage* output produced by the tool
>
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
>    > >    from a sample is considered to hint at possible contamination.
>    > >    Here, however, we are dealing with sequencing data from captured
>    > >    exomes, *i.e*, the reads are not representing random sequences from
>    > >    a genome, but rather an arbitrary selection.
>    > >    In fact, the samples at hand were prepared using Agilent's
>    > >    SureSelect V5 technology for exome enrichment, and bimodal GC
>    > >    content distributions can be seen as a hallmark of that capture
>    > >    method in several publications (see, for example, Fig. 4C in
>    > >    [Meienberg et al., 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4477645/)
>    > >    ).
>    > >
>    > {: .solution}
>    {: .question}
{: .hands_on}

## Read trimming and filtering
Although the raw reads used in this tutorial are of relatively good overall
quality already, we will apply read trimming and filtering to see if we can
improve things still a bit, but also to demonstrate the general concept.

> ### {% icon hands_on %} Hands-on: Read trimming and filtering of the normal tissue reads
> 1. **Trimmomatic** {% icon tool %} to trim and filter the normal tissue reads
>    - *"Single-end or paired-end reads?"*: `Paired-end (two separate
>      input files)`
>      
>      This makes the tool treat the forward and reverse reads simultaneously.
>
>      - {% icon param-file %} *"Input FASTQ file (R1/first of pair)"*: the
>        forward reads (r1) dataset of the normal tissue sample
>      - {% icon param-file %} *"Input FASTQ file (R2/second of pair)"*: the
>        reverse reads (r2) dataset of the normal tissue sample
>
>    - *"Perform initial ILLUMINACLIP step?"*: `Yes`
>       - *"Select standard adapter sequences or provide custom?"*: `Standard`
>          - *"Adapter sequences to use"*: `TruSeq3 (paired-ended, for MiSeq and
>          HiSeq)`
>       - *"Maximum mismatch count which will still allow a full match to be
>        performed"*: `2`
>       - *"How accurate the match between the two 'adapter ligated' reads must
>        be for PE palindrome read alignment"*: `30`
>       - *"How accurate the match between any adapter etc. sequence must be
>        against a read"*: `10`
>       - *"Minimum length of adapter that needs to be detected (PE specific/
>        palindrome mode)"*: `8`
>       - *"Always keep both reads (PE specific/palindrome mode)?"*: `Yes`
>
>       These parameters are used to cut ILLUMINA-specific adapter sequences
>       from the reads.
>
>    - In *"Trimmomatic Operation"*
>       - In *"1: Trimmomatic Operation"*
>           -  *"Select Trimmomatic operation to perform"*: `Cut the specified number
>       of bases from the start of the read (HEADCROP)`
>              - *"Number of bases to remove from the start of the read"*: `3`
>       - {% icon param-repeat %} "Insert Trimmomatic Operation"*
>       - In *"2: Trimmomatic Operation"*
>           -  *"Select Trimmomatic operation to perform"*: `Cut bases off the end of
>       a read, if below a threshold quality (TRAILING)`
>              - *"Minimum quality required to keep a base"*: `10`
>       - {% icon param-repeat %} "Insert Trimmomatic Operation"*
>       - In *"3: Trimmomatic Operation"*
>           - *"Select Trimmomatic operation to perform"*: `Drop reads below a
>       specified length (MINLEN)`
>              - *"Minimum quality required to keep a base"*: `25`
>
>       These three trimming and filtering operations will be applied to the
>       reads in the given order after adapter trimming.
>    
{: .hands_on}

Running this job will generate four output datasets:
- two for the trimmed forward and reverse reads that still have a proper
  mate in the other dataset
- two more datasets of orphaned forward and reverse reads, for which the
  corresponding mate got dropped because of insufficient length after
  trimming; when you inspect these two files, however, you should find
  that they are empty because none of our relatively high quality reads
  got trimmed that excessively. You can delete the two datasets to keep
  your history more compact.

> ### {% icon details %} More on handling unpaired reads in paired-end data
> Splitting out potential unpaired reads into separate datasets like this
> is important because read mappers, typically, expect reads in forward
> and reverse input datasets to be arranged in proper pairs, and reads
> in one of the datasets without a counterpart in the other would destroy
> that expected structure. Therefor, when your data is paired-end data
> always make sure you use Trimmomatic in paired-end mode, too!
{: .details}

> ### {% icon hands_on %} Hands-on: Read trimming and filtering of the tumor tissue reads
> 2. Trim and filter the **tumor tissue** reads following the same steps as above, just change the two input datasets to treat the tumor tissue reads with identical settings.
>
> 3. Check that the two unpaired reads datasets are empty, and delete them.   
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
> 1. **Map with BWA-MEM** {% icon tool %} to map the reads from the **normal tissue** sample to the reference genome
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a built-in genome index`
>        - *"Using reference genome"*: `Human: hg19` (or a similarly named option)
> 
>      > ### {% icon comment %} Using the imported `hg19` sequence
>      > If you have imported the `hg19` sequence as a fasta dataset into your
>      > history instead:
>      >   - *"Will you select a reference genome from your history or use a
>      >     built-in index?"*: `Use a genome from history and build index`
>      >      - {% icon param-file %} *"Use the following dataset as the reference sequence"*: your imported `hg19` fasta dataset.
>      {: .comment}
>
>    - *"Single or Paired-end reads"*: `Paired`
>       - {% icon param-file %} *"Select first set of reads"*: the trimmed 
>         forward reads (r1) dataset of the **normal tissue** sample; output of
>         **Trimmomatic** {% icon tool %}
>       - {% icon param-file %} *"Select second set of reads"*: the trimmed
>         reverse reads (r2) dataset of the **normal tissue** sample; output of
>         **Trimmomatic** {% icon tool %}
>    - *"Set read groups information?"*: `Set read groups (SAM/BAM specification)`
>      - *"Auto-assign"*: `No`
>        - *"Read group identifier (ID)"*: `231335` (this value being taken
>          from the original name of the normal tissue input files)
>      - *"Auto-assign"*: `No`
>        - *"Read group sample name (SM)"*: `Normal`
>
>    > ### {% icon comment %} More on read group identifiers and sample names
>    > In general, you are free to choose ID and SM values to your liking,
>    > but the ID should unambiguously identify the sequencing run that
>    > produced the reads, while the SM value should identify the
>    > biological sample.
>    {: .comment}
>
> 2. **Map with BWA-MEM** {% icon tool %} to map the reads from the **tumor tissue** sample,
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a built-in genome index`
>        - *"Using reference genome"*: `Human: hg19` (or a similarly named option)
>
>      Adjust these settings as before if you are using the imported reference
>      genome.
>    - *"Single or Paired-end reads"*: `Paired`
>       - {% icon param-file %} *"Select first set of reads"*: the trimmed
>         forward reads (r1) dataset of the **tumor tissue** sample; output of
>         **Trimmomatic** {% icon tool %}
>       - {% icon param-file %} *"Select second set of reads"*: the reverse
>         reads (r2) dataset of the **tumor tissue** sample; output of
>         **Trimmomatic** {% icon tool %}
>    - *"Set read groups information?"*: `Set read groups (SAM/BAM specification)`
>      - *"Auto-assign"*: `No`
>        - *"Read group identifier (ID)"*: `231336` (this value, again, being
>          taken from the original name of the tumor tissue input files)
>      - *"Auto-assign"*: `No`
>        - *"Read group sample name (SM)"*: `Tumor`
>
{: .hands_on}

# Mapped reads postprocessing

To ensure that we base our variant analysis only on unambiguous, high-quality
read mappings we will do some postprocessing next.

## Filtering on mapped reads properties

To produce new filtered BAM datasets with only those reads retained that
have a minimal mapping quality of 1 and are mapped in a proper pair:

> ### {% icon hands_on %} Hands-on: Filtering for mapping status and quality
>
> 1. **Filter SAM or BAM, output SAM or BAM** {% icon tool %} with the following
> parameters (leaving non-mentioned ones at their defaults):
>   - {% icon param-files %} *"SAM or BAM file to filter"*: mapped reads datasets from
>     the normal *and* the tumor tissue data, outputs of **Map with BWA-MEM** {% icon tool %}
>   - *"Minimum MAPQ quality score"*: `1`
>   - *"Filter on bitwise flag"*: `yes`
>     - *"Only output alignments with all of these flag bits set"*:
>        - {% icon param-check %} *"Read is mapped in a proper pair"*
>
{: .hands_on}

This will result in two new datasets, one for each of the normal and tumor data.

## Removing duplicate reads

> ### {% icon hands_on %} Hands-on: Remove duplicates
>
> 1. **RmDup** {% icon tool %} with the following parameters:
>   - {% icon param-files %} *"BAM file"*: filtered reads datasets from
>     the normal *and* the tumor tissue data; the outputs of
>     **Filter SAM or BAM**
>   - *"Is this paired-end or single end data"*: `BAM is paired-end`
>     - *"Treat as single-end"*: `No`
>
{: .hands_on}

Again, this will produce two new datasets, one for each of the normal and
tumor data.

## Left-align reads around indels

> ### {% icon hands_on %} Hands-on: Left-align
>
> 1. **BamLeftAlign** {% icon tool %} with the following parameters:
>    - *"Choose the source for the reference genome"*: `Locally cached`
>      - {% icon param-files %} *"BAM dataset to re-align"*: your
>        filtered and deduplicated reads datasets from the normal *and* the tumor
>        tissue data; the outputs of **RmDup**
>      - *"Using reference genome"*: `Human: hg19` (or a similarly
>        named choice)
>   - *"Maximum number of iterations"*: `5`
>
> > ### {% icon comment %} Using the imported `hg19` sequence
> > If you have imported the `hg19` sequence as a fasta dataset into your
> > history instead:
> >   - *"Choose the source for the reference genome"*: `History`
> >      - {% icon param-file %} *"Using reference file"*: your imported `hg19` fasta dataset
> {: .comment}
{: .hands_on}

As before, this will generate two new datasets, one for each of the normal and tumor data.

## Recalibrate read mapping qualities

> ### {% icon hands_on %} Hands-on: Recalibrate read quality scores
>
> 1. **CalMD** {% icon tool %} with the following parameters:
>    - {% icon param-files %} *"BAM file to recalculate"*: the left-aligned
>      datasets from the normal and the tumor tissue data; the outputs of
>      **BamLeftAlign** {% icon tool %}
>    - *"Choose the source for the reference genome"*: `Use a built-in genome`
>      - *"Using reference genome"*: `Human: hg19` (or a similarly named choice)
>
>      > ### {% icon comment %} Using the imported `hg19` sequence
>      > If you have imported the `hg19` sequence as a fasta dataset into your
>      > history instead:
>      >   - *"Choose the source for the reference genome"*: `Use a genome from the history`
>      >      - {% icon param-file %} *"Using reference file"*: your imported `hg19` fasta dataset.
>      {: .comment}
>
>    - *"Do you also want BAQ (Base Alignment Quality) scores to be calculated?"*: `No`
>
>      The *VarScan somatic* tool that we are going to use for calling variants at
>      the next step is typically used in combination with unadjusted base quality
>      scores because the general opinion is that the base quality downgrades
>      performed by *CalMD* and other tools from the *samtools* suite of tools
>      are too severe for *VarScan* to retain good sensitivity. We are sticking
>      with this practice in this tutorial.
>
>      > ### {% icon comment %} Using adjusted base quality scores
>      > If, for your own data, you would like to experiment with adjusted base
>      > quality scores, it is important to understand that *VarScan somatic* will
>      > only make use of the adjusted scores if they are incorporated directly into
>      > the read base qualities of a BAM input dataset, but not if they are written
>      > to the dataset separately. 
>      > 
>      > Hence, should you ever decide to use:
>      > - *"Do you also want BAQ (Base Alignment Quality) scores to be calculated?"*: `Yes, run BAQ calculation`
>      >
>      >   and you want this setting to affect downstream variant calling with
>      >   *VarScan somatic*, make sure you also set then:
>      >   - *"Use BAQ to cap read base qualities"*: `Yes`
>      >
>      > Please also note that BAQ scores are quite expensive to calculate so
>      > be prepared to see a substantial (up to 10x!) increase in job run time 
>      > when enabling it.
>      {: .comment}
>
>    - *"Additional options"*: `Advanced options`
>      - *"Change identical bases to '='"*: `No`
>      - *"Coefficient to cap mapping quality of poorly mapped reads"*: `50`
>
>        This last setting is the real reason why we use CalMD at this point.
>        It is an empirical, but well-established finding that the mapping
>        quality of reads mapped with *bwa* should be capped this way before
>        variant calling.
>
{: .hands_on}

This will, once more, produce two new datasets, one for each of the normal
and tumor data.


# Variant calling and classification

Having generated a high-quality set of mapped read pairs, we can proceed to
variant calling. The tool **VarScan somatic** is a dedicated solution for
somatic variant calling that:

- detects variant alleles in tumor/normal sample pairs
- calls sample genotypes at variant sites
- classifies variants into germline, somatic and LOH event variants using
  solid classical statistics even in the presence of non-pure samples like
  those obtained from biopsies

> ### {% icon hands_on %} Hands-on: Variant calling and classification
>
> 1. **VarScan somatic** {% icon tool %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in
>   genome?"*: `Use a built-in genome`
>      - *"reference genome"*: `Human: hg19` (or a similarly named choice)
>
>      > ### {% icon comment %} Using the imported `hg19` sequence
>      > If you have imported the `hg19` sequence as a fasta dataset into your
>      > history instead:
>      >   - *"Will you select a reference genome from your history or use a built-in
>      >     genome?"*: `Use a genome from the history`
>      >      - {% icon param-file %} *"reference genome"*: your imported `hg19` fasta dataset.
>      {: .comment}
>
>    - {% icon param-file %} *"aligned reads from normal sample"*: the mapped
>      and fully post-processed normal tissue dataset; one of the two outputs
>      of **CalMD** {% icon tool %}
>    - {% icon param-file %}*"aligned reads from tumor sample"*: the mapped
>      and fully post-processed tumor tissue dataset; the other output of
>      **CalMD** {% icon tool %}
>    - *"Estimated purity (non-tumor content) of normal sample"*: `1`
>    - *"Estimated purity (tumor content) of tumor sample"*: `0.5`
>    - *"Generate separate output datasets for SNP and indel calls?"*: `No`
>    - *"Settings for Variant Calling"*: `Customize settings`
>      - *"Minimum base quality"*: `28`
>
>        We have seen, at the quality control step, that our sequencing data is of
>        really good quality, and we have chosen not to downgrade base qualities at
>        the quality scores recalibration step above, so we can increase the base
>        quality required at any given position without throwing away too much of
>        our data. 
>
>      - *"Minimum mapping quality"*: `1`
>
>        During postprocessing, we have filtered our reads for ones with a mapping
>        quality of at least one, so requiring this quality also here does not
>        actually change the results, but it makes the requirement more explicit.
>
>      Leave all other settings in this section at their default values.
>    - *"Settings for Posterior Variant Filtering"*: `Use default values`
{: .hands_on}

# Variant annotation and reporting

For this tutorial we are going to use variant and gene annotations from many
different sources. Most of these are handled for us by the tools we are going
to use in this section, but we need to take care of importing the data from
four sources into Galaxy separately:

- variant annotations from [Cancer Hotspots](https://www.cancerhotspots.org)
- variant and gene information from the
  [Cancer Biomarkers database](https://www.cancergenomeinterpreter.org/biomarkers)
  of the Cancer Genome Interpreter (CGI) project
- variant and gene information from the [CIViC](https://civicdb.org) database
- variant annotations from [dbSNP](https://www.ncbi.nlm.nih.gov/snp)
- lists of genes annotated with the keywords *proto-oncogene* or *tumor
  suppressor* at [UniProt](https://www.uniprot.org/)

Each of these annotation sets has been released either in the public domain or
under a free data license, which allows you to use them as part of this
tutorial, but also for other purposes.

Starting from the data downloaded from each of these sites, we have generated
a set of new data files tailored to the requirements of the workflow of this
tutorial and have made them available through Zenodo, again under a free data
license.

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
>
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
>
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
> 1. **SnpEff eff** {% icon tool %}
>    - {% icon param-file %} *"Sequence changes (SNPs, MNPs, InDels)"*: the
>      output of **VarScan somatic** {% icon tool %}
>    - *"Input format"*: `VCF`
>    - *"Output format"*: `VCF (only if input is VCF)`
>    - *"Genome source"*: `Locally installed reference genome`
>       - *"Genome"*: `Homo sapiens: hg19` (or a similarly named option)
>
>      > ### {% icon comment %} Using the imported `hg19` SnpEff genome database
>      > If you have imported the `hg19` SnpEff genome database into your
>      > history instead:
>      >   - *"Genome source"*: `Downloaded snpEff database in your history`
>      >      - {% icon param-file %} *"SnpEff4.3 Genome Data"*: your imported `hg19` SnpEff dataset.
>      {: .comment}
>
>    - *"Produce Summary Stats"*: `No`
> 
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
> 1. **GEMINI load** {% icon tool %} with
>    - {% icon param-file %} *"VCF dataset to be loaded in the GEMINI database"*:
>      the output of **SnpEff eff** {% icon tool %}
>    - *"The variants in this input are"*: `annotated with snpEff`
>    - *"This input comes with genotype calls for its samples"*: `Yes`
>
>      Calling sample genotypes was part of what we used VarScan somatic for.
>
>    - *"Choose a gemini annotation source"*: select the latest available annotations snapshot (most likely, there will be only one)
>    - *"Sample and family information in PED format"*: `Nothing selected`
>    - *"Load the following optional content into the database"* 
>      - {% icon param-check %} *"GERP scores"*
>      - {% icon param-check %} *"CADD scores"*
>      - {% icon param-check %} *"Gene tables"*
>      - {% icon param-check %} *"Sample genotypes"*
>      - {% icon param-check %} *"variant INFO field"*
>
>      Be careful to leave **unchecked**:
>      - *"Genotype likelihoods (sample PLs)"*
>
>         VarScan somatic does not generate these values
>
>      - *"only variants that passed all filters"*
>
>        It is simple and more flexible to filter for this property later
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
> 1. **GEMINI annotate** {% icon tool %}:
>    - {% icon param-file %} *"GEMINI database"*: output of **GEMINI load**
>    - {% icon param-file %} *"Dataset to use as the annotation source"*: output of **VarScan somatic**
>    - *"Strict variant-identity matching of database and annotation records
>   (VCF format only)"*: `Yes`
>
>      This setting does not really matter here since you have built the GEMINI
>      database from the exact same list of variants that we are now retrieving
>      annotations from and because VarScan somatic does not call multiple alleles
>      at single sites. Matching on variant-identity is the behavior we would
>      like to see though, so we may as well be explicit about it.
>
>    - *"Type of information to add to the database variants"*: `Specific values
>      extracted from matching records in the annotation source (extract)`
>
>      - In *"1: Annotation extraction recipe"*:
>        - *"Elements to extract from the annotation source"*: `SS`
>        - *"Database column name to use for recording annotations"*:
>        `somatic_status`
>        - *"What type of data are you trying to extract?"*: `Integer numbers`
>        - *"If multiple annotations are found for the same variant,
>        store ..."*: `the first annotation found`
>
>        This is the recipe for extracting the VarScan-generated *SS* field and
>        adding it as a new column *somatic_status* to the GEMINI database.
>
>      - {% icon param-repeat %} *"Insert Annotation extraction recipe"*
>      - In *"2: Annotation extraction recipe"*:
>        - *"Elements to extract from the annotation source"*: `GPV`
>        - *"Database column name to use for recording annotations"*:
>       `germline_p`
>        - *"What type of data are you trying to extract?"*: `Numbers with
>        decimal precision`
>        - *"If multiple annotations are found for the same variant,
>        store ..."*: `the first annotation found`
>
>        This is the recipe for extracting the VarScan-generated *GPV* field
>        and adding it as a new column *germline_p* to the GEMINI database.
>
>      - {% icon param-repeat %} *"Insert Annotation extraction recipe"*
>      - In *"3: Annotation extraction recipe"*:
>        - *"Elements to extract from the annotation source"*: `SPV`
>        - *"Database column name to use for recording annotations"*:
>       `somatic_p`
>        - *"What type of data are you trying to extract?"*: `Numbers with
>        decimal precision`
>        - *"If multiple annotations are found for the same variant,
>        store ..."*: `the first annotation found`
>
>        This is the recipe for extracting the VarScan-generated *SPV* field
>        and adding it as a new column *somatic_p* to the GEMINI database.
>
{: .hands_on}

Next, we are going to add additional annotations beyond the ones directly
obtainable through *GEMINI* or from the input VCF dataset. Specifically we
want to add:

- more information from **dbSNP**

  As part of the database creation process, GEMINI already checks all
  variants whether they occur in dbSNP and, if so, stores their dbSNP IDs.
  Hence, we only need to extract some additional information of interest.

- information from **Cancer Hotspots**
- links to the **CIViC** database
- information from the **Cancer Genome Interpreter (CGI)**

> ### {% icon hands_on %} Hands-on: Adding further annotations
> 1. **GEMINI annotate** {% icon tool %} to add further annotations from **dbSNP**
>    - {% icon param-file %} *"GEMINI database"*: the output of the last
>      **GEMINI annotate** {% icon tool %} run
>    - {% icon param-file %} *"Dataset to use as the annotation source"*: the imported `dbsnp.b147.chr5_12_17.vcf`
>    - *"Strict variant-identity matching of database and annotation records
>      (VCF format only)"*: `Yes`
>   
>      dbSNP stores information about specific SNPs observed in human
>      populations and we would like to record if any exact same SNPs are among
>      our variants.
>
>    - *"Type of information to add to the database variants"*: `Specific
>      values extracted from matching records in the annotation source
>      (extract)`
>      - In *"1: Annotation extraction recipe"*:
>        - *"Elements to extract from the annotation source"*: `SAO`
>        - *"Database column name to use for recording annotations"*:
>        `rs_ss`
>        - *"What type of data are you trying to extract?"*: `Integer numbers`
>        - *"If multiple annotations are found for the same variant,
>        store ..."*: `the first annotation found`
>
>        This recipe extracts the dbSNP *SAO* field and adds it as *rs_ss* to
>        the GEMINI database.
>
>      - {% icon param-repeat %} *"Insert Annotation extraction recipe"*
>      - In *"2: Annotation extraction recipe"*:
>        - *"Elements to extract from the annotation source"*: `CFL`
>        - *"Database column name to use for recording annotations"*:
>       `rs_cfl`
>        - *"What type of data are you trying to extract?"*: `Numbers with
>        decimal precision`
>        - *"If multiple annotations are found for the same variant,
>        store ..."*: `the first annotation found`
>
>        This recipe extracts the dbSNP *CFL* field and adds it as *rs_cfl* to
>        the GEMINI database.
>
>      - {% icon param-repeat %} *"Insert Annotation extraction recipe"*
>      - In *"3: Annotation extraction recipe"*:
>        - *"Elements to extract from the annotation source"*: `ASP`
>        - *"Database column name to use for recording annotations"*:
>       `rs_asp`
>        - *"What type of data are you trying to extract?"*: `Integer numbers`
>        - *"If multiple annotations are found for the same variant,
>        store ..."*: `the first annotation found`
>
>        This recipe extracts the dbSNP *ASP* field and adds it as *rs_asp* to
>        the GEMINI database.
> 2. **GEMINI annotate** {% icon tool %} to add further annotations from **cancerhotspots**
>    - {% icon param-file %} *"GEMINI database"*: the output of the last
>      **GEMINI annotate** {% icon tool %} run
>    - {% icon param-file %} *"Dataset to use as the annotation source"*: the imported `cancerhotspots_v2.bed`
>    - *"Strict variant-identity matching of database and annotation records
>      (VCF format only)"*: `Yes` (with input in BED format this setting will
>      be ignored)
>
>      For the cancerhotspots data, we are simply going to record the best
>      q-value associated with any cancerhotspots variant overlapping one
>      of our variant sites.
>
>    - *"Type of information to add to the database variants"*: `Specific
>      values extracted from matching records in the annotation source
>      (extract)`
>      - In *"1: Annotation extraction recipe"*:
>        - *"Elements to extract from the annotation source"*: `5`
>
>          The q-values are stored in the fifth column of the BED dataset.
>
>        - *"Database column name to use for recording annotations"*:
>        `hs_qvalue`
>        - *"What type of data are you trying to extract?"*: `Numbers with
>         decimal precision`
>        - *"If multiple annotations are found for the same variant,
>        store ..."*: `the smallest of the (numeric) values`
>
>        This is the recipe for extracting the *q-values* of overlapping
>        cancerhotspots sites and adding them as *hs_qvalue* to the GEMINI
>        database.
>
> 3. **GEMINI annotate** {% icon tool %} to add links to **CIViC**
>    - {% icon param-file %} *"GEMINI database"*: the output of the last
>      **GEMINI annotate** {% icon tool %} run
>    - {% icon param-file %} *"Dataset to use as the annotation source"*: the imported `CIViC.bed`
>    - *"Strict variant-identity matching of database and annotation records
>      (VCF format only)"*: `Yes` (with input in BED format this setting will
>      be ignored)
>
>      For the CIViC data, we are going to record the hyperlink to any variant
>      found in the CIViC database that overlaps one of our variant sites.
>
>    - *"Type of information to add to the database variants"*: `Specific
>      values extracted from matching records in the annotation source
>      (extract)`
>      - In *"1: Annotation extraction recipe"*:
>        - *"Elements to extract from the annotation source"*: `4`
>
>          The hyperlinks are stored in the fourth column of the BED dataset.
>
>        - *"Database column name to use for recording annotations"*:
>        `overlapping_civic_url`
>        - *"What type of data are you trying to extract?"*: `Text (text)`
>        - *"If multiple annotations are found for the same variant,
>        store ..."*: `a comma-separated list of non-redundant (text) values`
>
>        This is the recipe for extracting the hyperlinks of overlapping CIViC
>        sites and adding them as a list of *overlapping_civic_urls* to the
>        GEMINI database.
>
> 3. **GEMINI annotate** {% icon tool %} to add information from the **Cancer Genome Interpreter (CGI)**
>    - {% icon param-file %} *"GEMINI database"*: the output of the last
>      **GEMINI annotate** {% icon tool %} run
>    - {% icon param-file %} *"Dataset to use as the annotation source"*: the imported `cgi_variant_positions.bed`
>    - *"Strict variant-identity matching of database and annotation records
>      (VCF format only)"*: `Yes` (with input in BED format this setting will
>      be ignored)
>
>      For the CGI data, we are going to record if any of our variant sites
>      is overlapped by a variant in the CGI biomarkers database.
>
>    - *"Type of information to add to the database variants"*: `Binary
>      indicator (1=found, 0=not found) of whether the variant had any match in
>      the annotation source (boolean)`
>      - *"Database column name to use for recording annotations"*: `in_cgidb`
>
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
somatic variants. Our strategy for retrieving them is to:

1. rely on the *somatic status* of the variants called by **VarScan somatic**
   {% icon tool %}
2. disregard questionable variants, for which either a non-negligible amount of
   supporting sequencing reads is also found in the normal tissue data, or
   which are only supported by a very small fraction of the reads from the
   tumor sample

> ### {% icon hands_on %} Hands-on: Querying the GEMINI database for somatic variants
> 1. **GEMINI query** {% icon tool %} with:
>    - {% icon param-file %} *"GEMINI database"*: the fully annotated database
>      created in the last **GEMINI annotate** {% icon tool %} step
>    - *"Build GEMINI query using"*: `Basic variant query constructor`
>      - {% icon param-repeat %} *"Insert Genotype filter expression"*
>        - *"Restrictions to apply to genotype values"*:
>          `gt_alt_freqs.NORMAL <= 0.05 AND gt_alt_freqs.TUMOR >= 0.10`
>
>        With this genotype-based filter, we retain only those variants that
>        are supported by less than 5% of the reads of the normal sample, but
>        by more than 10% of the reads of the tumor sample.
>
>      - *"Additional constraints expressed in SQL syntax"*: `somatic_status = 2`
>
>        Among the info stored in the GEMINI database is the somatic status
>        VarScan somatic has called for every variant (remember we
>        used GEMINI annotate to add it). With the condition
>        `somatic_status = 2` we retain only those variants passing the
>        genotype filter above **and** considered somatic variants by the
>        variant caller.
>
>      - In *"Output format options"*
>        - *"Type of report to generate"*: `tabular (GEMINI default)`
>          - *"Add a header of column names to the output"*: `Yes`
>          - *"Set of columns to include in the variant report table"*: `Custom
>          (report user-specified columns)`
>            - In *"Choose columns to include in the report"*:
>              - {% icon param-check %} *"chrom"*
>              - {% icon param-check %} *"start"*
>              - {% icon param-check %} *"ref"*
>              - {% icon param-check %} *"alt"*
>            - *"Additional columns (comma-separated)"*: `gene, aa_change, rs_ids,
>              hs_qvalue, cosmic_ids`
>
>        Here we specify, which columns (from the *variants* table of the
>        GEMINI database) we want to have included, in the specified order, in
>        a `tabular` variant report.
>
{: .hands_on}

> ### {% icon comment %} How am I supposed to know these column names?
> Obviously, you need to know the names of the columns (in the tables of the
> GEMINI database) to include them in the report, but how are you supposed to
> know them?
>
> The standard ones (added by *GEMINI load* when building the database) are
> listed in the
> [GEMINI documentation](https://gemini.readthedocs.io/en/latest/content/database_schema.html).
> The non-standard columns (the ones you added with *GEMINI annotate*) have
> the names you gave them, when you added them.
>
> Alternatively, to get the tables and column names of a specific database
> listed, you can use **GEMINI database info** {% icon tool %} like so:
>
> - *"GEMINI database"*: the database you want to explore
> - *"Information to retrieve from the database"*: `Names of database tables
> and their columns`
{: .comment}

What about more sophisticated filtering?

> ### {% icon hands_on %} Hands-on: More complex filter criteria
> 1. **GEMINI query** {% icon tool %} with the exact same settings as before, but:
>    - *"Additional constraints expressed in SQL syntax"*: `somatic_status = 2 AND somatic_p <= 0.05 AND (filter IS NULL OR rs_ids IS NOT NULL) AND rs_cfl != 1 AND rs_asp != 1`
>
>     This translates into "variants classified as somatic with a p-value <=
>     0.05, which haven't been flagged as likely false-positives or, if so, are
>     listed in dbSNP, but in there, are not flagged as being assembly-dependent
>     or -specific".
{: .hands_on}

> ### {% icon comment %} SQL keywords
> In the condition above, SQL keywords are given in uppercase. This is not a
> requirement, but it makes it easier to understand the syntax.
>
> You can check whether any cell in a data table is empty with `IS NULL`, and
> whether it contains *any* value with `IS NOT NULL`. To combine different
> filter criteria logically, you can use `AND` and `OR`, and parentheses to
> group conditions if required.
{: .comment}

If you have followed all steps up to here exactly, running this job should give
you a tabular dataset of 43 variants, and with the annotations in the report
it is relatively easy to pick out a few interesting ones.
Before we focus on the content of the report, however, we could enhance the
report format a bit more.

> ### {% icon hands_on %} Hands-on: SQL-based output formatting
> 1. **GEMINI query** {% icon tool %} with the exact same settings as in the
> last example, but:
>    - In *"Output format options"*
>      - *"Additional columns (comma-separated)"*: `type, gt_alt_freqs.TUMOR, gt_alt_freqs.NORMAL,
>        ifnull(nullif(round(max_aaf_all,2),-1.0),0) AS MAF, gene, impact_so,
>        aa_change, ifnull(round(cadd_scaled,2),'.') AS cadd_scaled,
>        round(gerp_bp_score,2) AS gerp_bp, ifnull(round(gerp_element_pval,2),'.')
>        AS gerp_element_pval, ifnull(round(hs_qvalue,2), '.') AS hs_qvalue,
>        in_omim, ifnull(clinvar_sig,'.') AS clinvar_sig,
>        ifnull(clinvar_disease_name,'.') AS clinvar_disease_name,
>        ifnull(rs_ids,'.') AS dbsnp_ids, rs_ss, ifnull(cosmic_ids,'.') AS
>        cosmic_ids, ifnull(overlapping_civic_url,'.') AS overlapping_civic_url,
>        in_cgidb`
{: .hands_on}

This last query adds a lot more annotations to the report, and it also
demonstrates the use of the `AS` keyword to rename columns and of some
[SQLite functions](https://sqlite.org/lang_corefunc.html) to clean up the
output.

> ### {% icon question %} Question
> Compare the new report to the previous one to see what has changed.
{: .question}

> ### {% icon details %} Limitations of genotype column queries
> In case you are wondering why the above query does not use rounding on
> the alternate allele frequencies of the samples, *i.e.*, on
> `gt_alt_freqs.TUMOR` and `gt_alt_freqs.NORMAL`, or why it does not rename
> these columns, that is because it would break the query.
>
> As a general rule, note that all columns in the variants table starting
> with `gt` (the genotype columns, calculated from the genotype fields of a
> VCF dataset) cannot be used like regular SQLite columns, but are parsed by
> GEMINI separately. That is why you cannot mix them with regular SQLite
> elements like functions and alias specifications with `AS`.
> You may also have noticed that `gt_alt_freqs.TUMOR` and
> `gt_alt_freqs.NORMAL` do not obey the column order specification of the
> query, but end up as the last columns of the tabular report. This is
> another artefact of GEMINI's special treatment of `gt` columns.
{: .details}

## Generating reports of genes affected by variants

As a final step, let us now try to generate a gene-centered report based on the
same somatic variants we just selected above.

Such a gene-centered report would include annotations that apply to a whole
gene affected by a variant rather than to the variant itself. Examples of such
annotations include known synonyms of an affected gene, its NCBI entrez number,
the ClinVar phenotype, if any, associated with the gene, a hyperlink to the
gene's page at CIViC.org, *etc.*.

Some of this information comes built-in into every GEMINI database, but it is
stored in a separate table called `gene_detailed`, while all information we
used and queried so far was from the `variants` table.

To access information from the `variants` and the `gene_detailed` table in the
same query we need to join the two tables. Such an operation is not possible
with the `Basic query constructor` we have used so far, but requires an
advanced mode for composing the query.

> ### {% icon hands_on %} Hands-on: Turning query results into gene-centered reports
> 1. **GEMINI query** {% icon tool %} in advanced mode by choosing
>   - *"Build GEMINI query using"*: `Advanced query constructor`
>   - *"The query to be issued to the database"*: `SELECT v.gene, v.chrom,
>     g.synonym, g.hgnc_id, g.entrez_id, g.rvis_pct, v.clinvar_gene_phenotype
>     FROM variants v, gene_detailed g WHERE v.chrom = g.chrom AND
>     v.gene = g.gene AND v.somatic_status = 2 AND v.somatic_p <= 0.05 AND
>     (v.filter IS NULL OR v.rs_ids IS NOT NULL) AND v.rs_cfl != 1 and
>     v.rs_asp != 1 GROUP BY g.gene`
>
>     > ### {% icon comment %} Elements of the SQL query
>     > In this full SQL query, the part between `SELECT` and `FROM` states which
>     > columns from which database tables we wish to retrieve, while the part
>     > between `FROM` and `WHERE` specifies the database tables that need to be
>     > consulted and gives them simpler aliases (`v` becomes an alias for
>     > the `variants` table, `g` for the `gene_detailed` table), which we can then
>     > use everywhere else in the query.
>     >
>     > The part following `WHERE` are the filter criteria we want to apply. Note
>     > that these criteria are almost the same as those we have used in our
>     > earlier somatic variants query, but since we are working with two tables
>     > instead of just one, we need to state which table the filter columns come
>     > from through table prefixes. Thus, `somatic_status` becomes
>     > `v.somatic_status`, *etc.*. In addition, we want to report, of course,
>     > corresponding information from the two tables, which is ensured by the
>     > additional criteria `v.chrom = g.chrom` and `v.gene = g.gene` (the SQL
>     > terminology for this is: we want to join the `variants` and the
>     > `gene_detailed` tables on their `chrom` and `gene` columns).
>     >
>     > The `GROUP BY` part, finally, specifies that we want to collapse records
>     > affecting the same gene into one.
>     {: .comment}
>
>   - *"Genotype filter expression"*: `gt_alt_freqs.NORMAL <= 0.05 AND
>     gt_alt_freqs.TUMOR >= 0.10`
>
>     This remains the same as in the previous somatic variants query.    
{: .hands_on}

## Adding additional annotations to the gene-centered report

Unfortunately, *GEMINI annotate* lets you add columns only to the variants
table of a GEMINI database, but there is no simple way to enhance the
`gene_detailed` table with additional annotations. That's why we are going to
add such extra annotations now to the tabular gene-centered report using more
general-purpose Galaxy tools.

Annotating the tabular gene report produced with GEMINI is a two-step
process, in which we first *join* the report and tabular annotation sources
into a larger tabular dataset, from which we then eliminate redundant and
unwanted columns, while *rearranging* the remaining ones.

**Step 1** consists of three separate *join* operations that sequentially pull in the annotations found in the three gene-based tabular datasets that you imported in the *Get Data* step of this section.

> ### {% icon hands_on %} Hands-on: Join 
> 
> 1. **Join two files** {% icon tool %} to add **UniProt cancer genes** information
>    - {% icon param-file %} *"1st file"*: the GEMINI-generated gene report from the previous step
>    - *"Column to use from 1st file"*: `Column: 1`
>    - {% icon param-file %}  *"2nd file"*: the imported `Uniprot_Cancer_Genes` dataset 
>    - *"Column to use from 2nd file"*: `Column: 1`
>    - *"Output lines appearing in"*: `Both 1st & 2nd file, plus unpairable
>      lines from 1st file. (-a 1)`
>
>      If a variant-affected gene is not listed as a Uniprot cancer gene, then,
>      obviously, we still want to have it included in the final report.
>
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
> 2. **Join two files** {% icon tool %} to add **CGI biomarkers** information
>    - {% icon param-file %} *"1st file"*: the partially annotated dataset from the previous
>    - *"Column to use from 1st file"*: `Column: 1`
>    - {% icon param-file %} *"2nd file"*: the imported `cgi_genes` dataset
>    - *"Column to use from 2nd file"*: `Column: 1`
>    - *"Output lines appearing in"*: `Both 1st & 2nd file, plus unpairable
>      lines from 1st file. (-a 1)`
>    - *"First line is a header line"*: `Yes`
>    - *"Ignore case"*: `No`
>    - *"Value to put in unpaired (empty) fields"*: `0`
>
>    Inspect the input and the result dataset to make sure you understand what
>    happened at this step.
>
> 3. **Join two files** {% icon tool %} to add gene information from **CIViC**
>    - {% icon param-file %} *"1st file"*: the partially annotated dataset from step 2
>    - *"Column to use from 1st file"*: `Column: 1`
>    - {% icon param-file %} *"2nd file"*: the imported `GeneSummaries` dataset
>    - *"Column to use from 2nd file"*: `Column: 3`
>
>      The gene column in the CIViC gene summaries annotation dataset is *not* the first one!
>   
>    - *"Output lines appearing in"*: `Both 1st & 2nd file, plus unpairable
>      lines from 1st file. (-a 1)`
>    - *"First line is a header line"*: `Yes`
>    - *"Ignore case"*: `No`
>    - *"Value to put in unpaired (empty) fields"*: `.`
>
>    Inspect the input and the result dataset to make sure you understand what
>    happened at this step.
>
{: .hands_on}

If you took a look at all output datasets as suggested, you will have noticed
that each of the *join* operations kept the gene columns from both input
datasets. In addition, we had no control over the order, in which columns got
added to the report, nor could we exclude columns.

In **Step 2** of the report preparation we are going to address all of these
issues and rearrange to get a fully annotated gene report.


> ### {% icon hands_on %} Hands-on: Rearrange to get a fully annotated gene report 
> 4. **Column arrange by header name** {% icon tool %} configured like
>    this:
>    - {% icon param-file %} *"file to rearrange"*: the final **Join** result dataset from step 3
>    - In *"Specify the first few columns by name"*
>      - In *"1: Specify the first few columns by name"*
>        - *"column"*: `gene`
>      - {% icon param-repeat %} *"Specify the first few columns by name"*   
>      - In *"2: Specify the first few columns by name"*   
>        - *"column"*: `chrom`
>      - {% icon param-repeat %} *"Specify the first few columns by name"*   
>      - In *"3: Specify the first few columns by name"*   
>        - *"column"*: `synonym`
>      - {% icon param-repeat %} *"Specify the first few columns by name"*   
>      - In *"4: Specify the first few columns by name"*   
>        - *"column"*: `hgnc_id`
>      - {% icon param-repeat %} *"Specify the first few columns by name"*   
>      - In *"5: Specify the first few columns by name"*   
>        - *"column"*: `entrez_id`
>      - {% icon param-repeat %} *"Specify the first few columns by name"*   
>      - In *"6: Specify the first few columns by name"*   
>        - *"column"*: `rvis_pct`
>      - {% icon param-repeat %} *"Specify the first few columns by name"*   
>      - In *"7: Specify the first few columns by name"*   
>        - *"column"*: `is_OG`
>      - {% icon param-repeat %} *"Specify the first few columns by name"*   
>      - In *"8: Specify the first few columns by name"*   
>        - *"column"*: `is_TS`
>      - {% icon param-repeat %} *"Specify the first few columns by name"*   
>      - In *"9: Specify the first few columns by name"*   
>        - *"column"*: `in_cgi_biomarkers`
>      - {% icon param-repeat %} *"Specify the first few columns by name"*   
>      - In *"10: Specify the first few columns by name"*   
>        - *"column"*: `clinvar_gene_phenotype`
>      - {% icon param-repeat %} *"Specify the first few columns by name"*   
>      - In *"11: Specify the first few columns by name"*   
>        - *"column"*: `gene_civic_url`
>      - {% icon param-repeat %} *"Specify the first few columns by name"*   
>      - In *"12: Specify the first few columns by name"*   
>        - *"column"*: `description`
>    - *"Discard unspecified columns"*: `Yes`
>
> > ### {% icon comment %} Alternative tool suggestion
> > If your Galaxy server does not offer the *Column arrange* tool, it will
> > almost certainly offer **Cut columns from a table** {% icon tool %},
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
