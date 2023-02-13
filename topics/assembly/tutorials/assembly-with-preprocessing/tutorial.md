---
layout: tutorial_hands_on

title: "Unicycler assembly of SARS-CoV-2 genome with preprocessing to remove human genome reads"
zenodo_link: "https://doi.org/10.5281/zenodo.3732358"
questions:
  - How can a genome of interest be assembled against a background of contaminating reads from other genomes?
  - How can sequencing data from public sources be turned into assembly-ready polished datasets?
objectives:
  - Obtain viral (SARS-CoV-2) sequencing data with contaminating human reads from public sources
  - Organize the data into collections and check its quality
  - Detect and remove human reads
  - Assemble retained reads and explore the results
time_estimation: "4h"  # plus additional time for (optional) NCBI SRA downloads
level: Intermediate
key_points:
  - Certain types of NGS samples can be heavily contaminated with sequences from other genomes
  - Reads from known/expected contaminating sources can be identified by mapping to the respective genomes
  - After mapping, use filtering tools to remove identified contaminating reads, and use conversion tools to convert remaining mapped reads back into raw sequenced reads expected by most downstream tools
requirements:
  -
    type: "internal"
    topic_name: assembly
    tutorials:
      - unicycler-assembly
  -
    type: "internal"
    topic_name: sequence-analysis
    tutorials:
      - mapping
  -
    type: "internal"
    topic_name: galaxy-interface
    tutorials:
      - collections
tags:
  - covid19
contributors:
  - wm75

---

# Introduction


In some research or clinical contexts it is not possible, or very hard, to
purify DNA/RNA for sequencing from just the specimen of interest.
Instead you will isolate DNA that is contaminated, sometimes heavily, with
DNA/RNA of a different origin.
This is the case for example with microbiome samples, which typically display
considerable contamination with host DNA, or with samples of body fluids for
pathogen detection. Such contamination can pose an issue with certain types of
analyses, in particular with genome assembly.

This tutorial guides you through the preprocessing of sequencing data of
bronchoalveolar lavage fluid (BALF) samples obtained from early COVID-19
patients in China. Since such samples are expected to be contaminated
signficantly with human sequenced reads, the goal is to enrich the data for
SARS-CoV-2 reads by identifying and discarding reads of human origin before
trying to assemble the viral genome sequence.

> <comment-title>The usegalaxy.* COVID-19 analysis project</comment-title>
> This tutorial uses the same data as, and recapitulates to a large extent, the
> [Pre-processing](https://covid19.galaxyproject.org/genomics/1-PreProcessing/)
> and [Assembly](https://covid19.galaxyproject.org/genomics/2-Assembly/) steps
> of the [Genomics](https://covid19.galaxyproject.org/genomics/) section of
> [covid19.galaxyproject.org](https://covid19.galaxyproject.org/).
>
{: .comment}

> <agenda-title></agenda-title>
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get data

We are going to create an assembly of the SARS-CoV-2 genome based on mixed
short-reads (Illumina) and long-reads (Nanopore) data from a total of six
different samples, all of which are publicly accessible through the NCBI
and EBI short-reads archives (SRAs).

Since automated data downloads from SRAs can be unreliable at times and could
get broken by changes to the download interface on the NCBI/EBI side, this
tutorial offers two ways to access the sequenced reads input data:

1. Direct download from the NCBI SRA based on accession numbers and using the
   dedicated **Faster Download** tool

   Use this method if it works and is fast enough for you, and if you are
   interested in learning to obtain short-reads data directly from NCBI, in
   general.

2. Download of the same data deposited as a copy at [Zenodo](https://zenodo.org/record/3732359)
   This method uses Galaxy's generic data import functionality, and should be
   very reliable and faster than the download from NCBI.
   It also showcases **rule-based** uploads and demonstrates how they can be
   used to download several datasets and to arrange them into easy to handle
   data structures at the same time.

   > <details-title>Rule-based uploads</details-title>
   > In this tutorial you will only use the features of Galaxy's rule-based
   > uploader that are required to get the input data ready for our analysis,
   > and we will not explain those features in much detail.
   >
   > If, after this first taste, you are interested in a thorough introduction
   > we recommend the advanced tutorial
   > [Collections: Rule Based Uploader](../../../galaxy-interface/tutorials/upload-rules/tutorial.html).
   {: .details}

   Use this method if the direct download from the NCBI SRA does not work, or
   is too slow for your time frame, or if you are interested in advanced use
   of Galaxy's data import functionality.

The corresponding two step-by-step instructions below have been crafted to
produce identically arranged data structures in your history so all subsequent
steps are independent of the data source you choose.

## Get data from NCBI SRA

> <hands-on-title>Data upload to Galaxy from NCBI SRA</hands-on-title>
>
> 1. Create a new history for this tutorial and give it a proper name
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>    {% snippet faqs/galaxy/histories_rename.md %}
>
> 2. Create a new dataset listing the SRA accession numbers of the Illumina paired-end input data for this tutorial:
>
>    ```
>    SRR10903401
>    SRR10903402
>    SRR10971381
>    ```
>
>    call it, *e.g.*, `Illumina accessions` and set its datatype to `tabular`.
>
>    {% snippet faqs/galaxy/datasets_create_new_file.md format="tabular" %}
>
> 3. Create another new dataset listing the SRA accession numbers of the Nanopore input data for this tutorial:
>
>    ```
>    SRR10948550
>    SRR10948474
>    SRR10902284
>    ```
>
>    call it, *e.g.*, `Nanopore accessions` and set its datatype to `tabular`.
>
> 4. Add `#illumina`/`#nanopore` tags to the datasets
>
>    > <comment-title>Name tags in the analysis</comment-title>
>    > We are going to treat the Illumina- and the Nanopore-sequenced data
>    > separately in this tutorial up to the actual genome assembly step.
>    >
>    > To make it easier to keep track of the two branches of the analysis, we
>    > recommend the use of Galaxy's dataset **name tags**.
>    > A name tag will automatically propagate to any new dataset derived
>    > from the tagged dataset.
>    {: .comment}
>    
>    You can create a name tag by attaching a tag starting with `#` to any
>    dataset.
>
>    Name tags are meant to help you identify the origin of datasets quickly.
>    Feel free to either use the suggested names above or choose ones you like.
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
> 4. Retrieve the Illumina reads data from NCBI:
>
>    Run **Faster Download and Extract Reads in FASTQ** {% icon tool %} with
>    the following parameter settings:
>    - *"select input type"*: `List of SRA accession, one per line`
>      - {% icon param-file %} *"sra accession list"*:
>        the `Illumina accessions` dataset created above
>    - in *"Advanced Options"*
>      - *"Select how to split the spots"*: `--split-3`
>
>    The tool run should generate four new items in your history - three
>    collections and one *log* dataset with a summary of what happened.
>
>    Since all three datasets that we tried to retrieve contain only
>    paired-end reads, only the `Pair-end data` collection is expected to
>    contain downloaded data. Click on the other two collections to verify that
>    they are empty, then delete them from your history (since the collections
>    do not contain any datasets, it does not matter, which delete option -
>    "Collection Only", "Delete Datasets" or "Permanently Delete Datasets" you
>    are choosing when prompted).
>
> 5. Retrieve the Nanopore reads data from NCBI:
>
>    Run **Faster Download and Extract Reads in FASTQ** {% icon tool %} with
>    the following parameter settings:
>    - *"select input type"*: `List of SRA accession, one per line`
>      - {% icon param-file %} *"sra accession list"*:
>        the `Nanopore accessions` dataset created above
>    - in *"Advanced Options"*
>      - *"Select how to split the spots"*: `--split-3`
>
>    As in the previous step, the tool run should generate four new items in
>    your history.
>
>    Since all three datasets that we tried to retrieve in this run contain
>    only single-end reads, only the `Single-end data` collection is expected to
>    contain downloaded data this time. Click on the other two collections to
>    verify that they are empty, then delete them from your history.
>
{: .hands_on}


## Get data from Zenodo

> <hands-on-title>Data upload to Galaxy from Zenodo</hands-on-title>
>
> 1. Create a new history for this tutorial and give it a proper name
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>    {% snippet faqs/galaxy/histories_rename.md %}
>
> 2. Import Illumina-sequenced reads data from [Zenodo](https://zenodo.org/record/3732359)
>
>    The Zenodo links for the data are these:
>    ```
>    https://zenodo.org/record/3732359/files/SRR10903401_r1.fq.gz
>    https://zenodo.org/record/3732359/files/SRR10903401_r2.fq.gz
>    https://zenodo.org/record/3732359/files/SRR10903402_r1.fq.gz
>    https://zenodo.org/record/3732359/files/SRR10903402_r2.fq.gz
>    https://zenodo.org/record/3732359/files/SRR10971381_r1.fq.gz
>    https://zenodo.org/record/3732359/files/SRR10971381_r2.fq.gz
>    ```
>
>    To upload these data to your Galaxy history in structured form:
>
>      - Copy the above list of links
>      - Open the Galaxy Upload Manager
>        ({% icon galaxy-upload %} on the top right of the tool panel)
>      - In the pop-up window, switch to the **Rule-based** tab and select
>        - *"Upload data as:"*: `Collection(s)`
>        - *"Load tabular data from:"*: `Pasted Table`
>        - Paste the copied links into the text field
>        - Click **Build**
>      - In the next screen, select
>        - {% icon param-repeat %} *"Column"*: `Using a Regular expression`
>          - *"From Column"*: `A`
>          - {% icon param-check %} *"Create columns matching expression groups."*
>          - *"Regular Expression"*: `.+/(SRR\d+)_r(\d).fq.gz`
>          - *"Number of Groups"*: `2`
>          - Click **Apply**
>        - {% icon param-repeat %} *"Rules"*: `Add / Modify Column Definitions`
>          - {% icon param-repeat %} *"Add Definition"*: `URL`
>            - *"URL"*: `A`
>          - {% icon param-repeat %} *"Add Definition"*: `List Identifier(s)`
>            - *"List Identifier(s)"*: `B`
>          - {% icon param-repeat %} *"Add Definition"*: `Paired-end Indicator`
>            - *"Paired-end Indicator"*: `C`
>          - Click **Apply**
>        - *"Type"*: `fastqsanger.gz`
>        - *"Name"*: `Illumina PE data` (or similar)
>        - *"Add nametag for name:"* {% icon param-check %}
>
>          > <comment-title>Name tags in the analysis</comment-title>
>          > We are going to treat the Illumina- and the Nanopore-sequenced data
>          > separately in this tutorial up to the actual genome assembly step.
>          >
>          > To make it easier to keep track of the two branches of the analysis, we
>          > recommend the use of Galaxy's dataset **name tags**.
>          > A name tag will automatically propagate to any new dataset derived
>          > from the tagged dataset.
>          {: .comment}
>
>          Checking this option tells Galaxy to reuse the collection name above
>          as a name tag on the collection.
>
>        - Click **Upload**
>
> 3. Import Nanopore-sequenced reads data from [Zenodo](https://zenodo.org/record/3732359)
>
>    ```
>    https://zenodo.org/record/3732359/files/SRR10902284_ONT.fq.gz
>    https://zenodo.org/record/3732359/files/SRR10948474_ONT.fq.gz
>    https://zenodo.org/record/3732359/files/SRR10948550_ONT.fq.gz
>    ```
>
>    Again, we want to upload this data to our history in structured form.
>    To do so:
>
>      - Copy the above list of links
>      - Open the Galaxy Upload Manager
>        ({% icon galaxy-upload %} on the top right of the tool panel)
>      - In the pop-up window, switch to the **Rule-based** tab and select
>        - *"Upload data as:"*: `Collection(s)`
>        - *"Load tabular data from:"*: `Pasted Table`
>        - Paste the copied links into the text field
>        - Click **Build**
>      - In the next screen, select
>        - {% icon param-repeat %} *"Column"*: `Using a Regular expression`
>          - *"From Column"*: `A`
>          - {% icon param-check %} *"Create columns matching expression groups."*
>          - *"Regular Expression"*: `.+/(SRR\d+)_ONT.fq.gz`
>          - *"Number of Groups"*: `1`
>          - Click **Apply**
>        - {% icon param-repeat %} *"Rules"*: `Add / Modify Column Definitions`
>          - {% icon param-repeat %} *"Add Definition"*: `URL`
>            - *"URL"*: `A`
>          - {% icon param-repeat %} *"Add Definition"*: `List Identifier(s)`
>            - *"List Identifier(s)"*: `B`
>          - Click **Apply**
>        - *"Type"*: `fastqsanger.gz`
>        - *"Name"*: `Nanopore data` (or similar)
>        - *"Add nametag for name:"* {% icon param-check %}    
>        - Click **Upload**
>
{: .hands_on}


# Read trimming and quality control

In the end, we do not want to base our mapping on low-quality reads that may
cause misassembly of fragments or introduce sequencing artefacts into the final
assembled sequence. Since quality control, read filtering and read trimming are
all quite fast and computationally cheap operations compared to the read
mapping that we will use to identify and eliminate human reads, it is best to
perform these steps up front.

Due to their very different nature, however, Illumina- and Nanopore-sequenced
reads should be treated rather differently.

## Trimming and filtering of Illumina reads

Galaxy offers a panel of different NGS reads trimming/filtering tools. Here,
we use **fastp** {% icon tool %}, which is straightforward to configure,
and when combined with **MultiQC** {% icon tool %}, enables nice and
easy-to-interpret visualizations of the effects of preprocessing, in particular
for multiple samples.

In the following, we configure the tool to retain reads only if at most 20% of
their bases have a Phred-scaled quality >= 20 and if there length in bases after
trimming of adapter sequences (which the tool auto-detects for us) is at least
50.

The JSON-formatted report produced by the tool, can serve as input to
**MultiQC** {% icon tool %} for a direct visual comparison of key quality
metrics for all samples before and after preprocessing.

> <hands-on-title>Reads preprocessing and quality reporting</hands-on-title>
> 1. **fastp** {% icon tool %} with the following parameters
>    - *"Single-end or paired reads"*: `Paired Collection`
>      - *"Select paired collection(s)"*: the collection of Illumina-sequenced
>        reads as produced in the *Get Data* section
>    - in *"Filter Options"*
>      - in *"Quality filtering options"*
>        - *"Disable quality filtering"*: `No`
>        - *"Qualified quality phred"*: `20`
>        - *"Unqualified percent limit"*: `20`
>      - in *"Length filtering options"*
>        - *"Disable length filtering"*: `No`
>        - *"Length required"*: `50`
>    - in *"Output Options"*
>      - *"Output HTML report"*: `No`
>      - *"Output JSON report"*: `Yes`
>
>    The tool run produces two collections - one with the actual preprocessed
>    reads of all input samples, another one with a JSON-formatted report of
>    the processing for every sample.
>
> 2. **MultiQC** {% icon tool %} with the following parameters
>    - {% icon param-repeat %} *"Results"*
>      - *"Which tool was used generate logs?"*: `fastp`
>        - {% icon param-collection %} *"Output of fastp"*: the collection of
>          JSON-formatted reports, second collection produced by **fastp**
>          {% icon tool %}
>
>    This tool run generates a single webpage output with the combined quality reports
>    for all samples before and after processing with **fastp** {% icon tool %}.
>
{: .hands_on}

## Quality control of Nanopore reads

Nanopore-sequenced reads differ greatly in length from one another and are, on
average, of relatively low quality (in particular when compared to
Illumina-sequenced reads). These properties make it challenging to preprocess
them with standard tools. Quality assessment tools calibrated to work well with
Illumina-sequenced reads are less useful for Nanopore-sequenced reads, too, for
the same reasons. Here we restrict ourselves to a simple quality check with
**NanoPlot** {% icon tool %}, a dedicated QC tool for Nanopore-sequenced reads.

> <hands-on-title>Checking the quality of Nanopore reads with NanoPlot</hands-on-title>
> 1. **NanoPlot** {% icon tool %} with the following parameters
>    - *"Select multifile mode"*: `batch`
>      - *"Type of the file(s) to work on"*: `fastq`
>        - {% icon param-collection %} *"files"*: the collection of
>          Nanopore-sequenced reads as produced in the *Get Data* section
>    - in *"Options for filtering or transforming input prior to plotting"*
>      - *"Logarithmic scaling of lengths in plots."*: `Yes`
>
>    > <question-title></question-title>
>    >
>    > 1. Looking at the three generated quality reports, which of the three
>    >    samples seems to be of better quality overall than the other two,
>    >    and what are some criteria that support this conclusion?
>    >
>    > > <solution-title></solution-title>
>    > > 1. Sample `SRR10948474` has the best overall quality.
>    > >
>    > >    It has both higher average read length and quality than the other
>    > >    two (see the *Summary statistics* table), has a distribution of
>    > >    read qualities that peaks around intermediate quality scores (not
>    > >    low ones as for the other samples), and contains some extra-long
>    > >    (though at least partly rather low-quality) reads. Those last two
>    > >    points become most obvious when looking at the *Read lengths vs
>    > >    Average read quality* plots.
>    >  {: .solution }
>    {: .question}
>
{: .hands_on}

# Subtraction of reads mapping to the human reference genome

## Mapping of Illumina reads

In this tutorial, we are using **Bowtie2** for mapping our short-reads data to
the human genome. *BWA-MEM* would be a good alternative for mapping the 150
nucleotides (see the QC report above) reads from our samples.

According to its authors, the *Minimap2* aligner, which we will be using for
mapping the Nanopore-sequenced data in the next step, is supposed to outcompete
*Bowtie2* and *BWA-MEM* in terms of speed even for Illumina-sequenced reads
of length > 100 nts, but we opt for the conservative approach of using a
widely-used, well-tested tool here.

> <hands-on-title>Mapping with Bowtie2</hands-on-title>
> 1. **Bowtie2** {% icon tool %} with the following parameters
>    - *"Is this single or paired library"*: `Paired-end Dataset Collection`
>       - *"FASTQ Paired Dataset"*: the collection of preprocessed
>         Illumina-sequenced reads, output of **fastp** {% icon tool %}
>       - *"Write unaligned reads (in fastq format) to separate file(s)"*: `No`
>
>         Activating this option may seem attractive since the unaligned reads
>         are what we are interested in, but filtering for those reads with
>         a dedicated tool in a separate step allows us to filter on the
>         properties of the read pairs instead of those of individual reads.
>
>       - *"Write aligned reads (in fastq format) to separate file(s)"*: `No`
>       - *"Do you want to set paired-end options?"*: `No`
>
>     - *"Will you select a reference genome from your history or use a built-in index?"*:
>       `Use a built-in genome index`
>       - *"Select reference genome"*: `Human (Homo sapiens): hg38 Full`
>     - *"Set read groups information?"*: `Do not set`
>     - *"Select analysis mode"*: `Default setting only`
>     - *"Do you want to tweak SAM/BAM Options?"*: `No`
>     - *"Save the bowtie2 mapping statistics to the history"*: `Yes`
>
>    The tool run should produce two collections of output datasets. One with
>    the actual mapped reads of the three samples and one with the corresponding
>    mapping statistics for each sample, which we want to have a brief look at
>    next.
>
> 2. Inspect the `mapping stats` of each sample by clicking on the corresponding collection, then on the {% icon galaxy-eye %} (eye) icon of each individual sample data
>
>    > <question-title></question-title>
>    >
>    > 1. What percentage of reads of each sample has been aligned to the `hg38` reference genome?
>    > 2. Which sample is the least contaminated with human reads?
>    > 3. Which sample contains the highest amount of SARS-CoV2 reads?
>    >
>    > > <solution-title></solution-title>
>    > > 1. The samples have between 13% and 21% of reads aligned to `hg38`.
>    > >    The information can be found on the last line of output for each
>    > >    sample.
>    > > 2. Sample `SRR10971381` is the least contaminated with just above 13%
>    > >    of human reads
>    > > 3. You cannot answer this question with this data. While `SRR10971381`
>    > >    shows the lowest relative contamination with human reads, that
>    > >    does not mean that all other reads are from SARS-CoV-2. They could
>    > >    come from other species (*e.g.*, bacteria or other viruses)
>    > >    contained in this BALF sample.
>    >  {: .solution }
>    {: .question}
>
{: .hands_on}

## Mapping of Nanopore reads

For the mapping of the Nanopore-sequenced data we are using the **Minimap2**
aligner, which is particularly efficient for mapping long reads.

> <hands-on-title>Nanopore reads mapping</hands-on-title>
>
> 1. **Map with minimap2** {% icon tool %} with the following parameters
>    - *"Will you select a reference genome from your history or use a built-in index?"*:
>      `Use a built-in genome index`
>      - *"Using reference genome"*: `Human Dec. 2013 (GRCh38/hg38) (hg38)`
>
>    - *"Single or Paired-end reads"*: `Single`
>      - {% icon param-collection %} *"Select fastq dataset"*: the collection
>        of Nanopore-sequenced reads as set up in the *Get Data* section
>
>    - *"Select a profile of preset options"*: `Oxford Nanopore read to reference mapping. ...`
>
>    This tool run produces one collection with the actual mapped reads for
>    each Nanopore-sequenced sample. Unlike **Bowtie2** it does not have an
>    option to output mapping statistics directly. However, we can generate
>    that information through an extra step.
>
> 2. **Samtools stats** {% icon tool %} with the following parameters
>    - {% icon param-collection %} *"BAM file"*: the collection of mapped
>      Nanopore-sequenced reads, output of **Map with minimap2 {% icon tool %} (step 1)
>    - *"Output"*: `Separate datasets for each statistic`
>      - *"Desired output files"*
>        - {% icon param-check %} *"Summary numbers"*
>
>          These simple summary stats correspond approximately to the
>          statistics generated by **Bowtie2** and are enough for our purpose.
>
>    > <comment-title>Mapping stats for Nanopore-sequenced long reads</comment-title>
>    >
>    > Since, unlike Illumina-generated reads, Nanopore-sequenced reads can
>    > have very different lengths, it makes limited sense to calculate a
>    > ratio of mapped to overall reads for them.
>    >
>    > Instead, `bases mapped` / `total length` should give a more reliable
>    > estimate of which fraction of the data is due to human genome sequence.
>    >
>    > Try to calculate this ratio for the three samples on your own!
>    {: .comment}
>
{: .hands_on}

## Human reads subtraction

At this point you should have two collections of mapped reads - one holding the
mapping results obtained with **Bowtie2** of the three Illumina-sequenced
samples, the other one holding the **minimap2** output for the three
Nanopore-sequenced samples.

Next, we are going to filter the data from both collections to retain only
those reads that were *not* mapped to the human genome, *i.e* those of
potential viral origin.

> <hands-on-title>Mapped reads filtering</hands-on-title>
>
> 1. **Samtools view** {% icon tool %} to filter the Illumina-sequenced reads mapped with Bowtie2:
>    - {% icon param-collection %} *"SAM/BAM/CRAM data set"*: the collection of
>      mapped Illumina-sequenced reads, output of **Bowtie2** {% icon tool %}
>    - *"What would you like to look at?"*: `A filtered/subsampled selection of reads`
>      - in *"Configure filters"*
>        - *"Require that these flags are set"*: `Read is unmapped` and
>          `Mate is unmapped`
>
>      - *"What would you like to have reported?"*: `All reads retained after filtering and subsampling`
>        - *"Produce extra dataset with dropped reads?"*: `No`
>        - *"Output format"*: `BAM (-b)`
>
> 2. **Samtools view** {% icon tool %} to filter the Nanopore-sequenced reads mapped with minimap2:
>    - {% icon param-collection %} *"SAM/BAM/CRAM data set"*: the collection of
>      mapped Nanopore-sequenced reads, output of **minimap2** {% icon tool %}
>    - *"What would you like to look at?"*: `A filtered/subsampled selection of reads`
>      - In *"Configure filters"*
>        - *"Require that these flags are set"*: `Read is unmapped`
>
>      - *"What would you like to have reported?"*: `All reads retained after filtering and subsampling`
>        - *"Produce extra dataset with dropped reads?"*: `No`
>        - *"Output format"*: `BAM (-b)`
>
> 3. (Optional) Remove the database `hg38` attribute from the output files
>
>    > <details-title>Why do this?</details-title>
>    > When we ran the **Bowtie2** {% icon tool %} and **minimap2**
>    > {% icon tool %} mappers before, these tools set the *database* attribute
>    > on their outputs to `hg38` to indicate that the mapped reads in these
>    > outputs have been mapped against this version of the human reference
>    > genome - an important piece of information for further work with those
>    > mapped reads.
>    >
>    > In step 1 and 2 above, however, we have eliminated all mapped reads so
>    > the `database: hg38` info is misleading from this point onward in the
>    > analysis. While not directly harmful, it is best practice to remove this
>    > metadata now.
>    {: .details}
>
>    For the outputs of step 1 and step 2 above, reset the database/build
>    (dbkey) to `unspecified (?)`.
>
>    {% snippet faqs/galaxy/datasets_change_dbkey.md dbkey="unspecified (?)" %}
>
{: .hands_on}


# Format conversion of remaining reads

## Conversion to fastq format

Assembly tools, typically, expect their input data to be fastq-formatted, but
what we have after mapping and filtering is data in BAM format. Hence, we need
to convert the retained Illumina- and Nanopore-sequenced reads back into their
original format before proceeding to assembly.

> <hands-on-title>BAM to fastq format conversion</hands-on-title>
>
> 1. **Samtools fastx** {% icon tool %} to convert the filtered Illumina-sequenced reads to fastq format
>    - {% icon param-collection %} *"BAM or SAM file to convert"*: the collection of
>      filtered Illumina-sequenced reads, output of first **Samtools view** {% icon tool %} run
>    - *"Output format"*: `compressed FASTQ`
>    - *"outputs"*: `READ1` and `READ2`
>
> 2. **Samtools fastx** {% icon tool %} to convert the filtered Nanopore-sequenced reads to fastq format
>    - {% icon param-collection %} *"BAM or SAM file to convert"*: the collection of
>      filtered Nanopore-sequenced reads, output of second **Samtools view** {% icon tool %} run
>    - *"Output format"*: `compressed FASTQ`
>    - *"outputs"*: `unspecific`
>
{: .hands_on}

## Optional: Rearrange the filtered data into the original nested data structure

If you compare the outputs of the last step to the input data we started out
with, you will notice that the Illumina-sequenced data is arranged differently
now than initially. It is arranged now into separate collections of forward and
reverse reads, whereas we started with a single nested collection of the data.

We can easily cast the data back into its original structure with one of
Galaxy's collection manipulation tools, but note that we will not use the
resulting nested collection for this tutorial because the **Unicycler** tool
for assembling the reads would not be able to handle the nested data correctly.

Thus, the following just serves as an illustration and is entirely optional.

> <hands-on-title>Arrange two list collections into a list of pairs</hands-on-title>
>
> 1. **Zip Collection** {% icon tool %} with the following parameters
>    - {% icon param-collection %} *"Input Dataset (Forward)"*: the collection
>      of filtered Illumina-sequenced forward reads in fastq format,
>      first output of first **Samtools fastx** {% icon tool %} run
>    - {% icon param-collection %} *"Input Dataset (Reverse)"*: the collection
>      of filtered Illumina-sequenced reverse reads in fastq format,
>      second output of first **Samtools fastx** {% icon tool %} run
>
{: .hands_on}

## Merging of reads with collection operations

To merge reads from several samples into a combined final assembly, we need to
pass the data to **Unicycler** {% icon tool %} in partially merged form. The
forward and reverse reads of paired-end data should be kept separate, and so
should short and long reads. However, the tool has no option to combine data
from individual samples, so we need to merge the forward, reverse, and the long
reads data, respectively, across samples. Conveniently for us, the outputs of
the earlier **Samtools fastx** {% icon tool %} runs have already returned the
data structured into three corresponding collections for us.

> <hands-on-title>Collapsing each collection into a single dataset</hands-on-title>
>
> 1. **Collapse Collection** {% icon tool %} of Illumina-sequenced *forward* reads
>    - {% icon param-collection %} *"Collection of files to collapse into single dataset"*:
>      the collection of filtered Illumina-sequenced forward reads in fastq format,
>      first output of first **Samtools fastx** {% icon tool %} run
>    - *"Keep one header line"*: `No`
>    - *"Prepend File name"*: `No`
>
> 2. **Collapse Collection** {% icon tool %} of Illumina-sequenced *reverse* reads
>    - {% icon param-collection %} *"Collection of files to collapse into single dataset"*:
>      the collection of filtered Illumina-sequenced reverse reads in fastq format,
>      second output of first **Samtools fastx** {% icon tool %} run
>    - *"Keep one header line"*: `No`
>    - *"Prepend File name"*: `No`
>
> 3. **Collapse Collection** {% icon tool %} of Nanopore-sequenced reads
>    - {% icon param-collection %} *"Collection of files to collapse into single dataset"*:
>      the collection of filtered Nanopore-sequenced reads in fastq format,
>      output of second **Samtools fastx** {% icon tool %} run
>    - *"Keep one header line"*: `No`
>    - *"Prepend File name"*: `No`
>
{: .hands_on}


# SARS-CoV-2 genome assembly

## Optional: Subsampling of reads

The actual assembly of the sequenced reads represents the real bottleneck in
this tutorial. Assembly of the full set of sequences can easily take 10 hours
and would best be conducted overnight.

If you do not have that much time, you should downsample the Illumina-sequenced
combined reads now. Which will reduce the time required to finish the
subsequent assembly step to approximately 1-2 hours.

> <comment-title>If you are in a hurry</comment-title>
> The downsampling parameters below have been chosen to have minimal impact on
> the assembly results. Further speed-ups are certainly possible, but will
> likely lead to poor assembly outcomes.
{: .comment}

> <hands-on-title>Subsampling of paired-end short-reads data</hands-on-title>
>
> 1. **seqtk_sample** {% icon tool %} with the following parameters
>    - {% icon param-files %} *"Input FASTA/Q file"*: The two datasets with the
>      combined Illumina-sequenced forward and reverse reads, outputs of the
>      first and the second run of **Collapse Collection**
>    - *"RNG seed"*: 4
>    - *"Subsample (decimal fraction or number)"*: 0.1
>
{: .hands_on}

## Create assembly

> <hands-on-title>Assembly of SARS-CoV2 genome</hands-on-title>
>
> 1. **Create assemblies with Unicycler** {% icon tool %} with the following parameters
>    - *"Paired or Single end data?"*: `Paired`
>      - {% icon param-file %} *"Select first set of reads"*: the combined
>        Illumina-sequenced forward reads from all samples, output of the
>        first **Collapse Collection** {% icon tool %} run (or first output of
>        **seqtk_sample** {% icon tool %} for a subsample-based assembly)
>      - {% icon param-file %} *"Select second set of reads"*: the combined
>        Illumina-sequenced reverse reads from all samples, output of the
>        second **Collapse Collection** {% icon tool %} run  (or second output
>        of **seqtk_sample** {% icon tool %} for a subsample-based assembly)
>
>    - {% icon param-file %} *"Select long reads. If there are no long reads, leave this empty"*:
>      Nanpore-sequenced reads from all samples, output of the third
>      **Collapse Collection** {% icon tool %} run
>    - *"Select Bridging mode"*: `Normal (moderate contig size and misassembly rate)`
>    - *"Exclude contigs from the FASTA file which are shorter than this length (bp)"*:
>      `100`
>    - *"The expected number of linear (i.e. non-circular) sequences in the assembly"*:
>      `1`
>
{: .hands_on}

## Explore assembly

The **Unicycler** {% icon tool %} run above should produce two output datasets:

- a final assembly in FASTA format
- an assembly graph

Of these, the assembly graph is more information-rich because it not only
contains the sequences of *all* assembled fragments (including the ones shorter
than the threshold length defined for inclusion of the fragments into the FASTA
output), but also indicates the relative average coverage of the fragments by
sequenced reads and how some of the fragments could potentially be bridged
after resolving ambiguities manually.

### Assembly inspection with Bandage

On the downside, the assembly graph format takes some getting used to before
you can make sense out of the information it provides.

This issue can be alleviated through the use of **Bandage**, a package for
exploring assembly graphs through summary reports and visualizations of their
contents.

> <hands-on-title>Assembly stats and visualization with Bandage</hands-on-title>
>
> 1. **Bandage Info** {% icon tool %} with the following parameters
>    - {% icon param-file %} *"Graphical Fragment Assembly"*: the assembly graph dataset produced by
>      **Unicycler**
>    - *"Output the information in a single tab-delimited line starting with the graph file"*:
>      `No`
>
> 2. **Bandage Image** {% icon tool %} with the following parameters
>    - {% icon param-file %} *"Graphical Fragment Assembly"*: the assembly graph dataset produced by
>      **Unicycler**
>    - *"Node name labels?"*: `Yes`
>    - *"Node length labels?"*: `Yes`
>
{: .hands_on}

Let us inspect the summary report produced by **Bandage Info** {% icon tool %}
first:

You may be rather disappointed by the large *percentage of of dead ends* in the
assembly graph (in general, lower is better here), and by the correspondingly
large *node count*. After all, should the viral sequence not be encoded on a
single small contig (a quick check at
[Wikipedia](https://en.wikipedia.org/wiki/Coronavirus#Genome) reveals that
coronaviruses have genomes in the size range of 30kb)?

On the other hand, there is the *Longest node* of 29768 bp of assembled
sequence, which is suspiciously close to the expected genome size, but a much
larger *Estimated sequence length*.

Next, take a look at the assembly graph visualization generated by **Bandage
Image** {% icon tool %} to see if that tells us more:

Indeed, this output shows that **Unicycler** {% icon tool %} managed to
assemble a good number of contigs of moderate size, then had trouble with a
number of really small fragments that it could only assemble with lots of
ambiguities (leading to that ugly clutter of nodes in the top row of the
image). Those small fragments will probably be hard to make sense of, but the
manageable list of moderate-size contigs (nodes 1-23, 25, 26) is encouraging.

Of these, node 1 is the longest node mentioned in the report with a size close
to our expectations.

Node 2 looks peculiar since **Unicycler** claims it is circular, while
Coronavirus genomes are known to be linear.

### Check origin of assembled sequences with BLAST

While we could view the actual contents of the assembly graph output of
**Unicycler** {% icon tool %} and extract node sequences of interest from it
(the longest node and that circular one could be a start), things are much
easier if we work with the FASTA output of **Unicycler** instead.

From the visualization with **Bandage Image** {% icon tool %} we know that the
separately assembled nodes are all longer than 1000 bp. We can extract those
sequences based on the length threshold in Galaxy, then BLAST all retained
sequences in one go.

> <hands-on-title>Filter FASTA sequences by their length</hands-on-title>
>
> 1. **Filter sequences by length** {% icon tool %} with the following parameters
>    - {% icon param-file %} *"Fasta file"*: the FASTA output produced by
>      **Unicycler**
>    - *"Minimal length"*: `1000`
>
>    This outputs a new FASTA datasets with only the sequences satisfying our
>    length threshold.
>
> > <tip-title>Apply length filters after instead of during assembly</tip-title>
> > You may have noted that in the **Unicycler** {% icon tool %} run we kept
> > the tool's *"Exclude contigs from the FASTA file which are shorter than
> > this length (bp)"* option at its default value of `100` instead of using
> > the 1,000 bp threshold there directly to save a step in the analysis.
> >
> > The reason we did this is that normally you will not know the exact length
> > threshold you want until *after* having explored the generated assembly.
> >
> > Length-filtering some FASTA sequences is a trivial process that takes very
> > little time, but you would not want to rerun an hours-long assembly job
> > just because you accidentally stripped some interesting assembled sequences
> > from the output.
> {: .tip}
>
{: .hands_on}

> <hands-on-title>NCBI BLAST of multiple contigs</hands-on-title>
>
> 1. View the output of **Filter sequences by length** {% icon tool %} by
>    clicking the {% icon galaxy-eye %} (eye) icon attached to that dataset.
> 2. Click into the middle panel, which should now display the content of the
>    dataset, select all sequences by pressing <kbd>Ctrl</kbd>+<kbd>A</kbd>,
>    then copy the selection to the clipboard with <kbd>Ctrl</kbd>+<kbd>C</kbd>
> 3. Head over to the
>    [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=MegaBlast&PROGRAM=blastn)
>    service, paste the copied content into the `Enter Query Sequence` text
>    box and click `BLAST` at the bottom of the form.
> 4. Wait for the *BLAST* run to finish.
> 5. On the results page, look for the drop-down menu next to `Results for`.
>    It lets you toggle the BLAST hits list further down to include only the
>    matches to individual sequences from your multi-sequence query.
>
{: .hands_on}

Now take a bit of time to explore the BLAST hits uncovered for some of your
assembled nodes. Pay attention, specifically, to the node with the longest
sequence (node #1) and the circular node #2, but also investigate the results
for a few others.

> <question-title></question-title>
>
> 1. Which genome is represented by node #1?
> 2. Which genome corresponds to node #2? Does this finding remind you of
>    something you have learnt before?
> 3. What do most other node sequences have in common? Do these additinal
>    findings make sense?
>
> > <solution-title></solution-title>
> > 1. The sequence of node #1 is the assembled SARS-CoV-2 sequence we are
> >    looking for. It is a perfect match to various SARS-CoV-2 genome
> >    sequences found in Genbank over the entire assembled length, and we have
> >    been able to assemble almost the entire genome even from the subsampled
> >    sequencing data.
> > 2. The circular sequence of node #2 corresponds to the 5,386 bp genome of
> >    bacteriophage phiX174. As explained in the more general
> >    [Unicycler assembly](../unicycler-assembly/tutorial.html) tutorial,
> >    this genome is often used as a spike-in in Illumina sequencing.
> >    Finding the complete sequence here is, thus,  another indication that
> >    our analysis worked and produced meaningful results.
> > 3. Almost all other assembled sequences appear to represent parts of
> >    bacterial genomes. The only exceptions are the node #10 sequence, for
> >    which no significant BLAST hits could be found, and the node #11
> >    sequence, which represents a small stretch of left-over human genomic
> >    DNA, which seems to have survived our subtraction approach.
> >
> >    What all the bacterial genomes have in common is that they represent
> >    genera of bacteria that are known to colonize the oral cavity and
> >    mucosa. Since all samples in this analysis are BALF samples the presence
> >    of DNA from such bacteria should not be surprising. In addition, some
> >    members of the identified genera are known as opportunistic pathogens.
> >    In particular, members of the genus *Prevotella* can infect the
> >    respiratory tract and contribute to inflammation under anaerobic
> >    conditions caused by primary infections. Hence, an alternative
> >    explanation for the presence of some of these sequences in the samples
> >    might be that the corresponding bacteria contributed to the clinical
> >    picture of some of the Covid-19 patients they were obtained from.
> {: .solution}
>
{: .question}

# Conclusion


The power of modern genome assembly tools is remarkable, and so is their
robustness in the face of data of metagenomic nature. Assembling reads derived
from a virus and a good handful of copurified bacteria back into separate
contigs is a challenging task, which Unicycler solved without major issues!

Nevertheless, good quality assemblies still rely on proper preprocessing and
filtering to reduce the number of misassembly events, ambiguous assemblies and
of incorporation of sequencing errors into the final assembly.
