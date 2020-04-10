---
layout: tutorial_hands_on

title: "Unicycler assembly after preprocessing to remove reads from a contaminating genome"
zenodo_link: "https://doi.org/10.5281/zenodo.3732358"
questions:
  - How can I assemble a genome of interest against a background of contaminating reads from another genome?
objectives:
  - Obtain viral (SARS-CoV-2) sequencing data with contaminating human reads from public sources
  - Detect and remove human reads
  - Regenerate fastq data of remaining reads for downstream analysis
time_estimation: "4h"
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
    topic_name: galaxy-data-manipulation
    tutorials:
      - collections
      - upload-rules
tags:
  - covid19
contributors:
  - wm75

---

# Introduction
{:.no_toc}

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

> ### {% icon comment %} The usegalaxy.* COVID-19 analysis project
> This tutorial uses the same data as, and recapitulates to a large extent, the
> [Pre-processing](https://covid19.galaxyproject.org/genomics/1-PreProcessing/)
> and [Assembly](https://covid19.galaxyproject.org/genomics/2-Assembly/) steps
> of the [Genomics](https://covid19.galaxyproject.org/genomics/) section of
> [covid19.galaxyproject.org](https://covid19.galaxyproject.org/).
>
{: .comment}

> ### Agenda
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
short-reads archive (SRA).

Since automated data downloads from the SRA can be unreliable at times and
could get broken by changes to the download interface on the NCBI side, this
tutorial offers two ways to access the sequenced reads input data:

1. Direct download from the NCBI SRA based on accession numbers and using the
   dedicated **Faster Download** tool
   
   Use this method if it works and is fast enough for you, and if you are
   interested in learning to obtain short-reads data directly from NCBI, in
   general.

2. Download of the same data deposited as a copy at [Zenodo](https://zenodo.org)

   This method uses Galaxy's generic data import functionality, and should be
   very reliable and faster than the download from NCBI.
   It also showcases **rule-based** uploads and demonstrates how they can be
   used to download several datasets and to arrange them into easy to handle
   data structures at the same time.
   
   Use this method if the direct download from the NCBI SRA does not work, or
   is too slow for your time frame, or if you are interested in advanced use
   of Galaxy's data import functionality.
   
The corresponding two step-by-step instructions below have been crafted to
produce identically arranged data structures in your history so all subsequent
steps are independent of the data source you choose.

## Get data from NCBI SRA

> ### {% icon hands_on %} Hands-on: Data upload to Galaxy from NCBI SRA
>
> 1. Create a new history for this tutorial and give it a proper name
>
>    {% include snippets/create_new_history.md %}
>    {% include snippets/rename_history.md %}
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
>    {% include snippets/create_new_file.md format="tabular" %}
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
>    > ### {% icon comment %} Name tags in the analysis
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
>    {% include snippets/add_tag.md %}
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

> ### {% icon hands_on %} Hands-on: Data upload to Galaxy from Zenodo
>
> 1. Create a new history for this tutorial and give it a proper name
>
>    {% include snippets/create_new_history.md %}
>    {% include snippets/rename_history.md %}
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
>        - *"Name"*: `Illumnia PE data` (or similar)
>        - *"Add nametag for name:"* {% icon param-check %}
>
>          > ### {% icon comment %} Name tags in the analysis
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

1. Trimming of Illumina reads with fastp / visualizing results with MultiQC

   **TO DO**

2. QC of Nanopore reads with NanoPlot and/or FastQC

   **TO DO**

# Subtraction of reads mapping to the human reference genome

## Map Illumina reads with bowtie2

In this tutorial, we are using **Bowtie2** for mapping our short-reads data to
the human genome. *BWA

> ### {% icon hands_on %} Hands-on: Mapping with Bowtie2
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
        `Use a built-in genome index`
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
>    > ### {% icon question %} Questions
>    >
>    > 1. What percentage of reads of each sample has been aligned to the `hg38` reference genome?
>    > 2. Which sample is the least contaminated with human reads?
>    > 3. Which sample contains the highest amount of SARS-CoV2 reads?
>    >
>    > > ### {% icon solution %} Solution
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

## Map Nanopore reads with minimap2

> ### {% icon hands_on %} Hands-on: Nanopore reads mapping
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
>    - *"Select analysis mode (sets default)"*: `Oxford Nanopore read to reference mapping. ...`
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
>    > ### {% icon comment %} Mapping stats for Nanopore-sequenced long reads
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
   
> ### {% icon hands_on %} Hands-on: Mapped reads filtering
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
>    > ### {% icon details %} Why do this?
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
>    {% include snippets/change_dbkey.md dbkey="unspecified (?)" %}
>
{: .hands_on}
   

# Format conversion of remaining reads

## Conversion to fastq with samtools fastx

> ### {% icon hands_on %} Hands-on: BAM to fastq format conversion
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

> ### {% icon hands_on %} Hands-on: Arrange to list collections into a list of pairs
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

   required by unicycler

> ### {% icon hands_on %} Hands-on: Collapsing each collection into a single dataset
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
subsequent assembly step to approximately 2 hours.

> ### {% icon comment %} If you are in a hurry
> The downsampling parameters below have been chosen to have minimal impact on
> the assembly results. Further speed-ups are certainly possible, but will
> likely lead to different assembly outcomes.
{: .comment}

> ### {% icon hands_on %} Hands-on: Subsampling of paired-end short-reads data
>
> 1. **seqtk_sample** {% icon tool %} with the following parameters
>    - {% icon param-files %} *"Input FASTA/Q file"*: The two datasets with the
>      combined Illumina-sequenced forward and reverse reads, outputs of the
>      first and the second run of **Collapse Collection**
>    - *"RNG seed"*: 4
>    - *"Subsample (decimal fraction or number)"*: 0.2
>
{: .hands_on}

## Create assembly

> ### {% icon hands_on %} Hands-on: Assembly of SARS-CoV2 genome
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
>      `10000`
>    - *"The expected number of linear (i.e. non-circular) sequences in the assembly"*:
>      `1`
>
{: .hands_on}

## Explore assembly

### Assembly inspection with Bandage

> ### {% icon hands_on %} Hands-on: Assembly stats and visualization with Bandage
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
{: .hands_on}

### Check origin of assembled sequences with BLAST

**TO DO** (point out 100% identity to published SARS-CoV-2 sequence of largest
assembled sequence. Discuss second assembled sequence (*Prevonella* species as
opportunistic pathogens)

# Conclusion
{:.no_toc}

**TO DO**
