---
layout: tutorial_hands_on

title: Detection of SARS-CoV-2 variants by genomic analysis of wastewater samples
level: Intermediate
zenodo_link: "https://zenodo.org/record/7468942"
questions:
- How to detect SARS-CoV-2 lineages abundances in wastewater samples in Galaxy?
- How to learn other species diversity in wastewater samples?
- Which tools and workflows can we use to identify SARS-CoV-2 lineages in wastewater samples in Galaxy?
objectives:
- Repeat SARS-CoV-2 data preparation
- Select and run workflow to define lineages abundances of SARS-CoV-2 in wastewater samples
- Run workflow to summarize and generate report for previously called allelic variants
- Interpret results
time_estimation: 3H
key_points:
- 2 specialized workflows are available for the identification of lineages abundances of SARS-CoV-2 from raw wastewater sequencing data depending on the exact type of input
- Data from batches of samples can be processed in parallel using collections
contributors:
- plushz
- bebatut
- wm75
tags:
- covid19
- wastewater
---


# Introduction

Millions of people have been affected by the COVID-19 pandemic after the first report of SARS-CoV-2 in Wuhan, China. During the COVID-19 pandemic, wastewater surveillance has received extensive public attention as a passive monitoring system that complements clinical and genomic surveillance. The detection and quantification of viral RNA in wastewater samples are already possible through several methods and protocols, and viral RNA concentrations in wastewater have been shown to correlate with reported cases. The wastewater surveillance methods are a good solution for enabling early, economical, and efficient detection so that public-health measures can be implemented as soon as they are necessary.

![Upper branch shows a clinical surveillance of sars-cov-2, from the infection moment to bioinformatics data analysis; lower branch, in turn, represents wastewater surveillance.](./images/ww-process-ww.png "Schematic diagram shows the process of detecting viruses by wastewater surveillance against clinical surveillance.")

Two sequencing platforms (e.g. Illumina and Oxford Nanopore) in combination with several established library preparation (e.g. ampliconic and metatranscriptomic) strategies are predominantly used to generate SARS-CoV-2 sequence data. However, data alone do not equal knowledge: they need to be analyzed. The Galaxy community has developed workflows to perform bioinformatics analysis.

> <details-title>Further reading</details-title>
> More information about the workflows, including benchmarking, can be found
> - on the Galaxy Covid-19 effort [website](https://galaxyproject.org/projects/covid19/)
{: .details}

This tutorial will teach you how to obtain and run these workflows appropriately for different types of input data, be it:

- Paired-end data derived from Illumina-based Metatranscriptomic experiments, or
- Paired-end data generated with Illumina-based Ampliconic (ARTIC) protocols

Both workflows are repurposed and improved Galaxy workflows for clinical data. To learn more about workflows for clinical data you can follow this [dedicated tutorial]({% link topics/variant-analysis/tutorials/sars-cov-2-variant-discovery/tutorial.md %}).

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Prepare Galaxy and data

Any analysis should get its own Galaxy history. So let's start by creating a new one:

> <hands-on-title>Prepare the Galaxy history</hands-on-title>
>
> 1. Create a new history for this analysis
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Rename the history
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
{: .hands_on}


## Get sequencing data

Before we can begin any Galaxy analysis, we need to upload the input data: FASTQ files with the sequenced wastewater samples cantained SARS-CoV-2. Several types of data are possible:

- Paired-end data derived from Illumina-based RNAseq experiments
- Paired-end data generated with Illumina-based Ampliconic (ARTIC) protocols

We encourage you to use your own data here (with at least 2 samples). If you do not have any datasets available, we provide some example datasets: 1) paired-end data generated with Illumina-based metatranscriptomic protocol from [ENA](https://www.ebi.ac.uk/ena/browser/view/PRJNA661613), European Nucleotide Archive; 2) [paired-end data generated synthetically with Illumina-based Ampliconic (ARTIC) protocols](https://github.com/suskraem/ww_benchmark).

There are several possibilities to upload the data depending on how many datasets you have and what their origin is:

- Import datasets

  - from your local file system,
  - from a given URL or
  - from a shared data library on the Galaxy server you are working on

  and organize the imported data as a dataset collection.

  > <comment-title>Collections</comment-title>
  >
  > A dataset collection is a way to represent an arbitrarily large collection of samples as a singular entity within a user's workspace. For an in-depth introduction to the concept you can follow this [dedicated tutorial]({% link topics/galaxy-interface/tutorials/collections/tutorial.md %}).
  >
  {: .comment}

- Import from [NCBI's Sequence Read Archive (SRA) at NCBI](https://www.ncbi.nlm.nih.gov/sra) with the help of a dedicated tool, which will organize the data into collections for you.

   > <comment-title>Getting data from SRA</comment-title>
   >
   > [A dedicated tutorial is available to explain how to find and import SARS-CoV-2 data from SRA]({% link topics/variant-analysis/tutorials/sars-cov-2/tutorial.md %}).
   >
   {: .comment}

> <hands-on-title>Import datasets</hands-on-title>
>
> 1. Import the datasets
>
>    - Option 1 [{% icon video %}](https://youtu.be/FFCDx1rMGAQ): Your own local data using **Upload Data** (recommended for 1-10 datasets). 
>
>      {% snippet faqs/galaxy/datasets_upload.md %}
>
>    - Option 2 [{% icon video %}](https://youtu.be/hC8KSuT_OP8): Your own local data using **FTP** (recommended for >10 datasets)
>
>      {% snippet faqs/galaxy/datasets_upload_ftp.md %}
>
>    - Option 3: From an external server via URL
>
>      {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>      For our example datasets, the datasets are stored on [Zenodo]({{ page.zenodo_link }}) and can be retrieved using the following URLs:
>      For metatranscriptomic-illumina dataset:
>
>      ```
>      {{ page.zenodo_link }}/files/SRR12596165_1.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596165_2.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596166_1.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596166_2.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596167_1.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596167_2.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596168_1.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596168_2.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596169_1.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596169_2.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596170_1.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596170_2.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596171_1.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596171_2.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596172_1.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596172_2.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596173_1.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596173_2.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596174_1.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596174_2.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596175_1.fastq.gz
>      {{ page.zenodo_link }}/files/SRR12596175_2.fastq.gz
>      ```
>      For ampliconic-illumina dataset:
>
>      ```
>      https://zenodo.org/record/7469383/files/sample1_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample1_R2.fastq.gz
>      https://zenodo.org/record/7469383/files/sample2_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample2_R2.fastq.gz
>      https://zenodo.org/record/7469383/files/sample3_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample3_R2.fastq.gz
>      https://zenodo.org/record/7469383/files/sample4_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample4_R2.fastq.gz
>      https://zenodo.org/record/7469383/files/sample5_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample5_R2.fastq.gz
>      https://zenodo.org/record/7469383/files/sample6_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample6_R2.fastq.gz
>      https://zenodo.org/record/7469383/files/sample7_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample7_R2.fastq.gz
>      https://zenodo.org/record/7469383/files/sample8_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample8_R2.fastq.gz
>      https://zenodo.org/record/7469383/files/sample9_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample9_R2.fastq.gz
>      https://zenodo.org/record/7469383/files/sample10_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample10_R2.fastq.gz
>      https://zenodo.org/record/7469383/files/sample11_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample11_R2.fastq.gz
>      https://zenodo.org/record/7469383/files/sample12_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample12_R2.fastq.gz
>      https://zenodo.org/record/7469383/files/sample13_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample13_R2.fastq.gz
>      https://zenodo.org/record/7469383/files/sample14_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample14_R2.fastq.gz
>      https://zenodo.org/record/7469383/files/sample15_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample15_R2.fastq.gz
>      https://zenodo.org/record/7469383/files/sample16_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample16_R2.fastq.gz
>      https://zenodo.org/record/7469383/files/sample17_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample17_R2.fastq.gz
>      https://zenodo.org/record/7469383/files/sample18_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample18_R2.fastq.gz
>      https://zenodo.org/record/7469383/files/sample19_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample19_R2.fastq.gz
>      https://zenodo.org/record/7469383/files/sample20_R1.fastq.gz
>      https://zenodo.org/record/7469383/files/sample20_R2.fastq.gz
>      ```
>
> 2. Create a collection to organize the data
>
>      {% snippet faqs/galaxy/collections_build_list_paired.md %}
>
>      For the example datasets:
>      - Since the datasets carry `_1` and `_2` or `_R1` and `_R2` in their names, Galaxy may already have detected a possible pairing scheme for the data, in which case the datasets will appear in green in the lower half (the paired section) of the dialog.
>
>        You could accept this default pairing, but as shown in the middle column of the paired section, this would include the `.fastqsanger` suffix in the pair names (even with `Remove file extensions?` checked Galaxy would only remove the last suffix, `.gz`, from the dataset names.
>
>        It is better to undo the default pairing and specify exactly what we want:
>        - at the top of the *paired section*: click `Unpair all`
>
>          This will move all input datasets into the *unpaired section* in the upper half of the dialog.
>        - set the text of *unpaired forward* to: `_1.fastq.gz` or `_R1.fastq.gz`
>        - set the text of *unpaired reverse* to: `_2.fastq.gz` or `_R2.fastq.gz`
>        - click: `Auto-pair`
>
>        All datasets should be moved to the *paired section* again, but the middle column should now show that only the sample accession numbers will be used as the pair names.
>
>      - Make sure *Hide original elements* is checked to obtain a cleaned-up history after building the collection.
>      - Click *Create Collection*
>
{: .hands_on}

> <comment-title>Learning to build collections automatically</comment-title>
>
> It is possible to build collections from tabular data containing URLs, sample sheets, list of accessions or identifiers, etc., directly during upload of the data. [A dedicated tutorial is available to explain the different possibilities]({% link topics/galaxy-interface/tutorials/upload-rules/tutorial.md %}).
>
{: .comment}

## Import auxiliary datasets

Besides the sequenced reads data, we need at least two additional datasets for calling variants and annotating them:

- the SARS-CoV-2 reference sequence [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=fasta) to align and compare our sequencing data against

- a tabular dataset defining aliases for viral gene product names, which will let us translate NCBI RefSeq Protein identifiers (used by the SnpEff annotation tool) to the commonly used names of coronavirus proteins and cleavage products.

Another three datasets are needed only for the analysis of ampliconic, e.g. ARTIC-amplified, input data:

- a BED file specifying the primers used during amplification and their binding sites on the viral genome
- a custom tabular file describing the amplicon grouping of the primers
- bedfile defining the amplicons that is used by Cojac: mutbamscan


> <hands-on-title>Import auxiliary datasets</hands-on-title>
>
> 1. Import the auxiliary datasets:
>    - the SARS-CoV-2 reference (`NC_045512.2_reference.fasta`)
>    - BED file containing ARTIC primer positions (`ARTIC_primer_BED.bed`)
>    - ARTIC primer amplicon grouping info (`ARTIC_primers_to_amplicon_assignments.bed`)
>    - bedfile defining the amplicons that is used by Cojac: mutbamscan (`BED_defining_amplicons_for_COJAC.bed`)
>
>    > <details-title>Not using ARTIC v4.1 amplified sequencing data?</details-title>
>    >
>    > The instructions here assume you will be analyzing the example samples
>    > suggested above, which have been amplified using version 4.1 of the ARTIC
>    > network's SARS-CoV-2 primer set. If you have decided to work through
>    > this tutorial using your own samples of interest, and if those samples
>    > have been amplified with a different primer set, you will have to upload
>    > your own datasets with primer and amplicon information at this point.
>    > If the primer set is from the ARTIC network, just not version 4.1, you
>    > should be able to obtain the primer BED file from
>    > [their SARS-CoV-2 github repo](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019).
>    {: .details}
>
>    Import these datasets from [Zenodo](https://zenodo.org/record/7469383)
>
>      ```
>      https://zenodo.org/record/7469383/files/NC_045512.2_FASTA_sequence_of_SARS-CoV-2.fasta
>      https://zenodo.org/record/7469383/files/ARTIC_primer_BED.bed
>      https://zenodo.org/record/7469383/files/ARTIC_primers_to_amplicon_assignments.bed
>      https://zenodo.org/record/7469383/files/BED_defining_amplicons_for_COJAC.bed
>      ```
>
>      {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    For the example datasets, you will need to import all 4 auxiliary datasets.
>
> 2. Check and manually correct assigned datatypes
>
>    If you have imported the auxiliary datasets via their Zenodo links, Galaxy
>    will have tried to autodetect the format of each imported dataset, but
>    will not always be right with its guess. It's your task now to check and
>    possibly correct the format assignments for each of the datasets!
>
>    - Expand the view of each of the uploaded auxiliary datasets and see if
>      Galaxy shows the following `format` values:
>      - for `NC_045512.2_FASTA_sequence_of_SARS-CoV-2.fasta`: `fasta`
>      - for `ARTIC_primer_BED.bed`: `bed`
>      - for `ARTIC_primers_to_amplicon_assignments.bed`: `bed` or `tabular`
>      - for `BED_defining_amplicons_for_COJAC.bed`: `bed`
>
>    - If any of the above assignments are not what they should be, then change
>      the datatype of the corresponding dataset now to the intended format.
>
>      {% snippet faqs/galaxy/datasets_change_datatype.md %}
>
>    If you have imported the auxiliary datasets into your history from a
>    shared data library or history, then the above steps are not necessary
>    (though checking the datatypes of imported data is good practice in
>    general) because the shared datasets have their format configured
>    correctly already.
>
{: .hands_on}


# From FASTQ to SARS-CoV-2 lineages abundances

Before identifying lineages aboundances, several steps have to be done to preprocess data.

Quality control step is often used because there is no perfect sequencing technology, and each
instrument will generate different types and amounts of errors, such as incorrect nucleotide calls.
Each sequencing platform has technical limitations that result in these incorrectly called bases.
Thus, it is important to identify and exclude error types that may affect downstream analysis
interpretation. As a result, sequence quality control is an essential first step in the analysis
process.

Another step, primer trimming, is a specific step for datasets generated with ARTIC protocol.
The auxiliary file is used for this step - a BED file specifying the primers used during amplification
and their binding sites on the viral genome. Primer trimmer uses primer positions supplied in
a BED file to soft clip primer sequences from an aligned and sorted BAM file. Following this,
the reads are trimmed based on a quality threshold. More specifically, some primer trimmers,
in order to do quality trimming, use a sliding window approach. The window slides from the
5’ end to the 3’ end and if at any point the average base quality in the window falls below
the threshold, the remaining read is softly clipped. If after trimming, the length of the read is
greater than the minimum length specified, the read is written to the new trimmed BAM file. It
should be noted, for datasets that were not generated with primer-based protocol like ARTIC,
this primer-trimming step is not applicable.

Moreover, adapter trimming step is processed. For instance, upon Illumina sequencing we receive
raw reads with adapters at 3’ end. The adapters contain the sequencing primer binding sites,
the index sequences, and the sites that allow library fragments to attach to the flow cell lawn.
This might influence a downstream analysis, thus, adapter trimming is required.

A decontamination step can then be included to remove reads from the human genome, since
viral sequence data from clinical samples commonly contain human contamination. Prior to
sharing, it needs to be removed for legal and ethical reasons as well as to speed up downstream
analysis {% cite hunt2022 %}.

The crucial step is mapping with reference SARS-CoV-2 sequence NC_045512.2 that is publicly available in NCBI database {% cite genome2020 %}. A mapping tool of choice can differ from one pipeline to another,
depending on read length, sequencing technology, and other factors.

Some pipeline steps are not always included in pipelines, such as removing duplicates. This step
can be important for Illumina sequencing reads. During the sequencing process with Illumina
sequencing technology, some duplicate reads/sequences can be produced, which can create bias
in downstream analyses. It is, therefore, possible to remove duplicates or mark them without
removing them. When removing duplicates, one should be certain that they are duplicates and
not repeated regions. It can therefore be reasonable to keep duplicates marked rather than
remove them, as this can be useful for downstream analysis.

Another step, which is not present everywhere, is helpful due to potential ambiguity, while indels
are not parsed when they overlap the beginning or end of alignment boundaries. Input insertions
and deletions must be homogenized with left realignment in order to gain a more homogeneous
distribution. Left realignment will place all indels in homopolymer and microsatellite repeats at
the same position, provided that doing so does not introduce mismatches between the read and
reference other than the indel {% cite garrison2012 %}. Basically, this step is considered to correct mapping errors
and prepare reads for further variant calling.

Additionally, realigned reads can be taken and checked for the quality of alignment using bioinformatics
tools (e.g., Qualimap {% cite qualimap %}, {% cite qualimap2 %}). Based on the features of the mapped reads, it analyzes SAM/BAM alignment data and provides a global picture of the data that can help detect biases
in sequencing and/or mapping of the data and ease decision-making for further analysis.
After mapping and other additional preparation steps, variant calling should be run where variants
from sequence data are identified.

![Here is simplified process of bioinformatics steps used to analyze sequenced data for sars-cov-2 surveillance. Tools can differ from one pipeline to another. But the main steps, in general, are more or less the same. Raw data are sequencing data. Then, primer trimming is a specific step for ampliconic datasets. The auxiliary file is used for this step - a BED file specifying the primers used during amplification. Variant calling should be run where variants from sequence data are identified. Variant calling step is followed by mutation annotation. The data is not changed; here, only format is changed to be more readable](./images/sars-surveillance-bioinf-last.png "Main steps to be done for bioinformatics of SARS-CoV-2 surveillance.")

We offer two workflow with their tools and parameters optimized for different types of input data as outlined in the following table:

Workflow version 	| Input data | Read aligner | Variant caller | Delineation tool
--- | --- | --- | --- | ---
Illumina metatranscriptomic PE | Paired-end data derived from RNAseq experiments | **bwa-mem** {% cite li_2010 %} | **lofreq** {% cite wilm_2012 %} | Freyja {% cite karthikeyan2022 %}
Illumina ampliconic PE | Paired-end data generated with Illumina-based Ampliconic (ARTIC) protocols | **bwa-mem** {% cite li_2010 %} | **lofreq** {% cite wilm_2012 %} | Freyja {% cite karthikeyan2022 %}, COJAC {% cite cojac2022 %}

> <details-title>About the workflows</details-title>
>
> - The Illumina metatranscriptomic workflow (Illumina metatranscriptomic PE) perform read mapping with **bwa-mem** and **bowtie2**, respectively, followed by sensitive allelic-variant (AV) calling across a wide range of AFs with **lofreq**. Computation of SARS-CoV-2 lineages abundances is performed with **Freyja**.
> - The workflow for Illumina-based ARTIC data (Illumina ARTIC) builds on the RNASeq workflow for paired-end data using the same steps for mapping (**bwa-mem**) and AV calling (**lofreq**), but adds extra logic operators for trimming ARTIC primer sequences off reads with the **ivar** package. In addition, this workflow uses **ivar** also to identify amplicons affected by ARTIC primer-binding site mutations and excludes reads derived from such “tainted” amplicons when calculating alternative allele frequences (AFs) of other AVs. Computation of SARS-CoV-2 lineages abundances is performed with **Freyja** and **COJAC** in two different branches.
>
{: .details}

> <hands-on-title>From FASTQ to SARS-CoV-2 lineages abundances</hands-on-title>
>
> 1. **Get the workflow** for your data into Galaxy 
>
>    - Option 1: Find workflows on the [WorkflowHub](https://workflowhub.eu) and run them directly on [usegalaxy.eu](https://usegalaxy.eu/)
>
>      - Open the workflow page on the WorkflowHub
>        - [Illumina ampliconic PE](https://workflowhub.eu/workflows/)
>        - [Illumina metatranscriptomic PE](https://workflowhub.eu/workflows/)
>
>      - Click on `Run on usegalaxy.eu` at the top right of the page
>      
>        The browser will open a new tab with Galaxy's workflow invocation interface.
>
>    - Option 2: Import the workflow via its github link
>
>      - Open the GitHub repository of your workflow
>        - [Illumina ampliconic PE](https://github.com/iwc-workflows/sars-cov-2-wastewater-pe-illumina-artic-variant-analysis)
>        - [Illumina metatranscriptomic PE](https://github.com/iwc-workflows/sars-cov-2-wastewater-pe-illumina-metatranscriptomic-variant-analysis)
>      - Open the `.ga` file
>      - Click on `Raw` at the top right of the file view
>      - Save the file or Copy the URL of the file
>      - Import the workflow to Galaxy
>
>        {% snippet faqs/galaxy/workflows_import.md %}
>
> 2. Run **SARS-CoV-2: lineages analysis on ...** {% icon workflow %} using the following parameters:
>
>    {% snippet faqs/galaxy/workflows_run.md %}
>
>    - *"Send results to a new history"*: `No`
>
>    - For **Illumina ampliconic PE** workflow (named **WW-SARS-CoV-2: lineages analysis on ARTIC PE data**),  *to use for example datasets*
>      - {% icon param-collection %} *"1: Paired Collection"*: paired collection created for the input datasets
>      - {% icon param-file %} *"2: NC_045512.2 FASTA sequence of SARS-CoV-2"*: `NC_045512.2_FASTA_sequence_of_SARS-CoV-2.fasta`
>      - {% icon param-file %} *"3: ARTIC primer BED"*: `ARTIC_primer_BED.bed`
>      - {% icon param-file %} *"4: ARTIC primers to amplicon assignments"*: `ARTIC_primers_to_amplicon_assignments.bed`
>      - {% icon param-file %} *"5: BED defining amplicons for COJAC"*: `BED_defining_amplicons_for_COJAC.bed`
>
>    - For **Illumina metatranscriptomic PE** workflow (named **WW-SARS-CoV-2: lineages analysis on Metatranscriptomics PE data**)
>      - {% icon param-collection %} *"1: Paired Collection"*: paired collection created for the input datasets
>      - {% icon param-file %} *"2: NC_045512.2 FASTA sequence of SARS-CoV-2"*: `NC_045512.2_FASTA_sequence_of_SARS-CoV-2.fasta`
>
{: .hands_on}

The execution of the workflow takes some time. It is possible to launch the next step even if it is not done, as long as all steps are successfully scheduled.


# Results interpretation

Once the jobs of previous workflows are done, we can look at the results.

**Freyja** produces the following reports

1. **Aggregated data**: This table provides a aggregated information about lineages abundances for all samples. Each line in the table represents one sample.

   Column | Field | Meaning
   --- | --- | ---
   1 | `Sample` | sample name
   2 | `summarized` | denotes a sum of all lineage abundances in a particular WHO designation (i.e. B.1.617.2 and AY.6 abundances are summed in the above example), otherwise they are grouped into "Other"
   3 | `lineages` | lists the identified lineages in descending order
   4 | `abundances` | contains the corresponding abundances estimates
   5 | `resid` | corresponds to the residual of the weighted least absolute devation problem used to estimate lineage abundances
   6 | `coverage` | provides the 10x coverage estimate (percent of sites with 10 or greater reads

   > <question-title></question-title>
   > 
   > 1. How many lineages were identified in sample SRR12596165 in metatranscriptomic-illumina dataset provided for this tutorial?
   > 2. Which proportion of B.1.533 lineage was identified in sample SRR12596170?
   >
   > > <solution-title></solution-title>
   > >
   > > 1. In the second line of the table we see 7 lineages listed in column "lineages": B.10 B.47 B.23 B.26 B.1.14 B.20 B
   > >
   > > 2. In the 9th line of the table that corresponds to SRR12596170 sample, B.1.533 lineage is in the first position in the column "lineages". Then we look at the next column "abundances" and see in the first position is the proportion of B.1.533 which is 0.15773800.
   > >
   > {: .solution}
   {: .question}

2. **Lineages abundances plot**: This plot provides a fractional abundance estimate for all aggregated samples. Each bar in the plot represents one sample. Different colors represent different lineages.

   > <question-title></question-title>
   > 
   > 1. Which WHO designated variant is prevalent in sample1 from ampliconic-illumina dataset provided for this tutorial?
   > 1. What other lineages are present in this sample?
   >
   > > <solution-title></solution-title>
   > >
   > > 1. Omicron
   > >
   > > 2. Delta and Other
   > >
   > {: .solution}
   {: .question}

**COJAC** produces the following report

1. **Mutation cooccurrence**: This table provides a aggregated information about lineages abundances and mutations for all samples grouped by amplicons. Each line in the table represents one sample.

   Column | Field | Meaning
   --- | --- | ---
   1 | `Sample` | sample name
   2 | `count` | total count of amplicons carrying the sites of interest
   3 | `mut_all` | amplicons carrying mutations on all site of interest (e.g.: variant mutations observed on all sites)
   4 | `mut_oneless` | amplicons where one mutation is missing (e.g.: only 2 out of 3 sites carried the variant mutation, 1 sites carries wild-type)
   5 | `frac` | fraction (mut_all/count) or empty if no counts
   6 | `cooc` | number of considered site (e.g.: 2 sites of interests) or empty if no counts

   > <question-title></question-title>
   > 
   > 1. Which fraction of BA.1 lineage is identified on amplicon A72 in sample1 from ampliconic-illumina dataset provided for this tutorial?
   > 2. How many amplicons carrying mutations on all site of interest are identified in sample 3?
   >
   > > <solution-title></solution-title>
   > >
   > > 1. 0.8181818182
   > >
   > > 2. 3
   > >
   > {: .solution}
   {: .question}



# Conclusion


In this tutorial, we used Galaxy workflows for the detection and interpretation of lineages abundances in SARS-CoV-2 wastewater samples.

The workflows can be freely used and immediately accessed from the three global Galaxy instances. Each is capable of supporting thousands of users running hundreds of thousands of analyses per month. 

It is also possible to automate the workflow runs using the command line as explained in [a dedicated tutorial]({% link topics/galaxy-interface/tutorials/workflow-automation/tutorial.md %}).
