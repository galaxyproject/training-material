---
layout: tutorial_hands_on

title: Using SamestrGal to identify shared strains in FMT-treated rCDI samples
level: Intermediate
zenodo_link: https://zenodo.org/records/20745835
questions:
- How can shared microbial strains between metagenomic samples be identified?
- What information is needed to determine whether two samples share the same strain rather than just the same species?
- How can shared strain patterns show donor engraftment and strain persistence in fecal microbiota transplantation (FMT)?
objectives:
- Run the SamestrGal workflow in Galaxy to detect shared microbial strains across metagenomic samples
- Explain the role of each tool in the SameStr workflow
- Interpret the shared strain outputs produced by SameStr Summarize
- Distinguish between strain engraftment and persistence in FMT-treated samples based on SamestrGal's output
time_estimation: 2H
key_points:
- SamestrGal enables shared strain detection between metagenomic samples
- SamestrGal can be run as a complete workflow or tool by tool to understand each analysis step
- Shared strains are identified by comparing the genetic variation of species-specific marker genes between samples
- In FMT-treated rCDI samples, shared strains between donor and post-FMT samples indicate engraftment, while shared strains between pre-FMT and post-FMT samples indicate persistence

contributions:
  authorship:
    - xens25

---


# Introduction

<!-- This is a comment. -->

*Clostridium difficile* is a pathogen found in the gut that can proliferate once antibiotics remove its competing bacteria, leading to recurrent *Clostridium difficile* infection (rCDI) {% cite cdifficile %}. Fecal microbiota transplantation (FMT) treats rCDI by restoring a donor's balanced gut microbiota in the patient {% cite zanellaterrier2014recurrent %}. Since bacterial strains within the same species can behave differently, confirming that a bacterium detected in the patient after treatment is the same strain that came from the donor, rather than a different strain of the same species that was already present, requires strain-level resolution rather than species-level identification alone {% cite smillie2018strainfinder %}.

SameStr is a bioinformatic tool that detects shared microbial strains across metagenomic samples by comparing the genetic variation observed at species-specific marker genes. In FMT studies, this allows tracking of strain persistence and engraftment between donor, pre-FMT, and post-FMT samples {% cite samestr_podlesny %}.

> <details-title>Key terms used in this tutorial</details-title>
>
> - **rCDI**: recurrent *Clostridium difficile* infection, caused by *C. difficile* spores that survive antibiotic treatment and later reactivate.
> - **FMT**: fecal microbiota transplantation, a treatment for rCDI that restores a donor's gut microbiota in the patient.
> - **Strain**: a classification within a bacterial species that exhibits distinct genotypic or phenotypic traits.
> - **Clade**: a group of organisms in a phylogenetic tree that share a common ancestor and all of its descendants.
> - **SNV profile**: the frequency of each nucleotide observed across the aligned reads at every covered position of a clade's marker genes in a sample, used by SameStr to compare samples.
> - **MVS (Maximum Variant Profile Similarity)**: the fraction of shared positions with enough coverage between two samples that contain at least one matching nucleotide.
> - **Persistence**: a strain detected in a patient before FMT (pre-FMT) that is still detected in the same patient after treatment (post-FMT).
> - **Engraftment**: a strain from the donor that is detected in the patient's sample after FMT (post-FMT).
{: .details}


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# SamestrGal Workflow Overview

SamestrGal uses metagenomic samples as input and runs a sequence of tools to detect shared strains between them. The reads are first trimmed and quality-filtered with KneadData, then aligned to species-specific marker genes with MetaPhlAn, which also produces a taxonomic profile for each sample. The resulting alignments are then processed by the SameStr tools, which convert them into a species-specific SNV profile for each clade detected in a sample. These SNV profiles are merged across samples for each clade. The Maximum Variant Profile Similarity (MVS) is then calculated between samples for each clade. Thresholds for the overlap and the MVS between samples are also set, and the results are reported in a table across all sample comparisons. Based on these thresholds, the table indicates whether each clade is classified as a shared strain.

![Diagram of the SamestrGal workflow](../../images/samestr-workflow.png "Overview of the SamestrGal workflow.")

The following sections show how to run the entire SameStr workflow with SamestrGal, and then go through each tool individually to explain what it does, and how it contributes to the entire workflow.


## Get data

The following steps use metagenomic shotgun sequencing data from FMT-treated rCDI patient stool samples, part of the FRICKE cohort {% cite samestr_podlesny %}. Three paired-end samples are used: a donor sample, a sample collected from the patient before treatment (pre-FMT), and a sample collected from the patient after treatment (post-FMT). The reads used in this tutorial have been reduced from the original dataset to minimize the time needed to complete the workflow.


> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial and name it (e.g. "SamestrGal tutorial")
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    https://zenodo.org/api/records/20745835/files/28C_R1.fastq.gz/content
>    https://zenodo.org/api/records/20745835/files/28C_R2.fastq.gz/content
>    https://zenodo.org/api/records/20745835/files/28A_R1.fastq.gz/content
>    https://zenodo.org/api/records/20745835/files/28B_R2.fastq.gz/content
>    https://zenodo.org/api/records/20745835/files/28A_R2.fastq.gz/content
>    https://zenodo.org/api/records/20745835/files/28B_R1.fastq.gz/content
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets to show the role of each sample in FMT:
>    - `28A_R1` / `28A_R2` → `Pre-FMT_R1` / `Pre-FMT_R2`
>    - `28B_R1` / `28B_R2` → `Donor_R1` / `Donor_R2`
>    - `28C_R1` / `28C_R2` → `Post-FMT_R1` / `Post-FMT_R2`
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
> 4. Check that the datatype of each dataset is set to `fastqsanger.gz`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="fastqsanger.gz" %}
>
> 5. Add a tag to each dataset (`#donor`, `#pre-fmt`, or `#post-fmt`) so they are easy to identify later in the workflow
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

# Running SamestrGal as a Complete Workflow

This section shows how to run the entire SamestrGal analysis in a single step, using the published workflow. The following section then goes through each tool individually to explain what it does.

> <hands-on-title> Run SamestrGal </hands-on-title>
>
> 1. Import the [SamestrGal workflow](https://usegalaxy.eu/published/workflow?id=8e13c87e751e88d8) into your account
>
>    {% snippet faqs/galaxy/workflows_import.md %}
>
> 2. Build a paired collection from the three samples uploaded in the previous step. Each sample's forward and reverse reads will be grouped together
>
>    {% snippet faqs/galaxy/collections_build_list_paired.md %}
>
> 3. Run the **SamestrGal** {% icon workflow %} workflow with the following parameters:
>    - {% icon param-collection %} *"Raw Reads"*: the paired collection created in step 2
>    - *"Select Host reference genome"*: `hg19`
>    - *"Select Sequencer for Trimmomatic"*: `NexteraPE`
>    - *"Run MetaPhlAn"*: `Yes`
>    - *"Select MetaPhlAn Database"*: `mpa_vJan25_CHOCOPhlAnSGB_202503-11062025`
>    - *"Select SameStr database"*: `mpa_vJan25_CHOCOPhlAnSGB_202503-11062025`
>    - In *"Alignment Filtering"*:
>        - *"Percent identity"*: `0.9`
>        - *"Minimum alignment length"*: `40`
>        - *"Minimum base quality"*: `20`
>        - *"Minimum alignment quality"*: `0`
>        - *"Minimum vertical coverage"*: `3`
>    - In *"SameStr Filter settings"*:
>        - *"Minimum samples per clade"*: `2`
>        - *"Marker truncation length"*: `20`
>        - *"Nucleotides minimum variant coverage"*: `2`
>        - *"Minimum variant coverage fraction"*: `0.1`
>        - *"Minimum position coverage"*: `1`
>        - *"Position coverage standard deviation cutoff"*: `3.0`
>        - *"Minimum horizontal coverage"*: `5000`
>        - *"Minimum sample coverage per position"*: `2`
>    - In *"SameStr Summarize settings"*:
>        - *"Minimum overlap for comparison"*: `5000`
>        - *"Minimum similarity for shared strains"*: `0.999`
>
>
>    {% snippet faqs/galaxy/workflows_run.md %}
>
>    > <details-title> Why these parameter values? </details-title>
>    >
>    > These values include some of the parameters mentioned in the original SameStr publication {% cite samestr_podlesny %}. For parameters not specified there, the tool's default values are used instead. A few important parameters are:
>    >
>    > - *"Select SameStr database"* must always match the profiler selected in *"Select MetaPhlAn Database"*. Mixing databases from different profilers will cause errors.
>    > - *"Minimum overlap for comparison"* (5000) and *"Minimum similarity for shared strains"* (0.999) are the two thresholds that ultimately decide whether a clade is classified as a shared strain.
>    >
>    > The remaining alignment and filtering parameters are explained individually, tool by tool, in the next section.
>    {: .details}
>
>
>    > <warning-title> Changing Parameters </warning-title>
>    >
>    > The tutorial dataset was reduced to the minimum size that still produces results with the parameter values listed above. Changing these parameters may therefore cause the workflow to produce no output. To experiment with different parameter values, download the complete samples from the ENA under accession PRJEB39023 and repeat the Data Upload step with the full dataset.
>    {: .warning}
>
{: .hands_on}

# Shared Strain Analysis

This section goes through each tool used by SamestrGal individually, explaining what it does and why it is configured this way. Whether you ran the complete workflow in the previous section or are running the tools separately, start by uploading the data and building a paired collection as described in step 2 of *Running SamestrGal as a Complete Workflow*.

## Pre-processing of reads with **KneadData**

KneadData trims reads based on quality parameters and performs host-read removal on the raw sequencing reads {% cite kneaddata_website %}. This step is performed to ensure that no sequencing-quality errors are introduced into the shared strain analysis. KneadData has several tools internally, which can be seen as different sections in the tool. 

Trimmomatic is the internal tool that performs quality-based trimming and removes adapters. Adapters are short synthetic DNA sequences attached to each read during library preparation so that it can bind to the sequencer. These are not part of the original biological sample, so they must be removed before alignment. Tandem Repeat Finder (TRF) is another of the internal tools, which removes repetitive DNA sequences, which are short sequences repeated consecutively along the genome that can otherwise align to multiple locations and introduce ambiguous mappings. FastQC reports is another of the internal tools that can be enabled to obtain additional outputs of quality metrics of the sequencing reads. Host removal is done because raw metagenomic reads can contain both microbial and human DNA, and only the microbial DNA is relevant for the shared strain analysis. TRF and FastQC are optional and are skipped in this tutorial, matching how KneadData was configured in the original SameStr publication {% cite samestr_podlesny %}.


> <hands-on-title> Trim and Remove Host Reads </hands-on-title>
>
> 1. {% tool [KneadData](toolshed.g2.bx.psu.edu/repos/iuc/kneaddata/kneaddata/0.12.1+galaxy1) %} with the following parameters:
>    - *"Read type"*: `Paired reads`
>        - *"Save all unmatched reads"*: `Yes`
>        - {% icon param-collection %} *"Paired reads collection"*: the paired collection created earlier
>    - In *"Alignment Tool"*:
>        - *"Select alignment tool for host removal"*: `Bowtie2 (default)`
>            - *"Select reference genome database"*: `hg19`
>    - In *"Quality Control with FastQC"*:
>        - *"FastQC quality reports"*: `Do not generate FastQC reports`
>    - In *"Trimmomatic"*:
>        - *"Run Trimmomatic quality/adapter trimming"*: `Enable Trimmomatic trimming`
>            - *"Available sequencers"*: `NexteraPE`
>    - In *"Tandem Repeat Finder (TRF)"*:
>        - *"Specify if you want to include the TRF step in the workflow."*: `Skip TRF Step`
>
> {% snippet faqs/galaxy/tools_select_collection.md %}
>
>    > <warning-title> If TRF is enabled </warning-title>
>    >
>    > Enabling TRF adds two extra output collections containing only the reads affected by repeat removal. These can be useful for checking how many reads remain after TRF, but are not meant to be carried to the next step. `paired_output` and `unmatched_paired` are still the correct outputs to use in the next step.
>    {: .warning}
>
> 2. {% tool [Merge collections](__MERGE_COLLECTION__) %} with the following parameters:
>    - In *"Input Collections"*:
>        - {% icon param-repeat %} *"Insert Input Collections"*
>            - {% icon param-file %} *"Input Collection"*: `paired_output` (output of **KneadData** {% icon tool %})
>        - {% icon param-repeat %} *"Insert Input Collections"*
>            - {% icon param-file %} *"Input Collection"*: `unmatched_paired` (output of **KneadData** {% icon tool %})
>    - In *"Advanced Options"*:
>        - *"How should conflicts (or potential conflicts) be handled?"*: `Keep first instance`
>
>    > <comment-title> Why merge these two collections? </comment-title>
>    >
>    > KneadData outputs two separate collections of reads: `paired_output`, the quality-filtered and host-decontaminated reads that kept their mate pair, and `unmatched_paired`, the orphaned reads. Orphaned reads are those that lost their paired read during quality filtering but still contain information themselves. Merging both collections ensures that all of the surviving microbial sequence data is used in the strain-level analysis.
>    {: .comment}
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Why is `hg19`, a human reference genome, selected as the reference genome database for host removal?
>
> > <solution-title></solution-title>
> >
> > 1. Because the samples used are human stool, so the sequencing reads include human DNA alongside the microbial DNA. This human DNA must be removed before the strain analysis.
> >
> {: .solution}
>
{: .question}

The result of these two steps is a single collection containing one set of reads per sample: all the reads that passed KneadData's quality filtering and host decontamination, both paired and orphaned. The reads are now ready to be used for taxonomic profiling in the next step.


## Taxonomic Profiling and Marker-Based Alignment with **MetaPhlAn**

MetaPhlAn performs both the taxonomic profiling and the marker-based alignment used by SameStr. It maps the pre-processed reads against a set of clade-specific markers derived from reference genomes {% cite blanco-miguez2023metaphlan4 %}. MetaPhlAn produces several outputs. The ones used in the subsequent workflow steps are the taxonomic profile, which lists the relative abundance of each clade detected in the sample, and the reads that aligned to the markers.

> <hands-on-title> Taxonomic profiling and Read Alignment </hands-on-title>
>
> 1. {% tool [MetaPhlAn](toolshed.g2.bx.psu.edu/repos/iuc/metaphlan/metaphlan/4.2.4+galaxy0) %} with the following parameters:
>    - In *"Inputs"*:
>        - *"Input(s)"*: `Fasta/FastQ file(s) with microbiota reads`
>            - *"Fasta/FastQ file(s) with microbiota reads"*: `Paired-end collection`
>                - {% icon param-collection %} *"Paired-end Fasta/FastQ file with microbiota reads"*: the paired collection produced by **Merge collections** {% icon tool %}
>        - *"Database with clade-specific marker genes"*: `Locally cached`
>            - *"Cached database with clade-specific marker genes"*: `mpa_vJan25_CHOCOPhlAnSGB_202503-11062025`
>    - In *"Analysis"*:
>        - *"Type of analysis to perform"*: `rel_ab: Profiling a microbiota in terms of relative abundances`
>            - *"Taxonomic level for the relative abundance output"*: `All taxonomic levels`
>    - *"Subsample"*: `No`
>
>
> 2. {% tool [Samtools view](toolshed.g2.bx.psu.edu/repos/iuc/samtools_view/samtools_view/1.22+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"SAM/BAM/CRAM data set"*: `sam_output_file` (output of **MetaPhlAn** {% icon tool %})
>    - *"What would you like to look at?"*: `All reads in the input dataset`
>        - *"What would you like to have reported?"*: `The actual reads`
>            - *"Output format"*: `BAM (-b)`
>    - *"Use a reference sequence"*: `No`
>
>    > <comment-title> Why convert to BAM? </comment-title>
>    >
>    > MetaPhlAn outputs the aligned reads in SAM format, but SameStr Convert requires BAM input. Samtools performs this conversion.
>    {: .comment}
>
{: .hands_on}

> <details-title> How will taxonomic profiling be used in the workflow? </details-title>
>
> The taxonomic profile identifies which clades are present in each sample and at what relative abundance. SameStr can only compare strains for a clade if that clade is detected in more than one sample, so the taxonomic profile determines, for each pair of samples, which clades are even candidates for shared strain detection in the following steps.
>
{: .details}

## Strain-Level Analysis with **SameStr**

SameStr is not a single tool but a chain of command-line tools, each carrying out one step of the shared-strain detection process. The chain starts with the marker alignments, which are converted into per-clade SNV profiles. These profiles are then merged across samples and filtered to remove low-quality positions, before the samples are compared to identify shared strains between them. The following subsections go through each of these tools in the order they are run.


### **SameStr Convert**


SameStr Convert takes the aligned reads and the taxonomic profile from MetaPhlAn and produces one SNV profile per detected clade, describing the nucleotides observed at each covered position of that clade's marker genes.


> <hands-on-title> Generate Per-Clade SNV Profiles </hands-on-title>
>
> 1. {% tool [SameStr Convert](toolshed.g2.bx.psu.edu/repos/iuc/samestr/samestr_convert/1.2025.111+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Aligned reads"*: the BAM collection produced by **Samtools view** {% icon tool %}
>    - {% icon param-collection %} *"Taxonomic profile"*: the taxonomic profile collection produced by **MetaPhlAn** {% icon tool %}
>    - *"SameStr database"*: `mpa_vJan25_CHOCOPhlAnSGB_202503-11062025`
>    - In *"Alignment Parameters"*:
>        - *"Percent identity"*: `0.9`
>        - *"Minimum alignment length"*: `40`
>        - *"Minimum base quality"*: `20`
>        - *"Minimum alignment quality"*: `0`
>        - *"Minimum vertical coverage"*: `3`
>
>    > <comment-title> What do these alignment parameters control? </comment-title>
>    >
>    > Each of these settings decides whether a given aligned read, or a given position within it, has enough quality to be counted in the SNV profile:
>    >
>    > - *"Percent identity"* is the minimum fraction of bases in a read that must match the reference marker sequence. Reads with more mismatches than this are discarded.
>    > - *"Minimum alignment length"* is the shortest amount of aligned base pairs a read can have and still be kept. Shorter alignments are considered too unreliable to use.
>    > - *"Minimum base quality"* is the lowest sequencing quality score a nucleotide call can have and still be counted. A score of 20 corresponds to a 99% probability that the base was called correctly.
>    > - *"Minimum alignment quality"* is the lowest mapping quality score a read can have and still be kept. Setting it to 0 means no read is excluded on this basis alone.
>    > - *"Minimum vertical coverage"* is the minimum number of reads that must cover a position for that position to be included in the SNV profile.
>    >
>    {: .comment}
>
{: .hands_on}

> <details-title> How is the quality score calculated? </details-title>
>
> Sequencing quality scores follow the Phred scale, where the score Q relates to the probability P of an incorrect base call by Q = -10·log₁₀(P). For a quality score of 20, this gives P = 10⁻² = 0.01, a 1% chance of an incorrect call.
>
{: .details}

> <hands-on-title> Flatten the SNV Profile Collection </hands-on-title>
>
> 1. {% tool [Flatten collection](__FLATTEN__) %} with the following parameters:
>    - {% icon param-collection %} *"Input Collection"*: the per-sample SNV profile collection produced by **SameStr Convert** {% icon tool %}
>    - *"Join collection identifiers using"*: `underscore (_)`
>
>    > <comment-title> Why flatten the collection? </comment-title>
>    >
>    > SameStr Convert organizes its output by sample, with each sample containing one npz file per detected clade. SameStr Merge instead needs a single flat list of npz files, since it groups them by clade internally based on their filenames. Flatten collection removes the per-sample nesting so all npz files can be passed to SameStr Merge together.
>    {: .comment}
>
{: .hands_on}


The Flatten collection tool joins each identifier with an underscore, which adds the sample name as an extra prefix, producing `sample_clade.sample` (e.g. `Pre-FMT_t__SGB14809.Pre-FMT`). This breaks the `clade.sample` format that SameStr Merge expects, so the following step removes this prefix and restores the original identifiers.

> <hands-on-title> Fix the Collection Identifiers </hands-on-title>
>
> 1. {% tool [Apply rules](__APPLY_RULES__) %} with the following parameters:
>    - {% icon param-collection %} *"Input Collection"*: the flattened SNV profile collection produced by **Flatten collection** {% icon tool %}
>    - In the Rule Builder:
>        - *"Add Column"* > *"Using a Regular Expression"*
>        - Select *"Create column from expression replacement"*
>            - *"Regular Expression"*: `^[^_]+_(.*)$`
>            - *"Replacement Expression"*: `\1`
>        - Set *"List Identifier(s)"* to column B
>
{: .hands_on}


> <comment-title> Expected output </comment-title>
>
> If you are following this tutorial with the same parameters and datasets, you should have a single flat list of 28 npz datasets, combining all samples together.
>
{: .comment}

### **SameStr Merge**

SameStr Merge combines the per-sample SNV profiles for the same clade into a single multi-sample profile. For each clade, it produces one merged npz profile combining the data from every sample that contains it, along with a sample name file listing which samples are included.

> <hands-on-title> Merge SNV Profiles Across Samples </hands-on-title>
>
> 1. {% tool [SameStr Merge](toolshed.g2.bx.psu.edu/repos/iuc/samestr/samestr_merge/1.2025.111+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input SNV profiles"*: the SNV profile collection produced by **Apply rules** {% icon tool %}
>    - *"SameStr database"*: `mpa_vJan25_CHOCOPhlAnSGB_202503-11062025`
>
{: .hands_on}


> <question-title></question-title>
>
> 1. The input to SameStr Merge was a list of 28 npz datasets. After running SameStr Merge, this list decreased to 22 datasets. Why do we have less datasets than before?
>
> > <solution-title></solution-title>
> >
> > 1. SameStr Merge combines the per-sample profiles for the same clade into a single profile, so any clade detected in more than one sample is merged into one entry instead of being counted once per sample. 
> >
> {: .solution}
>
{: .question}

### **SameStr Filter**

SameStr Filter removes clades, samples, and positions that do not meet a set of coverage and quality thresholds, before the profiles are used to detect shared strains.

> <hands-on-title> Filter the Merged SNV Profiles </hands-on-title>
>
> 1. {% tool [SameStr Filter](toolshed.g2.bx.psu.edu/repos/iuc/samestr/samestr_filter/1.2025.111+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"SNV profile to filter"*: the merged SNV profile collection produced by **SameStr Merge** {% icon tool %}
>    - {% icon param-collection %} *"Sample identifiers"*: the sample name file collection produced by **SameStr Merge** {% icon tool %}
>    - *"SameStr database"*: `mpa_vJan25_CHOCOPhlAnSGB_202503-11062025`
>    - In *"Clade and marker settings"*:
>        - *"Minimum samples per clade"*: `2`
>        - *"Marker truncation length"*: `20`
>    - In *"Sample Variant Filtering"*:
>        - *"Minimum variant coverage (nucleotides)"*: `2`
>        - *"Minimum variant coverage (fraction)"*: `0.1`
>    - In *"Sample Position Filtering"*:
>        - *"Minimum position coverage (nucleotides)"*: `1`
>        - *"Position coverage standard deviation cutoff"*: `3.0`
>    - In *"Sample Filtering"*:
>        - *"Minimum horizontal coverage (nucleotides)"*: `5000`
>    - In *"Global Position Filtering"*:
>        - *"Minimum sample coverage per position (count)"*: `2`
>
>    > <comment-title> What do these filtering parameters control? </comment-title>
>    >
>    > - *"Minimum samples per clade"* discards clades detected in fewer samples than this.
>    > - *"Marker truncation length"* trims this many nucleotides from each end of every marker gene, removing edge regions where alignments tend to be less reliable.
>    > - *"Minimum variant coverage"* (nucleotides and fraction) sets the minimum read count and minimum proportion of reads that must support a nucleotide call for it to be counted as a real variant, rather than sequencing noise.
>    > - *"Minimum position coverage"* and the *"standard deviation cutoff"* remove positions with too few reads, or with coverage that is too high or too low than a sample's average.
>    > - *"Minimum horizontal coverage"* discards a clade from a sample if too few of its marker positions are covered overall.
>    > - *"Minimum sample coverage per position"* keeps a position in the final profile only if it is covered in at least this many samples.
>    >
>    {: .comment}
>
{: .hands_on}

The resulting collections contain some empty datasets. In order to use the collections in the subsequent steps, the empty datasets need to be eliminated from the collections.

> <hands-on-title> Remove Empty SNV Profile Datasets </hands-on-title>
>
> 1. {% tool [Filter empty datasets](__FILTER_EMPTY_DATASETS__) %} with the following parameters:
>    - {% icon param-collection %} *"Input Collection"*: the filtered SNV profile collection produced by **SameStr Filter** {% icon tool %}
>
{: .hands_on}

> <hands-on-title> Remove Empty Sample Name Datasets </hands-on-title>
>
> 1. {% tool [Filter empty datasets](__FILTER_EMPTY_DATASETS__) %} with the following parameters:
>    - {% icon param-collection %} *"Input Collection"*: the filtered sample name file collection produced by **SameStr Filter** {% icon tool %}
>
{: .hands_on}


> <comment-title> Expected output </comment-title>
>
> If you are following this tutorial with the same parameters and datasets, you should now have 3 datasets in the filtered SNV profile collection, and 3 tabular datasets in the sample name collection.
>
{: .comment}



### **SameStr Stats**

SameStr Stats produces one of the final outputs of the workflow, since its results are not required by any subsequent tool. It calculates statistics on the coverage and nucleotide diversity of the filtered SNV profiles.

> <hands-on-title> Compute Alignment Statistics </hands-on-title>
>
> 1. {% tool [SameStr Stats](toolshed.g2.bx.psu.edu/repos/iuc/samestr/samestr_stats/1.2025.111+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input SNV profiles"*: the SNV profile collection produced by **Filter empty datasets** {% icon tool %}
>    - {% icon param-collection %} *"Sample identifiers"*: the sample name collection produced by **Filter empty datasets** {% icon tool %}
>    - *"SameStr database"*: `mpa_vJan25_CHOCOPhlAnSGB_202503-11062025`
>
{: .hands_on}



> <details-title> What can you learn from these statistics? </details-title>
>
> Two columns are especially useful for a first look at the results. *f_covered* shows what fraction of the clade's marker positions were actually covered in a given sample, which reflects how much of that clade was captured by sequencing. *f_mono* and *f_poly* show whether the covered positions are dominated by a single nucleotide (monomorphic) or contain multiple nucleotides (polymorphic). A high f_mono value suggests the sample is colonized by a single dominant strain for that clade, while a higher f_poly value would point to multiple co-existing strains.
>
{: .details}



### **SameStr Compare**

SameStr Compare calculates pairwise similarity between samples for each clade, using the filtered SNV profiles to determine whether samples share the same strain. SameStr Compare generates 3 output matrices for each clade: an overlap matrix giving the number of shared covered positions between each pair of samples, a fraction matrix giving the Maximum Variant Profile Similarity (MVS) between each pair, and a closest matrix identifying each sample's most similar match. 

> <hands-on-title> Compare Samples for Each Clade </hands-on-title>
>
> 1. {% tool [SameStr Compare](toolshed.g2.bx.psu.edu/repos/iuc/samestr/samestr_compare/1.2025.111+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input SNV profiles"*: the SNV profile collection produced by the *Remove Empty SNV Profile Datasets* step
>    - {% icon param-collection %} *"Sample identifiers"*: the sample name collection produced by the *Remove Empty Sample Name Datasets* step
>    - *"SameStr database"*: `mpa_vJan25_CHOCOPhlAnSGB_202503-11062025`
>
{: .hands_on}


> <details-title> Variant comparison settings </details-title>
>
> This step is run with the default setting, which uses all variants without distinguishing between minority and dominant alleles. Two other options exist: using only dominant variants at each position, or using all variants while also reporting a separate set of metrics based on dominant variants alone. A further option, disabled by default, additionally generates FASTA alignments of the dominant variants for each clade as an extra output.
>
{: .details}


### **SameStr Summarize**

SameStr Summarize is the last tool in the workflow, which produces the shared-strain results. It combines the overlap and similarity matrices from SameStr Compare with the taxonomic profiles from MetaPhlAn, and generates three summary tables: a taxon count table listing the clades detected per sample, a co-occurrence table showing which clades are shared between sample pairs, and a strain events table reporting, for each shared clade, whether the two samples count as sharing the same strain.


> <hands-on-title> Generate the Shared Strain Summary </hands-on-title>
>
> 1. {% tool [SameStr Summarize](toolshed.g2.bx.psu.edu/repos/iuc/samestr/samestr_summarize/1.2025.111+galaxy0) %} with the following parameters:
>    - *"SameStr database"*: `mpa_vJan25_CHOCOPhlAnSGB_202503-11062025`
>    - {% icon param-collection %} *"Overlap metrics (collection output of SameStr compare)"*: the overlap matrix collection produced by **SameStr Compare** {% icon tool %}
>    - {% icon param-collection %} *"Fractional similarities (collection output of SameStr compare)"*: the fraction matrix collection produced by **SameStr Compare** {% icon tool %}
>    - {% icon param-collection %} *"Taxonomic profiles"*: the taxonomic profile collection produced by **MetaPhlAn** {% icon tool %}
>    - In *"Summary threshold arguments"*:
>        - *"Minimum overlap for comparison"*: `5000`
>        - *"Minimum similarity for shared strains"*: `0.999`
>
>    > <comment-title> What do these two thresholds control? </comment-title>
>    >
>    > These are the two values that decide whether a clade counts as a shared strain between two samples or not.
>    >
>    > - *"Minimum overlap for comparison"* (5000) requires at least this many shared covered positions between two samples before their similarity is even considered, so the comparison isn't based on too little data to be reliable.
>    > - *"Minimum similarity for shared strains"* (0.999) is the MVS threshold a sample pair must have to be classified as sharing a strain.
>    >
>    {: .comment}
>
{: .hands_on}


## Interpreting the Results from SameStr Summarize

### Taxon Counts table

![Taxon counts table from SameStr Summarize](../../images/samestr-taxon-counts.png "Number of taxa detected in each sample at each taxonomic level.")

This table shows how many different taxa were detected in each sample at every taxonomic level, from kingdom down to clade. It gives a quick way to compare the overall taxonomic groups found in each sample before looking at any shared-strain results.


### Co-occurrence table

![Co-occurrence table from SameStr Summarize](../../images/samestr-cooccurrences.png "Shared taxa between each pair of samples, from kingdom down to strain level.")

This table is the fastest way to check whether clades are shared at all between samples. *shared_species* and *shared_clade* count how many taxa were detected in both samples, regardless of whether they turned out to be the same strain, while *shared_strain* counts only those clades confirmed as the same strain by SameStr Compare. The gap between *shared_clade* and *shared_strain* tells you how many shared taxa were only confirmed at the species/clade level instead of being detected as strains.

> <question-title></question-title>
>
> 1. Why is it less likely that Pre-FMT and Donor share anything, compared to the other two pairs?
> 2. Can you tell from this table whether any strain persisted in the patient? How would you tell?
>
> > <solution-title></solution-title>
> >
> > 1. The Pre-FMT sample reflects the patient's gut microbiome before treatment, while the Donor reflects a healthy and unrelated microbiome. These two have no direct biological relationship, so any shared strain between them would be coincidental.
> > 2. No. The Pre-FMT/Post-FMT row shows shared_strain = 0, meaning no clade was confirmed as the same strain between these two samples, so this data shows no evidence of persistence.
> >
> {: .solution}
>
{: .question}

### Strain Events table

![Strain events table from SameStr Summarize](../../images/samestr-strain-events.png "Per-clade pairwise comparisons showing overlap, similarity, and the resulting shared-strain classification.")

This table shows the exact clades analyzed between each pair of samples, which is what the co-occurrence table's totals are built from. Each clade is given one of three events: 

- *shared_strain*: when both the overlap and similarity thresholds are met
- *other_strain*: when the overlap threshold is met but similarity falls below 0.999
- *same_clade*: when the overlap threshold itself is not met, so the comparison never reaches strain-level resolution


> <question-title></question-title>
>
> 1. `t__SGB15299` and `t__SGB6362` both have a similarity of 1.0, a perfect match, yet neither is labeled `shared_strain`. Why not?
>
> > <solution-title></solution-title>
> >
> > 1. Their overlap (895 and 2783) is below the minimum overlap threshold of 5000, so the comparison could only be resolved at clade level, not strain level: we can confirm the same clade was found in both samples, but not whether it's the same strain.
> >
> {: .solution}
>
{: .question}


From the table in the image, we can see one shared strain (clade t__SGB4285) between the Donor and Post-FMT, meaning that this strain successfully engrafted.




# Conclusion

This tutorial showed how to run the SamestrGal workflow to detect shared microbial strains between metagenomic samples from an FMT-treated rCDI patient. The tools can also be run separately for pre-processing the reads, aligning them, and finally generating the shared strain results. This tutorial also showed where to find the shared strains between the samples from the workflow outputs, as well as how to distinguish between strain persistence and engraftment based on the results. The analysis found a shared strain between the Donor and the Post-FMT sample, which can be classified as engraftment in the FMT treatment.

