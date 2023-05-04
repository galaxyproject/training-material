---
layout: tutorial_hands_on

title: Mutation calling, viral genome reconstruction and lineage/clade assignment from SARS-CoV-2 sequencing data
subtopic: one-health
level: Intermediate
zenodo_link: "https://zenodo.org/record/5036687"
questions:
- How can we extract annotated allelic variants in SARS-Cov-2 sequences in Galaxy?
- Which tools and workflows can we use to identify SARS-CoV-2 lineages in Galaxy?
objectives:
- Repeat SARS-CoV-2 data preparation
- Select and run workflow to extract annotated allelic variants from FASTQ files
- Run workflow to summarize and generate report for previously called allelic variants
- Interpret summaries for annotated allelic variants
- Run workflow to extract consensus sequences
- Select and run tools to assign clades/lineages
time_estimation: 3H
key_points:
- 4 specialized, best-practice variant calling workflows are available for the identification of annotated allelic variants from raw sequencing data depending on the exact type of input
- Data from batches of samples can be processed in parallel using collections
- Annotated allelic variants can be used to build consensus sequences for and assign each sample to known viral clades/lineages
contributors:
- wm75
- bebatut
tags:
- covid19
- virology
---


# Introduction


Effectively monitoring global infectious disease crises, such as the COVID-19 pandemic, requires capacity to generate and analyze large volumes of sequencing data in near real time. These data have proven essential for monitoring the emergence and spread of new variants, and for understanding the evolutionary dynamics of the virus.

Two sequencing platforms (Illumina and Oxford Nanopore) in combination with several established library preparation (Ampliconic and metatranscriptomic) strategies are predominantly used to generate SARS-CoV-2 sequence data. However, data alone do not equal knowledge: they need to be analyzed. The Galaxy community has developed high-quality analysis workflows to support

- sensitive identification of SARS-CoV-2 allelic variants (AVs) starting with allele frequencies as low as 5% from deep sequencing reads
- generation of user-friendly reports for batches of results
- reliable and configurable consensus genome generation from called variants

> <details-title>Further reading</details-title>
> More information about the workflows, including benchmarking, can be found
> - on the Galaxy Covid-19 effort website: [covid19.galaxyproject.org](https://covid19.galaxyproject.org/)
> - as a BioRxiv preprint: [Global platform for SARS-CoV-2 analysis](https://www.biorxiv.org/content/10.1101/2021.03.25.437046v1)
{: .details}

This tutorial will teach you how to obtain, run and combine these workflows appropriately for different types of input data, be it:

- Single-end data derived from Illumina-based RNAseq experiments
- Paired-end data derived from Illumina-based RNAseq experiments
- Paired-end data generated with Illumina-based Ampliconic (ARTIC) protocols, or
- ONT FASTQ files generated with Oxford nanopore (ONT)-based Ampliconic (ARTIC) protocols

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Prepare Galaxy and data

The suggested input for this tutorial is a special batch of data that is of particular interest as it represents a turning point in the COVID-19 pandemic.
It is a subset (16 samples) of the first sequencing data reported from South Africa at the end of November 2021 for the then novel, fast-spreading SARS-CoV-2 variant that would later be named Omicron.
This data has been Illumina paired-end sequenced after amplification with the ARTIC v4 set of tiled-amplicon primers.

Alternatively, you can also follow this tutorial using your own SARS-CoV-2 sequencing data (you need at least two samples) as long as it is of one of the following types:

- Single-end data derived from Illumina-based whole-genome RNAseq experiments
- Paired-end data derived from Illumina-based whole-genome RNAseq experiments
- Paired-end data generated with Illumina-based tiled-amplicon (e.g. ARTIC) protocols
- ONT FASTQ files generated with Oxford nanopore (ONT)-based tiled-amplicon (e.g. ARTIC) protocols

If you are using your own *tiled-amplicon* data, you are also expected to know the primer scheme used at the amplification step.

## Prepare a new Galaxy history

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

If you are going to use your own sequencing data, there are several possibilities to upload the data depending on how many datasets you have and what their origin is:

- You can import data

  - from your local file system,
  - from a given URL or
  - from a shared data library on the Galaxy server you are working on

  In all of these cases you will also have to organize the imported data into a dataset collection.

- Alternatively, if your data is available from the [NCBI's Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra), you can import it with the help of a dedicated tool, which will organize the data into collections for you.

  > <comment-title>Getting data from SRA</comment-title>
  >
  > The simpler [SARS-CoV-2 sequencing data analysis tutorial]({% link topics/variant-analysis/tutorials/sars-cov-2/tutorial.md %}) uses and explains this alternative way of importing.
  >
  {: .comment}

For the suggested batch of early Omicron data we suggest downloading it via URLs from the
[European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home). In case your Galaxy server offers that same data through a shared data library, this represents a faster (data is already on the server) alternative, so you'll find instructions for this scenario as well below.

> <hands-on-title>Import the sequencing data</hands-on-title>
>
> - Option 1: Import from the ENA
>
>   ```
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/002/SRR17054502/SRR17054502_1.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/002/SRR17054502/SRR17054502_2.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/003/SRR17054503/SRR17054503_1.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/003/SRR17054503/SRR17054503_2.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/004/SRR17054504/SRR17054504_1.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/004/SRR17054504/SRR17054504_2.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/005/SRR17054505/SRR17054505_1.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/005/SRR17054505/SRR17054505_2.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/006/SRR17054506/SRR17054506_1.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/006/SRR17054506/SRR17054506_2.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/007/SRR17054507/SRR17054507_1.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/007/SRR17054507/SRR17054507_2.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/008/SRR17054508/SRR17054508_1.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/008/SRR17054508/SRR17054508_2.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/009/SRR17054509/SRR17054509_1.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/009/SRR17054509/SRR17054509_2.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/010/SRR17054510/SRR17054510_1.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/010/SRR17054510/SRR17054510_2.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/033/SRR17051933/SRR17051933_1.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/033/SRR17051933/SRR17051933_2.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/034/SRR17051934/SRR17051934_1.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/034/SRR17051934/SRR17051934_2.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/035/SRR17051935/SRR17051935_1.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/035/SRR17051935/SRR17051935_2.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/036/SRR17051936/SRR17051936_1.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/036/SRR17051936/SRR17051936_2.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/037/SRR17051937/SRR17051937_1.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/037/SRR17051937/SRR17051937_2.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/038/SRR17051938/SRR17051938_1.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/038/SRR17051938/SRR17051938_2.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/039/SRR17051939/SRR17051939_1.fastq.gz
>   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/039/SRR17051939/SRR17051939_2.fastq.gz
>   ```
>
>   1. Copy the links above
>   2. Open the {% tool [Upload](upload1) %} Manager
>   3. In the top row of tabs select **Collection**
>   4. Configure the drop-down select boxes on that tab like this:
>      - *"Collection Type"*: `Pair`
>      - *"File Type"*: `fastqsanger.gz`
>   5. Click on **Paste/Fetch data** and paste the links you copied into the empty text box
>   6. Press **Start**
>   7. Wait for the **Build** button to become enabled, click it, and, in the next dialogue, give a suitable **Name**, like `Sequencing data`, to the new collection.
>   8. Click on **Create collection**
>
> - Option 2: Import from a shared data library
>
>   {% snippet faqs/galaxy/datasets_import_from_data_library.md astype="as a Collection" collection_type="List of Pairs" collection_name="Sequencing data" tohistory="the history you created for this tutorial" path="GTN - Material / Variant analysis / Mutation calling, viral genome reconstruction and lineage/clade assignment from SARS-CoV-2 sequencing data / DOI: 10.5281/zenodo.5036686" %}
>
{: .hands_on}

## Import reference sequence and auxiliary datasets

Besides the sequenced reads data, we need at least two additional datasets for calling variants and annotating them:

- the SARS-CoV-2 reference sequence [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=fasta) to align and compare our sequencing data against

- a tabular dataset defining aliases for viral gene product names, which will let us translate NCBI RefSeq Protein identifiers (used by the SnpEff annotation tool) to the commonly used names of coronavirus proteins and cleavage products.

> <hands-on-title>Get reference sequence and feature mappings</hands-on-title>
>
> 1. Get the SARS-CoV-2 reference sequence
>
>    A convenient public download link for this sequence is best obtained from the ENA again, where the sequence is known under its [INSDC](https://www.insdc.org/) alias [MN908947.3](https://www.ebi.ac.uk/ena/browser/view/MN908947.3):
>    ```
>    https://www.ebi.ac.uk/ena/browser/api/fasta/MN908947.3?download=true
>    ```
>
>    1. {% tool [Upload](upload1) %} the reference to your history via the link above and make sure the dataset format is set to `fasta`.
>
>       {% snippet faqs/galaxy/datasets_import_via_link.md format="fasta" %}
>    2. {% tool [Replace Text in entire line](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_line/1.1.2) in entire line %} to simplify the reference sequence name
>       - {% icon param-file %} *"File to process"*: the uploaded reference sequence from the ENA
>       - In {% icon param-repeat %} *"1. Replacement"*:
>         - *"Find pattern"*: `^>.+`
>         - *"Replace with"*: `>NC_045512.2`
>
>       While the identifiers `MN908947.3` and `NC_045512.2` really refer to the same sequence, one of the tools we are going to use during the analysis (SnpEff) requires the NCBI RefSeq identifier.
>    3. When the Replace Text tool run is finished, **rename** the output dataset to make it clear that this is the SARS-CoV-2 reference dataset to use in the analysis.
>
>       {% snippet faqs/galaxy/datasets_rename.md name="SARS-CoV-2 reference" format="fasta" %}
> 2. {% tool [Upload](upload1) %} the mapping for translation product identifiers in tabular format
>
>    This mapping really consists of just a few lines of text. Each line lists the NCBI Protein Refseq identifier of a SARS-CoV-2 translation product (which the tool SnpEff knows about and will use for annotating mutation effects), followed by a more commonly used name for that product (which we would like to see in final mutation reports). The last line specifies a "mapping" for the **.** annotation, which SnpEff uses for mutations that do not affect any viral open-reading frame. We do not have a better name for it so we specify that we want to retain this annotation unaltered in reports.
>
>    ```
>    YP_009725297.1    leader
>    YP_009725298.1    nsp2
>    YP_009725299.1    nsp3
>    YP_009725300.1    nsp4
>    YP_009725301.1    3Cpro
>    YP_009725302.1    nsp6
>    YP_009725303.1    nsp7
>    YP_009725304.1    nsp8
>    YP_009725305.1    nsp9
>    YP_009725306.1    nsp10
>    YP_009725307.1    RdRp
>    YP_009725308.1    helicase
>    YP_009725309.1    ExoN
>    YP_009725310.1    endoR
>    YP_009725311.1    MethTr
>    YP_009725312.1    nsp11
>    GU280_gp02        S
>    GU280_gp03        orf3a
>    GU280_gp04        E
>    GU280_gp05        M
>    GU280_gp06        orf6
>    GU280_gp07        orf7a
>    GU280_gp08        orf7b
>    GU280_gp09        orf8
>    GU280_gp10        N
>    GU280_gp11        orf10
>    .                 .
>    ```
>
>    Two remarks on this content:
>    - Since the feature mapping dataset is expected to be in tabular format, but the above display uses spaces to separate columns, please make sure, you have **Convert spaces to tabs** checked when creating the dataset from the copied content!
>    - If you prefer other names for certain translation products than the ones defined above (*e.g.* you might be used to call the first peptide *nsp1* instead of *leader*), you are, of course, free to change those names in the pasted content before uploading it to Galaxy. What needs to be kept is only the Refseq identifiers in the first column as these are fixed by the SnpEff tool.
>
>    {% snippet faqs/galaxy/datasets_create_new_file.md name="SARS-CoV-2 feature mapping" format="tabular" convertspaces="true" %}
>
{: .hands_on}

Another two datasets are needed only for the analysis of ampliconic, e.g. ARTIC-amplified, input data:

- a BED file specifying the primers used during amplification and their binding sites on the viral genome
- a custom tabular file describing the amplicon grouping of the primers

> <details-title>Using your own tiled-amplicon data? Provide the correct primer scheme and amplicon info.</details-title>
>
> The instructions below assume that you are going to analyze viral samples amplified using **version 4 of the ARTIC network's SARS-CoV-2 set of primers**, which is the case for the suggested Omicron batch of data.
>
> If you have decided to analyze your own tiled-amplicon sequencing data in this tutorial, and if your samples have been amplified with a **different** set of primers, you are supposed, at the following step, to upload a primer scheme file and corresponding amplicon information that describes this set of primers.
>
> The Galaxy Project maintains a collection of such files for some commonly used sets of primers at [Zenodo](https://doi.org/10.5281/zenodo.4555734).
> If the files describing your set of primers are part of this collection, you can simply upload them using their Zenodo download URLs (analogous to what is shown for the ARTIC v4 primers below).
>
> For sets of primers *not* included in the collection, you will have to create those files yourself (or obtain them from other sources).
> The expected format for the primer scheme file is 6-column BED format, while the amplicon info file is a simple tabular format that lists all primer names (which must be identical to the ones used in the primer scheme file) that contribute to formation of the same amplicon on a single tab-separated line.
> If in doubt, download the ARTIC v4 files through the URLs provided below, and use them as a template for your own custom files.
>
{: .details}

> <hands-on-title>Get primer scheme and amplicon info</hands-on-title>
>
> 2. Get ARTIC v4 primer files
>
>    ```
>    https://zenodo.org/record/5888324/files/ARTIC_nCoV-2019_v4.bed
>    https://zenodo.org/record/5888324/files/ARTIC_amplicon_info_v4.tsv
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
{: .hands_on}

At this point you should have the following items in your history:

1. A collection of sequenced reads for all samples to be analyzed
   - This collection should be organized as a list with one element per sample that you are planning to analyze.

     If you are going to analyze the suggested batch of Omicron data, the list should have 16 elements.

   - Each of the list elements should itself be a paired collection of a *forward* and a *reverse* reads dataset, each in *fastqsanger.gz* format.

2. The SARS-CoV-2 reference as a single *fasta* dataset
3. The SARS-CoV-2 feature mappings as a single *tabular* dataset
4. Only required if you are analyzing *tiled-amplicon* data (which is the case for the suggested batch):
   - a primer scheme as a single *bed* or *bed6* dataset
   - amplicon information as a single *tabular* dataset

If any of these items are still missing, head back to the corresponding section(s) and upload them now.

If any of these items have not been assigned the correct format (expand the view of each dataset to reveal the format Galaxy has recorded for it), please fix them now.

{% snippet faqs/galaxy/datasets_change_datatype.md %}

If everything looks fine, you are ready to start the actual data analysis.

# Analysis

> The workflows default to requiring an AF ≥ 0.05 and AV-supporting reads of ≥ 10 (these and all other parameters can be easily changed by the user). For an AV to be listed in the reports, it must surpass these thresholds in at least one sample of the respective dataset. We estimate that for AV calls with an AF ≥ 0.05, our analyses have a false-positive rate of < 15% for both Illumina RNAseq and Illumina ARTIC data, while the true-positive rate of calling such low-frequency AVs is ~80% and approaches 100% for AVs with an AF ≥ 0.15. This estimate is based on an initial application of the Illumina RNAseq and Illumina ARTIC workflows to two samples for which data of both types had been obtained at the virology department of the University of Freiburg and the assumption that AVs supported by both sets of sequencing data are true AVs. The second threshold of 10 AV-supporting reads is applied to ensure that calculated AFs are sufficiently precise for all AVs.
>
> More details about the workflows can be found on the [Covid-19 project pages](https://galaxyproject.org/projects/covid19/workflows/) of the [Galaxy Community Hub](https://galaxyproject.org/)

## From FASTQ to annotated allelic variants

To identify the SARS-CoV-2 allelic variants (AVs), a first workflow converts the FASTQ files to annotated AVs through a series of steps that include quality control, trimming, mapping, deduplication, AV calling, and filtering.

Four versions of this workflow are available with their tools and parameters optimized for different types of input data as outlined in the following table:

Workflow version 	| Input data | Read aligner | Variant caller
--- | --- | --- | ---
Illumina RNAseq SE | Single-end data derived from RNAseq experiments | **bowtie2** {% cite langmead_2012 %} | **lofreq** {% cite wilm_2012 %}
Illumina RNAseq PE | Paired-end data derived from RNAseq experiments | **bwa-mem** {% cite li_2010 %} | **lofreq** {% cite wilm_2012 %}
Illumina ARTIC | Paired-end data generated with Illumina-based Ampliconic (ARTIC) protocols | **bwa-mem** {% cite li_2010 %} | **lofreq** {% cite wilm_2012 %}
ONT ARTIC | ONT FASTQ files generated with Oxford nanopore (ONT)-based Ampliconic (ARTIC) protocols | **minimap2** {% cite li_2018 %} | **medaka**

> <details-title>About the workflows</details-title>
>
> - The two Illumina RNASeq workflows (Illumina RNAseq SE and Illumina RNAseq PE) perform read mapping with **bwa-mem** and **bowtie2**, respectively, followed by sensitive allelic-variant (AV) calling across a wide range of AFs with **lofreq**.
> - The workflow for Illumina-based ARTIC data (Illumina ARTIC) builds on the RNASeq workflow for paired-end data using the same steps for mapping (**bwa-mem**) and AV calling (**lofreq**), but adds extra logic operators for trimming ARTIC primer sequences off reads with the **ivar** package. In addition, this workflow uses **ivar** also to identify amplicons affected by ARTIC primer-binding site mutations and excludes reads derived from such “tainted” amplicons when calculating alternative allele frequences (AFs) of other AVs.
> - The workflow for ONT-sequenced ARTIC data (ONT ARTIC) is modeled after the alignment/AV-calling steps of the [ARTIC pipeline](https://artic.readthedocs.io/). It performs, essentially, the same steps as that pipeline’s minion command, i.e. read mapping with **minimap2** and AV calling with **medaka**. Like the Illumina ARTIC workflow it uses **ivar** for primer trimming. Since ONT-sequenced reads have a much higher error rate than Illumina-sequenced reads and are therefore plagued more by false-positive AV calls, this workflow makes no attempt to handle amplicons affected by potential primer-binding site mutations.
>
> All four workflows use **SnpEff**, specifically its 4.5covid19 version, for AV annotation.
>
{: .details}

> <hands-on-title>From FASTQ to annotated AVs</hands-on-title>
>
> 1. **Get the workflow** for your data into Galaxy 
>
>    All workflows developed as part of the Galaxy Covid-19 project can be retrieved via the two popular workflow registries [Dockstore](https://dockstore.org/) and [WorkflowHub](https://workflowhub.eu/), and Galaxy makes this process really easy for you.
>
>    {% snippet faqs/galaxy/workflows_import_search.md search_query='organization:"iwc" name:"sars-cov-2"' box_type="none" %}
>
>    *IWC* (the intergalactic workflow commission) is the organization that the Galaxy Covid-19 project uses to publish its workflows, and the name restriction makes sure we are only getting workflows from that organization that deal with SARS-CoV-2 data analysis.
>
>    Depending on your input data you will need to select the appropriate workflow from the list of hits returned by the workflow registry server.
>    This would be:
>    - **sars-cov-2-pe-illumina-artic-variant-calling/COVID-19-PE-ARTIC-ILLUMINA**
>
>      if you are working with the suggested batch of samples for this tutorial, or if your own data is tiled-amplicon data sequenced on the Illumina platform in paired-end mode
>    - **sars-cov-2-ont-artic-variant-calling/COVID-19-ARTIC-ONT**
>
>      if you are working with tiled-amplicon data sequenced on the ONT platform
>    - **sars-cov-2-pe-illumina-wgs-variant-calling/COVID-19-PE-WGS-ILLUMINA**
>
>      if you are working with WGS data obtained on the Illumina platform in paired-end mode
>    - **sars-cov-2-se-illumina-wgs-variant-calling/COVID-19-SE-WGS-ILLUMINA**
>
>      if you are working with WGS data obtained on the Illumina platform in single-end mode
>
> 2. Run **COVID-19: variation analysis on ...** {% icon workflow %} using the following parameters:
>
>    {% snippet faqs/galaxy/workflows_run.md %}
>
>    Note: the {% icon galaxy-gear %} icon next to **Run Workflow** offers the option to *Send results to a new history*.
>    This is very useful if you are planning to analyze the data in your current history in multiple different ways, and you would like to have each analysis end up in its own dedicated history.
>    Here, however, we only want to do one analysis of our batch of data so we are fine with results of the workflow run getting added to the current history.
>
>    - For **Illumina ARTIC PE** workflow (named **COVID-19: variation analysis on ARTIC PE data**),  *to use for example datasets*
>      - {% icon param-file %} *"1: ARTIC primers to amplicon assignments"*: `ARTIC_amplicon_info_v3.tsv` or `ARTIC amplicon info v3`
>      - {% icon param-file %} *"2: ARTIC primer BED"*: `ARTIC_nCoV-2019_v3.bed` or `ARTIC nCoV-2019 v3`
>      - {% icon param-file %} *"3: FASTA sequence of SARS-CoV-2"*: `NC_045512.2_reference.fasta` or `NC_045512.2 reference sequence`
>      - {% icon param-collection %} *"4: Paired Collection (fastqsanger) - A paired collection of fastq datasets to call variant from"*: paired collection created for the input datasets
>
>    - For **Illumina RNAseq PE** workflow (named **COVID-19: variation analysis on WGS PE data**)
>      - {% icon param-collection %} *"1: Paired Collection (fastqsanger)"*: paired collection created for the input datasets
>      - {% icon param-file %} *"2: NC_045512.2 FASTA sequence of SARS-CoV-2"*: `NC_045512.2_reference.fasta` or `NC_045512.2 reference sequence`
>
>    - For **Illumina RNAseq SE** workflow (named **COVID-19: variation analysis on WGS SE data**)
>      - {% icon param-collection %} *"1: Input dataset collection"*: dataset collection created for the input datasets
>      - {% icon param-file %} *"2: NC_045512.2 FASTA sequence of SARS-CoV-2"*: `NC_045512.2_reference.fasta` or `NC_045512.2 reference sequence`
>
>    - For **ONT ARTIC** workflow (named **COVID-19: variation analysis of ARTIC ONT data**)
>      - {% icon param-file %} *"1: ARTIC primer BED"*: `ARTIC_nCoV-2019_v3.bed` or `ARTIC nCoV-2019 v3`
>      - {% icon param-file %} *"2: FASTA sequence of SARS-CoV-2"*: `NC_045512.2_reference.fasta` or `NC_045512.2 reference sequence`
>      - {% icon param-collection %} *"3: Collection of ONT-sequenced reads"*: dataset collection created for the input datasets
{: .hands_on}

The execution of the workflow takes some time. It is possible to launch the next step even if it is not done, as long as all steps are successfully scheduled.

## From annotated AVs per sample to AV summary

Once the jobs of previous workflows are done, we identified AVs for each sample. We can run a "Reporting workflow" on them to generate a final AV summary.

This workflow takes the collection of called (with lofreq) and annotated (with SnpEff) variants (one VCF dataset per input sample) that got generated as one of the outputs of any of the four variation analysis workflows above, and generates two tabular reports and an overview plot summarizing all the variant information for your batch of samples.

> <warning-title>Use the right collection of annotated variants!</warning-title>
> The variation analysis workflow should have generated *two* collections of annotated variants - one called `Final (SnpEff-) annotated variants`, the other one called `Final (SnpEff-) annotated variants with strand-bias soft filter applied`.
> 
> If you have analyzed ampliconic data with any of the **variation analysis of ARTIC** data workflows, then please consider the strand-bias soft-filtered collection experimental and proceed with the `Final (SnpEff-) annotated variants` collection as input to the next workflow.
>
> If you are working with WGS data using either the **variation analysis on WGS PE data** or the **variation analysis on WGS SE data** workflow, then (and only then) you should continue with the `Final (SnpEff-) annotated variants with strand-bias soft filter applied` collection to eliminate some likely false-postive variant calls.
>
{: .warning}

> <hands-on-title>From annotated AVs per sample to AV summary</hands-on-title>
>
> 1. **Get the workflow** into Galaxy
>    
>    - Option 1: Find workflows on the [WorkflowHub](https://workflowhub.eu) and run them directly on [usegalaxy.eu](https://usegalaxy.eu/)
>
>      Please note that this option currently works *only* with usegalaxy.eu!
>
>      - Open the [workflow page on WokflowHub](https://workflowhub.eu/workflows/109)
>      - Click on `Run on usegalaxy.eu` on the top right of the page
>      
>        The browser will open a new tab with Galaxy's workflow invocation interface.
>
>    - Option 2: Import the workflow from [Dockstore](https://dockstore.org/) using Galaxy's workflow search
>
>      {% snippet faqs/galaxy/workflows_import_search.md trs_server="Dockstore" search_query='organization:"iwc-workflows"' workflow_name="sars-cov-2-variation-reporting/COVID-19-VARIATION-REPORTING" %}
>
>    - Option 3: Import the workflow via its github repo link
>
>      - Open the [workflow GitHub repository](https://github.com/iwc-workflows/sars-cov-2-variation-reporting)
>      - Open the `.ga` file
>      - Click on `Raw` on the top right of the file
>      - Save the file or Copy the URL of the file
>      - Import the workflow to Galaxy
>
>        {% snippet faqs/galaxy/workflows_import.md %}
>
> 2. Run **COVID-19: variation analysis reporting** {% icon workflow %} using the following parameters:
>
>    {% snippet faqs/galaxy/workflows_run.md %}
>
>    - *"Send results to a new history"*: `No`
>    - *"1: AF Filter - Allele Frequency Filter"*: `0.05`
> 
>       This number is the minimum allele frequency required for variants to be included in the report.
>
>    - *"2: DP Filer"*: `1`
>
>       The minimum depth of all alignments required at a variant site;
>       the suggested value will, effectively, deactivate filtering on overall DP and will result in the DP_ALT Filter to be used as the only coverage-based filter.
>
>    - *"3: DP_ALT Filter"*: `10`
>
>       The minimum depth of alignments at a site that need to support the respective variant allele
>
>    - *"4: Variation data to report"*: `Final (SnpEff-) annotated variants`
>
>       The collection with variation data in VCF format: the output of the previous workflow
>
>
>    - *"4: gene products translations"*: `NC_045512.2_feature_mapping.tsv` or `NC_045512.2 feature mapping`
>
>       The custom tabular file mapping NCBI RefSeq Protein identifiers (as used by snpEff version 4.5covid19) to their commonly used names, part of the auxillary data; the names in the second column of this dataset are the ones that will appear in the reports generated by this workflow.
>
>    - *"5: Number of Clusters"*: `3`
>
>      The variant frequency plot generated by the workflow will separate the samples into this number of clusters.
>
{: .hands_on}

The three key results datasets produced by the Reporting workflow are:

1. **Combined Variant Report by Sample**: This table combines the key statistics for each AV call in each sample. Each line in the dataset represents one AV detected in one specific sample

   Column | Field | Meaning
   --- | --- | ---
   1 | `Sample` | SRA run ID
   2 | `POS` | Position in [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/1798174254)
   3 | `FILTER` | `Filter` field from VCF
   4 | `REF` |  Reference base
   5 | `ALT` | Alternative base
   6 | `DP` | Sequencing depth
   7 | `AF` | Alternative allele frequency
   8 | `SB` | Strand bias P-value from Fisher's exact test calculated by [`lofreq`](https://csb5.github.io/lofreq/)
   9 |`DP4` | Depth for Forward Ref Counts, Reverse Ref Counts, Forward Alt Counts, Reverse Alt Counts
   10 |`IMPACT` | Functional impact (from SNPEff)
   11 |`FUNCLASS` | Funclass for change (from SNPEff)
   12 |`EFFECT` | Effect of change (from SNPEff)
   13 |`GENE` | Gene name
   14 |`CODON` | Codon
   15 |`AA` | Amino acid
   16 |`TRID` | Short name for the gene
   17 |`min(AF)` | Minimum Alternative Allele Freq across all samples containing this change
   18 |`max(AF)` | Maximum Alternative Allele Freq across all samples containing this change
   19 |`countunique(change)` | Number of distinct types of changes at this site across all samples
   20 |`countunique(FUNCLASS)` | Number of distinct FUNCLASS values at this site across all samples
   21 |`change` | Change at this site in this sample

   > <question-title></question-title>
   > 
   > 1. How many AVs are found for all samples?
   > 2. How many AVs are found for the first sample in the document?
   > 3. How many AVs are found for each sample?
   >
   > > <solution-title></solution-title>
   > >
   > > 1. By expanding the dataset in the history, we have the number of lines in the file. 868 lines for the example datasets. The first line is the header of the table. Then 867 AVs.
   > >
   > > 2. We can filter the table to get only the AVs for the first sample {% tool [Filter data on any column using simple expressions](Filter1) %} with the following parameters:
   > >    - {% icon param-file %} *"Filter*": `Combined Variant Report by Sample`
   > >    - *"With following condition*": `c1=='ERR5931005'` (to adapt with the sample name)
   > >    - *"Number of header lines to skip*": `1`
   > >
   > >    We got then only the AVs for the selected sample (48 for ERR5931005).
   > >
   > > 3. To get the number of AVs for each sample, we can run {% tool [Group data](Grouping1) %} with the following parameters:
   > >    - {% icon param-file %} *"Select data"*: `Combined Variant Report by Sample`
   > >    - *"Group by column"*: `Column: 1`
   > >    - In *"Operation"*:
   > >      - In *"1: Operation"*:
   > >        - *"Type"*: `Count`
   > >        - *"On column"*: `Column: 2`
   > >    
   > >    With our example datasets, it seems that samples have between 42 and 56 AVs. 
   > {: .solution}
   {: .question}

2. **Combined Variant Report by Variant**: This table combines the information about each AV *across* samples.

   Column | Field | Meaning
   --- | --- | ---
   1 | `POS` | Position in [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/1798174254)
   2 | `REF` | Reference base
   3 | `ALT` | Alternative base
   4 | `IMPACT` | Functional impact (from SnpEff) 
   5 | `FUNCLASS` | Funclass for change (from SnpEff)
   6 | `EFFECT` | Effect of change (from SnpEff)
   7 | `GENE` | Gene 
   8 | `CODON` | Codon 
   9 | `AA` | Amino acid 
   10 |`TRID` | Short name for the gene (from the feature mapping dataset)
   11 |`countunique(Sample)` | Number of distinct samples containing this change 
   12 |`min(AF)` | Minimum Alternative Allele Freq across all samples containing this change 
   13 |`max(AF)` | Maximum Alternative Allele Freq across all samples containing this change 
   14 |`SAMPLES(above-thresholds)` | List of distinct samples where this change has frequency abobe threshold (5%)
   15 |`SAMPLES(all)` | List of distinct samples containing this change at any frequency (including below threshold) 
   16 |`AFs(all)` | List of all allele frequencies across all samples 
   17 |`change` |  Change 

   > <question-title></question-title>
   > 
   > 1. How many AVs are found?
   > 1. What are the different impacts of the AVs?
   > 2. How many variants are found for each impact?
   > 3. What are the different effects of HIGH impact?
   > 4. Are there any AVs impacting all samples?
   >
   > > <solution-title></solution-title>
   > >
   > > 1. By expanding the dataset in the history, we have the number of lines in the file. 184 lines for the example datasets. The first line is the header of the table. Then 183 AVs.
   > >
   > > 2. The different impacts of the AVs are HIGH, MODERATE and LOW.
   > >
   > > 2. To get the number of AVs for each impact levels, we can run {% tool [Group data](Grouping1) %} with the following parameters:
   > >    - {% icon param-file %} *"Select data"*: `Combined Variant Report by Variant`
   > >    - *"Group by column"*: `Column: 4`
   > >    - In *"Operation"*:
   > >      - In *"1: Operation"*:
   > >        - *"Type"*: `Count`
   > >        - *"On column"*: `Column: 1`
   > >    
   > >    With our example datasets, we find:
   > >    - 11 AVs with no predicted impact
   > >    - 52 LOW AVs
   > >    - 111 MODERATE AVs
   > >    - 9 HIGH AVs
   > >
   > > 3. We can filter the table to get only the AVs with HIGH impact by running {% tool [Filter data on any column using simple expressions](Filter1) %} with the following parameters:
   > >    - {% icon param-file %} *"Filter*": `Combined Variant Report by Variant`
   > >    - *"With following condition*": `c4=='HIGH'`
   > >    - *"Number of header lines to skip*": `1`
   > >
   > >    The different effects for the 9 HIGH AVs are STOP_GAINED and FRAME_SHIFT.
   > >
   > > 4. We can filter the table to get the AVs for which `countunique(Sample)` is equal the number of samples (18 in our example dataset): {% tool [Filter data on any column using simple expressions](Filter1) %} with the following parameters:
   > >    - {% icon param-file %} *"Filter*": `Combined Variant Report by Variant`
   > >    - *"With following condition*": `c11==18` (to adapt to the number of sample)
   > >    - *"Number of header lines to skip*": `1`
   > >
   > >    For our example datasets, 4 AVs are found in all samples
   > {: .solution}
   {: .question}

3. **Variant frequency plot**

   ![Variant frequency plot](../../images/sars-cov-2-variant-discovery/variant-frequency.svg)

   This plot represents AFs (cell color) for the different AVs (columns) and the different samples (rows). The AVs are grouped by genes (different colors on the 1st row). Information about their effect is also represented on the 2nd row. The samples are clustered following the tree displayed on the left.

   In the example datasets, the samples are clustered in 3 clusters (as we defined when running the workflow), that may represent different SARS-CoV-2 lineages as the AVs profiles are different.

## From AVs to consensus sequences

For the variant calls, we can now run a workflow which generates reliable consensus sequences according to transparent criteria that capture at least some of the complexity of variant calling:

- Each consensus sequence is guaranteed to capture all called, filter-passing variants as defined in the VCF of its sample that reach a user-defined consensus allele frequency threshold.
- Filter-failing variants and variants below a second user-defined minimal allele frequency threshold are ignored.
- Genomic positions of filter-passing variants with an allele frequency in between the two thresholds are hard-masked (with N) in the consensus sequence of their sample.
- Genomic positions with a coverage (calculated from the read alignments input) below another user-defined threshold are hard-masked, too, unless they are consensus variant sites.

The workflow takes a collection of VCFs and a collection of the corresponding aligned reads (for the purpose of calculating genome-wide coverage) such as produced by the first workflow we ran.

> <warning-title>Use the right collections as input!</warning-title>
> The variation analysis workflow should have generated *two* collections of annotated variants - one called `Final (SnpEff-) annotated variants`, the other one called `Final (SnpEff-) annotated variants with strand-bias soft filter applied`.
> 
> If you have analyzed ampliconic data with any of the **variation analysis of ARTIC** data workflows, then please consider the strand-bias soft-filtered collection experimental and proceed with the `Final (SnpEff-) annotated variants` collection as input to the next workflow.
>
> If you are working with WGS data using either the **variation analysis on WGS PE data** or the **variation analysis on WGS SE data** workflow, then (and only then) you should continue with the `Final (SnpEff-) annotated variants with strand-bias soft filter applied` collection to eliminate some likely false-postive variant calls.
> 
> As for the collection of aligned reads, you should choose the collection of BAM datasets that was used for variant calling (specifically for the first round of variant calling in the case of Illumina ampliconic data). This collection should be called `Fully processed reads for variant calling (primer-trimmed, realigned reads with added indelquals)`, if it was generated by any of the Illumina variation analysis workflows, or `BamLeftAlign on collection ...`, if it was produced by the ONT variation analysis workflow.
> 
{: .warning}


> <hands-on-title>From AVs to consensus sequences</hands-on-title>
>
> 1. **Get the workflow** into Galaxy
>
>    - Option 1: Find workflows on the [WorkflowHub](https://workflowhub.eu) and run them directly on [usegalaxy.eu](https://usegalaxy.eu/)
>
>      Please note that this option currently works *only* with usegalaxy.eu!
>
>      - Open the [workflow page on WokflowHub](https://workflowhub.eu/workflows/138)
>      - Click on `Run on usegalaxy.eu` on the top right of the page
>      
>        The browser will open a new tab with Galaxy's workflow invocation interface.
>
>    - Option 2: Import the workflow from [Dockstore](https://dockstore.org/) using Galaxy's workflow search
>
>      {% snippet faqs/galaxy/workflows_import_search.md trs_server="Dockstore" search_query='organization:"iwc-workflows"' workflow_name="sars-cov-2-consensus-from-variation/COVID-19-CONSENSUS-CONSTRUCTION" %}
>
>    - Option 2: Import the workflow via its github repo link
>
>      - Open the [workflow GitHub repository](https://github.com/iwc-workflows/sars-cov-2-consensus-from-variation)
>      - Open the `.ga` file
>      - Click on `Raw` on the top right of the file
>      - Save the file or Copy the URL of the file
>      - Import the workflow to Galaxy
>
>        {% snippet faqs/galaxy/workflows_import.md %}
>
> 2. Run **COVID-19: consensus construction** {% icon workflow %} using the following parameters:
>
>    {% snippet faqs/galaxy/workflows_run.md %}
>
>    - *"Send results to a new history"*: `No`
>    - *"1: Variant calls"*: `Final (SnpEff-) annotated variants`
>
>       The collection with variation data in VCF format: the output of the first workflow
>
>    - *"2: min-AF for consensus variants*: `0.7`
>
>       Only variant calls with an AF greater than this value will be considered consensus variants.
>
>    - *"3: min-AF for failed variants"*: `0.25`
>
>       Variant calls with an AF higher than this value, but lower than the AF threshold for consensus variants will be considered questionable and the respective sites be masked (with Ns) in the consensus sequence. Variants with an AF below this threshold will be ignored.
>
>    - *"4: aligned reads data for depth calculation"*: `Fully processed reads for variant calling`
>
>       Collection with fully processed BAMs generated by the first workflow. 
>
>       For ARTIC data, the BAMs should NOT have undergone processing with **ivar removereads**
>
>    - *"5: Depth-threshold for masking"*: `5`
>
>       Sites in the viral genome covered by less than this number of reads detection of variants is considered to become unreliable. Such sites will be masked (with Ns) in the consensus sequence unless there is a consensus variant call at the site.
>
>    - *"6: Reference genome*": `NC_045512.2_reference.fasta` or `NC_045512.2 reference sequence`
>
>       SARS-CoV-2 reference genome, part of the auxillary data.
>
{: .hands_on}

The main outputs of the workflow are:
- A collection of viral consensus sequences.
- A multisample FASTA of all these sequences.

The last one can be used as input for tools like **Pangolin** or **Nextclade**.

## From consensus sequences to clade/lineage assignments

To assign lineages to the different samples from their consensus sequences, two tools are available: **Pangolin** and **Nextclade**.

### With Pangolin

Pangolin (Phylogenetic Assignment of Named Global Outbreak LINeages) can be used to assign a SARS-CoV-2 genome sequence the most likely lineage based on the PANGO nomenclature system.

> <hands-on-title>From consensus sequences to clade assignations using Pangolin</hands-on-title>
>
> 1. {% tool [Pangolin](toolshed.g2.bx.psu.edu/repos/iuc/pangolin/pangolin/3.1.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input FASTA File(s)"*: `Multisample consensus FASTA`
>
> 2. Inspect the generated output
{: .hands_on}

Pangolin generates a table file with taxon name and lineage assigned. Each line corresponds to each sample in the input consensus FASTA file provided. The columns are:

Column | Field | Meaning
--- | --- | ---
1 | `taxon` | The name of an input query sequence, here the sample name.
2 | `lineage` | The most likely lineage assigned to a given sequence based on the inference engine used and the SARS-CoV-2 diversity designated. This assignment may be is sensitive to missing data at key sites. [Lineage Description List](https://cov-lineages.org/lineage_description_list.html)
3 | `conflict` | In the pangoLEARN decision tree model, a given sequence gets assigned to the most likely category based on known diversity. If a sequence can fit into more than one category, the conflict score will be greater than 0 and reflect the number of categories the sequence could fit into. If the conflict score is 0, this means that within the current decision tree there is only one category that the sequence could be assigned to.
4 | `ambiguity_score` | This score is a function of the quantity of missing data in a sequence. It represents the proportion of relevant sites in a sequence which were imputed to the reference values. A score of 1 indicates that no sites were imputed, while a score of 0 indicates that more sites were imputed than were not imputed. This score only includes sites which are used by the decision tree to classify a sequence.
5 | `scorpio_call` | If a query is assigned a constellation by scorpio this call is output in this column. The full set of constellations searched by default can be found at the constellations repository.
6 | `scorpio_support` | The support score is the proportion of defining variants which have the alternative allele in the sequence.
7 | `scorpio_conflict` | The conflict score is the proportion of defining variants which have the reference allele in the sequence. Ambiguous/other non-ref/alt bases at each of the variant positions contribute only to the denominators of these scores.
8 | `version` | A version number that represents both the pango-designation number and the inference engine used to assign the lineage.
9 | `pangolin_version` | The version of pangolin software running.
10 | `pangoLEARN_version` | The dated version of the pangoLEARN model installed.
11 | `pango_version` | The version of pango-designation lineages that this assignment is based on.
12 | `status` | Indicates whether the sequence passed the QC thresholds for minimum length and maximum N content.
13 | `note` | If any conflicts from the decision tree, this field will output the alternative assignments. If the sequence failed QC this field will describe why. If the sequence met the SNP thresholds for scorpio to call a constellation, it’ll describe the exact SNP counts of Alt, Ref and Amb (Alternative, Reference and Ambiguous) alleles for that call.

> <question-title></question-title>
> 
> How many different lineages have been found? How many samples for each lineage?
>
> > <solution-title></solution-title>
> >
> > To summarize the number of lineages and number of samples for  each lineage, we can run {% tool [Group data](Grouping1) %} with the following parameters:
> >    - {% icon param-file %} *"Select data"*: output of **pangolin**
> >    - *"Group by column"*: `Column: 2`
> >    - In *"Operation"*:
> >      - In *"1: Operation"*:
> >        - *"Type"*: `Count`
> >        - *"On column"*: `Column: 1`
> >    
> > For our example datasets, we obtain then:
> > - 13 samples B.1.1.7 / Alpha (B.1.1.7-like)
> > - 3 samples B.1.617.2 / Delta (B.1.617.2-like)
> > - 1 sample B.1.525	
> > - 1 sample P.1
> >
> {: .solution}
{: .question}

### With Nextclade

Nextclade assigns clades, calls mutations and performs sequence quality checks on SARS-CoV-2 genomes.

> <hands-on-title>From consensus sequences to clade assignations using Nextclade</hands-on-title>
>
> 1. {% tool [Nextclade](toolshed.g2.bx.psu.edu/repos/iuc/nextclade/nextclade/0.14.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"SARS-CoV-2 consensus sequences (FASTA)"*: `Multisample consensus FASTA`
>    - {% icon param-check %} *"Output options"*: `Tabular format report`
>
> 2. Inspect the generated output
{: .hands_on}

Column | Field | Meaning
--- | --- | ---
1 | `seqName` | Name of the sequence in the source data, here the sample name
2 | `clade` | The result of the clade assignment of a sequence, as defined by Nextstrain. Currently known clades are depicted in the schema below
3 | `qc.overallScore` | Overall QC score
4 | `qc.overallStatus` | Overall QC status
5 | `totalGaps` | Number of - characters (gaps)
6 | `totalInsertions` | Total length of insertions
7 | `totalMissing` | Number of N characters (missing data) 
8 | `totalMutations` | Number of mutations. Mutations are called relative to the reference sequence Wuhan-Hu-1
9 | `totalNonACGTNs` | Number of non-ACGTN characters
10 | `totalPcrPrimerChanges` | Total number of mutations affecting user-specified PCR primer binding sites
11 | `substitutions` | List of mutations
12 | `deletions` | List of deletions (positions are 1-based)
13 | `insertions` | Insertions relative to the reference Wuhan-Hu-1 (positions are 1-based)
14 | `missing` | Intervals consisting of N characters
15 | `nonACGTNs` | List of positions of non-ACGTN characters (for example ambiguous nucleotide codes)
16 | `pcrPrimerChanges` | Number of user-specified PCR primer binding sites affected by mutations
17 | `aaSubstitutions` | List of aminoacid changes
18 | `totalAminoacidSubstitutions` | Number of aminoacid changes
19 | `aaDeletions` | List of aminoacid deletions
20 | `totalAminoacidDeletions` | Number of aminoacid deletions
21 | `alignmentEnd` | Position of end of alignment
22 | `alignmentScore` | Alignment score
23 | `alignmentStart` | Position of beginning of alignment
24 | `qc.missingData.missingDataThreshold` | Threshold for flagging sequences based on number of sites with Ns
25 | `qc.missingData.score` | Score for missing data
26 | `qc.missingData.status` | Status on missing data
27 | `qc.missingData.totalMissing` | Number of sites with Ns
28 | `qc.mixedSites.mixedSitesThreshold` | Threshold for flagging sequences based on number of mutations relative to the reference sequence
29 | `qc.mixedSites.score` | Score for high divergence
30 | `qc.mixedSites.status` | Status for high divergence 
31 | `qc.mixedSites.totalMixedSites` | Number of sites with mutations
32 | `qc.privateMutations.cutoff` | Threshold for the number of non-ACGTN characters for flagging sequences
33 | `qc.privateMutations.excess` | Number of ambiguous nucleotides above the threshold
34 | `qc.privateMutations.score` | Score for ambiguous nucleotides
35 | `qc.privateMutations.status` | Status for ambiguous nucleotides
36 | `qc.privateMutations.total` | Number of ambiguous nucleotides
37 | `qc.snpClusters.clusteredSNPs` | Clusters with 6 or more differences in 100 bases
38 | `qc.snpClusters.score` | Score for clustered differences
39 | `qc.snpClusters.status` | Status for clustered differences
40 | `qc.snpClusters.totalSNPs` | Number of differences in clusters
41 | `errors` | Other errors (e.g. sequences in which some of the major genes fail to translate because of frame shifting insertions or deletions)

> <question-title></question-title>
> 
> How many different lineages have been found? How many samples for each lineage?
>
> ![Illustration of phylogenetic relationship of clades, as used in Nextclade](../../images/sars-cov-2-variant-discovery/clades.svg "Illustration of phylogenetic relationship of clades, as used in Nextclade (Source: <a href="https://github.com/nextstrain/ncov-clades-schema/#ncov-clade-schema">Nextstrain</a>)")
>
> > <solution-title></solution-title>
> >
> > To summarize the number of lineages and number of samples for  each lineage, we can run {% tool [Group data](Grouping1) %} with the following parameters:
> >    - {% icon param-file %} *"Select data"*: output of **Nextclade**
> >    - *"Group by column"*: `Column: 2`
> >    - In *"Operation"*:
> >      - In *"1: Operation"*:
> >        - *"Type"*: `Count`
> >        - *"On column"*: `Column: 1`
> >    
> > For our example datasets, we obtain then:
> > - 10 samples 20I (Alpha, V1)
> > - 4 samples 20B (ancestor of 20I)
> > - 3 samples 21A (Delta)
> > - 1 sample 21D (Eta)
> >
> {: .solution}
{: .question}

### Comparison between Pangolin and Nextclade clade assignments

We can compare **Pangolin** and **Nextclade** clade assignments by extracting interesting columns and joining them into a single dataset using sample ids.

> <hands-on-title>Comparison clade assignations</hands-on-title>
>
> 1. {% tool [Cut columns from a table](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1,c2`
>    - *"Delimited by"*: `Tab`
>    - {% icon param-file %} *"From"*: output of **Nextclade**
>
> 2. {% tool [Cut columns from a table](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1,c2,c5`
>    - *"Delimited by"*: `Tab`
>    - {% icon param-file %} *"From"*: output of **Pangolin**
>
> 3. {% tool [Join two Datasets](join1) %}
>    - {% icon param-file %} *"Join*": output of first **cut**
>    - *"using column"*: `Column: 1`
>    - {% icon param-file %} *"with*": output of second **cut**
>    - *"and column"*: `Column: 1`
>
> 4. Inspect the generated output
{: .hands_on}

We can see that **Pangolin** and **Nextclade** are globally coherent despite differences in lineage nomenclature.

# Conclusion


In this tutorial, we used a collection of Galaxy workflows for the detection and interpretation of sequence variants in SARS-CoV-2:

![Analysis flow in the tutorial](../../images/sars-cov-2-variant-discovery/schema.png "Analysis flow in the tutorial")

The workflows can be freely used and immediately accessed from the three global Galaxy instances. Each is capable of supporting thousands of users running hundreds of thousands of analyses per month. 

It is also possible to automate the workflow runs using the command line as explained in [a dedicated tutorial]({% link topics/galaxy-interface/tutorials/workflow-automation/tutorial.md %}).
