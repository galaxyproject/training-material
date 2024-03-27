---
layout: tutorial_hands_on

title: "Pathway analysis with MINERVA"
subtopic: visualisation
priority: 2

tags:
    - bulk
    - rna-seq
    - viz
level: Intermediate
zenodo_link: https://zenodo.org/records/10405036
questions:
    - TODO
objectives:
    - Perform an analysis using a workflow from WorkflowHub
    - Visualise and interpret the results with MINERVA
time_estimation: 1h
key_points:
    - TODO

contributions:
  authorship:
    - mjostaszewski
    - MattiHoch
    - PapXis
    - mbaardwijk
    - shiltemann
    - hexylena
  testing:
    - hexylena
  # Helped manage the plugin, deployment of EU where we tested this
  infrastructure:
    - mira-miracoli
    - kysrpex
    - sanjaysrikakulam
    - bgruening
  funding:
    - by-covid
---

This tutorial is a partial reproduction of {% cite Togami_2022 %} wherein they evaluated mRNA and miRNA in a selection of COVID-19 patients and healthy controls.
While that paper uses a closed source pipeline, we'll be reproducing the analysis with open source tools in Galaxy, using a workflow on WorkflowHub developed for the BY-COVID project.

> <comment-title>Full data</comment-title>
> The original data is available at NCBI under [BioProject PRJNA754796](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA754796)
{: .comment}
>
> <agenda-title></agenda-title>
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

There are several places you can jump to in this tutorial, using pre-calculated
data. We recommend you jump directly to the [Analysis](#Analysis) section, as
that precludes the slowest and most data intensive parts of this tutorial.
However, the entire process is documented in case you want to reproduce our
work.

## Study Design

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Create a new file with the following factor data:
>
>    ```text
>    Run	Group
>    SRR16681566	healthy
>    SRR16681565	healthy
>    SRR16681564	healthy
>    SRR16681563	healthy
>    SRR16681562	healthy
>    SRR16681561	healthy
>    SRR16681560	healthy
>    SRR16681559	healthy
>    SRR16681558	healthy
>    SRR16681557	healthy
>    SRR16681556	healthy
>    SRR16681555	healthy
>    SRR16681554	healthy
>    SRR16681553	healthy
>    SRR16681552	healthy
>    SRR16681551	healthy
>    SRR16681550	COVID
>    SRR16681549	COVID
>    SRR16681548	COVID
>    SRR16681547	COVID
>    SRR16681546	COVID
>    SRR16681545	COVID
>    SRR16681544	COVID
>    SRR16681543	COVID
>    SRR16681542	COVID
>    SRR16681541	COVID
>    SRR16681540	COVID
>    SRR16681539	COVID
>    SRR16681538	COVID
>    SRR16681537	COVID
>    SRR16681536	COVID
>    SRR16681535	COVID
>    SRR16681534	COVID
>    SRR16681533	COVID
>    SRR16681532	COVID
>    SRR16681531	COVID
>    SRR16681530	COVID
>    SRR16681529	COVID
>    SRR16681528	COVID
>    SRR16681527	COVID
>    SRR16681526	COVID
>    SRR16681525	COVID
>    SRR16681524	COVID
>    SRR16681523	COVID
>    SRR16681522	COVID
>    SRR16681521	COVID
>    SRR16681520	COVID
>    SRR15462530	COVID
>    SRR15462529	COVID
>    SRR15462528	COVID
>    SRR15462527	COVID
>    SRR15462526	COVID
>    SRR15462525	COVID
>    SRR15462524	COVID
>    SRR15462523	COVID
>    SRR15462522	COVID
>    SRR15462521	COVID
>    SRR15462520	healthy
>    SRR15462519	healthy
>    SRR15462518	healthy
>    SRR15462517	healthy
>    SRR15462516	healthy
>    ```
>
> - Cut the data
> - Select the SRR lines
> - FasterQ to download
>
{: .hands_on}

## Analysis

We have split this workflow into two parts, based only on how long the first portion of the workflow takes to execute. The rough runtime of the workflow portions when this was being developed can be broken down as follows:

Step                     | Time
---                      | ---
Data Download            | ~6h
Processing Counts        | ~8h
Analysis & Visualisation | 15m

These numbers were generated on UseGalaxy.eu and may not represent the most
efficient possible computation, as they are executed on a shared cluster that can, at times, be more or less busy.

As such we recommend you skip to [Limma](#limma) to progress to the efficient
portion. The data provided in the Zenodo record is from the entire analysis,
analysed with the Counts step that can be skipped:

### Data Download

We'll start by downloading our fastq files from the [GEO Dataset GSE182152](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182152)

> <hands-on-title>Download the data from GEO (ETA: 6 Hours)</hands-on-title>
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1`
>    - {% icon param-file %} *"Select lines from"*: `factordata`
>
> 1. {% tool [Select](Grep1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `out_file1` (output of **Cut** {% icon tool %})
>    - *"that"*: `NOT Matching`
>    - *"the pattern"*: `Run`
>
> 1. {% tool [Faster Download and Extract Reads in FASTQ](toolshed.g2.bx.psu.edu/repos/iuc/sra_tools/fasterq_dump/3.0.8+galaxy1) %} with the following parameters:
>    - *"select input type"*: `List of SRA accession, one per line`
>        - {% icon param-file %} *"sra accession list"*: `out_file1` (output of **Select** {% icon tool %})
>
{: .hands_on}

### Counts


With that done, we can start to analyse the data using HISAT2 and featureCounts

> <hands-on-title>Run the Workflow</hands-on-title>
>
> 1. **Import the workflow** into Galaxy
>
>    {% snippet faqs/galaxy/workflows_run_wfh.md title="mRNA-Seq BY-COVID Pipeline" wfhub_id="688" box_type="none" %}
>
{: .hands_on}

This workflow produces a handful of outputs: the featureCounts results, and a
MultiQC report. Looking at the report we see generally reasonable quality data.

### limma

> <hands-on-title>Only If You Skipped Here: Download the Counts Files</hands-on-title>
>
> 1. Open the Rule Builder
>    - *"Upload data as"*: `Collection(s)`
>    - *"Load tabular data from"*: `Pasted Table`
>    - **Paste** the following table:
>
>      ```
>      https://zenodo.org/records/10405036/files/gene_lengths.tabular	gene_lengths	Gene Lengths
>      https://zenodo.org/records/10405036/files/SRR15462516.featureCounts.tabular	SRR15462516	featureCounts
>      https://zenodo.org/records/10405036/files/SRR15462517.featureCounts.tabular	SRR15462517	featureCounts
>      https://zenodo.org/records/10405036/files/SRR15462518.featureCounts.tabular	SRR15462518	featureCounts
>      https://zenodo.org/records/10405036/files/SRR15462519.featureCounts.tabular	SRR15462519	featureCounts
>      https://zenodo.org/records/10405036/files/SRR15462520.featureCounts.tabular	SRR15462520	featureCounts
>      https://zenodo.org/records/10405036/files/SRR15462521.featureCounts.tabular	SRR15462521	featureCounts
>      https://zenodo.org/records/10405036/files/SRR15462522.featureCounts.tabular	SRR15462522	featureCounts
>      https://zenodo.org/records/10405036/files/SRR15462523.featureCounts.tabular	SRR15462523	featureCounts
>      https://zenodo.org/records/10405036/files/SRR15462524.featureCounts.tabular	SRR15462524	featureCounts
>      https://zenodo.org/records/10405036/files/SRR15462525.featureCounts.tabular	SRR15462525	featureCounts
>      https://zenodo.org/records/10405036/files/SRR15462526.featureCounts.tabular	SRR15462526	featureCounts
>      https://zenodo.org/records/10405036/files/SRR15462527.featureCounts.tabular	SRR15462527	featureCounts
>      https://zenodo.org/records/10405036/files/SRR15462528.featureCounts.tabular	SRR15462528	featureCounts
>      https://zenodo.org/records/10405036/files/SRR15462529.featureCounts.tabular	SRR15462529	featureCounts
>      https://zenodo.org/records/10405036/files/SRR15462530.featureCounts.tabular	SRR15462530	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681520.featureCounts.tabular	SRR16681520	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681521.featureCounts.tabular	SRR16681521	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681522.featureCounts.tabular	SRR16681522	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681523.featureCounts.tabular	SRR16681523	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681524.featureCounts.tabular	SRR16681524	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681525.featureCounts.tabular	SRR16681525	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681526.featureCounts.tabular	SRR16681526	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681527.featureCounts.tabular	SRR16681527	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681528.featureCounts.tabular	SRR16681528	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681529.featureCounts.tabular	SRR16681529	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681530.featureCounts.tabular	SRR16681530	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681531.featureCounts.tabular	SRR16681531	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681532.featureCounts.tabular	SRR16681532	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681533.featureCounts.tabular	SRR16681533	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681534.featureCounts.tabular	SRR16681534	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681535.featureCounts.tabular	SRR16681535	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681536.featureCounts.tabular	SRR16681536	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681537.featureCounts.tabular	SRR16681537	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681538.featureCounts.tabular	SRR16681538	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681539.featureCounts.tabular	SRR16681539	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681540.featureCounts.tabular	SRR16681540	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681541.featureCounts.tabular	SRR16681541	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681542.featureCounts.tabular	SRR16681542	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681543.featureCounts.tabular	SRR16681543	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681544.featureCounts.tabular	SRR16681544	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681545.featureCounts.tabular	SRR16681545	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681546.featureCounts.tabular	SRR16681546	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681547.featureCounts.tabular	SRR16681547	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681548.featureCounts.tabular	SRR16681548	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681549.featureCounts.tabular	SRR16681549	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681550.featureCounts.tabular	SRR16681550	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681551.featureCounts.tabular	SRR16681551	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681552.featureCounts.tabular	SRR16681552	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681553.featureCounts.tabular	SRR16681553	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681554.featureCounts.tabular	SRR16681554	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681555.featureCounts.tabular	SRR16681555	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681556.featureCounts.tabular	SRR16681556	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681557.featureCounts.tabular	SRR16681557	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681558.featureCounts.tabular	SRR16681558	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681559.featureCounts.tabular	SRR16681559	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681560.featureCounts.tabular	SRR16681560	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681561.featureCounts.tabular	SRR16681561	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681562.featureCounts.tabular	SRR16681562	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681563.featureCounts.tabular	SRR16681563	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681564.featureCounts.tabular	SRR16681564	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681565.featureCounts.tabular	SRR16681565	featureCounts
>      https://zenodo.org/records/10405036/files/SRR16681566.featureCounts.tabular	SRR16681566	featureCounts
>      ```
>
>    - Click `Build`
>
> 1. From **Rules** menu, select `Add / Modify Column Definitions`
>     - `Add Definition` → `Collection Name` → Select Column `C`
>     - `Add Definition` → `List Identifier(s)` → Select Column `B`
>     - `Add Definition` → `URL` → Column `A`
>
{: .hands_on}


> <hands-on-title>Analyse the Counts</hands-on-title>
>
> 1. Run the workflow with the Factor Data from the first Hands on, and the datasets from the workflow or Zenodo download, depending on your path:
>
>    {% snippet faqs/galaxy/workflows_run_wfh.md title="mRNA-Seq BY-COVID Pipeline" wfhub_id="689" box_type="none" %}
>
{: .hands_on}

You should have a few outputs, namely the `goseq` outputs, and a table ready for visualisation in MINERVA!

## MINERVA

TODO
