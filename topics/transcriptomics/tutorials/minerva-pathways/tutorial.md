---
layout: tutorial_hands_on

title: "Downstream analysis with MINERVA"
subtopic: visualisation
priority: 2

tags:
    - bulk
    - rna-seq
    - viz
level: Intermediate
zenodo_link:
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
  funding:
    - by-covid
  infrastructure:
    - mira-miracoli
    - kysrpex
    - sanjaysrikakulam
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

# Data upload

TODO

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


# Analysis

## limma

With the

> <hands-on-title>Run the Workflow</hands-on-title>
> 1. **Import the workflow** into Galaxy
>
>    {% snippet faqs/galaxy/workflows_run_wfh.md title="mRNA-Seq BY-COVID Pipeline" wfhub_id="685" %}
>
>    Provide the factor data, RNA Sequencing data, and annotations files to the appropriate fields.
>
{: .hands_on}

## MINERVA
