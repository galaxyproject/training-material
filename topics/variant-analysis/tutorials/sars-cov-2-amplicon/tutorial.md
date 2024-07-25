---
layout: tutorial_hands_on

title: SARS-CoV-2 ARTIC sequence analysis
zenodo_link: 'https://zenodo.org/record/4944064'
questions:
- What is genomic epidemiology and what can we discover from SARS-CoV-2 whole genome sequence (WGS) data?
- What do we mean by genetic variation, and how do we find variations in SARS-CoV-2 genomic data?
- What are SARS-CoV-2 clades and lineages, and how can we assign SARS-CoV-2 genomic samples to clades and lineage?
objectives:
- How to apply bioinformatics tools to identify variants from SARS-CoV-2 WGS data
- Know what Nextstrain clades and PANGO lineages are in the context of SARS-CoV-2
- Apply the nextclade and pangolin tools to assign SARS-CoV-2 lineages and clades
- Understand how to identify the clades and lineages that SARS-CoV-2 WGS samples belong to
- What quality control steps to apply when analysing SARS-CoV-2 WGS data
keywords:
- SARS-CoV-2, genome, variants
time_estimation: 3H
requirements:
- 
  type: "internal"
  topic_name: "introduction"
  tutorials:
  - galaxy-intro-101
-
  type: "internal"
  topic_name: sequence-analysis
  tutorials:
  - quality-control
  - mapping
-
  type: "internal"
  topic_name: galaxy-interface
  tutorials:
  - get-data
  - collections
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- pvanheus
- mudiboevans
- abright087 

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/introduction.md %}


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Get data

> ### {% icon hands_on %} Hands-on: Data upload 
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>https://zenodo.org/record/4944064/files/ERR4970105_1.fastq.gz
>https://zenodo.org/record/4944064/files/ERR4970105_2.fastq.gz
>https://zenodo.org/record/4944064/files/ERR4970106_1.fastq.gz
>https://zenodo.org/record/4944064/files/ERR4970106_2.fastq.gz
>https://zenodo.org/record/4944064/files/ERR4970107_1.fastq.gz
>https://zenodo.org/record/4944064/files/ERR4970107_2.fastq.gz
>https://zenodo.org/record/4944064/files/MN908947_3_Wuhan-Hu-1.fasta
>https://zenodo.org/record/4944064/files/SARS-CoV-2-ARTICv3.bed
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets (to get rid of the `http://zenodo.org...` prefix)
>
> 4. Check that the datatype
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to the sample ID
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
> 6. Convert the FASTQ datasets to a [List Paired collection]({{ site.baseurl }}/topics/galaxy-interface/tutorials/collections/tutorial.html).
>
{: .hands_on}

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/rbu.md %}

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/fastp.md %}

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/bwa_mem.md %}

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/samtools_stats.md %}

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/bam_filter.md %}

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/bamqc.md %}

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/ivar_trim.md %}

After mapping the reads, quality filtering the mapped reads and trimming primers we have now produce
the final dataset of aligned reads for each of our samples. The next steps that we are going to take
follow two paths: we need to identify the variation in the genome that is present in each sample,
and we need to construct a consensus genome from which we can infer the viral lineage of the sample (and make the genome 
available for deposit in a database or use in constructing a phylogeny).

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/variant_annotation.md %}

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/consensus_genome_characterization.md %}

{% include topics/variant-analysis/tutorials/sars-cov-2-amplicon/lineage_designation.md %}


# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.