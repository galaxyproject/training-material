---
layout: tutorial_hands_on

title: "Mapping"
zenodo_link: "https://doi.org/10.5281/zenodo.1324070"
questions:
  - What is mapping?
  - What two things are crucial for a correct mapping?
  - What is BAM?
objectives:
  - You will learn what mapping is
  - A genome browser is shown that helps you to understand your data
time_estimation: "1h"
key_points:
  - Know your data!
  - Mapping is not trivial
  - There are many mapping algorithms, it depends on your data which one to choose
requirements:
  -
    type: "internal"
    topic_name: sequence-analysis
    tutorials:
      - quality-control
follow_up_training:
  -
    type: "internal"
    topic_name: transcriptomics
    tutorials:
      - ref-based
  -
    type: "internal"
    topic_name: epigenetics
    tutorials:
      - formation_of_super-structures_on_xi
contributors:
  - joachimwolff
  - bebatut
  - hexylena
---

# Introduction
{:.no_toc}

Sequencing produces a collection of sequences without genomic context. We do not know to which part of the genome the sequences correspond to. Mapping the reads of an experiment to a reference genome is a key step in modern genomic data analysis. With the mapping the reads are assigned to a specific location in the genome and insights like the expression level of genes can be gained.

The short reads do not come with position information, so we do not know what part of the genome they came from. We need to use the sequence of the read itself to find the corresponding region in the reference sequence. But the reference sequence can be quite long (~3 billion bases for human), making it a daunting task to find a matching region. Since our reads are short, there may be several, equally likely places in the reference sequence from which they could have been read. This is especially true for repetitive regions.

In principle, we could do a BLAST analysis to figure out where the sequenced pieces fit best in the known genome. We would need to do that for each of the millions of reads in our sequencing data. Aligning millions of short sequences this way may, however, take a couple of weeks. And we do not care about the exact base to base correspondence (alignment). What we are interested in is "where these reads came from". This approach is called **mapping**.

In the following, we will process a dataset with the mapper **Bowtie2** and we will visualize the data with the program **IGV**.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Prepare the data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial and give it a proper name
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>    {% snippet faqs/galaxy/histories_rename.md %}
>
> 2. Import `wt_H3K4me3_read1.fastq.gz` and `wt_H3K4me3_read2.fastq.gz` from [Zenodo](https://zenodo.org/record/1324070) or from the data library (ask your instructor)
>
>    ```
>    https://zenodo.org/record/1324070/files/wt_H3K4me3_read1.fastq.gz
>    https://zenodo.org/record/1324070/files/wt_H3K4me3_read2.fastq.gz
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
>    As default, Galaxy takes the link as name, so rename them.
>
> 3. Rename the files to `reads_1` and `reads_2`
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
{: .hands_on}

We just imported in Galaxy FASTQ files corresponding to paired-end data as we could get directly from a sequencing facility.

During sequencing, errors are introduced, such as incorrect nucleotides being called. Sequencing errors might bias the analysis and can lead to a misinterpretation of the data. The first step for any type of sequencing data is always to check their quality.

There is a dedicated tutorial for [quality control]({% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %}) of sequencing data. We will not repeat the steps there. You should follow the [tutorial]({% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %}) and apply it to your data before going further.

# Map reads on a reference genome

Read mapping is the process to align the reads on a reference genomes. A mapper takes as input a reference genome and a set of reads. Its aim is to align each read in the set of reads on the reference genome, allowing mismatches, indels and clipping of some short fragments on the two ends of the reads:

![Explanation of mapping](../../images/mapping/mapping.png "Illustration of the mapping process. The input consists of a set of reads and a reference genome. In the middle, it gives the results of mapping: the locations of the reads on the reference genome. The first read is aligned at position 100 and the alignment has two mismatches. The second read is aligned at position 114. It is a local alignment with clippings on the left and right. The third read is aligned at position 123. It consists of a 2-base insertion and a 1-base deletion.")

We need a reference genome to map the reads on.

{% include topics/sequence-analysis/tutorials/mapping/ref_genome_explanation.md answer_3="This data comes from the ChIP-seq of mice, so we will use mm10 (*Mus musculus*)."%}

Currently, there are over 60 different mappers, and their number is growing. In this tutorial, we will use [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/), a fast and memory-efficient open-source tool particularly good at aligning sequencing reads of about 50 up to 1,000s of bases to relatively long genomes.

> ### {% icon hands_on %} Hands-on: Mapping with Bowtie2
> 1. {% tool [Bowtie2](toolshed.g2.bx.psu.edu/repos/devteam/bowtie2/bowtie2/2.4.2+galaxy0) %} with the following parameters
>    - *"Is this single or paired library"*: `Paired-end`
>       - {% icon param-file %} *"FASTA/Q file #1"*: `reads_1`
>       - {% icon param-file %} *"FASTA/Q file #2"*: `reads_2`
>       - *"Do you want to set paired-end options?"*: `No`
>
>           You should have a look at the parameters there, specially the mate orientation if you know it. They can improve the quality of the paired-end mapping.
>
>     - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a built-in genome index`
>       - *"Select reference genome"*: `Mouse (Mus musculus): mm10`
>     - *"Select analysis mode"*: `Default setting only`
>
>         You should have a look at the non default parameters and try to understand them. They can have an impact on the mapping and improving it.
>
>     - *"Save the bowtie2 mapping statistics to the history"*: `Yes`
>
> 2. Inspect the `mapping stats` file by clicking on the {% icon galaxy-eye %} (eye) icon
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What information is provided here?
>    > 2. How many reads have been mapped exactly 1 time?
>    > 3. How many reads have been mapped more than 1 time? How is it possible? What should we do with them?
>    > 4. How many pair of reads have not been mapped? What are the causes?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. The information given here is a quantity one. We can see how many sequences are aligned. It does not tell us something about the quality.
>    > > 2. ~90% reads have been aligned exactly 1 time
>    > > 3. ~7% reads have been aligned concordantly >1 times. These are called multi-mapped reads. It can happen because of repetitions in the reference genome (multiple copies of a gene for example), particularly when the reads are small. It is difficult to decide where these sequences come from and therefore most of the pipelines ignore them. Always check the statistics there to be sure of not discarding too much information in any downstream analyses.
>    > > 4. ~3% pair of reads have not been mapped because
>    > >     - both reads in the pair aligned but their positions do not concord with pair of reads (`aligned discordantly 1 time`)
>    > >     - reads of these pairs are multi-mapped (`aligned >1 times` in `pairs aligned 0 times concordantly or discordantly`)
>    > >     - one read of these pairs are mapped but not the paired read (`aligned exactly 1 time` in `pairs aligned 0 times concordantly or discordantly`)
>    > >     - the rest are not mapped at all
>    >  {: .solution }
>    {: .question}
>
{: .hands_on}

Checking the mapping statistics is an important step to do before continuing any analyses. There are several potential sources for errors in mapping, including (but not limited to):

- **Polymerase Chain Reaction (PCR) artifacts**: Many high-throughput sequencing (HTS) methods involve one or multiple PCR steps. PCR errors will show as mismatches in the alignment, and especially errors in early PCR rounds will show up in multiple reads, falsely suggesting genetic variation in the sample. A related error would be PCR duplicates, where the same read pair occurs multiple times, skewing coverage calculations in the alignment.
- **Sequencing errors**: The sequencing machine can make an erroneous call either for physical reasons (e.g. oil on an Illumina slide) or due to properties of the sequenced DNA (e.g., homopolymers). As sequencing errors are often random, they can be filtered out as singleton reads during variant calling.
- **Mapping errors**: The mapping algorithm can map a read to the wrong location in the reference. This often happens around repeats or other low-complexity regions.

So if the mapping statistics are not good, you should investigate the cause of these errors before going further in your analyses.

After that, you should have a look at the reads and inspect the BAM file where the read mappings are stored.

# Inspection of a BAM file

{% include topics/sequence-analysis/tutorials/mapping/bam_explanation.md mapper="Bowtie2" %}

The BAM file includes a lot of information about each read, particularly the quality of mapping.

> ### {% icon hands_on %} Hands-on: Summary of mapping quality
> 1. {% tool [Samtools Stats](toolshed.g2.bx.psu.edu/repos/devteam/samtools_stats/samtools_stats/2.0.2+galaxy2) %} with the following parameters
>    - {% icon param-file %} *"BAM file"*: `aligned reads` (output of **Bowtie2** {% icon tool %})
>    - *"Use reference sequence"*: `Locally cached/Use a built-in genome`
>      - *"Using genome"*: `Mouse (Mus musculus): mm10 Full`
>
> 2. Inspect the {% icon param-file %} `Stats` file
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What is the proportion of mismatches in the mapped reads when aligned to the reference genome?
>    > 2. What does the error rate represent?
>    > 3. What is the average quality? How is it represented?
>    > 4. What is the insert size average?
>    > 5. How many reads have a mapping quality score below 20?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. There are ~21,900 mismatches for ~4,753,900 bases mapped which on average produces ~0.005 mismatches per mapped bases.
>    > > 2. The error rate is the proportion of mismatches per mapped bases, so the ratio computed right before.
>    > > 3. The average quality is the mean quality score of the mapping. It is a Phred score like the one used in the FASTQ file for each nucleotide. But here the score is not per nucleotide, but per read and it represents the probability of mapping quality.
>    > > 4. The insert size is the distance between the two reads in the pairs.
>    > > 5. To get the info:
>    > >      1. **Filter BAM datasets on a variety of attributes** {% icon tool %} with a filter to keep only the reads with a mapping quality >= 20
>    > >      2. **Stats generate statistics for BAM dataset** {% icon tool %} on the output of **Filter**
>    > >
>    > >    Before filtering: 95,412 reads and after filtering: 89,664 reads.
>    >  {: .solution }
>    {: .question}
{: .hands_on}

# Visualization using a Genome Browser (IGV)

The Integrative Genomics Viewer (IGV) is a high-performance visualization tool for interactive exploration of large, integrated genomic datasets. It supports a wide variety of data types, including array-based and next-generation sequence data, and genomic annotations. In the following, we will use it to visualize the mapped reads.

{% include topics/sequence-analysis/tutorials/mapping/igv.md tool="Bowtie2" region_to_zoom="chr2:98,666,236-98,667,473" %}

# Visualization using a Genome Browser (JBrowse)

{% tool [JBrowse](toolshed.g2.bx.psu.edu/repos/iuc/jbrowse/jbrowse/1.16.11+galaxy0) %} is an alternative, web-based genome browser. Whereas IGV is a piece of software you must download and run, JBrowse instances are websites hosted online that provide an interface to browse genomics data. We'll use it to visualise the mapped reads.

{% include topics/sequence-analysis/tutorials/mapping/jbrowse.md tool="Bowtie2" region_to_zoom="chr2:98,666,236-98,667,473" %}


# Conclusion
{:.no_toc}

After quality control, mapping is an important step of most analyses of sequencing data (RNA-Seq, ChIP-Seq, etc) to determine where in the genome our reads originated from and use this information for downstream analyses.
