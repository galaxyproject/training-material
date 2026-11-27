---
layout: tutorial_hands_on

title: Bacterial hybrid genome assembly for AMR surveillance
zenodo_link: 'https://zenodo.org/records/24393063'
questions:
- How can we run Unicycler with Illumina and Oxford Nanopore sequencing data?
- What should we check when using hybrid assemblies for AMR surveillance?
objectives: 
- Follow the bioinformatics workflow from raw sequencing reads to an AMR profile.
- Improve the quality of sequencing data.
- Understand the outputs of tools to assess the sequencing read quality.
- Learn how to run Unicycler.
- Interpret the outputs of the assembly tool.
- Learn how to run AMR gene detection tools.
time_estimation: 6H
key_points:
- Compare read quality, read length, and sequencing yield before and after preprocessing to decide whether the reads are suitable for assembly.
- Illumina reads provide accurate base calls, and Oxford Nanopore reads can span repeats and help connect regions in the assembly graph.
- Use Kraken2 and Bracken to flag unexpected taxa. Use Bandage, QUAST, and CheckM2 to assess assembly structure, contiguity, completeness, and contamination.
- Compare contig IDs in the StarAMR report. Finding an AMR gene and a plasmid replicon on the same contig supports a plasmid location.
contributions:
  authorship:
    - trongthucnguyen
    - praissyzefi
    - jvmourao
---



Bacterial genome assembly reconstructs the DNA sequences of an isolate from overlapping sequencing reads. A **hybrid assembly** combines short and long reads from the same isolate. Illumina short reads provide accurate base calls, while Oxford Nanopore long reads can span regions that are difficult to reconstruct from short reads alone. In this tutorial, we use **Unicycler** {% cite Wick2017Unicycler %} to assemble the bacterial genome. Unicycler uses the short reads to build an initial assembly graph, then uses the long reads to help resolve connections in that graph.

A repeat is a DNA sequence that occurs in more than one location, sometimes on both the bacterial chromosome and a plasmid. A short read contained entirely within a repeat may match several copies equally well. A long read that spans the repeat and extends into unique sequence on both sides can help connect the correct flanking regions. This can produce a less fragmented assembly. When the reads provide enough support, the chromosome and individual plasmids can be resolved as separate sequences.

For **antimicrobial resistance (AMR)** surveillance, the location of a resistance gene can be as informative as its identity. When the same resistance gene is detected in multiple isolates, its surrounding sequence can help investigators distinguish between the spread of a resistant bacterial lineage and the movement of a resistance-carrying plasmid between lineages. A hybrid assembly can provide evidence that a gene is carried on a plasmid or on the chromosome. Establishing transmission requires comparing multiple isolates and considering the epidemiological evidence.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# General introduction

In this tutorial, we will use paired-end Illumina and Nanopore reads from *Escherichia coli* isolate **GPT24-02** ([BioSample SAMEA117852473](https://www.ncbi.nlm.nih.gov/biosample/49875565)) as an example. We will assess the quality of the raw reads and screen for potential contamination. We will then trim adapters, filter the reads, and build a hybrid assembly with Unicycler.

After assembly, we will inspect the assembly graph, assess assembly quality, and estimate genome completeness and contamination. Finally, we will screen the assembled genome for AMR determinants and examine what the contigs and assembly graph can tell us about where detected resistance genes are located.

# Let's get started

Your first hybrid assembly might feel like a big task, but you’ll build it one step at a time. Take your time and see what you discover as your reads come together into a genome.

## Importing the data into Galaxy

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new, empty analysis history for this tutorial. Name it `Hybrid assembly`.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import both Illumina mates and one ONT FASTQ from [Zenodo]({{ page.zenodo_link }}). Use the full ONT file for the subsampling exercise or the prepared subset to start at NanoPlot.
>
>    ***Illumina paired-end reads***
>
>    ```
>    https://zenodo.org/records/24393063/files/GCF_965641715.1_Illumina_R1.fastq.gz?download=1
>    https://zenodo.org/records/24393063/files/GCF_965641715.1_Illumina_R2.fastq.gz?download=1
>    ```
>    
>    ***Full Oxford Nanopore reads***
>
>    ```
>    https://zenodo.org/records/24393063/files/GCF_965641715.1_ONT.fastq.gz?download=1
>    ```
>
>    ***Pre-subsampled Oxford Nanopore reads***
>
>    ```
>    https://zenodo.org/records/24393063/files/GCF_965641715.1.ONT_subsample.fastq.gz?download=1
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Rename the datasets `Ecoli_Illumina_R1` and `Ecoli_Illumina_R2`. Name the ONT dataset `Ecoli_ONT_full` or `Ecoli_ONT_subset`, depending on the file you imported.
> 4. After the uploads finish, check that Galaxy recognized each dataset as FASTQ.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md %}
>
{: .hands_on}

> <details-title>Extra practice with Enterococcus faecium GPT24-06</details-title>
>
> Try the hybrid assembly workflow on a second isolate, *Enterococcus
> faecium* **GPT24-06**. Create a new Galaxy history named
> `GPT24-06 hybrid assembly`. Import both Illumina mates and one ONT
> FASTQ from [Zenodo]({{ page.zenodo_link }}).
>
> **Illumina paired-end reads**
>
> ```text
> https://zenodo.org/records/24393063/files/GCF_965641755.1_Illumina_R1.fastq.gz?download=1
> https://zenodo.org/records/24393063/files/GCF_965641755.1_Illumina_R2.fastq.gz?download=1
> ```
>
> **Pre-subsampled ONT reads**
>
> Use this **10% subset** to start directly at NanoPlot.
>
> ```text
> https://zenodo.org/records/24393063/files/GCF_965641755.1_ONT_subsample.fastq.gz?download=1
> ```
>
> **Full ONT reads**
>
> Use this file to practise subsampling. In `seqtk_sample`, set the
> fraction to `0.10` and the random seed to `42`, then assess the
> subset with NanoPlot.
>
> ```text
> https://zenodo.org/records/24393063/files/GCF_965641755.1_ONT.fastq.gz?download=1
> ```
>
> {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> After assembly, compare your GPT24-06 results with GPT24-02. What
> differences do you see in the assembly graph, contig count, and
> total assembly length?
>
{: .details}

## Subsample the Oxford Nanopore reads with seqtk_sample

The full *E. coli* GPT24-02 ONT dataset contains ~1.45 Gbp of sequence.
We will randomly retain approximately 20% of the reads to reduce the
amount of data processed during the training. We will use a fixed random
seed so that the same input and settings produce a reproducible subset.

> <hands-on-title>Create a 20% subset of the *E. coli* GPT24-02 ONT reads</hands-on-title>
>
> 1. Run {% tool [seqtk_sample](toolshed.g2.bx.psu.edu/repos/iuc/seqtk/seqtk_sample/1.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input FASTA/Q file"*: `Ecoli_ONT_full`
>    - *"RNG seed"*: `42`
>    - *"Subsample (decimal fraction or number)"*: `0.2`
>
> 2. Rename the output `Ecoli_ONT_subset`.
>
> 3. Use this FASTQ as the input for NanoPlot in the next section.
>
{: .hands_on}

> <tip-title>Using the pre-subsampled Zenodo file?</tip-title>
>
> If you imported `Ecoli_ONT_subset` from Zenodo, continue
> directly to NanoPlot. Sampling that file again at `0.2` would retain
> about 4% of the original reads.
>
{: .tip}

> <question-title>Check the yield after subsampling</question-title>
>
> The full ONT dataset contains 1,447,257,134 bases. If you retain 20% of
> the reads, approximately how many bases do you expect? Compare your
> estimate with the total bases reported by NanoPlot for the subset.
>
> > <solution-title></solution-title>
> >
> > We expect about **289.45 Mbp**. NanoPlot reported **289.06 Mbp** for the `Ecoli_ONT_subset`.
> > The values differ slightly because seqtk_sample selects whole reads
> > and ONT reads vary in length.
> >
> {: .solution}
>
{: .question}

# Overview of read quality control and preprocessing

## Introduction to read quality control

Before performing bacterial hybrid genome assembly, sequencing reads generated using Illumina and Oxford Nanopore Technologies (ONT) should be assessed and, where appropriate, preprocessed.

The two sequencing technologies have different read characteristics:

- **Illumina sequencing** produces relatively short reads, generally with high per-base accuracy. Quality issues may include low-quality bases, adapter contamination, and differences in read quality towards the ends of reads.
- **Oxford Nanopore sequencing** produces long reads that can span repetitive regions and help resolve genome structure. Read lengths and quality scores can vary considerably, and reads may contain adapters, low-quality regions, or chimeric sequences.

In this tutorial, we will use Galaxy to:

1. Assess raw Illumina reads using {% cite FastQC %} and summarize the reports using {% cite MultiQC %}.
2. Assess raw ONT reads using {% cite NanoPlot %}.
3. Preprocess paired-end Illumina reads using {% cite Fastp %}.
4. Preprocess ONT reads using {% cite Fastplong %} followed by {% cite Filtlong %}.
5. Reassess the processed reads and compare the before-and-after QC results.
6. Discuss whether the processed reads are appropriate for bacterial hybrid assembly.

> <details-title> With Illumina paired-end reads </details-title>
>
> We will examine R1 and R2 separately because the two read directions can have different quality profiles.
>
{: .details}

## Quality assessment of raw Illumina reads

### Assess Illumina read quality with **FastQC**

**FastQC** provides a quality-control report for individual FASTQ files. It summarizes several characteristics of the reads, including per-base sequence quality, sequence length distribution, GC content, adapter content, and overrepresented sequences.

> <hands-on-title> Run FastQC </hands-on-title>
>
> 1. Run {% tool [FastQC](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.74+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Raw read data from your current history"*: `Ecoli_Illumina_R1`.
>
> 2. Repeat the step with `Ecoli_Illumina_R2`.
>
> 3. Open the **Webpage** output from each FastQC run and compare the R1 and R2 reports.
>
{: .hands_on}

**Pay particular attention to the following FastQC modules.**

| FastQC module | What it shows | Why it matters |
| --- | --- | --- |
| Per base sequence quality | Distribution of quality scores at each position along the reads | Helps identify positions where quality drops, often toward the ends of reads. |
| Per sequence quality scores | Distribution of mean quality scores across reads | Shows whether a subset of reads has consistently low average quality. |
| Per base sequence content | Proportion of A, T, G, and C at each read position | Positional biases can reflect library preparation or sequencing effects. |
| Per sequence GC content | Distribution of GC content across reads | An unexpected peak or shape can prompt a check for contamination or a mixture of organisms. |
| Sequence length distribution | Distribution of read lengths | Shows whether read lengths match the expected sequencing settings or vary across the dataset. |
| Sequence duplication levels | How often identical sequences occur | High duplication can reflect library composition, abundant sequences, or technical duplication. |
| Overrepresented sequences | Sequences that occur unusually often | Helps flag adapters, primers, or other abundant sequences for investigation. |
| Adapter content | Frequency of detected adapter sequences along the reads | Shows whether adapters are detectable and where adapter read-through tends to appear. |

> <question-title>Compare the two read directions</question-title>
>
> 1. Do R1 and R2 have the same number of sequences?
> 2. If adapter content rises near the 3′ end, what might explain it?
>
> > <solution-title></solution-title>
> >
> > 1. For *E. coli* GPT24-02, expect 998,943 sequences in each file.
> > 2. Some DNA fragments can be shorter than the sequencing read.
> >    When sequencing continues past the end of a fragment, it can
> >    read into the adapter.
> >
> {: .solution}
>
{: .question}

### **MultiQC**: summarize Illumina QC reports

When working with multiple samples, inspecting individual FastQC reports
one by one can be time-consuming. **MultiQC** collects
results from supported tools and presents them in a combined report.

In this exercise, you will use MultiQC to compare the FastQC results for
`Ecoli_Illumina_R1` and `Ecoli_Illumina_R2` side by side. In a larger
surveillance project, the same approach can help you compare read quality
across isolates or sequencing batches.

> <hands-on-title>Combine the FastQC reports for R1 and R2</hands-on-title>
>
> 1. Run {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.9+galaxy1) %}.
>    - Under **Results**, select **FastQC** as the tool that generated the logs.
>    - Select the two FastQC **RawData** outputs from `Ecoli_Illumina_R1`
>      and `Ecoli_Illumina_R2`.
>
> 2. Open the **MultiQC HTML report**. Check that both mates are listed,
>    then compare their **sequence counts**, **per-base sequence quality**,
>    and **adapter content**. Use the report to compare R1 and R2 and identify whether the same patterns occur in both read directions.
>
{: .hands_on}

> <question-title>Compare the two Illumina mates</question-title>
>
> How many sequences does **MultiQC** report for each mate? Which plots would
> you inspect before deciding whether trimming is needed?
>
> > <solution-title></solution-title>
> >
> > Both mates have **998,943 sequences**. Check **Per Base Sequence
> > Quality** for a drop in quality near the read ends and **Adapter
> > Content** for adapter read-through.
> >
> {: .solution}
>
{: .question}

## Quality assessment of raw ONT reads

### **NanoPlot**: assess ONT reads before trimming

ONT sequencing generates reads with a wide range of lengths. Long reads can span repeats and help produce more contiguous bacterial genome assemblies. To understand what the reads can contribute to the assembly, examine both their length and quality.

**NanoPlot** summarizes read lengths and quality scores and shows how these two measurements relate to each other.

> <hands-on-title>Run NanoPlot ONT subset</hands-on-title>
>
> 1. Run {% tool [NanoPlot](toolshed.g2.bx.psu.edu/repos/iuc/nanoplot/nanoplot/1.48.0+galaxy0) %} with the following parameters:
>    - *"Select multifile mode"*: `batch`
>    - *"Type of the file(s) to work on"*: `fastq`
>    - {% icon param-file %} **Input FASTQ**: `Ecoli_ONT_subset`, either the output from **seqtk_sample** or the pre-subsampled file from Zenodo.
>    - In *"Options for filtering or transforming input prior to plotting"*:
>        - *"Logarithmic scaling of lengths in plots."*: `Yes`
>    - Leave NanoPlot's filtering and downsampling options unset.
>
> 2. Open the NanoPlot report. Record the read count, total bases,
>    read N50, and mean read quality. Keep the report for comparison
>    after trimming and filtering.
>
{: .hands_on}

**Interpret the NanoPlot report**

| Metric | What it shows | Why it matters |
| --- | --- | --- |
| Number of reads | Number of reads in the input FASTQ | Shows how many reads are available for the next step. |
| Total bases | Sum of all read lengths | Dividing this value by the estimated genome size gives an estimate of sequencing depth. |
| Read length distribution | How read lengths vary across the dataset | Shows the range of read lengths available to span repeats. |
| Read N50 | Read length at which reads of that length or longer account for at least half of all sequenced bases | Summarizes read length by weighting each read by the number of bases it contributes. |
| Mean read quality | Average of the per-read quality scores | Gives an overview of read quality across the dataset. |
| Quality versus read length | How read quality varies with read length | Helps you assess whether the longest reads also have usable quality. |

> <question-title>Estimate ONT depth before trimming</question-title>
>
> NanoPlot reported **289,062,596 bases** for the GPT24-02 subset.
> Assuming a **5.0 Mbp** genome, what is the estimated ONT depth?
>
> > <solution-title></solution-title>
> >
> > Divide **289,062,596 sequenced bases** by **5,000,000 bases**.
> > The estimated depth is **57.8×** before trimming and filtering.
> > Repeat the calculation with the total bases from the post-filtering
> > NanoPlot report to estimate the depth going into assembly.
> >
> {: .solution}
>
{: .question}

### Preprocess the paired-end Illumina reads with **fastp**

After assessing the raw Illumina data, we can perform preprocessing. **fastp** is a tool for quality control and preprocessing of FASTQ data that can detect and trim adapters, trim low-quality bases, filter reads according to selected criteria, and produce reports describing the processing.

> <details-title>How much trimming is appropriate?</details-title>
>
> Avoid applying aggressive trimming without a reason: excessive trimming
> can reduce read length and the amount of data available for assembly.
> The objective is not to remove as many reads as possible. It is to improve
> the usable read dataset while retaining sufficient high-quality data for
> assembly.
>
{: .details}

> <hands-on-title> Create a paired collection and run fastp </hands-on-title>
>
> 1. Select `Ecoli_Illumina_R1` and `Ecoli_Illumina_R2` in your Galaxy history.
>    Create a **List of Paired Datasets** containing one pair. Check that
>    `Ecoli_Illumina_R1` is **forward** and `Ecoli_Illumina_R2` is **reverse**.
>    Name the collection `Ecoli_Illumina_paired`.
>
>    {% snippet faqs/galaxy/collections_build_list_paired.md %}
>
> 2. Run {% tool [fastp](toolshed.g2.bx.psu.edu/repos/iuc/fastp/fastp/1.4.0+galaxy0) %} with the following parameters:
>    - *"Single-end or paired reads"*: `Paired Collection`
>        - {% icon param-collection %} *"Select paired collection(s)"*: `Ecoli_Illumina_paired`
>    - In *"Duplicated Reads Options"*:
>        - *"Enable duplicated reads analysis"*: `Disable (--dont_eval_duplication)`
>    - Leave the other settings at their defaults. Keep read merging off
>      so the output contains separate forward and reverse FASTQ files.
>
> 3. Open the fastp HTML report. Under **Summary**, compare **total reads**
>    and **total bases** before and after filtering.
>
> 4. Rename the cleaned paired collection `Ecoli_Illumina_cleaned` and keep it for the Unicycler step.
>
{: .hands_on}

> <question-title> Interpret the fastp report </question-title>
>
> 1. Each FastQC report lists 998,943 sequences. Why does fastp report
>    1,997,886 total reads? How many read pairs does that represent?
> 2. How many reads and bases did fastp remove? What percentage of bases
>    reached Q30 after filtering?
>
> > <solution-title></solution-title>
> >
> > 1. FastQC reports R1 and R2 separately. fastp counts both mates together,
> >    so 998,943 R1 reads plus 998,943 R2 reads give 1,997,886 reads,
> >    representing **998,943 read pairs**.
> > 2. fastp removed **zero reads and zero bases** in this run. The before
> >    and after totals were identical at 1,997,886 reads and 295,680,581
> >    bases. After filtering, **83.45% of bases** had a quality score
> >    of Q30 or higher, the same percentage as before filtering.
> >
> {: .solution}
>
{: .question}

### **FastQC** and **MultiQC**: assess processed Illumina reads

After preprocessing, run FastQC on the processed R1 and R2 files from
fastp. Then use MultiQC to summarize the two reports and compare them
with the results from the raw reads.

> <hands-on-title> Check Illumina read quality after fastp </hands-on-title>
>
> 1. Run {% tool [FastQC](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.74+galaxy1) %} on the cleaned **R1** FASTQ from **fastp**. Leave the other settings at their defaults.
> 2. Repeat for the cleaned **R2** FASTQ. If the files are inside the
>    cleaned paired collection `Ecoli_Illumina_cleaned`, expand the collection to find each mate.
> 3. Run **MultiQC** and select **FastQC** as the tool that generated the
>    logs. Use the two **RawData** outputs from the new FastQC runs.
> 4. Open the MultiQC HTML report. Compare processed R1 with raw R1, and
>    processed R2 with raw R2. Check **sequence counts**,
>    **per-base sequence quality**, and **adapter content**.
>
{: .hands_on}

### Preprocessing of ONT reads using **fastplong** and **Filtlong**

ONT reads are generally longer than Illumina reads, and their length and quality can vary substantially. We use two tools designed for long-read data:

1. **fastplong** performs long-read FASTQ preprocessing, including supported adapter handling and quality/length filtering.
2. **Filtlong** applies additional read selection based on read quality and length criteria.

### **fastplong**: preprocess ONT reads

**fastplong** is designed for preprocessing long-read FASTQ data. Depending on the selected options, it can handle adapters, filter reads, and process chimeric reads. For this exercise, we will use the raw ONT FASTQ file `Ecoli_ONT_subset` as input.

> <hands-on-title> Run fastplong </hands-on-title>
>
> 1. Run {% tool [fastplong](toolshed.g2.bx.psu.edu/repos/iuc/fastplong/fastplong/0.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input Long Reads"*: `Ecoli_ONT_subset`
>    - *"Disable adapter trimming"*: `No`
>    - *"Disable quality filtering"*: `No`
>    - *"Disable length filtering"*: `No`
>    - *"Sampling / selection mode"*: `No sampling - keep all reads passing the filters`
>
> 2. Open the fastplong HTML report. Record the number of reads and total
>    bases before and after preprocessing. Check the adapter section
>    for evidence of trimming.
>
{: .hands_on}


### **Filtlong**: select ONT reads for downstream assembly

After **fastplong**, we use **Filtlong** to apply additional read selection. **Filtlong** can prioritize reads based on quality and length and can be configured using minimum thresholds or a target amount of data.

For this exercise, we will apply a minimum read length of 2,000 bases. Choose a threshold based on read quality and the goal of the assembly. This threshold is not a universal requirement for bacterial genome assembly.

> <hands-on-title> Run filtlong </hands-on-title>
>
> 1. Run {% tool [Filtlong](toolshed.g2.bx.psu.edu/repos/iuc/filtlong/filtlong/0.3.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} **Input FASTQ**: the processed ONT FASTQ from **fastplong**
>    - *"Min. length"*: `2000`
>    - Leave **Keep percentage**, **Total bases**, and the quality
>      thresholds unset.
>
{: .hands_on}

### **NanoPlot**: assess processed ONT reads

> <hands-on-title>Assess the ONT reads after filtering</hands-on-title>
>
> 1. Run {% tool [NanoPlot](toolshed.g2.bx.psu.edu/repos/iuc/nanoplot/nanoplot/1.48.0+galaxy0) %} on the **filtered FASTQ** from **Filtlong**.
>    - *"Select multifile mode"*: `batch`
>    - *"Type of the file(s) to work on"*: `fastq`
>    - *"Logarithmic scaling of lengths in plots."*: `Yes`
>    - Leave NanoPlot's filtering and downsampling options unset.
>
> 2. Open the NanoPlot report and record the **read count**, **total
>    bases**, **read N50**, and **mean read quality**.
>
{: .hands_on}

### Taxonomic screening for potential contamination using **Kraken2** and **Bracken**

We expect the reads for *E. coli* GPT24-02 to come mainly from *Escherichia coli*.
**Kraken2** {% cite Wood2019Kraken2 %} assigns reads to taxa using a reference database. We will use
its report to check whether the expected species dominates and whether
other assignments warrant a closer look. **Bracken** {% cite Lu2017Bracken %} uses the **Kraken2**
report to estimate species abundance for the Illumina reads.

> <hands-on-title> Classify the raw Illumina reads with Kraken2 </hands-on-title>
>
> 1. Run {% tool [Kraken2](toolshed.g2.bx.psu.edu/repos/iuc/kraken2/kraken2/2.17.2+galaxy0) %} with the following settings:
>    - *"Single or paired reads"*: choose the paired option and select the raw **Ecoli_Illumina_paired**.
>    - *"Select a Kraken2 database"*: `PlusPF-16`
>    - *"Print a report with aggregate counts/clade to file"*: `Yes`.
>    - Leave the Kraken 1 MPA format and minimizer data options off.
>
> 2. Open the **Report** output and find the *Escherichia coli* row.
>    Record the most abundant assignments and the unclassified fraction.
>
{: .hands_on}

> <hands-on-title> Estimate species abundance with Bracken </hands-on-title>
>
> 1. Run {% tool [Bracken](toolshed.g2.bx.psu.edu/repos/iuc/bracken/est_abundance/3.1+galaxy0) %} with the following settings:
>    - *"Kraken report file"*: the **Report** output from Kraken2.
>    - *"Select a Kmer distribution"*: `PlusPF-16` from the same database build used in Kraken2. If Galaxy lists multiple read lengths for that build, select the `150 bp` distribution.
>    - *"Level"*: `Species`.
>
> 2. Open the Bracken abundance table. Identify the species with the
>    largest estimated read count and note any other species reported.
>
{: .hands_on}

> <question-title>What do Kraken2 and Bracken tell us?</question-title>
>
> 1. Why can **Kraken2** assign a read to a genus when you expect it to
>    come from a particular species?
> 2. What does **Bracken** do with the **Kraken2** report, and why should it
>    use a distribution built from the same database?
>
> > <solution-title></solution-title>
> >
> > 1. Related species can share sequences. When the evidence in a
> >    read supports several species within the same genus, **Kraken2**
> >    can assign the read to the genus.
> >
> > 2. **Bracken** uses the counts in the **Kraken2** report and a read-length
> >    distribution model to estimate abundance at the selected rank.
> >    The model must match the Kraken2 database because it describes
> >    how reads are expected to be classified against that database.
> >
> {: .solution}
>
{: .question}

# Create a hybrid assembly

## Bacterial hybrid assembly with **Unicycler**

**Unicycler** builds an initial assembly graph from the Illumina short reads.
Nanopore long reads can span repeats and help resolve connections in the
graph. Both read sets must come from the same isolate.


> <hands-on-title> Assemble *E. coli* GPT24-02 with Illumina and Nanopore reads </hands-on-title>
>
> 1. Run {% tool [Create assemblies with Unicycler](toolshed.g2.bx.psu.edu/repos/iuc/unicycler/unicycler/0.5.1+galaxy1) %} with the following parameters:
>    - *"Paired or Single end data?"*: `Paired Collection`
>        - {% icon param-collection %} *"Select a paired collection"*: `Ecoli_Illumina_cleaned` (output of **fastp**)
>    - {% icon param-file %} *"Select long reads. If there are no long reads, leave this empty"*: the **filtered GPT24-02 ONT FASTQ** from **Filtlong**
>    - *"Select Bridging mode"*: `Normal`
>    - Leave the other settings at their defaults.
>
> 2. When the job finishes, locate the **Assembly FASTA** and the
>    **Final Assembly Graph (GFA)**. If the outputs are in collections,
>    expand them to find the GPT24-02 datasets.
>
{: .hands_on}

> <question-title>Interpret the Unicycler inputs and outputs</question-title>
>
> 1. How can a Nanopore read help Unicycler resolve a repeat in the
>    short-read assembly graph?
> 2. Which **Unicycler** output will you use for **Bandage**, and which will
>    you use for **QUAST**, **CheckM2**, and **StarAMR**?
>
> > <solution-title></solution-title>
> >
> > 1. A long read that spans the repeat and reaches unique sequence
> >    on both sides can help identify which flanking regions are connected.
> > 2. Use the **Final Assembly Graph (GFA)** for Bandage to visualize
> >    connections between segments. Use the **Assembly FASTA** for
> >    **QUAST**, **CheckM2**, and **StarAMR**, which analyze the assembled sequences.
> >
> {: .solution}
>
{: .question}

## Visualize the assembly result with **Bandage Image**

An assembly graph shows sequence segments as nodes and connections between them as edges. Inspecting the graph with **Bandage** {% cite Wick2015Bandage %} can help you identify separate components, loops, and branches that are hard to see in a list of contigs. In a fully resolved assembly of an isolate with a circular chromosome and circular plasmids, the chromosome and each plasmid would appear as separate, unbranched loops.

> <hands-on-title> Draw the assembly graph with Bandage </hands-on-title>
>
> 1. Run {% tool [Bandage Image](toolshed.g2.bx.psu.edu/repos/iuc/bandage/bandage_image/2022.09+galaxy4) %} with the following parameter:
>    - {% icon param-file %} *"Graphical Fragment Assembly"*: **Final Assembly Graph** (`gfa`) from **Unicycler**
>    - Leave the other settings at their defaults.
>
> 2. When the job finishes, open the image output. Look for separate
>    components and note which ones form loops or contain branches.
>
{: .hands_on}

> <tip-title>Choose the graph output</tip-title>
>
> Use the **Final Assembly Graph** for **Bandage**. Keep the **Assembly FASTA**
> for **QUAST**, **CheckM2**, and **StarAMR**. Galaxy labels the two formats `gfa` and
> `fasta`, respectively.
>
{: .tip}

> <question-title> Interpret a circular component </question-title>
>
> If you see a small, separate loop in the graph, what might it represent?
> What additional evidence would support your interpretation?
>
> > <solution-title></solution-title>
> >
> > A small, separate loop is consistent with a circular sequence and
> > could represent a plasmid. A plasmid replicon detected on the
> > corresponding contig would support that interpretation. Reads spanning
> > the circular junction would support the inference of circularity.
> >
> {: .solution}
>
{: .question}

## Assess assembly contiguity with **QUAST**

**QUAST** {% cite Gurevich2013QUAST %} summarizes assembled sequences using metrics such as contig count,
total assembly length, largest contig, and N50. We will run **QUAST** without
a reference genome to describe the contiguity of the *E. coli* GPT24-02 assembly.

> <hands-on-title> Run QUAST on the assembly </hands-on-title>
>
> 1. Run {% tool [Quast](toolshed.g2.bx.psu.edu/repos/iuc/quast/quast/5.3.0+galaxy3) %} with the following parameters:
>    - *"Assembly mode?"*: `Individual assembly (1 contig file per sample)`
>        - *"Use customized names for the input files?"*: `No, use dataset names`
>            - {% icon param-file %} *"Contigs/scaffolds file"*: **Assembly FASTA** from **Unicycler**
>        - *"Reads options"*: `Disabled`
>    - *"Type of assembly"*: `Genome`
>        - *"Use a reference genome?"*: `No`
>    - *"Lower threshold for a contig length (in bp)"*: `500`
>    - In *"Genes"*:
>        - *"Tool for gene prediction"*: `Don't predict genes`
>    - *"Output files"*: `HTML reports` and `Tabular reports`
>
> 2. Open the HTML report. Under **Statistics without reference**, record
>    the **number of contigs**, **total length**, **largest contig**, and **N50**.
>
{: .hands_on}

> <question-title>Explain the N50 of the *E. coli* GPT24-02 assembly</question-title>
>
> 1. How many contigs are in the assembly, and what is its total length?
> 2. Why is the N50 equal to the length of the largest contig? What can
>    you learn from the Bandage graph that the N50 alone cannot tell you?
>
> > <solution-title></solution-title>
> >
> > 1. **QUAST** reports **9 contigs** and a total length of **5,489,945 bp**.
> >    All nine contigs meet the 500 bp cutoff.
> > 2. The largest contig is **5,159,512 bp**, about **94%** of the assembly.
> >    When contigs are ordered from longest to shortest, this first
> >    contig already accounts for at least half of the assembled bases,
> >    so its length is also the **N50**. The **Bandage** graph can show whether
> >    smaller components are separate, circular, or branched.
> >
> {: .solution}
>
{: .question}

## Estimate genome completeness and contamination with **CheckM2**

**QUAST** describes how the assembled bases are distributed among contigs.
**CheckM2** {% cite Chklovski2023CheckM2 %} uses features of the predicted genes to estimate genome
completeness and contamination. We use both reports to assess
different aspects of the same assembly.

> <hands-on-title>Run CheckM2 on the assembly</hands-on-title>
>
> 1. Run {% tool [CheckM2](toolshed.g2.bx.psu.edu/repos/iuc/checkm2/checkm2/1.1.0+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input MAG/SAG datasets"*: the full **Assembly FASTA** from **Unicycler**.
>    - *"Select reference genome"*: select the installed **CheckM2 DIAMOND database**.
>    - *"Treat input files as protein files"*: `No`.
>    - *"Model options"*: `None` to let **CheckM2** choose the model.
>
> 2. Open the **Quality Report** and record the **Completeness** and
>    **Contamination** estimates. Compare the genome size and contig
>    count with the **QUAST** report.
>
{: .hands_on}

> <tip-title>Submit all contigs together</tip-title>
>
> Use the same nine-contig FASTA you used for **QUAST**. **CheckM2** evaluates
> the contigs in that file together and returns one row for the assembly.
>
{: .tip}

> <question-title>Interpret the CheckM2 results</question-title>
>
> 1. What completeness and contamination did **CheckM2** estimate for *E. coli* GPT24-02?
> 2. How can the assembly have a completeness estimate of 100% while
>    still containing nine contigs?
>
> > <solution-title></solution-title>
> >
> > 1. **CheckM2** estimated **100.0% completeness** and **0.26% contamination**.
> >    The reported assembly size was **5,489,945 bp across nine contigs**,
> >    matching **QUAST**.
> > 2. Completeness estimates how much of the expected genome content is
> >    represented. That content can be distributed across several contigs.
> >    Use **QUAST** to describe contiguity and the **Bandage** graph to inspect
> >    how individual components are arranged.
> >
> {: .solution}
>
{: .question}

# Antimicrobial resistance gene screening with **staramr**

In the final step of this tutorial, we use **StarAMR**
{% cite Bharat2022Staramr %} to screen the assembly for resistance genes.
For each detected gene, we can use its contig ID to check whether it
shares a contig with features such as a plasmid replicon. Long reads
can span repeats and help assemble a resistance gene together with its
flanking sequences, making its genomic location easier to infer.

> <hands-on-title> Antimicrobial resistance gene screening </hands-on-title>
>
> 1. Run {% tool [StarAMR](toolshed.g2.bx.psu.edu/repos/iuc/staramr/staramr_search/0.12.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} **Input assembly**: the full **Assembly FASTA** from **Unicycler**
>    - *"Enable scanning for point mutations using the PointFinder database"*: `Escherichia coli (escherichia_coli)`
>    - In *"Advanced options"*:
>        - *"Exclude certain AMR genes from the results"*: `Exclude default list of AMR genes`
>
> 2. Open the **Detailed Summary**. For each acquired AMR gene reported
>    as **Resistance**, record the gene name and its **Contig** ID.
>    Then find the plasmid replicon hits reported as **Plasmid** and
>    compare their contig IDs with those of the AMR genes.
>
{: .hands_on}

> <question-title>Locate AMR genes in the *E. coli* GPT24-02 assembly</question-title>
>
> 1. Which AMR genes did StarAMR detect, and on which contigs?
> 2. Compare the contig IDs for the AMR genes with those for the plasmid
>    replicons. Which genes have evidence of a plasmid location?
>
> > <solution-title></solution-title>
> >
> > 1. StarAMR detected four acquired AMR genes. Contig **3** contains
> >    `aph(3'')-Ib`, `aph(6)-Id`, and `blaTEM-1C`. Contig **5**
> >    contains `mcr-1.1`.
> >
> > 2. The three AMR genes on contig **3** share that contig with the
> >    `IncI1-I(Alpha)` replicon. On contig **5**, `mcr-1.1` shares
> >    the contig with the `IncX4` replicon. These matches support a
> >    plasmid location for all four genes.
> >
> {: .solution}
>
{: .question}


# Bottom line

We used Illumina and Oxford Nanopore reads from *E. coli*
GPT24-02 to build a hybrid assembly with Unicycler. Comparing the reads
before and after preprocessing helped us see how filtering affected
the amount and length of sequence available for assembly. Bandage,
QUAST, and CheckM2 then helped us assess the assembly's structure,
contiguity, estimated completeness, and contamination.

For AMR surveillance, detecting a resistance gene is the first step.
The hybrid assembly also lets us examine the sequence around it.
In *E. coli* GPT24-02, StarAMR detected `mcr-1.1` and an `IncX4` plasmid replicon
on the same contig, supporting a plasmid location for the resistance
gene. Long reads can help preserve this genomic context by spanning
repeats that would otherwise fragment the assembly.
