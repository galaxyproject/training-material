---
layout: tutorial_hands_on

title: "Antibiotic resistance detection"
tags: [nanopore, plasmids]
zenodo_link: "https://doi.org/10.5281/zenodo.3247504"
questions:
  - "How do I assemble a genome with Nanopore data?"
  - "How do I get more information about the structure of the genomes?"
  - "How do I get more information about the antimicrobial resistance genes?"
objectives:
  - "Perform Quality control on your reads"
  - "Perform a genome assembly with Minimap2/Miniasm/Racon"
  - "Determine the structure of the genome(s)"
  - "Perform a scan for antimicrobial resistance genes with Staramr"
time_estimation: "3h"
key_points:
  - "Minimap2, Miniasm and Racon can be used for creating a fast assembly of Nanopore sequences"
  - "Nanopore sequencing is useful for reconstruction of genomes"
  - "Antimicrobial resistance genes are detectable after fast assembly"
contributors:
  - willemdek11
  - shiltemann
---

# Overview
{:.no_toc}

Pervasive use (and misuse) of antibiotics for human disease treatment, as well as for various agricultural purposes, has resulted in the evolution of multidrug resistant (MDR) pathogenic bacteria. The [Center for Disease Control estimates](https://www.cdc.gov/drugresistance/) that in the U.S. alone, every year at least 2 million people get an antibiotic-resistant infection, and at least 23,000 people die. Antibiotic resistance poses a major public health challenge, and its causes and mitigations are widely studied.

Plasmids are small DNA molecules within a cell which are physically separated from chromosomal DNA and can replicate independently.

![Depiction of plasmids](../../images/plasmid-metagenomics-nanopore/plasmids.png)

Plasmids are considered a major vector facilitating the transmission of drug resistant genes among bacteria via [horizontal transfer](https://en.wikipedia.org/wiki/Horizontal_gene_transfer) ({% cite Beatson2014 %}, {% cite Smillie2010 %}). Careful characterization of plasmids and other MDR mobile genetic elements is vital for understanding their evolution and transmission and adaptation to new hosts.

![Illustration of transformation of a bacteria to drug resistance](../../images/plasmid-metagenomics-nanopore/bacterial_transformation.png)

Due to the high prevalence of repeat sequences and inserts in plasmids, using traditional NGS short-read sequencing to assemble plasmid sequences is difficult and time-consuming. With the advent of third-generation single-molecule long-read sequencing technologies, full assembly of plasmid sequences is now possible.

In this tutorial we will recreate the analysis described in the paper by {% cite LiXie2018 %} entitled *Efficient generation of complete sequences of MDR-encoding plasmids by rapid assembly of MinION barcoding sequencing data*. We will use data sequenced by the [Nanopore](https://nanoporetech.com/) MinION sequencer.

The assembly is performed with **Minimap2** {% icon tool %} ({% cite Li2018 %}),
**Miniasm** {% icon tool %} ({% cite Li2016 %}) and **Racon** {% icon tool %} ({% cite Vaser2017 %}).
The downstream analysis is done with **Nanoplot** {% icon tool %} ({% cite DeCoster2018 %}),
**Bandage** {% icon tool %} ({% cite Wick2015 %}) and **PlasFlow** {% icon tool %} ({% cite Krawczyk2018 %}).

A schematic view of the workflow we will perform in this tutorial is given below:

![Workflow representation of this tutorial](../../images/plasmid-metagenomics-nanopore/Workflow.png)


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


{% include snippets/warning_results_may_vary.md %}

# Obtaining and preparing data

In this tutorial we use metagenomic Nanopore data, but similar pipelines can be used for other types of datasets or other long-read sequencing platforms.

> ### {% icon comment %} Background: Nanopore sequencing
>
> Nanopore sequencing has several properties that make it well-suited for our purposes
>
> 1. Long-read sequencing technology offers simplified and less ambiguous genome assembly
> 2. Long-read sequencing gives the ability to span repetitive genomic regions
> 3. Long-read sequencing makes it possible to identify large structural variation
>
> ![How nanopore sequencing works](../../images/plasmid-metagenomics-nanopore/sequence_method.jpg) <br><br>
>
> (Image credit: [Nanopore sequencing: The advantages of long reads for genome assembly](https://nanoporetech.com/sites/default/files/s3/white-papers/WGS_Assembly_white_paper.pdf?submissionGuid=40a7546b-9e51-42e7-bde9-b5ddef3c3512 ))
{: .comment}


## Understanding our input data

In this tutorial we are interested in determing the antimicrobial resistance genes.

As training data we use plasmids from a dataset (created by {% cite LiXie2018 %}) used for evaluation of the efficiency of MDR plasmid sequencing by MinION platform. In the experiment, 12 MDR plasmid-bearing strains were selected for plasmid extraction, including *E. coli, S. typhimurium*, *V. parahaemolyticus*, and *K. pneumoniae*.


> ### {% icon details %} More details about datasets
> Overnight cultures (100 mL) were harvested and subjected to plasmid extraction using the QIAGEN Plasmid Midi Kit. The extracted plasmids were dissolved in ultrapure distilled water, and concentrations were measured by Qubit 3.0 Fluorometer with a dsDNA BR Assay Kit. The plasmids were stored in –20°C until library preparation.
>
> Library preparation was performed using the Rapid Barcoding Sequencing Kit (SQK-RBK001) according to the standard protocol provided by the manufacturer (Oxford Nanopore). Briefly, 7.5-μL plasmid templates were combined with a 2.5-μL Fragmentation Mix Barcode (1 barcode for each sample). The mixtures were incubated at 30°C for 1 minute and at 75°C for 1 minute. The barcoded libraries were pooled together with designated ratios in 10 μL; 1 μL of RAD (Rapid 1D Adapter) was added to the pooled library and mixed gently; 0.2 μL of Blunt/TA Ligase Master Mix was added and incubated for 5 minutes at room temperature. The constructed library was loaded into the Flow Cell R9.4 (FLO-MIN106) on a MinION device and run with the SQK-RBK001_plus_Basecaller script of MinKNOW1.5.12 software. The run was stopped after 8 hours, and the flow cell was washed by a Wash Kit (EXP-WSH002) and stored in 4°C for later use.
>
> To obtain high-quality short read data, paired-end (2 × 150 bp) libraries were prepared by the focused acoustic shearing method with the NEBNext Ultra DNA Library Prep Kit and the Multiplex Oligos Kit for Illumina (NEB). The libraries were quantified by employing quantitative PCR with P5-P7 primers, and they were pooled together and sequenced on the NextSeq 500 platform according to the manufacturer's protocol (Illumina).
>
> Although a local basecaller script was used during the run, there was still a small amount of reads that were not basecalled due to the generation of raw data in a rapid mode. Albacore basecalling software (v1.0.3) was used to generate fast5 files harboring the 1D DNA sequence from fast5 files with only raw data in the tmp folder. Also, the read_fast5_basecaller.py script in Albacore was used to de-multiplex the 12 samples from basecalled fast5 files (except the files in fail folder) based on the 12 barcodes in SQK-RBK001. The Poretools toolkit was utilized to extract all the DNA sequences from fast5 to fasta format among the 12 samples, respectively (Poretools, RRID:SCR_015879).
>
{: .details}

> ### {% icon comment %} Dataset details
> Because of the large size of the original datasets (1.15 GB) you are given 1 of the 12 plasmids
> files.
> <br><br>
> This sequence file is 51 MB of nanopore sequences. The 10026 reads found in this file contain 49190798 nucleotides.
> As mentioned before nanopore sequences are long reads and this is confirmed by a mean read length of 4906.3.
{: .comment}

## Importing the data into Galaxy

Now that we know what our input data is, let's get it into our Galaxy history:

> ### {% icon hands_on %} Hands-on: Obtaining our data
>
> 1. Make sure you have an empty analysis history. Give it a name.
>
>    {% include snippets/create_new_history.md %}
>
> 2. **Import Sample Data** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3247504.svg)](https://doi.org/10.5281/zenodo.3247504)
>    ```
>    https://zenodo.org/record/3247504/files/RB01.fasta
>    https://zenodo.org/record/3247504/files/RB02.fasta
>    https://zenodo.org/record/3247504/files/RB03.fasta
>    https://zenodo.org/record/3247504/files/RB04.fasta
>    https://zenodo.org/record/3247504/files/RB05.fasta
>    https://zenodo.org/record/3247504/files/RB06.fasta
>    https://zenodo.org/record/3247504/files/RB07.fasta
>    https://zenodo.org/record/3247504/files/RB08.fasta
>    https://zenodo.org/record/3247504/files/RB09.fasta
>    https://zenodo.org/record/3247504/files/RB10.fasta
>    https://zenodo.org/record/3247504/files/RB11.fasta
>    https://zenodo.org/record/3247504/files/RB12.fasta
>    ```
>    {% include snippets/import_via_link.md %}
>
> 3. **Build a list collection**
>
>    {% include snippets/build_list_collection.md %}
>
{: .hands_on}

# Quality Control

## NanoPlot to explore data


The first thing we want to do is to get a feeling for our input data and its quality. This is done
using the NanoPlot tool. This will create several plots, a statisical report and an HTML
report page.

![NanoPlot example](../../images/plasmid-metagenomics-nanopore/NanoPlot.png)

> ### {% icon hands_on %} Hands-on: Plotting scripts for long read sequencing data
>
> 1. **NanoPlot** {% icon tool %} with the following parameters
>   - *"Type of the file(s) to work on"*: `fasta`
>   - *"files"*: `Data list` you just created
>
{: .hands_on}

The `Histogram Read Length` gives an overview of the read distribution of all files in the given data collection.
![NanoPlot Output](../../images/plasmid-metagenomics-nanopore/NanoPlot_output.png) <br><br>


> ### {% icon question %} Question
>
> What was the mean read length for this (RB01) sample?
>
> > ### {% icon solution %} Solution
> > 4906.3
> >
> > This can be determined by looking at the NanoStats or HTML output of NanoPlot.
> {: .solution }
{: .question}


For more information on the topic of quality control, please see our training materials
[here]({{site.baseurl}}/topics/sequence-analysis/)

# De Novo Assembly

## Pairwise alignment using Minimap2

In this experiment we used Nanopore sequencing; this means that sequencing results in long reads with overlap.
To find this overlap, Minimap2 is used. Minimap2 is a sequence alignment program that can be used for different
purposes, but in this case we'll use it to find overlaps between long reads with an error rate up to ~15%.
Typical other use cases for Minimap2 include: (1) mapping PacBio or Oxford Nanopore genomic reads to the human genome;
(2) splice-aware alignment of PacBio Iso-Seq or Nanopore cDNA or Direct RNA reads against a reference genome;
(3) aligning Illumina single- or paired-end reads; (4) assembly-to-assembly alignment; (5) full-genome alignment
between two closely related species with divergence below ~15%.

Minimap2 is faster and more accurate than mainstream long-read mappers such as BLASR, BWA-MEM, NGMLR and GMAP and
therefore widely used for Nanopore aligning. Detailed evaluations of Minimap2 are available from
the [Minimap2 paper](https://doi.org/10.1093/bioinformatics/bty191).

![Pairwise alignment](../../images/plasmid-metagenomics-nanopore/Minimap2.png)


> ### {% icon hands_on %} Hands-on: Pairwise sequence alignment
>
> 1. **Map with minimap2** {% icon tool %} with the following parameters
>   - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>   - *"Use the following data collection as the reference sequence"*: `Created dataset collection`
>   - *"Select analysis mode (sets default)"*: `Oxford Nanopore all-vs--all overlap mapping`
>   - *"Select an output format"*: `paf`
>
{: .hands_on}

This step maps the Nanopore sequence reads against itself to find overlaps. The result is a paf file.
PAF is a text format describing the approximate mapping positions between two
set of sequences. PAF is TAB-delimited with each line consisting of the
following predefined fields:

|Col|Type  |Description                               |
|--:|:----:|:-----------------------------------------|
|1  |string|Query sequence name                       |
|2  |int   |Query sequence length                     |
|3  |int   |Query start (0-based)                     |
|4  |int   |Query end (0-based)                       |
|5  |char  |Relative strand: "+" or "-"               |
|6  |string|Target sequence name                      |
|7  |int   |Target sequence length                    |
|8  |int   |Target start on original strand (0-based) |
|9  |int   |Target end on original strand (0-based)   |
|10 |int   |Number of residue matches                 |
|11 |int   |Alignment block length                    |
|12 |int   |Mapping quality (0-255; 255 for missing)  |

View the outcome of RB12, it should look something like this:
```
channel_100_69f2ea89-01c5-45f4-8e1b-55a09acdb3f5_template	4518	114	2613	+	channel_139_250c7e7b-f063-4313-8564-d3efbfa7e38d_template	3657	206	2732	273	2605	0	tp:A:S	cm:i:29	s1:i:240	dv:f:0.2016	rl:i:1516
channel_100_69f2ea89-01c5-45f4-8e1b-55a09acdb3f5_template	4518	148	1212	+	channel_313_35f447cb-7e4b-4c3d-977e-dc0de2717a4d_template	3776	2433	3450	218	1064	0	tp:A:S	cm:i:31	s1:i:210	dv:f:0.1291	rl:i:1516
channel_100_69f2ea89-01c5-45f4-8e1b-55a09acdb3f5_template	4518	251	1328	+	channel_313_a83f7257-52db-46e4-8e2a-1776500c7363_template	3699	2327	3382	208	1082	0	tp:A:S	cm:i:29	s1:i:203	dv:f:0.1363	rl:i:1516
```

## Ultrafast de novo assembly using Miniasm

The mapped reads are ready to be assembled with Miniasm. Miniasm is a very fast Overlap Layout Consensus based de novo assembler for noisy long reads.
It takes all-vs-all read self-mappings (typically by Minimap2) as input and outputs an assembly graph in the GFA format.
Different from mainstream assemblers, miniasm does not have a consensus step.
It simply concatenates pieces of read sequences to generate the final sequences.
The optimal case would be to recreate a complete chromosome or plasmid.
Thus the per-base error rate is similar to the raw input reads.
![Pairwise alignment](../../images/plasmid-metagenomics-nanopore/Miniasm.png)

> ### {% icon hands_on %} Hands-on: De novo assembly
>
> 1. **miniasm** {% icon tool %} with the following parameters
>   - *"Sequence Reads"*: `Input dataset collection`
>   - *"PAF file"*: `Output Minimap dataset collection` created by the Minimap2 tool
>
{: .hands_on}

The `Assembly Graph` output file gives information about the steps taken in the asssembly.

The output should look like:

```
S	utg000001l	GAAATCATCAGGCGTTTTTCACGATATGGACGGGAAGATGCGGAAATAGGCAGGAGGACATAGAAATGCCTGAGGGGTCTGGGATGGTGCGGGCAACGGATGTTATGGTAAATAAGCTTCCGTTGGTAAACCTGTAAGTTTTCAGGAACGAGACTCGTTTAGAACATCTAAAAAGCACATGAATGCTGCTATAGAAGCGACTGTTTGATGGTTCATGGTTTTATATTAAGGTAGATGAAAACTCAATACAGCGCTATTTAGAGACGCTTAAATGGCATTCTTGTATTTTAACGATGCAGAGTAACAGGCTCTAGCTAGTTATAGTACATCGAAAAGTTCATATAGGGACAAAACAAAGGGACTGAAATACATAGCCATAAGCTCGCTTCAAAGTCCTAACCACCAGCTATATTGCATGGGTTTGGTTAAGACATACGTCATGTCTTTGGGGCATATTATGGTGAAAAACGCGTTCATCTACACCAATATTTATACATTTCTCTTCAATAATGACGCCGTTTCCAGCCAACTTTGCTATATGGGTGACAATTTTCGCTTTTAAACGGAGATGTCTCCATGATTGAGCAGATCTAAAAATCTGCATTTTGTGCTGTGCAGGCATTACCATTCTTTCTATGCACAAACGTGCCGCGTCTTGCATTCTGCTTTGTATGCGGAGGTGCTCAATGAGGCTTTAACTACTCCCCTCCGTTAACGCCGTTCGTGTTTTAAACCGCCTTCTCGATAAAGAGTTTCACAGGTGCTGAAGAAATTCGCCTTTACGATGATCAGGGCGGGTGATTGGTACATCAACAACAGGCTTCTTGTTCTGTTTTGAGAGCATGAAGATCAGCCGCGATGGTAATATGTTCGTGTATTTGTCTCAAATACATTGATGACTATTACAACCGTCGTGAGTATCAGTTTAGTCAGCAGAGGCATGCTGTCGCTGTGCTCGTCAGCATAGAAAAGCCTCTACGGCCGCTCGTTTGCTTTTATCGCTAGGTGAGCTTGTCCCAGAAAGTATCCGTATGATTCTTGCGCTAATGATGATTAGGGGATATTATGTTTCTGGGGTTATTTTGTTGTTTTGCTGGCTTTCTCTACTATTAGAGAGCTTGGGGTTTCTTGGAAACGTTGAATGGGGTCCTCTAGCATGATTTGTTTGAGCTTCATCTGGTCTATGAAGCAGTCAGGAAGAGTCTTAAGCACGATGTTGGGTCTTAAACTAGACTAAAAATATAATATGCGAGCTTTTAACTCATAGACGTGTAGTTTTCACTCTTTTTGCTGGTGTTCTAAGCAAAAACTTCTGCATTACGCCTGTTTGTGGTTGCCTGTTATAAGATAACAGTAAAATCCATGCGTTTCTAAGGAGAAAATGAACAAAAGCTATCACTAAATCGTATTTCCAATAGAGGTTTCGTCCTCTTTGAAAGGGGTAGTTCATGAAAATTAGTTCATGGAAATGTTTATGATTTAGGTTGGAGTTATTAGGCGCAAAGAAGAAATATAAAAGGATAATTGAAGTTCTAAATCATTATGGTGATTGAGAGGCGTTATGAGCACAGCTAAGGAAAGTTAAGGAAAAAGGCACTCCAATGAGTAGTGCCAATGAATATTCTGGTCGTGATAGTATGAATCAAGATTACATCTAGGCAATCCATGGATATTATTATAGGTTTCTGTGGGCTGTTGGTTCGGTATTAAGTCAATTGTAGATGTAACCGAAGAAGAACTGAAGATTTAGGATACGAAGTATATAGGATAAGAGTTAGCGATCTTATCAGAAACTTTTTCTGAGTGGTGGGGATACATCTACTCCATTCAAAAGATATAATTCCTTGCGGGATTTTGGTGATTCCTTAAGATGTAGAGTACTACCCACACATACTAGCTGAAGCTGTAATTCATGAAATAAAGATAGAAAAATCTAGAAGGAAAATCAAAGGTTAAAGAAACGACTTTCCTTATGAATCAACTAAAAACATAAGATGAAGTGAGTTGCCCAAGTGGTCTAACATAATTTCTATTTATTAGGTGATAGGGCGAAACTGAAAGAAAAGAAATCTTCGAGATAGGGGTTTAAGCAAATCTGAAGTTGATATAATCATCCATCATGACAGGAGTCTGAGCACCCTTATGGCCAGCAAACAGCAAAAGCTATACTTGATGCTGACTATTTTGTTAAAATAATCAGGCCAGAGCTCACCTTAAATCTAAGTTAAACGCTTTCTTGGTCTATTCATGGGCTTAATGGCCTTACGCCAAACGCGCATGAGGAAGACATGTACGCAGCTCATTCTGCCTCTTTACAATCTGCATGTCTTTCCAGACAGGTAGGTGCTGCTATTTTAGATATAAAGGAAATTTAATTGCTGTTGGGCGCAATGATGTTCTAAATTTAACAGTGTGGCCTTTATAGTGCTGATGATGGTGAATGATCATCGCTGTGTCTACAAAGGTGGCAAATGTTACAATGAAGCAAGAAAATTAAATAAAGATAAATACAGCAAATATTGCTTTTCTGATAGGTTAAGAATTTGTCTGTTCAGTTAACTCAAGAACAAGCAGAAAAGGTTGCTGGAATGTCTACGAAGGTACTCCTGTATCATCAATCGACAGTACTCGCTCTATTCATGCTGAAATGGATGCTTCACAGCCCATCTTGCTAGACTCGAGATAGTGGTTTTGAAGATAAGGTTTTATATACAACTACGTTCCCCTGTCATAATTGTGCAGGGTATATTGTTGCTGTTGGTATCAATAAAGTTGTTTATTGAACATATGAAAAAGCAGCTTTAGAGTTGCTGATGATGCAATAACTAGAAGTTAACGAAAGTGGTAAAGTTTTATTTGAGTCCTTTGAAGGGGTTTCTACTGGAATATCATAGAAGTTCTTCTTTTCTACAGATGAAAGGATAGTAAGGTAATGCGGAGAAATATTCAACGAAGTATAAGAATCATATTGATATACAGTTTATTGATAGTTGGATTTTAATCGTGTGGCGGATATTTTCATTAACACAAAACAGATGCTGAAGCTCAAGTAGCCGTTCCATCTGGGACATAGTTAATCCGCTTCTGCCACCGATTCTCGACACCACCAACATCGGGTATGAATCTGTGACTCTGATGTTACAGAGCTTAATCTTGTGTACCAAAAACCACCATACCAACGGTGGTTTTCTCTGAGCTACTGCTCTTTGAGCCGAGGTAACTGGCTTGGAGGAGCTCAGTCACCAAAACTTGTTCTTTCAGTTTAGCCTTAACAGGCGCACAACTTAAGATGCTCCTCTAAATCAGTTACCAGTGGCTGCTGCCAGTGGCGCTTTGTCGTGCCTTCTGGGTTGGACTCAAGACGATAGTTACCGGGGCAGGCGCAGCGGTCGGACTGAACGGGGTTCGTGCGCTGGTCAGCTTGGGCGAACTGCCTGCGGAACGCTGAGTGTCAGGCGATGAGTAAAGTTTAAACAGCCATAACAGCAGGTGAAGAGCACCCCGGTGCCAAACCGAAGGCAGGAACAGAGCGCAGGGAGCTACCGGGGAAACGCAGGATCTTTATAGTCCTGTCCGGGTTTCGCCACCATGATTTGAGCGTCAGATTCTGTGATGCTTGTCAGGGGGCGGAGCCTATCGGAAAACCGGCAACGCGGCCTTTCTTGTTGCTTCTCCCATTCTCTCATTGCCTATGGCTTAATGTCTCTGTTCCCTCCGCTCGCCGCAGCGAACGACCGAGCGGGCGAGTCAGTGAGCGGAAGCGGAATATCTGCGGGCTTCTCTTTGGCACCGTACGCCATAGCGCATTTTAATACGATGCAGAATAGGGCGGGTACGCCGCAAAGTGACGTCACCTGACGTTCTGGATTACAAAGGTTAAAGCAGCGGACAGGAAATGTTTTGTGCCTAGCTATGCTATTCACAAGTAGCAGGACAGATGTGTTTTGGAGTACCATGACATAAAGACTCTGGTTAGCTCCCGCATGTTTGCTATCTTTGCTACCATTTGCTACTTTTGCTATCTCCAGGGTCTTCTGGTCTTTCAGTTGCTATCCTTTTGCTACTTTTGCTATCAAAATGCTACCTCTCCCTTCTTGCAATAAATGACCAAGGCACGTTAGCGATGTATACATCGCGCACCATAAGTGCCCCTTCCGAATCCCGTTCATCACTGTGTGATTAATCGAGGTTAAATCGACTACCAACGTAAAATCCTGCCATGCCTGACGGCTGACAACGCCCTCACGGTCGCGCAGTAAACTGCGCGCCCCTTTCACCGTTATCGGTGGGGTTCGTGGTGGCTTTCAGGGCGCTCGCCGGATTTTAACCGCTAAAATGAGCGATCCATGCGTTCGTG	LN:i:4357
a	utg000001l	0	channel_364_204a2254-2b6f-4f10-9ec5-6d40f0b870e4_template:101-4457	+	4357
```

## Remapping using Minimap2

The Assembly graph created can be used for mapping again with minimap2, but first the graph should be transformed to FASTA format

> ### {% icon hands_on %} Hands-on: Pairwise sequence alignment
>
> 1. **GFA to Fasta** {% icon tool %} with the following parameters
>   - *"Input GFA file"*: the `Assembly Graph` created by the Miniasm tool
>
> 2. **Map with minimap2** {% icon tool %} with the following parameters
>   - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>   - *"Use the following dataset as the reference sequence"*: `FASTA file` created by the GFA to Fasta tool
>   - *"Select fastq dataset"*: `RB01.fasta`
>   - *"Select an output format"*: `paf`
>
> > ### {% icon question %} Question
> >
> > How many reads are over after the use of Minimap2 and Miniasm? (Hint: run the Nanoplot tool on the output of GFA to Fasta)
> >
> > > ### {% icon solution %} Solution
> > > 456
> > >
> > > This can be determined by looking at the NanoStats output of NanoPlot.
> > {: .solution }
> {: .question}
{: .hands_on}

## Ultrafast consensus module using racon

The mapped reads can be improved even more using Racon to find a consensus sequence.
Racon is a standalone consensus module to correct raw contigs generated by rapid assembly methods which do not include a consensus step.
Racon generates genomic consensus which is of similar or better quality compared to the output generated by assembly methods which employ both error correction and consensus steps.
This while providing a speedup of several times compared to those methods.
It supports data produced by both Pacific Biosciences and Oxford Nanopore Technologies.
![Consensus Module](../../images/plasmid-metagenomics-nanopore/Racon.png)

> ### {% icon hands_on %} Hands-on: Consensus module
>
> 1. **Racon** {% icon tool %} with the following parameters
>   - *"Sequences"*: `RB01.fasta`
>   - *"Overlaps"*: the latest `PAF file` created by the Minimap2 tool
>   - *"Target sequences"*: the `FASTA file` created by the GFA to Fasta tool
>
{: .hands_on}

The `Racon` output file gives the final contigs.

The output should look like:

```
>utg000001l LN:i:4399 RC:i:43 XC:f:1.000000
GGGCAACGGATGTTATGGTAAATAGCTTCCGTTGGTAAGTATCTCAATGTCTGAAGTTTTCAGGGTTTGAGACTCGTTTAGAACATCTAAAAAAAAGCACATGAATGCTGCTATAGAAGCGACTGTTTGATGGTTCATGGTTTACCATTAAGGTAGATGAAAACTCAATACAGCGCTATTTAGGGCGCTTAAATGGCATTCTTATTATTTTTAACGATGCAGAGTAACAGGCTCTAGCTAGGTTATGGTACACTCGAAAAGTTCCCTTATAGGGACAAAACAAAGGGACTTTGAAGCTACGCAACCGCCAAGCTCGCTTCAAAGTCCCTGTACCGCAGCTATATTGCATGGGTTTGGTTAAGACATACGTCATGTCTTTGGGGCATATTATGGTGAAAACGCGTTCGCGTCTACACCAATATTTATACATTTCTCTTCAATATTGACGCCGTTTCCAGCCAACTTTTGCTATATGGGTGACAATTTTCGCTTTTAAAACGGAGATGTCTCCATGATTGAGCATGATCCTAAAATCTGCATTTTGTGCTGTGCAGGCATTACCATTCTTTCTATGCACAAACGTGCCGCGTCTTGCATTCTGCTTTTGTATGCGGAGGTGCTCAATGAGAGCTTTAACTACTCCCCTCGTTATGCTAACGCCGTTCGTGTTTTAAACCGCCTTCTCGATAAAGAGTTTCTGAGTGCTGAAGAAATTCGCCTTTACGATGATCAGGGCGGGTGATTGAGTACATCAACAACAGGCTTCTTGTTCCTGTTTTGAGCATGAAGATCAGCCGCGATGGTGAGTATGTTCGTGTCTGGTATTTGTCTCAAGAAGCCATTGATGACTATTACAACCGTCGTGAGTATCAGTTTAGTCAGCAGAGGCATGCTGTCGCTGTTGCTCGTCAGCATAGAAAAGCCTCTACGGCCGCTCGTTTGCTTTTATCGCTAGGTGAAACTGTCCCAGAAAGTATCCGTATGATTCTTGCTGCTAATGATGATTAGGAGGATATTATGTTTCTGGGGTTATTTGTTGTTTTGCTGGCTTTCTCTACCTATTAGAGAGCTTGGGGTTTCTTCACGGAAACGTTGAATGGGAGTGCCTCTCGCCATGATTTGTTTTAGGCTTCATCTGGTCTATGAAGCAGTCAGGAGAGTGCTTAAGCACGATGTTTAGGTCTTAAACTAGACTAAAAGCTCGCAATATGCGAGCTTTTATTTATTCATAGACTTGTGTAGTTTTCACTCTTTTTGCTGGTGTTCTAAGCAAAAACTTCTGCGCATTACGCCTGTTTGTGGTTGCCTGTTATAAAGATAACAGTAAAATCCATGCGTTTCGCAGGAGAAAATATGAGCAAAGAGCTATCACTAAATCGTATTTCCAATAGAGTTTCGTCCTCTTTGAAGAGGGGTAGTTCATGAAAGAGTGGTTCATGGAAATGTTTATGATTTGGAGTTGGAGTTATTAAGGCGCAAAAGAGAAGAAGCTATAAAAGGATAATTGAAGATTCTAAATCATTAAAATGGTGATTGAGAGGCGTTATGAGCACAGCTAGGAAGTTAAAGGAAAAAAAGCACTCCAATGAGTAGTGCCAATGAAGATATTTCTGGTCGTGATAGTATGAATCAGATTACATCTAGGCAATCCATGGATATTATTATAGGTTTCTGTGGGGCTGTTGGTTCCGGTATTAAGTCAATTGTAGATGTAACCGAAGAAGAACTGGAAGATTTAGGATACGAAGTATATAGGATAAGAGTTAGCGATCTTATGCAGAAACTTTTTTTCTGAGTGGTGTGGGGATACATCTACTCCATTCAAAGATATAATTCCTTGCAGGATTTTGGTGATTCCTTAAGATGTAAGTACTACCCACACATACTAGCTGAAGCTGTAATTCATGAAATAAAGATAGAAAAATCTAGAAGGAAAATCAAAAGAGTTAAAGAAACGAGCTTTCCTTATAGATCAACTAAAGCACCAAGATGAAATTGAGTTGCTTAGAATGGTCTATCAACATAATTTCTATTTATTAGGTGTTATTAGGAGCGAGACTGAAAGAAAAGAAATCTTCGAGATGAGGGTTTAAGCAAATCTGAAGTTGATATAATCATCCATCATGACAGGAAGTCTGAGCACCCTTATGGCCAGCAAACAGCAAAAGCTATACTTGATGCTGACTATTTTGTTAAAAATAATCAAGGCCAGAGCTCACCTTAAATCTAAAGTTAAACGCTTTCTTGGTCTAATTCATGGGCTTAATGGCCTTACGCCAAACGCGCATGAGAAGGGCATGTACGCAGCTCATTCTGCCTCTTTACAATCTGCATGTCTTTCCAGACAGGTAGGTGCTGCTATTTTAGATATAAAGGAAATTTAATTGCTGTTGGGCGCAATGATGTTCCTAAATTTAACGGTGGCCTTTATAGTGCTGATGATGGGTGAATGATCATCGCTGTGTCTACAAAGGTGGCAAATGTTACAATGAAGTAAGAAAATTAAAATAAAAGATAAAATACAGCAAATATTGCTTTCTGATAAGGTTAAGGAGTTGTCTGTTCAGTTAACTCAAGAACAAGCAGAAAAGGTTGCTGGAATAATCTACGAAGGTACTCCTGTATCATCAATCATTGAGTACTCTCGCTCTATTCATGCTGAAATGGATGCTATTACATCTCTTGCTAGACTCGGGAATGGTGGTTTTGAAGATAAGGTTTTATATACAACTACGTTCCCCTGTCATAATTGTGCGAGACATATTGTTGCTGTTGGTATCAATAAAGTTGTTTACATTGAACCTTATGAAAAAAGCTTAGCTTTAGAGTTGCACGATGATGCAATAACTGAAGTTAACGAAAGTGGTAAGGTTTTATTTGAGTCCTTTGAAGGGGTTTCACCAAGGAGATATCATAAGTTCTTCTTTTCTACAGATGAAAGGAAGGATAGTAAAGGTAATGCGGAGAAATATTCAACGAAGTATAAGAATCATATTGATATACAGTTTATTGATAGTTATGGATTTTGAAAATCGTGTGGCGGATATTTTCATTAGTAACACAAAACAGATGCTGAAGCTCAAGTAGCCGTTCCATCTGGGACGCCAGTTAATCCGCTTCTGCCACCGATTCCTCCGACGCCACCAACATCGGTATGAAGATCTATTGACTCTGATGTTACACAGAGCTTAATCTTGCAAAGAAGAAAAACCACCGCTACCAACGGTGGTTTTCAGGTTCTCTGAGCTACCAACTCTTTGAACCGAGGTAACTGGCTTGGAGGAGCTCAGTCACCAAAACTTGTTCTTTCAGTTTAGCCTTAACAGGCGCACAACTTCAAGACTAACTCCTCTAAATCAGTTACCAGTGGCTGCTGCCAGTGGCGCTTTGTCGTGCCTTCTCGGGTTGGACTCAAGACGATAGTTACCGGGGCAGGCGCAGCGGTCGGACTGAACGGGGGTTCGTGCATACAGTCCAGCTTGGAGCGAACTGCCTACCCGGAACTGAGTGTCAGGCGTGGAATGAGACAAACGCGGCCATAACAGCGGAATGACACACCGGTAAACCGAAAGGCAGGAACAGGAGAGCGCACGAGGGAGCTACCGGGGAAACGCAGGATCTTTATAGTCCTGTCGGGTTTCGCCACCACTGATTTGAGCGTCAGATTCTGTGATGCTTGTCAGGGGGCGGAGCCTATGGAAAACCGGCAACGCCGCGGCCTTTCTTGTTGCTTCTCCATTCTCTCATTGCCTATGGCTTAATGTCTCTGTTCCCTCCGCTCGCCGCAGCCGAACGACCGAGCGGAGCGAGTCAGTGAGCGAGGAAGCGGAATATCTGCGGGCTTCTCTTTGGCACCGTACGCGCCATAGCTGCATTTTTAATACGATGCAGGATAGGGCGAACGCCGCAAAGTGACGTCACCTGACGTTCTGGATTACAAAGGTTAAAGCAGCGGACAGGAAATGTTTTGTGCCCACAGCTATGCTGGCTTATTCACAAAACAGCGGGACAAGATGTGTTTTGGAGTACCGCCGGCCCGCTAAAGACTCTGGTTAGCTCCCGCATGTTTGCTATCTTTTTGCTACCATTTTGCTACTTTTTTGCTATCTCCAGAGTCTTCTGGTCTTTCGGTTGCTATCCTTTTGCTACTTTTTTGCTATCAAAATGCTACCTCTCCCTTCTTGCAATAAATGACCAAGGCACGTTAGCGATGTATACATCGCCGCGCCCATAAGTGCCCCTTCGAATCCCGTTCATCACTGTGCTGATTAATCGAGGTTAAATCCGACCCCTTAGCAATAAATCCTGCCATGCCTGACGGGCTGCCAACGCCTCACGGTCGCATTGAGTAAACTGCGCACCCTTTCACCGTTATCGGTGGGGTTCGTGGTGGCTTTCA
```

# Species and plasmids

## Prediction of plasmid sequences and classes using PlasFlow

To determine whether the contigs are chomosomal or plasmid DNA PlasFlow can be used. Furthermore, it
assigns the contigs to a bacterial class.

PlasFlow is a set of scripts used for prediction of plasmid sequences in metagenomic contigs.
It relies on the neural network models trained on full genome and plasmid sequences and is able to differentiate between plasmids and chromosomes with accuracy reaching 96%.

![Pairwise alignment](../../images/plasmid-metagenomics-nanopore/PlasFlow.png)


> ### {% icon hands_on %} Hands-on: Prediction of plasmid sequences
>
> 1. **PlasFlow** {% icon tool %} with the following parameters
>   - *"Sequence Reads"*: the `contig file` created by the Racon tool
>
> > ### {% icon question %} Question
> >
> > What is the classification of utg000144c? (Hint: Check the probability table created by PlasFlow)
> >
> > > ### {% icon solution %} Solution
> > > plasmid.Proteobacteria
> > >
> > > This can be determined by looking at the 5th column of the probability table.
> > {: .solution }
> {: .question}
{: .hands_on}


The most important output of PlasFlow is a tabular file containing all predictions (specified with `--output` option), consisting of several columns including:

```
contig_id 	contig_name 	contig_length 	id 	label 	...
```

where:

- `contig_id`is an internal id of sequence used for the classification
- `contig_name` is a name of contig used in the classification
- `contig_length` shows the length of a classified sequence
- `id` is an internal id of a produced label (classification)
- `label` is the actual classification
- `...` represents additional columns showing probabilities of assignment to each possible class

Additionaly, PlasFlow produces fasta files containing input sequences binned to plasmids, chromosomes and unclassified.

## Visualising de novo assembly graphs using Bandage

To determine whether the contigs are chomosomal or plasmid DNA Bandage can give a clear view of the assembly.

Bandage (a Bioinformatics Application for Navigating De novo Assembly Graphs Easily), is a program that creates visualisations of assembly graphs.
Sequence assembler programs (such as Miniasm, Velvet, SPAdes, Trinity and MEGAHIT) carry out assembly by building a graph, from which contigs are generated.
By visualisation of these assembly graphs, Bandage allows users to better understand, troubleshoot and improve their assemblies.
![Bandage gui](../../images/plasmid-metagenomics-nanopore/Bandage.png)

> ### {% icon hands_on %} Hands-on: Visualising de novo assembly graphs
>
> 1. **Bandage image** {% icon tool %} with the following parameters
>   - *"Graphical Fragment Assembly"*: the `Assembly graph` created by the Miniasm tool
>
{: .hands_on}

The Assembly graph image shows one hypothetical plasmid, where the other sequences seem to be chromosal DNA.
![Bandage output](../../images/plasmid-metagenomics-nanopore/Bandage_output.png)

# Antibiotic Resistance

## Scans genome contigs for antimicrobial resistance genes

To determine whether the contigs contain antimirocbial resistance genes (AMR) staramr can be used.
Staramr (*AMR) scans bacterial genome contigs against both the ResFinder and PointFinder databases (used by the ResFinder webservice)
and compiles a summary report of detected antimicrobial resistance genes.

![Pairwise alignment](../../images/plasmid-metagenomics-nanopore/StarAmr.png)


> ### {% icon hands_on %} Hands-on: Prediction of AMR genes
>
> 1. **staramr** {% icon tool %} with the following parameters
>   - *"genomes"* the `contig file` created by the Racon tool
>
> > ### {% icon question %} Question
> >
> > Which contig contains the Gene: dfrA17? (Hint: Check the resfinder.tsv created by staramr)
> >
> > > ### {% icon solution %} Solution
> > > utg000144c
> > >
> > > This can be determined by looking at the 2nd and 7th column of the resfinder.tsv.
> > {: .solution }
> {: .question}
{: .hands_on}

There are 5 different output files produced by `staramr`:

1. `summary.tsv`:  A summary of all detected AMR genes/mutations in each genome, one genome per line.
2. `resfinder.tsv`: A tabular file of each AMR gene and additional BLAST information from the **ResFinder** database, one gene per line.
3. `pointfinder.tsv`: A tabular file of each AMR point mutation and additional BLAST information from the **PointFinder** database, one gene per line.
4. `settings.txt`: The command-line, database versions, and other settings used to run `staramr`.
5. `results.xlsx`: An Excel spreadsheet containing the previous 4 files as separate worksheets.

The summary file is most important and gives all the resistance genes found.

# Conclusion
You have now seen how to perform an assembly on Nanopore sequencing data, with some extra analysis tools. You have worked your way through the following pipeline:
![Workflow representation of this tutorial](../../images/plasmid-metagenomics-nanopore/Workflow.png)
