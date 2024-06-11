---
layout: tutorial_hands_on

title: 'Very Short Introductions: QC'
zenodo_link: https://zenodo.org/records/10870107
questions:
  - 'How do I know my sequencing data is good?'
objectives:
  - 'Learn about QCing Illumina and Element data'
  - 'Learn about fastp'
  -  'Learn about multiqc'
time_estimation: 30M
key_points:
  - 'Galaxy is an excellent tool for quickly QCing your data'
subtopic: core
priority: 10
contributors:
  - nekrut
draft: true

---

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Introduction

In this tutorial we will look at assessing quality of data from two short read technologies: [Illumina](http://www.nature.com/doifinder/10.1038/nature07517) and [Element Biosciences](http://dx.doi.org/10.1038/s41587-023-01750-7). 

## FASTQ format

[FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) is not a very well defined format. In the beginning various manufacturers of sequencing instruments were free to interpret FASTQ as they saw fit, resulting in a multitude of FASTQ flavors. This variation stemmed primarily from different ways of encoding quality values as described [on the Wikipedia article for FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) (below you will find an explanation of quality scores and their meaning). Today, the [FASTQ Sanger](https://www.ncbi.nlm.nih.gov/pubmed/20015970) version of the format is considered to be the standard form of FASTQ. Galaxy is using FASTQ Sanger as the only legitimate input for downstream processing tools and provides [a number of utilities for converting FASTQ files](https://www.ncbi.nlm.nih.gov/pubmed/20562416) into this form (see **FASTQ Quality Control** section of Galaxy tools).

The FASTQ format looks like this:


```

@M02286:19:000000000-AA549:1:1101:12677:1273 1:N:0:23
CCTACGGGTGGCAGCAGTGAGGAATATTGGTCAATGGACGGAAGTCTGAACCAGCCAAGTAGCGTGCAG
+
ABC8C,:@F:CE8,B-,C,-6-9-C,CE9-CC--C-<-C++,,+;CE<,,CD,CEFC,@E9<FCFCF?9
@M02286:19:000000000-AA549:1:1101:15048:1299 1:N:0:23
CCTACGGGTGGCTGCAGTGAGGAATATTGGACAATGGTCGGAAGACTGATCCAGCCATGCCGCGTGCAG
+
ABC@CC77CFCEG;F9<F89<9--C,CE,--C-6C-,CE:++7:,CF<,CEF,CFGGD8FFCFCFEGCF
@M02286:19:000000000-AA549:1:1101:11116:1322 1:N:0:23
CCTACGGGAGGCAGCAGTAGGGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGT
+
AAC<CCF+@@>CC,C9,F9C9@9-CFFFE@7@:+CC8-C@:7,@EFE,6CF:+8F7EFEEF@EGGGEEE

```

Each sequencing read is represented by four lines:

1. `@` followed by read ID and optional information about sequencing run
2. sequenced bases
3. `+` (optionally followed by the read ID and some additional info)
4. quality scores for each base of the sequence encoded as [ASCII symbols](https://en.wikipedia.org/wiki/ASCII)

## Paired end data

It is common to prepare pair-end and mate-pair sequencing libraries. This is highly beneficial for a number of applications discussed in subsequent topics. For now let's just briefly discuss what these are and how they manifest themselves in FASTQ form.

![Paired-end and mate-pair reads](../../images/pe_mp.png "<b>Paired-end and mate-pair reads</b>. In paired end sequencing (left) the actual ends of rather short DNA molecules (less than 1kb) are determined, while for mate pair sequencing (right) the ends of long molecules are joined and prepared in special sequencing libraries. In these mate pair protocols, the ends of long, size-selected molecules are connected with an internal adapter sequence (i.e. linker, yellow) in a circularization reaction. The circular molecule is then processed using restriction enzymes or fragmentation. Fragments are enriched for the linker and outer library adapters are added around the two combined molecule ends. The internal adapter can then be used as a second priming site for an additional sequencing reaction in the same orientation or sequencing can be performed from the second adapter, from the reverse strand. (From <i>Understanding and improving high-throughput sequencing data production and analysis</i>, Ph.D. dissertation by <a href='https://ul.qucosa.de/api/qucosa%3A11231/attachment/ATT-0/'>Martin Kircher</a>)")


Thus in both cases (paired-end and mate-pair) a single physical piece of DNA (or RNA in the case of RNA-seq) is sequenced from two ends and so generates two reads. These can be represented as separate files (two FASTQ files with first and second reads) or a single file were reads for each end are interleaved. Here are examples:

### Two single files

File 1

```
 @M02286:19:000000000-AA549:1:1101:12677:1273 1:N:0:23
 CCTACGGGTGGCAGCAGTGAGGAATATTGGTCAATGGACGGAAGTCT
 +
 ABC8C,:@F:CE8,B-,C,-6-9-C,CE9-CC--C-<-C++,,+;CE
 @M02286:19:000000000-AA549:1:1101:15048:1299 1:N:0:23
 CCTACGGGTGGCTGCAGTGAGGAATATTGGACAATGGTCGGAAGACT
 +
 ABC@CC77CFCEG;F9<F89<9--C,CE,--C-6C-,CE:++7:,CF
```

File 2

```
@M02286:19:000000000-AA549:1:1101:12677:1273 2:N:0:23
CACTACCCGTGTATCTAATCCTGTTTGATACCCGCACCTTCGAGCTTA
+
--8A,CCE+,,;,<CC,,<CE@,CFD,,C,CFF+@+@CCEF,,,B+C,
@M02286:19:000000000-AA549:1:1101:15048:1299 2:N:0:23
CACTACCGGGGTATCTAATCCTGTTCGCTCCCCACGCTTTCGTCCATC
+
-6AC,EE@::CF7CFF<<FFGGDFFF,@FGGGG?F7FEGGGDEFF>FF
```

> <comment-title>Read order is important</comment-title>
> Note that read IDs are **identical** in two files and they are listed in **the same** order. In some cases read IDs in the first and second file may be appended with `/1` and `/2` tags, respectively.
{: .comment}

### Interleaved file

```
@1/1
AGGGATGTGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA
+
EGGEGGGDFGEEEAEECGDEGGFEEGEFGBEEDDECFEFDD@CDD<ED
@1/2
CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC
+
GHHHDFDFGFGEGFBGEGGEGEGGGHGFGHFHFHHHHHHHEF?EFEFF
@2/1
AGGGATGTGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA
+
HHHHHHEGFHEEFEEHEEHHGGEGGGGEFGFGGGGHHHHFBEEEEEFG
@2/2
CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC
+
HHHHHHHHHHHHHGHHHHHHGHHHHHHHHHHHFHHHFHHHHHHHHHHH
```

Here the first and the second reads are identified with `/1` and `/2` tags.

> <comment-title>FASTQ format is a loose standard</comment-title>
> FASTQ format is not strictly defined and its variations will always cause headache for you. See [this page](https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/) for more information.
{: .comment}

## What are base qualities?

As we've seen above, FASTQ datasets contain two types of information:

- *sequence of the read*
- *base qualities* for each nucleotide in the read.

The base qualities allow us to judge how trustworthy each base in a sequencing read is. The following excerpt from an excellent [tutorial](http://chagall.med.cornell.edu/RNASEQcourse/Intro2RNAseq.pdf) by Friederike Dündar, Luce Skrabanek, Paul Zumbo explains what base qualities are:

> <comment-title>From "Introduction to differential gene expression analysis using RNA-seq"</comment-title>
> Illumina sequencing is based on identifying the individual nucleotides by the fluorescence signal emitted upon their incorporation into the growing sequencing read. Once the fluorescence intensities are extracted and translated into the four letter code. The deduction of nucleotide sequences from the images acquired during sequencing is commonly referred to as base calling.
><br><br>
> Due to the imperfect nature of the sequencing process and limitations of the optical instruments, base calling will always have inherent uncertainty. This is the reason why FASTQ files store the DNA sequence of each read together with a position-specific quality score that represents the error probability, i.e., how likely it is that an individual base call may be incorrect. The score is called [Phred score](http://www.phrap.com/phred/), $$Q$$, which is proportional to the probability $$p$$ that a base call is incorrect, where $$Q = −10lg(p)$$. For example, a Phred score of 10 corresponds to one error in every ten base calls ($$Q = −10lg(0.1)$$), or 90% accuracy; a Phred score of 20 corresponds to one error in every 100 base calls, or 99% accuracy. A higher Phred score thus reflects higher confidence in the reported base.
><br><br>
> To assign each base a unique score identifier (instead of numbers of varying character length), Phred scores are typically represented as ASCII characters. At http://ascii-code.com/ you can see which characters are assigned to what number.
><br><br>
> For raw reads, the range of scores will depend on the sequencing technology and the base caller used (Illumina, for example, used a tool called Bustard, or, more recently, RTA). Unfortunately, Illumina has been anything but consistent in how they calculated and ASCII-encoded the Phred score (see below)! In addition, Illumina now allows Phred scores for base calls with as high as 45, while 41 used to be the maximum score until the HiSeq X. This may cause issues with downstream sapplications that expect an upper limit of 41.
{: .comment}

![Illumina quality score](../../images/illumina_qs.png "Comparison between different encodings")

Base call quality scores are represented with the Phred range. Different Illumina (formerly Solexa) versions
used different scores and ASCII offsets. Starting with Illumina format 1.8, the score now represents the standard
Sanger/Phred format that is also used by other sequencing platforms and the sequencing archives.

![FASTQ quality score](../../images/fastq_qs.png "The ASCII interpretation and ranges of the different Phred score notations used by Illumina and the original Sanger interpretation. Although the Sanger format allows a theoretical score of 93, raw sequencing reads typically do not exceed a Phred score of 60. In fact, most Illumina-based sequencing will result in maximum scores of 41 to 45 (image from <a href='https://en.wikipedia.org/wiki/FASTQ_format'>Wikipedia</a>).")


# Assessing quality

This tutorial provides two sample datasets: one generated with Illumina and the other with ElemBio (AVITI) machine. Pick any. We will compare qualities of the two platforms at the and of this tutorial. 

## Get data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}):
>
>    **Element** (AVITI). This datasets contains forward (r1) and reverse (r2) reads for two samples---`s1` and `s2`:
>
>    ```
>    https://zenodo.org/records/10870107/files/elem_s1_r1.fq.gz
>    https://zenodo.org/records/10870107/files/elem_s1_r2.fq.gz
>    https://zenodo.org/records/10870107/files/elem_s2_r1.fq.gz
>    https://zenodo.org/records/10870107/files/elem_s2_r2.fq.gz
>    ```
>
>    **Illumina**. This datasets contains forward (_1) and reverse (_2) reads for two samples---`bl` and `ch`:
>
>    ```
>    https://zenodo.org/record/5119008/files/M117-bl_1.fq.gz
>    https://zenodo.org/record/5119008/files/M117-bl_2.fq.gz
>    https://zenodo.org/record/5119008/files/M117-ch_1.fq.gz
>    https://zenodo.org/record/5119008/files/M117-ch_2.fq.gz
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Check that the datatype is set to `fastqsanger.gz` (Galaxy usually automatically assigns the right datatype, but check anyway).
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 4. Create a paired collection from your data (relevant part beging at 0:57 and and the case of this tutorial we have 4 datasets): 
>
>    <iframe width="560" height="315" src="https://www.youtube-nocookie.com/embed/uZUt9XIHUQo?si=F0qUj76L_lrM7R2j&amp;start=57" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
>
{: .hands_on}

## Assess quality and clean the reads with `fastp`

[`fastp`](https://github.com/OpenGene/fastp) is fast utility for computing read qualities and performing initial clean up of the data such as [adapter removal](https://github.com/OpenGene/fastp?tab=readme-ov-file#adapters).

> <hands-on-title> Running `fastp` </hands-on-title>
>
> 1. {% tool [fastp](toolshed.g2.bx.psu.edu/repos/iuc/fastp/fastp/0.23.2+galaxy0) %} with the following parameters:
>    - *"Single-end or paired reads"*: `Paired Collection`
>        - {% icon param-collection %} *"Select paired collection(s)"*: `output` (Input dataset collection)
>    - In *"Output Options"*:
>        - *"Output JSON report"*: `Yes`
>
{: .hands_on}

## Sub-step with `MultiQC`[]

[`multiqc`](https://multiqc.info/) is a report generation tool. It can process outputs on multiple tools including `fastp`:

> <hands-on-title> Generating quality report </hands-on-title>
>
> 1. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.11+galaxy1) %} with the following parameters:
>    - In *"Results"*:
>        - {% icon param-repeat %} *"Insert Results"*
>            - *"Which tool was used generate logs?"*: `fastp`
>                - {% icon param-file %} *"Output of fastp"*: `report_json` (output of **fastp** {% icon tool %})
>
{: .hands_on}

# Interpreting the results

The comparison between Element and Illumina QC plot is shown below:

![QC scores for Element and Illumina](../../images/qc_elem_ill.svg "Element (top) versus Illumina (bottom). Element base qualities are consistently over 40 for most bases within each read. Blue and black lines represent samples <tt>s1</tt> and <tt>s2</tt> in case of Element and <tt>bl</tt> and <tt>ch</tt> in case of Illumina.")

Both datasets are exceptionally good. However, the elemnt data has even higher quality.

# Conclusion

Now you can easily evaluate quality of your own data in Galaxy!
