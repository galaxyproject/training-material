---
layout: tutorial_hands_on

title: "NGS data logistics"
zenodo_link: ""
questions:
  - "How to manipulate and process NGS data"
objectives:
  - "Understand most common types of NGS-related datatypes"
  - "Learn about how Galaxy handles NGS data using Illumina data derived from patients infected with SARS-CoV-2"
time_estimation: "1H30M"
key_points:
  - "FASTQ Sanger version of the format is considered to be the standard form of FASTQ."
  - "Paired end data can be provided as two files or as an interleaved one."
  - "FastqQC is a tool allowing to check the quality of FASTQ datasets."
  - "The most common tools for mapping are Bowtie, BWA, BWA-MEM. You can use in-built genome to map against or upload one if it is missing."
  - "The standard format for storing aligned reads is SAM/BAM. The major toolsets to process these datasets are DeepTools, SAMtools, BAMtools and Picard."
  - "Data can be uploaded directly from a computer, from EBI SRA and from NCBI SRA, also using FTP or URL."
  - "One can retrieve NGS data from Sequence Read Archive"
  - "Galaxy can analyze massive amounts of data and make them suitable for secondary analysis"
subtopic: next-steps
contributors:
  - nekrut
  - mvdbeek
  - tnabtaf
  - blankenberg
---

In this section we will look at practical aspects of manipulation of next-generation sequencing data. We will start with the FASTQ format produced by most sequencing machines and will finish with the SAM/BAM format representing mapped reads. The cover image above shows a screen dump of a SAM dataset.

# Introduction to sequencing data

## FASTQ manipulation and quality control

[FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) is not a very well defined format. In the beginning various manufacturers of sequencing instruments were free to interpret FASTQ as they saw fit, resulting in a multitude of FASTQ flavors. This variation stemmed primarily from different ways of encoding quality values as described [here](https://en.wikipedia.org/wiki/FASTQ_format) (below you will find an explanation of quality scores and their meaning). Today, the [FASTQ Sanger](https://www.ncbi.nlm.nih.gov/pubmed/20015970) version of the format is considered to be the standard form of FASTQ. Galaxy is using FASTQ Sanger as the only legitimate input for downstream processing tools and provides [a number of utilities for converting FASTQ files](https://www.ncbi.nlm.nih.gov/pubmed/20562416) into this form (see **FASTQ Quality Control** section of Galaxy tools).

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

|   |
|----|
| ![Paired-end and mate-pair reads](../../images/pe_mp.png) |
|<small>**Paired-end and mate-pair reads**. In paired end sequencing (left) the actual ends of rather short DNA molecules (less than 1kb) are determined, while for mate pair sequencing (right) the ends of long molecules are joined and prepared in special sequencing libraries. In these mate pair protocols, the ends of long, size-selected molecules are connected with an internal adapter sequence (i.e. linker, yellow) in a circularization reaction. The circular molecule is then processed using restriction enzymes or fragmentation. Fragments are enriched for the linker and outer library adapters are added around the two combined molecule ends. The internal adapter can then be used as a second priming site for an additional sequencing reaction in the same orientation or sequencing can be performed from the second adapter, from the reverse strand. (From "Understanding and improving high-throughput sequencing data production and analysis", Ph.D. dissertation by [Martin Kircher](https://ul.qucosa.de/api/qucosa%3A11231/attachment/ATT-0/))</small>|


Thus in both cases (paired-end and mate-pair) a single physical piece of DNA (or RNA in the case of RNA-seq) is sequenced from two ends and so generates two reads. These can be represented as separate files (two FASTQ files with first and second reads) or a single file were reads for each end are interleaved. Here are examples:

#### Two single files

**File 1**

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

 **File 2**

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

#### Interleaved file

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

The base qualities allow us to judge how trustworthy each base in a sequencing read is. The following excerpt from an excellent [tutorial](http://chagall.med.cornell.edu/RNASEQcourse/Intro2RNAseq.pdf) by Friederike D&uuml;ndar, Luce Skrabanek, Paul Zumbo explains what base qualities are:

> <comment-title>From "Introduction to differential gene expression analysis using RNA-seq"</comment-title>
> Illumina sequencing is based on identifying the individual nucleotides by the fluorescence signal emitted upon their incorporation into the growing sequencing read. Once the fluorescence intensities are extracted and translated into the four letter code. The deduction of nucleotide sequences from the images acquired during sequencing is commonly referred to as base calling.
><br>
> Due to the imperfect nature of the sequencing process and limitations of the optical instruments, base calling will always have inherent uncertainty. This is the reason why FASTQ files store the DNA sequence of each read together with a position-specific quality score that represents the error probability, i.e., how likely it is that an individual base call may be incorrect. The score is called [Phred score](http://www.phrap.com/phred/), $$Q$$, which is proportional to the probability $$p$$ that a base call is incorrect, where $$Q = −10lg(p)$$. For example, a Phred score of 10 corresponds to one error in every ten base calls ($$Q = −10lg(0.1)$$), or 90% accuracy; a Phred score of 20 corresponds to one error in every 100 base calls, or 99% accuracy. A higher Phred score thus reflects higher confidence in the reported base.
><br>
> To assign each base a unique score identifier (instead of numbers of varying character length), Phred scores are typically represented as ASCII characters. At http://ascii-code.com/ you can see which characters are assigned to what number.
><br>
> For raw reads, the range of scores will depend on the sequencing technology and the base caller used (Illumina, for example, used a tool called Bustard, or, more recently, RTA). Unfortunately, Illumina has been anything but consistent in how they calculated and ASCII-encoded the Phred score (see below)! In addition, Illumina now allows Phred scores for base calls with as high as 45, while 41 used to be the maximum score until the HiSeq X. This may cause issues with downstream sapplications that expect an upper limit of 41.
{: .comment}

![Illumina quality score](../../images/illumina_qs.png)

Base call quality scores are represented with the Phred range. Different Illumina (formerly Solexa) versions
used different scores and ASCII offsets. Starting with Illumina format 1.8, the score now represents the standard
Sanger/Phred format that is also used by other sequencing platforms and the sequencing archives.

|                                                              |
|--------------------------------------------------------------|
| ![FASTQ quality score](../../images/fastq_qs.png) |
| <small>The ASCII interpretation and ranges of the different Phred score notations used by Illumina and the original Sanger interpretation. Although the Sanger format allows a theoretical score of 93, raw sequencing reads typically do not exceed a Phred score of 60. In fact, most Illumina-based sequencing will result in maximum scores of 41 to 45 (image from [Wikipedia](https://en.wikipedia.org/wiki/FASTQ_format))</small> |

## Assessing data quality

One of the first steps in the analysis of NGS data is seeing how good the data actually is. [FastqQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a fantastic tool allowing you to assess the quality of FASTQ datasets (and deciding whether to blame or not to blame whoever has done sequencing for you).

|                                        |                                    |
|:---------------------------------------|:-----------------------------------|
| ![Good quality in FastQC](../../images/good_fq.png)    | ![Bad quality in FastQC](../../images/bad_fq.png) |
|<small>**A.** Excellent quality</small> | <small>**B.** Hmmm...OK</small>    |

Here you can see FastQC base quality reports (the tools gives you many other types of data) for two datasets: **A** and **B**. The **A** dataset has long reads (250 bp) and very good quality profile with no qualities dropping below [phred score](http://www.phrap.com/phred/) of 30. The **B** dataset is significantly worse with ends of the reads dipping below phred score of 20. The **B** reads may need to be trimmed for further processing.

<!--
<div class="embed-responsive embed-responsive-16by9"><iframe src="https://player.vimeo.com/video/123453134?portrait=0" webkitallowfullscreen mozallowfullscreen allowfullscreen></iframe></div>
-->

## Mapping your data

Mapping of NGS reads against reference sequences is one of the key steps of the analysis. Now it is time to see how this is done in practice. Below is a list of key publications highlighting mainstream mapping tools:

- 2009 Bowtie 1 - [Langmead et al.](http://genomebiology.com/content/10/3/R25)
- 2012 Bowtie 2 - [Langmead and Salzberg](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3322381/)
- 2009 BWA - [Li and Durbin](https://academic.oup.com/bioinformatics/article/25/14/1754/225615/Fast-and-accurate-short-read-alignment-with)
- 2010 BWA - [Li and Durbin](https://academic.oup.com/bioinformatics/article/26/5/589/211735/Fast-and-accurate-long-read-alignment-with-Burrows)
- 2013 BWA-MEM - [Li](https://arxiv.org/abs/1303.3997)

### Mapping against a pre-computed genome index

Mappers usually compare reads against a reference sequence that has been transformed into a highly accessible data structure called genome index. Such indexes should be generated before mapping begins. Galaxy instances typically store indexes for a number of publicly available genome builds.

|                                                              |
|--------------------------------------------------------------|
| ![Cached genome](../../images/cached_genome.png)                    |
|<small>Mapping against a pre-computed index in Galaxy.</small>|

For example, the image above shows indexes for `hg38` version of the human genome. You can see that there are actually three choices: (1) `hg38`, (2) `hg38 canonical` and (3) `hg38 canonical female`. The `hg38` contains all chromosomes as well as all unplaced contigs. The `hg38 canonical` does not contain unplaced sequences and only consists of chromosomes 1 through 22, X, Y, and mitochondria. The
`hg38 canonical female` contains everything from the canonical set with the exception of chromosome Y.

### What if pre-computed index does not exist?

If Galaxy does not have a genome you need to map against, you can upload your genome sequence as a FASTA file and use it in the mapper directly as shown below (**Load reference genome** is set to `History`).

|                                                              |
|--------------------------------------------------------------|
| ![Uploaded genome](../../images/uploaded_genome.png) |
|<small>Mapping against a pre-computed index in Galaxy </small>|

In this case Galaxy will first create an index from this dataset and then run mapping analysis against it.

## SAM/BAM datasets

The [SAM/BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) format is an accepted standard for storing aligned reads (it can also store unaligned reads and some mappers such as BWA are accepting unaligned BAM as input). The binary form of the format (BAM) is compact and can be rapidly searched (if indexed). In Galaxy BAM datasets are always indexed (accompanies by a .bai file) and sorted in coordinate order. In the following duscussion I once again rely on [tutorial](http://chagall.med.cornell.edu/RNASEQcourse/Intro2RNAseq.pdf) by Friederike D&uuml;ndar, Luce Skrabanek, and Paul Zumbo.

The Sequence Alignment/Map (SAM) format is, in fact, a generic nucleotide alignment format that describes the alignment of sequencing reads (or query sequences) to a reference. The human readable, TABdelimited SAM files can be compressed into the Binary Alignment/Map format. These BAM files are bigger than simply gzipped SAM files, because they have been optimized for fast random access rather than size reduction. Position-sorted BAM files can be indexed so that all reads aligning to a locus can be efficiently retrieved without loading the entire file into memory.

As shown below, SAM files typically contain a short header section and a very long alignment section where each row represents a single read alignment. The following sections will explain the SAM format in a bit more detail. For the most comprehensive and updated information go to https://github.com/samtools/hts-specs.

|                                                              |
|--------------------------------------------------------------|
| ![BAM structure](../../images/bam_structure.png)   |
|<small>**Schematic representation of a SAM file**. Each line of the optional header section starts with “@”, followed by the appropriate abbreviation (e.g., SQ for sequence dictionary which lists all chromosomes names (SN) and their lengths (LN)). The vast majority of lines within a SAM file typically correspond to read alignments where each read is described by the 11 mandatory entries (black font) and a variable number of optional fields (grey font). From [tutorial](http://chagall.med.cornell.edu/RNASEQcourse/Intro2RNAseq.pdf) by Friederike D&uuml;ndar, Luce Skrabanek, and Paul Zumbo.</small>|

### SAM Header

The header section includes information about how the alignment was generated and stored. All lines in the header section are tab-delimited and begin with the “@” character, followed by tag:value pairs, where tag is a two-letter string that defines the content and the format of value. For example, the “@SQ” line in the header section contains the information about the names and lengths of the *reference sequences to which the reads were aligned. For a hypothetical organism with three chromosomes of length 1,000 bp, the SAM header should contain the following three lines:

```
@SQ SN:chr1 LN:1000
@SQ SN:chr2 LN:1000
@SQ SN:chr3 LN:1000
```

### SAM alignment section

The optional header section is followed by the alignment section where each line corresponds to one sequenced read. For each read, there are 11 mandatory fields that always appear in the same order:

```
<QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL>
```

If the corresponding information is unavailable or irrelevant, field values can be ‘0’ or ‘*’ (depending on the field, see below), but they cannot be missing! After the 11 mandatory fields, a variable number of optional fields can be present. Here’s an example of one single line of a real-life SAM file (you may need to scroll sideways):

```
ERR458493 .552967 16 chrI 140 255 12 M61232N37M2S * 0 0 CCACTCGTTCACCAGGGCCGGCGGGCTGATCACTTTATCGTGCATCTTGGC BB?HHJJIGHHJIGIIJJIJGIJIJJIIIGHBJJJJJJHHHHFFDDDA1+B NH:i:1 HI:i:1 AS:i:41 nM:i:2
```

The following table explains the format and content of each field. The `FLAG`, `CIGAR`, and the optional fields (marked in blue) are explained in more detail below. The number of optional fields can vary widely between different SAM files and even between reads within in the same file. The field types marked in blue are explained in more detail in the main text below.

![SAM fields](../../images/sam_fields.png)

### `FLAG` field

The FLAG field encodes various pieces of information about the individual read, which is particularly important for PE reads. It contains an integer that is generated from a sequence of bits (0, 1). This way, answers to multiple binary (Yes/No) questions can be compactly stored as a series of bits, where each of the single bits can be addressed and assigned separately.

The following table gives an overview of the different properties that can be encoded in the FLAG field. The developers of the SAM format and samtools tend to use the hexadecimal encoding as a means to refer to the different bits in their documentation. The value of the FLAG field in a given SAM file, however, will always be the decimal representation of the sum of the underlying binary values (as shown in Table below, row 2).

|                                                              |
|--------------------------------------------------------------|
| ![SAM flag](../../images/sam_flag.png) |
|<small>The `FLAG` field of SAM files stores information about the respective read alignment in one single decimal number. The decimal number is the sum of all the answers to the Yes/No questions associated with each binary bit. The hexadecimal representation is used to refer to the individual bits (questions). A bit is set if the corresponding state is true. For example, if a read is paired, `0x1` will be set, returning the decimal value of 1. Therefore, all `FLAG` values associated with paired reads must be uneven decimal numbers. Conversely, if the `0x1` bit is unset (= read is not paired), no assumptions can be made about `0x2`, `0x8`, `0x20`, `0x40` and `0x80` because they refer to paired reads. From [tutorial](http://chagall.med.cornell.edu/RNASEQcourse/Intro2RNAseq.pdf) by Friederike D&uuml;ndar, Luce Skrabanek, and Paul Zumbo</small>|

In a run with single reads, the flags you most commonly see are:

- 0: This read has been mapped to the forward strand. (None of the bit-wise flags have been set.)
- 4: The read is unmapped (`0x4` is set).
- 16: The read is mapped to the reverse strand (`0x10` is set)

(`0x100`, `0x200` and `0x400` are not used by most aligners/mappers, but could, in principle be set for single reads.) Some common `FLAG` values that you may see in a PE experiment include:


|                     |                                      |
----------------------|---------------------------------------
|**69** (= 1 + 4 + 64) | The read is paired, is the first read in the pair, and is unmapped.|
|**77** (= 1 + 4 + 8 + 64) | The read is paired, is the first read in the pair, both are unmapped.|
|**83** (= 1 + 2 + 16 + 64) | The read is paired, mapped in a proper pair, is the first read in the pair, and it is mapped to the reverse strand.|
|**99** (= 1 + 2 + 32 + 64) | The read is paired, mapped in a proper pair, is the first read in the pair, and its mate is mapped to the reverse strand.|
|**133** (= 1 + 4 + 128) | The read is paired, is the second read in the pair, and it is unmapped.|
|**137** (= 1 + 8 + 128) | The read is paired, is the second read in the pair, and it is mapped while its mate is not.|
|**141** (= 1 + 4 + 8 + 128) | The read is paired, is the second read in the pair, but both are unmapped.|
|**147** (= 1 + 2 + 16 + 128) | The read is paired, mapped in a proper pair, is the second read in the pair, and mapped to the reverse strand.|
|**163** (= 1 + 2 + 32 + 128) | The read is paired, mapped in a proper pair, is the second read in the pair, and its mate is mapped to the reverse strand.|

A useful website for quickly translating the FLAG integers into plain English explanations like the ones shown above is: https://broadinstitute.github.io/picard/explain-flags.html

### `CIGAR` string

`CIGAR` stands for *Concise Idiosyncratic Gapped Alignment Report*. This sixth field of a SAM file
contains a so-called CIGAR string indicating which operations were necessary to map the read to the reference sequence at that particular locus.

The following operations are defined in CIGAR format (also see figure below):

- **M** - Alignment (can be a sequence match or mismatch!)
- **I** - Insertion in the read compared to the reference
- **D** - Deletion in the read compared to the reference
- **N** - Skipped region from the reference. For mRNA-to-genome alignments, an N operation represents an intron. For other types of alignments, the interpretation of N is not defined.
- **S** - Soft clipping (clipped sequences are present in read); S may only have H operations between them and the ends of the string
- **H** - Hard clipping (clipped sequences are NOT present in the alignment record); can only be present as the first and/or last operation
- **P** - Padding (silent deletion from padded reference)
- **=** - Sequence match (not widely used)
- **X** - Sequence mismatch (not widely used)

The sum of lengths of the **M**, **I**, **S**, **=**, **X** operations must equal the length of the read. Here are some examples:

|                                 |
|---------------------------------|
|![CIGAR](../../images/cigar.png)|
|<small>From [tutorial](http://chagall.med.cornell.edu/RNASEQcourse/Intro2RNAseq.pdf) by Friederike D&uuml;ndar, Luce Skrabanek, and Paul Zumbo.</small>|

### Optional fields

Following the eleven mandatory SAM file fields, the optional fields are presented as key-value
pairs in the format of `<TAG>:<TYPE>:<VALUE>`, where `TYPE` is one of:

- `A` - Character
- `i` - Integer
- `f` - Float number
- `Z` - String
- `H` - Hex string

The information stored in these optional fields will vary widely depending on the mapper and new tags can be added freely. In addition, reads within the same SAM file may have different numbers of optional fields, depending on the program that generated the SAM file. Commonly used optional tags include:

- `AS:i` - Alignment score
- `BC:Z` - Barcode sequence
- `HI:i` - Match is i-th hit to the read
- `NH:i` - Number of reported alignments for the query sequence
- `NM:i` - Edit distance of the query to the reference
- `MD:Z` - String that contains the exact positions of mismatches (should complement the CIGAR string)
- `RG:Z` - Read group (should match the entry after ID if @RG is present in the header.

Thus, for example, we can use the NM:i:0 tag to select only those reads which map perfectly to the reference(i.e., have no mismatches). While the optional fields listed above are fairly standardized, tags that begin with `X`, `Y`, and `Z` are reserved for particularly free usage and will never be part of the official SAM file format specifications. `XS`, for example, is used by TopHat (an RNA-seq analysis tool we will discuss later) to encode the strand information (e.g., `XS:A:+`) while Bowtie2 and BWA use `XS:i:` for reads with multiple alignments to store the alignment score for the next-best-scoring alignment (e.g., `XS:i:30`).

### Read Groups

One of the key features of SAM/BAM format is the ability to label individual reads with readgroup tags. This allows pooling results of multiple experiments into a single BAM dataset. This significantly simplifies downstream logistics: instead of dealing with multiple datasets one can handle just one. Many downstream analysis tools such as variant callers are designed to recognize readgroup data and output results on per-readgroup basis.

One of the best descriptions of BAM readgroups is on [GATK support site](https://gatkforums.broadinstitute.org/discussion/1317/collected-faqs-about-bam-files). We have gratefully stolen two tables describing the most important readgroup tags - `ID`, `SM`, `LB`, and `PL` - from GATK forum and provide them here:

![Read groups](../../images/rg.png)

GATK forum also provides the following example:

![Read group example](../../images/rg_example.png)

### Manipulating SAM/BAM datasets

We support four major toolsets for processing of SAM/BAM datasets:

 * [DeepTools](https://deeptools.readthedocs.io) - a suite of user-friendly tools for the visualization, quality control and normalization of data from deep-sequencing DNA sequencing experiments.
 * [SAMtools](http://www.htslib.org/) - various utilities for manipulating alignments in the SAM/BAM format, including sorting, merging, indexing and generating alignments in a per-position format.
 * [BEDtools](https://bedtools.readthedocs.io/en/latest/) - a toolkit originally written for BED format was expanded for analysis of BAM and VCF datasets.
 * [Picard](https://broadinstitute.github.io/picard/) - a set of Java tools for manipulating high-throughput sequencing data (HTS) data and formats.

## The challenge of read duplicates

### PCR duplicates

Preparation of sequencing libraries (at least at the time of writing) for technologies such as Illumina (used in this example) involves PCR amplification. It is required to generate sufficient number of sequencing templates so that a reliable detection can be performed by base callers. Yet PCR has its biases, which are especially profound in cases of multitemplate PCR used for construction of sequencing libraries (Kanagawa et al. [2003](https://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&dopt=Abstract&list_uids=16233530)).

|                                              |
|----------------------------------------------|
| ![PCR duplicates](../../images/pcr-duplicates.png) |
|<small>Analyzing molecules aligning with the same outer coordinates, a mapping quality of at least 30 and a length of at least 30nt, resulted in an average coverage of 12.9 per PCR duplicate and an empirical coverage distribution similar to an exponential/power law distribution (left upper panel). This indicates that many molecules are only observed for deeper sequencing while other molecules are available at higher frequencies. Analyzing length (left middle panel) and GC content (left lower panel) patterns as well as the combination (right panel) shows higher PCR duplicate counts for a GC content between 30% to 70% as well as for shorter molecules compared to longer molecules. This effect may be due to an amplification bias from the polymerase or the cluster generation process necessary for Illumina sequencing. From Ph.D. dissertation of [Martin Kircher](https://ul.qucosa.de/api/qucosa%3A11231/attachment/ATT-0/)).</small>|

Duplicates can be identified based on their outer alignment coordinates or using sequence-based clustering. One of the common ways for identification of duplicate reads is the `MarkDuplicates` utility from [Picard](https://broadinstitute.github.io/picard/command-line-overview.html) package. It is designed to identify both PCR and optical duplicates:

<div class="well well-lg">

Duplicates are identified as read pairs having identical 5' positions (coordinate and strand) for both reads in a mate pair (and optionally, matching unique molecular identifier reads; see BARCODE_TAG option). Optical, or more broadly Sequencing, duplicates are duplicates that appear clustered together spatially during sequencing and can arise from optical/imagine-processing artifacts or from bio-chemical processes during clonal amplification and sequencing; they are identified using the READ_NAME_REGEX and the OPTICAL_DUPLICATE_PIXEL_DISTANCE options. The tool's main output is a new SAM or BAM file in which duplicates have been identified in the SAM flags field, or optionally removed (see REMOVE_DUPLICATE and REMOVE_SEQUENCING_DUPLICATES), and optionally marked with a duplicate type in the 'DT' optional attribute. In addition, it also outputs a metrics file containing the numbers of READ_PAIRS_EXAMINED, UNMAPPED_READS, UNPAIRED_READS, UNPAIRED_READ DUPLICATES, READ_PAIR_DUPLICATES, and READ_PAIR_OPTICAL_DUPLICATES. Usage example: java -jar picard.jar MarkDuplicates I=input.bam \ O=marked_duplicates.bam M=marked_dup_metrics.txt.`

</div>

### Sampling coincidence duplicates

However, one has to be careful when removing duplicates in cases when the sequencing targets are small (e.g., sequencing of bacterial, viral, or organellar genomes as well as amplicons). This is because when sequencing target is small reads will have the same coordinates by chance and not because of PCR amplification issues. The figure below illustrates the fine balance between estimates allele frequency, coverage, and variation in insert size:

|                                              |
|----------------------------------------------|
| ![Sampling bias](../../images/sampling-bias.png) |
| <small>The Variant Allele Frequency (VAF) bias determined by coverage and insert size variance. Reads are paired-end and read length is 76. The insert size distribution is modeled as a Gaussian distribution with mean at 200 and standard deviation shown on the x-axis. The true VAF is 0.05. The darkness at each position indicates the magnitude of the bias in the VAF. (From Zhou et al. [2013](https://bioinformatics.oxfordjournals.org/content/30/8/1073)).</small> |


# Getting NGS data to Galaxy

You can upload data in Galaxy using one of these ways:

## From your computer

This works well for small files because web browser do not like lengthy file transfers:

<iframe width="560" height="315" src="https://www.youtube.com/embed/FFCDx1rMGAQ" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

## Using FTP

FTP ([file transfer protocol](https://en.wikipedia.org/wiki/File_Transfer_Protocol)) allows transferring large collection of files:

<iframe width="560" height="315" src="https://www.youtube.com/embed/hC8KSuT_OP8" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

## From NCBI short read archive

Finally, datasets can be uploaded directly from NCBI's short read archive:

<iframe width="560" height="315" src="https://www.youtube.com/embed/Q4t-beYZ-do" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

> <comment-title></comment-title>
>
> We will use this last approach, getting data from NCBI SRA, in this tutorial.
>
{: .comment}

# Let's do it: From reads to variants

In primary analysis we start with raw sequencing data (e.g., fastq reads) and convert them into a dataset for secondary analysis. Such dataset can be a list of sequence variants, a collection of ChIP-seq peaks, a list of differentially expressed genes and so on.

SARS-CoV-2 is undoubtedly the most talked about biological system of the day. We will use SARS-CoV-2 sequencing data to demonstrate how Galaxy handles NGS data. In this case we will go from FASTQ data to a list of variants.

## Find necessary data in SRA

First we need to find a good dataset to play with. The [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) is the primary archive of *unassembled reads*  operated by the [US National Institutes of Health (NIH)](https://www.ncbi.nlm.nih.gov/).  SRA is a great place to get the sequencing data that underlie publications and studies. Let's do that:

> <hands-on-title>Get Metadata from NCBI SRA</hands-on-title>
>
> 1. Go to NCBI's SRA page by pointing your browser to https://www.ncbi.nlm.nih.gov/sra
> 2. In the search box enter `SARS-CoV-2 Patient Sequencing From Partners / MGH`:
> ![Find data](../../images/find_mgh_data.png) (Alternatively, you simply click on this [link](https://www.ncbi.nlm.nih.gov/sra/?term=SARS-CoV-2+Patient+Sequencing+From+Partners+%2F+MGH))
> 3. The web page will show a large number of SRA datasets (at the time of writing there were 2,223). This is data from a [study](https://science.sciencemag.org/content/early/2020/12/09/science.abe3261) describing analysis of SARS-CoV-2 in Boston area.
> 4. Download metadata describing these datasets by:
>   - clicking on **Send to:** dropdown
>   - Selecting `File`
>   - Changing **Format** to `RunInfo`
>   - Clicking **Create file**
> Here is how it should look like:
> ![Getting RunInfo from SRA](../../images/get_runinfo.png)
> 5. This would create a rather large `SraRunInfo.csv` file in your `Downloads` folder.
{: .hands_on}

Now that we have downloaded this file we can go to a Galaxy instance and start processing it.

> <comment-title></comment-title>
>
> Note that the file we just downloaded is **not** sequencing data itself. Rather, it is *metadata* describing properties of sequencing reads. We will filter this list down to just a few accessions that will be used in the remainder of this tutorial.
>
{: .comment}

## Process and filter `SraRunInfo.csv` file in Galaxy

> <hands-on-title>Upload `SraRunInfo.csv` file into Galaxy</hands-on-title>
>
> 1. Go to your Galaxy instance of choice such as one of the [usegalaxy.org](https://usegalaxy.org/), [usegalaxy.eu](https://usegalaxy.eu), [usegalaxy.org.au](https://usegalaxy.org.au) or any other. (This tutorial uses usegalaxy.org).
> 1. Click *Upload Data* button:
> ![Data upload button](../../images/upload_button.png)
> 1. In the dialog box that would appear click "*Choose local files*" button:
> ![Choose local files button](../../images/choose_local_files_button.png)
> 1. Find and select `SraRunInfo.csv` file from your computer
> 1. Click *Start* button
> 1. Close dialog by pressing **Close** button
> 1. You can now look at the content of this file by clicking {% icon galaxy-eye %} (eye) icon. You will see that this file contains a lot of information about individual SRA accessions. In this study every accession corresponds to an individual patient whose samples were sequenced.
{: .hands_on}

Galaxy can process all 2,000+ datasets, but to make this tutorial bearable we need to selected a smaller subset. In particular our previous experience with this data shows two interesting datasets `SRR11954102` and `SRR12733957`. So, let's pull them out.

{% snippet faqs/galaxy/analysis_cut.md %}

> <hands-on-title>Creating a subset of data</hands-on-title>
>
> 1. Find {% tool [Select lines that match an expression](Grep1) %} tool in **Filter and Sort** section of the tool panel.
>    > <tip-title>Finding tools</tip-title>
>    > Galaxy may have an overwhelming amount of tools installed. To find a specific tool type the tool name in the tool panel search box to find the tool.
>    {: .tip}
> 1. Make sure the `SraRunInfo.csv` dataset we just uploaded is listed in the {% icon param-file %} "*Select lines from*" field of the tool form.
> 1. In "*the pattern*" field enter the following expression &rarr; `SRR12733957|SRR11954102`. These are two accession we want to find separated by the pipe symbol `|`. The `|` means `or`: find lines containing `SRR12733957` **or** `SRR11954102`.
> 1. Click `Execute` button.
> 1. This will generate a file containing two lines (well ... one line is also used as the header, so it will appear the the file has three lines. It is OK.)
> 1. Cut the first column from the file using {% tool [Cut](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/1.1.0) %} tool, which you will find in **Text Manipulation** section of the tool pane.
> 1. Make sure the dataset produced by the previous step is selected in the "*File to cut*" field of the tool form.
> 1. Change "*Delimited by*" to `Comma`
> 1. In "*List of fields*" select `Column: 1`.
> 1. Hit `Execute`
> This will produce a text file with just two lines:
> ```
> SRR12733957
> SRR11954102
>```
{: .hands_on}

Now that we have identifiers of datasets we want we need to download the actual sequencing data. You can also watch [this video](#from-ncbi-short-read-archive).

## Download sequencing data

> <hands-on-title>Get data from SRA</hands-on-title>
>
> 1. Run {% tool [Faster Download and Extract Reads in FASTQ](toolshed.g2.bx.psu.edu/repos/iuc/sra_tools/fasterq_dump/2.10.9+galaxy0) %} with the following parameters:
>    - *"select input type"*: `List of SRA accession, one per line`
>        - The parameter {% icon param-file %} *"sra accession list"* should point the output of the {% icon tool %} "**Cut**" from the previous step.
>    - **Click** the `Execute` button. This will run the tool, which retrieves the sequence read datasets for the runs that were listed in the `SRA` dataset. It may take some time. So this may be a good time to take a break.
>
> 2. Several entries are created in your history panel when you submit this job:
>    - **`Pair-end data (fasterq-dump)`**: Contains Paired-end datasets (if available)
>    - **`Single-end data (fasterq-dump)`** Contains Single-end datasets (if available)
>    - **`Other data (fasterq-dump)`** Contains Unpaired datasets (if available)
>    - **`fasterq-dump log`** Contains Information about the tool execution
{: .hands_on}

The first three items are actually *collections* of datasets. *Collections* in Galaxy are logical groupings of datasets that reflect the semantic relationships between them in the experiment / analysis. In this case the tool creates separate collections for paired-end reads, single reads, and *other*. See the [Collections tutorial](https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/collections/tutorial.html) and watch [videos](https://youtube.com/playlist?list=PLNFLKDpdM3B9UaxWEXgziHXO3k-003FzE) (with names beginning with "Dataset Collections") for more information.

Explore the collections by first **clicking** on the collection name in the history panel. This takes you inside the collection and shows you the datasets in it.  You can then navigate back to the outer level of your history.

Once `fasterq` finishes transferring data (all boxes are green / done), we are ready to analyze it.


## Now what?

You can now analyze the retrieved data using any sequence analysis tools and workflows in Galaxy.  SRA holds backing data for every imaginable type of *-seq experiment.

If you ran this tutorial, but retrieved datasets that you were interested in, then see the rest of the GTN library for ideas on how to analyze in Galaxy.

However, if you retrieved the datasets used in this tutorial's examples above, then you are ready to run the SARS-CoV-2 variant analysis below.

In this part of the tutorial we will perform variant calling and basic analysis of the datasets downloaded above. We will start by downloading the Wuhan-Hu-1 SARS-CoV-2 reference sequence, then run adapter trimming, alignment and variant calling.

> <comment-title>The usegalaxy.* COVID-19 analysis project</comment-title>
> This tutorial uses a subset of the data and runs through the
> [Variation Analysis](https://covid19.galaxyproject.org/genomics/4-variation/)
> section of [covid19.galaxyproject.org](https://covid19.galaxyproject.org/).
> The data for [covid19.galaxyproject.org](https://covid19.galaxyproject.org/) is
> being updated continuously as new datasets are made public.
{: .comment}

## Get the reference genome data

The reference genome data for today is for SARS-CoV-2, "Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome", having the accession ID of NC_045512.2.

This data is available from directly from GenBank using the following [link](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz).

> <hands-on-title>Get the reference genome data</hands-on-title>
>
> 1. Import the following file into your history:
>
>    ```
>    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
{: .hands_on}

## Adapter trimming with **fastp**

Removing sequencing adapters improves alignments and variant calling. **fastp** {% icon tool %} can automatically detect widely used sequencing adapters.

> <hands-on-title>Running `fastp`</hands-on-title>
>
> Run {% tool [fastp](toolshed.g2.bx.psu.edu/repos/iuc/fastp/fastp/0.20.1+galaxy0) %} with the following parameters:
>    - *"Single-end or paired reads"*: `Paired Collection`
>        - {% icon param-file %} *"Select paired collection(s)"*: `list_paired` (output of **Faster Download and Extract Reads in FASTQ** {% icon tool %})
>    - In *"Output Options"*:
>        - *"Output JSON report"*: `Yes`
{: .hands_on}


## Alignment with  **Map with BWA-MEM**

**BWA-MEM** {% icon tool %} is a widely used sequence aligner for short-read sequencing datasets such as those we are analysing in this tutorial.

> <hands-on-title>Map sequencing reads to reference genome</hands-on-title>
>
> Run {% tool [BWA-MEM](toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.17.1) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `output` (Input dataset)
>    - *"Single or Paired-end reads"*: `Paired Collection`
>        - {% icon param-file %} *"Select a paired collection"*: `output_paired_coll` (output of **fastp** {% icon tool %})
>    - *"Set read groups information?"*: `Do not set`
>    - *"Select analysis mode"*: `1.Simple Illumina mode`
>
{: .hands_on}

## Remove duplicates with **MarkDuplicates**

**MarkDuplicates** {% icon tool %} removes duplicate sequences originating from library preparation artifacts and sequencing artifacts. It is important to remove these artefactual sequences to avoid artificial overrepresentation of single molecule.

> <hands-on-title>Remove duplicates</hands-on-title>
>
> Run {% tool [MarkDuplicates](toolshed.g2.bx.psu.edu/repos/devteam/picard/picard_MarkDuplicates/2.18.2.2)%} with the following parameters:
>    - {% icon param-file %} *"Select SAM/BAM dataset or dataset collection"*: `bam_output` (output of **Map with BWA-MEM** {% icon tool %})
>    - *"If true do not write duplicates to the output file instead of writing them with appropriate flags set"*: `Yes`
>
{: .hands_on}

<!-- Not used in this version since this step was not included in video

## Generate alignment statistics with **Samtools stats**

After the duplicate marking step above we can generate statistic about the alignment we have generated.

> <hands-on-title>Generate alignment statistics</hands-on-title>
>
> 1. **Samtools stats** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"BAM file"*: `outFile` (output of **MarkDuplicates** {% icon tool %})
>    - *"Set coverage distribution"*: `No`
>    - *"Output"*: `One single summary file`
>    - *"Filter by SAM flags"*: `Do not filter`
>    - *"Use a reference sequence"*: `No`
>    - *"Filter by regions"*: `No`
{: .hands_on}

-->

## **Realign reads** with lofreq viterbi

**Realign reads** {% icon tool %} corrects misalignments around insertions and deletions. This is required in order to accurately detect variants.

> <hands-on-title>Realign reads around indels</hands-on-title>
>
> Run {% tool [Realign reads](toolshed.g2.bx.psu.edu/repos/iuc/lofreq_viterbi/lofreq_viterbi/2.1.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Reads to realign"*: `outFile` (output of **MarkDuplicates** {% icon tool %})
>    - *"Choose the source for the reference genome"*: `History`
>        - {% icon param-file %} *"Reference"*: `output` (Input dataset)
>    - In *"Advanced options"*:
>        - *"How to handle base qualities of 2?"*: `Keep unchanged`
{: .hands_on}

## Add indel qualities with lofreq **Insert indel qualities**

This step adds indel qualities into our alignment file. This is necessary in order to call variants using **Call variants** with lofreq {% icon tool %}

> <hands-on-title>Add indel qualities</hands-on-title>
>
> Run {% tool [Insert indel qualities](toolshed.g2.bx.psu.edu/repos/iuc/lofreq_indelqual/lofreq_indelqual/2.1.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Reads"*: `realigned` (output of **Realign reads** {% icon tool %})
>    - *"Indel calculation approach"*: `Dindel`
>        - *"Choose the source for the reference genome"*: `History`
>            - {% icon param-file %} *"Reference"*: `output` (Input dataset)
>
{: .hands_on}

## Call Variants using lofreq **Call variants**

We are now ready to call variants.

> <hands-on-title>Call variants</hands-on-title>
>
> Run {% tool [Call variants](toolshed.g2.bx.psu.edu/repos/iuc/lofreq_call/lofreq_call/2.1.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input reads in BAM format"*: `output` (output of **Insert indel qualities** {% icon tool %})
>    - *"Choose the source for the reference genome"*: `History`
>        - {% icon param-file %} *"Reference"*: `output` (Input dataset)
>    - *"Call variants across"*: `Whole reference`
>    - *"Types of variants to call"*: `SNVs and indels`
>    - *"Variant calling parameters"*: `Configure settings`
>        - In *"Coverage"*:
>            - *"Minimal coverage"*: `50`
>        - In *"Base-calling quality"*:
>            - *"Minimum baseQ"*: `30`
>            - *"Minimum baseQ for alternate bases"*: `30`
>        - In *"Mapping quality"*:
>            - *"Minimum mapping quality"*: `20`
>    - *"Variant filter parameters"*: `Preset filtering on QUAL score + coverage + strand bias (lofreq call default)`
{: .hands_on}

The output of this step is a collection of VCF files that can be visualized in a genome browser.

## Annotate variant effects with **SnpEff eff:**

We will now annotate the variants we called in the previous step with the effect they have on the SARS-CoV-2 genome.

> <hands-on-title>Annotate variant effects</hands-on-title>
>
> Run {% tool [SnpEff](toolshed.g2.bx.psu.edu/repos/iuc/snpeff_sars_cov_2/snpeff_sars_cov_2/4.5covid19) %} with the following parameters:
>    - {% icon param-file %} *"Sequence changes (SNPs, MNPs, InDels)"*: `variants` (output of **Call variants** {% icon tool %})
>    - *"Output format"*: `VCF (only if input is VCF)`
>    - *"Create CSV report, useful for downstream analysis (-csvStats)"*: `Yes`
>    - *"Annotation options"*: ``
>    - *"Filter output"*: ``
>    - *"Filter out specific Effects"*: `No`
>
{: .hands_on}

The output of this step is a VCF file with added variant effects.

## Create table of variants using **SnpSift Extract Fields**

We will now select various effects from the VCF and create a tabular file that is easier to understand for humans.

> <hands-on-title>Create table of variants</hands-on-title>
>
> Run {% tool [SnpSift Extract Fields](toolshed.g2.bx.psu.edu/repos/iuc/snpsift/snpSift_extractFields/4.3+t.galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Variant input file in VCF format"*: `snpeff_output` (output of **SnpEff eff:** {% icon tool %})
>    - *"Fields to extract"*: `CHROM POS REF ALT QUAL DP AF SB DP4 EFF[*].IMPACT EFF[*].FUNCLASS EFF[*].EFFECT EFF[*].GENE EFF[*].CODON`
>    - *"One effect per line"*: `Yes`
>    - *"empty field text"*: `.`
>
{: .hands_on}

We can inspect the output files and see check if Variants in this file are also described in [an observable notebook that shows the geographic distribution of SARS-CoV-2 variant sequences](https://observablehq.com/@spond/distribution-of-sars-cov-2-sequences-that-have-a-particular)

Interesting variants include the C to T variant at position 14408 (14408C/T) in SRR11772204, 28144T/C in SRR11597145 and 25563G/T in SRR11667145.


## Summarize data with **MultiQC**

We will now summarize our analysis with MultiQC, which generates a beautiful report for our data.

> <hands-on-title>Summarize data</hands-on-title>
>
> Run {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.8+galaxy1) %} with the following parameters:
>    - In *"Results"*:
>        - {% icon param-repeat %} *"Insert Results"*
>            - *"Which tool was used generate logs?"*: `fastp`
>                - {% icon param-file %} *"Output of fastp"*: `report_json` (output of **fastp** {% icon tool %})
>        - {% icon param-repeat %} *"Insert Results"*
>            - *"Which tool was used generate logs?"*: `Picard`
>                - In *"Picard output"*:
>                    - {% icon param-repeat %} *"Insert Picard output"*
>                        - *"Type of Picard output?"*: `Markdups`
>                        - {% icon param-file %} *"Picard output"*: `metrics_file` (output of **MarkDuplicates** {% icon tool %})
{: .hands_on}

The above state allows us to judge the quality of the data. In this particular case data is not bad as quality values never drop below 30:

![MultiQC plot](../../images/multiqc.png)


## Collapse data into a single dataset

We now extracted meaningful fields from VCF datasets. But they still exist as a collection. To move towards secondary analysis we need to **collapse** this collection into a single dataset. For more information about collapsing collections see this video:

<iframe width="560" height="315" src="https://www.youtube.com/embed/ypuFZ1RKMIY" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

> <hands-on-title>Collapse a collection</hands-on-title>
>
> Run {% tool [Collapse Collection](toolshed.g2.bx.psu.edu/repos/nml/collapse_collections/collapse_dataset/4.0) %} with the following parameters:
>    - {% icon param-collection %} *"Collection of files to collapse into single dataset"*: `snpsift extract fields` (output of **SnpSift Extract Fields** {% icon tool %})
>    - "*Keep one header line*": `Yes`
>    - "*Prepend File name*": `Yes`
>    - "*Where to add dataset name*": `Same line and each line in dataset`
{: .hands_on}

You can see that this tool takes lines from all collection elements (in our case we have two), add element name as the first column, and pastes everything together. So if we have a collection as an input:

> <code-in-title>A collection with two items</code-in-title>
> A collection element named `SRR11954102`
>
>```
>NC_045512.2  84 PASS  C T  960.0  28  1.0       0 0,0,15,13 MODIFIER  NONE  INTERGENIC
>NC_045512.2 241 PASS  C T 2394.0  69  0.971014  0 0,0,39,29 MODIFIER  NONE  INTERGENIC
>```
>
>A collection element named `SRR12733957`:
>
>```
>NC_045512.2 241 PASS  C T 1954.0  63  0.888889  0 0,0,42,21 MODIFIER  NONE  INTERGENIC
>NC_045512.2 823 PASS  C T 1199.0  50  0.76      3 5,6,13,26 LOW       LOW   LOW
>```
{: .code-in}

We will have a single dataset as the output:

> <code-out-title>A single dataset</code-out-title>
>
>then the **Collapse Collection** {% icon tool %} will produce this:
>
>```
>SRR11954102 NC_045512.2  84 PASS  C T  960.0  28  1.0       0 0,0,15,13 MODIFIER  NONE  INTERGENIC
>SRR11954102 NC_045512.2 241 PASS  C T 2394.0  69  0.971014  0 0,0,39,29 MODIFIER  NONE  INTERGENIC
>SRR12733957 NC_045512.2 241 PASS  C T 1954.0  63  0.888889  0 0,0,42,21 MODIFIER  NONE  INTERGENIC
>SRR12733957 NC_045512.2 823 PASS  C T 1199.0  50  0.76      3 5,6,13,26 LOW       LOW   LOW
>```
{: .code-out}

you can see that added a column with dataset ID taken from collection element name.

## Anything interesting?

These data are now ready for downstream analysis. One of the interesting things about these data is that it contains some of the mutations identified in the [B.1.1.7](https://cov-lineages.org/global_report_B.1.1.7.html) lineage. For example, dataset `SRR12733957` contains a nonsense mutation (change for legitimate codon specifying an amino acid to a stop codon) in ORF8a:

```
SRR12733957 NC_045512.2 27972 PASS  C T 56.0  39  0.076923  0 7,29,0,3  HIGH  NONSENSE  STOP_GAINED ORF8  Caa/Taa
```

This is interesting because these datasets were collected well before B.1.1.7 became widely spread. Can you find more mutations here in these data?


<!--

# Let's do it: A peek at secondary analysis in Galaxy

After a dataset is converted into an intermediary result (such as variant table generated during the previous step) by primary analysis, this intermediary data needs to be processed further to produce biological insight. IN Galaxy this can be done by starting a Jupyter notebook directly on history datasets. Here is a short demo.

> <warning-title>Jupyter is in Beta</warning-title>
> Jupyter notebook functionality in Galaxy is currently in Beta. It means that we are still working on it and its implementation is constantly changing. Be aware of this fact.
{: .warning}

> <hands-on-title>Analyze data in Jupyter from Galaxy</hands-on-title>
>
> 1. Open {% icon tool %} **Interactive Jupyter Notebook** from section **Interactive Tools**
> 1. *"Do you already have a notebook?"*: `Start with a fresh notebook`
> 1. *"Include data into the environment"*: Select history dataset containing the output of the collection collapse step
> 1. Click `Execute`
> 1. It will create a dataset in the history that will be permanently yellow (running)
> 1. You will see the following banner:
>
>    ![Jupyter startup banner](../../images/jupyter_start.png)
>
> 1. Click on `User menu`. You may need to wait a bit (go get coffee quickly). In the end you will a `Jupyter Interactive Tool` link. Click on it.
> 1. This will open a new tab with a fully functional Jupyter environment:
> ![Jupyter running](../../images/jupyter_running.png)
> 1. Start a new notebook by clicking on Python 3 image. This will start a fresh notebook.
> 1. Type `!ls data/` in the first cell and you will see a file from Galaxy history.
> 1. Upgrade [Pandas](https://pandas.pydata.org/) and [Seaborn](https://seaborn.pydata.org/index.html): `!pip install -U pandas seaborn`
> 1. Read data into a dataframe:
>
>    ```
>    import pandas as pd
>    df = pd.read_csv('data/Collapsed_collection',sep='\t')
>    # Drop duplicates is required to get rid of SnpEff artifacts
>    df = df.drop_duplicates(ignore_index=True)
>    ```
>    {% icon warning %} Note that your file may be named differently! Use it exact name instead of `Collapsed_collection` above.
>
> 1. Filter variants with intermediate frequencies between 20% and 80%:
>
>    ```
>    df = df[ (df['AF']>=0.05) & (df['AF']<=0.8) ]
>    ```
>
> 1. Plot relationship between variant position and alternative allele frequency:
>
>    ```
>    import seaborn as sns
>    %matplotlib inline
>    sns.relplot(x='POS',y='AF',hue='EFF[*].FUNCLASS',data=df,height=5,aspect=2,s=100,style='Sample')
>    ```
>
>    This plot shows how variants are distributed across the SARS-CoV-2 genome:
>    ![Seaborn plot](../../images/sns_plot.png)
> 1. Play with it!
{: .hands_on}


# Let's do it: A peek at secondary analysis in Google Colab

We can also import data directly into Google Colab:

> <hands-on-title>Analyze data in Google Colab</hands-on-title>
>
> 1. Expand history dataset containing the output of the collection collapse step
> 1. Right click on {% icon galaxy-save %} (disk) icon and copy link to clipboard
> 1. Go to https://colab.research.google.com/ and create a new notebook
> 1. Type `!wget <link>` where you must replace `<link>` with the URL you copied into clipboard at the previous step
> 1. This will copy this dataset into your Colaboratory environment
> 1. See the dataset file by typing `!ls` in the Colab cell
> 1. You will see a file with a name that looks something like this: `'display?to_ext=tabular'`
> 1. Rename it into something palatable:
>```
>!mv 'display?to_ext=tabular' data.tsv`
>```
> 1. Read data into a [Pandas](https://pandas.pydata.org/) dataframe:
> ```
> import pandas as pd
> df = pd.read_csv('data.tsv',sep='\t')
> # Drop duplicates is required to get rid of SnpEff artifacts
> df = df.drop_duplicates(ignore_index=True)
>```
> 1. Filter variants with intermediate frequencies between 20% and 80%:
>```
>df = df[ (df['AF']>=0.05) & (df['AF']<=0.8) ]
>```
> 1. Plot relationship between variant position and alternative allele frequency:
>```
>import seaborn as sns
>sns.relplot(x='POS',y='AF',hue='EFF[*].FUNCLASS',data=df,height=5,aspect=2,s=100,style='Sample')
>```
> This plot shows how variants are distributed across the SARS-CoV-2 genome:
> ![Seaborn plot](../../images/sns_plot.png)
> 1. Play with it!
{: .hands_on}

-->

# Conclusion


Congratulations, you now know how to import sequence data from the SRA and how to run an example analysis on these datasets.
