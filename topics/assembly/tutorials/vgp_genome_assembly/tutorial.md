---
layout: tutorial_hands_on
title: VGP assembly pipeline
zenodo_link: 'https://zenodo.org/record/5887339'
level: Intermediate
tags:
 - pacbio
 - eukaryote
 - VGP
questions:
- "what combination of tools can produce the highest quality assembly of vertebrate genomes?"
- "How can we evaluate how good it is?"
objectives:
- "Learn the tools necessary to perform a de novo assembly of a vertebrate genome"
- "Evaluate the quality of the assembly"
time_estimation: '3h'
key_points:
- "The VGP pipeline allows to generate error-free, near gapless reference-quality genome assemblies"
- "The assembly can be divided in four main stages: genome profile analysis, HiFi long read phased assembly with hifiasm, Bionano hybrid scaffolding and Hi-C hybrid scaffolding"
contributors:
- delphine-l
- astrovsky01
- gallardoalba
- annasyme
- abueg
- pickettbd
- gf777
- msozzoni
---


# Introduction


Advances in sequencing technologies over the last few decades have revolutionized the field of genomics, allowing for a reduction in both the time and resources required to *de novo* genome assembly. Until recently, second-generation sequencing technologies (also known as Next Generation Sequencing or NGS) allowed to produce highly accurate but short (up to 800bp) reads, whose extension was not long enough to cope with the difficulties associated with repetitive regions. Today, so-called third-generation sequencing (TGS) technologies, usually known as single-molecule real-time (SMRT) sequencing, have become dominant in *de novo* assembly of large genomes. TGS can use native DNA without amplification, reducing sequencing error and bias ({% cite Hon2020 %}, {% cite Giani2020 %}). Very recently, Pacific Biosciences introduced High-Fidelity (HiFi) sequencing, which produces reads 10-25 kpb in length with a minimum accuracy of 99% (Q20). In this tutorial you will use HiFi reads in combination with data from additional sequencing technologies to generate a high-quality genome assembly.

Deciphering the structural organization of complex vertebrate genomes is currently one of the largest challenges in genomics ({% cite Frenkel2012 %}). Despite the significant progress made in recent years, a key question remains: what combination of data and tools can produce the highest quality assembly? In order to adequately answer this question, it is necessary to analyse two of the main factors that determine the difficulty of genome assembly processes: repetitive content and heterozygosity.

Repetitive elements can be grouped into two categories: interspersed repeats, such as transposable elements (TE) that occur at multiple loci throughout the genome, and tandem repeats (TR) that occur at a single locus ({% cite Trresen2019 %}). Repetitive elements are an important component of eukaryotic genomes, constituting over a third of the genome in the case of mammals ({% cite SoteroCaio2017 %}, {% cite Chalopin2015 %}). In the case of tandem repeats, various estimates suggest that they are present in at least one third of human protein sequences ({% cite Marcotte1999 %}). TE content is among the main factors contributing to the lack of continuity in the reconstruction of genomes, especially in the case of large ones, as TE content is highly correlated with genome size ({% cite SoteroCaio2017 %}). On the other hand, TR usually lead to local genome assembly collapse, especially when their length is close to that of the reads ({% cite Trresen2019 %}).

Heterozygosity is also an important factor in genome assembly. Haplotype phasing, that is, the identification of alleles that are co-located on the same chromosome, has become a fundamental problem in heterozygous and polyploid genome assemblies ({% cite Zhang2020 %}). When no reference sequence is available, the *state-of-the-art* strategy consists of constructing a string graph with vertices representing reads and edges representing consistent overlaps. In this kind of graph, after transitive reduction, heterozygous alleles in the string graph are represented by bubbles. When combined with Hi-C data, this approach allows complete diploid reconstruction ({% cite DominguezDelAngel2018 %}, {% cite Zhang2020 %}, {% cite Dida2021 %}).

The G10K consortium launched the Vertebrate Genomes Project (VGP), whose goal is generating high-quality, near-error-free, gap-free, chromosome-level, haplotype-phased, annotated reference genome assemblies for every vertebrate species ({% cite Rhie2021 %}). This tutorial will guide you step by step to assemble a high-quality genome using the VGP assembly pipeline.



> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# VGP assembly pipeline overview

This training is an adaptation of the VGP assembly pipeline 2.0 (fig. 1).

![Figure 1: VGP pipeline](../../images/vgp_assembly/VGP_Pipeline.png "VPG Pipeline 2.0. The pipeline starts with assembly of the HiFi reads into contigs, yielding the primary and alternate assemblies. Then, duplicated and erroneously assigned contigs will be removed by using purge_dups. Finally, Bionano optical maps and HiC data are used to generate a scaffolded primary assembly.")

With the aim of making it easier to understand, the training has been organized into four main sections: genome profile analysis, HiFi phased assembly with hifiasm, post-assembly processing and hybrid scaffolding.    
    
# Get data

In order to reduce computation time, we will assemble samples from the yeast _Saccharomyces cerevisiae_ S288C, a widely used laboratory strain isolated in the 1950s by Robert Mortimer. Using _S. cerevisae_, one of the most intensively studied eukaryotic model organisms, has the additional advantage of allowing us to evaluate the final result of our assembly with great precision. For this tutorial, we generated a set of synthetic HiFi reads corresponding to a theoretical diploid genome.

> <details-title>Generation of synthetic reads</details-title>
> The synthetic HiFi reads were generated by using the [S288C assembly](https://www.ncbi.nlm.nih.gov/genome/15?genome_assembly_id=22535) as reference genome. With this objective we used [HIsim](https://github.com/thegenemyers/HI.SIM), a HiFi shotgun simulator. The commands used are detailed below:
>
> ```bash
> ./FastK -T8 -Nsynthetic -k40 -t4 ./hifi_reads.fastq.gz
> ./FastK -T8 -Nsynthetic -k40 -p:synthetic.ktab .hifi_reads.fastq.gz
> ./Histex -G synthetic > synthetic.histogram
> ./GeneScopeFK.R -i synthetic.histogram -o ./synthetic -k 40 -p 1
> ./HImodel -v -osynthetic -g10:55 -e5 -T8 ./synthetic
> ./HIsim ./GCF_000146045.2_R64_genomic.fasta ./synthetic.model -p.2,.2 -c30 -r1 -U -f -h -e
> ```
> The selected mutation rate was 2%. Note that HIsim generates the synthetic reads in FASTA format. This is perfectly fine for illustrating the workflow, but you should be aware that usually you will work with HiFi reads in FASTQ format.
>
{: .details}
    
The first step is to get the datasets from Zenodo. The VGP assembly pipeline uses data generated by a variety of technologies, including PacBio HiFi reads, Bionano optical maps, and Hi-C chromatin interaction maps.
    
> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }})
>
>    - Open the file {% icon galaxy-upload %} __upload__ menu
>    - Click on **Rule-based** tab
>    - *"Upload data as"*: `Datasets`
>    - Copy the tabular data, paste it into the textbox and press <kbd>Build</kbd>
>
>       ```
>   Hi-C_dataset_F   https://zenodo.org/record/5550653/files/SRR7126301_1.fastq.gz?download=1   fastqsanger.gz    Hi-C
>   Hi-C_dataset_R   https://zenodo.org/record/5550653/files/SRR7126301_2.fastq.gz?download=1   fastqsanger.gz    Hi-C
>   Bionano_dataset    https://zenodo.org/record/5550653/files/bionano.cmap?download=1   cmap    Bionano
>       ```
>
>    - From **Rules** menu select `Add / Modify Column Definitions`
>       - Click `Add Definition` button and select `Name`: column `A`
>       - Click `Add Definition` button and select `URL`: column `B`
>       - Click `Add Definition` button and select `Type`: column `C`
>       - Click `Add Definition` button and select `Name Tag`: column `D`
>    - Click `Apply` and press <kbd>Upload</kbd>
>   
> 3. Import the remaining datasets from [Zenodo]({{ page.zenodo_link }})
>
>    - Open the file {% icon galaxy-upload %} __upload__ menu
>    - Click on **Rule-based** tab
>    - *"Upload data as"*: `Collections`
>    - Copy the tabular data, paste it into the textbox and press <kbd>Build</kbd>
>
>       ```
>   dataset_01    https://zenodo.org/record/6098306/files/HiFi_synthetic_50x_01.fasta?download=1  fasta    HiFi  HiFi_collection
>   dataset_02    https://zenodo.org/record/6098306/files/HiFi_synthetic_50x_02.fasta?download=1  fasta    HiFi  HiFi_collection
>   dataset_03    https://zenodo.org/record/6098306/files/HiFi_synthetic_50x_03.fasta?download=1  fasta    HiFi  HiFi_collection
>       ```
>
>    - From **Rules** menu select `Add / Modify Column Definitions`
>       - Click `Add Definition` button and select `List Identifier(s)`: column `A`
>       - Click `Add Definition` button and select `URL`: column `B`
>       - Click `Add Definition` button and select `Type`: column `C`
>       - Click `Add Definition` button and select `Group Tag`: column `D`
>       - Click `Add Definition` button and select `Collection Name`: column `E`
>    - Click `Apply` and press <kbd>Upload</kbd>
>
{: .hands_on}

### HiFi reads preprocessing with **Cutadapt**
    
Adapter trimming usually means trimming the adapter sequence off the ends of reads, which is where the adapter sequence is usually located in NGS reads. However, due to the nature of SMRT sequencing technology, adapters do not have a specific, predictable location in  HiFi reads. Additionally, the reads containing adapter sequence could be of generally lower quality compared to the rest of the reads. Thus, we will use Cutadapt not to trim, but to remove the entire read if a read is found to have an adapter inside of it.

> <comment-title>Background on PacBio HiFi reads</comment-title>
>
> PacBio HiFi reads rely on the Single Molecule Real-Time (SMRT) sequencing technology. SMRT is based on real-time imaging of fluorescently tagged nucleotides as they are added to a newly synthesized DNA strand. HiFi further combine multiple subreads from the same circular template to produce one highly accurate consensus sequence (fig. 2).
>
> ![Figure 2: PacBio HiFi technology](../../images/vgp_assembly/pacbio_technology.png "HiFi reads are produced by calling consensus from subreads generated by multiple passes of the enzyme around a circularized template. This results in a HiFi read that is both long and accurate.")
>
> This technology allows to generate long-read sequencing data with read lengths in the range of 10-25 kb and minimum read consensus accuracy  greater than 99% (Q20).
>
{: .comment}

> <hands-on-title>Primer removal with Cutadapt</hands-on-title>
>
> 1. {% tool [Cutadapt](toolshed.g2.bx.psu.edu/repos/lparsons/cutadapt/cutadapt/3.4) %} with the following parameters:
>    - *"Single-end or Paired-end reads?"*: `Single-end`
>        - {% icon param-collection %} *"FASTQ/A file"*: `HiFi_collection`
>        - In *"Read 1 Options"*:
>            - In *"5' or 3' (Anywhere) Adapters"*:
>                - {% icon param-repeat %} *"Insert 5' or 3' (Anywhere) Adapters"*
>                    - *"Source"*: `Enter custom sequence`
>                        - *"Enter custom 5' or 3' adapter name"*: `First adapter`
>                        - *"Enter custom 5' or 3' adapter sequence"*: `ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT`
>                - {% icon param-repeat %} *"Insert 5' or 3' (Anywhere) Adapters"*
>                    - *"Source"*: `Enter custom sequence`
>                        - *"Enter custom 5' or 3' adapter name"*: `Second adapter`
>                        - *"Enter custom 5' or 3' adapter sequence"*: `ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT`
>    - In *"Adapter Options"*:
>        - *"Maximum error rate"*: `0.1`
>        - *"Minimum overlap length"*: `35`
>        - *"Look for adapters in the reverse complement"*: `Yes`
>    - In *"Filter Options"*:
>        - *"Discard Trimmed Reads"*: `Yes`
>
>    > <tip-title>Select collection dataset</tip-title>
>    >
>    > 1. Click on {% icon param-collection %} **Dataset collection** in front of the input parameter you want to supply the collection to.
>    > 2. Select the collection you want to use from the list
>    >
>    {: .tip}
>
> 2. Rename the output file as `HiFi_collection (trim)`.
>
> {% snippet faqs/galaxy/datasets_rename.md %}
>
{: .hands_on}


# Genome profile analysis

[{% icon exchange %} Switch to short version]({% link topics/assembly/tutorials/vgp_workflow_training/tutorial.md %}#genome-profile-analysis)

Before starting a *de novo* genome assembly project, it is useful to collect metrics on the properties of the genome under consideration, such as the expected genome size. Traditionally, DNA flow cytometry was considered the golden standard for estimating the genome size. Nowadays, experimental methods have been replaced by computational approaches ({% cite wang2020estimation %}). One of the widely used genome profiling methods is based on the analysis of k-mer frequencies. It allows one to provide information not only about the genomic complexity, such as the genome size and levels of heterozygosity and repeat content, but also about the data quality.

> <comment-title>K-mer size, sequencing coverage and genome size</comment-title>
>
>*K*-mers are unique substrings of length k contained within a DNA sequence. For example, the DNA sequence *TCGATCACA* can be decomposed into five unique *k*-mers that have five bases long: *TCGAT*, *CGATC*, *GATCA*, *ATCAC* and *TCACA*. A sequence of length L will have  L-k+1 *k*-mers. On the other hand, the number of possible *k*-mers can be calculated as  n<sup>k</sup>, where n is the number of possible monomers and k is the k-mer size.
>
>
>---------| -------------|-----------------------
>  Bases  |  K-mer size  |  Total possible k-mers
>---------| -------------|-----------------------
>    4    |       1      |            4          
>    4    |       2      |           16          
>    4    |       3      |           64          
>    4    |      ...     |          ...          
>    4    |      10      |    1.048.576         
>---------|--------------|-----------------------
>
> Thus, the k-mer size is a key parameter, which must be large enough to map  uniquely to the genome, but not too large, since it can lead to wasting computational resources. In the case of the human genome, *k*-mers of 31 bases in length lead to 96.96% of unique *k*-mers.
>
> Each unique k-mer can be assigned a value for coverage based on the number of times it occurs in a sequence, whose distribution will approximate a Poisson distribution, with the peak corresponding to the average genome sequencing depth. From the genome coverage, the genome size can be easily computed.
{: .comment}


## Generation of k-mer spectra with **Meryl**

Meryl will allow us to generate the *k*-mer profile by decomposing the sequencing data into *k*-length substrings, counting the occurrence of each *k*-mer and determining its frequency. The original version of Meryl was developed for the Celera Assembler. The current Meryl version comprises three main modules: one for generating *k*-mer databases, one for filtering and combining databases, and one for searching databases. *K*-mers are stored in lexicographical order in the database, similar to words in a dictionary ({% cite Rhie2020 %}).

> <comment-title><i>k</i>-mer size estimation</comment-title>
>
>  Given an estimated genome size (G) and a tolerable collision rate (p), an appropriate k can be computed as k = log4 (G(1 − p)/p).
>
{: .comment}

> <hands-on-title>Generate <i>k</i>-mers count distribution</hands-on-title>
>
> 1. Run {% tool [Meryl](toolshed.g2.bx.psu.edu/repos/iuc/meryl/meryl/1.3+galaxy2) %} with the following parameters:
>    - *"Operation type selector"*: `Count operations`
>        - *"Count operations"*: `Count: count the occurrences of canonical k-mers`
>        - {% icon param-collection %} *"Input sequences"*: `HiFi_collection (trim)`
>        - *"k-mer size selector"*: `Set a k-mer size`
>            - "*k-mer size*": `32`
>
>    > <comment-title>Selection of <i>k</i>-mer size</comment-title>
>    >
>    > We used 32 as *k*-mer size, as this length has demonstrated to be sufficiently long that most *k*-mers are not repetitive and is short enough to be more robust to sequencing errors. For very large (haploid size > 10 Gb) and/or very repetitive genomes, larger *k*-mer length is recommended to increase the number of unique *k*-mers.
>    {: .comment}
>
> 2. Rename it `Collection meryldb`
>
> 3. Run {% tool [Meryl](toolshed.g2.bx.psu.edu/repos/iuc/meryl/meryl/1.3+galaxy1) %} again with the following parameters:
>    - *"Operation type selector"*: `Operations on sets of k-mers`
>        - *"Operations on sets of k-mers"*: `Union-sum: return k-mers that occur in any input, set the count to the sum of the counts`
>        - {% icon param-file %} *"Input meryldb"*: `Collection meryldb`
>
> 4. Rename it as `Merged meryldb`    
>
> 5. Run {% tool [Meryl](toolshed.g2.bx.psu.edu/repos/iuc/meryl/meryl/1.3+galaxy0) %} for the third time with the following parameters:
>    - *"Operation type selector"*: `Generate histogram dataset`
>        - {% icon param-file %} *"Input meryldb"*: `Merged meryldb`
>
> 6. Finally, rename it as `Meryldb histogram`.
>
{: .hands_on}


## Genome profiling with **GenomeScope2**

The next step is to infer the genome properties from the *k*-mer histogram generated by Meryl, for which we will use GenomeScope2. Genomescope2 relies on a nonlinear least-squares optimization to fit a mixture of negative binomial distributions, generating estimated values for genome size, repetitiveness, and heterozygosity rates ({% cite RanalloBenavidez2020 %}).

> <hands-on-title>Estimate genome properties</hands-on-title>
>
> 1. {% tool [GenomeScope](toolshed.g2.bx.psu.edu/repos/iuc/genomescope/genomescope/2.0) %} with the following parameters:
>    - {% icon param-file %} *"Input histogram file"*: `Meryldb histogram`
>    - *"k-mer length used to calculate k-mer spectra"*: `32`
>
>   - In "*Output options*": mark `Summary of the analysis`
>   - In "*Advanced options*":
>       - *"Create testing.tsv file with model parameters"*: `Yes`
>
>    {: .comment}
>
{: .hands_on}

Genomescope will generate six outputs:

- Plots
    - *Linear plot*: *k*-mer spectra and fitted models: frequency (y-axis) versus coverage.
    - *Log plot*: logarithmic transformation of the previous plot.
    - Transformed linear plot: *k*-mer spectra and fitted models: frequency times coverage (y-axis) versus coverage (x-axis). This transformation increases the heights of higher-order peaks, overcoming the effect of high heterozygosity.
    - Transformed log plot: logarithmic transformation of the previous plot.
- Model: this file includes a detailed report about the model fitting.
- Summary: it includes the properties inferred from the model, such as genome haploid length and the percentage of heterozygosity.

Now, let's analyze the *k*-mer profiles, fitted models and estimated parameters (fig. 3).

![Figure 3: Genomescope plot](../../images/vgp_assembly/genomescope_plot.png "GenomeScope2 32-mer profile. The first peak located at coverage 32x corresponds to the heterozygous peak. The second peak at coverage 50x, corresponds to the homozygous peak. Estimate of the heterozygous portion is 0.576%. The plot also includes informatin about the inferred total genome length (len), genome unique length percent (uniq), overall heterozygosity rate (het), mean k-mer coverage for heterozygous bases (kcov), read error rate (err), average rate of read duplications (dup) and k-mer size (k).")

This distribution is the result of the Poisson process underlying the generation of sequencing reads. As we can see, the k-mer profile follows a bimodal distribution, indicative of a diploid genome. The distribution is consistent with the theoretical diploid model (model fit > 93%). Low frequency *k*-mers are the result of sequencing errors. GenomeScope2 estimated a haploid genome size is around 11.7 Mb, a value reasonably close to *Saccharomyces* genome size. Additionally, it revealed that the variation across the genomic sequences is 0.576%.
    
Before proceeding to the next section, we need to carry out some text manipulation operations on the output generated by GenomeScope2 to make it compatible with downstream tools. The goal is to extract some parameters which at a later stage will be used by purge_dups.

The first relevant parameter is the `estimated genome size`.

> <hands-on-title>Get estimated genome size</hands-on-title>
>
> 1. {% tool [Replace parts of text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_find_and_replace/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `summary` (output of **GenomeScope** {% icon tool %})
>    - *"Find pattern"*: `bp`
>    - *"Replace with*": leave this field empty (as it is)
>    - *"Replace all occurrences of the pattern"*: `Yes`
>    - *"Find and Replace text in"*: `entire line`
>
> 2. {% tool [Replace parts of text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_find_and_replace/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: output file of **Replace** {% icon tool %}
>    - *"Find pattern"*: `,`
>    - *"Replace with*": leave this field empty (as it is)
>    - *"Replace all occurrences of the pattern"*: `Yes`
>    - *"Find and Replace text in"*: `entire line`
>
> 3. {% tool [Search in textfiles (grep)](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_grep_tool/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: output file of the previous step.
>    - *"Type of regex"*: `Basic`
>    - *"Regular Expression"*: `Haploid`
>
> 4. {% tool [Convert delimiters to TAB](Convert characters1) %} with the following parameters:
>    - {% icon param-file %} *"in Dataset"*: output of **Search in textfiles** {% icon tool %}
>
> 5. {% tool [Advanced Cut](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/1.1.0) %} with the following parameters:
>    - {% icon param-file %} *"File to cut"*: output of **Convert delimiters to TAB** {% icon tool %}
>    - *"Cut by"*: `fields`
>        - *"List of Fields"*: `Column: 5`
>
> 6. Rename the output as `Estimated genome size`.
>
> > <question-title></question-title>
> >
> > What is the estimated genome size?
> >
> > > <solution-title></solution-title>
> > >
> > > The estimated genome size is 11.747.076 bp.
> > >
> > {: .solution}
> >
> {: .question}
>
{: .hands_on}

Now let's parse the `upper bound for the read depth estimation` parameter. This parameter will be used later by purge_dups as high read depth cutoff to identify *junk contigs*.

> <hands-on-title>Get maximum read depth</hands-on-title>
>
> 1. {% tool [Compute an expression on every row](toolshed.g2.bx.psu.edu/repos/devteam/column_maker/Add_a_column1/1.6) %} with the following parameters:
>    - *"Add expression"*: `1.5*c3`
>    - {% icon param-file %} *"as a new column to"*: `model_params` (output of **GenomeScope** {% icon tool %})
>    - *"Round result?"*: `Yes`
>    - *"Input has a header line with column names?"*: `No`
>
> 2. {% tool [Compute an expression on every row](toolshed.g2.bx.psu.edu/repos/devteam/column_maker/Add_a_column1/1.6) %} with the following parameters:
>    - *"Add expression"*: `3*c7`
>    - {% icon param-file %} *"as a new column to"*: output of the previous step.
>    - *"Round result?"*: `Yes`
>    - *"Input has a header line with column names?"*: `No`
>
> 3. Rename it as `Parsing temporal output`
>
> 4. {% tool [Advanced Cut](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/1.1.0) %} with the following parameters:
>    - {% icon param-file %} *"File to cut"*: `Parsing temporal output`
>    - *"Cut by"*: `fields`
>        - *"List of Fields"*: `Column 8`
>
> 6. Rename the output as `Maximum depth`
>
> > <question-title></question-title>
> >
> > 1. What is the estimated maximum depth?
> > 2. What does this parameter represent?
> >
> > > <solution-title></solution-title>
> > >
> > > 1. The estimated maximum depth is 114 reads.
> > > 2. The maximum depth indicates the maximum number of sequencing reads that align to specific positions in the genome.
> > >
> > {: .solution}
> >
> {: .question}
>
{: .hands_on}

Finally, let's parse the `transition between haploid and diploid coverage depths` parameter. This will be used as threshold to discriminate between haploid and diploid depths.

> <hands-on-title>Get transition parameter        </hands-on-title>
> 1. {% tool [Advanced Cut](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/1.1.0) %} with the following parameters:
>    - {% icon param-file %} *"File to cut"*: `Parsing temporal output`
>    - *"Cut by"*: `fields`
>        - *"List of Fields"*: `Column 7`
>
> 2. Rename the output as `Transition parameter`
>
> > <question-title></question-title>
> >
> > What is the estimated transition parameter?
> >
> > > <solution-title></solution-title>
> > >
> > > The estimated transition parameter is 38 reads.
> > >
> > {: .solution}
> >
> {: .question}
>
{: .hands_on}


# HiFi phased assembly with hifiasm

[{% icon exchange %} Switch to short version]({% link topics/assembly/tutorials/vgp_workflow_training/tutorial.md %}#hifi-phased-assembly-with-hifiasm)

Once we have finished the genome profiling stage, we can start the genome assembly with hifiasm,  a fast open-source *de novo* assembler specifically developed for PacBio HiFi reads.

## Genome assembly with **hifiasm**

One of the key advantages of hifiasm is that it allows us to resolve near-identical, but not exactly identical sequences, such as repeats and segmental duplications ({% cite Cheng2021 %}).

> <comment-title>Hifiasm algorithm details</comment-title>
>
> By default hifiasm performs three rounds of haplotype-aware error correction to correct sequence errors but keeping heterozygous alleles. A position on the target read to be corrected is considered informative if there are two different nucleotides at that position in the alignment, and each allele is supported by at least three reads.
>
> ![Figure 4: Hifiasm algorithm overview](../../images/vgp_assembly/hifiasm_algorithm.png "Hifiasm algorithm overview. Orange and blue bars represent the reads with heterozygous alleles carrying local phasing information, while green bars come from the homozygous regions without any heterozygous alleles.")
>
> Then, hifiasm builds a phased assembly string graph with local phasing information from the corrected reads. Only the reads coming from the same haplotype are connected in the phased assembly graph. After transitive reduction, a pair of heterozygous alleles is represented by a _bubble_ in the string graph. If there is no additional data, hifiasm arbitrarily selects one side of each bubble and outputs a primary assembly. In the case of a heterozygous genome, the primary assembly generated at this step may still retain haplotigs from the alternate allele.
>
>
{: .comment}

> <hands-on-title>Phased assembly with <b>hifiasm</b></hands-on-title>
>
> 1. {% tool [Hifiasm](toolshed.g2.bx.psu.edu/repos/bgruening/hifiasm/hifiasm/0.14+galaxy0) %} with the following parameters:
>    - *"Assembly mode"*: `Standard`
>        - {% icon param-file %} *"Input reads"*: `HiFi_collection (trim)` (output of **Cutadapt** {% icon tool %})
>    - *"Options for purging duplicates"*: `Specify`
>       - *"Purge level*": `Light`
>       - *"Coverage upper bound"*: `114` (maximum depth previously obtained)
>    - *"Options for Hi-C partition"*: `Specify`
>       - *"Hi-C R1 reads"*: `Hi-C_dataset_F`
>       - *"Hi-C R2 reads"*: `Hi-C_dataset_R`
>
> 2. After the tool has finished running rename its outputs as follows:
>   - Rename the `Hi-C hap1 balanced contig graph` as `Primary contigs graph` and add a `#primary` tag
>   - Rename the `Hi-C hap2 balanced contig graph` as `Alternate contigs graph` and  add a `#alternate` tag
>
{: .hands_on}

Hifiasm generates four outputs in Graphical Fragment Assembly (GFA) format; this format is designed to represent genome variation, splice graphs in genes, and even overlaps between reads.

We have obtained the fully phased contig graphs of the primary and alternate haplotypes, but the output format of hifiasm must be converted to FASTA format for the subsequent steps.


> <hands-on-title>convert GFA to FASTA</hands-on-title>
>
> 1. {% tool [GFA to FASTA](toolshed.g2.bx.psu.edu/repos/iuc/gfa_to_fa/gfa_to_fa/0.1.2) %} with the following parameters:
>    - {% icon param-files %} *"Input GFA file"*: select `Primary contigs graph` and the `Alternate contigs graph` datasets
>
>    > <tip-title>Select multiple datasets</tip-title>
>    > 1. Click on {% icon param-files %} **Multiple datasets**
>    > 2. Select several files by keeping the <kbd>Ctrl</kbd> (or <kbd>COMMAND</kbd>) key pressed and clicking on the files of interest
>    {: .tip}
>
> 2. Rename the outputs as `Primary contigs FASTA` and `Alternate contigs FASTA`
>
{: .hands_on}

## Initial assembly evaluation

The VGP assembly pipeline contains several built-in QC steps, including QUAST, BUSCO (Benchmarking Universal Single-Copy Orthologs), Merqury, and Pretext. QUAST will generate summary statistics, BUSCO will search for universal single-copy ortholog genes, Merqury will evaluate assembly copy-numbers using *k*-mers, and Pretext will be used to evaluate the assembly contiguity.

> <comment-title>QUAST statistics</comment-title>
>
> QUAST will provide us with the following statistics:
>
> - No. of contigs: The total number of contigs in the assembly.
> - Largest contig: The length of the largest contig in the assembly.
> - Total length: The total number of bases in the assembly.
> - Nx: The largest contig length, *L*, such that contigs of length >= *L* account for at least *x*% of the bases of the assembly.
> - NGx: The contig length such that using equal or longer length contigs produces *x*% of the length of the reference genome, rather than *x*% of the assembly length.
> - GC content: the percentage of nitrogenous bases which are either guanine or cytosine.
>
{: .comment}

QUAST default pipeline utilizes Minimap2. Functional elements prediction modules include GeneMarkS, GeneMark-ES, GlimmerHMM, Barrnap, and BUSCO.
    
> <hands-on-title>assembly evaluation with QUAST</hands-on-title>
>
> 1. {% tool [Quast](toolshed.g2.bx.psu.edu/repos/iuc/quast/quast/5.0.2+galaxy1) %} with the following parameters:
>    - *"Use customized names for the input files?"*: `Yes, specify custom names`
>    - In *"1. Contigs/scaffolds"*:
>        - {% icon param-file %} *"Contigs/scaffolds file"*: `Primary contigs FASTA`
>        - *"Name"*: `Primary assembly`
>    - Click in *"Insert Contigs/scaffolds"*
>    - In *"2. Contigs/scaffolds"*:
>        - {% icon param-file %} *"Contigs/scaffolds file"*: `Alternate contigs FASTA`
>        - *"Name"*: `Alternate assembly`
>    - *"Reads options"*: `Pacbio SMRT reads`
>        - {% icon param-collection %} *"FASTQ file"*: `HiFi collection (trim)`
>    - *"Type of assembly"*: `Genome`
>        - *"Use a reference genome?"*: `No`
>            - *"Estimated reference genome size (in bp) for computing NGx statistics"*: `11747076` (previously estimated)
>        - *"Type of organism"*: `Eukaryote: use of GeneMark-ES for gene finding, Barrnap for ribosomal RNA genes prediction, BUSCO for conserved orthologs finding (--eukaryote)`
>    - *"Is genome large (>100Mpb)?"*: `No`
>
>    > <comment-title></comment-title>
>    >
>    > Remember that for this training we are using _S. cerevisiae_, a reduced genome. In the case of assembling a vertebrate genome, you must select `yes` in the previous option.
>    {: .comment}
>
> 2. Rename the HTML report as `QUAST initial report`
>
{: .hands_on}

    
Let's have a look at the HTML report generated by QUAST. 

![Figure 5: QUAST initial plot](../../images/vgp_assembly/QUAST_initial.png "QUAST HTML report. Statistics of the primary and alternate assembly (a). Cumulative length plot. It shows the growth of contig lengths. On the x-axis, contigs are ordered from the largest to smallest. The y-axis gives the size of the x largest contigs in the assembly (b).")

According to the report, both assemblies are quite similar; the primary assembly includes 18 contigs, whose cumulative length is around 12.1Mbp. On the other hand, the second alternate assembly includes 17 contigs, whose total lenght is 11.3Mbp. As we can see in the figure 5b, the alternate assembly is slighly smaller than the primary assembly.

> <question-title></question-title>
>
> 1. What is the longest contig in the primary assembly? And in the alternate one?
> 2. What is the N50 of the primary assembly?
> 3. Which percentage of reads mapped to each assembly?
>
> > <solution-title></solution-title>
> >
> > 1. The longest contig in the primary assembly is 1.531.728 bp, and 1.532.843 bp in the alternate assembly.
> > 2. The N50 of the primary assembly is 813.039 bp.
> > 3. According to the report, 100% of reads mapped to both the primary assembly and 97.98% to the alternate assembly.
> >
> {: .solution}
>
{: .question}

Next, we will use BUSCO, which will provide quantitative assessment of the completeness of a genome assembly in terms of expected gene content. It relies on the analysis of genes that should be present only once in a complete assembly or gene set, while allowing for rare gene duplications or losses ({% cite Simo2015 %}).

> <hands-on-title>assessing assembly completeness with BUSCO</hands-on-title>
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.0.0+galaxy0) %} with the following parameters:
>    - {% icon param-files %} *"Sequences to analyze"*: `Primary contigs FASTA` and `Alternate contigs FASTA`
>    - *"Mode"*: `Genome assemblies (DNA)`
>        - *"Use Augustus instead of Metaeuk"*: `Use Metaeuk`
>    - *"Auto-detect or select lineage?"*: `Select lineage`
>       - *"Lineage"*: `Saccharomycetes`
>    - *"Which outputs should be generated"*: `short summary text`
>
>    > <comment-title></comment-title>
>    >
>    > Remember to modify the lineage option if you are working with vertebrate genomes.
>    {: .comment}
>
> 2. Rename the summary as `BUSCO initial report primary` and `BUSCO initial report alternate`.
>
{: .hands_on}

BUSCO generates two outputs by default: a complete report for each BUSCO gene (fig. 6), and a short summary.

![Figure 6 : BUSCO](../../images/vgp_assembly/BUSCO_full_table.png "BUSCO full table. It contains the complete results in a tabular format with scores and lengths of BUSCO matches, and coordinates.")

As we can see in the figure, the detailed report includes valuable information, such as the status of the gene identified in the input sequence, its position (start/end coordinates and strand), the BUSCO score, the gene length and a keyword description.
    
> <question-title></question-title>
>
> 1. How many complete BUSCO genes have been identified in the primary assembly? You can find that information in the short summary.
> 2. How many BUSCOs genes are absent?
>
> > <solution-title></solution-title>
> >
> > 1. According to the report, our assembly contains the complete sequence of 2046 complete BUSCO genes (95.7%).
> > 2. 29 BUSCO genes are missing.
> >
> {: .solution}
>
{: .question}

Despite BUSCO being robust for species that have been widely studied, it can be inaccurate when the newly assembled genome belongs to a taxonomic group that is not well represented in [OrthoDB](https://www.orthodb.org/). Merqury provides a complementary approach for assessing genome assembly quality metrics in a reference-free manner via *k*-mer copy number analysis.

> <hands-on-title><i>k</i>-mer based evaluation with Merqury</hands-on-title>
>
> 1. {% tool [Merqury](toolshed.g2.bx.psu.edu/repos/iuc/merqury/merqury/1.3) %} with the following parameters:
>    - *"Evaluation mode"*: `Default mode`
>        - {% icon param-file %} *"k-mer counts database"*: `Merged meryldb`
>        - *"Number of assemblies"*: `Two assemblies
>            - {% icon param-file %} *"First genome assembly"*: `Primary contigs FASTA`
>            - {% icon param-file %} *"Second genome assembly"*: `Alternate contigs FASTA`    
>
{: .hands_on}


By default, Merqury generates three collections as output: stats, plots and QV stats. The "stats" collection contains the completeness statistics, while the "QV stats" collection contains the quality value statistics. Let's have a look at the primary assembly copy number (CN) spectrum plot, known as the *spectra-cn* plot (fig. 7).

![Figure 7: Merqury plot](../../images/vgp_assembly/merqury_cn_plot.png "Merqury CN plot. This plot tracks the multiplicity of each k-mer found in the Hi-Fi read set and colors it by the number of times it is found in a given assembly. Merqury connects the midpoint of each histogram bin with a line, giving the illusion of a smooth curve."){:width="80%"}

The black region in the left side corresponds to k-mers found only in the read set; it is usually indicative of sequencing error in the read set, although it can also be indicative of missing sequences in the assembly. The read area represents one-copy k-mers in the genome, while blue area represents two-copy k-mers originating from homozygous sequence or haplotype-specific duplications. From that figure we can state that the sequencing coverage is around 50x.

# Post-assembly processing

[{% icon exchange %} Switch to short version]({% link topics/assembly/tutorials/vgp_workflow_training/tutorial.md %}#post-assembly-processing)

An ideal haploid representation would consist of one allelic copy of all heterozygous regions in the two haplomes, as well as all hemizygous regions from both haplomes ({% cite Guan2019 %}). However, in highly heterozygous genomes, assembly algorithms are frequently not able to identify the highly divergent allelic sequences as belonging to the same region, resulting in the assembly of those regions as separate contigs. This can lead to issues in downstream analysis, such as scaffolding, gene annotation and read mapping in general ({% cite Small2007 %}, {% cite Guan2019 %}, {% cite Roach2018 %}). In order to solve this problem, we are going to use purge_dups; this tool will allow us to identify and reassign allelic contigs.

## Remove haplotypic duplication with **purge_dups**
        
This stage consists of three substages: read-depth analysis, generation of all versus all self-alignment and resolution of haplotigs and overlaps (fig. 8).

![Figure 8: Post-processing with purge_dups](../../images/vgp_assembly/purge_dupspipeline.png "Purge_dups pipeline. Adapted from github.com/dfguan/purge_dups. Purge_dups is integrated in a multi-step pipeline consisting in three main substages. Red indicates the steps which require to use Minimap2.")

### Read-depth analysis


Initially, we need to collapse our HiFi trimmed reads collection into a single dataset.
    
> <hands-on-title>Collapse the collection</hands-on-title>
>
> 1. {% tool [Collapse Collection](toolshed.g2.bx.psu.edu/repos/nml/collapse_collections/collapse_dataset/4.2) %} with the following parameters:
>    - {% icon param-collection %} *"Collection of files to collapse into single dataset"*:`HiFi_collection (trim)`
> 2. Rename the output as `HiFi reads collapsed`
{: .hands_on}

Now, we will map the reads against the primary assembly by using Minimap2 ({% cite Li2018 %}), an alignment program designed to map long sequences.

> <hands-on-title>Map the reads to contigs with <b>Minimap2</b></hands-on-title>
>
> 1. {% tool [Map with minimap2](toolshed.g2.bx.psu.edu/repos/iuc/minimap2/minimap2/2.17+galaxy4) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `Primary contigs FASTA`
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-collection %} *"Select fastq dataset"*: `HiFi reads collapsed`
>        - *"Select a profile of preset options"*: `Long assembly to reference mapping (-k19 -w19 -A1 -B19 -O39,81 -E3,1 -s200 -z200 --min-occ-floor=100). Typically, the alignment will not extend to regions with 5% or higher sequence divergence. Only use this preset if the average divergence is far below 5%. (asm5)`
>    - In *"Set advanced output options"*:
>        - *"Select an output format"*: `paf`
>
> 2. Rename the output as `Reads mapped to contigs`
{: .hands_on}

Finally, we will use the `Reads mapped to contigs` pairwise mapping format (PAF) file for calculating some statistics required in a later stage. In this step, purge_dups (listed as **Purge overlaps** in Galaxy tool panel) initially produces a read-depth histogram from base-level coverages. This information is used for estimating the coverage cutoffs, taking into account that collapsed haplotype contigs will lead to reads from both alleles mapping to those contigs, whereas if the alleles have assembled as separate contigs, then the reads will be split over the two contigs, resulting in half the read-depth ({% cite Roach2018 %}). 

> <hands-on-title>Read-depth analisys</hands-on-title>
> 1. {% tool [Purge overlaps](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy3) %} with the following parameters:
>    - *"Function mode"*: `Calculate coverage cutoff, base-level read depth and create read depth histogram for PacBio data (calcuts+pbcstats)`
>        - {% icon param-file %} *"PAF input file"*: `Reads mapped to contigs`
>        - In *"Calcuts options"*:
>            - *"Transition between haploid and diploid"*: 38
>            - *"Upper bound for read depth"*: `114` (the previously estimated maximum depth)
>            - *"Ploidity"*: `Diploid`
>
> 2. Rename the outputs as `PBCSTAT base coverage primary`, `Histogram plot primary` and `Calcuts cutoff primary`.
{: .hands_on}

Purge overlaps (purge_dups) generates three outputs:

- PBCSTAT base coverage: it contains the base-level coverage information.
- Calcuts-cutoff: it includes the thresholds calculated by purge_dups.
- Histogram plot.


    
### Generation of all versus all self-alignment

Now, we will segment the draft assembly into contigs by cutting at blocks of *N*s, and use minimap2 to generate an all by all self-alignment.

> <hands-on-title>purge_dups pipeline    </hands-on-title>
> 1. {% tool [Purge overlaps](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy2) %} with the following parameters:
>    - *"Function mode"*: `split assembly FASTA file by 'N's (split_fa)`
>        - {% icon param-file %} *"Assembly FASTA file"*: `Primary contigs FASTA`
>
> 2. Rename the output as `Split FASTA`
>
> 3. {% tool [Map with minimap2](toolshed.g2.bx.psu.edu/repos/iuc/minimap2/minimap2/2.17+galaxy4) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `Split FASTA`
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-file %} *"Select fastq dataset"*: `Split FASTA`
>        - *"Select a profile of preset options"*: `Construct a self-homology map - use the same genome as query and reference (-DP -k19 -w 19 -m200) (self-homology)`
>    - In *"Set advanced output options"*:
>        - *"Select an output format"*: `PAF`
>
> 4. Rename the output as `Self-homology map primary`
{: .hands_on}


### Resolution of haplotigs and overlaps        

During the final step of the purge_dups pipeline, it will use the self alignments and the cutoffs for identifying the haplotypic duplications.

> <comment-title>Purge overlaps (purge_dups) algorithm details</comment-title>
>
> In order to identify the haplotypic duplications, purge_dups uses the  base-level coverage information to flag the contigs according the following criteria:
> - If more than 80% bases of a contig are above the high read depth cutoff or below the noise cutoff, it is discarded.
> - If more than 80% bases are in the diploid depth interval, it is labeled as a primary contig, otherwise it is considered further as a possible haplotig.
>
> Contigs that were flagged for further analysis according to read-depth are then evaluated to attempt to identify synteny with its allelic companion contig. In this step, purge_dups uses the information contained in the self alignments:
> - If the alignment score is larger than the cutoff *s* (default 70), the contig is marked for reassignment as haplotig. Contigs marked for reassignment with a maximum match score greater than the cutoff *m* (default 200) are further flagged as repetitive regions.
>
> - Otherwise contigs are considered as a candidate primary contig.
>
> Once all matches associated with haplotigs have been removed from the self-alignment set, purge_dups ties consistent matches between the remaining candidates to find collinear matches, filtering all the matches whose score is less than the minimum chaining score *l*.
>
> Finally, purge_dups calculates the average coverage of the matching intervals for each overlap, and mark an unambiguous overlap as heterozygous when the average coverage on both contigs is less than the read-depth cutoff, removing the sequences corresponding to the matching interval in the shorter contig.
>
{: .comment}

> <hands-on-title>Resolution of haplotigs and overlaps</hands-on-title>
>    
> 1. {% tool [Purge overlaps](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy5) %} with the following parameters:
>    - *"Select the purge_dups function"*: `Purge haplotigs and overlaps for an assembly (purge_dups)`
>        - {% icon param-file %} *"PAF input file"*: `Self-homology map primary`
>        - {% icon param-file %} *"Base-level coverage file"*: `PBCSTAT base coverage primary`
>        - {% icon param-file %} *"Cutoffs file"*: `calcuts cutoff primary`
>
> 2. Rename the output as `purge_dups BED`
>
> 3. {% tool [Purge overlaps](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy2) %} with the following parameters:
>    - *"Select the purge_dups function"*: `Obtain sequences after purging (get_seqs)`
>        - {% icon param-file %} *"Assembly FASTA file"*: `Primary contigs FASTA`
>        - {% icon param-file %} *"BED input file"*: `purge_dups BED` (output of the previous step)
>
> 4. Rename the output `get_seq purged sequences` as `Primary contigs purged` and the `get_seq haplotype` file as `Alternate haplotype contigs`
>
{: .hands_on}


### Process the alternative assembly

Now we should repeat the same procedure with the alternate contigs generated by hifiasm.  In that case, we should start by merging the `Alternate haplotype contigs` generated in the previous step and the `Alternate contigs FASTA` file.
    
> <hands-on-title>Merge the purged sequences and the Alternate contigs</hands-on-title>
>
> 1. {% tool [Concatenate datasets](cat1) %} with the following parameters:
>    - {% icon param-file %} *"Concatenate Dataset"*: `Alternate contigs FASTA`
>    - In *"Dataset"*:
>        - {% icon param-repeat %} *"Insert Dataset"*
>            - {% icon param-file %} *"Select"*: `Alternate haplotype contigs`
>
>    > <comment-title></comment-title>
>    >
>    > Remember that the `Alternate haplotype contigs` file contains those contigs that were considered to be haplotypic duplications of the primary contigs.
>    {: .comment}
>
> 2. Rename the output as `Alternate contigs full`
>
{: .hands_on}

Once we have merged the files, we should run the purge_dups pipeline again, but using the `Alternate contigs full` file as input.

> <hands-on-title>Process the alternate assembly with <i>purge_dups</i></hands-on-title>
>
> 1. {% tool [Map with minimap2](toolshed.g2.bx.psu.edu/repos/iuc/minimap2/minimap2/2.17+galaxy4) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `Alternate contigs full`
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-collection %} *"Select fastq dataset"*: `HiFi reads collapsed`
>        - *"Select a profile of preset options"*: `Long assembly to reference mapping (-k19 -w19 -A1 -B19 -O39,81 -E3,1 -s200 -z200 --min-occ-floor=100). Typically, the alignment will not extend to regions with 5% or higher sequence divergence. Only use this preset if the average divergence is far below 5%. (asm5)` (**Note** `asm5` at the end!)
>    - In *"Set advanced output options"*:
>        - *"Select an output format"*: `paf`
>
> 2. Rename the output as `Reads mapped to contigs alternate`
>
> 3. {% tool [Purge overlaps](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy3) %} with the following parameters:
>    - *"Function mode"*: `Calculate coverage cutoff, base-level read depth and create read depth histogram for PacBio data (calcuts+pbcstats)`
>        - {% icon param-file %} *"PAF input file"*: `Reads mapped to contigs alternate`
>        - In *"Calcuts options"*:
>            - *"Upper bound for read depth"*: `114`
>            - *"Ploidity"*: `Diploid`
>
> 3. Rename the outputs as `PBCSTAT base coverage alternate`, `Histogram plot alternate` and `Calcuts cutoff alternate`.
>
> 4. {% tool [Purge overlaps](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy2) %} with the following parameters:
>    - *"Function mode"*: `split assembly FASTA file by 'N's (split_fa)`
>        - {% icon param-file %} *"Assembly FASTA file"*: `Alternate contigs full`
>
> 5. Rename the output as `Split FASTA alternate`
>
> 6. {% tool [Map with minimap2](toolshed.g2.bx.psu.edu/repos/iuc/minimap2/minimap2/2.17+galaxy4) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `Split FASTA alternate`
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-file %} *"Select fastq dataset"*: `Split FASTA alternate`
>        - *"Select a profile of preset options"*: `Construct a self-homology map - use the same genome as query and reference (-DP -k19 -w 19 -m200) (self-homology)`
>    - In *"Set advanced output options"*:
>        - *"Select an output format"*: `PAF`
>
> 7. Rename the output as `Self-homology map alternate`
>        
> 8. {% tool [Purge overlaps](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy5) %} with the following parameters:
>    - *"Select the purge_dups function"*: `Purge haplotigs and overlaps for an assembly (purge_dups)`
>        - {% icon param-file %} *"PAF input file"*: `Self-homology map alternate`
>        - {% icon param-file %} *"Base-level coverage file"*: `PBCSTAT base coverage alternate`
>        - {% icon param-file %} *"Cutoffs file"*: `calcuts cutoff alternate`
>
>
> 9. Rename the output as `purge_dups BED alternate`
> 
> 10. {% tool [Purge overlaps](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy2) %} with the following parameters:
>    - *"Select the purge_dups function"*: `Obtain sequences after purging (get_seqs)`
>        - {% icon param-file %} *"Assembly FASTA file"*: `Alternate contigs full`
>        - {% icon param-file %} *"BED input file"*: `purge_dups BED alternate`
>
> 11. Rename the outputs as `Alternate contigs purged` and `Alternate haplotype contigs`.
>
{: .hands_on}


# Hybrid scaffolding

At this point, we have obtained the primary and alternate assemblies, each of which consists in a collection of contigs (contiguous sequences assembled from overlapping reads). Next, the contigs will be assembled into scaffolds, i.e., sequences of contigs interspaced with gaps. For this purpose, we will carry out a hybrid scaffolding by taking advantage of two additional technologies: Bionano optical maps and Hi-C data. 

## Hybrid scaffolding using Bionano data

[{% icon exchange %} Switch to short version]({% link topics/assembly/tutorials/vgp_workflow_training/tutorial.md %}#hybrid-scaffolding-with-bionano-optical-maps)

In this step, the linkage information provided by optical maps is integrated with primary assembly sequences, and the overlaps are used to orient and order the contigs, resolve chimeric joins, and estimate the length of gaps between adjacent contigs. One of the advantages of optical maps is that they can easily span genomic regions that are difficult to resolve using DNA sequencing technologies ({% cite Savara2021 %}, {% cite Yuan2020 %}).

> <comment-title>Background on Bionano optical maps</comment-title>
>
> Bionano technology relies on the isolation of kilobase-long DNA fragments, which are labeled at specific sequence motifs with a fluorescent dye, resulting in a unique fluorescent pattern for each genome. DNA molecules are stretched into nanoscale channels and imaged with a high-resolution camera, allowing us to build optical maps that include the physical locations of labels rather than base-level information ({% cite Lam2012 %}, {% cite Giani2020 %}, {% cite Savara2021 %}).
>
> ![Figure 9: Bionano optical maps](../../images/vgp_assembly/bionano.png "Bionano optical maps. Source: https://bionanogenomics.com")
>
> The average optical map molecule length, around 225 kbp, is substantially larger than the PacBio HiFi reads, with read lengths averaging 10-25 kbp.
>
{: .comment}

The *Bionano Hybrid Scaffold* tool automates the scaffolding process, which includes five main steps:

1. Generate *in silico* maps for sequence assembly.
2. Align *in silico* sequence maps against Bionano genome maps to identify and resolve potential conflicts.
3. Merge the non-conflicting maps into hybrid scaffolds.
4. Align sequence maps to the hybrid scaffolds
5. Generate AGP and FASTA files for the scaffolds.


> <hands-on-title>Bionano hybrid scaffolding</hands-on-title>
>
> 1. {% tool [Bionano Hybrid Scaffold](toolshed.g2.bx.psu.edu/repos/bgruening/bionano_scaffold/bionano_scaffold/3.6.1+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"NGS FASTA"*: `Primary contigs purged`
>    - {% icon param-file %} *"BioNano CMAP"*: `Bionano_dataset`
>    - *"Configuration mode"*: `VGP mode`
>    - *"Genome maps conflict filter"*: `Cut contig at conflict`
>    - *"Sequences conflict filter"*: `Cut contig at conflict`
>
>    > <comment-title></comment-title>
>    >
>    > If your data are not associated with VGP, make sure that the configuration mode fits with your samples.
>    {: .comment}
>
> 2. {% tool [Concatenate datasets](cat1) %} with the following parameters:
>    - {% icon param-file %} *"Concatenate Dataset"*: `NGScontigs scaffold NCBI trimmed` (output of **Bionano Hybrid Scaffold** {% icon tool %})
>    - In *"Dataset"*:
>        - {% icon param-repeat %} *"Insert Dataset"*
>            - {% icon param-file %} *"Select"*: `NGScontigs not scaffolded trimmed` (output of **Bionano Hybrid Scaffold** {% icon tool %})
>
> 3. Rename the output as `Primary assembly bionano`
{: .hands_on}

## Second round of assembly evaluation

Once we have run purge_dups and Bionano, we can evaluate the assembly again, and compare the results.
    
> <hands-on-title>Bionano assembly evaluation with QUAST and BUSCO</hands-on-title>
>
> 1. {% tool [Quast](toolshed.g2.bx.psu.edu/repos/iuc/quast/quast/5.0.2+galaxy1) %} with the following parameters:
>    - *"Use customized names for the input files?"*: `Yes, specify custom names`
>        - {% icon param-file %} *"Contigs/scaffolds file"*: `Primary assembly bionano`
>        - *"Name"*: `Primary assembly bionano`
>    - *"Reads options"*: `Pacbio SMRT reads`
>        - {% icon param-collection %} *"FASTQ file"*: `HiFi reads collapsed`
>    - *"Type of assembly"*: `Genome`
>        - *"Use a reference genome?"*: `No`
>            - *"Estimated reference genome size (in bp) for computing NGx statistics"*: `11747076` (previously estimated)
>        - *"Type of organism"*: `Eukaryote (--eukaryote): use of GeneMark-ES for gene finding, Barrnap for ribosomal RNA genes prediction, BUSCO for conserved orthologs finding`
>    - In *"Genes"*:
>        - *"Tool for gene prediction"*: `Don't predict genes`
>
> 2. Rename the HTML report as `QUAST second report`
>
> 3. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.0.0+galaxy0) %} with the following parameters:
>    - {% icon param-files %} *"Sequences to analyze"*: `Primary assembly Bionano`
>    - *"Mode"*: `Genome assemblies (DNA)`
>        - *"Use Augustus instead of Metaeuk"*: `Use Metaeuk`
>    - *"Auto-detect or select lineage?"*: `Select lineage`
>       - *"Lineage"*: `Saccharomycetes`
>    - In *"Advanced Options"*:
>        - *"Which outputs should be generated"*: `short summary text` and `summary image`
>
> 4. Rename the summary as `BUSCO second report`
>
{: .hands_on}

As we can observe in the cumulative plot (fig. 10a), the total length of the assembly (12.160.926 bp) is slightly larger than the expected genome size. With respect to the NG50 statistic (fig. 10b), the value is 922.430 bp, which is significantly higher than the value obtained during the first evaluation stage (813.039 bp).
        
![Figure 10: QUAST and BUSCO plots](../../images/vgp_assembly/QUAST_cummulative.png ". Cumulative length plot (a). NGx plot. The y-axis represents the NGx values in Mbp, and the x-axis is the percentage of the genome (b). Assembly evaluation after runnig Bionano. BUSCO genes are defined as "Complete (C) and single copy (S)" when are found once in the single-copy ortholog database, "Complete (C) and duplicated (D)" when single-copy ortholog genes which were found more than once, "Fragmented (F)" when genes are matching just partially to a single-copy ortholog DB, and "Missing (M)" when genes which are expected but were not detected (c).")
    
Finally, the BUSCO's summary image (fig. 10c) shows that most of the universal single-copy orthologs are present in our assembly.

> <question-title></question-title>
>
> 1. How many scaffolds are in the primary assembly after the hybrid scaffolding?
> 2. What is the size of the largest scaffold? Has improved with respect to the previous evaluation?
> 3. What is the percertage of completeness on the core set genes in BUSCO? Has increased the completeness?
>
> > <solution-title></solution-title>
> >
> > 1. The number of contigs is 17.
> > 2. The largest contig is 1.531.728 bp long. This value hasn't changed.
> > 3. The percentage of complete BUSCOs is 95.7%. Yes, it has increased, since in the previous evaluation the completeness percentage was 88.7%.
> >
> {: .solution}
>
{: .question}
    
## Hybrid scaffolding based on Hi-C mapping data

[{% icon exchange %} Switch to short version]({% link topics/assembly/tutorials/vgp_workflow_training/tutorial.md %}#hybrid-scaffolding-with-hi-c-data)

Hi-C is a sequencing-based molecular assay designed to identify regions of frequent physical interaction in the genome by measuring the contact frequency between all pairs of loci, allowing us to provide an insight into the three-dimensional organization of a genome  ({% cite Dixon2012 %}, {% cite LiebermanAiden2009 %}). In this final stage, we will exploit the fact that the contact frequency between a pair of loci strongly correlates with the one-dimensional distance between them with the objective of linking the Bionano scaffolds to a chromosome scale.

> <comment-title>Background about Hi-C data</comment-title>
>
> The high-throughput chromosome conformation capture (Hi-C) technology is based on the capture of the chromatin in three-dimensional space. During Hi-C library preparation, DNA is crosslinked in its 3D conformation. Then, the DNA is digested using restriction enzymes, and the digested ends are filled with biotinylated nucleotides (fig. 11). The biotinylated nucleotides enable the specific purification of the ligation junctions, preventing the sequencing of DNA molecules that do not contain such junctions which are thus mostly uninformative ({% cite Lajoie2015 %}).
>
> ![Figure 11: Hi-C protocol](../../images/vgp_assembly/hi-c_protocol.png "Hi-C protocol. Adapted from Rao et al. 2014")
>
> Next, the blunt ends of the spatially proximal digested end are ligated. Each DNA fragment is then sequenced from each end of this artificial junction, generating read pairs. This provides contact information that can be used to reconstruct the proximity of genomic sequences belonging to the same chromosome ({% cite Giani2020 %}). Hi-C data are in the form of two-dimensional matrices (contact maps) whose entries quantify the intensity of the physical interaction between genome regions.
>
{: .comment}


### Pre-processing Hi-C data

Despite Hi-C generating paired-end reads, we need to map each read separately. This is because most aligners assume that the distance between paired-end reads fit a known distribution, but in Hi-C data the insert size of the ligation product can vary between one base pair to hundreds of megabases ({% cite Lajoie2015 %}).

> <hands-on-title>Mapping Hi-C reads</hands-on-title>
>
> 1. {% tool [BWA-MEM2](toolshed.g2.bx.psu.edu/repos/iuc/bwa_mem2/bwa_mem2/2.2.1+galaxy0) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `Primary assembly bionano`
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-file %} *"Select fastq dataset"*: `Hi-C_dataset_F`
>    - *"Set read groups information?"*: `Do not set`
>    - *"Select analysis mode"*: `1.Simple Illumina mode`
>    - *"BAM sorting mode"*: `Sort by read names  (i.e., the QNAME field) `
>
> 2. Rename the output as `BAM forward`
>
> 3. {% tool [BWA-MEM2](toolshed.g2.bx.psu.edu/repos/iuc/bwa_mem2/bwa_mem2/2.2.1+galaxy0) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `Primary assembly bionano`
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-file %} *"Select fastq dataset"*: `Hi-C_dataset_R`
>    - *"Set read groups information?"*: `Do not set`
>    - *"Select analysis mode"*: `1.Simple Illumina mode`
>    - *"BAM sorting mode"*: `Sort by read names  (i.e., the QNAME field) `
>
> 4. Rename the output as `BAM reverse`
>
> 5. {% tool [Filter and merge](toolshed.g2.bx.psu.edu/repos/iuc/bellerophon/bellerophon/1.0+galaxy0) %} chimeric reads from Arima Genomics with the following parameters:
>    - {% icon param-file %} *"First set of reads"*: `BAM forward`
>    - {% icon param-file %} *"Second set of reads"*: `BAM  reverse`
>
> 6. Rename it as `BAM Hi-C reads`
{: .hands_on}

Finally, we need to convert the BAM file to BED format and sort it.


### Generate initial Hi-C contact map

After mapping the Hi-C reads, the next step is to generate an initial Hi-C contact map, which will allow us to compare the Hi-C contact maps before and after using the Hi-C for scaffolding.

> <comment-title>Biological basis of Hi-C contact maps</comment-title>
>
> Hi-C contact maps reflect the interaction frequency between genomic loci. In order to understand the Hi-C contact maps, it is necessary to take into account two factors: the higher interaction frequency between loci that reside in the same chromosome (_i.e._, in cis), and the distance-dependent decay of interaction frequency ({% cite Lajoie2015 %}).
>
> The higher interaction between cis regions can be explained, at least in part, by the territorial organization of chromosomes in interphase (chromosome territories), and in a genome-wide contact map, this pattern appears as blocks of high interaction centered along the diagonal and matching individual chromosomes (fig. 12) ({% cite Cremer2010 %}, {% cite Lajoie2015 %}).
>
> ![Figure 12: Hi-C map](../../images/vgp_assembly/hic_map.png "An example of a Hi-C map. Genomic regions are arranged along the x and y axes, and contacts are colored on the matrix like a heat map; here darker color indicates greater interaction frequency.") {:width="10%"}
>   
> On the other hand, the distance-dependent decay may be due to random movement of the chromosomes, and in the contact map appears as a gradual decrease of the interaction frequency the farther away from the diagonal it moves ({% cite Lajoie2015 %}).
>
>
{: .comment}


> <hands-on-title>Generate a contact map with <b>PretextMap</b> and <b>Pretext Snapshot</b></hands-on-title>
>
> 1. {% tool [PretextMap](toolshed.g2.bx.psu.edu/repos/iuc/pretext_map/pretext_map/0.1.6+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input dataset in SAM or BAM format"*: `BAM Hi-C reads`
>    - *"Sort by"*: `Don't sort`
>
> 3. Rename the output as `PretextMap output`
>
> 2. {% tool [Pretext Snapshot](toolshed.g2.bx.psu.edu/repos/iuc/pretext_snapshot/pretext_snapshot/0.0.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input Pretext map file"*: `PretextMap output`
>    - *"Output image format"*: `png`
>    - *"Show grid?"*: `Yes`
{: .hands_on}

Let's have a look at the Hi-C contact maps generated by Pretext Snapshot.

![Figure 13: Pretext optical map](../../images/vgp_assembly/hic_map_pretext.png "Hi-C map generated by Pretext. Primary assembly full contact map generated in this training (a) Hi-C map representative of a typical missasembly (b).")

In the contact generated from the Bionano-scaffolded assembly can be identified 17 scaffolds, representing each of the haploid chromosomes of our genome (fig. 13.a). The fact that all the contact signals are found around the diagonal suggest that the contigs were scaffolded in the right order. However, during the assembly of complex genomes, it is common to find in the contact maps indicators of errors during the scaffolding process, as shown in the figure 13b. In that case, a contig belonging to the second chromosome has been misplaced as part of the fourth chromosome. We can also note that the final portion of the second chromosome should be placed at the beginning, as the off-diagonal contact signal suggests.

Once we have evaluated the quality of the scaffolded genome assembly, the next step consists in integrating the information contained in the HiC reads into our assembly, so that any errors identified can be resolved. For this purpose we will use SALSA2 ({% cite Ghurye2019 %}).  
     
### SALSA2 scaffolding

SALSA2 is an open source software that makes use of Hi-C to linearly orient and order assembled contigs along entire chromosomes ({% cite Ghurye2019 %}). One of the advantages of SALSA2 with respect to most existing Hi-C scaffolding tools is that it doesn't require the estimated number of chromosomes.

> <comment-title>SALSA2 algorithm overview</comment-title>
>
> Initially SALSA2 uses the physical coverage of Hi-C pairs to identify suspicious regions and break the sequence at the likely point of mis-assembly. Then, a hybrid scaffold graph is constructed using edges from the Hi-C reads, scoring the edges according to a *best buddy* scheme (fig. 14a).
>
> ![Figure 14: SALSA2 algorithm](../../images/vgp_assembly/salsa2_algorithm.png "Overview of the SALSA2 algorithm. Solid edges indicate the linkages between different contigs and dotted edges indicate the links between the ends of the same contig. B and E denote the start and end of contigs, respectively. Adapted from Ghurye et al. 2019.")
>
> From this graph scaffolds are iteratively constructed using a greedy weighted maximum matching. After each iteration, a mis-join detection step is performed to check if any of the joins made during this round are incorrect. Incorrect joins are broken and the edges blacklisted during subsequent iterations. This process continues until the majority of joins made in the prior iteration are incorrect. This provides a natural stopping condition, when accurate Hi-C links have been exhausted ({% cite Ghurye2019 %}).
>
{: .comment}

Before launching SALSA2, we need to carry out some modifications on our datasets.

> <hands-on-title>BAM to BED conversion</hands-on-title>
>
> 1. {% tool [bedtools BAM to BED](toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_bamtobed/2.30.0+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Convert the following BAM file to BED"*: `BAM Hi-C reads`
>    - *"What type of BED output would you like"*: `Create a full, 12-column "blocked" BED file`
>
> 2. Rename the output as `BED unsorted`
>
> 3. {% tool [Sort](sort1) %} with the following parameters:
>    - {% icon param-file %} *"Sort Dataset"*: `BED unsorted`
>    - *"on column"*: `Column: 4`
>    - *"with flavor"*: `Alphabetical sort`
>    - *"everything in"*: `Ascending order`
>
> 4. Rename the output as `BED sorted`
>
> 5. {% tool [Replace](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_find_and_replace/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `Primary assembly bionano`
>    - *"Find pattern"*: `:`
>    - *"Replace all occurrences of the pattern"*: `Yes`
>    - *"Find and Replace text in"*: `entire line`
>
> 6. Rename the output as `Primary assembly bionano edited`
{: .hands_on}

Now we can launch SALSA2 in order to generate the hybrid scaffolding based on the Hi-C data.

> <hands-on-title>Salsa scaffolding</hands-on-title>
>
>
> 1. {% tool [SALSA](toolshed.g2.bx.psu.edu/repos/iuc/salsa/salsa/2.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Initial assembly file"*: `Primary assembly bionano edited`
>    - {% icon param-file %} *"Bed alignment"*: `BED sorted`
>    - *"Restriction enzyme sequence(s)"*: `CTTAAG`
>
> 2. Rename the output as `SALSA2 scaffold FASTA` and `SALSA2 scaffold AGP`
>
{: .hands_on}

### Evaluate the final genome assembly with Pretext

Finally, we should repeat the procedure described previously for generating the contact maps, but in that case, we will use the scaffold generated by SALSA2. 

> <hands-on-title>Mapping reads against the scaffold</hands-on-title>
>
> 1. {% tool [BWA-MEM2](toolshed.g2.bx.psu.edu/repos/iuc/bwa_mem2/bwa_mem2/2.2.1+galaxy0) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `SALSA2 scaffold FASTA`
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-file %} *"Select fastq dataset"*: `Hi-C_dataset_F`
>    - *"Set read groups information?"*: `Do not set`
>    - *"Select analysis mode"*: `1.Simple Illumina mode`
>    - *"BAM sorting mode"*: `Sort by read names  (i.e., the QNAME field) `
>
> 2. Rename the output as `BAM forward SALSA2`
>
> 3. {% tool [BWA-MEM2](toolshed.g2.bx.psu.edu/repos/iuc/bwa_mem2/bwa_mem2/2.2.1+galaxy0) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `SALSA2 scaffold FASTA`
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-file %} *"Select fastq dataset"*: `Hi-C_dataset_R`
>    - *"Set read groups information?"*: `Do not set`
>    - *"Select analysis mode"*: `1.Simple Illumina mode`
>    - *"BAM sorting mode"*: `Sort by read names  (i.e., the QNAME field) `
>
> 4. Rename the output as `BAM reverse SALSA2`
>
> 5. {% tool [Filter and merge](toolshed.g2.bx.psu.edu/repos/iuc/bellerophon/bellerophon/1.0+galaxy0) %} chimeric reads from Arima Genomics with the following parameters:
>    - {% icon param-file %} *"First set of reads"*: `BAM forward SALSA2`
>    - {% icon param-file %} *"Second set of reads"*: `BAM reverse SALSA2`
>
> 6. Rename the output as `BAM Hi-C reads SALSA2`
>
> 7. {% tool [PretextMap](toolshed.g2.bx.psu.edu/repos/iuc/pretext_map/pretext_map/0.1.6+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input dataset in SAM or BAM format"*: `BAM Hi-C reads SALSA2`
>    - *"Sort by"*: `Don't sort`
>
> 8. Rename the output as `PretextMap output SALSA2`
>
> 9. {% tool [Pretext Snapshot](toolshed.g2.bx.psu.edu/repos/iuc/pretext_snapshot/pretext_snapshot/0.0.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input Pretext map file"*: `PretextMap output SALSA2`
>    - *"Output image format"*: `png`
>    - *"Show grid?"*: `Yes`
>
{: .hands_on}

In order to evaluate the Hi-C hybrid scaffolding, we are going to compare the contact maps before and after running SALSA2 (fig. 15). 
  
![Figure 15: Pretext final contact map](../../images/vgp_assembly/hi-c_pretext_final.png "Hi-C map generated by Pretext after the hybrid scaffolding based on Hi-C data. The red circles indicate the  differences between the contact map generated after (a) and before (b) Hi-C hybrid scaffolding.")

Among the most notable differences that can be identified between the contact maps, it can be highlighted the regions marked with red circles, where inversion can be identified.
        
# Conclusion

To sum up, it is worthwhile to compare the final assembly with the [_S. cerevisiae_ S288C reference genome](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_assembly_stats.txt).

![Table 1: Final stats](../../images/vgp_assembly/stats_conclusion.png "Comparison between the final assembly generating in this training and the reference genome. Contiguity plot using the reference genome size (a). Assemby statistics (b).")

With respect to the total sequence length, we can conclude that the size of our genome assembly is almost identical to the reference genome (fig.16a,b). Regarding the number of scaffolds, the obtained value is similar to the reference assembly, which consist in 16 chromosomes plus the mitochondrial DNA, which consists of 85,779 bp. The remaining statistics exhibit very similar values (fig. 16b).

![Figure 16: Comparison reference genome](../../images/vgp_assembly/hi-c_pretext_conclusion.png "Comparison between contact maps generated by using the final assembly (a) and the reference genome (b).")

If we compare the contact map of our assembled genome (fig. 17a) with the reference assembly (fig. 17b), we can see that the two are essentially identical. This means that we have achieved an almost perfect assembly at the chromosome level.
        

