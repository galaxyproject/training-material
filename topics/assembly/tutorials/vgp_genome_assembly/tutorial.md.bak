---
layout: tutorial_hands_on

title: VGP assembly pipeline
zenodo_link: ''
enable: false
questions:
- what combination of tools can produce the highest quality assembly of vertebrate genomes?
- How can we evaluate how good it is? 
objectives:
- Learn the tools necessary to perform a de novo assembly of a vertebrate genome
- Evaluate the quality of the assembly
time_estimation: '2h'
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- delphine-l
- astrovsky01
- gallardoalba

---


# Introduction
{:.no_toc}

Advances in sequencing technologies over the last few decades have revolutionised the field of genomics, allowing for a reduction in both the time and resources required to *de novo* genome assembly. Until recently, second-generation sequencing technologies (also known as Next Generation Sequencing or NGS) allowed to produce highly accurate but short (up to 800bp), whose extension was not long enough to cope with the difficulties associated with repetitive regions. Today, so-called third-generation sequencing (TGS) technologies, usually known as single-molecule real-time (SMRT) sequencing, have become dominant in de novo assembly of large genomes. TGS can use native DNA without amplication, reducing sequencing error and bias ({% cite Hon2020 %}, {% cite Giani2020 %}). Very recently, Pacific Biosciences introducedHigh-Fidelity (HiFi) sequencing, which produces reads 10-20 kpb in length with a minimum accuracy of 99% (Q20). In this tutorial you will use HiFi reads in combination with data from additional sequencing technologies to generate a high-quality reference genome assembly.

Deciphering the structural organisation of complex vertebrate genomes is currently one of the most challenges in genomics ({% cite Frenkel2012 %}). Despite the significant progress made in recent years, a key question remains: what combination of data and tools can produce the highest quality assembly? In order to adequately answer it, it is necessary to analyse two of the main factors that determine the difficulty of genome assembly processes: repetitive content and heterozigosity.

Repetitive elements can be grouped into two categories: interspersed repeats, such as transposable elements (TE) that occur at multiple loci throughout the genome, and tandem repeats (TR), that occur at a single locus ({% cite Trresen2019 %}). Repetitive elements are an important component of eukariotyc genomes, constituting over a third of the genome in the case of mammals ({% cite SoteroCaio2017 %}, {% cite Chalopin2015 %}). In the case of tamdem repeats, various estimates suggest that they are present in at least one third of human protein sequences ({% cite Marcotte1999 %}). TE content is among the main factors contributing to the lack of continuity in the reconstruction of genomes, especially in the case of large ones, as its content is highly correlated with genome size ({% cite SoteroCaio2017 %}). On the other hand, TR usually lead to local genome assembly collapse, especially when their length is close to that of the reads ({% cite Trresen2019 %}).

Heterozygosity is also an important factor in genome assembly. Haplotype phasing, that is, the identification of alleles that are co-located on the same chromosome, has become a fundamental problem in heterozygous and polyploid genome assemblies ({% cite Zhang2020 %}). A common strategy to overcome these difficulties is to remap genomes to a single haplotype, which represents the whole genome. This approach is useful for highly inbred samples that are nearly homozygous, but when applied to highly heterozygous genomes, such as aquatic organism, it missses potential differences in sequence, structure, and gene presence, usually leading to ambiguties and redundancies in the initial contig-level assemblies ({% cite DominguezDelAngel2018 %}, {% cite Zhang2020 %}).

To address these problems, the G10K consortium launched the Vertebrate Genomes Project (VGP), whose goal is generating high-quality, near-error-free, gap-free, chromosome-level, haplotype-phased, annotated reference genome assembly for each of the vertebrate species currently present on planet Earth ({% cite Rhie2021 %}). The protocol proposed in this tutorial, the VGP assembly pipeline, is the result of years of study and analysis of the available tools and data sources.


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# VGP assembly pipeline overview

The figure 1 represents the VGP assembly pipeline. 

![fig1:VGP pipeline](../../images/vgp_assembly/VGP_Pipeline.png "VPG Pipeline 2.0")

In order to facilitate the development of the workflow, we will structure it in four main sections:

- Genome profile analysis
- HiFi long read phased assembly with Hifiasm
- Hybrid scaffolding based on phased assembly and Bionano data
- Hybrid scaffolding based on a phased assembly and Hi-C mapping data

**TODO:** suggest including here something about how the Galaxy workflow has additional steps (e.g. parse parameter value), so that it can run automatically, but if run step by step in this tutorial, then some of those steps are just manually entered value (e.g. genome size parameter)

## Background on datasets

In order to reduce processing times, we will use samples from _Saccharomyces cerevisiae_, one of the most intensively studied eukaryotic model organisms in molecular and cell biology. This organisms can be haploid or diploid, depending the stage of its life cycle. Both cell types are stable and can reproduce asexually by mitosis.

The VGP assembly pipeline requires datasets generated by three different technologies: PacBio HiFi reads, Bionano optical maps, and Hi-C chromatin interaction maps.

PacBio HiFi reads rely on the Single Molecule Real-Time (SMRT) sequencing technology. It is based on real-time imaging of fluorescently tagged nucleotides as they are synthesized along individual DNA template molecules, combining multiple subreads of the same circular template using a statistical model to produce one highly accurate consensus sequence, along with base quality values (figure 2). This technology allows to generate long-read sequencing dataseets with read lengths averaging 10-25 kb and accuracies greater than 99.5%.

![fig2:PacBio sequencing technolgoy](../../images/vgp_assembly/pacbio_hifi.png "PacBio HiFi sequencing")

Optical genome mapping is a method for detecting structural variants. The generation of Bionano optical maps starts with high molecular weight DNA, which is labeled at specific sequences motif with a fluorescent dye, resulting in a unique fluorescence pattern for each individual genome. The comparison of the labelled fragments among different samples enables the detection of structural variants. Optical maps are integrated with the primary assemby sequence in order to identify and correct potential chimeric joints, and estimate the gap sizes.

The high-throughput chromosome conformation capture (Hi-C) technology is based on the capture of the chromatin conformation, enabling the identification of topological domains. Hi-C chromatin interaction maps methods first crosslink the chromatin in its 3D conformation. The crosslinked DNA is digested using restriction enzymes, and the digested ends are filled with biotinylated nucleotides. Next, the blunt ends  of spatially proximal digested end are ligated, preserving the chromosome interaction regions. Finally, the DNA is purified to assure that only fragments originating from ligation events are sequenced. 

# Get data

> ### {% icon hands_on %} Hands-on: Data upload
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
>       - Clich `Add Definition` button and select `Name Tag`: column `D`
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
>   SRR13577846_1    https://zenodo.org/record/5550653/files/SRR13577846_1.30x.wgaps.fastq.gz?download=1  fastqsanger.gz    HiFi  HiFi_collection
>   SRR13577846_2    https://zenodo.org/record/5550653/files/SRR13577846_2.30x.wgaps.fastq.gz?download=1  fastqsanger.gz    HiFi  HiFi_collection
>   SRR13577846_3    https://zenodo.org/record/5550653/files/SRR13577846_3.30x.wgaps.fastq.gz?download=1  fastqsanger.gz    HiFi  HiFi_collection
>       ```
>
>    - From **Rules** menu select `Add / Modify Column Definitions`
>       - Click `Add Definition` button and select `List Identifier(s)`: column `A`
>       - Click `Add Definition` button and select `URL`: column `B`
>       - Click `Add Definition` button and select `Type`: column `C`
>       - Clich `Add Definition` button and select `Group Tag`: column `D`
>       - Clich `Add Definition` button and select `Collection Name`: column `E`
>    - Click `Apply` and press <kbd>Upload</kbd>
>
{: .hands_on}


# Data quality assessment

To begin our analysis we will carry out the evaluation and pre-processing of our data. In order to identify potential anomalies in the data, we are going to use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), an open-source tool that provides a simple way to quality control raw sequence data.

> ### {% icon hands_on %} Hands-on: Quality check
>
> 1. Run **FastQC** {% icon tool %} with the following parameters
>    - {% icon param-collection %} *"Raw read data from your current history"*: `HiFi_collection`
>
> 2. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.8+galaxy1) %} with the following parameters:
>    - In *"Results"*:
>      - *"Which tool was used generate logs?"*: `FastQC`
>      - {% icon param-collection %} *"Dataset collection"*: select the `FastQC on collection:Raw Data` dataset.
>    - In *"Report title"*: `HiFi quality report`
> 3. Click on the {% icon galaxy-eye %} (eye) icon and inspect the generated HTML file
>
{: .hands_on}

![fig3:HiFi Quality report](../../images/vgp_assembly/quality_plot.png "PacBio HiFi qualiry report")

As we can see, the mean Phred score is over 80 in all the samples, which means that the base call accuracy is around 99.999999%!


> ### {% icon comment %} Comments
> For more information on the topic of quality control, please see our training materials [here](https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/quality-control/tutorial.html).
{: .comment}

According the quality report, less that 0.1% of the reads include adaptor sequences. Despide of this, we will trim the residual adaptors sequences by using Cutadapt, in order to avoid potential reads which could interfer in the subsequent steps.

> ### {% icon hands_on %} Hands-on: Optional step: primer removal
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
>        - *"Match times"*: `3`
>        - *"Maximum error rate"*: `0.1`
>        - *"Minimum overlap length"*: `35`
>        - *"Look for adapters in the reverse complement"*: `True`
>    - In *"Filter Options"*:
>        - *"Discard Trimmed Reads"*: `Yes`
>
> 2. Rename the output file as `HiFi_collection (trim)`. To rename an output file, click on the result, and then click again on the title to change it. After closing, you may need to refresh the history to see the name change.
>
{: .hands_on}

# Genome profile analysis

An important step before starting a de novo genome assembly project is to proceed with the analysis of the genome profile. Determining these characteristics in advance has the potential to reveal whether an analysis is not reflecting the full complexity of the genome, for example, if the number of variants is underestimated or a significant fraction of the genome is not assembled ({% cite Vurture2017 %}).

Traditionally DNA flow citometry was considered the golden standart for estimating the genome size, one of the most important factors to determine the required coverage level. However, nowadays experimental methods have been replaced by computational approaches {% cite wang2020estimation %}. One of the most widely used procedures for undertaking genomic profiling is the analyis of k-mer frequencies. It allows to provide information not only about the genomic complexity, such as the genome size, levels of heterozygosity and repeat content, but also about the data quality. In addition, k-mer spectra analysis can be used in a reference-free manner for assessing genome assembly quality metrics ({% cite Rhie2020 %}).

> ### {% icon details %} K-mer size, sequencing coverage and genome size
>
>K-mers are unique substrings of length k contained within a DNA sequence. For example, the DNA sequence *TCGATCACA* can be decomposed into six unique k-mers that have five bases long: *TCGAT*, *CGATC*, *GATCA*, *ATCAC* and *TCACA*. A sequence of length L will have  L-k+1 k-mers. On the other hand, the number of possible k-mers can be calculated as  n<sup>k</sup>, where n is number of possible monomers and k is the k-mer size.
>
>
>---------| -------------|-----------------------
>  Bases  |  K-mer size  |  Total possible k-mers
>---------| -------------|-----------------------
>    4    |       1      |            4          
>    4    |       2      |           16          
>    4    |       3      |           64          
>    4    |       4      |          256          
>    4    |      ...     |          ...          
>    4    |      10      |    1.048.576         
>---------|--------------|-----------------------
>
> Thus, the k-mer size is a key parameter, which must be large enough to map  uniquely to the genome, but not too large, since it can lead to wasting computational resources. In the case of the human genome, k-mers of 31 bases in length lead to 96.96% of unique k-mers.
>
>Each unique k-mer can be assigned a value for coverage based on the number of times it occurs in a sequence, whose distribution will approximate a Poisson distribution, with the peak corresponding to the average genome sequencing depth. From the genome coverage, the genome size can be easily computed.
{: .details}

    
In section we will use two basic tools to computationally estimate the genome features: Meryl and GenomeScope.

## Generation of k-mer spectra with **Meryl**

Meryl will allow us to perform the k-mer profiling by decomposing the sequencing data into k-lenght substrings and determining its frequency. The original version was developed for use in the Celera Assembler, and it comprises three modules: one for generating k-mer databases, one for filtering and combining databases, and one for searching databases. The k-mer database is stored in sorted order, similar to words in a dictionary ({% cite Rhie2020 %}).

> ### {% icon comment %} K-mer size estimation
>
>  Given an estimated genome size (G) and a tolerable collision rate (p), an appropriate k can be computed as k = log4 (G(1 − p)/p).
>
{: .comment}

> ### {% icon hands_on %} Hands-on: Generate k-mers count distribution
>
> 1. {% tool [Meryl](toolshed.g2.bx.psu.edu/repos/iuc/meryl/meryl/1.3+galaxy2) %} with the following parameters:
>    - *"Operation type selector"*: `Count operations`
>        - *"Count operations"*: `Count: count the ocurrences of canonical k-mers`
>        - {% icon param-collection %} *"Input sequences"*: `HiFi_collection (trim)`
>        - *"K-mer size selector"*: `Set a k-mer size`
>            - "*K-mer size*": `21`
>
>    > ### {% icon comment %} Election of k-mer size
>    >
>    > We used 21 as k-mer size, as this length has demonstrated to be sufficiently long that most k-mers are not repetitive and is short enough that the analysis will be more robust to sequencing errors. For extremely large (haploid size over 10 Gb) and/or very repetitive genomes, it is recommended to use larger k-mer lengths to increase the number of unique k-mers. 
>    {: .comment}
>
> 2. Rename it `Collection meryldb`
>
> 3. {% tool [Meryl](toolshed.g2.bx.psu.edu/repos/iuc/meryl/meryl/1.3+galaxy1) %} with the following parameters:
>    - *"Operation type selector"*: `Operations on sets of k-mers`
>        - *"Operations on sets of k-mers"*: `Union-sum: return k-mers that occur in any input, set the count to the sum of the counts`
>        - {% icon param-file %} *"Input meryldb"*: `Collection meryldb`
>
> 4. Rename it as `Merged meryldb`    
>
> 5. {% tool [Meryl](toolshed.g2.bx.psu.edu/repos/iuc/meryl/meryl/1.3+galaxy0) %} with the following parameters:
>    - *"Operation type selector"*: `Generate histogram dataset`
>        - {% icon param-file %} *"Input meryldb"*: `Merged meryldb`
>
> 6. Finally, rename it as `Meryldb histogram`.
>
{: .hands_on}


## Genome profiling with **GenomeScope2**

The next step is to computationally infer the genome properties from the k-mer count distribution generated by Meryl, for which we'll use GenomeScope2. It relies in a nonlinear least-squares optimization to fit a mixture of negative binomial distributions, generating estimated values for genome size, repetitiveness, and heterozygosity rates ({% cite RanalloBenavidez2020 %}).

> ### {% icon hands_on %} Hands-on: Estimate genome properties
>
> 1. {% tool [GenomeScope](toolshed.g2.bx.psu.edu/repos/iuc/genomescope/genomescope/2.0) %} with the following parameters:
>    - {% icon param-file %} *"Input histogram file"*: `Meryldb histogram`
>    - *"K-mer length used to calculate k-mer spectra"*: `21`
>
>   - In "*Output options*": mark `Summary of the analysis`
>   - In "*Advanced options*":
>       - *"Create testing.tsv file with model parameters"*: `true`
>
>    {: .comment}
>
{: .hands_on}

Genomescope will generate six outputs:
    
- Plots
    - *Linear plot*: K-mer spectra and fitted models: frequency (y-axis) versus coverage.
    - *Log plot*: logarithmic transformation of the previous plot.
    - Transformed linear plot: K-mer spectra and fitted models: frequency times coverage (y-axis) versus coverage (x-axis). It allows to increases the heights of higher-order peaks, overcoming the effect of high heterozygosity.
    - Transformed log plot: logarithmic transformation of the previous plot.
- Model: this file includes a detailed report about the model fitting.
- Summary: it includes the properties infered from the model, such as genome haploid length and the percentage of heterozygosity.

Now, let's analyze the k-mer profiles, fitted models and estimated parameters:

![fig3:Genomescope plot](../../images/vgp_assembly/genomescope_plot.png "Genomescope2 plot")

As we can see, there is an unique peak centered around 28, which is the coverage with the highest number of different 21-mers. According the normal-like k-mer spectra, we can infer that it is a haploid genome. The large number of unique k-mers on the left side with frequence around one is due to error during the sequencing process.

Before jumping to the next section, we need to carry out some operation on the output generated by Genomescope2. The goal is to generate some parameters which at a later stage will be used by **purge_dups** ({% cite Guan2019 %}). Lets start with the `estimated genome size`.

> ### {% icon hands_on %} Hands-on: Get estimated genome size
>
> 1. {% tool [Replace](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_find_and_replace/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `summary` (output of **GenomeScope** {% icon tool %})
>    - *"Find pattern"*: `bp`
>    - *"Replace all occurences of the pattern"*: `Yes`
>    - *"Find and Replace text in"*: `entire line`
>
> 2. {% tool [Replace](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_find_and_replace/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: output file of **Replace** {% icon tool %})
>    - *"Find pattern"*: `,`
>    - *"Replace all occurences of the pattern"*: `Yes`
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
> 6. {% tool [Parse parameter value](param_value_from_file) %}(param_value_from_file) with the following parameters:
>    - {% icon param-file %} *"Input file containing parameter to parse out of"*: output of **Advanced Cut** {% icon tool %}
>    - *"Select type of parameter to parse"*: `Integer`
>
> 7. Rename the output as `Estimated genome size`.
>
> > ### {% icon question %} Questions
> >
> > What is the estimated genome size?
> >
> > > ### {% icon solution %} Solution
> > >
> > > The estimated genome size is 12664060 bp.
> > >
> > {: .solution}
> >
> {: .question}
> 
{: .hands_on}

Now let's parse the `upper bound for the read depth estimation` parameter.
       
> ### {% icon hands_on %} Hands-on: Get maximum read depth
>
> 1. {% tool [Compute an expression on every row](toolshed.g2.bx.psu.edu/repos/devteam/column_maker/Add_a_column1/1.6) %} with the following parameters:
>    - *"Add expression"*: `1.5*c3`
>    - {% icon param-file %} *"as a new column to"*: `model_params` (output of **GenomeScope** {% icon tool %})
>    - *"Round result?"*: `Yes`
>    - *"Input has a header line with column names?"*: `No`
>
> 2. {% tool [Compute an expression on every row](toolshed.g2.bx.psu.edu/repos/devteam/column_maker/Add_a_column1/1.6) %} with the following parameters:
>    - *"Add expression"*: `3*c7`
>    - {% icon param-file %} *"as a new column to"*: output of **Compute** {% icon tool %})
>    - *"Round result?"*: `Yes`
>    - *"Input has a header line with column names?"*: `No`
>
> 3. Rename it as `Parsing temporal output`
>
> 4. {% tool [Advanced Cut](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/1.1.0) %} with the following parameters:
>    - {% icon param-file %} *"File to cut"*: `Parsing temporal output` (output of **Compute** {% icon tool %})
>    - *"Cut by"*: `fields`
>        - *"List of Fields"*: `Column 8`
>
> 5. {% tool [Parse parameter value](param_value_from_file) %} with the following parameters:
>    - {% icon param-file %} *"Input file containing parameter to parse out of"*: output of **Advanced Cut** {% icon tool %}
>    - *"Select type of parameter to parse"*: `Integer`
>
> 6. Rename it as `Maximum depth`
>
> > ### {% icon question %} Questions
> >
> > What is the estimated maximum depth?
> >
> > > ### {% icon solution %} Solution
> > >
> > > The estimated maximum depth is  63 reads.
> > >
> > {: .solution}
> >
> {: .question}
>
{: .hands_on}

Finally, let's parse the `transition between haploid and diploid coverage depths` parameter.

> ### {% icon hands_on %} Hands-on: Get transition parameter        
> 1. {% tool [Advanced Cut](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/1.1.0) %} with the following parameters:
>    - {% icon param-file %} *"File to cut"*: `Parsing temporal output` (output of **Compute** {% icon tool %})
>    - *"Cut by"*: `fields`
>        - *"List of Fields"*: `Column 7`
>
> 2. {% tool [Parse parameter value](param_value_from_file) %} with the following parameters:
>    - {% icon param-file %} *"Input file containing parameter to parse out of"*: output of **Advanced Cut** {% icon tool %}
>    - *"Select type of parameter to parse"*: `Integer`
>
> 3. Rename it as `Transition parameter`
>
> > ### {% icon question %} Questions
> >
> > What is the estimated transition parameter?
> >
> > > ### {% icon solution %} Solution
> > >
> > > The estimated transition parameter is  21 reads.
> > >
> > {: .solution}
> >
> {: .question}
>
{: .hands_on}


# HiFi long read phased assembly with hifiasm

Once we have finished the genome profiling stage, we can start the genome assembly with **hifiasm**,  a fast open-source de novo assembler specifically developed for PacBio HiFi reads.

## Genome assembly with **hifiasm**

One of the key focus of hifiams is to different copies of a segmental duplication involving a single segregating site, allowing to resolve near-identical, but not exactly identical, repeats and segmental duplications ({% cite Cheng2021 %}).

> ### {% icon comment %} Hifiasm algorithm details
>
>By default hifiasm performs three rounds of haplotype-aware error correction to correct sequence errors but keeping heterozygous alleles. A position on the target read to be corrected is considered informative if there are tow different nucleotides at the position of the alignment, and each type is supported by at least tree reads.
>
> ![fig4:Hifiasm algorithm overview](../../images/vgp_assembly/hifiasm_algorithm.png "Hifiasm algorithm overview. Orange and blue bars represent the reads with heterozygous alleles carrying local phasing information, while green bars come from the homozygous regions without any heterozygous alleles.")
>
>Then, hifiasm builds a phased assembly string graph with local phasing information from the corrected reads. Only the reads coming from the same haplotype are connected in the phased assembly graph. After transitive reduction, a pair of heterozygous alleles will be represented by a _bubble_ in the string graph. If there are no additional data, hifiasm arbitrarily selects one side of each bubble and outputs a primary assembly. For a heterozygous genome, the primary assembly generated at this step may still contain haplotigs from more than one homologous haplotype.
>
>
{: .comment}

> ### {% icon hands_on %} Hands-on: Phased assembly with **hifiasm**
>
> 1. {% tool [Hifiasm](toolshed.g2.bx.psu.edu/repos/bgruening/hifiasm/hifiasm/0.14+galaxy0) %} with the following parameters:
>    - *"Assembly mode"*: `Standard`
>        - {% icon param-file %} *"Input reads"*: `HiFi_collection (trim)` (output of **Cutadapt** {% icon tool %})
>    - *"Options for purging duplicates"*: `Specify`
>       - *"Coverage upper bound"*: `63` (maximum depth previously obtained)
>    - *"Options for Hi-C partition"*: `Specify`
>       - *"Hi-C R1 reads"*: `Hi-C_dataset_F`
>       - *"Hi-C R2 reads"*: `Hi-C_dataset_R`
>
> 2. Rename the `Hi-C hap1 contig graph` as `Primary contig graph` and add a `#primary` tag
> 3. Rename the `Hi-C hap2 contig graph` as `Alternate contig graph` and  add a `#alternate` tag
>
{: .hands_on}

Hifiasm generates four outputs in GFA format; this format is specially designed to capture sequence graphs as the product of an assembly, a representation of variation in genomes, splice graphs in genes, or even overlap between reads from long-read sequencing technology.

## Convert GFA format to FASTA with **GFA to FASTA** 

We have obtained the fully phased contig graphs of the primary and alternate haplotypes, but the output format of **hifiasm** is not adequate for the subsequent steps, so we will convert them into fasta format.

> ### {% icon hands_on %} Hands-on: convert GFA to FASTA
>
> 1. {% tool [GFA to FASTA](toolshed.g2.bx.psu.edu/repos/iuc/gfa_to_fa/gfa_to_fa/0.1.2) %} with the following parameters:
>    - {% icon param-files %} *"Input GFA file"*: select `Primary contig graph` and the `Alternate contig graph` datasets
>
> 2. Rename the outputs as `Primary contig FASTA` and `Alternate contig FASTA`
>
{: .hands_on}

## Initial assembly evaluation

Once generated the draft assembly, it is a good idea to evaluate its quality. 

> ### {% icon hands_on %} Hands-on: assembly evaluation with Quast
>
> 1. {% tool [Quast](toolshed.g2.bx.psu.edu/repos/iuc/quast/quast/5.0.2+galaxy1) %} with the following parameters:
>    - *"Use customized names for the input files?"*: `Yes, specify custom names`
>    - In *"1. Contigs/scaffolds"*:
>        - {% icon param-file %} *"Contigs/scaffolds file"*: `Primary contig FASTA`
>        - *"Name"*: `Primary assembly`
>    - Click in *"Insert Contigs/scaffolds"*
>    - In *"2. Contigs/scaffolds"*:
>        - {% icon param-file %} *"Contigs/scaffolds file"*: `Alternate contig FASTA`
>        - *"Name"*: `Alternate assembly`
>    - *"Reads options"*: `Pacbio SMRT reads`
>        - {% icon param-collection %} *"FASTQ file"*: `HiFi collection (trim)`
>    - *"Type of assembly"*: `Genome`
>        - *"Use a reference genome?"*: `No`
>            - *"Estimated reference genome size (in bp) for computing NGx statistics"*: `12664060` (previously estimated)
>        - *"Type of organism"*: `Eukaryote: use of GeneMark-ES for gene finding, Barrnap for ribosomal RNA genes prediction, BUSCO for conserved orthologs finding (--eukaryote)`
>    - *"Is genome large (>100Mpb)?"*: `No`
>
>    > ### {% icon comment %} Comment
>    >
>    > Remember that for this training we are using _S. cerevisiae_, a reduced genome. In the case of assembling a vertebrate genome, you must select `yes` in the previous option.
>    {: .comment}
>
> 2. Rename the HTML report as `QUAST initial report`
>
{: .hands_on}

Let's have a look at the HTML report.

![fig5:QUAST plot](../../images/vgp_assembly/QUAST_initial.png "Quast initial report.")

> ### {% icon question %} Questions
>
> 1. What is the longest contig in the primary assembly? And in the alternate one?
> 2. What is the N50 of the primary assembly?
> 3. Which percentage of reads mapped to each assembly? 
>
> > ### {% icon solution %} Solution
> >
> > 1. The longest contig in the primary assembly is 914.549 bp, and 15.845 bp in the alternate assembly.
> > 2. The N50 of the primary assembly is 425.706 bp.
> > 3. According the report, 100% of reads mapped to the primary assembly, but only around 57% mapped to the alternate assembly.
> > 
> {: .solution}
>
{: .question}

> ### {% icon hands_on %} Hands-on: assessing assembly completness with BUSCO
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.0.0+galaxy0) %} with the following parameters:
>    - {% icon param-files %} *"Sequences to analyse"*: `Primary contig FASTA` and `Alternate contig FASTA`
>    - *"Mode"*: `Genome assemblies (DNA)`
>        - *"Use Augustus instead of Metaeuk"*: `Use Metaeuk`
>    - *"Auto-detect or select lineage?"*: `Select lineage`
>       - *"Lineage"*: `Saccharomycetes`
>    - *"Which outputs should be generated"*: `short summary text`
>
>    > ### {% icon comment %} Comment
>    >
>    > Remember to modify the lineage option if you are working with vertebrate genomes.
>    {: .comment}
>
> 2. Rename the summary as `BUSCO initial report`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. Which percentage of Benchmarking Universal Single-Copy Orthologs (BUSCO) genes have been identified?
> 2. How many BUSCOs gene are absent?
>
> > ### {% icon solution %} Solution
> >
> > 1. According the report, our assembly contains the complete sequence of  99.3% of BUSCO genes.
> > 2. 8 BUSCO genes are missing.
> > 
> {: .solution}
>
{: .question}


> ### {% icon hands_on %} Hands-on: K-mer based evaluation with Merqury
>
> 1. {% tool [Merqury](toolshed.g2.bx.psu.edu/repos/iuc/merqury/merqury/1.3) %} with the following parameters:
>    - *"Evaluation mode"*: `Default mode`
>        - {% icon param-file %} *"K-mer counts database"*: `Merged meryldb`
>        - *"Number of assemblies"*: `Two assemblies
>            - {% icon param-file %} *"First genome assembly"*: `Primary contig FASTA`
>            - {% icon param-file %} *"Second genome assembly"*: `Alternate contig FASTA`    
>
{: .hands_on}
    
    
# Post-assembly processing

An ideal haploid representation would consist of one allelic copy of all heterozygous regions in the two haplomes (haplotype contigs), as well as all hemizygous regions from both haplomes. However, the allelic relationship between haplotypes still present a problem for *de novo* genome assembly, specially in high heterozygous genomes; sequence divergence between pair of allelic sequences can lead to assemble there regions as separate contigs, rather than the expected single haplotype-fused contig. It can result in assemblies signicantly larger than the haploid genome size, which can lead to interferences in downstream stages, such as scaffolding and gene annotation ({% cite Guan2019 %}, {% cite Roach2018 %}). 

Usually, allelic relationships are inferred at the post-assembly stage. Despite the haplotig processing requites multiple steps, the approach used in this tutorial can be summaryzed in two steps: firstly we will identify the syntenic contigs by using the mapped read coverage and **minimap2** ({% cite Li2018 %}) alignments. Then, we will resolve the haplotigs and overlaps in the primary assembly by using **purge_dups**.

## Remove haplotypic duplication with **purge_dups**

This step includes 11 steps, summarized in the following scheme:

![fig4:Post-processing step](../../images/vgp_assembly/purge_dupspipeline.png "Purge_dups pipeline")


> ### {% icon hands_on %} Hands-on: purge_dups pipeline
>
> 1. {% tool [Collapse Collection](toolshed.g2.bx.psu.edu/repos/nml/collapse_collections/collapse_dataset/4.2) %} with the following parameters:
>    - {% icon param-collection %} *"Collection of files to collapse into single dataset"*:`HiFi_collection (trim)`
> 
> 2. Rename de output as `HiFi reads collapsed`
>
> 3. {% tool [Map with minimap2](toolshed.g2.bx.psu.edu/repos/iuc/minimap2/minimap2/2.17+galaxy4) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `Primary contig FASTA`
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-collection %} *"Select fastq dataset"*: `HiFi_collection (trim)` (output of **Cutadapt** {% icon tool %})
>        - *"Select a profile of preset options"*: `Long assembly to reference mapping (-k19 -w19 -A1 -B19 -O39,81 -E3,1 -s200 -z200 --min-occ-floor=100). Typically, the alignment will not extend to regions with 5% or higher sequence divergence. Only use this preset if the average divergence is far below 5%. (asm5)`
>    - In *"Set advanced output options"*:
>        - *"Select an output format"*: `paf`
>
> 4. Rename the output as `Reads mapped to contigs`
> 
> 5. {% tool [purge_dups](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy3) %} with the following parameters:
>    - *"Function mode"*: `Calculate coverage cutoff, base-level read depth and create read depth histogram for PacBio data (calcuts+pbcstats)`
>        - {% icon param-file %} *"PAF input file"*: `Reads mapped to contigs`
>        - In *"Calcuts options"*:
>            - *"Upper bound for read depth"*: `63` (the previously estimated maximum depth)
>            - *"Ploidity"*: `Haploid`
>
>    > ### {% icon comment %} Comment
>    >
>    > In the case you are working with a diploid orgasm, you should select `diploid` in the ploidity option.
>    > It will generate three outputs: the base-level coverage file (PBCSTAT base coverage), the cutoff file (calcuts cutoff) and a histogram plot.
>    {: .comment}
>
> 6. {% tool [purge_dups](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy2) %} with the following parameters:
>    - *"Function mode"*: `split assembly FASTA file by 'N's (split_fa)`
>        - {% icon param-file %} *"Assembly FASTA file"*: `Primary contig FASTA`
>
> 7. Rename the output as `Split FASTA`
>
> 8. {% tool [Map with minimap2](toolshed.g2.bx.psu.edu/repos/iuc/minimap2/minimap2/2.17+galaxy4) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `Split FASTA`
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-file %} *"Select fastq dataset"*: `Split FASTA`
>        - *"Select a profile of preset options"*: `Construct a self-homology map - use the same genome as query and reference (-DP -k19 -w 19 -m200) (self-homology)`
>    - In *"Set advanced output options"*:
>        - *"Select an output format"*: `PAF`
> 
> 9. Rename the output as `Self-homology map`
>
> 10. {% tool [purge_dups](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy5) %} with the following parameters:
>    - *"Select the purge_dups function"*: `Purge haplotigs and overlaps for an assembly (purge_dups)`
>        - {% icon param-file %} *"PAF input file"*: `Self-homology map`
>        - {% icon param-file %} *"Base-level coverage file"*: `PBCSTAT base coverage` (output of the fifth step)
>        - {% icon param-file %} *"Cutoffs file"*: `calcuts cutoff` (output of the fifth step)
>
> 11. {% tool [purge_dups](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy2) %} with the following parameters:
>    - *"Select the purge_dups function"*: `Obtain sequences after purging (get_seqs)`
>        - {% icon param-file %} *"Assembly FASTA file"*: `Primary contig FASTA`
>        - {% icon param-file %} *"BED input file"*: `purge_dups BED` (output of the previous step)
>
>
{: .hands_on}

## Second assembly evaluation assembly evaluation

Once we have purged the duplications, let's evaluate the assembly again. 

> ### {% icon hands_on %} Hands-on: assembly evaluation with Quast
>
> 1. {% tool [Quast](toolshed.g2.bx.psu.edu/repos/iuc/quast/quast/5.0.2+galaxy1) %} with the following parameters:
>    - *"Use customized names for the input files?"*: `Yes, specify custom names`
>    - In *"1. Contigs/scaffolds"*:
>        - {% icon param-file %} *"Contigs/scaffolds file"*: `Primary contig FASTA`
>        - *"Name"*: `Primary assembly`
>    - Click in *"Insert Contigs/scaffolds"*
>    - In *"2. Contigs/scaffolds"*:
>        - {% icon param-file %} *"Contigs/scaffolds file"*: `Alternate contig FASTA`
>        - *"Name"*: `Alternate assembly`
>    - *"Reads options"*: `Pacbio SMRT reads`
>        - {% icon param-collection %} *"FASTQ file"*: `HiFi collection (trim)`
>    - *"Type of assembly"*: `Genome`
>        - *"Use a reference genome?"*: `No`
>            - *"Estimated reference genome size (in bp) for computing NGx statistics"*: `12664060` (previously estimated)
>        - *"Type of organism"*: `Eukaryote: use of GeneMark-ES for gene finding, Barrnap for ribosomal RNA genes prediction, BUSCO for conserved orthologs finding (--eukaryote)`
>    - *"Is genome large (>100Mpb)?"*: `No`
>
>
> 2. Rename the HTML report as `QUAST second report`
>
{: .hands_on}


> ### {% icon hands_on %} Hands-on: assessing assembly completness with BUSCO
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.0.0+galaxy0) %} with the following parameters:
>    - {% icon param-files %} *"Sequences to analyse"*: `Primary contig FASTA` and `Alternate contig FASTA`
>    - *"Mode"*: `Genome assemblies (DNA)`
>        - *"Use Augustus instead of Metaeuk"*: `Use Metaeuk`
>    - *"Auto-detect or select lineage?"*: `Select lineage`
>       - *"Lineage"*: `Saccharomycetes`
>    - In *"Advanced Options"*:
>        - *"Which outputs should be generated"*: `short summary text`
>
>    > ### {% icon comment %} Comment
>    >
>    > Remember to modify the lineage option if you are working with vertebrate genomes.
>    {: .comment}
>
> 2. Rename the summary as `BUSCO initial report`
>
{: .hands_on}
    
----

<!--

Bibliography https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3409271/

We segment the input draft assembly into contigs by cutting at blocks ‘N’s, and use minimap2 to generate an all by all self-alignment.

We next recognize and remove haplotigs in essentially the same way as purge_haplotigs, and remove all matches associated with haplotigs from the self-alignment set.

Finally we chain consistent matches in the remainder to find overlaps, then calculate the average coverage of the matching intervals for each overlap, and mark an unambiguous overlap as heterozygous when the average coverage on both contigs is less than the read depth cutoff found in step 1, removing the sequence corresponding to the matching interval in the shorter contig.

purge_dups can significantly improve genome assemblies by removing overlaps and haplotigs caused by sequence divergence in heterozygous regions. This both removes false duplications in primary draft assemblies while retaining completeness and sequence integrity, and can improve scaffolding. 

Along with sequence similarity, purge_dups and purge_haplotigs take into account the coverage depth obtained by mapping short or long reads to the contigs. Coverage depth represents the number of reads covering a position in a contig (computed after mapping reads on the assembly). The contigs are then aligned to select duplicates accurately and remove them. While purge_dups sets its coverage thresholds automatically, purge_haplotigs requires user-provided values.

-->

## Sub-step with **Concatenate datasets**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Concatenate datasets](cat1) %} with the following parameters:
>    - {% icon param-file %} *"Concatenate Dataset"*: `get_seqs_hap` (output of **Purge overlaps** {% icon tool %})
>    - In *"Dataset"*:
>        - {% icon param-repeat %} *"Insert Dataset"*
>            - {% icon param-file %} *"Select"*: `out_fa` (output of **GFA to FASTA** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}


# Hybrid scaffolding based on phased assembly and Bionano data


## Sub-step with **Bionano Hybrid Scaffold**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Bionano Hybrid Scaffold](toolshed.g2.bx.psu.edu/repos/bgruening/bionano_scaffold/bionano_scaffold/3.6.1+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"NGS FASTA"*: `output` (Input dataset)
>    - {% icon param-file %} *"BioNano CMAP"*: `output` (Input dataset)
>    - *"Configuration mode"*: `VGP mode`
>    - {% icon param-file %} *"Conflict resolution file"*: `output` (Input dataset)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Concatenate datasets**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Concatenate datasets](cat1) %} with the following parameters:
>    - {% icon param-file %} *"Concatenate Dataset"*: `ngs_contigs_scaffold_trimmed` (output of **Bionano Hybrid Scaffold** {% icon tool %})
>    - In *"Dataset"*:
>        - {% icon param-repeat %} *"Insert Dataset"*
>            - {% icon param-file %} *"Select"*: `ngs_contigs_not_scaffolded_trimmed` (output of **Bionano Hybrid Scaffold** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Parse parameter value**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Parse parameter value](param_value_from_file) %} with the following parameters:
>    - {% icon param-file %} *"Input file containing parameter to parse out of"*: `output` (Input dataset)
>    - *"Select type of parameter to parse"*: `Integer`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

## Sub-step with **Quast**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Quast](toolshed.g2.bx.psu.edu/repos/iuc/quast/quast/5.0.2+galaxy1) %} with the following parameters:
>    - *"Use customized names for the input files?"*: `No, use dataset names`
>        - {% icon param-file %} *"Contigs/scaffolds file"*: `out_file1` (output of **Concatenate datasets** {% icon tool %})
>    - *"Type of assembly"*: `Genome`
>        - *"Use a reference genome?"*: `No`
>            - *"Estimated reference genome size (in bp) for computing NGx statistics"*: `{'id': 7, 'output_name': 'integer_param'}`
>        - *"Type of organism"*: `Eukaryote (--eukaryote): use of GeneMark-ES for gene finding, Barrnap for ribosomal RNA genes prediction, BUSCO for conserved orthologs finding`
>    - In *"Genes"*:
>        - *"Tool for gene prediction"*: `Don't predict genes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

# Hybrid scaffolding based on a phased assembly and HiC mapping data

***TODO***: need to re-name a lot of the inputs and outputs here. They have been auto-generated from the workflow but I think we want people to be able to run this step by step. I've taken out some of the steps that are "parse parameter value" etc. 

In this section we map HiC reads to scaffold the genome assembly. In HiC sequencing, parts of the genome that are close together are artificially joined. A DNA fragment is then sequenced from each end of this artificial junction, giving a read pair. If reads from this read pair map to two contigs, it indicates that those two contigs are close together in the genome. A good short video showing the HiC process is here: https://youtu.be/-MxEw3IXUWU

**TODO**: add image here of HiC

Inputs required for this section:

* An assembly FASTA file. This can be the output of the phased assembly section, and/or the output of the Bionano scaffolding section. 
* HiC reads, one set of forward reads and one set of reverse reads. If there is more than one set of Hi-C pair-read datasets, concatenate all the forward reads into one file, and the reverse reads into another file, in the same order.
* Genome size estimate: we can get this from an earlier step using GenomeScope. This is the haploid length. 

A summary of the 5 steps in this section: 

* Map the HiC reads to the assembly 
* View a contact map of HiC reads against the assembly, before scaffolding
* Scaffold the assembly with HiC reads: using the assembly file, and the mapped HiC reads
* Evaluate the scaffolding results: use Busco and Quast 
* View a contact map of the HiC reads after scaffolding

Outputs from this section:

* A scaffolded assembly FASTA file
* contact maps of HiC reads pre- and post scaffolding
* post-scaffolding reports from Busco and Quast 


A simplified image of the workflow for this section (not showing Quast and Busco):


![Hic-wf-summary](../../images/vgp_assembly/hic-wf-summary.png "HiC workflow")




## 1. Map the HiC reads to the assembly

We will do this separately for the forward and reverse set of HiC reads. We have to do this separately because these are not standard paired-end reads. 

**Map the forward HiC reads:**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Map with BWA-MEM](toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.17.2) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `output` (Input dataset)
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-file %} *"Select fastq dataset"*: `output` (Input dataset)
>    - *"Set read groups information?"*: `Do not set`
>    - *"Select analysis mode"*: `1.Simple Illumina mode`
>    - *"BAM sorting mode"*: `Sort by read names  (i.e., the QNAME field) `
>
{: .hands_on}

**Map the reverse HiC reads:**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Map with BWA-MEM](toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.17.2) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `output` (Input dataset)
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-file %} *"Select fastq dataset"*: `output` (Input dataset)
>    - *"Set read groups information?"*: `Do not set`
>    - *"Select analysis mode"*: `1.Simple Illumina mode`
>    - *"BAM sorting mode"*: `Sort by read names  (i.e., the QNAME field) `
>
{: .hands_on}


**Merge the mapped reads:**

Now we will merge these two BAM files:

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Filter and merge](toolshed.g2.bx.psu.edu/repos/iuc/bellerophon/bellerophon/1.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"First set of reads"*: `bam_output` (output of **Map with BWA-MEM** {% icon tool %})
>    - {% icon param-file %} *"Second set of reads"*: `bam_output` (output of **Map with BWA-MEM** {% icon tool %})
>
{: .hands_on}

**Convert the mapped BAM file to a BED file:**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [bedtools BAM to BED](toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_bamtobed/2.30.0+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Convert the following BAM file to BED"*: `outfile` (output of **Filter and merge** {% icon tool %})
>    - *"What type of BED output would you like"*: `Create a full, 12-column "blocked" BED file`
>
{: .hands_on}

**Sort the BED file:**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Sort](sort1) %} with the following parameters:
>    - {% icon param-file %} *"Sort Dataset"*: `output` (output of **bedtools BAM to BED** {% icon tool %})
>    - *"on column"*: `c4`
>    - *"with flavor"*: `Alphabetical sort`
>    - *"everything in"*: `Ascending order`
>
{: .hands_on}


## 2. View a contact map of the mapped HiC reads

Most of the paired reads from HiC will map to the same (or nearby) contigs. On a graph, with ordered contigs on each axis, a lot of the contacts will be along the diagonal (mapping to self), or nearby (around that diagonal line). But some may be in odd places - for example, showing a lot of reads mapped to both contig 4 and contig 19. We will now generate a contact map of the assembly before it is scaffolded, to compare to the contact map after scaffolding.

**Generate a contact map:**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [PretextMap](toolshed.g2.bx.psu.edu/repos/iuc/pretext_map/pretext_map/0.1.6+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input dataset in SAM or BAM format"*: `outfile` (output of **Filter and merge** {% icon tool %})
>    - *"Sort by"*: `Don't sort`
>
{: .hands_on}

**Convert the map to an image:**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Pretext Snapshot](toolshed.g2.bx.psu.edu/repos/iuc/pretext_snapshot/pretext_snapshot/0.0.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input Pretext map file"*: `pretext_map_out` (output of **PretextMap** {% icon tool %})
>    - *"Output image format"*: `png`
>    - *"Show grid?"*: `Yes`
>
{: .hands_on}

***TODO***: explain the output here. What does it mean. What does this show about our data/assembly so far (e.g. do the contigs look fairly well ordered, or not). 


## 3. Salsa scaffolding

Files required: The assembly file (optional: and the assembly graph), the sorted BED file, and the restriction enzyme sequence from the HiC sequencing. If you are using VGP GenomeArk data, you can get this information from the same file as the HiC reads, in a file called re_bases.txt.

**Prepare the assembly file:**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Replace](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_find_and_replace/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `output` (Input dataset)
>    - *"Find pattern"*: `:`
>    - *"Replace all occurences of the pattern"*: `Yes`
>    - *"Find and Replace text in"*: `entire line`
>
{: .hands_on}


**SALSA scaffolding:**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [SALSA](toolshed.g2.bx.psu.edu/repos/iuc/salsa/salsa/2.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Initial assembly file"*: `outfile` (output of **Replace** {% icon tool %})
>    - {% icon param-file %} *"Bed alignment"*: `out_file1` (output of **Sort** {% icon tool %})
>    - {% icon param-file %} *"Sequence graphs"*: `output` (Input dataset)
>    - *"Restriction enzyme sequence(s)"*: add the enzyme sequence(s) here
>
{: .hands_on}


## 4. Evaluate the Salsa scaffolding results

The scaffolded assembly fasta file can then be analysed in Busco and Quast.

**Busco:**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.2.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `scaffolds_fasta` (output of **SALSA** {% icon tool %})
>    - *"Mode"*: `Genome assemblies (DNA)`
>        - *"Use Augustus instead of Metaeuk"*: `Use Metaeuk`
>    - *"Lineage"*: ``
>    - In *"Advanced Options"*:
>        - *"Which outputs should be generated"*: ``
>
{: .hands_on}


There are four outputs: short summary, summary as an image, and two tables (full results and missing buscos). 

***TODO***: explain what these outputs mean; are the results "good" ?

**Quast:**

Inputs required for Quast: scaffolded assembly file from Salsa, and estimated genome size. The estimated genome size is obtained from an earlier step with GenomeScope.

Run Quast:

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Quast](toolshed.g2.bx.psu.edu/repos/iuc/quast/quast/5.0.2+galaxy1) %} with the following parameters:
>    - *"Use customized names for the input files?"*: `No, use dataset names`
>        - {% icon param-file %} *"Contigs/scaffolds file"*: `scaffolds_fasta` (output of **SALSA** {% icon tool %})
>    - *"Type of assembly"*: `Genome`
>        - *"Use a reference genome?"*: `No`
>            - *"Estimated reference genome size (in bp) for computing NGx statistics"*: `enter estimated genome size`
>        - *"Type of organism"*: `Eukaryote (--eukaryote): use of GeneMark-ES for gene finding, Barrnap for ribosomal RNA genes prediction, BUSCO for conserved orthologs finding`
>    - *"Is genome large (> 100 Mbp)?"*: `Yes`
>    - In *"Genes"*:
>        - *"Tool for gene prediction"*: `Don't predict genes`
>
{: .hands_on}


There are four outputs: the Quast report in three formats, and a log file. 

***TODO***: explain what these outputs mean; are the results "good" ?


## 5. Generate a post-scaffolding contact map

There are five steps: 

* Map the forward HiC reads to the scaffolded assembly
* Map the reverse HiC reads to the scaffolded assembly
* Combine these bam files into a single file
* Generate a contact map
* Conver the map to an image

**Map the forward HiC reads:**


> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Map with BWA-MEM](toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.17.2) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `scaffolds_fasta` (output of **SALSA** {% icon tool %})
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-file %} *"Select fastq dataset"*: `output` (Input dataset)
>    - *"Set read groups information?"*: `Do not set`
>    - *"Select analysis mode"*: `1.Simple Illumina mode`
>
{: .hands_on}

**Map the reverse HiC reads:**


> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Map with BWA-MEM](toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.17.2) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `scaffolds_fasta` (output of **SALSA** {% icon tool %})
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-file %} *"Select fastq dataset"*: `output` (Input dataset)
>    - *"Set read groups information?"*: `Do not set`
>    - *"Select analysis mode"*: `1.Simple Illumina mode`
>
{: .hands_on}

**Merge the mapped reads:**


> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Filter and merge](toolshed.g2.bx.psu.edu/repos/iuc/bellerophon/bellerophon/1.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"First set of reads"*: `bam_output` (output of **Map with BWA-MEM** {% icon tool %})
>    - {% icon param-file %} *"Second set of reads"*: `bam_output` (output of **Map with BWA-MEM** {% icon tool %})
>
>
{: .hands_on}

**Generate a contact map:**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [PretextMap](toolshed.g2.bx.psu.edu/repos/iuc/pretext_map/pretext_map/0.1.6+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input dataset in SAM or BAM format"*: `outfile` (output of **Filter and merge** {% icon tool %})
>    - *"Sort by"*: `Don't sort`
>
>
{: .hands_on}

**Convert the map to an image:**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Pretext Snapshot](toolshed.g2.bx.psu.edu/repos/iuc/pretext_snapshot/pretext_snapshot/0.0.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input Pretext map file"*: `pretext_map_out` (output of **PretextMap** {% icon tool %})
>    - *"Output image format"*: `png`
>    - *"Show grid?"*: `Yes`
>
>
{: .hands_on}


***TODO***: explain the output here. What does the pretext map show. How does it compare to the pre-scaffolding map. 

***TODO***: overall, explain what the scaffolding section results mean. What are the next possible steps.


***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}



***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}



# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
