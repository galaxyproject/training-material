---
layout: tutorial_hands_on

title: VGP assembly pipeline
zenodo_link: ''
enable: false
questions:
- "what combination of tools can produce the highest quality assembly of vertebrate genomes?"
- "How can we evaluate how good it is?"
objectives:
- "Learn the tools necessary to perform a de novo assembly of a vertebrate genome"
- "Evaluate the quality of the assembly"
time_estimation: '2h'
key_points:
- "The VGP pipeline allows to generate error-free, near gapless reference-quality genome assemblies"
- "The assembly can be divided in four main stages: genome profile analysis, HiFi long read phased assembly with hifiasm, Bionano hybrid scaffolding and Hi-C hybrid scaffolding"
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

Heterozygosity is also an important factor in genome assembly. Haplotype phasing, that is, the identification of alleles that are co-located on the same chromosome, has become a fundamental problem in heterozygous and polyploid genome assemblies ({% cite Zhang2020 %}). When there's not a reference sequence available, the *state-of-the--art* strategy consist in constructing a string graph with vertexes representing reads and edges representing consistent overlaps. In this kind of graph, after transitive reduction, heterozigous alleles in the string graph are represented by bubbles. When combined with Hi-C data, this approach allows complete diploid reconstruction ({% cite DominguezDelAngel2018 %}, {% cite Zhang2020 %}, {% cite Dida2021 %}).

The G10K consortium launched the Vertebrate Genomes Project (VGP), whose goal is generating high-quality, near-error-free, gap-free, chromosome-level, haplotype-phased, annotated reference genome assembly for each of the vertebrate species ({% cite Rhie2021 %}). Thhis tutorial will guide you step by step to assemble a high-quality referece genome by using the VGP assembly pipeline.


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

The tutorial is structured in four main sections:

- Genome profile analysis
- HiFi phased assembly with Hifiasm
- Hybrid scaffolding based using Bionano data
- Hi-C scaffolding

## Run the VGP workflows automatically

The pipeline proposed in this training is an adaption of the [current workflow versions](https://github.com/Delphine-L/iwc/tree/VGP/workflows/VGP-assembly-v2), whose purpuse is to explain each of the stages in with VGP assembly pipeline is structured. If you desire to run the *state-of-art* VGP pipelines, just follow the following instructions:

> ### {% icon hands_on %} Hands-on: Run workflow
>
> 1. Dowload the worflow files (whose extesion is *ga*) from this [GitHub repository](https://github.com/Delphine-L/iwc/tree/VGP/workflows/VGP-assembly-v2).
> 2. Click on Workflow on the top menu bar of Galaxy. You will see a list of all your workflows.
> 3. Click on the upload icon {% icon galaxy-upload %} at the top-right of the screen.
> ![Import Workflow](../../images/vgp_assembly/import_workflow.png "Import workflow from a file or URL.")
> 4. Provide your workflow:
>    - Option 1: Upload the workflow file in the box labelled “Archived Workflow File”
>    - Option 2: Paste the URL of the workflow into the box labelled “Archived Workflow URL”
> 5. Click the **import workflow** button.
>
{: .hands_on}

> ### {% icon comment %} Comments
> The Galaxy workflows include additional steps (e.g. parse parameter value) required for running it automatically, but which are not necessary when we run the pipeline step by step manually.
{: .comment}
    
## Background on datasets

To reduce compute times, we will use samples from the yeast _Saccharomyces cerevisiae_, one of the most intensively studied eukaryotic model organisms in molecular and cell biology. Yeast can be haploid or diploid, depending the stage of its life cycle. Both cell types are stable and can reproduce asexually by mitosis.

The VGP assembly pipeline uses data generated by a variery of technologies, including PacBio HiFi reads, Bionano optical maps, and Hi-C chromatin interaction maps.

PacBio HiFi reads rely on the Single Molecule Real-Time (SMRT) sequencing technology. SMRT is based on real-time imaging of fluorescently tagged nucleotides as they are added to a newly synthesized DNA strand. HiFi further combine multiple subreads from the same circular template to produce one highly accurate consensus sequence (fig. 3). This technology allows to generate long-read sequencing data with read lengths in the range of 10-25 kb and minimum read consensus accuracy  greater than 99% (Q20).

<!--
![fig2:PacBio sequencing technolgoy](../../images/vgp_assembly/pacbio_hifi.png "PacBio HiFi sequencing")
-->
    
Bionano technology relies on the isolation of kilobase-long DNA fragments, which are labeled at specific sequence motifs with a fluorescent dye, resulting in a unique fluorescent pattern for each genome. DNA molecules are stretched into nanoscale channels and imaged with a high-resolution camera, allowing to build optical maps based on the distance between motif-specific patterns ({% cite Lam2012 %}, {% cite Giani2020 %}).

Linkage information provided by optical maps is integrated with primary assembly sequences, and the overlaps are used to orient and order the contigs, resolve chimeric joins, and estimate the length of gaps between adjacent contigs (fig. 4).
    

![Bionano optical maps](../../images/vgp_assembly/bionano.png "Bionano optical maps")

    
Finally, the thigh-throughput chromosome conformation capture (Hi-C) technology is based on the capture of the chromatin three-dimensional. During Hi-C library preparation, DNA is crosslinked the chromatin in its 3D conformation. The crosslinked DNA is digested using restriction enzymes, and the digested ends are filled with biotinylated nucleotides. Next, the blunt ends of spatially proximal digested end are ligated. This provides contact information that can be used to reconstruct the proximity of genomic sequences belonging to the same chromosome ({% cite Giani2020 %}).

# Get data

Now we can start with the pipeline. The first step is to get the datasets from Zenodo:
    
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

    

## Primer removal from HiFi reads

Now, we will trim the residual adaptors sequences by using Cutadapt, in order remove the reads  that could interfere with the assembly process.

> ### {% icon hands_on %} Hands-on: Primer removal
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

Before starting a de novo genome assembly project, it useful to collect metrics on the properties of the genome under consideration, such as the expected genome size. Traditionally DNA flow citometry was considered the golden standard for estimating the genome size. Nowadays experimental methods have been replaced by computational approaches {% cite wang2020estimation %}. One of the widely used genome profilling methods is based on the analysis of k-mer frequencies. It allows to provide information not only about the genomic complexity, such as the genome size, levels of heterozygosity and repeat content, but also about the data quality.

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

Meryl will allow us to generate the *k*-mer profile by decomposing the sequencing data into *k*-lenght substrings, counting the ocurrence of each *k*-mer and determining its frequency. The original version of Meryl was developed for the Celera Assembler. The current Meryl version compises three main modules: one for generating *k*-mer databases, one for filtering and combining databases, and one for searching databases. *K*-mers are stored in lexicographical order in the database, similar to words in a dictionary ({% cite Rhie2020 %}).

> ### {% icon comment %} *k*-mer size estimation
>
>  Given an estimated genome size (G) and a tolerable collision rate (p), an appropriate k can be computed as k = log4 (G(1 − p)/p).
>
{: .comment}

> ### {% icon hands_on %} Hands-on: Generate *k*-mers count distribution
>
> 1. {% tool [Meryl](toolshed.g2.bx.psu.edu/repos/iuc/meryl/meryl/1.3+galaxy2) %} with the following parameters:
>    - *"Operation type selector"*: `Count operations`
>        - *"Count operations"*: `Count: count the ocurrences of canonical *k*-mers`
>        - {% icon param-collection %} *"Input sequences"*: `HiFi_collection (trim)`
>        - *"*k*-mer size selector"*: `Set a *k*-mer size`
>            - "**k*-mer size*": `21`
>
>    > ### {% icon comment %} Election of *k*-mer size
>    >
>    > We used 21 as *k*-mer size, as this length has demonstrated to be sufficiently long that most *k*-mers are not repetitive and is short enough to be more robust to sequencing errors. For very large (haploid size > 10 Gb) and/or very repetitive genomes, larger *k*-mer length is recommended to increase the number of unique *k*-mers. 
>    {: .comment}
>
> 2. Rename it `Collection meryldb`
>
> 3. {% tool [Meryl](toolshed.g2.bx.psu.edu/repos/iuc/meryl/meryl/1.3+galaxy1) %} with the following parameters:
>    - *"Operation type selector"*: `Operations on sets of *k*-mers`
>        - *"Operations on sets of *k*-mers"*: `Union-sum: return *k*-mers that occur in any input, set the count to the sum of the counts`
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

The next step is to infer the genome properties from the *k*-mer histogram generated by Meryl, for which we willl use GenomeScope2. Genomescope2 relies on a nonlinear least-squares optimization to fit a mixture of negative binomial distributions, generating estimated values for genome size, repetitiveness, and heterozygosity rates ({% cite RanalloBenavidez2020 %}).

> ### {% icon hands_on %} Hands-on: Estimate genome properties
>
> 1. {% tool [GenomeScope](toolshed.g2.bx.psu.edu/repos/iuc/genomescope/genomescope/2.0) %} with the following parameters:
>    - {% icon param-file %} *"Input histogram file"*: `Meryldb histogram`
>    - *"*k*-mer length used to calculate *k*-mer spectra"*: `21`
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
    - *Linear plot*: *k*-mer spectra and fitted models: frequency (y-axis) versus coverage.
    - *Log plot*: logarithmic transformation of the previous plot.
    - Transformed linear plot: *k*-mer spectra and fitted models: frequency times coverage (y-axis) versus coverage (x-axis). It allows to increases the heights of higher-order peaks, overcoming the effect of high heterozygosity.
    - Transformed log plot: logarithmic transformation of the previous plot.
- Model: this file includes a detailed report about the model fitting.
- Summary: it includes the properties infered from the model, such as genome haploid length and the percentage of heterozygosity.

Now, let's analyze the *k*-mer profiles, fitted models and estimated parameters:

![fig3:Genomescope plot](../../images/vgp_assembly/genomescope_plot.png "Genomescope2 plot")

This distribution is the result of the Poisson process underlying the generation of sequencing reads. As we can see, there is an unique peak centered around 28x, the modal *k*-mer coverage. The absence of a secondary peak at half diploid coverage is suggestive of the haploid nature of this genome, but could generally result also from very low heterozygosity. Low frequency *k*-mers are the result of sequencing errors.

Before proceeding to the next section, we need to carry out some operations on the output generated by GenomeScope2. The goal is to extract some parameters which at a later stage will be used by **purge_dups** ({% cite Guan2019 %}). The first relevant parameter is the `estimated genome size`.

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
> 6. Rename the output as `Estimated genome size`.
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
>    - {% icon param-file %} *"File to cut"*: `Parsing temporal output`
>    - *"Cut by"*: `fields`
>        - *"List of Fields"*: `Column 7`
>
> 2. Rename the output as `Transition parameter`
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


# HiFi phased assembly with hifiasm

Once we have done genome profiling stage, we can start the genome assembly with **hifiasm**,  a fast open-source de novo assembler specifically developed for PacBio HiFi reads.

## Genome assembly with **hifiasm**

One of the key advantages of hifiasm is that it allows to resolve near-identical, but not exacly identical sequences, such as repeats and segmental duplications ({% cite Cheng2021 %}).

> ### {% icon comment %} Hifiasm algorithm details
>
>By default hifiasm performs three rounds of haplotype-aware error correction to correct sequence errors but keeping heterozygous alleles. A position on the target read to be corrected is considered informative if there are two different nucleotides at that position in the alignment, and each allele is supported by at least tree reads.
>
> ![fig4:Hifiasm algorithm overview](../../images/vgp_assembly/hifiasm_algorithm.png "Hifiasm algorithm overview. Orange and blue bars represent the reads with heterozygous alleles carrying local phasing information, while green bars come from the homozygous regions without any heterozygous alleles.")
>
>Then, hifiasm builds a phased assembly string graph with local phasing information from the corrected reads. Only the reads coming from the same haplotype are connected in the phased assembly graph. After transitive reduction, a pair of heterozygous alleles is represented by a _bubble_ in the string graph. If there is no additional data, hifiasm arbitrarily selects one side of each bubble and outputs a primary assembly. In the case of a heterozygous genome, the primary assembly generated at this step may still retain haplotigs from the alternate allele.
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

Hifiasm generates four outputs in GFA format; this format is designed to represent genome variation, splice graphs in genes, and even overlaps between reads.


    
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

The VGP assembly pipeline contains several built-in QC steps, including QUAST, BUSCO, Merqury and Pretext. QUAST will generate summary statistics, BUSCO will search for universal single-copy ortholog genes, Merqury will evaluate assembly copy-numbers using k-mers, and Pretext will be used to evaluate the assembly contiguity.

> ### {% icon comment %} QUAST statistics
>
> QUAST will provide us with the following statistics:
> 
> - No. of contigs: The total number of contigs in the assembly.
> - Largest contig: The length of the largest contig in the assembly.
> - Total length: The total number of bases in the assembly.
> - Nx: The largest contig length, L, such that using contigs of length >= L accounts for at least x% of the bases of the assembly.
> - NGx: The contig length such that using equal or longer length contigs produces x% of the length of the reference genome, rather than x% of the assembly length.
> - GC content: the percentage of nitrogenous bases which are either guanine or cytosine.
>
{: .comment}
    
> ### {% icon hands_on %} Hands-on: assembly evaluation with QUAST
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

Let's have a look at the report generated by QUAST.

![fig5:QUAST plot](../../images/vgp_assembly/QUAST_initial.png "QUAST initial report.")

According with the plot, the primary assembly includes 40 contigs, whose accumulative length is around 12 Mbp. On the other hand, the accumulative length of the secondary assembly is much smaller, despite including more than 120 contigs. This is expected, since we are working with an haploid organism.

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

Next, we will use BUSCO, which will allow us to quantitative assessment of the completeness of a genome assembly in terms of expected gene content. It relies in the analysis of genes that should be present only once in a complete assembly or gene set, while allowing for rare gene duplications or losses ({% cite Simo2015 %}).

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

BUSCO generates two outputs, a tabular file which contains the complete results with scores, descriptions, lengths and coordinates of matches (fig. 8), and a short summary.

![fig8:BUSCO plot](../../images/vgp_assembly/BUSCO_full_table.png "BUSCO full table")

As we can see in the report, the results are simplified into four categories: *complete and single-copy*, *complete and duplicated*, *fragmented* and  *Missing BUSCOs*.
    
> ### {% icon question %} Questions
>
> 1. How many complete BUSCO genes have been identified?
> 2. How many BUSCOs gene are absent?
>
> > ### {% icon solution %} Solution
> >
> > 1. According the report, our assembly contains the complete sequence of  2121 complete BUSCO genes.
> > 2. 8 BUSCO genes are missing.
> > 
> {: .solution}
>
{: .question}

Despite BUSCO is robust for species that have been widely studied, it can be inaccurate when the newly assembled genome belong to taxas that are not well represented in the [OrtologDB](https://www.orthodb.org/). Merqury provides a complementary approach for assessing genome assembly quality metrics in a reference-free manner via k-mer copy number analysis.
    
> ### {% icon hands_on %} Hands-on: *k*-mer based evaluation with Merqury
>
> 1. {% tool [Merqury](toolshed.g2.bx.psu.edu/repos/iuc/merqury/merqury/1.3) %} with the following parameters:
>    - *"Evaluation mode"*: `Default mode`
>        - {% icon param-file %} *"k-mer counts database"*: `Merged meryldb`
>        - *"Number of assemblies"*: `Two assemblies
>            - {% icon param-file %} *"First genome assembly"*: `Primary contig FASTA`
>            - {% icon param-file %} *"Second genome assembly"*: `Alternate contig FASTA`    
>
{: .hands_on}

By default, Merqury generates three collections as output: stats, plots and QV stats.
    
The copy number spectum plot, known as *spectra-cn* plot (fig. 9a), tracks the multiplicity of each k-mer found in the Hi-Fi read set and colors it by the number of times it is found in a given assembly . Merqury connects the midpoint of each histogram bin with a line, giving the illusion of a smooth curve. Remember that the histogram dataset was generated previously by Meryl.

![fi9:Merqury plot](../../images/vgp_assembly/merqury_cn_plot.png "Merqury CN plot")

The black region in the left side corresponds to k-mers found only in the read set; it is usually indicative of sequencing error in the read set, although also can be indicative of missing sequences in the assembly. The read area represents one-copy k-mers in the genome, while blue area represents two-copy k-mers originating. from homozygous sequence or haplotype-specific duplications. From that figure we can state that the sequencing coverage is around 30x. 

The figure 1b corresponds to the assembly spectrum plot. This plot is helpful for measuring diploid assembly completeness as it shows the fraction of k-mers specific to both the primary and alternate assemblies ({% cite Rhie2020 %}).

# Post-assembly processing

An ideal haploid representation would consist of one allelic copy of all heterozygous regions in the two haplomes, as well as all hemizygous regions from both haplomes ({% cite Guan2019 %}). However, in high heterozygous genomes, assembly algorithms are frequently not able identify the highly divergent allelic sequences, resulting in the assembly of those regions as separate contigs, which can lead to issues in downstream analysis, such as scaffolding, gene annotation and read mapping in general ({% cite Small2007 %}, {% cite Guan2019 %}, {% cite Roach2018 %}). In order to solve this problem, we are going to use purge_dups; this tool will allow us to identify and reassign allelic contigs.

## Remove haplotypic duplication with **purge_dups**
        
This stage can consist in three substages: coverage related calculation, generation of all versus all self-alignment and resolution of haplotigs and overlaps (fig. 10).

![fig4:Post-processing step](../../images/vgp_assembly/purge_dupspipeline.png "Purge_dups pipeline. Adapted from github.com/dfguan/purge_dups.")

### Read-depth analysis

Initially we need to collapse the collection into a single dataset.
    
> ### {% icon hands_on %} Hands-on: Collapse the collection
>
> 1. {% tool [Collapse Collection](toolshed.g2.bx.psu.edu/repos/nml/collapse_collections/collapse_dataset/4.2) %} with the following parameters:
>    - {% icon param-collection %} *"Collection of files to collapse into single dataset"*:`HiFi_collection (trim)`
> 2. Rename de output as `HiFi reads collapsed`
{: .hands_on}

Now, we will map the reads against the primary assembly by using Minimap2 ({% cite Li2018 %}), an alignment program designed to map long sequences. 
   
> ### {% icon hands_on %} Hands-on: Map the reads to contigs with **Minimap2**
>
> 1. {% tool [Map with minimap2](toolshed.g2.bx.psu.edu/repos/iuc/minimap2/minimap2/2.17+galaxy4) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `Primary contig FASTA`
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-collection %} *"Select fastq dataset"*: `HiFi reads collapsed`
>        - *"Select a profile of preset options"*: `Long assembly to reference mapping (-k19 -w19 -A1 -B19 -O39,81 -E3,1 -s200 -z200 --min-occ-floor=100). Typically, the alignment will not extend to regions with 5% or higher sequence divergence. Only use this preset if the average divergence is far below 5%. (asm5)`
>    - In *"Set advanced output options"*:
>        - *"Select an output format"*: `paf`
>
> 2. Rename the output as `Reads mapped to contigs`
{: .hands_on}

Finally, we will use the `Reads mapped to contigs` PAF file for calculating some statistics required in a later stage. In this step, purge_dups initially produces a read-depth histogram from base-level coverages. This information is used for estimating the coverage cutoffs, taking on the fact that for collapsed haplotype contigs the reads from both alleles will map, whereas if the alleles have assembled as separate contigs the reads will be split over the two contigs, resulting in half the read-depth ({% cite Roach2018 %}). 
    
> ### {% icon hands_on %} Hands-on: Read-depth analsys
> 1. {% tool [purge_dups](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy3) %} with the following parameters:
>    - *"Function mode"*: `Calculate coverage cutoff, base-level read depth and create read depth histogram for PacBio data (calcuts+pbcstats)`
>        - {% icon param-file %} *"PAF input file"*: `Reads mapped to contigs`
>        - In *"Calcuts options"*:
>            - *"Upper bound for read depth"*: `63` (the previously estimated maximum depth)
>            - *"Ploidity"*: `Haploid`
>
>    > ### {% icon comment %} Comment
>    >
>    > In case you are working with a diploid organism, you should select `diploid` in the ploidy option.
>    {: .comment}
{: .hands_on}

Purge_dups generates three outputs:

- PBCSTAT base coverage: it contains the base-level coverage information.
- Calcuts-cutoff: it includes the thresholds calculated by purge_dups.
- Histagram plot.

As we can see in the histogram plot, the read depth distribution follows a simetric and unimodal distribution (fig. 11). This result is expected, since we are assemblying an haploid genome.
    
![fi9:purge_dups plot](../../images/vgp_assembly/purge_dups_plot.png "purge_dups histogram plot")

> ### {% icon details %} Expected results in diploid organisms
>
> For diploid organisms, it is frequent to observe bimodal distributions when allelic sequences from the same genomic region are assembled separately.
>
> ![fi9:purge_dups plot](../../images/vgp_assembly/purge_dups_plot2.png "purge_dups histogram generated in diploid genomes.")
>
{: .details}


### Generation of all versus all self-alignmnet

Now, we will segment the draft assembly into contigs by cutting at blocks of ‘N’s, and use minimap2 to generate an all by all self-alignment.

> ### {% icon hands_on %} Hands-on: purge_dups pipeline    
> 1. {% tool [purge_dups](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy2) %} with the following parameters:
>    - *"Function mode"*: `split assembly FASTA file by 'N's (split_fa)`
>        - {% icon param-file %} *"Assembly FASTA file"*: `Primary contig FASTA`
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
> 4. Rename the output as `Self-homology map`
{: .hands_on}

        
### Resolution of haplotigs and overlaps        

During the final step of the purge_dups pipeline, it will use the self alignments and the cutoffs for identifying the haplotypic duplications.

> ### {% icon comment %} Purge_dups algorithm details
>
> In order to identify the haplotypic duplications, purge_dups uses the  base-level coverage information to flag the contigs according the following criteria:
> - If more than 80% bases of a contig are above the high read depth cutoff or below the noise cutoff, it is discarted.
> - If more than 80% bases are in the diploid depth interval, it is labelled as a primary contig, otherwise it is considered further as a possible haplotig.
>
> Contigs that were flagged for further analysis according to read-depth are then evaluated to attempt to identify synteny with its allelic companion contig. In this step, purge_dups uses the information contained in the self alignments:
> - If the alignment score is larger than the cutoff *s* (default 70), the contig is marked for reassignment as haplotig. Contigs marked for reassignment with a maximum match score greater than the cutoff *m* (default 200) are further flagged as repetitive regions.
> 
> - Otherwise contigs are considered as a candidate primary contig.
>
> Once all matches associated with haplotigs have been removed from the self-alignment set, purge_dups ties consistent matches between the remaining candidates to find collinear matches, filtering all the matches whose score is less than the minimum chaining score *l*.
>
> Finally, purge-dups calculates the average coverage of the matching intervals for each overlap, and mark an unambiguous overlap as heterozygous when the average coverage on both contigs is less than the read-depth cutoff, removing the sequences corresponding to the matching interval in the shorter contig.
>
{: .comment}

> ### {% icon hands_on %} Hands-on: Resolution of haplotigs and overlaps
>    
> 1. {% tool [purge_dups](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy5) %} with the following parameters:
>    - *"Select the purge_dups function"*: `Purge haplotigs and overlaps for an assembly (purge_dups)`
>        - {% icon param-file %} *"PAF input file"*: `Self-homology map`
>        - {% icon param-file %} *"Base-level coverage file"*: `PBCSTAT base coverage` (output of the fifth step)
>        - {% icon param-file %} *"Cutoffs file"*: `calcuts cutoff` (output of the fifth step)
>
> 2. {% tool [purge_dups](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy2) %} with the following parameters:
>    - *"Select the purge_dups function"*: `Obtain sequences after purging (get_seqs)`
>        - {% icon param-file %} *"Assembly FASTA file"*: `Primary contig FASTA`
>        - {% icon param-file %} *"BED input file"*: `purge_dups BED` (output of the previous step)
> 
> 3. Rename the output `get_seq purged sequences` as `Primary contigs purged` and the `get_seq haplotype` file as `Alternate haplotype contigs`.
>
{: .hands_on}


### Process the alternative assembly: optional

Now we should repeat the same procedure with the alternate contigs generated by hifiasm.  In that case, we should start by merging the merging the `Purged sequences` generated in the previos step and the `Alternate contigs FASTA` file.
    
> ### {% icon hands_on %} Hands-on: Merge the purged sequences and the Alternate contigs
>
> 1. {% tool [Concatenate datasets](cat1) %} with the following parameters:
>    - {% icon param-file %} *"Concatenate Dataset"*: `Alternate contigs FASTA`
>    - In *"Dataset"*:
>        - {% icon param-repeat %} *"Insert Dataset"*
>            - {% icon param-file %} *"Select"*: `Purged sequences`
>
>    > ### {% icon comment %} Comment
>    >
>    > Remember that the `Purged sequences` file contains those contigs that were considered to be haplotypic duplications of the primary contigs.
>    {: .comment}
>
> 2. Rename the output as `Alternate contigs full`
>
{: .hands_on}

Once we have merged the files, we should run the pipeline again, but using the `Alternate contigs full` file as input.

> ### {% icon hands_on %} Hands-on: Process the alternate assembly with *purge_dups*
>
> 1. {% tool [Map with minimap2](toolshed.g2.bx.psu.edu/repos/iuc/minimap2/minimap2/2.17+galaxy4) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `Alternate contigs full`
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-collection %} *"Select fastq dataset"*: `HiFi reads collapsed`
>        - *"Select a profile of preset options"*: `Long assembly to reference mapping (-k19 -w19 -A1 -B19 -O39,81 -E3,1 -s200 -z200 --min-occ-floor=100). Typically, the alignment will not extend to regions with 5% or higher sequence divergence. Only use this preset if the average divergence is far below 5%. (asm5)`
>    - In *"Set advanced output options"*:
>        - *"Select an output format"*: `paf`
>
> 2. Rename the output as `Reads mapped to contigs alternate`
>
> 3. {% tool [purge_dups](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy3) %} with the following parameters:
>    - *"Function mode"*: `Calculate coverage cutoff, base-level read depth and create read depth histogram for PacBio data (calcuts+pbcstats)`
>        - {% icon param-file %} *"PAF input file"*: `Reads mapped to contigs alternate`
>        - In *"Calcuts options"*:
>            - *"Upper bound for read depth"*: `63`
>            - *"Ploidity"*: `Haploid`
>
>    > ### {% icon comment %} Comment
>    >
>    > In case you are working with a diploid organism, you should select `diploid`.
>    {: .comment}
>
> 4. {% tool [purge_dups](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy2) %} with the following parameters:
>    - *"Function mode"*: `split assembly FASTA file by 'N's (split_fa)`
>        - {% icon param-file %} *"Assembly FASTA file"*: `Alternate contig full`
>
> 2. Rename the output as `Split FASTA alternate`
>
> 3. {% tool [Map with minimap2](toolshed.g2.bx.psu.edu/repos/iuc/minimap2/minimap2/2.17+galaxy4) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `Split FASTA alternate`
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-file %} *"Select fastq dataset"*: `Split FASTA alternate`
>        - *"Select a profile of preset options"*: `Construct a self-homology map - use the same genome as query and reference (-DP -k19 -w 19 -m200) (self-homology)`
>    - In *"Set advanced output options"*:
>        - *"Select an output format"*: `PAF`
> 
> 4. Rename the output as `Self-homology map alternate`
>        
> 5. {% tool [purge_dups](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy5) %} with the following parameters:
>    - *"Select the purge_dups function"*: `Purge haplotigs and overlaps for an assembly (purge_dups)`
>        - {% icon param-file %} *"PAF input file"*: `Self-homology map`
>        - {% icon param-file %} *"Base-level coverage file"*: `PBCSTAT base coverage alternate`
>        - {% icon param-file %} *"Cutoffs file"*: `calcuts cutoff alternate`
>
> 6. {% tool [purge_dups](toolshed.g2.bx.psu.edu/repos/iuc/purge_dups/purge_dups/1.2.5+galaxy2) %} with the following parameters:
>    - *"Select the purge_dups function"*: `Obtain sequences after purging (get_seqs)`
>        - {% icon param-file %} *"Assembly FASTA file"*: `Alternate contig full`
>        - {% icon param-file %} *"BED input file"*: `purge_dups BED alternate`
> 
> 3. Rename the outputs as `Purged sequences` and `Alternate haplotype`.
>
{: .hands_on}


## Second round of assembly evaluation

Once we have run purge_dups, we can evaluate assembly again, and compare the results before and after purging. In that case, we will only evaluate the primary assembly.

> ### {% icon hands_on %} Hands-on: assembly evaluation with QUAST
>
> 1. {% tool [Quast](toolshed.g2.bx.psu.edu/repos/iuc/quast/quast/5.0.2+galaxy1) %} with the following parameters:
>    - *"Use customized names for the input files?"*: `Yes, specify custom names`
>    - In *"1. Contigs/scaffolds"*:
>        - {% icon param-file %} *"Contigs/scaffolds file"*: `Primary contigs purged`
>        - *"Name"*: `Primary assembly`
>    - *"Reads options"*: `Pacbio SMRT reads`
>        - {% icon param-collection %} *"FASTQ file"*: `HiFi reads collapsed`
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

In the figure 13 we can see the accumulative length plot corresponding to the primary assembly after being processed with purge_dups.

![fig4:Post-processing step](../../images/vgp_assembly/QUAST_cummulative.png "QUAST accumutive plot")

> ### {% icon question %} Questions
>
> 1. How many contigs include the primary assembly after being processed with purge_dups? And before?
> 2. Which percentage of reads mapped to the primary assembly?
>
> > ### {% icon solution %} Solution
> >
> > 1. The primary assembly includes 37 contigs after being processed. The original assembly included 40 contigs.
> > 2. According the report, 100% of reads mapped to the primary assembly. 
> > 
> {: .solution}
>
{: .question}

Now, let's evaluate the assemby with BUSCO.

> ### {% icon hands_on %} Hands-on: assessing assembly completness with BUSCO
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.0.0+galaxy0) %} with the following parameters:
>    - {% icon param-files %} *"Sequences to analyse"*: `Primary contigs purged`
>    - *"Mode"*: `Genome assemblies (DNA)`
>        - *"Use Augustus instead of Metaeuk"*: `Use Metaeuk`
>    - *"Auto-detect or select lineage?"*: `Select lineage`
>       - *"Lineage"*: `Saccharomycetes`
>    - In *"Advanced Options"*:
>        - *"Which outputs should be generated"*: `short summary text` and `summary image`
>
> 2. Rename the summary as `BUSCO second report`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> How many complete BUSCO genes have been identified?
>
> > ### {% icon solution %} Solution
> >
> > According the report, the primary assembly contains 2021 complete BUSCO genes.
> > 
> {: .solution}
>
{: .question}

# Hybrid scaffolding based using Bionano data


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

In this section we map Hi-C reads to scaffold the genome assembly. During Hi-C library preparation, the DNA is fragmented and re-ligated, leading to region of the genome that are closer in the 3D space being prefentially joined. Each DNA fragment is then sequenced from each end of this artificial junction, generating read pairs, with most of the contacts in the kpb range, bu with many contacts in the Mbp range as well. This information can be used to reconstruct the order and orientation of the contigs/scaffolds generated in the previous steps.
    

## Pre-processing Hi-C data

Even though Hi-C generated paired-end reads, we need to map each read separately. The reason is that most aligners assumes that the ends of a single continuous genomic fragment are being sequenced, and the distance between these two ends
fits a known distribution, but in Hi-C data,  the insert size of the ligation product can vary between
1bp to hundreds of megabases ({% cite Lajoie2015 %}).

> ### {% icon hands_on %} Hands-on: Mapping Hi-C reads
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
> 2. Rename the output as `BAM forward`
>
> 3. {% tool [Map with BWA-MEM](toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.17.2) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `output` (Input dataset)
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-file %} *"Select fastq dataset"*: `output` (Input dataset)
>    - *"Set read groups information?"*: `Do not set`
>    - *"Select analysis mode"*: `1.Simple Illumina mode`
>    - *"BAM sorting mode"*: `Sort by read names  (i.e., the QNAME field) `
>
> 4. Rename the output as `BAM reverse`
>
> 5. {% tool [Filter and merge](toolshed.g2.bx.psu.edu/repos/iuc/bellerophon/bellerophon/1.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"First set of reads"*: `BAM forward` 
>    - {% icon param-file %} *"Second set of reads"*: `Bam reverse`
>
> 6. Rename it as `BAM Hi-C reads`
{: .hands_on}

Finally, we need to convert the BAM file to BED format, and sorting it.

> ### {% icon hands_on %} Hands-on: BAM to BED conversion
>
> 1. {% tool [bedtools BAM to BED](toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_bamtobed/2.30.0+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Convert the following BAM file to BED"*: `outfile` (output of **Filter and merge** {% icon tool %})
>    - *"What type of BED output would you like"*: `Create a full, 12-column "blocked" BED file`
>
> 2. {% tool [Sort](sort1) %} with the following parameters:
>    - {% icon param-file %} *"Sort Dataset"*: `output` (output of **bedtools BAM to BED** {% icon tool %})
>    - *"on column"*: `c4`
>    - *"with flavor"*: `Alphabetical sort`
>    - *"everything in"*: `Ascending order`
>
{: .hands_on}


## Generate Hi-C contact map

Most of the paired reads from HiC will map to the same  contigs. In a typic Hi-C contact map, contigs are ordered by size and are evaluated against each other in a triangular matrix. Most contacts should appear close to the diagonal (mapping to self). Off-diagonal signal is indicative of chromosomal interactions in non-chromosome level assemblies, or of assembly mis-joints. We will now generate a Hi-C contact map before scaffolding the assembly to compare the Hi-C contact map after scaffolding.

> ### {% icon hands_on %} Hands-on: Generate a contact map with PretextMap
>
> 1. {% tool [PretextMap](toolshed.g2.bx.psu.edu/repos/iuc/pretext_map/pretext_map/0.1.6+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input dataset in SAM or BAM format"*: `outfile` (output of **Filter and merge** {% icon tool %})
>    - *"Sort by"*: `Don't sort`
>
> 2. {% tool [Pretext Snapshot](toolshed.g2.bx.psu.edu/repos/iuc/pretext_snapshot/pretext_snapshot/0.0.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input Pretext map file"*: `pretext_map_out` (output of **PretextMap** {% icon tool %})
>    - *"Output image format"*: `png`
>    - *"Show grid?"*: `Yes`
>
{: .hands_on}

***TODO***: explain the output here. What does it mean. What does this show about our data/assembly so far (e.g. do the contigs look fairly well ordered, or not). 


## Salsa scaffolding

> ### {% icon hands_on %} Hands-on: Salsa scaffolding
>
> 1. {% tool [Replace](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_find_and_replace/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `output` (Input dataset)
>    - *"Find pattern"*: `:`
>    - *"Replace all occurences of the pattern"*: `Yes`
>    - *"Find and Replace text in"*: `entire line`
>
> 2. {% tool [SALSA](toolshed.g2.bx.psu.edu/repos/iuc/salsa/salsa/2.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Initial assembly file"*: `outfile` (output of **Replace** {% icon tool %})
>    - {% icon param-file %} *"Bed alignment"*: `out_file1` (output of **Sort** {% icon tool %})
>    - {% icon param-file %} *"Sequence graphs"*: `output` (Input dataset)
>    - *"Restriction enzyme sequence(s)"*: add the enzyme sequence(s) here
>
{: .hands_on}


## Evaluate the Salsa scaffolding results

Now, the scaffolded assembly will be evaluated using BUSCO and QUAST.

> ### {% icon hands_on %} Hands-on: Evaluation with BUSCO
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

> ### {% icon hands_on %} Hands-on: Evaluation with QUAST
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


## Generate a Hi-C contact map after Salsa2 scaffolding

Now, we repeat the producedure described previously for generating the optical maps, but in that case, we will use the scaffold generated by Salsa2.
    
> ### {% icon hands_on %} Hands-on: Mapping reads against the scaffold
>
> 1. {% tool [Map with BWA-MEM](toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.17.2) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `scaffolds_fasta` (output of **SALSA** {% icon tool %})
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-file %} *"Select fastq dataset"*: `output` (Input dataset)
>    - *"Set read groups information?"*: `Do not set`
>    - *"Select analysis mode"*: `1.Simple Illumina mode`
>
> 2. {% tool [Map with BWA-MEM](toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.17.2) %} with the following parameters:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>        - {% icon param-file %} *"Use the following dataset as the reference sequence"*: `scaffolds_fasta` (output of **SALSA** {% icon tool %})
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-file %} *"Select fastq dataset"*: `output` (Input dataset)
>    - *"Set read groups information?"*: `Do not set`
>    - *"Select analysis mode"*: `1.Simple Illumina mode`
>
> 3. {% tool [Filter and merge](toolshed.g2.bx.psu.edu/repos/iuc/bellerophon/bellerophon/1.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"First set of reads"*: `bam_output` (output of **Map with BWA-MEM** {% icon tool %})
>    - {% icon param-file %} *"Second set of reads"*: `bam_output` (output of **Map with BWA-MEM** {% icon tool %})
>
>
> 4. {% tool [PretextMap](toolshed.g2.bx.psu.edu/repos/iuc/pretext_map/pretext_map/0.1.6+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input dataset in SAM or BAM format"*: `outfile` (output of **Filter and merge** {% icon tool %})
>    - *"Sort by"*: `Don't sort`
>
> 5. {% tool [Pretext Snapshot](toolshed.g2.bx.psu.edu/repos/iuc/pretext_snapshot/pretext_snapshot/0.0.3+galaxy0) %} with the following parameters:
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
