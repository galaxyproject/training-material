---
layout: tutorial_hands_on

title: Genome-wide isoform switching analysis
zenodo_link: ''
questions:
- Which isoform switching events are preemitant in hepatoblastoma?
- Which tools do we need to use in order to perform global isoform profiling?
objectives:
- Perform isoform analyis in order to evaluate hepatoblastoma isoform expression profile
time_estimation: 3H
key_points:
- Isoform expression analysis
contributors:
- gallardoalba
---

# Introduction

Hepatoblastoma (HB) is the most common malignant pediatric liver tumor and one of the fastest-rising cancers in children younger than 5 years (incidence has tripled in the last 30 years) ({% cite Nagae2021 %}, {% cite Zhang2021 %}). The origin of HB is largely unknown; nearly all cases of hepatoblastoma occur in children with no previous known family history of hepatoblastoma ({% cite Tomlinson2012 %}). The molecular analysis has shown that mutations of the Wnt/β-catenin cascade (a key regulator of cell fate and proliferation during liver development and regeneration) occur in the vast majority of human HB samples, and almost exclusively affect the CTNNB1 gene (encoding β-catenin), suggesting that β-catenin pathway activation is the driver event in HB ({% cite Bell2017 %}). 

Epigenetic profile analysis indicates HB tumors are characterized by genome-wide RNA editing and DNA methylation dysregulation ({% cite CarrilloReixach2020 %}). Diverse studies have demonstrated the influence of aberrant methylation in hepatoblastoma biology by affecting genes involved in signaling and tumor suppression as well as its clinical relevance ({% cite Zhang2021 %}). Although DNA methylation was originally thought to only affect transcription, emerging evidence shows that it also regulates alternative splicing, an evolutionarily conserved mechanism that increases transcriptome and proteome diversity by allowing the generation of multiple mRNA products from a single gene ({% cite LevMaor2015 %}). So far, seven basic types of alternative splicing have been identified, including exon skipping, alternative 5′-splice site, alternative 3′-splice site, mutually exclusive exons, intron retention, alternative promoter, and alternative polyadenylation (fig. 1).

![figX:Isoform usage](../../images/differential_isoform/isoformSwitcher_splicing_patterns.png "Splicing patterns. The observed splice patterns (left colum) of two isoforms compared as indicated by the color of the splice patterns. The corresponding classification of the event (middle column) and the abreviation used (right column).")

Discovered over 40 years ago, alternative splicing formed a large part of the puzzle explaining how proteomic complexity can be achieved with a limited set of genes  ({% cite Alt1980 %}). The majority of eukaryote genes have multiple transcriptional isoforms, and recent data indicate that each transcript of protein-coding genes contain 11 exons and produce 5.4 mRNAs on average ({% cite Piovesan2016 %}). In humans,  approximately 95% of multi-exon genes show evidence of alternative splicing (AS) and approximately 60% of genes have at least one alternative transcription start site, some of which exert antagonistic functions ({% cite Carninci2006 %}, {% cite Miura2012 %}). AS regulation is essential for providing cells and tissues their specific features, and for their response to environmental changes ({% cite Wang2008 %}, {% cite Kalsotra2011 %}). Differential usage of isoforms in different conditions, often referred to as isoform switching, can have substantial biological impact, caused by the difference in the functional potential of the two isoforms ({% cite VittingSeerup2017 %}). Isoform switches are implicated in many diseases and are especially prominent in cancer; this fact has  has motivated genome-wide screens for isoform switches with predicted functional consequences ({% cite VittingSeerup2019 %}).

In this tutorial, we aim to perform a genome-wide analysis of the isoform switching phaenomena in hepatoblasmoma, which offers improved resolution over gene expression, with the objective of identify genes of clinical relevance.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Background on data

The datasets consist of ten FASTQ files, generated through the Illumina HiSeq 4000 sequencing system. The samples were obtaned by strand-specific RNA sequencing on hepatoblastoma paired samples. The procotol used for extracting the samples includes the depletion of rRNAs by subtractive hybridization, a general strategy for mRNA enrichment in RNA-seq samples. The original datasets are available in the NCBI SRA database, with the accession number [PRJNA416439](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA416439). For this tutorial, subsets from the original data were generated in order to reduce the analysis run time.

## Get data

The first step of our analysis consists of retrieving the RNA-seq datasets from Zenodo and organizing them into collections.

> <hands-on-title>Retrieve miRNA-Seq and mRNA-Seq datasets</hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from Zenodo:
>
>    - Open the file {% icon galaxy-upload %} __upload__ menu
>    - Click on __Rule-based__ tab
>    - *"Upload data as"*: `Collection(s)`
>    - Copy the following tabular data, paste it into the textbox and press <kbd>Build</kbd>
>
>      ```
>      SRR6263336	Health tissue	https://zenodo.org/record/7649293/files/sub_SRR6263336_forward.fastqsanger.gz	fastqsanger.gz	1
>      SRR6263336	Health tissue	https://zenodo.org/record/7649293/files/sub_SRR6263336_reverse.fastqsanger.gz	fastqsanger.gz	2
>      SRR6263339	Health tissue	https://zenodo.org/record/7649293/files/sub_SRR6263339_forward.fastqsanger.gz	fastqsanger.gz	1
>      SRR6263339	Health tissue	https://zenodo.org/record/7649293/files/sub_SRR6263339_reverse.fastqsanger.gz	fastqsanger.gz	2
>      SRR6263341	Health tissue	https://zenodo.org/record/7649293/files/sub_SRR6263341_forward.fastqsanger.gz	fastqsanger.gz	1
>      SRR6263341	Health tissue	https://zenodo.org/record/7649293/files/sub_SRR6263341_reverse.fastqsanger.gz	fastqsanger.gz	2
>      SRR6263343	Health tissue	https://zenodo.org/record/7649293/files/sub_SRR6263343_forward.fastqsanger.gz	fastqsanger.gz	1
>      SRR6263343	Health tissue	https://zenodo.org/record/7649293/files/sub_SRR6263343_reverse.fastqsanger.gz	fastqsanger.gz	2
>      SRR6263345	Health tissue	https://zenodo.org/record/7649293/files/sub_SRR6263345_forward.fastqsanger.gz	fastqsanger.gz	1
>      SRR6263345	Health tissue	https://zenodo.org/record/7649293/files/sub_SRR6263345_reverse.fastqsanger.gz	fastqsanger.gz	2
>      SRR6263371	Hepatoblastoma tissue	https://zenodo.org/record/7649293/files/sub_SRR6263371_forward.fastqsanger.gz	fastqsanger.gz	1
>      SRR6263371	Hepatoblastoma tissue	https://zenodo.org/record/7649293/files/sub_SRR6263371_reverse.fastqsanger.gz	fastqsanger.gz	2
>      SRR6263373	Hepatoblastoma tissue	https://zenodo.org/record/7649293/files/sub_SRR6263373_forward.fastqsanger.gz	fastqsanger.gz	1
>      SRR6263373	Hepatoblastoma tissue	https://zenodo.org/record/7649293/files/sub_SRR6263373_reverse.fastqsanger.gz	fastqsanger.gz	2
>      SRR6263375	Hepatoblastoma tissue	https://zenodo.org/record/7649293/files/sub_SRR6263375_forward.fastqsanger.gz	fastqsanger.gz	1
>      SRR6263375	Hepatoblastoma tissue	https://zenodo.org/record/7649293/files/sub_SRR6263375_reverse.fastqsanger.gz	fastqsanger.gz	2
>      SRR6263379	Hepatoblastoma tissue	https://zenodo.org/record/7649293/files/sub_SRR6263379_forward.fastqsanger.gz	fastqsanger.gz	1
>      SRR6263379	Hepatoblastoma tissue	https://zenodo.org/record/7649293/files/sub_SRR6263379_reverse.fastqsanger.gz	fastqsanger.gz	2
>      SRR6263381	Hepatoblastoma tissue	https://zenodo.org/record/7649293/files/sub_SRR6263381_forward.fastqsanger.gz	fastqsanger.gz	1
>      SRR6263381	Hepatoblastoma tissue	https://zenodo.org/record/7649293/files/sub_SRR6263381_reverse.fastqsanger.gz	fastqsanger.gz	2
>      ```
>
>    - From **Rules** menu select `Add / Modify Column Definitions`
>       - Click `Add Definition` button and select `List Identifier(s)`: column `A`
>
>         > <tip-title>Can't find <i>List Identifier</i>?</tip-title>
>         > Then you've chosen to upload as a 'dataset' and not a 'collection'. Close the upload menu, and restart the process, making sure you check *Upload data as*: **Collection(s)**
>         {: .tip}
>
>       - Click `Add Definition` button and select `Collection Name`: column `B`
>       - Click `Add Definition` button and select `URL`: column `C`
>       - Click `Add Definition` button and select `Type`: column `D`
>       - Click `Add Definition` button and select `Pair-end Indicator`: column `E`
>
>    - Click `Apply` and press <kbd>Upload</kbd>
>
> {% snippet faqs/galaxy/datasets_add_tag.md type="name" %}
{: .hands_on}

Next we will retrieve the remaining datasets.

> <hands-on-title>Retrieve the additional datasets</hands-on-title>
>
> 1. Import the files from Zenodo:
>
>    - Open the file {% icon galaxy-upload %} __upload__ menu
>    - *"Upload data as"*: `Datasets`
>    - Once again, copy the tabular data, paste it into the textbox and press <kbd>Build</kbd>
>
>      ```
>      CPAT_header.tab  https://zenodo.org/record/7656544/files/CPAT_header.tab
>      active_site.dat.gz	https://zenodo.org/record/7649293/files/active_site.dat.gz
>      gencode.v42.lncRNA_transcripts.fa.gz	https://zenodo.org/record/7649293/files/gencode.v42.lncRNA_transcripts.fa.gz
>      gencode.v42.pc_transcripts.fa.gz	https://zenodo.org/record/7649293/files/gencode.v42.pc_transcripts.fa.gz
>      gencode.v43.annotation.gtf.gz	https://zenodo.org/record/7649293/files/gencode.v43.annotation.gtf.gz
>      gencode.v43.transcripts.fa.gz	https://zenodo.org/record/7649293/files/gencode.v43.transcripts.fa.gz
>      GRCh38.p13.genome.fa.gz	https://zenodo.org/record/7649293/files/GRCh38.p13.genome.fa.gz
>      Pfam-A.hmm.dat.gz	https://zenodo.org/record/7649293/files/Pfam-A.hmm.dat.gz
>      Pfam-A.hmm.gz	https://zenodo.org/record/7649293/files/Pfam-A.hmm.gz
>      ```
>
>    - From **Rules** menu select `Add / Modify Column Definitions`
>       - Click `Add Definition` button and select `Name`: column `A`
>       - Click `Add Definition` button and select `URL`: column `B`
>    - Click `Apply` and press <kbd>Upload</kbd>
>
>
{: .hands_on}

> <details-title>Dataset subsampling</details-title>
>
> As indicated above, for this tutorial the depth of the samples was reduced in order to speed up the time needed to carry out the analysis. This was done as follows:
>
> > <hands-on-title>Dataset subsampling</hands-on-title>
> >
> > 1. {% tool [Sub-sample sequences files ](toolshed.g2.bx.psu.edu/repos/peterjc/sample_seqs/sample_seqs/0.2.5) %} with the following parameters:
> >    - {% icon param-files %} *"Multiple datasets"*: Each of the fastq files
> >    - *"Subsampling approach"*: `Take every N-th sequence (or pair e.g. every fifth sequence)`
> >    - *"N"*: `100`
> {: .hands_on}
>
> In this way, we will only take 1% of reads at a random sampling rate.
{: .details}

# Quality assessment

Once we have got the datasets, we can start with the analysis. The first step is to perform the quality assessment. Since this step is deeply covered in the tutorial [Quality control](https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/quality-control/tutorial.html), we won't describe this section in detail.


## Inicial quality evaluation

For the initial quality evaluation we will interlace the paired-read file, so instead of having the forward and the reverse reads in seperate files, each right read will be stored after the corresponding paired left read into a single file; this step is required because **FastQC** is not able to process separated paired-end reads. As probably you know, **FastQC** allows to generate a  QC report that allows to have easily an overview of basic quality control metrics. Then we will combine the individual report into a single one by making use of **MultiQC**.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [FASTQ interlacer](toolshed.g2.bx.psu.edu/repos/devteam/fastq_paired_end_interlacer/fastq_paired_end_interlacer/1.2.0.1+galaxy0) %} with the following parameters:
>    - *"Type of paired-end datasets"*: `1 paired dataset collection`
>        - {% icon param-collection %} *"Paired-end reads collection"*: `output` (Input dataset collection)
>
> 2. Repeat the previos step
>
> 3. {% tool [FastQC](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.73+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Raw read data from your current history"*: `outfile_pairs_from_coll` (output of **FASTQ interlacer** {% icon tool %})
>
> 4. Repeat the previous step
>
>
> 5. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.11+galaxy1) %} with the following parameters:
>    - In *"Results"*:
>        - {% icon param-repeat %} *"Insert Results"*
>            - *"Which tool was used generate logs?"*: `FastQC`
>                - In *"FastQC output"*:
>                    - {% icon param-repeat %} *"Insert FastQC output"*
>                        - {% icon param-file %} *"FastQC output"*: `text_file` (output of **FastQC** {% icon tool %})
>    - *"Report title"*: `Raw data QC`
>
{: .hands_on}

Let's evaluate the per base sequence quality and the adapter content.

![figX:MultiQC reports](../../images/differential_isoform/fastqc_per_base_sequence_quality.png "MultiQC aggregated reports. Mean quality scores (a). Adapter content (b).")

As we can appreciate in the figure 2.a, the per base quality of all reads seems to be very good, with values over 30 in all cases. With respect to the adapter content, the adapter content is over 20% in most samples; it means that we need to pre-process the reads before pretending to start the isoform analysis, because otherwise the results could be affected as a result of adapter contaminations.

## Read pre-processing with fastp

In order to remove the adaptors we will make use of **fastp**, which is able to detect the adapter sequence by performing a per-read overlap analysis, so we won't even need to specify the adapter sequences. 

> <hands-on-title> Task description </hands-on-title>
>
> {% tool [fastp](toolshed.g2.bx.psu.edu/repos/iuc/fastp/fastp/0.23.2+galaxy0) %} with the following parameters:
>    - *"Single-end or paired reads"*: `Paired Collection`
>        - {% icon param-collection %} *"Select paired collection(s)"*: `output` (Input dataset collection)
>        - In *"Global trimming options"*:
>            - *"Trim front for input 1"*: `10`
>    - In *"Overrepresented Sequence Analysis"*:
>        - *"Enable overrepresented analysis"*: `Yes`
>        - *"Overrepresentation sampling"*: `50`
>    - In *"Filter Options"*:
>        - In *"Quality filtering options"*:
>            - *"Qualified quality phred"*: `20`
>
{: .hands_on}

# RNA-seq mapping and isoform quantification 

The following section can be considered as the hard-core part of the training, the reason is  not because of it’s complexity (not all the details about the computational procedures will be presented, just those elements required for a basic understanding), but because it allows to characterize isoform quantification approach as genome-guided-based method. 

> <comment-title>Transcriptome-reconstruction approaches</comment-title>
>
> Transcriptome reconstruction provides critical information for isoform quantification. The different methods for estimating transcript/isoform abundance can be classified, roughly, depending on two main requierements: reference sequence and alignment. reference-guided transcriptome assembly strategy requires to aligning sequencing reads to a reference genome first, and then assembling overlapping alignments into transcripts. In contrast, *de novo* transcriptome assembly methods directly reconstructs overlapping reads into transcripts by utilizing the redundancy of sequencing reads themselves ({% cite Lu2013 %}).
>
{: .comment}

In that section makes use of three main tools: **RNA STAR**, considered a state-of-the-art mapping tool for RNA-seq data, **RSeQC**, a package that allows comprehensively evaluate different aspects of the RNA-seq data, and **Stringtie**, which uses a genome-guided transcriptome assembly approach along with concepts from de novo genome assembly to perform transcript assembly and quantification. 

## RNA-seq mapping with **RNA STAR**

**RNA STAR** is a splice-aware RNA-seq alignment tool that allows to identify canonical and non-canonical splice junctions by making use of sequential maximum mappable seed search in uncompressed suffix arrays followed by seed clustering and stitching procedure ({% cite Krianovi2017 %}). One advantage of **RNA STAR** with respect to other tools is that it includes a feature called *two-pass mode*, a framework in which splice junctions are separately discovered and quantified, allowing robustly and accurately identify  splice junction patterns for differential splicing analysis and variant discovery.

> <comment-title>Intron spanning in RNA-seq analysis</comment-title>
>
> RNA-seq mappers need to face the challenge associated with intron spanning of mature mRNA molecules, from which introns have been removed by splicing (single short read might align to two locations that are separated by 10 kbp or more). The complexity of this operation can be better understood if we take in account that for a typical human RNA-seq data set using 100-bp reads, more than 35% of the reads will span multiple exons.
>
{: .comment}

During two-pass mode splice junctions are discovered in a first alignment pass with high stringency, and are used as annotation in a second pass to permit lower stringency alignment, and therefore higher sensitivity (fig. 3). Two-pass alignment enables sequence reads to span novel splice junctions by fewer nucleotides, confering greater read depth and providing significantly more accurate quantification of novel splice junctions that one-pass alignment ({% cite Veeneman2015 %}).

![figX:Sequence quality](../../images/differential_isoform/RNASTAR_twopass_mode.png "Two-pass alignment flowchart. Center and right, stepwise progression of two-pass alignment. First, the genome is indexed with gene annotation. Next, novel splice junctions are discovered from RNA sequencing data at a relatively high stringency (12 nt minimum spanning length). Third, these discovered splice junctions, and expressed annotated splice junctions are used to re-index the genome. Finally, alignment is performed a second time, quantifying novel and annotated splice junctions using the same, relatively lower stringency (3 nt minimum spanning length), producing splice junction expression (source: Veeneman et al., 2016)")

The choice of **RNA STAR** as mapper is also determined by the sequencing technology; it has been demonstrated adequate for short-read sequencing data, but when using long-read data, such as PacBio or ONT reads, it is recommended to use **GMAP** as aligment tool ({% cite Krianovi2017 %}). So, let's perform the mapping step.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [RNA STAR](toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.8a+galaxy1) %} with the following parameters:
>    - *"Single-end or paired-end reads"*: `Paired-end (as collection)`
>        - {% icon param-file %} *"RNA-Seq FASTQ/FASTA paired reads"*: `output_paired_coll` (output of **fastp** {% icon tool %})
>    - *"Custom or built-in reference genome"*: `Use reference genome from history and create temporary index`
>        - {% icon param-file %} *"Select a reference genome"*: `output` (Input dataset)
>        - *"Build index with or without known splice junctions annotation"*: `build index with gene-model`
>            - {% icon param-file %} *"Gene model (gff3,gtf) file for splice junctions"*: `output` (Input dataset)
>    - *"Use 2-pass mapping for more sensitive novel splice junction discovery"*: `Yes, perform single-sample 2-pass mapping of all reads`
>
{: .hands_on}

Before moving to the transcriptome assembly and quantification step, we are going to use **RSeQC** in order to obtain some RNA-seq-specific quality control metrics.

## RNA-seq specific quality control metrics with **RSeQC**

RNA-seq-specific quality control metrics, such as sequencing depth, read distribution and coverage uniformity, are essential to ensure that the RNA-seq data are adequate for transcriptome reconstruction and alternative splicing analysis. For example, the use of RNA-seq with unsaturated sequencing depth gives imprecise estimations  and fails to detect low abundance splice junctions, limiting the precision of many analyses ({% cite Wang2012 %}).

In this section we will make use of of the **RSeQC** toolkit in order to generate the RNA-seq-specific quality control metrics. But before starting, we need to convert the annotation GTF file into BED12 format, which will be required in subsequent steps.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Convert GTF to BED12](toolshed.g2.bx.psu.edu/repos/iuc/gtftobed12/gtftobed12/357) %} with the following parameters:
>    - {% icon param-file %} *"GTF File to convert"*: `output` (Input dataset)
>    - *"Advanced options"*: `Set advanced options`
>        - *"Ignore groups without exons"*: `Yes`
>
{: .hands_on}

We are going to use the following RSeQC modules:

- **Infer Experiment**: inference of RNA-seq configuration
- **Gene Body Coverage**: compute read coverage over gene bodies
- **Junction Saturation**: check junction saturation
- **Junction Annotation**: compares detected splice junctions to a reference gene model
- **Read Distribution**: calculates how mapped reads are distributed over genome features

Once all required outputs have been generated, we will integrate them by using **MultiQC** in order to interpretate the results.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Infer Experiment](toolshed.g2.bx.psu.edu/repos/nilesh/rseqc/rseqc_infer_experiment/5.0.1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input BAM file"*: `mapped_reads` (output of **RNA STAR** {% icon tool %})
>    - {% icon param-file %} *"Reference gene model"*: `bed_file` (output of **Convert GTF to BED12** {% icon tool %})
>
> 2. {% tool [Gene body coverage (BAM)](toolshed.g2.bx.psu.edu/repos/nilesh/rseqc/rseqc_geneBody_coverage/5.0.1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input BAM file"*: `mapped_reads` (output of **RNA STAR** {% icon tool %})
>    - {% icon param-file %} *"Reference gene model"*: `bed_file` (output of **Convert GTF to BED12** {% icon tool %})
>
> 3. {% tool [Junction Saturation](toolshed.g2.bx.psu.edu/repos/nilesh/rseqc/rseqc_junction_saturation/5.0.1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input BAM/SAM file"*: `mapped_reads` (output of **RNA STAR** {% icon tool %})
>    - {% icon param-file %} *"Reference gene model"*: `bed_file` (output of **Convert GTF to BED12** {% icon tool %})
>    - *"Sampling bounds and frequency"*: `Default sampling bounds and frequency`
>    - *"Output R-Script"*: `Yes`
>
> 4. {% tool [Junction Annotation](toolshed.g2.bx.psu.edu/repos/nilesh/rseqc/rseqc_junction_annotation/5.0.1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input BAM/SAM file"*: `mapped_reads` (output of **RNA STAR** {% icon tool %})
>    - {% icon param-file %} *"Reference gene model"*: `bed_file` (output of **Convert GTF to BED12** {% icon tool %})
>
> 5. {% tool [Read Distribution](toolshed.g2.bx.psu.edu/repos/nilesh/rseqc/rseqc_read_distribution/5.0.1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input BAM/SAM file"*: `mapped_reads` (output of **RNA STAR** {% icon tool %})
>    - {% icon param-file %} *"Reference gene model"*: `bed_file` (output of **Convert GTF to BED12** {% icon tool %})
>
{: .hands_on}

Now, we will integrate the outputs into **MultiQC**.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.11+galaxy1) %} with the following parameters:
>    - In *"Results"*:
>        - {% icon param-repeat %} *"Insert Results"*
>            - *"Which tool was used generate logs?"*: `STAR`
>                - In *"STAR output"*:
>                    - {% icon param-repeat %} *"Insert STAR output"*
>                        - *"Type of STAR output?"*: `Log`
>                            - {% icon param-file %} *"STAR log output"*: `output_log` (output of **RNA STAR** {% icon tool %})
>                    - {% icon param-repeat %} *"Insert STAR output"*
>                        - *"Type of STAR output?"*: `Log`
>        - {% icon param-repeat %} *"Insert Results"*
>            - *"Which tool was used generate logs?"*: `RSeQC`
>                - In *"RSeQC output"*:
>                    - {% icon param-repeat %} *"Insert RSeQC output"*
>                        - *"Type of RSeQC output?"*: `Infer experiment`
>                            - {% icon param-file %} *"RSeQC infer experiment: configuration output"*: `output` (output of **Infer Experiment** {% icon tool %})
>                    - {% icon param-repeat %} *"Insert RSeQC output"*
>                        - *"Type of RSeQC output?"*: `Infer experiment`
>                    - {% icon param-repeat %} *"Insert RSeQC output"*
>                        - *"Type of RSeQC output?"*: `Read distribution`
>                    - {% icon param-repeat %} *"Insert RSeQC output"*
>                        - *"Type of RSeQC output?"*: `Read distribution`
>                    - {% icon param-repeat %} *"Insert RSeQC output"*
>                        - *"Type of RSeQC output?"*: `Junction saturation`
>                    - {% icon param-repeat %} *"Insert RSeQC output"*
>                        - *"Type of RSeQC output?"*: `Junction saturation`
>                    - {% icon param-repeat %} *"Insert RSeQC output"*
>                        - *"Type of RSeQC output?"*: `Junction annotation`
>                    - {% icon param-repeat %} *"Insert RSeQC output"*
>                        - *"Type of RSeQC output?"*: `Junction annotation`
>                    - {% icon param-repeat %} *"Insert RSeQC output"*
>                        - *"Type of RSeQC output?"*: `Gene body coverage`
>                    - {% icon param-repeat %} *"Insert RSeQC output"*
>                        - *"Type of RSeQC output?"*: `Gene body coverage`
>
{: .hands_on}

Now we can have a look at the **RSeQC** results; we can start by evaluating the plot corresponding to the **Infer Experiment** module, which allows to speculate the experimental design (whether sequencing is strand-specific, and if so, how reads are stranded) by sampling a subset of reads from the BAM file and comparing their genome coordinates and strands with those of the reference gene model ({% cite Wang2012 %}).  

![figX:RSeQC infer experiment](../../images/differential_isoform/rseqc_infer_experiment.png "RSeQC Infer Experiment plot. It counts the percentage of reads and read pairs that match the strandedness of overlapping transcripts. It can be used to infer whether RNA-seq library preps are stranded (sense or antisense).")

As can be appreciated in the image, the proportion of reads assigned as *sense* is similar to the ones assigned as *antisense*, which indicates that in that case our RNA-seq data is non-strand specific. 

Now, let's evaluate the results generated by the **Gene Body Coverage** module. It scales all transcripts to 100 nt and calculates the number of reads covering each nucleotide position. The plot generated from this information illustrates the coverage profile along the gene body, defined as the entire gene from the transcription start site to the end of the transcript (fig. 5).

![figX:RSeQC gene body coverage](../../images/differential_isoform/rseqc_gene_body_coverage_plot.png "RSeQC gene body coverage plot. It calculates read coverage over gene bodies. This is used to check if reads coverage is uniform and if there is any 5' or 3' bias.")

The gene body coverage pattern is highly influenced by the RNA-seq protocol, and it is useful for identifying artifacts such as 3' skew in libraries. For example, a skew towards increased 3' coverage can happen in degraded samples prepared with poly-A selection. According the figure 5, there're not bias in our reads as a result of sequencing technical problems.

Other important metric for alternative splicing analysis is the one provided by the **Junction Saturation** module, which allows to determine if the current sequencing depth is sufficient to perform alternative splicing analyses; splice junctions are detected for each re-sampled subset of reads, and the number of detected splice junctions increases as the resample percentage increases before finally reaching a fixed value. Using an unsaturated sequencing depth would miss many rare splice junctions ({% cite Wang2012 %}). In that case, we will focus only in the **known junctions** plot (figure 6).

![figX:RSeQC junction saturation](../../images/differential_isoform/rseqc_junction_saturation.png "RSeQC junction saturation.")

As we can appreciate in the plot, the known junctions tend to stabilize around 160.000, which indicates that the read sequencing depth is good enough for performing the alternative splicing analysis.

<!-- Needs to be done -->

It separates all detected splice junctions into ‘known’, ‘complete novel’ and ‘partial novel’ by comparing them with the reference gene model {% cite Wang2012 %}.

It compare detected splice junctions to reference gene model. splicing annotation is performed in two levels: splice event level and splice junction level.

- Splice read: An RNA read, especially long read, can be spliced more than once, therefore, 100 spliced reads can produce >= 100 splicing events.
- Splice junction: multiple splicing events spanning the same intron can be consolidated into one splicing junction.

Detected junctions were divided to 3 exclusive categories:

- Annotated (known): The junction is part of the gene model. Both splice sites, 5’ splice site (5’SS) and 3’splice site (3’SS) are annotated by reference gene model.
- Complete_novel: Both 5’SS and 3’SS are novel.
- Partial_novel: One of the splice site (5’SS or 3’SS) is novel, and the other splice site is annotated.

### Genome features analysis with **RSeQC Read Distribution**

<!-- Needs to be edited -->

Provided a BAM/SAM file and reference gene model, this module will calculate how mapped reads were distributed over genome feature (like CDS exon, 5’UTR exon, 3’ UTR exon, Intron, Intergenic regions). When genome features are overlapped (e.g. a region could be annotated as both exon and intron by two different transcripts) , they are prioritize as: CDS exons > UTR exons > Introns > Intergenic regions, for example, if a read was mapped to both CDS exon and intron, it will be assigned to CDS exons.




![figX:STAR alignment](../../images/differential_isoform/star_alignment.png "RNA star alignment")





![figX:RSeQC junction annotation](../../images/differential_isoform/rseqc_junction_annotation_junctions.png "RSeQC junction annotation")



![figX:RSeQC read distribution](../../images/differential_isoform/rseqc_read_distribution_plot.png "RSeQC read distribution")


## Transcriptome assembly and quantification with **StringTie**

<!-- Needs to be edited -->

Using a mapping of reads to the reference genome, genome-guided transcript assemblers cluster the reads and build graph models representing all possible isoforms for each gene. One such model is a splice graph, in which nodes represent exons or parts of exons, and paths through the graph represent possible splice variants.

StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. It uses a novel network flow algorithm as well as an optional de novo assembly step to assemble and quantitate full-length transcripts representing multiple splice variants for each gene locus. Its input can include not only alignments of short reads that can also be used by other transcript assemblers, but also alignments of longer sequences that have been assembled from those reads {% cite Pertea2015 %}.

StringTie assembles transcripts and estimates their expression levels simultaneously. StringTie first groups the reads into clusters, then creates a splice graph for each cluster from which it identifies transcripts, and then for each transcript it creates a separate flow network to estimate its expression level using a maximum flow algorithm {% cite Pertea2015 %}.

![figX:Stringtie algorithm](../../images/differential_isoform/stringtie_algorithm.png "Transcript assembly pipeline for StringTie. It begin with a set of RNA-seq reads that have been mapped to the genome. StringTie iteratively extracts the heaviest path from a splice graph, constructs a flow network, computes maximum flow to estimate abundance, and then updates the splice graph by removing reads that were assigned by the flow algorithm. This process repeats until all reads have been assigned.")

The main reason underlying the greater accuracy of StringTie most likely derives from its optimization criteria. By balancing the coverage (or flow) of each transcript across each assembly, it incorporates depth of coverage constraints into the assembly algorithm itself. When assembling a whole genome, coverage is a crucial parameter that must be used to constrain the algorithm; otherwise an assembler may incorrectly collapse repetitive sequences. Similarly, when assembling a transcript, each exon within an isoform should have similar coverage, and ignoring this parameter may produce sets of transcripts that are parsimonious but wrong {% cite Pertea2015 %}.

StringTie is a transcript assembler that uses the optimization technique of maximum flow in a specially constructed flow network to determine gene expression levels, and does so while simultaneously assembling each isoform of a gene. And unlike other transcript assemblers, it incorporates alignment to both a genome and a de novo assembly of reads {% cite Pertea2015 %}.

StringTie takes as input a SAM, BAM or CRAM file sorted by coordinate (genomic location). Any SAM record with a spliced alignment (i.e. having a read alignment across at least one junction) should have the XS tag (or the ts tag, see below) which indicates the transcription strand, the genomic strand from which the RNA that produced the read originated . 

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [StringTie](toolshed.g2.bx.psu.edu/repos/iuc/stringtie/stringtie/2.2.1+galaxy1) %} with the following parameters:
>    - *"Input options"*: `Short reads`
>        - {% icon param-file %} *"Input short mapped reads"*: `mapped_reads` (output of **RNA STAR** {% icon tool %})
>    - *"Use a reference file to guide assembly?"*: `Use reference GTF/GFF3`
>        - *"Reference file"*: `Use a file from history`
>            - {% icon param-file %} *"GTF/GFF3 dataset to guide assembly"*: `output` (Input dataset)
>        - *"Use Reference transcripts only?"*: `Yes`
>        - *"Output files for differential expression?"*: `Ballgown`
>
{: .hands_on}

Stringtie generates different table files; in our case, we are only interested in the collection called **transcript-level expression measurements** (t_tab.ctab files). Each file includes one row per transcript, with the following columns:

- t_id: numeric transcript id
- chr, strand, start, end: genomic location of the transcript
- t_name: generated transcript id
- num_exons: number of exons comprising the transcript
- length: transcript length, including both exons and introns
- gene_id: gene the transcript belongs to
- gene_name: HUGO gene name for the transcript, if known
- cov: per-base coverage for the transcript (available for each sample)
- FPKM: Estimated FPKM for the transcript (available for each sample)

# Genome-wide isoform switch analysis with **IsoformSwitchAnalyzeR**

<!-- Needs to be edited! -->

IsoformSwitchAnalyzeR enables analysis of changes in genome-wide patterns of alternative splicing and isoform switch consequences.

A genome-wide analysis is both useful for getting an overview of the extent of isoform switching as well as discovering general patterns. IsoformSwitchAnalyzeR supports this by providing four different summaries/analyses for both the analysis of alternative splicing and isoform switches with predicted consequences. All functions provide a visual overview as well as a data.frame with the summary statistics. The four analysis types supported are:

- Global summary statistics, implemented in the extractConsequenceSummary() and extractSplicingSummary() functions, which summarizes the number of switches with predicted consequences and the number of splicing events occurring in the different comparisons, respectively.
- Analysis of splicing/consequence enrichment, implemented in the extractConsequenceEnrichment() and extractSplicingEnrichment() functions, which analyzes whether a particular consequence/splice type occurs more frequently than the opposite event (e.g. domain loss vs domain gain) in a giving comparison (e.g. WT->KO1).
- Comparison of enrichment, implemented in extractConsequenceEnrichmentComparison() and extractSplicingEnrichmentComparison() functions, which compares the enrichment of a particular consequence/splice type between comparisons (e.g. compares the changes in WT->KO1 vs WT->KO2)
- Analysis of genome-wide changes in isoform usage, implemented in extractConsequenceGenomeWide() and extractSplicingGenomeWide() functions, which analyses the genome-wide changes in isoform usage for all isoforms with particular opposite pattern events. This type of analysis is particular interesting if the expected difference between conditions is large, since such effects could result in genome-wide changes. The analysis works by simultaneously analyzing all isoforms with a specific feature (e.g. intron retention) for changes in isoform usage.

## Import data into **IsoformSwitchAnalyzeR**


> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [IsoformSwitchAnalyzeR](toolshed.g2.bx.psu.edu/repos/iuc/isoformswitchanalyzer/isoformswitchanalyzer/1.20.0+galaxy0) %} with the following parameters:
>    - *"Tool function mode"*: `Import data`
>        - In *"1: Factor level"*:
>            - *"Specify a factor level, typical values could be 'tumor' or 'treated'"*: `Cancer`
>            - {% icon param-file %} *"Transcript-level expression measurements"*: `transcript_expression` (output of **StringTie** {% icon tool %})
>        - In *"2: Factor level"*:
>            - *"Specify a factor level, typical values could be 'tumor' or 'treated'"*: `Health`
>            - {% icon param-file %} *"Transcript-level expression measurements"*: `transcript_expression` (output of **StringTie** {% icon tool %})
>        - *"Quantification data source"*: `StringTie`
>            - *"Average read length"*: `140`
>        - {% icon param-file %} *"Genome annotation (GTF)"*: `output` (Input dataset)
>        - {% icon param-file %} *"Transcriptome"*: `output` (Input dataset)
>
{: .hands_on}


## Filtering and isoform switching identification

<!-- Needs to be edited! -->

Once you have a switchAnalyzeRlist, there is a good chance that it contains a lot of genes/isoforms that are irrelevant for an analysis of isoform switches. Examples of such could be single isoform genes or non-expressed isoforms. These extra genes/isoforms will make the downstream analysis take (much) longer than necessary. Therefore we have implemented a pre-filtering step to remove these features before continuing with the analysis. Importantly, filtering can enhance the reliability of the downstream analysis as described in detail below.

By using preFilter() it is possible to remove genes and isoforms from all aspects of the switchAnalyzeRlist by filtering on:

- Multi-isoform genes
- Gene expression
- Isoform expression
- Isoform Fraction (isoform usage)
- Unwanted isoform classes
- Unwanted gene biotypes
- Genes without differential isoform usage
- Removal of single isoform genes is the default setting in preFilter() since these genes, per definition, cannot have changes in isoform usage.

Filtering on isoform expression allows removal of non-used isoforms that only appear in the switchAnalyzeRlist because they were in the isoform/gene annotation used. Furthermore, the expression filtering allows removal of lowly expressed isoforms where the expression levels might be untrustworthy. 


<!-- Needs to be edited! -->

Two major challenges in testing differential isoform usage have been controlling false discovery rates (FDR) and applying effect size cutoffs in experimental setups with confounding effects. Recent studies such as Love at al highlights DEXSeq (developed by Anders et al., see What To Cite — please remember to cite it) as being a good solution as it controls FDR quite well. We have therefore implemented a DEXSeq based test as the default in IsoformSwitchAnalyzeR. This test furthermore utilizes limma to produce effect sizes corrected for confounding effects.

An important argument in isoformSwitchTestDEXSeq is the ‘reduceToSwitchingGenes’. When TRUE this argument will cause the function to reduce/subset the switchAnalyzeRlist to the genes which each contains at least one differential used isoform, as indicated by the alpha and dIFcutoff cutoffs. This option ensures the rest of the workflow runs significantly faster since isoforms from genes without isoform switching are not analyzed.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [IsoformSwitchAnalyzeR](toolshed.g2.bx.psu.edu/repos/iuc/isoformswitchanalyzer/isoformswitchanalyzer/1.20.0+galaxy0) %} with the following parameters:
>    - *"Tool function mode"*: `Analysis part one: Extract isoform switches and their sequences`
>        - {% icon param-file %} *"IsoformSwitchAnalyzeR R object"*: `switchList` (output of **IsoformSwitchAnalyzeR** {% icon tool %})
>
{: .hands_on}

### Importing external sequence analysis and genome-wide analysis

### Protein domain identification with **PfamScan**

<!-- Needs to be edited! -->

PfamScan is used to search a FASTA sequence against a library of Pfam HMM.

Pfam is a database of protein families and domains that is widely used to analyse novel genomes, metagenomes and to guide experimental work on particular proteins and systems (1,2). Each Pfam family has a seed alignment that contains a representative set of sequences for the entry. A profile hidden Markov model (HMM) is automatically built from the seed alignment and searched against a sequence database called pfamseq using the HMMER software (http://hmmer.org/). All sequence regions that satisfy a family-specific curated threshold, also known as the gathering threshold, are aligned to the profile HMM to create the full alignment {% cite Mistry2020 %}.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [PfamScan](toolshed.g2.bx.psu.edu/repos/bgruening/pfamscan/pfamscan/1.6+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Protein sequences FASTA file"*: `isoformAA` (output of **IsoformSwitchAnalyzeR** {% icon tool %})
>    - {% icon param-file %} *"Pfam-A HMM library"*: `output` (Input dataset)
>    - {% icon param-file %} *"Pfam-A HMM Stockholm file"*: `output` (Input dataset)
>    - *"Predict active site residues"*: `Enabled`
>        - {% icon param-file %} *"Active sites file"*: `output` (Input dataset)
>
{: .hands_on}

### RNA coding probablity prediction with **CPAT**

<!-- Needs to be edited! -->

CPAT is a bioinformatics tool to predict RNA’s coding probability based on the RNA sequence characteristics. To achieve this goal, CPAT calculates scores of these 4 linguistic features from a set of known protein-coding genes and another set of non-coding genes.

- ORF size
- ORF coverage
- Fickett TESTCODE
- Hexamer usage bias

CPAT will then builds a logistic regression model using these 4 features as predictor variables and the “protein-coding status” as the response variable. After evaluating the performance and determining the probability cutoff, the model can be used to predict new RNA sequences.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [CPAT](toolshed.g2.bx.psu.edu/repos/bgruening/cpat/cpat/3.0.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Query nucletide sequences"*: `isoformNT` (output of **IsoformSwitchAnalyzeR** {% icon tool %})
>    - {% icon param-file %} *"Reference genome"*: `output` (Input dataset)
>    - {% icon param-file %} *"Coding sequences file"*: `output` (Input dataset)
>    - {% icon param-file %} *"Non coding sequeces file"*: `output` (Input dataset)
>
{: .hands_on}

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Remove beginning](https://usegalaxy.eu/root?tool_id=Remove beginning1) %} with the following parameters:
>    - *"Remove first"*: `1`
>    - {% icon param-file %} *"From"*: `CPAT on data...`
>
> 2. {% tool [Text reformatting](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_awk_tool/1.1.2) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `Remove beginning on data...`
>    - *"AWK Program"*: `{print i++"\t"$1"\t"$3"\t"$8"\t"$9"\t"$10"\t"$11"\t""-"}`
>
> 3. {% tool [Concatenate datasets](https://usegalaxy.eu/root?tool_id=cat1) %} with the following parameters:
>    - {% icon param-file %} *"Concatenate Dataset"*: `CPAT_header.tab`
>    - In *"Dataset"*:
>       - Click in "*Insert Dataset*"
>    - In *"1: Dataset"*:
>       - {% icon param-file %} *"Select"*: `Text reformatting on...`
>
{: .hands_on}

If an isoform has a significant change in its contribution to gene expression, there must per definition be reciprocal changes in one (or more) isoforms in the opposite direction, compensating for the change in the first isoform. We utilize this by extracting the isoforms that are significantly differentially used and compare them to the isoforms that are compensating. Using all the information gathered through the workflow described above, the annotation of the isoform(s) used more (positive dIF) can be compared to the isoform(s) used less (negative dIF) and by systematically identify differences annotation we can identify potential function consequences of the isoform switch.

Specifically, IsoformSwitchAnalyzeR contains a function analyzeSwitchConsequences() which extracts the isoforms with significant changes in their isoform usage (defined by the alpha and dIFcutoff parameters, see Identifying Isoform Switches for details) and the isoform, with a large opposite change in isoform usage (also controlled via the dIFcutoff parameters) that compensate for the changes. Note that if an isoform-level test was not used, the gene is require to be significant (defined by the alpha parameter); but, isoforms are then selected purely based on their changes in dIF values.

These isoforms are then divided into the isoforms that increase their contribution to gene expression (positive dIF values larger than dIFcutoff) and the isoforms that decrease their contribution (negative dIF values smaller than -dIFcutoff). The isoforms with increased contribution are then (in a pairwise manner) compared to the isoform with decreasing contribution. In each of these comparisons the isoforms compared are analyzed for differences in their annotation (controlled by the consequencesToAnalyze parameter). Currently 22 different features of the isoforms can be compared, which include features such as intron retention, coding potential, NMD status, protein domains and the sequence similarity of the amino acid sequence of the annotated ORFs. 

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [IsoformSwitchAnalyzeR](toolshed.g2.bx.psu.edu/repos/iuc/isoformswitchanalyzer/isoformswitchanalyzer/1.20.0+galaxy0) %} with the following parameters:
>    - *"Tool function mode"*: `Analysis part two: Plot all isoform switches and their annotation`
>        - {% icon param-file %} *"IsoformSwitchAnalyzeR R object"*: `switchList` (output of **IsoformSwitchAnalyzeR** {% icon tool %})
>        - *"Analysis mode"*: `Full analysis`
>        - *"Include prediction of coding potential information"*: `CPAT`
>            - {% icon param-file %} *"CPAT result file"*: `orf_seqs_prob_best` (output of **CPAT** {% icon tool %})
>        - *"Include Pfam information"*: `Enabled`
>            - {% icon param-file %} *"Include Pfam results (sequence analysis of protein domains)"*: `output` (output of **PfamScan** {% icon tool %})
>        - *"Include SignalP results"*: `Disabled`
>        - *"Include prediction of intrinsically disordered Regions (IDR) information"*: `Disabled`
>
{: .hands_on}

![figX:Consequences isoform](../../images/differential_isoform/isoformSwitchAnalyzer_consequences_isoform.png "Consequences isoform")

![figX:Consequences feature](../../images/differential_isoform/isoformSwitchAnalyzer_consequences_features.png "Consequences features")

![figX:Splicing event](../../images/differential_isoform/isoformSwitchAnalyzer_splicing_event.png "Splicing event")

![figX:Summary](../../images/differential_isoform/isoformSwitchAnalyzer_summary.png "Summary")

![figX:Isoform usage](../../images/differential_isoform/isoformSwitchAnalyzer_isoform_usage.png "Isoform usage")

### Analysis of individual isoform switching

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [IsoformSwitchAnalyzeR](toolshed.g2.bx.psu.edu/repos/iuc/isoformswitchanalyzer/isoformswitchanalyzer/1.20.0+galaxy0) %} with the following parameters:
>    - *"Tool function mode"*: `Analysis part two: Plot all isoform switches and their annotation`
>        - {% icon param-file %} *"IsoformSwitchAnalyzeR R object"*: `switchList` (output of **IsoformSwitchAnalyzeR** {% icon tool %})
>        - *"Analysis mode"*: `Analyze specific gene`
>        - *"Gene name"*: `Example name`
>        - *"Include prediction of coding potential information"*: `CPAT`
>            - {% icon param-file %} *"CPAT result file"*: `orf_seqs_prob_best` (output of **CPAT** {% icon tool %})
>        - *"Include Pfam information"*: `Enabled`
>            - {% icon param-file %} *"Include Pfam results (sequence analysis of protein domains)"*: `output` (output of **PfamScan** {% icon tool %})
>        - *"Include SignalP results"*: `Disabled`
>        - *"Include prediction of intrinsically disordered Regions (IDR) information"*: `Disabled`
>
{: .hands_on}


![figX:Isoform usage](../../images/differential_isoform/isoformSwitchAnalyzer_gene.png "PHLPP2 isoform expression profile")


# Conclusion

Despite the large amount of RNA-seq data and computational methods available, isoform-based expression analysis is rare. This means that the potential of existing RNA-seq data is untapped, and as a consequence, our general understanding of differential isoform usage is poor. The few efforts at analyzing individual isoform switches have typically dealt with isoforms by describing their frequent occurrence rather than trying to systematically predict their consequence. Overall, this is unsatisfying, as isoform usage is important in disease and especially cancer, where many individual isoform switches have been described {% cite VittingSeerup2017 %}.

Here we present methods for the statistical identification and analysis of isoform switches with predicted functional consequences {% cite VittingSeerup2017 %}.