---
layout: tutorial_hands_on

title: Differential isoform expression
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

<!-- Needs to be edited! -->

Hepatoblastoma (HB) is the most common malignant pediatric liver tumor and one of the fastest-rising cancers in children younger than 5 years {% cite Nagae2021 %}.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Experimental design

Text about experimental design.

# Background on data

Text about data background.

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
>      SRR11611349	Control miRNA	https://zenodo.org/record/4710649/files/SRR11611349_MIRNASEQ_CTL.fastqsanger.gz	fastqsanger.gz
>      SRR11611350	Control miRNA	https://zenodo.org/record/4710649/files/SRR11611350_MIRNASEQ_CTL.fastqsanger.gz	fastqsanger.gz
>      SRR11611351	Control miRNA	https://zenodo.org/record/4710649/files/SRR11611351.MIRNASEQ_CTLfastqsanger.gz	fastqsanger.gz
>      SRR11611352	BR treated miRNA	https://zenodo.org/record/4710649/files/SRR11611352_MIRNASEQ_BL.fastqsanger.gz	fastqsanger.gz
>      SRR11611353	BR treated miRNA	https://zenodo.org/record/4710649/files/SRR11611353_MIRNASEQ_BL.fastqsanger.gz	fastqsanger.gz
>      SRR11611354	BR treated miRNA	https://zenodo.org/record/4710649/files/SRR11611354_MIRNASEQ_BL.fastqsanger.gz	fastqsanger.gz
>      SRR1019436	Control mRNA	https://zenodo.org/record/4710649/files/SRR1019436_RNASEQ_CTL.fastqsanger.gz	fastqsanger.gz
>      SRR1019437	Control mRNA	https://zenodo.org/record/4710649/files/SRR1019437_RNASEQ_CTL.fastqsanger.gz	fastqsanger.gz
>      SRR1019438	BR treated mRNA	https://zenodo.org/record/4710649/files/SRR1019438_RNASEQ_BL.fastqsanger.gz	fastqsanger.gz
>      SRR1019439	BR treated mRNA	https://zenodo.org/record/4710649/files/SRR1019439_RNASEQ_BL.fastqsanger.gz	fastqsanger.gz
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
>      annotation_AtRTD2.gtf	https://zenodo.org/record/4710649/files/annotation_AtRTD2_19April2016.gtf.gz
>      transcriptome.fasta	https://zenodo.org/record/4710649/files/transcriptome_AtRTD2_12April2016.fasta.gz
>      star_miRNA_seq.fasta	https://zenodo.org/record/4710649/files/star_miRNA_seq.fasta
>      mature_miRNA_AT.fasta	https://zenodo.org/record/4710649/files/mature_miRNA_AT.fasta
>      miRNA_stem-loop_seq.fasta	https://zenodo.org/record/4710649/files/miRNA_stem-loop_seq.fasta
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

# Initial quality control

Quality assessment text.

## Sub-step with **FASTQ interlacer**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [FASTQ interlacer](toolshed.g2.bx.psu.edu/repos/devteam/fastq_paired_end_interlacer/fastq_paired_end_interlacer/1.2.0.1+galaxy0) %} with the following parameters:
>    - *"Type of paired-end datasets"*: `1 paired dataset collection`
>        - {% icon param-collection %} *"Paired-end reads collection"*: `output` (Input dataset collection)
>
{: .hands_on}

Repeat with the other collection.

## Sub-step with **FastQC**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [FastQC](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.73+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Raw read data from your current history"*: `outfile_pairs_from_coll` (output of **FASTQ interlacer** {% icon tool %})
>
{: .hands_on}

Repeat with the other collection.

## Sub-step with **MultiQC**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.11+galaxy1) %} with the following parameters:
>    - In *"Results"*:
>        - {% icon param-repeat %} *"Insert Results"*
>            - *"Which tool was used generate logs?"*: `FastQC`
>                - In *"FastQC output"*:
>                    - {% icon param-repeat %} *"Insert FastQC output"*
>                        - {% icon param-file %} *"FastQC output"*: `text_file` (output of **FastQC** {% icon tool %})
>    - *"Report title"*: `Raw data QC`
>
{: .hands_on}

# Quality assessment

## Sub-step with **fastp**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [fastp](toolshed.g2.bx.psu.edu/repos/iuc/fastp/fastp/0.23.2+galaxy0) %} with the following parameters:
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
>    - In *"Output Options"*:
>        - *"Output HTML report"*: `Yes`
>        - *"Output JSON report"*: `Yes`
>
{: .hands_on}

Repeat with the other collection

## Sub-step with **MultiQC**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.11+galaxy1) %} with the following parameters:
>    - In *"Results"*:
>        - {% icon param-repeat %} *"Insert Results"*
>            - *"Which tool was used generate logs?"*: `fastp`
>                - {% icon param-file %} *"Output of fastp"*: `report_json` (output of **fastp** {% icon tool %})
>        - {% icon param-repeat %} *"Insert Results"*
>            - *"Which tool was used generate logs?"*: `fastp`
>                - {% icon param-file %} *"Output of fastp"*: `report_json` (output of **fastp** {% icon tool %})
>
{: .hands_on}


# RNA-seq assembly

## Sub-step with **RNA STAR**

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
>    - *"Per gene/transcript output"*: `No per gene or transcript output`
>    - In *"BAM output format specification"*:
>        - *"Read alignment tags to include in the BAM output"*: ``
>    - In *"Output filter criteria"*:
>        - *"Would you like to set additional output filters?"*: `No`
>    - In *"Algorithmic settings"*:
>        - *"Configure seed, alignment and limits options"*: `Use Defaults`
>
{: .hands_on}

## Post-alignment RNA-seq analysis


## Sub-step with **Convert GTF to BED12**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Convert GTF to BED12](toolshed.g2.bx.psu.edu/repos/iuc/gtftobed12/gtftobed12/357) %} with the following parameters:
>    - {% icon param-file %} *"GTF File to convert"*: `output` (Input dataset)
>    - *"Advanced options"*: `Set advanced options`
>        - *"Ignore groups without exons"*: `Yes`
>
{: .hands_on}

## Sub-step with **Junction Annotation**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Junction Annotation](toolshed.g2.bx.psu.edu/repos/nilesh/rseqc/rseqc_junction_annotation/5.0.1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input BAM/SAM file"*: `mapped_reads` (output of **RNA STAR** {% icon tool %})
>    - {% icon param-file %} *"Reference gene model"*: `bed_file` (output of **Convert GTF to BED12** {% icon tool %})
>
{: .hands_on}

## Sub-step with **Junction Saturation**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Junction Saturation](toolshed.g2.bx.psu.edu/repos/nilesh/rseqc/rseqc_junction_saturation/5.0.1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input BAM/SAM file"*: `mapped_reads` (output of **RNA STAR** {% icon tool %})
>    - {% icon param-file %} *"Reference gene model"*: `bed_file` (output of **Convert GTF to BED12** {% icon tool %})
>    - *"Sampling bounds and frequency"*: `Default sampling bounds and frequency`
>    - *"Output R-Script"*: `Yes`
>
{: .hands_on}

## Sub-step with **Read Distribution**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Read Distribution](toolshed.g2.bx.psu.edu/repos/nilesh/rseqc/rseqc_read_distribution/5.0.1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input BAM/SAM file"*: `mapped_reads` (output of **RNA STAR** {% icon tool %})
>    - {% icon param-file %} *"Reference gene model"*: `bed_file` (output of **Convert GTF to BED12** {% icon tool %})
>
{: .hands_on}

## Sub-step with **Infer Experiment**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Infer Experiment](toolshed.g2.bx.psu.edu/repos/nilesh/rseqc/rseqc_infer_experiment/5.0.1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input BAM file"*: `mapped_reads` (output of **RNA STAR** {% icon tool %})
>    - {% icon param-file %} *"Reference gene model"*: `bed_file` (output of **Convert GTF to BED12** {% icon tool %})
>
{: .hands_on}

## Sub-step with **MultiQC**

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
>
{: .hands_on}

## Transcriptome assembly

## Sub-step with **StringTie**

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

## Isoform analysis

## Sub-step with **IsoformSwitchAnalyzeR**

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

## Sub-step with **IsoformSwitchAnalyzeR**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [IsoformSwitchAnalyzeR](toolshed.g2.bx.psu.edu/repos/iuc/isoformswitchanalyzer/isoformswitchanalyzer/1.20.0+galaxy0) %} with the following parameters:
>    - *"Tool function mode"*: `Analysis part one: Extract isoform switches and their sequences`
>        - {% icon param-file %} *"IsoformSwitchAnalyzeR R object"*: `switchList` (output of **IsoformSwitchAnalyzeR** {% icon tool %})
>
{: .hands_on}

## Sub-step with **PfamScan**

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

## Sub-step with **CPAT**

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
>    - {% icon param-file %} *"Concatenate Dataset"*: `Text reformatting on...`
>    - In *"Dataset"*:
>       - Click in "*Insert Dataset*"
>    - In *"1: Dataset"*:
>       - {% icon param-file %} *"Select"*: `CPAT_header.tab`
>
{: .hands_on}

## Sub-step with **IsoformSwitchAnalyzeR**

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

## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.