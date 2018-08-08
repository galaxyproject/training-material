---
layout: tutorial_hands_on
topic_name: transcriptomics
tutorial_name: clipseq
---

# Introduction
{:.no_toc}

The eCLIP data provided here is a subset of the eCLIP data of RBFOX2 from a study published by *Nostrand et al.* (2016, http://dx.doi.org/10.1038/nmeth.3810). The dataset contains the first biological replicate of RBFOX2 CLIP-seq and the input control experiment (fastq files). The data was changed and downsampled to reduce data processing time, thus the datasets does not correspond to the original data pulled from *Nostrand et al.* (2016, http://dx.doi.org/10.1038/nmeth.3810). Also included is a text file (.txt) encompassing the chromosome sizes of hg19 obtained from UCSC (http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes) and a genome annotation (.gtf) file taken from Ensembl (http://ftp.ensemblorg.ebi.ac.uk/pub/release-74/gtf/homo_sapiens/).

**Table 1**: Metadata for CLIP-seq experiments in this tutorial. PE: paired-end.

| Cellular state | Datatype | Description | Replicate | ENCODE Accession | Library type | Read length | Stranded? |
| ---            | ---      | :-:     | :-:       | ---           | :-:          | :-:         | :-:       |
| HepG2            | eCLIP | RBFOX2   | 1         | ENCSR987FTF     | PE           | 175-300          | Yes        |
| HepG2            | eCLIP | input   | 1         | ENCSR799EKA     | PE           | 175-300          | Yes        |


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Finding Binding Motifs for RBFOX2

RBFOX2 is a relevant development and tissue-specific splicing factor with a very conserved motif: `TGCATG`. We therefore want to process and validate the data to find this conserved motif and in the process identify the function of RBFOX2 as well as describe the function of the targeted RNA.

## Step 1: Get data

> ### {% icon hands_on %} Hands-on: Data upload
> 1. Create and name a new history for this tutorial.
> 2. Import the following files from [Zenodo](https://zenodo.org/record/1327423) or from a data
>    library named `TODO` if available (ask your instructor)
>
>    ```
>    https://zenodo.org/api/files/102d29d5-2180-490b-be7c-bb0e4ca7b109/hg19_chr_sizes.txt
>    https://zenodo.org/api/files/102d29d5-2180-490b-be7c-bb0e4ca7b109/Homo_sapiens.GRCh37.74.gtf
>    https://zenodo.org/api/files/102d29d5-2180-490b-be7c-bb0e4ca7b109/RBFOX2-204-CLIP_S1_R1_RBFOX2.fastq
>    https://zenodo.org/api/files/102d29d5-2180-490b-be7c-bb0e4ca7b109/RBFOX2-204-CLIP_S1_R2_RBFOX2.fastq
>    https://zenodo.org/api/files/102d29d5-2180-490b-be7c-bb0e4ca7b109/RBFOX2-204-INPUT_S2_R1.fastq
>    https://zenodo.org/api/files/102d29d5-2180-490b-be7c-bb0e4ca7b109/RBFOX2-204-INPUT_S2_R2.fastq
>    ```
>
>    > ### {% icon tip %} Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * Press **Start**
>    >
>    > By default, Galaxy uses the url as the name, so please rename them to something more pleasing.
>    {: .tip}
>
>    ![upload](../../images/upload_data_page.png "Data can be imported directly with links.")
>
>   ![data](../../images/clipseq_data_uploaded.png "Imported datasets will appear in the history panel.")
>    > ### {% icon tip %} Tip: Importing data from a data library
>    >
>    > * Go into "Shared data" (top panel) then "Data libraries"
>    > * Click on "Training data" and then "Transcriptomics"
>    > * Select interesting file
>    > * Click on "Import selected datasets into history"
>    > * Import in a new history
>    {: .tip}
>
{: .hands_on}

# Step 2: Quality Control

As for any NGS data analysis, CLIP-seq data must be quality controlled before being aligned to a reference genome. For more detailed information on NGS quality control, check out the tutorial [here]({{site.baseurl}}/topics/sequence-analysis).

## Report with **FastQC**

> ### {% icon hands_on %} Hands-on: Quality control with FastQC
>
> 1. **FastQC** {% icon tool %}: Run the tool **FastQC** on each FASTQ file to assess the quality of the raw data. An explanation of the results can be found on the [FastQC web page](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
>
>    > ### {% icon tip %} Tip: Running a tool on multiple data files
>    >
>    > You can run this tool - and many other tools - on all the FASTQ files at once!
>    > To do this, first select the "Multiple datasets" icon (two stacked pages) under the "Input FASTQ file" heading in the **FASTQC** Tool Form, then shift+click to select multiple FASTQ files.
>    {: .tip}
>
>   Check the **Sequence Duplication Levels** plot.   
>
>   ![fastqbefore](../../images/clipseq_duplication_level_1.png "Sequence duplication levels <b>before</b> de-duplication.")
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What does the y-axis represent in Figure 3?
>    > 2. What is the meaning of the read and blue line?
>    > 3. What does the headline of Figure 3 tell you?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. PCR duplication occur naturally in any NGS experiment during the PCR amplification of the genetic material. CLIP-Seq is prone to many PCR duplicates because of the sparse material that is obtained during a CLIP-Seq experiment resulting in many occasions in high PCR cycles. The y-axis in Figure 3 represents the portion of reads with the specific duplication level. An exact sequence match is needed to detect a duplicated read. More information can be found here: http://www.cureffi.org/2012/12/11/how-pcr-duplicates-arise-in-next-generation-sequencing/ .
>    > > 2. The blue line shows the duplication levels distribution of the full sequence set. The red line depicts an ideal curve after a de-duplication step (duplicates filtered out). Spikes in the red line come from different duplication levels in the original data (blue line).
>    > > 3. The headline states an expected value of reads that would remain after duplicated reads would be filtered out. A high percentage suggest, that no de-duplication step is needed. This should also correspond with the red line.
>    > {: .solution }
>    {: .question}
{: .hands_on}

# Step 2: Removal of Adapters, Barcodes and Unique Molecular Identifiers (UMIs)

It is often necessary to remove adapter and barcodes sequences as well as UMIs. <br/>
**Adapters** (or primers) are needed for PCR amplification and sequencing in a standard NGS protocol. Unfortunately, it might happen during the sequencing that the machine does not stop at the read end and sequences through the adapter as well. That is why, we need to check if our reads contain those sequences which we are then cutting out.<br/>
**Barcodes** on the other hand are especially designed for a read library and intentionally sequenced. Sometimes experiments are sequenced at the same time which is called **multiplexing**. Multiplexing allows for a better data normalization and comparison. The barcodes are then used to divide the un-multiplexed data set into the individual read libraries. (**Note: Our data is already de-multiplexed, i.e., we do not have to take barcode seqeunces into account.**)<br/>
**UMIs** are similar to barcodes but these sequences are unique for each read. UMIs were introduced since iCLIP to deal with the high duplication levels of a CLIP experiment. Because each read contain an UMI, PCR duplicates of that read also contain the same UMI, which makes it possible to fuse all reads with the same UMI.

## Removal of Adapter Sequences with **Cutadapt**

In this task we are going to remove two 3' and two 5' adapters from the reads (Note: The eCLIP protocol uses more adapter sequences, for more information take a look [here](http://dx.doi.org/10.1038/nmeth.3810)). Because **Cutadapt** can only process one site of the read pair, we have to trigger **Cutadapt** twice.

> ### {% icon hands_on %} Hands-on: Adapter Removal
>
> 1. **Cutadapt (v. 1.6)** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Fastq file to trim"*: `R1`
>    - *"Track Paired Reads"*: `Yes`
>    - {% icon param-file %} *"Paired fastq file (NOT trimmed)"*: `R2`
>    - In *"3' Adapters"*:
>        - Click on *"Insert 3' Adapters"*:
>        - In *"1: 3' Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 3' adapter sequence"*: `AACTTGTAGATCGGA`
>        - Click on *"Insert 3' Adapters"*:
>        - In *"2: 3' Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 3' adapter sequence"*: `AGGACCAAGATCGGA`
>    - In *"5' (Front) Adapters"*:
>        - Click on *"Insert 5' (Front) Adapters"*:
>        - In *"1: 5' (Front) Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 5' adapter sequence"*: `CTTCCGATCTACAAGTT`
>        - Click on *"Insert 5' (Front) Adapters"*:
>        - In *"2: 5' (Front) Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 5' adapter sequence"*: `CTTCCGATCTTGGTCCT`
>    - *"Minimum overlap length"*: `"5"`
>    - *"Output filtering options"*: `Set Filters`
>        - *"Minimum length"*: `10`
>    - *"Additional output options"*: `Default`
>    - *"Additional modifications to reads"*: `Set Modification Options`
>        - *"Cut bases from reads before adapter trimming"*: `-5`
>
>    > ### {% icon comment %} Why do we remove 5 bp from the first read?
>    >
>    > In eCLIP it can happen that the sequencing goes over the first read into the UMI, which is at the 3' end of the first read. The UMI is 5 bp long in our data. To make sure our first read in the read pair does not contain the UMI, we simply remove the last 5 bp from it. The UMI that we actually need for the de-duplication is located on the 5' end of our second read in the read pair.
> {: .comment}
>
> 1. **Cutadapt (v. 1.6)** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Fastq file to trim"*: `R2`
>    - *"Track Paired Reads"*: `Yes`
>    - {% icon param-file %} *"Paired fastq file (NOT trimmed)"*: `R1`
>    - In *"3' Adapters"*:
>        - Click on *"Insert 3' Adapters"*:
>        - In *"1: 3' Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 3' adapter sequence"*: `AACTTGTAGATCGGA`
>        - Click on *"Insert 3' Adapters"*:
>        - In *"2: 3' Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 3' adapter sequence"*: `AGGACCAAGATCGGA`
>    - In *"5' (Front) Adapters"*:
>        - Click on *"Insert 5' (Front) Adapters"*:
>        - In *"1: 5' (Front) Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 5' adapter sequence"*: `CTTCCGATCTACAAGTT`
>        - Click on *"Insert 5' (Front) Adapters"*:
>        - In *"2: 5' (Front) Adapters"*:
>            - *"Source"*: `Enter custom sequence`
>                - *"Enter custom 5' adapter sequence"*: `CTTCCGATCTTGGTCCT`
>    - *"Minimum overlap length"*: `"5"`
>    - *"Output filtering options"*: `Set Filters`
>        - *"Minimum length"*: `10`
>    - *"Additional output options"*: `Default`
>    - *"Additional modifications to reads"*: `No Read Modifications`
>
>    > ### {% icon comment %} Do the Same thing for the input control data set.
>    >
>    > If you processed the RBFOX2 fastq dataset then do the same thing for input control data set or *vice verca*.
>    {: .comment}
>
{: .hands_on}

## Removal of UMIs with **UMI-tools extract**

In this task we are going to remove the 5 bp UMI in the 5' end of the second read.

> ### {% icon hands_on %} Hands-on: UMI Removal
>
> 1. **UMI-tools extract** {% icon tool %} with the following parameters:
>    - *"Library type"*: `Paired-end`
>     - {% icon param-file %} *"Reads in FASTQ format"*: `R1 from Cutadapt output`
>     - {% icon param-file %} *"Reads in FASTQ format"*: `R2 from Cutadapt output`
>     - *"Barcode on both reads?"*: `Barcode on first read only`
>    - *"Use Known Barcodes?"*: `No`
>    - *"Method to extract barcodes"*: `String`
>    - *"Barcode pattern for first read"*: `"NNNNN"`
>    - *"Is the barcode at the 5' end?"*: `Yes`
>    - *"Output log?"*: `Yes`
>    - *"Enable quality filter?"*: `No`
>
>    > ### {% icon comment %} Do the Same thing for the input control data set.
>    >
>    > If you processed the RBFOX2 fastq dataset then do the same thing for input control data set or *vice verca*.
>    {: .comment}
>
{: .hands_on}

# Step 3: Aligning Reads to a Reference Genome

To determine where DNA fragments originated in the genome, the sequenced reads must be aligned to a reference genome. This is equivalent to solving a jigsaw puzzle, but unfortunately, not all pieces are unique. In principle, you could do a BLAST analysis to figure out where the sequenced pieces fit best in the known genome. Aligning millions of short sequences this way, however, this can take a couple of weeks.
Nowadays, there are many read alignment programs, `STAR` is one of them that works well with CLIP-Seq data, for more information read  [here](doi:10.1093/bioinformatics/bts635). STAR is able to use genome as well as transcriptome data. This ability is handy, since CLIP-Seq generetas transcriptome data, thus, we have to take RNA processing steps like splicing events into account.

## Aligning with **RNA STAR**

> ### {% icon hands_on %} Hands-on: Alignment
>
> 1. **RNA STAR** {% icon tool %} with the following parameters:
>    - *"Single-end or paired-end reads"*: `Paired-end (as individual datasets)`
>    - *"Custom or built-in reference genome"*: `Use a built-in index`
>        - *"Reference genome with or without an annotation"*: `use genome reference with builtin gene-model`
>            - *"Select reference genome"*: `Homo sapiens (hg19+GRCh37.75)`
>    - *"Count number of reads per gene"*: `No`
>    - *"Would you like to set output parameters (formatting and filtering)?"*: `No`
>    - *"Other parameters (seed, alignment, limits and chimeric alignment)"*: `Extended parameter list`
>        - In *"Alignment parameters"*:
>            - *"Use end-to-end read alignments, with no soft-clipping?"*: `Yes`
>        - *"Would you like to set chimeric alignment parameters?"*: `No`
>
>    > ### {% icon comment %} Soft-Clipping vs Hard-Clipping
>    >
>    > Clipping is a way to deal with low quality bases during the alignment step. In **Soft-Clipping** the bases at the 5' and 3 end of the read are not part of the alignment. In **Hard-Clipping** the bases at the 5' and 3' end of the read are not part of the alignment **and** will be completely removed from the read sequence in the BAM file.
>    {: .comment}
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. Can you think of a reason why we disabled the soft-clipping?
>
> > ### {% icon solution %} Solution
> >
> > 1. In eCLIP the crosslinking position should be at the beginning of the second read. If we would enable soft-clipping, we would add potential bases with low quality at the end of our second reads that would blur our crosslinking position and we would lose precision to detect potential binding regions of RBFOX2.
> >
> {: .solution}
>
{: .question}

# Step 4: De-Duplication

More information on **UMI-tools** can be found [here](10.1101/gr.209601.116).

## Sub-step with **UMI-tools deduplicate**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **UMI-tools deduplicate** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Reads to deduplicate in SAM or BAM format"*: `mapped_reads` (output of **RNA STAR** {% icon tool %})
>    - *""*: ``
>    - *"Separator between read id and UMI."*: `"_"`
>    - *"Tag which contains UMI."*: `""`
>    - *"Method used to identify PCR duplicates within reads."*: ``
>    - *"Edit distance threshold"*: `"1"`
>    - *"BAM is paired end"*: `Yes`
>    - *"Spliced reads are unique"*: `Yes`
>    - *"Soft clip threshold"*: `"4"`
>    - *"Use the read length as as a criterion when deduping"*: `Yes`
>    - *"Consider all alignments to a single contig together"*: `Yes`
>    - *"Only consider a random selection of the reads"*: `"1.0"`
>    - *"Only consider a single chromosome"*: `Yes`
>    - *"Deduplicate per contig"*: `Yes`
>    - *"Deduplicate per gene"*: `Yes`
>    - *"Deduplicate by this gene tag"*: `""`
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

## Sub-step with **FastQC**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **FastQC** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Short read data from your current history"*: `output` (output of **UMI-tools deduplicate** {% icon tool %})
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

# Step 5: Second Quality Control

## Sub-step with **multiBamSummary**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **multiBamSummary** {% icon tool %} with the following parameters:
>    - *"Sample order matters"*: `No`
>    - *"Choose computation mode"*: `Bins`
>    - *"Region of the genome to limit the operation to"*: `""`
>    - *"Show advanced options"*: `no`
>    - *"Save raw counts (coverages) to file"*: `Yes`
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

## Sub-step with **plotFingerprint**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **plotFingerprint** {% icon tool %} with the following parameters:
>    - *"Sample order matters"*: `No`
>    - *"Region of the genome to limit the operation to"*: `""`
>    - *"Show advanced options"*: `no`
>    - *"Show advanced output settings"*: `no`
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

## Sub-step with **plotCorrelation**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **plotCorrelation** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Matrix file from the multiBamSummary tool"*: `outFile` (output of **multiBamSummary** {% icon tool %})
>    - *"Correlation method"*: ``
>    - *"Plotting type"*: `Heatmap`
>    - *"Skip zeros"*: `Yes`
>    - *"Image file format"*: ``
>    - *"Remove regions with very large counts"*: `Yes`
>    - *"Save the matrix of values underlying the heatmap"*: `Yes`
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

# Step 6: Peakcalling

## Sub-step with **PEAKachu**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PEAKachu** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Experiment Libraries"*: `output` (output of **UMI-tools deduplicate** {% icon tool %})
>    - {% icon param-file %} *"Control Libraries"*: `output` (output of **UMI-tools deduplicate** {% icon tool %})
>    - *"Pairwise Replicates"*: `Yes`
>    - *"Paired End"*: `Yes`
>    - *"Maximum Insert Size"*: `"200"`
>    - *"Features"*: `""`
>    - *"Sub-Features"*: `""`
>    - *"Select Mode"*: `Adaptive`
>        - *"Normalisation Method."*: `DESeq2`
>    - *"Mad Multiplier"*: `"0.0"`
>    - *"Fold Change Threshold"*: `"2.0"`
>    - *"Adjusted p-value Threshold"*: `"0.05"`
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

## Sub-step with **Extract alignment ends**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Extract alignment ends** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Alignments in SAM or BAM format"*: `output` (output of **UMI-tools deduplicate** {% icon tool %})
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

## Sub-step with **Text reformatting**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Text reformatting** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `peak_tables` (output of **PEAKachu** {% icon tool %})
>    - *"AWK Program"*: `"NR>1{\nif ($3 < $4) {\n   print $1,$3,$4,\"clip_peak_\"NR-1,$9,$5;\n}\nelse {\n   print $1,$4,$3,\"clip_peak_\"NR-1,$9,$5;\n}\n}"`
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

## Sub-step with **SortBED**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **SortBED** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sort the following BED file"*: `alignment_ends` (output of **Extract alignment ends** {% icon tool %})
>    - *"Sort by"*: ``
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

## Sub-step with **Get crosslinked nucleotides**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Get crosslinked nucleotides** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Alignments in BED format"*: `alignment_ends` (output of **Extract alignment ends** {% icon tool %})
>    - *"Set position one nt downstream of 3'-end as crosslinked nucleotide"*: `Yes`
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

## Sub-step with **SlopBed**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **SlopBed** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"BED/VCF/GFF file"*: `outfile` (output of **Text reformatting** {% icon tool %})
>    - *"Genome file"*: `Genome file from your history`
>    - *"Define -l and -r as a fraction of the feature’s length"*: `Yes`
>    - *"Define -l and -r based on strand"*: `Yes`
>    - *"Choose what you want to do"*: `Increase the BED/GFF/VCF entry by the same number base pairs in each direction.`
>        - *"Number of base pairs"*: `20`
>    - *"Print the header from the A file prior to results"*: `Yes`
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

## Sub-step with **Create a BedGraph of genome coverage**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Create a BedGraph of genome coverage** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"The BAM or BED file from which coverage should be computed"*: `output` (output of **SortBED** {% icon tool %})
>    - *"Report regions with zero coverage"*: `Yes`
>    - *"Treat split/spliced BAM or BED12 entries as distinct BED intervals when computing coverage."*: `Yes`
>    - *"Calculate coverage based on"*: ``
>    - *"Scale the coverage by a constant factor"*: `""`
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

## Sub-step with **Text reformatting**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Text reformatting** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `crosslinking_coordinates` (output of **Get crosslinked nucleotides** {% icon tool %})
>    - *"AWK Program"*: `"$2!=0"`
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

## Sub-step with **Text reformatting**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Text reformatting** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `output` (output of **SlopBed** {% icon tool %})
>    - *"AWK Program"*: `"{print $0}"`
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

## Sub-step with **Extract Genomic DNA**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Extract Genomic DNA** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Fetch sequences for intervals in"*: `output` (output of **SlopBed** {% icon tool %})
>    - *"Interpret features when possible"*: ``
>    - *"Choose the source for the reference genome"*: `locally cached`
>        - *"Using reference genome"*: ``
>    - *"Select output format"*: ``
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

## Sub-step with **Wig/BedGraph-to-bigWig**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Wig/BedGraph-to-bigWig** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Convert"*: `output` (output of **Create a BedGraph of genome coverage** {% icon tool %})
>    - *"Converter settings to use"*: `Default`
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

## Sub-step with **SortBED**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **SortBED** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sort the following BED file"*: `outfile` (output of **Text reformatting** {% icon tool %})
>    - *"Sort by"*: ``
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

## Sub-step with **RNA Centric Annotation System**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **RNA Centric Annotation System** {% icon tool %} with the following parameters:
>    - *"Genome Version"*: ``
>    - {% icon param-file %} *"Target regions in BED format"*: `outfile` (output of **Text reformatting** {% icon tool %})
>    - {% icon param-file %} *"Reference annotation in ENSEMBL GTF format"*: `output` (Input dataset)
>    - *"Run annotation."*: `Yes`
>    - *"Run GO term enrichment"*: `Yes`
>    - *"Run gene set enrichment"*: `No`
>    - *"Run motif search"*: `Yes`
>    - *"Downsampling (N)"*: `"0"`
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

## Sub-step with **MEME-ChIP**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **MEME-ChIP** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Primary sequences"*: `output` (output of **Extract Genomic DNA** {% icon tool %})
>    - *"Sequence alphabet"*: ``
>    - *"Options Configuration"*: `Advanced`
>        - *"Limit of sequences to pass to MEME"*: `100`
>        - *"Should subsampling be random?"*: `Yes`
>            - *"Seed for the randomized selection of sequences"*: `123`
>        - *"maximum size of a sequence before it is cut down to a centered section"*: `0`
>        - *"Search given strand only"*: `Yes`
>        - *"What is the expected motif site distribution?"*: `Zero or one occurances per sequence`
>        - *"Minimum motif width"*: `5`
>        - *"Maximum motif width"*: `20`
>        - *"Maximum number of motifs to find"*: `20`
>        - *"Stop DREME searching after finding this many motifs"*: `5`
>    - *"I certify that I am not using this tool for commercial purposes."*: `Yes`
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

## Sub-step with **Create a BedGraph of genome coverage**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Create a BedGraph of genome coverage** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"The BAM or BED file from which coverage should be computed"*: `output` (output of **SortBED** {% icon tool %})
>    - *"Report regions with zero coverage"*: `Yes`
>    - *"Treat split/spliced BAM or BED12 entries as distinct BED intervals when computing coverage."*: `Yes`
>    - *"Calculate coverage based on"*: ``
>    - *"Scale the coverage by a constant factor"*: `""`
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

## Sub-step with **Wig/BedGraph-to-bigWig**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Wig/BedGraph-to-bigWig** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Convert"*: `output` (output of **Create a BedGraph of genome coverage** {% icon tool %})
>    - *"Converter settings to use"*: `Default`
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


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
