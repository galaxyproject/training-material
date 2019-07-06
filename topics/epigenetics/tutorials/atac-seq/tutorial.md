---
layout: tutorial_hands_on

title: ATAC-Seq data analysis
zenodo_link: "https://zenodo.org/record/3270536"
questions:
- Which DNA regions are accesible in the human lymphoblastoid cell line GM12878?
- How to analyse and visualise ATAC-Seq data?
objectives:
- Apply appropriate analysis and quality control steps for ATAC-Seq
- Visualise coverage and peaks of specific regions
- Generate heatmaps
time_estimation: 'No idea'
key_points:
- ATAC-Seq can be used to identify accessible gene promoters and enhancers
- ATAC-seq analysis uses similar to ChIP-Seq but different parameters
contributors:
- lldelisle
- mblue9
- heylf
- glorybsm

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

General introduction about the topic and then an introduction of the
tutorial (the questions and the objectives). It is nice also to have a
scheme to sum up the pipeline used during the tutorial. The idea is to
give to trainees insight into the content of the tutorial and the (theoretical
and technical) key concepts they will learn.

**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

![ATAC-Seq](../../images/atac-seq/atac-seq.jpeg "Buenrostro et al. 2013 Nat Methods")

The original dataset had 2 x 200 million reads. This would be too long to be processed in a training session. So, we downsampled the original dataset to 200 000 reads but added about 200 000 reads pairs that will map to chr22 to have a good profile on this chromosome similar to a 2 x 20 million reads original fastq.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Preprocessing

We first need to download the dataset that we downsampled as well as other annotations files. Then, to increase the number of reads that will map to the assembly (here Human genome version 38), we need to preprocess the reads.

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
>
>    {% include snippets/create_new_history.md %}
>
> 2. Import the files from [Zenodo](https://zenodo.org/record/3270536) or from the shared data library
>
>    ```
>    https://zenodo.org/record/3270536/files/SRR891268_R1.fastq.gz
>    https://zenodo.org/record/3270536/files/SRR891268_R2.fastq.gz
>    ```
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. By default, galaxy will give as name the full path. Rename the datasets to keep only the file names.
>
>    {% include snippets/rename_dataset.md %}
>
> 4. Check that the datatype of the 2 fastq files is fasqsanger.gz
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> > ### {% icon details %} FASTQ format
> > If you are not familiar with FASTQ format, see the [Quality Control tutorial]({{ site.baseurl }}{% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %})
{: .details}
{: .hands_on}

***TODO***: Add the bed files links.

***TODO***: *Consider adding a detail box to expand the theory*

## Quality control

The first step is to check the quality of the reads and the presence of the Nextera adapters.
> ### {% icon details %} Presence of adapters
> When we do ATAC-Seq, two transposase could cut the DNA about 40 bp apart. This can be smaller than the sequencing length so we expect to have Nextera adapters at the end these reads. 
{: .details}
We can do this with FastQC. The FastQC web page Adapter Content section shows the presence of Nextera Transposase Sequence in the reads. We will remove the adapters with Cutadapt.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **FastQC** {% icon tool %} with the default parameters:
>       - *"Short read data from your current history"*: Choose here either only the `SRR891268_R1.fastq.gz` file with {% icon param-file %} or use {% icon param-files %} **Multiple datasets** to choose both `SRR891268_R1.fastq.gz` and `SRR891268_R2.fastq.gz`.
>
>    {% include snippets/select_multiple_datasets.md %}
>
> 2. Inspect the webpage output of **FastQC** {% icon tool %} for the `SRR891268_R1` sample. Check which are the adaptors found at the end of the reads.
>
>    > ### {% icon question %} Questions
>    >
>    > What is the read length?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > The read length is 50 bp.
>    > >
>    > {: .solution}
>    >
>    >
>    > How many reads are in the fastq?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > There are 285247 reads.
>    > >
>    > {: .solution}
>    >
>    >
>    > Which are the steps which have a warning?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > 1) Per base sequence content
>    > > > ### {% icon details %} Tn5 sequence bias
>    > > > It is well known that the tn5 have a strong sequence bias at the insersion site. You can have more information about it reading [Green et al. 2012]({% link https://mobilednajournal.biomedcentral.com/articles/10.1186/1759-8753-3-3 %})
{: .details}
>    > > 2) Sequence Duplication Levels
>    > > 3) Overrepresented sequences
>    > >
>    > {: .solution}
>    >
>    {: .question}

This is what you should expect from the adapter section:
![FastQC screenshot on the adapter content section](../../images/atac-seq/Screenshot_fastqcBeforecutadapt.png "FastQC screenshot on the adapter content section")

## Trimming reads

To trim the adapters we provide the Nextera adapter sequences to Cutadapt. Thesse adapters are shown in the image below.

![Nextera library with the sequence of adapters](../../images/atac-seq/nexteraLibraryPicture.svg "Nextera library with the sequence of adapters")

The forward and reverse adapters are slightly different.We will also trim low quality bases at the ends of reads (quality less than 20). We will only keep reads that are 20 bases or more after trimming so we don't have reads that are too short to be mapped.

> ### {% icon hands_on %} Hands-on: Trim reads
>
> 1. **Cutadapt** {% icon tool %} with the following parameters:
>    - *"Single-end or Paired-end reads?"*: `Paired-end`
>        - {% icon param-file %} *"FASTQ/A file #1"*: select `SRR891268_R1.fastq.gz`
>        - {% icon param-file %} *"FASTQ/A file #2"*: select `SRR891268_R2.fastq.gz`
>        - In *"Read 1 Options"*:
>            - In *"3' (End) Adapters"*:
>                - {% icon param-repeat %} *"Insert 3' (End) Adapters"*
>                    - *"Source"*: `Enter custom sequence`
>                        - {% icon param-text %} *"Enter custom 3' adapter name (Optional if Multiple output is 'No')"*: `Nextera R1`
>                        - {% icon param-text %} *"Enter custom 3' adapter sequence"*: `CTGTCTCTTATACACATCTCCGAGCCCACGAGAC`
>        - In *"Read 2 Options"*:
>            - In *"3' (End) Adapters"*:
>                - {% icon param-repeat %} *"Insert 3' (End) Adapters"*
>                    - *"Source"*: `Enter custom sequence`
>                        - {% icon param-text %} *"Enter custom 3' adapter name (Optional)"*: `Nextera R2`
>                        - {% icon param-text %} *"Enter custom 3' adapter sequence"*: `CTGTCTCTTATACACATCTGACGCTGCCGACGA`
>    - In *"Filter Options"*:
>        - {% icon param-text %} *"Minimum length"*: `20`
>    - In *"Read Modification Options"*:
>        - {% icon param-text %} *"Quality cutoff"*: `20`
>    - In *"Output Options"*:
>        - *"Report"*: `Yes`
>
> 2. Click on the {% icon galaxy-eye %} (eye) icon of the report and read the first lines.
{: .hands_on}

This is what you should see:
![Summary of cutadapt](../../images/atac-seq/Screenshot_cutadaptSummary.png "Summary of cutadapt")

> ### {% icon question %} Questions
>
> 1. Which portion of reads contains adapters?
> 2. How many reads were still longer than 20 after the trimming?
>
> > ### {% icon solution %} Solution
> >
> > 1. 15%
> > 2. 284,864 (99.9%)
> >
> {: .solution}
>
{: .question}

 If we run FastQC again we should see under Adapter Content that the Nextera adapters are no longer present. 

# Mapping

## Mapping reads to reference genome

Next we map the trimmed reads to the human reference genome. We could use BWA or Bowtie to do this, here we will use Bowtie2. We use the settings...


> ### {% icon hands_on %} Hands-on: Mapping reads to reference genome
>
> 1. **Bowtie2** {% icon tool %} with the following parameters:
>    - *"Is this single or paired library"*: `Paired-end`
>        - *"Do you want to set paired-end options?"*: `Yes`
>            - *"Set the maximum fragment length for valid paired-end alignments"*: `1000`
>            - *"Allow mate dovetailing"*: `Yes`
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a built-in genome index`
>        - *"Select reference genome"*: `Human Dec. 2013 (GRCh38/hg38) (hg38)`
>    - *"Set read groups information?"*: `Do not set`
>    - *"Select analysis mode"*: `1: Default setting only`
>        - *"Do you want to use presets?"*: `Very sensitive end-to-end (--very-sensitive)`
>    - *"Do you want to tweak SAM/BAM Options?"*: `No`
>    - *"Save the bowtie2 mapping statistics to the history"*: `Yes`
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

# Filtering mapped reads

## Filter uninformative reads

We apply some filters to the reads after mapping. ATAC-seq datasets can have a lot of reads mapping to mitchondria and we remove these reads. We also remove reads with low mapping quality and that are not properly paired.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Filter** {% icon tool %} with the following parameters:
>    - In *"Condition"*:
>        - {% icon param-repeat %} *"Insert Condition"*
>            - In *"Filter"*:
>                - {% icon param-repeat %} *"Insert Filter"*
>                    - *"Select BAM property to filter on"*: `mapQuality`
>                        - *"Filter on read mapping quality (phred scale)"*: `>=30`
>                - {% icon param-repeat %} *"Insert Filter"*
>                    - *"Select BAM property to filter on"*: `isProperPair`
>                        - *"Select properly paired reads"*: `Yes`
>                - {% icon param-repeat %} *"Insert Filter"*
>                    - *"Select BAM property to filter on"*: `reference`
>                        - *"Filter on the reference name for the read"*: `!chrM`
>    - *"Would you like to set rules?"*: `Yes`
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

## Filter duplicate reads

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **MarkDuplicates** {% icon tool %} with the following parameters:
>    - In *"Comment"*:
>        - {% icon param-repeat %} *"Insert Comment"*
>            - *"Add this comment to BAM dataset"*: `[`
>        - {% icon param-repeat %} *"Insert Comment"*
>            - *"Add this comment to BAM dataset"*: `]`
>    - *"If true do not write duplicates to the output file instead of writing them with appropriate flags set"*: `Yes`
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

## Check fragment sizes

We check the fragment sizes (the sizes of the pieces of DNA that were sequenced) with Picard CollectInsertSizeMetrics.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **CollectInsertSizeMetrics** {% icon tool %} with the following parameters:
>    - *"Load reference genome from"*: `Local cache`
>        - *"Using reference genome"*: `Human Dec. 2013 (GRCh38/hg38) (hg38)`
>    - *"The level(s) at which to accumulate metrics"*: ``
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

# Peak calling

## Prepare input reads

We convert the BAM file to BED format because when shifting reads with MACS2 will only consider one of the read pairs.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Convert, Merge, Randomize** {% icon tool %} with the following parameters:
>    - *"Select BAM manipulation"*: `Convert`
>        - *""*: `BED`
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

## Call peaks

We call peaks with a peak caller such as MACS2 or Genrich. Here we will use MACS2.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **MACS2 callpeak** {% icon tool %} with the following parameters:
>    - *"Are you pooling Treatment Files?"*: `No`
>    - *"Do you have a Control File?"*: `No`
>    - *"Format of Input Files"*: `Single-end BED`
>    - *"Effective genome size"*: `H. sapiens (2.7e9)`
>    - *"Build Model"*: `Do not build the shifting model (--nomodel)`
>        - *"Set shift size"*: `-100`
>    - *"Peak detection based on"*: `q-value`
>    - *"Additional Outputs"*: ``
>    - In *"Advanced Options"*:
>        - *"Composite broad regions"*: `No broad regions`
>            - *"Use a more sophisticated signal processing approach to find subpeak summits in each enriched peak region"*: `Yes`
>        - *"How many duplicate tags at the exact same location are allowed?"*: `all`
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

# Visualisation of coverage

## Visualise region of interest

> ### {% icon hands_on %} Hands-on: Convert peaks to BED3 format
>
> 1. **Cut** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"File to cut"*: `output_narrowpeaks` (output of **MACS2 callpeak** {% icon tool %})
>    - *"Cut by"*: `fields`
>        - *"List of Fields"*: `c1-3`
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

## Visualise region with **pyGenomeTracks**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **pyGenomeTracks** {% icon tool %} with the following parameters:
>    - *"Region of the genome to limit the operation"*: `chr22:20,477,597-20,545,090`
>    - In *"Include tracks in your plot"*:
>        - {% icon param-repeat %} *"Insert Include tracks in your plot"*
>            - *"Choose style of the track"*: `Bedgraph track `
>                - *"Plot title"*: `coverage from macs2 (extended +/-100bp)`
>                - {% icon param-file %} *"Track file bedgraph format"*: `output_treat_pileup` (output of **MACS2 callpeak** {% icon tool %})
>                - *"Color of track"*: `#c0504d`
>                - *"height"*: `5.0`
>                - *"Show visualization of data range"*: `Yes`
>        - {% icon param-repeat %} *"Insert Include tracks in your plot"*
>            - *"Choose style of the track"*: `Gene track / Bed track`
>                - *"Plot title"*: `Peaks from macs2 (extended +/-100bp)`
>                - {% icon param-file %} *"Track file bed format"*: `output` (output of **Cut** {% icon tool %})
>                - *"Color of track"*: `#c0504d`
>                - *"Plot labels"*: `Yes`
>        - {% icon param-repeat %} *"Insert Include tracks in your plot"*
>            - *"Choose style of the track"*: `Gene track / Bed track`
>                - *"Plot title"*: `Genes`
>                - {% icon param-file %} *"Track file bed format"*: `output` (output of **bedtools SortBED** {% icon tool %})
>                - *"Type"*: `genes`
>        - {% icon param-repeat %} *"Insert Include tracks in your plot"*
>            - *"Choose style of the track"*: `Gene track / Bed track`
>                - *"Plot title"*: `CTCF peaks`
>                - {% icon param-file %} *"Track file bed format"*: `output` (Input dataset)
>                - *"Color of track"*: `#4bacc6`
>                - *"Plot labels"*: `Yes`
>    - *"Configure x-axis"*: `Yes`
>        - *"Where to place the x-axis"*: `Bottom`
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

## Create heatmap of genes

### Convert BedGraph to bigWig

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Wig/BedGraph-to-bigWig** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Convert"*: `output_treat_pileup` (output of **MACS2 callpeak** {% icon tool %})
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

## Generate computeMatrix

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **computeMatrix** {% icon tool %} with the following parameters:
>    - In *"Select regions"*:
>        - {% icon param-repeat %} *"Insert Select regions"*
>            - {% icon param-file %} *"Regions to plot"*: `output` (output of **bedtools SortBED** {% icon tool %})
>    - *"Sample order matters"*: `No`
>        - {% icon param-file %} *"Score file"*: `out_file1` (output of **Wig/BedGraph-to-bigWig** {% icon tool %})
>    - *"computeMatrix has two main output options"*: `reference-point`
>    - *"Show advanced output settings"*: `no`
>    - *"Show advanced options"*: `no`
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


### Run **plotHeatmap**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **plotHeatmap** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Matrix file from the computeMatrix tool"*: `outFileName` (output of **computeMatrix** {% icon tool %})
>    - *"Show advanced output settings"*: `no`
>    - *"Show advanced options"*: `no`
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


# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.