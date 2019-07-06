---
layout: tutorial_hands_on

title: ATAC-Seq data analysis
zenodo_link: https://zenodo.org/record/3270536
questions:
- Which biological questions are addressed by the tutorial?
- Which bioinformatics techniques are important to know for this type of data?
objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentences describing what a learner should be able to do once they
  have completed the tutorial
- You can use Bloom's Taxonomy to write effective learning objectives
time_estimation: ''
key_points:
- The take-home messages
- They will appear at the end of the tutorial
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

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Preprocessing

Give some background about what the trainees will be doing in the section.

Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
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
> 3. Rename the datasets
> 4. Check that the datatype is fasqsanger.gz
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

# Title of the section usually corresponding to a big step in the analysis

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


## Quality control

The first step is to check the quality of the reads and the presence of the Nextera adapters. We can do this with FastQC. The FastQC web page Adapter Content section shows the presence of Nextera Transposase Sequence in the reads. We can remove the adapters with Cutadapt.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **FastQC** {% icon tool %} with the following parameters:
>       - {% icon param-files %} *"Short read data from your current history"*: input datasets selected with **Multiple datasets**
>
>    {% include snippets/select_multiple_datasets.md %}
>
> 2. Inspect the webpage output of **FastQC** {% icon tool %} for the `SRR891268_R1` sample
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
>    {: .question}


## Sub-step with **Cutadapt**

To trim the adapters we provide the Nextera adapter sequences to Cutadapt. Thesse adapters are shown in the image below.

![Nextera adapters](../../images/atac-seq/nexteraLibraryPicture.svg "Nextera adapters")

The forward and reverse adapters are slightly different. If we run FastQC again we should see under Adapter Content that the Nextera adapters are no longer present. We will also trim low quality bases at the ends of reads (quality less than 20). We will only keep reads that are 20 bases or more after trimming so we don't have reads that are too short.

> ### {% icon hands_on %} Hands-on: Trim reads
>
> 1. **Cutadapt** {% icon tool %} with the following parameters:
>    - *"Single-end or Paired-end reads?"*: `Paired-end`
>        - In *"Read 1 Options"*:
>            - In *"3' (End) Adapters"*:
>                - {% icon param-repeat %} *"Insert 3' (End) Adapters"*
>                    - *"Source"*: `Enter custom sequence`
>                        - *"Enter custom 3' adapter name (Optional if Multiple output is 'No')"*: `Nextera R1`
>                        - *"Enter custom 3' adapter sequence"*: `CTGTCTCTTATACACATCTCCGAGCCCACGAGAC`
>        - In *"Read 2 Options"*:
>            - In *"3' (End) Adapters"*:
>                - {% icon param-repeat %} *"Insert 3' (End) Adapters"*
>                    - *"Source"*: `Enter custom sequence`
>                        - *"Enter custom 3' adapter name (Optional)"*: `Nextera R2`
>                        - *"Enter custom 3' adapter sequence"*: `CTGTCTCTTATACACATCTGACGCTGCCGACGA`
>    - In *"Filter Options"*:
>        - *"Minimum length"*: `20`
>    - In *"Read Modification Options"*:
>        - *"Quality cutoff"*: `20`
>    - In *"Output Options"*:
>        - *"Report"*: `Yes`
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

## Sub-step with **Convert GTF to BED12**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Convert GTF to BED12** {% icon tool %} with the following parameters:
>    - *"Advanced options"*: `Use default options`
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

# Mapping

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

## Sub-step with **bedtools SortBED**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **bedtools SortBED** {% icon tool %} with the following parameters:
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

## Sub-step with **Filter**

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

## Sub-step with **MarkDuplicates**

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

## Sub-step with **CollectInsertSizeMetrics**

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

## Sub-step with **Convert, Merge, Randomize**

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

## Sub-step with **MACS2 callpeak**

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

## Sub-step with **Wig/BedGraph-to-bigWig**

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

## Sub-step with **Cut**

> ### {% icon hands_on %} Hands-on: Task description
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

## Sub-step with **computeMatrix**

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

## Sub-step with **pyGenomeTracks**

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

## Sub-step with **plotHeatmap**

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