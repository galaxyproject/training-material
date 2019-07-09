---
layout: tutorial_hands_on

title: Metatranscriptomics analysis using microbiome RNA-seq data
zenodo_link: https://zenodo.org/record/3269404
questions:
- "How to analyze metatranscriptomics data?"
- "What information can be extracted of metatranscriptomics data?"
- "How to assign taxa and function to the identified sequences?"
objectives:
- "Choosing the best approach to analyze metatranscriptomics data"
- "Exposure to functional microbiome characterization using metatranscriptomic results"
- "Learn where metatranscriptomics fits in 'multi-omic' analysis of microbiomes"
- "Visualisation of a community structure"
time_estimation: 3H
key_points:
- "With shotgun data, we can extract information about the studied community structure and also the functions realised by the community"
- "Metatranscriptomics data analyses are complex and time-consuming"
contributors:
- pratikdjagtap
- subinamehta
- jraysajulga
- bebatut
- emmaleith

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

General introduction about the topic and then an introduction of the
tutorial (the questions and the objectives). It is nice also to have a
scheme to sum up the pipeline used during the tutorial. The idea is to
give to trainees insight into the content of the tutorial and the (theoretical
and technical) key concepts they will learn.

## ASaiM Workflow

![ASaiM diagram](../../images/asaim-wf.png "The Auvergne Sequence analysis of intestinal Microbiota (ASaiM)")

**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Title for your first section

Give some background about what the trainees will be doing in the section.

Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

# Metatranscriptomics for characterizing microbial communities

###  Multi-omics

A multi-omic approach to studying microbiomes can cover the various taxa and functions in a system:

![meta-momics diagram](../../images/meta-omics.png "Metagenomics is packed with information about the present taxonomies in a microbiome, but do not tell much about important functions. That is where metatranscriptomics and metaproteomics play a big part.")

Microbiomes play a critical role in host health, disease, and the environment. Microbiome dynamics, especially at the taxonomy and functional level requires studying the DNA content (metagenomics), RNA expression (metatranscriptomics), protein expression (metaproteomics) or small molecules (metabolomics). Metatranscriptomics analysis enables understanding of how the microbiome responds to the environment by studying the functional analysis of genes expressed by the microbiome. It can also estimate the taxonomic composition of the microbial population.

Batut ​et al.​(​10.1093/gigascience/giy057​) developed ASaiM (Auvergne Sequence analysis of intestinal Microbiota) an open-source Galaxy-based workflow that enables microbiome analyses. ASaiM workflow offers a streamlined Galaxy workflow for users to explore metagenomic/metatranscriptomic data in a reproducible and transparent environment. In this tutorial, we demonstrate the use of an updated version of the ASaiM workflow. This workflow takes in paired-end datasets of raw shotgun sequences (in FastQ format) as an input; preprocess it; extract taxonomic and functional information and combines them to offer insights into taxonomic contribution to a function or functions expressed by a particular taxonomy.

In this ASaiM metatranscriptomics tutorial, the fastqsanger files are used as input dataset. We have used a trimmed version of the original fastq file in this tutorial for the purpose of saving time and resources.


## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>    https://zenodo.org/api/files/84d7d6c9-2b7b-4569-87f6-dabc5ee42bc2/FORWARD_T4A_F.fastqsanger
>    https://zenodo.org/api/files/84d7d6c9-2b7b-4569-87f6-dabc5ee42bc2/REVERSE_T4A_R.fastqsanger
>    ```
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
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


## **FASTQ Groomer**

> ### {% icon hands_on %} Hands-on: The FASTQ Groomer tool provides the user with several formatting options related to the input FASTQ file. In the Default option, the output will be sanger formatted (i.e. file.fastqsanger). If you choose the Advanced option, the user has the capability to choose between any format and force quality scoring i.e. if the score falls outside the target score range, it will be converted to the closest available value.
>
> 1. **FASTQ Groomer** {% icon tool %} with the following parameters:
>    - {% icon param-select %} *"File to groom”*: `T4A_F`
>    - {% icon param-select %} *"Input FASTQ quality scores type"*:`Sanger & Illumina 1.8+`
>    - *"Advanced Options"*: `Hide Advanced Options`
>   
> 2. **Execute**
>
>    > ### {% icon comment %} Comment
>    >
>    > Perform this action on both the files (**Forward and Reverse(T4_R)**)
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. What is the reason behind using FASTQ Groomer?
> 2. Are there any other tools which can perform similar function?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}


## **Cutadapt**

> ### {% icon hands_on %} Hands-on: Cutadapt helps finding and removing adapter sequences, primers, poly-A tails or other unwanted sequences from the input FASTQ files. It trims the input reads by finding the adapter or primer sequences in an error-tolerant way. It also has additional features of modifying and filtering reads.
>
> 1. **Cutadapt** {% icon tool %} with the following parameters:
>    - *"Single-end or Paired-end reads?"*: `Paired-end`
>  - {% icon param-select %} *”FASTQ/A file #1 ”*: `FASTQ groomer on data 1` 
>  - {% icon param-select %} *”FASTQ/A file #2 ”*:`FASTQ groomer on data 2`
>  - {% icon param-select %} *”Adapter options ”*: `Default`
>  - {% icon param-select %} *”Filter options ”*: `Default`
>  - {% icon param-select %} *”Read Modifications options ”*: `Default`
>  - {% icon param-select %} *”Output options ”*: `Default`
>
> 2. **Execute**
>
>    > ### {% icon comment %} Comment
>    >
>    > Keep all the other parameters default
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Are there any other tool which can perform similar functions?
> 2. Reason behind selecting Cutadapt??
>
> > ### {% icon solution %} Solution
> >
> > 1.Yes, there are many tools such as Trimmomatic, Trim Galore, Clip, trim putative adapter sequences. etc.
> > 2. The reason for choosing Cutadapt for this function is because it is error tolerant and the version is pretty stable.
> {: .solution}
>
{: .question}

## **Filter with SortMeRNA**

> ### {% icon hands_on %} Hands-on: SortMeRNA is a tool used to filter ribosomal RNA fragments from metatranscriptomic data produced by next-generation sequencers. It is capable of handling large RNA databases and sorting out all fragments matching to the database with high accuracy and specificity. It can take paired-end reads as well as single-end reads.

>
> 1. **Filter with SortMeRNA** {% icon tool %} with the following parameters:
>
>   - {% icon param-select %} *”Sequencing type”*: `Reads are paired`
>   - {% icon param-file %} *”Forward reads”*: `Read 1: trimmed`
>   - {% icon param-file %} *”Reverse reads”*: `Read 2: trimmed`
>   - {% icon param-check %} *”If one of the paired-end reads aligns and the other one does not”*: `Output both reads to rejected file (--paired_out)`
>   - {% icon param-select %} *”Which strands to search”*: `Search both strands`
>   - {% icon param-select %} *”Databases to query”*: `Public pre-indexed ribosomal databases`
>   - {% icon param-check %} *”rRNA databases”*: ‘Select all’
>   - {% icon param-select %} *”Include aligned reads in FASTA/FASTQ format?”*: `Yes (--fastx)`
>   - {% icon param-check %} *”Include rejected reads file?”*: `Yes`
>   - {% icon param-check %} *”Generate statistics file”*: `No`
>   - {% icon param-select %} *”Alignment report”*: `Do not report alignments`
>   - {% icon param-text %} *”E-value threshold”*: `1`
>   - {% icon param-text %} *”SW score for a match”*: `2`
>   - {% icon param-text %} *”SW penalty for a mismatch”*: `-3`
>   - {% icon param-text %} *”SW penalty for introducing a gap”*: `5`
>   - {% icon param-text %} *”SW penalty for extending a gap”*: `2`
>   - {% icon param-text %} *”SW penalty for ambiguous letters (N’s)”*: `-3`
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > We used the Galaxy Version 2.1b.6 for our research as it was the most stable and latest version.
>    {: .comment}
>
{: .hands_on}




## **FASTQ interlacer**

> ### {% icon hands_on %} Hands-on: This tool joins paired-end FASTQ reads from two separate files, one with the forward reads and one with the reverse reads, into a single files where forward alternate with their reverse reads. The join is performed using sequence identifiers (headers), allowing the two files to contain differing ordering. If a sequence identifier does not appear in both files, it provides a separate file.

>    > ### {% icon comment %} Comment
>    >
>    > There are two FASTQ interlacers in this workflow. The first one is for MetaPhlAn2 input and the other for HUMANn2 input.
>    {: .comment}
>
{: .hands_on}
>
> 1. **FASTQ interlacer** {% icon tool %} with the following parameters: For MetaPhlAn2
>    - *"Type of paired-end datasets"*: `2 separate datasets`
>
>    ***TODO***: *Check parameter descriptions*
>   - {% icon param-select %} *”Type of paired-end datasets”*: ‘2 separate datasets’
>   - {% icon param-file %} *”Left-hand mates”*: `Cutadapt on data 4 and data 3: Read 1 Output`
>   - {% icon param-file %} *”Right-hand mates”*:`Cutadapt on data 4 and data 3: Read 2 Output`
>
> 2. **Execute**

>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > The FASTQ joiner tool was not used in this workflow as the format of the output was not suitable for our analysis.
>    {: .comment}
>
{: .hands_on}


## **FASTQ interlacer**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **FASTQ interlacer** {% icon tool %} with the following parameters: For HUMANn2
>    - *"Type of paired-end datasets"*: `2 separate datasets`
>
>    ***TODO***: *Check parameter descriptions*
>   - {% icon param-select %} *”Type of paired-end datasets”*: ‘2 separate datasets’
>   - {% icon param-file %} *”Left-hand mates”*: `Filter with sortmeRNA: Unaligned Forward reads`
>   - {% icon param-file %} *”Right-hand mates”*: `Filter with sortmeRNA: Unaligned Reverse reads`
>


## **MetaPhlAn2**

> ### {% icon hands_on %} Hands-on: MetaPhlAn2 is a computational tool used to profile the structure and the composition of microbial communities (Bacteria, Archaea, Eukaryotes and Viruses) from metagenomic shotgun sequencing data with species level resolution. MetaPhlAn2 takes as input a sequence file, in this case the FASTQ interlacer output. The output is a tab-separated file with the predicted taxon relative abundances.
>
> 1. **MetaPhlAn2** {% icon tool %} with the following parameters:
>   - {% icon param-tool %} *”Input file”*: `FASTQ interlacer output`
>   - {% icon param-select %} *”Database with clade-specific marker genes”*: `Locally cached`
>   - {% icon param-select %} *”Cached database with clade-specific marker genes”*: `MetaPhlAn2 clade-specific marker genes`
>   - {% icon param-select %} *”Type of analysis to perform”*: `Profiling a metagenomes in terms of relative abundances`
>   - {% icon param-select %} *”Taxonomic level for the relative abundance output”*: `All taxonomic levels`
>   - {% icon param-text %} *”Minimum total nucleotide length for the markers in a clade for estimating the abundance without considering sub-clade abundances”*: `2000`
>   - {% icon param-text %} *”Sam records for aligned reads with the longest subalignment length smaller than this threshold will be discarded”*: `0`
>   - {% icon param-check %} *”Profile viral organisms?”*: `Yes`
>   - {% icon param-check %} *”Profile eukaryotic organisms?”*: `Yes`
>   - {% icon param-check %} *”Profile bacteria organisms?”*: `Yes`
>   - {% icon param-check %} *”Profile archea organisms?”*: `Yes`
>   - {% icon param-text %} *”Quantile value for the robust average”*: `0.1`
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

## **HUMAnN2**

> ### {% icon hands_on %} Hands-on: HUMAnN2 is a pipeline developed for efficiently and accurately profiling the presence/absence and abundance of microbial pathways in a community from metagenomic or metatranscriptomic sequencing data. For creating the input for HUMAnN2, we interlaced the output from SortmeRNA into a single FASTQ file. This interlaced file and the MetaPhlAn2 output (taxonomic profile) is used as the HUMAnN2 input. For functional profiling, it uses multiple databases that are locally cached nucleotide or protein databases.
>
> 1. **HUMAnN2** {% icon tool %} with the following parameters:
>    - "Input sequence file" to the imported sequence file
>    - "Use of a custom taxonomic profile" to `Yes`
>    - "Taxonomic profile file" to `Community profile` output of `MetaPhlAn2`
>    - "Nucleotide database" to `Locally cached`
>    - "Nucleotide database" to `Full`
>    - "Protein database" to `Locally cached`
>    - "Protein database" to `Full UniRef50`
>    - "Search for uniref50 or uniref90 gene families?" to `uniref50`
>    - "Database to use for pathway computations" to `MetaCyc`
>    - "Advanced Options"
>    - "Remove stratification from output" to `Yes`
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

##  **Export to GraPhlAn**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Export to GraPhlAn** {% icon tool %} with the following parameters:
>   - {% icon param-file %} *”Input file”*: `MetaPhlAn2 Community file`
>   - {% icon param-select %} *”Use a LEfSe output file as input?”*: `No`
>   - {% icon param-text %} *”List which levels should be annotated in the tree”*: `blank`
>   - {% icon param-text %} *”List which levels should use the external legend for the annotation”*: `blank`
>   - {% icon param-text %} *”List which levels should be highlight with a shaded background”*: `blank`
>   - {% icon param-text %} *”List which of the clades that should be highlight with a shaded background”*: `blank`
>   - {% icon param-text %} *”List of color to use for the shaded background”*: `blank`
>   - {% icon param-text %} *”Title of the GraPhlAn plot”*: `blank`
>   - {% icon param-text %} *”Title font size”*: `15`
>   - {% icon param-text %} *”Default size  for clades that are not found as biomarkers”*: `10`
>   - {% icon param-text %} *”Minimum value of clades that are biomarkers”*: `20`
>   - {% icon param-text %} *”Maximum value of clades that are biomarkers”*: `200`
>   - {% icon param-text %} *”Default font size”*: `10`
>   - {% icon param-text %} *”Minimum font size”*: `8`
>   - {% icon param-text %} *”Maximum font size”*: `12`
>   - {% icon param-text %} *”Font size for the annotation legend”*: `10`
>   - {% icon param-text %} *”Minimum abundance value for a clade to be annotated”*: `20.0`
>   - {% icon param-text %} *”Number of clades to highlight”*: `blank`
>   - {% icon param-text %} *”Minimum number of biomarkers to extract”*: `blank`
>   - {% icon param-text %} *”Row number containing the names of the features”*: `0`
>   - {% icon param-text %} *”Row containing the names of the samples”*: `0`
>   - {% icon param-text %} *”Row number to use as metadata”*: `blank`
>   - {% icon param-text %} *”Row number to skip from the input file”*: `blank`
>   - {% icon param-text %} *”Percentile of sample value distribution for sample selection”*: `blank`
>   - {% icon param-text %} *”Percentile of feature value distribution for sample selection”*: ‘’
>   - {% icon param-text %} *”Number of top samples to select”*: `blank`
>   - {% icon param-text %} *”Number of top features to select”*: `blank`
>
>    > ### {% icon comment %} Comment
>    >
>    > It is a conversion tool to produce both annotation and tree file that can serve as an input to GraPhlAn.
>    {: .comment}
>
{: .hands_on}

**Format MetaPhlAn2 output** {% icon tool %} with the following parameters:
>   - {% icon param-file %} *”Input file (MetaPhlAN2 output)”*: ‘MetaPhlAn2 Community File’

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

## **Krona pie chart**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Krona pie chart** {% icon tool %} with the following parameters:
>   - {% icon param-select %} *”What is the type of your input data”*: `Tabular`
>   - {% icon param-file %} *”Input file”*: `Format MetaPhlAn2 Krona output`
>   - {% icon param-text %} *”Provide a name for the basal rank”*: `Root`
>   - {% icon param-check %} *”Combine data from multiple datasets”*: `No`
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

## Sub-step with **Group abundances**

> ### {% icon hands_on %} Hands-on: This tool groups and converts UniRef50 gene family abundances generated with HUMAnN2 into GO slim terms. Custom databases can be provided: a custom Gene Ontology file with the core ontology, a custom slim Gene Ontology with description of slim GO terms to use, and a custom mapping file for mapping UniRef50 gene families to GO terms.

>
> 1. **Group abundances** {% icon tool %} with the following parameters:
>   - {% icon param-file %} *”HUMAnN2 output with UniRef 50 gene family abundance”*: `HUMAnN2 Gene families and their abundance file`
>   - {% icon param-check %} *”Use a custom Gene Ontology file?”*: `No`
>   - {% icon param-check %} *”Use a custom slim Gene Ontology file?”*: `No`
>   - {% icon param-check %} *”Use a custom correspondence between UniReg50 and precise GO?”*: `No`  
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

## Sub-step with **Create a genus level gene families file**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Create a genus level gene families file** {% icon tool %} with the following parameters:
>   - {% icon param-file %} *”Gene families input table”*: `HUMAnN2 Gene families and their abundance file`

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

## **Combine MetaPhlAn2 and HUMAnN2 outputs**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Combine MetaPhlAn2 and HUMAnN2 outputs** {% icon tool %} with the following parameters:
>   - {% icon param-file %} *”Input file corresponding to MetaPhlAN2 output”*: `MetaPhlAn2 Community File`
>   - {% icon param-file %} *”Input file corresponding HUMAnN2 output”*: `HUMAnN2 Gene families and their abundance file`
>   - {% icon param-select %} *”Type of characteristics in HUMAnN2 file”*: `Gene families`
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

## Sub-step with **Combine MetaPhlAn2 and HUMAnN2 outputs**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Combine MetaPhlAn2 and HUMAnN2 outputs** {% icon tool %} with the following parameters:
>   - {% icon param-file %} *”Input file corresponding to MetaPhlAN2 output”*: `MetaPhlAn2 Community File`
>   - {% icon param-file %} *”Input file corresponding HUMAnN2 output”*: `HUMAnN2 Pathways and their abundance file`
>   - {% icon param-select %} *”Type of characteristics in HUMAnN2 file”*: `Gene families`
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

##**Generation, personalization and annotation of tree**

> ### {% icon hands_on %} Hands-on: This tool takes the Export to GraphlAn output ( tree and annotation ) and creates a phylogenetic tree (format: phyloxml).
>
> 1. **Generation, personalization and annotation of tree** {% icon tool %} with the following parameters:
>   - {% icon param-file %} *”Input tree”*: `Export to GraPhlAn Tree output`
>   - {% icon param-file %} *”Annotation file”*: `Export to GraPhlAn Annotation output`

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

## **GraPhlAn**

> ### {% icon hands_on %} Hands-on: GraPhlAn is a software tool for producing high-quality circular representations of taxonomic and phylogenetic trees. It takes the input from the Generation, personalization and annotation of tree for Graphlan tool to create these .png images.
>
> 1. **GraPhlAn** {% icon tool %} with the following parameters:
>   - {% icon param-file %} *”Input tree”*: `Generation, personalization and annotation of tree output`
>   - {% icon param-select %} *”Output format”*: `PNG`
>   - {% icon param-text %} *”Dpi of the output image”*: `blank`
>   - {% icon param-text %} *”Size of the output image”*: `7`
>   - {% icon param-text %} *”Distance between the most external graphical element and the border of the image”*: `blank`
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
