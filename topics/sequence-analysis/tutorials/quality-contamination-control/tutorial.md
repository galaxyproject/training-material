---
layout: tutorial_hands_on
title: Quality and contamination control in bacterial isolate using Illumina MiSeq
  Data
zenodo_link: https://zenodo.org/record/10669812
questions:
- How to check the quality of MiSeq data?
- What are the species in bacterial isolate sequencing data?
objectives:
- Run tools to evaluate sequencing data on quality and quantity
- Evaluate the output of quality control tools
- Improve the quality of sequencing data
- Run a series of tool to identify species in bacterial isolate sequencing data
- Visualize the species abundance
time_estimation: 2H
key_points:
- Conduct quality control on every dataset before performing any other bioinformatics
  analysis
- Review the quality metrics and, if necessary, improve the quality of your data
- Check the impact of the quality control
- Detect witch microorganisms are present and extract the species level
- Check possible contamination in a bacterial isolate
tags:
- illumina
- bacteria
- microgalaxy
level: Introductory
edam_ontology:
- topic_3673
- topic_0622
- topic_3301
- topic_3697
contributions:
  authorship:
  - bebatut
  - clsiguret
  funding:
  - abromics
recordings:
- youtube_id: Cx31r5emUJk
  length: 34M
  galaxy_version: 24.0.4.dev0
  date: '2024-08-01'
  speakers: [clsiguret]
  captioners: [clsiguret]
  bot-timestamp: 1722527209


---


Sequencing (determining of DNA/RNA nucleotide sequence) is used all over the world for all kinds of analysis. The product of these sequencers are reads, which are sequences of detected nucleotides. Depending on the technique these have specific lengths (30-500bp) or using Oxford Nanopore Technologies sequencing have much longer variable lengths.

{% snippet faqs/galaxy/sequencing_illumina_miseq.md %}

Contemporary sequencing technologies are capable of generating an enormous volume of sequence reads in a single run. Nonetheless, each technology has its imperfections, leading to various types and frequencies of errors, such as miscalled nucleotides. These inaccuracies stem from the inherent technical constraints of each sequencing platform. When sequencing bacterial isolates, it is crucial to verify the quality of the data but also to check the expected species or strains are found in the data or if there is any contamination.

To illustrate the process, we take raw data of a bacterial genome (KUN1163 sample) from sequencing data produced in "Complete Genome Sequences of Eight Methicillin-Resistant *Staphylococcus aureus* Strains Isolated from Patients in Japan" ({% cite Hikichi_2019 %}). 

> Methicillin-resistant *Staphylococcus aureus* (MRSA) is a major pathogen
> causing nosocomial infections, and the clinical manifestations of MRSA
> range from asymptomatic colonization of the nasal mucosa to soft tissue
> infection to fulminant invasive disease. Here, we report the complete
> genome sequences of eight MRSA strains isolated from patients in Japan.
{: .quote cite="{% cite_url Hikichi_2019 %}"}

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Galaxy and data preparation

Any analysis should get its own Galaxy history. So let's start by creating a new one and get the data (forward and reverse quality-controlled sequences) into it.

> <hands-on-title>Prepare Galaxy and data</hands-on-title>
>
> 1. Create a new history for this analysis
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Rename the history
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
{: .hands_on}

Now, we need to import the data: 2 FASTQ files containing the reads from the sequencer.

> <hands-on-title>Import datasets</hands-on-title>
> 1. {% tool [Import](upload1) %} the files from [Zenodo]({{ page.zenodo_link }}) or from Galaxy shared data libraries:
>
>    ```
>    {{ page.zenodo_link }}/files/DRR187559_1.fastqsanger.bz2
>    {{ page.zenodo_link }}/files/DRR187559_2.fastqsanger.bz2
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
>
> 2. Rename the datasets to remove `.fastqsanger.bz2` and keep only the sequence run ID (`DRR187559_1` and `DRR187559_2`)
>
>    {% snippet faqs/galaxy/datasets_rename.md name="DRR187559_1" %}
>
> 3. Tag both datasets `#unfiltered`
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
> 4. **View** {% icon galaxy-eye %} the renamed file
>
{: .hands_on}

The datasets are both FASTQ files.

> <question-title></question-title>
>
> 1. What are the 4 main features of each read in a FASTQ file.
> 2. What does the `_1` and `_2` mean in your filenames?
>
> > <solution-title></solution-title>
> > 1. The following:
> >
> >    -   A `@` followed by a name and sometimes information of the read
> >    -   A nucleotide sequence
> >    -   A `+` (optional followed by the name)
> >    -   The quality score per base of nucleotide sequence (Each symbol
> >        represents a quality score, which will be explained later)
> >
> > 2. Forward and reverse reads, by convention, are labelled `_1` and `_2`, but they might also be `_f`/`_r` or `_r1`/`_r2`.
> {: .solution}
{: .question}

{% include _includes/cyoa-choices.html option1="Step-by-step" option2="Workflow" default="Step-by-step" text="Do you want to go step-by-step or using a workflow?" disambiguation="workflow"%}

<div class="Workflow" markdown="1">

In this section we will run a Galaxy workflow that performs the following tasks with the following tools:
1. Assess the reads quality before preprocessing it using [__FastQC__](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
2. Trimming and filtering reads by length and quality using **Fastp** ({% cite Chen2018 %}).
3. Find witch microorgasnims are present using [__Kraken2__](https://ccb.jhu.edu/software/kraken2/) ({% cite Wood2014 %}).
4. Extract the species level with **Bracken** (Bayesian Reestimation of Abundance after Classification with Kraken) ({% cite Lu.2017 %}).
5. Detect minority organisms or contamination using **Recentrifuge** ({% cite marti2019recentrifuge %}).

We will run all these steps using a single workflow, then discuss each step and the results in more detail.

> <hands-on-title>Pre-Processing</hands-on-title>
>
> 1. **Import the workflow** into Galaxy
>    - Copy the URL (e.g. via right-click) of [this workflow]({{ site.baseurl }}{{ page.dir }}workflows/Quality_and_contamination_control_in_bacterial_isolate_using_Illumina_MiSeq_Data.ga) or download it to your computer.
>    - Import the workflow into Galaxy
>
>    {% snippet faqs/galaxy/workflows_import.md %}
>
> 2. Run **Workflow : Quality and contamination control in bacterial isolate using Illumina MiSeq Data** {% icon workflow %} using the following parameters
>    - *"DRR187559_1"*: `DRR187559_1`, which is the forward read data.
>
>    - *"DRR187559_2"*: `DRR187559_2`, which is the reverse read data.
>
>    {% snippet faqs/galaxy/workflows_run.md %}
>
{: .hands_on}

The workflow will take some time. Once completed, results will be available in your history. While waiting, read the next sections for details on each workflow step and the corresponding outputs.
</div>


# Read quality control and improvement

During sequencing, errors are introduced, such as incorrect nucleotides being called. These are due to the technical limitations of each sequencing platform. Sequencing errors might bias the analysis and can lead to a misinterpretation of the data. Adapters may also be present if the reads are longer than the fragments sequenced and trimming these may improve the number of reads mapped. **Sequence quality control is therefore an essential first step in any analysis.**

Assessing the quality by hand would be too much work. That's why tools like
[NanoPlot](https://github.com/wdecoster/NanoPlot) or
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) are made, as they  generate a summary and plots of the data statistics. NanoPlot is
mainly used for long-read data, like ONT and PACBIO and FastQC for short read,
like Illumina and Sanger. You can read more in our dedicated [Quality Control
Tutorial]({% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %}).

Before doing any analysis, the first questions we should ask about the input
reads include:

- What is the coverage of my genome?
- How good are my reads?
- Do I need to ask/perform for a new sequencing run?
- Is it suitable for the analysis I need to do?

## Quality control

<div class="Step-by-step" markdown="1">
> <hands-on-title>Quality Control</hands-on-title>
>
> 1. {% tool [FastQC](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.74+galaxy0) %} with the following parameters:
>    - {% icon param-files %} *"Short read data from your current history"*: both `DRR187559_1` and `DRR187559_2`
>
>    {% snippet faqs/galaxy/tools_select_multiple_datasets.md %} 
>
> 2. Inspect the webpage outputs
>
{: .hands_on}
</div>

<div class="Workflow" markdown="1">

> <hands-on-title>Quality Control</hands-on-title>
> Inspect the webpage outputs of **FastQC**
>
{: .hands_on}

</div>

**FastQC** combines quality statistics from all separate reads and combines them in plots. An important plot is the Per base sequence quality. 

DRR187559_1 | DRR187559_2
----------- | -----------
![FastQC plot showing reads that mostly stay in the green](../../../assembly/images/mrsa/fastqc-1.png) | ![Same as previous plot, but the beginning of the reads are slightly better quality](../../../assembly/images/mrsa/fastqc-2.png)

Here you have the reads sequence length on the x-axes against the quality score (Phred-score) on the y-axis. The y-axis is divided in three sections: 
- Green = good quality, 
- Orange = mediocre quality, and 
- Red = bad quality.

For each position, a boxplot is drawn with:

- the median value, represented by the central red line
- the inter-quartile range (25-75%), represented by the yellow box
- the 10% and 90% values in the upper and lower whiskers
- the mean quality, represented by the blue line

> <question-title></question-title>
>
> How does the mean quality score change along the sequence?
>
> > <solution-title></solution-title>
> > The mean quality score (blue line) decreases at the sequences end. It is common for the mean quality to drop towards the end of the sequences, as the sequencers are incorporating more incorrect nucleotides at the end. For Illumina data it is normal that the first few bases are of some lower quality and how longer the reads get the worse the quality becomes. This is often due to signal decay or phasing during the sequencing run.
> >
> {: .solution }
{: .question}

## Quality improvement

Depending on the analysis it could be possible that a certain quality or length
is needed. In this case we are going to trim the data using **fastp** ({% cite Chen2018 %}):

- Trim the start and end of the reads if those fall below a quality score of 20

  Different trimming tools have different algorithms for deciding when to cut but trimmomatic will cut based on the quality score of one base alone. Trimmomatic starts from each end, and as long as the base is below 20, it will be cut until it reaches one greater than 20. A sliding window trimming will be performed where if the average quality of 4 bases drops below 20, the read will be truncated there. 
  
- Filter for reads to keep only reads with at least 30 bases: Anything shorter will complicate the assembly


<div class="Step-by-step" markdown="1">
> <hands-on-title>Quality improvement</hands-on-title>
>
> 1. {% tool [fastp](toolshed.g2.bx.psu.edu/repos/iuc/fastp/fastp/0.23.4+galaxy0) %} with the following parameters:
>    - *"Single-end or paired reads"*: `Paired`
>        - {% icon param-file %} *"Input 1"*: `DRR187559_1`
>        - {% icon param-file %} *"Input 2"*: `DRR187559_2`
>    - In *"Filter Options"*:
>        - In *"Length filtering Options"*:
>            - *Length required*: `30`
>    - In *"Read Modification Options"*:
>        - In *"Per read cuitting by quality options"*:
>            - *Cut by quality in front (5')*: `Yes`
>            - *Cut by quality in front (3')*: `Yes`
>            - *Cutting window size*: `4`
>            - *Cutting mean quality*: `20`
>    - In *"Output Options"*:
>        - *"Output JSON report"*: `Yes`
>
> 2. Edit the tags of the **fastp** FASTQ outputs to 
>    1. Remove the `#unfiltered` tag
>    2. Add a new tag `#filtered`
{: .hands_on}
</div>

**fastp** generates also a report, similar to **FastQC**, useful to compare the impact of the trimming and filtering.

> <question-title></question-title>
>
> Looking at **fastp** HTML report
>
> 1. How did the average read length change before and after filtering?
> 3. Did trimming improve the mean quality scores?
> 4. Did trimming affect the GC content?
> 5. Is this data ok to assemble? Do we need to re-sequence it?
>
> > <solution-title></solution-title>
> >
> > 1. Read lengths went down more significantly:
> >    - Before filtering: 190bp, 221bp
> >    - After filtering: 189bp, 219bp
> > 3. It increased the percentage of Q20 and Q30 bases (bases with quality score above 20 and 30 respectively)
> > 4. No, it did not. If it had, that would be unexpected.
> > 5. This data looks OK. The number of short reads in R1 is not optimal but assembly should partially work but not the entire, closed genome.
> >
> {: .solution}
{: .question}

# Identification of expected species and detection of contamination

When working with bacterial isolates, it is crucial to verify whether the expected species or strains are present in the data and to identify any potential contamination. Ensuring the presence of the intended species is essential for the accuracy and reliability of the research, as deviations could lead to erroneous conclusions. Additionally, detecting contamination is vital to maintain the integrity of the samples and to avoid misleading results that could compromise subsequent analyses and applications.

## Taxonomic profiling

To find out which microorganisms are present, we will compare the filtered reads of the sample to a reference database, i.e. sequences of known microorganisms stored in a database, using **Kraken2** ({% cite wood2019improved %}).

{% snippet topics/microbiome/faqs/kraken.md %}

For this tutorial, we will use the PlusPF database which contains the RefSeq Standard (archaea, bacteria, viral, plasmid, human, UniVec_Core), protozoa and fungi data.

<div class="Step-by-step" markdown="1">
> <hands-on-title> Assign taxonomic labels with Kraken2</hands-on-title>
>
> 1. {% tool [Kraken2](toolshed.g2.bx.psu.edu/repos/iuc/kraken2/kraken2/2.1.3+galaxy1) %} with the following parameters:
>    - *"Single or paired reads"*: `Paired`
>        - {% icon param-file %} *"Forward strand"*: **fastp** `Read 1 output`
>        - {% icon param-file %} *"Reverse strand"*: **fastp** `Read 2 output`
>    - *"Minimum Base Quality"*: `10`
>    - In *"Create Report"*:
>        - *"Print a report with aggregrate counts/clade to file"*: `Yes`
>    - *"Select a Kraken2 database"*: `PlusPF-16`
>
{: .hands_on}
</div>

**Kraken2** generates 2 outputs:

- **Classification**: tabular files with one line for each sequence classified by Kraken and 5 columns:

   1. `C`/`U`: a one letter indicating if the sequence classified or unclassified
   2. Sequence ID as in the input file
   3. NCBI taxonomy ID assigned to the sequence, or 0 if unclassified
   4. Length of sequence in bp (`read1|read2` for paired reads)
   5. A space-delimited list indicating the lowest common ancestor (LCA) mapping of each k-mer in the sequence     
      
      For example, `562:13 561:4 A:31 0:1 562:3` would indicate that:
      1. The first 13 k-mers mapped to taxonomy ID #562
      2. The next 4 k-mers mapped to taxonomy ID #561
      3. The next 31 k-mers contained an ambiguous nucleotide
      4. The next k-mer was not in the database
      5. The last 3 k-mers mapped to taxonomy ID #562     
      
      `|:|` indicates end of first read, start of second read for paired reads

   ```
   Column 1	Column 2	Column 3	Column 4	Column 5
   C	DRR187559.1	1280	164|85	0:1 1279:1 0:41 1279:10 0:5 1280:5 1279:1 0:1 1279:1 0:7 1279:5 0:2 1279:6 0:12 1279:5 0:19 1279:2 0:6 |:| 0:39 1279:2 0:6 1279:3 0:1
   C	DRR187559.2	1280	70|198	0:2 1279:5 0:29 |:| 0:52 1279:5 0:13 1279:3 0:23 1279:2 0:45 1280:1 0:9 1280:3 A:8
   C	DRR187559.3	1279	106|73	0:4 1279:3 0:36 1279:4 0:10 1279:5 0:3 1279:5 0:2 |:| 0:39
   C	DRR187559.4	1279	121|189	1279:6 0:17 1279:4 0:28 1279:2 0:30 |:| 0:7 1279:5 0:19 1279:5 0:20 1279:1 0:8 1279:1 0:1 1279:6 0:25 1279:2 0:44 A:11
   C	DRR187559.5	1279	68|150	1279:2 0:20 1279:3 0:9 |:| 0:10 1279:3 0:24 1279:2 0:9 1279:5 0:21 1279:5 0:20 1279:5 0:9 1279:1 0:2
   C	DRR187559.6	1280	137|246	0:2 1280:5 0:28 1279:1 0:28 1279:2 0:8 1279:2 0:23 1279:1 0:3 |:| 1279:1 0:2 1279:3 0:61 1279:2 0:14 1279:2 0:97 A:30 
   ```

   > <question-title></question-title>
   >
   > 1. Is the first sequence in the file classified or unclassified?
   > 2. What is the taxonomy ID assigned to the first classified sequence?
   > 3. What is the corresponding taxon?
   >
   > > <solution-title></solution-title>
   > > 1. classified
   > > 2. 1280
   > > 3. 1280 corresponds to [Staphylococcus aureus](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/tree/?taxon=1280) .
   > {: .solution}
   {: .question}

- **Report**: tabular files with one line per taxon and 6 columns or fields

   1. Percentage of fragments covered by the clade rooted at this taxon
   2. Number of fragments covered by the clade rooted at this taxon
   3. Number of fragments assigned directly to this taxon
   4. A rank code, indicating
      - (U)nclassified
      - (R)oot
      - (D)omain
      - (K)ingdom
      - (P)hylum
      - (C)lass
      - (O)rder
      - (F)amily
      - (G)enus, or
      - (S)pecies

      Taxa that are not at any of these 10 ranks have a rank code that is formed by using the rank code of the closest ancestor rank with a number indicating the distance from that rank. E.g., `G2` is a rank code indicating a taxon is between genus and species and the grandparent taxon is at the genus rank.

   5. NCBI taxonomic ID number
   6. Indented scientific name

   ```
   Column 1	Column 2	Column 3	Column 4	Column 5	Column 6
   0.24 	1065 	1065 	U 	0 	unclassified
   99.76 	450716 	14873 	R 	1 	root
   96.44 	435695 	2 	R1 	131567 	cellular organisms
   96.44 	435675 	3889 	D 	2 	Bacteria
   95.56 	431709 	78 	D1 	1783272 	Terrabacteria group
   95.53 	431578 	163 	P 	1239 	Firmicutes
   95.49 	431390 	4625 	C 	91061 	Bacilli
   94.38 	426383 	1436 	O 	1385 	Bacillales
   94.04 	424874 	2689 	F 	90964 	Staphylococcaceae
   93.44 	422124 	234829 	G 	1279 	Staphylococcus 
   ```

   > <question-title></question-title>
   >
   > 1. How many taxa have been found?
   > 2. What are the percentage on unclassified?
   > 3. What are the domains found?
   >
   > > <solution-title></solution-title>
   > >
   > > 1. 627, as the number of lines
   > > 2. 0.24%
   > > 3. Only Bacteria
   > {: .solution}
   >
   {: .question}

## Species identification

In Kraken output, there are quite a lot of identified taxa with different levels. To obtain a more accurate and detailed understanding at the species level, we will use __Bracken__. 

__Bracken__ refines the Kraken results by re-estimating the abundances of species in metagenomic samples, providing a more precise and reliable identification of species, which is crucial for downstream analysis and interpretation.

__Bracken__ (Bayesian Reestimation of Abundance after Classification with Kraken) is a "simple and worthwile addition to Kraken for better abundance estimates" ({% cite Ye.2019 %}). Instead of only using proportions of classified reads, it takes a probabilistic approach to generate final abundance profiles. It works by re-distributing reads in the taxonomic tree: "Reads assigned to nodes above the species level are distributed down to the species nodes, while reads assigned at the strain level are re-distributed upward to their parent species" ({% cite Lu.2017 %}).

<div class="Step-by-step" markdown="1">
> <hands-on-title>Extract species with Bracken</hands-on-title>
>
> 1. {% tool [Bracken](toolshed.g2.bx.psu.edu/repos/iuc/bracken/est_abundance/2.9+galaxy0) %} with the following parameters:
>     - {% icon param-collection %} *"Kraken report file"*: **Report** output of **Kraken**
>     - *"Select a kmer distribution"*: `PlusPF-16`, same as for Kraken
>
>        It is important to choose the same database that you also chose for Kraken2
>
>     - *"Level"*: `Species`
>     - *"Produce Kraken-Style Bracken report"*: `Yes`
>
{: .hands_on}
</div>

**Bracken** generates 2 outputs:

- **Kraken style report**: tabular files with one line per taxon and 6 columns or fields. Same configuration as the **Report** output of **Kraken**:

   1. Percentage of fragments covered by the clade rooted at this taxon
   2. Number of fragments covered by the clade rooted at this taxon
   3. Number of fragments assigned directly to this taxon
   4. A rank code, indicating
      - (U)nclassified
      - (R)oot
      - (D)omain
      - (K)ingdom
      - (P)hylum
      - (C)lass
      - (O)rder
      - (F)amily
      - (G)enus, or
      - (S)pecies
   5. NCBI taxonomic ID number
   6. Indented scientific name

   ```
   Column 1	Column 2	Column 3	Column 4	Column 5	Column 6
   100.00 	450408 	0 	R 	1 	root
   99.14 	446538 	0 	R1 	131567 	cellular organisms
   99.14 	446524 	0 	D 	2 	Bacteria
   99.14 	446524 	0 	D1 	1783272 	Terrabacteria group
   99.13 	446491 	0 	P 	1239 	Firmicutes
   99.13 	446491 	0 	C 	91061 	Bacilli
   99.06 	446152 	0 	O 	1385 	Bacillales
   99.06 	446152 	0 	F 	90964 	Staphylococcaceae
   99.04 	446101 	0 	G 	1279 	Staphylococcus
   95.62 	430661 	430661 	S 	1280 	Staphylococcus aureus 
   ```

   > <question-title></question-title>
   >
   > 1. How many taxa have been found?
   > 2. What is the family found?
   >
   > > <solution-title></solution-title>
   > >
   > > 1. 119, as the number of lines
   > > 2. Staphylococcaceae, with 99.06%!
   > {: .solution}
   >
   {: .question}

- **Report**: tabular files with one line per taxon and 7 columns or fields

   1. Taxon name
   2. Taxonomy ID
   3. Level ID (S=Species, G=Genus, O=Order, F=Family, P=Phylum, K=Kingdom)
   4. Kraken assigned reads
   5. Added reads with abundance re-estimation
   6. Total reads after abundance re-estimation
   7. Fraction of total reads

   ```
   Column 1	Column 2	Column 3	Column 4	Column 5	Column 6	Column 7
   Staphylococcus aureus 	1280 	S 	182995 	247666 	430661 	0.95616
   Staphylococcus roterodami 	2699836 	S 	698 	48 	746 	0.00166
   Staphylococcus epidermidis 	1282 	S 	587 	13 	600 	0.00133
   Staphylococcus lugdunensis 	28035 	S 	511 	3 	514 	0.00114
   Staphylococcus schweitzeri 	1654388 	S 	409 	26 	435 	0.00097
   Staphylococcus simiae 	308354 	S 	398 	2 	400 	0.00089 
   ```

> <question-title></question-title>
>
> 1. How many species have been found?
> 2. Which the species has been the most identified? And in which proportion?
> 3. What are the other species?
>
> > <solution-title></solution-title>
> > 1. 51 (52 lines including 1 line with header)
> > 2. *Staphylococcus aureus* with 95.6% of the reads.
> > 3. Most of the other species are from *Staphylococcus* genus, so same as *Staphylococcus aureus*. The other species in really low proportion.
> {: .solution}
{: .question}

As expected *Staphylococcus aureus* represents most of the reads in the data.

## Contamination detection

To explore **Kraken** report and specially to detect more reliably minority organisms or contamination, we will use **Recentrifuge** ({% cite marti2019recentrifuge %}). 

**Recentrifuge** enhances analysis by reanalyzing metagenomic classifications with interactive charts that highlight confidence levels. It supports 48 taxonomic ranks of the [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy), including strains, and generates plots for shared and exclusive taxa, facilitating robust comparative analysis. 

**Recentrifuge** includes a novel contamination removal algorithm, useful when negative control samples are available, ensuring data integrity with control-subtracted plots. It also excels in detecting minority organisms in complex datasets, crucial for sensitive applications such as clinical and environmental studies.

<div class="Step-by-step" markdown="1">
> <hands-on-title> Identify contamination </hands-on-title>
>
> 1. {% tool [Recentrifuge](toolshed.g2.bx.psu.edu/repos/iuc/recentrifuge/recentrifuge/1.12.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select taxonomy file tabular formated"*: **Classification** output of **Krancken2** {% icon tool %}
>    - *"Type of input file (Centrifuge, CLARK, Generic, Kraken, LMAT)"*: `Kraken`
>    - In *"Database type"*:
>        - *"Cached database whith taxa ID"*: `NCBI-2023-06-27`
>    - In *"Output options"*:
>        - *"Type of extra output to be generated (default on CSV)"*: `TSV`
>        - *"Remove the log file"*: `Yes`
>    - In *" Fine tuning of algorithm parameters"*:
>        - *"Strain level instead of species as the resolution limit for the robust contamination removal algorithm; use with caution, this is an experimental feature"*: `Yes`
>
{: .hands_on}
</div>

**Recentrifuge** generates 3 outputs:

- A **statistic table** with general statistics about assignations

   > <question-title></question-title>
   >
   > 1. How many sequences have been used?
   > 2. How many sequences have been classified?
   >
   > > <solution-title></solution-title>
   > > 1. 451,780
   > > 2. 450,715
   > {: .solution}
   {: .question}

- A **data table** with a report for each taxa

   > <question-title></question-title>
   >
   > 1. How many taxa have been kept?
   > 2. What is the lowest level in the data?
   >
   > > <solution-title></solution-title>
   > > 1. 187 (190 lines including 3 header lines)
   > > 2. The lowest level is strain.
   > {: .solution}
   {: .question}

- A **HTML report** with a Krona chart

   <iframe id="recentrifuge" src="recentrifuge_report.html" frameBorder="0" width="100%" height="900px"> ![Krona chart with multi-layered pie chart representing the community profile with in the center the higher taxonomy levels (i.e. domain) and on the exterior the more detailed ones (i.e. species)](./images/recentrifuge.png) </iframe>

   > <question-title></question-title>
   >
   > 1. What is the percentage of classified sequences?
   > 2. When clicking on *Staphylococcus aureus*, what can we say about the strains?
   > 3. Is there any contamination?
   >
   > > <solution-title></solution-title>
   > >
   > > 1. 99.8%
   > > 2. 99% of sequences assigned to *Staphylococcus aureus* are not assigned to any strains, probably because they are too similar to several strains. *Staphylococcus aureus* subsp. aureus JKD6159 is the strain with the most classified sequences, but only 0.3% of the sequences assigned to *Staphylococcus aureus*.
   > > 3. There is no sign of a possible contamination. Most sequences are classified to taxon on the *Staphylococcus aureus* taxonomy. Only 3% of the sequences are not classified to *Staphylococcus*.
   > >
   > {: .solution}
   >
   {: .question}

Once we have identified contamination, if any is present, the next step is to **remove the contaminated sequences** from the dataset to ensure the integrity of the remaining data. This can be done using bioinformatics tools designed to filter out unwanted sequences, such as **BBduk** ({% cite Bushnell2017 %}).

**BBduk** {% icon tool %} is a member of the [BBTools](https://sourceforge.net/projects/bbmap/) package, where 'Duk' stands for *Decontamination Using Kmers*. **BBduk** filters or trims reads for adapters and contaminants using k-mers, effectively removing unwanted sequences and improving the quality of the dataset.

Additionally, it's important to document and **report the contamination findings** to maintain transparency and guide any necessary adjustments in sample collection or processing protocols.

# Conclusion

In this tutorial, we inspected the quality of the bacterial isolate sequencing data and checked the expected species and potential contamination. Prepared short reads **can be used in downstream analysis**, like [Genome Assembly]({% link topics/assembly/tutorials/mrsa-illumina/tutorial.md %}). Even after the assembly, you can identify contamination using a tool like **CheckM** ({% cite Parks_2015 %}). **CheckM** {% icon tool %} can be used for screening the 'cleanness' of bacterial assemblies.

To learn more about data quality control, you can follow this tutorial: [Quality Control]({% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %}).
