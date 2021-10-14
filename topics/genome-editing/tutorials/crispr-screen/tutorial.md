---
layout: tutorial_hands_on
enable: false

title: CRISPR screen analysis
zenodo_link: https://zenodo.org/record/5570011
questions:
- What are the steps to process CRISPR screen data?
- How to identify differentially enriched guides across multiple experimental conditions?
objectives:
- Check quality of raw reads
- Trim sequencing adapters
- Count guide sequences in samples
- Evaluate the quality of count results
- Test for differential enrichment of guides across conditions
time_estimation: 2H
key_points:
- CRISPR screen data can be analysed using MAGeCK and standard read quality tools
contributors:
- mblue9
- kenjifujihara
requirements:
  -
    type: "internal"
    topic_name: introduction
    tutorials:
      - galaxy-intro-short
  -
    type: "internal"
    topic_name: sequence-analysis
    tutorials:
      - quality-control
  -
    type: "internal"
    topic_name: galaxy-interface
    tutorials:
      - collections
      - upload-rules

---


# Introduction
{:.no_toc}

**C**lustered **R**egularly **I**nterspaced **S**hort **P**alindromic **R**epeats (**CRISPR**)-Cas9 is a groundbreaking technology of recent years. It enables editing of the genome and resulted in a Nobel Prize for Emmanuelle Charpentier and Jennifer Doudna in 2020. 

The CRISPR repeat sequences guide the Cas9 enzyme to introduce breaks in DNA. With the CRISPR-Cas9 editing system, a 20 base target region in the genome is added into the CRISPR guide and the Cas9 enzyme then cuts this region.

CRISPR screens provide a high-throughput way to identify genes and pathways that enable cells to survive. In a CRISPR screen all the genes in the genome can be targeted (genome-wide screen) or just a selection (boutique screen). Knockout or activation screens can be performed. 

![Illustration of CRISPR Screen Method](../../images/crispr-screen/crispr-screen.jpg "CRISPR knockout and activation method (from <a href='#Joung2016'> Joung <i>et al.</i> 2016</a>)")

Here we will demonstrate analysing CRISPR screen using data from {% cite Fujihara2020 %}.


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Preparing the reads

## Data upload

We will use fastq files containing 1% of reads from the original samples to demonstrate the read processing steps. 

> ### {% icon hands_on %} Hands-on: Retrieve CRISPR screen fastq datasets
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
>      T0-Control https://zenodo.org/api/files/efc64b27-6db2-4931-ba8c-f05393f520e3/T0-Control.fastq.gz
>      T8-APR-246 https://zenodo.org/api/files/efc64b27-6db2-4931-ba8c-f05393f520e3/T8-APR-246.fastq.gz
>      T8-Vehicle https://zenodo.org/api/files/efc64b27-6db2-4931-ba8c-f05393f520e3/T8-Vehicle.fastq.gz
>      ```
>  
>    ![Rule-based Uploader](../../images/crispr-screen/crispr_rule_uploader.png)
>
>    - From **Rules** menu select `Add / Modify Column Definitions`
>       - Click `Add Definition` button and select `List Identifier(s)`: column `A`
>
>         > ### {% icon tip %} Can't find *List Identifier*?
>         > Then you've chosen to upload as a 'dataset' and not a 'collection'. Close the upload menu, and restart the process, making sure you check *Upload data as*: **Collection(s)**
>         {: .tip}
>
>       - Click `Add Definition` button and select `URL`: column `B`
>
>    - Click `Apply` 
>    - In the Name: box type `fastqs` and press <kbd>Upload</kbd>
>      
>    ![Rule-based Editor](../../images/crispr-screen/crispr_rule_editor.png)
>    
{: .hands_on}

## Raw reads QC

First we'll check the quality of the raw read sequences with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and aggregate the reports from the multiple samples with [MultiQC](https://multiqc.info/) ({% cite ewels2016multiqc %}). We will check if the base quality is good and for presence of adapters. For more details on quality control and what the FastQC plots mean see the ["Quality control" tutorial]({% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %}).


> ### {% icon hands_on %} Hands-on: Quality control
>
> 1. {% tool [FastQC](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.72+galaxy1) %} with the following parameters:
>    - {% icon param-collection %} *"Short read data from your current history"*: input datasets selected with **Multiple datasets**
>
>    {% snippet faqs/galaxy/tools_select_multiple_datasets.md %}
>
> 2. Inspect the webpage output of **FastQC** for the T0 sample
>
>    > ### {% icon question %} Questions
>    >
>    > What is the read length?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > The read length is 75 bp.
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
> 3. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.9+galaxy1) %} with the following parameters to aggregate the FastQC reports:
>     - In *"Results"*
>       - *"Which tool was used generate logs?"*: `FastQC`
>       - In *"FastQC output"*
>         - *"Type of FastQC output?"*: `Raw data`
>         - {% icon param-collection %} *"FastQC output"*: `Raw data` files (output of **FastQC**)
>
> 4. Inspect the webpage output from MultiQC for each FASTQ
>
>    > ### {% icon question %} Questions
>    >
>    > What do you think of the quality of the sequences?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > The quality seems good for the 3 files.
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
{: .hands_on}

We need to trim the adapters to leave just the 20bp guide sequences.

## Trim adapters

We'll trim these sequences using [Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) ({% cite marcel2011cutadapt %}) and its linked adapter format `MY_5PRIME_ADAPTER...MY_3PRIME_ADAPTER`, as discussed [here](https://github.com/marcelm/cutadapt/issues/261#issue-261019127).

> ### {% icon hands_on %} Hands-on: Trim adapters
>
> 1. {% tool [Cutadapt](toolshed.g2.bx.psu.edu/repos/lparsons/cutadapt/cutadapt/3.4+galaxy1) %} with the following parameters:
>    - *"Single-end or Paired-end reads?"*: `Single-end`
>        - {% icon param-collection %} *"FASTQ/A file #1"*: all fastq.gz files 
>        - In *"Read 1 Options"*:
>            - In *"3' (End) Adapters"*:
>                - {% icon param-repeat %} *"Insert 3' (End) Adapters"*
>                    - *"Source"*: `Enter custom sequence`
>                        - *"Enter custom 3' adapter sequence"*: `TGTGGAAAGGACGAAACACCG...GTTTTAGAGCTAGAAATAGCAAG`
>    - In *"Filter Options"*:
>        - *"Minimum length (R1)"*: `20`
>    - In *"Read Modification Options"*:
>        - *"Quality cutoff"*: `20`
>    - *"Outputs selector"*: 
>        - *"Report"*: tick
>
> 2. Inspect the report output from Cutadapt for the T8-APR sample
>
>    > ### {% icon question %} Questions
>    >
>    > What % of reads had adapters?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > 99.9% 
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
{: .hands_on}


# Counting

For the rest of the CRISPR screen analysis, counting and testing, we'll use MAGeCK ({% cite Li2014 %}, {% cite Li2015 %}).

To count how many guides we have for each gene, we need a library file that tells us which guide sequence belongs to which gene. The guides used here are from the Brunello library so we use that file. The file must be tab-separated and contain no spaces within the target names. If necessary, there are tools in Galaxy that can format the file removing spaces and converting commas to tabs.

> ### {% icon hands_on %} Hands-on: Count guides per gene
> 1. Import the library file 
>    ```
>    https://zenodo.org/api/files/efc64b27-6db2-4931-ba8c-f05393f520e3/brunello.tsv
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 2. {% tool [MAGeCK count](toolshed.g2.bx.psu.edu/repos/iuc/mageck_count/mageck_count/0.5.9.2.4) %} with the following parameters:
>    - *"Reads Files or Count Table?"*: `Separate Reads files`
>        - {% icon param-collection %} *"Sample reads"*: the `Read 1 Output` (outputs of **Cutadapt**)
>    - {% icon param-file %} *"sgRNA library file"*: the `brunello.tsv` file
>    - In *"Output Options"*:
>        - *"Output Count Summary file"*: `Yes`
>        - *"Output plots"*: `Yes`
>
> 3. Inspect the Count Summary file
>
>    > ### {% icon question %} Questions
>    >
>    > What percent of reads mapped?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > More than 80% reads mapped in each sample 
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
{: .hands_on}

# Testing

We have been using 1% of reads from the samples in the original dataset to save time as FASTQ files are large. As counts files are small, here we will import and use the MAGeCK counts file generated using all the reads for the samples.

We want to compare the drug treated sample (T8-APR-246) to the control (T8-Vehicle). We could specify them using their names, which must match the names used in the columns of the counts file, but hypens aren't allowed. We can also specify by their positions in the counts file with the first sample column being 0.

> ### {% icon hands_on %} Hands-on: Test for enrichment
> 1. Import the count file from the full dataset [Zenodo]({{ page.zenodo_link }}) or the Shared Data library (if available):
>    ```
>    https://zenodo.org/api/files/efc64b27-6db2-4931-ba8c-f05393f520e3/mageck_counts_full.tsv
>    ```
>
> 2. {% tool [MAGeCKs test](toolshed.g2.bx.psu.edu/repos/iuc/mageck_test/mageck_test/0.5.9.2.1) %} with the following parameters:
>    - {% icon param-file %} *"Counts file"*: the `mageck_counts_full.tsv` file
>    - *"Specify Treated samples or Control"*: `Treated samples`
>        - *"Treated Sample Labels (or Indexes)"*: `0`
>    - *"Control Sample Labels (or Indexes)"*: `1`
>    - In *"Output Options"*:
>        - *"Output normalized counts file"*: `Yes`
>        - *"Output plots"*: `Yes`
>    - In *"Advanced Options"*:
>        - *"Method for normalization"*: `Total`
>
> 3. Inspect the PDF Report output.
>
>    > ### {% icon question %} Questions
>    >
>    > What are the top 3 genes?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > ESD, MTHFD1L and SHMT2 which are part of the glutathione pathway. That was the main pathway found to be altered in the published paper for this dataset.
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
{: .hands_on}


# Conclusion
{:.no_toc}

CRISPR Screen reads can be assessed for quality using standard sequencing tools such as FASTQC, MultiQC and trimmed of adapters using Cutadapt. The detection of enriched guides can be performed using MAGeCK.
