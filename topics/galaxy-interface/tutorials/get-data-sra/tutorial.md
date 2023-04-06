---
layout: tutorial_hands_on

title: "Importing data from NCBI's Sequence Read Archive (SRA) to Galaxy"
zenodo_link: 'http://zenodo.org/record/3906454'
questions:
- Learn how to get data from the Sequence Read Archive in Galaxy.
objectives:
- Understand how Galaxy and the Sequence Read Archive interact.
- Be able to go from Galaxy to the Short Reach Archive, query SRA, use the SRA Run Selector to send selected metadata to Galaxy, and then import sequence data from SRA into Galaxy.
time_estimation: '15m'
key_points:
- Sequence data in the SRA can be directly imported into Galaxy
contributors:
- mvdbeek
- tnabtaf
- blankenberg
- nekrut
- hexylena
tags:
- ncbi
- sra
- get-data
subtopic: upload

---

The aim of this tutorial is to introduce you to the processing of next generation sequencing data in Galaxy.  This tutorial uses a COVID-19 variant calling from Illumina data, but it isn't about variant calling *per se*.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Find necessary data in SRA

First we need to find a good dataset to play with. The [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) is the primary archive of *unassembled reads*  operated by the [US National Institutes of Health (NIH)](https://www.ncbi.nlm.nih.gov/).  SRA is a great place to get the sequencing data that underlie publications and studies. Let's do that:

> <hands-on-title>Task description</hands-on-title>
>
> 1. Go to NCBI's SRA page by pointing your browser to https://www.ncbi.nlm.nih.gov/sra
> 2. In the search box enter `SARS-CoV-2 Patient Sequencing From Partners / MGH`:
>
>    ![The above search query in NCBI SRA's interface.]({% link topics/variant-analysis/images/find_mgh_data.png %}) (Alternatively, you [simply access the data via link](https://www.ncbi.nlm.nih.gov/sra/?term=SARS-CoV-2+Patient+Sequencing+From+Partners+%2F+MGH))
>
> 3. Note that some of the datasets found say "ARTICv3 amplicon sequencing". This is a sequencing technique that requires addition analysis steps not discussed in this tutorial. The data that we will analyse (datasets mentioned below) uses a technique called "metagenomic sequencing".
> 4. The web page will show a large number of SRA datasets (at the time of writing there were 3,927). This is data from a [study](https://science.sciencemag.org/content/early/2020/12/09/science.abe3261) describing analysis of SARS-CoV-2 in Boston area.
> 5. Download metadata describing these datasets by:
>   - clicking on **Send to:** dropdown
>   - Selecting `File`
>   - Changing **Format** to `RunInfo`
>   - Clicking **Create file**
>   Here is how it should look like:
>
>   ![A dropdown in NCBI with aforementioned selections]({% link topics/variant-analysis/images/get_runinfo.png %})
>
> 6. This would create a rather large `SraRunInfo.csv` file in your `Downloads` folder.
{: .hands_on}

Now that we have downloaded this file we can go to a Galaxy instance and start processing it.

> <comment-title></comment-title>
>
> Note that the file we just downloaded is **not** sequencing data itself. Rather, it is *metadata* describing properties of sequencing reads. We will filter this list down to just a few accessions that will be used in the remainder of this tutorial.
>
{: .comment}

## Process and filter `SraRunInfo.csv` file in Galaxy

> <hands-on-title>Upload `SraRunInfo.csv` file into Galaxy</hands-on-title>
>
> 1. Go to your Galaxy instance of choice such as one of the [usegalaxy.org](https://usegalaxy.org/), [usegalaxy.eu](https://usegalaxy.eu), [usegalaxy.org.au](https://usegalaxy.org.au) or any other. (This tutorial uses usegalaxy.org).
> 1. Click *Upload Data* button in Galaxy
> 1. In the dialog box that would appear click "*Choose local files*" button:
> 1. Find and select `SraRunInfo.csv` file from your computer
> 1. Click *Start* button
> 1. Close dialog by pressing **Close** button
> 1. You can now look at the content of this file by clicking {% icon galaxy-eye %} (eye) icon. You will see that this file contains a lot of information about individual SRA accessions. In this study every accession corresponds to an individual patient whose samples were sequenced.
{: .hands_on}

Galaxy can process all 2,000+ datasets but to make this tutorial bearable we need to selected a smaller subset. In particular our previous experience with this data shows two interesting datasets `SRR11954102` and `SRR12733957`. So, let's pull them out.

{% snippet faqs/galaxy/analysis_cut.md %}

> <hands-on-title>Creating a subset of data</hands-on-title>
>
> 1. Find {% icon tool %} "**Select lines that match an expression**" tool in **Filter and Sort** section of the tool panel.
>    > <tip-title>Finding tools</tip-title>
>    > Galaxy may have an overwhelming amount of tools installed. To find a specific tool type the tool name in the tool panel search box to find the tool.
>    {: .tip}
> 1. Make sure the `SraRunInfo.csv` dataset we just uploaded is listed in the {% icon param-file %} "*Select lines from*" field of the tool form.
> 1. In "*the pattern*" field enter the following expression &rarr; `SRR12733957|SRR11954102`. These are two accession we want to find separated by the pipe symbol `|`. The `|` means `or`: find lines containing `SRR12733957` **or** `SRR11954102`.
> 1. Click `Run Tool` button.
> 1. This will generate a file containing two lines (well ... one line is also used as the header, so it will appear the the file has three lines. It is OK.)
> 1. Cut the first column from the file using {% icon tool %} "**Cut**" tool, which you will find in **Text Manipulation** section of the tool pane.
> 1. Make sure the dataset produced by the previous step is selected in the "*File to cut*" field of the tool form.
> 1. Change "*Delimited by*" to `Comma`
> 1. In "*List of fields*" select `Column: 1`.
> 1. Hit `Execute`
>
> This will produce a text file with just two lines:
> ```
> SRR12733957
> SRR11954102
> ```
{: .hands_on}

Now that we have identifiers of datasets we want we need to download the actual sequencing data.

## Download sequencing data with **Faster Download and Extract Reads in FASTQ**

> <hands-on-title>Task description</hands-on-title>
>
> 1. **Faster Download and Extract Reads in FASTQ** {% icon tool %} with the following parameters:
>    - *"select input type"*: `List of SRA accession, one per line`
>        - The parameter {% icon param-file %} *"sra accession list"* should point the output of the {% icon tool %} "**Cut**" from the previous step.
>    - **Click** the `Run Tool` button. This will run the tool, which retrieves the sequence read datasets for the runs that were listed in the `SRA` dataset. It may take some time. So this may be a good time to do get coffee.
>
> 2. Several entries are created in your history panel when you submit this job:
>    - **`Pair-end data (fasterq-dump)`**: Contains Paired-end datasets (if available)
>    - **`Single-end data (fasterq-dump)`** Contains Single-end datasets (if available)
>    - **`Other data (fasterq-dump)`** Contains Unpaired datasets (if available)
>    - **`fasterq-dump log`** Contains Information about the tool execution
{: .hands_on}

The first three items are actually *collections* of datasets.  *Collections* in Galaxy are logical groupings of datasets that reflect the semantic relationships between them in the experiment / analysis.  In this case the tool creates a separate collection each for paired-end reads, single reads, and *other*.
See the Collections tutorials for more.

Explore the collections by first **clicking** on the collection name in the history panel.  This takes you inside the collection and shows you the datasets in it.  You can then navigate back to the outer level of your history.

Once `fasterq` finishes transferring data (all boxes are green / done), you are ready to analyse it!
