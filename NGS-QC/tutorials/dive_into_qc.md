---
layout: tutorial_hands_on
topic_name: NGS-QC
tutorial_name: dive_into_qc
---

# Introduction

During the sequencing process, some errors may be introduced like incorporation of ambiguous nucleotides. Analyzing poor data wastes CPU and people time. 

The quality control of the sequences right after sequencing is then an essential step to ensure that the raw data looks good before processing them. It reduces the biases in data that may compromised the downstream analyses with low-quality sequences, sequence artifacts, ... The process is the same for any type of sequencing data.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. [Sequence dataset importing](#sequence-dataset-importing)
> 2. [Quality checking of the sequences](#quality-checking-of-the-sequences)
> 3. [Improvement of the quality of the sequences](#improvement-of-the-quality-of-the-sequences)
> {: .agenda}

# Sequence dataset importing

> ### :pencil2: Hands-on: Data upload
>
> 1. Create a new history
> 2. Import the FASTQ file: [`GSM461178_untreat_paired_subset_1`](https://zenodo.org/record/61771/files/GSM461178_untreat_paired_subset_1.fastq)
>
>    > ### :bulb: Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * Press **Start**    
>    {: .tip}
>
>    > ### :bulb: Tip: Changing the file type `fastq` to `fastqsanger` once the data file is in your history
>    >
>    > * Click on the pencil button displayed in your dataset in the history
>    > * Choose **Datatype** on the top
>    > * Select `fastqsanger`
>    > * Press **save**
>    {: .tip}
>
>    > ### :nut_and_bolt: Comments
>    > 
>    > Rename the dataset to "First dataset"
>    {: .comment}
> As default, Galaxy takes the link as name.
{: .hands_on}

# Quality checking of the sequences

To estimate sequence quality and treatments to do on the data, many indicators can be checked:

- Quality score of the sequences with
    - Per base sequence quality
    - Per sequence quality scores
    - Per tile sequence quality
- Sequence content with 
    - Per base sequence content
    - Per sequence GC content
    - Per base N content
- Sequence length with the sequence length distribution
- Duplicated sequences
- Tag sequences with 
    - Adapter contamination
    - Kmer Content

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is an open-source tool provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It generates quality graphics and estimates numerous quality informations and threshold. For each studied indicators, FastQC providing a quick overview to tell in which areas there may be problems.

> ### :pencil2: Hands-on: Run FastQC
>
> 1. **FastQC** :wrench:: Run FastQC on the imported FastQ file with default parameters
> 2. Inspect the FastQC report on the webpage
>
>    > ### :bulb: Tip: Inspecting the content of a file in Galaxy
>    >
>    > * Click on the eye ("View data") on the right of the file name in the history
>    > * Inspect the content of the file on the middle 
>    {: .tip}
>
>    > ### :question: Questions
>    >
>    > 1. How good are the quality scores?
>    > 2. Why is there warning for the per base sequence content and the per sequence GC content graphs?
>    > 3. What must be done to improve the sequences?
>    > 
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>The scores of the sequences are quite good: no warning from FastQC, even if we can see a slight decrease of the quality at the end of sequences </li>
>    >    <li>In the beginning of sequences, the sequence content per base is not really good and the percentage are not equal. For the GC content, the distribution is bit shifted on the left and too high</li>
>    >    <li>We can trim a bit the end of the sequences, but not too much as the sequences are already small</li>
>    >    </ol>
>    >    </details>
>    {: .question}
{: .hands_on}

# Improvement of the quality of the sequences

Based on previous quality graphs, sequences must be treated to obtain good dataset and then the bias in downstream analysis.

In general, quality treatments are:

- Filtering of sequences
    - with small mean quality score
    - too small
    - with too many N bases
    - based on their GC content
    - ...
- Cutting/Trimming sequences
    - from low quality score parts
    - tails
    - ...

To improve the quality of the sequences, we use [Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) tool. It automates quality and adapter trimming as well as quality control. 

> ### :pencil2: Hands-on
>
> 1. **Trim Galore** :wrench:: Run Trim Galore on the imported data file
>
>    > ### :question: Questions
>    >
>    > Which parameters must be applied to follow the previous recommendations?
>    > 
>    > <details>
>    > <summary>Click to view answers</summary>
>    > We use the default ones:
>    > <ul>
>    > <li>Automatic detection and trimming of adapter sequences</li>
>    > <li>Trimming low-quality ends (below 20) from reads in addition to adapter removal</li>
>    > <li>Removing reads shorter than 20 bp</li>
>    > </ul>
>    > </details>
>    {: .question}
>
> 2. **FastQC** :wrench:: Re-run FastQC on the quality controlled data file and inspect the new FastQC report
>
>    > ### :question: Questions
>    > 
>    > 1. How many sequences have been removed?
>    > 2. Has the quality of the sequences been improved?
>    > 3. Can you explain why the per base sequence content is not good now?
>    > 
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>Before Trim Galore, the dataset was made of 100,000 sequences. After Trim Galore, there is 99,653 sequences.</li>
>    >    <li>The per base quality score is better but other indicators are bad now. The sequence length distribution is not clear as before because sequences have different size after the trimming</li>
>    >    <li>The per base sequence content is red now. The trimming of the end of some sequences may have biased </li>
>    >    </ol>
>    >    </details>
>    {: .question}
{: .hands_on}

The quality of the previous dataset was pretty good from beginning. The quality treatment improved the quality score but to the cost of other parameters.

# Control the quality of a second dataset

Now, we would like to see the impact to quality control and treatment on a bad dataset.

> ### :pencil2: Hands-on
>
> 1. Create a new history
> 2. Import the FASTQ file: [`GSM461182_untreat_single_subset`](https://zenodo.org/record/61771/files/GSM461182_untreat_single_subset.fastq)
> 3. **FastQC** :wrench:: Run FastQC on this new dataset
>
>    > ### :question: Questions
>    >
>    > 1. How good is this dataset?
>    > 2. What must be done to improve the sequences?
>    > 
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>There is red warning for the per base sequence quality (pretty bad along the sequence but worst at the end of sequences), the per base sequence content (bad at the beginning of the sequences), the per sequence GC content. </li>
>    >    <li>The end of sequences must be cut. And we can also also remove the beginning of the sequences (the first 13 bases) to improve the per base sequence content</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 3. **Trim Galore** :wrench:: Run Trim Galore on the new dataset to fit the previous choices
> 4. **FastQC** :wrench:: Re-run FastQC to check the impact of Trim Galore
>
>    > ### :question: Questions
>    >
>    > 1. How many sequences have been removed?
>    > 2. Has the quality of the sequences been improved?
>    > 3. Can you explain why the per base sequence content is not good now?
>    > 
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>Before Trim Galore, the dataset was made of 100,000 sequences. After Trim Galore, there is 97,644 sequences.</li>
>    >    <li>The per base quality score is better (not red anymore). But the per base sequence content is still red even if it is a bit better.</li>
>    >    </ol>
>    >    </details>
>    {: .question}
{: .hands_on}

# Conclusion

In this tutorial, we have controlled the quality of two datasets to ensure that the raw data looks good before analysing them with tools to extract RNA-Seq, ChIP-Seq or any other type of information. The approach of quality control is similar for any type of sequencing data:

![](../images/dive_into_qc_workflow.png)

