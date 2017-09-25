---
layout: tutorial_hands_on
topic_name: sequence-analysis
tutorial_name: quality-control
---

# Introduction
{:.no_toc}

During sequencing, errors might be introduced, such as the incorporation of ambiguous nucleotides. These are due to the technical limitations of each sequencing platform. Sequencing errors might bias the analysis, ultimately leading to a misinterpretation of the data.

Sequence quality control is therefore an essential step to take right after receiving the raw sequencing data. This ensures proper, analysis regardless the sequencing platform used to obtain the data.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Importing sequences

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history
> 2. Import the FASTQ file: [`GSM461178_untreat_paired_subset_1`](https://zenodo.org/record/61771/files/GSM461178_untreat_paired_subset_1.fastq)
>
>    > ### {% icon tip %} Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * Press **Start**    
>    {: .tip}
>
>    > ### {% icon tip %} Tip: Change the file type `fastq` to `fastqsanger` once the data file is in your history
>    >
>    > * Click on the pencil button displayed in your data file in the history
>    > * Choose **Datatype** on the top
>    > * Select `fastqsanger`
>    > * Press **save**
>    {: .tip}
>
>    > ### {% icon comment %} Comments
>    >
>    > Rename the dataset to "First dataset"
>    {: .comment}
> By default, when data is imported via its link, Galaxy names it with its URL.
{: .hands_on}

# Quality check

To estimate sequence quality and how to further filter raw data, different indicators can be checked:

- Quality score of the sequences with
    - Per-base sequence quality
    - Per-sequence quality scores
    - Per-tile sequence quality
- Sequence content with
    - Per-base sequence content
    - Per-sequence GC content
    - Per-base N content
- Sequence length with the sequence length distribution
- Duplicated sequences
- Tag sequences with
    - Adapter contamination
    - Kmer Content

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is an open-source tool that provides a simple way to quality-control raw sequence data coming from high throughput sequencing pipelines. It gets rid of low quality score reads, and generates graphics and estimates providing a quick overview about which data might be the source of bias in the analysis.

> ### {% icon hands_on %} Hands-on: Quality check
>
> 1. **FastQC** {% icon tool %}: Run FastQC on the imported FastQ file with default parameters
> 2. Inspect the FastQC report on its webpage output
>
>    > ### {% icon tip %} Tip: Inspecting the content of a file in Galaxy
>    >
>    > * Click on the eye ("View data") on the right of the file name in the history
>    > * Inspect the content of the file on the middle
>    {: .tip}
>
>    > ### {% icon question %} Questions
>    >
>    > 1. How good are the quality scores?
>    > 2. Why is there a warning for the per-base sequence content and the per-sequence GC content graphs?
>    > 3. What needs to be done to improve the sequences?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>The sequence scores are quite good: no warnings from FastQC, even if we can see a slight decrease of the quality at the end of the reads</li>
>    >    <li>In the beginning of sequences, the sequence content per base is not really good and the percentages are not equal. For the GC content, the distribution is slightly shifted on the left, and too high</li>
>    >    <li>We can trim the end of the sequences a little, but not too much as the sequences are already small</li>
>    >    </ol>
>    >    </details>
>    {: .question}
{: .hands_on}

# Improvement of sequence quality

Based on the informations provided by the quality graphs, the sequences must to be treated to avoid bias in downstream analyis.

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

To improve the overall sequence quality, we use the [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) tool. This tool enhances sequence quality by automating adapter trimming as well as quality control.

> ### {% icon hands_on %} Hands-on: Quality treatment and re-checking
>
> 1. **Trim Galore** {% icon tool %}: Run Trim Galore on the imported data file
>
>    > ### {% icon question %} Questions
>    >
>    > Which parameters must be applied to follow the previous recommendations?
>    >
>    > <details>
>    > <summary>Click to view answers</summary>
>    > We use the default ones:
>    > <ul>
>    > <li>​
If you already know the which adapter sequences were used during the library preparation, please use them. Otherwise, use the option for automatic detection and trimming of adapter sequences</li>
>    > <li>Trimming low-quality ends (below 20) from reads in addition to adapter removal</li>
>    > <li>Option for required number bases overlap with adapter sequence can be tweaked. The default value "1" is too stringent that on average 25% of reads will be trimmed. In order to reduce these falsely trimmed bases, please set it to 5 bases.</li>
>    > <li>Removing reads shorter than 20 bp</li>
>    > </ul>
>    > </details>
>    {: .question}
>
> 2. **FastQC** {% icon tool %}: Re-run FastQC on the quality controlled data file and inspect the new FastQC report
>
>    > ### {% icon question %} Questions
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

> ### {% icon hands_on %} Hands-on: Quality control and treatment
>
> 1. Create a new history
> 2. Import the FASTQ file: [`GSM461182_untreat_single_subset`](https://zenodo.org/record/61771/files/GSM461182_untreat_single_subset.fastq)
> 3. **FastQC** {% icon tool %}: Run FastQC on this new dataset
>
>    > ### {% icon question %} Questions
>    >
>    > 1. How good is this dataset?
>    > 2. What must be done to improve the sequences?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>There is red warning for the per base sequence quality (pretty bad along the sequence but worst at the end of sequences), the per base sequence content (bad at the beginning of the sequences), the per sequence GC content. </li>
>    >    <li>The end of sequences must be cut.</li>
>    >    <li>Generally, 5' ends of the reads are not of bad quality unless there is something went wrong. The problem with this particular data set is that it was sequenced using the old Illumina sequencing machine. The machine calibrates while reading fragments that are in the beginning of the flowcell. Unfortunately, the first 100k reads which we selected for the analysis are generated during the calibration. But with the latest sequencing machines, usually we do not see this problem. If you used latest sequencing machine and still see bad quality bases in the beginning of the reads, please investigate and not just trim them.</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 3. **Trim Galore** {% icon tool %}: Run Trim Galore on the new dataset to fit the previous choices
> 4. **FastQC** {% icon tool %}: Re-run FastQC to check the impact of Trim Galore
>
>    > ### {% icon question %} Questions
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
{:.no_toc}

In this tutorial, we have controlled the quality of two datasets to ensure that the raw data looks good before analysing them with tools to extract RNA-Seq, ChIP-Seq or any other type of information. The approach of quality control is similar for any type of sequencing data:

![The quality control tutorial workflow](../../images/dive_into_qc_workflow.png)
