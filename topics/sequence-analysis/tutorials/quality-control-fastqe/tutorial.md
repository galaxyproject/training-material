---
layout: tutorial_hands_on

title: "Quality Control with FASTQE 🧬😎"
zenodo_link: "https://doi.org/10.5281/zenodo.61771"
questions:
  - How to perform quality control of NGS raw data (FASTQ) with FASTQE?
  - What are the quality parameters to check for a dataset?
  - How to improve the quality of a dataset?
objectives:
  - Manipulate FASTQ files
  - Assess quality from a FASTQ file
  - Use FASTQE tool
  - Understand FASTQE output
  - Use tools for quality correction
  - Use a tool to aggregate FASTQE output
  - Process single-end and paired-end data
follow_up_training:
  -
    type: "internal"
    topic_name: sequence-analysis
    tutorials:
      - mapping
time_estimation: "1H30M"
level: Introductory
key_points:
  - Run quality control on every dataset before running any other bioinformatics analysis
  - Take care of the parameters used to improve the sequence quality
  - Re-run FASTQE to check the impact of the quality control
  - For paired-end reads analyze the forward and reverse reads together
contributors:
  - mblue9
---

# Introduction
{:.no_toc}

During sequencing, the nucleotide bases in a DNA or RNA sample (library) are determined by the sequencer. For each fragment in the library, a short sequence is generated, also called a **read**, which is simply a succession of nucleotides.

Modern sequencing technologies can generate a massive number of sequence reads in a single experiment. However, no sequencing technology is perfect, and each instrument will generate different types and amount of errors, such as incorrect nucleotides being called. These wrongly called bases are due to the technical limitations of each sequencing platform.

Therefore, it is necessary to understand, identify and exclude error-types that may impact the interpretation of downstream analysis.
Sequence quality control is therefore an essential first step in your analysis. Catching errors early saves time later on.

This FASTQE tutorial is adapted from the [Galaxy Quality Control]({% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %}) tutorial.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Inspect a raw sequence file

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial and give it a proper name
>
>    {% include snippets/create_new_history.md %}
>    {% include snippets/rename_history.md %}
>
> 2. Import `GSM461178_untreat_paired_subset_1.fastq` from [Zenodo](https://zenodo.org/record/61771) or from the data library (ask your instructor)
>
>    ```
>    https://zenodo.org/record/61771/files/GSM461178_untreat_paired_subset_1.fastq
>    ```
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
>    As default, Galaxy takes the link as name, so rename them.
>
> 4. Rename the file to `reads_1`
>
>    {% include snippets/rename_dataset.md %}
>
{: .hands_on}

We just imported a file into Galaxy. This file is similar to the data we could get directly from a sequencing facility: a [FASTQ file](https://en.wikipedia.org/wiki/FASTQ_format).

> ### {% icon hands_on %} Hands-on: Inspect the FASTQ file
>
> 1. Inspect the file by clicking on the {% icon galaxy-eye %} (eye) icon
>
{: .hands_on}

Although it looks complicated (and maybe it is), the FASTQ format is easy to understand with a little decoding.

Each read, representing a fragment of the library, is encoded by 4 lines:

Line  | Description
--- | ---
1 | Always begins with `@` followed by the information about the read
2 | The actual nucleic sequence
3 | Always begins with a `+` and contains sometimes the same info in line 1
4 | Has a string of characters which represent the quality scores associated with each base of the nucleic sequence; must have the same number of characters as line 2

So for example, the first sequence in our file is:

```
@SRR031716.1 HWI-EAS299_4_30M2BAAXX:3:1:944:1798 length=37
GTGGATATGGATATCCAAATTATATTTGCATAATTTG
+SRR031716.1 HWI-EAS299_4_30M2BAAXX:3:1:944:1798 length=37
IIIIIIIIIIIIIIIIIIIIIIIIIIIII8IIIIIII
```

It means that the fragment named `SRR031716.1` corresponds to the DNA sequence `GTGGATATGGATATCCAAATTATATTTGCATAATTTG` and this sequence has been sequenced with a quality `IIIIIIIIIIIIIIIIIIIIIIIIIIIII8IIIIIII`.

But what does this quality score mean?

The quality score for each sequence is a string of characters, one for each base of the nucleic sequence, used to characterize the probability of mis-identification of each base. The score is encoded using the ASCII character table (with [some historical differences](https://en.wikipedia.org/wiki/FASTQ_format#Encoding)):

![Encoding of the quality score with ASCII characters for different Phred encoding](../../images/quality_score_encoding.png)

So there is an ASCII character associated with each nucleotide, representing its [Phred quality score](https://en.wikipedia.org/wiki/Phred_quality_score), the probability of an incorrect base call:

Phred Quality Score | Probability of incorrect base call | Base call accuracy
--- | --- | ---
10 | 1 in 10 | 90%
20 | 1 in 100 | 99%
30 | 1 in 1000 | 99.9%
40 | 1 in 10,000 | 99.99%
50 | 1 in 100,000 | 99.999%
60 | 1 in 1,000,000 | 99.9999%

> ### {% icon question %} Questions
>
> 1. Which ASCII character corresponds to the worst Phred score for Illumina 1.8+?
> 2. What is the Phred quality score of the 3rd nucleotide of the 1st sequence?
> 3. What is the accuracy of this 3rd nucleotide?
>
> > ### {% icon solution %} Solution
> > 1. The worst Phred score is the smallest one, so 0. For Illumina, it corresponds to the `!` character.
> > 2. The 3rd nucleotide of the 1st sequence has a ASCII character `I`, which correspond to a score of 40.
> > 3. The corresponding nucleotide `G` has an accuracy of 99.99%
> >
> {: .solution }
{: .question}

When looking at the file in Galaxy, it looks like most the nucleotides have a high score (`I` corresponding to a score 40). Is it true for all sequences? And along the full sequence length?

# Assess the Read Quality

To estimate sequence quality along all sequences, we now use [FASTQE](https://fastqe.com/). It is an open-source tool that provides a simple way to quality control raw sequence data and print them as emoji. You can use it to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.

> ### {% icon hands_on %} Hands-on: Quality check
>
> 1. Run **FASTQE** {% icon tool %} with the following parameters
>    - {% icon param-files %} *"FastQ data"*: `reads_1`
>
> 2. Inspect the generated HTML file
>
{: .hands_on}

You can see what score each emoji corresponds to [here](https://github.com/fastqe/fastqe)

> ### {% icon question %} Questions
>
> What is the lowest score in the datasest?
>
> > ### {% icon solution %} Solution
> > The lowest score in this dataset is 1 ❌.
> {: .solution }
{: .question}

## Per base sequence quality

Rather than looking at quality scores for each individual read, FASTQE looks at quality collectively across all reads within a sample and calculates the mean, maximum and minimum for each nucleotide position along the length of the reads.

TODO: add image of output


It is normal with all Illumina sequencers for the quality score to start out lower over the first 5-7 bases and to then rise. The quality of reads on most platforms will drop at the end of the read. This is often due to signal decay or phasing during the sequencing run. The recent developments in chemistry applied to sequencing has improved this somewhat, but reads are now longer than ever.


> ### {% icon details %} Signal decay and phasing
>
> - Signal decay
>
>  The fluorescent signal intensity decays with each cycle of the sequencing process. Due to the degrading fluorophores, a proportion of the strands in the cluster are not being elongated. The proportion of the signal being emitted continues to decrease with each cycle, yielding to a decrease of quality scores at the 3' end of the read.
>
> - Phasing
>
>  The signal starts to blur with the increase of number of cycles because the cluster looses synchronicity. As the cycles progress, some strands get random failures of nucleotides to incorporate due to:
>
>  - Incomplete removal of the 3' terminators and fluorophores
>  - Incorporation of nucleotides without effective 3' terminators
>
>  This leads to a decrease in quality scores at the 3' end of the read.
{: .details}


> ### {% icon question %} Questions
>
> 1. How is the mean score changing along the sequence?
>
> > ### {% icon solution %} Solution
> > 1. The mean score over the sequence is dropping at the end of the sequences. This is very common: the sequencers are incorporating more incorrect nucleotides at the end. But the overall score stays good: over 28.
> >
> {: .solution }
{: .question}

When the median quality is below a FASTQE score of ~⚠️ (Phred score of 20), we should consider trimming away bad quality bases from the sequence. We will explain that process in the next section.



# Filter and Trim

The quality of the sequences drops at the end of the sequences. This could cause bias in downstream analyses with these potentially incorrectly called nucleotides. Sequences must be treated to reduce bias in downstream analysis. In general, quality treatments include:

1. Cutting/Trimming/masking sequences
    - from low quality score regions
    - beginning/end of sequence
    - removing adapters
2. Filtering of sequences
    - with low mean quality score
    - too short
    - with too many ambiguous (N) bases

To accomplish this task we will use [Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html), a tool that enhances sequence quality by automating adapter trimming as well as quality control.

> ### {% icon hands_on %} Hands-on: Improvement of sequence quality
>
> 1. Run **Cutadapt** {% icon tool %} with the following parameters
>    - *"Single-end or Paired-end reads?"*: `Single-end`
>       - {% icon param-file %} *"Reads in FASTQ format"*: `reads_1` (Input dataset)
>
>          > ### {% icon tip %} Tip: Files not selectable?
>          > If your FASTQ files cannot be selected, you might check whether their format is FASTQ with Sanger-scaled quality values (`fastqsanger`). You can edit the data type by clicking on the pencil symbol.
>          {: .tip}
>
>       - In *Read 1 Options*
>
>          > ### {% icon comment %} Known adapters
>          > If you see or know which adapter sequences were used during library preparation, provide their sequences there.
>          {: .comment}
>
>    - In *"Filter Options"*
>       - *"Minimum length"*: `20`
>
>           It will remove reads that are shorter than 20 bp, after trimming of adapters and bad regions.
>
>    - In *"Read Modification Options"*
>       - *"Quality cutoff"*: `20`
>
>           After adapter removal (if any), we choose to remove ends ( 5' and/or 3') with low-quality, here below 20 in quality).
>
>    - In *"Output Options"*
>       - *"Report"*: `Yes`
>
> 2. Inspect the generated txt file (`Report`)
>
>    > ### {% icon question %} Questions
>    >
>    > 1. How many reads have been found with adapters?
>    > 2. How many basepairs have been removed from the reads because of bad quality?
>    > 3. How many sequence pairs have been removed because they were too short?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. 0 reads with adapters
>    > > 2. 44,164 bp (1.2%) (`Quality-trimmed:`)
>    > > 3. 322 sequences
>    > {: .solution }
>    {: .question}
>
> 2. (Optional) **FASTQE** {% icon tool %}: Re-run **FASTQE** on the quality-controlled data, and inspect the new FASTQE report
>
>    > ### {% icon question %} Questions
>    >
>    > 1. How many sequences have been removed?
>    > 2. Has sequence quality been improved?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. Before Cutadapt, the dataset comprised 100,000 sequences. After Cutadapt, there are 99,678 sequences
>    > > 2. The per-base quality score looks better, but other indicators show bad values now. The sequence length distribution is not clear anymore because sequences have different size after the trimming operation
>    > {: .solution }
>    {: .question}
{: .hands_on}

The quality of the previous dataset was pretty good from the beginning and we improved it with with trimming and filtering step (in a reasonable way to not lose too much information)

> ### {% icon comment %} Bad quality sequences
> If the quality of the reads is not good, we should always first check what is wrong and think about it: it may come from the type of sequencing or what we sequenced (high quantity of overrepresented sequences in transcriptomics data, biased percentage of bases in HiC data).
>
> You can also ask the sequencing facility about it, especially if the quality is really bad: the quality treatments can not solve everything. If too many bad quality bases are cut away, the corresponding reads then will be filtered out and you loose them.
{: .comment}

> ### {% icon details %} Trimming with Cutadapt
>
> One of the biggest advantage of Cutadapt compared to other trimming tools (e.g. TrimGalore!) is that it has a good [documentation](https://cutadapt.readthedocs.io) explaining how the tool works in detail.
>
> Cutadapt quality trimming algorithm consists of three simple steps:
>
> 1. Subtract the chosen threshold value from the quality value of each position
> 2. Compute a partial sum of these differences from the end of the sequence to each position
>    (as long as the partial sum is negative)
> 3. Cut at the minimum value of the partial sum
>
> In the following example, we assume that the 3’ end is to be quality-trimmed with a threshold of 10 and we have the following quality values
>
> ```
> 42 40 26 27 8 7 11 4 2 3
> ```
>
> 1. Subtract the threshold
>
>     ```
>     32 30 16 17 -2 -3 1 -6 -8 -7
>     ```
>
> 2. Add up the numbers, starting from the 3' end (partial sums) and stop early if the sum is greater than zero
>
>     ```
>     (70) (38) 8 -8 -25 -23 -20, -21 -15 -7
>     ```
>
>     The numbers in parentheses are not computed (because 8 is greater than zero), but shown here for completeness.
>
> 3. Choose the position of the minimum (`-25`) as the trimming position
>
> Therefore, the read is trimmed to the first four bases, which have quality values
>
> ```
> 42 40 26 27
> ```
>
> Note that thereby also positions with a quality value larger than the chosen threshold are removed if they are embedded in regions with lower quality (the partial sum is decreasing if the quality values are smaller than the threshold). The advantage of this procedure is that it is robust against a small number of positions with a quality higher than the threshold.
>
>
> Alternatives to this procedure would be:
>
> * Cut after the first position with a quality smaller than the threshold
> * Sliding window approach
>
>     The sliding window approach checks that the average quality of each sequence window of specified length is larger than the threshold. Note that in contrast to cutadapt's approach, this approach has one more parameter and the robustness depends of the length of the window (in combination with the quality threshold). Both approaches are implemented in Trimmomatic.
{: .details}

# Process paired-end data

With paired-end sequencing, the fragments are sequenced from both sides. This approach results in two reads per fragment, with the first read in forward orientation and the second read in reverse-complement orientation. With this technique, we have the advantage to get more information about each DNA fragment compared to reads sequenced by only single-end sequencing:

```
    ------>                       [single-end]

    ----------------------------- [fragment]

    ------>               <------ [paired-end]
```
The distance between both reads is known and therefore is additional information that can improve read mapping.

Paired-end sequencing generates 2 FASTQ files:
- One file with the sequences corresponding to **forward** orientation of all the fragments
- One file with the sequences corresponding to **reverse** orientation of all the fragments

Usually we recognize these two files which belong to one sample by the name which has the same identifier for the reads but a different extension, e.g. `sampleA_R1.fastq` for the forward reads and `sampleA_R2.fastq` for the reverse reads. It can also be `_f` or `_1` for the forward reads and `_r` or `_2` for the reverse reads.

The data we analyzed in the previous step was not single-end data but the forward reads of paired-end data. We will now do the quality control on the reverse reads.

> ### {% icon hands_on %} Hands-on: Assessing the quality of paired-end reads
>
> 1. Import the reverse read `GSM461178_untreat_paired_subset_2.fastq` from [Zenodo](https://zenodo.org/record/61771) or from the data library (ask your instructor)
>
>    ```
>    https://zenodo.org/record/61771/files/GSM461178_untreat_paired_subset_2.fastq
>    ```
>
> 2. Rename the file to `reads_2`
> 3. **FASTQE** {% icon tool %} with the reverse reads
> 4. Inspect the webpage output from FASTQE
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. What do you think about the quality of the sequences?
> 2. What should we do?
>
> > ### {% icon solution %} Solution
> >
> > 1. The quality of the sequences seems worse for the reverse reads than for the forward reads:
> >     - Per Sequence Quality Scores: distribution more on the left, i.e. a lower mean quality of the sequences
> >     - Per base sequence quality: less smooth curve and stronger decrease at the end with a mean value below 28
> >     - Per Base Sequence Content: stronger bias at the beginning and no clear distinction between C-G and A-T groups
> >
> >    The other indicators (adapters, duplication levels, etc) are similar.
> >
> > 2. We should trim the end of the sequences and filter them with **Cutadapt** {% icon tool %}
> >
> {: .solution}
{: .question}

With paired-end reads the average quality scores for forward reads will almost always be higher than for reverse reads.

After trimming, reverse reads will be shorter because of their quality and then will be eliminated during the filtering step. If one of the reverse reads is removed, its corresponding forward read should be removed too. Otherwise we will get different number of reads in both files and in different order, and order is important for the next steps. Therefore **it is important to treat the forward and reverse reads together for trimming and filtering**.

> ### {% icon hands_on %} Hands-on: Improving the quality of paired-end data
> 1. **Cutadapt** {% icon tool %} with the following parameters
>    - *"Single-end or Paired-end reads?"*: `Paired-end`
>       - {% icon param-file %} *"FASTQ/A file #1"*: `reads_1` (Input dataset)
>       - {% icon param-file %} *"FASTQ/A file #2"*: `reads_2` (Input dataset)
>
>          The order is important here!
>
>       - In *Read 1 Options* or *Read 2 Options*
>
>         As before, no adapters were found in these datasets. When you process your own data and you know which adapter sequences were used during library preparation, you should provide their sequences here.
>
>    - In *"Filter Options"*
>       - *"Minimum length"*: `20`
>    - In *"Read Modification Options"*
>       - *"Quality cutoff"*: `20`
>    - In *"Output Options"*
>       - *"Report"*: `Yes`
>
> 2. Inspect the generated txt file (`Report`)
>
>    > ### {% icon question %} Questions
>    >
>    > 1. How many basepairs has been removed from the reads because of bad quality?
>    > 2. How many sequence pairs have been removed because they were too short?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. 44,164 bp (`Quality-trimmed:`) for the forward reads and 138,638 bp for the reverse reads.
>    > > 2. 1,376 sequences have been removed because at least one read was shorter than the length cutoff (322 when only the forward reads were analyzed).
>    > {: .solution }
>    {: .question}
>
{: .hands_on}

In addition to the report, Cutadapt generates 2 files:
- Read 1 with the trimmed and filtered forward reads
- Read 2 with the trimmed and filtered reverse reads

These datasets can be used for the downstream analysis, e.g. mapping.

> ### {% icon question %} Questions
>
> 1. What kind of alignment is used for finding adapters in reads?
> 2. What is the criterion to choose the best adapter alignment?
>
> > ### {% icon solution %} Solution
> >
> > 1. Semi-global alignment, i.e., only the overlapping part of the read and the adapter sequence is used for scoring.
> > 2. An alignment with maximum overlap is computed that has the smallest number of mismatches and indels.
> >
> {: .solution}
{: .question}

# Conclusion
{:.no_toc}

In this tutorial we checked the quality of two FASTQ files to ensure that their data looks good before inferring any further information. This step is the usual first step for analyses such as RNA-Seq, ChIP-Seq, or any other OMIC analysis relying on NGS data. Quality control steps are similar for any type of sequencing data:

- Quality assessment with a tool like **FASTQE** {% icon tool %}
- Trimming and filtering with a tool like **Cutadapt** {% icon tool %}
