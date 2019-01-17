---
layout: tutorial_hands_on

title: "Quality Control"
zenodo_link: "https://doi.org/10.5281/zenodo.61771"
questions:
  - How to perform quality control of NGS raw data (FASTQ)?
  - What are the quality parameters to check for a dataset?
  - How to improve the quality of a dataset?
objectives:
  - Manipulate FASTQ files
  - Assess quality from a FASTQ file
  - Use FastQC tool
  - Understand FastQC output
  - Use tools for quality correction
  - Use a tool to aggregate FastQC output
  - Process single-end and paired-end data
follow_up_training:
  -
    type: "internal"
    topic_name: sequence-analysis
    tutorials: 
      - mapping
time_estimation: "1H"
key_points:
  - Run quality control on every dataset before running any other bioinformatics analysis
  - Take care of the parameters used to improve the sequence quality
  - Re-run FastQC to check the impact of the quality control
  - For paired-end reads analyze the forward and reverse reads together
contributors:
  - bebatut
---

# Introduction
{:.no_toc}

During sequencing, the nucleotide bases in a DNA or RNA sample (library) are determined by the sequencer. For each fragment in the library, a short sequence is generated, also called a **read**, which is simply a succession of nucleotides

Modern sequencing technologies can generate a massive number of sequence reads in a single experiment. However, no sequencing technology is perfect, and each instrument will generate different types and amounts of error, such as incorrect nucleotides being called. These are due to the technical limitations of each sequencing platform.

Therefore, it is necessary to understand, identify and exclude error-types that may impact the interpretation of downstream analysis.
Sequence quality control is therefore an essential first step in your analysis. Catching errors early saves time later on.

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

We just imported a file into Galaxy. This file is similar the data we could get directly from a sequencing facility: a [FASTQ file](https://en.wikipedia.org/wiki/FASTQ_format).

> ### {% icon hands_on %} Hands-on: Inspect the FASTQ file
>
> 1. Inspect the file by clicking on the {% icon galaxy-eye %} (eye) icon
>
{: .hands_on}

Although it looks complicated (and maybe it is), the FASTQ format is easy to understand with a little decoding.

Each read, representing a fragment of the library, is encoded by 4 lines:

Line  | Description
--- | ---
1 | Always begins with `@` and then information about the read
2 | The actual nucleic sequence
3 | Always begins with a `+` and sometimes the same info in line 1
4 | Has a string of characters which represent the quality scores associated with each base of the nucleic sequence; must have the same number of characters as line 2

So for example, the first sequence in our file is:

```
@SRR031716.1 HWI-EAS299_4_30M2BAAXX:3:1:944:1798 length=37
GTGGATATGGATATCCAAATTATATTTGCATAATTTG
+SRR031716.1 HWI-EAS299_4_30M2BAAXX:3:1:944:1798 length=37
IIIIIIIIIIIIIIIIIIIIIIIIIIIII8IIIIIII
```

It means that the fragment named `SRR031716.1` corresponds to the DNA sequence `GTGGATATGGATATCCAAATTATATTTGCATAATTTG` and this sequence has been sequenced with a quality `IIIIIIIIIIIIIIIIIIIIIIIIIIIII8IIIIIII`.

But what does this quality mean?

The quality for each sequence is a string of characters, one for each base of the nucleic sequence, used to characterize the probability of mis-indentification of each base. The score is encoded using the ASCII character table (with [some historical differences](https://en.wikipedia.org/wiki/FASTQ_format#Encoding)):

```
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................  Sanger
 ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................  Solexa
 ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................  Illumina 1.3+
 .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................  Illumina 1.5+
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................  Illumina 1.8+
 !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
 |                         |    |        |                              |                     |
33                        59   64       73                            104                   126
 0........................26...31.......40
                          -5....0........9.............................40
                                0........9.............................40
                                   3.....9.............................40
 0.2......................26...31........41
```

So each nucleotide is assigned an ASCII character, representing its [Phred quality score](https://en.wikipedia.org/wiki/Phred_quality_score), the probability of an incorrect base call:

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
> 1. What is the Phred quality score of the 3rd nucleotide of the 1st sequence?
> 2. What is the accuracy of this 3rd nucleotide?
>
> > ### {% icon solution %} Solution
> > 1. The worst Phred score is the smallest one, so 0. For Illumina, it corresponds to the `!` character.
> > 2. The 3rd nucleotide of the 1st sequence has a ASCII character `I`, which correspond to a score of 40.
> > 3. The corresponding nucleotide `G` has an accuracy of 99.99%
> >
> {: .solution }
{: .question}

> ### {% icon comment %} Comment
> The current lllumina (1.8+) uses Sanger format (Phred+33). If you are working with older datasets you may encounter the older scoring schemes. **FastQC** {% icon tool %}, the tool we will use in the next step, can be used to try to determine what type of quality encoding is used (through assessing the range of Phred values seen in the FASTQ).
{: .comment}

When looking at the file in Galaxy, it looks like most the nucleotides have a high score (`I` corresponding to a score 40). Is it true for all sequences? And along the full sequence length?

# Assess the Read Quality

To estimate sequence quality along all sequences, we now use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). It is an open-source tool that provides a simple way to quality control raw sequence data. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.

> ### {% icon hands_on %} Hands-on: Quality check
>
> 1. Run **FastQC** {% icon tool %} with the following parameters
>    - {% icon param-files %} *"Short read data from your current history"*: `reads_1`
>
> 2. Inspect the generated HTML files
>
{: .hands_on}

Rather than looking at quality scores for each individual read, FastQC looks at quality collectively across all reads within a sample:

![Per base sequence quality](../../images/quality-control/per_base_sequence_quality.png "Per base sequence quality")

On the x-axis are the base position in the read, and on the y-axis quality scores. In this example, the sample contains reads that are 37 bp long. For each position, there is a box plotted to illustrate the distribution of values (with the whiskers indicating the 90th and 10th percentile scores). The blue line represents the mean score at the postion. The plot background is color-coded to identify good (green), acceptable (yellow), and bad (red) quality scores.

> ### {% icon question %} Questions
>
> 1. Which Phred encoding is used in the FASTQ file for these sequences?
> 2. How is the mean score changing along the sequence?
> 3. Is this tendency seen in all sequences?
>
> > ### {% icon solution %} Solution
> > 1. The Phred scores are encoded using `Sanger / Illumina 1.9` (`Encoding` in the top table)
> > 2. The mean score over the sequence is dropping at the end of the sequences. It is usual: the sequencers are incorporating more errored nucleotides at the end. But the score stays good: over 28.
> > 3. The box plots are getting wider at the end of the sequences. It means a lot of sequences have their score dropping at the end of the sequence. At after 31 nucleotides, more than 10% of the sequences have scores below 20.
> >
> {: .solution }
{: .question}

FastQC produces other diagnostic plots to assess sample quality.

> ### {% icon details %} More details about the other FastQC plots
>
> The explanation is adapted from [FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/))
>
> - **Per sequence quality scores**
>
>     ![Per sequence quality scores](../../images/quality-control/per_sequence_quality_scores.png "Per sequence quality scores")
>
>     The per sequence quality score reports if a subset of the sequences have universally low quality values. It is often the case that a subset of sequences will have universally poor quality, often because they are poorly imaged (on the edge of the field of view etc), however these should represent only a small percentage of the total sequences.
>
> - **Per tile sequence quality**
>
>     ![Per tile sequence quality](../../images/quality-control/per_tile_sequence_quality.png "Per tile sequence quality")
>
>      This graph will only appear for Illumina library which retains its original sequence identifiers. Encoded in these is the flowcell tile from which each read came. The graph enables you to look at the quality scores from each tile across all of your bases to see if there was a loss in quality associated with only one part of the flowcell.
>
>      The plot shows the deviation from the average quality for each tile. The colours are on a cold to hot scale, with cold colours being positions where the quality was at, or above, the average for that base in the run, and hotter colours indicate that a tile had worse qualities than other tiles for that base. In the example above you can see that certain tiles show consistently poor quality. A good plot should be blue all over.
>
> - **Per base sequence content**
>
>     ![Per base sequence content](../../images/quality-control/per_base_sequence_content.png "Per base sequence content")
>
>     Per Base Sequence Content plots out the proportion of each of the four DNA bases (T, C, A, G) at each position along the reads.
>
>     In a random library we would expect that there would be little to no difference between the four bases, so the lines in this plot should run parallel with each other. The relative amount of each base should reflect the overall amount of these bases in the genome, but in any case they should not be hugely imbalanced from each other.
>
>     It's worth noting that some types of library will always produce biased sequence composition, normally at the start of the read. Libraries produced by priming using random hexamers (including nearly all RNA-Seq libraries), and those which were fragmented using transposases, will contain an intrinsic bias in the positions at which reads start. This bias does not involve a specific sequence, but instead provides enrichment of a number of different K-mers at the 5' end of the reads. Whilst this is a true technical bias, it isn't something which can be corrected by trimming and in most cases doesn't seem to adversely affect the downstream analysis. It will, however, produce a warning or error in this module.
>
> - **Per sequence GC content**
>
>     ![Per sequence GC content](../../images/quality-control/per_sequence_gc_content.png "Per sequence GC content")
>
>     This plot displays the GC content across the whole length of each sequence in a file and compares it to a modelled normal distribution of GC content.
>
>     In a normal random library we would expect to see a roughly normal distribution of GC content, where the central peak corresponds to the overall GC content of the underlying genome. Since we don't know the GC content of the genome the modal GC content is calculated from the observed data and used to build a reference distribution.
>
>     An unusually-shaped distribution could indicate a contaminated library or some other kind of biased subset. A normal distribution which is shifted indicates some systematic bias, which is independent of base position. If there is a systematic bias which creates a shifted normal distribution then this won't be flagged as an error by the module since it doesn't know what your genome's GC content should be.
>
> - **Per base N content**
>
>     ![Per base N content](../../images/quality-control/per_base_n_content.png "Per base N content")
>
>     If a sequencer is unable to make a base call with sufficient confidence then it will normally substitute an N rather than a conventional base call.
>
>     It plots out the percentage of base calls at each position for which an N was called.
>
>     It's not unusual to see a very low proportion of Ns appearing in a sequence, especially near the end of a sequence. However, if this proportion rises above a few percent it suggests that the analysis pipeline was unable to interpret the data well enough to make valid base calls.
>
> - **Sequence length distribution**
>
>     ![Sequence length distribution](../../images/quality-control/sequence_length_distribution.png "Sequence > length distribution")
>
>     Some high throughput sequencers generate sequence fragments of uniform length, but others can contain reads of wildly varying lengths. Even within uniform length libraries some pipelines will trim sequences to remove poor quality base calls from the end.
>
>     This graph shows the distribution of fragment sizes in the file which was analysed. In many cases this will produce a simple graph showing a peak only at one size, but for variable length FASTQ files this will show the relative amounts of each different size of sequence fragment.
>
> - **Sequence Duplication Levels**
>
>     ![Sequence Duplication Levels](../../images/quality-control/sequence_duplication_levels.png "Sequence Duplication Levels")
>
>     In a diverse library most sequences will occur only once in the final set. A low level of duplication may indicate a very high level of coverage of the target sequence, but a high level of duplication is more likely to indicate some kind of enrichment bias (eg PCR over amplification).
>
>     FastQC counts the degree of duplication for every sequence in a library and creates a plot showing the relative number of sequences with different degrees of duplication.
>
>     There are two lines on the plot. The blue line takes the full sequence set and shows how its duplication levels are distributed. In the red plot the sequences are de-duplicated and the proportions shown are the proportions of the deduplicated set which come from different duplication levels in the original data.
>
>     In a properly diverse library most sequences should fall into the far left of the plot in both the red and blue lines. A general level of enrichment, indicating broad oversequencing in the library will tend to flatten the lines, lowering the low end and generally raising other categories. More specific enrichments of subsets, or the presence of low complexity contaminants will tend to produce spikes towards the right of the plot. These high duplication peaks will most often appear in the blue trace as they make up a high proportion of the original library, but usually disappear in the red trace as they make up an insignificant proportion of the deduplicated set. If peaks persist in the blue trace then this suggests that there are a large number of different highly duplicated sequences which might indicate either a contaminant set or a very severe technical duplication.
>
> - **Overrepresented sequences**
>
>      A normal high-throughput library will contain a diverse set of sequences, with no individual sequence making up a tiny fraction of the whole. Finding that a single sequence is very overrepresented in the set either means that it is highly biologically significant, or indicates that the library is contaminated, or not as diverse as you expected.
>
>      FastQC lists all of the sequence which make up more than 0.1% of the total. For each overrepresented sequence FastQC will look for matches in a database of common contaminants and will report the best hit it finds. Hits must be at least 20bp in length and have no more than 1 mismatch. Finding a hit doesn't necessarily mean that this is the source of the contamination, but may point you in the right direction. It's also worth pointing out that many adapter sequences are very similar to each other so you may get a hit reported which isn't technically correct, but which has very similar sequence to the actual match.
>
> - **Adapter Content**
>
>     ![Adapter Content](../../images/quality-control/adapter_content.png "Adapter Content")
>
>     FastQC does a generic analysis of all of the Kmers in the library to find those which do not have even coverage through the length of your reads. This can find a number of different sources of bias in the library which can include the presence of read-through adapter sequences building up on the end of your sequences.
>
>     The presence of any overrepresented sequences in your library (such as adapter dimers) causes the Kmer plot to be dominated by the Kmers these sequences contain, and then it's not always easy to see if there are other biases present in which you might be interested.
>
>     The plot shows a cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position. Once a sequence has been seen in a read it is counted as being present right through to the end of the read so the percentages you see will only increase as the read length goes on.
>
{: .details}


> ### {% icon question %} Questions
>
> Why is there a warning for the per-base sequence content and the per-sequence GC content graphs?
>
> > ### {% icon solution %} Solution
> > In the beginning of sequences, the sequence content per base is not really good and the percentages are not equal. For the GC content, the distribution is slightly shifted on the left, and too high.
> {: .solution }
{: .question}


# Trim and filter

Based on the information provided by the quality graphs, the quality of the sequences drops at the end of the sequences. It could cause bias in downstream analyis with these potentially wrong nucleotides.

The sequences should be treated to reduce bias in downstream analyis. In general, quality treatments include:

- Filtering of sequences
    - with low mean quality score
    - too short
    - with too many ambiguous (N) bases
    - based on their GC content
- Cutting/Trimming/masking sequences
    - from low quality score regions
    - beginning/end of sequence
    - removing adapters

To accomplish this task we use [Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) tool that enhances sequence quality by automating adapter trimming as well as quality control.

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
>          If you know which adapter sequences were used during library preparation, provide their sequences there.
>
>    - In *"Adapter Options"*
>       - *"Minimum overlap length"*: `3`
>
>    - In *"Filter Options"*
>       - *"Minimum length"*: `20`
>
>           To remove reads shorter than 20 bp
>
>    - In *"Read Modification Options"*
>       - *"Quality cutoff"*: `20`
>
>           To trim low-quality 3' ends (below 20) from reads in addition to adapter removal
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
> 2. (Optional) **FastQC** {% icon tool %}: Re-run **FastQC** on the quality-controlled data, and inspect the new FastQC report
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
> If the quality of the reads is not good, we should always first check what is wrong and think about it: it may come from the type of sequencing or what we sequenced (high quantity of overrepresented sequences in transcriptomics data, biaised percentage of bases in HiC data).
>
> You can also ask the sequencing facility about it, especially if the quality is really bad: the quality treatments can not solve everything and you may lose too much information.
{: .comment}

# Process paired-end data

With paired-end sequencing, the fragments are sequenced from both sides. This approach results in two reads per fragment, with the first read in forward orientation and the second read in reverse-complement orientation. The distance between both reads is known and therefore is additional information that can improve read mapping.

With paired-end sequencing, more of each DNA fragment can be covered than with single-end sequencing (only forward orientation sequenced):

```
    ------>
      ------>
        ------>
          ------>

    ----------------------------- [fragment]

    ------>         <------
      ------>         <------
        ------>         <------
          ------>         <------
```

Paired-end sequencing generates 2 FASTQ files:
- One file with the sequences corresponding to foward orientation of all the fragments
- One file with the sequences corresponding to reverse orientation of all the fragments

The data we analyzed in the previous step was not single-end data but the forward reads of paired-end data. We will now do the quality control on both forward and reverse reads together.

> ### {% icon question %} Questions
>
> Why is performing quality control on both forward and reverse reads together important?
>
> > ### {% icon solution %} Solution
> > During the filtering, some reads are eliminated because of their length. If one of the reverse reads is removed, its corresponding forward read should be removed too. So the forward and reverse reads should be processed together.
> {: .solution }
{: .question}

Let's first have a look at the quality of our reads!

> ### {% icon hands_on %} Hands-on: Assessing the quality of paired-end data
>
> 1. Import the reverse read `GSM461178_untreat_paired_subset_2.fastq` from [Zenodo](https://zenodo.org/record/61771) or from the data library (ask your instructor)
>
>    ```
>    https://zenodo.org/record/61771/files/GSM461178_untreat_paired_subset_2.fastq
>    ```
>
> 2. Rename the file to `reads_2`
> 3. **FastQC** {% icon tool %} with the reverse reads
>
>     FastQC has already been run on the forward reads
>
> 4. **MultiQC** {% icon tool %} with the following parameters to aggregate the FastQC reports
>      - In *"Results"*
>        - *"Which tool was used generate logs?"*: `FastQC`
>        - In *"FastQC output"*
>           - *"Type of FastQC output?"*: `Raw data`
>           - {% icon param-files %} *"FastQC output"*: `Raw data` files (output of both **FastQC** {% icon tool %})
>
>    {% include snippets/select_multiple_datasets.md %}
>
> 5. Inspect the webpage output from MultiQC
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What do you think of the quality of the sequences?
>    > 2. What should we do?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > 1. The quality of the sequences seems worse for the reverse reads than for the forward reads: lower mean quality of the sequences and stronger decrease at the end (mean value below 28). Lower-quality reverse reads .
>    > > 2. We should trim the end of the sequences and filter them with **Cutadapt** {% icon tool %}
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

It is usual that the quality of the sequences is worse for the reverse than for the forward reads. Therefore it is important to treat the forward and reverse reads together.

> ### {% icon hands_on %} Hands-on: Assessing the quality of paired-end data
> 1. **Cutadapt** {% icon tool %} with the following parameters
>    - *"Single-end or Paired-end reads?"*: `Paired-end`
>       - {% icon param-file %} *"FASTQ/A file #1"*: `reads_1` (Input dataset)
>       - {% icon param-file %} *"FASTQ/A file #2"*: `reads_2` (Input dataset)
>
>          The order is important here!
>
>       - In *Read 1 Options* or *Read 2 Options*
>
>          If you know which adapter sequences were used during library preparation, provide their sequences there.
>
>    - In *"Adapter Options"*
>       - *"Minimum overlap length"*: `3`
>    - In *"Filter Options"*
>       - *"Minimum length"*: `20`
>    - In *"Filter Options"*
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
>    > > 2. 1,376 sequences have been removed because at least one read was shorter than the length cutoff (334 when only the forward reads were analyzed).
>    > {: .solution }
>    {: .question}
{: .hands_on}


> ### {% icon details %} Algorithmic details
>
> One of the biggest advantage of cutadapt over other trimming tools is that it has a nice [documentation](https://cutadapt.readthedocs.io) explaining how the tool works in detail.
>
> Cutadapt quality trimming algorithm consists of three simple steps:
>
> 1. substract the chosen threshold value from the quality value of each position
> 2. compute a partial sum of these differences from the end of the sequence to each position
>    (as long as the partial sum is negative)
> 3. cut at the minimum value of the partial sum
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
> 2. Add up the numbers, starting from the end (partial sums) and stop early if the sum is greater than zero
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
> Note that thereby also positions with a quality value larger than the chosen threshold are removed if they are embedded in regions with lower quality (note that the partial sum is decreasing if the quality values are smaller than the threshold). The advantage of this procedure is that it is robust against a small number of positions with a quality higher than the threshold.
>
>
> Alternatives to this procedure would be
>
> * to cut all positions with a quality smaller than the threshold
> * a sliding window approach
>
> The sliding window approach checks that the average quality of each sequence window of specified length is larger than the threshold. Note that in contrast to cutadapt's approach, this approach has one more parameter and the robustness depends of the length of the windows (in combination with the quality threshold). Both approaches are implemented in Trimmomatic.
{: .details}

> ### {% icon question %} Questions
>
> 1. What kind of alignment is used for finding adapters in reads?
> 2. What is the criterion to choose the best adapter alignment?
>
> > ### {% icon solution %} Solution
> >
> > 1. Semi global alignment, i.e., only the overlapping part of the read and the adapter sequence is used for scoring.
> > 2. An alignment with maximum overlap is computed that has a smallest number of mismatches and indels.
> >
> {: .solution}
{: .question}


# Conclusion
{:.no_toc}

In this tutorial we checked the quality of two FASTQ files to ensure that their data looks good before inferring any further information. This step is the usual first step for analyses such as RNA-Seq, ChIP-Seq, or any other OMIC analysis relying on NGS data. Quality control steps are similar for any type of sequencing data:

- Quality assessment with a tool like **FastQC** {% icon tool %}
- Trimming and filtering with a tool like **Cutadapt** {% icon tool %}
