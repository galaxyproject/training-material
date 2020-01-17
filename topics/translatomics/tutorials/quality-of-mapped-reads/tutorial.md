---
layout: tutorial_hands_on

title: Quality of Ribo-Seq data
zenodo_link: ''
questions:
- How is the quality of the mapping results?
- What is the triplet nucleotide periodicity?
- How is the distribution of read lengths?
objectives:
- Learn quality control metrics of Ribo-Seq post-alignment data
- Understand the characteristic of Ribo-Seq data
time_estimation: '1h'
key_points:
- Triplet nucleotide periodicity is an important characteristic for Ribo-Seq data
- Read length of Ribo-Seq data is usually concentrated around 30 nt.
contributors:
- ldyang14
- IceApink
---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

Although we got the corresponding file, the BAM format file, from the mapping results, we have no way of knowing specifics of mapping results, for example, the proportion of mapped reads, read distribution, etc. Not only can we learn the quality of mapping, but also we can make a reasonable explanation based on it for the subsequent abnormal analysis results. Therefore, we should check the quality of mapping results firstly to lay the foundation for the downstream analysis.

There are lots of tools that can be used to check the quality of mapping results, such as Samtools, [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [RSeQC](http://rseqc.sourceforge.net/), [Ribo-seQC](https://github.com/ohlerlab/RiboseQC). You will learn about comprehensive information of your mapping results with the help of these tools. Below, we are going to introduce how to check the quality of mapping results through these tools in the Galaxy.

Here, we use results produced by full data of samples to show information more visually.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# General statistics of alignment

## Mapping ratio

How many reads were mapped to the genome successfully? If there is a very low mapping ratio in one or more samples, we should think about the reason leading to this. Of course the following analysis should be suspended until the reason for the abnormal mapping result was found and solved.

> ### {% icon hands_on %} Hands-on: Aggregate the HISAT2 summary files with **MultiQC**
>
> 1. **MultiQC** {% icon tool %} with the following parameters:
>    - In *"Results"*
>      - {% icon param-select %} *"Which tool was used generate logs?"*: `HISAT2`
>      - {% icon param-collection %} *"Output of HISAT2"*: `Mapping summary` (output of **HISAT2** {% icon tool %})
> 2. Inspect the `Webpage` output from MultiQC
	{: .hands_on}

![Mapping results of Hisat2](../../images/quality-of-mapped-reads/hisat2_mapstats_mqc.png "Mapping results of Hisat2")

We can see that the number of successfully mapped reads is relatively high from the above figure. However, the percentage of multiple mapped reads is very high, and some even exceed 60%. Usually, such a strange proportion is quite abnormal, but this phenomenon is common in the Ribo-Seq data. Because reads from Ribo-Seq easily mapped to repeat regions on the genome due to the relatively short length. Besides, reads may contain the contaminant from the rRNA, elevating the ratio of unaligned.

## Other statistics of alignment

We can acquire more detailed information besides mapping ratio from the alignment through `samtools stats`. Then, we aggregate all results into one sheet using MultiQC.

> ### {% icon hands_on %} Hands-on: Calculate QC metric using samtools stats
>
> - **Samtools stats** {% icon tool %} with following parameters:
>   - {% icon param-collection %} *"BAM File"*: `aligned reads (BAM)`
>   - {% icon param-select %} *"Output"*: `One single summary file`
> - **MultiQC** {% icon tool %} with following parameters:
>   - In *"Results"*:
>     - {% icon param-select %} *"Which tool was used generate logs?"*: `Samtools`
>       - In *"Samtools output"*:
>         - {% icon param-select %} *"Type of Samtools output?"*: `stat`
>           - {% icon param-collection %} *"Samtools flagstat output"*: `Samtools stats on collection xxx`
>
{: .hands_on}

The picture below indicates a part of aggregated statistical information about `samtools stats`. Each blue point represents a sample, when the mouse arrow stay on it, the sample name and the number of corresponding mapped reads will be displayed in the top-left corner.

![Samtools general statistics](../../images/quality-of-mapped-reads/samtools_stats_1.png "Samtools general statistics")

![Samtools alignment metric](../../images/quality-of-mapped-reads/samtools_stats_2.png "Samtools alignment metric")

# Quality control using Ribo-seQC

The above is generally inspect items for quality control. However, due to the characteristic of ribosomes slide on the mRNA as 3nt in each step, certain quality control metrics other than above are also able to as the evaluation items for quality control of Ribo-Seq data.

> ### {% icon hands_on %} Hands-on: Check triplet nucleotide periodicity using Ribo-seqC
>
> - Because Ribo-SeqC needs twobit format of fasta, we should transfer fasta file to twobit file with faToTwoBit.
>
> > ### {% icon hands_on %} Hands-on: Transfer .fasta to .2bit
> >
> > - **faToTwoBit** {% icon tool %} with following parameters:
> >   - {% icon param-file %} *"fasta"*: `hg38_ucsc.fasta`
> >
> > > ### {% icon comment %} Comment
> > >
> > > If you want to know more about twobit file, you can read [twoBit](https://genome.ucsc.edu/goldenpath/help/twoBit.html).
> > >
> > {: .comment}
> >
> {: .hands_on}
>
> - **Ribo-seQC** {% icon tool %} with following parameters:
>   - In *"Prepare annotation files"*:
>     - In *"Inputs"*:
>       - {% icon param-file %} *"gtf"*: `gencode.v32.annotation.gtf`
>       - {% icon param-file %} *"fa in twobit format"*: `hg38_ucsc.2bit`
>       - {% icon param-file %} *"BSgenome"*: `Hsapiens.UCSC.hg38`
>   - In *"Main analysis"*:
>     - ~~***TODO*** {% icon param-file %} *"annotation"*: `build-in`~~
>     - {% icon param-collection %} *"BAM"*: `aligned reads (BAM)`
>
{: .hands_on}

You will obtain plenty of information about mapped results when the process above was completed. Hence, we introduce some parts of the results to display and check the mapped quality.

## Distribution of read lengths

Sequencing reads of Ribo-Seq data are from fragments of ribosome-enclosed, so the distribution of read lengths will mainly concentrate on a specific length.

![Distribution of read lengths](../../images/quality-of-mapped-reads/riboseqc1_read_length_distribution.png "Distribution of read lengths")

> ### {% icon question %} Questions
>
> Do you think the quality of read length distribution is good or bad?
>
> > ### {% icon solution %} Solution
> >
> > The quality of read length distribution is good. Because there is a significant peak at 32 nt and the distribution is ridge type.
> >
> {: .solution}
>
{: .question}

## Triplet nucleotide periodicity

As mentioned in the [Introduction to translatomics](), the most significant feature of Ribo-Seq data is triplet nucleotide periodicity, in addition, this feature is also the criterion to judge the quality of Ribo-Seq data. If we can't observe this feature, we should reflect on reasons leading to these results. For example, whether there is an error during the library preparation. If the triplet nucleotide periodicity can not be observed after we rule out all of the points that we may make an error, we should consider to drop out this data.

![Triplet nucleotide periodicity](../../images/quality-of-mapped-reads/riboseqc3_3nt.png "Triplet nucleotide periodicity")

We can observe triplet nucleotide periodicity from the figure above, so the quality of the sample data is not bad. Therefore, we can execute the subsequent analysis.


# Conclusion

{:.no_toc}

You will have a more comprehensive understanding of your data through this tutorial, and you will know how to check the quality of mapped results to find out the reason leading to abnormal indicators.
