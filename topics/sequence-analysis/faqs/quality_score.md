---
title: Quality Scores
area: format
box_type: details
layout: faq
contributors: [bebatut, nakucher, hexylena]
---

But what does this quality score mean?

The quality score for each sequence is a string of characters, one for each base of the nucleotide sequence, used to characterize the probability of misidentification of each base. The score is encoded using the ASCII character table (with [some historical differences](https://en.wikipedia.org/wiki/FASTQ_format#Encoding)):

To save space, the sequencer records an [ASCII character](http://drive5.com/usearch/manual/quality_score.html) to represent scores 0-42. For example 10 corresponds to "+" and 40 corresponds to "I". FastQC knows how to translate this. This is often called "Phred" scoring.

![Encoding of the quality score with ASCII characters for different Phred encoding. The ascii code sequence is shown at the top with symbols for 33 to 64, upper case letters, more symbols, and then lowercase letters. Sanger maps from 33 to 73 while solexa is shifted, starting at 59 and going to 104. Illumina 1.3 starts at 54 and goes to 104, Illumina 1.5 is shifted three scores to the right but still ends at 104. Illumina 1.8+ goes back to the Sanger except one single score wider. Illumina]({{site.baseurl}}/topics/sequence-analysis/faqs/images/fastq-quality-encoding.png)

So there is an ASCII character associated with each nucleotide, representing its [Phred quality score](https://en.wikipedia.org/wiki/Phred_quality_score), the probability of an incorrect base call:

Phred Quality Score | Probability of incorrect base call | Base call accuracy
--- | --- | ---
10 | 1 in 10 | 90%
20 | 1 in 100 | 99%
30 | 1 in 1000 | 99.9%
40 | 1 in 10,000 | 99.99%
50 | 1 in 100,000 | 99.999%
60 | 1 in 1,000,000 | 99.9999%


What does 0-42 represent? These numbers, when plugged into a formula, tell us the probability of an error for that base. This is the formula, where Q is our quality score (0-42) and P is the probability of an error:

```
Q = -10 log10(P)
```

Using this formula, we can calculate that a quality score of 40 means only 0.00010 probability of an error!
