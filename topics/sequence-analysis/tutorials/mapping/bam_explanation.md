A BAM ([Binary Alignment Map](https://en.wikipedia.org/wiki/SAM_(file_format))) file is a compressed binary file storing the read sequences, whether they have been aligned to a reference sequence (e.g. a chromosome), and if so, the position on the reference sequence at which they have been aligned.

> <hands-on-title>Inspect a BAM/SAM file</hands-on-title>
>
> 1. Inspect the {% icon param-file %} output of **{{ include.mapper }}** {% icon tool %}
>
{: .hands_on}

A BAM file (or a SAM file, the non-compressed version) consists of:

- A header section (the lines starting with `@`) containing metadata particularly the chromosome names and lengths (lines starting with the `@SQ` symbol)
- An alignment section consisting of a table with 11 mandatory fields, as well as a variable number of optional fields:

    Col | Field | Type | Brief Description
    --- | --- | --- | ---
    1 | QNAME | String | Query template NAME
    2 | FLAG | Integer | Bitwise FLAG
    3 | RNAME | String | References sequence NAME
    4 | POS | Integer | 1- based leftmost mapping POSition
    5 | MAPQ | Integer | MAPping Quality
    6 | CIGAR | String | CIGAR String
    7 | RNEXT | String | Ref. name of the mate/next read
    8 | PNEXT | Integer | Position of the mate/next read
    9 | TLEN | Integer | Observed Template LENgth
    10 | SEQ | String | Segment SEQuence
    11 | QUAL | String | ASCII of Phred-scaled base QUALity+33 

> <question-title></question-title>
>
> 1. Which information do you find in a SAM/BAM file?
> 2. What is the additional information compared to a FASTQ file?
>
> > <solution-title></solution-title>
> > 1. Sequences and quality information, like a FASTQ
> > 2. Mapping information, Location of the read on the chromosome, Mapping quality, etc
> {: .solution }
{: .question}
