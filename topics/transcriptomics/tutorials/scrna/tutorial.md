---
layout: tutorial_hands_on
topic_name: transcriptomics
tutorial_name: scrna
---

# Introduction
{:.no_toc}

The advent of single-cell RNA sequencing has provided the means to explore samples at the individual cell level, enabling a greater understanding of the development and function of such samples by the characteristics of their constituent cells. These cells can be classified into different cell subpopulations or "clusters" by the variation of gene expression in each cell, where cells in the same cluster exhibit the same levels of differential expression in the same set of related genes.

By identifying significant genes in each cluster, cell types can be inferred and similarity between clusters can be determined based on their proximity to one another. If the cells all come from the same sample tissue, the cell types will likely correspond to the different stages of cell differentiation expected of that tissue, wherein a lineage/heirarchy could potentially be derived.

This tutorial is in part inspired by aspects of the [Hemberg workflow](https://hemberg-lab.github.io/scRNA.seq.course/) at the Sangar insititute, as well as the [CGATOxford](https://github.com/CGATOxford/UMI-tools) workflow. The barcoding follows the CelSeq2 protocol and uses the lane configuration as that utilized by the Freiburg MPI Grün lab.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Background

Most scRNA sequencing techniques use pooled-sequencing approaches to generate a higher throughput of data by performing to amplification and sequencing upon multiple cells in the same "pool". From a bioinformatics standpoint, this means that the output FASTQ data from the sequencer is batch-specific and contains all the sequences from multiple cells.

The size of these FASTQ files are typically in the gigabyte range and are somewhat impractical for training purposes, so we will expediate the analysis by using a smaller subset of actual batch data. 

## CelSeq2

The CelSeq2 protocol is a paired-end protocol, meaning that the primers contain strand-specific information. In this case; Read1 contains the barcoding information followed by the polyT tail of the messenger RNA, and Read2 contains the actual sequence. Here, Read1 is regarded as the 'forward' strand and Read2 as the 'reverse' strand, though this is more a convention when dealing with paired-end data rather than an indication of the actual strand orientation.

 ![CelSeq2 Scheme](../../images/celseq2_schema.png "Read1 encapsulates the barcodes, and Read2 the mRNA sequence" )


## Barcodes

Barcodes are small random oligonucleotides that are inserted into the captured sequence at a specific point, and provide two pieces of information about the sequence:
  
 1. Which cell the sequence came from
 2. Which transcript the sequence came from

When the sequence is mapped against a reference genome, we can then see which gene locus it aligns to and qualitavely assert that, together with the two pieces of information above, the sequence depicts a transcript from a specific a gene that originated from a specific cell.

> ### {% icon question %} Question
>
> Why is it important to know which cell a sequence came from?
>
> > ### {% icon solution %} Solution
> >
> > If our sequence codes for a *GeneX* which is a gene of interest, we may want to know which cells express GeneX more than others.
> >
> > If CellA has 10 times more GeneX sequences than CellB, then we know that CellA and CellB differ at GeneX - which might suggest a causative source of variation for any change in function between CellA and CellB (or cells in the same cluster as CellA, and cells in the same cluster as CellB)
> >
> {: .solution}
{: .question}


> ### {% icon question %} Question
>
> Barcoding the cell makes sense, but why do we need to barcode the transcript too?  
> Can we not infer which gene the sequence originates from by simply mapping it against the reference genome?
> 
> > ### {% icon solution %} Solution
> >
> > *Yes* and *no*! 
> >
> > **Yes**: We can indeed align our sequence against a reference genome and obtain the name of the gene it aligns against. This sequence will then contribute to the 'count' of sequences that gene has, and increase the expression of that gene.
> > 
> > **No**: We do not know whether these 'counts' are *unique*. Many of these counts could be duplicates as a result of the amplification process. To explain further, we must look at UMIs and their role in the analysis.
> > 
> {: .solution}
{: .question}

> ### {% icon details %} UMIs: Mitigating duplicate transcript counts:
>
> One of the major issues with sequencing is that the read fragments require amplification before they can be sequenced. A gene with a single mRNA transcript will not be detected by most sequencers, so it needs to be duplicated 100-1000x times for the sequencer to 'see' it.
>
> Amplification is an imprecise process however, since some reads are amplified more than others, and subsequent amplification can lead to these over-amplified reads being over-amplified even more, leading to an exponential bias of some reads over others.
>
>  ![Amplification Bias](../../images/scrna_amplif_errors.png "A cell with two reads from different trascripts being amplified unevenly" )
>
> Consider the above example where two reads from different transcripts are amplified unevenly. The resulting frequency table would look like so:
>
>  |  | Reads in Cell 1 |
>  |--|------------------|
>  | Gene Red | 5 |
>  | Gene Blue | 0 |
>  
> But the truth is entirely different (i.e. Gene Red should have 1 count, and Gene Blue should also have 1 count).  
> How do we correct for this bias?
> 
> ### UMIs to the rescue
> 
> **Unique Molecular Identifiers** (or *UMIs*) constitute the second portion of a barcode, where their role is to *uniquely* count reads such that amplicons of the same read are only counted once, e.g:
> 
>  ![Amplification Bias with UMIs](../../images/scrna_amplif_errors_umis.png "A cell with four unique transcripts, two from Gene Red and two from Gene Blue")
> 
> Here, we see two unique transcripts from Gene Red and two unique transcripts from Gene Blue, each given a (coloured) UMI. After amplificaiton, Gene Red has more reads than Gene Blue. If we were to construct a frequency table as before to count the reads, we would have:
> 
>  |  | Reads in Cell 1 |
>  |--|-----------------|
>  | Gene Red | 6 |
>  | Gene Blue | 3 |
>  
>  Which is wrong, because it shows that Red has twice the expression that Blue does. However, we can reconstitute the true count by considering the UMI information:
>  
>  |  | UMI colour  | Reads in Cell 1 |
>  |--|-------------|-----------------|
>  | Gene Red | Pink | 2 |
>  |          | Blue | 4 |
>  | Gene Blue | Brown | 1 |
>  |           | Green | 2 |
>  
> We can then make the decision to ignore the frequencies of these UMIs, and simply count how many *unique* UMIs we see in each gene:
> 
>  |  | Set of UMIs in Gene | UMIs in Cell 1 |
>  |--|---------------------|----------------|
>  | Gene Red | {Pink, Blue} | 2 |
>  | Gene Blue | {Brown, Green} | 2 |
>  
>  Which provides us with the true count of the number of true transcripts for each gene:
>  
>  |  | UMIs in Cell 1 |
>  |--|----------------|
>  | Gene Red | 2 |
>  | Gene Blue | 2 |
{: .details}


> ### {% icon question %} Question
>
> 1. Are UMIs specific to genes? i.e. Can the same UMI map to different genes?
> 2. Can the same UMI map to different mRNA molecules of the same gene?
>
> > ### {% icon solution %} Solution
> >
> > 1. The same UMI can tag transcripts of different genes. UMIs are just 'added randomness' that help reduce amplification bias, but they are not unique to any particular gene.
> > 2. Yes, UMIs are not precise but work on a probabalistic level. In most cases, two transcripts of the same gene will be tagged by different UMIs. In rarer (but still prevalent) cases, the same UMI will capture different transcripts of the same gene.


are added to each cell
> >
> > If CellA has 10 times more GeneX sequences than CellB, then we know that CellA and CellB differ at GeneX - which might suggest a causative source of variation for any change in function between CellA and CellB (or cells in the same cluster as CellA, and cells in the same cluster as CellB)
> >
> {: .solution}
{: .question}



### Barcoding Format

We now know the role of UMIs and cell barcodes, but how do we handle them in the analysis?

Let us look at 4 example sequences:

<!-- 

These are reads that all map to ENSDARG00000019692. In [Cell, UMI] format:

(J00182:75:HTKJNBBXX:2:1114:12469:11073|J00182:75:HTKJNBBXX:2:2222:13301:35690|J00182:75:HTKJNBBXX:2:1203:25022:13763|J00182:75:HTKJNBBXX:2:1115:8501:46961)

            Cell  , UMI
1: 46961 -- ACCAGA, GGAAGA
2: 13763 -- GGTAAC, GTCCCA -> same umi, same cell
3: 35690 -- GGTAAC, GTCCCA -> same umi, same cell
4: 11073 -- GGTAAC, CGGCGT -> diff umi, same cell

-->

The Forward reads:

    @J00182:75:HTKJNBBXX:2:1115:8501:46961 1:N:0:ATCACG
    GGAAGAACCAGATTTTTTTTTTTTTTTTTT
    +
    AAFFFJJJJJJJFFFJJJJJJJJJJJJJJJ
    
    @J00182:75:HTKJNBBXX:2:1203:25022:13763 1:N:0:ATCACG
    GTCCCAGGTAACTTTTTTTTTTTTTTTTTT
    +
    AAFFFJJJJJJJJFFJJJJJJJJJFJ<FF-
    
    @J00182:75:HTKJNBBXX:2:2222:13301:35690 1:N:0:ATCACG
    GTCCCAGGTAACTTTTTTTTTTTTTTTTTT
    +
    AAFFFJJJJJJJ<AFJJJJJFFJJFJJJFF
    
    @J00182:75:HTKJNBBXX:2:1114:12469:11073 1:N:0:ATCACG
    CGGCGTGGTAACTTTTTTTTTTTTTTTTCC
    +
    AAFFFJJJJJJJFAFFJJJJJJJJF---<F

The Reverse reads:

    @J00182:75:HTKJNBBXX:2:1115:8501:46961 2:N:0:ATCACG
    GACCTCTGATCTTTACGAAAGGCCAACGCGTTTTCAGTCTGGACACGGTTCAGCTCCTGTTCATTATTCA
    +
    A<<A-777F<AA<AJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
    
    @J00182:75:HTKJNBBXX:2:1203:25022:13763 2:N:0:ATCACG
    GCCACCTAATTTCCGTCATCACACTCCTCTCCGTTTTCAACTTGCACAATGCTGTCTCCGCAGAATCCCT
    +
    ---<----<A---77-7A-FJ<JJFFJJ<JJAJ7<-FAFFJJFF<FFJJFFAJFA-AFFFJFFFFFJAJJ
    
    @J00182:75:HTKJNBBXX:2:2222:13301:35690 2:N:0:ATCACG
    CAATCCTCTCCGTTATCAACTTGCACAATGCTGTCTCCGCAGAATCCCTCCGGATCAGGATCGCTCTCCA
    +
    <<A-77--77F<----7AFJ-A--FJJJFAJF-AFAJAJ<JFJ<JJJFFJJJFJJJJJAAFJJJFJJJF-
    
    @J00182:75:HTKJNBBXX:2:1114:12469:11073 2:N:0:ATCACG
    ATCCACTTATTGCAAAGCAGAGGACATTGAGTCTCACCTTTTGTCCAGGTCTTCCAATTTCACCCTGCAA
    +
    A-77AA-7FF<7FFJFFFJJJJJJJJJJJJJ-AFJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJ

As shown in image 1, the barcode is on the forward read and sequence is on the reverse read. The barcode consists of the 6bp UMI followed by the 6bp cell barcode (followed by the poly-T tail). The sequence is 70 bp long as contains the information that we wish to map. 

Mini Exercise:
Q: How do we know the barcoding format of reads?
A: Usually you are explicitly told, as well being given a list of cell barcodes to use. Most cell barcodes have an edit distance >= 2 between them in order to correct against 1bp read errors.
Q: How do we verify the barcode format of our reads?
A: Look at the FASTQC
 We can see that the distribution of the first 6bp has far more variability than the following 6bp which seems to have less variation. From this we can deduce where the boundaries of our barcodes are.
 SQ: Which side is the cell barcode? Which side the UMI?
 We expect to see more genes than cells, with an absolute maximum of 4**6 UMIs in a given cell. So the left side is the UMIs since there appears to be more of an even distribution, and the right side is the CB since there the variation seems more spiky.
 

### Uniting Barcodes with Sequence

In a sense, we have a disparity in our data: the reverse reads contain the sequences we wish to map, but not the barcodes; the forward reads contain the barcode, but not the sequence. For the forward and reverse reads given above, the information that we really want from both can be summarized in this table:

| Read | Cell | UMI | Sequence |
|------|------|-----|----------|
| @J00182:75:HTKJNBBXX:2:1115:8501:46961  | ACCAGA | GGAAGA | GACCTCTGATCTTTACGAAAGGCCAACGCGTTTTCAGTCTGGACACGGTTCAGCTCCTGTTCATTATTCA |
| @J00182:75:HTKJNBBXX:2:1203:25022:13763 | GGTAAC | GTCCCA | GCCACCTAATTTCCGTCATCACACTCCTCTCCGTTTTCAACTTGCACAATGCTGTCTCCGCAGAATCCCT |
| @J00182:75:HTKJNBBXX:2:2222:13301:35690 | GGTAAC | GTCCCA | CAATCCTCTCCGTTATCAACTTGCACAATGCTGTCTCCGCAGAATCCCTCCGGATCAGGATCGCTCTCCA |
| @J00182:75:HTKJNBBXX:2:1114:12469:11073 | GGTAAC | CGGCGT | ATCCACTTATTGCAAAGCAGAGGACATTGAGTCTCACCTTTTGTCCAGGTCTTCCAATTTCACCCTGCAA |

Mini Exercise:
Q: Which of these reads are duplicates?
Q: Which of these reads come from the same cell?
Q: How many cells in total are depicted here? How many UMIs?

Exercise: These four reads all actually map to same gene. Which gene is it?
 - Take 1 sequence, use BLAT


Q: How should we couple these two source of information into a single location without impacting the data content?

A: We take the barcode information from the forward reads, and stick it into the *header* of the reverse reads. That way we can align our sequence to the reference and still keep the barcode information assosciated with the reads.

i.e. 

    @J00182:75:HTKJNBBXX:2:1115:8501:46961 2:N:0:ATCACG ACCAGA_GGAAGA
    GACCTCTGATCTTTACGAAAGGCCAACGCGTTTTCAGTCTGGACACGGTTCAGCTCCTGTTCATTATTCA
    +
    A<<A-777F<AA<AJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
    
    @J00182:75:HTKJNBBXX:2:1203:25022:13763 2:N:0:ATCACG GGTAAC_GTCCCA
    GCCACCTAATTTCCGTCATCACACTCCTCTCCGTTTTCAACTTGCACAATGCTGTCTCCGCAGAATCCCT
    +
    ---<----<A---77-7A-FJ<JJFFJJ<JJAJ7<-FAFFJJFF<FFJJFFAJFA-AFFFJFFFFFJAJJ
    
    @J00182:75:HTKJNBBXX:2:2222:13301:35690 2:N:0:ATCACG GGTAAC_GTCCCA
    CAATCCTCTCCGTTATCAACTTGCACAATGCTGTCTCCGCAGAATCCCTCCGGATCAGGATCGCTCTCCA
    +
    <<A-77--77F<----7AFJ-A--FJJJFAJF-AFAJAJ<JFJ<JJJFFJJJFJJJJJAAFJJJFJJJF-
    
    @J00182:75:HTKJNBBXX:2:1114:12469:11073 2:N:0:ATCACG GGTAAC_CGGCGT
    ATCCACTTATTGCAAAGCAGAGGACATTGAGTCTCACCTTTTGTCCAGGTCTTCCAATTTCACCCTGCAA
    +
    A-77AA-7FF<7FFJFFFJJJJJJJJJJJJJ-AFJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJ

where we have added CellBarcode_UMI to the headers of each read. We can now throw away our forward reads, as they have no more useful information within them and proceed to the mapping stage.


# Barcode Extraction

TODO!

# Mapping

Mapping is a relatively straightforward process:

 1. Select your genome
 2. Select your GTF file
 3. Run MultiQC on the resulting STAR log

The FASTQ data was generated from sequencing Zebra working with Zebrafish data, so to perform the alignment we will need to gather all data relevant to that genome. We will use the latest version (DanRerv10).

> ### {% icon Hands-on %} Performing the Alignment
>
> 1. Obtain the GTF file and import it into our history.
>   - This file contains all the gene, exon, intron, and other regions of interest that we will use to annotate our reads, should our reads align to any of the regions specified in this file.
>   - *"Shared Data"* → *"Data Libraries"* → *"Genomes + Annotations"* → *"Annotations"* → "Dario_Rerio_v10.gtf (danRer10)" → Click
>   - *"to History"*
> 
> 2. Select **RNA-STAR** {%icon tool %} with the following parameters:
>   - *"Single End"*: Our input FASTQ has sequencing data from only a single read
>   - *"RNA-Seq FASTQ/FASTA file"*: Of the FASTQ files output by UMI_tools extract, select the one with the sequencing reads intact.
>   - *"Custom or built-in reference"*: `Use a built-in index`
>   - *"Reference genome with or without an annotation"*: `use genome reference without builtin gene-model`
>   - *"Select reference genome"*: The FASTQ data was generated from sequencing zebrafish, so select `DanRer10`
>   - *"Gene model (gff3,gtf) file for splice junctions"*: Select the GTF file we imported into our history
>   - Execute

This should take a minute or two depending on your position in the queue. Once your output files are green, proceed to the next step.

> ### {% icon Hands-on %} Performing the QC on the Alignment
> 
> Let us examine how well our alignment went.
> 
> 1. Select **MultiQC** {%icon tool %} with the following parameters:
>  - *"Results"* → *"1: Results"* → *"Which tool was used to generate logs?"*: `STAR`
>  - *"STAR output"* → *"1: STAR output"* → *"Type of STAR output?"*: `Log`
>  - *"STAR log output"* : Select the file that ends in " log"
>  - Execute
>
> 2. Once green, click on the "MultiQC on data : Webpage" eye symbol.
>
> > ### {% icon question %} Question
> > 
> > 1. What percentage of our reads are uniquely mapped? How many millions of reads is this percentage?
> > 2. What percentage of our reads are mapped to more than one locus?
> > 3. Is our overall mapping 'good' ?
> >
> >
> > > ### {% icon solution %} Solution
> > > 
> > > 1. 73.5% or 8.3 million reads were successfully mapped
> > > 2. 11.3% are multiply mapped, and 2.2% were mapped to too many loci
> > >   - Multiply mapped means that a read was aligned to more than one gene
> > >   - Mapped to too many loci means that a read was aligned to 10 or more loci, and should be ignored.
> > > 3. It depends on how good the sequencing protocol is, and how many reads in total were mapped.
> > >   - 90% is amazing, reserved for bulk RNA-seq which typically has high coverage
> > >   - 70% is weak for bulk RNA-seq, but good for single-cell RNA-seq
> > >   - 6 million mapped reads should be enough to generate a downstream analysis from.
> > >
> > {: .solution}
> {: .question}
{: .hands_on}

# Quantification

> ### {% icon comment %} Recap of previous stages
>
> 1. *Barcode Extraction*:
>  Here we used `umi_tools extract` on our input forward and reverse FASTQ files, and extracted the umi and cell barcode from the forward read *sequence*, and placed it into the *header* of both forward and reverse FASTQ files. i.e. FASTQ files → Modified FASTQ files
> 2. *Mapping*
>  We took the sequencing data from the reverse FASTQ file (with modified headers) and aligned it to the Zebrafish genome, using annotations presented in the GTF file for that genome. i.e. Modified FASTQ file (reverse) → BAM file
>

We now have a BAM file of our aligned reads, with cell and UMI barcodes embedded in the read headers. We also have the chromosome and base-pair positions of where these reads are aligned.

> ### {% icon Hands-on %} We can confirm this by peeking into the BAM file
>
> 1. Click on the eye symbol of the BAM output from STAR.
> 2. There are many header lines that begin with '@' which we are not interested in. 
> 3. Do a Ctrl+F search for `@co` and then look at the lines directly below it.
> 
> One such read is given as so:
>
>    J00182:75:HTKJNBBXX:2:1121:9729:45889_GACGAA_GTGGTC	16	chr1	2030	3	70M	*	0	0	AGAGGTTCCAATATTCCCATGAAATTGAGATTTTGTAAAAGAGTGAAGTGTGGTTACTTTCACTGAGAGG	JJJJJJJJJJJJJJJJJJJJJJJJFJJJJJAJJJJJJJJJJFJFJFFJJJJJJJJJJJJFF7AJA-77<A	NH:i:2 HI:i:1 AS:i:64 nM:i:2
>   
> 

The fields of the BAM file can be better explained at section 1.4 of [the SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf), but we will summarize the main fields of interest here:

 * `J00182..._GACGAA_GCGGTC`: The *readname* appended by an underscore '_', the cell barcode, another '_', and then the UMI barcode.
 * `16`: The FLAG value, which can be explained in the SAM specification, or more interactively [here](https://broadinstitute.github.io/picard/explain-flags.html).
> ### {% icon question %} What does the FLAG value of 16 tell us about this read?
> > ### {% icon solution %} Solution
>
> The read aligns to the reverse strand
 * `chr1` `2030`: The position and base-pair of alignment of the first base of the sequence.
 * We next have a series of quality fields, as well as the `sequence



In a sense, we can say that the forward read is effectively not useful at all for mapping, it only holds information as to where each read came from via the barcodes. How can we extract this information and 



We can decode this for each of our 4 reads


As before, we 




The read primers are added to each cell before amplification and sequencing

The use of 


When 




FASTQ data output from scRNA sequencers are batch-specific, meaning that sequences from individual cells are not demultiplexed into individual FASTQ files, but are

Introduction about this part

## Subpart 1

Short introduction about this subpart.

<!--
{% icon hands_on %} will render the hands_on icon as specified in
_config.yml in the root of this repository.
-->

> ### {% icon hands_on %} Hands-on: Data upload and organization
>
> 1. Step1
> 2. Step2
>
>    > ### {% icon comment %} Comments
>    > A comment
>    {: .comment}
>
>    > ### {% icon tip %}Tip: A tip
>    >
>    > * Step1
>    > * Step2
>    {: .tip}
{: .hands_on}

## Subpart 2

Short introduction about this subpart.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Step1
> 2. Step2
>
>    > ### {% icon question %} Question
>    >
>    > Question?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > Answer to question
>    > >
>    > {: .solution}
>    >
>    {: .question}
{: .hands_on}

Some blabla
> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Step1
> 2. **My Tool** {% icon tool %} with the following parameters
>   - *"param1"*: the file `myfile`
>   - *"param2"*: `42`
>   - *"param3"*: `Yes`
>
> 3. **My Tool** {% icon tool %} with the following parameters
>   - {% icon param-text %} *"My text parameter"*: `my value`
>   - {% icon param-file %} *"My input file"*: `my file`
>   - {% icon param-files %} *"My multiple file input or collection"*: `my collection`
>   - {% icon param-select %} *"My select menu"*: `my choice`
>   - {% icon param-check %} *"My check box"*: `yes`
>
>    > ### {% icon question %} Questions
>    >
>    > 1. Question1?
>    > 2. Question2?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > 1. Answer for question1
>    > > 2. Answer for question2
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
> 3. Step3
{: .hands_on}

# Part 2

Short introduction about this subpart.

> ### {% icon comment %} Comment
>
> Do you want to learn more about the principles behind mapping? Follow our [training](../../NGS-mapping)
{: .comment}


> ### {% icon details %} More details on the ....
>
> Add more details in Markdown. By default the box is collapsed. And is expanded when clicked
>
{: .details}

# Conclusion
{:.no_toc}

Conclusion about the technical key points. And then relation between the techniques and the biological question to end with a global view.
