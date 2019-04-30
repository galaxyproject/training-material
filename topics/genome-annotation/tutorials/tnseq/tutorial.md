---
layout: tutorial_hands_on

title: "Essential genes detection with Transposon insertion sequencing"
zenodo_link: "https://doi.org/10.5281/zenodo.940733"
tags:
  - bacteria
  - tnseq
  - essential genes
questions:
  - "What is Transposon insertion Sequencing?"
  - "How to get TA sites Coverage ? "
  - "How to predict essential genes ?"
objectives:
  - "Understand the read structure of TnSeq analyses"
  - "Predict Essential genes with Transit"
  - "Compare gene essentiality in control sample and an different experimental condition"
time_estimation: "7H"
key_points:
contributors:
  - delphine-l
---

# Introduction
{:.no_toc}

In microbiology, identifying links between genotype and phenotype is key to understand bacteria growth and virulence mechanisms, and to identify targets for drugs and vaccines. These analysis are limitated by the lack of bacterial genome annotations (eg 30% of genes for S. pneumoniae are of unknown function) and by the fact that genotypes often arose from complex composant interactions.


## Overview of Transposon insertion Sequencing
{:.no_toc}

Transposon insertion sequencing is a technique used to functionally annotate bacterial genomes. The genome is saturated by transposon insertions, and the insertion of a transposon being disruptive for the region, the analysis of insertion frequency provides information on how the bacteria fitness change due to this disruption (see [Transposon insertion sequencing method](#figure-1)) :
 - An insertion mutant with lower fitness decrease frequency in the population
 - An increased fitness lead to increased frequency in the population




 ![Illustration of tnseq Method](../../images/tnseq/principle_tnseq.png "<b>Transposon insertion sequencing method</b> - <b>a. Data production</b> The initial population genomes are mutated so that the genome is saturated with transposon insertions.  A library is <i>saturated</i> if in the genomes across the whole population of bacteria, each potential insertion site has at least one insertion. The population is then divided into several media containing different growth conditions. After growth, the regions flanking the insertion are amplified and sequenced, allowing to determine the location of the insertion. <b>b. Analysis</b> After alignement to the reference genome, the resulting data will show a discrete repartition of reads on each TA site. If a gene present several insertions, like the two leftmost genes in <i>Condition A</i>, it means that its disruption has little or no impact to the bacterial growth. On the other hand, when a gene shows no insertions at all, like the rightmost gene in <i>Condition A</i>, is means that any disruption in this gene killed the bacteria, meaning its a gene essential to bacteria survival. If the library is sufficiently saturated, there is a clear threshold between essential and non-essential genes when you analyze the insertion rate per gene. (From <a href='#Chao2016'> Chao <i>et al.</i> 2016</a>)")


Two type of transposon insertion methods exist:
 -  Gene disruption, where we analyze only the disruptions. (The object of this tutorial)
 -  Regulatory element insertion, where different promoters are inserted by the transposon, and we analyze the change in gene expression in addition to the disruption. (Will be the subject or another tutorial)

## <a name="BuildLibrary">Building a TnSeq library</a>
{:.no_toc}

In this tutorial, we are using mariner transposon, that target TA sequences, in ordered to target the whole genome uniformely. Different types of transposon can be used depending of the goal of you analysis :
- Randomly pooled tranposon :
    - Mariner-based transposons, which target TA dinucleotides
        - The TA are distributed relatively evenly along genome, which allows to impact statistically every gene
        - There is in average more than 30 insertions site per kb
        - The local variations means less loci and less statistical power
        - Advantages: low insertion bias, easy to build saturated libraries
    - Tn5-based vectors, which insert at random sites
        - Require no target sequences
        - It has a preference for high GC content, causing insertion bias
        - Useful for specie where it is difficult to build mariner based transposons
- Defined sequence transposon
    - Can be used to study interactions in pathways of interest
    - More precise targeting (small genes, pathways) for specific analyses

Independently of your choice of transposon you need to be careful about your library complexity . A large complexity means that there is multiple insertion in every potential locus. The higher density of insertion you have, the greater precision  you have in identifying limits of regions of interest. If the density of the library is too low, some genes might not be disrupted by chance and mistaken for essential. The advantage of a target specific transposon, like the mariner, in opposition of a Tn5-based transposon inserting randomly, it that the limited number of insertion sites makes it easier to build high complexity libraries.

After you selected the type of transposon corresponding to your goals, you need to modify it to allow insertion site amplification and sequencing so that you get a library representative of the tranposon insertion.
Biases could get introduced by the process due to uneven fragment sizes. To avoid this problem, we can introduce a Type I restriction site to cleave DNA downstream of transposon, and get uniform fragment sizes

It has been shown that a minimum length of 16 bp is necessary for precise mapping on the genome {% cite Kwon2015 %} . We can therefore use the MmeI restriction site (21pb) but not BsmFI (11 to 12 bp). It is not important in that case to have longer reads as we do not care to have a good coverage on the entire genome, the only information we need is the Ta site affected by an insertion, to do that we only need the location of the start of the reads.

In this tutorial, the transposon is a mariner Himar1 with the structure described in the figure [Structure of the tranposon constructs](#TranspStructure) {% cite Santiago2015 %}.

![Structure of the transposon containing several parcodes and adapters](../../images/tnseq/tranposon_structure.png "<b>Structure of the tranposon constructs</b> - The transposon construct is a mariner transposon with two specific region used to specifically sequence the region upstream of the insertion. The transposon inserts at TA site at the ITR junctions. These ITR junctions have been modified to include a Mme1 restriction site. Using MmeI enzyme to determine the size of the reads allow to have a homogeneous read size and therefore avoid a bias in the representation of the insertions. It also includes a NotI restriction site (cut 21 bp upstream from the restriction site). These two site are the 5' and 3' limits to the genomic DNA we want to sequence. <b>A. Sequence flanking genomic regions</b> After digestion by NotI restriction enzyme, the fragments are attached to biotinylated adaptors that link to NotI restriction site. The attached fragment are then digested by MMeI at a site upstream , where an Illumina primer is then linked. The sequencing is then done, adding Illumina adaptors and an additional barcode to the read for multiplexed sequencing. <b>B. Removing incorrect fragments</b> An insertion can sometimes be composed of one or more copies of the transposon (multimer). There is therefore a risk to select plasmid backbone sequence. To solve this problem, an additional NotI has been add in the backbone to create different length construct, that can later be filtrated (. Different promoters are added to the construct along with an additional 3 bp barcode to analyze differential expression impact, but this will be the subject of another tutorial. (From <a href='#Santiago2015'> Santiago <i>et al.</i> 2015</a>)")


Because of this complex tranposon structure, the reads obtained after sequencing contain a large portion of tranposon sequence for a 16-17 bp genomic sequence. This will necessitate several step of pre-processing to extract this genomic sequence.




## Tnseq analysis
{:.no_toc}

Once we extracted the genomic sequences from the initial reads, we need to locate each of them on the genome to link them to a TA site. To do that we need to map them to a reference genome, link them to a specific insertion site, and then count the number of insertion for each TA site.

Once we have the count of insertion at every insertion site, there is several methods existing to identify essential genes of regions. They can be divided in two major categories (see figure [Methods of TnSeq Analyses](#AnalysesMethods) {% cite Chao2016 %}) :
1. Annotation dependent : The read counts and/or insertion frequency are calculated across defined regions (genes, promoters ... )
2. Annotation independent : The read counts and/or disrupted sites are considered across the whole genome, independently of defined structures.




![Different types of TnSeq Analyses](../../images/tnseq/type_of_analyses.png "<b>Methods of TnSeq Analyses</b> - <b>Annotation dependent method</b> The total read count an/or percentage of disrupted site are computed per annotated regions. The values are then compared to the rest of the genome to classify the genes into the categories <i>essential</i> or <i>non-essential</i>. <b>Annotation independent method</b> The total read count and/or disrupted sites are computed independently of annotated regions. One of these methods is using a sliding window. Each window is then classified into the categories <i>essential</i> or <i>non-essential</i>. After the windows have been classified, they are linked annotations, and the genes/regions can be classified as <i>essential</i>, <i>non-essential</i>, or <i>domain essential</i> according to the classification of the windows they cover. The same classification can be done using HMM based methods instead of sliding windows. In that case, each insertion site will be predicted as <i>essential</i> or <i>non essential</i>. (From <a href='#Chao2016'> Chao <i>et al.</i> 2016</a>)")

The objectives of this tutorial will then be to remove non genomic sequences from the reads, align them to a reference genome, and use the location of genes to determine a list of essential genes.
> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

Let's start with uploading the data.


> ### {% icon hands_on %} Hands-on: Import the data
>
> 1. Create a new history for this tutorial and give it a proper name
>
>    {% include snippets/create_new_history.md %}
>    {% include snippets/rename_history.md %}
>
> 2. Import `Tnseq-Tutorial-reads.fastqsanger.gz`, `staph_aur.fasta`, `construct_barcodes.fasta`, `condition_barcodes.fasta`, and `staph_aur.gff3` from [Zenodo](https://zenodo.org/record/2579335)
>
>
>    {% include snippets/import_via_link.md %}
>
>    As default, Galaxy takes the link as name, so rename them.
>
> 4. Rename the files ...
>
>    {% include snippets/rename_dataset.md %}
>
> 3. Inspect the files by clicking on the {% icon galaxy-eye %} (eye) icon (**View data**)
>
{: .hands_on}


# Removing all non genomic sequences from the sequenced reads

## Data Structure

The experimental design of transposon insertion sequencing produces raw reads containing a lot of adapters and foreign sequences that has been used to insert and target the transposon. In order to obtain the core reads that contain only genomic sequence, we have a number of steps to do to remove them and divide the reads per experimental condition and type of transposon.

![Workflow of data pre-processing](../../images/tnseq/preprocessing.png "<b>Data pre-processing</b> - The pre-processing of the data will be done through several steps of Cutadapt software, first we sill separate the reads of each experimental condition based on a 8 bp barcode at the beginning of the read (Illumina demultiplexing). The tail of each set of read is then removed. It immediately follows the 3bp barcode specific to transposon constructs, and contains illumina adapter sequence and downstream. To be sure all our reads have been trimmed correctly we filter out the reads too large. We then separate the reads per transposon construct and then remove the remaining transposon sequence containing MmeI.")


## Separating reads from different experimental conditions

 First we divide the initial data set by experimental conditions thanks to a 8bp barcode added by the Illumina multiplexing protocol.
 Take a look at the fasta file containing the barcodes for each condition (`condition_barcodes.fasta`) :


 Barcode data:
 ```
 >control
 ^CTCAGAAG
 >condition
 ^GACGTCAT
 ```

 The '^' at the beginning of the sequence means we want to anchor the barcode at the beginning of the read. To know more about the symbols used by cutadapt, see cutadapt [manual](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types)

> ### {% icon hands_on %} Hands-on:  Barcode splitting with Cutadapt
>
> 1. **Cutadapt** {% icon tool %} Split the dataset based on barcodes with:
>    - Set *"Single-end or Paired-end reads?"* to `Single-end`
>    - Set *"FASTQ/A file"* to the `fastq.gz file` containing the training set of reads.
>    - **Click** on `Insert 5' (Front) Adapters`
>         - Set *"Source"* to `File From History`
>         - Set *"Choose file containing 5' adapters"* to the `condition barcodes` file in the history
>     ![Use cutadapt with front adapters](../../images/tnseq/frontadaptercutadapt.png)
>
>    - **Click** on `Adapter Options`
>         - Set *"Maximum error rate"* to `0.15` to allow 1 mismatch
>         - Set *"Match times"* to `3` in case the barcode attached several times
>     ![Adapter options in Cutadapt](../../images/tnseq/optionaptercutadapt.png)
>    - **Click** on `Output Options`
>         - Set *"Report"* to `yes`
>         - Set *"Multiple output"* to `yes` to separate the reads into one file per condition
>     ![Output options in Cutadapt](../../images/tnseq/cutadaptOutput.png)
>    - **Click** on `Execute`
>
>
>    > ### {% icon question %} Questions
>    >
>    >  What would change if our barcodes were at the end of the reads ?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > >  We would have used the option `Insert 3' (End) Adapters` and anchored them at the end of the read with the symbol `$` in our fasta file containing the barcodes.
>    > >
>    > {: .solution }
>    >
> {: .question}
{: .hands_on}

The output is a collection of the different conditions dataset, here control and condition, and a report text file.
You can look at the report and see that 100% of the reads has been trimmed.


## Removing Adapter sequence

Now that we have divided the set per condition, we are going to trim the tail containing the illumina adapter. We want to remove the adapter and everything downstream, we are therefore gonna use the end adapter option of cutadapt and not anchor the sequence anywhere.
To remove the reads that might not have been trimmed because of too many mismatches or other reasons, we will then filter the reads by size. Because we know the approximate size of the remaining sequences, we can filter the reads based on this estimation.


> ### {% icon hands_on %} Hands-on:  Remove Adapter with Cutadapt
>
> {% icon tool %} Select the **Cutadapt** tool in the tool bar and run with the following parameters to remove adapters:
>    - Set *"Single-end or Paired-end reads?"* to `Single-end`
>    - Set *"FASTQ/A file"* to the `Cutadapt on data... Output` collection output of the previous step.
>    - **Click** on `Insert 3' (End) Adapters`
>        - Set *"Source"* to `Enter custom Sequence`
>        - Set *"Enter custom 3' adapter sequence"* to `CGTTATGGCACGC`
>    - **Click** on `Adapter Options`
>        - Set *"Match times"* to `3` in case the barcode attached several times
>    - **Click** on `Output Options`
>        - Set *"Report"* to `yes`
>    - **Click** on `Execute`
>
>
> {% icon tool %} Run  **Cutadapt**  again with the following parameters to filter reads based on length:
>    - Set *"Single-end or Paired-end reads?"* to `Single-end`
>    - Set *"FASTQ/A file"* to the `Cutadapt on data... Output` collection output of the previous step.
>    - **Click** on `Filter Options`
>        - Set *"Minimum length"* to `64`
>        - Set *"Maximum length"* to `70`
>    - **Click** on `Output Options`
>        - Set *"Report"* to `yes`
>    - **Click** on `Execute`
>
>
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What are the outputs at this step?
>    >
>    > 2. What percentage of the reads contained the adapter?
>    >
>    > 3. How many reads where discarded after filtering?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > 1. The outputs are two collections : one containing the reads in two conditions, and one containing the Cutadapt reports for each condition.
>    > >
>    > > 2. More than 99% of the reads contained the adapter in both conditions. (Line `Reads with adapters` in Report output)
>    > >
>    > > 3. Less than 2% of the reads were discarded in both conditions. (Lines `Reads that were too short` and `Reads that were too long` in Report output)
>    > >
>    > {: .solution }
> {: .question}
>
{: .hands_on}

We can see that is both samples the reads have pass the filtering at more than 98%. If the reads passing is very low, it means the previous trimming steps is incomplete or faulty.

## Separate reads from different transposon constructs

We have now removed the transposon sequences that was outside of the 3 bp barcode specific to the type of construct. The constructs used in this experiment contain different strengths and directions of promoters. It means that in addition of disrupting a gene at the location of the insertion, it will modify the expression of either upstream or downstream regions. The analysis of this modification will be studied in another training material, but for now we will consider it does not impact the essentiality analysis and we will use the different constructs as replicates.
We therefore need to separate the reads based on the construct specific 3bp barcodes.


> ### {% icon question %} Questions
>
> Take a look at the fasta file containing the construct barcodes. What does the "$" means in the barcode sequence file ?
>
>   > ### {% icon solution %} Solution
>   >
>   > It means the barcode is anchored at the end of the reads.
>   >
> {: .solution }
>
{: .question}


> ### {% icon hands_on %} Hands-on:  Barcode split with Cutadapt
>
> {% icon tool %} Select the **Cutadapt** tool in the tool bar and run with the following parameters:
>    - Set *"Single-end or Paired-end reads?"* to `Single-end`
>    - Set *"FASTQ/A file"* to the  `Cutadapt on data... Output` collection output of the previous step.
>    - **Click** on `Insert 3' (End) Adapters`
>        - Set *"Source"* to `File From History`
>        - Set *"Choose file containing 3' adapters"* to the `construct barcodes` file in our history
>    - **Click** on `Adapter Options`
>        - Set *"Match times"* to `3` in case the barcode attached several times
>    - **Click** on `Output Options`
>        - Set *"Report"* to `yes`
>        - Set *"Multiple output"* to `yes` to separate the reads into one file per condition
>    - **Click** on `Execute`
>
>    > ### {% icon question %} Questions
>    >
>    > Are the reads equally divided between constructs ?
>    >
>    > > ### {% icon solution %} Solution
>    > >  When you look at the reports, you can see that most of the reads have been assigned to the *blunt* construct, this is because the blunt construct is the control and does not contain any promoters. This means that there is less negative selective pressure on blunt than the other ones, that have affected flanking region in addition to the disrupted gene at the insertion site.  This won't be a problem here as the tool we use for the essentiality prediction consider the sum of reads in the replicates.
>    > >
>    > {: .solution }
> {: .question}
>
{: .hands_on}

You can notice that the output of this split is a *nested collection*, a collection of collection.

## Removing remaining transposon sequence.

The last remaining transposon sequence is the linker containing the MmeI restriction site.


> ### {% icon hands_on %} Hands-on:  Remove Linker with Cutadapt
>
> {% icon tool %} Select the **Cutadapt** tool in the tool bar and run with the following parameters:
>    - Set *"Single-end or Paired-end reads?"* to `Single-end`
>    - Set *"FASTQ/A file"* to the `Cutadapt on data... Output` nested collection output of the previous step.
>    - **Click** on `Insert 3' (End) Adapters`
>        - Set *"Source"* to `Enter custom Sequence`
>        - Set *"Enter custom 3' adapter sequence"* to `ACAGGTTGGATGATAAGTCCCCGGTCTATATTGAGAGTAACTACATTT`
>    - **Click** on `Adapter Options`
>        - Set *"Maximum error rate"* to `0.15`
>    - **Click** on `Output Options`
>        - Set *"Report"* to `yes`
>    - **Click** on `Execute`
>
{: .hands_on}

Verify that the majority of the read have been trimmed. Now that we isolated the genomic sequences from the initial reads. We want to align them to count how many insertion have been retained at each TA sites.

# Counting the number of insertion per TA sites

## Aligning the reads to a reference genome

The first step is to map our read to the reference genome. We are going to use the tool Bowtie. We could also use Bowtie2 with an end-to-end option, but Bowtie is more suitable for very short reads like ours (16-17 bp).

> ### {% icon hands_on %} Hands-on:  Map reads with Bowtie
>
> {% icon tool %} Select the **Map with Bowtie for Illumina** tool in the tool bar and run with the following parameters:
>   - Set *"Will you select a reference genome from your history or use a built-in index?"* to `Use one from the history`
>   - Set *"Select the reference genome"* to the `staph_aur.fasta` file.
>   - Set *"Is this library mate-paired?"* to `Single-end`
>   -  Set *"FASTQ file"* to `Cutadapt on...Output` that you got at the end of the preprocessing section
>   - Set *"Bowtie settings to use"* to `Full parameters list` to change parameters so that bowtie stricly enforces no mismatches.
>       - Set *"Skip the first n reads (-s)"* to `0`
>       - Set *"Maximum number of mismatches permitted in the seed (-n)"* to `0`
>       - Set *"Seed length (-l)"* to `17`
>       ![Full parameters in Bowtie](../../images/tnseq/bowtie1.png)
>       - Set *"Whether or not to try as hard as possible to find valid alignments when they exist (-y)"* to `Try Hard`
>       - Set *"Whether or not to make Bowtie guarantee that reported singleton alignments are 'best' in terms of stratum and in terms of the quality values at the mismatched positions (--best)"* to `Use best`
>       ![Full parameters in Bowtie](../../images/tnseq/bowtie2.png)
>   - **Click** on `Execute`
>
>
> {% icon tool %} *Optional* : Rename your collection for better clarity:
>   - Open the collection by clicking on it in the history panel
>   -  **Click** on the name of the collection `Bowtie...`
>   - Change the title to  `Mapped reads`
>   - Hit Enter key
>
>
>    > ### {% icon question %} Questions
>    >
>    > Why are we strictly enforcing no mismatch mapping?
>    >
>    > > ### {% icon solution %} Solution
>    > > Our reads being very short, the smallest size allowing precise mapping (see [Introduction](#BuildLibrary) ), allowing even one mismatch would risk having reads mapping in wrong positions.
>    > >
>    > {: .solution }
> {: .question}
>
{: .hands_on}



## Getting coverage of the genome

Now that we have mapped the reads on the reference genome, we are going to calculate the coverage of the genome to later cross them with our TA sites position.
In our case, the reads cover the flanking region on one side of the TA site where the transposon inserted. That means we do not want to have the coverage across the whole reads, as it could cover several TA sites, but only the coverage at the end of the read. (See Figure [Mapping read and TA site coverage](#figure-5))


![A read align to the genome with its 3' end covering half of the TA site](../../images/tnseq/Map_cov.png "<b>Mapping read and TA site coverage</b> - The sequenced read cover the 5' region flanking the site of insertion. To assign the read to its correct insertion, we need to compute the coverage at the 3' end of the read.")


> ### {% icon hands_on %} Hands-on:  Compute genome coverage
>
> {% icon tool %} Select the **bamCoverage** tool in the tool bar and run with the following parameters:
>   - Set *"BAM/CRAM file"* to `Mapped Reads` collection
>   - Set *"Bin size in bases"* to `1`
>   - Set *"Scaling/Normalization method"* to `Do not normalize or scale`
>   - Set *"Coverage file format"* to `bedgraph`
>   - Set *"Show advanced options"* to `yes`
>       - Set *"Ignore missing data?"* to `yes` to get only region with read counts
>       - Set *"Offset inside each alignment to use for the signal location."* to `-1` to read the signal of the coverage only at the 3' end of the read.
>   - **Click** on `Execute`
>
{: .hands_on}

## Getting TA sites positions

In order to get the coverage on each TA site we need to prepare a file containing the position of each TA site. As you can see on the figure "[Mapping read and TA site coverage](#MapCoverage)" the read can cover one side or the other of the TA depending on the direction of insertion of the transposon. Depending on the direction of insertion the coverage will be counted on the leftmost position of the TA site or on the rightmost. As we are not considering strand separately in this analyses, we will consider both count as attached to the leftmost base of the TA site. To do that we will create two list of TA site positions, listing the 5' end of each TA site for forward and reverse strand, and then merge them to get a global count per TA site.
We first need to create a fasta file containing the motif :
```
>TA_site
TA
```

> ### {% icon hands_on %} Hands-on:  Get TA sites coordinates
>
> {% icon tool %} Select the **wordmatch** tool in the tool bar and run with the following parameters:
>   - Set *Sequence 1"* to `TA site` file you just created
>   - Set *"Sequence 2"* to the genome sequence file `staph_aur.fasta`
>   - Set *"Word size"* to `2`
>   - Set *"Output Alignment File Format"* to `match`
>   - Set *"Output Feature 1 File Format"* to `GFF`
>   - Set *"Output Feature 2 File Format"* to `GFF`
>   - **Click** on `Execute`
>
{: .hands_on}

Wordmatch tools provides you with three outputs:
    - A gff file containing the locations where sequence 1 maps on sequence 2
    - A gff file containing the locations where sequence 2 maps on sequence 1
    - A report file providing data on each alignment

In our case we are only interested in the locations of TA sites mapping on the genomic sequence. We will therefore use the file looking like that :



| Seqname |Source|Feature|Start|End|Score|Strand|Frame|Group|
|:--------------:|:--------------:|:--------------:|:--------------:|:--------------:|:--------------:|:--------------:|:--------------:|:--------------:|
| NC_007795.1 |wordmatch|misc_feature|5|6|1.000|+|.|Sequence "NC_007795.1.1"; note "TA_site"|
||

The only information we need here are the positions of the 5' end of each TA site for each strand.

> ### {% icon hands_on %} Hands-on:  Get forward and reverse strand positions
>
> {% icon tool %} Select the **Cut columns from a table (cut)** tool in the tool bar and run with the following parameters:
>   - Set *File to cut"* to `Wordmatch on...` file containing the TA sites positions
>   - Set *"Operation"* to `Keep`
>   - Set *"Delimited by"* to `tab`
>   - Set *"Cut by"* to `fields`
>   - Set *"List of Fields"* to `Column: 1` and `Column: 4`  to get forward strand coordinate
>   - **Click** on `Execute`
>   - Add `#forward` tag to the output
>
> {% icon tool %} Select the **Cut columns from a table (cut)** tool in the tool bar and run with the following parameters:
>   - Set *File to cut"* to `Wordmatch on...` file containing the TA sites positions
>   - Set *"Operation"* to `Keep`
>   - Set *"Delimited by"* to `tab`
>   - Set *"Cut by"* to `fields`
>   - Set *"List of Fields"* to  `Column: 1` and  `Column: 5`  to get reverse strand coordinate
>   - **Click** on `Execute`
>   - Add `#reverse` tag to the output
>
{: .hands_on}

The coordinates provided by wordmatch are based on count 1. Meaning the first nucleotide is counted as number 1. However, bamCoverage count the first nucleotide as nucleotide number 0. To be able to compare the two results, we need to shift the TA site positions by 1.

> ### {% icon hands_on %} Hands-on:  Shift TA sites positions
>
> {% icon tool %} To shift the positions of TA sites, select the **Compute an expression on every row** tool in the tool bar and run with the following parameters:
>   - Set *Add expression"* to `c2-1` to shift the position by 1
>   - Set *"as a new column to"* to `Multiple datasets` and select the two files with the tags `forward` and `reverse`
>   - Set *"Round result?"* to `yes`
>   - **Click** on `Execute`
>
> {% icon tool %} To keep only the new coordinates, select the **Cut columns from a table (cut)** tool in the tool bar and run with the following parameters:
>   - Set *File to cut"* to `Multiple datasets` and select the two files output of the previous step
>   - Set *"Operation"* to `Keep`
>   - Set *"Delimited by"* to `tab`
>   - Set *"Cut by"* to `fields`
>   - Set *"List of Fields"* to `Column: 3`  to get new coordinates
>   - **Click** on `Execute`
>
{: .hands_on}



## Merging overall coverage and TA sites to get the coverage of each TA sites

Now that have the files containing the coordinates of our TA sites for both strands. We will cross them with the coverage files to get the coverage on both sides of each TA site.

> ### {% icon hands_on %} Hands-on:  Get coverage of TA sites
>
> {% icon tool %} Select the **Join two files** tool in the tool bar and run with the following parameters:
>   - Set *1st file"* to `bamCoverage on ...` collection
>   - Set *"Column to use from 1st file"* to `Column 2` to select the position of the end of the read
>   - Set *"2nd File"* to the file with the `forward` tags
>   - Set *"Column to use from 2nd file"* to `Column 1`
>   - Set *"Output lines appearing in"*  to `Both 1st & 2nd file, plus unpairable lines from 2st file. (-a 2)` to get all TA sites.
>   - Set *"First line is a header line"* to `No`
>   - Set *"Value to put in unpaired (empty) fields"* to `0`: assign a count of 0 to TA sites not covered in  the bamCoverage file.
>   - **Click** on `Execute`
>   - Add `#forward` tag to the output collections
>
>
> {% icon tool %} Repeat the operation with the `reverse` file
>
>
> {% icon tool %} Select the **Cut columns from a table (cut)** tool in the tool bar and run with the following parameters:
>   - Set *File to cut"* to the collection `Join on collection...` with the `forward` tag
>   - Set *"Operation"* to `Keep`
>   - Set *"Delimited by"* to `tab`
>   - Set *"Cut by"* to `fields`
>   - Set *"List of Fields"* to `Column: 1`  and `Column: 4` to keep only position and coverage
>   - **Click** on `Execute`
>
> {% icon tool %} Repeat the operation with the `reverse` collection
>
{: .hands_on}

We now have a read count for each nucleotides of the TA sites. The insertions counted in the forward strand files are in a different direction than those counted on the reverse strand file. In our case, we are only interested in studying the gene disruption, and therefore we do not want to consider the directions separately. We will add the forward and reverse count together to get a total count per TA site. In order to do that we need to assign the count of both strand to the same nucleotide. We will do that by shifting by one position the count on the reverse strand.

> ### {% icon hands_on %} Hands-on:  Get total count per TA site
>
> {% icon tool %} To shift the positions of reverse strand counts, select the **Compute an expression on every row** tool in the tool bar and run with the following parameters:
>   - Set *Add expression"* to `c1-1` to shift the position by 1
>   - Set *"as a new column to"* to the collection `Cut on...` with the tags `reverse`
>   - Set *"Round result?"* to `yes`
>   - **Click** on `Execute`
>
> {% icon tool %} To select the new coordinates, select the **Cut columns from a table** tool in the tool bar and run with the following parameters:
>   - Set *"Cut columns"* to `c3,c2` to keep only new position and counts in that order
>   - Set *"Delimited by"* to `tab`
>   - Set *From"* to the output of the previous step collection `Compute on collection...`
>   - **Click** on `Execute`
>
> {% icon tool %} Before adding the counts we need to sort the files. Select the **Sort data in ascending or descending order** tool in the *"Text Manipulation"* section of the tool bar and run with the following parameters:
>   - Set *"Sort Query"* to the newest collection with the `forward` tag
>   - Set *"Number of header lines"* to `0`
>   - Set *"on column"* to `Column: 1` to sort on positions
>   - Set *"in"* to `Ascending Order`
>   - Set *"Flavor"* to `Fast numeric sort (-n)`
>   - **Click** on `Execute`
>
> {% icon tool %} Repeat the operation with the newest `reverse` file
>
> {% icon tool %} To merge the two collections, select the **Join two files** tool in the tool bar and run with the following parameters:
>   - Set *1st file"* to `Sort on ...` collection with the `forward` tag
>   - Set *"Column to use from 1st file"* to `Column 1` to join on positions
>   - Set *"2nd File"* o `Sort on ...` collection with the `reverse` tag
>   - Set *"Column to use from 2nd file"* to `Column 1`
>   - Set *"Output lines appearing in"*  to `Both 1st & 2nd file` to get all TA sites.
>   - Set *"First line is a header line"* to `No`
>   - Set *"Value to put in unpaired (empty) fields"* to `.`
>   - **Click** on `Execute`
>
> Now that we have joined our files we want to create a new column with the addition of bth counts
>
> {% icon tool %} Select the **Compute an expression on every row** tool in the tool bar and run with the following parameters:
>   - Set *Add expression"* to `c2+c3` to get the total count
>   - Set *"as a new column to"* to the collection `Join on ...` with the tags `reverse` and `forward`
>   - Set *"Round result?"* to `yes`
>   - **Click** on `Execute`
>
>  {% icon tool %} To get a file with only the position and total count, select the **Cut columns from a table (cut)** tool in the tool bar and run with the following parameters:
>   - Set *File to cut"* to the collection `Compute...`
>   - Set *"Operation"* to `Keep`
>   - Set *"Delimited by"* to `tab`
>   - Set *"Cut by"* to `fields`
>   - Set *"List of Fields"* to `Column: 1`  and `Column: 4` to keep only position and total counts
>   - **Click** on `Execute`
>
> {% icon tool %} Select the **Sort data in ascending or descending order** tool in the *"Text Manipulation"* section of the tool bar and run with the following parameters:
>   - Set *Sort Query"* to the collection `Cut columns on...`
>   - Set *"Number of header lines"* to `0`
>   - Set *"on column"* to `Column: 1` to sort on positions
>   - Set *"in"* to `Ascending Order`
>   - Set *"Flavor"* to `Fast numeric sort (-n)`
>   - **Click** on `Execute`
>
>
{: .hands_on}

# Predicting Essential Genes with Transit

## Predict the essentiality of genes

Now that we have the counts of insertions per TA site, we can use them to predict gene esssentiality with Transit. In order to do that, we need to create a an annotation file in the `prot_table` format, specifique to the Transit tool. You can create this file from a gff3 from GenBank like the one you uploaded earlier.

Verify that the format of your file is `gff3` and not `gff`. If that is not the case you can change the datatype.

{% include snippets/change_datatype.md %}

> ### {% icon hands_on %} Hands-on : Create annotation file in prot_table format
>
>  {% icon tool %} Select the **Convert GFF3 to prot_table for TRANSIT** tool in the TRANSIT section of the tool bar and run with the following parameters:
>   - Set *GenBank GFF file"* to the `gff3` file
>   - **Click** on `Execute`
>
>
{: .hands_on}


Now that we have prepared the annotation file, we can use the count per TA site to predict essential genes using Transit tool.  We will modify a couple parameters from the default settings :
-   We want to ignore counts lower than 2, in order to reduce the number of sites presenting a single read, which could be artefactual.
-   We want to ignore insertion near the extremities of the genes. A disrupted site that would be very close to the border of the gene may not actually disturb the gene function, and therefore not be an actual signal of disruption.

> ### {% icon hands_on %} Hands-on : Predict gene essentiality with Transit
>
>  {% icon tool %} Select the **Gumbel** tool in the TRANSIT section of the tool bar and run with the following parameters:
>   - Set *Input .wig files"* to the collection `Sort...`
>   - Set *Input annotation* to the `Convert...` file generated at the previous step
>   - Set *Smallest read-count to consider* to 2, to ignore single count insertions
>   - Set *Ignore TAs occuring at given fraction of the N terminus.* to 0.1, to ignore 10% of the insertion at the N-terminus extremity of the gene
>   - Set *Ignore TAs occuring at given fraction of the C terminus.* to 0.1, to ignore 10% of the insertion at the C-terminus extremity of the gene
>   - **Click** on `Execute`
>
>
{: .hands_on}

The output of Transit is a tabulated file containing the following columns (you can find more information on [Transit manual page](https://transit.readthedocs.io/en/latest/transit_methods.html)):
-   Gene ID
-   Name of the gene
-   Gene description
-   Number of Transposon Insertions Observed within the Gene
-   Total Number of TA sites within the Gene
-   Length of the Maximum Run of Non-Insertions observed (in number of TA sites)
-   Span of nucleotides for the Maximum Run of Non-Insertions
-   Posterior Probability of Essentiality
-   Essentiality call for the gene : E=Essential, U=Uncertain, NE=Non-Essential, S=too short

We can obtain the list of genes predicted as essential by filtering on the essentiality call.

> ### {% icon hands_on %} Hands-on : Get table of essential genes
>
>  {% icon tool %} Select the **Filter data on any column using simple expressions** tool and run with the following parameters:
>   - Set *With following condition"* to `c9=='E'` to select essential genes`
>   - **Click** on `Execute`
>
>
{: .hands_on}

## Compare the essential genes between two conditions

Now let's compare the results between out two conditions :


> ### {% icon hands_on %} Hands-on : Separate Files from the collection
>  {% icon tool %} Select the **Extract Dataset from a list** tool and run with the following parameters:
>   - Set *Input List"* to the `Filter on...`  result
>   - Set *How should a dataset be selected?"* to `Select by index`
>   - Set *Element index:"* to `0` to select the first dataset
>   - **Click** on `Execute`
>
>  {% icon tool %} Select the **Extract Dataset from a list** tool and run with the following parameters:
>   - Set *Input List"* to the `Filter on...`  result
>   - Set *How should a dataset be selected?"* to `Select by index`
>   - Set *Element index:"* to `1` to select the second dataset
>   - **Click** on `Execute`
{: .hands_on}


> ### {% icon hands_on %} Hands-on : Get gene essential in both conditions
>
>  {% icon tool %} Select the **Join two files** tool and run with the following parameters:
>   - Set *1st file"* to `Control`
>   - Set *Column to use from 1st file"* to `Column : 1` to compare on gene ID
>   - Set *2nd file"* to `Condition`
>   - Set *Column to use from 2nd file"* to `Column : 1` to compare on gene ID
>   - Set *Output lines appearing in"* to `Both 1st and 2nd files` to get common essential genes
>   - **Click** on `Execute`
>
>
{: .hands_on}

> ### {% icon hands_on %} Hands-on : Get gene essential only in control
>
>  {% icon tool %} Select the **Join two files** tool and run with the following parameters:
>   - Set *1st file"* to `Control`
>   - Set *Column to use from 1st file"* to `Column : 1` to compare on gene ID
>   - Set *2nd file"* to `Condition`
>   - Set *Column to use from 2nd file"* to `Column : 1` to compare on gene ID
>   - Set *Output lines appearing in"* to `1st but not 2nd`
>   - **Click** on `Execute`
>
>
{: .hands_on}


> ### {% icon hands_on %} Hands-on : Get gene essential only in condition
>
>  {% icon tool %} Select the **Join two files** tool and run with the following parameters:
>   - Set *1st file"* to `Control`
>   - Set *Column to use from 1st file"* to `Column : 1` to compare on gene ID
>   - Set *2nd file"* to `Condition`
>   - Set *Column to use from 2nd file"* to `Column : 1` to compare on gene ID
>   - Set *Output lines appearing in"* to `2nd but not 1st`
>   - **Click** on `Execute`
>
>
{: .hands_on}
