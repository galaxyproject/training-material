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

In this tutorial, the transposon is a mariner Himar1 with the structure described in the figure [Structure of the tranposon constructs](#figure-2) {% cite Santiago2015 %}.

![Structure of the transposon containing several parcodes and adapters](../../images/tnseq/tranposon_structure.png "<b>Structure of the tranposon constructs</b> - The transposon construct is a mariner transposon with two specific region used to specifically sequence the region upstream of the insertion. The transposon inserts at TA site at the ITR junctions. These ITR junctions have been modified to include a Mme1 restriction site. Using MmeI enzyme to determine the size of the reads allow to have a homogeneous read size and therefore avoid a bias in the representation of the insertions. It also includes a NotI restriction site (cut 21 bp upstream from the restriction site). These two site are the 5' and 3' limits to the genomic DNA we want to sequence. <b>A. Sequence flanking genomic regions</b> After digestion by NotI restriction enzyme, the fragments are attached to biotinylated adaptors that link to NotI restriction site. The attached fragment are then digested by MMeI at a site upstream , where an Illumina primer is then linked. The sequencing is then done, adding Illumina adaptors and an additional barcode to the read for multiplexed sequencing. <b>B. Removing incorrect fragments</b> An insertion can sometimes be composed of one or more copies of the transposon (multimer). There is therefore a risk to select plasmid backbone sequence. To solve this problem, an additional NotI has been add in the backbone to create different length construct, that can later be filtrated (. Different promoters are added to the construct along with an additional 3 bp barcode to analyze differential expression impact, but this will be the subject of another tutorial. (From <a href='#Santiago2015'> Santiago <i>et al.</i> 2015</a>)")


Because of this complex tranposon structure, the reads obtained after sequencing contain a large portion of tranposon sequence for a 16-17 bp genomic sequence. This will necessitate several step of pre-processing to extract this genomic sequence.




## Tnseq analysis
{:.no_toc}

Once we extracted the genomic sequences from the initial reads, we need to locate each of them on the genome to link them to a TA site. To do that we need to map them to a reference genome, link them to a specific insertion site, and then count the number of insertion for each TA site.

Once we have the count of insertion at every insertion site, there is several methods existing to identify essential genes of regions. They can be divided in two major categories (see figure [Methods of TnSeq Analyses](#figure-3) {% cite Chao2016 %}) :
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
>  1. **Cutadapt** {% icon tool %} Split the dataset based on barcodes with:
>    - *"Single-end or Paired-end reads?"* : `Single-end`
>    - *"FASTQ/A file"* : `fastq.gz file` containing the training set of reads.
>    - Unfold `Insert 5' (Front) Adapters`
>    - *"Source"* : `File From History`
>    - *"Choose file containing 5' adapters"* : `condition barcodes` file
>     ![Use cutadapt with front adapters](../../images/tnseq/frontadaptercutadapt.png)
>    - Unfold `Adapter Options`
>    -  *"Maximum error rate"* : `0.15` to allow 1 mismatch
>    -  *"Match times"* : `3` in case the barcode attached several times
>     ![Adapter options in Cutadapt](../../images/tnseq/optionaptercutadapt.png)
>    - Unfold `Output Options`
>    -  *"Report"* : `yes`
>    -  *"Multiple output"* : `yes` to separate the reads into one file per condition
>     ![Output options in Cutadapt](../../images/tnseq/cutadaptOutput.png)
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
>  1. **Cutadapt** {% icon tool %} Remove adapters with:
>    - *"Single-end or Paired-end reads?"* : `Single-end`
>    - *"FASTQ/A file"* : `Cutadapt on data... Output` collection output of the previous step.
>    - Unfold `Insert 3' (End) Adapters`
>        - *"Source"* : `Enter custom Sequence`
>        - *"Enter custom 3' adapter sequence"* : `CGTTATGGCACGC`
>    - Unfold `Adapter Options`
>        - *"Match times"* : `3` in case the barcode attached several times
>    - Unfold `Output Options`
>        - *"Report"* : `yes`
>    - Unfold `Insert 3' (End) Adapters`
>        - *"Source"* : `Enter custom Sequence`
>        - *"Enter custom 3' adapter sequence"* : `CGTTATGGCACGC`
>    - Unfold `Adapter Options`
>        - *"Match times"* : `3` in case the barcode attached several times
>    - Unfold `Output Options`
>        - *"Report"* : `yes`
>
>
>  2. **Cutadapt** {% icon tool %} Filter reads based on length with:
>    - Set *"Single-end or Paired-end reads?"* : `Single-end`
>    - Set *"FASTQ/A file"* : `Cutadapt on data... Output` collection output of the previous step.
>    - Unfold `Filter Options`
>        - *"Minimum length"* : `64`
>        - *"Maximum length"* : `70`
>    - Unfold `Output Options`
>        - *"Report"* to `yes`
>
>    {% include snippets/select_collection.md %}
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
>  1. **Cutadapt** {% icon tool %} with:
>    - *"Single-end or Paired-end reads?"* to `Single-end`
>    - *"FASTQ/A file"* : `Cutadapt on data... Output` collection output of the previous step.
>    - Unfold `Insert 3' (End) Adapters`
>       - *"Source"* : `File From History`
>       - *"Choose file containing 3' adapters"* : `construct barcodes` file
>    - Unfold `Adapter Options`
>       - *"Match times"* : `3` in case the barcode attached several times
>    - Unfold `Output Options`
>       - *"Report"* : `yes`
>       - *"Multiple output"* : `yes` to separate the reads into one file per condition
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
>  1. **Cutadapt** {% icon tool %} with:
>    - *"Single-end or Paired-end reads?"* to `Single-end`
>    - *"FASTQ/A file"* : `Cutadapt on data... Output` nested collection output of the previous step.
>    - Unfold`Insert 3' (End) Adapters`
>        - *"Source"* : `Enter custom Sequence`
>        - *"Enter custom 3' adapter sequence"* : `ACAGGTTGGATGATAAGTCCCCGGTCTATATTGAGAGTAACTACATTT`
>    - Unfold `Adapter Options`
>        - *"Maximum error rate"* : `0.15`
>    - Unfold `Output Options`
>        - *"Report"* to `yes`
>
{: .hands_on}

Verify that the majority of the read have been trimmed. Now that we isolated the genomic sequences from the initial reads. We want to align them to count how many insertion have been retained at each TA sites.

# Counting the number of insertion per TA sites

## Aligning the reads to a reference genome

The first step is to map our read to the reference genome. We are going to use the tool Bowtie. We could also use Bowtie2 with an end-to-end option, but Bowtie is more suitable for very short reads like ours (16-17 bp).

> ### {% icon hands_on %} Hands-on:  Map reads with Bowtie
>
>  1. **Map with Bowtie for Illumina** {% icon tool %} with:
>   - *"Will you select a reference genome from your history or use a built-in index?"* : `Use one from the history`
>   - *"Select the reference genome"* : `staph_aur.fasta` file.
>   - *"Is this library mate-paired?"* : `Single-end`
>   - *"FASTQ file"* : `Cutadapt on...Output` that you got at the end of the preprocessing section
>   - *"Bowtie settings to use"* : `Full parameters list` to change parameters so that bowtie strictly enforces no mismatches.
>       - *"Skip the first n reads (-s)"* : `0`
>       - *"Maximum number of mismatches permitted in the seed (-n)"* : `0`
>       - *"Seed length (-l)"* : `17`
>       ![Full parameters in Bowtie](../../images/tnseq/bowtie1.png)
>       - *"Whether or not to try as hard as possible to find valid alignments when they exist (-y)"* : `Try Hard`
>       - *"Whether or not to make Bowtie guarantee that reported singleton alignments are 'best' in terms of stratum and in terms of the quality values at the mismatched positions (--best)"* : `Use best`
>       ![Full parameters in Bowtie](../../images/tnseq/bowtie2.png)
>
>
>  2. *Optional* : Rename your collection for better clarity
>
>
>    > ### {% icon tip %} Tip: Renaming a collection
>    >
>    > 1. Open the collection by clicking on it in the history panel
>    > 2. Click on the name of the collection
>    > 3. Enter the new name
>    > 4. Hit the `Enter` key
>    {: .tip}
>
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
>  1. **bamCoverage** {% icon tool %} with:
>   - *"BAM/CRAM file"* : `Mapped Reads` collection
>   - *"Bin size in bases"* : `1`
>   - *"Scaling/Normalization method"* : `Do not normalize or scale`
>   - *"Coverage file format"* : `bedgraph`
>   - *"Show advanced options"* : `yes`
>       - *"Ignore missing data?"* : `yes` to get only region with read counts
>       - *"Offset inside each alignment to use for the signal location."* : `-1` to read the signal of the coverage only at the 3' end of the read.
>
{: .hands_on}

## Getting TA sites positions

In order to get the coverage on each TA site we need to prepare a file containing the position of each TA site. As you can see on the figure "[Mapping read and TA site coverage](#figure-5)" the read can cover one side or the other of the TA depending on the direction of insertion of the transposon. Depending on the direction of insertion the coverage will be counted on the leftmost position of the TA site or on the rightmost. As we are not considering strand separately in this analyses, we will consider both count as attached to the leftmost base of the TA site. To do that we will create two list of TA site positions, listing the 5' end of each TA site for forward and reverse strand, and then merge them to get a global count per TA site.
We first need to create a fasta file containing the motif :
```
>TA_site
TA
```

> ### {% icon hands_on %} Hands-on:  Get TA sites coordinates
>
>  1. **wordmatch** {% icon tool %} with:
>   - *Sequence 1"* : `TA site` file you just created
>   - *"Sequence 2"* : `staph_aur.fasta`
>   - *"Word size"* : `2`
>   - *"Output Alignment File Format"* : `match`
>   - *"Output Feature 1 File Format"* : `GFF`
>   - *"Output Feature 2 File Format"* : `GFF`
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
>  1. **Cut columns from a table (cut)** {% icon tool %} with:
>   - *File to cut"* to `Wordmatch on...` file containing the TA sites positions
>   - *"Operation"* to `Keep`
>   - *"Delimited by"* to `tab`
>   - *"Cut by"* to `fields`
>   - *"List of Fields"* to `Column: 1` and `Column: 4`  to get forward strand coordinate
>   
>  2. Add the `#forward` tag to the output
>
>
>  3. **Cut columns from a table (cut)** {% icon tool %} with:
>   - *File to cut"* : `Wordmatch on...` file containing the TA sites positions
>   - *"Operation"* : `Keep`
>   - *"Delimited by"* : `tab`
>   - *"Cut by"* : `fields`
>   - *"List of Fields"* : `Column: 1` and  `Column: 5`  to get reverse strand coordinate
>
>  4. Add the `#reverse` tag to the output
>
>    > ### {% icon tip %} Tip: Adding a tag to a collection
>    > * Click on the collection
>    > * Add a tag starting with `#` in the `Add tags` field
>    >
>    >     Tags starting with `#` will be automatically propagated to the outputs of tools using this dataset.
>    >
>    > * Hit the `Enter` key
>    > * Check that the tag is appearing below the dataset name
>    >
>    {: .tip}
>
>
>
{: .hands_on}


The coordinates provided by wordmatch are based on count 1. Meaning the first nucleotide is counted as number 1. However, bamCoverage count the first nucleotide as nucleotide number 0. To be able to compare the two results, we need to shift the TA site positions by 1.

> ### {% icon hands_on %} Hands-on:  Shift TA sites positions
>
>  1. **Compute an expression on every row** {% icon tool %} to shift the positions of TA sites with:
>   - *Add expression"* : `c2-1` to shift the position by 1
>   - *"as a new column to"* : `Multiple datasets` and select the two files with the tags `forward` and `reverse`
>   - *"Round result?"* : `yes`
>
>  2. **Cut columns from a table (cut)** {% icon tool %} to keep only the new coordinates with:
>   - *File to cut"* : `Multiple datasets` and select the two files output of the previous step
>   - *"Operation"* : `Keep`
>   - *"Delimited by"* : `tab`
>   - *"Cut by"* : `fields`
>   - *"List of Fields"* : `Column: 3`  to get new coordinates
>
>
>    {% include snippets/select_multiple_datasets.md %}
>
{: .hands_on}



## Merging overall coverage and TA sites to get the coverage of each TA sites

Now that have the files containing the coordinates of our TA sites for both strands. We will cross them with the coverage files to get the coverage on both sides of each TA site.

> ### {% icon hands_on %} Hands-on:  Get coverage of TA sites
>
>  1. **Join two file** {% icon tool %} with:
>   - *1st file"* : `bamCoverage on ...` collection
>   - *"Column to use from 1st file"* : `Column 2` to select the position of the end of the read
>   - *"2nd File"* : the file with the `forward` tags
>   - *"Column to use from 2nd file"* : `Column 1`
>   - *"Output lines appearing in"*  : `Both 1st & 2nd file, plus unpairable lines from 2st file. (-a 2)` to get all TA sites.
>   - *"First line is a header line"* : `No`
>   - *"Value to put in unpaired (empty) fields"* : `0` to assign a count of 0 to TA sites not covered in the bamCoverage file.
>
>  2. Add `#forward` tag to the output collections
>
>  3. Repeat the operation with the `reverse` file
>
>  4. **Cut columns from a table (cut)** {% icon tool %} with:
>   - *File to cut"* : `Join on collection...` collection with the `forward` tag
>   - *"Operation"* : `Keep`
>   - *"Delimited by"* : `tab`
>   - *"Cut by"* : `fields`
>   - *"List of Fields"* : `Column: 1` and `Column: 4` to keep only position and coverage
>
>  5. Repeat the operation with the `reverse` collection
>
{: .hands_on}

We now have a read count for each nucleotides of the TA sites. The insertions counted in the forward strand files are in a different direction than those counted on the reverse strand file. In our case, we are only interested in studying the gene disruption, and therefore we do not want to consider the directions separately. We will add the forward and reverse count together to get a total count per TA site. In order to do that we need to assign the count of both strand to the same nucleotide. We will do that by shifting by one position the count on the reverse strand.

> ### {% icon hands_on %} Hands-on:  Get total count per TA site
>
>  1. **Compute an expression on every row** {% icon tool %} to shift the positions of reverse strand counts with:
>   - *Add expression"* : `c1-1` to shift the position by 1
>   - *"as a new column to"* : `Cut on...` collection with the tags `reverse`
>   - *"Round result?"* : `yes`
>
>  2. **Cut columns from a table** {% icon tool %} with:
>   - *"Cut columns"* : `c3,c2` to keep only new position and counts in that order
>   - *"Delimited by"* : `tab`
>   - *From"* : `Compute on collection...`
>
>  3. **Cut columns from a table** {% icon tool %} in the *"Text Manipulation"* section of the tool bar to sort the files with:
>   - *"Sort Query"* : the newest collection with the `forward` tag
>   - *"Number of header lines"* : `0`
>   - *"on column"* : `Column: 1` to sort on positions
>   - *"in"* : `Ascending Order`
>   - *"Flavor"* : `Fast numeric sort (-n)`
>
>  4. Repeat the operation with the newest `reverse` file
>
>  6. **Join two files** {% icon tool %} to merge the two collections and get the forward and reverse count in two different columns with:
>   - *1st file"* : `Sort on ...` collection with the `forward` tag
>   - *"Column to use from 1st file"* : `Column 1` to join on positions
>   - *"2nd File"* : `Sort on ...` collection with the `reverse` tag
>   - *"Column to use from 2nd file"* : `Column 1`
>   - *"Output lines appearing in"*  : `Both 1st & 2nd file` to get all TA sites.
>   - *"First line is a header line"* : `No`
>   - *"Value to put in unpaired (empty) fields"* : `.`
>
>  7. **Compute an expression on every row** {% icon tool %} to get the total read count with:
>   - *Add expression"* : `c2+c3` to get the total count
>   - *"as a new column to"* : the collection `Join on ...` with the tags `reverse` and `forward`
>   - *"Round result?"* : `yes`
>
>  8. **Cut columns from a table (cut)** {% icon tool %} to get a file with only the position and total count with:
>   - *File to cut"* : the collection `Compute...`
>   - *"Operation"* : `Keep`
>   - *"Delimited by"* : `tab`
>   - *"Cut by"* : `fields`
>   - *"List of Fields"* : `Column: 1` and `Column: 4` to keep only position and total counts
>
>  9. **Sort data in ascending or descending order** {% icon tool %} in the *"Text Manipulation"* section of the tool bar with:
>   - *Sort Query"* : the collection `Cut columns on...`
>   - *"Number of header lines"* : `0`
>   - *"on column"* : `Column: 1` to sort on positions
>   - *"in"* : `Ascending Order`
>   - *"Flavor"* : `Fast numeric sort (-n)`
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
>  1. **Convert GFF3 to prot_table for TRANSIT** {% icon tool %} with:
>   - *GenBank GFF file"* : the `gff3` file
>
>
{: .hands_on}


Now that we have prepared the annotation file, we can use the count per TA site to predict essential genes using Transit tool.  We will modify a couple parameters from the default settings :
-   We want to ignore counts lower than 2, in order to reduce the number of sites presenting a single read, which could be artefactual.
-   We want to ignore insertion near the extremities of the genes. A disrupted site that would be very close to the border of the gene may not actually disturb the gene function, and therefore not be an actual signal of disruption.

> ### {% icon hands_on %} Hands-on : Predict gene essentiality with Transit
>
>  1. **Gumbel** {% icon tool %} with:
>   - *Input .wig files"* : the collection `Sort...`
>   - *Input annotation* : the `Convert...` file generated at the previous step
>   - *Smallest read-count to consider* : `2`, to ignore single count insertions
>   - *Ignore TAs occuring at given fraction of the N terminus.* : `0.1`, to ignore 10% of the insertion at the N-terminus extremity of the gene
>   - *Ignore TAs occuring at given fraction of the C terminus.* : `0.1`, to ignore 10% of the insertion at the C-terminus extremity of the gene
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
>  1. **Filter data on any column using simple expressions** {% icon tool %} with:
>   - *With following condition"* : `c9=='E'` to select essential genes`
>
>
{: .hands_on}

## Compare the essential genes between two conditions

Now let's compare the results between out two conditions :


> ### {% icon hands_on %} Hands-on : Separate Files from the collection
>
>  1. **Extract Dataset from a list** {% icon tool %} with:
>   - *Input List"* : the `Filter on...`  result
>   - *How should a dataset be selected?"* : `Select by index`
>   - *Element index:"* : `0` to select the first dataset
>
>  2. **Extract Dataset from a list** {% icon tool %} with:
>   - *Input List"* : the `Filter on...`  result
>   - *How should a dataset be selected?"* : `Select by index`
>   - *Element index:"* : `1` to select the second dataset
>
{: .hands_on}


> ### {% icon hands_on %} Hands-on : Get gene essential in both conditions
>
>  1. **Join two files** {% icon tool %} with:
>   - *1st file"* : `Control`
>   - *Column to use from 1st file"* : `Column : 1` to compare on gene ID
>   - *2nd file"* : `Condition`
>   - *Column to use from 2nd file"* : `Column : 1` to compare on gene ID
>   - *Output lines appearing in"* : `Both 1st and 2nd files` to get common essential genes
>
>
{: .hands_on}

> ### {% icon hands_on %} Hands-on : Get gene essential only in control
>
>  1. **Join two files** {% icon tool %} with:
>   - *1st file"* : `Control`
>   - *Column to use from 1st file"* : `Column : 1` to compare on gene ID
>   - *2nd file"* : `Condition`
>   - *Column to use from 2nd file"* : `Column : 1` to compare on gene ID
>   - *Output lines appearing in"* : `1st but not 2nd`
>
>
{: .hands_on}


> ### {% icon hands_on %} Hands-on : Get gene essential only in condition
>
>  1. **Join two files** {% icon tool %} with:
>   - *1st file"* : `Control`
>   - *Column to use from 1st file"* : `Column : 1` to compare on gene ID
>   - *2nd file"* : `Condition`
>   - *Column to use from 2nd file"* : `Column : 1` to compare on gene ID
>   - *Output lines appearing in"* : `2nd but not 1st`
>
>
{: .hands_on}
