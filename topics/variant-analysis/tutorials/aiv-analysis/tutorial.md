---
layout: tutorial_hands_on

title: Avian influenza viral strain analysis from gene segment sequencing data
enable: false
level: Intermediate
questions:
- With reassortment of gene segments being a common event in avian influenza virus (AIV) evolution, does it make sense to use a reference-based mapping approach for constructing consensus genome sequences for AIV samples?
- Is it possible to reuse existing tools and workflows developed for the analysis of sequencing data from other viruses?
- How can we obtain meaningful phylogenetic insight from AIV consensus sequences? 
objectives:
- Determine how reassortment impacts reference-based mapping approaches
- Use a collection of per-segment reference sequences to construct a hybrid reference genome that is sufficiently close to a sequenced sample to be useful as a reference for mapping
- Construct a sample consensus genome from mapped reads
- Generate per-segment phylogenetic trees of AIV consensus sequences
time_estimation: 4H
key_points:
- Reassortment of gene segments makes reference-based mapping of influenza sequencing data challenging
- An alternative to *de-novo* assembly can be mapping to a dynamically chosen reference genome
- Variant calling and consensus genome construction can follow workflows used also for other viral sequence data 
- Standard phylogenetic tools can be used to find relationships between influenza samples but should be used on a per-segment basis
contributors:
- wm75
tags:
- virology
- public health
---


# Introduction

Of the four species of influenza viruses (Influenza A-D), *Influenza A* is the most virulent in human hosts and subtypes of it have been responsible for all historic flu pandemics.

The different influenza species infect distinct animal hosts (though all of them can infect humans and pigs) and the most important natural reservoir for *Influenza A* are birds (in particular, wild aquatic birds), in which it causes **Avian influenza**.

Importantly, all flu pandemics of the last century have been triggered by **reassortment** events between avian and human *Influenza A* strains (although at least the 2009 swine flu pandemic involved an additional mixing event with an *Influenza A* strain from pigs).

Reassortment events are a key trait of *Influenza A* made possible by its wide range of natural hosts combined with two molecular characteristics of the species:

1. the linear negative-sense single-stranded RNA genomes of influenza viruses (and of all viruses of the *Orthomyxoviridae* family) are segmented, i.e. consist of several (eight in the case of *Influenza A* and *B*) distinct pieces of RNA, which typically encode just one, sometimes two viral proteins.
2. the mutation rate in the genome is high in influenza viruses, and particularly high in *Influenza A*, since their RNA polymerase lacks exonuclease activity necessary for proof-reading.

Together these characteristics enable *Influenza A* to evolve relatively rapidly in non-human hosts, then re-adapt to humans through exchange of segments (reassortment) in a host co-infected with a human and, *e.g.*, an avian strain. If such a reassortment introduces new variants of the two most antigenic proteins of influenza, HA (hemagglutinin) and NA (neuraminidase), this antigenic shift provides the reassorted strain with immune escape potential, which, if it happens in a genetic background compatible with efficient transmission in humans, can trigger an unusually strong wave of influenza or even a new pandemic.

In order to estimate the likelihood of an epidemic or pandemic influenza event in the near future it is, therefor, important to monitor closely the genome composition of circulating Avian Influenza strains in wild and domestic birds and the huge technological advances over the last decade make it possible to use next-generation sequencing for this purpose.

At the same time, the segmented nature of their genomes combined with high genetic variability, especially in the HA and NA genes, requires rather specialized bioinformatics workflows for the analysis of data from influenza viruses compared to other viruses with similar genome size (like, *e.g.* SARS-CoV-2):

The viral surface proteins HA and (to a lesser extent) NA are the main targets of the host antibody response and are, thus, under constant selection pressure to mutate into forms capable of evading an existing host immune response. As a consequence, these segments have evolved into a much richer panel of sequence variants than the other segments, to the point that sequences of the HA segment can, at present, be classified into 18 distinct subfamilies, H1-H18, while there are 11 recognized subfamilies of NA, and *Influenza A* strains are subtyped (as, for example, H5N1, H3N2, etc.) according to the combination of HA- and NA-encoding segments they are carrying.

Importantly, the sequence diversity of HA and (again to a lesser extent) of NA segments is big enough to prevent a naive approach of mapping sequenced reads to one specific agreed-upon *Influenza A* reference sequence. While this would work for the other six segments, mapping software would regularly fail to find enough plausible mappings for sequenced reads of HA and NA origin to continue analysis with. This is why, in this tutorial we are going to explore an alternative approach, which is also mapping-based but chooses a suitable reference for each segment dynamically based on the input sequencing data.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Prepare analysis history and data

Any analysis should get its own Galaxy history. So let's start by creating a new one:

> <hands-on-title>Prepare the Galaxy history</hands-on-title>
>
> 1. Create a new history for this analysis
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Rename the history
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
{: .hands_on}

## Get reference and sequenced samples data

In this tutorial you are going to work on a single avian influenza sample sequenced in paired-end mode on the Illumina platform, *i.e.* we are going to download two datasets of sequenced reads for that sample.

In addition, we are going to base the analysis on a small collection of multiple reference sequences for each influenza gene segment. We prepared this collection for you from public INSAFlu data.

If you have your own curated collection of reference sequences, you should be able to use it for following this tutorial without any problem. Note, however, that the analysis results referred to in many of the questions will be different if you are exchanging the reference collection.
If you want to learn how you can create a Galaxy collection from your own references that is structured like the default we are suggesting here, you may want to follow the dedicated tutorial on [Using dataset collections]({% link topics/galaxy-interface/tutorials/collections/tutorial.md %}).

> <hands-on-title>Get the data</hands-on-title>
>
> 1. {% tool [Upload](upload1) %} the forward and reverse reads of your sequenced sample to your history and turn them into a *Paired Collection*.
>    1. In the *"Download from web or upload from disk"* dialogue, switch to the `Collection` tab,
>       then configure the drop-down select boxes on that tab like this:
>       - *"Collection Type"*: `Pair`
>       - *"File Type"*: `fastqsanger.gz`
>
>    2. Click on `Paste/Fetch data` and copy the following links over to the empty text box
>
>       ```
>       https://usegalaxy.eu/api/datasets/4838ba20a6d86765171e4f201205515c/display?to_ext=fastqsanger.gz
>       https://usegalaxy.eu/api/datasets/4838ba20a6d867657d09280baae4e1ec/display?to_ext=fastqsanger.gz
>       ```
>
>       and press `Start`.
>
>    3. Wait for the `Build` button to become enabled, click it, and, in the next dialogue, give a suitable **Name** to the new collection.
>
> 2. Import the reference data collection into your history
>    1. Open the [shared history](https://usegalaxy.eu/u/wolfgang-maier/h/influenza-resources) with the pre-prepared INSAFlu reference data
>    2. Click on {% icon new-history %} *"Import"* in the upper right corner of the page
>    3. Choose a name for the new history (or go with the default), confirm by clicking Import and wait for the page to reload
>    4. In the history panel (which should now show the newly imported history),
>       click on the {% icon galaxy-gear %} icon in the bar just above the history's datasets,
>       select `Copy Datasets`, and wait for the center panel to refresh
>    5. In the center panel,
>       - make sure the imported history is selected as the *"Source History"*
>       - select the history you created for this tutorial as the *"Destination History"*
>       - under *"Source History"* select the `References per segment (INSAFlu)` collection for copying
>       - click *"Copy History Items"*
>    6. Switch to the tutorial history
>
{: .hands_on}

## Inspect the reference data

> <hands-on-title>What is in the reference collection?</hands-on-title>
>
> 1. Step inside the uploaded reference collection by clicking on it.
> 2. Expand individual elements of the collection by clicking on them to reveal more details about them.
> 3. View the content of any element by clicking its {% icon galaxy-eye %}.
>
>    You can then scroll through the data and use your browser's search function to look for particular sequences.
>
>    > <question-title>Questions</question-title>
>    >
>    > 1. How many elements are in the collection? What do they represent?
>    > 2. How many sequences are stored in each element?
>    > 3. Are all subtypes of HA and NA represented in the sequences?
>    > 4. Are there sequences of non-A Influenza species?
>    >
>    > > <solution-title>Answers</solution-title>
>    > >
>    > > 1. The collection consists of eight elements and each element provides reference sequences for one genome segment.
>    > > 2. Each element has data for 56 sequences (Galaxy knows how to interpret the fasta format of the sequencing data and counts sequences for you; you can see the count for any of the elements by expanding it).
>    > > 3. The collection has reference sequence data for H1-H18 (with the exception of H11) and for N1-N11. For many of those subtypes, however, there is only one reference, i.e., within-subtype variation is not captured well by the collection.
>    > > 4. There are 6 *Influenza B* references in the collection. Since *Influenza B* is not normally seen in birds and we are going to analyze an avian influenza sample here, these can be considered controls: if we get an *Influenza B* assignment for the sample at any step in our analysis, we know something is very suspicious. This also means that we only have 50 *informative* references in the collection.
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

{% icon galaxy-info %} There are lots of text processing and filtering tools for Galaxy, of which this tutorial is going to introduce just a small subset. You could use them to extract statistics about which subtypes are found how many times in the collection of references instead of relying on scrolling and searching. If you are interested and have the time, you may want to try the [Data Manipulation Olympics]({% link topics/introduction/tutorials/data-manipulation-olympics/tutorial.md %}) to learn more about available such tools and how to combine them.

# Quality control

As a very first step, we would like to make sure that we base our analysis only on the high-quality parts of the sequenced reads.

NGS reads often have lower base calling quality near their ends, so here we are going to trim low-quality stretches of bases from both ends. In addition, we are going to discard reads shorter than 30 bases after trimming. This is reasonable since our next step will involve extracting all possible 21-mers from the reads and there are very few possibilities for a read shorter than 30 bases.

The tool **fastp** lets us perform these tasks and obtain a nice quality report for our reads before and after processing in one go, but many other options exist to perform sequenced reads quality control and trimming/filtering in Galaxy and the dedicated tutorial on [quality control]({% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %}) introduces more of them.

> <hands-on-title>QC and read trimming/filtering with fastp</hands-on-title>
>
> 1. {% tool [fastp](toolshed.g2.bx.psu.edu/repos/iuc/fastp/fastp/0.23.2+galaxy0) %}
>    - *"Single-end or paired reads"*: `Paired Collection`
>    - *"Select paired collection(s)"*: the uploaded paired collection of sequenced reads
>    - In *"Filter Options"* under *"Length filtering options"*:
>      - *"Length required"*: `30`
>    - In *"Read Modification Options"* under *"Per read cutting by quality options"*:
>      - *"Cut by quality in front (5')"*: `Yes`
>      - *"Cut by quality in tail (3')"*: `Yes`
>      - *"Cutting mean quality"*: `30`
>
>    Running the tool will result in two outputs. One is a new paired collection with the processed reads, the other one is a report of initial quality, the processing actions performed and their effect on key quality metrics.
>
> 2. {% icon galaxy-eye %} **Inspect** the report and try to answer the following questions:
>
>    > <question-title>Questions</question-title>
>    >
>    > 1. Which set of reads (forward or reverse reads) did profit most quality-wise from our low-quality base trimming?
>    > 2. What percentage of reads got discarded completely with our settings, and what percentage of bases?
>    > 3. Why, according to the Filtering result, have some reads been discarded because of too many Ns?
>    > 
>    > > <solution-title>Answers</solution-title>
>    > >
>    > > 1. The forward reads were of worse quality in particular near their 3'-ends than the reverse reads and, consequently, were affected more strongly by trimming (see the different "quality" plots in the report).
>    > > 2. Less than 2% of reads got discarded completely (the "Filtering result" section of the report says 98.03% of all reads were retained). In contrast, around 10% of all bases got discarded (compare total bases before and after filtering) so most reads got trimmed but not discarded.
>    > > 3. **fastp** is a tool with many hidden default actions in its various sections. Inside the *Filter Options* check *Quality filtering options* and the various defaults mentioned in the parameters help therein. One of these defaults is to discard reads with more than 5 Ns. These defaults are often reasonable, but it's good to be aware of them.
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

# Per-segment subtyping and hybrid reference construction

Our goal at this step is to find best matches between our sequencing data and the reference sequences of each genome segment. This will:
- give us preliminary subtyping results with regard to the HA and NA segments for the sequenced sample
- suggest the best reference to map the sequencing data to for each segment.

We are then going to combine the best reference of each segment into a hybrid reference genome to use for mapping our sequenced reads against.

To identify the best matching reference segments, we are going to run the tool **VAPOR** asking it to report the stats for any hits it identifies.

> <hands-on-title>Exploring best-matching reference segment scores</hands-on-title>
>
> 1. {% tool [VAPOR](toolshed.g2.bx.psu.edu/repos/iuc/vapor/vapor/1.0.2+galaxy3) %}
>    - {% icon param-collection %} *"Reference sequences"*: `References per segment (INSAFlu)`
>    - *"Type of sequencing data"*: `Paired-end as collection`
>    - *"Paired collection of sequenced reads"*: quality-trimmed reads; output of **fastp**
>    - *"Desired output"*: `Return scores of best matches`
>    - *"Limit number of reported matches to"*: `0`
>    - In *"Optional arguments"*:
>      - *"Read kmer filtering threshold"*: `0.1`
>      - *"Minimum k-mer proportion"*: `0.0`
>
> 2. {% icon galaxy-eye %} Explore the output collection produced by the tool
>
>    > <question-title>Questions</question-title>
>    >
>    > 1. Is there a difference between the results for segment 4 (encoding HA), segment 6 (encoding NA) and those for the remaining segments?
>    > 2. According to the tool, what is the likely subtype with regard to HA and NA of the sample?
>    > 
>    > > <solution-title>Answers</solution-title>
>    > >
>    > > 1. The values for `% of query bases in reads` and for `Total score` are drastically lower for the HA and somewhat lower for the NA segment than for any of the other segments. This could be due to the structure of the collection of references we have used, but also fits very well with the higher selection pressure on NA and, in particular on HA, as the main antigenic proteins of the virus, which leads to higher variability in these segments than in the others.
>    > > 2. The best match (assigned an almost 2x higher score than the second-best match) found for the HA segment is from an H4N6 strain. The top two matches with regard to the NA segment (again assigned ~ 2x higher scores than the next best matches) are from strains with the N6 subtype of NA.
>    > >
>    > >    The most likely subtype of the sample, thus, appears to be H4N6.
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

Now that we have established that things *may* make sense, we can use the output of **VAPOR** to extract the actual sequence of the top hit for each reference segment. We then concatenate these best matches into a hybrid reference genome for mapping.

> <hands-on-title>Obtaining sequences of top hits identified by VAPOR</hands-on-title>
>
> 1. {% tool [Replace parts of text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_find_and_replace/1.1.4) %} to extract just the name of the sequence from each line of VAPOR's output
>    - {% icon param-collection %} *"File to process"*: collection of score outputs of **VAPOR**
>    - In {% icon param-repeat %} *"1. Find and Replace"*:
>      - *"Find pattern"*: `^.+\t>(.+)$`
>      - *"Replace with"*: `$1`
>      - *"Find-Pattern is a regular expression"*: `Yes`
>
> 2. {% tool [Select first lines from a dataset](Show beginning1) %} to get the first line, i.e. the best match from VAPOR's output
>    - *"Select first"*: `1` lines
>    - {% icon param-collection %} *"from"*: output of **Replace**
>
> 3. {% tool [seqtk_subseq](toolshed.g2.bx.psu.edu/repos/iuc/seqtk/seqtk_subseq/1.3.1) %} to extract the reference sequences based on their names reported by VAPOR
>    - {% icon param-collection %} *"Input FASTA/Q file"*: `References per segment (INSAFlu)`
>    - *"Select source of sequence choices"*: `FASTA/Q ID list`
>    - {% icon param-collection %} *"Input ID list"*: the collection output of **Select first**
>
> 4. {% tool [Collapse Collection](toolshed.g2.bx.psu.edu/repos/nml/collapse_collections/collapse_dataset/5.1.0) %} to combine the best-matching sequence of each segment into one dataset
>    - {% icon param-collection %} *"Collection of files to collapse into single dataset"*: the selected sequences collection produced by **seqtk_subseq**
>
{: .hands_on}

At this point, you may want to {% icon galaxy-eye %} inspect the output of the last step to see if it is an eight-segments reference genome as expected, and if the segments correspond to the top hits found by VAPOR.

# Mapping to a hybrid reference

If things went well, the hybrid reference we just obtained should be close enough across all segments to our sample to allow successful mapping of reads. Before we start the mapping we may want to truncate the segment names in our hybrid reference genome though because currently these names still reflect the full origin of the segment sequences, but from now on we are fine with just the segment ID. We can use the **Replace** tool from before again to truncate the names.

> <hands-on-title>Shortening sequence titles</hands-on-title>
>
> 1. {% tool [Replace parts of text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_find_and_replace/1.1.4) %}
>    - {% icon param-file %} *"File to process"*: the hybrid reference genome; output of **Collapse Collection**
>    - In {% icon param-repeat %} *"1. Find and Replace"*:
>      - *"Find pattern"*: `^>([^|]+).+$`
>      - *"Replace with"*: `>$1`
>      - *"Find-Pattern is a regular expression"*: `Yes`
>
{: .hands_on}

Having polished the titles of the segments in our hybrid reference genome we are finally ready for mapping, which we will carry out with **BWA-MEM**, clean up a bit with **Samtools view** and produce a quality report for with **QualiMap BamQC**.

> <hands-on-title>Read mapping and quality control</hands-on-title>
>
> 1. {% tool [Map with BWA-MEM](toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.17.2) %}
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>    - {% icon param-file %} *"Use the following dataset as the reference sequence"*: hybrid reference genome with shortened names; output of **Replace**
>    - *"Single or Paired-end reads"*: `Paired Collection`
>    - *"Select a paired collection"*: quality-trimmed reads; output of **fastp** in QC and Trimming
>
> 2. {% tool [Samtools view](toolshed.g2.bx.psu.edu/repos/iuc/samtools_view/samtools_view/1.15.1+galaxy0) %}
>    - {% icon param-file %} *"SAM/BAM/CRAM data set"*: mapped reads BAM dataset; output of **BWA-MEM**
>    - *"What would you like to look at?"*: `A filtered/subsampled selection of reads`
>    - In *"Configure filters"*:
>      - *"Filter by quality"*: `20`
>      - *"Require that these flags are set"*: `Read is paired` and `Read is mapped in a proper pair`
>
> 3. {% tool [QualiMap BamQC](toolshed.g2.bx.psu.edu/repos/iuc/qualimap_bamqc/qualimap_bamqc/2.2.2d+galaxy3) %}
>    - {% icon param-file %} *"Mapped reads input dataset"*: filtered mapped reads BAM dataset; output of **Samtools view**
>    - *"Skip duplicate reads"*: `Unselect all`
>    - In *"Settings affecting specific plots"*:
>      - *"Number of bins to use in across-reference plots"*: `40`
>
> 4. {% icon galaxy-eye %} Study the report generated with QualiMap
>
>    > <question-title>Questions</question-title>
>    >
>    > 1. What is the coverage of each segment by the sequenced reads, and is it uniform?
>    > 2. Look for a plot showing read mapping quality across the reference. What can you conclude?
>    >
>    > > <solution-title>Answers</solution-title>
>    > >
>    > > 1. Coverage is rather different for the different segments. Looking at "Mean coverage" and its "Standard deviation" in the "Chromosome stats" table and at the "Coverage across reference" plot, coverage of the HA segment seems to be most critical since it approaches zero in some regions.
>    > > 2. Mapping quality is almost constant at or very near 60 (which happens to be the maximum mapping quality value emitted by BWA-MEM) across all segments with the exception of HA.
>    > >
>    > > These two observations combined with the VAPOR statistics for HA show again that even the best matching reference for this segment is not exactly close to the sequence of our sample. The read mapper had difficulties placing the sequenced reads on that sub-optimal reference and it looks as if we might have very little sequence information for some HA regions.
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

# Consensus sequence construction

From the polished mapping of reads to our custom reference we can now construct the consensus sequence of our sample.

Unfortunately, the tool we are going to use for this, **ivar consensus**, is not capable of working with more than one reference name at a time, but because this is influenza data we have mappings to 8 different segments described in our data. So we need to take a little detour and split the mapped reads data into a collection of datasets each containing the mappings for just one segment first again, then perform the consensus construction for all of them in parallel.

> <hands-on-title>Splitting mapped reads by genome segment</hands-on-title>
>
> 1. {% tool [Split BAM by Reference](toolshed.g2.bx.psu.edu/repos/iuc/bamtools_split_ref/bamtools_split_ref/2.5.1+galaxy0) %}
>    - {% icon param-file %} *"BAM dataset to split by reference"*: filtered mapped reads; output of **Samtools view**
>    - *"Select references (chromosomes and contigs) you would like to restrict bam to"*: `Unselect all`
>
{: .hands_on}

The output from this step has the desired collection structure, but the names of the collection elements are not the nicest. Ideally, we would just reuse the segment names, which are already provided in our mapping reference genome. So lets extract these names again and use them as new element labels

> <hands-on-title>Relabeling collection elements</hands-on-title>
>
> 1. {% tool [Select lines that match an expression](Grep1) %}
>    - {% icon param-file %} *"Select lines from"*: the mapping reference; output of **Replace** in Mapping to a hybrid reference
>    - *"that"*: `Matching`
>    - *"the pattern"*: `^>.+`
>
> 2. {% tool [Replace parts of text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_find_and_replace/1.1.4) %}
>    - {% icon param-file %} *"File to process"*: output of **Select**
>    - In {% icon param-repeat %} *"1. Find and Replace"*:
>      - *"Find pattern"*: `^>(.+)$`
>      - *"Replace with"*: `$1`
>      - *"Find-Pattern is a regular expression"*: `Yes`
>
> 3. {% tool [Relabel identifiers](__RELABEL_FROM_FILE__) %}
>    - {% icon param-file %} *"Input collection"*: the split mappings; output of **Split BAM by Reference**
>    - *"How should the new labels be specified?"*: `Using lines in a simple text file.`
>    - {% icon param-file %} *"New Identifiers"*: output of **Replace**
>    - *"Ensure strict mapping"*: `Yes`
>
{: .hands_on}

And with that we are ready for consensus sequence generation!

To accept any base suggested by the mapped sequenced reads as the consensus base for the corresponding genome position, we ask for the following requirements to be fulfilled:

- at least ten sequenced reads have to provide information about the base in question
- at a minimum, 70% of these reads have to agree on the base at this position.

To avoid getting misled too much by sequencing errors, we are also going to ignore bases with a base calling quality less than 20 in the above counts (i.e., we are going to base our decisions only on bases in sequenced reads that the basecaller of the sequencer was reasonably sure about.

Now what if we cannot obtain a consensus base for a position with the above criteria? In such cases of uncertainty we want to insert an N (i.e. an unknown base) to express that we either did not have enough information about the position or that this information was ambiguous.

{% icon galaxy-info %} All of the above limits for consensus base calling are arbitrary to some degree, and depend somewhat on the quality of the sequencing data. With very high overall coverage, for example, it is possible to increase the coverage threshold, but if you increase that threshold too much, you may end up with a consensus sequence consisting mostly of Ns.

> <hands-on-title>Per-segment consensus construction</hands-on-title>
>
> 1. {% tool [ivar consensus](toolshed.g2.bx.psu.edu/repos/iuc/ivar_consensus/ivar_consensus/1.3.1+galaxy0) %}
>    - {% icon param-collection %} *"Bam file"*: the relabeled collection of mapped reads; output of **Relabel identifiers**
>    - *"Minimum quality score threshold to count base"*: `20`
>    - *"Minimum frequency threshold"*: `0.7`
>    - *"Minimum depth to call consensus"*: `10`
>    - *"Exclude regions with smaller depth than the minimum threshold"*: `No`
>    - *"Use N instead of - for regions with less than minimum coverage"*: `Yes`
>
>    The output is a consensus sequence in FASTA format, one per segment, with
>    the names just providing a bit too much detail for our purpose.
>
> 2. {% tool [Replace parts of text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_find_and_replace/1.1.4) %}
>    - {% icon param-collection %} *"File to process"*: consensus sequences; output of **ivar consensus**
>    - In {% icon param-repeat %} *"1. Find and Replace"*:
>      - *"Find pattern"*: `^>Consensus_(.*)_threshold_.*`
>      - *"Replace with"*: `>$1`
>      - *"Find-Pattern is a regular expression"*: `Yes`
>
> 3. {% icon galaxy-eye %} Inspect each consensus sequence generated for the different segments
>
>    > <question-title>Question</question-title>
>    >
>    > Does everything look ok?
>    >
>    > > <solution-title>Answer</solution-title>
>    > >
>    > > As expected from the findings so far, the consensus sequence for the HA segment has stretches of Ns in it, which likely reflect the mapping issues and associated loss of coverage caused by our insufficiently sized collection of references.
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

# Placing segments on a phylogenetic tree

The next logical step after obtaining the consensus sequences of segments of our sample is to explore how those sequences are related to the sequences in our reference collection.
To do so, we are going to build multiple sequence alignments (MSAs) per segment of all reference sequences, add our sample segments to those alignments, and build phylogenetic trees, one per segment. We are going to use two rather standard tools, **MAFFT** and **IQTree**, for generating MSAs and trees, respectively.

> <hands-on-title>Exploring phylogeny</hands-on-title>
>
> 1. {% tool [MAFFT](toolshed.g2.bx.psu.edu/repos/rnateam/mafft/rbc_mafft/7.508+galaxy0) %}
>    - {% icon param-collection %} *"Sequences to align"*: `References per segment (INSAFlu)`
>    - *"Data type"*: `Nucleic Acids`
>    - *"Matrix selection"*: `No matrix`
>
>    The result is a collection of MSAs, each representing all reference sequences of one genome segment.
>
> 2. {% tool [MAFFT add](toolshed.g2.bx.psu.edu/repos/rnateam/mafft/rbc_mafft_add/7.508+galaxy0) %}
>    - {% icon param-collection %} *"Sequences to add to the alignment"*: consensus sequences with simplified names; output of **Replace**
>    - {% icon param-collection %} *"Alignment"*: collection of MSAs; output of MAFFT
>    - *"What do you want to add to the alignment"*: `A single sequence`
>    - *"Keep alignment length"*: `No`
>
>    The result is a new collection with each of our sample consensus sequences added to the respective segment MSA.
>
> 3. {% tool [IQTree](toolshed.g2.bx.psu.edu/repos/iuc/iqtree/iqtree/2.1.2+galaxy2) %}
>    - {% icon param-collection %} *"Specify input alignment file in PHYLIP, FASTA, NEXUS, CLUSTAL or MSF format."*: output of **MAFFT add**
>    - *"Specify sequence type ..."*: `DNA`
>
> 4. {% icon galaxy-eye %} Explore each of the final trees produced by IQTree for the different segments
>
>    > <question-title>Question</question-title>
>    >
>    > What are your conclusions about the sample in general and its HA and NA segments in particular?
>    >
>    > > <solution-title>Answers</solution-title>
>    > >
>    > > - For most of its segments the sample resembles relatively recent (from the last decade) Eurasian reference sequences.
>    > > - For HA and NA the sample clusters with the few available samples of the corresponding subtype.
>    > > - None of the references closest to the sample with respect to HA and NA are close to the recent Eurasian reference cluster for their remaining segments.
>    > > - A plausible explanation is that the H4 and N6 segments of the sample have been brought into the recent Eurasian background through a reassortment event. Caveat: interpretations like this can be heavily influenced by the size of the reference collection!
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

# Conclusion

Analysis workflows for influenza whole-genome sequencing data need to take into account the specific characteristics of the viral genome. Due to their higher natural variability this is especially true for avian influenza samples and for the HA- and NA-encoding segments of the genome.

Nevertheless, it looks possible, with carefully chosen reference segment sequences and bioinformatic tools, to avoid a computationally expensive de-novo assembly approach and to use mapping against a dynamically compiled reference genome instead.

The rather small reference segment collection suggested for this tutorial consists of 56 different samples, of which only a single one has the H4N6 subtype of the sample analyzed here. Still it allowed us to perform subtyping of the sample, to construct complete consensus sequences for 7 segments including NA, and to draw valuable conclusions about the origin of the sample. It is conceivable that a larger collection of references chosen to capture several strains from each HA subtype could solve the remaining issue of the incomplete HA consensus sequence.

