---
layout: tutorial_hands_on

title: "Reproducing Critical Assessment of Metagenome Interpretation assembly challenge on marine dataset with Galaxy, including tool updates and benchmarking analysis of results"
zenodo_link: ""
questions:
  - "How to analyze metagenomics data?"
  - "What information can be extracted of metagenomics data?"
  - "What is the difference between amplicon and shotgun data?"
  - "What are the difference in the analyses of amplicon and shotgun data?"
objectives:
  - "Be familiar with CAMI challenge"
  - "Be able to select one challenge and benchmarking datasets"
  - "Be able to upload benchmarking datasets into Galaxy"
  - "Be familiar with Galaxy computational resources"
  - "Be familiar with existing tools in Galaxy"
  - "Be familiar with databases / reference genomes that are available in Galaxy"
  - "Be able to add/update the tool in Galaxy"
  - "Be able to choose assembling tool based on dataset metadata"
  - "Perform assembly with Flye, Megahit, Abyss, MetaSPAdes tools"
  - "Produce a benchmarking analysis of these assemblies with Quast, Bowtie2, Map with Minimap2, Samtools, and MultiQC tools"
  - "Create plots with Python to compare results"
time_estimation: "2H30M"
key_points:
  - "With amplicon data, we can extract information about the studied community structure"
  - "With shotgun data, we can extract information about the studied community structure and also the functions realised by the community"
  - "The tools used to analyze amplicon and shotgun data are different, except for the visualisation"
  - "Metagenomics data analyses are complex and time-consuming"
contributors:
  - shiltemann
  - bebatut
---

# Introduction
{:.no_toc}

In metagenomics, information about micro-organisms in an environment can be extracted with two main techniques:

- Amplicon sequencing, which sequences only the rRNA or ribosomal DNA of organisms
- Shotgun sequencing, which sequences full genomes of the micro-organisms in the environment

In this tutorial, we will introduce the two main types of analyses with their general principles and differences. For a more in-depth look at these analyses, we recommend our detailed tutorials on each analysis.

We will use two datasets (one amplicon and one shotgun) from the same [project on the Argentinean agricultural pampean soils](https://www.ebi.ac.uk/metagenomics/projects/SRP016633). In this project, three different geographic regions that are under different types of land uses and two soil types (bulk and rhizospheric) were analyzed using shotgun and amplicon sequencing. We will focus on data from the Argentina Anguil and Pampas Bulk Soil (the original study included one more geographical regions, [see](https://doi.org/10.1186/2049-2618-1-21)).

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Amplicon data

Amplicon sequencing is a highly targeted approach for analyzing genetic variation in specific genomic regions.
In the metagenomics fields, amplicon sequencing refers to capture and sequence of rRNA data in a sample.
It can be 16S for bacteria or archea or 18S for eukaryotes.

> ### {% icon comment %} Background: The 16S ribosomal RNA gene
> ![The 16S ribosomal RNA gene](../../images/16S_gene.png "The 16S ribosomal RNA gene. Credit: https://www.slideshare.net/beiko/ccbc-tutorial-beiko")
>
> The 16S rRNA gene has several properties that make it ideally suited to our purposes
>
> 1. Present in all living organisms
> 2. Highly conserved + highly variable regions
> 3. Huge reference databases
>
> ![Variable regions](../../images/16S_variableregions.jpg "Variable regions of the 16S rRNA")
>
> The highly conserved regions make it easy to target the gene across different organisms, while the highly variable regions allow us to distinguish between different species.
>
{: .comment}

With amplicon data, we can determine the micro-organisms from which the sequences in our sample are coming from. This is called taxonomic assignation.
We try to assign sequences to taxons and then classify or extract the taxonomy in our sample.

In this analysis, we will use the [mothur tool suite](https://mothur.org), but only a small portion of its tools and possibilities.
To learn more in detail about how to use this tool, check out the full [mothur tutorial](../mothur-miseq-sop/tutorial.html).


## Quality Control

The first step in any analysis should be to check and improve the quality of our data.


> ### {% icon comment %} Comment
>
> For more information on the topic of quality control, please see our training materials [here]({% link topics/sequence-analysis/index.md %}).
{: .comment}


First, let's get a feel for our data:

> ### {% icon hands_on %} Hands-on: Summarize data
>
> 1. **Summary.seqs** {% icon tool %} with the following parameters
>   - "fasta" parameter to the fasta from `Unique.seqs`
>   - "count" to count table from `Count.seqs`
>   - "output logfile?" to `yes`
>
{: .hands_on}

The `summary` output files give information per read. The `logfile` outputs also contain some summary
statistics:

```
              Start    End        NBases     Ambigs   Polymer  NumSeqs
Minimum:      1        80         80         0        3        1
2.5%-tile:    1        104        104        0        3        501
25%-tile:     1        242        242        0        4        5001
Median:       1        245        245        0        4        10001
75%-tile:     1        245        245        0        4        15001
97.5%-tile:   1        247        247        0        6        19501
Maximum:      1        275        275        2        31       20000
Mean:         1        237.519    237.519    0.00495  4.24965
# of unique seqs:   19502
total # of seqs:    20000
```

This tells us that we have a total of 19,502 unique sequences, representing 20,000 total sequences that vary in length between 80 and 275 bases. Also, note that at least some of our sequences had some ambiguous base calls.
Furthermore, at least one read had a homopolymer stretch of 31 bases, this is likely an error so we would like to filter such reads out as well.

If you are thinking that 20,000 is an oddly round number, you are correct; we downsampled the original datasets to 10,000 reads per sample for this tutorial to reduce the amount of time the analysis steps will take.

We can filter our dataset on length, base quality, and maximum homopolymer length using the `Screen.seqs` tool

The following tool will remove any sequences with ambiguous bases (`maxambig` parameter), homopolymer stretches of 9 or more bases (`maxhomop` parameter) and any reads longer than 275 bp or shorter than 225 bp.

> ### {% icon hands_on %} Hands-on: Filter reads based on quality and length
>
> 1. **Screen.seqs** {% icon tool %} with the following parameters
>   - "fasta" to the fasta file from `Unique.seqs`
>   - "minlength" parameter to `225`
>   - "maxlength" parameter to `275`
>   - "maxambig" parameter to `0`
>   - "maxhomop" parameter to `8`
>   - "count" to the count file from `Count.seqs`
>
{: .hands_on}

> ### {% icon question %} Question
>
> How many reads were removed in this screening step? (Hint: run the `Summary.seqs` tool again)
>
> > ### {% icon solution %} Solution
> > 1,804.
> >
> > This can be determined by looking at the number of lines in bad.accnos output of screen.seqs step or by comparing the total number of seqs between of the summary.seqs log before and after this screening step
> {: .solution }
{: .question}


# Conclusion
{:.no_toc}

We can summarize the analyses with amplicon and shotgun metagenomic data:

![Scheme to sum up the analysis](../../images/general-tutorial-scheme.png)

Both analyses are quite complex! However, in this tutorial, we only showed simple cases of metagenomics data analysis with subset of real data.

Check our other tutorials to learn more in detail of how to analyze metagenomics data.
