---
layout: tutorial_hands_on
topic_name: metagenomics
tutorial_name: general-tutorial
---

# Introduction
{:.no_toc}

In metagenomics, information about micro-organisms in an environment can be extracted with two main techniques:

- Amplicon sequencing (or 16S rRNA/rDNA), which sequence only on the rRNA/rDNA of organisms
- Shotgun sequencing, which sequence full genomes of the micro-organisms in the environment

In this tutorial, we will introduce the two types of analyses with the general principles behind and the differences. To go deeper in such analyses, we recommend to check our detailed tutorials on each analysis.

For that, we will use two datasets (one amplicon and one shotgun) from the same environment: the Argentina Anguil Bulk Soil, studied in a [project on the Argentinean agricultural pampean soils](https://www.ebi.ac.uk/metagenomics/projects/SRP016633). In this project, three different types of land uses and two soil types (bulk and rhizospheric) were analyzed using shotgun and amplicon sequencing. We will focus on the Argentina Anguil Bulk Soil.

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

> ### {% icon tip %} Background: The 16S ribosomal RNA gene
> ![The 16S ribosomal RNA gene](../../images/16S_gene.png) <br><br>
>
> The 16S rRNA gene has several properties that make it ideally suited for our purposes
>
> 1. Present in all living organisms
> 2. Single copy (no recombination)
> 3. Highly conserved + highly variable regions
> 4. Huge reference databases
>
> ![Variable regions](../../images/16S_variableregions.jpg)
>
> The highly conserved regions make it easy to target the gene across different organisms,
> while the highly variable regions allow us to distinguish between different species.
>
> (slide credit [https://www.slideshare.net/beiko/ccbc-tutorial-beiko ](http://www.slideshare.net/beiko/ccbc-tutorial-beiko ))
{: .tip}

With amplicon data, we can extract from which micro-organisms the sequences in our sample are coming from. This is called taxonomic assignation.
We try to assign sequences to taxons and then classify or extract the taxonomy in our sample.

In this analysis, we will use [mothur tool suite](https://mothur.org), but only a small portion of its tools and possibilities.
To learn more in detail how to use, check out the full [mothur tutorial](../mothur-miseq-sop/tutorial.html).

## Importing the data

Our datasets comes from a soil samples in two different Argentinian locations, with capture and sequencing of the 16S rDNA V4 region
using 454 GS FLX Titanium. The original data are available at EBI Metagenomics under the following run numbers:

- Pampa soil: [SRR531818](https://www.ebi.ac.uk/metagenomics/projects/SRP016633/samples/SRS353016/runs/SRR531818/results/versions/2.0)
- Anguil soil: [SRR651839](https://www.ebi.ac.uk/metagenomics/projects/SRP016633/samples/SRS386929/runs/SRR651839/results/versions/2.0)


> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Import from [Zenodo](https://zenodo.org/record/815875) or from the data library (in "Analyses of metagenomics data") the files
>    - `SRR531818_pampa.fasta`
>    - `SRR651839_anguil.fasta`
>
>    > ### {% icon tip %} Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * Press **Start**
>    {: .tip}
>
>    > ### {% icon tip %} Tip: Importing data from a data library
>    >
>    > * Go into "Shared data" (top panel) then "Data libraries"
>    > * Click on "Training data" and then "Analyses of metagenomics data"
>    > * Select interesting file
>    > * Click on "Import selected datasets into history"
>    > * Import in a new history
>    {: .tip}
>
>    As default, Galaxy takes the link as name, so rename them.
>
{: .hands_on}

<!--

Anguil Soil: https://www.ebi.ac.uk/metagenomics/projects/SRP016633/samples/SRS386929
Pampa Soil https://www.ebi.ac.uk/metagenomics/projects/SRP016633/samples/SRS353016/runs/SRR531818/results/versions/2.0

Project's data: https://www.ebi.ac.uk/metagenomics/projects/SRP016633/samples/SRS353016/runs/SRR531818/results/versions/2.0
Project's pipeline: https://www.ebi.ac.uk/metagenomics/pipelines/2.0
Project's QC results: https://www.ebi.ac.uk/metagenomics/projects/SRP016633/samples/SRS386929/runs/SRR651839/results/versions/2.0
-->

### Preparing datasets

We will perform a multisample analysis with mothur, in order to do so, we will merge all reads into a single file,
and create a *group file*, indicating which reads belong to which samples.

> ### {% icon hands_on %} Hands-on: Prepare multisample analysis
>
> 1. **Merge.files** {% icon tool %} with the following parameters
>   - "Merge" to `fasta files`
>   - "Inputs" to the two sample fasta files
>
> 2. **Make.group** {% icon tool %} with the following parameters
>   - "Method" to `Manually specify fasta files and group names`
>   - "Additional": Add two elements to this repeat: Pampa sample fasta file, with group name `pampa` and Anguil sample fasta file, with group name `anguil`
>
{: .hands_on}

> ### {% icon tip %} Tip
>
> Because we only have a small number of samples, we used the manual specification. If you have hundreds of samples this would quickly become bothersome. The solution? use a collection! To read more about collections in Galaxy, please see [this]() tutorial
{: .tip}

Have a look at the group file. It is a very simple file, it contains two columns, first contains the read names, second the group (sample) name, in our case `pampa` or `anguil`.


### Optimization of files for computation

Because we are sequencing many of the same organisms, we anticipate that many of our sequences are
duplicates of each other. Because it's computationally wasteful to align the same thing a bazillion
times, we'll unique our sequences using the `Unique.seqs` command:

> ### {% icon hands_on %} Hands-on: Remove duplicate sequences
>
> 1. **Unique.seqs** {% icon tool %} with the following parameters
>   - "fasta" to the merged fasta file
>
>    > ### {% icon question %} Question
>    >
>    > How many sequences were unique? How many duplicates were removed?
>    >
>    >    <details>
>    >    <summary>Click to view answer</summary>
>    >    9,502 unique sequences and 498 duplicates. <br>
>    >    This can be determined from the number of lines in the fasta (or names) output, compared to the
>    >    number of lines in the fasta file before this step.
>    >    </details>
>    {: .question}
>
{: .hands_on}

This `Unique.seqs` tool outputs two files, one is a fasta file containing only the unique sequences, and a *names files*.
The names file consists of two columns, the first contains the sequence names for each of the unique sequences, and the second column contains all other sequence names that are identical to the representative sequence in the first column.

```
name          representatives
read_name1    read_name2,read_name,read_name5,read_name11
read_name4    read_name6,read_name,read_name10
read_name7    read_name8
...
```

> ### {% icon hands_on %} Hands-on: Count sequences
>
> 1. **Count.seqs** {% icon tool %} with the following parameters
>   - "name" to the name file from `Unique.seqs`
>   - "Use a group file" to `yes`
>   - "group" to the group file from `Make.group`
{: .hands_on}

The `Count.seqs` file keeps track of the number of sequences represented by each unique representative across multiple samples. We will pass this file to many of the following tools to be used or updated as needed.

## Quality Control

The first step in any analysis should be to check and improve the quality of our data.


> ### {% icon comment %} Comment
>
> For more information on the topic of quality control, please see our training materials [here]({{site.url}}/topics/sequence-analysis/).
{: .comment}


First, let's get a feel of our data:

> ### {% icon hands_on %} Hands-on: Summarize data
>
> 1. **Summary.seqs** {% icon tool %} with the following parameters
>   - "fasta" parameter to the fasta from `Unique.seqs`
>   - "count" to count table from `Count.seqs`
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

If you are thinking that 20,000 is an oddly round number, you are correct, we downsampled the original datasets to 10,000 reads per sample for this tutorial to reduce the amount of time the analysis steps will take.

We can filter our dataset on length, base quality, and maximum homopolymer length using the `Screen.seqs` tool

The following tool will remove any sequences with ambiguous bases and anything longer than 275 bp.

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
> > ### {% icon question %} Question
> >
> > How many reads were removed in this screening step? (Hint: run the `Summary.seqs` tool again)
> >
> >    <details>
> >    <summary>Click to view answer</summary>
> >    1,804. <br>
> >    This can be determined by looking at the number of lines in bad.accnos output of screen.seqs step or by comparing the total number of seqs between of the summary.seqs log before and after this screening step
> >    </details>
> {: .question}
{: .hands_on}

## Sequence Alignment

Aligning our sequences to a reference helps improve OTU assignment [[Schloss et. al.](https://www.ncbi.nlm.nih.gov/pubmed/23018771)], so we will now align our sequences to the Silva reference database.

> ### {% icon hands_on %} Hands-on: Align sequences
>
> 1. Import the `silva.v4.fasta` file in your history
> 2. **Align.seqs** {% icon tool %} with the following parameters
>   - "fasta" to the `good.fasta` output from `Screen.seqs`
>   - "reference" to the `silva.v4.fasta` reference file
>   - "flip" to `Yes`
>
>    This step may take a few minutes, please be patient.
>
> 3. **Summary.seqs** {% icon tool %} with the following parameters
>   - "fasta" parameter to the aligned output from `Align.seqs`
>   - "count" parameter to count_table output from `Screen.seqs`
>
{: .hands_on}

To get an idea of the quality of the alignment, we can view the log output from the summary step:

```
        Start   End NBases  Ambigs  Polymer NumSeqs
Minimum:    2391    10674   9   0   2   1
2.5%-tile:  3080    12071   234 0   4   455
25%-tile:   3080    13424   244 0   4   4545
Median:     3080    13424   245 0   4   9090
75%-tile:   3080    13424   245 0   4   13634
97.5%-tile: 3082    13424   246 0   6   17724
Maximum:    13396   13425   267 0   7   18178
Mean:   3080.6  13380   244.212 0   4.27946
# of unique seqs:   17698
total # of seqs:    18178
```

> ### {% icon question %} Questions
>
> 1. How many sequences have been aligned?
> 2. Between which positions most of the reads are aligned to this references?
>
>    <details>
>    <summary>Click to view answers</summary>
>    <ol type="1">
>    <li>17,698 are aligned</li>
>    <li>From this we can see that most of our reads align nicely to positions `3080-13424` on this reference.
This corresponds exactly to the V4 target region of the 16S gene.</li>
>    </ol>
>    </details>
{: .question}

To make sure that everything overlaps the same region we'll re-run `Screen.seqs` to get sequences that start at or before position 3,080 and end at or after position 13,424.

> ### {% icon hands_on %} Hands-on: Remove poorly aligned sequences
>
> 1. **Screen.seqs** {% icon tool %} with the following parameters
>   - "fasta" to the aligned fasta file
>   - "start" to `3080`
>   - "end" to `13424`
>   - "count" to the group file created by the previous run of `Screen.seqs`
>
> > ### {% icon question %} Question
> >
> >  How many sequences were removed in this step?
> > <details>
> >   <summary> Click to view answer</summary>
> >   4,579 sequences were removed. This is the number of lines in the bad.accnos output.
> > </details>
> {: .question}
{: .hands_on}

Now we know our sequences overlap the same alignment coordinates, we want to make sure they *only* overlap that region. So we'll filter the sequences to remove the overhangs at both ends. In addition, there are many columns in the alignment that only contain gap characters (*i.e.* "."). These can be pulled out without losing any information. We'll do all this with `Filter.seqs`:

> ### {% icon hands_on %} Hands-on: Filter sequences
>
> 1. **Filter.seqs** {% icon tool %} with the following parameters
>   - "fasta"" to `good.fasta` output from `Screen.seqs`
>   - "vertical" to Yes
>   - "trump" to `.`
{: .hands_on}


## Extraction of taxonomic information

The main questions when analyzing amplicon data are: Which micro-organisms are present in an environmental samples? And in which proportion? What is the structure of the community of the micro-organisms?

The idea is to take the sequences and assign them to a taxon. To do that, we group (or cluster) sequences based on their similarity to define Operational Taxonomic Units (OTUs); groups of similar sequences that can be treated as a single "genus" or "species" (depending on the clustering threshold)

> ### {% icon tip %} Background: Operational Taxonomic Units (OTUs)
>
> In 16S metagenomics approaches, OTUs are clusters of similar sequence variants of the 16S rDNA marker gene sequence. Each of these clusters is intended to represent a taxonomic unit of a bacteria species or genus depending on the sequence similarity threshold. Typically, OTU cluster are defined by a 97% identity threshold of the 16S gene sequence variants at genus level. 98% or 99% identity is suggested for species separation.
>
> ![OTU and cluster with 97% identity threshold](../../images/otu.png)
>
> ![OTU graph](../../images/OTU_graph.png)
>
> (Image credit: Danzeisen et al. 2013, 10.7717/peerj.237)
{: .tip}



The first thing we want to do is to further de-noise our sequences, by pre-clustering the sequences using the `Pre.cluster` command, allowing for up to 2 differences between sequences. This command will split the sequences by group and then sort them by abundance and go from most abundant to least and identify sequences that differ no more than 2 nucleotides from on another. If this is the case, then they get merged. We generally recommend allowing 1 difference for every 100 basepairs of sequence:

> ### {% icon hands_on %} Hands-on: Perform preliminary clustering of sequences and remove undesired sequences
>
> 1. **Pre.cluster** {% icon tool %} with the following parameters
>   - "fasta" to the fasta output from the last `Filter.seqs` run
>   - "name file or count table" to the count table from the last `Screen.seqs` step
>   - "diffs" to 2
>
>   > ### {% icon question %} Question
>   >
>   >  How many unique sequences are we left with after this clustering of highly similar sequences?
>   > <details>
>   >   <summary> Click to view answer</summary>
>   >   10,386 <br>
>   >   This is the number of lines in the fasta output
>   > </details>
>   {: .question}
>
{: .hands_on}

<!-- optional additional QC: chimera.uchime -->
We would like to classify the sequences using a training set.

> ### {% icon hands_on %} Hands-on: Classify the sequences
>
> 1. Import the `trainset16_022016.pds.fasta` and `trainset16_022016.pds.tax` in your history
> 2. **Classify.seqs** {% icon tool %} with the following parameters
>   - "fasta" to the fasta output from `Pre.cluster`
>   - "reference" to `trainset16_022016.pds.fasta` from your history
>   - "taxonomy" to `trainset16_022016.pds.tax` from your history
>   - "cutoff" to 80
>   - "count" to the count table from `Pre.cluster`
>
> This step may take a couple of minutes, now may be a good time to grab a cup of tea :coffee:
>
{: .hands_on}

Have a look at the taxonomy output. You will see that every read now has a classification.

The next step is then to use this information to know the abundance of the different found taxons.

> ### {% icon hands_on %} Hands-on: Cluster our data into OTUs
>
> 1. **Cluster.split** {% icon tool %} with the following parameters
>   - "Split by" to `Classification using fasta`
>   - "fasta" to the fasta output from `Pre.cluster`
>   - "taxonomy" to the taxonomy output from `Classify.seqs`
>   - "count" to the count table output from `Pre.cluster`
>   - "cutoff" to `0.15`
>
>     Next we want to know how many sequences are in each OTU from each group and we can do this using the `Make.shared` command. Here we tell mothur that we're really only interested in the 0.03 cutoff level:
>
> 2. **Make.shared** {% icon tool %} with the following parameters
>   - "Select input type" to `OTU list`
>   - "list" to list output from `Cluster.split`
>   - "count" to the count table from `Pre.cluster`
>   - "label" to `0.03`
>
>     We probably also want to know the taxonomy for each of our OTUs. We can get the consensus taxonomy for each OTU using the `Classify.otu` command:
>
> 3. **Classify.otu** {% icon tool %} with the following parameters
>   - "list" to output from `Cluster.split`
>   - "count" to the count table from `Pre.cluster`
>   - "taxonomy" to the taxonomy output from `Classify.seqs`
>   - "label" to `0.03`
{: .hands_on}

> ### {% icon question %} Questions
>
> How many OTUs with taxonomic assignation are found for the Anguil sample? And for the Pampa sample?
>
>    <details>
>    <summary>Click to view answers</summary>
>    2,212 for Anguil and 2,490 for Pampa ("tax.summary" output of `Classify.otus`)
>    </details>
{: .question}

## Visualization

We have now determined our OTUs and classified them, but looking at a long text file is not very informative.
Let's visualize our data using Krona:

> ### {% icon hands_on %} Hands-on: Krona
>
> 1. **Visualize with Krona** {% icon tool %} with the following parameters
>   - "Input file" to taxonomy output from `Classify.otu` (collection)
>   - Set **Is this output from mothur?** to `Yes`
>
{: .hands_on}

The result is an HTML file with an interactive visualization, for instance try clicking
on one of the rings in the image or playing around with some of the settings.

![Krona output](../../images/krona.png)

This produced a single plot for both your samples, but what if you want to compare
the two samples?

> ### {% icon hands_on %} Hands-on: Per-sample Krona plots
>
> 1. **Classify.otu** {% icon tool %}
>
>    Hit the rerun button on the `Classify.otu` job in your history and see if you can find settings that will give you per-sample taxonomy data
>
> 2. **Visualize with Krona** {% icon tool %}
>
>    Now use this new output collection to create per-sample Krona plots
>
{: .hands_on}


To further explore the community structure, we can visualize it with dedicated tools such as Phinch:

> ### {% icon hands_on %} Hands-on: Visualization of the community structure with Phinch
>
> 1. **Make.biom** {% icon tool %} with the following parameters
>   - "shared" to Make.shared
>   - "constaxonomy" to taxonomy output from the first run of `classify.otu` (collection)
>
>    The Galaxy project runs an instance of Phinch, and if you look at the output BIOM file, you will see a link to view the file at Phinch:
>
>    ![Link to Phinch](../../../../shared/images/viewatphinch.png)
>
> 2. Click on the icon
>
>    It will lead you to the Phinch website, which will automatically load in your file, and where you can several interactive visualisations:
>
>     ![Phinch website interface](../../../../shared/images/phinch_overviewpage.png)
{: .hands_on}

Once we have information about the community structure (OTUs with taxonomic structure), we can do more analysis on it: estimation of the diversity of micro-organism, comparison fo diversity between samples, analysis of populations, ... We will not go into detail of such analyses here but you follow our tutorials on amplicon data analyses to learn about them.

# Shotgun metagenomics data

In the previous section, we see how to analyze amplicon data to extract the community structure. Such information can also be extracted from shotgun metagenomic data.

In shotgun data, full genomes of the micro-organisms in the environment are sequenced (not only the 16S or 18S). We can then have access to the rRNA (only a small part of the genomes), but also to the genes of the micro-organisms. Using this information, we can try to answer to questions "What are the micro-organisms doing?" in addition to the question "What micro-organisms are present?".

In this second part, we will use a metagenomic sample of the Pampas Soil ([SRR606451](https://www.ebi.ac.uk/metagenomics/projects/SRP016633/samples/SRS372043/runs/SRR606451/results/versions/2.0)).

## Data upload

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history
> 2. Import the "SRR606451_pampa" Fasta file from [Zenodo](http://zenodo.org/record/815875) or from the data library (in "Analyses of metagenomics data")
{: .hands_on}

## Extraction of taxonomic information

As for amplicon data, we can extract taxonomic and community structure information from shotgun data. Different approaches can be used:

- Same approaches as for amplicon data with identification and classification of OTUs

    Such approaches imply a first step of sequence sorting to extract only the 16S and 18S sequences on which the same tools as for amplicon data. However, rRNA sequences represent a low proportion (< 1%) of the shotgun sequences so such an approach is not the most statistically supported

- Assignation of taxonomy on the whole sequences using databases with marker genes

In this tutorial, we use the second approach with MetaPhlAn2. This tools is using a database of ~1M unique clade-specific marker genes (not only the rRNA genes) identified from ~17,000 reference (bacterial, archeal, viral and eukaryotic) genomes.

> ### {% icon hands_on %} Hands-on: Taxonomic assignation with MetaPhlAn2
>
> 1. **MetaPhlAN2** {% icon tool %} with
>    - "Input file" to the imported file`
>    - "MetaPhlAn2 clade-specific marker genes" to `locally cached`
>    - "Cached database with clade-specific marker genes" to `MetaPhlAn2 clade-specific marker genes`
>
> This step may take a couple of minutes :coffee:
{: .hands_on}

3 files are generated:

- A tabular file with the community structure

    ```
    #SampleID   Metaphlan2_Analysis
    k__Bacteria 100.0
    k__Bacteria|p__Proteobacteria   86.20712
    k__Bacteria|p__Actinobacteria   13.79288
    k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria    86.20712
    k__Bacteria|p__Actinobacteria|c__Actinobacteria 13.79288
    ```

    Each line contains a taxa and its relative abundance found for our sample. The file starts with high level taxa (kingdom: `k__`) and go to more precise taxa.


- A BIOM file with the same information as the previous file but in BIOM format

    It can be used then by mothur and other tools requiring community structure information in BIOM format

- A SAM file with the results of the mapping of the sequences on the reference database

> ### {% icon question %} Questions
>
> 1. What is the most precise level we have access to with MetaPhlAn2?
> 2. What are the two orders found in our sample?
> 3. What is the most abundant family in our sample?
>
>    <details>
>    <summary>Click to view answers</summary>
>    <ol type="1">
>    <li>We have access to species level</li>
>    <li>Pseudomonadales and Solirubrobacterales are found in our sample</li>
>    <li>The most abundant family is Pseudomonadaceae with 86.21 % of the assigned sequences</li>
>    </ol>
>    </details>
{: .question}

Even if the output of MetaPhlAn2 is bit easier to parse than the BIOM file, we want to visualize and explore the community structure with KRONA

> ### {% icon hands_on %} Hands-on: Interactive visualization with KRONA
>
> 1. **Format MetaPhlAn2 output for Krona** {% icon tool %} with
>    - "Input file" to `Community profile` output of `MetaPhlAn2`
>
> 2. **KRONA pie chart** {% icon tool %}: with
>    - "What is the type of your input data" as `MetaPhlan`
>    - "Input file" to the output of `Format MetaPhlAn2`
>
{: .hands_on}

## Extraction of functional information

We would like now to answer the question "What are the micro-organisms doing?" or "Which functions are done by the micro-organisms in the environment?".

In the shotgun data, we have access to the sequences from the full genome, with gene sequences then. We use that to identify the genes, associate them to a function, build pathways, etc to investigate the functional part of the community.

> ### {% icon hands_on %} Hands-on: Metabolism function identification
>
> 1. **HUMAnN2** {% icon tool %} with
>    - "Input sequence file" to the imported sequence file
>    - "Use of a custom taxonomic profile" to `Yes`
>    - "Taxonomic profile file" to `Community profile` output of `MetaPhlAn2`
>    - "Nucleotide database" to `Locally cached`
>    - "Nucleotide database" to `Full`
>    - "Protein database" to `Locally cached`
>    - "Protein database" to `Full UniRef50`
>    - "Search for uniref50 or uniref90 gene families?" to `uniref50`
>    - "Database to use for pathway computations" to `MetaCyc`
>    - "Advanced Options"
>    - "Remove stratification from output" to `Yes`
>
>    This step is long so we generated the output for you!
>
> 2. Import the 3 files whose the name is starting with "humann2"
{: .hands_on}

HUMAnN2 generates 3 files

- A file with the abundance of gene families

    Gene family abundance is reported in RPK (reads per kilobase) units to normalize for gene length. It reflects the relative gene (or transcript) copy number in the community.

    "UNMAPPED" value is the total number of reads which remain unmapped after both alignment steps (nucleotide and translated search). Since other gene features in the table are quantified in RPK units, "UNMAPPED" can be interpreted as a single unknown gene of length 1 kilobase recruiting all reads that failed to map to known sequences.

- A file with the coverage of pathways

    Pathway coverage provides an alternative description of the presence (1) and absence (0) of pathways in a community, independent of their quantitative abundance.

- A file with the abundance of pathways

> ### {% icon question %} Questions
>
> How many gene families and pathways have been identified?
>
>    <details>
>    <summary>Click to view answers</summary>
>    44 gene families but no pathways are identified
>    </details>
{: .question}

The RPK for the gene families are quite difficult to interpret in term of relative abundance. We decide then to normalize the values

> ### {% icon hands_on %} Hands-on: Normalize the gene family abundances
>
> 1. **Renormalize a HUMAnN2 generated table** {% icon tool %} with
>    - "Gene/pathway table" to the gene family table generated with `HUMAnN2`
>    - "Normalization scheme" to `Relative abundance`
>    - "Normalization level" to `Normalization of all levels by community total`
>
>  > ### {% icon question %} Questions
>  >
>  > 1. Which percentage of sequences has not be assigned to a gene family?
>  > 2. What is the most abundant gene family?
>  >
>  >    <details>
>  >    <summary>Click to view answers</summary>
>  >    <ol type="1">
>  >    <li>55% of the sequences has not be assigned to a gene family</li>
>  >    <li>The most abundant gene family with 25% of sequences is a putative secreted protein</li>
>  >    </ol>
>  >    </details>
>  {: .question}
{: .hands_on}

With the HUMAnN2 output, we have access to UniRef50 gene families. However, the names can remains cryptic and sometimes we would like a more general view about the functions. HUMAnN proposes a tool to regroup the gene families into different meta-groups: GO (Gene Ontology), EC, etc.

> ### {% icon hands_on %} Hands-on: Regroup the gene families into GO terms
>
> 1. **Regroup a HUMAnN2 generated table by features** {% icon tool %} with
>    - "Gene/pathway table" to the gene family table generated with `HUMAnN2`
>    - "How to combine grouped features?" to `Sum`
>    - "Use built-in grouping options?" to `Yes`
>    - "Gene family type" to `UniRef50 gene families`
>    - "Grouping options" to `UniRef50 gene families into GO`
>
> 2. **Renormalize a HUMAnN2 generated table** {% icon tool %} with
>    - "Gene/pathway table" to the gene family table generated with `Regroup`
>    - "Normalization scheme" to `Relative abundance`
>    - "Normalization level" to `Normalization of all levels by community total`
>
>  > ### {% icon question %} Questions
>  >
>  > 1. What is the most abundant GO term?
>  > 2. What is related to in [Gene Ontology](http://www.geneontology.org/)?
>  >
>  >    <details>
>  >    <summary>Click to view answers</summary>
>  >    <ol type="1">
>  >    <li>GO:0007275 (found after sorting) with 2.68%</li>
>  >    <li>It seems to correspond to "multicellular organism development"(?)</li>
>  >    </ol>
>  >    </details>
>  {: .question}
{: .hands_on}

With the previous analyses, we investigate "Which micro-organims are present in my sample?" and "What function are done by the micro-organisms in my sample?". We can go further in these analyses (for example with combination of functional and taxonomic results). We did not detail that in this tutorial but you can found more analyses in our tutorials on shotgun metagenomic data analyses.

# Conclusion
{:.no_toc}

We can summarize the analyses with amplicon and shotgun metagenomic data:

![Scheme to sum up the analysis](../../images/general-tutorial-scheme.png)

Both analyses are quite complex! However, in this tutorial, we only showed simple cases of metagenomics data analysis with subset of real data.

Check our other tutorials to learn more in details how to analyze metagenomics data.
