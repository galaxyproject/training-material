---
layout: tutorial_hands_on
topic_name: metagenomics
tutorial_name: general-tutorial
---

# Introduction

In metagenomics, information about micro-organisms in an environment can be extracted with two main techniques:

- Amplicon sequencing (or 16S rRNA/rDNA), which sequence only on the rRNA/rDNA of organisms
- Whole-genome sequencing (WGS), which sequence full genomes of the micro-organisms in the environment

In this tutorial, we will introduce the two types of analyses with the general principles behind and the differences. To go deeper in such analyses, we recommend to check our detailed tutorials on each analysis.

For that, we will use two datasets (one amplicon and one WGS) from the same environment: the Argentina Anguil Bulk Soil, studied in a [project on the Argentinean agricultural pampean soils](https://www.ebi.ac.uk/metagenomics/projects/SRP016633). In this project, three different types of land uses and two soil types (bulk and rhizospheric) were analyzed using WGS and amplicon sequencing. We will focus on the Argentina Anguil Bulk Soil.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. [Amplicon data](#amplicon-data)
> 2. [Whole-genome sequencing data](#whole-genome-sequencing-data)
{: .agenda}

# Amplicon data

Amplicon sequencing is a highly targeted approach for analyzing genetic variation in specific genomic regions.
In the metagenomics fields, amplicon sequencing refers to capture and sequence of rRNA data in a sample.
It can be 16S for bacteria or archea or 18S for eukaryotes.

> ### :book: Background: The 16S ribosomal RNA gene
> ![](../../images/16S_gene.png) <br><br>
>
> The 16S rRNA gene has several properties that make it ideally suited for our purposes
>
> 1. Present in all living organisms
> 2. Single copy (no recombination)
> 3. Highly conserved + highly variable regions
> 4. Huge reference databases
>
> ![](../../images/16S_variableregions.jpg)
>
> The highly conserved regions make it easy to target the gene across different organisms,
> while the highly variable regions allow us to distinguish between different species.
>
> (slide credit [http://www.slideshare.net/beiko/ccbc-tutorial-beiko ](http://www.slideshare.net/beiko/ccbc-tutorial-beiko ))
{: .tip}

With amplicon data, we can extract from which micro-organisms the sequences in our sample are coming from. This is called taxonomic assignation.
We try to assign sequences to taxons and then classify or extract the taxonomy in our sample.

Our datasets comes from a soil samples in two different Argentinian locations, with capture and sequencing of the 16S rDNA V4 region
using 454 GS FLX Titanium. The original data are available at EBI Metagenomics under the following run numbers:

Pampa soil: [SRR531818](https://www.ebi.ac.uk/metagenomics/projects/SRP016633/samples/SRS353016/runs/SRR531818/results/versions/2.0) and
Anguil soil: [SRR651839](https://www.ebi.ac.uk/metagenomics/projects/SRP016633/samples/SRS386929/runs/SRR651839/results/versions/2.0)

In this analysis, we will use [Mothur tool suite](http://mothur.org), but only a small portion of its tools and possibilities.
To learn more in detail how to use, check out the full [Mothur tutorial](../mothur-miseq-sop/tutorial.html).

## Importing the data

> ### :pencil2: Hands-on: Data upload
>
> 1. Import the FASTQ file from [Zenodo]() or from the data library (in "Analyses of metagenomics data" the "..." file)
>
>    > ### :bulb: Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * Press **Start**
>    {: .tip}
>
>    > ### :bulb: Tip: Importing data from a data library
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

> ### :pencil2: Hands-on: prepare multisample analysis
>
> - **merge.files** :wrench: with the following parameters
>   - "Merge" to `fasta files`
>   - "Inputs" to the two sample fasta files
>
> - **make.group** :wrench: with the following parameters
>   - "Method" to Manually
>   - "Additional" Add two elements to this repeat
>     - Pampa sample fasta file, with group name `pampa`
>     - Anguil sample fasta file, with group name `anguil`
>
{: .hands_on}

> ### :bulb: Tip
>
> Because we only have a small number of samples, we used the manual specification. If you have hundreds of samples
> this would quickly become bothersome. The solution? use a collection! To read more about collections in Galaxy
> please see [this]() tutorial
{: .tip}

Have a look at the group file. It is a very simple file, it contains two columns, first contains the read names,
second the group (sample) name, in our case `pampa` or `anguil`.


### Optimize files for computation
Because we are sequencing many of the same organisms, we anticipate that many of our sequences are
duplicates of each other. Because it's computationally wasteful to align the same thing a bazillion
times, we'll unique our sequences using the `unique.seqs` command:

> ### :pencil2: Hands-on: Remove duplicate sequences
>
> - **Unique.seqs** :wrench: with the following parameters
>   - "fasta" to the `good.fasta` output from Screen.seqs
>
>
> > ### :question: Question
> >
> > How many sequences were unique? how many duplicates were removed?
> >
> >    <details>
> >    <summary>Click to view answer</summary>
> >    16,426 unique sequences and 112,446 duplicates. <br>
> >    This can be determined from the number of lines in the fasta (or names) output, compared to the
> >    number of lines in the fasta file before this step. The log file also contains a line showing the
> >    total number of sequences before and the command: <br><br>
> >    mothur > unique.seqs(fasta=fasta.dat) <br>
> >    128872	16426
> >    </details>
> {: .question}
{: .hands_on}

This tool outputs two files, one is a fasta file containing only the unique sequences, and a *names files*.
The names file consists of two columns, the first contains the sequence names for each of the unique
sequences, and the second column contains all other sequence names that are identical to the representative
sequence in the first column.

```
name          representatives
read_name1    read_name2,read_name,read_name5,read_name11
read_name4    read_name6,read_name,read_name10
read_name7    read_name8
...
```


## Quality Control

The first step in any analysis should be to check and improve the quality of our data.
For more information on the topic of quality control, please see our training materials
[here](https://galaxyproject.github.io/training-material/NGS-QC/)

First, let's get a feel of our data:

> ### :pencil2: Hands-on: Summarize data
>
> - **Summary.seqs** :wrench: with the following parameters
>   - "fasta" parameter to the fasta file you just imported.
>   - We do not need to supply a names or count file
>
{: .hands_on}

The `summary` output files give information per read. The `logfile` outputs also contain some summary
statistics:

```
            Start  End      NBases    Ambigs   Polymer  NumSeqs
Minimum:    1      80       80        0        3        1
2.5%-tile:  1      104      104       0        3        501
25%-tile:   1      242      242       0        4        5001
Median:     1      245      245       0        4        10001
75%-tile:   1      245      245       0        4        15001
97.5%-tile: 1      247      247       0        6        19501
Maximum:    1      275      275       2        31        20000
Mean:       1      237.519  237.519   0.00495  4.24965
# of Seqs:      20000
```

This tells us that we have a total of 20,000 sequences that vary in length between 80 and 275 bases.
Also, note that at least some of our sequences had some ambiguous base calls. Furthermore, at least one
read had a homopolymer stretch of 31 bases, this is likely an error so we would like to filter such reads
out as well.

If you are thinking that 20,000 is an oddly round number, you are correct, we downsampled the original
datasets to 10,000 reads each for this tutorial to reduce the amount of time the analysis steps will take.

We can filter our dataset on length, base quality, and maximum homopolymer length using the `screen.seqs` tool

The following tool will remove any sequences with ambiguous bases and anything longer than 275 bp.

> ### :pencil2: Hands-on: Filter reads based on quality and length
>
> - **Screen.seqs** :wrench: with the following parameters
>   - "fasta" to the merged fasta file
>   - "group" the group file created in the make.contigs step
>   - "minlength" parameter to `225`
>   - "maxlength" parameter to `275`
>   - "maxambig" parameter to `0`
>   - "maxhomop" parameter to `8`
>   - "group" to the group file we created
>
> > ### :question: Question
> >
> > How many reads were removed in this screening step? (Hint: run the summary.seqs tool again)
> >
> >    <details>
> >    <summary>Click to view answer</summary>
> >    1,822. <br>
> >    This can be determined by looking at the number of lines in bad.accnos output of screen.seqs step
> >    or by comparing the total number of seqs between of the summary.seqs log before and after this screening
> >    step
> >    </details>
> {: .question}
{: .hands_on}



## Sequence Alignment

Aligning our sequences to a reference helps improve OTU assignment [[Schloss et. al.](https://www.ncbi.nlm.nih.gov/pubmed/23018771)],
so we will now align our sequences to the Silva reference database.

> ### :pencil2: Hands-on: Align sequences
>
> - **Align.seqs** :wrench: with the following parameters
>   - "fasta" to the `good.fasta` output from screen.seqs
>   - "reference" to the `silva.v4.fasta` reference file from your history
>   - "flip" to `Yes`
>
> This step may take a few minutes, please be patient.
>
> - **Summary.seqs** :wrench: with the following parameters
>   - "fasta" parameter to the aligned output from previous step
>   - "count" parameter to `count_table` output from Count.seqs
>
{: .hands_on}

View the log output from the summary step.

```
              Start      End        NBases   Ambigs   Polymer  NumSeqs
Minimum:      1044       1065       10       0        2        1
2.5%-tile:    14974      23965      234      0        4        455
25%-tile:     14974      25318      244      0        4        4545
Median:       14974      25318      245      0        4        9090
75%-tile:     14974      25318      245      0        4        13634
97.5%-tile:   14976      25318      247      0        6        17724
Maximum:      15630      26169      274      0        7        18178
Mean:         14973.3    25277.7    244.314  0        4.2799
# of Seqs:      18178


		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	138	152	10	0	2	1
2.5%-tile:	2308	4062	234	0	4	455
25%-tile:	2308	4090	244	0	4	4545
Median: 	2308	4090	245	0	4	9090
75%-tile:	2308	4090	245	0	4	13634
97.5%-tile:	2310	4090	247	0	6	17724
Maximum:	2323	4152	274	0	7	18178
Mean:	2307.95	4087.47	244.314	0	4.2799
# of Seqs:	18178
```

2308 4090

From this we can see that most of our reads align nicely to positions `14974-25318` on this reference.
This corresponds exactly to the V4 target region of the 16S gene. If our data was not so nice, we could
now remove any sequences not mapped to our target region using `screen.seqs`


To make sure that everything overlaps the same region we'll re-run screen.seqs to get sequences that
start at or before position 1968 and end at or after position 11550. We'll also set the maximum
homopolymer length to 8 since there's nothing in the database with a stretch of 9 or more of the same
base in a row (this also could have been done in the first execution of screen.seqs above).

> ### :pencil2: Hands-on: Remove poorly aligned sequences
>
> - **Screen.seqs** :wrench: with the following parameters
>   - "fasta" to the aligned fasta file
>   - "start" to 2308
>   - "end" to 4090
>   - "group" to the group file created by the previous screen.seqs step
>
> **Note:** we supply the count table so that it can be updated for the sequences we're removing.
>
> > ### :question: Question
> >
> >  How many sequences were removed in this step?
> > <details>
> >   <summary> Click to view answer</summary>
> >   128 sequences were removed. This is the number of lines in the bad.accnos output.
> > </details>
> {: .question}
{: .hands_on}


Now we know our sequences overlap the same alignment coordinates, we want to make sure they *only* overlap
that region. So we'll filter the sequences to remove the overhangs at both ends. In addition, there are many
columns in the alignment that only contain gap characters (i.e. "."). These can be pulled out without
losing any information. We'll do all this with filter.seqs:

> ### :pencil2: Hands-on: Filter sequences
>
> - **Filter.seqs** :wrench: with the following parameters
>   - "fasta"" to good.fasta output from Sreen.seqs
>   - "vertical" to Yes
>   - "trump" to `.`
{: .hands_on}



## Cluster sequences

The next thing we want to do to further de-noise our sequences, is to pre-cluster the sequences using the
`pre.cluster` command, allowing for up to 2 differences between sequences. This command will split the
sequences by group and then sort them by abundance and go from most abundant to least and identify
sequences that differ no more than 2 nucleotides from on another. If this is the case, then they get
merged. We generally recommend allowing 1 difference for every 100 basepairs of sequence:

> ### :pencil2: Hands-on: Perform preliminary clustering of sequences
>
> - **Pre.cluster** :wrench: with the following parameters
>   - "fasta" to the fasta output from the last Unique.seqs run
>   - "name file or count table" to the count table from the last Unique.seqs
>   - "diffs" to 2
>
> > ### :question: Question
> >
> >  How many unique sequences are we left with after this clustering of highly similar sequences?
> > <details>
> >   <summary> Click to view answer</summary>
> >   5672. <br>
> >   This is the number of lines in the fasta output
> > </details>
> {: .question}
{: .hands_on}


<!-- optional additional QC: chimera.uchime -->


> ### :pencil2: Hands-on: Remove undesired sequences
>
> - **Classify.seqs** :wrench: with the following parameters
>   - "fasta" to the fasta output from Pre.cluster
>   - "reference" to `trainset16_022016.pds.fasta` from your history
>   - "taxonomy" to `trainset16_022016.pds.tax` from your history
>   - "cutoff" to 80
>
> Have a look at the taxonomy output. You will see that every read now has a classification.
>
{: .hands_on}

> ### :pencil2: Hands-on: Cluster our data into OTUs
>
> - **Cluster.split** :wrench: with the following parameters
>   - "Split by" to `Classification using fasta`
>   - "fasta" to the fasta output from Remove.groups
>   - "taxonomy" to the taxonomy output from Remove.groups
>   - "taxlevel" to `4`
>   - "count" to the count table output from Remove.groups
>   - "cutoff" to `0.15`
>
> Next we want to know how many sequences are in each OTU from each group and we can do this using the
> `Make.shared` command. Here we tell Mothur that we're really only interested in the 0.03 cutoff level:
>
> - **Make.shared** :wrench: with the following parameters
>   - "Select input type" to `OTU list`
>   - "list" to list output from Cluster.split
>   - "count" to the count table from Remove.groups
>   - "label" to `0.03`
>
> We probably also want to know the taxonomy for each of our OTUs. We can get the consensus taxonomy for each
> OTU using the `Classify.otu` command:
>
> - **Classify.otu** :wrench: with the following parameters
>   - "list" to output from Cluster.split
>   - "count" to the count table from Remove.groups
>   - "taxonomy" to the taxonomy output from Remove.groups
>   - "label" to `0.03`
>
{: .hands_on}



## Extraction of taxonomic information

The main questions when analyzing amplicon data are: Which micro-organisms are present in an environmental samples? And in which proportion? What is the structure of the community of the micro-organisms?

The idea is to take the sequences and assign them to a taxon. To do that, we group (or cluster) sequences based on their similarity to define Operational Taxonomic Units, groups of similar sequences that can be treated as a single "genus" or "species" (depending on the clustering threshold)

![](../../images/otu.png)

> ### :pencil2: Hands-on: Extraction of OTUs with Mothur
>
> 1. Step1
> 2. Step2
>
>    > ### :nut_and_bolt: Comments
>    > A comment
>    {: .comment}
>
>    > ### :bulb: Tip: A tip
>    >
>    > * Step1
>    > * Step2
>    {: .tip}
{: .hands_on}

Once the sequences are clustered into OTUs, one sequence of each OTU is selected as a representative sequence for the OTU. The taxonomic assignation (genus, species, ...) for this sequence is searched and then assigned to all sequences of the OTU.

> ### :pencil2: Hands-on: Taxonomic assignation of the OTUs
>
> 1. Step1
> 2. Step2
>
>    > ### :nut_and_bolt: Comments
>    > A comment
>    {: .comment}
>
>    > ### :bulb: Tip: A tip
>    >
>    > * Step1
>    > * Step2
>    {: .tip}
{: .hands_on}

With the taxonomic assignation for each OTU, we can now extract for each genus (or other taxonomic level) how many OTUs (with how many sequences) are assigned to this genus (or other taxonomic level): extracting the community structure (taxon and their abundance) for the sample.

To explore the community structure, we can visualize it with dedicated tools such as Phinch:

> ### :pencil2: Hands-on: Visualization of the community structure with Phinch
>
> 1. Step1
> 2. Step2
>
>    > ### :nut_and_bolt: Comments
>    > A comment
>    {: .comment}
>
>    > ### :bulb: Tip: A tip
>    >
>    > * Step1
>    > * Step2
>    {: .tip}
{: .hands_on}

Once we have information about the community structure (OTUs with taxonomic structure), we can do more analysis on it: estimation of the diversity of micro-organism, comparison fo diversity between samples, analysis of populations, ... We will not detail such analyses here but you follow our tutorials on amplicon data analyses to learn about them.

# Whole-genome sequencing data

In the previous section, we see how to analyze amplicon data to extract the community structure. Such information can also be extracted from whole-genome sequencing (WGS) metagenomic data.

In WGS data, full genomes of the micro-organisms in the environment are sequenced (not only the 16S or 18S). We can then have access to the rRNA (only a small part of the genomes), but also to the genes of the micro-organisms. Using this information, we can try to answer to questions "What are the micro-organisms doing?" in addition to the question "What micro-organisms are present?".

## Extraction of taxonomic information

As for amplicon data, we can extract taxonomic and community structure information from WGS data. Different approaches can be used:

- Same approaches as for amplicon data with identification and classification of OTUs

    Such approaches imply a first step of sequence sorting to extract only the 16S and 18S sequences on which the same tools as for amplicon data. However, rRNA sequences represent a low proportion (< 1%) of the WGS sequences so such approache is not the most statistically supported

- Assignation of taxonomy on the whole sequences using databases with marker genes

In this tutorial, we use the second approach with MetaPhlAn2

> ### :pencil2: Hands-on: Taxonomic assignation with MetaPhlAn2
>
> 1. **MetaPhlAN2** :wrench:: Run **MetaPhlAN2** on the dereplicated sequences
>
>    > ### :question: Questions
>    >
>    > 1. What does the main output file contain?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li></li>
>    >    <li></li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
{: .hands_on}

Even if the output of MetaPhlAn2 is bit easier to parse than the BIOM file, we want to visualize and explore the community structure. We use an interactive tool called KRONA

> ### :pencil2: Hands-on: Interactive visualization with KRONA
>
> 1. **Format MetaPhlAn2 output for Krona** :wrench:: Run **Format MetaPhlAn2 output for Krona** to format MetaPhlAn2 output for KRONA
> 2. **KRONA** :wrench:: Run **KRONA** on the formatted MetaPhlAn2 output
>
>    > ### :question: Questions
>    >
>    > 1. What are the main species found for the bacteria?
>    > 2. Is the proportion of ... similar to the one found with EBI Metagenomics?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li></li>
>    >    <li></li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
{: .hands_on}

*One sentence to conclude the taxonomic analysis*

## Extraction of functional information

We would like now to answer the question "What are the micro-organisms doing?" or "Which functions are done by the micro-organisms in the environment?".

In the WGS data, we have access to the gene sequences. We use that to identify the genes, associate them to a function, build pathways, etc to investigate the functional part of the community.

> ### :pencil2: Hands-on: Metabolism function identification
>
> 1. **HUMAnN2** :wrench:: Run **HUMAnN2** on non rRNA sequences (SortMeRNA output) to extract the gene families and pathways in the sample
>
>    > ### :question: Questions
>    >
>    > 1. Which gene families is the most found one?
>    > 2. Which pathway is the most found one?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li></li>
>    >    <li></li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 2. Inspection of HUMAnN2 results?? Extraction of interesting info?
{: .hands_on}

*One sentence to conclude the functional analysis*

With the previous analyses, we investigate "Which micro-organims are present in my sample?" and "What function are done by the micro-organisms in my sample?". We can go further in these analyses (for example with combination of functional and taxonomic results). We did not detail that in this tutorial but you can found more analyses in our tutorials on whole-genome sequencing data analyses.

# Conclusion

*scheme to sum up the analyses*

*Conclusion about the technical key points. And then relation between the technics and the biological question to end with a global view.*
