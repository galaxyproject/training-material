Metagenomics with Mothur: MiSeq SOP
==================================

:grey_question: ***Questions***

- *What is the effect of normal variation in the gut microbiome on host health?*

:dart: ***Objectives***

- *Learn how to analyse 16S rRNA sequences in Galaxy using the Mothur toolsuite*
- *Learn how to use dataset collections to process a large number of samples at once*

:heavy_check_mark: ***Requirements***

- *Basic knowledge of Galaxy (e.g. Galaxy Introduction module)*

:hourglass: ***Time estimation*** *1d/3h/6h*

# Introduction

This tutorial will demonstrate how to perform the *standard operating procedure (SOP)* for the analysis of 16S rRNA gene sequences generated using Illumina's MiSeq platform, with the [Mothur toolsuite](http://www.mothur.org/wiki) within Galaxy. This SOP was developed by the Schloss lab and described [here](http://www.mothur.org/wiki/MiSeq_SOP) on the Mothur wiki.

<!-- TODO: add citation
Kozich JJ, Westcott SL, Baxter NT, Highlander SK, Schloss PD. (2013): Development of a dual-index sequencing strategy and curation pipeline for analyzing amplicon sequence data on the MiSeq Illumina sequencing platform. Applied and Environmental Microbiology. 79(17):5112-20. -->

## Overview
In this tutorial we will perform the following steps:

1. Obtaining and preparing input data
2. Quality Control
  - Reducing sequencing and PCR errors
  - Processing improved sequences
  - Assessing error rate
  - Preparing for analysis  
3. OTU-based analysis
  - alpha diversity
  - beta diversity
  - population-level analysis
4. Vizualisations
  - Phinch
  - Krona

Each of the Mothur tools in Galaxy contains a link to the mothur wiki in the help section. Here you can find more details about all the inputs, outputs and parameters for the tool.

:information_source: Your results may deviate slightly from the ones presented in this tutorial due to differing tool or reference data versions or stochastic processes in the algorithms.

# Part 1: Obtaining and preparing data

## Understanding our input data
In this tutorial we are interested in understanding the effect of normal variation in the gut microbiome on host health. To that end fresh feces from mice were collected on a daily basis for 365 days post weaning. During the first 150 days post weaning (dpw), nothing was done to our mice except allow them to eat, get fat, and be merry. We were curious whether the rapid change in weight observed during the first 10 dpw affected the stability microbiome compared to the microbiome observed between days 140 and 150. We will address this question in this tutorial using a combination of OTU, phylotype, and phylogenetic methods.

To make this tutorial easier to execute, we are providing only part of the data - you are given the flow files for one animal at 10 time points (5 early and 5 late). In addition, to sequencing samples from mice fecal material, we resequenced a mock community composed of genomic DNA from 21 bacterial strains. We will use the 10 fecal samples to look at how to analyze microbial communities and the mock community to measure the error rate and its effect on other analyses.

:information_source: **Dataset details**  
Because of the large size of the original dataset (3.9 GB) you are given 21 of the 362 pairs of fastq files. For example, you will see two files: `F3D0_S188_L001_R1_001.fastq` and `F3D0_S188_L001_R2_001.fastq`. These two files correspond to Female 3 on Day 0 (i.e. the day of weaning). The first and all those with R1 correspond to read 1 while the second and all those with R2 correspond to the second or reverse read. These sequences are 250 bp and overlap in the V4 region of the 16S rRNA gene; this region is about 253 bp long. So looking at the datasets you will see 22 fastq files representing 10 time points from Female 3 and 1 mock community. You will also see `HMP_MOCK.v35.fasta` which contains the sequences used in the mock community that were sequenced in fasta format.

:pencil2: ***Hands on!***

Now that we know what our input data is, let's get it into our history:

1. Create a **new history** and name it "Mothur MiSeq SOP"

2. **Import** the training data **to your history**. There are two ways to do this. The easiest is using the data available from a *shared data library*, if this is not possible you can download the data yourself and upload it to your Galaxy instance.
 - From data library:
    - Navigate to the shared data library named *Galaxy training: Metagenomics with Mothur - MiSeq SOP* and import all fastq files you encounter there.
  - From your computer:
    - obtain data directly from [here](http://www.mothur.org/w/images/d/d6/MiSeqSOPData.zip)  <!-- TODO: zenodo link-->
    - unzip it
    - upload all fastq files to your history (40 files)

3. Create a **paired collection**. Since we have paired end data, each sample consist of two seperate fastq files, one containing the forward reads, and one containing the reverse reads. We can recognize the pairing from the file names, which will differ only by `_R1` or `_R2` in the filename. We can tell Galaxy about this paired naming convention, so that our tools will know which files belong together.

 - Click on the **checkmark icon** at top of your history.

    ![](../../shared/images/history_menu_buttons2.png)

 - Select all the fastq files (40 in total), then click on **for all selected..** and select **Build list of dataset pairs** from the dropdown menu.

 - In the next dialog window you can create the list of pairs. By default Galaxy will look for pairs of files that differ only by a `_1` and `_2` part in their names. In our case however, these should be `_R1` and `_R2`. Please change these values accordingly. You should now see a list of pairs suggested by Galaxy,

     ![](../images/create_collection.png)

 - Examine the pairings, if it looks good, you can click on **auto-pair** to create the suggested pairs.

     ![](../images/create_collection2.png)

    The middle segment is the name for each pair. You can change these names by clicking on them. These names will be used as sample names in the downstream analysis so make sure they are informative.

 - Once you are happy with your pairings, enter a name for your new collection at the bottom right of the screen.

 - Click the **Create List** button. A new dataset collection item will now appear in your history.

## Reference Data

Apart from our input data we will also need some additional reference data.

:pencil2: ***Hands on!***

1. Go back to the data library and locate the following files or download them from Zenodo and upload them to your history
  - `silva.v4.fasta`
  - `HMP_MOCK.v35.fasta`
  - `trainset9_032012.pds.fasta`
  - `trainset9_032012.pds.tax`

# Step 2: Quality Control

## Reducing sequencing and PCR errors

The first thing we want to do is combine our two sets of reads for each sample and then to combine the data from all of the samples. This is done using the make.contigs command, which requires the paired collection as input. This command will extract the sequence and quality score data from your fastq files, create the reverse complement of the reverse read and then join the reads into contigs.

:information_source: **Algorithm details**  
We have a very simple algorithm to do this. First, we align the pairs of sequences. Next, we look across the alignment and identify any positions where the two reads disagree. If one sequence has a base and the other has a gap, the quality score of the base must be over 25 to be considered real. If both sequences have a base at that position, then we require one of the bases to have a quality score 6 or more points better than the other. If it is less than 6 points better, then we set the consensus base to an N.

:pencil2: ***Hands on!***

**1: *Combine forward and reverse reads***  
  - **Tool:** Make.contigs
  - **Parameters:**
    - Set `Fastq pair` to the collection you created in the previous step
    - Leave all other parameters to the default settings
  - Execute

  :bulb: **Collections as input**  
  To provide a collection as input for a tool, click on the `Dataset collection` button in front of the parameter. The dropdown menu will now list collections as possible inputs.

   ![](../../shared/images/tools_collection_input.png)

The output from this tool will be 6 new collections, one for each type of output file (e.g. one collection with the `trim.contig.fasta` files for each pair, and another for all `scrap.contig.fasta` files). Before we continue with the analysis, we would like to combine this data into a single file. To this end we will merge all the trimmed fasta files and create a group files

**2: *Merge fasta outputs***
  - **Tool:** Merge.files
  - **Parameters:**
    - Set **merge** to *fasta files*
    - Set **inputs** to the `trim.contigs.fasta` collection output of the previous step

This will simply concatenate all the fasta files into a single file. In order to retain the knowledge of which reads came from which samples, we also create a group file.

:bulb: ***Collections as input***  
To provide a collection as input for a tool, click on the `Dataset collection` button in front of the parameter. The dropdown menu will now list collections as possible inputs.

 ![](../../shared/images/tools_collection_input.png)


**3: *Make group file***
  - **Tool:** Make.group
  - **Parameters:**
    - Set **Method** to *Automatically from collection*
    - Set **fasta collection** to the `trim.contigs.fasta` output of the make.contigs step

The *group file* consists of two columns, the first is the sequence read name, and the second is the group (sample) it belongs to. The sample names were generated from the dataset names of the pairs in our collection.

```
seq1  sample1
seq2  sample1
seq3  sample2
..
```

**4: *Summarize data***  
To get more information about the resulting fasta files, we can use the `Summary.seqs` tool.

  - **Tool:** Summary.seqs
  - **Parameters:**
    - Set `fasta` parameter to the merged fasta file from step 2
    - We do not need to supply a names or count file

The `summary` output files give information per read. The `logfile` outputs also contain some summary statistics:

```
             Start    End        NBases     Ambigs   Polymer  NumSeqs
Minimum:     1        248        248        0        3        1
2.5%-tile:   1        252        252        0        3        3810
25%-tile:    1        252        252        0        4        38091
Median:      1        252        252        0        4        76181
75%-tile:    1        253        253        0        5        114271
97.5%-tile:  1        253        253        6        6        148552
Maximum:     1        502        502        249      243      152360
Mean:        1        252.811    252.811    0.70063  4.44854
# of Seqs:   152360
```

This tells us that we have 152360 sequences that for the most part vary between 248 and 253 bases. Interestingly, the longest read in the dataset is 502 bp. Be suspicious of this. Recall that the reads are supposed to be 251 bp each. This read clearly didn't assemble well (or at all). Also, note that at least 2.5% of our sequences had some ambiguous base calls. We'll take care of these issues in the next step when we run `screen.seqs`.

**5: Filter reads based on quality and length**  

The following tool will remove any sequences with ambiguous bases and anything longer than 275 bp.

- **Tool:** Screen.seqs
- **Parameters:**
  - Set **fasta** parameter to the merged fasta file from step 2
  - Set **group** parameter to the group file created in step 3
  - Set **maxlength** parameter to `275`
  - Set **maxambig** parameter to `0`

After the tool has finished, run **summary.seqs** tool on the filtered fasta (`good.fasta`) and examine the logfile output.

:question: How many reads were removed in this screening step?
<!-- answer: 23,488. Answer can be found by looking at number of lines in bad.accnos output of screen.seqs or comparing total number of seqs between the two summary.seqs logfiles -->

## Processing improved sequences
We anticipate that many of our sequences are duplicates of each other. Because it's computationally wasteful to align the same thing a bazillion times, we'll unique our sequences using the unique.seqs command:

:pencil2: ***Hands on!***

**1: Remove duplicate sequences**  

- **Tool:** Unique.seqs
- **Parameters:**
  - Set **fasta** to the `good.fasta` output from Screen.seqs

This tool outputs two files, one is the fasta file containing only unique sequences, and a *names files*. The names file consists of two columns, the first contains the sequence names for each of the unique sequences, and the second column contains all other sequence names that are identical to the sequence in the first column.

```
name    representatives
seq1    seq2,seq3,seq5,seq11
seq4    seq6,seq9,seq10
seq7    seq8
...
```

:question: How many sequences were unique? how many duplicates were removed?
<!-- answer: 16,426 unique sequences, 112,446 duplicates-->

**2: Generate count table**

To reduce file sizes further and streamline analysis, we can now summarize the data in a *count table*.

- **Tool:** Count.seqs
- **Parameters:**
  - Set **name** to the `names` output from Unique.seqs
  - Set **Use a Group file** to `yes`
  - Set **group** to the group file output by the Screen.seqs tools
  - Set **groups** to all (or leave blank)

The count_table will look something like this:

```
Representative_Sequence                      total   F3D0   F3D1  F3D141  F3D142  ...
M00967_43_000000000-A3JHG_1_1101_14069_1827  4402    370    29    257     142
M00967_43_000000000-A3JHG_1_1101_18044_1900  28      1      0     1       0
M00967_43_000000000-A3JHG_1_1101_13234_1983  10522   425    281   340     205
...
```

The first column contains the representative sequence, and the subsequent columns contain the number of duplicates of this sequence observed per sample.

**3: Align sequences**

- **Tool:** Align.seqs
- **Parameters:**
  - Set **fasta** to the fasta output from Unique.seqs
  - Set **reference** to the `silva.v4.fasta` reference file
    - If your Galaxy is preconfigured with this reference data you will be able to find it in dropdown menu.
    - If not, set **Select Reference Template From** to `Your History` and select the appropriate file from your history.

This can take a minute or two. After the job completes, let's create another summary to see what is going on

- **Tool:** Summary.seqs
- **Parameters:**
  - Set **fasta** parameter to the aligned output from previous step
  - Set **count** parameter to count_table we made earlier

```
            Start    End      NBases  Ambigs   Polymer  NumSeqs
Minimum:    1250     10693    250     0        3        1
2.5%-tile:  1968     11550    252     0        3        3222
25%-tile:   1968     11550    252     0        4        32219
Median:     1968     11550    252     0        4        64437
75%-tile:   1968     11550    253     0        5        96655
97.5%-tile: 1968     11550    253     0        6        125651
Maximum:    1982     13400    270     0        12       128872
Mean:       1967.99  11550    252.462 0        4.36693
# of unique seqs:   16426
total # of seqs:    128872
```

So what does this mean? You'll see that the bulk of the sequences start at position 1968 and end at position 11550. Some sequences start at position 1250 or 1982 and end at 10693 or 13400. These deviants from the mode positions are likely due to an insertion or deletion at the terminal ends of the alignments. Sometimes you'll see sequences that start and end at the same position indicating a very poor alignment, which is generally due to non-specific amplification. To make sure that everything overlaps the same region we'll re-run screen.seqs to get sequences that start at or before position 1968 and end at or after position 11550. We'll also set the maximum homopolymer length to 8 since there's nothing in the database with a stretch of 9 or more of the same base in a row (this really could have been done in the first execution of screen.seqs above).

- **Tool:** Screen.seqs
- **Parameters:**
  - Set **fasta** parameter to the aligned fasta file
  - Set **start** parameter to 1968
  - Set **end** parameter to 11550
  - Set **maxhomop** to 8
  - Set **count** to our count_table

**Note:** we supply the count table so that it can be updated for the sequences we're removing.

:question: How many sequences were removed in this step?
<!--answer: 128 -->

Now we know our sequences overlap the same alignment coordinates, we want to make sure they only overlap that region. So we'll filter the sequences to remove the overhangs at both ends. Since we've done paired-end sequencing, this shouldn't be much of an issue, but whatever. In addition, there are many columns in the alignment that only contain gap characters (i.e. "."). These can be pulled out without losing any information. We'll do all this with filter.seqs:

- **Tool:** Filter.seqs
- **Parameters:**
  - Set **fasta** to good.fasta output from the last Sreen.seqs run
  - Set **vertical** to Yes
  - Set **trump** to `.`

In the log file we see the following information:

```
Length of filtered alignment: 376
Number of columns removed: 13049
Length of the original alignment: 13425
Number of sequences used to construct filter: 16298
```

This means that our initial alignment was 13425 columns wide and that we were able to remove 13049 terminal gap characters using `trump=.` and vertical gap characters using `vertical=yes`. The final alignment length is 376 columns. Because we've perhaps created some redundancy across our sequences by trimming the ends, we can re-run unique.seqs:

- **Tool:** Unique.seqs
- **Parameters:**
  - Set **fasta** to the `filtered fasta` output from Filter.seqs
  - Set **name file or count table** to the count table from the last Screen.seqs


This identified 3 duplicate sequences that we've now merged with previous unique sequences. The next thing we want to do to further de-noise our sequences is to pre-cluster the sequences using the pre.cluster command allowing for up to 2 differences between sequences. This command will split the sequences by group and then sort them by abundance and go from most abundant to least and identify sequences that are within 2 nt of each other. If they are then they get merged. We generally favor allowing 1 difference for every 100 bp of sequence:

- **Tool:** Pre.cluster
- **Parameters:**
  - Set **fasta** to the fasta output from the last Unique.seqs run
  - Set **name file or count table** to the count table from the last Unique.seqs
  - Set **diffs** to 2

We now have 5672 unique sequences. At this point we have removed as much sequencing error as we can and it is time to turn our attention to removing chimeras. We'll do this using the UCHIME algorithm that is called within mothur using the chimera.uchime command. Again, this command will split the data by sample and check for chimeras. Our preferred way of doing this is to use the abundant sequences as our reference. In addition, if a sequence is flagged as chimeric in one sample, the the default (dereplicate=No) is to remove it from all samples. Our experience suggests that this is a bit aggressive since we've seen rare sequences get flagged as chimeric when they're the most abundant sequence in another sample. This is how we do it:

- **Tool:** Chimera.uchime
- **Parameters:**
  - Set **fasta** to the fasta output from Pre.cluster
  - Set **Select Reference Template from** to `Self`
    - Set **count** to the count table from the last Pre.cluster
  - Set **dereplicate** to Yes

Running chimera.uchime with the count file will remove the chimeric sequences from the count file. But you still need to remove those sequences from the fasta file. We do this using remove.seqs:

- **Tool:** Remove.seqs
- **Parameters:**
  - Set **accnos** to the uchime.accnos file from Chimera.uchime
  - Set **fasta** to the fasta output from Pre.cluster


If we run summary.seqs again with latest fasta and count table we see:


```
            Start   End    NBases    Ambigs    Polymer    NumSeqs
Minimum:    1       376    249        0        3          1
2.5%-tile:  1       376    252        0        3          2984
25%-tile:   1       376    252        0        4          29833
Median:     1       376    252        0        4          59666
75%-tile:   1       376    253        0        5          89498
97.5%-tile: 1       376    253        0        6          116347
Maximum:    1       376    256        0        8          119330
Mean:       1       376    252.467    0        4.37281
# of unique seqs:    2628
total # of seqs:    119330
```

Note that we went from 128,655 to 119,330 sequences for a reduction of 7.3%; this is a reasonable number of sequences to be flagged as chimeric. As a final quality control step, we need to see if there are any "undesirables" in our dataset. Sometimes when we pick a primer set they will amplify other stuff that gets to this point in the pipeline such as 18S rRNA gene fragments or 16S rRNA from Archaea, chloroplasts, and mitochondria. There's also just the random stuff that we want to get rid of. Now you may say, "But wait I want that stuff". Fine. But, the primers we use, are only supposed to amplify members of the Bacteria and if they're hitting Eukaryota or Archaea, then its a mistake. Also, realize that chloroplasts and mitochondria have no functional role in a microbial community. But I digress. Let's go ahead and classify those sequences using the Bayesian classifier with the **classify.seqs** command:

- **Tool:** Classify.seqs
- **Parameters:**
  - Set **fasta** to the fasta output from Remove.seqs
  - Set **reference** to `trainset9032012.pds.fasta` from your history
  - Set **taxonomy** to `trainset9032012.pds.tax` from your history
  - Set **count** to the count table file from Chimera.uchime
  - Set **cutoff** to 80

Have a look at the taxonomy output. You will see that every read now has a classification.

Now that everything is classified we want to remove our undesirables. We do this with the remove.lineage command:

- **Tool:** Remove.lineage
- **Parameters:**
  - Set **taxonomy** to the taxonomy output from Classify.seqs
  - Set **taxon** to `Chloroplast-Mitochondria-unknown-Archaea-Eukaryota` in the textbox under *Manually select taxons for filtering*
  - Set **fasta** to the fasta output from Remove.seqs
  - Set **count** to count table from Chimera.uchime

Also of note is that *unknown* only pops up as a classification if the classifier cannot classify your sequence to one of the domains. If you run summary.seqs you'll see that we now have 2608 unique sequences and a total of 119168 total sequences. This means about 200 of our sequences were in these various groups. At this point we have curated our data as far as possible and we're ready to see what our error rate is.


## Assessing error rates

Measuring the error rate of your sequences is something you can only do if you have co-sequenced a mock community. This is something we include for every 95 samples we sequence. You should too because it will help you gauge your error rates and allow you to see how well your curation is going and whether something is wrong with your sequencing set up. First we want to pull the sequences out that were from our "Mock" sample using the get.groups command:

- **Tool:** Get.groups
- **Parameters:**
  - Set **group file or count table** to the count table from Remove.lineage
  - Set **groups** to `Mock`
  - Set **fasta** to fasta output from Remove.lineage

In the log file we see the following

```
Selected 67 sequences from your fasta file.
Selected 4060 sequences from your count file
```
This tells us that we had 67 unique sequences and a total of 4,060 total sequences in our Mock sample. We can now use the **seq.error** command to measure the error rates:

- **Tool:** Get.groups
- **Parameters:**
  - Set **fasta** to the fasta from Get.groups
  - Set **reference** to `HMP_MOCK.v35.fasta` file from your history
  - Set **count** to the count table from Get.groups

In the log file we see somthing like this:

```
It took 0 to read 32 sequences.
Overall error rate:    6.5108e-05
Errors    Sequences
0    3998
1    3
2    0
3    2
4    1
```

That rocks, eh? Our error rate is 0.0065%. We can now cluster the sequences into OTUs to see how many spurious OTUs we have:

- **Tool:** Dist.seqs
- **Parameters:**
  - Set **fasta** to the fasta from Get.groups
  - Set **cutoff** to `0.20`

- **Tool:** Cluster
- **Parameters:**
  - Set **column** to the dist output from Dist.seqs
  - Set **count** to the count table from Get.groups

- **Tool:** Make.shared
- **Parameters:**
    - Set **list** to the OTU list from Cluster
    - Set **count** to the count table from Get.groups
    - Set **label** to `0.03`

- **Tool:** Rarefaction.single
- **Parameters:**
  - Set **shared** to the shared file from Make.shared

Open the rarefaction output (in the output collection). You'll see that for 4060 sequences, we'd have 34 OTUs from the Mock community. This number of course includes some stealthy chimeras that escaped our detection methods. If we used 3000 sequences, we would have about 31 OTUs. In a perfect world with no chimeras and no sequencing errors, we'd have 20 OTUs. This is not a perfect world. But this is pretty darn good!

Now that we have assessed our error rates we are ready for some real analysis

## Preparing for analysis
We're almost to the point where you can have some fun with your data (I'm already having fun, aren't you?). We'd like to do two things- assign sequences to OTUs and phylotypes. First, we want to remove the Mock sample from our dataset using the remove.groups command:

- **Tool:** Remove.groups
- **Parameters:**
  - Set **Select input type** to `fasta , name, taxonomy, or list with a group file or count table`
  - Set **count table**, **fasta**, and **taxonomy** to the respective outputs from Remove.lineage
  - Set **groups** to `Mock`

Now we have a couple of options for clustering sequences into OTUs. For a small dataset like this, we could do the traditional approach using dist.seqs and cluster as we did with the Mock sample.

The alternative is to use the cluster.split command. In this approach, we use the taxonomic information to split the sequences into bins and then cluster within each bin. The Schloss lab have published results showing that if you split at the level of Order or Family, and cluster to a 0.03 cutoff, you'll get just as good of clustering as you would with the "traditional" approach. The advantage of the cluster.split approach is that it should be faster, use less memory, and can be run on multiple processors. In an ideal world we would prefer the traditional route because "Trad is rad", but we also think that kind of humor is funny.... In this command we use `taxlevel=4`, which corresponds to the level of *Order*. This is the approach that we generally use in the Schloss lab.

- **Tool:** Cluster.split
- **Parameters:**  
  - Set **Split by** to `classification`
  - Set **fasta** to the fasta output from Remove.groups
  - Set **taxonomy** to the taxonomy output from Remove.groups
  - Set **taxlevel** to `4`
  - Set **count** to the count table output from Remove.groups
  - Set **cutoff** to `0.15`

Next we want to know how many sequences are in each OTU from each group and we can do this using the `Make.shared` command. Here we tell mothur that we're really only interested in the 0.03 cutoff level:

- **Tool:** Make.shared
- **Parameters:**  
  - Set **Select input type** to `OTU list`
  - Set **list** to list output from Cluster.split
  - Set **count** to the count table from Remove.groups
  - Set **label** to `0.03`

We probably also want to know the taxonomy for each of our OTUs. We can get the consensus taxonomy for each OTU using the `Classify.otu` command:

- **Tool:** Classify.otu
- **Parameters:**  
  - Set **list** to output from Cluster.split
  - Set **count** to the count table from Remove.groups
  - Set **taxonomy** to the taxonomy output from Remove.groups
  - Set **label** to `0.03`

Opening the taxonomy output for level 0.03 shows something like:

```
OTU       Size    Taxonomy
Otu001    17      Bacteria(100);"Verrucomicrobia"(100);Verrucomicrobiae(100);Verrucomicrobiales(100);Verrucomicrobiaceae(100);Akkermansia(100);
Otu002    1       Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);Aeromonadales(100);Aeromonadaceae(100);Aeromonas(100);
Otu003    6       Bacteria(100);"Proteobacteria"(100);Betaproteobacteria(100);Neisseriales(100);Neisseriaceae(100);Neisseria(100);
Otu004    1       Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);unclassified(100);unclassified(100);unclassified(100);
Otu005    1       Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);Xanthomonadales(100);Xanthomonadaceae(100);Stenotrophomonas(100);
Otu006    598     Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Ruminococcaceae(100);unclassified(100);
Otu007    513     Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Lachnospiraceae(100);unclassified(100);
Otu008    1442    Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Lachnospiraceae(100);unclassified(100);
```

This is telling you that Otu001 was observed 17 times in your samples and that all of the sequences (100%) were classified as being members of the Akkermansia.

In this tutorial we will continue with this otu-based approach, for the phylotype and phylogenic approaches, please see the extra credit section at the end of this tutorial.

# Step 3: OTU-based Analysis

Let's do something more interesting and actually analyze our data. We'll focus on the OTU-based dataset. The analysis using the phylotype-based analysis is essentially the same. Also, remember that our initial question had to do with the stability and change in community structure in these samples when comparing early and late samples. Keep in mind that the group names have either a F or M (sex of animal) followed by a number (number of animal) followed by a D and a three digit number (number of days post weaning).

:pencil2: ***Hands on!***

What we now want to do is see how many sequences we have in each sample. We'll do this with the `Count.groups` command:

- **Tool:** Count.groups
- **Parameters:**
  - Set **shared** to the shared file from Make.shared

Take a look at the output. We see that our smallest sample had 2440 sequences in it. That is a reasonable number. Despite what some say, subsampling and rarefying your data is an important thing to do. We'll generate a subsampled file for our analyses with the `Sub.sample` command:

- **Tool:** Sub.sample
- **Parameters:**
  - Set **Select type of data to subsample** to `OTU Shared`
  - Set **shared** to output from Make.shared from the OTU section above
  - Set **size** to `2440`


## Alpha diversity

Let's start our analysis by analyzing the alpha diversity of the samples. First we will generate rarefaction curves describing the number of OTUs observed as a function of sampling effort. We'll do this with the `Rarefaction.single` command:

- **Tool:** Rarefaction.single
- **Parameters:**
  - Set **shared** to shared file from Make.shared

<!-- TODO: plot rarefaction curves -->
Examine the rarefaction curve output.

```
numsampled    0.03-F3D0    lci-F3D0    hci-F3D0    0.03-F3D1   ...
1              1.0000       1.0000      1.0000      1.0000
100           41.6560      35.0000     48.0000     45.0560
200           59.0330      51.0000     67.0000     61.5740
300           70.5640      62.0000     79.0000     71.4700
400           78.8320      71.0000     87.0000     78.4730
500           85.3650      77.0000     94.0000     83.9990
...
```

This file displays the number of OTUs identified per amount of sequences used (numsampled). What we would like to see is the number of additional OTUs identified when adding more sequences reaching a plateau. Then we know we have covered our full diversity. This information would be easier to interpret in the form of a graph. Let's plot the rarefaction curve for a couple of our sequences:

<!-- the following is because plotting tool will not detect columns in file in collections yet -->
First let's make our life a little bit easier. As we only have one dataset in our collection anyways, we can collapse it:
- **Tool:** Collapse Collection
- **Parameters:**
  - Set **Collection of files to collapse to a single dataset** to the rarefaction curve collection

Now we are ready to plot our rarefaction curves

- **Tool:** Plotting tool
- **Parameters:**
  - Set **Plot Title** to `Rarefaction`
  - Set **Label for x axis** to `Number of Sequences`
  - Set **Label for y axis** to `Number of OTUs`
  - Click on Insert Series,
      - Set **Dataset** to the collapsed rarefaction curve collection
      - Set **Column for x axis** to `Column 1`
      - Set **Column for y-axis** to `Column 2`
   - Add as many of the samples as you want, ignore the columns starting with lci or hci and use only those like `0.03-<samplename>`



Alas, rarefaction is not a measure of richness, but a measure of diversity. If you consider two communities with the same richness, but different evenness then after sampling a large number of individuals their rarefaction curves will asymptote to the same value. Since they have different evennesses the shapes of the curves will differ. Therefore, selecting a number of individuals to cutoff the rarefaction curve isn't allowing a researcher to compare samples based on richness, but their diversity. Finally, let's get a table containing the number of sequences, the sample coverage, the number of observed OTUs, and the Inverse Simpson diversity estimate using the `Summary.single` command. To standardize everything, let's randomly select 2441 sequences from each sample 1000 times and calculate the average:

- **Tool:** Summary.single
- **Parameters:**
  - Set **share** to shared file from Make.shared
  - Set **calc** to nseqs,coverage,sobs,invsimpson
  - Set **size** to 2440

These data will be outputted to a table called the summary file. Interestingly, the sample coverages were all above 97%, indicating that we did a pretty good job of sampling the communities. Plotting the richness or diversity of the samples would show that there was little difference between the different animals or between the early and late time points. You could follow this up with a repeated-measures ANOVA and find that there was no significant difference based on sex or early vs. late.

<!-- TODO: plot this table too? -->

## Beta diversity

Let's calculate the similarity of the membership and structure found in the various samples. We'll do this with the `Dist.shared` command that will allow us to rarefy our data to a common number of sequences.

- **Tool:** Dist.shared
- **Parameters:**
  - Set **shared** to the shared file from Make.shared
  - Set **calc** to thetayc,jclass
  - Set **size** to 2440

## Population-level Analysis

# Step 4: Visualisations


:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

# Conclusion

Conclusion about the technical key points. And then relation between the technics and the biological question to end with a global view.

:grey_exclamation: ***Key Points***

- *Simple sentence to sum up the first key point of the tutorial (Take home message)*
- *Second key point*
- *Third key point*
- *...*

# :clap: Thank you




# Extra Credit (Optional)

In this tutorial we focussed on the OTU based approach, but the original SOP also mentions phylotype- and phylogeny-based methods. These are described below for the interested reader.

## Phylotype-based analysis

For some analyses you may desire to bin your sequences in to phylotypes according to their taxonomic classification. We can do this using the phylotype command:

- **Tool:** Phylotype
- **Parameters:**  
  - Set **taxonomy** to output from Remove.groups

The cutoff numbering is a bit different for phylotype compared to Cluster and Cluster.split. Here you see 1 through 6 listed as options; these correspond to Genus through Kingdom levels, respectively. So if you want the genus-level shared file we'll do the following:

- **Tool:** Make.shared
- **Parameters:**  
  - Set **Select input type** to `OTU list`
  - Set **list** to list output from Phylotype
  - Set **count** to the count table from Remove.groups
  - Set **label** to `1`

We also want to know who these OTUs are and can run classify.otu on our phylotypes:

- **Tool:** Classify.otu
- **Parameters:**  
  - Set **list** to output from Phylotype
  - Set **count** to output from Remove.groups
  - Set **taxonomy** to output from Remove.groups
  - Set **label** to `1`

the subsequent analysis is virtually the same as the OTU based apporach.

## Phylogeny-based analysis

If you are interested in using methods that depend on a phylogenetic tree such as calculating phylogenetic diversity or the unifrac commands, you'll need to generate a tree.

- **Tool:** Dist.seqs
- **Parameters:**
  - Set **fasta** to fasta output from Remove.groups
  - Set **output** to `Phylip formatted Lower Triangle Matrix`

- **Tool:** Clearcut
- **Parameters:**
  - Set **Distance Matrix** to `Phylip Distance Matrix`
  - Set **phylip** to the distance matrix from Dist.seqs
