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
3. OTU-based analysis
4. Phylotype-based analysis
5. Phylogeny-based analysis
6. Vizualisation with Phinch

Each of the Mothur tools in Galaxy contains a link to the mothur wiki in the help section. Here you can find more details about all the inputs, outputs and parameters for the tool.

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
      - upload all fastq files to your history
3. Create a **paired collection**. Since we have paired end data, each sample consist of two seperate fastq files, one containing the forward reads, and one containing the reverse reads. We can recognize the pairing from the file names, they will differ only by `_R1` or `_R2` in the filename. We can tell Galaxy about this pairing naming convention, so that our tools will know which files belong together.

 - Click on the **checkmark icon** at top of your history.

    ![](../../shared/images/history_menu_buttons2.png)

 - Select all fastq files, then click on **for all selected..** and select **Build list of dataset pairs** from the dropdown menu.
 - In the next dialog window you can create the list of pairs. By default Galaxy will look for pairs of files that differ only by a `_1` and `_2` part in their names. In our case however, these should be `_R1` and `_R2`. Please change these values accordingly. You should now see a list of pairs detected by galaxy, examine the file names, if it looks good, you can click on **auto-pair** to create the suggested pairs.

     ![](../images/create_collection.png)

 - You can change the name for each of your pairs here as well. These names will be used as sample names in the downstream analysis.

 - Once you are happy with your pairings, enter a name for your new collection at the bottom right of the screen.

 - Click the **Create List** button. A new dataset collection item should now appear in your history.

## Reference Data

Apart from our input data we will also need some additional reference data.

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step


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
    - Set **merge** to *fasta*
    - Set **inputs** to the `trim.contigs.fasta` collection output of the previous step

This will simply concatenate all the fasta files into a single file. In order to retain the knowledge of which reads came from which samples, we create a group file.

**3: *Make group file***
  - **Tool:** Make.group
  - **Parameters:**
    - Set **method** to *automatically from collection*
    - Set **fasta collection** to the `trim.contigs.fasta` output of the make.contigs step

The *group file* consists of two columns, the first is the sequence read name, and the second is the group (sample) it belongs to. The sample names were generated from the dataset names of the pairs in our collection.

```
seq1  sample1
seq2  sample1
seq3  sample2
..
```

:information_source: For those familiar with using mothur from the commandline, steps 2 and 3 generated files equivalent to the `stability.groups` and `stability.trim.fasta` files normally output by make.contigs when providing a stability file.

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

After the tool has finished, run summary.seqs tool on the filtered fasta (`good.fasta`) and examine the logfile output.

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
name   representatives
seq1    seq2,seq3,seq5,seq11
seq4    seq6,seq9,seq10
seq7    seq8
...
```

:question: How many sequences were unique? how many duplicates were removed?

**2: Generate count table**

To reduce file sizes further and streamline analysis, we can now summarize the data in a *count table*.

- **Tool:** Count.seqs
- **Parameters:**
  - Set **name** to the `names` output from Unique.seqs
  - Set **Use a Group file** to `yes`
  - Set **group** to the group file output by the Make.groups tools
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

This should finish fairly quickly, let's create another summary to see what is going on

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

So what does this mean? You'll see that the bulk of the sequences start at position 1968 and end at position 11550. Some sequences start at position 1250 or 1982 and end at 10693 or 13400. These deviants from the mode positions are likely due to an insertion or deletion at the terminal ends of the aliignments. Sometimes you'll see sequences that start and end at the same position indicating a very poor alignment, which is generally due to non-specific amplification. To make sure that everything overlaps the same region we'll re-run screen.seqs to get sequences that start at or before position 1968 and end at or after position 11550. We'll also set the maximum homopolymer length to 8 since there's nothing in the database with a stretch of 9 or more of the same base in a row (this really could have been done in the first execution of screen.seqs above).

- **Tool:** Screen.seqs
- **Parameters:**
  - Set **fasta** parameter to the aligned fasta
  - Set **start** parameter to 1968
  - Set **end** parameter to 11550
  - Set **maxhomop** to 8
  - Set **summary** to the most recent summary file (from after align.seqs)
  - Set **count** to our count_table


**Note:** we supply the count table so that it can be updated for the sequences we're removing and we're also using the summary file so we don't have to figure out again all the start and stop positions.


## Assessing error rates
## Preparing for analysis

# OTU-based analysis
## alpha diversity
## beta diversity
# Phylotype-based analysis
# Phylogeny-based analysis

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
