---
layout: tutorial_hands_on

title: Identifying Mycorrhizal Fungi from ITS2 sequencing using LotuS2

zenodo_link: https://zenodo.org/records/13710591

abbreviations:
  SPUN: Society for the Protection of Underground Networks
  MF: Mycorrhizal fungi
  EcMF: ectomycorrhizal fungi
  AMF: arbuscular mycorrhizal fungi
  eDNA: environmental DNA
  ITS2: Internal Transcribed Spacer 2
  SSU: Small Subunit (RNA)
  OTU: Operating Taxonomic Unit
  ASV: Amplicon Sequence Variant
  sdm: simple demultiplexer
  TSV: tab-separated-values


questions:
- What is the fungal community composition in a given soil sample?

objectives:
- Understand the files needed for running LotuS2
- Learn to use datasets from Zenodo in a Galaxy tool
- Learn to upload data files to Galaxy
- Learn to run the Galaxy LotuS2 tool and what the parameters mean
- Understand the output files from LotuS2
- Understand the structure of the mapping file needed by LotuS2 to link sample metadata to a pair of fastq files

time_estimation: 3H

key_points:
- LotuS is a metagenomics tool for identifying species and {OTU}s and {ASV}s from sequencing data
- Galaxy is an easy way to run LotuS2 in the cloud for bioinformatics beginners

contributors:
- sujaikumar
- bethanmanley

---

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Introduction

This tutorial is for you if you are a mycologist using metabarcoding data to understand the fungal composition of soil samples. In particular, this tutorial and its downstream applications will be of most interest to mycorrhizal scientists. {SPUN} uses Illumina sequencing of environmental {eDNA} from soil samples to identify mycorrhizal fungi. At {SPUN}, we use both {ITS2} and {SSU} amplicon sequencing to identify ectomycorrhizal and arbuscular mycorrhizal fungi, respectively. This tutorial focuses on the {ITS2} sequencing used to identify {MF}.

We do this by using a specific primer pair that amplifies just the {ITS2} region of the ribosomal DNA present in a soil sample. This tutorial covers data for an {ITS2} amplicon produced using the ITS3F and ITS4R primer set and sequenced using Illumina sequencing, which generates paired-end data. In this case, the example data has been generated using a NextSeq 2000, using 2x300bp chemistry.

This tutorial concentrates on the bioinformatics part of the process, i.e. the steps we need to do **after** we get data back from the sequencer. Starting from some example fastq files, we demonstrate how to upload data, run LotuS2, and examine the output files. We end with a brief description of what we can do with these output files.

We use LotuS2 at {SPUN} because we have found that it performs best out of all the tools we tried for identifying {MF} from sequencing data {% cite ozkurt2022 %}.

## Pre-requisites for this tutorial

If you have never used Galaxy before, we highly recommend doing their [interactive tour](https://usegalaxy.org/tours/core.galaxy_ui) first (takes a few minutes). Or, if you have more time (1h 30 min recommended), you can do the [Galaxy Basics for Everyone]({% link topics/introduction/tutorials/galaxy-intro-101-everyone/tutorial.md %}) tutorial.

We recommend signing up for a user account at one of the UseGalaxy servers (see the [first part of the Galaxy Basics tutorial]({% link topics/introduction/tutorials/galaxy-intro-101-everyone/tutorial.md %})

However, if you do not have any extra time, or just want to get started with LotuS2, you can start following the steps below by following each step carefully and using the *Tips* for each *Hands-on* exercise. You can start without even creating a user account, but then you won't be able to save your histories or come back to your analysis.

> <hands-on-title>Launch galaxy without a user account</hands-on-title>
>
> 1. Open your favorite browser (Chrome/Chromium, Safari, or Firefox, but not Internet Explorer/Edge!)
> 2. Browse to [https://usegalaxy.eu](https://usegalaxy.eu)
> 3. We recommend keeping this tutorial open in a separate browser window side-by-side if you have space on your desktop.
{: .hands_on}

# Understanding the files needed for running LotuS2

At {SPUN}, we run LotuS2 to identify mycorrhizal fungi in a set of samples, using the following input files:

1. DNA sequence files for each sample in **gzipped FASTQ** format from the Illumina MiSeq sequencer
2. A mapping file in a **tab-separated-values** format which specifies which FASTQ files correspond to which samples
3. A {sdm} file in **text** format which specifies how the sequence FASTQ files should be quality filtered and demultiplexed.

For this tutorial we have already provided a few example files at {{ page.zenodo_link }}. You can click on this link to see which files are available (there is no need to download them).

These files corresond to the three types of input files above, as shown:

- DNA sequence files:
    1. C_ITS2_S160_R1_001.fastq.gz
    1. C_ITS2_S160_R2_001.fastq.gz
    1. N5_ITS2_S140_R1_001.fastq.gz
    1. N5_ITS2_S140_R2_001.fastq.gz
    1. Pcov3_ITS2_S151_R1_001.fastq.gz
    1. Pcov3_ITS2_S151_R2_001.fastq.gz
- A mapping file:
    1. Colombia_ITS2_Mapping.tsv
- An {sdm} file
    1. sdm_miSeq_ITS.txt

For this tutorial, you do not need to download these files. Galaxy allows you to fetch data from a remote location directly into a tool or workflow without first downloading the files.

In the next section we will get the data in to Galaxy, and after that we will look at the files to see what they look like.

# Get Data

This section describes three options that will allow you to download and use {SPUN}'s example data files. Select one option and follow the instructions.

> <hands-on-title>Option 1: Data upload - Import history</hands-on-title>
>
> 1. Import history from: [tutorial input history](https://usegalaxy.eu/u/sujai_spun_earth/h/identifying-mf-from-its2-sequencing-using-lotus2---tutorial-input)
>
>    {% snippet faqs/galaxy/histories_import.md %}
>
> 2. The history should be visible in a pane on the right of your window called ‘History’. Here, you can see the files and rename them.
> 3. **Rename** {% icon galaxy-pencil %} the history to your name of choice. This should be something that helps you to remember the project, such as "SPUN Colombia 24"
>
{: .hands_on}

> <hands-on-title>Option 2: Data upload - Add to history</hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the 6 sequencing read files, 1 mapping file, and 1 sdm parameters file from [Zenodo]({{ page.zenodo_link }})
>
>    ```
>    {{ page.zenodo_link }}/files/C_ITS2_S160_R1_001.fastq.gz
>    {{ page.zenodo_link }}/files/C_ITS2_S160_R2_001.fastq.gz
>    {{ page.zenodo_link }}/files/N5_ITS2_S140_R1_001.fastq.gz
>    {{ page.zenodo_link }}/files/N5_ITS2_S140_R2_001.fastq.gz
>    {{ page.zenodo_link }}/files/Pcov3_ITS2_S151_R1_001.fastq.gz
>    {{ page.zenodo_link }}/files/Pcov3_ITS2_S151_R2_001.fastq.gz
>    {{ page.zenodo_link }}/files/Colombia_ITS2_Mapping.tsv
>    {{ page.zenodo_link }}/files/sdm_miSeq_ITS.txt
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>
{: .hands_on}

> <hands-on-title>Option 3: Data upload - Upload from your own computer</hands-on-title>
>
> 1. This option will take the longest time, and should not be used if you have a slow internet connection.
> 2. However, this option will most closely mimic how you will use this tool in your real analysis, i.e. showing you how to upload data from your own computer
> 3. First download all 8 files from [Zenodo]({{ page.zenodo_link }}) to your own local computer (desktop or laptop). You can do do this by scrolling down to the "Files" list section, and clicking the "Download All" link at the top of this section. This will download a 55 {MB} zip file to your local computer
> 4. Unzip the file - the folder should have 8 files in it.
> 5. Create a new history in Galaxy
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 6. Upload these 8 files to the new history
>
>    {% snippet faqs/galaxy/datasets_upload.md %}
>
{: .hands_on}

## Examine the input files

We will briefly look at each type of file to see that it has uploaded correctly. This doesn’t need to be done for every file each time you use LotuS2 if you have a large number of files, but it is good practice to check some files to see that they have uploaded in the correct format.

> <hands-on-title>Inspect the FASTQ file</hands-on-title>
>
> * Click on the file name of one of the fastq.gz files in the history
>    - The filename should expand to show you some information such as the size of the file, and the type of the file
> * Click on the {% icon galaxy-eye %} (eye) icon next to a fastq.gz file in the history
> * You should see 4 lines for each sequence:
     1. A header beginning with `>` followed by sequence identifiers
     2. The DNA sequence of the read, made up of ATGCN
     3. A line with just `+` on it (indicates that the next line has sequence quality values)
     4. A line with sequence quality values for each nucleotide, in ASCII format
>
> You do not need to know more about FASTQ files for this tutorial, or about sequence quality in ASCII format, but if you want to learn more, you can do the [Sequence Analysis: Quality Control tutorial](https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/quality-control/tutorial.html).
>
{: .hands_on}


> <hands-on-title>Inspect the mapping {TSV} file</hands-on-title>
>
> * Click on the {% icon galaxy-eye %} (eye) icon next to *Colombia_ITS2_Mapping.tsv* in your history
> * Each row is a sample
> * The columns tell you information about each sample:
>    - #SampleID is a short name/ID for the sample. Here, "C" stands for "Control" sample, and the other two are called samples "N5" and "Pcov3"
>    - fastqFile lists the sequence file names that correspond to each sample. We have two fastq files for paired-end sequencing: first the forward file is listed, followed by a "," (comma), and then the reverse file name is listed, with no spaces between these.
>    - ForwardPrimer: Lists the forward primer sequence used during PCR to produce the sequenced amplicon (ITS3 forward primer here).
>    - ReversePrimer: Lists the forward primer sequence used during PCR to produce the sequenced amplicon (ITS4 reverse primer here).
> * The remaining columns (Latitude, Longitude, Country, Vegetation, etc) are metadata for each sample. LotuS2 does not **need** these columns, but it will copy them to the final R Phyloseq object that it creates. These columns will be needed for doing ecological analyses later. You can add metadata to the phyloseq object that is created through LotuS2 later in R, as long as the sample ID in your metadata sheet matches the sample ID in your mapping file.
>
{: .hands_on}

> <comment-title>Note</comment-title>
> - In this case, the ForwardPrimer and ReversePrimer are the same across all samples. This is typical as we usually sequence the same region in all samples using the same primer pair.
> - LotuS2 allows you to specify primer set as a parameter, and if this is done, the primer sequence does not need to be written in the mapping file. In this case, we are including primer sequences in the mapping file. When using different primer sets, if you are following these instructions for your own data, you may need to change the primer sequences to reflect the primer set used for your own sequencing.
{: .comment}

> <hands-on-title>Inspect the {sdm} options file</hands-on-title>
>
> * Click on the {% icon galaxy-eye %} (eye) icon next to *sdm_miSeq_ITS.txt* in your history
> * This file specifies some of the parameters needed for processing ITS data with Illumina MiSeq paired-end sequencing using the software 'SDM' used by LotuS2. It can be also used for our files that were generated using an Illlumina NextSeq 2000.
>
{: .hands_on}

> <question-title>Check your understanding of the inputs</question-title>
>
> 1. What is the reverse primer for sample "N5"?
> 2. What is the DNA concentration in ng per µl for sample "Pcov3"?
> 3. When LotuS2 processes ITS paired-end sequences using these inputs, what is the minimum sequence length (minSeqLength) it will consider, below which it will discard the sequence?
> 4. Which of the 6 fastq files has the most data? (i.e., is biggest in size?)
>
> > <solution-title></solution-title>
> >
> >1. In *Colombia_ITS2_Mapping.tsv*: **TCCTCCGCTTATTGATATGC**
> >2. In *Colombia_ITS2_Mapping.tsv*, under column *DNA_concentration_ng_ul*: **20.4**
> >3. In *sdm_miSeq_ITS.txt*, next to *minSeqLength*: **110**
> >4. Click on the names of each of the fastq.gz files in ths history. The biggest is N5_ITS2_S140_R2_001.fastq.gz: **19.0 MB**
> >
> {: .solution}
>
{: .question}


## Create a "list of pairs" from the fastq.gz files

- This is a very important step when you first upload the data, because the LotuS2 tool needs our fastq.gz files in pairs.
- Galaxy has a very handy feature that automatically detects pairs of files from the filenames

> <hands-on-title>Create pairs of fastq files</hands-on-title>
>
> 1. Select ONLY the 6 fastq.gz files in your history (do not select `sdm_miSeq_ITS.txt` and `Colombia_ITS2_Mapping.tsv`)
> 2. Choose the "Build List of Dataset Pairs" option for these 6 files ![Screenshot showing how to select 6 out of 8 items, and then clicking the top right dropdown to choose the "Build List of Dataset Pairs" option](images/history-build-list-of-dataset-pairs.png)
> 3. Galaxy will automatically pair the files for you on the next screen
> 4. Rename the new collection as "Colombia ITS2 fastq pairs"
> 5. By default these two options are checked: "Hide original elements" and "Remove file extensions". You can leave them checked or unchecked. If you leave them checked, the 6 individual fastq.gz files will disappear from your history and be replaced by one collection with 3 pairs of fastq files with the new name "Colombia ITS2 fastq pairs"
>
{: .hands_on}

# Run LotuS2

When we run the LotuS2 tool on our data, it runs many steps in the background:

1. demultiplexing and filtering raw fastq sequences
2. denoising, removing chimeric sequences and clustering sequences into very high quality {OTU}s/{ASV}s
3. determining taxonomic origin of each OTU using specialized and general purpose databases and statistical algorithms
4. constructing OTU, genus, family, class, order and phylum abundance tables in .txt or .biom format
5. reconstructing the OTU phylogenetic tree
6. generating phyloseq objects for downstream analysis

As LotuS2 is a very powerful, general-purpose tool used in many metabarcoding projects for bacteria, fungi, and eukaryotes, it provides many different parameters (options for running the software) specified for each special use case.

In the next subsection we show how to run LotuS2 in Galaxy and how to set the parameters needed for a fungal dataset.

## Run LotuS2 with the example fungal dataset

> <hands-on-title>Run LotuS2</hands-on-title>
>
> In the panel on the left, select "Tools", search for "LotuS2" and then run the tool using the parameters below. Leave the rest of the parameters at their default settings. 
>
> Make sure you have the right version (2.32+galaxy0). You can check the version by clicking the {% icon tool-versions %} (blocks) icon.
>
> 1. {% tool [LotuS2](.bx.psu.edu/repos/earlhaminst/lotus2/lotus2/2.32+galaxy0) %} with the following parameters:
>    - *"Single- or Paired-end data?"*: `Paired-end list`
>        - In *"List of paired reads"*: choose the paired-list you created in the previous section: `Colombia ITS2 fastq pairs` (or whatever name you gave to the collection)
>    - In *"Mapping file (optional)"*: `Colombia_ITS2_mapping.tsv`
>    - Forward (and Reverse) Primer: Leave blank
>    - *"Clustering algorithm"*: `VSEARCH`
>    - *"Taxonomy aligner for taxonomic profiling of OTUs"*: `Lambda, LCA against custom reference database`
>        - *"Taxonomy reference database"*: `Use a built-in taxonomy database`
>            - *"Using reference database"*: `ITS fungi specific (UNITE)`
>    - In *"Other Clustering Options"*:
>        - *"Minimum size of dereplicated raw reads (optional)"*: `10:1,5:2,3:3`
>    - In *"Other Taxonomy Options"*:
>        - *"Amplicon type"*: `ITS2`
>        - *"Tax group"*: `fungal 18S/23S/ITS annotation`
>    - In *"Other options"*:
>        - *"SDM option file"*: `sdm_miSeq_ITS.txt`
>
> 2. Once all the parameters above are entered, click the **"Run Tool"** button at the end.
> 3. You should see a green box on the next page saying "Started tool LotuS2 and successfully added 1 job to the queue". The box lists the 6 outputs that the tool produces.
> 4. This step can take 5-10 minutes to run or longer depending on how many other jobs are running on the usegalaxy server. While we are waiting for it to finish, you can do the next step on "Creating your own mapping.tsv" file
> 5. You will know when the tool has finished, because all the outputs in the history will turn green.
> 
> 
>    > <comment-title>Notes</comment-title>
>    >
>    > - Remember to choose _Paired-end list_ in the sequencing read data section. Galaxy will pick up the Paired-end list available in the History, which will have the name you gave it in the _Create a list of pairs_ step
>    > - In *"Forward (and Reverse) Primer"*: Leave blank, as we have already provided them in the mapping.tsv file. If you are carrying out an analysis using your own data, you may add here the primer sequences used in your analyses in stead of in the mapping file, if you prefer.
>    > - In *"Other Clustering Options"*: *"Minimum size of dereplicated raw reads (optional)"*: we put `10:1,5:2,3:3`. Each "X:Y" pair means "A unique dereplicated read must be seen at least X times in at least Y samples". So, if a sequence read is only found in 1 sample, it must be present in 10 copies. If a sequence read is found in only 2 samples, it must be found in 5 copies in each, etc. This is so that sequence errors are not taken as novel biological sequences, reducing the occurrences of false positive OTUs or ASVs
>    {: .comment}
>
{: .hands_on}

> <question-title>Examine the parameter options</question-title>
>
> 1. What other clustering algorithms can you use inside LotuS2?
> 2. What other amplicon types can LotuS2 classify?
>
> > <solution-title></solution-title>
> >
> >1. In the "LotuS2" tool's parameter page, under *"Clustering Algorithms"*, click on the drop down. You should see these tools:
> >    - DADA2
> >    - swarm
> >    - CD-HIT
> >    - VSEARCH
> >
> >    SPUN uses VSEARCH for our {ITS2} amplicon for the identification of fungi from soils.
> >2. In *"Other Taxonomy Options"*, click on the drop down under *"Amplicon type"*:
> >    - Default
> >    - LSU
> >    - SSU
> >    - ITS
> >    - ITS1
> >    - ITS2
> >
> >    Here, we were using general fungal primers that amplify the {ITS2} region. In another tutorial we will see how to use the {SSU} region for specific amplification of {AMF}
> >
> {: .solution}
>
{: .question}

# Examine the outputs

The LotuS2 Galaxy tool creates 6 output files that you should see in your history, with names like this (the exact number after each data might be different):
1. LotuS2 on data 8, data 6, and others: main log file
2. LotuS2 on data 8, data 6, and others: mapping file
3. LotuS2 on data 8, data 6, and others: Newick-formatted phylogenetic trees between sequences
4. LotuS2 on data 8, data 6, and others: FASTA-formatted extended OTU seed sequences
5. LotuS2 on data 8, data 6, and others: OTU abundance matrix
6. LotuS2 on data 8, data 6, and others: Complete LotuS2 output

> <hands-on-title>Examine outputs</hands-on-title>
> Click on the little {% icon galaxy-eye %} (eye) icon next to each output in the history to examine it. An explanation of each file is given below.
>
> 1. **main log file**: If we had run LotuS on the command line, the contents of the `output/LotuSLogS/LotuS_run.log` would be this main LotuS2 run log file.
> 
>     It has information on all the parameters actually passed to the tool, and information on the start time of each step in the program, plus the output of each step.
>   
>     You should always check this file first to see if LotuS2 completed correctly, and approximately how many reads were used/classified. If this number is much lower than what you expected, that might indicate a problem with the run and the parameters.
> 
> 2. **mapping file**: You can ignore this file as it is a repeat of the mapping tsv file that we used as an input.
> 
> 3. **Newick-formatted phylogenetic trees between sequences**: A phylogenetic tree created from the ITS2 sequences. If you want to see what the tree looks like, you can copy the contents of this file and paste it at this [online tree viewer](http://etetoolkit.org/)treeview/
> 
> 4. **FASTA-formatted extended OTU seed sequences**: OTU sequences created by the LotuS2 program after clustering near-identical reads
> 
> 5. **OTU abundance matrix**: A columnar file with tab-separated-values. Each row is an OTU. The first column has the OTU name, and the remaining columns have the OTU abundance (i.e. how many reads were seen for that OTU) in each sample
> 
> 6. **Complete LotuS2 output**: If you try to view this file you will see some unreadable binary characters on the screen. That's because this is a zip file with the complete LotuS2 output folder in one zip folder. You should download this zip file and unzip it on your local computer if you want to see everything inside this folder. ![screenshot of how to download the zip file](images/history-download-zip.png)
> 
>     One of the most useful files in this output folder is the `phyloseq.Rdata` file which you can load in R to do further ecological analysis.
{: .hands_on}

> <question-title>Test your understanding of the outputs</question-title>
>
> If for some reason, your LotuS2 run did not complete, you can look at this [example LotuS2 run](https://usegalaxy.eu/u/sujai_spun_earth/h/identifying-mf-from-its2-sequencing-using-lotus2---tutorial-example-run) that we generated using the steps above.
>
> 1. How many total reads were in the OTU abundance matrix? (hint: main log file)
> 2. What percentage of reads are assigned at the phylum level and at the genus level? (hint: main log file)
> 3. How long did the entire LotuS2 pipeline take to run? (hint: main log has time stamps at the left of each step in the format hh:mm:ss)
> 4. Which OTU was most abundant in all samples? (hint: OTU abundance matrix). How many reads were present in each sample?
> 5. What was the sequence of the least abundant OTU? (hint: FASTA-formatted extended OTU seed sequences)  
>
> > <solution-title>Answers</solution-title>
> >
> > For the first 3 questions, scroll to near the bottom of the *"main log file"*. Under "Calculating Taxonomic Abundance Tables" the answers in our [example LotuS2 run](https://usegalaxy.eu/u/sujai_spun_earth/h/identifying-mf-from-its2-sequencing-using-lotus2---tutorial-example-run) were:
> >
> > 1. Total reads in matrix: 132123
> > 2. Phylum	95%, Genus	65%
> > 3. 00:02:34 (i.e. 2 minutes 34 seconds)
> > 4. In the *"OTU abundance matrix"*, the OTUs are arranged in order of greatest-to-least abundance so OTU1 was the most abundant with 0 reads in sample **C** (Control), 25,659 reads in sample **N5**, and 0 reads in sample **Pcov3**
> > 5. In the *"OTU abundance matrix"*, the two least abundant sequences with 6 reads total each were OTU256 and OTU257. We can look up the sequences of these files in the *"FASTA-formatted extended OTU seed sequences"* at the bottom of the file:
> >
> > ```
> > >OTU256
> > GAAATGCGATACGTAATGCGAATTGCAGAACTCAGTGAATCATCAAATCTTTGAACGCAAATTGCGCTTTCGGGATAGGCCTGAAAGCACGTTTCTTTGAGTATCGATTCAACCAACGCTCGATTGCGTGCTTTTTATTTTTCAAAATAAAAGCTTGCATTCGGTTGTAATGAGTTTTCTTTCTCTTTGAAAGTGACTTGAAGAATCGACAGTGAATGAACGATTTTCAAATCGAAACATCTGTCGGCATAAGCAATGCGGTTAAACTATTGCGTTGTGAGCTGATAGGATGTAGGCGATGTAAGAAATCGTTGGATCTTGTAACTGTTCGCAAGTGACAAGAATGACAAAATTTGATTAAGATCTCAAATGAAGCGAGGATACCCGCCGAACTTAA
> > >OTU257
> > GAAATGCGATAAGTAGTGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCAAATGGCGCTCCCTGGCATCCCGGGAAGCATGCCTGTTTGAGAACCGTAGATAGGCCAAGCTCACCAGGTTCCCCGTATCATCATCCAGCCCGGACGATGTTTACCACTCTGGTGGTGATGGGTCACTGGCGCTTTTTCGCTGCCCGTGACCTGCAACCTTTTGTCCCCGGAAGTCGACGGACCCGGATTCCCATGTCAGGTTGACACGGAGAGGGCATCACCCCCCAATCCGATGGTTTTCTGGGGAAAGACGCGGCCATCCTTCCCGGGGCGCCCCTTCAACaaaacaaaaaaaaacaaaaaaaaaCTTGGTCTCAAATCAGGCAAGAGAACCCGCTGAACTTAA
> > ```
> >
> {: .solution}
>
{: .question}


# Run LotuS2 on your own data

To run the LotuS2 workflow on your own fungal data, you only need the sequencing fastq files and some information about what primers were used.

In the example dataset, we used the ITS3/ITS4 primer pair to sequence the ITS2 reqion, so our Forward primer was `GCATCGATGAAGAACGCAGC` and the Reverse primer was `TCCTCCGCTTATTGATATGC`

When you want to process your own sequencing files, you can specify the Forward and Reverse primer in the Galaxy LotuS2 Parameters directly in the options if all your samples have the same forward and reverse primers:
  - *"Forward primer used to amplify DNA region (optional) - optional"*: `GCATCGATGAAGAACGCAGC`
  - *"Reverse primer used to amplify DNA region (optional) - optional"*: `TCCTCCGCTTATTGATATGC`

We advise adding metadata to the phyloseq object after running LotuS2, so that a large mapping file with all metadata does not have to be created. See future tutorials on adding metadata from a .csv file to a phyloseq object in R.

## Create your own mapping tsv file

Mapping files for LotuS2 have to be in tab-separated-values format. You can create this file in any spreadsheet sofware, eg: Microsoft Excel or Google Sheets. The important thing to note is that when you save or export the file, you should select the TSV format (called "Tab separated values" in Google Sheets or "Tab-delimited text .txt" in Excel)

LotuS2 allows many columns in the mapping tsv file according to the [specification](https://lotus2.earlham.ac.uk//main.php?site=documentation#mappingfile). However, for simplicity, we recommend this format for sequencing data using the SPUN protocol:

- Only 1 header row that begins with `#SampleID`
- One row per sample
- The only essential columns are `#SampleID`, `fastqFile`
- If you have the same ForwardPrimer and ReversePrimer in all your samples, then you can provide them in the LotuS2 Galaxy tool and delete the ForwardPrimer and ReversePrimer columns
- If you have any additional metadata variables that you want to analyse downstream using R, then you can provide them as columns
- We recommend using "SequencingRun" and labelling with the name of the sequencing run, which is arbitrary. Something like "Scripps run 1". This is because if samples were run across multiple different sequencing runs it is important for LotuS2 to factor this into analyses. If samples were all run on the same run, write the same for every sample.

For more information, see the [LotuS2 Mapping File Format Documentation](https://lotus2.earlham.ac.uk/main.php?site=documentation#mappingfile)

In the exercise below you will create your own mapping tsv file for a new run.

> <hands-on-title>Create and use your own mapping tsv file</hands-on-title>
>
> Steps:
>
> 1. Import history from: [mapping example history](https://usegalaxy.eu/u/sujai_spun_earth/h/identifying-mf-from-its2-sequencing-using-lotus2--create-mapping)
>
>    {% snippet faqs/galaxy/histories_import.md %}
>
> 2. Create a table in your favourite spereadsheet software (eg *Google Sheets* or *Microsoft Excel*) with the following columns (you can copy paste from the example below. We have filled in the first two columns for you #SampleID and fastqFile):
>    | #SampleID | fastqFile | ForwardPrimer | ReversePrimer | Vegetation |
>    | N11 | N11_ITS2_S134_R1_001.fastq.gz,N11_ITS2_S134_R2_001.fastq.gz | | | |
>    | N16 | N16_ITS2_S138_R1_001.fastq.gz,N16_ITS2_S138_R2_001.fastq.gz | | | |
>
> 3. Complete the table:
>   - Add the correct ForwardPrimer `GCATCGATGAAGAACGCAGC` and the correct ReversePrimer `TCCTCCGCTTATTGATATGC` to each row
>   - Add the Vegetation `Trigonobalanus excelsa forest` to each row
>
> 4. Export or Save the spreadsheet as a new file in tab-separated-values format (Google Sheets: _File_: _Download_: _Tab separated values (.tsv)_; Excel: _File_: _Save as_: choose Tab delimited text in the format)
>
> 5. Name the file something like "mapping.tsv" or "mapping.txt"
>
> 6. Upload this file to your history as a local file
>
>    {% snippet faqs/galaxy/datasets_upload.md %}
>
> 7. Now you should have all the files you need in this new history to run LotuS2:
>    - 4 paired-end fastq files for 2 samples
>    - 1 mapping tsv file specifying sample IDs, fastqFiles, primers, and metadata for these 2 samples
>    - 1 sdm_miSeq_ITS.txt
>
> 8. Rerun LotuS2 following the [steps above](#run-lotus2)
>
> 9. If your LotuS2 run finishes - then CONGRATULATIONS! You now know how to create your own mapping tsv file and run LotuS2 on a new set of fastq.gz files
>
{: .hands_on}

# Conclusion

- LotuS2 is a powerful tool that can be used to identify species present in a DNA sequencing sample
- You have learnt how to run LotuS2 on example data
- You have learnt how to create a mapping tsv file for a new set of data
