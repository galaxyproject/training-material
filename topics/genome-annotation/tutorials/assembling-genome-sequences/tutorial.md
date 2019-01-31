---
layout: tutorial_hands_on
topic_name: genome-annotation
tutorial_name: assembling-genome-sequences
---

# Introduction

In this tutorial users will learn how to assemble phage genomes that have been sequenced at the Center for Phage Technology. This starts from retrieving the raw read data into the CPT Galaxy, trimming those reads, and proceeds through assembly with SPAdes and initial analysis of the contig.

> ### Agenda
>
> 1. Import Sequencing Data into a New History in Galaxy
> 2. Running the Quality Control Report and Trimming the Reads
> 3. Assemble the Contig Using SPAdes
> {:toc}
>
{: .agenda}

# Import Sequencing Data into a New History in Galaxy

The user performing the sequencing will upload the final BaseSpace data into a shared Data Library in Galaxy.

> ### {% icon hands_on %}
> 1. Under the "Shared Data" drop-down menu at the top of the Galaxy home page.
>
> ![](../../images/assembling-genome-sequences-screenshots/1_shared_data.png)
>
> 2. Click on "CPT Sequencing." On that page, the different sequencing runs are sorted by Year-Month.
>
> 3. Identify which index contains the R1 and R2 reads in question (the forward and reverse reads, respectively). Your data file should be a .fastq, fastqsanger, or .fastqsolexa/equivalent.
>
> 4. If there are multiple samples across multiple indexes, make sure each index gets its own, new history. This can be done by clicking the "to History" button at the top of the current index.
>
> ![](../../images/assembling-genome-sequences-screenshots/2_to_history.png)
>
> 5. Once the data has been imported and it is ready to assemble, return to the Galaxy Homepage where all the tools necessary can be accessed.
{: .hands_on}

# Running the Quality Control Report and Trimming the Reads

To learn about the Galaxy tool used for quality control analysis, read the [FastQC Manual](http://bficores.colorado.edu/biofrontiers-core-facility-workshops/workshops-in-the-series/short-read-2016-course-materials/day-4-sequencing-qc/day-4-files-2016/fastqc-manual/view) and watch [this quick video](https://www.youtube.com/watch?v=bz93ReOv87Y) that explains each analysis module.

> ### {% icon hands_on %} FastQC to Trimming
> 1. Using the Search bar in the **Tools** column on the left side of the Galaxy interface, search "FastQC." Find the *FastQC Read Quality reports* result.
>
>![](../../images/assembling-genome-sequences-screenshots/4_search_fastqc.png)
> 
> 2. Run the [FastQC tool](https://cpt.tamu.edu/galaxy/root?tool_id=toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.63) on both R1 and R2 reads separately. This tool results in two output entries in the history.
>
> ![](../../images/assembling-genome-sequences-screenshots/3_fastqc_tool.png)
>
> 3. Once the tool has run, click the eye {% icon solution %} symbol on the resulting dataset under the History column to view the **Webpage** generated (*not* the RawData output).
>
> 4. Identify Trimming Parameters by observing two modules, beginning with the Per Base Sequence Content Module. The lines should be parallel and fairly consistent throughout. Take note of the base # where the lines diverge too much, as these will be the points to trim at. Below is an example in which one might trim after base 245.
>
> ![](../../images/assembling-genome-sequences-screenshots/5_per_base.png)
>
> 5. The second module to analyze is the Per base Sequence Quality Module. Ideally, the quality score for reads that are trimmed will remain high (at least above 20). To better understand what is and isn't a good quality score, watch [this short video.](https://www.youtube.com/watch?v=bz93ReOv87Y) Below is an example of a good report where the first 20 bases and the bases after 250 might be trimmed.
>
> ![](../../images/assembling-genome-sequences-screenshots/6_per_base_quality.png)
>
> 6. Using the Galaxy [Trim sequences](https://cpt.tamu.edu/galaxy/root?tool_id=toolshed.g2.bx.psu.edu/repos/devteam/fastx_trimmer/cshl_fastx_trimmer/1.0.0) tool, set the base parameters, and execute for both R1 and R2.
>
> ![](../../images/assembling-genome-sequences-screenshots/9-search_trim_sequences.png)
>
> ![](../../images/assembling-genome-sequences-screenshots/8_trim_sequences.png)
>
>    > ### {% icon comment %} Note that...
>    > The Trim sequences tool cannot use a .fastq file as input. If the data is in the .fastq format, it can be converted into the appropriate format in Galaxy. Under the history entry containing the .fastq file, find the dataset for the .fastq file, and click the {% icon hands_on %} icon (**Edit Attributes**). When the attributes appear in the center panel, select "Datatype" at the top. Use the drop-down menu to select .fastqsanger (or equivalent), then click save.
>    > ![](../../images/assembling-genome-sequences-screenshots/7_change_file_type.png)
> {: .comment}
{: .hands_on}

# Assemble the Contig Using SPAdes

Search for SPAdes in the search bar underneath the Tools column on the left side of the Galaxy interface.

![](../../images/assembling-genome-sequences-screenshots/10_search_spades.png)

Selecting the Galaxy [spades tool](https://cpt.tamu.edu/galaxy/root?tool_id=toolshed.g2.bx.psu.edu/repos/lionelguy/spades/spades/1.0) and adjust the following parameters:

1. K-mers
> * All values must be *odd*, less than 128, listed in *ascending* order, and *smaller* than the read length.
> * Smaller values = more stringent (less contigs assembled, more accurate contigs)
> * Larger values = less stringent (more possible contigs, possible error in assembly)

2. Library type is **Paired-end/Single** reads

3. File format: can use one or both of the reads sets

**Unpaired/Single reads** assembly can be run using trimmed forward (R1) *or* backward (R2) reads.

![](../../images/assembling-genome-sequences-screenshots/11_unpaired_reads.png)

**Separate input files** assembly uses *both* forward and backwards reads and it usually gives the best output.

> ### {% icon comment %} Note that...
> This does not always yield the best output. For this reason, it is recommended that three instances of spades are run to start.
{: .comment}

![](../../images/assembling-genome-sequences-screenshots/12_separate_input_files.png)

Once all of the parameters have been set, **execute** the spades tool by clicking the Execute button at the bottom of the tool. This could take a few hours; however, multiple spades assemblies can be run at once. After the tools have finished running, there will be 5 outputs from each SPAdes run in the history column.

![](../../images/assembling-genome-sequences-screenshots/13_spades_datasets.png)

The contigs from each SPAdes iteration can be organized based on size by running the [Sort tool.](https://cpt.tamu.edu/galaxy/root?tool_id=sort1)

![](../../images/assembling-genome-sequences-screenshots/14_sort_tool.png)

![](../../images/assembling-genome-sequences-screenshots/15_sort_parameters.png)

Choosing "Spades scaffold stats" as the **Sort Dataset** option and "Column: 2" for **on column** will yield and outcome that looks something like this:

![](../../images/assembling-genome-sequences-screenshots/16_spades_scaffold.png)

> ### {% icon comment %} Note that...
> "SPAdes scaffold stats" and "SPAdes contig stats" should have very similar sorting results, though there can be a few differences between the two. If there is trouble finding a contig using one dataset, one possible troubleshooting solution is to check the other dataset.
{: .comment}

> ### {% icon tip %} What to Look For
> * Most of the time, there is a general idea of the expected size of the genome that has been sequenced, based on restriction digest or PFGE results. Look for contigs in that size range.
>
> * Look for contigs with a coverage much higher than the rest of the list. Often this is something the user is looking for.
> 
> * Look for contigs that are the same size in the assembly that came out of using RL alone, R2 along, and both R1 and R2 together>
>
> * If unsure, extract a few candidate nodes and try blasting them as outlined below.
>
> * In the end, it is unsure if the correct contig have been chosen (or that the contig is even complete) until a confirmation PCR is attempted, and close the genome with  closure PCR/Sanger sequencing.
{: . tip}

To extract the contig, run the **Fasta Extract Sequence** tool; this pulls out the FASTA file associated with a specific node.

![](../../images/assembling-genome-sequences-screenshots/17_fasta_extraction.png)

![](../../images/assembling-genome-sequences-screenshots/18_fasta_tool_parameters.png)

Choose "SPAdes contigs (fasta)" as the data file. Under "Sequence ID (or partial)," type in the node you want to extract. **Be sure to type in the underscore after the node number!** Example: NODE_21_. If the underscore is left out, it will extract *the wrong node*

Now, [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) that node! The will help identify which genome from the index this sequence is most likely associated with (ideally, do PCR confirmations to be 100% sure).

> * Use BLASTn + megablast, as it will yield the most closely related organisms. For a broader net, choose a different algorithm.
>    > * Doing the analysis on the home BLAST website will give a quick answer, but this is not a result that can be saved. This is recommended if an immediate answer is desired, but the following *must also* be done.
> * Doing the BLAST in Galaxy will yield a permanent link that can be stored for reference. Choose the output as BLAST XML (later, the [blast2html tool](https://cpt.tamu.edu/galaxy/root?tool_id=toolshed.g2.bx.psu.edu/repos/jankanis/blast2html/blast2html/0.0.14) must be used to convert to html) or html. After the result is ready, right-click on the eye {% icon solution %} icon and choose "open in a new tab." The hyperlink can be copied and subsequently pasted into a tracking sheet where confirmation and closure information will also be compiled.

Extracted sequences appear in the history as such:

![](../../images/assembling-genome-sequences-screenshots/19_extracted_sequence_dataset.png)

> 1. First, click on the pencil ![](../../images/assembling-genome-sequences-screenshots/20_pencil_icon.png) icon to bring upon the attributes of that dataset. Rename the dataset by editing the "Name" section of the attributes, and then select "Save." Example name: NODE_X:RawName ("raw" because this is the raw, unclosed contig).

> 2. Next, rename it using the ["Fasta Sequence Renamer" tool.](https://cpt.tamu.edu/galaxy/root?tool_id=edu.tamu.cpt.fasta.rename) Do NOT put "raw" in this name, because this is the identifier that Galaxy will use to identify the contig, no matter how it is edited.

![](../../images/assembling-genome-sequences-screenshots/21_fasta_sequence_renamer.png)

![](../../images/assembling-genome-sequences-screenshots/21_fasta_renamer_parameters.png)

Run the [PhageTerm tool](https://cpt.tamu.edu/galaxy/root?tool_id=PhageTerm) to generate a report that suggests the type of genome ends.

> * Choose the input files based on which dataset gave the contig for the phage genome. For the FASTQ **mandatory input** use the better set (usually R1). For the **optional input**, use the other dataset (usually R2).
> * Name the output file with the phage name.
> * Execute.
> * When complete, open the output called report.

> ### {% icon tip %} More Information
> * For the protocol on confirmation and closure, read [this document](https://docs.google.com/document/d/1uoSG9o2J2oWQZbHIhVRVGyIwPkbryujdeSBB2HQXTSQ/edit).
> * For the protocol on polishing a genome after annotation to prepare for depositing in GenBank, read [this document](https://docs.google.com/document/d/1aXE01fphROysxygPMrYcWdkl4iU05H6jKvhtOalTh_A/edit#heading=h.lwwpw4ay44cb).
{: .tip}


