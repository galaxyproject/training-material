---
layout: tutorial_hands_on
topic_name: your_topic
tutorial_name: your_tutorial_name
---

# Introduction
{:.no_toc}

<!-- This is a comment. -->

Here we will demonstrate genome analyses strategies for understanding structural differences between a newly assembled genome and a set of published, annotated genomes. 

> ### Agenda
>
> In this tutorial we begin with a new genome assembly just produced in the [Unicycler tutorial]({{site.baseurl}}/topics/assembly/tutorials/unicycler-assembly/tutorial.html). This is an assembly of *E. coli* C, which we will be comparing to assemblies of all other complete genes of *E. coli*.
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Finding and loading all complete *E. coli* genomes

*E. coli* is one of the most studied organisms. Naturally, there are hundreds of complete genomes. Here we will shows how to uploaded all (!) complete *E. coli* genomes as once. 

## Preparing the data

Our initial objective is to compare our assembly against all complete *E. coli* genomes to identify the most related ones and to identify any interesting genome alterations. In order to do this we need to align our assembly against all other genomes. And in order to do that we need to first obtain all these other genomes. 

[NCBI](https://www.ncbi.nlm.nih.gov/) is the resource that would store all complete *E. coli* genomes. Specifically, they can be found [here](https://www.ncbi.nlm.nih.gov/genome/genomes/167). As we will see this list contains over 500 genomes and so uploading them by hand will likely result in carpal tunnel syndrome, which we want to prevent. Galaxy has several features that are specifically designed for uploading and managing large sets of similar types of data. The following two **Hands-on** section show how they can be used to import all completed *E. coli* genomes into Galaxy. 

<!--
{% icon hands_on %} will render the hands_on icon as specified in
_config.yml in the root of this repository.
-->

> ### {% icon hands_on %} Hands-on: Preparing a list of all complete *E. coli* genomes
>
>Open [the NCBI list of of *E. coli* genomes](https://www.ncbi.nlm.nih.gov/genome/genomes/167) in a new window and position two browser windows (one the tutorial and the one you just opened) side by side. Then follow the steps in the following video. 
>
>---------------------
>
><div style="padding:56.25% 0 0 0;position:relative;"><iframe src="https://player.vimeo.com/video/271328293?title=0&byline=0&portrait=0" style="position:absolute;top:0;left:0;width:100%;height:100%;" frameborder="0" webkitallowfullscreen mozallowfullscreen allowfullscreen></iframe></div><script src="https://player.vimeo.com/api/player.js"></script>
>
{: .hands_on}

## Getting complete *E. coli* genomes into Galaxy

Now that the list is formatted as a table in a spreadsheet it is time to upload it into Galaxy. There is a problem though &uarr; the URLs (web addresses) in the list do not actually point to sequence files that we would need to perform alignments. Instead they point to directories. For example, this URL:

```
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/008/865/GCA_000008865.1_ASM886v1					
```

points to a directory (rather than a file) containing the following files:

![GenBank assembly files for an E. coli strain](../../images/genbank_dir.png "A list of files for an E. coli assembly. For further analyses we need datasets ending with '_genomic.fna.gz'.")

So to download sequence files we need to edit URLs by adding filenames to them. For example, in the case of URL shown above we need to add `/GCA_000008865.1_ASM886v1_genomic.fna.gz` to the end to get this:

```
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/008/865/GCA_000008865.1_ASM886v1/GCA_000008865.1_ASM886v1_genomic.fna.gz
```

this can be done as a two step process where we first copy the end part of the existing URL (`/GCA_000008865.1_ASM886v1_genomic.fna.gz`) and then add a fixed string `_genomic.fna.gz` to the end of it. Doing by hand is crazy and trying it in a spreadsheet is complicated. Fortunately Galaxy's new rule-based unloader helps with that as shown in the next **Hands-on** section.

Short introduction about this subpart.

> ### {% icon hands_on %} Hands-on: Data upload
>
>Here we copy data from the spreadsheet described in the previous section into Galaxy's rule-based uploader to download several hunder complete genomes into a Collection. Follow the steps in the video below.
>
>----------------------
>
><div style="padding:56.25% 0 0 0;position:relative;"><iframe src="https://player.vimeo.com/video/271336444?title=0&byline=0&portrait=0" style="position:absolute;top:0;left:0;width:100%;height:100%;" frameborder="0" webkitallowfullscreen mozallowfullscreen allowfullscreen></iframe></div><script src="https://player.vimeo.com/api/player.js"></script>
{: .hands_on}

## Preparing assembly

Before doing any analyses we need to upload assembly produced in [Unicycler tutorial]({{site.baseurl}}/topics/assembly/tutorials/unicycler-assembly/tutorial.html) from Zenodo. 

 > ### {% icon hands_on %} Uploading *E. coli* assembly into Galaxy
 >
 > 1. Open upload tool (Upload icon on the top of the left pane)
 > 2. Click **Paste/Fetch data** button (Bottom of the interface box)
 > 3. Paste `https://zenodo.org/record/1251125/files/Ecoli_C_assembly.fna` into the box.
 > 4. Set **Type** to `fasta`
 > 5. Click **Start**
{: .hands_on}

The assembly we just uploaded has two issues that need to be addressed before proceeding with our analysis:

 1. It contains two sequences: the one of *E. coli* C genome (the one we really needed) and another representing phage phiX174 (a by product of Illumina sequencing were it is used a spike in DNA). 
 2. Sequences have unwieldy names like `>1 length=4576293 depth=1.00x circular=true`. We need to rename it to something more meaningful.

 Let's fix these two problems.

 > ### {% icon tip %} Tip: Finding tools mentioned in this tutorial
 >Galaxy instances contain hundreds of tools. As a result it is sometimes hard to find tools mentioned in tutorials such as this one. 
 >
 >To help with this challenge Galaxy has a search box at the top of the left panel. Use this box to find the tools mentioned here.
 >![](../../images/tool_search.png "Use search box to find tools!")
 {: .tip}

> ### {% icon hands_on %} Hands-on: Fixing assembly
>
> 1. First we will use **Filter sequences by length** {% icon tool %} tool to remove phiX 174 genome. Set parameters as follows:
>   - **Fasta file** to the dataset you've just uploaded. (It will have a name something like `https://zenodo.org/record/1251125/files/Ecoli_C_assembly.fna` ).
>   - **Minimal length** to `10000` (this was phiX174, which is around 5,000 bp, will be filtered out)
>   - **Maxumum length** you do not need to change.
> 2. Second we will use **Text transformation with sed** {% icon tool %} tool to convert name of the assembly from `>1 length=4576293 depth=1.00x circular=true` to `>Ecoli_C`. Set parameters as follows:
>   - **File to process** should be set to the output of the previous step
>   - Inside **SED program** box enter the following expression (so called [Regular Expression](https://en.wikipedia.org/wiki/Regular_expression): `s/^\>1.*$/\>Ecoli_C/`
>
>	> ### {% icon tip %} Highlight: SED editor and Regular Expressions
>	>The expression `s/^>1.*$/>Ecoli_C/` contains several pieces that you need to understand. Let's write it top-to-bottom and explain:
>	>
>	> - `s` - tell SED to *Substitute*
>	> - `/` - opens a section of the commands telling SED *what* to substitute. 
>	> - `^` - tell SED to start looking at *the beginning* of each line
>	> - `>` - is the first character we want to match. Remember that name of the sequence in FASTA files starts with `>`
>	> - `1` - is the number present is our old name (`>1 length=4576293 depth=1.00x circular=true` to `>Ecoli_C`)
>	> - `.` - dot has a special meaning. It signifies *any* character. 
>	> - `*` - is a *quantifier*. From [Wikipedia](https://en.wikipedia.org/wiki/Regular_expression): "The asterisk indicates zero or more occurrences of the preceding element. For example, ab*c matches `ac`, `abc`, `abbc`, `abbbc`, and so on."
>	> - `$` - signifies *the end* of a line
> 	> - `/` - is *the end* of the *what to substitute* section. It also serves as the beginning of *what to substitute WITH* section
> 	> - `>` - is the required element of the FASTA sequence name
>	> - `Ecoli_C` is the *name* we want the sequence to have
>	> - `/` - is the end of SED command
>	>
>	>So in short we are replacing `>1 length=4576293 depth=1.00x circular=true` with `>Ecoli_C`. The *Regular expression* `^\>1.*$` is used here to represent `>1 length=4576293 depth=1.00x circular=true`.<br>
>	>Detailed description of regular expressions is outside of the scope of this tutorial, but there are other great resources. Start with [Software Carpentry Regular Expressions tutorial](http://v4.software-carpentry.org/regexp/index.html)!
>	{: .tip}
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What is the meaning of `^` character is SED expression?
>    > 2. Where do you go to learn more about regular expressions?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>It tells SED to start matching from the beginning of the string.</li>
>    >    <li>Software carpentry at <a href="https://software-carpentry.org">https://software-carpentry.org</a></li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
{: .hands_on}

# Generating alignments

Now everything is loaded and ready to go. We will now align our assembly against each of the *E. coli* genomes we have uploaded into the collection. To do this we will use [LASTZ](https://lastz.github.io/lastz/) an aligner designed to align long sequences. 

> ### {% icon hands_on %} Hands-on: Running LASTZ
> 1. Open **LASTZ** interface
> 2. Change **Select TARGET sequence(s) to align against** to `from your history`
> 3. In **Select a reference dataset** click on the folder icon (![](../../images/folder-o.png)) and the collection containing all *E. coli* genomes we uploaded below. 
> 4. In **Select QUERY sequence(s)** choose our assembly which was prepared in the previous step.
> 5. Find section of LASTZ interface called **Chaining** and expand it.
> 6. Set **Perform chaining of HSPs with no penalties** to `Yes`
> 7. Find section of LASTZ interface called **Output** and expand it.
> 8. Set **Specify the output format** to `blastn`
> 9. Run LASTZ by clicking **Execute** button
{: .hands_on}

Note that because we started LASTZ on *a collection* of *E. coli* genomes is will output alignment information as *a collection* as well. Collection is simply a way to represent large sets of similar data in compact way within Galaxy's interface.

> ### {% icon warning %} It will take a while!
> Please understand that alignment is not an instantaneous process: allow several hours for these jobs to clear.
{: .warning-box}

# Finding closely related assemblies

## Understanding LASTZ output

LASTZ produced data in so-called `blastn` format, which looks like this:

```
         1       2     3   4  5 6       7       8    9   10      11    12
-------------------------------------------------------------------------
BA000007.2 Ecoli_C 66.81 232 51 6 3668174 3668397 5936 6149 3.2e-40 162.7
BA000007.2 Ecoli_C 57.77 206 38 8  643802  643962 5945 6146 1.6e-18  90.6
BA000007.2 Ecoli_C 67.03 185 32 6 4849373 4849528 5965 6149 2.9e-28 122.9
BA000007.2 Ecoli_C 63.06 157 33 3 1874604 1874735 5991 6147 5.8e-26 115.3
```

where columns are:

 1. `qseqid` - query (e.g., gene) sequence id
 2.	`sseqid` - subject (e.g., reference genome) sequence id
 3.	`pident` - percentage of identical matches
 4.	`length` - alignment length
 5.	`mismatch` - number of mismatches
 6.	`gapopen` - number of gap openings
 7.	`qstart` - start of alignment in query
 8.	`qend` - end of alignment in query
 9.	`sstart` - start of alignment in subject
 10. `send` - end of alignment in subject
 11. `evalue` - [expect value](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ#expect)
 12. `bitscore`	- [bit score](https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html)

 The alignment information produced by LASTZ is a collection. In this collection each element contains alignment data between each of the *E. coli* genomes and our assembly:

![LASTZ collection](../../images/lastz_collection.png "LASTZ produced a collection where each element corresponds to an alignment between an <i>E. coli</i> genome and our assembly. Here one of the elements is expanded.").

## Collapsing collection

Collections are a wonderful way to organize large sets of data and parallelize data processing like we did here. However, at this point we need to combine all data into one dataset. Follow the steps below to accomplish this:

> ### {% icon hands_on %} Hands-on: Combining collection into a single dataset
> 1. Open **Collapse Collection** tool.
> 2. In **Collection of files to collapse** click on the folder icon (![](../../images/folder-o.png)) and select the output of **LASTZ** produced on previous step. 
> 3. Leave all other options as they are - no changes needed. 
> 4. Click **Execute**
{: .hands_on}

This will produce one gigantic table (over 12 million lines) containing combined LASTZ output for all genomes. 

## Getting taste of the alignment data

To make further analyses we need to obtain some understanding about the nature of the alignment data. To do this let's select a random subsample of the large dataset we've generated above. This is necessary because processing the entire dataset will take time and will not give us a better insight anyway. So first we will select 10,000 lines from the alignment data:

> ### {% icon hands_on %} Hands-on: Selecting random subset of data
> 1. Open **Select random lines from a file** tool.
> 2. Set **Randomly select** to `10000`
> 3. In **from** choose the dataset generated on previous step (it will called something like `Collapse Collection on data...` )
> 4. Click **Execute**
{: .hands_on}

Now we can visualize this dataset to discover generalities:

> ### {% icon hands_on %} Hands-on: Graphing alignment data
> 1. Expand the random subset of alignment data generated on the previous step by clicking on it.
> 2. You will see "chart" button (![](../../images/bar-chart-o.png)). Click on it.
> 3. In the center pane you will see a list of visualizations. Select **Scatter plot (NVD3)**
> 4. Click **Select data** button (![](../../images/disks.png))
> 5. Set **Values for x-axis** to `Column: 3` (alignment identity)
> 6. Set **Values for y-axis** to `Column: 4` (alignment length)
> 7. You can also click on configuration button (![](../../images/chart_cog.png)) and specify axis labels etc.
{: .hands_on}

The relationship between the alignment identity and alignment length looks like this (remember that this only a subsample of the data):

![Identity versus length](../../images/id_vs_len.png "Alignment identity (%) versus length (bp). This graph is truncated at teh top")

you can see that most alignments are short and have relatively low identity. Thus we can filter the original dataset by identity and length. Judging form this graph we can selected alignment longer than 10,000 bp with identity above 90%. 

> ### {% icon hands_on %} Hands-on: Filtering data
> 1. Select **Filter data on any column using simple expressions** tool.
> 2. In **Filter** select the full dataset. **IMPORTANT** you need to select the full dataset, not the down-sampled one, but the one generated by collection collapsing operation. 
> 3. In **With following condition** enter our filtering criteria: `c3 >= 90 and c4 >= 10000` (here `c` stands for *column*).
> 4. Click **Execute**
{: .hands_on}

## Aggregating data

Remember, our objective is to find genomes that are most similar to our. Given the alignment data in the table we just created we can define similarity as follows:

------
*Genomes that have the smallest number of alignment blocks but the highest overall alignment length are most similar to our assembly. This essentially means that they have longest uninterrupted region of high similarity to our assembly.*

-------

However, to extract this information from our data we need to aggregate it. In other words, for each *E. coli* genome we need to calculate the total number of alignment blocks, their combined length, and average identity. The following section explain how to do this:

> ### {% icon hands_on %} Hands-on: Aggregating the data
> 1. Select **Datamash (operations on tabular data)**
> 2. In **Input tabular dataset** select results of the filtering performed in the previous step (it will have a name with `Filter on data...` in it. 
> 3. In **Group by fields** enter `1`. This is because column 1 contains name of the *E. coli* genome we mapped against. 
> 4. Set **Sort input** to `Yes`.
> 5. On **Operation to perform on each group** set **Type** to `Count` and **On column** to `Column: 1`.
> 6. Click **Insert operation to perform on each group** button twice to add two more input boxes.
> 7. In second **Operation to perform on each group** set **Type** to `Mean` and **On column** to `Column: 3`.
> 8. In third **Operation to perform on each group** set **Type** to `Sum` and **On column** to `Column: 4`.
> 9. Click **Execute**
{: .hands_on}

## Finding closest relatives

Dataset generated above lists each *E. coli* genome accession only once and will have aggregate information of the number of alignment blocks, mean identity, and total length. Let's graph these data:

> ### {% icon hands_on %} Hands-on: Graphing aggregated data
> 1. Expand the aggregated data generated on the previous step by clicking on it.
> 2. You will see "chart" button (![](../../images/bar-chart-o.png)). Click on it.
> 3. In the center pane you will see a list of visualizations. Select **Scatter plot (NVD3)**
> 4. Click **Select data** button (![](../../images/disks.png))
> 5. Set **Data point labels** to `Column: 1` (Accession number of each *E. coli* genome)
> 5. Set **Values for x-axis** to `Column: 2` (# of alignment blocks)
> 6. Set **Values for y-axis** to `Column: 4` (Total alignment length)
> 7. You can also click on configuration button (![](../../images/chart_cog.png)) and specify axis labels etc.
{: .hands_on}

The relationship between the number of alignment blocks and total alignment length looks like this:

![Identity versus length](../../images/best_genomes_chart.png "Number of alignment blocks versus total alignment length (bp).")

A group of three dots in the upper left corner of this scatter plot represents genomes that are most similar to our assembly: they have the small number of alignment blocks and high total alignment length. Mousing over these three dots (if you set **Data point labels** correctly in the previous step) will reveal they accession numbers: `LT906474.1`, `CP024090.1`, and `CP020543.1`. 

> ### {% icon warning %} Things change
> It is possible that when you will be repeating these steps the set of sequences in NCBI will change and you will obtain different accession numbers. Keep this in mind.
{: .warning-box}

Let's find table entries corresponding to these:

> ### {% icon hands_on %} Hands-on: Extracting into about best hits
> 1. Open **Select lines that match an expression** tool. 
> 2. Set **Select lines from** to previous dataset (named like `Datamash on data...`)
> 3. In **the pattern** enter `LT906474|CP024090|CP020543`. Here `|` means `or`.
> 4. Click **Execute**
{: .hands_on}

This will generate a short table like this:

```
CP020543.1 11 99.91	4487098
CP024090.1 12 99.91	4540487
LT906474.1  8 99.94	4575223
```

From this it appears that `LT906474.1` is closes to our assembly as it has eight alignment blocks, longest total alignment length (4,575,223) and highest mean identity (99.94%). 

# Comparing genome architectures

Now that we know the three genomes most closely related to ours, let's take a closer look at them. First we will re-download sequence and annotation data. 

## Getting sequences and annotations

> ### {% icon hands_on %} Hands-on: Uploading sequences and annotations
> Using the three accession listed above we will fetch necessary data from NCBI. Follow the steps in the video below:
>
>----------
><div style="padding:56.25% 0 0 0;position:relative;"><iframe src="https://player.vimeo.com/video/272379016?title=0&byline=0&portrait=0" style="position:absolute;top:0;left:0;width:100%;height:100%;" frameborder="0" webkitallowfullscreen mozallowfullscreen allowfullscreen></iframe></div><script src="https://player.vimeo.com/api/player.js"></script>
>-----------
>
>At the end of this you should have two collections: one containing genomic sequences and another containing annotations. 
{: .hands_on}

## Visualizing rearrangements

Now we will perform alignments between our assembly and the three most closely related genomes to get a detailed look at any possible genome architecture changes. We will again use LASTZ:

> ### {% icon hands_on %} Hands-on: Aligning again
> 1. Open **LASTZ** interface
> 2. Change **Select TARGET sequence(s) to align against** to `from your history`
> 3. In **Select a reference dataset** click on the folder icon (![](../../images/folder-o.png)) and select the collection of the three genomes (in the video above we called it `DNA`). 
> 4. In **Select QUERY sequence(s)** choose our assembly which was prepared in the beginning (it has a name `Text transformation on data...`).
> 5. Find section of LASTZ interface called **Chaining** and expand it.
> 6. Set **Perform chaining of HSPs with no penalties** to `Yes`
> 7. Find section of LASTZ interface called **Output** and expand it.
> 8. Set **Specify the output format** to `Customized general`
> 9. Within **Select which fields to include** select the following:
>	* `score` - alignment score
>	* `name1` - name of the *target* sequence
>	* `strand` - strand for the *target* sequence
>	* `zstart` - 0-based start of alignment in *target*
>	* `end1` - end of alignment in *target*
>	* `length1` - length of alignment in *target*
>	* `name2` - name of *query* sequence
>	* `strand2` - strand for the *query* sequence
>	* `zstart2` - 0-based start of alignment in *query*
>	* `end2` - end of alignment in *query*
>	* `identity` - alignment identity 
>	* `number` - alignment number 
> 9. In **Create a dotplot representation of alignments?** select `Yes`
> 10. Run LASTZ by clicking **Execute** button
{: .hands_on}

Because we chose to produce Dot Plots as well LASTZ will generate two collections: one containing alignment data and the other containing DotPlots in PNG format:

![Dot Plots](../../images/three_dot_plots.png "Dot Plot representations of alignments between three <i>E. coli</i> genomes and our assembly. Query (Y-axis) is indicated above each dot plot. Query (X-axis) is our assembly. Red circle indicates a region deleted in our assembly.")

A quick conclusion that can be drawn here is that there is a large inversion in CP020543 and deletion in our assembly. If you are not sure how to interpret Dot Plots here is a great explanation by [Michael Schatz](http://schatz-lab.org/):

![Interpreting Dot Plots](../../images/dotplot.png "A quick reference to interpreting Dot Plots. Our case is identical to <i>Insertion into Reference</i> shown in the upper left.")

For a moment let's leave LASTZ result and create a browser that would allows us to display our results. 

## Producing a Genome Browser for this experiment

Dot plots we've produced above are great, but they are static. It would be wonderful to load these data into a genome browser where one can zoom in and out as well as add tracks such as those containing genes. To create a browser we need a genome and a set of tracks. Tracks are features such as genes or SNPs with start and end positions corresponding to a coordinate system provided by the genome. Thus the first thing to do is to create a *genome* that would represent our experiment. We can create such a genome by simply combining the three genomes of closely related strains with our assembly in a single dataset - a hybrid genome.

> ### {% icon hands_on %} Hands-on: Creating a single FASTA dataset with all genomes
> First step will be collapsing the collection containing the three genomes into a single file. To do this:
> 1. Open **Collapse Collection** tool.
> 2. In **Collection of files to collapse** click on the folder icon (![](../../images/folder-o.png)) and select the collection of the three genomes (in the video above we called it `best hits`). 
> 3. Leave all other options as they are - no changes needed. 
> 4. Click **Execute**
{: .hands_on}

This will produce a single FASTA dataset containing the three genomes. There is one problem though. If we look at the data in this file, we will see that FASTA headers look like this:
```
>CP020543.1 Escherichia coli C, complete genome
```
This is a problem because a browser will "think" that this particular genome is called `CP020543.1 Escherichia coli C, complete genome` while in the alignment files produced by LASTZ the same genome will be listed as simply `CP020543.1`. Because these two seemingly identical things are technically different it will notbe possible to render alignment results (or any other annotation) on a browser. To solve this issue we simply need to remove ` Escherichia coli C, complete genome` from `>CP020543.1 Escherichia coli C, complete genome` and convert it into `>CP020543.1`. For this we will use **sed** tool we already used above to [prepare assembly files](#preparing-assembly):

> ### {% icon hands_on %} Hands-on: Cleaning sequence names
> 1. Open **Text transformation with sed** {% icon tool %} tool
> 2. Set **File to process** to the output of the previous step (a FASTA file produced by collection collapse)
> 3. Inside **SED program** box enter the following expression (so called [Regular Expression](https://en.wikipedia.org/wiki/Regular_expression): `s/\ Esc.*$//`. Here we are matching from space (`\ `) separating `CP020543.1` and `Escherichia coli C, complete genome` and substituting this with nothing. 
> 4. Click **Execute**
{: .hands_on}

To make sure that everything completed correctly let's grab FASTA headers from all sequences in the dataset produced by the last tool:

> ### {% icon hands_on %} Hands-on: "Grepping" FASTA headers
> 1. Open **Search in textfiles (grep)** tool.
> 2. In **Select lines from** select output of the previous step (called `Text transformation on data...`)
> 3. In **Regular Expression** box enter `^>`. This tells to return all line that begin (`^` signifies the beginning) with `>`.
> 4. Click **Execute**
{: .hands_on}

If everything went well we will see something like this:
 ```
 >CP020543.1
 >CP024090.1
 >LT906474.1
 ```

Finally, we need to add our own assembly to the FASTA dataset containing the three genomes. This can be done by a simple concatenation:

> ### {% icon hands_on %} Hands-on: Concatenate FASTA files
> 1. Open **Concatenate datasets tail-to-head (cat)**
> 2. In **Datasets to concatenate** select output of **sed** tool we performed one step ago (*before* last step; it is called `Text transformation on...`)
> 3. Click **Insert Dataset** button
> 4. Select our assembly (its name also begins with `Text transformation on...` but is located earlier in the history)
{: .hands_on}

The resulting dataset contains four sequences: three genomes plus our assembly. Let's start a browser using these sequences:

> ### {% icon hands_on %} Hands-on: Starting a custom IGV browser
> 1. Go to [IGV web page](http://software.broadinstitute.org/software/igv/download) and launch a browser appropriate for your platform. Wait for it to start. It will display human genome, but we will change that. 
> 2. Go back to your Galaxy session and expand the dataset generated during the last step.
> 3. Click on `local` link in **display with IGV local**
> 4. Wait a bit and IGV will refresh displaying "chromosomes" of our *hybrid* genome:
>
>-------------
>
> ![Empty IGV](../../images/igv_empty.png "IGV instance displaying <i>Hybrid</i> genome without tracks")
{: .hands_on}

## Preparing and displaying alignments

[Above](#-hands-on-aligning-again) we computed alignments using LASTZ. Because we ran LASTZ on a collection containing genomic sequences, LASTZ produced a collection as well (actually two collections: one containing alignments an the other with dot plots). To display alignments in the browser we need to do several things:

 1. Fix unwanted `%` signs in LASTZ output
 2. Create names for alignment blocks
 3. Convert LASTZ output into [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format
 4. Create a single BED track containing alignments against all four genomes.

To begin, let's look at the LASTZ output:

```
       1          2 3      4      5      6       7 8      9     10            11    12  13
------------------------------------------------------------------------------------------
10141727 CP020543.1 +     48 106157 106109 Ecoli_C +      0 106109 106107/106109 100.0%  1
    5465 CP020543.1 + 121267 121367    100 Ecoli_C + 109317 109418    76/100     76.0%   2
    4870 CP020543.1 + 159368 159512    144 Ecoli_C + 128706 128828    95/115     82.6%   3
```

One immediate problem is `%` character in column 12 (alignment identity). We need to remove it. For this we will use **SED** tool that should be familiar to us from [previous hands-on exercises](#-hands-on-cleaning-sequence-names):

> ### {% icon hands_on %} Hands-on: Removing `%` character from LASTZ output
>
> 1. Open **Text transformation with sed** {% icon tool %} tool
> 2. In **File to process** click on the folder icon (![](../../images/folder-o.png)) and select the output of LASTZ (called `LASTZ on collection ...: mapped reads` ).
> 3. Inside **SED program** box enter the following expression (so called [Regular Expression](https://en.wikipedia.org/wiki/Regular_expression): `s/\%//`. Here we are matching percent character `%` (it is pre-pended with `\` because it is a special character, but we want `sed` to interpret it literally, as the percentage sign) and substituting this with nothing. 
> 4. Click **Execute**
{: .hands_on}

As a result LASTZ output will look like this (no `%` signs):

```
      1          2 3      4      5      6       7 8      9     10            11    12  13
------------------------------------------------------------------------------------------
10141727 CP020543.1 +     48 106157 106109 Ecoli_C +      0 106109 106107/106109 100.0  1
    5465 CP020543.1 + 121267 121367    100 Ecoli_C + 109317 109418    76/100     76.0   2
    4870 CP020543.1 + 159368 159512    144 Ecoli_C + 128706 128828    95/115     82.6   3
```

One of the fields chosen by us for [LASTZ run](#-hands-on-aligning-again) is `number`. This is an incrementing number given by LASTZ to every alignment block so it can be uniquely identified. The problem is that by running LASTZ on a collection on three genomes it generated number for each output independently starting with `1` each time. So these alignment identified are unique within each individual run but are redundant for multiple runs. We can fix that by pre-pending each alignment identified (column 13) with name of the target sequence (column 2). This would create alignment identified that are truly unique. For example, in case of LASTZ output shown above alignment identifier `1` will become `CP020543.11`, `2` will become `CP020543.12` and so on. Here is how we will do that:

> ### {% icon hands_on %} Hands-on: Creating unique alignment identifiers
> 
> 1. Open **Merge Columns together** tool.
> 2. In **Select data** choose the output of the previous step (called `Text transformation of collection ...`)
> 3. In **Merge column** select `Column: 2` (this is Targe sequence name)
> 4. In **with column** select `Column: 13` (this is the alignment block identified created by LASTZ)
> 5. Click **Execute**
{: .hands_on}

The output will look like this:

```
     1          2 3      4      5      6       7 8      9     10            11    12  13           14
-----------------------------------------------------------------------------------------------------
10141727 CP020543.1 +     48 106157 106109 Ecoli_C +      0 106109 106107/106109 100.0  1 CP020543.11 
    5465 CP020543.1 + 121267 121367    100 Ecoli_C + 109317 109418    76/100     76.0   2 CP020543.12 
    4870 CP020543.1 + 159368 159512    144 Ecoli_C + 128706 128828    95/115     82.6   3 CP020543.13 
```

the tool added a new column (Column 14) containing a merge between the target name and alignment id. 

Our next goals is to convert this into a format that will be acceptable to the genome browser [created above](#producing-a-genome-browser-for-this-experiment)


The columns were chosen by us [above](#-hands-on-aligning-again) and represent coordinates of alignment blocks identified by LASTZ. Because these are coordinates they can be easily visualized in the browser we just created. There is only problem though, to visualize these alignments we need to convert them into a format that the browser would understand. One of these formats is [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). In one of its simplest forms it has six columns:

 1. Chromosome ID
 2. Start
 3. End
 4. Name of the feature
 5. Score
 6. Strand (`+`, `-`, or `.` for no strand data)



Short introduction about this subpart.

> ### {% icon comment %} Comment
>
> Do you want to learn more about the principles behind mapping? Follow our [training](../../NGS-mapping)
{: .comment}

# Conclusion
{:.no_toc}

Conclusion about the technical key points. And then relation between the technical and the biological question to end with a global view.
