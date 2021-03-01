---
layout: tutorial_hands_on

title: Generating a cell matrix using Alevin
subtopic: single-cell
priority: 9

zenodo_link: ''
questions:
- I have some single cell fastq files I want to analyse. Where do I start?
objectives:
- Repeat matrix generation for any droplet-based single cell sequencing data
- Apply data combination and metadata editing for particular experimental designs
- Interpret quality control (QC) plots to make informed decisions on cell thresholds
- Find relevant information in GTF files for the particulars of their study, and include this in data matrix metadata
time_estimation: 3H
key_points:
- Create a scanpy-accessible AnnData object from fastq files, including relevant cell and gene metadata
- Combine multiple samples and label according to study design
requirements:
-
    type: "internal"
    topic_name: transcriptomics
    tutorials:
        - scrna-intro
        - scrna-umis
tags:
- single-cell
- 10x
contributors:
- nomadscientist
- pinin4fjords

---

# Introduction
{:.no_toc}

<!-- This is a comment. -->

This tutorial will take you from raw fastq files to a cell x gene data matrix in AnnData format. What's a data matrix, and what's AnnData format? Well you'll find out! Importantly, this is the first step in processing single cell data in order to start analysing it. Currently you have a bunch of strings of ATGGGCTT etc. in your sequencing files, and what you need to know is how many cells you have and what genes appear in those cells. In the second part of this tutorial, we will also look at combining fastq files and adding in metadata (for instance, SEX or GENOTYPE) for analysis later on. These steps are the most computationally heavy in the single cell world, as you're starting with 100s of millions of reads, each 4 lines of text. Later on in analysis this data becomes simple gene counts such as 'Cell A has 4 GAPDHs', which is a lot easier to store! Because of this data overload, we have downsampled the fastq files to speed up the analysis a bit. Saying that, you're still having to map loads of reads to the massive murine genome, so get yourself a cup of coffee and prepare to analyse!

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# 1: Generating a matrix

In this section, we will show you the principles of the initial phase of single-cell RNA-seq analysis: generating expression measures in a matrix. We'll concentrate on droplet-based (rather than plate-based) methodology, since this is the process with most differences with respect to conventional approaches developed for bulk RNA-seq.

Droplet-based data consists of three components: cell barcodes, unique molecular identifiers (UMIs) and cDNA reads. To generate cell-wise quantifications we need to:

 * Process cell barcodes, working out which ones correspond to 'real' cells, which to sequencing artefacts, and possibly correct any barcodes likely to be the product of sequencing errors by comparison to more frequent sequences.
 * Map biological sequences to the reference genome or transcriptome.
 * 'De-duplicate' using the UMIs.

This used to be a complex process involving multiple algorithms, or was performed with technology-specific methods (such as 10X's 'Cellranger' tool)  but is now much simpler thanks to the advent of a few new methods. When selecting methodology for your own work you should consider:

 * [STARsolo](https://github.com/alexdobin/STAR) - a dscRNA-seq-specific variant of the popular genome alignment method STAR. Produces results very close to those of Cellranger (which itself uses STAR under the hood).
 * [Kallisto/ bustools](https://www.kallistobus.tools/) - developed by the originators of the transcriptome quantification method, Kallisto.
 * [Alevin](https://salmon.readthedocs.io/en/latest/alevin.html) - another transcriptome method developed by the authors of the Salmon tool.

We're going to use Alevin for demonstration purposes, but we do not endorse one method over another.

## 1.1 Get Data

We've provided you with some example data to play with, a small subset of the reads in a mouse dataset of fetal growth restriction (see the study in Single Cell Expression Atlas [here](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-6945/results/tsne) and the project submission [here](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6945/)). This is a study using the Drop-seq chemistry, however this tutorial is almost identical to a 10x chemistry. We will point out the one tool parameter change you will need to run 10x samples.

Down-sampled reads and some associated annotation can be downloaded from Zenodo below, or you can import this EXAMPLE INPUT HISTORY. How did I downsample these fastq files? Check out this history to find out! FIXME
Additionally, to map your reads, you will need a transcriptome to align against (a FASTA) as well as the gene information for each transcript (a gtf) file. You can download these for your species of interest from Ensembl [here](https://www.ensembl.org/info/data/ftp/index.html). Additionally, these files are available in the above history as well as the Zenodo links below. Keep in mind, these are big files, so they may take a bit to import!

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}


> ### {% icon question %} Questions
>
> Have a look at the files you now have in your history.
> 1. Which of the FASTQ files do you think contains the barcode sequences?
> 2. Given the chemistry this study should have, are the barcode/UMI reads the correct length?
> 3. What is the 'N701' referring to?
>
> > ### {% icon solution %} Solution
> >
> > 1. Read 1 (SLX-7632.TAAGGCGA.N701.s_1.r_1.fq-400k) contains the cell barcode and UMI because it is significantly shorter (indeed, 20 bp!) compared to the longer, r_2 transcript read. For ease, rename these files N701-Read1 and N701-Read2.
> > 2. You can see Read 1 is only 20 bp long, which for original Drop-Seq is 12 bp for cell barcode and 8 bp for UMI. This is correct! Be warned - 10x Chromium (and many technologies) change their chemistry over time, so particularly when you are accessing public data, you want to check and make sure you have your numbers correct!
> > 3. 'N701' is referring to an index read. This sample was run alongside 6 other samples, each denoted by an Illumina Nextera Index (N70X). Later, this will tell you batch information. If you look at the 'Experimental Design' file, you'll see that the N701 sample was from a male wildtype neonatal thymus.
> {: .solution}
{: .question}


## 1.2 Generate a transcript to gene map

Gene-level, rather than transcript-level, quantification is standard in scRNA-seq, which means that that the expression level of alternatively spliced RNA molecules are combined to create gene-level values. Droplet-based scRNA-seq techniques only sample one end each transcript, so lack the full-molecule coverage that would be required to accurately quantify different transcript isoforms.  

To generate gene-level quantifications based on transcriptome quantification, Alevin and similar tools require a conversion between transcript and gene identifiers. We can derive a transcript-gene conversion from the gene annotations available in genome resources such as Ensembl. The transcripts in such a list need to match the ones we will use later to build a binary transcriptome index. If you were using spike-ins, you'd need to add these to the transcriptome and the transcript-gene mapping.

In your example data you will see the murine reference annotation as retrieved from Ensembl in GTF format. This annotation contains gene, exon, transcript and all sorts of other information on the sequences. We will use these to generate the transcript/ gene mapping by passing that information to a tool that extracts just the transcript identifiers we need.

> ### {% icon question %} Questions
>
> Which of the 'attributes' in the last column of the GTF files contains the transcript and gene identifiers?
>
>
>   > ### {% icon tip %} Hint
>   >
>   > The file is organised such that the last column (headed 'Group') contains a wealth of information in the format: attribute1 "information associated with attribute 1";attribute2 "information associated with attribute 2" etc.
>   {: .tip}
>
> > ### {% icon solution %} Solution
> > *gene_id* and *transcript_id* are each followed by "ensembl gene_id" and "ensembl transcript_id"
> {: .solution}
{: .question}

It's now time to parse the GTF file using the [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html) package in R. This parsing will give us a conversion table with a list of transcript identifiers and their corresponding gene identifiers for counting. Additionally, because we will be generating our own binary index (more later!), we also need to input our FASTA so that it can be filtered to only contain transcriptome information found in the GTF.

> ### {% icon hands_on %} Hands-on: Generate a filtered FASTA and transcript-gene map
>
> 1. {% tool [GTF2GeneList](https://humancellatlas.usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/ebi-gxa/gtf2gene_list/_ensembl_gtf2gene_list/1.42.1+galaxy6) %} with the following parameters:
>    - {% icon param-file %} *"Ensembl GTF file"*: `GTF file in the history` {% icon galaxy-history %}
>    - *"Feature type for which to derive annotation"*: `transcript` (Your sequences are transcript sequencing, so this is your starting point)
>    - *"Field to place first in output table"*: `transcript_id` (This is accessing the column you identified above!)
>    - *"Suppress header line in output?"*: `Yes` (The next tool (Alevin) does not expect a header)
>    - *"Comma-separated list of field names to extract from the GTF (default: use all fields)"*: `transcript_id,gene_id` (This calls the first column to be the transcript_id, and the second the gene_id. Thus, your key can turn transcripts into genes)
>    - *"Append version to transcript identifiers?"*: `Yes` (The Ensembl FASTA files usually have these, and since we need the FASTA transcriptome and the GTF gene information to work together, we need to append these!)
>    - *"Flag mitochondrial features?"*: `No`
>    - *"Filter a FASTA-format cDNA file to match annotations?"*: `Yes`
>    - {% icon param-file %} *"FASTA-format cDNA/transcript file"*: `FASTA file in your history` {% icon galaxy-history %}
>    - *"Annotation field to match with sequences"*: `transcript_id`
> 2. Rename {% icon galaxy-pencil %} the annotation table to `Map`
>
> 3. Rename {% icon galaxy-pencil %} the uncompressed filtered FASTA file to `Filtered FASTA`
{: .hands_on}

## 1.3 Generate a transcriptome index & quantify!

Alevin collapses the steps involved in dealing with dscRNA-seq into a single process. Such tools need to compare the sequences in your sample to a reference containing all the likely transcript sequences (a 'transcriptome'). This will contain the biological transcript sequences known for a given species, and perhaps also technical sequences such as 'spike ins' if you have those.

> ### {% icon details %} How does Alevin work?
>
> To be able to search a transcriptome quickly, Alevin needs to convert the text (FASTA) format sequences into something it can search quickly, called an 'index'. The index is in a binary rather than human-readable format, but allows fast lookup by Alevin. Because the types of biological and technical sequences we need to include in the index can vary between experiments, and because we often want to use the most up-to-date reference sequences from Ensembl or NCBI, we can end up re-making the indices quite often. Making these indices is time-consuming! Have a look at the uncompressed FASTA to see what it starts with.
>
{: .details}

We now have:

* Barcode/ UMI reads
* cDNA reads
* transcript/ gene mapping
* filtered FASTA

We can now run Alevin. In some public instances, Alevin won't show up if you search for it. Instead, you have to click the Single Cell tab at the left and scroll down to the Alevin tool. Tip: If you click the tools from the tutorial option within Galaxy, you'll always have the correct version of the tool! In this case, it is: (Galaxy Version 0.14.1.2+galaxy1) - it should be default. If not, click 'Versions' and choose that version.   

![Clicking the tool](../../images/wab-tutorial-in-galaxy.png "Tutorial option at the top right in Galaxy")

![Accessing tutorial option within Galaxy](../../images/wab-tutorial-option-filler.png "Click the tool in the tutorial with Galaxy")

> ### {% icon hands_on %} Hands-on: Running Alevin
>
> 1. {% tool [Alevin](https://humancellatlas.usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/bgruening/alevin/alevin/1.3.0+galaxy2) %}
>
>   > ### {% icon question %} Questions
>   >
>   > Try to fill in the parameters of Alevin using what you know!
>   >
>   >   > ### {% icon tip %} Tip: Strandedness?
>   >   >
>   >   > The Salmon documentation on 'Fragment Library Types' and running the Alevin command (https://salmon.readthedocs.io/en/latest/library_type.html and https://salmon.readthedocs.io/en/latest/alevin.html) will help here, although keep in mind the image there is drawn with the RNA 5' on top, whereas in this scRNA-seq protocol, the polyA is captured by its 3' tail and thus effectively the bottom or reverse strand...)
>   >   {: .tip}
>   >
>   >   > ### {% icon solution %} Solution  
>   >   >    - *"Select a reference transcriptome from your history or use a built-in index?"*: `Use one from the history`
>   >   >       - You are going to generate the binary index using your filtered FASTA!
>   >   >    - {% icon param-file %} *"Transcripts FASTA file"*: `Filtered FASTA`
>   >   >    - *"Single or paired-end reads?"*: `Paired-end`
>   >   >    - {% icon param-file %} *"file1"*: `N701-Read1`
>   >   >    - {% icon param-file %} *"file2"*: `N701-Read2`
>   >   >    - *"Relative orientation of reads within a pair"*: `Mates are oriented towards each other (IR)`
>   >   >    - *"Specify the strandedness of the reads"*: `read comes from the reverse strand (SR)`
>   >   >    - *"Protocol"*: `DropSeq Single Cell protocol`
>   >   >    - {% icon param-file %} *"Transcript to gene map file"*: `Map`
>   >   >    - *"Retrieve all output files"*: `Yes`
>    - In *"Optional commands"*:
>        - *"dumpFeatures"*: `Yes`
>        - *"dumpMTX"*: `Yes`
>   >   {: .solution}
>   {: .question}
{: .hands_on}

> ### {% icon comment %} What if I'm running a 10x sample?
>
> The main parameter that needs changing for a 10X Chromium sample is the 'Protocol' parameter of Alevin. Just select the correct 10x Chemistry there instead. Additionally, under 'Optional Commands', you can input a 'whitelist file'. This is the file of all the cell barcodes there are in 10x beads. You can get this file here 'https://zenodo.org/record/3457880/files/3M-february-2018.txt.gz' although you ought to check for any changes in chemistry or updates in the whitelist file directly with 10x themselves.
{: .comment}

This tool will take a while to run. Alevin produces many file outputs, not all of which we'll use. You can refer to the [Alevin documentation](https://salmon.readthedocs.io/en/latest/alevin.html) if you're curious what they all are, but we're most interested in is:

    * the matrix itself (quants_mat.mtx.gz - the count by gene and cell)
    * the row (cell/ barcode) identifiers (quants_mat_rows.txt) and
    * the column (gene) labels (quants_mat_cols.txt).

This is the matrix market (MTX) format.

> ### {% icon question %} Questions
>
> After you've run Alevin, {% icon galaxy-eye %} look through all the different files. Can you find:
> 1. The Mapping Rate?
> 2. How many cells are present in the matrix output?
>
> > ### {% icon solution %} Solution
> >
> > 1. Inspect {% icon galaxy-eye %} the file {% icon param-file %} meta_info.json. You can see the mapping rate is a paltry `24.75%`. This is a terrible mapping rate. Why might this be? Remember this was downsampled, and specifically by taking only the last 400,000 reads of the fastq file. The overall mapping rate of the file is more like 50%, which is still quite poor, but for early Drop-Seq samples and single-cell data in general, you might expect a slightly poorer mapping rate. 10x samples are much better these days! This is real data, not test data, after all!
> > 2. Inspect {% icon galaxy-eye %} the file {% icon param-file %} 'quants_mat_rows.txt', and you can see it has `2163` lines. The rows refer to the cells in the cell x gene matrix. According to this (rough) estimate, your sample has 2163 cells in it!
> >
> {: .solution}
>
{: .question}

> ### {% icon warning %} Warning!: Choose the appropriate input going forward!
> Make certain to use **quants_mat.mtx.gz**  and NOT **quants_tier.mtx.gz** going forward.
{: .warning}

{% icon congratulations %} Congratulations - you've made an expression matrix! We could almost stop here. But it's sensible to do some basic QC, and one of the things we can do is look at a barcode rank plot.

# 2: Basic QC

The question we're looking to answer here, is: "do we have mostly a have a single cell per droplet"? That's what experimenters are normally aiming for, but it's not entirely straightforward to get exactly one cell per droplet. Sometimes almost no cells make it into droplets, other times we have too many cells in each droplet. At a minimum, we should easily be able to distinguish droplets with cells from those without.   

> ### {% icon hands_on %} Hands-on: Generate a raw barcode QC plot
>
> 1. {% tool [Droplet barcode rank plot](https://humancellatlas.usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/ebi-gxa/droplet_barcode_plot/_dropletBarcodePlot/1.6.1+galaxy2) %} with the following parameters:
>    - *"Input MTX-format matrix?"*: `No`
>    - {% icon param-file %} *"A two-column tab-delimited file, with barcodes in the first column and frequencies in the second"*: `raw_cb_frequencies.txt`
>    - *"Label to place in plot title"*: `Barcode rank plot (raw barcode frequencies)`
>
> 2. Rename {% icon galaxy-pencil %} the image output `Barcode Plot - raw barcode frequencies`
{: .hands_on}

![raw droplet barcode plots-400k](../../images/wab-raw_barcodes-400k.png "400k subsample raw")

Now, the image generated here (400k) isn't the most informative - but you are dealing with a fraction of the reads! If you run the total sample (so identical steps above, but with significantly more time!) you'd get the image below.

![raw droplet barcode plots-total](../../images/wab-raw_barcodes-total.png "Total sample - 32,579,453 reads - raw")

This is our own formulation of the barcode plot based on a [discussion](https://github.com/COMBINE-lab/salmon/issues/362#issuecomment-490160480) we had with community members. The left hand plots with the smooth lines are the main plots, showing the UMI counts for individual cell barcodes ranked from high to low. We expect a sharp drop-off between cell-containing droplets and ones that are empty or contain only cell debris. Now, this data is not an ideal dataset, so for perspective, in an ideal world with a very clean 10x run, data will look a bit more like the following (see the study in Single Cell Expression Atlas [here](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-6653/results/tsne) and the project submission [here](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6653/)).

![raw droplet barcode plots - lung atlas](../../images/wab-lung-atlas-barcodes-raw.png "Pretty data - raw")

In that plot, you can see the clearer 'knee' bend, showing the cut-off between empty droplets and cell-containing droplets.

The right hand plots are density plots from the first one, and the thresholds are generated either using [dropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html) or by the method described in that discussion. We could use any of these thresholds to select cells, assuming that anything with fewer counts is not a valid cell. By default, Alevin does something similar, and we can learn something about that by plotting just the barcodes Alevin retains.

> ### {% icon hands_on %} Hands-on: Generate Alevin's barcode QC plot
>
> 1. {% tool [Droplet barcode rank plot](https://humancellatlas.usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/ebi-gxa/droplet_barcode_plot/_dropletBarcodePlot/1.6.1+galaxy2) %} with the following parameters:
>    - *"Input MTX-format matrix?"*: `Yes`
>    - *"Matrix-market format matrix file, with cells by column (overrides --barcode-frequencies if supplied)"*: `quants_mat.mtx`
>    - *"For use with --mtx-matrix: force interpretation of matrix to assume cells are by row, rather than by column (default)"*: `Yes`
>    - *"Label to place in plot title"*: `Barcode rank plot (Alevin-processed)`
>    - *"Number of bins used in barcode count frequency distribution"*: `50`
>    - *"Above-baseline multiplier to calculate roryk threshold"*: `1.5`
>
> 2. Rename {% icon galaxy-pencil %} the image output `Barcode Plot - Alevin processed barcodes`
{: .hands_on}

![raw droplet barcode plots - 400k](../../images/wab-alevin-barcodes-400k.png "400k subsample - Alevin processed")

And the full sample looks like:

![raw droplet barcode plots - total](../../images/wab-alevin-barcodes-total.png "Total sample - 32,579,453 reads - Alevin processed")

And to round this off, here's the lung atlas plot.

![raw droplet barcode plots - total](../../images/wab-alevin-barcodes-lung.png "Pretty data - Alevin processed")

You should see a completely vertical drop-off where Alevin has trunctated the distribution (after excluding any cell barcode that had <10 UMI, Alevin then chose a threshold based off the curve and removed all barcodes with fewer UMIs).

In experiments with relatively simple characteristics, this 'knee detection' method works relatively well. But some populations present difficulties due to sub-populations of small cells that cannot be distinguished from empty droplets based purely on counts by barcode. Some libraries produce multiple 'knees' for multiple sub-populations. The [emptyDrops](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y) method has become a popular way of dealing with this. emptyDrops still retains barcodes with very high counts, but also adds in barcodes that can be statistically distinguished from the ambient profiles, even if total counts are similar. In order to ultimately run emptyDrops (or indeed, whatever tool you like that accomplishes biologically relevant thresholding), we first need to re-run Alevin, but prevent it from applying its own less than ideal thresholds.

To use emptyDrops effectively, we need to go back and re-run Alevin, stopping it from applying it's own thresholds. Click the re-run icon {% icon galaxy-refresh %} on any Alevin output in your history, because almost every parameter is the same as before, except you need to change the following:

## 2.1 Generate an unprocessed matrix in a usable format

> ### {% icon hands_on %} Hands-on: Stopping Alevin from thresholding
> 1. {% tool [Alevin](https://humancellatlas.usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/bgruening/alevin/alevin/1.3.0+galaxy2) %}
>   - *"Optional commands"*
>   - *"keepCBFraction"*: '1' - i.e. keep them all!
>   - *"freqThreshold"*: '3' - This will only remove cell barcodes with a frequency of less than 3, a low bar to pass but useful way of avoiding processing a bunch of almost certainly empty barcodes
{: .hands_on}

> ### {% icon question %} Question
>
> How many cells are in the output now?
>
> > ### {% icon solution %} Solution
> >
> > 1. `22539` cells are in the quants_mat_rows now! Far more than the Alevin-filtered `2163`. This needs some serious filtering with EmptyDrops!
> >
> {: .solution}
>
{: .question}

Alevin outputs MTX format, which we can pass to the dropletUtils package and run emptyDrops. Unfortunately the matrix is in the wrong orientation for tools expecting files like those produced by 10X software (which dropletUtils does). We need to 'transform' the matrix such that cells are in columns and genes are in rows.

> ### {% icon warning %} Be careful!
> Don't mix up files from the different Alevin runs! Use the later run, which has higher numbers in the history!
{: .warning}

> ### {% icon hands_on %} Hands-on: Transform matrix
>
> 1. {% tool [salmonKallistoMtxTo10x](https://humancellatlas.usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/ebi-gxa/salmon_kallisto_mtx_to_10x/_salmon_kallisto_mtx_to_10x/0.0.1+galaxy5) %} with the following parameters:
>    - {% icon param-file %} *".mtx-format matrix"*: `quants_mat.mtx.gz` (output of **Alevin** {% icon tool %})
>    - {% icon param-file %} *"Tab-delimited genes file"*: `quants_mat_cols.txt` (output of **Alevin** {% icon tool %})
>    - {% icon param-file %} *"Tab-delimited barcodes file"*: `quants_mat_rows.txt` (output of **Alevin** {% icon tool %})
>
> 2. Rename {% icon galaxy-pencil %} 'salmonKallistoMtxTo10x....:genes' to `Gene table`
> 3. Rename {% icon galaxy-pencil %} 'salmonKallistoMtxTo10x....:barcodes' to `Barcode table`
> 4. Rename {% icon galaxy-pencil %} 'salmonKallistoMtxTo10x....:matrix' to `Matrix table`
{: .hands_on}

The output is a matrix in the correct orientation for the rest of our tools. However, our matrix is looking a bit sparse - for instance, click on `Gene table`. I don't know about you, but I'd struggle to have a good biological discussion using only Ensembl gene_ids! What I'd really like is the more understandable 'GAPDH' or other gene acronym, as well as information on mitochondrial genes so that I can assess if my cells were stressed out or not. In order to prepare our data for emptyDrops, we're going to combine this information into an object, and it's easiest to add in that information now.

## 2.2 Adding in Gene metadata

> ### {% icon question %} Question
>
> Where can we find this gene information?
>
> > ### {% icon solution %} Solution
> >
> > Our old friend the GTF file!
> >
> {: .solution}
>
{: .question}

> ### {% icon question %} Question
>
> Which of the 'attributes' in the last column of that file contains the gene acronym?
>
> > ### {% icon solution %} Solution
> >
> > gene_name
> >
> {: .solution}
>
{: .question}

We're now going to re-run {% icon galaxy-refresh %} the tool that extracts information from our GTF file.

> ### {% icon hands_on %} Hands-on: Generate gene information
>
> 1. {% tool [GTF2GeneList](https://humancellatlas.usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/ebi-gxa/gtf2gene_list/_ensembl_gtf2gene_list/1.42.1+galaxy6) %} with the following parameters:
>    - *"Feature type for which to derive annotation"*: `gene`
>    - *"Field to place first in output table"*: `gene_id`
>    - *"Suppress header line in output?"*: `Yes`
>    - *"Comma-separated list of field names to extract from the GTF (default: use all fields)"*: `gene_id,gene_name,mito`
>    - *"Append version to transcript identifiers?"*: `Yes`
>    - *"Flag mitochondrial features?"*: `Yes` - note, this will auto-fill a bunch of acronyms for searching in the GTF for mitochondrial associated genes. This is good!
>    - *"Filter a FASTA-format cDNA file to match annotations?"*: `No` - we don't need to, we're done with the FASTA!
> 2. Check that the output file type is `tabular`. If not, change the file type by clicking the 'Edit attributes'{% icon galaxy-pencil %} on the dataset in the history (as if you were renaming the file.) Then click `Datatypes` and type in `tabular`. Click `Change datatype`.)
> 2. Rename {% icon galaxy-pencil %} the annotation table to `Gene Information`
{: .hands_on}

Inspect {% icon galaxy-eye %} the **Gene Information** object in the history. Now you have made a new key for gene_id, with gene name and a column of mitochondrial information (false = not mitochondrial, true = mitochondrial). We need to add this information into the salmonKallistoMtxTo10x output 'Gene table'. But we need to keep 'Gene table' in the same order, since it is referenced in the 'Matrix table' by row.

> ### {% icon hands_on %} Hands-on: Combine MTX Gene Table with Gene Information
>
> 1. {% tool [Join two Datasets](https://humancellatlas.usegalaxy.eu/root?tool_id=join1) %} with the following parameters:
>    - *"Join"*: `Gene Table`
>    - *"Using column"*: `Column: 1`
>    - *"with"*: `Gene Information`
>    - *"and column"*: `Column: 1`
>    - *"Keep lines of first input that do not join with second input"*: `Yes`
>    - *"Keep lines of first input that are incomplete"*: `Yes`
>    - *"Fill empty columns"*: `No`
>    - *"Keep the header lines"*: `No`
>
>
> If you inspect {% icon galaxy-eye %} the object, you'll see we have joined these tables and now have quite a few gene_id repeats. Let's take those out, while keeping the order of the original 'Gene Table'.
>
>
> 2.  {% tool [Cut columns from a table](https://humancellatlas.usegalaxy.eu/root?tool_id=Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1,c4,c5`
>    - *"Delimited by"*: `Tab`
>    - *"From"*: output of **Join two Datasets** {% icon tool %}
>
> 3. Rename output `Annotated Gene Table`
{: .hands_on}

Inspect {% icon galaxy-eye %} your `Annotated Gene Table`. That's more like it! You now have `gene_id`, `gene_name`, and `mito`. Now let's get back to your journey to emptyDrops and sophisticated thresholding of empty droplets!

# 3: emptyDrops

emptyDrops works with a specific form of R object called a SingleCellExperiment. We need to convert our transformed MTX files into that form, using the DropletUtils Read10x tool:

> ### {% icon hands_on %} Hands-on: Converting to SingleCellExperiment format
>
> 1. {% tool [DropletUtils Read10x](https://humancellatlas.usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/ebi-gxa/dropletutils_read_10x/dropletutils_read_10x/1.0.3+galaxy2){% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Expression matrix in sparse matrix format (.mtx)"*: `Matrix table`
>    - {% icon param-file %} *"Gene Table"*: `Annotated Gene Table`
>    - {% icon param-file %} *"Barcode/cell table"*: `Barcode table`
>    - *"Should metadata file be added?"*: `No`
>
> 2. Rename {% icon galaxy-pencil %} output: `SCE Object`
{: .hands_on}

Fantastic! Now that our matrix is combined into an object, specifically the SingleCellExperiment format, we can now run emptyDrops! Let's get rid of those background droplets containing no cells!

> ### {% icon hands_on %} Hands-on: Emptydrops
>
> 1. {% tool [DropletUtils emptyDrops](https://humancellatlas.usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/ebi-gxa/dropletutils_empty_drops/dropletutils_empty_drops/1.0.3+galaxy1){% icon tool %} with the following parameters:
>    - {% icon param-file %} *"SingleCellExperiment rdata object"*: `SCE Object`
>    - *"Should barcodes estimated to have no cells be removed from the output object?"*: `Yes`
>
> 2. Rename {% icon galaxy-pencil %} `serialised SingleCellExperiment` output as `Stringent-Object`
>
> 3. Rename {% icon galaxy-pencil %} `tabular output` as `Stringent-Tabular Output`
{: .hands_on}

> ### {% icon question %} Question
>
> How many cell barcodes remain after the emptyDrops treatment? Why might that be?
>
>   > ### {% icon tip %} Hint
>   > If you click on the `Stringent-Object` in the {% icon galaxy-history %} history, the text in that window says `22 barcodes`. Why is this so low??
>   > Consider...is this a complete set of data?
>   {: .tip}
>
>
> > ### {% icon solution %} Solution
> >
> > Remember this is a subsampled dataset. If you look carefully at the parameters of emptyDrops, you'll see it set a minimum threshold at 100 UMI. If you look at the barcode plots above for the 400k read sample, you'll see this is far too stringent for this subsampled data! To satisfy your curiosity, this minimum threshold would yield `3011` barcodes for the total sample.
> >
> {: .solution}
>
{: .question}

Let's go back and tweak parameters, re-running the tool with a looser threshold minimum.

> ### {% icon details %} Decision-time - UMI count lower bound
> If you are working in a group, you can now divvy up a decision here with one control and the rest varied numbers so that you can compare results throughout the tutorials.
> - Variable: **UMI count lower bound**
> - Control:  `5`
> - Everyone else: Consider the droplet barcode rank plots and choose (different) appropriate lower bounds.
{: .details}

> ### {% icon hands_on %} Hands-on: emptyDrops - do-over!
>
> 1. {% tool [DropletUtils emptyDrops](https://humancellatlas.usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/ebi-gxa/dropletutils_empty_drops/dropletutils_empty_drops/1.0.3+galaxy1){% icon tool %} with the following parameters:
>    - {% icon param-file %} *"SingleCellExperiment rdata object"*: `SCE Object`
>    - *"UMI count lower bound"*: `5` - you can input different numbers here and see what happens!
>    - *"Should barcodes estimated to have no cells be removed from the output object?"*: `Yes`
>
> 2. Rename {% icon galaxy-pencil %} `serialised SingleCellExperiment` output as `<yournumberhere>UMI-Object`
>
> 3. Rename {% icon galaxy-pencil %} `tabular output` as `<yournumberhere>UMI-Tabular Output`
{: .hands_on}

You should now have `111` barcodes! You now have an annotated expression matrix ready to go for further processing and analysis! Well done! However, the next tutorials we will link to use a tool called Scanpy. You need to convert this SingleCellExperiment object into a format called `annData`, which is a variant of a file format called `hdf5`.

> ### {% icon hands_on %} Hands-on: Converting to AnnData format
>
> 1. {% tool [SCEasy convert](https://humancellatlas.usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/ebi-gxa/sceasy_convert/sceasy_convert/0.0.5+galaxy1){% icon tool %} with the following parameters:
>    - *"Direction of conversion"*: `SingleCellExperiment to AnnData`
>    - {% icon param-file %} *"Input object in SingleCellExperiment RDS format"*: `<yournumberhere>UMI-Object`
>    - *"Name of the assay to be transferred as main layer"*: `counts`
>
> 2. Rename {% icon galaxy-pencil %} output `N701-400k-AnnData`
>
{: .hands_on}

{% icon congratulations %} Congrats! Your object is ready to for the scanpy pipeline! However, it may be that you want to combine this object with others like it, for instance, maybe you ran 5 samples, and you are starting with 10 fastq files...

# 4: Combining fastq files

This sample was originally one of seven. So to run the other 12 fastq files, I strongly suggest you use the workflow (LINK!) to run each of the other pairs! It's going to take a while, so go and have a cup of tea... Or, for ease, I have downsampled the other 6 pairs of fastq files associated with this data, ran them myself (history here), and plopped them in a new clean history for you to import alongside the result you have already created above. You can access it here or get data with zenodo.


## 4.1 Data
> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial (if you're not importing the history above)
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

Inspect the {% icon galaxy-eye %} 'Experimental Design' text file. This shows you how each N70X corresponds to a sample, and whether that sample was from a male or female. This will be important metadata to add to our sample, which we will add very similarly to how you added the gene_name and mito metadata above!

## 4.2 Concatenating Objects
> ### {% icon hands_on %} Hands-on: Concatenating AnnData objects
>
> 1. {% tool [Manipulate AnnData](https://humancellatlas.usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.6.22.post1+galaxy4){% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: 'FIXME'
>    - *"Function to manipulate the object"*: 'Concatenate along the observations axis'
>    - {% icon param-file %} *"Annotated data matrix to add"*: 'Select all the other matrix files FIXME'
>    - *"Join method"*: 'Intersection of variables'
>    - *"Key to add the batch annovation to obs"*: `batch'
>    - *"Separator to join the existing index names with the batch category"*: '-'
> 2. Rename the object 'Combined AnnData object'
{: .hands_on}

Now let's look at what we've done! Unfortunately, AnnData objects are quite complicated, so the {% icon galaxy-eye %} won't help us too much here. Instead, we're going to use a tool to look into our object from now on.

> ### {% icon hands_on %} Hands-on: Inspecting AnnData Objects
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.6.22.post1+galaxy4) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Combined AnnData object'
>    - *"What to inspect?"*: `General information about the object`
> 2. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.6.22.post1+galaxy4) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Combined AnnData object'
>    - *"What to inspect?"*: `Key-indexed observations annotation (obs)`
> 3. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.6.22.post1+galaxy4) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Combined AnnData object'
>    - *"What to inspect?"*: `Key-indexed annotation of variables/features (var)`
{: .hands_on}

Now have a look at the three {% icon tool %} **Inspect AnnData** outputs.

> ### {% icon question %} Question
>
> 1. How many cells do you have now?
> 2. Where is 'batch' information stored?
>
> > ### {% icon solution %} Solution
> >
> > 1. FIXME
> > 2. Under 'Key-indexed observations annotation (obs)'. Batch refers to the order in which the matrices were added.
> > FIXME show image of how they are added
> {: .solution}
>
{: .question}

# 5: Adding metadata

Importantly, the batch numbers provided here are not necessarily in the same order as your samples N701,2,3,4 etc. Rather, the first input you had will be considered 'batch 0', and then the lowest object you added in the 'Annotated data matrix to add' will be 'batch 1', etc. I've set this up to run so it's in order (to save my brain some processing time!), so for how I ran it, it looks like this:

 FIXME

If you ran yours in a different order, update the batching parameters below to suit!
The two critical pieces of metadata in this experiment are **sex** and **genotype**. I will later want to color my cell plots by these parameters, so I want to add them in now!

> ### {% icon hands_on %} Hands-on: Labelling sex
>
> 1. {% tool [Replace Text in a specific column](https://humancellatlas.usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `Key-indexed observations annotation (obs) (output of **Inspect AnnData: Key-indexed observations annotation (obs)** {% icon tool %})
>    - *"1. Replacement"*
>          - *"in column"*: 'Column: 3'
>          - *"Find pattern"*: 'FIXME'
>          - *"Replace with"*: 'male'
>    - Select **Insert Replacement**
>    - *"2. Replacement"*
>          - *"in column"*: 'Column: 3'
>          - *"Find pattern"*: 'FIXME'
>          - *"Replace with"*: 'female'
>    - Select **Insert Replacement**
>    - *"3. Replacement"*
>          - *"in column"*: 'Column: 3'
>          - *"Find pattern"*: 'batch'
>          - *"Replace with"*: 'sex'
> Now we want only the column containing the sex information - we will ultimately add this into the cell annotation in the AnnData object.
> 2. {% tool [Cut columns from a table](https://humancellatlas.usegalaxy.eu/root?tool_id=Cut1) %} with the following parameters:
>    - *"Cut columns"*: 'c3'
>    - *"Delimited by"*: 'Tab'
>    - % icon param-file %} *"From"*: (output of **Replace text** {% icon tool %})
> 3. Rename output 'Sex metadata'
{: .hands_on}

That was so fun, let's do it all again but for genotype!

> ### {% icon hands_on %} Hands-on: Labelling genotype
>
> 1. {% tool [Replace Text in a specific column](https://humancellatlas.usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `Key-indexed observations annotation (obs) (output of **Inspect AnnData: Key-indexed observations annotation (obs)** {% icon tool %})
>    - *"1. Replacement"*
>          - *"in column"*: 'Column: 3'
>          - *"Find pattern"*: 'FIXME'
>          - *"Replace with"*: 'wildtype'
>    - Select **Insert Replacement**
>    - *"2. Replacement"*
>          - *"in column"*: 'Column: 3'
>          - *"Find pattern"*: 'FIXME'
>          - *"Replace with"*: 'knockout'
>    - Select **Insert Replacement**
>    - *"3. Replacement"*
>          - *"in column"*: 'Column: 3'
>          - *"Find pattern"*: 'batch'
>          - *"Replace with"*: 'genotype'
> Now we want only the column containing the genotype information - we will ultimately add this into the cell annotation in the AnnData object.
> 2. {% tool [Cut columns from a table](https://humancellatlas.usegalaxy.eu/root?tool_id=Cut1) %} with the following parameters:
>    - *"Cut columns"*: 'c3'
>    - *"Delimited by"*: 'Tab'
>    - % icon param-file %} *"From"*: (output of **Replace text** {% icon tool %})
> 3. Rename output 'Genotype metadata'
{: .hands_on}

You might want to do this with all sorts of different metadata - which labs handled the samples, which days they were run, etc. Once you've added created all your metadata columns, we can add them together before plugging them into the AnnData object itself.

> ### {% icon hands_on %} Hands-on: Combining metadata columns
>
> 1. {% tool [Paste two files side by side](https://humancellatlas.usegalaxy.eu/root?tool_id=Paste1) %} with the following parameters:
>    - {% icon param-file %} *"Paste"*: `Genotype metadata'
>    - {% icon param-file %} *"and"*: 'Sex metadata'
>    - *"Delimit by"*: 'Tab'
> 2. Rename - 'Cell Metadata'
{: .hands_on}

Let's add it to the AnnData object!

> ### {% icon hands_on %} Hands-on: Adding metadata to AnnData object
>
> 1. {% tool [Manipulate AnnData](https://humancellatlas.usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.6.22.post1+galaxy4) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Combined AnnData object'
>    - *"Function to manipulate the object"*: `Add new annotation(s) for observations or variables`
>    - *"What to annotate?"*: 'Observations (obs)'
>    - {% icon param-file %} *"Table with new annotations"*: 'Cell Metadata'
{: .hands_on}

Woohoo! We're there! You can run an **Inspect AnnData** to check now, but I want to clean up this AnnData object just a bit more first. It would be a lot nicer if 'batch' meant something, rather than 'the order in which the Manipulate AnnData tool added my datasets'.

> ### {% icon hands_on %} Hands-on: Labelling batches
>
> 1. {% tool [Manipulate AnnData](https://humancellatlas.usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.6.22.post1+galaxy4) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: (output of **Manipulate AnnData - Add new annotations** {% icon tool %})
>    - *"Function to manipulate the object"*: `Rename categories of annotation`
>    - *"Key for observations or variables annotation"*: 'batch'
>    - *"Comma-separated list of new categories"*: FIXME 'N701,N707,N706,N705,N704,N703,N702'
> 2. Rename the output 'Batched Object'
{: .hands_on}

Huzzah! We are JUST about there. However, while we've been focussing on our cell metadata (sample, batch, genotype, etc.) to relabel the 'obs' in our object.. remember when I mentioned mitochondria ages ago? And how often in single cell samples, mitochondrial RNA is often an indicator of stress during dissociation? We should probably do something with our column of true/false in the gene annotation that tells us information about the cells.
> ### {% icon hands_on %} Hands-on: Calculating mitochondrial RNA in cells
>
> 1. {% tool [AnnData Operations](https://humancellatlas.usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/ebi-gxa/anndata_ops/anndata_ops/0.0.3+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input object in hdf5 AnnData format"*: 'Batched object'
>    - *"Function of output object"*: `AnnData format`
>    - *"Copy AnnData to .raw"*: 'No'
>    - *"Gene symbols field in AnnData"*: 'NA.'
>    - *"Flag genes that start with these names"*: 'Insert Flag genes that start with these names'
>    - *"Starts with"*: 'True'
>    - *"Var name"*: 'mito'
>    - *"Number of top genes"*: '50'
> 2. Rename output 'Pre-processed object'
{: .hands_on}

Well done!  I strongly suggest have a play with the **Inspect AnnData** {% icon tool %} on your final 'Pre-processed object' to see the wealth of information that has been added. You are now ready to move along to the next tutorial FIXME

## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
