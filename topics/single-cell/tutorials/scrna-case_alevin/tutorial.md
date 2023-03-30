---
layout: tutorial_hands_on

title: "Generating a single cell matrix using Alevin"
subtopic: single-cell-CS
priority: 1

zenodo_link: 'https://zenodo.org/record/4574153'

redirect_from:
- /topics/transcriptomics/tutorials/droplet-quantification-preprocessing/tutorial
- /topics/transcriptomics/tutorials/scrna-case_alevin/tutorial

questions:
  - I have some single cell FASTQ files I want to analyse. Where do I start?

objectives:
  - Generate a cellxgene matrix for droplet-based single cell sequencing data
  - Interpret quality control (QC) plots to make informed decisions on cell thresholds
  - Find relevant information in GTF files for the particulars of their study, and include this in data matrix metadata

time_estimation: 2H

key_points:
  - Create a scanpy-accessible AnnData object from FASTQ files, including relevant gene metadata

tags:
  - single-cell
  - 10x
  - paper-replication
  - español
  - transcriptomics

contributions:
  authorship:
    - nomadscientist
    - pinin4fjords

  editing:
    - hexylena

  testing:
    - wee-snufkin

requirements:
  - type: "internal"
    topic_name: single-cell
    tutorials:
        - scrna-intro
        - scrna-umis

translations:
  - es

gitter: Galaxy-Training-Network/galaxy-single-cell

follow_up_training:
  -
    type: "internal"
    topic_name: single-cell
    tutorials:
        - scrna-case_alevin-combine-datasets
---

# Introduction


<!-- This is a comment. -->

This tutorial will take you from raw FASTQ files to a cell x gene data matrix in AnnData format. What's a data matrix, and what's AnnData format? Well you'll find out! Importantly, this is the first step in processing single cell data in order to start analysing it. Currently you have a bunch of strings of `ATGGGCTT` etc. in your sequencing files, and what you need to know is how many cells you have and what genes appear in those cells. These steps are the most computationally heavy in the single cell world, as you're starting with 100s of millions of reads, each with 4 lines of text. Later on in analysis, this data becomes simple gene counts such as 'Cell A has 4 GAPDHs', which is a lot easier to store! Because of this data overload, we have downsampled the FASTQ files to speed up the analysis a bit. Saying that, you're still having to map loads of reads to the massive murine genome, so get yourself a cup of coffee and prepare to analyse!

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Generating a matrix

In this section, we will show you the principles of the initial phase of single-cell RNA-seq analysis: generating expression measures in a matrix. We'll concentrate on droplet-based (rather than plate-based) methodology, since this is the process with most differences with respect to conventional approaches developed for bulk RNA-seq.

Droplet-based data consists of three components: cell barcodes, unique molecular identifiers (UMIs) and cDNA reads. To generate cell-wise quantifications we need to:

 * Process cell barcodes, working out which ones correspond to 'real' cells, which to sequencing artefacts, and possibly correct any barcodes likely to be the product of sequencing errors by comparison to more frequent sequences.
 * Map biological sequences to the reference genome or transcriptome.
 * 'De-duplicate' using the UMIs.

This used to be a complex process involving multiple algorithms, or was performed with technology-specific methods (such as 10X's 'Cellranger' tool)  but is now much simpler thanks to the advent of a few new methods. When selecting methodology for your own work you should consider:

 * [STARsolo](https://github.com/alexdobin/STAR) - a droplet-based scRNA-seq-specific variant of the popular genome alignment method STAR. Produces results very close to those of Cellranger (which itself uses STAR under the hood).
 * [Kallisto/ bustools](https://www.kallistobus.tools/) - developed by the originators of the transcriptome quantification method, Kallisto.
 * [Alevin](https://salmon.readthedocs.io/en/latest/alevin.html) - another transcriptome analysis method developed by the authors of the Salmon tool.

We're going to use Alevin {% cite article-Alevin %} for demonstration purposes, but we do not endorse one method over another.

## Get Data

We've provided you with some example data to play with, a small subset of the reads in a mouse dataset of fetal growth restriction {% cite Bacon2018 %} (see the [study in Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-6945/results/tsne) and the [project submission](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6945/)). This is a study using the Drop-seq chemistry, however this tutorial is almost identical to a 10x chemistry. We will point out the one tool parameter change you will need to run 10x samples. This data is not carefully curated, standard tutorial data - it's real, it's messy, it desperately needs filtering, it has background RNA running around, and most of all it will give you a chance to practice your analysis as if this data were yours.

Down-sampled reads and some associated annotation can be imported below. How did I downsample these FASTQ files? Check out [this history](https://humancellatlas.usegalaxy.eu/u/wendi.bacon.training/h/pre-processing-with-alevin---part-1---how-to-downsample) to find out!

Additionally, to map your reads, you will need a transcriptome to align against (a FASTA) as well as the gene information for each transcript (a gtf) file. You can download these for your species of interest [from Ensembl](https://www.ensembl.org/info/data/ftp/index.html). These files are included in the data import step below. Keep in mind, these are big files, so the fastest way to get these into your Galaxy account is through importing them by history.

> <hands-on-title>Option 1: Data upload - Import history</hands-on-title>
>
> 1. Import history from: [example input history](https://humancellatlas.usegalaxy.eu/u/wendi.bacon.training/h/cs1pre-processing-with-alevin---input-1)
>
>
>    {% snippet faqs/galaxy/histories_import.md %}
>
> 2. **Rename** {% icon galaxy-pencil %} the the history to your name of choice.
>
{: .hands_on}

> <hands-on-title>Option 2: Data upload - Add to history</hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the Experimental Design table, sequencing reads 1 & 2, the GTF and fasta files from [Zenodo]({{ page.zenodo_link }})
>
>    ```
>    {{ page.zenodo_link }}/files/Experimental_Design.tabular
>    {{ page.zenodo_link }}/files/Mus_musculus.GRCm38.100.gtf.gff
>    {{ page.zenodo_link }}/files/Mus_musculus.GRCm38.cdna.all.fa.fasta
>    {{ page.zenodo_link }}/files/SLX-7632.TAAGGCGA.N701.s_1.r_1.fq-400k.fastq
>    {{ page.zenodo_link }}/files/SLX-7632.TAAGGCGA.N701.s_1.r_2.fq-400k.fastq
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Rename {% icon galaxy-pencil %} the datasets
>
{: .hands_on}

> <question-title></question-title>
>
> Have a look at the files you now have in your history.
> 1. Which of the FASTQ files do you think contains the barcode sequences?
> 2. Given the chemistry this study should have, are the barcode/UMI reads the correct length?
> 3. What is the 'N701' referring to?
>
> > <solution-title></solution-title>
> >
> > 1. Read 1 (SLX-7632.TAAGGCGA.N701.s_1.r_1.fq-400k) contains the cell barcode and UMI because it is significantly shorter (indeed, 20 bp!) compared to the longer, r_2 transcript read. For ease, rename these files N701-Read1 and N701-Read2.
> > 2. You can see Read 1 is only 20 bp long, which for original Drop-Seq is 12 bp for cell barcode and 8 bp for UMI. This is correct! Be warned - 10x Chromium (and many technologies) change their chemistry over time, so particularly when you are accessing public data, you want to check and make sure you have your numbers correct!
> > 3. 'N701' is referring to an index read. This sample was run alongside 6 other samples, each denoted by an Illumina Nextera Index (N70X). Later, this will tell you batch information. If you look at the 'Experimental Design' file, you'll see that the N701 sample was from a male wildtype neonatal thymus.
> {: .solution}
{: .question}

# Important tips for easier analysis

{% snippet faqs/galaxy/tutorial_mode.md %}

> <comment-title></comment-title>
> - The Galaxy tool search panel sometimes doesn't find the tools we need from the thousands available.
> - You'll have a much easier time selecting tools from the panel (if you aren't using tutorial mode!) if you are on the [https://humancellatlas.usegalaxy.eu](https://humancellatlas.usegalaxy.eu)
{: .comment}

## Generate a transcript to gene map

Gene-level, rather than transcript-level, quantification is standard in scRNA-seq, which means that the expression level of alternatively spliced RNA molecules are combined to create gene-level values. Droplet-based scRNA-seq techniques only sample one end each transcript, so lack the full-molecule coverage that would be required to accurately quantify different transcript isoforms.

To generate gene-level quantifications based on transcriptome quantification, Alevin and similar tools require a conversion between transcript and gene identifiers. We can derive a transcript-gene conversion from the gene annotations available in genome resources such as Ensembl. The transcripts in such a list need to match the ones we will use later to build a binary transcriptome index. If you were using spike-ins, you'd need to add these to the transcriptome and the transcript-gene mapping.

In your example data you will see the murine reference annotation as retrieved from Ensembl in GTF format. This annotation contains gene, exon, transcript and all sorts of other information on the sequences. We will use these to generate the transcript/ gene mapping by passing that information to a tool that extracts just the transcript identifiers we need.

> <question-title></question-title>
>
> Which of the 'attributes' in the last column of the GTF files contains the transcript and gene identifiers?
>
>
>   > <tip-title>Hint</tip-title>
>   >
>   > The file is organised such that the last column (headed 'Group') contains a wealth of information in the format: attribute1 "information associated with attribute 1";attribute2 "information associated with attribute 2" etc.
>   {: .tip}
>
> > <solution-title></solution-title>
> > *gene_id* and *transcript_id* are each followed by "ensembl gene_id" and "ensembl transcript_id"
> {: .solution}
{: .question}

It's now time to parse the GTF file using the [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html) package in R. This parsing will give us a conversion table with a list of transcript identifiers and their corresponding gene identifiers for counting. Additionally, because we will be generating our own binary index (more later!), we also need to input our FASTA so that it can be filtered to only contain transcriptome information found in the GTF.

> <hands-on-title>Generate a filtered FASTA and transcript-gene map</hands-on-title>
>
> 1. {% tool [GTF2GeneList](toolshed.g2.bx.psu.edu/repos/ebi-gxa/gtf2gene_list/_ensembl_gtf2gene_list/1.52.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Ensembl GTF file"*: `GTF file in the history` {% icon galaxy-history %}
>    - *"Feature type for which to derive annotation"*: `transcript` (Your sequences are transcript sequencing, so this is your starting point)
>    - *"Field to place first in output table"*: `transcript_id` (This is accessing the column you identified above!)
>    - *"Suppress header line in output?"*: `Yes` (The next tool (Alevin) does not expect a header)
>    - *"Comma-separated list of field names to extract from the GTF (default: use all fields)"*: `transcript_id,gene_id` (This calls the first column to be the transcript_id, and the second the gene_id. Thus, your key can turn transcripts into genes)
>    - *"Append version to transcript identifiers?"*: `Yes` (The Ensembl FASTA files usually have these, and since we need the FASTA transcriptome and the GTF gene information to work together, we need to append these!)
>    - *"Flag mitochondrial features?"*: `No`
>    - *"Provide a cDNA file for extracting annotations and/ or possible filtering?"*: `Yes`
>    - {% icon param-file %} *"FASTA-format cDNA/transcript file"*: `FASTA file in your history` {% icon galaxy-history %}
>    - *"Annotation field to match with sequences"*: `transcript_id`
>    - *"Filter the cDNA file to match the annotations?"*: `Yes`
>
> 2. Rename {% icon galaxy-pencil %} the annotation table to `Map`
>
> 3. Rename {% icon galaxy-pencil %} the uncompressed filtered FASTA file to `Filtered FASTA`
{: .hands_on}

## Generate a transcriptome index & quantify!

Alevin collapses the steps involved in dealing with dscRNA-seq into a single process. Such tools need to compare the sequences in your sample to a reference containing all the likely transcript sequences (a 'transcriptome'). This will contain the biological transcript sequences known for a given species, and perhaps also technical sequences such as 'spike ins' if you have those.

> <details-title>How does Alevin work?</details-title>
>
> To be able to search a transcriptome quickly, Alevin needs to convert the text (FASTA) format sequences into something it can search quickly, called an 'index'. The index is in a binary rather than human-readable format, but allows fast lookup by Alevin. Because the types of biological and technical sequences we need to include in the index can vary between experiments, and because we often want to use the most up-to-date reference sequences from Ensembl or NCBI, we can end up re-making the indices quite often. Making these indices is time-consuming! Have a look at the uncompressed FASTA to see what it starts with.
>
{: .details}

We now have:

* Barcode/ UMI reads
* cDNA reads
* transcript/ gene mapping
* filtered FASTA

We can now run Alevin. In some public instances, Alevin won't show up if you search for it. Instead, you may have to click the Single Cell tab at the left and scroll down to the Alevin tool. Alternatively, use Tutorial Mode as described above and you'll easily navigate to all the tools, and their versions will all be the tried and tested ones of this tutorial. It's often a good idea to check your tool versions. In this case the tutorial was built with Alevin Galaxy Version 1.5.1+galaxy0. To identify which version of a tool you are using. Select {% icon tool-versions %} 'Versions' and choose the appropriate version.


> <hands-on-title>Running Alevin</hands-on-title>
>
> 1. {% tool [Alevin](toolshed.g2.bx.psu.edu/repos/bgruening/alevin/alevin/1.5.1+galaxy0) %}
>
>     > <question-title></question-title>
>     >
>     > Try to fill in the parameters of Alevin using what you know!
>     >
>     >   > <tip-title>Strandedness?</tip-title>
>     >   >
>     >   > The Salmon documentation on 'Fragment Library Types' and running the Alevin command is here: [salmon.readthedocs.io/en/latest/library_type.html](https://salmon.readthedocs.io/en/latest/library_type.html) and [salmon.readthedocs.io/en/latest/alevin.html](https://salmon.readthedocs.io/en/latest/alevin.html). These links will help here, although keep in mind the image there is drawn with the RNA 5' on top, whereas in this scRNA-seq protocol, the polyA is captured by its 3' tail and thus effectively the bottom or reverse strand...)
>     >   {: .tip}
>     >
>     >   > <solution-title></solution-title>
>     >   >    - *"Select a reference transcriptome from your history or use a built-in index?"*: `Use one from the history`
>     >   >       - You are going to generate the binary index using your filtered FASTA!
>     >   >    - {% icon param-file %} *"Transcripts FASTA file"*: `Filtered FASTA`
>     >   >    - *"Single or paired-end reads?"*: `Paired-end`
>     >   >    - {% icon param-file %} *"file1"*: `N701-Read1`
>     >   >    - {% icon param-file %} *"file2"*: `N701-Read2`
>     >   >    - *"Specify the strandedness of the reads"*: `Infer automatically (A)`
>     >   >    - *"Protocol"*: `DropSeq Single Cell protocol`
>     >   >    - {% icon param-file %} *"Transcript to gene map file"*: `Map`
>     >   >    - *"Retrieve all output files"*: `Yes`
>     >   >    - In *"Optional commands"*:
>     >   >        - *"dumpFeatures"*: `Yes`
>     >   >        - *"dumpMTX"*: `Yes`
>     >   {: .solution}
>     {: .question}
{: .hands_on}

> <comment-title>What if I'm running a 10x sample?</comment-title>
>
> The main parameter that needs changing for a 10X Chromium sample is the 'Protocol' parameter of Alevin. Just select the correct 10x Chemistry there instead.
{: .comment}

This tool will take a while to run. Alevin produces many file outputs, not all of which we'll use. You can refer to the [Alevin documentation](https://salmon.readthedocs.io/en/latest/alevin.html) if you're curious what they all are, but we're most interested in is:

* the matrix itself (quants_mat.mtx.gz - the count by gene and cell)
* the row (cell/ barcode) identifiers (quants_mat_rows.txt) and
* the column (gene) labels (quants_mat_cols.txt).

This is the matrix market (MTX) format.

> <question-title></question-title>
>
> After you've run Alevin, {% icon galaxy-eye %} look through all the different files. Can you find:
> 1. The Mapping Rate?
> 2. How many cells are present in the matrix output?
>
> > <solution-title></solution-title>
> >
> > 1. Inspect {% icon galaxy-eye %} the file {% icon param-file %} meta_info.json. You can see the mapping rate is a paltry `25.45%`. This is a terrible mapping rate. Why might this be? Remember this was downsampled, and specifically by taking only the last 400,000 reads of the FASTQ file. The overall mapping rate of the file is more like 50%, which is still quite poor, but for early Drop-Seq samples and single-cell data in general, you might expect a slightly poorer mapping rate. 10x samples are much better these days! This is real data, not test data, after all!
> > 2. Inspect {% icon galaxy-eye %} the file {% icon param-file %} 'quants_mat_rows.txt', and you can see it has `2163` lines. The rows refer to the cells in the cell x gene matrix. According to this (rough) estimate, your sample has 2163 cells in it!
> >
> {: .solution}
>
{: .question}

> <warning-title>Choose the appropriate input going forward!</warning-title>
> Make certain to use **quants_mat.mtx.gz**  and NOT **quants_tier.mtx.gz** going forward.
{: .warning}

{% icon congratulations %} Congratulations - you've made an expression matrix! We could almost stop here. But it's sensible to do some basic QC, and one of the things we can do is look at a barcode rank plot.

# Basic QC

The question we're looking to answer here, is: "do we mostly have a single cell per droplet"? That's what experimenters are normally aiming for, but it's not entirely straightforward to get exactly one cell per droplet. Sometimes almost no cells make it into droplets, other times we have too many cells in each droplet. At a minimum, we should easily be able to distinguish droplets with cells from those without.

> <hands-on-title>Generate a raw barcode QC plot</hands-on-title>
>
> 1. {% tool [Droplet barcode rank plot](toolshed.g2.bx.psu.edu/repos/ebi-gxa/droplet_barcode_plot/_dropletBarcodePlot/1.6.1+galaxy2) %} with the following parameters:
>    - *"Input MTX-format matrix?"*: `No`
>    - {% icon param-file %} *"A two-column tab-delimited file, with barcodes in the first column and frequencies in the second"*: `raw_cb_frequencies.txt`
>    - *"Label to place in plot title"*: `Barcode rank plot (raw barcode frequencies)`
>
> 2. Rename {% icon galaxy-pencil %} the image output `Barcode Plot - raw barcode frequencies`
{: .hands_on}

![raw droplet barcode plots-400k](../../images/scrna-casestudy/wab-raw_barcodes-400k.png "400k subsample raw")

Now, the image generated here (400k) isn't the most informative - but you are dealing with a fraction of the reads! If you run the total sample (so identical steps above, but with significantly more time!) you'd get the image below.

![raw droplet barcode plots-total](../../images/scrna-casestudy/wab-raw_barcodes-total.png "Total sample - 32,579,453 reads - raw")

This is our own formulation of the barcode plot based on a [discussion](https://github.com/COMBINE-lab/salmon/issues/362#issuecomment-490160480) we had with community members. The left hand plots with the smooth lines are the main plots, showing the UMI counts for individual cell barcodes ranked from high to low. We expect a sharp drop-off between cell-containing droplets and ones that are empty or contain only cell debris. Now, this data is not an ideal dataset, so for perspective, in an ideal world with a very clean 10x run, data will look a bit more like the following taken from the lung atlas (see the [study in Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-6653/results/tsne) and the [project submission](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6653/)).

![raw droplet barcode plots - lung atlas](../../images/scrna-casestudy/wab-lung-atlas-barcodes-raw.png "Pretty data - raw")

In that plot, you can see the clearer 'knee' bend, showing the cut-off between empty droplets and cell-containing droplets.

The right hand plots are density plots from the first one, and the thresholds are generated either using [dropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html) or by the method described in the discussion mentioned above. We could use any of these thresholds to select cells, assuming that anything with fewer counts is not a valid cell. By default, Alevin does something similar, and we can learn something about that by plotting just the barcodes Alevin retains.

> <hands-on-title>Generate Alevin's barcode QC plot</hands-on-title>
>
> 1. {% tool [Droplet barcode rank plot](toolshed.g2.bx.psu.edu/repos/ebi-gxa/droplet_barcode_plot/_dropletBarcodePlot/1.6.1+galaxy2) %} with the following parameters:
>    - *"Input MTX-format matrix?"*: `Yes`
>    - *"Matrix-market format matrix file, with cells by column (overrides --barcode-frequencies if supplied)"*: `quants_mat.mtx`
>    - *"For use with --mtx-matrix: force interpretation of matrix to assume cells are by row, rather than by column (default)"*: `Yes`
>    - *"Label to place in plot title"*: `Barcode rank plot (Alevin-processed)`
>
> 2. Rename {% icon galaxy-pencil %} the image output `Barcode Plot - Alevin processed barcodes`
{: .hands_on}

![raw droplet barcode plots - 400k](../../images/scrna-casestudy/wab-alevin-barcodes-400k.png "400k subsample - Alevin processed")

And the full sample looks like:

![raw droplet barcode plots - total](../../images/scrna-casestudy/wab-alevin-barcodes-total.png "Total sample - 32,579,453 reads - Alevin processed")

And to round this off, here's the lung atlas plot.

![raw droplet barcode plots - total](../../images/scrna-casestudy/wab-alevin-barcodes-lung.png "Pretty data - Alevin processed")

You should see a completely vertical drop-off where Alevin has trunctated the distribution (after excluding any cell barcode that had <10 UMI, Alevin then chose a threshold based off the curve and removed all barcodes with fewer UMIs).

In experiments with relatively simple characteristics, this 'knee detection' method works relatively well. But some populations (such as our sample!) present difficulties. For instance, sub-populations of small cells may not be distinguished from empty droplets based purely on counts by barcode. Some libraries produce multiple 'knees' for multiple sub-populations. The [emptyDrops](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y) method has become a popular way of dealing with this. emptyDrops still retains barcodes with very high counts, but also adds in barcodes that can be statistically distinguished from the ambient profiles, even if total counts are similar. In order to ultimately run emptyDrops (or indeed, whatever tool you like that accomplishes biologically relevant thresholding), we first need to re-run Alevin, but prevent it from applying its own less than ideal thresholds.

To use emptyDrops effectively, we need to go back and re-run Alevin, stopping it from applying it's own thresholds. Click the re-run icon {% icon galaxy-refresh %} on any Alevin output in your history, because almost every parameter is the same as before, except you need to change the following:


## Generate an unprocessed matrix in a usable format

> <hands-on-title>Stopping Alevin from thresholding</hands-on-title>
> 1. {% tool [Alevin](toolshed.g2.bx.psu.edu/repos/bgruening/alevin/alevin/1.5.1+galaxy0) %} (Click re-run on the last Alevin output)
>    - *"Optional commands"*
>    - *"keepCBFraction"*: '1' - i.e. keep them all!
>    - *"freqThreshold"*: '3' - This will only remove cell barcodes with a frequency of less than 3, a low bar to pass but useful way of avoiding processing a bunch of almost certainly empty barcodes
>
> {% snippet faqs/galaxy/tools_rerun.md %}
{: .hands_on}

> <question-title></question-title>
>
> How many cells are in the output now?
>
> > <solution-title></solution-title>
> >
> > 1. `22952` cells are in the quants_mat_rows now! Far more than the Alevin-filtered `2163`. This needs some serious filtering with EmptyDrops!
> >
> {: .solution}
>
{: .question}

Alevin outputs MTX format, which we can pass to the dropletUtils package and run emptyDrops. Unfortunately the matrix is in the wrong orientation for tools expecting files like those produced by 10X software (which dropletUtils does). We need to 'transform' the matrix such that cells are in columns and genes are in rows.

> <warning-title>Be careful!</warning-title>
> Don't mix up files from the different Alevin runs! Use the later run, which has higher numbers in the history!
{: .warning}

> <hands-on-title>Transform matrix</hands-on-title>
>
> 1. {% tool [salmonKallistoMtxTo10x](toolshed.g2.bx.psu.edu/repos/ebi-gxa/salmon_kallisto_mtx_to_10x/_salmon_kallisto_mtx_to_10x/0.0.1+galaxy6) %} with the following parameters:
>    - {% icon param-file %} *".mtx-format matrix"*: `quants_mat.mtx.gz` (output of **Alevin** {% icon tool %})
>    - {% icon param-file %} *"Tab-delimited genes file"*: `quants_mat_cols.txt` (output of **Alevin** {% icon tool %})
>    - {% icon param-file %} *"Tab-delimited barcodes file"*: `quants_mat_rows.txt` (output of **Alevin** {% icon tool %})
>
> 2. Rename {% icon galaxy-pencil %} 'salmonKallistoMtxTo10x....:genes' to `Gene table`
> 3. Rename {% icon galaxy-pencil %} 'salmonKallistoMtxTo10x....:barcodes' to `Barcode table`
> 4. Rename {% icon galaxy-pencil %} 'salmonKallistoMtxTo10x....:matrix' to `Matrix table`
{: .hands_on}

The output is a matrix in the correct orientation for the rest of our tools. However, our matrix is looking a bit sparse - for instance, click on `Gene table`. I don't know about you, but I'd struggle to have a good biological discussion using only Ensembl gene_ids! What I'd really like is the more understandable 'GAPDH' or other gene acronym, as well as information on mitochondrial genes so that I can assess if my cells were stressed out or not. In order to prepare our data for emptyDrops, we're going to combine this information into an object, and it's easiest to add in that information now.

## Adding in Gene metadata

> <question-title></question-title>
>
> Where can we find this gene information?
>
> > <solution-title></solution-title>
> >
> > Our old friend the GTF file!
> >
> {: .solution}
>
{: .question}

> <question-title></question-title>
>
> Which of the 'attributes' in the last column of that file contains the gene acronym?
>
> > <solution-title></solution-title>
> >
> > gene_name
> >
> {: .solution}
>
{: .question}

We're now going to re-run {% icon galaxy-refresh %} the tool that extracts information from our GTF file.

> <hands-on-title>Generate gene information</hands-on-title>
>
> 1. {% tool [GTF2GeneList](toolshed.g2.bx.psu.edu/repos/ebi-gxa/gtf2gene_list/_ensembl_gtf2gene_list/1.52.0+galaxy0) %} with the following parameters:
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

> <hands-on-title>Combine MTX Gene Table with Gene Information</hands-on-title>
>
> 1. {% tool [Join two Datasets](join1) %} with the following parameters:
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
>    If you inspect {% icon galaxy-eye %} the object, you'll see we have joined these tables and now have quite a few gene_id repeats. Let's take those out, while keeping the order of the original 'Gene Table'.
>
>
> 2.  {% tool [Cut columns from a table](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1,c4,c5`
>    - *"Delimited by"*: `Tab`
>    - *"From"*: output of **Join two Datasets** {% icon tool %}
>
> 3. Rename output `Annotated Gene Table`
{: .hands_on}

Inspect {% icon galaxy-eye %} your `Annotated Gene Table`. That's more like it! You now have `gene_id`, `gene_name`, and `mito`. Now let's get back to your journey to emptyDrops and sophisticated thresholding of empty droplets!

# emptyDrops

emptyDrops {% cite article-emptyDrops %} works with a specific form of R object called a SingleCellExperiment. We need to convert our transformed MTX files into that form, using the DropletUtils Read10x tool:

> <hands-on-title>Converting to SingleCellExperiment format</hands-on-title>
>
> 1. {% tool [DropletUtils Read10x](toolshed.g2.bx.psu.edu/repos/ebi-gxa/dropletutils_read_10x/dropletutils_read_10x/1.0.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Expression matrix in sparse matrix format (.mtx)"*: `Matrix table`
>    - {% icon param-file %} *"Gene Table"*: `Annotated Gene Table`
>    - {% icon param-file %} *"Barcode/cell table"*: `Barcode table`
>    - *"Should metadata file be added?"*: `No`
>
> 2. Rename {% icon galaxy-pencil %} output: `SCE Object`
{: .hands_on}

Fantastic! Now that our matrix is combined into an object, specifically the SingleCellExperiment format, we can now run emptyDrops! Let's get rid of those background droplets containing no cells!

> <hands-on-title>Emptydrops</hands-on-title>
>
> 1. {% tool [DropletUtils emptyDrops](toolshed.g2.bx.psu.edu/repos/ebi-gxa/dropletutils_empty_drops/dropletutils_empty_drops/1.0.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"SingleCellExperiment rdata object"*: `SCE Object`
>    - *"Should barcodes estimated to have no cells be removed from the output object?"*: `Yes`
>
> 2. Rename {% icon galaxy-pencil %} `serialised SingleCellExperiment` output as `Emptied-Object`
>
> 3. Rename {% icon galaxy-pencil %} `tabular output` as `Emptied-Tabular Output`
{: .hands_on}

> <question-title></question-title>
>
> How many cell barcodes remain after the emptyDrops treatment? Why might that be?
>
>   > <tip-title>Hint</tip-title>
>   > If you click on the `Emptied-Object` in the {% icon galaxy-history %} history, the text in that window says `37 barcodes` or something similar to that - there is an element of random in the algorithm, so yours might differ slightly. Why is this so low?? And why might the number be different?
>   > Consider...is this a complete set of data?
>   {: .tip}
>
>
> > <solution-title></solution-title>
> >
> > Remember this is a subsampled dataset. If you look carefully at the parameters of emptyDrops, you'll see it set a minimum threshold at 100 UMI. If you look at the barcode plots above for the 400k read sample, you'll see this is far too stringent for this subsampled data! To satisfy your curiosity, this minimum threshold would yield `4332` barcodes for the total sample. Also, the number may vary slightly as the output depends on a large number of random iterations.
> >
> {: .solution}
>
{: .question}

We will nevertheless proceed with your majestic annotated expression matrix of 37 cells, ready to go for further processing and analysis! However, the next tutorials we will link to use a tool suite called Scanpy {% cite Wolf2018 %}. You need to convert this SingleCellExperiment object into a format called `annData`, which is a variant of a file format called `hdf5`.

> <hands-on-title>Converting to AnnData format</hands-on-title>
>
> 1. {% tool [SCEasy convert](toolshed.g2.bx.psu.edu/repos/ebi-gxa/sceasy_convert/sceasy_convert/0.0.5+galaxy1) %} with the following parameters:
>    - *"Direction of conversion"*: `SingleCellExperiment to AnnData`
>    - {% icon param-file %} *"Input object in SingleCellExperiment RDS format"*: `Emptied-Object`
>    - *"Name of the assay to be transferred as main layer"*: `counts`
>
> 2. Rename {% icon galaxy-pencil %} output `N701-400k-AnnData`
>
{: .hands_on}

{% icon congratulations %} Congrats! Your object is ready to for the scanpy pipeline! You can can check your work against the [example history](https://humancellatlas.usegalaxy.eu/u/wendi.bacon.training/h/cs1pre-processing-with-alevin---answer-key), or compare how the subsampled datasets you've generated compare with the [total sample](https://humancellatlas.usegalaxy.eu/u/wendi.bacon.training/h/cs1pre-processing-with-alevin---total-n701-answer-key)

However, it may be that you want to combine this object with others like it, for instance, maybe you ran 5 samples, and you are starting with 10 FASTQ files...

# Analysing multiple FASTQ files

This sample was originally one of seven. So to run the other [12 downsampled FASTQ files](https://humancellatlas.usegalaxy.eu/u/wendi.bacon.training/h/alevin-tutorial---all-samples---400k), you can use a [workflow](https://humancellatlas.usegalaxy.eu/u/wendi.bacon.training/w/cs1generating-a-single-cell-matrix-using-alevin) and run them all at once! All these samples are going to take a while, so go and have several cups of tea... Or, better yet, I have [run them myself](https://humancellatlas.usegalaxy.eu/u/wendi.bacon.training/h/cs1generating-a-single-cell-matrix-using-alevin---all-downsample-fastq-processed-to-anndata). To combine the resultant files into a single matrix, you can look at the next tutorial in this case study: [Combining datasets after pre-processing]({% link topics/single-cell/tutorials/scrna-case_alevin-combine-datasets/tutorial.md %})

# Mitochondrial flagging

We have assumed you will be combining multiple files - but if that's not the case, you'll need to perform this step to turn your column of `true` and `false` labelling the mitochondrial genes into some metrics telling you the % of mitochondrial genes in each cell. You can follow that step here: [Mitochondrial calculations](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/scrna-case_alevin-combine-datsets#mitochondrial-reads).

# Conclusion


![Workflow Part 1](../../images/scrna-casestudy/wab-alevin-part1workflow.png "Workflow  - Steps 1-3")

We have:

 * Examined raw read data, annotations and necessary input files for quantification.
 * Run Alevin in two different parameterisations, both allowing Alevin to make its own calls on what constitutes empty droplets, and applying emptyDrops instead.
 * Deployed barcode rank plots as a way of quickly assessing the signal present in droplet datasets.
 * Applied the necessary conversion to pass these data to downstream processes.

 To discuss with like-minded scientists, join our Gitter channel for all things Galaxy-single cell!
 [![Gitter](https://badges.gitter.im/Galaxy-Training-Network/galaxy-single-cell.svg)](https://gitter.im/Galaxy-Training-Network/galaxy-single-cell?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
