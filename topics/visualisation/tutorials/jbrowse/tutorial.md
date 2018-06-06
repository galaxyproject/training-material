---
layout: tutorial_hands_on
topic_name: visualisation
tutorial_name: jbrowse
---

# Introduction
{:.no_toc}

> JBrowse is a fast, embeddable genome browser built completely with JavaScript
> and HTML5, with optional run-once data formatting tools written in Perl.

The Galaxy tool accepts data in many formats:

- Intervals/Feature Tracks (GFF/GFF3, BED, GenBank)
- BAM Pileups
- Blast XML results
- Wig/BigWig
- VCF files

and executes the "run-once data formatting tools" mentioned in its description. The JBrowse tool has an incredibly extensive number of options, more than anyone needs most of the time. We'll go through them in detail but feel free to skip the sections that don't apply to the data types you use. Not everyone has Blast XML results to visualise.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Preparation

## Tool Installation

This tutorial covers versions 0.7.0.3 or greater of the JBrowse tool.

## Data Upload

The data for today is a subset of a real dataset from a Staphylococcus aureus bacteria.
We have a closed genome sequence and an annotation for our "wildtype" strain.
We have used a whole genome shotgun approach to produce a set of short sequence reads on an Illumina DNA sequencing instrument for our mutant strain.

- The reads are paired-end
- Each read is on average 150 bases
- The reads would cover the original wildtype genome to a depth of 19x

The files we will be using are:

- `mutant_R1.fastq` & `mutant_R2.fastq` - the read files in fastq format.
- `wildtype.fna` - The sequence of the reference strain in fasta format.
- `wildtype.gbk` - The reference strain with gene and other annotations in genbank format.
- `wildtype.gff` - The reference strain with gene and other annotations in gff3 format.

This data is available at Zenodo using the following [link](https://doi.org/10.5281/zenodo.582600).
> ### {% icon hands_on %} Hands-on: Get the data
>
> 1.  Import all of the following files into a new history:
>     - [mutant_R1.fastq](https://zenodo.org/record/582600/files/mutant_R1.fastq)
>     - [mutant_R2.fastq](https://zenodo.org/record/582600/files/mutant_R2.fastq)
>     - [wildtype.fna](https://zenodo.org/record/582600/files/wildtype.fna)  
>     - [wildtype.gbk](https://zenodo.org/record/582600/files/wildtype.gbk)
>     - [wildtype.gff](https://zenodo.org/record/582600/files/wildtype.gff)
> 
>     ```
>     https://zenodo.org/record/582600/files/mutant_R1.fastq
>     https://zenodo.org/record/582600/files/mutant_R2.fastq
>     https://zenodo.org/record/582600/files/wildtype.fna
>     https://zenodo.org/record/582600/files/wildtype.gbk
>     https://zenodo.org/record/582600/files/wildtype.gff
>     ```  
>
>     > ### {% icon tip %} Tip: Importing data via links
>     >
>     > * Copy the link location
>     > * Open the Galaxy Upload Manager
>     > * Select **Paste/Fetch Data**
>     > * Paste the link into the text field
>     > * Press **Start**
>     {: .tip}
>
{: .hands_on}


# Simple Gene Tracks

We will start by adding a single track containing the genes from the wildtype.gff file.

> ### {% icon hands_on %} Hands-on: Build the JBrowse
>
> 1. *Reference Genome to display*: `Use a genome from History`
> 2. *Select the reference genome*: `wildtype.fna`
> 3. *Insert Track Group*
>    1. *Insert Annotation Track*
>         1. *Track Type*: `GFF/GFF3/BED/GBK Features`
>         2. *GFF/GFF3/BED/GBK Data*: `wildtype.gff`
> 4. Execute
> 5. View the contents of the file
>
> > ![Screenshot of JBrowse](../../images/jbrowse-gff-track.png "Screenshot of JBrowse")
{: .hands_on}

If you are not familiar with the operation of JBrowse there are some important points:

- "Tracks" are shown on the left. Clicking the checkboxes will make the tracks visible or invisible
- You can use your mouse scrollwheel to move around the genome view area, or you can click and drag to move.
- Double clicking will zoom in on the genome, or you can use the magnifying glass icons to zoom in our out

> ### {% icon tip %} Tip: Naming Tracks
>
> * The JBrowse tool takes track names directly from file names
> * If you want to rename tracks: **Click** on the **pencil** icon, edit the **Name** and click **Save**.
> * You can now re-run the JBrowse tool and it will produce a new JBrowse instance with corrected names.
{: .tip}

# Complex Gene Tracks

All of the track types in the JBrowse tool support a wide array of features. We've only looked at a simple track with default options, however there are more tools available to us to help create user-friendly JBrowse instances that can embed rich data.

> ### {% icon hands_on %} Hands-on: Build the JBrowse
>
> 1. *Reference Genome to display*: `Use a genome from History`
> 2. *Select the reference genome*: `wildtype.fna`
> 3. *Insert Track Group*
>    1. *Insert Annotation Track*
>         1. *Track Type*: `GFF/GFF3/BED/GBK Features`
>         2. *GFF/GFF3/BED/GBK Data*: `wildtype.gff`
>         3. *Index this Track*: `Yes`
>         4. *JBrowse Styling Options*
>             1. *JBrowse style.label*: `name`
>             2. *JBrowse style.description*: `product`
>         5. *JBrowse Contextual Menu Options*
>             1. *Menu Action*: `iframeDialog`
>             2. *Menu Label*: `{Locus_tag} on NCBI`
>             3. *Menu title*: `NCBI Protein {id}`
>             4. *Menu url*: `https://www.ncbi.nlm.nih.gov/protein/{Locus_tag}`
>             5. *Menu Icon*: `Search`
>         6. *Track Visibility*: `On for new users`
> 4. Execute
> 5. View the contents of the file
>
> > ### {% icon tip %} Tip: Templating in Contextual Menu Options
> >
> > There are more valid values for templating than mentioned in the help. When
> > you click on a feature in the JBrowse instance, it will present all of the
> > properties of the feature. Any of the top level properties can be used
> > directly in your templating
> {: .tip}
{: .hands_on}

Contextual menus can be used to link to more than just NCBI.

- The links can go anywhere such as web search services (e.g. Google) or genomics web services (e.g. EBI)
- Some sites use the IFrame action to link genes to local services where users are expected to submit annotation notes or data.

Unfortunately this object does not have access to the sequence, so this makes
building a link which could initiate a BLAST (or other) sequence search not
possible.


# Blast Results
