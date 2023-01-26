---
layout: tutorial_hands_on

title: "Viewing Cancer Alignments in a Genome Browser"
questions:
objectives:
  - Very quickly navigate around the genome
  - Visualise HTS read alignments
  - Examine SNP calls and structural re-arrangements by eye
  - Tell the difference between germline and somatic variants
key_points:
time_estimation: 2h
contributions:
  authorship: [soranamorrissey, rdeborjas]
  editing: [shiltemann]
  funding: [bioinformatics-ca,erasmusplus]


---

<!-- TODO add contributors. From bioinf.ca: "This lab is based on the HTS IGV lab originally by
Sorana Morrissy and was updated and modified by Heather Gibling for the Cancer Analysis workshop. " -->

{% snippet topics/cancer-analysis/faqs/attribution.md %}

# Introduction

This tutorial will introduce you to the Genome browsers; powerful tools for viewing many kinds of
genomic data, including data for DNA sequencing, RNA sequencing, microarrays, epigenetics, and
copy number alteration. For this tutorial we will use JBrowse, since it has a nice integration
with Galaxy, but many other Genome browsers exist, with [IGV](https://igv.org/app/) perhaps being the most widely used. TODO: link some examples/comparisons

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Data upload

Before we view our data in the Genome browser, let's upload it to Galaxy

> <hands-on-title> Upload data </hands-on-title>
>
> 1. Make sure you have an empty analysis history. Give it a name.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import data from Zenodo
>
>    ```
>     TODO: tumor.bam, normal.bam, genes.bed, dpsnp.bed
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
{: .hands_on}


These files are *alignment* files. This means the raw sequence reads have been mapped to the human reference genome. To learn more about this process of mapping, see our [dedicated tutorial]({% link topics/sequence-analysis/tutorials/mapping/tutorial.md %})


To view alignment files in a genome browser, you usually provide the BAM file, along with an *index* file, which allows the tools to quickly navigate the usually very large alignment files.
However, in our case, Galaxy automatically creates the index file whenever you upload a BAM file, and will supply it to the tools that need it behind the scenes.

TODO: note about test data we are using, scaled down for tutorial reasons, fact that it's paired-end data



# Getting familiar with JBrowse

Let's start by launching JBrowse and loading our data so that we can get a feel for JBrowe and its interface. JBrowse works just like any other tool in Galaxy, you provide your data, configure some settings, en execute. The output will be a copy of the JBrowse genome browser with your data loaded.

> <hands-on-title> Build JBrowse </hands-on-title>
>
> 1. {% tool [JBrowse genome browser](toolshed.g2.bx.psu.edu/repos/iuc/jbrowse/jbrowse/1.16.11+galaxy1) %} with the following parameters:
>    - *"Select a reference genome"*: `hg19`
>    - {% icon param-repeat %} *"Insert Track Group"*
>      - {% icon param-repeat %} *"Insert Annotation Track"*
>      - {% icon param-select %} *"Track Type"*: `BAM Pileups`
>      - {% icon param-files %} *"BAM Track Data"* (hold CTRL to select multiple)
>        - `normal.bam`
>        - `tumor.bam`
>    - {% icon param-toggle %} Autogenerate SNP Track: `Yes`
>
{: .hands_on}


## The JBrowse interface

Let's start by opening JBrowse, just to get a feel for what it looks like.

> <hands-on-title> A first look at JBrowse </hands-on-title>
>
> 1. Click on the eye icon {% icon galaxy-eye %} to view the data in JBrowse
>
>    ![A screenshot of the JBrowse interface, described by caption](./images/jbrowse-screenshot.png "The JBrowse Interface. At the top is the navigation bar; here you can zoom in and out, and provide a location on the gennome to jump to. On the left is the list of available data tracks, these can be checked and unchecked to show and hide them from view. In the main panel the data is shown. By default this is the reference data track and the GCContent track. This panel can be dragged left and right to move along the genome.")
>
> 2. If you need a bit more space on your screen, consider resizing Galaxy's side panels
>
>    {% snippet faqs/galaxy/interface_side_panels_resize.md %}
>
> 3. Play around in JBrowse
>
>    > <question-title> What do you see? </question-title>
>    >
>    > 1. Which part of the genome are you currently viewing?
>    > 3. Which data tracks do you see in the main panel?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > 1. Check the location bar at the top of the screen. You will see something like `chr1:10342..10677 (366b)`, This means you are currently viewing chromosome 1, from base 10342 to base 10677.
>    > >    By default, you will usually be viewing a part of chromosome 1, unless you specified a different region when starting JBrowse.
>    > > 2. The name of each data track is shown in the upper left corner of the track. By default, you should see two tracks labelled "Reference Sequence" and "GCContentXY". Tracks can be shown and hidden using the checkboxes in the left-hand panel.
>    > {: .solution}
>    {: .question}
>
>
> 4. Click on the **"Zoom in"** {% icon zoom-in %} button until you can see the individual bases in the reference genome track.
>    ![](./images/refbases.png)
>
>    > <question-title> What do you see? </question-title>
>    >
>    > 1. What are all the lettered rows you see?
>    > 2. What do the colours mean?
>    >
>    > > <solution-title></solution-title>
>    > > 1. The two center rows show the nucleotides of the forward and reverse strands. The three top rows show the potential amino acids, at different phasings. The bottom 3 rows show the same for the reverse strand.
>    > >
>    > > 2. Each nucleotide has a unique color. If you zoom out so far that you can not see the
>    > >  individual nucleotide letters anymore, you can still tell their identity by the color.
>    > >   Similarly, on the amino acid rows, stop codons are colored red and indicated with an asterisk (*). Start codons are coloured green.
>    > {: .solution}
>    {: .question}
>
> 5. What is the second track that is visible?
>
>    > <question-title> What do you see? </question-title>
>    >
>    > 1. What is the second track that is visible here?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > This track contains the GC content for the reference genome
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}





## Navigation

There are various ways to navigate around the genome, you can:
  - Zoom in or out of the current region
  - Drag the main panel left and right to move to adjacent streches of the genome
  - Jump to specific regions by entering a location in

Later in this tutorial we will add a genes track, so that you may also enter gene names in the location bar to jump directly to your gene of interest.

Let's get a feel for this


> <hands-on-title> Navigating around the genome </hands-on-title>
>
> 1. If Jbrowse is no longer open, click on the eye icon {% icon galaxy-eye %} to open JBrowse again
>
> 2. **Navigate** around the genome by dragging the main panel left and right to move to adjacent regions. You can also zoom out if you want to view bigger regions at a time.
>
> 3. You can also specify a specific region to jump to. Paste `chr1:10200-10800` into the location bar and hit the "Go" button (or press <kbd>Enter</kbd>)
>
>
>    > <question-title> What do you see? </question-title>
>    >
>    > 1. What do you notice about the reference genome here?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > Because of the colour coding that JBrowse applies, it is easy to see the repeats that are present at the beginning of this region
>    > >
>    > > ![](./images/ref-repeats.png)
>    > {: .solution}
>    {: .question}
>
{: .hands_on}




## Visualising read alignments

Now that we have a feel for JBrowse, let's view some of our data!

> <hands-on-title> Viewing Read Alignments </hands-on-title>
>
> 1. Check the box next to `normal.bam` {% icon param-check %} on the left-hand side
>
> 2. Navigate to `chr9:130,620,912-130,621,487` by copying the location into the location bar and hitting <kbd>Enter</kbd>
>
>    ![Screenshot of the JBrowse interface with the data track loaded](./images/jbrowse-data-track.png)
>
>    > <question-title> What do you see? </question-title>
>    >
>    > 1. What do you see in the data track?
>    > 2. What do the red and blue colours mean? (Hint: you can click on the reads to get more information about the read)
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > 1. The `normal.bam` tracks shows the reads at the location they mapped to. Each horizontal bar signifies one read.
>    > > 2. The red reads mapped to the `+` strand, the blue reads to the `-` strand. By clicking on a read you can get more information,
>    > >    including sequence, strand, mapping score, and more.
>    > >
>    > {: .solution}
>    {: .question}
>
>    > <details-title>Read colours</details-title>
>    >
>    > You may have seen some reads not coloured red or blue, these all have special meanings:
>    >
>    > | colour | meaning|
>    > |--------|--------|
>    > | standard red | forward strand |
>    > | standard blue | reverse strand |
>    > | hard red | forward strand, missing mate |
>    > | hard blue | reverse strand, missing mate |
>    > | light red | forward strand, not proper |
>    > | light blue | reverse strand, not proper |
>    > | black | forward strand, different chr |
>    > | gray | reverse strand, different chr |
>    >
>    > More information about alignment tracks can be found in the [JBrowse documentation](https://jbrowse.org/docs/alignments.html)
>    {: .details}
>
> 3. **Zoom in** {% icon zoom-in %} a bit until you can see individual mismatches on the reads
>
>    ![](./images/snvs.png)
>
>    > <question-title> What do you see? </question-title>
>    >
>    > 1. What do these mismatches mean?
>    > 2. Can you find a Single Nucleotide Variant (SNV)? Is it homozygous or heterozygous?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > 1. Any time a nucleotide in the read does not match the reference base it is mapped to, it is signified by a coloured block (colour depends on the nucleotide).
>    > >    This could signify a read error, a variant, or a mistake in the reference genome.
>    > > 2. Whenever a mismatch occurs in all reads overlapping the position, it is most likely a variant. If it occurs in roughly 50% of reads, it is a heterozygous SNV, if it occurs in roughly 100%
>    > >    it is a homosyzous SNV. In the screenshot above, the red line of `T` nucleotides represents an homozygous SNV at that location. The other mismatches are likely sequencing artifacts
>    > >    because they only occur in a single read.
>    > >
>    > {: .solution}
>    {: .question}
>
> 4. **Click on a read** to get more detailed information
>
>    ![](./images/read-info.png)
>
>    > <question-title> What do you see? </question-title>
>    >
>    > 1. How long is the read?
>    > 2. What can you say about the base quality?
>    > 3. How about the mapping quality?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > 1. This information can be found in the `Length` field. For the example in the screenshot the length is 101 base pairs
>    > > 2. The `Sequence and Quality` field shows the basecalling quality for each nucleotide in the read.
>    > > 3. The `CIGAR` string says something about how well the read mapped. In the screenshot this is `101M`, meaning "101 (mis)matches", so this read mapped without any insertions or deletions. And because we can tell from the colouring of the read there were no mismatches, we know these were all matches, so this read aligned perfectly.
>    > >
>    > >    {% snippet topics/sequence-analysis/faqs/cigar.md %}
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}


Depending on what we might like to do (e.g. quality control, SNP calling, CNV finding), or the type of data we have (RNASeq, paird/unpaired, etc) we might like
to customize how JBrowse displays our data. There are many options to control how alignment data is shown to us. Let's play around with some of those settings:

> <hands-on-title>Customizing Alignment display</hands-on-title>
>
> 1. You can customize the way alignments are displayed in JBrowse
>    - Hover over the `normal.bam` track name in the top-left of the track
>    - Click on the dropdown button to the right of the track name to access the track settings
>
>    ![screenshot of the track settings menu accessed by clicking on the track label](./images/track-settings.png){: width="50%"}
>
>
> 2. **Experiment** with the various settings. Think about which would be best for specific tasks (e.g. quality control, SNP calling, CNV finding).
>
>    - You can change how many reads are visible under `Display mode`
>      - try out the `normal`, `collapsed` and `compact` modes
>
>      ![](./images/jbrowse-display-mode.png){: width="75%"}
>
>    - Try the different settings for `Track visualisation types`
>
>      ![](./images/track-visualisation-types.png){: width="75%"}
>
>    - Try different settings for  `Colouring options`
>
>      ![](./images/jbrowse-colouring-options.png){: width="75%"}
>
>
{: .hands_on}




## Adding reference tracks

When we are viewing our data in a genome browser, it is often useful to include knowledge about the area we are viewing, such as genes or other features overlapping the position, or known variants or polymorphisms. To do this, we can load additional tracks into JBrowse that contain this data.


> <hands-on-title>Customizing Alignment display</hands-on-title>
>
> 1. Hit the **Re-run** {% icon galaxy-refresh %} button on the previous JBrowse output
>    - Scroll down past where we configured the BAM/Pileup tracks
>    - {% icon param-repeat %} Insert Annotation Track
>      - *"Track Type"*: `GFF/GFF3/BED Features`
>      - *"GFF/GFF3/BED Track Data"*: Select the following files (hold <kbd>CTRL</kbd> to select multiple files):
>           - `genes.bed`
>           - `dpsnp.bed`
>
> 2. View {% icon galaxy-eye %} the new JBrowse instance
>    - Enable the `normal.bam`, `genes.bed` and `dpSNP.bed` tracks
>    - Navigate to `chr9:130630233` <!-- chr9:130610720-130610929 -->
>
>      ![](./images/reftracks.png){: width="75%"}
>
>      > <question-title> What do you see? </question-title>
>      >
>      > 1. Is this a SNP or an SNV? Is it homozygous or heterozygous?
>      > 2. Does it impact a gene?
>      >
>      > > <solution-title></solution-title>
>      > >
>      > > 1. It is present in dbSNP (`rs4226`), so it is a common polymorphism in the human popluation. It is present in 100% of our reads, so it is homozygous.
>      > > 2. In the genes track you can see the `AK1` gene overlaps this region
>      > >
>      > {: .solution}
>      {: .question}
>
>
{: .hands-on}



# Inspecting small variants in the normal sample

## Heterozygous and Homozygous SNPs

> <hands-on-title></hands-on-title>
>
> 1. Navigate to `chr9:130607258-130607619`
>    - You should see two potential SNPs
>
>    ![](./images/two-snps.png)
>
> 2. Enable the track {% icon param-check %} `normal.bam - SNPs/Coverage` in the left-hand panel. Also enable the `dbSNP` track.
>    - make sure you can see both the SNPs/Coverage track as well as the alignment track. If your screen is too small to see both at once, set the alignment track (`normal.bam`) display to `compact`
>
>    ![](./images/two-snps-with-coverage.png)
>
>    > <question-title></question-title>
>    >
>    > 1. Which variant is heterozygous and which is homozygous?
>    > 2. What are the variant allele frequencies for each SNP? To find out, hover over the SNP in the `SNPs/Coverage` track.
>    >
>    > > <solution-title></solution-title>
>    > > 1. Left is heterozygous (~half the reads are marked as mismatches at this location) and right is homoszygous (all are marked as mismatches)
>    > > 2. 46% (left) and 100% (right)
>    > {: .solution}
>    {: .question}
>
> 2. Make sure the reads are coloured by read strand (default colouring option). So that red reads are in the forward orientation, and blue reads are in the reverse orientation.
>
>    > <question-title></question-title>
>    > Do these look like true SNPs? What evidence is there for this?
>    >
>    > > <solution-title></solution-title>
>    > > Yes, the allele frequencies are close to what we would expect to see for heterozygous and homozygous SNPS, and all of the mismatched bases are of high quality, and there is no strand bias. Additionally, they both line up with known SNPs in the dbSNP track.
>    > {: .solution}
>    {: .question}
>
> 2. Look at the other mismatched bases in the region between the two SNPs.
>
>    > <question-title></question-title>
>    >  1. Are these sequencing errors, SNPs, or SNVs?
>    >  2. Can a normal sample have somatice SNVs?
>    >
>    > > <solution-title></solution-title>
>    > > 1. Sequencing errors and SNVs. The low quality mismatches are most likely sequencing errors. The higher quality mismatches that only occur in one or two reads are most likely somatic SNVs. They are not SNPs because they are not known germline mutations
>    > > 2. Yes! There is always a chance there will be a mutation or error during genome replication.
>    > {: .solution}
>    {: .question}
>
{: .hands_on}




## Homozygous deletion

## GC coverage




# Inspecting small somatic variants in the tumor sample

## Somatic SNV

## Somatic SNP with change in heterozygosity

## Somatic indel next to SNP with change in heterozygosity


# Inspecting structural variants in NA12878

