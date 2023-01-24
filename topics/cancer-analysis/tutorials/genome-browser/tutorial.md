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
>     TODO: tumor.bam, normal.bam
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
{: .hands_on}


These files are *alignment* files. This means the raw sequence reads have been mapped to the human reference genome. To learn more about this process of mapping, see our [dedicated tutorial]({% link topics/sequence-analysis/tutorials/mapping/tutorial.md %})


To view alignment files in a genome browser, you usually provide the BAM file, along with an *index* file, which allows the tools to quickly navigate the usually very large alignment files.
However, in our case, Galaxy automatically creates the index file whenever you upload a BAM file, and will supply it to the tools that need it behind the scenes.

TODO: note about test data we are using, scaled down for tutorial reasons



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

> <hands-on-title> Navigating around the genome </hands-on-title>
>
> 1. Check the box next to `normal.bam` {% icon param-check %} on the left-hand side
>
> 2. Navigate to `chr9:130,620,912-130,621,487` by copying the location into the location bar and hitting <kbd>Enter</kbd>
>
{: .hands_on}




## Adding reference tracks

genes

dbsnp


# Inspecting small variants in the normal sample

# Inspecting small somatic variants in the tumor sample

# Inspecting structural variants in NA12878

