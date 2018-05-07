---
layout: tutorial_hands_on
topic_name: visualisation
tutorial_name: jbrowse
---

# Genomic Data Visualisation With JBrowse

> JBrowse is a fast, embeddable genome browser built completely with JavaScript
> and HTML5, with optional run-once data formatting tools written in Perl.

The Galaxy tool accepts data in many formats:

- Intervals/Feature Tracks (GFF/GFF3, BED, GenBank)
- BAM Pileups
- Blast XML results
- Wig/BigWig
- VCF files

and executes the "run-once data formatting tools" mentioned in the description. The JBrowse tool has an incredibly extensive number of options, more than anyone needs most of the time. We'll go through them in detail but feel free to skip the sections that don't apply to the data types you use. Not everyone has Blast XML results to visualise.

## General Options

### Reference Genome to Display

JBrowse features a standard reference genome / genome from history selector used by most other tools. [TODO]

Sequence files with multiple sequences can be provided, JBrowse will allow you to select between the sequences during viewing. Up to 30 will be shown from the dropdown selector within JBrowse (this is a [known issue](https://github.com/GMOD/jbrowse/issues/349).)

![An image showing the reference sequenc eselector in JBrowse](../../images/jbrowse-genome-selector.png)

### Produce Standalone Instance

The JBrowse web app loads data from a "data directory," a directory of files processed in the way JBrowse expects them to be. JBrowse natively displays some data types (bam, tabix, support for others will come in the future) while it requires converting other formats to JSON that it can load (fasta, gff3, etc.)

In orer to have a complete, viewable JBrowse instance we need to have an entire copy of the JBrowse source code (~12 Mb.) This is a lot of extra data to duplicate every time we want to build a JBrowse instance, so sometimes it is useful to not include the JBrowse source code, and instead only include the processed data. When the JBrowse tool is used in concert with the Galaxy-Apollo integration, the data directory is extremely useful as you won't be visualising the results in JBrowse, but instead sending them to Apollo which only needs the processed data.

### Genetic Code

Multiple genetic codes are supported in JBrowse, if your organism uses something other than the standard genetic code, then setting this correctly will affect display of start and stop codons in the genomic sequence.

![Reference sequence display in JBrowse showing start and stop codons highlighted.](../../images/jbrowse-genome-sequence.png)

### JBrowse-in-Galaxy Action

JBrowse allows you to update previously built instances. This is a less common operation for normal JBrowse usage (summarisation of an analysis workflow,) but much more common when JBrowse is used with the Galaxy-Apollo integration. In Apollo users can interactively make annotation based on the JBrowse evidence tracks, before pulling these new features back to Galaxy for re-analysis.

## Track Groups and Tracks

Tracks form the basis of data in JBrowse. Tracks can be organised into nested groups which can help users mentally organise the large amounts of data that can be displayed. For genomic annotation this often took the form of splitting data into Structural and Functional evidence tracks, and then below each of these categories splitting out the data further, e.g. in "Structural / Naive ORF Calls" and "Structural / Model based ORF Calls."

![Groups of tracks in JBrowse, each group has a header with a name and contains one or more tracks.](../../images/jbrowse-sections.png)

The `/` character allows nesting of categories. In newer versions of JBrowse (0.7.0+) you can also use `#date#` to insert the current date in the track group name. For analyses that might change over time (e.g. blasting against an always-updating NR database), it can be very helpful to organise these runs by date.

### GFF/GFF3/BED/GBK Features

These are your standard feature tracks. They usually show genes, mRNAs and other features of interest along a genomic region. This input allows selecting multiple tracks, doing so will cause them to all be styled identically. There is no separate option to rename the tracks, the Galaxy datasets will be brought into the JBrowse instance with the same name that they had in Galaxy.

#### `match`/`match_part` data

For hierarchical GFF/GFF3 data (i.e. with parent, child relationships,) the data may be match/match part data. If your features have a top level feature like `match` with a child `match_part` feature, then this option applies to those datasets. Sometimes those parent matches are under a different sequence ontology term, so you will need to supply the correct parent term there. Implementing sequence ontology traversal which would remove the need to specify this is a non-trivial feature and has not been implemented.

#### Track indexing

JBrowse provides a search functionality which will search over all indexed features. The indexing covers feature just the `Name` attributes. This feature is made optional due to bugs in the past with the indexing feature on JBrowse's side.

#### Track Type

JBrowse supports two rendering modes for tracks, HTML features and Canvas features. From an end user perspective they are mostly identical.



