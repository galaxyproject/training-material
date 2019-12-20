---
layout: tutorial_hands_on

title: "Visualisation with Circos"
questions:
  - "What can the Circos Galaxy tool be used for?"
  - "How can I visualise common genomic datasets using Circos?"
objectives:
  - "Create a number of Circos plots using the Galaxy tool"
  - "Familiarise yourself with the various different track types"
time_estimation: "2h"
key_points:
  - "Circos is an effective tool to make circular visualisation of high-dimensional datasets"
  - "Circos is often used for genomics, but can also be used for other types of data"
contributors:
  - shiltemann
  - erasche
---

# Introduction
{:.no_toc}


Circos {% cite krzywinski2009circos %} is a software package for visualizing data in a circular layout. This makes Circos ideal for exploring relationships between objects or positions. Circos plots have appeared in thousands of scientific publications. Although originally designed for visualizing genomic data, it can create figures from data in any field.

![Panel of example Circos images](../../images/circos-sample-panel.png)

In this tutorial you will learn how to create such publication-ready Circos plots within Galaxy, and hopefully you can draw inspiration from these for developing your own plots.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Background

Circos supports various different plot types, such as histograms, scatter plots and heat maps. Each Circos plot may contain multiple tracks containing different sub-plots, making it ideal for visualisation of high-dimensional data.


## Circos is an Iterative Process

Publication quality circos plots are *rarely* produced on the first try. Developing a quality Circos plot involves a lot of trial and error to find the best way to convey specific pieces of your data to your audience. Usually you will build up the circos plot one track at a time, and play around with different parameter settings until the plot looks exactly like you want it to:

![Gif of a circos plots through the development process](../../images/circos2.gif "Circos plots are an iterative process, requiring many iterative steps, each improving your plot and getting you closer to a final image.")


Circos is an extremely flexible but also very complex tool. The Galaxy Circos tool covers the most commonly used Circos features, but in order to avoid becoming too complex, it does not expose every single configuration option available in Circos. However, the Galaxy Circos tool allows you to download the full set of configuration files it uses, allowing you to manually tweak the plot further.


> ### {% icon comment %} Comment: Circos tutorials
>
> To learn more about using Circos outside of Galaxy (e.g. for tweaking the
> Circos configuration output by the Galaxy tool), there are a wide range
> of tutorials available from the [Circos website](https://circos.ca)
{: .comment}


# Circos Basics

## Ideogram

## Data Tracks

[track types, radius to choose position of track within plot, input format]

### 2D data track

### Link Track

### Text Labels

# Example 1: Cancer Genomics

In this section, we will recreate a Circos plot of the VCaP cancer cell line presented in {% cite alves2013gene %}. In this study, data from various sources were combined into a single integrative Circos plot.

![VCaP cancer Circos plot](../../images/circos/vcap.png){: width="60%"}


This plot has 4 tracks
 - Structural variants (2 tracks, data obtained from whole-genome NGS sequencing data)
 - B-allele Frequency (obtained from SNP-array data)
 - Copy Number (obtained from SNP array data)


In this section we will reproduce this Circos plot step by step.

## Ideogram

As the first step to this Circos plot, let's configure the ideogram (set of chromosomes to draw).

> ### {% icon hands_on %} Hands-on: Set ideogram configuration
>
> 1. **Circos** {% icon tool %} with the following parameters:
>    - In *"Reference Genome and Cytogenetic Bands"*:
>        - *"Reference Genome"*: `Karyotype`
>            - {% icon param-file %} *"Karyotype Configuration"*: `karyotype_human_hg18.txt`
>    - In *"Ideogram"*:
>        - *"Radius"*: `0.85`
>        - *"Thickness"*: `45`
>        - In *"Labels"*:
>            - *"Label Font Size"*: `64`
>        - In *"Cytogenic Bands"*:
>            - {% icon param-file %} *"Cytogenic Bands"*: `karyotype-bands.txt`
>            - *"Convert bands from BED format to circos karyotype band format"*: `No`
>            - *"Fill Bands"*: `2`
>            - *"Band Stroke Thickness"*: `1`
{: .hands_on}

You should now have a plot that looks like this:

![](../../images/circos/cancer_ideogram.png){: width="50%"}

We will use this as the basis for our plot, and add data tracks one at a time.

## Structural Variations

### Background

Structural variants (SVs) are large-scale genomic rearrangements. SVs involve large segments of DNA (>50 bp) that are deleted, duplicated, translocated or inverted.

![Overview of types of SVs](../../images/circos/sv-types.jpg){: width="50%"}


One of the first observations of SVs in the human genome is known as the [Philadelphia Chromosome](https://en.wikipedia.org/wiki/Philadelphia_chromosome), a SV observed in leukemia. In this mutation, a translocation of genetic material  occurs between chromosomes 9 and 22, resulting in a fusion between genes *BCR* and *ABL1*, causing the production of a hybrid protein, impairing various signalling pathways and causing the cell to divide uncontrollably.

![the Philedelphia chromosome](../../images/circos/sv.jpg){: width="50%"}

In cancer analyses it is therefore often useful to examine SVs and look for potential fusion genes that may affect cell function.

### Plotting SVs

We will now plot these SVs using a link track type, and colouring the links differently depending on whether the SVs are intrachromosomal (within a single chromosome) or interchromosomal (between different chromosomes).

Unfortunately, there is no standard file format for SV data, with different SV callers outputting different formats. Our first step will be to transform our input dataset to the Circos format for link tracks.

**SV File Format:**

```
#ASSEMBLY_ID	GS000008107-ASM
#SOFTWARE_VERSION	2.0.2.22
#GENERATED_BY cgatools
#GENERATED_AT	2012-Mar-07 19:33:32.930656
#FORMAT_VERSION	2.0
#GENOME_REFERENCE	NCBI build 36
#SAMPLE	GS00669-DNA_D02
#TYPE	JUNCTIONS
#DBSNP_BUILD	dbSNP build 130
#GENE_ANNOTATIONS	NCBI build 36.3
>Id	LeftChr	LeftPosition	LeftStrand	LeftLength	RightChr	RightPosition	RightStrand	RightLength	StrandConsistent	Interchromosomal	Distance	DiscordantMatePairAlignments	JunctionSequenceResolved	TransitionSequence	TransitionLength	LeftRepeatClassification	RightRepeatClassification	LeftGenes	RightGenes	XRef	DeletedTransposableElement	KnownUnderrepresentedRepeat	FrequencyInBaselineGenomeSet	AssembledSequence	EventId	Type	RelatedJunctions
3872	chr1	815629	-	155	chr1	5649523	-	822	Y	N	4833894	13	Y		0	AluS:SINE:Alu;SegDup;Self chain	AluSx:SINE:Alu;Self chain						0.92	gcaaccaaaatgactattctttctaccctcCTAGTTCAGACATAGCCTGAGACTTTTTTTTTTTTGAGATGAAGTCTCACTCTGTCACTCAGGCTGGAGTGCAGTGGCATGGTCTCGGCTCATTGCAATCTCTACCTCCCGGGTTCAAGTGATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGGCTACAGGCGTGCACCACCACACCTGGCTAATTTTCATATTATTAGTAGAGATGGGGTTTCACCATGTTGGTCAGACTGGTCTTGAACTCCTGACCTCAGGTGATCTGCCCGCCTCTGCCTCCCAAAATGCTGAGATTACAGATGTGAGCCACTGtgcccggccgcctgagacattttggacgac	490	complex	6269;6270;8575;8577;8578
8577	chr1	816163	+	449	chr1	5650075	+	435	Y	N	4833912	21	Y		0	SegDup;Self chain	Self chain						1.00	gactgcagggggcaggagctctctggctggGCCTTGATCGTGTTCAAGCCACAACCACAGACCTAGGCGTGGTCCCTCAGCCACCTTGTAGCCTTGGCTTGCAACATCTCGACATGGAAACCAAAATGCAGCAGGGCCAATGTGATCTGAAGTTTCCTGAAAAGTTTCTCAGACCCcctcttttaccccttgtgcaacctgcacac	490	complex	3872;6269;6270;8575;8578
8578	chr1	816176	+	363	chr1	5650768	+	228	Y	N	4834592	12	Y	CTTTTGTAACCTGCACACAGTGACCTGTATTCTAGAGGGTCCACACAGAGCTGCCATTCCTTCTGCCAGACCCTGCGGGACTCAGGCATTCTGGAGGCTTCCTGCCCTACAAAGGCAGCCAGACTCCCGCCATGCATCCCTGCACCAGCGGCTCACGGCCAGCTCCCTCACCTGCACCAGCGGCTCACGGCTAGCTCCCTCATCTGCATTCCAGTGGCTCATGGCCAGCTCCCTCACCTGCACCAGCGGCTCTCGGCCGGCTCCCTCCCCTGCACTCCAGCGGC	284	Self chain	Self chain;Tandem period 100;Tandem period 34;Tandem period 66						1.00	aaaccaaaacgcagcagagcccatgtgatcTGAAAGTTCCTGAAAAGTTGCCCAGACCCCCTCTTGTACCCCTTTTGTAACCTGCACACAGTGACCTGTATTCTAGAGGGTCCACACAGAGCTGCCATTCCTTCTGCCAGACCCTGCGGGACTCAGGCATTCTGGAGGCTTCCTGCCCTACAAAGGCAGCCAGACTCCCGCCATGCATCCCTGCACCAGCGGCTCACGGCCAGCTCCCTCACCTGCACCAGCGGCTCACGGCTAGCTCCCTCATCTGCATTCCAGTGGCTCATGGCCAGCTCCCTCACCTGCACCAGCGGCTCTCGGCCGGCTCCCTCCCCTGCACTCCAGCGGCtcaccgccggctccctcacctgcactccag	490	complex	3872;6269;6270;8575;8577
```


**Circos Input Format:**

```
chromosome - start - end - chromosome - start - end
```

So in order to convert this to Circos format, we need to
- Remove header lines (lines starting with `#`)
- Select the columns containing the chromosomes and positions of the breaks (junctions)

> ### {% icon hands_on %} Hands-on: Prepare input data
>
> 1. **Select** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `VCaP highConfidenceJunctions.tsv`
>    - *"that"*: `NOT Matching`
>    - *"the pattern"*: `^[#><]`
>
> 2. **Cut** {% icon tool %} with the following parameters:
>    - *"Cut columns"*: `c2,c3,c3,c6,c7,c7`
>    - {% icon param-file %} *"From"*: `out_file1` (output of **Select** {% icon tool %})
{: .hands_on}

Now that we have the correct format, we can

> ### {% icon hands_on %} Hands-on: Add Circos link track
>
> 1. Hit **Rerun** {% icon galaxy-refresh %} on the previous Circos {% icon tool %} run (where we set up the ideogram)
>
> 2. Add a new Link Track for the SV data, colouring by SV type:
>    - In *"Link Tracks"*:
>       - In *"Link Data"*:
>           - {% icon param-repeat %} *"Insert Link Data"*
>               - *"Inside Radius"*: `0.95`
>               - {% icon param-file %} *"Link Data Source"*: `out_file1` (output of **Cut** {% icon tool %})
>               - *"Link Type"*: `basic`
>               - *"Thickness"*: `3.0`
>               - *"Bezier Radius"*: `0.5`
>               - In *"Rules"*:
>                   - In *"Rule"*:
>                       - {% icon param-repeat %} *"Insert Rule"*
>                           - In *"Conditions to Apply"*:
>                               - {% icon param-repeat %} *"Insert Conditions to Apply"*
>                                   - *"Condition"*: `Interchromosomal`
>                           - In *"Actions to Apply"*:
>                               - {% icon param-repeat %} *"Insert Actions to Apply"*
>                                   - *"Action"*: `Change Link Colour`
>                                   - *"Link Color"*: `red` (Select from the colour wheel)
>
{: .hands_on}


Your output should look something like this:

![The plot with an SV track](../../images/circos/cancer_svs1.png){: width="60%"}

> ### {% icon question %} Questions
>
> 1. Are there more interchromosomal or intrachromosomal SVs?
> 2. Which chromosome appears to have the most SVs?
>
> > ### {% icon solution %} Solution
> >
> > 1. Interchromosomal SVs (between different chromosomes) are coloured red in this
> >    plot, while SVs within a single chromosome are coloured black. You can now easily see that there are more intrachromosomal SVs than interchromosomal.
> > 2. Chromosome 5 appears to have a lot more SVs than the other chromosomes (it looks almost completely black!)
> >
> {: .solution}
{: .question}

We see from this image that chromosome 5 has an unusually large number of SVs, let's look at that chromosome more closely:


> ### {% icon hands_on %} Hands-on: Plot only Chromosome 5
>
> 1. Hit **Rerun** {% icon galaxy-refresh %} on the previous Circos {% icon tool %} run
>
> 2. **Change** the following tool parameters:
>    - In *"Ideogram"*:
>       - *"Limit/Filter Chromosomes"*: `chr5`
>
{: .hands_on}

You should see a plot like:

![Circos plot of chromosome 5 SVs](../../images/circos/cancer_svs_chr5.png){: width="50%"}


> ### {% icon question %} Questions
>
> 1. Are there indeed significantly more SVs on chromosome 5 than on the other chromosomes? (hint: plot some of the other chromosomes as well)
> 2. Are the SVs equally distributed over chromosome 5? Can you think of an explanation for this?
>
> > ### {% icon solution %} Solution
> >
> > 1. Yes, plotting for example only chromosome 1 and comparing this with the chromosome 5 plot, reveals that chr5 has abnormally high number of SVs
> >
> >    ![Circos plot of chromosome 5 SVs](../../images/circos/cancer_svs_chr5.png){: width="30%"}
> >    ![Circos plot of chromosome 1 SVs](../../images/circos/cancer_svs_chr1.png){: width="30%"}
> >
> > 2. No, only part of chromosome 5 appears to be affected. It turns out that this region is exactly one arm of the chromosome.
> >    This could be caused by a phenomenon known as *chromothripsis*
> >
> {: .solution}
{: .question}



## Copy Number Variation

> ### {% icon hands_on %} Hands-on: Prepare the B-allele frequency table
>
> 1. **Select** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `B-allele frequency.tsv`
>    - *"that"*: `NOT Matching`
>    - *"the pattern"*: `^Chromosome`
>
> 1. **Select random lines** {% icon tool %} with the following parameters:
>    - *"Randomly select"*: `25000`
>    - {% icon param-file %} *"from"*: `out_file1` (output of **Select** {% icon tool %})
>    - *"Set a random seed"*: `Don't set seed`
>
{: .hands_on}





> ### {% icon hands_on %} Hands-on: Create copy number track
> 1. **Select** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `VCaP copy number.tsv`
>    - *"that"*: `NOT Matching`
>    - *"the pattern"*: `^Chromosome`
>
> 1. **Cut** {% icon tool %} with the following parameters:
>    - *"Cut columns"*: `c1,c2,c3,c4`
>    - {% icon param-file %} *"From"*: `out_file1` (output of **Select** {% icon tool %})
>
> 1. **Select random lines** {% icon tool %} with the following parameters:
>    - *"Randomly select"*: `25000`
>    - {% icon param-file %} *"from"*: `out_file1` (output of **Cut** {% icon tool %})
>    - *"Set a random seed"*: `Don't set seed`
{: .hands_on}


> ### {% icon hands_on %} Hands-on: Prepare VCF file
> 1. **SnpSift Extract Fields** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Variant input file in VCF format"*: `VCaP ListVariants.vcf`
>    - *"Fields to extract"*: `CHROM POS POS CHROM POS POS`
>
> 1. **Paste** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Paste"*: `output` (output of **SnpSift Extract Fields** {% icon tool %})
>    - {% icon param-file %} *"and"*: `output` (output of **SnpSift Extract Fields** {% icon tool %})
>
> 1. **Circos: Link Density Track** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Links file"*: `out_file1` (output of **Paste** {% icon tool %})
>    - *"Normalize"*: `Yes`
{: .hands_on}


> ### {% icon hands_on %} Hands-on: Plot the data!
> 1. **Circos** {% icon tool %} with the following parameters:
>    - In *"Reference Genome and Cytogenetic Bands"*:
>        - *"Reference Genome"*: `Karyotype`
>            - {% icon param-file %} *"Karyotype Configuration"*: `output` (Input dataset)
>    - In *"Plot Options"*:
>        - *"Plot Format"*: `Color`
>    - In *"Ideogram Configuration (Genome/Chromosomes)"*:
>        - *"Radius"*: `0.85`
>        - *"Thickness"*: `45.0`
>        - In *"Labels"*:
>            - *"Label Font Size"*: `64`
>        - In *"Cytogenic Bands"*:
>            - *"Fill Bands"*: `2`
>            - *"Band Stroke Thickness"*: `1`
>    - In *"Ticks"*:
>        - *"Show Ticks"*: `Yes`
>        - In *"Tick Group"*:
>            - {% icon param-repeat %} *"Insert Tick Group"*
>                - *"Tick Spacing"*: `10u`
>                - *"Show Tick Labels"*: `Yes`
>    - In *"2D Data"*:
>        - In *"2D Data Plot"*:
>            - {% icon param-repeat %} *"Insert 2D Data Plot"*
>                - *"Outside Radius"*: `0.95`
>                - *"Minimum / maximum options"*: `Supply min/max values`
>                    - *"Minimum value"*: `-1.0`
>                    - *"Maximum value"*: `1.0`
>                - *"Plot Format"*: `Scatter`
>                    - In *"Plot Format Specific Options"*:
>                        - *"Glyph Size"*: `4`
>                        - *"Color"*: {% color_picker #7f7f7f %} (gray)
>                        - *"Stroke Thickness"*: `0`
>                - In *"Rules"*:
>                    - In *"Rule"*:
>                        - {% icon param-repeat %} *"Insert Rule"*
>                            - In *"Conditions to Apply"*:
>                                - {% icon param-repeat %} *"Insert Conditions to Apply"*
>                                    - *"Condition"*: `Apply based on point value`
>                                        - *"Points above this value"*: `0.15`
>                            - In *"Actions to Apply"*:
>                                - {% icon param-repeat %} *"Insert Actions to Apply"*
>                                    - *"Action"*: `Change Fill Color for all points`
>                                        - *"Fill Color"*: {% color_picker #00b050 %} (green)
>                        - {% icon param-repeat %} *"Insert Rule"*
>                            - In *"Conditions to Apply"*:
>                                - {% icon param-repeat %} *"Insert Conditions to Apply"*
>                                    - *"Condition"*: `Apply based on point value`
>                                        - *"Points below this value"*: `-0.15`
>                            - In *"Actions to Apply"*:
>                                - {% icon param-repeat %} *"Insert Actions to Apply"*
>                                    - *"Action"*: `Change Fill Color for all points`
>                                        - *"Fill Color"*: {% color_picker #ff0000 %} (red)
>                - In *"Axes"*:
>                    - In *"Axis"*:
>                        - {% icon param-repeat %} *"Insert Axis"*
>                            - *"Inside Radius (y0)"*: `-1.0`
>                            - *"Radial-relative values"*: `Yes`
>                            - *"Spacing"*: `0.5`
>            - {% icon param-repeat %} *"Insert 2D Data Plot"*
>                - *"Outside Radius"*: `0.75`
>                - *"Inside Radius"*: `0.6`
>                - *"Minimum / maximum options"*: `Supply min/max values`
>                    - *"Minimum value"*: `0.0`
>                    - *"Maximum value"*: `1.0`
>                - *"Plot Format"*: `Scatter`
>                    - {% icon param-file %} *"Scatter Plot Data Source"*: `out_file1` (output of **Select random lines** {% icon tool %})
>                    - In *"Plot Format Specific Options"*:
>                        - *"Glyph Size"*: `4`
>                        - *"Color"*: {% color_picker #7f7f7f %} (gray)
>                        - *"Stroke Thickness"*: `0`
>                - In *"Axes"*:
>                    - In *"Axis"*:
>                        - {% icon param-repeat %} *"Insert Axis"*
>                            - *"Radial-relative values"*: `Yes`
>                            - *"Spacing"*: `0.25`
>            - {% icon param-repeat %} *"Insert 2D Data Plot"*
>                - *"Outside Radius"*: `0.55`
>                - *"Inside Radius"*: `0.4`
>                - *"Minimum / maximum options"*: `Plot all values`
>                - *"Plot Format"*: `Histogram`
>                    - In *"Plot Format Specific Options"*:
>                        - *"Transparency"*: `0.0`
>                        - *"Stroke Color"*: {% color_picker #00b0f0 %} (sky blue)
>                        - *"Fill underneath the histogram"*: `Yes`
>                - In *"Axes"*:
>                    - In *"Axis"*:
>                        - {% icon param-repeat %} *"Insert Axis"*
>                            - *"Radial-relative values"*: `Yes`
>                            - *"Spacing"*: `0.25`
>                - In *"Backgrounds"*:
>                    - In *"Background"*:
>                        - {% icon param-repeat %} *"Insert Background"*
>                            - *"Radial-relative values"*: `Yes`
>                            - *"Color"*: {% color_picker #f2f2f2 %} (light grey)
>    - In *"Link Tracks"*:
>        - In *"Link Data"*:
>            - {% icon param-repeat %} *"Insert Link Data"*
>                - *"Inside Radius"*: `0.38`
>                - {% icon param-file %} *"Link Data Source"*: `out_file1` (output of **Cut** {% icon tool %})
>                - *"Link Type"*: `basic`
>                - *"Thickness"*: `1.0`
>                - *"Bezier Radius"*: `0.25`
>                - In *"Advanced Settings"*:
>                    - *"Bezier Radius Purity"*: `1.0`
>                    - *"Perturb links?"*: `no`
>                - In *"Rules"*:
>                    - In *"Rule"*:
>                        - {% icon param-repeat %} *"Insert Rule"*
>                            - In *"Conditions to Apply"*:
>                                - {% icon param-repeat %} *"Insert Conditions to Apply"*
>                                    - *"Condition"*: `Interchromosomal`
>                            - In *"Actions to Apply"*:
>                                - {% icon param-repeat %} *"Insert Actions to Apply"*
>                                    - *"Action"*: `Change Visibility`
>            - {% icon param-repeat %} *"Insert Link Data"*
>                - *"Inside Radius"*: `0.3`
>                - *"Link Type"*: `basic`
>                - *"Link Color"*: {% color_picker #ff0000 %} (red)
>                - *"Thickness"*: `2.0`
>                - *"Bezier Radius"*: `0.0`
>                - In *"Advanced Settings"*:
>                    - *"Perturb links?"*: `no`
>                - In *"Rules"*:
>                    - In *"Rule"*:
>                        - {% icon param-repeat %} *"Insert Rule"*
>                            - In *"Conditions to Apply"*:
>                                - {% icon param-repeat %} *"Insert Conditions to Apply"*
>                                    - *"Condition"*: `Intrachromosomal`
>                            - In *"Actions to Apply"*:
>                                - {% icon param-repeat %} *"Insert Actions to Apply"*
>                                    - *"Action"*: `Change Visibility`
>
{: .hands_on}

# Nature Cover ENCODE Diagram

Here we will reproduce the output of the [Circos tutorial](http://circos.ca/documentation/tutorials/recipes/nature_cover_encode/lesson) for producing an image like that which was used on Nature's Cover:

![Nature Cover ENCODE](../../images/nature_cover_encode.png "Nature Cover")

The official Circos tutorial goes into detail on how to use rules and variables and automagic counting to help automate the production of such an image.

```
<plot>
r1   = dims(ideogram,radius_inner)
         - conf(plot_padding)*eval(remap(counter(plot),0,conf(num_plots),1,0.9))
         - eval((conf(plot_width)+conf(plot_padding))*counter(plot)*eval(remap(counter(plot),0,conf(num_plots),1,0.9)))
r0   = conf(.,r1)
         - conf(plot_width)*eval(remap(counter(plot),0,conf(num_plots),1,0.9))
post_increment_counter = plot:1
<<include rules.conf>>
</plot>
```

Due to the Circos Galaxy tool being a simplified interface and more user-friendly interface, we cannot make use of these advanced features like using Perl's math or `eval`uating code snippets to generate our image. <!-- TODO(hxr): fix the wording, this is gross -->

## Data Formats

<!-- TODO(hxr): document it all -->

The Circos Galaxy tool mostly accepts `tabular` files. These always have at least three columns:

1. Chromosome
2. Start
3. End

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    - [Chromosome File](../../files/chrom.tab)
>    - [Highlights](../../files/highlights.tab)
>
>    {% include snippets/import_via_link.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% include snippets/change_datatype.md datatype="tabular" %}
>
{: .hands_on}

- what is required for a plot (karyotype)
- outputs
- plot options
- ideogram
- ticks
- highlights
- 2d
- links


> ### {% icon hands_on %} Hands-on: Circos
>
> We will now run the Circos tool. In terms of tool interface, it is one of the most complex extant Galaxy tools.
>
> > ### {% icon tip %} Tip: Interface Complexity
> > The interface looks deceptively simple when all of the sections are collapsed, but as you start adding tracks it can be easy to get lost and become overwhelmed, so just go slowly. Do not worry if your plot does not look exactly like the expected output.
> {: .tip}
>
> 1. **Circos** {% icon tool %} with the following parameters:
>    - In *"Reference Genome and Cytogenetic Bands"*:
>        - *"Reference Genome"*: `Karyotype`
>            - {% icon param-file %} *"Karyotype Configuration"*: `chrom.tab`
>    - In *"Plot Options"*:
>        - *"Plot Format"*: `Color`
>            - *"Background Color"*: {% color_picker #000000 %}
>    - In *"Ideogram Configuration (Genome/Chromosomes)"*:
>        - *"Thickness"*: `0.0`
>        - In *"Labels"*:
>            - *"Show Label"*: `Yes`
>    - In *"Highlights"*:
>        - In *"Highlight"*:
>            - Click on *"Insert Highlight"*:
>            - In *"1: Highlight"*:
>                - *"Outside Radius"*: `0.99`
>                - *"Inside Radius"*: `0.9`
>                - {% icon param-file %} *"Highlight Data Source"*: `highlights.tab`
>            - Click on *"Insert Highlight"*:
>            - In *"2: Highlight"*:
>                - *"Outside Radius"*: `0.89`
>                - *"Inside Radius"*: `0.8`
>                - {% icon param-file %} *"Highlight Data Source"*: `highlights.tab`
>                - In *"Rules"*:
>                    - In *"Rule"*:
>                        - Click on *"Insert Rule"*:
>                        - In *"1: Rule"*:
>                            - In *"Conditions to Apply"*:
>                                - Click on *"Insert Conditions to Apply"*:
>                                - In *"1: Conditions to Apply"*:
>                                    - *"Condition"*: `Randomly`
>                                        - *"Percentage of bins"*: `0.1`
>                            - In *"Actions to Apply"*:
>                                - Click on *"Insert Actions to Apply"*:
>                                - In *"1: Actions to Apply"*:
>                                    - *"Action"*: `Change Fill Color for all points`
>                                        - *"Fill Color"*: {% color_picker #8064a2 %} (light purple)
>                            - *"Continue flow"*: `Yes`
>                        - Click on *"Insert Rule"*:
>                        - In *"2: Rule"*:
>                            - In *"Conditions to Apply"*:
>                                - Click on *"Insert Conditions to Apply"*:
>                                - In *"1: Conditions to Apply"*:
>                                    - *"Condition"*: `Randomly`
>                                        - *"Percentage of bins"*: `0.1`
>                            - In *"Actions to Apply"*:
>                                - Click on *"Insert Actions to Apply"*:
>                                - In *"1: Actions to Apply"*:
>                                    - *"Action"*: `Change Fill Color for all points`
>                                        - *"Fill Color"*: {% color_picker #ffff00 %} (yellow)
>                            - *"Continue flow"*: `Yes`
>            - Click on *"Insert Highlight"*:
>            - In *"3: Highlight"*:
>                - *"Outside Radius"*: `0.79`
>                - *"Inside Radius"*: `0.7`
>                - {% icon param-file %} *"Highlight Data Source"*: `highlights.tab`
>                - In *"Rules"*:
>                    - In *"Rule"*:
>                        - Click on *"Insert Rule"*:
>                        - In *"1: Rule"*:
>                            - In *"Conditions to Apply"*:
>                                - Click on *"Insert Conditions to Apply"*:
>                                - In *"1: Conditions to Apply"*:
>                                    - *"Condition"*: `Randomly`
>                                        - *"Percentage of bins"*: `0.2`
>                            - In *"Actions to Apply"*:
>                                - Click on *"Insert Actions to Apply"*:
>                                - In *"1: Actions to Apply"*:
>                                    - *"Action"*: `Change Fill Color for all points`
>                                        - *"Fill Color"*: {% color_picker #8064a2 %} (light purple)
>                            - *"Continue flow"*: `Yes`
>                        - Click on *"Insert Rule"*:
>                        - In *"2: Rule"*:
>                            - In *"Conditions to Apply"*:
>                                - Click on *"Insert Conditions to Apply"*:
>                                - In *"1: Conditions to Apply"*:
>                                    - *"Condition"*: `Randomly`
>                                        - *"Percentage of bins"*: `0.2`
>                            - In *"Actions to Apply"*:
>                                - Click on *"Insert Actions to Apply"*:
>                                - In *"1: Actions to Apply"*:
>                                    - *"Action"*: `Change Fill Color for all points`
>                                        - *"Fill Color"*: {% color_picker #ffff00 %} (yellow)
>
>
> 2. View the output PNG file
>
{: .hands_on}

When this has complete, your output should look similar to the following;

![Circos simplified Nature ENCODE Cover](../../images/circos-encode.png "Simplified Nature Cover")

<!--

What question could we ask about this plot?
***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}


-->

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
