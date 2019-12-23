---
layout: tutorial_hands_on

title: Genomic Data Visualisation with JBrowse
zenodo_link: https://doi.org/10.5281/zenodo.3591856
questions:
objectives:
time_estimation: "1h"
key_points:
contributors:
  - erasche
  - shiltemann
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

and executes the "run-once data formatting tools" mentioned in its description. The JBrowse tool has an incredibly extensive number of options, more than anyone needs most of the time. We'll go through them in detail but feel free to skip the sections that don't apply to the data types you use. Not everyone has Blast results to visualise.

This tutorial covers version 1.16.5+ of the JBrowse tool, earlier versions will have different behaviour and tool layout.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data Upload

> ### {% icon hands_on %} Hands-on: Getting the data
>
> 1. Create and name a new history for this tutorial.
>
>    {% include snippets/create_new_history.md %}
>
> 2. Import the datasets we will visualize:
>
>    ```
>    https://zenodo.org/record/3591856/files/blastp%20genes.gff3
>    https://zenodo.org/record/3591856/files/blastp%20vs%20swissprot.xml
>    https://zenodo.org/record/3591856/files/dna%20sequencing.bam
>    https://zenodo.org/record/3591856/files/dna%20sequencing%20coverage.bw
>    https://zenodo.org/record/3591856/files/genes%20(de%20novo).gff3
>    https://zenodo.org/record/3591856/files/genes%20(NCBI).gff3
>    https://zenodo.org/record/3591856/files/genome.fa
>    https://zenodo.org/record/3591856/files/RNA-Seq%20coverage%201.bw
>    https://zenodo.org/record/3591856/files/RNA-Seq%20coverage%202.bw
>    https://zenodo.org/record/3591856/files/variants.vcf
>    ```
>
>    {% include snippets/import_via_link.md %}
>
{: .hands_on}

The data for today is a subset of real datasets from *E. coli MG1655 strain K-12*

# Simple Gene Tracks

We will start by adding a couple of gene call tracks. In our case these are genes and gene predictions, but they don't have to be. In the general case this can be any interesting region from an analysis, where a tool has pointed out some region for further inspection for some reason, then this data can be visualised with the "gene" track type.

> ### {% icon hands_on %} Hands-on: Build the JBrowse
>
> 1. **JBrowse** {% icon tool %} with the following parameters:
>    - *"Reference genome to display"*: `Use a genome from history`
>        - {% icon param-file %} *"Select the reference genome"*: `genome.fa`
>    - *"Genetic Code"*: `11. The Bacterial, Archaeal and Plant Plastid Code`
>    - In *"Track Group"*:
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - *"Track Category"*: `Genes`
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `GFF/GFF3/BED Features`
>                        - {% icon param-file %} *"GFF/GFF3/BED Track Data"*: `genes (de novo).gff3`
> 4. Execute
>
> 5. View the contents of the file
>
>    ![Screenshot of JBrowse](../../images/jbrowse-gff-track.png "Screenshot of JBrowse")
>
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

> ### {% icon hands_on %} Hands-on: Build a JBrowse for viewing Genes
>
> 1. **JBrowse** {% icon tool %} with the following parameters:
>    - *"Reference genome to display"*: `Use a genome from history`
>        - {% icon param-file %} *"Select the reference genome"*: `genome.fa`
>    - *"Genetic Code"*: `11. The Bacterial, Archaeal and Plant Plastid Code`
>    - In *"Track Group"*:
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - *"Track Category"*: `Genes`
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `GFF/GFF3/BED Features`
>                        - {% icon param-file %} *"GFF/GFF3/BED Track Data"*: `genes (de novo).gff3`
>                        - *"JBrowse Track Type [Advanced]"*: `Canvas Features`
>                        - In *"JBrowse Feature Score Scaling & Coloring Options [Advanced]"*:
>                            - *"Color Score Algorithm"*: `Based on score`
>                                - *"How should minimum and maximum values be determined for the scores of the features"*: `Manually specify minimum and maximum expected scores for the feature track`
>                                    - *"Minimum expected score"*: `0`
>                                    - *"Maximum expected score"*: `1`
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `GFF/GFF3/BED Features`
>                        - {% icon param-file %} *"GFF/GFF3/BED Track Data"*: `genes (NCBI).gff3`
>                        - *"JBrowse Track Type [Advanced]"*: `Canvas Features`
>                        - In *"JBrowse Contextual Menu options [Advanced]"*:
>                            - In *"Track Menu"*:
>                                - {% icon param-repeat %} *"Insert Track Menu"*
>                                    - *"Menu label"*: `See protein at NCBI`
>                                    - *"Menu url"*: `https://www.ncbi.nlm.nih.gov/gene?term={locus_tag}%5BGene%20Name%5D`
>                                    - *"Menu icon"*: `Database`
>
> > ### {% icon tip %} Tip: Templating in Contextual Menu Options
> >
> > There are more valid values for templating than mentioned in the help. When
> > you click on a feature in the JBrowse instance, it will present all of the
> > properties of the feature. Any of the top level properties can be used
> > directly in your templating
> {: .tip}
>
> 2. View the output
>
> 3. Turn on both tracks of data.
>
> 4. Navigate to `21,200`, either manually, or by copying and pasting the location block: `NC_000913.3:18351..24780`
>
> 5. Right click on the `yaaY` gene, and click the "See protein at NCBI" menu option.
>
>    This menu option is dynamic, try it with a few other features from the `genes (NCBI).gff3` track. These features have a `locus_tag` and the menu button we added will open a URL to an NCBI search page for the value of this `locus_tag` attribute.
{: .hands_on}

Contextual menus can be used to link to more than just NCBI.

- The links can go anywhere such as web search services (e.g. Google) or genomics web services (e.g. EBI)
- Some sites use the IFrame action to link genes to local services where users are expected to submit annotation notes or data.

# Sequencing, Coverage, and Variation

This is the next major category of data that people wish to visualize, sequencing, coverage, and variation. The sequencing data can be of any type, as long as the results are formatted as BAM files.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **JBrowse** {% icon tool %} with the following parameters:
>    - *"Reference genome to display"*: `Use a genome from history`
>        - {% icon param-file %} *"Select the reference genome"*: `genome.fa`
>    - *"Genetic Code"*: `11. The Bacterial, Archaeal and Plant Plastid Code`
>    - In *"Track Group"*:
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - *"Track Category"*: `Coverage`
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `BigWig XY`
>                        - {% icon param-file %} *"BigWig Track Data"*: `dna sequencing coverage.bw`
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `BigWig XY`
>                        - {% icon param-file %} *"BigWig Track Data"*: `rna-seq coverage/1.bw`
>                        - *"Use XYPlot"*: `Yes`
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `BigWig XY`
>                        - {% icon param-file %} *"BigWig Track Data"*: `rna-seq coverage/2.bw`
>                        - *"Use XYPlot"*: `Yes`
>                        - *"Show variance band"*: `Yes`
>                        - *"Track Scaling"*: `Autoscale (global)`
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - *"Track Category"*: `Sequencing & Variation`
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `BAM Pileups`
>                        - {% icon param-file %} *"BAM Track Data"*: `sequencing.bam`
>                        - *"Autogenerate SNP Track"*: `Yes`
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `VCF SNPs`
>                        - {% icon param-file %} *"SNP Track Data"*: `variants.vcf`
>
> 2. Execute and then explore the resulting data.
>
{: .hands_on}


# Blast Results

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **JBrowse** {% icon tool %} with the following parameters:
>    - *"Reference genome to display"*: `Use a genome from history`
>        - {% icon param-file %} *"Select the reference genome"*: `output` (Input dataset)
>    - *"Genetic Code"*: `11. The Bacterial, Archaeal and Plant Plastid Code`
>    - In *"Track Group"*:
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `Blast XML`
>                        - {% icon param-file %} *"BlastXML Track Data"*: `blastp vs swissprot.xml`
>                        - {% icon param-file %} *"Features used in Blast Search"*: `blastp genes.gff3`
>                        - *"Minimum Gap Size"*: `5`
>                        - *"Is this a protein blast search?"*: `Yes`
>                        - In *"JBrowse Feature Score Scaling & Coloring Options [Advanced]"*:
>                            - *"Color Score Algorithm"*: `Based on score`
>                                - *"JBrowse style.color function's score scaling"*: `Blast scaling`
>
> 2. Execute and then explore the resulting data.
>
{: .hands_on}
