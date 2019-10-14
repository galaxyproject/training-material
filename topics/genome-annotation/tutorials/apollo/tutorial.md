---
layout: tutorial_hands_on

title: Apollo
zenodo_link: https://zenodo.org/record/3270822
tags:
- intro
- eukaryote
questions:
- How to visualize your genome after automated annotations have been performed?
- How to manually annotate genome after automated annotations have been performed?
- How to evaluate and visualize annotated genomic features?
objectives:
- Load genome into Galaxy
- View annotations in JBrowse
- Learn how to load JBrowse data into Apollo
- Learn how to manually refine genome annotations within Apollo
- Export refined genome annotations
time_estimation: 3H
key_points:
- Apollo allows a group to view and manually refine predicged genome annotations
- Use Apollo to edit annotations within your group.
- Export manual annotations as GFF3.
contributors:
- nathandunn
- erasche

requirements:
- type: "internal"
  topic_name: galaxy-data-manipulation
  tutorials:
    - upload-rules
#- type: "internal"
#  topic_name genome-annotation
#  tutorials:
#    -
---


# Introduction
{:.no_toc}

After automatically annotating your genome using [Prokka](../annotation-with-prokka/tutorial.html) or [Maker](../annotation-with-maker/tutorial.html), its important to visualize and then manually refine any additional data.
This process is most often done as part of a group, smaller organisms may be annotated individually though.

[Apollo](https://github.com/gmod/apollo) {% cite Dunn2019 %} provides a platform to do this, it is a web-based, collabortive genome annotation editor. Think of it as "Google Docs" for genome annotation, multiple users can work together to annotate a genome.

This demo was inspired by the [Apollo User's Guide](http://genomearchitect.github.io/users-guide/), which provides additional guidance. This tutorial was further developed by the [Galaxy Genome Annotation](https://galaxy-genome-annotation.github.io/) group who developed the Galaxy-Apollo bridge.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


![Picture of A. mel](https://upload.wikimedia.org/wikipedia/commons/thumb/f/f7/Bee_on_Lavender_Blossom_2.jpg/800px-Bee_on_Lavender_Blossom_2.jpg "<i>Apis mellifera</i> By Martin Falbisoner - <span class="int-own-work" lang="en">Own work</span>, <a href="https://creativecommons.org/licenses/by-sa/4.0" title="Creative Commons Attribution-Share Alike 4.0">CC BY-SA 4.0</a>, <a href="https://commons.wikimedia.org/w/index.php?curid=41709017">Link</a>"){: style="float:right;max-width:25%"}

# Data upload

To annotate a genome using Apollo, we need the reference genome sequence in FASTA format, and any evidence tracks we want to refine into our annotations. "Evidence tracks" can refer to something like:

- A set of prior gene predictions or other genomic feature predictions
- The output of a bioinformatic analysis like BLAST or InterProScan
- Sequencing reads from RNA-Seq or another HTS analysis
- If you are not doing a *de novo* annotation, then a previous released <abbr title="Official Gene Set">OGS</abbr>

In this tutorial we have obtained some data from [*Apis mellifera*](https://en.wikipedia.org/wiki/Apis_mellifera), this data comes from published scaffolds from the Apis mellifera assembly Amel_4.5 and Official Gene Set 3.2, available from [BeeBase](http://hymenopteragenome.org/beebase/?q=download_sequences).

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Click the upload icon {% icon galaxy-upload %}
>
> 2. Switch to the "Rule-based" tab
>
> 3. Copy & Paste the following table into the Rule-based uploader textbox:
>
>    ```
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/Amel_4.5_scaffolds.fa.gz	Scaffolds	fasta.gz
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/amel_OGSv3.2.gff3.gz	OGS v3.2	gff3
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/forager_Amel4.5_accepted_hits.bam	forager_Amel4.5_accepted_hits	bam
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/forager.bw	forager coverage	bigwig
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/nurse_Amel4.5_accepted_hits.bam	nurse_Amel4.5_accepted_hits	bam
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/nurse.bw	nurse coverage	bigwig
>    ```
> 4. Click **Build**
>
> 5. From **Rules** menu select `Add / Modify Column Definitions`
>
>    - Click `Add Definition` button and select `Name`
>      - *"Name"*: `B`
>    - Repeat this again and select `URL` instead.
>      - *"URL"*: `A`
>    - Repeat this again and select `URL` instead.
>      - *"Type"*: `C`
>    - Click `Apply`
>
> 6. At the bottom of the dialog on the left, set your Genome to `A. mellifera 04 Nov 2010 (Amel_4.5/apiMel4)`
>
> 7. Click **Upload**
>
{: .hands_on}

# Using Apollo for Annotation

Refining genomes happens in multiple steps:

- Create a JBrowse instance from the reference genome FASTA file and evidence tracks
- Import this data into Apollo
- Refine the annotations
- Export the refined genome annotations

In this tutorial we will focus more on the practical portions than the theoretical part of genome annotation, that will be covered in other tutorials. When you've completed this tutorial you should be comfortable manipulating genomic data in Galaxy and Apollo.

> ### {% icon details %} Why bother?
>
> Automated annotation programs continue to improve, however a simple score may not provide evidence necessary to confirm an accurate prediction.
> Therefore, it is necessary to both visually inspect the results and manually fix any issues with the predictions.
>
> Additionally, many times assemblies are less than perfect or read depth and quality may be insufficient.
{: .details}

## Build the JBrowse Instance



> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **JBrowse** {% icon tool %} with the following parameters:
>    - *"Reference genome to display"*: `Use a genome from history`
>        - {% icon param-file %} *"Select the reference genome"*: `output` (Input dataset)
>    - *"JBrowse-in-Galaxy Action"*: `New JBrowse Instance`
>    - In *"Track Group"*:
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - *"Track Category"*: `OGS 3.2`
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `GFF/GFF3/BED/GBK Features`
>                        - {% icon param-file %} *"GFF/GFF3/BED Track Data"*: `output` (Input dataset)
>                        - *"This is match/match_part data"*: `Yes`
>                        - *"JBrowse Track Type [Advanced]"*: `Neat HTML Features`
>                        - In *"JBrowse Feature Score Scaling & Coloring Options [Advanced]"*:
>                            - *"Color Score Algorithm"*: `Ignore score`
>                                - *"Color Selection"*: `Automatically selected`
>                        - *"Track Visibility"*: `On for new users`
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - *"Track Category"*: `Reads`
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `BAM Pileups`
>                        - {% icon param-file %} *"BAM Track Data"*: `output` (Input dataset)
>                        - *"Autogenerate SNP Track"*: `Yes`
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `BigWig XY`
>                        - *"Track Scaling"*: `Autoscale (local)`
>                        - In *"JBrowse Color Options [Advanced]"*:
>                            - *"Color Specification"*: `Automatically selected`
>                            - *"Bicolor Pivot"*: `Zero`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} JBrowse is highly configurable
>    >
>    > JBrowse is highly configurable and can be run standalone.
>    > A refined genome can be published separate from Apollo or used as further evidence.
>    {: .comment}

{: .hands_on}

You will now see three new tracks displaying all the evidences used by Maker to generate consensus gene models.


## Sub-step with **Create or Update Organism**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Create or Update Organism** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"JBrowse HTML Output"*: `output` (output of **JBrowse** {% icon tool %})
>    - *"Organism Common Name Source"*: `Direct Entry`
>        - *"Organism Common Name"*: `Honeybee`
>    - *"Genus"*: `Apis`
>    - *"Species"*: `mellifera`
>    - *"Is Organism Public"*: `Yes`
>    - *"Grant access to a user group"*: ``
>    - *"Remove old data directory"*: `Yes`
>
{: .hands_on}

## Sub-step with **Annotate**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Annotate** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Apollo Organism Listing"*: `output` (output of **Create or Update Organism** {% icon tool %})
>
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. Why DO we bother?
>
> > ### {% icon solution %} Solution
> >
> > 1. Annotation depth may be insufficient and assemblies may be incorrect.
> >
> {: .solution}
>
{: .question}

# Create refinements

View data at a particular position

## Search for a gene

Enter the gene XXX into the gene box.

Zoom in to the proper region.

## Create structural edits

Drag annotation from evidence to HTML region.

Conversely, if you right-click you can any type of genome feature annotation.

(SCREEN SHOT TO ADD of imperfection)
![Alternative text](../../images/image_name "Legend of the image")

Additional isoforms may be dragged up from the evidence.

## Edit structure

### Update exon position
Once isoforms have been created, the edges may be dragged to best match the biological evidence.

CDS's are automatically updated.

Conversely by selecting "choose the annotation" the individual code view is selected.

You will also notice that overlapping isoforms are highlighted.

(SCREEN SHOT of isoforms and dragging exons)
![Alternative text](../../images/image_name "Legend of the image")

Genes will automatically be predicted based on CDS overlap.  This can be unassigned by deselecting on the right-click menu.

By right-clicking on the refined genome feature the details of the genome feature can be retrieved quite readily.

(SCREEN SHOT of right-click menu)
![Alternative text](../../images/image_name "Legend of the image")

### View structured data

Selecting the features allows us to view the gene directly.

(SCREEN SHOT of feature menu)
![Alternative text](../../images/image_name "Legend of the image")


### Edit structured data

There are various structured data options from the figure.

All structured data


#### Editing and reverting history



(SCREEN SHOT of history menu)
![Alternative text](../../images/image_name "Legend of the image")


### Edit functional data

There is various functional data.


#### Edit names, etc.


#### Add comments, keys, and values



#### Show GO Annotations


## Executing the workflow


# Export refinements



# Add users to help with the refinement



***TODO***: *Add to workflow*


# Conclusion
{:.no_toc}

Congratulations, you finished this tutorial! You learned how to manually refine predicted eukaryotic genomes using Apollo and export them to other forms.

By using Apollo and JBrowse you can inspect and refine genome annotations with other researchers.
When refinement is sufficient an updated or new version of the genome may be exported as GFF3 as well as published as a new JBrowse directory for inspection.

# What's next?

After generating your refined genome, you'll want to merge it back into the official gene sets.

If a de novo set, you can export it as GFF3 and load it into a tool like [Tripal](http://tripal.info) to provide visualization.

