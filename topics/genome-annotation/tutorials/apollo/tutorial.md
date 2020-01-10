---
layout: tutorial_hands_on

title: Apollo
zenodo_link: https://zenodo.org/record/3270822
tags:
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
time_estimation: 3h
key_points:
- Apollo allows a group to view and manually refine predicged genome annotations
- Use Apollo to edit annotations within your group.
- Export manual annotations as GFF3.
contributors:
- abretaud
- erasche
- nathandunn

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

> {% icon warning %} Only works on UseGalaxy.eu
> Currently this tutorial requires an Apollo server to be deployed by the administrator. This will currently only work on UseGalaxy.eu, hopefully this list will expand in the future.
{: .warning}

# Introduction
{:.no_toc}

After automatically annotating your genome using [Prokka](../annotation-with-prokka/tutorial.html) or [Maker](../annotation-with-maker/tutorial.html), it is important to visualize your results so you can understand what your organism looks like, and then to manually refine these annotations along with any additional data you might have. This process is most often done as part of a group, smaller organisms may be annotated individually though.

[Apollo](https://github.com/gmod/apollo) {% cite Dunn2019 %} provides a platform to do this, it is a web-based, collaborative genome annotation editor. Think of it as "Google Docs" for genome annotation, multiple users can work together simultaneously to curate evidences and annotate a genome.

This demo is inspired by the [Apollo User's Guide](https://genomearchitect.readthedocs.io/en/latest/UsersGuide.html), which provides additional guidance. 

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data upload

To annotate a genome using Apollo, we need the reference genome sequence in FASTA format, and any evidence tracks we want to refine into our annotations. "Evidence tracks" can be any data like:

- A set of prior gene predictions or other genomic feature predictions
- The output of a bioinformatics analysis like BLAST or InterProScan
- Sequencing reads from RNA-Seq or another HTS analysis
- If you are not doing a *de novo* annotation, then a previous released <abbr title="Official Gene Set">OGS</abbr>

In this tutorial we have obtained some data from NCBI related to [*Escherichia coli K12 str. MG1655*](https://ecoliwiki.org/colipedia/index.php/Category:Strain:MG1655), and we will visualise this data and use it to make some annotations in order to familiarise you with the process.

> ### {% icon comment %} Real Data: Unreal Circumstances
> While the data for this tutorial is sourced from publicly available databases, and is all related to different experiments on *E. coli K12*, this is not necessarily the data *you* might use to annotate your genomes. You probably know best what data you should be using in your own circumstances, for the specific features on which you are focused.
{: .comment}

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 0. Create a new history and give it a good name
>
>    {% include snippets/create_new_history.md %}
>    {% include snippets/rename_history.md %}
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

Refining genome annotations happens in multiple steps:

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
>        - {% icon param-file %} *"Select the reference genome"*: Select the genome fasta file
>    - *"Genetic Code"*: `11. The Bacterial, Archael`
>    - In *"Track Group"*:
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - *"Track Category"*: `Gene Calls`
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `GFF/GFF3/BED/GBK Features`
>                        - {% icon param-files %} *"GFF/GFF3/BED Track Data"*: `Augustus` and `NCBI AnnotWriter Genes`
>                        - *"JBrowse Track Type [Advanced]"*: `Canvas Features`
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - *"Track Category"*: `Sequencing`
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `BAM Pileups`
>                        - {% icon param-files %} *"BAM Track Data"*: Both BWA-MEM Mappings
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `BigWig XY`
>                        - {% icon param-files %} *"BAM Track Data"*: Both of the BWA-MEM Coverage files (**not** the `(as bigwig)` files)
>                        - *"Use XYPlot"*: `Yes`
>                        - *"Show Variance Band"*: `Yes`
>                        - *"Track Scaling"*: `Autoscale (local)`
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - *"Track Category"*: `RNA-Seq`
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `BAM Pileups`
>                        - {% icon param-files %} *"BAM Track Data"*: Both TopHat Mappings
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `BigWig XY`
>                        - {% icon param-files %} *"BAM Track Data"*: Both of the `TopHat ... Coverage` files
>                        - *"Use XYPlot"*: `Yes`
>                        - *"Show Variance Band"*: `Yes`
>                        - *"Track Scaling"*: `Autoscale (local)`
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - *"Track Category"*: `Variation`
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `VCF SNPs`
>                        - {% icon param-files %} *"SNP Track Data"*: Both FreeBayes files
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - *"Track Category"*: `Similarity`
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `GFF/GFF3/BED Features`
>                        - {% icon param-file %} *"GFF/GFF3/BED Track Data"*: `O104:H4 LASTZ Alignments`
>                        - *"JBrowse Track Type [Advanced]"*: `Canvas Features`
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `Blast XML`
>                        - {% icon param-file %} *"Blast XML Track Data"*: The `blastp` results from swissprot
>                        - {% icon param-file %} *"Features used in Blast Search"*: The `NCBI AnnotWriter Genes` file
>                        - *"Minimum Gap Size"*: `3`
>                        - *"Is this a protein blast search?"*: `Yes`
>
>    > ### {% icon comment %} JBrowse is highly configurable
>    >
>    > JBrowse is highly configurable, we have set a very basic configuration but there are many more advanced features available to you, if you need them. You can choose precisely how data is displayed, and even what menu options are available when users click on features. If your features have some external identifiers like an NCBI Gene ID, you can even configure JBrowse that when the user clicks on the feature, it should show the gene page for that feature in a new tab. These sort of features are incredibly helpful for building very rich experiences.
>    >
>    > A static genome browser like this (just JBrowse, not in Apollo) is very useful for summarising results of a genomics workflow, where the next step is simply interpretation and not annotation.
>    >
>    > Currently we have built a standalone genome browser (data + the html page and user interface and javascript), but it's possible to just compile the data directory if you intend to send this data to Apollo, and don't need to view the static data in Galaxy.
>    {: .comment}
{: .hands_on}

This tool will take some time to run dependent on data size. All of the inputs need to be pre-processed by JBrowse into a form that it can render and visualise easily. Once this is complete, you can click on the {% icon galaxy-eye %} eyeball to view the JBrowse instance. This is a static view into the data, JBrowse does not let you make any annotations or save any changes. We will convert it into a dynamic view where we can make persistent annotations and share these with our colleagues.

## Sending data to Apollo

> ### {% icon hands_on %} Import to Apollo
>
> 1. **Create or Update Organism** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"JBrowse HTML Output"*: output of **JBrowse** {% icon tool %}
>    - *"Organism Common Name Source"*: `Direct Entry`
>        - *"Organism Common Name"*: `E. coli K12`
>    - *"Genus"*: `Escherichia`
>    - *"Species"*: `coli`
>
> 2. **Annotate** {% icon tool %}:
>    - {% icon param-file %} *"Apollo Organism Listing"*: output of **Create or Update Organism** {% icon tool %}
>
> 3. View {% icon galaxy-eye %} the output of the Annotate tool, when it is ready.
>
{: .hands_on}

Viewing the output will open a view into Apollo in the main panel. Here you can interact with your genome and make annotations. This "Annotate" output is a quick link to that specific genome, and while Apollo allows you to manage and annotate multiple genomes, this dataset will always take you back to that specific genome. You can additionally access the Apollo server outside of Galaxy. While the URL will be different for each Galaxy server that supports Apollo, UseGalaxy.eu's Apollo server it available at [https://usegalaxy.eu/apollo](https://usegalaxy.eu/apollo).
{: .warnin}


# Apollo


![Apollo in Galaxy](../../images/apollo/ui.png "This is the Apollo interface, when viewed from inside Galaxy. The Annotate tool produces a special output which acts like a link to the Apollo server. It is <b>not</b> sufficient to download this file to save a copy of your annotations, unlike other datasets in Galaxy. In the main panel Apollo is displayed, with two panels of its own: the Annotation Window (left), and the Annotator Panel (right). The annotation window is the main view into our genome and here we will see evidence, reconcile it, and make our annotations. The right panel allows us to show or hide individual evidence tracks, switch between organisms, and navigate around the genome.")

From the Apollo user manual:

> The major steps of manual annotation using Apollo can be summarized as follows:
>
> 1. Locate a chromosomal region of interest.
> 2. Determine whether a feature in an existing evidence track provides a reasonable gene model to start annotating.
> 3. Drag the selected feature to the ‘User Annotation’ area, creating an initial gene model.
> 4. Use editing functions to edit the gene model if necessary.
> 5. Check your edited gene model for consistency with existing homologs by exporting the FASTA formatted sequence and searching a protein sequence database, such as UniProt or the NCBI Non Redundant (NR) database, and by conducting preliminary functional assignments using the Gene Ontology (GO) database.
>
{: .quote}

The first four steps are generally the process of structural annotation (the process of identifying the correct gene model), and the last includes functional annotation (the process of assigning a putative function to a gene in your annotations).

## Structural Annotation

Let's start by looking at the tracks available to us, and then turning on the gene call tracks so we can start exploring our data.

> ### {% icon hands_on %} Visualize the Gene Calls
>
> 1. In the right hand panel at the top click on **Tracks** to open the track listing
>
>    ![Track menu](../../images/apollo/tracks.png)
>
> 2. In the **Gene Calls** group, click the checkbox to the right.
>
>    You can either activate tracks in bulk, by clicking on the checkbox to the right of the group name ("Gene Calls"), or by clicking on the group name to expand the section, and then selecting individual tracks.
>
> 3. Zoom to the first 10kb of the genome.
>
>    1. In the left hand Annotation Window, at the top navigation bar you will find a textbox which shows the current location on the genome.
>    2. Edit this and enter `1..10000`
>    3. Press *Go* or use <kbd>Enter</kbd> on your keyboard.
>
>    ![JBrowse Location Bar](../../images/apollo/location.png)
>
{: .hands_on}

We can now see two evidence tracks: one is the set of [genes from NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=511145&lvl=3&lin=f&keep=1&srchmode=1&unlock), the other is the output of [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus) {% cite Stanke_2008 %}. In a *de novo* annotation project, we probably will only have the outputs of various gene callers, and potentially some expression evidence like RNA-Seq.

## Create refinements

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
