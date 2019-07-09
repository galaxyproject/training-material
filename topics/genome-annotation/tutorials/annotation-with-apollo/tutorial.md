---
layout: tutorial_hands_on

title: Refining Manual Genome Annotations with Apollo
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
time_estimation: 3H
key_points:
- Apollo allows a group to view and manually refine predicged genome annotations
- Use Apollo to edit annotations within your group.
- Export manual annotations as GFF3.
contributors:
- nathandunn
---


# Introduction
{:.no_toc}

After automatically editing annotations using 
[Prokker](https://training.galaxyproject.org/training-material/topics/genome-annotation/tutorials/annotation-with-prokka/tutorial.html) or 
[Maker](https://training.galaxyproject.org/training-material/topics/genome-annotation/tutorials/annotation-with-maker/tutorial.html),
 its important to visualize and then manually refine any additional data. 

This process is most often done as part of a group.  

This demo is inspired by the [Apollo User's Guide](http://genomearchitect.github.io/users-guide/), which provides additional guidance. 

***TODO***: *Add links to the Apollo User's Gropu and GGA*

{% raw %} `{% cite Dunn2019 %}`{% endraw %}


**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data upload

To annotate a genome using Apollo, you simply need the reference genome sequence in FASTA format.

We will also need to provide evidence for our refined annotations.  

This will result in several types of evidence:

- If not *de novo*, then this could include a previous released official genome annotation.
- A set of prior gene predictions or other genomic feature types if available (provided for this example) typically as GFF3.  
- Individual read files in BAM format if available.
- Additional visual data to help indicate things like read density such as BigWig files.


## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo](https://doi.org/10.5281/zenodo.3270822) or from the shared data library
>
>    ```
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/Amel_4.5_scaffolds.fa.gz
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/amel_OGSv3.2.gff3.gz
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/forager_Amel4.5_accepted_hits.bam
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/forager_Amel4.5_accepted_hits.bam.bai
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/forager.bw
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/nurse_Amel4.5_accepted_hits.bam
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/nurse_Amel4.5_accepted_hits.bam.bai
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/nurse.bw
>    ```
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

# Setup Apollo for Annotation

Refining genomes happens in multiple steps:

- register a user, injecting the user from Galaxy into Apollo
- create a viewable genome or organism from the reference genome FASTA file
- add genomic evidence
- refine the genome (more on that)
- export the refined genome annotations


(SCREEN SHOT TO ADD of imporection)
![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

> ### {% icon details %} Why bother?
>
> Automated annotation programs continue to improve, however a simple score may not provide the visual evidence necessary to confirm an accurate prediction.
> Therefore, it is necessary to both visually inspect as well as fix any issues with the predicted genomes.
> 
> Additionally, many times assemblies are less than perfect or read depth and quality may be insufficient.
> 

{: .details}



## Sub-step with **JBrowse**

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
>    - *"Species"*: `Mellifera`
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

