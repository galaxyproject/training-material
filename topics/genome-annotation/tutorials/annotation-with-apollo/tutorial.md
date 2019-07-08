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

# Title for your first section

Give some background about what the trainees will be doing in the section.
Remember that many people reading your materials will likely be novices,
so make sure to explain all the relevant concepts.

## Title for a subsection
Section and subsection titles will be displayed in the tutorial index on the left side of
the page, so try to make them informative and concise!

# Hands-on Sections
Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/Amel_4.5_scaffolds.fa.gz
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/amel_OGSv3.2_cds.fa.gz
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/amel_OGSv3.2.gff3.gz
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/amel_OGSv3.2_pep.fa.gz
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/forager_Amel4.5_accepted_hits.bam
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/forager_Amel4.5_accepted_hits.bam.bai
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/forager.bw
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/nurse_Amel4.5_accepted_hits.bam
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/nurse_Amel4.5_accepted_hits.bam.bai
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/yeast_chr1%2B2.gff3
>    https://zenodo.org/api/files/55133323-b15b-45b4-98c9-dda00288b53f/yeast.fa
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
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

# Title of the section usually corresponding to a big step in the analysis

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


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
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

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
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

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

## Sub-step with **Annotate**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Annotate** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Apollo Organism Listing"*: `output` (output of **Create or Update Organism** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

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


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
