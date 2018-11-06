---
layout: tutorial_hands_on

title: Setting up molecular systems
zenodo_link: ''
questions:
- How to get started modelling a protein? 
objectives:
- learn about the Protein Data Bank
- learn how to setup up a model protein system (with CHARMM-GUI)
- learn how to upload the system to Galaxy 
requirements:
  -
    title: "Introduction to Computational chemistry"
    type: "internal"
    link: "/computational-chemistry/slides/introduction.html"
  -
    title: "Setting up molecular systems (slides)"
    type: "internal"
    link: "/computational-chemistry/tutorials/setting-up-molecular-systems/slides.html"
time_estimation: 30m
key_points:
  - "The PDB is a key resource for finding protein structures."
  - "Using CHARMM-GUI is one way to prepare a protein."
  - "To get data into Galaxy you can upload a file from your computer or paste in a web address."
contributors:
  - chrisbarnettster

---

# Introduction
{:.no_toc}

In this tutorial, we'll cover the basics of molecular modelling by setting up a protein and uploading this to Galaxy.

To start we'll look at the PDB and find the entry for an enzyme involved in Chagas Disease. The enzyme is 1S0I, a trans-sialidase.

> ### {% icon tip %} Background: What is the PDB (Protein Data Bank) and format?
>
> The Protein Data Bank (PDB) format contains atomic coordinates of biomolecules and provides a standard representation for macromolecular structure data derived from X-ray diffraction and NMR studies. 
> For example the `PDB`-file for a trans-sialadise with its substrate (PDB: [1S0I](https://www.rcsb.org/pdb/explore/explore.do?structureId=1s0i)):
>
> More resources:
>
>  -  Multiple structures are stored and can be queried at [https://www.rcsb.org/](https://www.rcsb.org/)
>  - Documentation describing the PDB file format is available from the wwPDB at [http://www.wwpdb.org/documentation/file-format.php](http://www.wwpdb.org/documentation/file-format.php).
{: .tip}


> ### {% icon tip %} Background: What is Chagas Disease?
>
>A tropical parasitic disease caused by Trypanosoma cruzi and spread by insect vectors (mainly the kissing bug).
>The trans-sialidase has been of interest as T. Cruzi incorporates sialic acid using this enzyme instead of a specific tranferase enzyme.
>
> More resources:
>
  - [https://en.wikipedia.org/wiki/Chagas_disease](https://en.wikipedia.org/wiki/Chagas_disease) 
  - [https://www.cdc.gov/parasites/chagas/gen_info/detailed.html](https://www.cdc.gov/parasites/chagas/gen_info/detailed.html)
  - [Trypanosoma cruzi Trans-sialidase: structural features and biological implications.](https://www.ncbi.nlm.nih.gov/pubmed/24264246)
  - [Profiling Transition-State Configurations on the Trypanosoma cruzi trans-Sialidase Free-Energy Reaction Surfaces](https://pubs.acs.org/doi/full/10.1021/jp506824r)
{: .tip}



> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Trans-Sialidase PDB

In this section we'll access the PDB, download the correct structure, import it and view in Galaxy.

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from the [PDB](https://files.rcsb.org/download/1S0I.pdb) or from the shared data library
>
>    ```
>    
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


## Sub-step with **My Tool**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **My Tool** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: File
>    - *"Parameter"*: `a value`
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

{% icon trophy %} Well done! You have started modelling a protein and uploaded it into Galaxy.
Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
