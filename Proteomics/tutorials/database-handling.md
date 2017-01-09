---
layout: tutorial_hands_on
topic_name: Database-Handling
tutorial_name: Database-Handling
---

This tutorial was based upon parts of the Galaxy-P 101 tutorial (https://usegalaxyp.readthedocs.io/en/latest/sections/galaxyp_101.html).

# Introduction

Identifying peptides in proteomic datasets is commonly done by using search engines that compare the MS2 spectra of peptides to theoretical spectra. These theoretical spectra are generated from a FASTA database containing proteins that are expected in the measured sample. Typically, those FASTA databases will contain all proteins of the organism the sample derived from.



The most common method of FDR calculation is by adding known non-existing proteins, so-called "decoys" to the database. This can be done by reverting or shuffling the real protein entries and adding the fake entries to the real database. Some search engines depend on decoy databases, others can perform this step automatically.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. [Loading a Search Database](#Loading-a-Search-Database)
> 2. [Merging Databases](#Merging-Databases) 
> 3. [Creating Decoy Databases](#Creating-Decoy-Databases)

# Loading a Search Database

There are a many options for how you can upload your search database (FASTA file with protein sequences). Three among these are:

*   Protein Database Downloader.
*   Use website link for the database.
*   Upload database from the data library.

In this tutorial, we will explore using Protein Database Downloader for database search. First we download the proteome of the organism of interest. In this tutorial, we will use the human proteome. 

> ### :pencil2: Hands-on: Loading a Search Database
>
> 1. Create a new history for this Database Handling exercise.
> 2. Go to Tools –> Get Data. Then click on Protein Database Downloader.
> 3. Select in the drop-down menues Taxonomy: Homo sapiens (Human) and reviewed: UniprotKB/Swiss-Prot (reviewed only).
> 4. Click on “Execute”.
>	>### :bulb: Tip: Types of uniprot databases
>	> Uniprot offers several types of databases. You may choose to download only reviewed (UniProtKB/Swissprot) databases, only unreviewed (UniProtKB/TREMBL) or both (UniProtKB).
>	> You may also include protein isoforms.

# Contaminant databases

In proteomic samples, some protein contaminants are very common, stemming from the experimenter or contaminated cell culture. Common contaminants are therefore added to the database. This has two benefits: 
1. Contamination can be observed, heavily contaminated samples can be excluded from analysis.
2. Contaminant peptides are not misassigned to similar peptides in the database reducing the risk of identifying false positives.

A widely used database for common contaminants is the **c**ommon **R**epository of **A**dventitious **P**roteins (cRAP). When using samples generated in cell cultures, it is furthermore recommended to include Mycoplasma proteomes in the search database. Mycoplasma infections are very common in cell culture and often go unnoticed ([Drexler and Uphoff 2002](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3463982/)).

> ### :pencil2: Hands-on: Contaminant databases
> 1. Open Protein Database Downloader.
> 2. Select Download from: cRAP (contaminants).
> 3. Click on "Execute".

## [Using the Database Merge Tool](#using-the-database-merge-tool)

In order to merge these two FASTA files in your history, go to Tools –> FASTA Manipulation –> Merge FASTA Databases.

![../_images/101_8.png](../_images/101_8.png)

Click on “Add new Input FASTA File(s)” – twice – so that your set up for merging files looks like this:

![../_images/101_9.png](../_images/101_9.png)

Click on Execute and you will be able to see your third History Item in Galaxy’s right pane.

![../_images/101_10.png](../_images/101_10.png)</div>

<div class="section" id="creating-a-decoy-database">

## Creating a Decoy Database[¶](#creating-a-decoy-database "Permalink to this headline")

Our next step is to create a decoy database out of the merged file (3rd history item on the list).

For this, go to Tools -> FASTA Manipulation –> ‘Create Decoy Database’.

![../_images/101_83.png](../_images/101_83.png)

Ensure that history item 3 shows up in the “FASTA Input:” box and that the box for “Include original entries in output database:” is checked. Change Decoy prefix to [<span class="problematic" id="id2">REV_</span>](#id1).

Your parameters for creating decoy database (reverse) tool should look like this.

![../_images/101_12.png](../_images/101_12.png)

Click on “Execute” to generate the fourth item in your history list. This is the FASTA database that we will be using for our search.

![../_images/101_13.png](../_images/101_13.png)

Now we will rename the history items to “Human UniProt”, “cRAP”, “Merged Human UniProt cRAP” and “Target_Decoy_Human_Contaminants” by clicking on the Pencil icon adjacent to each item. Also we will rename history to “Galaxy-P 101” (or whatever you want) by clicking on “Unnamed history” so everything looks like this:

![../_images/101_14.png](../_images/101_14.png)

Please feel free to explore tabs – Convert Format, Datatype, Permissions while you are editing the attributes. This is especially important while troubleshooting for steps that fail wherein a datatype has not been set properly or needs to be changed for subsequent steps.

![../_images/101_15.png](../_images/101_15.png)</div>

