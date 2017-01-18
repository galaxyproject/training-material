---
layout: tutorial_hands_on
topic_name: Database-Handling
tutorial_name: Database-Handling
---

This tutorial was based upon parts of the Galaxy-P 101 tutorial (https://usegalaxyp.readthedocs.io/en/latest/sections/galaxyp_101.html).

# Introduction

Identifying peptides in proteomic datasets is commonly done by using search engines that compare the MS2 spectra of peptides to theoretical spectra. These theoretical spectra are generated from a FASTA database containing proteins that are expected in the measured sample. Typically, those FASTA databases will contain all proteins of the organism the sample derived from.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. [Loading a Search Database](#Loading-a-Search-Database)
> 2. [Merging Databases](#Merging-Databases) 
> 3. [Creating Decoy Databases](#Decoy-database)

<a name="Loading-a-Search-Database"/></a>
# Loading a Search Database

There are a many options for how you can upload your search database (FASTA file with protein sequences). Three among these are:

*   Protein Database Downloader.
*   Use website link for the database.
*   Upload database from the data library.

In this tutorial, we will explore using Protein Database Downloader for database search. First we download the proteome of the organism of interest. In this tutorial, we will use the human proteome. 

> ### :pencil2: Hands-on: Loading a Search Database
>
> 1. Create a new history for this Database Handling exercise.
> 2. Open **Protein Database Downloader** :wrench: 
> 3. Select in the drop-down menues `Taxonomy`: "Homo sapiens (Human)" and `reviewed`: "UniprotKB/Swiss-Prot (reviewed only)".
> 4. Click on “Execute”. There will be a new dataset "Protein database" in your history, now.
> 5. Rename the "Protein database" to "Main database".
>
>	> ### :bulb: Tip: Types of uniprot databases
>	> Uniprot offers several types of databases. You may choose to download only reviewed (UniProtKB/Swissprot) databases, only unreviewed (UniProtKB/TREMBL) or both (UniProtKB). In well researched organisms, e.g. Homo sapiens or D. melanogaster, reviewed (Swissprot) databases are always kept up-to-date and may lead to cleaner search results. In other organisms, it might be wiser to include the unreviewed (TREMBL) database not to miss important proteins.
>	>
>	> You may also include protein isoforms.
>
>	> ### :question: Question
>	> What is the difference between a "reference proteome set" and a "complete proteome set"?

# Contaminant databases

In proteomic samples, some protein contaminants are very common, stemming from the experimenter or contaminated cell culture. Common contaminants are therefore added to the database. This has two benefits: 
1. Contamination can be observed, heavily contaminated samples can be excluded from analysis.
2. Contaminant peptides cannot be misassigned to similar peptides in the database reducing the risk of identifying false positives.

A widely used database for common contaminants is the **c**ommon **R**epository of **A**dventitious **P**roteins (cRAP). When using samples generated in cell cultures, it is furthermore recommended to include Mycoplasma proteomes in the search database. Mycoplasma infections are very common in cell culture and often go unnoticed ([Drexler and Uphoff, Cytotechnology, 2002](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3463982/)).

> ### :pencil2: Hands-on: Contaminant databases
> 1. Open **Protein Database Downloader** :wrench:. 
> 2. Select `Download from`: "cRAP (contaminants)" and execute.
> 3. Rename the new database to "crap database".
>
> ### :question: Question
> The cRAP database contains some human proteins. What does it mean if you identify those typical contaminants in a human sample? What does it mean in a non-human sample?
>    <details>
>    <summary>Click to view answers</summary>
>    	<ol type="1">
>    		<li> In samples stemming from a human source, identified human contaminants do not necessarily mean a contaminated sample. The proteins may as well stem from the original sample. Be careful with the interpretation. </li>
			<li> In samples from non-human sources, identified human contaminants do mean contamination by the experimenter. </li>
>   	 </ol>
>    </details>

> ### :pencil2: Optional Hands-On: Mycoplasma databases
> 90 - 95 % of mycoplasma infection in cell culture stem from the following species: M. orale, M. hyorhinis, M. arginini, M. fermentans, M. hominis and A. laidlawii ([Drexler and Uphoff, Cytotechnology, 2002](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3463982/)).
>
> 1. Use **Protein Database Downloader** :wrench: to download the six Mycoplasma databases. We will merge them to the main database in the next part of the tutorial.
> 
> ### :nut_and_bolt: Comment
> The reviewed mycoplasma databases do not contain all known proteins. It is better to include also the TREMBL database. Mycoplasma proteomes are very small, so even downloading TREMBL will not bloat your main database unneccessarily.

<a name="Merging-Databases"/></a>
# Merging databases

Depending on the search engine you are using you might need to merge all fasta entries in a single database. Before doing so, you should add an identifier to all cRAP entries to distinguish them from the proteins of interest.

> ### :pencil2: Hands-on: Merging databases
> First we will add the tag "CONTAMINANT" to each entry in the cRAP database.
>
> 1. Run **FASTA-to-Tabular** :wrench: on your crap database.
> 2. Run **Add column** :wrench: on the new output. In the field `Add this value` enter "CONTAMINANT" and execute.
> 3. Run **Tabular-to-FASTA** :wrench:. Use column 1 and column 3 as Title columns and column 2 as sequence column.
> 4. Rename the output to "Tagged cRAP database".
>	
> Now we can merge the databases:
>
> 1. Run **FASTA Merge Files and Filter Unique Sequences** :wrench: on the main database and the tagged cRAP database.
>
>	> ### Optional: Merging mycoplasma databases
>	> At this step you may also merge the mycoplasma protein databases that you downloaded earlier on. Simply enter them as additional inputs in **FASTA Merge Files and Filter Unique Sequences** :wrench:. You can enter any number of databases when you click on `Insert Input FASTA file(s)`.

<a name="Decoy-database"/></a>
# Creating a Decoy Database

The most common method of peptide and protein FDR calculation is by adding known non-existing proteins, so-called "decoys" to the database. This can be done by reverting or shuffling the real protein entries and adding those fake entries to the database. Some peptide search engines depend on premade decoy databases, others can perform this step automatically.

> ### :pencil2: Hands-on: Creating a Decoy Database
> 1. Run **DecoyDatabase**  :wrench: on the merged database.
> 2. Set the flag (-append) to `Yes` and execute.
>
>	> ### :bulb: Tip: Decoy tags
>	> The string you enter as a decoy tag will be added as a prefix or suffix (your choice) to the description of each decoy protein entry. Thus you can see from which entry in the real database the decoy was computed.
>	>
>	> ### :nut_and_bolt: Comment
>	> **DecoyDatabase**  :wrench: may also take several databases as input which are then automatically merged into one database.

# Concluding remarks
To keep your databases up-to-date, or if you need several databases for different organisms, it would make sense to create a workflow out of the Hands-On sections (to learn about workflows see [this tutorial](../../Introduction/tutorials/workflows.md)). You might also want to combine the mycoplasma databases to a single file, which you then easily can add to each of your main databases.

Often you may not want to use the most recent database for reasons of reproducibility. If so, you can transfer the final database of this tutorial into other histories to work with it.

Further reading about construction of the optimal database: ([Kumar et al., Methods in molecular biology, 2017](https://www.ncbi.nlm.nih.gov/pubmed/27975281)).