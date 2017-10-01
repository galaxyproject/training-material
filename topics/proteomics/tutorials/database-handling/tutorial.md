---
layout: tutorial_hands_on
topic_name: proteomics
tutorial_name: database-handling
---

# Introduction
{:.no_toc}

Identifying peptides in proteomic datasets is commonly done by using search engines that compare the MS2 spectra of peptides to theoretical spectra. These theoretical spectra are generated from a FASTA database containing proteins that are expected in the measured sample. Typically, those FASTA databases will contain all proteins of the organism the sample derived from.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Loading a Search Database

There are a many ways how you can upload your search database (FASTA file with protein sequences). Three among these are:

*   Using **Protein Database Downloader** {% icon tool %} .
*   Using a direct weblink to the database.
*   Uploading a database from the data library.

In this tutorial, we will explore **Protein Database Downloader** {% icon tool %} for generating a searchable protein database. First we download the proteome of the organism of interest. In this tutorial, we will use a database of the human proteome.

> ### {% icon hands_on %} Hands-on: Loading a Search Database
>
> 1. Create a new history for this Database Handling exercise.
> 2. Open **Protein Database Downloader** {% icon tool %}
> 3. Select in the drop-down menues `Taxonomy`: "Homo sapiens (Human)" and `reviewed`: "UniprotKB/Swiss-Prot (reviewed only)".
> 4. Click on `Execute`. There will be a new dataset named `Protein database` in your history, now.
> 5. Rename the `Protein database` to `Main database`.
>
>  > ### {% icon tip %} Tip: Types of uniprot databases
>  > Uniprot offers several types of databases. You may choose to download only reviewed (UniProtKB/Swissprot) databases, only unreviewed (UniProtKB/TREMBL) or both (UniProtKB). In well researched organisms, e.g. Homo sapiens or D. melanogaster, reviewed (Swissprot) databases are always kept up-to-date and may lead to cleaner search results. In other organisms, it might be wiser to include the unreviewed (TREMBL) database not to miss important proteins.
>  >
>  > You may also include protein isoforms by setting the tick box `Include isoform data` to `Yes`.
>  {: .tip}
>
>  > ### {% icon question %} Question
>  > What is the difference between a "reference proteome set" and a "complete proteome set"?
>  >
>  >  <details>
>  >  <summary>Click to view answer!</summary>
>  >  <ol type="1">
>  >  <li> A UniProt complete proteome consists of the set of proteins thought to be expressed by an organism whose genome has been completely sequenced. A reference proteome is the complete proteome of a representative, well-studied model organism or an organism of interest for biomedical research. Reference proteomes constitute a representative cross-section of the taxonomic diversity to be found within UniProtKB. They include the proteomes of well-studied model organisms and other proteomes of interest for biomedical and biotechnological research. Species of particular importance may be represented by numerous reference proteomes for specific ecotypes or strains of interest. [Link to source](http://www.uniprot.org/help/reference_proteome)</li>
>  >  </ol>
>  >  </details>
>  {: .question}
{: .hands_on}


# Contaminant databases

In proteomic samples, some protein contaminants are very common, stemming from the experimenter or contaminated cell culture. Common contaminants are therefore added to the database. This has two benefits:
1. Contamination can be observed, heavily contaminated samples can be excluded from analysis.
2. Contaminant peptides cannot be misassigned to similar peptides in the database reducing the risk of identifying false positives.

A widely used database for common contaminants is the **c**ommon **R**epository of **A**dventitious **P**roteins (cRAP). When using samples generated in cell cultures, it is furthermore recommended to include Mycoplasma proteomes in the search database. Mycoplasma infections are very common in cell culture and often go unnoticed ([Drexler and Uphoff, Cytotechnology, 2002](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3463982/)).

> ### {% icon hands_on %} Hands-on: Contaminant databases
> 1. Open **Protein Database Downloader** {% icon tool %}.
> 2. Select `Download from`: "cRAP (contaminants)" and execute.
> 3. Rename the new database to "crap database".
>
> To be able to distinguish contaminants from proteins of interest, you should add a tag to each contaminant protein.
>
> 1. Run **FASTA-to-Tabular** {% icon tool %} on your crap database.
> 2. Run **Add column** {% icon tool %} on the new output. In the field `Add this value` enter "CONTAMINANT" and execute.
> 3. Run **Tabular-to-FASTA** {% icon tool %}. Use column 1 and column 3 as Title columns and column 2 as sequence column.
> 4. Rename the **Tabular-to-FASTA** {% icon tool %} output to "Tagged cRAP database".
>
>
>  > ### {% icon question %} Question
>  > 1. The cRAP database contains some human proteins. What does it mean if you identify those typical contaminants in a human sample?
>  > 2. What does it mean in a non-human sample?
>  >
>  >  <details>
>  >  <summary>Click to view answers!</summary>
>  >  <ol type="1">
>  >  <li> In samples stemming from a human source, identified human contaminants do not necessarily mean a contaminated sample. The proteins may as well stem from the original sample. Be careful with the interpretation. </li>
>  >  <li> In samples from non-human sources, identified human contaminants do mean contamination by the experimenter. </li>
>  >  </ol>
>  >  </details>
>  {: .question}
{: .hands_on}


> ### {% icon hands_on %} Optional Hands-On: Mycoplasma databases
> 90 - 95 % of mycoplasma infection in cell culture stem from the following species: M. orale, M. hyorhinis, M. arginini, M. fermentans, M. hominis and A. laidlawii ([Drexler and Uphoff, Cytotechnology, 2002](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3463982/)).
>
> 1. Use **Protein Database Downloader** {% icon tool %} to download the six mycoplasma databases. We will merge them to the main database in the next part of the tutorial.
> 2. Run **FASTA Merge Files and Filter Unique Sequences** {% icon tool %} to combine all mycoplasma databases into a single one.
> 3. Tag each entry in the combined database with the string "MYCOPLASMA_CONTAMINANT" by using **FASTA-to-Tabular** {% icon tool %}, **Add column** :wrench: and **Tabular-to-FASTA** :wrench:, as explained [above](#contaminant-databases).
> 4. Rename the **Tabular-to-FASTA** {% icon tool %} output to "Tagged Mycoplasma database".
>
>  > ### {% icon comment %} Comment
>  > The reviewed mycoplasma databases do not contain all known proteins. It is better to include also the TREMBL database. Mycoplasma proteomes are relatively small, so even downloading TREMBL will not bloat your main database unneccessarily.
>  {: .comment}
{: .hands_on}


# Merging databases

Depending on the search engine you are using you might need to merge all fasta entries (i.e. proteins of interest and contaminants) in a single database. Make sure to merge the tagged versions of your contaminant databases.

> ### {% icon hands_on %} Hands-on: Merging databases
>
> 1. Run **FASTA Merge Files and Filter Unique Sequences** {% icon tool %} on the main database and the tagged cRAP database.
>
>	> ### Optional: Merging mycoplasma databases
>	> At this step you may also merge the mycoplasma protein databases that you downloaded earlier on. Simply enter them as additional inputs in **FASTA Merge Files and Filter Unique Sequences** {% icon tool %}. You can enter any number of databases when you click on `Insert Input FASTA file(s)`.
{: .hands_on}


# Creating a Decoy Database

The most common method of peptide and protein FDR calculation is by adding known non-existing proteins, so-called "decoys" to the database. This can be done by reverting or shuffling the real protein entries and adding those fake entries to the database. Some peptide search engines depend on premade decoy databases, others can perform this step automatically.

> ### {% icon hands_on %} Hands-on: Creating a Decoy Database
> 1. Run **DecoyDatabase**  {% icon tool %} on the merged database.
> 2. Set the flag (-append) to `Yes` and execute.
>
>  > ### {% icon tip %} Tip: Decoy tags
>  > The string you enter as a decoy tag will be added as a prefix or suffix (your choice) to the description of each decoy protein entry. Thus you can see from which entry in the real database the decoy was computed.
>  {: .tip}
>
>  > ### {% icon comment %} Comment
>  > **DecoyDatabase**  {% icon tool %} may also take several databases as input which are then automatically merged into one database.
>  {: .comment}
{: .hands_on}


# Concluding remarks
{:.no_toc}

To keep your databases up-to-date, or if you need several databases for different organisms, it would make sense to create a workflow out of the Hands-On sections (to learn about workflows see [this tutorial]({{site.baseurl}}/topics/introduction/tutorials/galaxy-intro-101/tutorial.html)). You might also want to combine the mycoplasma databases to a single file, which you then easily can add to each of your main databases.

Often you may not want to use the most recent database for reasons of reproducibility. If so, you can transfer the final database of this tutorial into other histories to work with it.

Further reading about construction of the optimal database: ([Kumar et al., Methods in molecular biology, 2017](https://www.ncbi.nlm.nih.gov/pubmed/27975281)).

This tutorial is based upon parts of the GalaxyP-101 tutorial (https://usegalaxyp.readthedocs.io/en/latest/sections/galaxyp_101.html).
