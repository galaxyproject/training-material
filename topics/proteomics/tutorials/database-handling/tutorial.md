---
layout: tutorial_hands_on

title: "Protein FASTA Database Handling"
edam_ontology: "topic_0121"
zenodo_link: ""
level: Introductory
questions:
  - "How to download protein FASTA databases of a certain organism?"
  - "How to download a contaminant database?"
  - "How to create a decoy database?"
  - "How to combine databases?"
objectives:
  - "Creation of a protein FASTA database ready for use with database search algorithms."
time_estimation: "30m"
key_points:
  - "There are several types of Uniprot databases."
  - "Search databases should always include possible contaminants."
  - "For analyzing cell culture or organic samples, search databases should include mycoplasma databases."
  - "Some peptide search engines depend on decoys to calculate the FDR."
contributors:
  - stortebecker
  - bgruening
---

# Introduction
{:.no_toc}

In mass spectrometry based proteomics experiments, peptides are assigned to experimentally acquired tandem mass spectra (MS2) by a method called peptide-spectral matching. Peptide spectral matching is commonly achieved by using search algorithms to match the acquired MS2 spectra to theoretical spectra. The theoretical spectra are generated from an in silico digestion and fragmentation of proteins in the FASTA database. Ideally, the protein FASTA databases will contain all proteins of the organism under investigation.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Uploading a protein search Database

There are a many ways how you can upload your protein search database (FASTA file with protein sequences). Three of these ways are:

*   Using **Protein Database Downloader** {% icon tool %} .
*   Using a direct weblink to the database.
*   Uploading a database from the data library.

In this tutorial, we will explore **Protein Database Downloader** {% icon tool %} for generating a protein search database. For this we will download the proteome of an organism of interest. In this tutorial, we will use a database of the human proteome.

> ### {% icon hands_on %} Hands-on: Uploading a protein search database
>
> 1. Create a new history for this Database Handling exercise.
>
>    {% include snippets/create_new_history.md %}
>
> 2. Open **Protein Database Downloader** {% icon tool %}
> 3. Select in the drop-down menues `Taxonomy`: "Homo sapiens (Human)" and `reviewed`: "UniprotKB/Swiss-Prot (reviewed only)".
> 4. Click on `Execute`. There will be a new dataset named `Protein database` in your history, now.
> 5. Rename the `Protein database` to `Main database`.
>
>  > ### {% icon comment %} Types of uniprot databases
>  > Uniprot offers several types of databases. You may choose to download only reviewed (UniProtKB/Swissprot) databases, only unreviewed (UniProtKB/TREMBL) or both (UniProtKB). In model organisms which are well-researched, e.g. _Homo sapiens_ or _D. melanogaster_, reviewed (Swissprot) databases contain curated proteins and may lead to smaller databases and cleaner search results. However, if the researcher is interested in identifying proteins that are unreviewed, it might be wiser to include the TrEMBL database.
>  >
>  > You may also include protein isoforms by setting the tick box `Include isoform data` to `Yes`.
>  {: .comment}
>
>  > ### {% icon question %} Question
>  > What is the difference between a "reference proteome set" and a "complete proteome set"?
>  >
>  > > ### {% icon solution %} Solution
>  > > 1.  A UniProt complete proteome consists of the set of proteins thought to be expressed by an organism whose genome has been completely sequenced. A reference proteome is the complete proteome of a representative, well-studied model organism or an organism of interest for biomedical and biotechnological research. Reference proteomes constitute a representative cross-section of the taxonomic diversity to be found within UniProtKB. Species of particular importance may be represented by numerous reference proteomes for specific ecotypes or strains of interest. [Link to source](http://www.uniprot.org/help/reference_proteome)
>  > {: .solution }
>  {: .question}
{: .hands_on}


# Contaminant databases

In proteomic samples, some protein contaminants are commonly present. These protein contaminants are introduced into the sample during sample preparation, either from contaminated samples (e.g. mycoplasma in cell culture), chemicals or the experimenter in person. In order to avoid misidentification of spectra derived from contaminants, protein sequences of common laboratory contaminants are added to the database. This has two benefits:
1. The degree of contamination can be observed, heavily contaminated samples can be excluded from analysis.
2. Contaminant peptides cannot be misassigned to similar peptides in the database reducing the risk of identifying false positives.

A widely used database for common contaminants is the **c**ommon **R**epository of **A**dventitious **P**roteins (cRAP). When using samples generated in cell cultures, it is furthermore recommended to include _Mycoplasma_ proteomes in the search database. _Mycoplasma_ infections are very common in cell culture and often go unnoticed ([Drexler and Uphoff, Cytotechnology, 2002](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3463982/)).

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
>  > > ### {% icon solution %} Solution
>  > > 1. In samples stemming from a human source, identified human contaminants do not necessarily mean a contaminated sample. The proteins may as well originate from the research study sample. Users are advised to use discretion when interpreting the data.
>  > > 2. In samples from non-human sources, identified human contaminants do mean contamination by the experimenter.
>  > {: .solution }
>  {: .question}
{: .hands_on}


> ### {% icon hands_on %} Optional Hands-On: _Mycoplasma_ databases
> 90 - 95 % of mycoplasma infection in cell culture originates from the following species: _M. orale, M. hyorhinis, M. arginini, M. fermentans, M. hominis_ and _A. laidlawii_ ([Drexler and Uphoff, Cytotechnology, 2002](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3463982/)).
>
> 1. Use **Protein Database Downloader** {% icon tool %} to download the six mycoplasma databases. We will merge them to the main database in the next part of the tutorial.
> 2. Run **FASTA Merge Files and Filter Unique Sequences** {% icon tool %} to combine all mycoplasma databases into a single one.
> 3. Tag each entry in the combined database with the string "MYCOPLASMA_CONTAMINANT" by using **FASTA-to-Tabular** {% icon tool %}, **Add column** {% icon tool %} and **Tabular-to-FASTA** {% icon tool %}, as explained [above](#contaminant-databases).
> 4. Rename the **Tabular-to-FASTA** {% icon tool %} output to "Tagged Mycoplasma database".
>
>  > ### {% icon comment %} Comment
>  > The reviewed mycoplasma databases do not contain all known proteins. It is better to also include the TREMBL database. _Mycoplasma_ proteomes are relatively small, so even downloading TrEMBL sequences will not incraese the size of your main database by much.
>  {: .comment}
{: .hands_on}


# Merging databases

Depending on the search algorithm in use,  you might need to merge all FASTA entries (i.e. proteins of interest and contaminants) in a single database. Make sure to merge the tagged versions of your contaminant databases.

> ### {% icon hands_on %} Hands-on: Merging databases
>
> 1. Run **FASTA Merge Files and Filter Unique Sequences** {% icon tool %} on the main database and the tagged cRAP database.
    - Set `How are sequences judged to be unique?` to 'Accession Only'.
>
> 2. **Optional**: Merging mycoplasma databases
>
>    At this step you may also merge the mycoplasma protein databases that you downloaded earlier on. Simply enter them as additional inputs in **FASTA Merge Files and Filter Unique Sequences** {% icon tool %}. You can enter any number of protein databases when you click on `Insert Input FASTA file(s)`.
{: .hands_on}


# Creating a Decoy Database

The most common method of peptide and protein False Discovery Rate (FDR) calculation is by adding protein sequences that are not expected to be present in the sample. These are also called decoy protein sequences. This can be done by generating reverse sequences of the target protein entries and appending these protein entries to the protein database. Some search algoritmms use premade target-decoy protein sequences while others can generate a target-decoy protein sequence database from a target protein sequence database before using them for peptide spectral matching.

> ### {% icon hands_on %} Hands-on: Creating a Decoy Database
> 1. Run **DecoyDatabase**  {% icon tool %} on the merged database.
>
>  > ### {% icon comment %} Decoy tags
>  > The string you enter as a decoy tag will be added as a prefix or suffix (your choice) to the description of each decoy protein entry. Thus you can see from which entry in the target database the decoy was computed.
>  {: .comment}
>
>  > ### {% icon comment %} Comment
>  > **DecoyDatabase**  {% icon tool %} may also take several databases as input which are then automatically merged into one database.
>  {: .comment}
{: .hands_on}


# Concluding remarks
{:.no_toc}

In order to keep your protein databases up-to-date, it is recommended to create a workflow out of the hands-on sections (to learn about workflows see [this tutorial]({{site.baseurl}}/topics/introduction/tutorials/galaxy-intro-101/tutorial.html)). You might also want to combine the mycoplasma databases to a single file, which you then easily can add to each of your main databases.

Often you may not want to use the most recent database for reasons of reproducibility. If so, you can transfer the final database of this tutorial into other histories to work with it.

Further reading about construction of the optimal database: ([Kumar et al., Methods in molecular biology, 2017](https://www.ncbi.nlm.nih.gov/pubmed/27975281)).

This tutorial is based upon parts of the GalaxyP-101 tutorial (https://usegalaxyp.readthedocs.io/en/latest/sections/galaxyp_101.html).
