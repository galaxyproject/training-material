# Introduction

In metagenomics, information about micro-organisms in an environment can be extracted with two main techniques:

- Amplicon sequencing (or 16S rRNA/rDNA), which sequence only on the rRNA/rDNA of organisms
- Whole-genome sequencing (WGS), which sequence full genomes of the micro-organisms in the environment

In this tutorial, we will introduce the two types of analyses with the general principles behind and the differences. To go deeper in such analyses, we recommend to check our detailed tutorials on each analysis.

For that, we will use two datasets (one amplicon and one WGS) from the same environment: the Argentina Anguil Bulk Soil, studied in a [project on the Argentinean agricultural pampean soils](https://www.ebi.ac.uk/metagenomics/projects/SRP016633). In this project, three different types of land uses and two soil types (bulk and rhizospheric) were analyzed using WGS and amplicon sequencing. We will focus on the Argentina Anguil Bulk Soil.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. [Amplicon data](#amplicon-data)
> 2. [Whole-genome sequencing data](#whole-genome-sequencing-data)
> {: .agenda}

# Amplicon data

Amplicon sequencing is a highly targeted approach for analyzing genetic variation in specific genomic regions. In the metagenomics fields, amplicon sequencing refers to capture and sequence of rRNA data in a sample. It can be 16S for bacteria or archea or 18S for eukaryotes. 

> ### :book: Background: The 16S ribosomal RNA gene
> ![](../../images/16S_gene.png) <br><br>
>
> The 16S rRNA gene has several properties that make it ideally suited for our purposes
>
> 1. Present in all living organisms
> 2. Single copy (no recombination)
> 3. Highly conserved + highly variable regions
> 4. Huge reference databases
>
> ![](../../images/16S_variableregions.jpg)
>
> The highly conserved regions make it easy to target the gene across different organisms,
> while the highly variable regions allow us to distinguish between different species.
>
> (slide credit [http://www.slideshare.net/beiko/ccbc-tutorial-beiko ](http://www.slideshare.net/beiko/ccbc-tutorial-beiko ))
{: .tip}

With amplicon data, we can extract from which micro-organisms the sequences in our sample are coming from. This is called taxonomic assignation. We try to assign sequences to taxons and then classify or extract the taxonomy in our sample.

Our dataset comes from a sample of Anguil soil with capture and sequencing of the 16S rDNA V4 region using 454 GS FLX Titanium. The original data are available at EBI Metagenomics under run number [SRR651839](https://www.ebi.ac.uk/metagenomics/projects/SRP016633/samples/SRS386929/runs/SRR651839/results/versions/2.0).

In this analysis, we will use [Mothur tool suite](http://mothur.org), but only a small portion of its tools and possibilities. To learn more in detail how to use, we recommend to check our [Mothur tutorial](). 

## Importing the data

> ### :pencil2: Hands-on: Data upload
>
> 1. Import the FASTQ file from [Zenodo]() or from the data library (in "Analyses of metagenomics data" the "..." file)
>
>    > ### :bulb: Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * Press **Start**
>    {: .tip}
>
>    > ### :bulb: Tip: Importing data from a data library
>    >
>    > * Go into "Shared data" (top panel) then "Data libraries"
>    > * Click on "Training data" and then "Analyses of metagenomics data"
>    > * Select interesting file
>    > * Click on "Import selected datasets into history"
>    > * Import in a new history
>    {: .tip}
>
>    As default, Galaxy takes the link as name, so rename them.
>
{: .hands_on}

## Extraction of taxonomic information

The main questions when analyzing amplicon data are: Which micro-organisms are present in an environmental samples? And in which proportion? What is the structure of the community of the micro-organisms?

The idea is to take the sequences and assign them to a taxon. To do that, we group (or cluster) sequences based on their similarity to define Operational Taxonomic Units, groups of similar sequences that can be treated as a single "genus" or "species" (depending on the clustering threshold)

![](../../images/otu.png)

> ### :pencil2: Hands-on: Extraction of OTUs with Mothur
>
> 1. Step1
> 2. Step2
>
>    > ### :nut_and_bolt: Comments
>    > A comment
>    {: .comment}
>
>    > ### :bulb: Tip: A tip
>    >
>    > * Step1
>    > * Step2
>    {: .tip}
{: .hands_on}

Once the sequences are clustered into OTUs, one sequence of each OTU is selected as a representative sequence for the OTU. The taxonomic assignation (genus, species, ...) for this sequence is searched and then assigned to all sequences of the OTU. 

> ### :pencil2: Hands-on: Taxonomic assignation of the OTUs
>
> 1. Step1
> 2. Step2
>
>    > ### :nut_and_bolt: Comments
>    > A comment
>    {: .comment}
>
>    > ### :bulb: Tip: A tip
>    >
>    > * Step1
>    > * Step2
>    {: .tip}
{: .hands_on}

With the taxonomic assignation for each OTU, we can now extract for each genus (or other taxonomic level) how many OTUs (with how many sequences) are assigned to this genus (or other taxonomic level): extracting the community structure (taxon and their abundance) for the sample. 

To explore the community structure, we can visualize it with dedicated tools such as Phinch:

> ### :pencil2: Hands-on: Visualization of the community structure with Phinch
>
> 1. Step1
> 2. Step2
>
>    > ### :nut_and_bolt: Comments
>    > A comment
>    {: .comment}
>
>    > ### :bulb: Tip: A tip
>    >
>    > * Step1
>    > * Step2
>    {: .tip}
{: .hands_on}

Once we have information about the community structure (OTUs with taxonomic structure), we can do more analysis on it: estimation of the diversity of micro-organism, comparison fo diversity between samples, analysis of populations, ... We will not detail such analyses here but you follow our tutorials on amplicon data analyses to learn about them.

# Whole-genome sequencing data

Difference with amplicon: whole genomes (not only 16S or 18S), more information

Why amplicon is still used?

## Extraction of taxonomic information

Same question as for the amplicon data

### Sequence sorting

After we can use similar approaches as for amplicon data

To show the few number of rRNA sequences in a WGS datasets, impact on the results

### Taxonomic analysis

Difference with OTUs?

> ### :pencil2: Hands-on: Data upload
>
> 1. Step1
> 2. Step2
>
>    > ### :question: Question
>    >
>    > Question?
>    >
>    > <details>
>    > <summary>Click to view answers</summary>
>    > Answer to question
>    > </details>
>    {: .question}
{: .hands_on}

And visualization?

## Extraction of functional information 

Question: Which functions are done by the micro-organisms in the environment?

> ### :pencil2: Hands-on: Data upload
>
> 1. Step1
> 2. Step2
>
>    > ### :question: Questions
>    >
>    > 1. Question1?
>    > 2. Question2?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>Answer for question1</li>
>    >    <li>Answer for question2</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 3. Step3
{: .hands_on}

To go further on these analyses, you can follow our tutorials on whole-genome sequencing data analyses.

# Conclusion

Conclusion about the technical key points. And then relation between the technics and the biological question to end with a global view.
