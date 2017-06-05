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

In the previous section, we see how to analyze amplicon data to extract the community structure. Such information can also be extracted from whole-genome sequencing (WGS) metagenomic data. 

In WGS data, full genomes of the micro-organisms in the environment are sequenced (not only the 16S or 18S). We can then have access to the rRNA (only a small part of the genomes), but also to the genes of the micro-organisms. Using this information, we can try to answer to questions "What are the micro-organisms doing?" in addition to the question "What micro-organisms are present?".

## Extraction of taxonomic information

As for amplicon data, we can extract taxonomic and community structure information from WGS data. Different approaches can be used:

- Same approaches as for amplicon data with identification and classification of OTUs

    Such approaches imply a first step of sequence sorting to extract only the 16S and 18S sequences on which the same tools as for amplicon data. However, rRNA sequences represent a low proportion (< 1%) of the WGS sequences so such approache is not the most statistically supported

- Assignation of taxonomy on the whole sequences using databases with marker genes

In this tutorial, we use the second approach with MetaPhlAn2

> ### :pencil2: Hands-on: Taxonomic assignation with MetaPhlAn2
>
> 1. **MetaPhlAN2** :wrench:: Run **MetaPhlAN2** on the dereplicated sequences
>
>    > ### :question: Questions
>    >
>    > 1. What does the main output file contain?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li></li>
>    >    <li></li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
{: .hands_on}

Even if the output of MetaPhlAn2 is bit easier to parse than the BIOM file, we want to visualize and explore the community structure. We use an interactive tool called KRONA

> ### :pencil2: Hands-on: Interactive visualization with KRONA
>
> 1. **Format MetaPhlAn2 output for Krona** :wrench:: Run **Format MetaPhlAn2 output for Krona** to format MetaPhlAn2 output for KRONA
> 2. **KRONA** :wrench:: Run **KRONA** on the formatted MetaPhlAn2 output
>
>    > ### :question: Questions
>    >
>    > 1. What are the main species found for the bacteria?
>    > 2. Is the proportion of ... similar to the one found with EBI Metagenomics?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li></li>
>    >    <li></li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
{: .hands_on} 

*One sentence to conclude the taxonomic analysis*

## Extraction of functional information 

We would like now to answer the question "What are the micro-organisms doing?" or "Which functions are done by the micro-organisms in the environment?". 

In the WGS data, we have access to the gene sequences. We use that to identify the genes, associate them to a function, build pathways, etc to investigate the functional part of the community.

> ### :pencil2: Hands-on: Metabolism function identification
>
> 1. **HUMAnN2** :wrench:: Run **HUMAnN2** on non rRNA sequences (SortMeRNA output) to extract the gene families and pathways in the sample
> 
>    > ### :question: Questions
>    >
>    > 1. Which gene families is the most found one?
>    > 2. Which pathway is the most found one?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li></li>
>    >    <li></li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 2. Inspection of HUMAnN2 results?? Extraction of interesting info?
{: .hands_on}

*One sentence to conclude the functional analysis*

With the previous analyses, we investigate "Which micro-organims are present in my sample?" and "What function are done by the micro-organisms in my sample?". We can go further in these analyses (for example with combination of functional and taxonomic results). We did not detail that in this tutorial but you can found more analyses in our tutorials on whole-genome sequencing data analyses.

# Conclusion

*scheme to sum up the analyses*

*Conclusion about the technical key points. And then relation between the technics and the biological question to end with a global view.*
