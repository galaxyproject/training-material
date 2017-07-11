gi---
layout: tutorial_hands_on
topic_name: metagenomics
tutorial_name: shotgun
---

# Introduction
{:.no_toc}

In metagenomics, information about micro-organisms in an environment can be extracted with two main techniques:

- Amplicon sequencing, which sequence only on the rRNA/rDNA of organisms
- Shotgun sequencing, which sequence full genomes of the micro-organisms in the environment

Data generated from these two techniques must be treated differently. In this tutorial, we will focus on the analysis of whole-genome sequencing. 

> ### :nut_and_bolt: Comments
> If you want to learn how to analyze amplicon data, please check our dedicated tutorials
{: .comment}

From both amplicon and shotgun metagenomics raw data, we can extract information about which micro-organisms are present in the studied environment. But, contrary to amplicon, shotgun metagenomics data contain also full genome information about the micro-organisms. It is then possible to identify genes associated to functions, to reconstruct metabolic pathways and then determine which functions are done by the micro-organisms in the studied environment. It is even possible to go further and determine which micro-organisms are involved in a given function or pathways.

> ### Agenda
>
> However, extraction of useful information from raw shotgun metagenomics sequences is a complex process
with numerous bioinformatics steps and tools to use. These steps can be get together in 4 main steps we will deal with in the following tutorial:
>
> 1. TOC
> {:toc}
>
{: .agenda}


In this tutorial, we will work on a sample from Arctic Ocean (at 451 m), sequenced with Illumina MiSeq. This sample has already been analyzed with the [EBI Metagenomics' pipeline](https://www.ebi.ac.uk/metagenomics/pipelines/3.0), which uses slightly different tools than the ones we will use here. But, we could compare our results with the ones obtained with EBI Metagenomics.

# Pretreatments

Before any extraction of information about the community, raw sequences have to be pre-processed with quality control of the raw sequences and sequence sorting. But, first, we need in get our data in Galaxy.

## Data upload

The original data are available at EBI Metagenomics under run number [ERR1855251](https://www.ebi.ac.uk/metagenomics/projects/ERP015773/samples/ERS1569001/runs/ERR1855251/results/versions/3.0).

> ### :pencil2: Hands-on: Data upload
>
> 1. Import the FASTQ file pair from [Zenodo]() or from the data library
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
>    > * Click on "Training data" and then "WGS input data"
>    > * Select both files
>    > * Click on "Import selected datasets into history"
>    > * Import in a new history
>    {: .tip}
>
>    As default, Galaxy takes the link as name, so rename them.
>
{: .hands_on}

## Quality control and treatment

For quality control, we use similar tools as described in [the Quality Control tutorial](../../NGS-QC/tutorials/dive_into_qc): [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).

> ### :pencil2: Hands-on: Quality control
>
> 1. **FastQC** :wrench:: on both FastQ files to control the quality of the reads
> 2. **MulitQC** :wrench:: with
>    - "Software name" to `FastQC`
>    - "Result file" to the raw data generated with FastQC
>
>    > ### :question: Questions
>    >
>    > 1. What can we say about the quality of the sequences in both sequence files?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>The quality of the sequence decrease a lot at the end of sequences for both datasets. We will need to trim them</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 2. **Trim Galore** :wrench:: with
>    - "Is this library paired- or single-end?" to `Paired-end`
>    - "Reads in FASTQ format" to the input datasets with first the forward (ending with `_1`) and then the reverse (ending with `_2`)
>    - "Trim Galore! advanced settings" to `Full parameter list`
>    - "Trim low-quality ends from reads in addition to adapter removal" to `20`
>    - "Discard reads that became shorter than length N" to `60`
>    - "Generate a report file" to `Yes`
>
> 3. **MulitQC** :wrench:: with
>    - "Software name" to `Cutadapt`
>    - "Result file" to the report file generated with Trim Galore!
>
>    > ### :question: Questions
>    >
>    > 1. How much of the sequences have been trimmed?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>6.7% and 18.3% of the bases for the forward reads and reverse reads (respectively)</li>
>    >    </ol>
>    >    </details>
>    {: .question}
{: .hands_on}

One sequence file in Fasta is expected for the next steps. We need then to assemble the paired sequences and convert them to Fasta.

> ### :pencil2: Hands-on: Dereplication
>
> 1. **fastq-join** :wrench:: with
>    - "Dataset type" to `Paired-end`
>    - "Read 1 Fastq" to the forward trimmed reads 
>    - "Read 2 Fastq" to the reverse trimmed reads
> 2. **FASTQ to FASTA** :wrench:: with
>    - "FASTQ Library to convert" to the joined FastQ file
>    - "Discard sequences with unknown (N) bases" to `no`
>    - "Rename sequence names in output file" to `no`
{: .hands_on}

## Dereplication

During sequencing, one sequence must have been added to the dataset in multiple exact copy. Removing such duplicates reduce the size of the dataset without loosing information, with the dereplication (identification of unique sequences in a dataset).

> ### :pencil2: Hands-on: Dereplication
>
> 1. **VSearch dereplication** :wrench:: with
>    - "Select your FASTA file" to the Fasta file generated previously
>    - "Strand specific clustering" to `Both strand`
>
>    > ### :question: Questions
>    >
>    > 1. How many sequences are removed with the dereplication?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>42,363 on 43,145</li>
>    >    <li></li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
{: .hands_on}

## Sequence sorting

With shotgun metagenomics data, full genome information can be accessed: information corresponding to CDS of the micro-organisms, sequences corresponding to ribosomal sequences (rDNA or rRNA) of the micro-organisms, ... Useful functional information are present in sequences corresponding to CDS, and some taxonomic information in sequences corresponding to ribosomomal sequences (like the amplicon). To reduce the dataset size for the extraction of functional information, we can remove rRNA/rDNA sequences from the original dataset. 

This task is also useful to inspect the rRNA/rDNA sequences. And as in EBI Metagenomics' pipeline, these sequences can be used for taxonomic analyses as any amplicon data

> ### :nut_and_bolt: Comments
> If you want to learn how to analyze amplicon data, please check our dedicated tutorials
{: .comment}

For this task, we use SortMeRNA ([Kopylova et al, 2012](https://academic.oup.com/bioinformatics/article-abstract/28/24/3211/246053)). This tool filter RNA sequences based on local sequence alignment (BLAST) against 8 rRNA databases (2 Rfam databases for 5.8S and 5S eukarya sequences and 6 SILVA datasets for 16S (archea and bacteria), 18S (eukarya), 23S (archea and bacteria) and 28S (eukarya) sequences.

> ### :pencil2: Hands-on: Sequence sorting
>
> 1. **SortMeRNA** :wrench:: with
>    - "Querying sequences" to the dereplicated dataset
>    - "Sequencing type" to `Reads are not paired`
>    - "Which strands to search" to `Search both strands`
>    - "Databases to query" to `Public pre-indexed ribosomal databases`
>    - "rRNA databases" to select all
>    - "Include aligned reads in FASTA/FASTQ format?" to `Yes`
>    - "Include rejected reads file?" to `Yes`
>    - "Include alignments in SAM format?" to `No`
>    - "Include alignments in BLAST-like format?" to `No`
>
>    > ### :question: Questions
>    >
>    > 1. Which percentage of the original data are assigned to rRNA/rDNA sequences?
>    > 2. How can you explain with low percentage?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>1,988 over 42,363 are aligned on rRNA databases so 4.7%</li>
>    >    <li>Shotgun metagenomics data with few rRNA genes</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 2. (Optional) **SortMeRNA** :wrench:: Run SortMeRNA on the assigned rRNA sequences with the selection of 16S rRNA databases
> 3. (Optional) **SortMeRNA** :wrench:: Run SortMeRNA on the assigned rRNA sequences with the selection of 18S rRNA databases
{: .hands_on}

# Extraction of taxonomic information

The first important information to extract from any metagenomics sample is which micro-organisms are present in the sequenced sample and in which proportion are they present. So we want to extract information about the structure of the community of micro-organisms in the environment.

To identify the community structure, several approaches can be used. With amplicon or rRNA data, the sequences are clustered into Operational Taxonomic Units (OTU) and one representative sequence of each OTU is assigned to the most plausible microbial lineage. This approach is possible because of rRNA data: data that evolved quite slowly compared to other part of genomes, rRNA sequences data (particularly 16S and 18S) are then well conserved and good taxonomic markers.

However, for shotgun data, applying such approaches implies using a really small proportion of sequences for the taxonomic assignation, which can induce statistical bias. Other approaches have been developed to cope with shotgun data. For example, MetaPhlAn2 ([Truong et al, 2015](https://www.nature.com/nmeth/journal/v12/n10/full/nmeth.3589.html)) uses a database of ~1M unique clade-specific marker genes (not only the rRNA genes) identified from ~17,000 reference (bacterial, archeal, viral and eukaryotic) genomes.

## Taxonomic assignation

> ### :pencil2: Hands-on: Taxonomic assignation
>
> 1. **MetaPhlAN2** :wrench:: with
>    - "Input file" to the dereplicated sequences
>    - "Database with clade-specific marker genes" to `Locally cached`
>    - "Cached database with clade-specific marker genes" to `MetaPhlAn2 clade-specific marker genes`
>    - "Type of analysis to perform" to `Profiling a metagenomes in terms of relative abundances`
>    - "Taxonomic level for the relative abundance output" to `All taxonomic levels`
>
{: .hands_on}

3 files are generated by MetaPhlAn2:

- A tabular file with the community structure

    ```
    #SampleID   Metaphlan2_Analysis
    k__Bacteria 100.0
    k__Bacteria|p__Proteobacteria   86.20712
    k__Bacteria|p__Actinobacteria   13.79288
    k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria    86.20712
    k__Bacteria|p__Actinobacteria|c__Actinobacteria 13.79288
    ```

    Each line contains a taxa and its relative abundance found for our sample. The file starts with high level taxa (kingdom: `k__`) and go to more precise taxa.


- A BIOM file with the same information as the previous file but in BIOM format

    It can be used then by mothur and other tools requiring community structure information in BIOM format

- A SAM file with the results of the mapping of the sequences on the reference database

> ### :question: Questions
>
> 1. What is the most precise level we have access to with MetaPhlAn2?
> 2. What are the two orders found in our sample?
> 3. What is the most abundant family in our sample?
>
>    <details>
>    <summary>Click to view answers</summary>
>    <ol type="1">
>    <li>We have access to species level</li>
>    <li>Pseudomonadales and Solirubrobacterales are found in our sample</li>
>    <li>The most abundant family is Pseudomonadaceae with 86.21 % of the assigned sequences</li>
>    </ol>
>    </details>
{: .question}

## Community structure visualization

The generated files remains difficult to inspect. Visualization would help. Here, we will look at 3 tools that can be used to visualize the community structure in our sample.

Krona ([Ondov et al, 2011](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-385)) is a visualization tool for intuitive exploration of relative abundances of taxonomic classifications. It produces an interactive HTML file

> ### :pencil2: Hands-on: Interactive visualization with KRONA
>
> 1. **Format MetaPhlAn2 output for Krona** :wrench:: with
>    - "Input file" to `Community profile` output of `MetaPhlAn2`
>
> 2. **KRONA pie chart** :wrench:: with
>    - "What is the type of your input data" as `MetaPhlan`
>    - "Input file" to the output of `Format MetaPhlAn2`
{: .hands_on}

> ### :question: Questions
>
> 1. What are the main species found for the bacteria?
>
>    <details>
>    <summary>Click to view answers</summary>
>    <ol type="1">
>    <li></li>
>    </ol>
>    </details>
{: .question}

[GraPhlAn](https://bitbucket.org/nsegata/graphlan/wiki/Home>) is a tool for producing circular static representation of taxonomic analyses, easily exportable.

> ### :pencil2: Hands-on: Static visualization with GraPhlAn
>
> 1. **Export to GraPhlAn** :wrench:: with
>    - "Input file" to `Community profile` output of `MetaPhlAn2`
>    - "Use a LEfSe output file as input?" to `No`
>    - "List which levels should be annotated in the tree" to `5`
>    - "List which levels should use the external legend for the annotation" to `6,7`
>    - "Title font size" to `15`
>    - "Default size for clades not found as biomarkers" to `10`
>    - "Minimum value of clades that are biomarkers" to `0`
>    - "Maximum value of clades that are biomarkers" to `250`
>    - "Font size for the annotation legend" to `11`
>    - "Minimun abundance value for a clade to be annotated" to `0`
>    - "Number of clades to highlight" to `100`
>
>    > ### :nut_and_bolt: Comments
>    > We decide to display the maximum of clade (100, here). If you want more or less, you can modulate the number of clades to highlight. And if you want to change displayed annotations, you can change levels to annotate.
>    {: .comment}
>
> 2. **Generation, personalization and annotation of tree for GraPhlAn** :wrench:: with
>    - "Input tree" to the tree generated by the previous tool
>    - "Annotation file" to the annotation file generated by the previous tool
> 3. **GraPhlAn** :wrench:: with
>    - "Input tree" to the tree generated in the last step
>    - "Output format" to `PNG`
>
{: .hands_on}

With our dataset, we obtain a nice graphical representation of taxonomic diversity inside our sample, with circle radius being proportional to relative abundance of the corresponding clade.

> ### :question: Questions
>
> 1. What are the main species found for the bacteria?
> 2. Is the main species the same as the one observed with KRONA? 
>
>    <details>
>    <summary>Click to view answers</summary>
>    <ol type="1">
>    <li></li>
>    </ol>
>    </details>
>
{: .question}

[Phinch](http://phinch.org/) is an open-source framework for visualizing biological data. 

Galaxy runs an instance of Phinch which is directly accessible for any BIOM file.

> ### :pencil2: Hands-on: Static visualization with GraPhlAn
>
> 1. Click on the BIOM file generated by MetaPhlAn2
> 2. Click on the "Visualize" icon and then on "Phinch"
>   
>    It will lead you to the Phinch website, which will automatically load in your file, and where you can several interactive visualisations:
>
{: .hands_on}

# Functional analyses

Investigation of the structure composition gives an insight on "What organisms are present in our sample". We now want to know "What are they doing in that environment?". With shotgun data, we have full genome information with sequence of genes. To determine which functions are done by the micro-organisms in the studied environment, we need then to identify genes, associate them to functions, combine such information to reconstruct metabolic pathways, ...

The first step is then to identify sequences and affiliate them to a known genes, using available database. [HUMAnN2](http://huttenhower.sph.harvard.edu/humann2) is a tool to profile the presence/absence and abundance of gene families and microbial pathways in a community from metagenomic or metatranscriptomic sequencing data.

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
{: .hands_on}

The 3 generated files give detailed insights into gene families and pathways. This is interesting when we want to look to a particular pathway or to check abundance of a given gene families. However, when we want a broad overview of metabolic processes in a community, we need tools to regroup gene families or pathways into global categories.

To get global categories from HUMAnN2 outputs, we decide to use the [Gene Ontology](http://geneontology.org/>) to describe the gene families in terms of their associated biological processes, cellular components and molecular functions.

> ### :pencil2: Hands-on: Metabolism function identification
>
> 1. **Group humann2 uniref50 abundances to Gene Ontology (GO) slim terms** :wrench:: Run **Group humann2 uniref50 abundances to Gene Ontology (GO) slim terms** on the gene family abundance generated with HUMAnN2
> 
>    > ### :question: Questions
>    >
>    > 1. Which molecular function is the most abundant one?
>    > 2. Which biological process is the most abundant one?
>    > 3. Which cellular component is the most abundant one?
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
> 2. Generate a barplot for each of the generated output by clicking on the visualization button on each dataset
{: .hands_on}

# Combination of taxonomic and functional results

With the previous analyses, we can now give some answers to the questions "Which micro-organims are present in my sample?" and "What function are done by the micro-organisms in my sample?". One remaining question stays unanswered: "Which micro-organisms are implied in the realization of a given function?".

To answer this question, we need to relate generated taxonomic and functional results. 

> ### :pencil2: Hands-on: Combination of taxonomic and functional results
>
> 1. **Combine... ** :wrench:: Run **Combine** on the gene family abundance generated with HUMAnN2 and the MetaPhlAn2 output
> 
>    > ### :question: Questions
>    >
>    > 1. 
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
> 2. 
{: .hands_on}

# Conclusion
{:.no_toc}

Conclusion about the technical key points. And then relation between the technics and the biological question to end with a global view.
