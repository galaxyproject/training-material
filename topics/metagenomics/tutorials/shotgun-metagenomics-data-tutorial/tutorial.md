---
layout: tutorial_hands_on

title: "Analysis of shotgun metagenomics data"
zenodo_link: "https://doi.org/10.5281/zenodo.815875"
questions:
  - "How to analyze metagenomics shotgun data?"
  - "What information can be extracted of metagenomics data?"
  - "What are the difference in the analyses of amplicon and shotgun data?"
objectives:
  - "Selection of tools to analyze shotgun data"
  - "Extracting taxonomic and functional information"
  - "Visualisation of a community structure"
time_estimation: "1H"
key_points:
  - "With amplicon data, we can extract information about the studied community structure"
  - "With shotgun data, we can extract information about the studied community structure and also the functions realised by the community"
  - "The tools used to analyze amplicon and shotgun data are different, except for the visualisation"
  - "Metagenomics data analyses are complex and time-consuming"
contributors:
  - shiltemann
  - bebatut
  - EngyNasr
---

# Introduction
{:.no_toc}

In metagenomics, information about micro-organisms in an environment can be extracted with two main techniques:

- [Amplicon sequencing]({% link topics/metagenomics/tutorials/mothur-miseq-sop/tutorial.md %}), which sequences only the rRNA or ribosomal DNA of organisms
- __Shotgun sequencing__, which sequences full genomes of the micro-organisms in the environment

In this tutorial, we will introduce the second type of analysis, shotgun data, with its general principles. For a more in-depth look at these analyses, we recommend our detailed tutorials on each analysis, Amplicon and Shotgun.

We will use the dataset from [project on the Argentinean agricultural pampean soils](https://www.ebi.ac.uk/metagenomics/projects/SRP016633). In this project, three different geographic regions that are under different types of land uses and two soil types (bulk and rhizospheric) were analyzed using shotgun and amplicon sequencing. We will focus on data from the Argentina Anguil and Pampas Bulk Soil (the original study included one more geographical regions, [see](https://doi.org/10.1186/2049-2618-1-21)).

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet faqs/galaxy/analysis_results_may_vary.md %}

# Shotgun metagenomics data

In [16S Microbial Analysis with mothur tutorial]({% link topics/metagenomics/tutorials/mothur-miseq-sop/tutorial.md %}), we saw how to analyze amplicon data to extract the community structure. Such information can also be extracted from shotgun metagenomic data.


In shotgun data analysis, full genomes of the micro-organisms in the environment are sequenced (not only the 16S or 18S). We can then have access to the rRNA (only a small part of the genomes), but also to the other genes of the micro-organisms. Using this information, we can try to answer questions such as "What are the micro-organisms doing?" in addition to the question "What micro-organisms are present?".

{% snippet topics/metagenomics/faqs/sequencing_shotgun.md %}

We will use a metagenomic sample of the Pampas Soil ([SRR606451](https://www.ebi.ac.uk/metagenomics/projects/SRP016633/samples/SRS372043/runs/SRR606451/results/versions/2.0)).

## Data upload

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the `SRR606451_pampa` Fasta file from [Zenodo](http://zenodo.org/record/815875) or from the data library (in "Analyses of metagenomics data")
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
>    ```
>    https://zenodo.org/record/815875/files/SRR606451_pampa.fasta
>    ```
>
{: .hands_on}

## Quality control

During sequencing, errors are introduced, such as incorrect nucleotides being called. These are due to the technical limitations of each sequencing platform. Sequencing errors might bias the analysis and can lead to a misinterpretation of the data.

Sequence quality control is therefore an essential first step in your analysis.
In this tutorial we use similar tools as described in the tutorial ["Quality control"]({% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %}):
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) generates a web report that will aid you in assessing the quality of your data
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) for trimming and filtering

> ### {% icon hands_on %} Hands-on: Quality control
>
> 1. {% tool [FastQC](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.73+galaxy0) %} with the following parameters:
>       - {% icon param-files %} *"Raw read data from your current history"*: imported file: `SRR606451_pampa.fasta`
>
>
> 2. Inspect the webpage output of **FastQC** {% icon tool %} for the `SRR606451_pampa.fasta` dataset
>
>    > ### {% icon question %} Questions
>    >
>    > What is the read length?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > The read length is  bp.
>    > >
>    > {: .solution}
>    >
>    {: .question}
>
{: .hands_on}

For more information about how to interpret the plots generated by FastQC, please see [this section]({% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %}#assess-quality-with-fastqc) in our dedicated Quality Control Tutorial.

Even though our data is already of pretty high quality, we can improve it even more by:
1. Trimming reads to remove bases that were sequenced with low certainty (= low-quality bases) at the ends of the reads
2. Removing reads of overall bad quality.
3. Removing reads that are too short to be informative in downstream analysis

> ### {% icon question %} Questions
>
> What are the possible tools to perform such functions?
>
>   > ### {% icon solution %} Solution
>   >
>   > There are many tools such as **Cutadapt**, **Trimmomatic**, **Trim Galore**, **Clip**, **trim putative adapter** sequences. etc. We choose here **Cutadapt** because it is error tolerant, it is fast and the version is pretty stable.
>   {: .solution}
>
{: .question}

There are several tools out there that can perform these steps, but in this analysis we use **Cutadapt** ({% cite martin2011cutadapt %}).

**Cutadapt** also helps find and remove adapter sequences, primers, poly-A tails and/or other unwanted sequences from the input FASTQ files. It trims the input reads by finding the adapter or primer sequences in an error-tolerant way. Additional features include modifying and filtering reads.


> ### {% icon hands_on %} Hands-on: Read trimming and filtering
>
> 1. {% tool [Cutadapt](toolshed.g2.bx.psu.edu/repos/lparsons/cutadapt/cutadapt/3.4+galaxy2) %} with the following parameters to trim low quality sequences:
>    - *"Single-end or Paired-end reads?"*: `Single-end`
>       - {% icon param-files %} *"FASTQ/A file"*: imported file: `SRR606451_pampa.fasta`
>
>    - In *"Filter Options"*
>       - *"Minimum length"*: `150`
>    - In *"Read Modification Options"*
>       - *"Quality cutoff"*: `20`
>    - In *"Output selector"*
>       - *"Report"*: select `Report: Cutadapt's per-adapter statistics.`
>
>      {% include topics/sequence-analysis/tutorials/quality-control/trimming_question.md %}
>
> 2. **Rename** {% icon galaxy-pencil %}
>     - `Read 1 output` to `QC controlled reads`
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
{: .hands_on}

**Cutadapt** {% icon tool %} outputs a report file containing some information about the trimming and filtering it performed.

> ### {% icon question %} Questions
>
> Inspect the output report from **Cutadapt** {% icon tool %}.
>
> 1. What % reads contain adapter?
> 2. What % reads have been trimmed because of bad quality?
> 3. What % reads have been removed because they were too short?
>
>   > ### {% icon solution %} Solution
>   >
>   > 
>   {: .solution}
>
{: .question}

## Ribosomal RNA fragments filtering

The ribosomal RNA **rRNAs** are useful for the taxonomic assignment (i.e. which organisms are found) but they do not provide any functional information, (i.e. which genes are expressed) To make the downstream functional annotation faster, we will sort the rRNA sequences using **SortMeRNA** ({% cite kopylova2012sortmerna %}). It can handle large RNA databases and sort out all fragments matching to the database with high accuracy and specificity:

![SortMeRNA](../../images/metatranscriptomics/sortmerna.png)

> ### {% icon hands_on %} Hands-on: Ribosomal RNA fragments filtering
>
> 1. {% tool [Filter with SortMeRNA](toolshed.g2.bx.psu.edu/repos/rnateam/sortmerna/bg_sortmerna/2.1b.6) %} with the following parameters:
>    - *"Sequencing type"*: `Reads are not paired`
>       - {% icon param-file %} *"Querying sequences"*: `QC controlled reads` (output of **Cutadapt** {% icon tool %})
>    - *"Databases to query"*: `Public pre-indexed ribosomal databases`
>       - *"rRNA databases"*: {% icon param-check %} Select all
>          - {% icon param-check %} `rfam-5s-database-id98`
>          - {% icon param-check %} `silva-arc-23s-id98`
>          - {% icon param-check %} `silva-euk-28s-id98`
>          - {% icon param-check %} `silva-bac-23s-id98`
>          - {% icon param-check %} `silva-euk-18s-id95`
>          - {% icon param-check %} `silva-bac-16s-id90`
>          - {% icon param-check %} `rfam-5.8s-database-id98`
>          - {% icon param-check %} `silva-arc-16s-id95`
>    - *"Include aligned reads in FASTA/FASTQ format?"*: `Yes (--fastx)`
>       - *"Include rejected reads file?"*: `Yes`
>    - *"Generate statistics file"*: `Yes`
>
> 2. Expand the aligned and unaligned forward reads datasets in the history
>
>    > ### {% icon question %} Questions
>    >
>    > How many sequences have been identified as rRNA and non rRNA?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > 
>    > >
>    > {: .solution }
>    {: .question}
>
{: .hands_on }

**SortMeRNA** {% icon tool %} removes any reads identified as rRNA from our dataset, and outputs a log file with more information about this filtering.

## Extraction of taxonomic information

As for amplicon data, we can extract taxonomic and community structure information from shotgun data. Different approaches can be used:

- Same approach as for amplicon data with identification and classification of OTUs

    Such an approach requires a first step of sequence sorting to extract only the 16S and 18S sequences, then using the same tools as for amplicon data. However, rRNA sequences represent a low proportion (< 1%) of the shotgun sequences so such an approach is not the most statistically supported

- Assignation of taxonomy on the whole sequences using databases with marker genes

In this tutorial, we use the second approach with [MetaPhlAn](https://elifesciences.org/articles/65088). This tools is using a database of ~1.1M unique clade-specific marker genes (not only the rRNA genes) identified from ~100,000 reference (bacterial, archeal, viral and eukaryotic) genomes.

> ### {% icon hands_on %} Hands-on: Taxonomic assignation with MetaPhlAn
>
> 1. {% tool [MetaPhlAn](toolshed.g2.bx.psu.edu/repos/iuc/metaphlan/metaphlan/3.0.13+galaxy0) %} with
>    - {% icon param-file %} *"1: Single-end Fasta/FastQ file with microbiota reads"*: imported file: `SRR606451_pampa.fasta`
>    - In *"Inputs Options"*
>       - *"Database with clade-specific marker genes"*: `locally cached`
>       - *"Cached database with clade-specific marker genes"*:`MetaPhlAn clade-specific marker genes`
>
> This step may take a couple of minutes.
{: .hands_on}

3 files are generated:

- A tabular file with the community structure

    ```
    #SampleID   Metaphlan_Analysis
    k__Bacteria 100.0
    k__Bacteria|p__Proteobacteria   86.20712
    k__Bacteria|p__Actinobacteria   13.79288
    k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria    86.20712
    k__Bacteria|p__Actinobacteria|c__Actinobacteria 13.79288
    ```

    Each line contains a taxa and its relative abundance found for our sample. The file starts with high level taxa (kingdom: `k__`) and go to more precise taxa.

> ### {% icon question %} Questions
>
> 1. What are the different taxa levels ordered with their level of precision (from highest to lowest)?
>
> > ### {% icon solution %} Solution
> > species, genus, family, order, class, phylum, kingdom, domain
> {: .solution }
{: .question}    


- A BIOM file with the same information as the previous file but in BIOM format

    It can be used then by mothur and other tools requiring community structure information in BIOM format

- A SAM file with the results of the mapping of the sequences on the reference database

> ### {% icon question %} Questions
>
> 1. What is the most precise level we have access to with MetaPhlAn?
> 2. What are the two orders found in our sample?
> 3. What is the most abundant family in our sample?
>
> > ### {% icon solution %} Solution
> > 1. We have access to species level
> > 2. Pseudomonadales and Solirubrobacterales are found in our sample
> > 3. The most abundant family is Pseudomonadaceae with 86.21 % of the assigned sequences
> {: .solution }
{: .question}

Even if the output of MetaPhlAn is bit easier to parse than the BIOM file, we want to visualize and explore the community structure with KRONA

> ### {% icon hands_on %} Hands-on: Interactive visualization with KRONA
>
> 1. {% tool [Format MetaPhlAn2 output for Krona](toolshed.g2.bx.psu.edu/repos/iuc/metaphlan2krona/metaphlan2krona/2.6.0.0) %} with
>    - *"Input file"*: `Predicted taxon relative abundances` output of `MetaPhlAn`
>
> 2. {% tool [KRONA pie chart](toolshed.g2.bx.psu.edu/repos/crs4/taxonomy_krona_chart/taxonomy_krona_chart/2.7.1) %} with
>    - *"What is the type of your input data"* as `Tabular`
>    - *"Input file"*: the output of `Format MetaPhlAn`
>
{: .hands_on}

## Extraction of functional information

We would now like to answer the question "What are the micro-organisms doing?" or "Which functions are performed by the micro-organisms in the environment?".

In the shotgun data, we have access to the gene sequences from the full genome. We use that to identify the genes, associate them with a function, build pathways, etc., to investigate the functional part of the community.

> ### {% icon hands_on %} Hands-on: Metabolism function identification
>
> 1. {% tool [HUMAnN](toolshed.g2.bx.psu.edu/repos/iuc/humann/humann/3.0.0+galaxy1) %} with
>    - {% icon param-file %} *"1: Input sequence file: Pre-computed mappings of reads to database sequences"*: `imported file: SRR606451_pampa.fasta`
>    - {% icon param-file %} *"2: Taxonomic profile file"*: `Predicted taxon relative abundances` output of `MetaPhlAn`
>    - *"Input(s) options"*: `Pre-computed mappings of reads database sequences`
>    - In *"Steps options"*: `Bypass the taxonomic profiling step and creates a custom ChocoPhlAn database of the species provided afterwards`
>       - *"Nucleotide database"*: `Locally cached`
>          - *"Nucleotide database"*: `Full ChocoPhlAn for HUManN`
>       - *"Protein database"*: `Locally cached`
>          - *"Protein database"*: `Full UniRef50 for HUManN`
>    - In *"Gene and pathway quantification options"*
>       - *"Database to use for pathway computations"*: `MetaCyc`
>    - In *"Outputs options"*
>       - *"Remove stratification from output"*: `Yes`
>
>    This step is long so we generated the output for you!
>
> 2. Import the 3 files whose the name is starting with "humann"
>
>    ```
>    https://zenodo.org/record/815875/files/humann2_gene_families_abundance.tsv
>    https://zenodo.org/record/815875/files/humann2_pathways_abundance.tsv
>    https://zenodo.org/record/815875/files/humann2_pathways_coverage.tsv
>    ```
{: .hands_on}

HUMAnN generates 3 files

- A file with the abundance of gene families

    Gene family abundance is reported in RPK (reads per kilobase) units to normalize for gene length. It reflects the relative gene (or transcript) copy number in the community.

    The "UNMAPPED" value is the total number of reads which remain unmapped after both alignment steps (nucleotide and translated search). Since other gene features in the table are quantified in RPK units, "UNMAPPED" can be interpreted as a single unknown gene of length 1 kilobase recruiting all reads that failed to map to known sequences.

- A file with the coverage of pathways

    Pathway coverage provides an alternative description of the presence (1) and absence (0) of pathways in a community, independent of their quantitative abundance.

- A file with the abundance of pathways

> ### {% icon question %} Questions
>
> How many gene families and pathways have been identified?
>
> > ### {% icon solution %} Solution
> > 44 gene families but no pathways are identified
> {: .solution }
{: .question}

The RPK for the gene families are quite difficult to interpret in term of relative abundance. We decide then to normalize the values

> ### {% icon hands_on %} Hands-on: Normalize the gene family abundances
>
> 1. {% tool [Renormalize a HUMAnN generated table](toolshed.g2.bx.psu.edu/repos/iuc/humann2_renorm_table/humann2_renorm_table/0.11.1.1) %} with
>    - "Gene/pathway table" to the gene family table generated with `HUMAnN`
>    - "Normalization scheme" to `Relative abundance`
>    - "Normalization level" to `Normalization of all levels by community total`
>
>  > ### {% icon question %} Questions
>  >
>  > 1. What percentage of sequences has not be assigned to a gene family?
>  > 2. What is the most abundant gene family?
>  >
>  > > ### {% icon solution %} Solution
>  > > 1. 55% of the sequences have not be assigned to a gene family
>  > > 2. The most abundant gene family with 25% of sequences is a putative secreted protein
>  > {: .solution }
>  {: .question}
{: .hands_on}

With the previous analyses, we investigate "Which micro-organims are present in my sample?" and "What function are performed by the micro-organisms in my sample?". We can go further in these analyses (for example, with a combination of functional and taxonomic results). We did not detail that in this tutorial but you can find more analyses in our tutorials on shotgun metagenomic data analyses.

# Conclusion
{:.no_toc}

We can summarize the analyses with amplicon and shotgun metagenomic data:

![Scheme to sum up the analysis](../../images/general-tutorial-scheme.png)

Both analyses are quite complex! However, in this tutorial, we only showed simple cases of metagenomics data analysis with subset of real data.

Check our other tutorials to learn more in detail of how to analyze metagenomics data.
