---
layout: tutorial_hands_on
title: Taxonomic Assignment of Metagenomic Data
zenodo_link: still needs to be created
questions:
- Which species (or genera, families, ...) are present in my sample?

objectives:

- explain what taxonomic assignment is
- explain how taxonomic assignment works
- apply Kraken2 to assign taxonomic labels
- apply Krona to visualize results of assignment and understand the output
- identify taxonomic classification tool that fits best depending on their data
level: Introductory
key_points:
- To do
time_estimation: 45M
contributors:
- Sophia120199
---
> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}
# Introduction
 
The term **"microbiome"** describes "a characteristic microbial community occupying a reasonably well-defined habitat which has distinct physio-chemical properties. The term thus not only refers to the microorganisms involved but also encompasses their theatre of activity" ({% cite whipps1988fungi %}). 

Microbiome data can be gathered from different environments such as soil, water or the human gut. The biological interest lies in general in the question how the microbiome present at a specific site influences this environment. **Metagenomic samples** contain DNA from different organisms at a specific site, where the sample was collected. Metagenomic data can be used to find out which organisms coexist in that niche and which genes are present in the different organisms. 

Metatranscriptomic samples include the transcribed gene products, thus RNA, that therefore allow to not only study the presence of genes but additionally their expression in the given environment. The following tutorial will focus on metagenomics data, but the principle is the same for metatranscriptomics data.

The investigation of microbes present at a specific site and their relative abundance is also called **“microbial community profiling”**.
Basic for this is to find out which microbes are present in the sample. This can be achieved for all known microbes, where the DNA sequence specific for a certain species is known. The process of sorting single DNA reads derived from sequencing to a species (or other ranks) is called **taxonomic assignment**. 

When talking about taxonomic assignment or taxonomic classification, most of the time we actually talk about two methods, that in practice are often used interchangeably: while **taxonomic binning** is the classification of individual sequence reads to reference taxa, **taxonomic profiling** refers to investigating relative abundances of taxa within a dataset but not to the classification of individual reads.
Tools for taxonomic classification can be divided into three groups. Nevertheless, all of them require a pre-computed database based on previously sequenced microbial DNA or protein sequences.
1. **DNA-to-DNA** classification tools compare sequencing reads with genomic databases of DNA sequences (Bracken, Kraken, Kraken2, MegaBLAST)
2. **DNA-to-Protein** classification tools compare sequencing reads with genomic databases of protein sequences (more computationally intensive because of analysis of all six frames of potential DNA-to amino acid translation) (DIAMOND)
3. **Marker based** classification tools use a reference database that only includes a subset of gene sequences (e.g. 16S rRNA sequence), which is quick, but introduces bias (MetaPhlAn2)
 
> ### {% icon details %} More details on taxonomy
>
> In general, taxonomy is the study of sorting organisms into different groups within a larger system according to similarities and differences. The groups are named, defined, classified and hierarchically ordered. The principal ranks from top to bottom (*with examples for the human being*) are domain (*eukarya*), kingdom (*animalia*), phylum (*chordata*), class (*mammalia*), order (*primates*), family (*hominidae*), genus (*homo*), and species (*sapiens*). From this classification, one can generate a tree of life, also known as a phylogenetic tree. It is a rooted tree that describes the relationship of all life on earth. At the root sits the “last universal common ancestor” and the three main branches (in taxonomy also called domains) are bacteria, archaea and eukaryotes. Most important for this is the idea that all life on earth is derived from a common ancestor and therefore when comparing two species, you will -sooner or later- find a common ancestor for all of them.
> 
> {: .details}
 
When we talk about metagenomic data here, what we start with is sequences derived from DNA fragments that could be isolated from the sample of interest. Ideally, from all microbes present in the sample, we would also find DNA. The underlying idea of taxonomic assignment is to compare the DNA sequences found in the sample (reads) to DNA sequences of a database. When a read matches a database DNA sequence of a known microbe, we can derive a list with microbes present in the sample.
The comparison of reads to database sequences can be done in different ways, leading to three different types of taxonomic assignment: k-mer based, gene-based and genome-based analysis.
 
- For the **k-mer based** analysis, databases as well as the samples DNA are broken into k-mers about 30 bp length for comparison. From all the genomes in the database, where a specific k-mer is found, a lowest common ancestor (LCA) tree is derived and the abundance of k-mers within the tree is counted. This is the basis for a root-to-leaf path calculation, where the path with the highest score is used for classification of the sample. By counting the abundance of k-mers, also an estimation of relative abundance of taxa is possible. The major advantage of k-mer based analysis is the low compute cost. Major disadvantages are the low detection accuracy, that the unclassified percentage is unknown and that there is no gene detection, no SNVs detection and no genomic comparison possible. An example for a k-mer based analysis tool is Kraken2, which will be used in this tutorial
 
![Kraken functionality](../../taxonomic-assignment/images/Kraken_algorithm.png "Kraken functionality.")  {% cite Wood.2014 %}
 
- For the **gene based** analysis, reads are aligned to reference genes about 1 kbp length. Next, marker genes are used to estimate species abundance. Furthermore, genes can be analyzed in isolation for presence or absence in a specific condition.
The major advantage is the detection of the pangenome (entire set of genes within a species). Major disadvantages are the high compute cost, low detection accuracy and that the unclassified percentage is unknown. At least intragenic SNVs can be detected and low-resolution genomic comparison is possible.
 
- For the **genome based** analysis, read pairs of 150 bp length are aligned to reference genomes of about 3 Mbp length. Considering the coverage and breadth, genomes are used to measure genome abundance. Furthermore, genes can be analyzed in genomic context. Advantages of this method are the high detection accuracy, that the unclassified percentage is known, that all SNVs can be detected and that high-resolution genomic comparisons are possible. This method takes medium compute cost.
 
After this theoretical introduction, let's now get hands on analyzing an actual dataset! 


 
# Background on data
 
The dataset we will use for this tutorial comes from an oasis in the mexican desert called Cuatro Ciénegas, that is studied because of its special environmental conditions {% cite Okie.2020 %}. The researchers collected samples directly from the pond (= control mesocosm, called JC1A in this tutorial) and fertilized some of the samples later on to receive nutrient enrichment (=Lagunita Fertilized Pond, called JP4D in this tutorial). In this way, they investigated the impact of nutrient enchriment on the microbial community.
The datafiles are named according to the first four characters of the filenames.
It is a collection of paired-end data with R1 being the forward reads and R2 being the reverse reads. Additionally, the reads have been trimmed using [__cutadapt__ ](https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/quality-control/tutorial.html#trim-and-filter---short-reads)

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this exercise
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the FASTQ file pairs from [Zenodo]({{ page.zenodo_link }}) or a data library:
>    - `JP4D_R1.fastq.gz`
>    - `JP4D_R2.fastq.gz`
    -`JC1A_R1.fastq.gz`
    - `JC1A_R2.fastq.gz`
>    ``
>
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}

>
> 3. Create a paired collection.
>
>    {% snippet faqs/galaxy/collections_build_list_paired.md %}
>
{: .hands_on}


 
# Hands on: k-mer based taxonomic assignment of metagenomic data

Our input data is the DNA reads of microbes present at Cuatro Ciénegas. To find out which microbes are present, we will compare the reads of the sample to sequences of known microbes stored in a database using the algorithm Kraken2, which is k-mer based. For this tutorial, we will use the Standard plus protozoa & fungi (2021) database. The Standard database includes archaea, bacteria, viral, plasmid, human, and UniVec_Core data. 
 
> ### {% icon hands_on %} Hands-on: Assign taxonomic labels with Kraken2
>
> 1. {% tool [Kraken2](toolshed.g2.bx.psu.edu/repos/iuc/kraken2/kraken2/2.0.8_beta+galaxy0) %} with the following parameters:
>    - *"Single or paired reads"*: `Paired Collection`
>    - *"Print scientific names instead of just taxids"*: `No`
>    - *"Confidence"*: `0.1`
>    - In *"Create Report"*:
>        - *"Print a report with aggregrate counts/clade to file"*: `Yes`
>        - *"Format report output like Kraken 1's kraken-mpa-report"*: `Yes`
>    - *"Select a Kraken2 database"*: `Standard plus protozoa & fungi (2021)`
>
>    > ### {% icon comment %} Comment
>    >
>    > A confidence score of 0.1 means that at least 10% of the k-mers should match entries in the database. This value can be reduced if a less restrictive taxonomic assignation is desired.
>    {: .comment}
>
{: .hands_on}
 
 
Kraken2 will create two output files called "Classification" and "Report".
 
Let's have a look at the classification file. It has 5 columns:
 
1. C/U: classified/unclassified
2. Sequence ID
3. Taxonomy ID
4. Length of sequence (read1|read2 for paired reads)
5. indicates LCA mapping of each k-mer in the sequence |:| indicates end of first read, start of second read for paired reads --> example: "n k-mers assigned to taxon xxx"
 
 
![Kraken2 Classification Output](../../taxonomic-assignment/images/Kraken2_classification_screenshot.png "Kraken2 Classification Output.")

> ### {% icon question %} Questions
>Let's have a look at the first line of the Kraken2 classification output.
>1. Is it classified or unclassified?
>2. How long is the sequenece?
>3. What is the taxonomy ID?
> 
>
>
> > ### {% icon solution %} Solution
> >1. classified
> >2. 142 bp
> >3. 398580
> > 
>
> {: .solution}
>
{: .question}

 
 
Let's also have a look at the report file. It has 2 columns:
1. taxon name grouped into d_domain, p_phylum, c_class, o_order, f_family, g_genus, s_species
2. number of reads assigned to specific taxon
 
![Kraken2 Report Output](../../taxonomic-assignment/images/Kraken2_report_screenshot.png "Kraken2 Report Output.")

> ### {% icon question %} Questions
>
> 1. What family does Paracoccus sp. MC1862 belong to? 
> 2. How many reads were assigned to Paracoccus sp. MC1862
>
>
> > ### {% icon solution %} Solution
> >
> > 1. Rhodobacteraceae
> > 2. 917
>
> {: .solution}
>
{: .question}

 
As both files contain a lot of information, we will use __Krona__ {% cite Ondov.2011 %}to visualize the data.
 
# Visualization of taxonomic assignment
 
Once we have assigned the corresponding taxa to each sequence, the next step is to properly visualize the data, for which we will use the __Krona pie chart__ tool ({% cite Ondov.2011 %}). But first, we need to convert the output generated by Kraken2 so it can be used as an input from the Krona tool.
 
## Convert output from Kraken2 so it can be used for Krona
 
__Convert Kraken__ tool is designed to translate results of the Kraken metagenomic classifier (see citations below) to the full representation of NCBI taxonomy. It does so by using Taxonomic ID field provided by Kraken. The output of this tool can be directly visualized by the Krona tool.
 
> ### {% icon hands_on %} Hands-on: Convert Kraken2 Output
>
> 1. {% tool [Convert Kraken](xxx) %} with the following parameters:
>    - *"Choose dataset to convert"*: Datset collection: Classification Output of Kraken2
>    - *"Select a taxonomy database"*: `2022-03-08`
>    - *"Read name"*: `column:2`
>    - *"Taxonomy ID field"*: `column:3`
>
{: .hands_on}
 
 
## Visualize the taxonomical classification with Krona
 
__Krona__ allows hierarchical data to be explored with zooming, multi-layered pie charts. With this tool, we can easily visualize the composition of the bacterial communities and compare how the populations of microorganisms are modified according to the conditions of the environment.
 
> ### {% icon hands_on %} Hands-on: Visualize data with Krona
>
> 1. {% tool [Krona pie chart](xxx) %} with the following parameters:
>    - *"Type of input data"*: `taxonomy`
>    - *"Input file"*: Dataset collection: Output file of Convert Kraken
>    - *"Combine data from multiple data sets?"*: yes
>
{: .hands_on}
 
 
Let's take a look at the [result](https://usegalaxy.eu/datasets/4838ba20a6d86765e92bccb62d7f6daa/display/?preview=True&dataset=0&node=0&collapse=true&color=false&depth=8&font=11&key=true). Using the search bar we can check if certain taxa are present.
 
 
> ### {% icon question %} Questions
>
> 1. How many percent of the bacteria consists of the genus "paracoccus"?
> 2. Is there any *Escherichia coli* present? If yes, how many reads were found?
> 3. Where might the eukaryotic DNA come from?
>
>
> > ### {% icon solution %} Solution
> >
> > 1. 5 %.
> > 2. It is present and 162 reads were found.
> > 3. It's probably human contamination.
>
> {: .solution}
>
{: .question}

> <br/><br/>
> <center><iframe id="krona" src="krona.html" frameBorder="0" width="90%" height="600px"> ![Krona at bacteria level](../../images/metatranscriptomics/krona_bacteria.png) </iframe></center>
 
 
 
# Discussion: Choosing the right tool

When it comes to taxonomic assignment while analyzing metagenomic data, in this tutorial presented Kraken2 is not the only tool available. Several papers do benchmarking of different tools ({% cite Meyer.2022 %},{% cite Sczyrba.2017 %},{% cite Ye.2019 %}) and their results are presented in the following section, with focus on tools that are available in Galaxy.
When it comes to taxonomic profiling, thus investigating the abundance of specific taxa, the biggest problem is the abundance bias. It is introduced during isolation of DNA (which might work for some organisms better then for others) and by PCR duplicates during PCR amplification.
 
When benchmarking different classification tools, several metrics are used to compare their performance:
1. **Precision**: proportion of true positive species identified in the sample divided by number of total species identified by the method
2. **Recall**: proportion of true positive species divided by the number of distinct species actually in the sample
3. Precision-recall curve: each point represents the precision and recall scores at a specific abundance threshold → **area under the precision-recall curve (AUPR)**
 
4. **L2 distance**: representation of abundance profiles → how accurately the abundance of each species or genera in the resulting classification reflects the abundance of each species in the original biological sample (“ground truth”)
 

## Profiling tools
 
Profilers, which are tools that investigate relative abundances of taxa within a dataset, fall into three groups depending on their performance:
1. Profilers, that correctly predict relative abundances
2. Precise profilers (suitable, when many false positives would increase cost and effort in downstream analysis)
3. Profilers with high recall (suitable for pathogen detection, when the failure of detecting an organism can have severe negative consequences)
 
However, some characteristics are common to all profilers:
- Most profilers only perform well until the family level
- Drastic decrease in performance between family and genus level, while little change between order and family level
- Fidelity of abundance estimates decreases notably when viruses and plasmids were present
- Taxonomic profilers vs profiles from taxonomic binning:
Precision and recall of the taxonomic binners were comparable to that of the profilers;
abundance estimation at higher ranks was more problematic for the binners
 
MetaPhlAn 2.0 belongs to the group of precise profilers. On the basis of the average of precision and recall, over all samples and taxonomic ranks, MetaPhlAn 2.0 performed second best of all 10 profilers tested.


> | Tool                | Best method across metrics for* :   | additional features           | available in Galaxy   |
> | -----------------   |-------------------------------------|
> | mOTUs 1.1.1         | -                                   |most memory efficient          | no                    |
> | mOTUs 2.5.1         |marine; plant-associated             | -                             | no                    |
> | mOTUs v.cami1       |strain-madness                       | -                             | no                    |
> | MetaPhlAn 2.9.21    |plant-associated                     | -                             | yes
> | MetaPhlAn 2.9.22    |marine; strain-madness               | -                             | yes
> | DUDes v.cami1       |strain-madness                       | -                             | no
> | FOCUS 1.5           | -                                   |fastest; most memory efficient | no                    |
> | Bracken 2.2         | -                                   |fastest                        | yes (version 2.7)     |
> | Bracken 2.6         |plant-associated                     | -                             | yes (version 2.7)     |
{: .matrix}

*metagenome benchmark datasets created by {% cite Meyer.2022 %} representing a marine, a high strain diversity environment (‘strain-madness’) and a plant-associated environment including fungal genomes and host plant material

 
 
