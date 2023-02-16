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
time_estimation: 1H
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
 
The dataset we will use for this tutorial comes from an oasis in the mexican desert called Cuatro Ciénegas {% cite Okie.2020 %}. The researchers were interested in genomic traits that affect the rates and costs of biochemical information processing within cells. They performed a whole-ecosystem experiment, thus fertilizing the pond to achieve nutrient enriched conditions. The microbe samples collected from the Lagunita Fertilized Pond is calles JP4D in this tutorial. As control they used samples from a control mesocosm, called JC1A in this tutorial. You will realize that the datasets differ in size, but according to the authors this doesn't matter for their analysis of genomic traits. Also, they underline that differences between the two samples reflect trait-mediated ecological dynamics instead of microevolutionary changes as the duration of the experiment was only 32 days. This means that depending on available nutrients, specific lineages within the pond grow more successfully than others because of their genomic traits.  

The datafiles are named according to the first four characters of the filenames.
It is a collection of paired-end data with R1 being the forward reads and R2 being the reverse reads. Additionally, the reads have been trimmed using [__cutadapt__ ](https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/quality-control/tutorial.html#trim-and-filter---short-reads)

So let's get started with uploading the datasets!

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this exercise
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the FASTQ file pairs from [Zenodo]({{ page.zenodo_link }}) or a data library: xxx
>    - `JP4D_R1.fastq.gz`
>    - `JP4D_R2.fastq.gz`
    -`JC1A_R1.fastq.gz`
    - `JC1A_R2.fastq.gz`
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
         - *"Format report output like Kraken 1's kraken-mpa-report"*: `No`
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
5. indicates LCA mapping of each k-mer in the sequence |:| indicates end of first read, start of second read for paired reads
 
 
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

 
As both files contain a lot of information, we will use __Krona__ {% cite Ondov.2011 %} to visualize the data.
 
# Hands on: estimating species abundance

A "simple and worthwile addition to Kraken for better abundance estimates" {% cite Ye.2019 %} is called __Bracken__ (Bayesian Reestimation of Abundance after Classification with Kraken). Instead of only using proportions of classified reads, it takes a probabilistic approach to generate final abundance profiles. It works by re-distributing reads in the taxonomic tree: "Reads assigned to nodes above the species level are distributed down to the species nodes, while reads assigned at the strain level are re-distributed upward to their parent species." {% cite Lu.2017 %}

> ### {% icon hands_on %} Hands-on: estimate species abundance with Bracken (Name in Galaxy: Estimate Abundance at Taxonomic Level)
>
> 1. {% tool [Bracken](https://toolshed.g2.bx.psu.edu/view/iuc/bracken/79450f7fd718) %} with the following parameters:
      - *"Kraken report file"*: Report: Kraken2 on dataset collection  (! do not use the classification file)
      - *"Select a kmer distribution"*: `PlusPF (2021-05-17)` (! it is important to choose the same database that you also chose for Kraken2)
>     - *"Produce Kraken-Style Bracken report"*: `yes` 
>
{: .hands_on}

# Visualization of taxonomic assignment
 
Once we have assigned the corresponding taxa to each sequence, the next step is to properly visualize the data, for which we will use the __Krona pie chart__ tool ({% cite Ondov.2011 %}). But first, we need to convert the output generated by Kraken2 so it can be used as an input from the Krona tool.
 
## Convert output from Kraken2 or Bracken so it can be used for Krona
 
__Krakentools: Convert kraken report file__ tool {% cite Lu.2017 %} is designed to translate results of the Kraken metagenomic classifier (see citations below) to the full representation of NCBI taxonomy. It does so by using Taxonomic ID field provided by Kraken. The output of this tool can be directly visualized by the Krona tool.
 
> ### {% icon hands_on %} Hands-on: Krakentools: Convert kraken report file
>
> {% tool [Krakentools: Convert kraken report file](toolshed.g2.bx.psu.edu/repos/iuc/krakentools_kreport2krona/krakentools_kreport2krona/1.2+galaxy0) %} with the following parameters:
>    - *"Kraken report file"*:      Report: Kraken2 on dataset collection 
                                    (or: Estimate Abundance at Taxonomic Level on collection: Kraken style report for Bracken)

>
{: .hands_on}
 
 
## Visualize the taxonomical classification with Krona
 
__Krona__ allows hierarchical data to be explored with zooming, multi-layered pie charts. With this tool, we can easily visualize the composition of the bacterial communities and compare how the populations of microorganisms are modified according to the conditions of the environment.
 
> ### {% icon hands_on %} Hands-on: Visualize data with Krona
>
> 1. {% tool [Krona pie chart](toolshed.g2.bx.psu.edu/repos/crs4/taxonomy_krona_chart/taxonomy_krona_chart/2.7.1) %} with the following parameters:
>    - *"Type of input data"*: `tabular`
>    - *"Input file"*: Dataset collection: Output file of Krakentools: Convert kraken report file
>
{: .hands_on}
 
 
Let's take a look at the [result for Kraken2](https://usegalaxy.eu/datasets/4838ba20a6d86765e1d95e6747f5e42d/display/?preview=True&dataset=0&node=0&collapse=true&color=false&depth=8&font=11&key=true)/([result for Bracken](https://usegalaxy.eu/datasets/4838ba20a6d867652a9ece4a5ae1276d/display/?preview=True&dataset=0&node=0&collapse=true&color=false&depth=8&font=11&key=true).) Using the search bar we can check if certain taxa are present.
 
 
> ### {% icon question %} Questions
>
>Choose the sample JC1A: 
> 1. How many percent of the bacteria consists of the genus "paracoccus"?
> 2. Is there any *Escherichia coli* present? If yes, how many reads were found?
> 3. Where might the eukaryotic DNA come from?
>
>
> > ### {% icon solution %} Solution
> >
> > 1. 3 %.
> > 2. It is present and 73 reads were found.
> > 3. It's probably human contamination.
>
> {: .solution}
>
{: .question}

> <br/><br/>
> <center><iframe id="krona" src="krona.html" frameBorder="0" width="90%" height="600px"> ![Krona at bacteria level](../../images/metatranscriptomics/krona_bacteria.png) </iframe></center>
 
 
## Visualize the taxonomical classification with Phinch

__Phinch__ ({% cite Bik.2014 %}) offers the user to visualize large biological datasets like our taxonomic classification. Taxonomy Bar Charts, Bubble Charts, Sankey Diagrams, Donut Partitions and Attributes Column Chart can be generated using this tool. Additionally, several different samples can easily be compared.
As a first step, we need to convert the Kraken output file into a kraken-biom file to make it accessible for Phinch. Fot this, we need to add a metadata file, provided here. xxx
When generating a metadata file for your own data, you can take this as an example and find the general guidelines [here](http://qiime.org/documentation/file_formats.html#mapping-file-overview) 

> ### {% icon hands_on %} Hands-on: Phinch
>1. Use {% tool [Kraken-biom] (toolshed.g2.bx.psu.edu/repos/iuc/kraken_biom/kraken_biom/1.2.0+galaxy1) %} to convert Kraken2 report into the correct format for phinch with the following parameters.
>    - *"Input"*: `Kraken2 on dataset collection`
>    - *"Sample Metadata file"*: `metadata.txt`
>    - *"Output format"*: `JSON`
>
> 2. **View** the file in Phinch
>   -  Use {% tool [Phinch Visualisation] (interactive_tool_phinch) %} with the following paramters:
>    - *"Input"*: `Kraken-biom output file`

Important note: don't wait for Galaxy to finish the job! Your results are available directly, as Phinch is an interactive tool. To see them, click on User in the blue part right at the top of the Galaxy page, there you can find "Active interactive tools". When you go there, you will find Phinch running and can visit the website with your results by following the link "Phinch Visualisation of Kraken-biom output file" which is provided as the name of the job. 

{: .hands_on} 

When you follow the link to the Phinch webpage, you first see an overview of your samples. Here, you have the possibility to further filter your data, for example by date or location, depending on which information you provided in your metadata file. 
Next, you click on ‘proceed to gallery’ to see an overview of all visualization options. 

Let’s have a look at the **taxonomy bar chart**. Here, you can see the abundance of different taxa depicted in different colors in your samples. On top of the chart you can select which rank is supposed to be shown in the chart. You can also change the display options to for example switch between value und percentage. If you hover over a specific taxon, you get exact information of the taxon’s name and the taxonomy occurrence in the sample. Furthermore, there is a search bar that allows you to search for specific taxa (within the rank that you chose to depict) and gives you the option to hide specific taxons. 

Next, lets go back to the gallery and choose the **bubble chart**. Here, you can find the distribution of taxa across the whole dataset at the rank that you can choose above the chart. When hovering over the bubbles, you get additional information concerning the taxon. Clicking on one bubble gives you the direct comparison of the distribution of this taxon in the different samples.  To oder the bubbles according to their size you can choose the ‘list’ option shown right next to the taxonomy level.

Another displaying option is the **Sankey diagram**, that is depicting the abundance of taxonomies as a flow chart. Again, you can choose the taxonomy level that you want to show in your diagram. When clicking on one bar of the diagram, this part is enlarged for better view.

The **donut partition** summarizes the microbial community according to non-numerical attributes. In the drop-down menu at the top right corner, you can switch between the different attributes provided in the metadata file. In our case, you can for example choose the ‘environmental medium’ to see the difference between sediment and water (It doesn’t really make a lot of sense in our very simple dataset, as this will show the same result as sorting them by sample 0 and 1, but if attributes change across different samples this might be an interesting visualization option). When clicking on one part of the donut you will also find the distribution of the taxon across the samples. On the right hand side you can additionally choose if you’d like to have dynamic y axis or prefer standard y axis to compare different donuts with each other. 

The **attributes column chart** summarizes the microbial community according to numerical attributes. In the drop-down menu at the top right corner, you can switch between the different attributes provided in the metadata file. In our case, you can for example choose the ‘geographic location’ to (again, it doesn’t really make a lot of sense in our very simple dataset, as this will show the same result as sorting them by sample 0 and 1, but if attributes change across different samples this might be an interesting visualization option).

> ### {% icon question %} Questions
>Let's have a look at our data and try to answer some questions.
>1. How many percent of the sample reads are bacteria and how many are eukaryota? 
>2. How do you filter the data from eukaryotic reads?
> 
>
>
> > ### {% icon solution %} Solution
> >1. Sample 0: 75,65 % bacteria; 24,51 % eukaryota
      Sample 1: 92,70 % bacteria; 6,87 % eukaryota

      (go to taxonomy bar chart, choose kingdom and hover over the bars to find "taxonomy occurence in this sample")

> >2. xxx Search function of phinch shows confusing results when clicking on “hide”
> > xxx phnich states invalid date format but when I change the date format in the metadata file from 20120706 to 2012-07-06, kraken-biom doesnt work any more
> > 
>
> {: .solution}
>
{: .question}


## Visualize the taxonomical classification with Pavian

__Pavian__ (pathogen visualization and more) ({% cite Breitwieser.2020 %}) is an interactive visualization tool for metagenomic data. It was developed for the clinical metagenomic problem to find a disease-causing pathogen in a patient sample, but it is useful to analyze and visualize any kind of metagenomics data. 

> ### {% icon hands_on %} Hands-on: Pavian

>   -  Use {% tool [Pavian] (interactive_tool_pavian) %} with the following paramters:
>    - *"Kraken and MetaPhlAn-style reports"*: `Kraken report file`
                                                (or: Estimate Abundance at Taxonomic Level on collection: Kraken style report for Bracken)

Important note: don't wait for Galaxy to finish the job! Your results are available directly, as Pavian is an interactive tool. To see them, click on User in the blue part right at the top of the Galaxy page, there you can find "Active interactive tools". When you go there, you will find Pavian running and can visit the website with your results by following the link "Pavian Visualisation" which is provided as the name of the job. 

{: .hands_on} 


 Following the provided link, you will find the Pavian webpage. When choosing ‘Use data on server’ as data source, you should find your selected report files from galaxy. If you wish you can choose several sample files at once to be able to compare them in Pavian. After selecting ‘Read selected directories’, Pavian will show you a classification summary of all samples. 

Select ‘Comparison’ in the sidebar to find a table with taxa as rows and samples as columns. Above the table you can change several displaying options like which rank you would like to see or if you want the Z-scores to be displayed. Furthermore, the data can be filtered to exclude for example human reads as a host genome. 

Select ‘Sample’ in the sidebar to zoom into one sample with a Sankey diagram, which shows different taxons as lines with their width based on the quantity of the taxon. The diagram is shown for one sample but when hovering over nodes, the abundance of that specific taxon compared in all samples shown in a bar plot will appear on the right-hand side. 
One additional feature of Pavian is the ‘Alignment Viewer’ to investigate for example if reads of a taxon cover the genome. Genome coverage can provide an indication whether an assignment is artificial or not. 

> ### {% icon question %} Questions
>Let's have a look at our data and try to answer some questions.
>1. What is the taxonomy occurence in percent of the Family Rhodobacteraceae in the two samples?
>2. Regarding the question the authors had, when collecting the data: What could it biologically mean, if one bacteria family increases their occurence in fertilized pond (JP4D)?
> 
>
>
> > ### {% icon solution %} Solution
> >1. Sample JC1A: 12,82 %
      Sample JP4D: 36,06 %

> >2. It could mean that they have a survival advantage in the new environment. According to the authors this correlates with specific genomic traits that enable them to cope better with high nutrient availability. 
> > 
>
> {: .solution}
>
{: .question}

For further information on how to use Pavian please have a look into the [Pavian Walktrough](https://github.com/fbreitwieser/pavian/blob/master/vignettes/pavian-walkthrough.pdf) provided by F. Breitwieser.


# Discussion: Choosing the right tool

When it comes to taxonomic assignment while analyzing metagenomic data, in this tutorial presented Kraken2 is not the only tool available. Several papers do benchmarking of different tools ({% cite Meyer.2022 %},{% cite Sczyrba.2017 %},{% cite Ye.2019 %}) and their results are presented in the following section, with focus on tools that are available in Galaxy.
The benchmarking papers present different methods for comparing the available tools: The CAMI challenge is based on results of different labs that each used the CAMI dataset to perform their analysis on and send it back to the authors. In contrast, {% cite Ye.2019 %} performed all the analysis themselves. Additionally, the datasets used for both benchmarking approaches differ: in the CAMI dataset only ~30%-40% of reads are simulated from known taxa while the rest of the reads are from novel taxa, plasmids or simulated evolved strains. In contrast, {% cite Ye.2019 %} used International Metagenomics and Microbiome Standards Alliance (IMMSA) datasets, wherein the taxa are described better.
 
When benchmarking different classification tools, several metrics are used to compare their performance:
1. **Precision**: proportion of true positive species identified in the sample divided by number of total species identified by the method
2. **Recall**: proportion of true positive species divided by the number of distinct species actually in the sample
3. Precision-recall curve: each point represents the precision and recall scores at a specific abundance threshold → **area under the precision-recall curve (AUPR)**
 
4. **L2 distance**: representation of abundance profiles → how accurately the abundance of each species or genera in the resulting classification reflects the abundance of each species in the original biological sample (“ground truth”)

When it comes to taxonomic profiling, thus investigating the abundance of specific taxa, the biggest problem is the abundance bias. It is introduced during isolation of DNA (which might work for some organisms better then for others) and by PCR duplicates during PCR amplification.

## Profiling tools
 
Profilers, which are tools that investigate relative abundances of taxa within a dataset, fall into three groups depending on their performance:
1. Profilers, that correctly predict relative abundances
2. Precise profilers (suitable, when many false positives would increase cost and effort in downstream analysis)
3. Profilers with high recall (suitable for pathogen detection, when the failure of detecting an organism can have severe negative consequences)
 
However, some characteristics are common to all profilers:
- Most profilers only perform well until the family level
- Drastic decrease in performance between family and genus level, while little change between order and family level
- poorer performance of all profilers on CAMI datasets compared to International Metagenomics and Microbiome Standards Alliance (IMMSA) 
- Fidelity of abundance estimates decreases notably when viruses and plasmids were present
- high numbers of false positive calls at low abundance 
- Taxonomic profilers vs profiles from taxonomic binning:
Precision and recall of the taxonomic binners were comparable to that of the profilers;
abundance estimation at higher ranks was more problematic for the binners
 


| Tool      | version              | available in Galaxy | In CAMI challenge, best method across metrics for\* : | additional features            | {% cite Ye.2019 %} benchmarking                                                                                                                                     |
| --------- | -------------------- | ------------------- | ----------------------------------------------------- | ------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| mOTUs     |  1.1.1               | no                  | \-                                                    | most memory efficient          |                                                                                                                                                            |
| mOTUs     |  2.5.1               | no                  | marine; plant-associated                              | \-                             |                                                                                                                                                            |
| mOTUs     | v.cami1              | no                  | strain-madness                                        | \-                             |                                                                                                                                                            |
| MetaPhlAn |  2.9.21              | yes                 | plant-associated                                      | \-                             | recommended for low computational requirements (< 2 Gb of memory)                                                                                          |
| MetaPhlAn |  2.9.22              | yes                 | marine; strain-madness                                | \-                             |
| DUDes     | v.cami1              | no                  | strain-madness                                        | \-                             |                                                                                                                                                            |
| FOCUS 1.5 |  1.5                 | no                  | \-                                                    | fastest; most memory efficient |                                                                                                                                                            |
| Bracken   |  2.2                 | yes (version 2.7)   | \-                                                    | fastest                        |  
| Bracken   |  2.6                 | yes (version 2.7)   | plant-associated                                                               
| Kraken    |  2.0.8 beta (GSA,Sr) | yes                 | marine                                                | fastest; most memory efficient | provide good performance metrics<br>very fast on large numbers of samples<br>allow custom databases<br>when high amounts of memory (>100 Gb) are available |

*metagenome benchmark datasets created by {% cite Meyer.2022 %} representing a marine, a high strain diversity environment (‘strain-madness’) and a plant-associated environment including fungal genomes and host plant material

 
 
