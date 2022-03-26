---
layout: tutorial_hands_on

title: "Reproducing Critical Assessment of Metagenome Interpretation assembly challenge on marine dataset with Galaxy"
zenodo_link: ""
questions:
  - "How to reproduce CAMI challenge?"
  - "How to add/update tool in Galaxy?"
  - "How to benchmark results of CAMI challenge?"
objectives:
  - "Be familiar with CAMI challenge"
  - "Be able to select one challenge and benchmarking datasets"
  - "Be able to upload benchmarking datasets into Galaxy"
  - "Be familiar with Galaxy computational resources"
  - "Be familiar with existing tools in Galaxy"
  - "Be familiar with databases / reference genomes that are available in Galaxy"
  - "Be able to add/update the tool in Galaxy"
  - "Be able to choose assembling tool based on dataset metadata"
  - "Perform assembly with Flye, Megahit, Abyss, MetaSPAdes tools"
  - "Produce a benchmarking analysis of these assemblies with Quast, Bowtie2, Map with Minimap2, Samtools, and MultiQC tools"
  - "Create plots with Python to compare results"
time_estimation: "40H00M"
key_points:
  - "Galaxy could be a good solution to reproduce CAMI"
  - "CAMI2 consists of 4 challenges: assembly, genome binning, taxonomix binning, taxonomic profiling"
  - "CAMI2 provides 5 datasets: pathogen detection, marine, strain madness, plant-associated"
  - "Short reads and long reads can be assembled with different tools"
  - "Co-assembly and individual assembly should be performed for dataset"
contributors:
  - PlushZ
---


## Introduction


### Global context 

Metagenomics analysis involves complex computational methods like assembly, genome binning, taxonomic binning, taxonomic profiling and other. It is only when these initial data processing methods make sense that downstream analyses and information extraction are meaningful. Even though vast progress has been made over the last few years, none of these approaches can fully recover the complex information encoded in metagenomes. These approaches all rely on simplifying assumptions that can have serious limitations.

It is usually difficult for microbiome data analysts to know which tools to use for each step. They must rely on evaluations made by the developers of new or improved methods. However, these evaluations are rarely comparable: there is no general standard for the evaluation of computational methods in metagenomics. As a result, users may not be properly informed and compute predictions may be misinterpreted.


### The CAMI Challenge

Essentially, the critical assessment concept involves evaluating a theory, situation, statement, or something else with the goal of supporting its dominant paradigms or disproving them as well as suggesting a better alternative plan. Any view or conclusion needs to be backed up by credible evidence in order to be considered critical. In metagenomic research, critical assessment of the data interpretation is especially important, because merely accepting the data as truth would not suffice.

The “Critical Assessment of Metagenome Interpretation” (CAMI) was established in 2012 [[1]](https://www.nature.com/articles/nmeth.4458). The objective of CAMI is to evaluate metagenomics methods independently, comprehensively, and without bias. The initiative provides users with comprehensive information about the performance of methods in all relevant scenarios. Therefore, it assists users in the selection and application of methods as well as their proper interpretation.

During the 1st CAMI challenge, extensive metagenome benchmark data sets were generated from newly sequenced genomes of around 700 microbial isolates and approximately 600 circular elements that were distinct from strains, species, genera or orders represented by public genomes. Four challenges were suggested (assembly, binning, taxonomic profiling, taxonomic binning). Overall, 16 teams worldwide submitted 215 submissions to the CAMI1 challenge, consisting of 25 programs and 36 biobox implementations, with permission to publish.

The 2nd CAMI started in 2019 {%cite noauthor_critical_nodate-1}. It includes 4 different challenges: 



* An assembly challenge
* A genome binning challenge
* A taxonomic binning challenge
* A taxonomic profiling challenge
* Clinical pathogen prediction challenge

Each challenge uses a similar set of benchmark datasets reproducing different environments and different sequencing technologies  long and short reads). The six benchmark datasets—reflecting a range of complexities—were created from 1,680 microbial genomes and 599 circular elements of viruses and plasmid:



* 2 ‘toy’ datasets created from public data and provided before a challenge, 
* 3 ‘challenge’ datasets derived exclusively from genomic data that were not publicly available at the time
    * Marine
        * Simulated short read and long read shotgun metagenome data from samples at different seafloor locations of a marine environment
    * Plant-associated
        * Simulated short read and long read shotgun metagenome data from samples taken from a plant rhizosphere environment
    * strain madness
        * Simulated long read and short read shotgun metagenome data, including a large amount of strain-level variation
* a pathogen detection challenge dataset was offered, based on a clinical metagenome sample of blood from a critically ill patient with hemorrhagic fever of unknown cause.

During the challenges, Datasets for the challenge were available for download only to participants via the [CAMI portal](https://data.cami-challenge.org/). After the challenges, all CAMI benchmark datasets were made available with digital object identifiers (DOIs) ([Table 2](https://www.nature.com/articles/s41596-020-00480-3/tables/2)) and the genomic data are now in public sequence repositories such as the National Center for Biotechnology Information (NCBI) to be used for further benchmarking in the field.

Different tools were run with different parameters on the various datasets. Participants ran their preferred tools with their preferred set of parameters and  submitted their results along with either a Docker container containing the complete workflow, a bioconda script or a software repository with detailed installation instructions, specifying all parameter settings and reference databases used. A total of 5,002 submissions were received from 30 external teams and CAMI developers for the four challenge datasets. Following that, the CAMI developers evaluated the results using standardised metrics and then made sense of the different results described in [Meyer et al, 2020](https://www.biorxiv.org/content/10.1101/2021.07.12.451567v1)

 Galaxy could be considered as a good platform for such a challenge as CAMI.

In addition to providing the computational resources, Galaxy provides full reproducibility, which makes it a great solution for CAMI and similar types of challenges. Everything required to complete such challenges is collected in one location:



* Datasets which can be uploaded to data libraries and shared with public
* Databases and tools of different versions can be available to use
* History of tool runs with documented parameters of tool runs, inputs, outputs, command line, job parameters (CPU, memory usage, runtime, etc.)
* Galaxy workflows can be extracted and shared for the reproducibility
* Etc.

In Galaxy, users can choose to share their history with certain users or to publish it so that everyone can access the history.


### Motivation

In this tutorial we show how Galaxy could be used as a platform to support the next CAMI challenges or for other similar critical assessment challenges, but also what were the issues we got and possible solutions (on Galaxy side but for organisers of such challenges).

## Select the challenge to reproduce

In this section we  choose the challenge to reproduce. As the entire challenge is not a one-day task and takes months to complete, we should focus on one challenge out of four - assembly, genome binning, taxonomic binning, or taxonomic profiling. The choice we make depends on our interests and objectives.


### Description of CAMI challenges


#### Assembly challenge

Assembly is one of the important components of metagenomic analysis. The challenge is called to assess the quality and performance of tools for assembling strain-resolved genomes. 

Three “challenge” datasets were offered:: 



* marine, 
* strain madness, and 
* plant associated data. 

Two types of assemblies were assessed: single-sample assembly, i.e. assembly of reads from only one sample at a time, and co-assembly, i.e. assembly of all reads from several samples together. For the sake of comparison, gold standard assemblies, i.e. the correct result which can be used to compare with and benchmark other methods’ results, were created with all regions covered by at least one read and were provided only after the challenges.

Participants provided 155 submissions and 20 assembler versions of 10 assembling tools to be evaluated during this challenge. The different tools performed differently while measuring various metrics. By CAMI2 developers there are provided plots ([Fig 1](https://www.biorxiv.org/content/biorxiv/early/2021/07/12/2021.07.12.451567/F1.large.jpg?width=800&height=600&carousel=1)) of benchmarking analyses of different tools depending on different statistics.



<p id="gdcalert1" ><span style="color: red; font-weight: bold">>>>>>  gd2md-html alert: inline image link here (to images/image1.png). Store image on your image server and adjust path/filename/extension if necessary. </span><br>(<a href="#">Back to top</a>)(<a href="#gdcalert2">Next alert</a>)<br><span style="color: red; font-weight: bold">>>>>> </span></p>


![alt_text](images/image1.png "image_tooltip")


Fig 1

Overall, GATB[[3]](https://pubmed.ncbi.nlm.nih.gov/24990603/) ranked best on the strain madness data across most metrics (mismatches, misassemblies, strain recall, and strain precision), while HipMer[[4]](https://ieeexplore.ieee.org/abstract/document/8665813) [[5]](https://ieeexplore.ieee.org/abstract/document/7832788) [[6]](https://www.nature.com/articles/s41598-020-67416-5) on the plant-associated data. Compared to the first CAMI challenge A-STAR[ [7]](https://www.a-star.edu.sg/) considerably improved in genome fraction for the strain madness data. The other leader SPAdes [[8]](https://pubmed.ncbi.nlm.nih.gov/22506599/) [[9]](https://genome.cshlp.org/content/27/5/824?ijkey=fef5a6b24fb0024b39509480f85afaa81724d3e8&keytype2=tf_ipsecsha) was not introduced in the CAMI1 challenge [[1]](https://www.nature.com/articles/nmeth.4458) and performed well in CAMI2. For the type of assembly, the different tools performed differently. Single-sample assembly was done well by all assemblers. HipMer [[4]](https://ieeexplore.ieee.org/abstract/document/8665813) [[5]](https://ieeexplore.ieee.org/abstract/document/7832788) [[6]](https://www.nature.com/articles/s41598-020-67416-5) performed better on single samples as well as for pooled samples.

#### Genome binning challenge

Genome binning  clusters or classifies sequences or contigs into bins representing genomes. 

It takes as input the output of assembly. Two types of short read assemblies were given as inputs for CAMI2 [[2]](https://www.biorxiv.org/content/10.1101/2021.07.12.451567v1): MEGAHIT [[10]](https://www.biorxiv.org/lookup/external-ref?access_num=25609793&link_type=MED&atom=%2Fbiorxiv%2Fearly%2F2021%2F07%2F12%2F2021.07.12.451567.atom) assembly and Gold Standard Assembly. Both types for strain madness, marine, and plant-associated datasets, so 6 possible datasets in total.

Participants submitted 95 results and 10 binning tools (overall, 18 binner versions). The results were evaluated based on .

 



* the average purity of bins, 
    * The average purity is the fraction of correctly assigned base pairs for all assignments to a given bin averaged over all predicted genome bins, where unmapped genomes are not considered. 
    * This value reflects how trustworthy the bin assignments are on average. 
* completeness of genomes recovered,
* Adjusted Rand Index which quantifies binning performance for the overall dataset.
    * According to the [paper](https://www.biorxiv.org/lookup/google-scholar?link_type=googlescholar&gs_type=article&author[0]=F.+Meyer&title=AMBER:+Assessment+of+Metagenome+BinnERs&publication_year=2018&journal=Gigascience&volume=7) [[11]](https://academic.oup.com/gigascience/article/7/6/giy069/5034950?login=true), genome binners generate groups or clusters of reads and contigs for a given dataset. Instead of calculating performance metrics established with a bin-to-genome mapping, the quality of a clustering can be evaluated by measuring the similarity between the obtained and correct cluster partitions of the dataset, corresponding here to the predicted genome bins and the gold standard contig or read genome assignments, respectively. Two contigs or reads of the same genome that are placed in the same predicted genome bin are considered true positives TP. Two contigs or reads of different genomes that are placed in different bins are considered true negatives TN. The Rand index ranges from 0 to 1 and is the number of true pairs, TP + TN, divided by the total number of pairs. However, for a random clustering of the dataset, the Rand index would be larger than 0. The adjusted Rand index (ARI) corrects for this by subtracting the expected value for the Rand index and normalising the resulting value, such that the values still range from 0 to 1.

The performance of different genome binning tools differed depending on metrics, datasets, and assembly type. Overall, the best trade-off performances were given by MetaBinner[[12]](https://www.biorxiv.org/content/10.1101/2021.07.25.453671v1) on marine and strain madness assemblies, and CONCOCT [[13]](https://www.biorxiv.org/lookup/external-ref?access_num=10.1038/nmeth.3103&link_type=DOI) on plant-associated assemblies.


#### Taxonomic binning challenge

Taxonomic binners group sequences into bins labelled with a taxonomic identifier.. For the challenge,  input datasets were marine, strain madness, and plant-associated of 2 types: reads or gold standard assemblies.

547 results for 7 tools were submitted and evaluated using   3 metrics : 



* average purity, 
* completeness of bins, and 
* the accuracy per sample (the fraction of contigs, or base pairs, that have been assigned by a method to the correct taxa for a taxonomic rank).

Performances for all datasets decreased from genus to species for all metrics, most notably for completeness by 22.2%.

Across all datasets and all metrics the best performance was shown by MEGAN [[14]](https://www.biorxiv.org/lookup/external-ref?access_num=10.1371/journal.pcbi.1004957&link_type=DOI), closely followed by Kraken [[15]](https://www.biorxiv.org/lookup/google-scholar?link_type=googlescholar&gs_type=article&author[0]=D.%20E.+Wood&author[1]=J.+Lu&author[2]=B+Langmead&title=Improved+metagenomic+analysis+with+Kraken+2&publication_year=2019&journal=Genome+Biol&volume=20) v.2.0.8-beta. However, in terms of certain metrics different binners performed better.

There was a post-processing data approach used to improve the results. This approach means filtering of the 1% smallest predicted bins per taxonomic level. Filtering increased purity whereas reduced completeness. Accuracy was not affected. The winner MEGAN was not changed though.


#### Taxonomic profiling challenge

Using metagenome samples, taxonomic profilers identify and quantify microbial community taxa. There were predicted taxonomic identities and relative abundances of microbial community members for the 64 samples of the mouse gut dataset. There were 4195 results for 13 tools evaluated. The majority of results were from short-read samples and a few from long-read.

The performance of the tools were evaluated using: 



* purity of identified taxa, 
* taxon abundance, 
    * For quantifying relative abundance estimates, the L1 norm and weighted UniFrac error are determined. The L1 norm assesses relative abundance estimates of taxa at a taxonomic rank, on the basis of the sum of the absolute differences between the true and predicted abundances across all taxa. The weighted UniFrac error uses a taxonomic tree storing the predicted abundances at the appropriate nodes for eight major taxonomic ranks. The UniFrac error is the total amount of predicted abundances that must be moved along the edges of the tree to cause them to overlap with the true relative abundances. 
* alpha diversity of profiling results.

There was a notable performance drop observed from genus to species level rank for all datasets: marina, strain madness, and plant-associated.

mOTUs [[16] ](https://www.biorxiv.org/lookup/external-ref?access_num=10.1038/s41467-019-08844-4&link_type=DOI)v.2.5.1 performed the best on the marine data at both genus and species. For strain madness data, leaders for genus rank were MetaPhlAn [[17]](https://www.nature.com/articles/nmeth.2066) [[18]](https://elifesciences.org/articles/65088) v.2.9.2, MetaPhyler [[19]](https://www.biorxiv.org/lookup/google-scholar?link_type=googlescholar&gs_type=article&author[0]=B.+Liu&author[1]=T.+Gibbons&author[2]=M.+Ghodsi&author[3]=T.+Treangen&author[4]=M+Pop&title=Accurate+and+fast+estimation+of+taxonomic+profiles+from+metagenomic+shotgun+sequences&publication_year=2011&journal=BMC+Genomics&volume=12) v.1.2571, and mOTUs of the CAMI1. For plant-associated data Bracken v2.6 performed best for completeness across ranks.

Abundances across ranks and submissions were on average predicted better for strain madness than marine data. On the marine data, mOTUs v.2.5.1 had the lowest L1 error while mOTUs v.2.5.1 and MetaPhlAn v.2.9.22 both had the lowest UniFrac error.

On the strain madness data, mOTUs cami1 performed best in L1 norm error whilst MetaPhlAn v.2.9.22 had the lowest UniFrac error.

On the plant-associated data, Bracken v.2.6 had the lowest L1 norm error across ranks and had the lowest UniFrac as “sourmash gather 3.3.2 k31” did.


### Compare tools used in CAMI2 with what is already exist in Galaxy

We compared tools submitted for every challenge to understand the start point of CAMI2 reproduction in the Galaxy. We collected metadata for each tool:



* source code, 
* Publication that describes the tool
* Tool version
* Availability in Galaxy
* owner/ developer in Galaxy	
* Galaxy wrapper directory	
* Tool version in Galaxy
* Relevance in Galaxy	
* Availability in conda-forge / bioconda	
* Notes on their performance	
* Ranking according to CAMI2 results	

DETAIL BOX (The result table - [gsheet](https://docs.google.com/spreadsheets/d/1bKDz8uTu1nsuqmwLNLDM_vQ05lAnHedxbj1YVKlkQBc/edit#gid=0))

After this comparison it turned out that the Assembly challenge is more interesting for us to reproduce on Galaxy platform. One reason for this choice is that highly ranked tools in CAMI2 (like Megahit, Flye, MetaSPAdes, etc.) are already available in Galaxy, so that we can focus more on benchmarking analysis of assembly results than on adding missing tools into Galaxy. 


## Select and get the data to benchmark

Once we have decided to reproduce the assembly challenge, we can move forward to determine which dataset is more representative for our purposes. There is no need to use all 3 benchmark datasets (marine, strain madness, plant-associated) or 2 toy datasets in order to demonstrate the utility of Galaxy for future CAMI and other similar tasks.


### Datasets description


##### Toy datasets

As in the first CAMI challenge metagenome “toy” (or “practice”) benchmark datasets were created from public genomes and provided together with the standard of truth before the challenges, to enable contest participants to familiarise themselves with data types and formats. These included



*  a 49-sample simulated metagenome data from Human Microbiome data to represent five different body sites of the human host, namely gastrointestinal tract, oral cavity, airways, skin and urogenital tract.

    DETAIL BOX - 

*  a 64-sample simulated metagenome data from the guts of different mice, vendors and positions in the gut. (**[https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMISIM_MOUSEGUT](https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMISIM_MOUSEGUT)**)
    * According to the tutorial [[20]](https://www.nature.com/articles/s41596-020-00480-3) the mouse gut metagenome toy dataset was generated with CAMISIM v.0.2 [[21]](https://pubmed.ncbi.nlm.nih.gov/30736849/) using a microbial community genome abundance distribution modelled from 791 public prokaryotic genomes marked as at least ‘scaffolds’ in the NCBI RefSeq [[22]](https://pubmed.ncbi.nlm.nih.gov/17130148/). They comprise 8 phyla, 18 classes, 26 orders, 50 families, 157 genera, and 549 species. The community genome abundance distribution matches as closely as possible the 16S taxonomic profiles for 64 mouse gut samples. As such, this dataset enables us to assess how well sequenced community members can be characterised with different techniques from the metagenomes of similar communities. On average, within each of the 64 samples, 91.8 genomes are represented. Both long- and short-read metagenomic sequencing data are available. The runtime to generate this data was approximately 3 weeks using eight CPU cores of a computer with an AMD Opteron 6378 CPU and 968 GB of main memory.

Both include  5 Gb long paired-end reads (Pacific Biosciences, variable length with a mean of 3000 bp) and 5 Gb short paired-end reads (Illumina HiSeq2000, 150 bp).


##### Challenge datasets

For the assembly challenge itself, 4 “challenge” datasets were offered: 



* clinical pathogen detection,
* marine, 
* strain madness, 
* plant associated data. 


###### Clinical Pathogen detection

The clinical pathogen detection dataset, a short-read metagenomic sequencing dataset of a blood sample from a patient with hemorrhagic fever of unknown cause, was excluded because it is provided only for the last challenge: identify a causal pathogen together with further pathogens.

DETAIL Box - [clinical pathogen detection](https://docs.google.com/spreadsheets/d/1bKDz8uTu1nsuqmwLNLDM_vQ05lAnHedxbj1YVKlkQBc/edit#gid=1317732159) 


###### Marine

Detail box - [The dataset represents a marine environment](https://docs.google.com/spreadsheets/d/1bKDz8uTu1nsuqmwLNLDM_vQ05lAnHedxbj1YVKlkQBc/edit#gid=1420005524)

The dataset represents a marine environment. It was created with CAMISIM from BIOM profiles of a deep-sea environment, using 155 newly sequenced marine isolate genomes from this environment and 622 public genomes with matching taxonomic provenance from MarRef, a manually curated database with completely sequenced marine genomes. Additionally, 200 newly sequenced circular elements including plasmids and viruses were added.

The dataset (100Gb) consists of 2x10 samples (5 Gb each) paired-end short (Illumina, with the length of 150bp) and long-read (Nanopore, with the length of 7408 bp) sequences. 


###### Strain madness

The [dataset](https://docs.google.com/spreadsheets/d/1bKDz8uTu1nsuqmwLNLDM_vQ05lAnHedxbj1YVKlkQBc/edit#gid=548769407) represents a very high strain diversity environment (“strain madness”) generated with CAMISIM, using 408 newly sequenced genomes. Dataset (400 Gb) includes 2x100 samples (2 Gb each) paired-rend short (Illumina, the length of 150bp) and long(Nanopore with the length of 7408 bp) read sequences. 


###### Plant-associated

**_Links to raw data: _**

The [dataset](https://docs.google.com/spreadsheets/d/1bKDz8uTu1nsuqmwLNLDM_vQ05lAnHedxbj1YVKlkQBc/edit#gid=934346728) represents a plant-associated environment (that included fungal genomes and host plant material). It includes 894 genomes:



* 224 are from the proGenomes terrestrial representative genomes, 
* 216 are newly sequenced genomes from an _Arabidopsis thaliana_ root rhizosphere,
* 55 are fungal genomes associated with the rhizosphere, 
* 398 are plasmids or circular elements and one _Arabidopsis thaliana_ genome. 

15.3% (137) of these genomes have at least one closely related genome present. 90% of metagenome sequence data originate from bacterial genomes, 9% are fungal genome sequences, and 1% is from _A. thaliana_. 399 circular elements including plasmids and viruses were added.

To evaluate the assembly quality of single-sample versus cross-assembly strategies, 23 new genomes from eight clusters of closely related genomes were selected and added to the dataset in certain samples with predetermined abundances.

Dataset (315Gb) includes 3x21 (5 Gb each paired-end read) paired-end short (with the length of 150bp) and long (with the length of 7408 bp or 3000 bp (sd 1000 bp)) read sequences


#### Dataset choice

We analysed all datasets and compared them with each other in terms of number of samples, dataset size, average read length.

Detail box -  The whole table is [here](https://docs.google.com/spreadsheets/d/1bKDz8uTu1nsuqmwLNLDM_vQ05lAnHedxbj1YVKlkQBc/edit#gid=650574223). 



1. After close consideration of all datasets we decide to choose a marine dataset:It contains paired-end reads from 10 samples, both short (sequenced with Illumina method) and long (sequenced with nanopore method). ITmeans that we can try different tools on either one sample or pooled sample collection and also rub tools working with long reads and short reads, and tools working with both short and long at the same time. 
2. Overall size of all marine data is 180 Gb. This size is big enough and allows us to test tools and measure such metrics as runtime and memory usage. Other datasets (strain madness and plant-associated) are too large, and since we have a limit in time of this tutorial we decided to give preference to smaller but highly representable marine dataset.


### Get datasets

So the assembly challenge is our choice of challenge that we are interested in. And marine dataset is our choice of dataset that we want to launch an assembly challenge on. Let’s start with the first step - uploading data to Galaxy.


#### Download data from CAMI2 resources.

Since data is of large size we cannot just download them to our computer by clicking the button download. Additionally, data cannot be downloaded at once, datasets are divided by samples - one file per one sample. There are also different files for long and short reads of the same sample.

It is possible to download datasets using DOI links from [Table 2](https://www.nature.com/articles/s41596-020-00480-3/tables/2) CAMI Tutorial [[20]](https://www.nature.com/articles/s41596-020-00480-3)

Hands-on



1. Use the table below to find unique identifier (DOI) of dataset you need

Detail box inside hands-on:


<table>
  <tr>
   <td><strong>CAMI benchmark dataset</strong>
   </td>
   <td><strong>DOI</strong>
   </td>
  </tr>
  <tr>
   <td>CAMI I: low, medium, high complexity, and ‘toy’ datasets
   </td>
   <td>10.5524/100344
   </td>
  </tr>
  <tr>
   <td>CAMI II: human microbiome project and mouse gut toy datasets
   </td>
   <td>10.4126/FRL01-006425518 and 10.4126/FRL01-006421672
   </td>
  </tr>
  <tr>
   <td>CAMI II: marine, strain madness, rhizosphere, and pathogen detection challenge datasets
   </td>
   <td>10.4126/FRL01-006425521
   </td>
  </tr>
</table>




2. Use [https://search.datacite.org/](https://search.datacite.org/) to search for data via DOI, and click the link that was found for the DOI. 
3. On the page of dataset there is the link provided to download dataset -  for marine dataset all samples available via [https://frl.publisso.de/data/frl:6425521/marine/](https://frl.publisso.de/data/frl:6425521/marine/)
4. Download datasets on your computer

Detail box:

Another way to download datasets is to use CAMI Client. [Here ](https://www.microbiome-cosi.org/cami/resources/cami-client) you can find instructions on how to use it.

Hands-on inside Detail box



1. Install Java 8 on your computer.
2. Create an account on[ CAMI portal](https://data.cami-challenge.org/)
3. Download the [camiClient.jar](https://data.cami-challenge.org/camiClient.jar)
4. Get text-file with links for dataset from [CAMI portal](https://data.cami-challenge.org/). You can download this file in the section “2nd CAMI Challenge Marine Dataset”
5. You can download data with CAMI client using the following command line in the terminal:

    ```
java -jar camiClient.jar -d <linkfileLocation> <targetDirectory> -p PATTERN
```




#### Upload data into Galaxy

There are different ways to upload data into the Galaxy. In case of a not very large dataset it’s more convenient to upload data directly from your computer to Galaxy. In this tutorial we use a large dataset, thus, we use FTP server for this purpose.

Hands-on



1. Create new history
2. Follow instructions [here](https://galaxyproject.eu/ftp/) to upload datasets into FTP Galaxy server
3. Press “upload data” -> “choose remote files” -> “FTP directory” and upload your data to your history from FTP server


#### Preprocessing data

We got the marine data in our history, with one file per sample (10 samples for long reads and 10 samples for short reads). We need to prepare our data.



1. We should organise our data into collections. For deeper learning about dataset collections follow this tutorial[ Using dataset collections tutoria](https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/collections/tutorial.html)l
2. These files consist of interleaved paired-end reads, i.e. forward and reverse reads are all together in one single file and not 2 as often. This can be a problem for some assemblers as not every assembling tool can take interleaved paired-end reads as input.For co-assembly of all samples, some assemblers expect one file with all samples  combined together.

Hands-on



1. Take files for short reads for all samples
2. Create a collection with all samples- We now have a collection with all files. Let's now combine all sample files into one file
3. Combine all samples from collection into one file -  follow the section “[Collapse data into a single dataset](https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/collections/tutorial.html#collapse-data-into-a-single-dataset)” in [Using dataset collections tutorial](https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/collections/tutorial.html)

    We can now split forward and reverse reads into 2 files

4.  FASTQ splitter with following parameters:
    1. **FASTQ reads:** output of Collapse collection tool
5. Check if the names of the output generated by Fastq splitter include the information about forward and reverse, rename them otherwise
6. Repeat 1-5 steps for long reads. Eventually we have 2 collections (one for short reads, one for long reads)

    Comment box:


    While working on this project we realised that both outputs (forward and reverse) had the same name - “FASTQ splitter on data X”, where X - is the internal number of data in this Galaxy history. Since we have run FASTQ splitter on one dataset and got 2 outputs, as a result filenames of both outputs were the same. Why is it not good? It can be mistaken by assembler (or any other tool where we want to use FASTQ splitter output files) because it will not know which file is reverse and which is forward.We updated the Galaxy tool to fix this issue



## Select and prepare the tools to run

Based on [CAMI2 paper](https://www.biorxiv.org/content/10.1101/2021.07.12.451567v1.full) [2] and [CAMI2 tutorial paper](https://www.nature.com/articles/s41596-020-00480-3#Bib1) [20], we can compare tools on a set of metrics to select the one to use for an analysis but also here to run the challenge. Tools performance ranking.

**_Performance of assembly tools on the marine dataset_**:


<table>
  <tr>
   <td>
   </td>
   <td>Genome fraction (%)
   </td>
   <td>Mismatches per 100 kbp
   </td>
   <td>Misassemblies
   </td>
   <td>NGA50
   </td>
   <td>Strain recall
   </td>
   <td>Strain precision
   </td>
  </tr>
  <tr>
   <td>GSA (gold standard assembly)
   </td>
   <td>76.9
   </td>
   <td>0
   </td>
   <td>0
   </td>
   <td>682,777
   </td>
   <td>The upper bound for strain recall was 54.9%, as provided by the marine gold standard assembly
   </td>
   <td>100
   </td>
  </tr>
  <tr>
   <td>ABySS [<a href="https://genome.cshlp.org/content/19/6/1117">23</a>]
   </td>
   <td>There were selected 50 unique, public genomes present as a single contig in the gold standard and with annotated 16S sequences. The mean fraction were very accurate for ABySS and HipMer (&lt;1% divergence). 
   </td>
   <td>
   </td>
   <td>ABySS created the fewest misassemblies for the marine
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>achieved 100% precision.
<p>
For the unique genome, the highest strain precision (100%)
   </td>
  </tr>
  <tr>
   <td>Ray <a href="https://www.biorxiv.org/lookup/google-scholar?link_type=googlescholar&gs_type=article&author[0]=S.+Boisvert&author[1]=F.+Raymond&author[2]=E.+Godzaridis&author[3]=F.+Laviolette&author[4]=J+Corbeil&title=Ray+Meta:+scalable+de+novo+metagenome+assembly+and+profiling&publication_year=2012&journal=Genome+Biol&volume=13">[24]</a>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>achieved 100% precision.
<p>
For the unique genome, the highest strain precision (100%)
   </td>
  </tr>
  <tr>
   <td>A-STAR
   </td>
   <td>A-STAR excelled in terms of genome fraction on marine and strain madness data sets. A-STAR improved the genome fraction to 44.1% on the marine dataset.
<p>
On marine common genomes, A-STAR (26.7%) achieved the highest genome fractions.
<p>
For unique genome, A-STAR provided the most complete assemblies (55.3% genome fraction).
<p>
There were selected 50 unique, public genomes present as a single contig in the gold standard and with annotated 16S sequences. A-STAR partially recovered 102 (78%) of 131 16S gold standard sequences.
   </td>
   <td>On marine and strain-madness datasets A-STAR created more mismatches than others. For marine dataset: 773 mm/100 kb.
   </td>
   <td>On marine and strain-madness datasets A-STAR created more misassemblies than others.
   </td>
   <td>
   </td>
   <td>On marine common genomes, A-STAR had 7.5% recall - the second highest.
   </td>
   <td>On marine common genomes, A-STAR had  69.4% precision  - the second highest
   </td>
  </tr>
  <tr>
   <td>OPERA-MS <a href="https://www.biorxiv.org/lookup/external-ref?access_num=10.1038/s41587-019-0191-2&link_type=DOI">[25]</a>
   </td>
   <td>There were selected 50 unique, public genomes present as a single contig in the gold standard and with annotated 16S sequences. The hybrid assembler OPERA-MS recovered one of the most complete 16S sequences (mean recovered gene fraction 47.1%)
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>The most contiguous assemblies were provided by the hybrid assembler OPERA-MS for the marine data, with an average NGA50 of 28,244 across genomes.
<p>
For the unique genome, OPERA-MS has an exceptional average NGA50 (187,083, 75% of the gold standard NGA50).
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>MEGAHIT
   </td>
   <td>41.1
<p>
On marine common genomes, MEGAHIT (24.6%) achieved the second highest (after A-STAR) genome fractions.
<p>
There were selected 50 unique, public genomes present as a single contig in the gold standard and with annotated 16S sequences. The mean fraction of genes recovered by short-read assembler MEGAHIT is 36.9% 
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>For short-read assembly, MEGAHIT had the highest contiguity on the marine (NGA50 of 26,599)
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>HipMer (ranked first across metrics on the marine data)
   </td>
   <td>high genome fraction
<p>
There were selected 50 unique, public genomes present as a single contig in the gold standard and with annotated 16S sequences. The mean fraction of genes recovered by short-read assembler HipMer is 29.6% 
   </td>
   <td>HipMer had the fewest mm per 100 kb on the marine data set with 96
   </td>
   <td>
   </td>
   <td>high NGA50
   </td>
   <td>HipMer had the highest strain recall (14.4% on marine).
<p>
For the unique genome, the HipMer assembly had the highest strain recall (20.4%).
   </td>
   <td>achieved 100% precision.
<p>
For the unique genome, the highest strain precision (100%)
   </td>
  </tr>
  <tr>
   <td>GATB
   </td>
   <td>There were selected 50 unique, public genomes present as a single contig in the gold standard and with annotated 16S sequences. The hybrid assembler GATB recovered the most complete 16S sequences (mean recovered gene fraction 60.1%)
   </td>
   <td>The fewest mismatches (173) for hybrid assemblers were introduced by GATB
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>10.8% on marine
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>SPAdes
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>The SPAdes hybrid submission had NGA50 of 43,014
   </td>
   <td>On marine common genomes, SPAdes had the highest strain recall (8.7%)
   </td>
   <td>On marine common genomes, SPAdes had 96.7% strain precision
   </td>
  </tr>
  <tr>
   <td>Flye [<a href="https://www.biorxiv.org/lookup/google-scholar?link_type=googlescholar&gs_type=article&author[0]=M.+Kolmogorov&title=metaFlye:+scalable+long-read+metagenome+assembly+using+repeat+graphs&publication_year=2020&journal=Nat.+Methods&volume=17&pages=1103-1110">26</a>]
   </td>
   <td colspan="6" >Flye performed less well than other assemblers across most metrics on the marine data
   </td>
  </tr>
</table>


Detail box:

**_Performance of the assembly tools on the strain madness dataset_**:


<table>
  <tr>
   <td>
   </td>
   <td>Genome fraction (%)
   </td>
   <td>Mismatches per 100 kbp
   </td>
   <td>Misassemblies
   </td>
   <td>NGA50
   </td>
   <td>Strain recall
   </td>
   <td>Strain precision
   </td>
  </tr>
  <tr>
   <td>GSA
   </td>
   <td>90.8
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>155,979
   </td>
   <td>The upper bound for strain recall was 67.4%, as provided by the strain madness gold standard assembly
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>ABySS
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>Ray
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>A-STAR
   </td>
   <td>A-STAR excelled in terms of genome fraction on marine and strain madness data sets.
<p>
For strain madness common genomes >75% genome fraction
   </td>
   <td>On marine and strain-madness datasets A-STAR created more mismatches than others.
<p>
For strain madness common genomes &lt;0.5% mm
   </td>
   <td>On marine and strain-madness datasets A-STAR created more misassemblies than others.
   </td>
   <td>A-STAR had the highest contiguity for the strain madness data (13,008)
   </td>
   <td>For strain madness common genomes, A-STAR recovered the most, with 1.5% recall
   </td>
   <td>For strain madness common genomes - 23.1% strain precision
   </td>
  </tr>
  <tr>
   <td>OPERA-MS
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>MEGAHIT
   </td>
   <td>10.4
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>For short-read assembly, MEGAHIT had the highest contiguity on the strain madness data (NGA50 of 4,793)
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>HipMer (in the second place)
   </td>
   <td>For strain madness common genomes, HipMer recovered a lower genome fraction (4.1% versus 30.4% for A-STAR)
   </td>
   <td>For strain madness common genomes, HipMer created fewer mismatches (0.1%)
   </td>
   <td>For strain madness common genomes, HipMer created fewer misassemblies per genome (0.5)
   </td>
   <td>
   </td>
   <td>HipMer had the highest strain recall (3.2% on strain madness)
<p>
For strain madness common genomes, 0.8% strain recall
   </td>
   <td>For strain madness common genomes, 100% strain precision
   </td>
  </tr>
  <tr>
   <td>GATB ( ranked best)
   </td>
   <td>
   </td>
   <td>GATB had the fewest mm per 100 kb on the strain madness data set with 98 
   </td>
   <td>GATB created the fewest misassemblies for the strain madness data
   </td>
   <td>
   </td>
   <td>2.9% on strain madness
   </td>
   <td>achieved 100% precision
   </td>
  </tr>
  <tr>
   <td>SPAdes
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>Flye
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
</table>


Detail box:

**_Performance of the assembly tools on the plant-associated dataset_**:


<table>
  <tr>
   <td>
   </td>
   <td>Genome fraction (%)
   </td>
   <td>Mismatches per 100 kbp
   </td>
   <td>Misassemblies
   </td>
   <td>NGA50
   </td>
   <td>Strain recall
   </td>
   <td>Strain precision
   </td>
  </tr>
  <tr>
   <td>GSA
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>ABySS
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>Ray
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>A-STAR
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>OPERA-MS
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>High number of mismatches
   </td>
   <td>
   </td>
   <td>Extremely small strain recall
   </td>
   <td>OPERA-MS has the lowest precision
   </td>
  </tr>
  <tr>
   <td>MEGAHIT
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>HipMer (performed best across metrics)
   </td>
   <td>HipMer, HipMer2, and Flye have the best coverage.
   </td>
   <td>HipMer has less mismatches per 100 kbp compared to other assemblers. Followed by HipMer2
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>HipMer and HipMer2 have the highest strain recall
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>GATB
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>SPAdes
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>SPAdes has the lowest number of misassemblies
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>Flye (the second place)
   </td>
   <td>HipMer, HipMer2, and Flye have the best coverage.
   </td>
   <td>
   </td>
   <td>Flye_hybrid has higher number of misassemblies than others
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>Flye, HipMer, amd HipMer2 have the highest precision 
   </td>
  </tr>
</table>


Using the described metrics, the different tools were evaluated in the CAMI paper and aggregated in tables (Supplementary Tables 3-7) from CAMI tutorial [20]. 

In these tables there are also ranking scores of the tools shown for every statistic as well as overall ranking scores. Overall, ranking scores for every dataset are computed as a sum of all ranking scores across metrics. The average ranking score of both datasets are calculated as weighted average sum of ranking for both datasets. We created a table showing all ranking results from previous tables:

Detail box - [table ](https://docs.google.com/spreadsheets/d/1bKDz8uTu1nsuqmwLNLDM_vQ05lAnHedxbj1YVKlkQBc/edit#gid=455354696)



* We summed up ranking scores for all methods, common genomes, and unique genomes to extract the following  priority list for tools

    **Marine **dataset

    * <span style="text-decoration:underline;">With</span> tool versions
        1. HipMer
        2. metaSPAdes_v3.13.1
        3. metaSPAdes_v3.13.0
        4. ABySS
        5. Ray-Meta
        6. Megahit_v1.1.2
        7. SPAdes_v3.14-dev
* <span style="text-decoration:underline;">Without</span> tool versions
1. HipMer
2. **metaSPAdes**
3. **ABySS**
4. **Ray-Meta**
5. **Megahit**

**Strain madness** dataset



    * <span style="text-decoration:underline;">With</span> tool versions
        8. HipMer
        9. Megahit_v1.1.2
        10. SPAdes_v3.14-dev
        11. OPERA-MS
        12. Megahit_V1.2.7
* <span style="text-decoration:underline;">Without</span> tool versions
1. HipMer
2. **Megahit**
3. **SPAdes**
4. **OPERA**

**Plant-associated** dataset

There are no certain ranking tables among Supplementary tables for plant-associated dataset. However,  in the [paper](https://www.biorxiv.org/content/10.1101/2021.07.12.451567v1.full) [2] there is information related to tools performance on plant-associated dataset. We created the priority list of tools for the plant-associated dataset. 



1. (Meta)HipMer 
2. **(meta)Flye **
3. **(meta)SPAdes **

Since in this tutorial we have decided to focus on marine dataset it would be reasonable to reproduce CAMI2 assembly challenge using HipMer, metaSPAdes, ABySS, Ray-Meta, Megahit assemblers which performed better. As we know from our comparison Galaxy and CAMI2 analysis, metaSPAdes, ABySS, Megahit tools are available in Galaxy while Ray-Meta and HipMer are not.

Detail box - 

Moreover, we also selected Flye assembler, even though Flye is not in the Top-5 assemblers in CAMI2. Flye is a long-read assembler that works well for this task.. We do that to show the assembler launched on long-reads and not only on short-reads and hybrid.


### Tool Update

Due to the limited time for this project, we have decided to skip adding new tools into Galaxy. So, in this tutorial we do not run Hipmer and Ray-Meta.

But, in order to show the possibility of adding, updating and fixing tools and the process of how to do that.

All tools in Galaxy are kept up-to-date with their corresponding command-line tool. This is why in this tutorial we checked if the tool we want to use is up-to-date before using it. We then updated Megahit from 1.1.3.5 version to 1.2.9 (the latest version for that moment).

The other case it might be necessary to update the tool when we get an error on tool run or we find out that there is a bug. It can be an issue/bug in the tool itself, then we have to contact the tool developers. In other cases, it can be an issue/bug with the tool wrapper. If this is the case, you may be able to fix it by yourself. 

This was the case with the metaSPAdes tool. When we launched the metaSPAdes tool it failed with the error:

ln: failed to create symbolic link 'paired_reads1/FASTQ_splitter_on_data_X.fastq': File exists

The reason for this error was the same filenames with forward and reverse reads. We figure out that this error came from the FASTQ splitter tool, so we updated the FASTQ splitter. But we also updated spades tool to avoid similar cases in the future.

Updating or fixing  tools has the same concept as adding tools but is less time consuming. Therefore, in this tutorial we will limit ourselves to fixing a tool, but the concept will be the same for adding or updating a new tool.

Now we are going to update the tool in Galaxy by updating its wrapper. A wrapper is a XML which serves as a layer between the command-line tool itself and the Galaxy interface. This wrapper helps Galaxy users to have a visual interface and selection of the parameters of the command-line.  Please read our dedicated tutorial to learn more about Galaxy tools and the detailed instruction of how to build or update the tool is in [Planemo documentation](https://planemo.readthedocs.io/en/latest/writing_standalone.html).

To update or fix a tool, we need to change its XML. Tools in a Galaxy server are installed by the administrators of this Galaxy server using the Galaxy ToolShed. The Galaxy ToolShed is like an AppStore with the XML for all possible Galaxy tools. So to update a tool, we need to change its XML in the ToolShed. \
But usually, XML of tools are developed collaboratively inside GitHub repositories and then pushed to the ToolShed. So we need to find the development repository of the Galaxy tool to update it. 

Hands-on: Finding the development repository of a tool to update



    2. Open the fastqc splitter in Galaxy
    3. Click on “See in the ToolShed” on the top right of the tool form
        1. You will be redirected to the tool on the ToolShed
    4. Follow the link in “**Development repository”**
        2. **You have now the tool wrapper, i.e. its XML**
        3. Once we have the development repository setup, we can update the tool. We recommend using Planemo, the Galaxy software development kit, to help the process. More information is in [Planemo documentation](https://planemo.readthedocs.io/en/latest/writing_standalone.html).

Hands-on: Update the tools



1. Update lines corresponding to output files: add “Forward” or “Reverse” to the filenames respectively.
2. Update tests for this tool
3. Lint and test the tool using `planemo lint` and `planemo test`
4. Push the changes to a branch of your fork of the GitHub repository
5. Create a Pull Request
6. Modify the tool given the suggested changes
7. In case you add new tool - add the tool to yaml file (for usegalaxy.eu server - [tools_iuc.yaml](https://github.com/usegalaxy-eu/usegalaxy-eu-tools/blob/master/tools_iuc.yaml))
8. Once the Pull request is merged to the main branch at GitHub, the tool will be added (or updated with increased version) in Tool Shed (for usegalaxy.eu - [Galaxy | Tool Shed](https://toolshed.g2.bx.psu.edu/)) on weekend
9. Go to the Tool Shed after weekend and check the tool changes
10. Go to usegalaxy.eu (or to other Galaxy server if you used one) and check that tool was added/updated


### Run the tools

We run tools both on individual sample and on all pooled samples as it was done in CAMI2 challenge as well. 

It is important to use both approaches: individual assembly and co-assembly and benchmark them in order to get the intuition of more applicable technique for the certain dataset. In case of individual assembly we have the possibility to de-replicate genome after binning. It will not be done in this tutorial as the purpose here is only assembly CAMI2 challenge. However, you can get to know both these practices in the [Assembly tutorial](https://training.galaxyproject.org/training-material/topics/assembly/) (in the section about de-replication).

Co-assembly is more commonly used than individual assembly and then de-replication after binning. Nevertheless, de-replication approach should be taken into account while choosing the assembly method.

The runtime and memory usage are obviously different for not only these two scenarios, but also for the tools. Some assemblers are fast (like Megahit), some are slow (like Flye). 

While choosing parameters to run every tool with we oriented on parameters which have been used in CAMI2 challenge by its participants. However, we didn’t consider it as absolute truth as in metagenomic analysis it can be not very clear which tools and which parameters would be better for the best assembly performance. Then we had to be critical to what had been done in CAMI2. It was the basis of this tutorial/project, even though we made some corrections. Additionally, we take into account recommendations from tool documentation.

A Galaxy history was created for each tool and type of assembly (single or co-assembly). We named them such that it is easy to identify the tool and type from the mame, but we also Created a spreadsheet to list all histories along with their metadata (tool, etc).


#### Long reads


##### Flye

 MetaFlye (metagenome mode of Flye) is a scalable long-read metagenome assembly using repeat graphs [[26](https://www.biorxiv.org/lookup/google-scholar?link_type=googlescholar&gs_type=article&author[0]=M.+Kolmogorov&title=metaFlye:+scalable+long-read+metagenome+assembly+using+repeat+graphs&publication_year=2020&journal=Nat.+Methods&volume=17&pages=1103-1110)]. It is supposed to use long reads as input.


######  Hands-on: Co-assembly of long-read samples with Flye



1. Create a new history with the name “CAMI2 Flye v29 -meta -nano-raw pooled interleaved long collection”
2. Run Flye with parameters:
    1. **Input reads**: the collection of long reads from all samples
    2. **Mode**:	--nano-raw

            Comment box: Flye was reported in CAMI2 using parameter –raw-bio while datasets were sequenced with nanopore, thus we corrected parameters and used –nano-raw instead.

    3. **Perform metagenomic assembly**:	True


######  Hands-on: Individual assembly of long-read sample with Flye



3. Create a new history with the name “CAMI2 Flye v29 -meta -nano-raw pooled interleaved long collection”
4. Run Flye with parameters:
    4. **Input reads**: long reads from sample1
    5. **Mode**:	--nano-raw
    6. **Perform metagenomic assembly**:	True

We just launched Flye, first, for collection of long reads from all 10 samples to produce co-assembly. In the second Hands-on we launched Flye on individual sample. In the following hands-on for other tools for short-reads and hybrid we will need to run assemblers for both options as well.


#### Short reads


##### Hands-on : Co-assembly of short-reads samples with Megahit, Abyss, and MetaSPAdes

Comment box:

In the next steps we use the Megahit assembler. MEGAHIT is an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph [10]. As input it uses metagenomics samples as paired-end fastq files. As metagenomics input data short reads should be used.

1. Create a new history with the name “Co-assembly of short-reads samples with Megahit, Abyss, and MetaSPAdes”

2. Run MEGAHIT with parameters:



    7. **Select your input option**: 	interleaved
    8. **Interleaved-paired-end file(s**): the collection of short reads from all samples

3. Rename output to “Output - CAMI2 MEGAHIT v129 pooled interleaved short collection”

4. Run MEGAHIT with parameters:



    9. **Select your input option	**paired
    10. **Mate 1 input reads	**Output of FASTQ splitter: Forward
    11. **Mate 2 input reads	**Output of FASTQ splitter: Reverse
    12. **K-mer specification method	**klim_method
    13. **Minimum kmer size	**21
    14. **Maximum kmer size	**91
    15. **Increment of kmer size of each iteration	**12
5. Rename output to “Output - CAMI2 MEGAHIT v129 paired-end k-min21 k-max91 k-step12 pooled deinterleaved short collapsed”
1. Run MEGAHIT with parameters:
    16. **Select your input option**: 	interleaved
    17. **Interleaved-paired-end file(s**): the collection of short reads from all samples
    18. **K-mer specification method	**klim_method
    19. **Minimum kmer size	**21
    20. **Maximum kmer size	**91
    21. **Increment of kmer size of each iteration	**10
6. Rename output to “Output - CAMI2 MEGAHIT v129 k-min21 k-max91 k-step10 pooled interleaved short collection”

    Comment box


    In the next steps we use Abyss assembler. ABySS is a parallel assembler for short read sequence data [23]. It uses de novo short read assembly algorithms. Short reads should be used as input data for Abyss.

7. Run Abyss with parameters:
        1. **Type of paired-end datasets**:	paired_il
        2. **Interleaved paired-end reads**: collapsed collection of short reads from all samples	
        3. **K-mer length (in bp):**	41
1. Rename output to “Output - CAMI2 Abyss v234 k41 pooled interleaved short collapsed”
8. Run Abyss with parameters:
        4. **Type of paired-end datasets**:	paired_il
        5. **Interleaved paired-end reads**: collapsed collection of short reads from all samples	
        6. **K-mer length (in bp):**	96

10. Rename output to “Output - CAMI2 Abyss v234 k96 pooled interleaved short collapsed”

Comment box

In the next steps we use metaSPAdes assembler. MetaSPAdes is a versatile metagenomic assembler [9]. As input for metaspades it can accept short reads. However, there is an option to use additionally long reads besides short reads to produce hybrid input.



9. Run MetaSPAdes with parameters:
        1. **Pair-end reads input format**	paired_interlaced
        2. **FASTQ file(s): interlaced:**	collapsed collection of short reads from all samples	
        3. **Select k-mer detection option**	manual
        4. **K-mer size values**	21,33,55,77
10. Rename output to “Output - CAMI2 MetaSPAdes v3_15_3 k21-33-55-77 pooled interleaved short collapsed”
11. Run MetaSPAdes with parameters:
        5. **Pair-end reads input format**	: paired_interlaced
        6. **FASTQ file(s): interlaced:**	collapsed collection of short reads from all samples	
        7. **Arf - Nanopore reads**: collapsed collection of long reads from all samples
        8. **Select k-mer detection option**	manual
        9. **K-mer size values**	21,33,55,77
12. Rename output to “Output - CAMI2 MetaSPAdes v3_15_3 k21-33-55-77 pooled interleaved short collapsed additional nanopore long collapsed”


##### Hands-on : Individual assembly of short-reads samples with Megahit, Abyss, and MetaSPAdes

Comment box:

In the next steps we use the Megahit assembler. MEGAHIT is an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph [10]. As input it uses metagenomics samples as paired-end fastq files. As metagenomics input data short reads should be used.

1. Create a new history with the name “Individual assembly of short-reads samples with Megahit, Abyss, and MetaSPAdes”

2. Run MEGAHIT with parameters:



    22. **Select your input option**: 	interleaved
    23. **Interleaved-paired-end file(s**): dataset of short reads from sample1

3. Rename output to “Output - CAMI2 MEGAHIT v129 interleaved short sample1”

4. Run MEGAHIT with parameters:



    24. **Select your input option	**paired
    25. **Mate 1 input reads	**Output of FASTQ splitter: Forward
    26. **Mate 2 input reads	**Output of FASTQ splitter: Reverse
    27. **K-mer specification method	**klim_method
    28. **Minimum kmer size	**21
    29. **Maximum kmer size	**91
    30. **Increment of kmer size of each iteration	**12
13. Rename output to “Output - CAMI2 MEGAHIT v129 paired-end k-min21 k-max91 k-step12 deinterleaved short sample1”
2. Run MEGAHIT with parameters:
    31. **Select your input option**: 	interleaved
    32. **Interleaved-paired-end file(s**): dataset of short reads from sample1
    33. **K-mer specification method	**klim_method
    34. **Minimum kmer size	**21
    35. **Maximum kmer size	**91
    36. **Increment of kmer size of each iteration	**10
14. Rename output to “Output - CAMI2 MEGAHIT v129 k-min21 k-max91 k-step10 interleaved short sample1”

    Comment box


    In the next steps we use Abyss assembler. ABySS is a parallel assembler for short read sequence data [23]. It uses de novo short read assembly algorithms. Short reads should be used as input data for Abyss.

15. Run Abyss with parameters:
        7. **Type of paired-end datasets**:	paired_il
        8. **Interleaved paired-end reads**: dataset of short reads from sample1	
        9. **K-mer length (in bp):**	41
2. Rename output to “Output - CAMI2 Abyss v234 k41 interleaved short sample1”
16. Run Abyss with parameters:
        10. **Type of paired-end datasets**:	paired_il
        11. **Interleaved paired-end reads**: dataset of short reads from sample1	
        12. **K-mer length (in bp):**	96

10. Rename output to “Output - CAMI2 Abyss v234 k96 interleaved short sample1”

Comment box

In the next steps we use metaSPAdes assembler. MetaSPAdes is a versatile metagenomic assembler [9]. As input for metaspades it can accept short reads. However, there is an option to use additionally long reads besides short reads to produce hybrid input.



17. Run MetaSPAdes with parameters:
        10. **Pair-end reads input format**	paired_interlaced
        11. **FASTQ file(s): interlaced:**	dataset of short reads from sample1	
        12. **Select k-mer detection option**	manual
        13. **K-mer size values**	21,33,55,77
18. Rename output to “Output - CAMI2 MetaSPAdes v3_15_3 k21-33-55-77 interleaved short sample1”
19. Run MetaSPAdes with parameters:
        14. **Pair-end reads input format**	: paired_interlaced
        15. **FASTQ file(s): interlaced:**	dataset of short reads from sample1	
        16. **Arf - Nanopore reads**: collapsed collection of long reads from all samples
        17. **Select k-mer detection option**	manual
        18. **K-mer size values**	21,33,55,77
20. Rename output to “Output - CAMI2 MetaSPAdes v3_15_3 k21-33-55-77 interleaved short sample1 additional nanopore long sample1”

Comment box:

We got an error “Out of memory” for the first launches of metaspades. The memory limit in Galaxy was set to 200 Gb. After that,[ the memory for this tool was increased up to 250 GB](https://github.com/usegalaxy-eu/infrastructure-playbook/pull/358), which is the recommended default value according to [SPAdes documentation](http://cab.spbu.ru/files/release3.12.0/manual.html#sec3.2) (see “Advanced options” section). Only Galaxy administrators can solve this kind of problem with memory. If you need the memory increase for the specific tool you need to report an issue via Galaxy by clicking “View and report this error” (small “bug” icon).

However, it was not enough. According to [Supplementary Table 2](https://www.biorxiv.org/content/10.1101/2021.07.12.451567v1.supplementary-material) CAMI2 participants used the memory limit of 500 Gb. Since the memory limit in Galaxy for every tool is set once and it is used for every launch of this tool no matter how much memory it is needed. It seems to be unreasonable to set 500 Gb instead of the recommended 250, because then any spades job would request that amount of memory.

Other issues we encountered during this project are more general. Protracted technical works occured on the Galaxy platform, all tools were stopped and launched from scratch. That means we ran out of time. Before these technical works tools have been running for around one month and they were re-launched. However, we increased the CPU which helped us to get results sooner.


### Description of metrics

Several metrics are used to assess the quality of the assemblies.



* **Genome fraction (%): **% of reference bases covered by assembled contigs obtained by similarity-based mapping. 
    * The total number of aligned bases in the reference, divided by the genome size. A base in the reference genome is counted as aligned if at least one contig has at least one alignment to this base. Contigs from repeat regions may map to multiple places, and thus may be counted multiple times in this quantity.
* **Mismatches per 100 kbp**: number of mismatched bases in the contig-reference alignment (average per 100 kb). 
    * The average number of mismatches per 100 000 aligned bases. This metric does not distinguish between _single-nucleotide polymorphisms,_ which are true differences in the assembled genome versus the reference genome, and _single-nucleotide errors,_ which are due to errors in reads or errors in the assembly algorithm.
    * **Duplication ratio**: total number of aligned bases / genome fraction * reference length.
* **Number of misassemblies**: number of contigs which:
* contain a gap of more than 1kb;
* contain inserts of more than 1kb; or
* align to different genomes
* **The number of misassemblies**: number of positions in the assembled contigs where the left flanking sequence aligns over 1 kb away from the right flanking sequence on the reference, or they overlap by >1 kb, or the flanking sequences align on opposite strands or different chromosomes.
* **NGA50**:  metric for measuring the contiguity of an assembly. 
    * For each reference genome, all contigs aligned to it are sorted by size and the NGA50 for that genome is defined as the length of the contig cumulatively surpassing 50% genome fraction. If a genome is not covered to 50%, NGA50 is undefined. Since we report the average NGA50 over all genomes, it was set to 0 for genomes with less than 50% genome fraction.
    * There is a similar metric NG50 that defines the contig length such that using equal or longer length contigs produces 50% of the length of the reference genome. NGA50 is NG50 such that the lengths of aligned blocks are counted instead of contig lengths. It means that first we perform the sequence alignment of contigs to the reference genome. Next, contigs with misassemblies are split into aligned blocks and are aligned independently to distinct parts of the genome.
* **Strain recall: **raction of high-quality (more than 90% genome fraction and less than 100 mismatches per 100 kb) genome assemblies recovered for all ground truth genomes. Strain recall measures how many genomes are recovered with high genome fraction and few mismatches (mm).
* **Strain precision: **fraction of high-quality assemblies among all high genome fraction (more than 90%) assemblies. 
    * Strain precision assesses how accurately reference genomes are recovered, based on the fraction of correctly assembled high-quality, near-complete genomes (>90% genome fraction, &lt;0.1% mm) divided by the overall number of assembled, near-complete genomes (>90% genome fraction).


## Assess the assembly outputs

In CAMI2, different tools are used to evaluate the quality of the tools for the different challenges.

Detail box:The list of assessment tools with metadata is [here](https://docs.google.com/spreadsheets/d/1bKDz8uTu1nsuqmwLNLDM_vQ05lAnHedxbj1YVKlkQBc/edit#gid=15494775)

In this tutorial we focus on the assembly challenge and on the metrics presented before

Before running Quast as in CAMI to extract general statistics, we want to know the percentage of reads that were used to build the assemblie


### Extraction of % reads used in assemblies

To extract this information, we map the input reads on the assemblies.


#### Short-reads

We use Bowtie2 for mapping short-reads raw data to the assembly we got after usage of different short-reads assemblers (Megahit, MetaSPAdes, Abyss) in order to get the percentage of reads used for assembly.

Hands-on



1. Create a new history named "Short reads co-assemblies"
2. Drag and drop /results from the different reshistories with short read co-assemblies
3. Run Bowtie2 with parameters
    1. **Is this single or paired library**:	paired_interleaved
    2. **Interleaved FASTQ file**: raw dataset collection
    3. **Will you select a reference genome from your history or use a built-in index**?	history
    4. **Select reference genome**: Megahit/MetaSPAdes/Abyss output (If the assembler, like Abyss, has Contigs and Scaffolds as output, use Contigs preferably. Contigs do not contain gaps represented with multiple-X letters.)
    5. **Save the bowtie2 mapping statistics to the history**: True
    6. Inspect the generated output' . 

Hands-on



4. Create a new history named "Short reads individual assemblies"
5. Drag and drop /results from the different histories with short read individual assemblies
6. Run Bowtie2 with parameters
    7. **Is this single or paired library**:	paired_interleaved
    8. **Interleaved FASTQ file**: raw dataset collection
    9. **Will you select a reference genome from your history or use a built-in index**?	history
    10. **Select reference genome**: Megahit/MetaSPAdes/Abyss output (If the assembler, like Abyss, has Contigs and Scaffolds as output, use Contigs preferably. Contigs do not contain gaps represented with multiple-X letters.)
    11. **Save the bowtie2 mapping statistics to the history**: True
    12. Inspect the generated output .

 


#### Long reads

Bowtie2 should not be used to align long reads to resulting assemblies because, as it is stated in <span style="text-decoration:underline;">Bowtie2 documentation</span>, “Bowtie2 is geared toward aligning relatively short sequencing reads to long genomes. That said, it handles arbitrarily small reference sequences (e.g. amplicons) and very long reads (i.e. upwards of 10s or 100s of kilobases), though it is slower in those settings. It is optimised for the read lengths and error modes yielded by typical Illumina sequencers”.

For long-reads instead, _Minimap2_ aligner can be used for mapping the Nanopore-sequenced data.

 Hands-on



1. Create a new history named "Long reads co-assemblies"
2. Drag and drop /results from the different histories with long read Co-assemblies
3. Run Map with minimap2 with parameters:
    1. **Will you select a reference genome from your history or use a built-in index? **Use a genome from history and build index
    2. **Use the following dataset as the reference sequence: **output of Flye
    3. **Single or Paired-end reads: **paired interleaved
    4. **Select fastq dataset: **raw dataset collection

        This tool run produces one collection with the actual mapped reads for each Nanopore-sequenced sample. Unlike Bowtie2 it does not have an option to output mapping statistics directly. However, we can generate that information through an extra step.

4. Run Samtools stats tool with the following parameters:
    5. **_BAM file_**: the collection of mapped Nanopore-sequenced reads, output of Map with minimap2 tool
    1. _“Output”_: One single summary file
5. Inspect the generated output

Hands-on



6. Create a new history named "Long reads individual assemblies"
7. Drag and drop /results from the different histories with long read individual  assemblies
8. Run Map with minimap2 with parameters:
    6. **Will you select a reference genome from your history or use a built-in index? **Use a genome from history and build index
    7. **Use the following dataset as the reference sequence: **output of Flye
    8. **Single or Paired-end reads: **paired interleaved
    9. **Select fastq dataset: **raw dataset collection
9. Run Samtools stats tool with the following parameters:
    10. **_BAM file_**: the collection of mapped Nanopore-sequenced reads, output of Map with minimap2 tool
    2. _“Output”_: One single summary file
10. Inspect the generated output


### Extraction of general metrics

 Assemblies were evaluated with metaQUAST (metagenomics mode of QUAST) version 5.0.2.


#### Quast

Hands-on: Metaquast without provided reference genome



1. Move to "Short reads co-assemblies" history
7. Run Quast with parameters:
    13. **Contigs/scaffolds file**: outputs of the assemblies
    14. **Type of assembly**:	metagenome
    15. **Output files**:	
        1. HTML report 
        2. PDF report 
        3. Tabular reports 
        4. Log file
8. Move to "short-read individual assembly" history
9. Run Quast
10. Move to "long-read co-assembly" history
11. Run Quest
12. Move to "long-read individual assemblies" history
13. Run Quast  

Hands-on: Metaquast with provided reference genome



1. Move to "Short reads co-assemblies" history
14. Run Quast with parameters:
    16. **Contigs/scaffolds file**: outputs of the assemblies
    17. **Type of assembly**:	metagenome
    18. Reference genome: history
        5. GSA file from CAMI2 data

Comment box:

Working on this tutorial we encountered some issues. For example, in CAMI2 Quast was launched with Gold Standard Assembly as a reference genome. In Galaxy it took a long time and the tool is still running. Withsout GSA as a reference genome we got QUAST results faster



    19. **Output files**:	
        6. HTML report 
        7. PDF report 
        8. Tabular reports 
        9. Log file
15. Move to "short-read individual assembly" history
16. Run Quast
17. Move to "long-read co-assembly" history
18. Run Quest
19. Move to "long-read individual assemblies" history
20. Run Quast 

Comment box:

While running Quast, some of the generated outputs (statistics and log files) were empty while the html report was not. After investigation, it seems that the outputs were not correctly retrieved in the Galaxy wrapper. We then fixed the tool.

 


#### Aggregate all metrics

To generate a nice report with all metrics combined we use MultiQC. This tool is a good choice when you want to combine results and make them visually represented with different graphs. MultiQC generates a webpage combining reports.

Hands-on



3. Create new history “Benchmarking analysis of reproduced CAMI2”
4. Drag&drop Quast, bowtie2 and Samtools stat outputs from all histories for different assemblers
5. Run MultiQC tool with the following parameters:
    1. Which tool used to generate report: Quast
        13. Data: all assemblers outputs of Quast
    2. (for long reads) Add report
        14. Which tool used to generate report: samtools
        15. Data: Output of Samtool stats 
    3.  (for short reads) Add report for Bowtie 2 
        16. Which tool used to generate report: Bowtie 2
        17. Data: all assemblers outputs of Bowtie2 (for short reads)


### Extract computational metrics from Galaxy - memory used, runtime

In CAMI2 there were compared computational characteristics of different tools such as runtime and memory usage. There were created the following plots - [Fig 5](https://www.biorxiv.org/content/biorxiv/early/2021/07/12/2021.07.12.451567/F5.large.jpg?width=800&height=600&carousel=1)



<p id="gdcalert2" ><span style="color: red; font-weight: bold">>>>>>  gd2md-html alert: inline image link here (to images/image2.png). Store image on your image server and adjust path/filename/extension if necessary. </span><br>(<a href="#">Back to top</a>)(<a href="#gdcalert3">Next alert</a>)<br><span style="color: red; font-weight: bold">>>>>> </span></p>


![alt_text](images/image2.png "image_tooltip")


Fig 5


        … this section will be filled after the tools finally finished


    To export job metrics (memory usage and runtime) from Galaxy in csv format use [biobland ](https://bioblend.readthedocs.io/en/latest/api_docs/galaxy/all.html#bioblend.galaxy.jobs.JobsClient.get_metrics)


### Galaxy results

We launched a Jupyter notebook to analyse results we got and reproduce graphs from CAMI2 paper for 4 tools (Flye, Megahit, Abyss, MetaSPAdes) and additional GSA in order to compare our results with Gold Standard Assembly. Further steps for benchmarking analysis could be to plot both individual assembly quast results and co-assembly quast results.

The Github repository with jupyter notebook is [here](https://github.com/PlushZ/cami2-galaxy-benchmarking).


# Conclusion
{:.no_toc}


* Attach Multiqc 4 webpages:
    * Short reads (Megahit, Metaspades, abyss):
        * Co-assembly
        * Individual assembly
    * Long reads (Flye):
        * Co-assembly
        * Individual assembly
* Attach supplementary tables to report (samples’ metadata, tools description, etc.)
* Attach Galaxy history Workflow 
    * [Galaxy histories](https://docs.google.com/spreadsheets/d/1bKDz8uTu1nsuqmwLNLDM_vQ05lAnHedxbj1YVKlkQBc/edit#gid=915388998) - one hist per one tool+set of parameters + all histories were created for pooled (all 10 samples together) and just for sample1 - to get any result to test quast, bowtie, multiqc and jupyter 
    * Add the workflow from galaxy
    * List of tools used (collection, fastq splitter, assembly tools, quast, bowtie, minimap2, multiqc) 
