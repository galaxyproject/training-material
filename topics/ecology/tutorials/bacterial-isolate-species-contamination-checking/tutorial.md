---
layout: tutorial_hands_on

title: Checking expected species and contamination in bacterial isolate
zenodo_link: 'https://zenodo.org/record/10572227'
questions:
- What are the species in bacterial isolate sequencing data?
objectives:
- Run a series of tool to identify species in bacterial isolate sequencing data
- Visualize the species abundance
time_estimation: 1H
key_points:
- Kraken assigns taxons to sequences
- Bracken extracts species from Kraken assignations 
- Recentrifuge extracts stats and creates visualization from Kraken report
tags:
- illumina
- bacteria
- microgalaxy
level: Introductory
edam_ontology:
- topic_3673 # Whole genome sequencing
- topic_0622 # Genomics
- topic_3301 # Microbiology
- topic_3697 # Microbial ecology
contributions:
  authorship:
  - bebatut
  editing:
  - clsiguret
  funding:
  - abromics
---

Sequencing (determining of DNA/RNA nucleotide sequence) is used all over the world for all kinds of analysis. The product of these sequencers are reads, which are sequences of detected nucleotides. 

When sequencing bacterial isolates, it might be good to check if the expected species or strains can be identified in the data or if there is any contamination.

To illustrate the process, we take quality-controlled data of a bacterial genome (KUN1163 sample) from sequencing data produced in "Complete Genome Sequences of Eight Methicillin-Resistant *Staphylococcus aureus* Strains Isolated from Patients in Japan" ({% cite Hikichi_2019 %}). 

> Methicillin-resistant *Staphylococcus aureus* (MRSA) is a major pathogen
> causing nosocomial infections, and the clinical manifestations of MRSA
> range from asymptomatic colonization of the nasal mucosa to soft tissue
> infection to fulminant invasive disease. Here, we report the complete
> genome sequences of eight MRSA strains isolated from patients in Japan.
{: .quote cite="{% cite_url Hikichi_2019 %}"}

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Galaxy and data preparation

Any analysis should get its own Galaxy history. So let's start by creating a new one and get the data (forward and reverse quality-controlled sequences) into it.

> <hands-on-title>Prepare Galaxy and data</hands-on-title>
>
> 1. Create a new history for this analysis
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Rename the history
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
> 3. {% tool [Import](upload1) %} the quality-controlled sequences from [Zenodo]({{ page.zenodo_link }}) or from Galaxy shared data libraries:
>
>    ```
>    {{ page.zenodo_link }}/files/DRR187559_after_fastp_1.fastq.gz
>    {{ page.zenodo_link }}/files/DRR187559_after_fastp_2.fastq.gz
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
>
{: .hands_on}

# Taxonomic profiling

To find out which microorganisms are present, we will compare the reads of the sample to a reference database, i.e. sequences of known microorganisms stored in a database, using **Kraken2** ({% cite wood2019improved %}).

{% snippet topics/microbiome/faqs/kraken.md %}

For this tutorial, we will use the PlusPF database which contains the Standard (archaea, bacteria, viral, plasmid, human, UniVec_Core), protozoa and fungi data.

> <hands-on-title> Assign taxonomic labels with Kraken2</hands-on-title>
>
> 1. {% tool [Kraken2](toolshed.g2.bx.psu.edu/repos/iuc/kraken2/kraken2/2.1.1+galaxy1) %} with the following parameters:
>    - *"Single or paired reads"*: `Paired`
>        - {% icon param-file %} *"Forward strand"*: `DRR187559_after_fastp_1`
>        - {% icon param-file %} *"Reverse strand"*: `DRR187559_after_fastp_2`
>    - *"Minimum Base Quality"*: `10`
>    - In *"Create Report"*:
>        - *"Print a report with aggregrate counts/clade to file"*: `Yes`
>    - *"Select a Kraken2 database"*: `PlusPF-16`
>
{: .hands_on}

- **Classification**: tabular files with one line for each sequence classified by Kraken and 5 columns:

   1. `C`/`U`: a one letter indicating if the sequence classified or unclassified
   2. Sequence ID as in the input file
   3. NCBI taxonomy ID assigned to the sequence, or 0 if unclassified
   4. Length of sequence in bp (`read1|read2` for paired reads)
   5. A space-delimited list indicating the lowest common ancestor (LCA) mapping of each k-mer in the sequence     
      
      For example, `562:13 561:4 A:31 0:1 562:3` would indicate that:
      1. The first 13 k-mers mapped to taxonomy ID #562
      2. The next 4 k-mers mapped to taxonomy ID #561
      3. The next 31 k-mers contained an ambiguous nucleotide
      4. The next k-mer was not in the database
      5. The last 3 k-mers mapped to taxonomy ID #562     
      
      `|:|` indicates end of first read, start of second read for paired reads

   ```
   Column 1	Column 2	Column 3	Column 4	Column 5
   C 	DRR187559.1 	1280 	164|85 	0:1 1279:1 0:41 1279:10 0:5 1280:5 1279:1 0:1 1279:1 0:7 1279:5 0:2 1279:6 0:12 1279:5 0:19 1279:2 0:6 |:| 0:39 1279:2 0:6 1279:3 0:1
   C 	DRR187559.2 	1280 	70|198 	0:2 1279:5 0:29 |:| 0:52 1279:5 0:13 1279:3 0:23 1279:2 0:45 1280:1 0:9 1280:3 A:8
   C 	DRR187559.3 	1279 	106|73 	0:4 1279:3 0:36 1279:4 0:10 1279:5 0:3 1279:5 0:2 |:| 0:39
   C 	DRR187559.4 	1279 	121|189 	1279:6 0:17 1279:4 0:28 1279:2 0:30 |:| 0:7 1279:5 0:19 1279:5 0:20 1279:1 0:8 1279:1 0:1 1279:6 0:25 1279:2 0:44 A:11
   C 	DRR187559.5 	1279 	68|150 	1279:2 0:20 1279:3 0:9 |:| 0:10 1279:3 0:24 1279:2 0:9 1279:5 0:21 1279:5 0:20 1279:5 0:9 1279:1 0:2
   C 	DRR187559.6 	1280 	137|246 	0:2 1280:5 0:28 1279:1 0:28 1279:2 0:8 1279:2 0:23 1279:1 0:3 |:| 1279:1 0:2 1279:3 0:61 1279:2 0:14 1279:2 0:97 A:30 
   ```

   > <question-title></question-title>
   >
   > 1. Is the first sequence in the file classified or unclassified?
   > 2. What is the taxonomy ID assigned to the first classified sequence?
   > 3. What is the corresponding taxon?
   >
   > > <solution-title></solution-title>
   > > 1. classified
   > > 2. 1280
   > > 3. 1280 corresponds to [Staphylococcus aureus](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/tree/?taxon=1280) .
   > {: .solution}
   {: .question}

- **Report**: tabular files with one line per taxon and 6 columns or fields

   1. Percentage of fragments covered by the clade rooted at this taxon
   2. Number of fragments covered by the clade rooted at this taxon
   3. Number of fragments assigned directly to this taxon
   4. A rank code, indicating
      - (U)nclassified
      - (R)oot
      - (D)omain
      - (K)ingdom
      - (P)hylum
      - (C)lass
      - (O)rder
      - (F)amily
      - (G)enus, or
      - (S)pecies

      Taxa that are not at any of these 10 ranks have a rank code that is formed by using the rank code of the closest ancestor rank with a number indicating the distance from that rank. E.g., `G2` is a rank code indicating a taxon is between genus and species and the grandparent taxon is at the genus rank.

   5. NCBI taxonomic ID number
   6. Indented scientific name

   ```
   Column 1	Column 2	Column 3	Column 4	Column 5	Column 6
   0.24 	1065 	1065 	U 	0 	unclassified
   99.76 	450716 	14873 	R 	1 	root
   96.44 	435695 	2 	R1 	131567 	cellular organisms
   96.43 	435675 	3889 	D 	2 	Bacteria
   95.56 	431709 	78 	D1 	1783272 	Terrabacteria group
   95.53 	431578 	163 	P 	1239 	Firmicutes
   95.49 	431390 	4625 	C 	91061 	Bacilli
   94.38 	426383 	1436 	O 	1385 	Bacillales
   94.04 	424874 	2689 	F 	90964 	Staphylococcaceae
   93.44 	422124 	234829 	G 	1279 	Staphylococcus 
   ```

   > <question-title></question-title>
   >
   > 1. How many taxa have been found?
   > 2. What are the percentage on unclassified?
   > 3. What are the domains found?
   >
   > > <solution-title></solution-title>
   > >
   > > 1. 621, as the number of lines
   > > 2. 0.24%
   > > 3. Only Bacteria
   > {: .solution}
   >
   {: .question}

# Species extraction

In Kraken output, there are quite a lot of identified taxa with different levels. To extract the species level, we will use __Bracken__. 

__Bracken__ (Bayesian Reestimation of Abundance after Classification with Kraken) is a "simple and worthwile addition to Kraken for better abundance estimates" ({% cite Ye.2019 %}). Instead of only using proportions of classified reads, it takes a probabilistic approach to generate final abundance profiles. It works by re-distributing reads in the taxonomic tree: "Reads assigned to nodes above the species level are distributed down to the species nodes, while reads assigned at the strain level are re-distributed upward to their parent species" ({% cite Lu.2017 %}).

> <hands-on-title>Extract species with Bracken</hands-on-title>
>
> 1. {% tool [Bracken](toolshed.g2.bx.psu.edu/repos/iuc/bracken/est_abundance/2.9+galaxy0) %} with the following parameters:
>     - {% icon param-collection %} *"Kraken report file"*: **Report** output of **Kraken**
>     - *"Select a kmer distribution"*: `PlusPF`, same as for Kraken
>
>        It is important to choose the same database that you also chose for Kraken2
>
>     - *"Level"*: `Species`
>
{: .hands_on}

**Bracken** generates one output as a table with 7 columns:

- Taxon name
- Taxonomy ID
- Level ID (S=Species, G=Genus, O=Order, F=Family, P=Phylum, K=Kingdom)
- Kraken assigned reads
- Added reads with abundance re-estimation
- Total reads after abundance re-estimation
- Fraction of total reads

> <question-title></question-title>
>
> 1. How many species have been found?
> 2. Which the species has been the most identified? And in which proportion?
> 3. What are the other species?
>
> > <solution-title></solution-title>
> > 1. 51 (52 lines including 1 line with header)
> > 2. *Staphylococcus aureus* with 95.5% of the reads.
> > 3. Most of the other species are from *Staphylococcus* genus, so same as *Staphylococcus aureus*. The other species in really low proportion.
> {: .solution}
{: .question}

As expected *Staphylococcus aureus* represents most of the reads in the data.

# Contamination identification

To explore **Kraken** report and specially to detect more reliably minority organisms or contamination, we will use **Recentrifuge** ({% cite marti2019recentrifuge %}).

> <hands-on-title> Identify contamination </hands-on-title>
>
> 1. {% tool [Recentrifuge](toolshed.g2.bx.psu.edu/repos/iuc/recentrifuge/recentrifuge/1.12.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select taxonomy file tabular formated"*: report of **Kraken2** {% icon tool %}
>    - *"Type of input file (Centrifuge, CLARK, Generic, Kraken, LMAT)"*: `Kraken`
>    - In *"Database type"*:
>        - *"Cached database whith taxa ID"*: latest
>    - In *"Output options"*:
>        - *"Type of extra output to be generated (default on CSV)"*: `TSV`
>        - *"Remove the log file"*: `Yes`
>    - In *" Fine tuning of algorithm parameters"*:
>        - *"Strain level instead of species as the resolution limit for the robust contamination removal algorithm; use with caution, this is an experimental feature"*: `Yes`
>
{: .hands_on}

**Recentrifuge** generates 3 outputs:

- A statistic table with general statistics about assignations

   > <question-title></question-title>
   >
   > 1. How many sequences have been used?
   > 2. How many sequences have been classified?
   >
   > > <solution-title></solution-title>
   > > 1. 451,780
   > > 2. 450,715
   > {: .solution}
   {: .question}

- A data table with a report for each taxa

   > <question-title></question-title>
   >
   > 1. How many taxa have been kept?
   > 2. What is the lowest level in the data?
   >
   > > <solution-title></solution-title>
   > > 1. 185 (189 lines including 3 header lines)
   > > 2. The lowest level is strain.
   > {: .solution}
   {: .question}

- A HTML report with a Krona chart

   <iframe id="recentrifuge" src="recentrifuge_report.html" frameBorder="0" width="100%" height="900px"> ![Krona chart with multi-layered pie chart representing the community profile with in the center the higher taxonomy levels (i.e. domain) and on the exterior the more detailed ones (i.e. species)](./images/recentrifuge.png) </iframe>

   > <question-title></question-title>
   >
   > 1. What is the percentage of classified sequences?
   > 2. When clicking on *Staphylococcus aureus*, what can we say about the strains?
   > 3. Is there any contamination?
   >
   > > <solution-title></solution-title>
   > >
   > > 1. 99.8%
   > > 2. 99% of sequences assigned to *Staphylococcus aureus* are not assigned to any strains, probably because they are too similar to several strains. *Staphylococcus aureus* subsp. aureus JKD6159 is the strain with the most classified sequences, but only 0.3% of the sequences assigned to *Staphylococcus aureus*.
   > > 3. There is no sign of a possible contamination. Most sequences are classified to taxon on the *Staphylococcus aureus* taxonomy. Only 3% of the sequences are not classified to *Staphylococcus*.
   > >
   > {: .solution}
   >
   {: .question}

# Conclusion

In this tutorial, we checked bacterial isolate sequencing data for expected species and potential contamination.
