---
layout: tutorial_hands_on

title: Genome Assembly of MRSA using Illumina MiSeq Data
zenodo_link: 'https://zenodo.org/record/4534098'
questions:
- How to check the quality of the MiSeq data?
- How to perform an assembly with MiSeq data?
- What species do I have and what is its sequence type?
- Which resistance genes are on my genome?
- Where are the genes located on my genome?
objectives:
- Assess your data on quality and quantity
- Assemble a genome
- Assess your assembly quality
- Assess presence of antimicrobial resistance genes
- Perform a species identification and MLST typing
- Search for resistance genes on the assembly
- Find a gene on your genome using Prokka + JBrowse
time_estimation: 2h
key_points:
- Illumina assemblies can be done with Shovill, but the assemblies are not completely resolved
- Annotation with Prokka is very easy
tags:
- illumina
- assembly
- amr

contributions:
  authorship:
  - bazante1
  editing:
  - hexylena
  - bazante1
  - shiltemann
  - miaomiaozhou88
  funding:
  - avans-atgm

follow_up_training:
- type: "internal"
  topic_name: visualisation
  tutorials:
  - jbrowse
- type: "internal"
  topic_name: galaxy-interface
  tutorials:
  - history-to-workflow
---


# Introduction


In this training you're going to make an assembly of data produced by
"Complete Genome Sequences of Eight Methicillin-Resistant
*Staphylococcus aureus* Strains Isolated from Patients in
Japan" from {% cite Hikichi_2019 %} which describes:

> Methicillin-resistant *Staphylococcus aureus* (MRSA) is a major pathogen
causing nosocomial infections, and the clinical manifestations of MRSA
range from asymptomatic colonization of the nasal mucosa to soft tissue
infection to fulminant invasive disease. Here, we report the complete
genome sequences of eight MRSA strains isolated from patients in Japan.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Background

Sequencing (determining of DNA/RNA nucleotide sequence) is used all over
the world for all kinds of analysis. The product of these sequencers are
reads, which are sequences of detected nucleotides. Depending on the
technique these have specific lengths (30-500bp) or using Oxford
Nanopore Technologies sequencing have much longer variable lengths.

{% snippet faqs/galaxy/sequencing_illumina_miseq.md %}

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. {% tool [Import](upload1) %} the files from [Zenodo]({{ page.zenodo_link }}) or from the shared data library
>
>    ```
>    {{ page.zenodo_link }}/files/DRR187559_1.fastqsanger.bz2
>    {{ page.zenodo_link }}/files/DRR187559_2.fastqsanger.bz2
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Convert the datatype of this output to uncompress it
>
>    {% snippet faqs/galaxy/datasets_convert_datatype.md conversion="Convert compressed to uncompressed" %}
>
> 4. Rename the uncompressed dataset to just the sequence run ID: `DRR187559_1` and `DRR187559_2`
>
>    {% snippet faqs/galaxy/datasets_rename.md name="DRR187559_1" %}
>
> 5. Tag the dataset `#unfiltered`
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
> 6. **View** {% icon galaxy-eye %} the renamed file
>
>    > <question-title></question-title>
>    >
>    > 1. What are the 4 main features of each read in a fastq file.
>    > 2. What does the `_1` and `_2` mean in your filenames?
>    >
>    > > <solution-title></solution-title>
>    > > 1. The following:
>    > >
>    > >    -   A `@` followed by a name and sometimes information of the read
>    > >    -   A nucleotide sequence
>    > >    -   A `+` (optional followed by the name)
>    > >    -   The quality score per base of nucleotide sequence (Each symbol
>    > >        represents a quality score, which will be explained later)
>    > >
>    > > 2. Forward and reverse reads, by convention, are labelled `_1` and `_2`, but they might also be `_f`/`_r` or `_r1`/`_r2`.
>    > {: .solution}
>    {: .question}
{: .hands_on}

# Quality Control

When assessing the fastq files all bases had their own quality (or Phred score)
represented by symbols. You can read more in our dedicated [Quality Control
Tutorial]({% link topics/sequence-analysis/tutorials/quality-control/tutorial.md %}).

To assess the quality by hand would be too much work. That's why tools like
[NanoPlot](https://github.com/wdecoster/NanoPlot) or
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) are made,
which will generate a summary and plots of the data statistics. NanoPlot is
mainly used for long-read data, like ONT and PACBIO and FastQC for short read,
like Illumina and Sanger.


Before doing any assembly, the first questions you should ask about your input
reads include:

- What is the coverage of my genome?
- How good are my reads?
- Do I need to ask/perform for a new sequencing run?
- Is it suitable for the analysis I need to do?

> <hands-on-title>QC & Filtering</hands-on-title>
>
> 1. {% tool [FastQC](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.72+galaxy1) %} with the following parameters:
>    - {% icon param-files %} *"Short read data from your current history"*: `DRR187559_1` and `DRR187559_2`
>
> 2. Examine the output "WebPage" files
>
{: .hands_on}

DRR187559_1 | DRR187559_2
----------- | -----------
![FastQC plot showing reads that mostly stay in the green](../../images/mrsa/fastqc-1.png) | ![Same as previous plot, but the beginning of the reads are slightly better quality](../../images/mrsa/fastqc-2.png)

FastQC combines all quality statistics from all separate reads and combines them in plots. An important plot is the Per base sequence quality. Here you have the reads sequence length on the x-axes against the quality score (Phred-score) on the y-axis. The y-axis is divided in three sections: Green = good quality, Orange = mediocre quality and red = bad quality. For Illumina data it is normal that the first few bases are of some lower quality and how longer the reads get the worse the quality becomes. This is often due to signal decay or phasing during the sequencing run.

For each position, a boxplot is drawn with:

- the median value, represented by the central red line
- the inter-quartile range (25-75%), represented by the yellow box
- the 10% and 90% values in the upper and lower whiskers
- the mean quality, represented by the blue line

Depending on the analysis it could be possible that a certain quality or length
is needed. In this case we’re going to trim the data using Trimmomatic. Here we
will trim the start (LEADING) and end (TRAILING)of the reads if those fall
below a quality score of 20. Different trimming tools have different algorithms
for deciding when to cut but trimmomatic will cut based on the quality csore of
one base alone. Trimmomatic starts from each end, and as long as the base is
below 20, it will be cut until it reaches one greater than 20. A sliding window
trimming will be performed where if the average quality of 4 bases drops below
20, the read will be truncated there. And last we will also filter for reads
which are at least 30 bases long, anything shorter will complicate the assembly
and we will drop.

> <hands-on-title>Assessing the data quality of the trimmed reads and comparing to the input reads</hands-on-title>
>
> 1. {% tool [Trimmomatic](toolshed.g2.bx.psu.edu/repos/pjbriggs/trimmomatic/trimmomatic/0.38.1) %}  with the following parameters:
>    - *"Single-end or paired-end reads?"*: `Paired-end (two separate input files)`
>        - {% icon param-file %} *"Input FASTQ file (R1/first of pair)"*: `DRR187559_1 uncompressed`
>        - {% icon param-file %} *"Input FASTQ file (R2/second of pair)"*: `DRR187559_2 uncompressed`
>    - *"Perform initial ILLUMINACLIP step?"*: `Yes`
>    - In *"Trimmomatic Operation"*:
>        - {% icon param-repeat %} *"Insert Trimmomatic Operation"*
>            - *"Select Trimmomatic operation to perform"*: `Cut bases off the start of a read, if below a threshold quality (LEADING)`
>                - *"Minimum quality required to keep a base"*: `20`
>        - {% icon param-repeat %} *"Insert Trimmomatic Operation"*
>            - *"Select Trimmomatic operation to perform"*: `Cut bases off the end of a read, if below a threshold quality (TRAILING)`
>                - *"Minimum quality required to keep a base"*: `20`
>        - {% icon param-repeat %} *"Insert Trimmomatic Operation"*
>            - *"Select Trimmomatic operation to perform"*: `Sliding window trimming (SLIDINGWINDOW)`
>                - *"Number of bases to average across"*: `4`
>                - *"Average quality required"*: `20`
>        - {% icon param-repeat %} *"Insert Trimmomatic Operation"*
>            - *"Select Trimmomatic operation to perform"*: `Drop reads below a specified length (MINLEN)`
>                - *"Minimum length of reads to be kept"*: `30`
>
>    This produces two sets of output files, the R1/R2 paired, and R1/R2 unpaired.
>
>    - The first first two are the paired reads, the reads in each of those files match up.
>    - The second two are the unpaired reads where one of the two reads was discarded.
>
> 2. Edit the tags of the trimmomatic outputs. Remove the `#unfiltered` tag and add a new tag, `#filtered`
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
> 3. {% tool [FastQC](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.72+galaxy1) %} with the following parameters:
>    - {% icon param-files %} *"Short read data from your current history"*: `Trimmomatic on DRR187559_1 (R1 Paired)` and `Trimmomatic on DRR187559_2 (R2 Paired)`
>
> 4. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.9) %} with the following parameters:
>    - In *"Results"*:
>        - {% icon param-repeat %} *"Insert Results"*
>            - *"Which tool was used generate logs?"*: `FastQC`
>                - In *"FastQC output"*:
>                    - {% icon param-repeat %} *"Insert FastQC output"*
>                        - {% icon param-files %} *"FastQC output"*: select **all** of the `FastQC on data ... RawData` files. You should have 4; 2 FastQCs from before trimming, 2 from after trimming.
>
{: .hands_on}

MultiQC is a tool to combine multiple outputs in one clear and easy to read overview. It generates plots to easily visualize the different quality statistics of all your fastq files.

> <question-title></question-title>
>
> Looking at your `MultiQC on data ....: Webpage` output, answer the following questions:
>
> 1. How did the duplicate reads change from before to after?
> 2. How did the average read length change?
> 3. Did trimming improve the mean quality scores?
> 4. Did trimming affect the GC content?
> 5. Is this data ok to assemble? Do we need to re-sequence it?
>
> > <solution-title></solution-title>
> >
> > 1. There were slight decreases in duplicate reads (19.8% →  18%, 21.2% →  20.5%)
> > 2. Read lengths went down more significantly; 191 to 152 bp, and 221 to 192 bp. These are still probably fine for assembly but you'll see that MultiQC has marked the 152 bp result in yellow.
> > 3. Yes, definitely. It is most visible at the end of the reads, they stay completely in the "green zone", unlike the untrimmed data which falls into the yellow.
> > 4. No, it did not. If it had, that would be unexpected.
> > 5. This data looks OK. The number of shorter reads in R1 is not optimal but it should assemble something. Probably not the entire, closed genomes but something.
> >
> {: .solution}
{: .question}


# Assembly

When the quality of the reads is determined and the data is filtered and/or trimmed (like we did with trimmomatic)
an assembly can be made.

There are many tools that create assembly for short-read data, but in this
tutorial Shovill will be used. Shovill is a SPAdes based genome assembler,
improved to work faster and only for smaller (bacterial) genomes.

{% snippet faqs/galaxy/analysis_results_may_vary.md %}

> <hands-on-title>Assembly using Shovill</hands-on-title>
>
> 1. {% tool [Shovill](toolshed.g2.bx.psu.edu/repos/iuc/shovill/shovill/1.1.0+galaxy0) %} with the following parameters:
>    - *"Input reads type, collection or single library"*: `Paired End`
>        - {% icon param-file %} *"Forward reads (R1)"*: `Trimmomatic on DRR187559_1 uncompressed (R1 paired)`
>        - {% icon param-file %} *"Reverse reads (R2)"*: `Trimmomatic on DRR187559_2 uncompressed (R2 paired)`
>    - In *"Advanced options"*:
>        - *"Estimated genome size"*: `2914567`
>        - *"Minimum contig length"*: `500`
>
> 2. View {% icon galaxy-eye %} the `Contigs` file:
>
>    > <question-title></question-title>
>    > 1. How big is the first contig?
>    > 2. What is the coverage of your biggest (first) contig?
>    >
>    > > <solution-title></solution-title>
>    > > Both of these can be found in the header line of the Contigs produced by Shovill
>    > >
>    > > *NOTE: The results can differ from this example, because Shovill can differ a bit per assembly*
>    > >
>    > > 1. 419691 bp
>    > > 2. 14.1
>    > {: .solution}
>    {: .question}
>
> 3. {% tool [Bandage Image](toolshed.g2.bx.psu.edu/repos/iuc/bandage/bandage_image/0.8.1+galaxy2) %} with the following parameters:
>
>    - {% icon param-file %} *"Graphical Fragment Assembly"*: `Shovill on ... Contig Graph`
>
> 4. View {% icon galaxy-eye %} the assembly graph image
>
>    ![Bandage output showing a mess of a genome graph](../../images/mrsa/bandage-illumina.jpg)
{: .hands_on}


> <question-title></question-title>
> First read [this page in the Bandage wiki](https://github.com/rrwick/Bandage/wiki/Simple-example) to help understand what the graph means.
>
> 1. What do you think of this assembly? Is it useful? Is it good enough?
>
> > <solution-title></solution-title>
> >
> > 1. This is a very messy assembly, with a lot of potential paths through the sequence. We cannot feel confident in the output .fasta file (as it is much smaller than the expected 2.9Mbp). In real life we might consider doing a hybrid assembly with Nanopore or other long read data to help resolve these issues.
> >
> {: .solution}
{: .question}


## Genome Assembly Evaluation

[Quast](http://quast.bioinf.spbau.ru/) ({% cite Gurevich2013 %})
is a tool providing quality metrics for assemblies, and can also be used
to compare multiple assemblies. The tool can also take an optional
reference file as input, and will provide complementary metrics. QUAST
stands for QUality ASsessment Tool. With later updates gene annotation
also possible with QUAST.

> <hands-on-title>Quality Control of assembly using Quast</hands-on-title>
>
> 1. {% tool [Quast](toolshed.g2.bx.psu.edu/repos/iuc/quast/quast/5.0.2+galaxy1) %} with the following parameters:
>    - *"Use customized names for the input files?"*: `No, use dataset names`
>        - {% icon param-file %} *"Contigs/scaffolds file"*: `Shovill on ... Contigs`
>
> 2. View {% icon galaxy-eye %} the HTML report from QUAST
>
>    The Quast tool outputs assembly metrics as an html file with metrics and
>    graphs. The image below looks exceptionally boring. This is a good
>    thing, because each corner means one contig. A contig is the contiguous
>    sequence made by combining all separate reads in the assembly
>
>    ![Image showing the HTML output of quast including a table over conting information and a cumulative length graph with the contigs.](../../images/mrsa/quast-illumina.png)
>
{: .hands_on}


This fasta file contains 32 contigs, meaning the chromosome is separated over multiple contigs. These contigs can also contain (parts of) plasmids.

> <question-title></question-title>
>
> 1. What is you GC content?
>
> > <solution-title></solution-title>
> >
> > 1. The GC content for our assembly was 32.77% (for comparison [MRSA Isolate HC1335 Strain](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5174738/) GC% is 32.89%). This means the length and GC% of the assembly could be good!
> >
> {: .solution}
{: .question}

# Identification of AMR Genes

Because we are working with a MRSA we are curious to see which resistance genes are located on the genome or on the plasmid. To determine whether the contigs contain antimicrobial resistance (AMR) genes [staramr](https://github.com/phac-nml/staramr) can be used  **Staramr** scans bacterial genome contigs against both the **ResFinder** ({% cite Zankari2012 %}), **PointFinder** ({% cite Zankari2017 %}), and **PlasmidFinder** ({% cite Carattoli2014 %}) databases (used by the ResFinder webservice) and compiles a summary report of detected antimicrobial resistance genes.

> <hands-on-title>Run staramr</hands-on-title>
>
> 1. {% tool [staramr](toolshed.g2.bx.psu.edu/repos/nml/staramr/staramr_search/0.7.1+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"genomes"*: `Shovill on data ... Contigs`
>
>    There are 7 different output files produced by **staramr** tool:
>
>    File                  | Contents
>    --------------------- | ---
>    summary.tsv           | A summary of all detected AMR genes/mutations in each genome, one genome per line.
>    detailed_summary.tsv  | A detailed summary of all detected AMR genes/mutations of each genome, one line per feature and multiple lines per genome.
>    resfinder.tsv         | A tabular file of each AMR gene and additional BLAST information from the **ResFinder** database, one gene per line.
>    Plasmidfinder.tsv     | A tabular file of each plastid sequences with additional BLAST information from the **PlasmidFinder** database, one sequence per line.
>    settings.txt          | The command-line, database versions, and other settings used to run staramr.
>    mlst.tsv              | A tabular file of the found loci per genome with its specified MLST scheme.
>    results.xlsx          | An Excel spreadsheet containing the previous 4 files as separate worksheets.
>
> 2. View {% icon galaxy-eye %} the detailed_summary.tsv file
>    - In this example the ST-typing could not be obtained. Multi-locus sequence type (MLST) is based on specific locus/alleles, which is sometimes hard to determine with error rich sequence data (like NanoPore).
>    - For the plasmid and resistance results the identity, overlap, length and the location on the contig can be found here.
>    - Multiple rep sequences are located on the second contig. (See "plasmid typing for gram-positive bacteria" {% cite Lozano_2012 %} for more information)
>    - Multiple resistance genes can be found on both contig 1 and contig 2.
>    - In the last column there are "Accession" numbers. These are references to NCBI, and you can search for these numbers there. E.g. [M13771](https://www.ncbi.nlm.nih.gov/nuccore/M13771)
>
{: .hands_on}


## CARD database

To get more information about these antibiotic resistant genes, you can check the [CARD database](https://card.mcmaster.ca) (**C**omprehensive **A**ntibiotic **R**esistance **D**atabase) ({% cite Jia2016 %})

![Screenshot of mecA sequence in CARD database with lots of metadata](../../images/mrsa/card.png "Screenshot of the CARD database interface. CARD gives information about the antibiotic resistance genes, as well as links to relevant publications.")

CARD can be very helpful to check all the resistance genes and check if
it is logical to find the resistance gene in a specific bacteria.

> <question-title></question-title>
>
> 1. To what family does [mecA](https://card.mcmaster.ca/ontology/36911) belong?
> 2. Do you expect to find this gene in this MRSA strain and why?
> 3. Is the accession number of the entry related to the accession reported by staramr?
>
> > <solution-title></solution-title>
> >
> > 1. [Methicillin resistant PBP2](https://card.mcmaster.ca/ontology/37589)
> > 2. The strain we use is a Methicillin(multi) resistant Staphylococcus aureus. As `mecA` has a perfect resistome mach with *S. aureus*, and the AMR Gene Family is methicillin resistant PBP2, we expect to see mecA in MRSA.
> > 3. No, these are completely unrelated. Unfortunately this is a **very** common issue in bioinformatics. Everyone builds their own numbering system for entries in their database (usually calling them 'accessions'), and then someone else needs to build a service to link these databases.
> >
> {: .solution}
{: .question}

## Gene annotation using Prokka

[Prokka](https://github.com/tseemann/prokka/blob/master/README.md) is a
tool software tool to rapidly annotate bacterial, archaeal and viral
genomes. Prokka will be used on your own made genome (assembly). Prokka
will try to annotate the bacteria based on related species and starting
codons can be chosen or default of the species can be used.

[JBrowse](https://jbrowse.org/docs/tutorial.html) is used to visualize
your genome file and merge multiple outputs.\
In this case you will use your assembly as your reference and the output
from prokka as an information track.


> <hands-on-title>Annotating the Genome</hands-on-title>
>
> 1. {% tool [Prokka](toolshed.g2.bx.psu.edu/repos/crs4/prokka/prokka/1.14.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Contigs to annotate"*: `Shovill on data ... Contigs`
>    - *"Genus name (--genus)"*: `staphylococcus `
>    - *"Species name (--species)"*: `aureus`
>    - *"Kingdom (--kingdom)"*: `Bacteria`
>    - *"Additional outputs"*: Select only the "Annotation in GFF3 format contaianing both sequences and annotations"
>
> 2. {% tool [Select lines that match an expression](Grep1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `staramr on data .. detailed_summary.tsv`
>    - *"that"*: Matching
>    - *"the pattern"*: `[0-9]+\.[0-9]+\t`
>
>    This will select lines with a decimal value (###.##) followed by a tab character, the column separator in Galaxy. As a result, any lines without an identity value will be filtered out.
>
> 3. {% tool [Table to GFF3](toolshed.g2.bx.psu.edu/repos/iuc/tbl2gff3/tbl2gff3/1.2) %}
>    - {% icon param-file %} *"Table"*: the output of the above **Select lines** {% icon tool %} step.
>    - *"Record ID column or value"*: `8`
>    - *"Feature start column or value"*: `9`
>    - *"Feature end column or value"*: `10`
>    - *"Feature score column or value"*: `5`
>    - *"Feature source column or value"*: `3`
>    - {% icon param-repeat %} *"Insert Qualifiers"*
>        - *"Name"*: `Name`
>        - *"Qualifier value column or raw text"*: `2`
>    - {% icon param-repeat %} *"Insert Qualifiers"*
>        - *"Name"*: `phenotype`
>        - *"Qualifier value column or raw text"*: `4`
>    - {% icon param-repeat %} *"Insert Qualifiers"*
>        - *"Name"*: `accession`
>        - *"Qualifier value column or raw text"*: `11`
>
> 4. {% tool [Bowtie2](toolshed.g2.bx.psu.edu/repos/devteam/bowtie2/bowtie2/2.3.4.3+galaxy0) %} with the following parameters:
>    - *"Is this single or paired library"*: `Paired-end`
>        - {% icon param-file %} *"FASTA/Q file #1"*: `Trimmomatic on DRR187559_1 uncompressed (R1 paired)` (output of **Trimmomatic** {% icon tool %})
>        - {% icon param-file %} *"FASTA/Q file #2"*: `Trimmomatic on DRR187559_2 uncompressed (R2 paired)` (output of **Trimmomatic** {% icon tool %})
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from the history and build index`
>        - {% icon param-file %} *"Select reference genome"*: `contigs` (output of **Shovill** {% icon tool %})
>    - *"Select analysis mode"*: `1: Default setting only` - You should have a look at the non-default parameters and try to understand them, but they don't need to be changed currently.
>    - *"Save the bowtie2 mapping statistics to the history"*: `Yes`
>
> 5. {% tool [JBrowse](toolshed.g2.bx.psu.edu/repos/iuc/jbrowse/jbrowse/1.16.9+galaxy0) %} with the following parameters:
>    - *"Reference genome to display"*: `Use a genome from history`
>        - {% icon param-file %} *"Select the reference genome"*: `Shovill on ...: Contigs` (output of **Shovill assembly** {% icon tool %})
>    - *"Genetic Code"*: `11. The Bacterial, Archaeal and Plant Plastid Code`
>    - In *"Track Group"*:
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - *"Track Category"*: `Prokka`
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `GFF/GFF3/BED Features`
>                        - {% icon param-file %} *"GFF/GFF3/BED Track Data"*: `Prokka on data ...: gff` (output of **Prokka** {% icon tool %})
>                        - *"JBrowse Track Type [Advanced]"*: `Neat Canvas Features`
>                        - *"Track Visibility"*: `On for new users`
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - *"Track Category"*: `AMR`
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `GFF/GFF3/BED Features`
>                        - {% icon param-file %} *"GFF/GFF3/BED Track Data"*: `Table to GFF3 on ...`, the output of the table to gff3 step
>                        - *"JBrowse Track Type [Advanced]"*: `Neat Canvas Features`
>                        - *"Track Visibility"*: `On for new users`
>        - {% icon param-repeat %} *"Insert Track Group"*
>            - *"Track Category"*: `Sequencing`
>            - In *"Annotation Track"*:
>                - {% icon param-repeat %} *"Insert Annotation Track"*
>                    - *"Track Type"*: `BAM Pileups`
>                        - {% icon param-file %} *"BAM Track Data"*: Bowtie2's output
>                        - *"Autogenerate SNP Track"*: `Yes`
>
> 3. View the output of JBrowse
>
{: .hands_on}

In the output of the JBrowse you can view the mapped reads and the found genes against the reference genome. With the search tools you can easily find genes of interest. JBrowse can handle many inputs and can be very useful. Using the bowtie2 mapping output, low coverage regions can be detected. This SNP detection can also give a clear view of where the data was less reliable or where variations were located.

If it takes too long to build the JBrowse instance, you can view an embedded one here. (**Warning**: feature name search will not work.)

{% snippet topics/visualisation/faqs/visualizations_jbrowse.html datadir="data" loc="contig00018:5488..27391" tracks="DNA,14e421a8469880793f2a8d774224e10d_0,6851a9d3f5d62263e4e4b34f1513980c_0" %}

# Conclusion

