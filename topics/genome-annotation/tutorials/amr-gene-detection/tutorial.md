---
layout: tutorial_hands_on

title: Identification of AMR Genes
zenodo_link: 'https://zenodo.org/record/4534098'
questions:
- What species do I have and what is its sequence type?
- Which resistance genes are on my genome?
- Where are the genes located on my genome?
objectives:
- Assess presence of antimicrobial resistance genes
- Perform a species identification and MLST typing
- Search for resistance genes on the assembly
- Find a gene on your genome using Prokka + JBrowse
time_estimation: 2h
key_points:
- Annotation with Prokka is very easy
tags:
- illumina
- amr
- one-health

contributions:
  authorship:
  - bazante1
  - bebatut
  editing:
  - hexylena
  - bazante1
  - shiltemann
  - miaomiaozhou88
  funding:
  - avans-atgm
  - abromics

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

In this training you're going to make an assembly of data produced by
"Complete Genome Sequences of Eight Methicillin-Resistant
*Staphylococcus aureus* Strains Isolated from Patients in
Japan" from {% cite Hikichi_2019 %} which describes:

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