---
layout: tutorial_hands_on

title: Metabarcoding/eDNA through Obitools
zenodo_link: https://zenodo.org/record/5932108/files/wolf_tutorial.zip?download=1
questions:
- how to analyze DNA metabarcoding / eDNA data produced on Illumina sequencers using the OBITools?
objectives:
- Deal with paired-end data to create consensus sequences
- Clean, filter and anlayse data to obtain strong results
time_estimation: 1H
key_points:
- From raw reads you can process, clean and filter data to obtain a list of species from environmental DNA (eDNA) samples.
contributors:
- colineroyaux
- onorvez
- obitools_team
- yvanlebras

---


# Introduction
{:.no_toc}

The data used in this tutorial correspond to those described in the [official OBITools tutorial](https://pythonhosted.org/OBITools/wolves.html) and show how to analyse four wolf scats, using the protocol published in Shehzad et al. (2012) for assessing carnivore diet. After extracting DNA from the faeces, the DNA amplifications were carried out using the primers TTAGATACCCCACTATGC and TAGAACAGGCTCCTCTAG amplifiying the 12S-V5 region (Riaz et al. 2011), together with a wolf blocking oligonucleotide.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Tutorial based on the official OBITools one
Based on this [OBITools official tutorial](https://pythonhosted.org/OBITools/wolves.html), you will learn here how to analyze DNA metabarcoding data produced on Illumina sequencers using:
 * the OBITools on Galaxy
 * some classical Galaxy tools

The data used in this tutorial correspond to the analysis of four wolf scats, using the protocol published in Shehzad et al. (2012) for assessing carnivore diet. After extracting DNA from the faeces, the DNA amplifications were carried out using the primers TTAGATACCCCACTATGC and TAGAACAGGCTCCTCTAG amplifiying the 12S-V5 region (Riaz et al. 2011), together with a wolf blocking oligonucleotide.

It is always a good idea to have a look at the intermediate results or to evaluate the best parameter for each step. Some commands are designed for that purpose, for example you can use :

- obicount to count the number of sequence records in a file
- obihead and obitail to view the first or last sequence records of a file
- obistat to get some basic statistics (count, mean, standard deviation) on the attributes (key=value combinations) in the header of each sequence record (see The extended OBITools fasta format in the fasta format description)
- any Galaxy tools corresponding to classical unix command such as less, awk, sort, wc to check your files.

# Manage input data
The data needed to run the tutorial are the following:

- fastq files resulting of a GA IIx (Illumina) paired-end (2 x 108 bp) sequencing assay of DNA extracted and amplified from four wolf faeces:
    - wolf_F.fastq
    - wolf_R.fastq
- the file describing the primers and tags used for all samples sequenced:
    - wolf_diet_ngsfilter.txt The tags correspond to short and specific sequences added on the 5’ end of each primer to distinguish the different samples
- the file containing the reference database in a fasta format:
    -db_v05_r117.fasta This reference database has been extracted from the release 117 of EMBL using ecoPCR

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the zip archive containing input files from [Zenodo](https://zenodo.org/record/5932108/files/wolf_tutorial.zip?download=1) 
>
>    ```
>    https://zenodo.org/record/5932108/files/wolf_tutorial.zip?download=1
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the dataset, here a zip archive, if needed
> 4. Check that the datatype is `zip`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>

## **Unzip** the downloaded archive

> ### {% icon hands_on %} Hands-on: Unzip the downladed .zip archive and prepare unzipped files to be used by OBITools
>
> 1. {% tool [Unzip](toolshed.g2.bx.psu.edu/repos/imgteam/unzip/unzip/0.2) %} with the following parameters:
>    - *"Extract single file"*: `All files`
>
>    > ### {% icon comment %} Comment
>    >
>    > To work properly, this unzip Galaxy tool is waiting "simple" archive as input, this means without sub directory.
>    {: .comment}
>
> 2. Add to each datafile a tag and/or modify names (*optional*)
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
> 3. Unhide all dataset from the resulting data collection so you can use these files independently.
>
> 4. Modify datatype from txt to tabular for the `wolf_diet_ngsfilter` dataset
{: .hands_on}


> ### {% icon question %} Questions
>
> 1. Why do we need to unhide manually datasets from the data collection?
>
> > ### {% icon solution %} Solution
> >
> > 1. Data collection is a functionality often used to deal with multiple datasets on the same format who can be analysed in batch mode. Here, the data collection is populated with heterogenous datafiles, coming from an archive. We thus need to treat separately each dataset of the collection, and to do so, we need to unhide corresponding datasets from the history, as datasets insides collections ar "just" like "symbolic link" to "classical" history datasets hidden by default.
> >
> {: .solution}
>
{: .question}




# Use OBITools

OBITools is a set of programs specifically designed for analyzing NGS data in a DNA metabarcoding context, taking into account taxonomic information. It is distributed as an open source software available on the following website: http://metabarcoding.org/obitools.

Citation: Boyer F., Mercier C., Bonin A., Taberlet P., Coissac E. (2016) OBITools: [a Unix-inspired software package for DNA metabarcoding](https://pubmed.ncbi.nlm.nih.gov/25959493/). Molecular Ecology Resources.

The OBITools commands consider a sequence record as an entity composed of five distinct elements. Two of them are mandatory, the identifier (id) and the DNA or protein sequence itself. The id is a single word composed of characters, digits, and other symbols like dots or underscores excluding spaces. Formally, the ids should be unique within a dataset and should identify each sequence record unambiguously, but only a few OBITools actually rely on this property. The sequence is an ordered set of characters corresponding to nucleotides or amino-acids according to the International Union of Pure and Applied Chemistry (IUPAC) nomenclature (Cornish-Bowden 1985). The three other elements composing a sequence record are optional. They consist in a sequence definition, a quality vector, and a set of attributes. The sequence definition is a free text describing the sequence briefly. The quality vector associates a quality score to each nucleotide or amino-acid. Usually this quality score is the result of the base-calling process by the sequencer. The last element is a set of attributes qualifying the sequence, each attribute being described by a key=value pair. The set of attributes is the central concept of the OBITools system. When an OBITools command is run on the sequence records included in a dataset, the result of the computation often consist in the addition of new attributes completing the annotation of each sequence record. This strategy of sequence annotation allows the OBITools to return their results as a new sequence record file that can be used as the input of another OBITools program, ultimately creating complex pipelines.


## Sub-step with **illuminapairedend**

> ### {% icon hands_on %} Hands-on: Recover consensus sequences from overlapping forward and reverse reads.
> 
> When using the result of a paired-end sequencing assay with supposedly overlapping forward and reverse reads, the first step is to recover the assembled sequence.
>
> The forward and reverse reads of the same fragment are at the same line position in the two fastq files obtained after sequencing. Based on these two files, the assembly of the forward and reverse reads is done with the illuminapairedend utility that aligns the two reads and returns the reconstructed sequence.
>
> 1. {% tool [illuminapairedend](toolshed.g2.bx.psu.edu/repos/iuc/obi_grep/obi_illuminapairedend/1.2.13) %} with the following parameters:
>    - *"Read from file"*: `wolf_F` for the 3p file
>    - *"Read from file"*: `wolfRF` for the 5p file
>    - *"minimum score for keeping aligment"*: `40.0`
>
>    > ### {% icon comment %} Comment
>    >
>    > Sequence records corresponding to the same read pair must be in the same order in the two files !
>    > 
>    > If the alignment score is below the defined score, here 40, the forward and reverse reads are not aligned but concatenated, and the value of the mode attribute in the sequence header is set to joined instead of alignment
>    {: .comment}
>
{: .hands_on}



## Sub-step with **obigrep**

> ### {% icon hands_on %} Hands-on: Remove unaligned sequence records
> We here use the value of the mode attribute in the sequence header to discard sequences not "joined" (see explanation about this mode on the previous step)
>
> 1. {% tool [obigrep](toolshed.g2.bx.psu.edu/repos/iuc/obi_grep/obi_grep/1.2.13) %} with the following parameters:
>    - *"Choose the sequence record selection option"*: `predicat`
>        - *"Python boolean expression to be evaluated for each sequence record."*: `mode!=joined`
>
>
>    > ### {% icon comment %} Comment
>    >
>    > The obigrep command is in some way analog to the standard Unix grep command. It selects a subset of sequence records from a sequence file.
>    > 
>    > A sequence record is a complex object composed of an identifier, a set of attributes (key=value), a definition, and the sequence itself.
>    > 
>    > Instead of working text line by text line as the standard Unix tool, selection is done sequence record by sequence record. A large set of options allows refining selection on any of the sequence record elements.
>    > 
>    > Moreover obigrep allows specifying simultaneously several conditions (that take the value TRUE or FALSE) and only the sequence records that fulfill all the conditions (all conditions are TRUE) are selected.
>    > 
>    > Sequence record selection options : 
>    > * sequence : Regular expression pattern to be tested against the sequence itself. ex: GAATTC
>    > * definition : Regular expression pattern to be tested against the definition of the sequence record. ex: [Cc]hloroplast
>    > * identifier : Regular expression pattern to be tested against the identifier of the sequence record. ex: ^GH
>    > * idlist : points to a text file containing the list of sequence record identifiers to be selected.
>    > * attribute : Regular expression pattern matched against the attributes of the sequence record. the value of this attribute is of the form : key:regular_pattern. ex:'family_name:Asteraceae'
>    > * hasattribute : Selects sequence records having an attribute whose key = KEY.
>    > * predicat : Python boolean expression to be evaluated for each sequence record. The attribute keys defined for each sequence record can be used in the expression as variable names. An extra variable named ‘sequence’ refers to the sequence record itself. ex: mode!="joined"
>    > * lmax : Keeps sequence records whose sequence length is equal or shorter than lmax. ex : 100
>    > * lmin : Selects sequence records whose sequence length is equal or longer than lmin. ex : 100
>    {: .comment}
>
{: .hands_on}


> ### {% icon question %} Questions
>
> 1. How do you verify the operation is successfull?
> 2. How many sequences are kept? Discarded?
>
> > ### {% icon solution %} Solution
> >
> > 1. you can search in the input file content the presence of `mode=joined` and same on the output file (just clicking the eye to visualize the content of each file and typing CTRL+C for example to search `mode=joined` in the file, or using a regex Galaxy tool for example). You can also at least look at the size of the output file, if smaller than input file, this is a first good indication.
> > 2. You can use a Galaxy tool like `Line/Word/Character count of a dataset` to count the number of lines of each dataset (input and output of obigrep) and divided by 4 (as in a FastQ file, each sequence is represented by a block of 4 lines). 45 276 sequences for input file. 44 717 for output file. Thus 559 sequences discarded.
> >
> {: .solution}
>
{: .question}

## Sub-step with **NGSfilter**

> ### {% icon hands_on %} Hands-on: Assigns sequence records to the corresponding experiment/sample based on DNA tags and primers
>
> 1. {% tool [NGSfilter](toolshed.g2.bx.psu.edu/repos/iuc/obi_ngsfilter/obi_ngsfilter/1.2.13) %} with the following parameters:
>    - *"Parameter file"*: `wolf_diet_ngsfilter`
>    - *"Read from file"*: `obigrep output`
>    - *"Number of errors allowed for matching primers"*: `2`
>    - *"Output data type"*: `fastq`
>
>
>    > ### {% icon comment %} Comment
>    >
>    > Each sequence record is assigned to its corresponding sample and marker using the data provided in a text file (here wolf_diet_ngsfilter.txt). This text file contains one line per sample, with the name of the experiment (several experiments can be included in the same file), the name of the tags (for example: aattaac if the same tag has been used on each extremity of the PCR products, or aattaac:gaagtag if the tags were different), the sequence of the forward primer, the sequence of the reverse primer, the letter T or F for sample identification using the forward primer and tag only or using both primers and both tags, respectively.
>    {: .comment}
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. How many sequences are not assigned?
>
> > ### {% icon solution %} Solution
> >
> > 1. 1391
> >
> {: .solution}
>
{: .question}

## Sub-step with **obiuniq**

> ### {% icon hands_on %} Hands-on: Groups together sequence records
>
> 1. {% tool [obiuniq](toolshed.g2.bx.psu.edu/repos/iuc/obi_uniq/obi_uniq/1.2.13) %} with the following parameters:
>    - *"Input sequences file"*: `Trimmed and annotated file by NGSfilter`
>    
>    - *"Attribute to merge"*: `sample`
>    
>    - *"Use specific option"*: `merge`
>
>
>    > ### {% icon comment %} Comment
>    >
>    > The same DNA molecule can be sequenced several times. In order to reduce both file size and computations time, and to get easier interpretable results, it is convenient to work with unique sequences instead of reads. To dereplicate such reads into unique sequences, we use the obiuniq command.
>    > Definition: Dereplicate reads into unique sequences
>    > * compare all the reads in a data set to each other
>    > * group strictly identical reads together
>    > * output the sequence for each group and its count in the original dataset (in this way, all duplicated reads are removed)
>    > *Definition adapted from Seguritan and Rohwer (2001)*
>    {: .comment}
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. How many sequences you had and how many you finally obtain?
>
> > ### {% icon solution %} Solution
> >
> > 1. From 43 326 to 3 962
> >
> {: .solution}
>
{: .question}

## Sub-step with **FastQC**

> ### {% icon hands_on %} Hands-on: Check quality of your data before and after analysis
>
> 1. {% tool [FastQC](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.73+galaxy0) %} with the following parameters:
>    - *"Raw read data from your current history"*: `wolf-F`, `wolf-R` and `obigrep output file`
>
>    > ### {% icon comment %} Comment
>    >
>    > To select more than one input dataset and execute the tool in parallel on multiple files, you have to select the `Multiple datasets` mode.
>    > 
>    {: .comment}
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. Do you know why Per base sequence quality is ok with obigrep output file and not original F and R files?
>
> > ### {% icon solution %} Solution
> >
> > 1. Raw F and R files have classical decrease of sequencing quality with the length of sequences. On obigrep output dataset, Forward and Reverse sequences are merged so nucleotides with low sequencing quality from the end of forward reads are replaced by the high sequencing quality from the start of reverse corresponding reads.
> >
> {: .solution}
>
{: .question}

## Sub-step with **obiannotate**

> ### {% icon hands_on %} Hands-on: Adds/Edits sequence record annotations
>
> 1. {% tool [obiannotate](toolshed.g2.bx.psu.edu/repos/iuc/obi_annotate/obi_annotate/1.2.13) %} with the following parameters:
>    - In *"Keep only attribute with key"*:
>        - *"key"*: `count`
>        - *"if you want to specify a second key"*: `merged_sample`

>
>    > ### {% icon comment %} Comment
>    >
>    > obiannotate is the command that allows adding/modifying/removing annotation attributes attached to sequence records.
>    > Once such attributes are added, they can be used by the other OBITools commands for filtering purposes or for statistics computing.
>    > 
>    > Here, the goal is to keep only `count` and `merged_sample` key=value attributes! 
>    {: .comment}
>
{: .hands_on}


## Sub-step with **obistat**

> ### {% icon hands_on %} Hands-on: Computes basic statistics for attribute values
>
> 1. {% tool [obistat](toolshed.g2.bx.psu.edu/repos/iuc/obi_stat/obi_stat/1.2.13) %} with the following parameters:
>    - In *"Category attribute"*:
>        - {% icon param-repeat %} *"Insert Category attribute"*
>            - *"How would you specify the category attribute key?"*: `simply by a key of an attribute`
>                - *"Attribute used to categorize the sequence records"*: `count`
>    - *"Use a specific option"*: `no`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > stats computes basic statistics for attribute values of sequence records. The sequence records can be categorized or not using one or several -c options. By default, only the number of sequence records and the total count are computed for each category. Additional statistics can be computed for attribute values in each category, like:
>    > * minimum value (-m option)
>    > * maximum value (-M option)
>    > * mean value (-a option)
>    > * variance (-v option)
>    > * standard deviation (-s option)
>    > 
>    > The result is a contingency table with the different categories in rows, and the computed statistics in columns.
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Can you use this result to say how many sequences occuring only once? You would need to use Galaxy tools like `Sort data in ascending or descending order` and ` Select first lines from a dataset` to answer the question
>
> > ### {% icon solution %} Solution
> >
> > 1. 3131 sequences are occuring once.
> >
> {: .solution}
>
{: .question}

## Sub-step with **obigrep**

> ### {% icon hands_on %} Hands-on: Keep only the sequences having a count greater or equal to 10 and a length shorter than 80 bp
>
> 1. {% tool [obigrep](toolshed.g2.bx.psu.edu/repos/iuc/obi_grep/obi_grep/1.2.13) %} with the following parameters:
>    - *"Choose the sequence record selection option"*: `predicat`
>        - *"Python boolean expression to be evaluated for each sequence record."*: `count>=10`
>
> 2. {% tool [obigrep](toolshed.g2.bx.psu.edu/repos/iuc/obi_grep/obi_grep/1.2.13) %} with the following parameters:
>    - *"Choose the sequence record selection option"*: `lmin`
>        - *"lmin"*: `80`
>
>    > ### {% icon comment %} Comment
>    > Based on the previous observation, we set the cut-off for keeping sequences for further analysis to a count of 10
>    > Based on previous knowledge we also remove sequences with a length shorter than 80 bp (option -l) as we know that the amplified 12S-V5 barcode for vertebrates must have a length around 100bp
>    {: .comment}
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. How many sequences are kept following the "count" filter?
> 2. How many sequences are kept following the "length" filter?
>
> > ### {% icon solution %} Solution
> >
> > 1. 178
> > 2. 175
> >
> {: .solution}
>
{: .question}

## Sub-step with **obiclean**

> ### {% icon hands_on %} Hands-on: Clean the sequences for PCR/sequencing errors (sequence variants)
>
> 1. {% tool [obiclean](toolshed.g2.bx.psu.edu/repos/iuc/obi_clean/obi_clean/1.2.13) %} with the following parameters:
>    - *"Input sequences file"*: `obigrep output file`
>    - *"Maximum numbers of differences between two variant sequences (default: 1)"*: `1`
>    -  *"Threshold ratio between counts (rare/abundant counts) of two sequence records so that the less abundant one is a variant of the more abundant (default: 1, i.e. all less abundant sequences are variants)"*: `0.05`
>    - *"Do you want to select only sequences with the head status in a least one sample?"*: `Yes`
>
>    > ### {% icon comment %} Comment
>    >
>    > As a final denoising step, using the obiclean program, we keep the head sequences that are sequences with no variants with a count greater than 5% of their own count
>    {: .comment}
>
{: .hands_on}


## Sub-step with **NCBI BLAST+ blastn**

> ### {% icon hands_on %} Hands-on: Search nucleotide database with nucleotide query sequence(s) from OBITools treatments
>
> 1. {% tool [NCBI BLAST+ blastn](toolshed.g2.bx.psu.edu/repos/devteam/ncbi_blast_plus/ncbi_blastn_wrapper/2.10.1+galaxy0) %} with the following parameters:
>    - *"Subject database/sequences"*: `FASTA file from your history (see warning note below)`
>    - *"Set expectation value cutoff"*: `0.0001`
>    - *"Output format"*: `Tabular (extended 25 columns)`
>    - *"Advanced Options"*: `Show Advanced Options`
>    - *"Maximum hits to consider/show": `1`
>
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}


## Sub-step with **Filter sequences by ID**

> ### {% icon hands_on %} Hands-on: Filter Blast results
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - *"Sequence file to be filtered"*: `db_v05_r117`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}





## Sub-step with **Filter**

> ### {% icon hands_on %} Hands-on: Filter Blast results
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - *"With following condition"*: `c3>99.99`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
