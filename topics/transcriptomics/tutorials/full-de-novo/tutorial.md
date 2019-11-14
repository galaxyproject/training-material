---
layout: tutorial_hands_on

title: De novo transcriptome assembly, annotation, and differential expression analysis
zenodo_link: 'https://zenodo.org/record/3541678'
questions:
- Which biological questions are addressed by the tutorial?
- Which bioinformatics techniques are important to know for this type of data?
objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentences describing what a learner should be able to do once they
  have completed the tutorial
- You can use Bloom's Taxonomy to write effective learning objectives
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- abretaud
- r1corre
- lecorguille

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

As a result of the development of novel sequencing technologies, the years between 2008 and 2012 saw a large drop in the cost of sequencing. Per megabase and genome, the cost dropped to 1/100,000th and 1/10,000th of the price, respectively. Prior to this, only transcriptomes of organisms that were of broad interest and utility to scientific research were sequenced; however, these developed in 2010s high-throughput sequencing (also called next-generation sequencing) technologies are both cost- and labor- effective, and the range of organisms studied via these methods is expanding.

Examining non-model organisms can provide novel insights into the mechanisms underlying the "diversity of fascinating morphological innovations" that have enabled the abundance of life on planet Earth. In animals and plants, the "innovations" that cannot be examined in common model organisms include mimicry, mutualism, parasitism, and asexual reproduction. De novo transcriptome assembly is often the preferred method to studying non-model organisms, since it is cheaper and easier than building a genome, and reference-based methods are not possible without an existing genome. The transcriptomes of these organisms can thus reveal novel proteins and their isoforms that are implicated in such unique biological phenomena.

[(source)](https://en.wikipedia.org/wiki/De_novo_transcriptome_assembly)

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Cleaning

Known sequencing biases:
- Unknown nucleotides (Ns)
- Bad quality nucleotides
- Hexamers biases (Illumina. Now corrected ?)

Why do we need to correct those?
- To remove a lot of sequencing errors (detrimental to the vast majority of assemblers)
- Because most de-bruijn graph based assemblers can’t handle unknown nucleotides

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
>
>    {% include snippets/create_new_history.md %}
>
> 2. Import the 12 `fq.gz` into a `List of Pairs` collection named `fastq_raw`
>    - Option 1: from a shared data library (ask your instructor)
>    - Option 2: from Zenodo using the URLs given below
>
>      [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3541678.svg)](https://doi.org/10.5281/zenodo.3541678)
>
>    ```
>    https://zenodo.org/record/3541678/files/A1_left.fq.gz
>    https://zenodo.org/record/3541678/files/A1_right.fq.gz
>    https://zenodo.org/record/3541678/files/A2_left.fq.gz
>    https://zenodo.org/record/3541678/files/A2_right.fq.gz
>    https://zenodo.org/record/3541678/files/A3_left.fq.gz
>    https://zenodo.org/record/3541678/files/A3_right.fq.gz
>    https://zenodo.org/record/3541678/files/B1_left.fq.gz
>    https://zenodo.org/record/3541678/files/B1_right.fq.gz
>    https://zenodo.org/record/3541678/files/B2_left.fq.gz
>    https://zenodo.org/record/3541678/files/B2_right.fq.gz
>    https://zenodo.org/record/3541678/files/B3_left.fq.gz
>    https://zenodo.org/record/3541678/files/B3_right.fq.gz
>    ```
>
>    {% include snippets/import_via_link.md collection=true collection_type="List of Pairs" collection_name="fastq_raw" %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

## Quality control with **FastQC**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **FastQC** {% icon tool %} with the following parameters:
>   - *"Short read data from your current history"*: `fastq_raw` (collection)
>
{: .hands_on}

<!-- ## Quality control with **MultiQC** - step 2/2

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **MultiQC** {% icon tool %} with the following parameters:
>    - In *"Results"*:
>        - {% icon param-repeat %} *"Insert Results"*
>            - *"Which tool was used generate logs?"*: `FastQC`
>                - In *"FastQC output"*:
>                    - *"Type of FastQC output?"*: `Raw data`
>                    - *"FastQC output"*: `data XX, data XX, and others (flattened)`
>
>    > ### {% icon comment %} Comment
>    >
>    > We agree that it's not comfortable. The wrapper of MultiQC must be improved
>    {: .comment}
>
{: .hands_on} -->

## Clean with **Trimmomatic**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Trimmomatic** {% icon tool %} with the following parameters:
>    - *"Single-end or paired-end reads?"*: `Paired-end (as collection)`
>    - *"Select FASTQ dataset collection with R1/R2 pair"*: `fastq_raw`
>    - *"Perform initial ILLUMINACLIP step?"*: `Yes`
>    - *"Adapter sequences to use"*: `TruSeq3 (additional seqs) (paired-ended, for MiSeq and HiSeq)`
>    - In *"Trimmomatic Operation"*:
>        - {% icon param-repeat %} *"Insert Trimmomatic Operation"*
>            - *"Select Trimmomatic operation to perform"*: `Cut bases off end of a read, if below a threshold quality (TRAILING)`
>        - {% icon param-repeat %} *"Insert Trimmomatic Operation"*
>            - *"Select Trimmomatic operation to perform"*: `Cut bases off start of a read, if below a threshold quality (LEADING)`
>        - {% icon param-repeat %} *"Insert Trimmomatic Operation"*
>            - *"Select Trimmomatic operation to perform"*: `Sliding window trimming (SLIDINGWINDOW)`
>        - {% icon param-repeat %} *"Insert Trimmomatic Operation"*
>            - *"Select Trimmomatic operation to perform"*: `Drop reads with average quality lower than a specific level (AVGQUAL)`
>                - *"Minimum length of reads to be kept"*: `25`
>        - {% icon param-repeat %} *"Insert Trimmomatic Operation"*
>            - *"Select Trimmomatic operation to perform"*: `Drop reads below a specified length (MINLEN)`
>                - *"Minimum length of reads to be kept"*: `50`
>    - *"Output trimmomatic log messages?"*: `Yes`
> 2. **Rename** the Dataset Collection
>    - `Trimmomatic on collection XX: paired` -> `fastqc_cleaned`
>
>    > ### {% icon comment %} Comment
>    >
>    > You can check the Trimmomatic log files to get the number of read before and after the cleaning
>    > ```
>    > Input Read Pairs: 10000
>    > Both Surviving: 8804 (88.04%)
>    > Forward Only Surviving: 491 (4.91%)
>    > Reverse Only Surviving: 456 (4.56%) Dropped: 249 (2.49%)
>    > ```
>    {: .comment}
>
>    {% include snippets/rename_collection.md %}
>
{: .hands_on}

## Quality control after cleaning with **FastQC**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **FastQC** {% icon tool %} with the following parameters:
>   - *"Short read data from your current history"*: `fastqc_cleaned` (collection)
>
{: .hands_on}
>
{: .question}

# Assembly

## Assembly with **Trinity**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Trinity** {% icon tool %} with the following parameters:
>    - *"Are you pooling sequence datasets?"*: `Yes`
>        - *"Paired or Single-end data?"*: `Paired-end collection`
>            - *"Strand specific data"*: `No`
>    - *"Run in silico normalization of reads"*: `No`
>    - In *"Additional Options"*:
>        - *"Use the genome guided mode?"*: `No`
> 2. **Rename** the Trinity output
>    - `Trinity on data 52, data 51, and others: Assembled Transcripts` -> `transcriptome_raw.fasta`
>
{: .hands_on}

# Assembly cleanning

## Checking of the assembly with **Trinity Statistics**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Trinity Statistics** {% icon tool %} with the following parameters:
>    - *"Trinity assembly"*: `transcriptome_raw.fasta`
>
>
{: .hands_on}

## Remapping on the raw transcriptome using **Align reads and estimate abundance**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Align reads and estimate abundance** {% icon tool %} with the following parameters:
>    - *"Paired or Single-end data?"*: `Paired`
>        - *"Left/Forward strand reads"* -> `Multiple datasets`
>            - Click on the *Folder* button at the right
>                - *Type to Search*: `left`
>                - Select the 6 `Trimmomatic on ..._left.fq.gz`
>        - *"Right/Reverse strand reads"* -> `Multiple datasets`
>            - Click on the *Folder* button at the right
>                - *Type to Search*: `right`
>                - Select the 6 `Trimmomatic on ..._left.fq.gz`
>        - *"Strand specific data"*: `Yes`
>    - *"Abundance estimation method"*: `Salmon`
>    - In *"Additional Options"*:
>        - *"Trinity assembly?"*: `Yes`
>
>    > ### {% icon comment %} Comment
>    >
>    > If you check at the Standard Error messages of your outputs. You can get the `Mapping rate`
>    > 1. Click on one dataset
>    > 2. Click on the little **i** icon
>    > 3. Click on *Tool Standard Error:	stderr*
>    > ```
>    > [2019-11-14 15:44:21.500] [jointLog] [info] Mapping rate = 44.4358%
>    > ```
>    {: .comment}
>
> ### {% icon comment %} Comment
>    >
>    > At this stage, you can now delete some useless datasets
>    > - `Trimmomatic on collection XX: unpaired`
>    > - `Align reads and estimate abundance on *: genes counts`
>    > Note that the dataset are just hidden. You can delete them permanently and make some room in the history options (the little wheel icon)
>    {: .comment}
>
> 2. **Rename** the 6 `* isoforms counts` :(
>    - Check in the information panel (**i** icon) the lineage of your file (ex: `A1_left.fq.gz` ... )
>    - Rename the datasets: `A1`, `A2`, `A3`, `B1`, `B2`, `B3`.
>
>    {% include snippets/rename_dataset.md %}
>
{: .hands_on}

## Merge the mapping tables and compute a TMM normalization with **Build expression matrix**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Build expression matrix** {% icon tool %} with the following parameters:
>    - *"Abundance estimates"*: `A1`, `A2`, `A3`, `B1`, `B2`, `B3`
>    - *"Abundance estimation method"*: `Salmon`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> What are the three tables?
>
> > ### {% icon solution %} Solution
> >
> > 1. `estimated RNA-Seq fragment isoform counts (raw counts)``
> > 2. `matrix of isoform TPM expression values (not cross-sample normalized)`
> > 3. `matrix of TMM-normalized expression values`
> >
> {: .solution}
>
{: .question}

# IN PROGRESS

## Sub-step with **RNASeq samples quality check**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **RNASeq samples quality check** {% icon tool %} with the following parameters:
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

## Sub-step with **Filter low expression transcripts**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Filter low expression transcripts** {% icon tool %} with the following parameters:
>    - *"Minimum expression level required across any sample"*: `1.0`
>    - *"Isoform filtering method"*: `Keep all isoforms above a minimum percent of dominant expression`
>        - *"Minimum percent of dominant isoform expression"*: `1`
>    - In *"Additional Options"*:
>        - *"Trinity assembly?"*: `Yes`
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

## Sub-step with **Compute contig Ex90N50 statistic and Ex90 transcript count**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Compute contig Ex90N50 statistic and Ex90 transcript count** {% icon tool %} with the following parameters:
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

## Sub-step with **Align reads and estimate abundance**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Align reads and estimate abundance** {% icon tool %} with the following parameters:
>    - *"Paired or Single-end data?"*: `Paired`
>        - *"Strand specific data"*: `Yes`
>    - *"Abundance estimation method"*: `Salmon`
>    - In *"Additional Options"*:
>        - *"Trinity assembly?"*: `Yes`
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

## Sub-step with **TransDecoder**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **TransDecoder** {% icon tool %} with the following parameters:
>    - In *"Training Options"*:
>        - *"Select the training method"*: `Train with the top longest ORFs`
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

## Sub-step with **Build expression matrix**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Build expression matrix** {% icon tool %} with the following parameters:
>    - *"Abundance estimation method"*: `Salmon`
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

## Sub-step with **Trinotate**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Trinotate** {% icon tool %} with the following parameters:
>    - *"Let Galaxy downloading the Trinotate Pre-generated Resource SQLite database"*: `Yes`
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

## Sub-step with **Differential expression analysis**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Differential expression analysis** {% icon tool %} with the following parameters:
>    - *"Differential analysis method"*: `DESeq2`
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

## Sub-step with **Extract and cluster differentially expressed transcripts**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Extract and cluster differentially expressed transcripts** {% icon tool %} with the following parameters:
>    - In *"Additional Options"*:
>        - *"Run GO enrichment analysis"*: `No`
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

## Sub-step with **Partition genes into expression clusters**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Partition genes into expression clusters** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"RData file"*: `rdata` (output of **Extract and cluster differentially expressed transcripts** {% icon tool %})
>    - *"Method for partitioning genes into clusters"*: `Cut tree based on x percent of max(height) of tree`
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
