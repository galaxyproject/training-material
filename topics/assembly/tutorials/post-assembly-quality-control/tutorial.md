---
layout: tutorial_hands_on

title: Post Assembly Quality Control
zenodo_link: ''
questions:
- what combination of tools can control the quality of an initial assembly?
- how to evaluate the quality and the completeness of the assemblies?
objectives:
- apply the post-assembly-QC-workflow using the necessary tools
- evaluate the quality of the post-assembly
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- contributor1
- contributor2

---


# Introduction

<!-- This is a comment. -->

An important part in genome assembly is quality control. During the whole process of DNA sequencing and assembling the genome, errors like mismatches and gaps can occur. This can affect the accuracy and completeness of the genome. Quality control helps to identify and correct these errors by evaluating the quality of the sequences and by detecting and removing potential contaminations.
A high-quality end-product as a result also ensures to avoid errors in downstream analyses.

Since there are many different ways how errors can occur there are also many different tools to identify and remove potential problems. The difficult part is to choose between them and to know when it's time to move on. It is important because time and resources play a big role in genome assembly.
In this tutorial you will learn how to use the tools for the post-assembly quality control
workflow. It's part of a post-assembly pipeline from ERGA to ensure high quality assemblies in appropriate time and resources.





General introduction about the topic and then an introduction of the
tutorial (the questions and the objectives). It is nice also to have a
scheme to sum up the pipeline used during the tutorial. The idea is to
give to trainees insight into the content of the tutorial and the (theoretical
and technical) key concepts they will learn.

You may want to cite some publications; this can be done by adding citations to the
bibliography file (`tutorial.bib` file next to your `tutorial.md` file). These citations
must be in bibtex format. If you have the DOI for the paper you wish to cite, you can
get the corresponding bibtex entry using [doi2bib.org](https://doi2bib.org).

With the example you will find in the `tutorial.bib` file, you can add a citation to
this article here in your tutorial like this:
{% raw %} `{% cite Batut2018 %}`{% endraw %}.
This will be rendered like this: {% cite Batut2018 %}, and links to a
[bibliography section](#bibliography) which will automatically be created at the end of the
tutorial.


**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Title for your first section

Give some background about what the trainees will be doing in the section.
Remember that many people reading your materials will likely be novices,
so make sure to explain all the relevant concepts.

## Title for a subsection
Section and subsection titles will be displayed in the tutorial index on the left side of
the page, so try to make them informative and concise!

# Hands-on Sections
Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

# Get data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

# Providing analysis information/statistics using k-mers:

It is common to analyse assemblies with the help of k-mer counting. During an assembly, the DNA fragments are broken down into k-mers. Then they are compared to identify regions of overlap. By aligning overlapping k-mers it's possible to piece the original DNA sequence together and generate a complete genome.
K-mers are also useful in genome analysis. The frequency and distribution of k-mers can be used to estimate the genome size and identify repetitive sequences.


## Sub-step with **Meryl**

Meryl is a k-mer counter. It is a powerful tool for counting k-mers in large-scale genomic datasets. Meryl uses a sorting-based approach that sorts the k-mers in lexicographical order.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Meryl](toolshed.g2.bx.psu.edu/repos/iuc/meryl/meryl/1.3+galaxy6) %} with the following parameters:
>    - *"Operation type selector"*: `Count operations`
>        - {% icon param-file %} *"Input sequences"*: `output` (Input dataset)
>        - *"K-mer size selector"*: `Estimate the best k-mer size`
>            - *"Genome size"*: `{'id': 4, 'output_name': 'output'}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Meryl**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Meryl](toolshed.g2.bx.psu.edu/repos/iuc/meryl/meryl/1.3+galaxy6) %} with the following parameters:
>    - *"Operation type selector"*: `Generate histogram dataset`
>        - {% icon param-file %} *"Input meryldb"*: `read_db` (output of **Meryl** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Merqury**

Merqury is designed for evaluating the completeness and accuracy of genome assemblies using short-read sequencing data. The quality of assemblies which are generated by using third-generation sequencing technologies can be reviewed and assessed by the tool. Merqury works by comparing the assembly to a high-quality reference genome. It uses k-mer-based methods to estimate and evaluate the completeness of the assembly and identify errors.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Merqury](toolshed.g2.bx.psu.edu/repos/iuc/merqury/merqury/1.3+galaxy2) %} with the following parameters:
>    - *"Evaluation mode"*: `Default mode`
>        - {% icon param-file %} *"K-mer counts database"*: `read_db` (output of **Meryl** {% icon tool %})
>        - *"Number of assemblies"*: `One assembly (pseudo-haplotype or mixed-haplotype)`
>            - {% icon param-file %} *"Genome assembly"*: `output` (Input dataset)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **GenomeScope**

Genomescope is used for analysing genomes using short-read sequencing data. It's a tool to analyse assemblies by estimating genome heterozygosity, repeat content, and size from sequencing reads using a kmer-based statistical approach.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [GenomeScope](toolshed.g2.bx.psu.edu/repos/iuc/genomescope/genomescope/2.0+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Input histogram file"*: `read_db_hist` (output of **Meryl** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

# Assembly decontamination

Extracted DNA from an organism always contains DNA from other species. This is why most assemblies need to go through a decontamination process to remove the non-target DNA for a higher-quality end product.
Contamination can also have  a significant impact on downstream analysis and interpretation of the genome assembly. For example, host contamination can lead to misinterpretation of gene content and functional annotation. Another point is that environmental contamination can affect the accuracy of taxonomic classification. Therefore, it is important to identify and minimise contamination in genome assemblies using appropriate quality control measures.


It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> <details-title> More details about the theory </details-title>
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


## Sub-step with **BlobToolKit**

Blobtoolkit is a tool for decontamination, analysing and visualising assemblies. It's a great tool to review the quality of an assembled genome. In our case optimal suited for the post-assembly quality control.
Blobtoolkit can help to identify contaminants by creating a dataset with taxonomic information and then comparing the assembly with the provided data from known contigs or scaffolds.


> <hands-on-title> Creating the BlobDir dataset </hands-on-title>
>
> 1. {% tool [BlobToolKit](toolshed.g2.bx.psu.edu/repos/bgruening/blobtoolkit/blobtoolkit/3.4.0+galaxy0) %} with the following parameters:
>    - *"Select mode"*: `Create a BlobToolKit dataset`
>        - {% icon param-file %} *"Genome assembly file"*: `output` (Input dataset)
>        - {% icon param-file %} *"Metadata file"*: `output` (Input dataset)
>        - *"NCBI taxonomy ID"*: `{'id': 2, 'output_name': 'output'}`
>        - {% icon param-file %} *"NCBI taxdump directory"*: `output` (Input dataset)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > To work with Blobtoolkit we need to create a new dataset structure called BlobDir. Therefore the minimum requirement is a fasta file which contains the sequence of our assembly. A list of sequence identifiers and some statistics like length, GC proportion and undefined bases will then be generated.
To get a more meaningful analysis and therefore more useful information about our assembly in order to get rid of contaminations, it is better to provide as much data as we can get. In our case we will also provide a Metadata file, NCBI taxonomy ID and NCBI taxdump directory.
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}


## Sub-step with **HISAT2**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [HISAT2](toolshed.g2.bx.psu.edu/repos/iuc/hisat2/hisat2/2.2.1+galaxy1) %} with the following parameters:
>    - *"Source for the reference genome"*: `Use a genome from history`
>        - {% icon param-file %} *"Select the reference genome"*: `output` (Input dataset)
>    - *"Is this a single or paired library"*: `Single-end`
>        - {% icon param-collection %} *"FASTA/Q file"*: `output` (Input dataset collection)
>    - In *"Advanced Options"*:
>        - *"Input options"*: `Use default values`
>        - *"Alignment options"*: `Use default values`
>        - *"Scoring options"*: `Use default values`
>        - *"Spliced alignment options"*: `Use default values`
>        - *"Reporting options"*: `Use default values`
>        - *"Output options"*: `Use default values`
>        - *"SAM options"*: `Use default values`
>        - *"Other options"*: `Use default values`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Busco**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.4.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `output` (Input dataset)
>    - *"Mode"*: `Genome assemblies (DNA)`
>        - *"Use Augustus instead of Metaeuk"*: `Use Metaeuk`
>    - *"Auto-detect or select lineage?"*: `Auto-detect`
>    - *"Which outputs should be generated"*: ``
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **BlobToolKit**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [BlobToolKit](toolshed.g2.bx.psu.edu/repos/bgruening/blobtoolkit/blobtoolkit/3.4.0+galaxy0) %} with the following parameters:
>    - *"Select mode"*: `Add data to a BlobToolKit dataset`
>        - {% icon param-file %} *"Blobdir.tgz file"*: `blobdir` (output of **BlobToolKit** {% icon tool %})
>        - {% icon param-file %} *"BUSCO full table file"*: `busco_table` (output of **Busco** {% icon tool %})
>        - *"BLAST/Diamond hits"*: `Disabled`
>        - {% icon param-file %} *"BAM/SAM/CRAM read alignment file"*: `output_alignments` (output of **HISAT2** {% icon tool %})
>        - *"Genetic text file"*: `Disabled`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **BlobToolKit**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [BlobToolKit](toolshed.g2.bx.psu.edu/repos/bgruening/blobtoolkit/blobtoolkit/3.4.0+galaxy0) %} with the following parameters:
>    - *"Select mode"*: `Generate plots`
>        - {% icon param-file %} *"Blobdir file"*: `blobdir` (output of **BlobToolKit** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}



## Sub-step with **gfastats**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [gfastats](toolshed.g2.bx.psu.edu/repos/bgruening/gfastats/gfastats/1.2.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `output` (Input dataset)
>    - *"Specify target sequences"*: `Disabled`
>    - *"Tool mode"*: `Summary statistics generation`
>        - *"Report mode"*: `Genome assembly statistics (--nstar-report)`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
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

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.