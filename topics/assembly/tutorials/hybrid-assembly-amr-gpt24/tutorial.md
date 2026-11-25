---
layout: tutorial_hands_on

title: Bacterial hybrid genome assembly for AMR surveillance
zenodo_link: ''
questions:
- How to run Unicycler with Illumina and Oxford Nanopore sequecing data?
- Which points should we pay attention when running hybrid assembly in an ARM surveillance research?
objectives: 
- Familiarize yourself with bioinformatic analysis process from sequencing raw read to antimicrobial resistance profile
- Improve the quality of sequencing data
- Understand the outputs of tools to assess the sequencing read quality
- Learn how to run Unicycler                           
- Understand the outputs of assembly tool
- Learn how to run AMR gene detection tools
time_estimation: 6H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- trongthucnguyen
- praissyzefi
- contributor3

---



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

<!-- This is a comment. -->

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

## Get data

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

# Title of the section usually corresponding to a big step in the analysis

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


## Sub-step with **seqtk_sample**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [seqtk_sample](toolshed.g2.bx.psu.edu/repos/iuc/seqtk/seqtk_sample/1.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input FASTA/Q file"*: `output` (Input dataset)
>    - *"RNG seed"*: `42`
>    - *"Subsample (decimal fraction or number)"*: `0.2`
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

## Sub-step with **FastQC**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [FastQC](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.74+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Raw read data from your current history"*: `output` (Input dataset)
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

## Sub-step with **fastp**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [fastp](toolshed.g2.bx.psu.edu/repos/iuc/fastp/fastp/1.4.0+galaxy0) %} with the following parameters:
>    - *"Single-end or paired reads"*: `Paired Collection`
>        - {% icon param-collection %} *"Select paired collection(s)"*: `output` (Input dataset collection)
>    - In *"Duplicated Reads Options"*:
>        - *"Enable duplicated reads analysis"*: `Disable (--dont_eval_duplication)`
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

## Sub-step with **FastQC**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [FastQC](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.74+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Raw read data from your current history"*: `output` (Input dataset)
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

## Sub-step with **NanoPlot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NanoPlot](toolshed.g2.bx.psu.edu/repos/iuc/nanoplot/nanoplot/1.48.0+galaxy0) %} with the following parameters:
>    - *"Select multifile mode"*: `batch`
>        - *"Type of the file(s) to work on"*: `fastq`
>            - {% icon param-file %} *""*: `default` (output of **seqtk_sample** {% icon tool %})
>    - In *"Options for filtering or transforming input prior to plotting"*:
>        - *"Logarithmic scaling of lengths in plots."*: `Yes`
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

## Sub-step with **Fastplong**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Fastplong](toolshed.g2.bx.psu.edu/repos/iuc/fastplong/fastplong/0.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input Long Reads"*: `default` (output of **seqtk_sample** {% icon tool %})
>    - In *"Adapter Trimming Options"*:
>        - *"Disable adapter trimming"*: `NO`
>            - *"Chimeric read handling"*: `Split chimeric reads into segments (default)`
>    - In *"Filter Options"*:
>        - In *"Quality filtering options"*:
>            - *"Disable quality filtering"*: `NO`
>        - In *"Length filtering options"*:
>            - *"Disable length filtering"*: `NO`
>        - In *"Low complexity filtering options"*:
>            - *"Enable low complexity filter"*: `NO`
>    - In *"Read Modification Options"*:
>        - In *"Polyx trimming on tail"*:
>            - *"PolyX tail trimming"*: `NO`
>        - In *"Per read cutting by quality options"*:
>            - *"Activate trimming by quality"*: `No quality trimming`
>        - In *"Low quality region handling"*:
>            - *"Low quality region handling"*: `None`
>    - In *"Read sampling and selection"*:
>        - *"Sampling / selection mode"*: `No sampling - keep all reads passing the filters`
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

## Sub-step with **filtlong**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [filtlong](toolshed.g2.bx.psu.edu/repos/iuc/filtlong/filtlong/0.3.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input FASTQ"*: `out` (output of **Fastplong** {% icon tool %})
>    - In *"Output thresholds"*:
>        - *"Min. length"*: `2000`
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

## Sub-step with **NanoPlot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NanoPlot](toolshed.g2.bx.psu.edu/repos/iuc/nanoplot/nanoplot/1.48.0+galaxy0) %} with the following parameters:
>    - *"Select multifile mode"*: `batch`
>        - *"Type of the file(s) to work on"*: `fastq`
>            - {% icon param-file %} *""*: `outfile` (output of **filtlong** {% icon tool %})
>    - In *"Options for filtering or transforming input prior to plotting"*:
>        - *"Logarithmic scaling of lengths in plots."*: `Yes`
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

## Sub-step with **Create assemblies with Unicycler**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Create assemblies with Unicycler](toolshed.g2.bx.psu.edu/repos/iuc/unicycler/unicycler/0.5.1+galaxy1) %} with the following parameters:
>    - *"Paired or Single end data?"*: `Paired Collection`
>        - {% icon param-file %} *"Select a paired collection"*: `output_paired_coll` (output of **fastp** {% icon tool %})
>    - {% icon param-file %} *"Select long reads. If there are no long reads, leave this empty"*: `outfile` (output of **filtlong** {% icon tool %})
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

## Sub-step with **Bandage Image**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Bandage Image](toolshed.g2.bx.psu.edu/repos/iuc/bandage/bandage_image/2022.09+galaxy4) %} with the following parameters:
>    - {% icon param-file %} *"Graphical Fragment Assembly"*: `assembly_graph` (output of **Create assemblies with Unicycler** {% icon tool %})
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

## Sub-step with **Quast**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Quast](toolshed.g2.bx.psu.edu/repos/iuc/quast/quast/5.3.0+galaxy3) %} with the following parameters:
>    - *"Assembly mode?"*: `Individual assembly (1 contig file per sample)`
>        - *"Use customized names for the input files?"*: `No, use dataset names`
>            - {% icon param-file %} *"Contigs/scaffolds file"*: `assembly` (output of **Create assemblies with Unicycler** {% icon tool %})
>        - *"Reads options"*: `Disabled`
>    - *"Type of assembly"*: `Genome`
>        - *"Use a reference genome?"*: `No`
>    - In *"Genes"*:
>        - *"Tool for gene prediction"*: `Don't predict genes`
>    - *"Output files"*: ``
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

## Sub-step with **staramr**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [staramr](toolshed.g2.bx.psu.edu/repos/iuc/staramr/staramr_search/0.12.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *""*: `assembly` (output of **Create assemblies with Unicycler** {% icon tool %})
>    - *"Enable scanning for point mutations using the PointFinder database"*: `Escherichia coli (escherichia_coli)`
>    - In *"Advanced options"*:
>        - *"Exclude certain AMR genes from the results"*: `Exclude default list of AMR genes`
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

## Sub-step with **checkm2**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [checkm2](toolshed.g2.bx.psu.edu/repos/iuc/checkm2/checkm2/1.1.0+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input MAG/SAG datasets"*: `assembly` (output of **Create assemblies with Unicycler** {% icon tool %})
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
