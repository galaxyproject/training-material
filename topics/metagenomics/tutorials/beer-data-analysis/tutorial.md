---
layout: tutorial_hands_on

title: 'Beer data analysis'
zenodo_link: ''
questions:
- What can be observed from the final chart?
objectives:
- To finds yeast strains contained in a sequenced beer sample.
- Learn how to fill parameters for Kraken2, 
- Create visualizations using Krona pie chart

time_estimation: 1H
key_points:
- Inputing correct values for the parameters of proper tools are important

contributors:
- Polina
- Siyu

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

What is a beer microbiome? There are collections of small living
creatures. These small creatures are called bacteria and they are
everywhere. In our gut, in the soil, even on vending machines. Some of these bacteria are actually very good for us. And some others can make us very ill. Bacteria come in different shapes and sizes but they have the same components. One crucial component is the DNA, the blueprint of life. The DNA encodes the shape and size and many other things unique for a bacterial species. Because of the encoding information the DNA can be used to identify what kind of bacteria the DNA is from. Therefore, within a sample form soil, our gut or beer one can specify what kind of species are inside the sample. Follow this tour to learn more about this kind of analysis.

### {% icon comment %} Background

>   ![The beerDeCoded process](../../images/beerprocess.png "The beerDeCoded process contains 3 consistent steps. The first step is DNA extraction from beer. Then, this DNA can be sequenced. That means that we can obtain the sequence of nucleotides for this specific DNA. Finally, we have to analyze received data in order to know which organisms this DNA is from.")
{: .comment}

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Getting started

First of all, this turorial will get you hannds on the basic Galaxy tasks, including creating a history and importing data.

## Create a history

Let's start by creating a new history.

> ### {% icon hands_on %} Hands-on: Create history
>
> 1. Creat an empty analysis history.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. **Rename your history** to be related to the beer data analysis project. 

>    {% snippet faqs/galaxy/histories_rename.md %}
>
{: .hands_on}
 

## Get data

 We will import fastq into the history we just created.

> ### {% icon hands_on %} Hands-on: Upload your dataset
>
> 1. Import the files from [Zenodo]({{ page.zenodo_link }}) or upload your own data
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_upload.md %}
>
> 2. Rename the datasets

{: .hands_on}

# Working with the dataset

Now that your data is ready, let's use some tools

## Assign taxonomic classifications 

One of the key steps in metagenomic data analysis is to identify the taxon to which the individual reads belong. Taxonomic classification tools are based on microbial genome databases to identify the origin of each sequence.

> ### {% icon hands_on %} Hands-on: Kraken2

To perform the taxonomic classification we will use Kraken2. 

> 1.  You can use <b>'tool search'</b> to locate tools. Search for <b>'Kraken2'</b> and select it. Tools may take a couple of moments to load.

> 2. {% tool [Kraken2](toolshed.g2.bx.psu.edu/repos/iuc/kraken2/kraken2/2.0.8_beta+galaxy0) %} with the following parameters:
>    - *"Single or paired reads"*: `Single`
>        - {% icon param-file %} *"Input sequences"*: `output` (Input dataset)
>    - *"Print scientific names instead of just taxids"*: `Yes`
>    - In *"Create Report"*:
>        - *"Print a report with aggregrate counts/clade to file"*: `Yes`
>        - *"Format report output like Kraken 1's kraken-mpa-report"*: `Yes`
> Leave the rest parameters as defaults.
>    > ### {% icon comment %} About Kraken2
>    >
![Kraken2](../../images/kraken2.png "This tool uses the minimizer method to sample the k-mers (all the read’s subsequences of length k) in a deterministic fashion in order to reduce memory constumption and processing time. In addition, it masks low-complexity sequences from reference sequences by using dustmasker.")
{: .comment}
>
{: .hands_on}



## Sub-step with **Reverse**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Reverse](toolshed.g2.bx.psu.edu/repos/iuc/datamash_reverse/datamash_reverse/1.1.0) %} with the following parameters:
>    - {% icon param-file %} *"Input tabular dataset"*: `report_output` (output of **Kraken2** {% icon tool %})
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

## Sub-step with **Replace Text**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_line/1.1.2) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `out_file` (output of **Reverse** {% icon tool %})
>    - In *"Replacement"*:
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `\|`
>            - *"Replace with:"*: `\t`
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

## Sub-step with **Krona pie chart**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Krona pie chart](toolshed.g2.bx.psu.edu/repos/crs4/taxonomy_krona_chart/taxonomy_krona_chart/2.6.1.1) %} with the following parameters:
>    - *"What is the type of your input data"*: `Tabular`
>        - {% icon param-file %} *"Input file"*: `outfile` (output of **Replace Text** {% icon tool %})
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