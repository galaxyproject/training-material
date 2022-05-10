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
![Kraken2](../../images/kraken2.jpg "This tool uses the minimizer method to sample the k-mers (all the read’s subsequences of length k) in a deterministic fashion in order to reduce memory constumption and processing time. In addition, it masks low-complexity sequences from reference sequences by using dustmasker.")
{: .comment}

 Congratulations! you have created two files. It contains Classification and Report file. It will remain stored in your history. Click the <b>'eye'</b> icon to view the data once the history item turns green.
>
{: .hands_on}

# Analyze taxonomic assigment

 Once we have assigned the corresponding taxa to the sequence, the next step is to properly visualize the data, for which we will use the Krona pie chart tool (<a href="https://doi.org/10.1186/1471-2105-12-385">Ondov et al. 2011</a>).

##Adjust dataset format

But before that, we need to adjust the format of the data output from Kraken2.

Search and select the <b>'Reverse'</b> tool.

> ### {% icon hands_on %} Hands-on: Reverse
>
> 1. {% tool [Reverse](toolshed.g2.bx.psu.edu/repos/iuc/datamash_reverse/datamash_reverse/1.1.0) %} with the following parameters:
>    - {% icon param-file %} *"Input tabular dataset"*: `report_output` (output of **Kraken2** {% icon tool %})
>
{: .hands_on}

Then search and select the <b>'Replace'</b> tool.

> ### {% icon hands_on %} Hands-on: Replace Text
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_line/1.1.2) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `out_file` (output of **Reverse** {% icon tool %})
>    - In *"Replacement"*:
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `\|`
>            - *"Replace with:"*: `\t`
>
{: .hands_on}


## Visualize the taxonomical classification  

<b>Krona</b> allows hierarchical data to be explored with zooming, multi-layered pie charts. With this tool, we can easily visualize the composition of the bacterial communities and compare how the populations of microorganisms are modified according to the conditions of the environment.

> ### {% icon hands_on %} Hands-on: Krona pie chart
>
> 1. {% tool [Krona pie chart](toolshed.g2.bx.psu.edu/repos/crs4/taxonomy_krona_chart/taxonomy_krona_chart/2.6.1.1) %} with the following parameters:
>    - *"What is the type of your input data"*: `Tabular`
>        - {% icon param-file %} *"Input file"*: `outfile` (output of **Replace Text** {% icon tool %})
>
  2. Let’s take a look at the result by clicking eye icon. Using the search bar we can check if certain taxa are present.


# Conclusion
{:.no_toc}
You have reached the end of the tour. Thank you for going through our tutorial.