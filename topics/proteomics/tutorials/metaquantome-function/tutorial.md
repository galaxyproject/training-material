---
layout: tutorial_hands_on

title: metaQuantome-Function
zenodo_link: ""
questions:
 - "How do I look at differentially expressed proteins?"
 - "How do I get top 5 or 10 function present in the sample?"
 - "How can I visualize my results?"
objectives:
  - "A functional and quantitational analysis of metaproteomic mass spectrometry data."
time_estimation: "1h"
key_points:
  - "With function and Quantitation input perform statistical analysis"
contributors:
  - subinamehta
  - timothygriffin
  - pratikdjagtap
  - emmaleith
  - mariecrane
  - pravs3683

follow_up_training:

    -
        type: "internal"
        topic_name: proteomics
        tutorials:
            - metaproteomics
            
        type: "internal"
        topic_name: proteomics
        tutorials:
            - metaquantome-data-creation
            
            


---


# Introduction
{:.no_toc}

In this tutorial we will learn about metaQuantome workflow which leverages peptide level quantitative information to analyze the functional expression within the microbial community in multiple conditions.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# **Pretreatments**

The first step in a tutorial is to get the data from the zenodo link provided and making sure that it is in the correct format.


## *Get data*

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial and give it a meaningful name
>
>    {% include snippets/create_new_history.md %}
>    {% include snippets/rename_history.md %}
>
> 2. Import the files: a Functional File, an Intensity file from [Zenodo]({{ page.zenodo_link }})
>    or from the shared data library (`GTN - Material` -> `{{ page.topic_name }}` -> `{{ page.title }}`)
>
>    ```
>
>    ```
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
>
> 3. Rename the datasets (If needed)
> 4. Check that the datatype ( Make sure they are in the correct formats).
>  Functional File (format=`tabular`), a Intensity file (format=`tabular`)
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to the name of the input data (optional).
>
{: .hands_on}

# **Download metaQuantome Databases**

> ### {% icon hands_on %} Hands-on: Run metaQuantome databases
>
> 1. {% tool [metaQuantome: databases](https://toolshed.g2.bx.psu.edu/view/galaxyp/metaquantome_db/7c45aa1c9464) %} with the following parameters:
>    - {% icon param-select %} *"Databases to Download"*: `Select All` (NCBI Taxonomy Database, Gene Ontology Database and Enzyme Commission (EC) database)
>
>
>    > ### {% icon comment %} Comment
>    >
>    > metaQuantome uses freely available bioinformatic databases to expand your set of direct annotations. 
>    > For most cases, all 3 databases can be downloaded (the default).

The databases are:

1. NCBI taxonomy database. This contains a list of all currently identified taxa and the relationships between them.

2. Gene Ontology (GO) term database. metaQuantome uses the OBO format of the database. Specifically, two files are used: the go-basic.obo file, which is a simplified version of the GO database that is guaranteed to be acyclic, and the metagenomics slim GO, which is a subset of the full GO that is useful for microbiome research. More details are available at http://geneontology.org/docs/download-ontology/

3. ENZYME database with Enzyme Classification (EC) numbers. This database classifies enzymes and organizes the relationships between them.

This module downloads the most recent releases of the specified databases and stores them in a single file, which can then be accessed by the rest of the metaQuantome modules. For reference, the taxonomy database is the largest (~500 Mb), while the GO and EC databases are smaller: ~34 Mb and ~10Mb, respectively.

Also, note that the databases will be stored in the history so that the date of download can be referenced later. Thus, the databases will not be modified, except for the NCBI database, which is updated every time metaQuantome: expand is run (this is a limitation of the Python package ete3, used within metaQuantome).
>    {: .comment}
>
{: .hands_on}


> ### {% icon question %} Questions
>
> 1. What is metaQuanome databases necessary?
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

# **Create metaQuantome sample file**


# **Run metaQuantome**


## **metaQuantome: expand**

> ### {% icon hands_on %} Hands-on: Run metaQuantome expand
>
> 1. {% tool [metaQuantome: expand](toolshed.g2.bx.psu.edu/repos/galaxyp/metaquantome_expand/metaquantome_expand/2.0.0-0) %} with the following parameters:
>    - {% icon param-file %} *"Database Archive File"*: `outfile` (output of **metaQuantome: database** {% icon tool %})
>    - {% icon param-file %} *"Samples file"*: `samples_file` (output of **metaQuantome: create samples file** {% icon tool %})
>    - *"Mode"*: `Functional analysis`
>        - {% icon param-file %} *"Functional file"*: `output` (Input dataset)
>        - *"Ontology"*: `Gene Ontology (GO terms)`
>        - *"Functional column name"*: `go_term`
>    - {% icon param-file %} *"Intensity file"*: `output` (Input dataset)
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

##  **metaQuantome: filter**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [metaQuantome: filter](toolshed.g2.bx.psu.edu/repos/galaxyp/metaquantome_filter/metaquantome_filter/2.0.0-0) %} with the following parameters:
>    - *"Mode"*: `Functional analysis`
>        - *"Ontology"*: `Gene Ontology (GO terms)`
>    - {% icon param-file %} *"Samples file"*: `samples_file` (output of **metaQuantome: create samples file** {% icon tool %})
>    - {% icon param-file %} *"metaquantome expand file"*: `outfile` (output of **metaQuantome: expand** {% icon tool %})
>    - *"min_peptides"*: `1`
>    - *"min_pep_nsamp"*: `1`
>    - *"min_children_non_leaf"*: `2`
>    - *"min_child_nsamp"*: `1`
>    - *"qthreshold"*: `2`
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

##  **metaQuantome: stat**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [metaQuantome: stat](toolshed.g2.bx.psu.edu/repos/galaxyp/metaquantome_stat/metaquantome_stat/2.0.0-0) %} with the following parameters:
>    - {% icon param-file %} *"Samples file"*: `samples_file` (output of **metaQuantome: create samples file** {% icon tool %})
>    - {% icon param-file %} *"metaquantome expanded or filtered file"*: `outfile` (output of **metaQuantome: filter** {% icon tool %})
>    - *"Mode"*: `Functional analysis`
>        - *"Ontology"*: `Gene Ontology (GO terms)`
>    - *"Name of the control sample group of control (used as denominator for fold change evaluation)."*: `T2`
>    - *"Statistical test to perform"*: `standard t-test (parametric)`
>    - *"Perform paired tests"*: `Yes`
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

## **metaQuantome: visualize**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [metaQuantome: visualize](toolshed.g2.bx.psu.edu/repos/galaxyp/metaquantome_viz/metaquantome_viz/2.0.0-0) %} with the following parameters:
>    - {% icon param-file %} *"Tabular file from metaQuantome stats or metaQuantome filter"*: `outfile` (output of **metaQuantome: stat** {% icon tool %})
>    - {% icon param-file %} *"Samples file"*: `outfile` (output of **metaQuantome: stat** {% icon tool %})
>    - *"Mode"*: `Functional analysis`
>        - *"Plot type"*: `Bar Chart (bar)`
>            - *"Mean intensity column name"*: `T2_mean`
>            - *"GO category to restrict to"*: `Molecular Function`
>            - *"Color for the bar fill"*: `Orange`
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

##  **metaQuantome: visualize**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [metaQuantome: visualize](toolshed.g2.bx.psu.edu/repos/galaxyp/metaquantome_viz/metaquantome_viz/2.0.0-0) %} with the following parameters:
>    - {% icon param-file %} *"Tabular file from metaQuantome stats or metaQuantome filter"*: `outfile` (output of **metaQuantome: stat** {% icon tool %})
>    - {% icon param-file %} *"Samples file"*: `samples_file` (output of **metaQuantome: create samples file** {% icon tool %})
>    - *"Mode"*: `Functional analysis`
>        - *"Plot type"*: `Bar Chart (bar)`
>            - *"Mean intensity column name"*: `T4_mean`
>            - *"GO category to restrict to"*: `Molecular Function`
>            - *"Color for the bar fill"*: `Orange`
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

## **metaQuantome: visualize**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [metaQuantome: visualize](toolshed.g2.bx.psu.edu/repos/galaxyp/metaquantome_viz/metaquantome_viz/2.0.0-0) %} with the following parameters:
>    - {% icon param-file %} *"Tabular file from metaQuantome stats or metaQuantome filter"*: `outfile` (output of **metaQuantome: stat** {% icon tool %})
>    - {% icon param-file %} *"Samples file"*: `samples_file` (output of **metaQuantome: create samples file** {% icon tool %})
>    - *"Mode"*: `Functional analysis`
>        - *"Plot type"*: `Bar Chart (bar)`
>            - *"Mean intensity column name"*: `T7_mean`
>            - *"GO category to restrict to"*: `Molecular Function`
>            - *"Color for the bar fill"*: `Orange`
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

## **metaQuantome: visualize**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [metaQuantome: visualize](toolshed.g2.bx.psu.edu/repos/galaxyp/metaquantome_viz/metaquantome_viz/2.0.0-0) %} with the following parameters:
>    - {% icon param-file %} *"Tabular file from metaQuantome stats or metaQuantome filter"*: `outfile` (output of **metaQuantome: stat** {% icon tool %})
>    - {% icon param-file %} *"Samples file"*: `samples_file` (output of **metaQuantome: create samples file** {% icon tool %})
>    - *"Mode"*: `Functional analysis`
>        - *"Plot type"*: `Volcano Plot (volcano)`
>            - *"Name of the fold change column in the stat dataframe"*: `log2fc_T7_over_T2`
>            - *"Name of the Corrected p-value column in the stat dataframe"*: `corrected_p_T7_over_T2`
>            - *"Name of the annotation column"*: `id`
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

## **metaQuantome: visualize**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [metaQuantome: visualize](toolshed.g2.bx.psu.edu/repos/galaxyp/metaquantome_viz/metaquantome_viz/2.0.0-0) %} with the following parameters:
>    - {% icon param-file %} *"Tabular file from metaQuantome stats or metaQuantome filter"*: `outfile` (output of **metaQuantome: stat** {% icon tool %})
>    - {% icon param-file %} *"Samples file"*: `samples_file` (output of **metaQuantome: create samples file** {% icon tool %})
>    - *"Mode"*: `Functional analysis`
>        - *"Plot type"*: `Volcano Plot (volcano)`
>            - *"Name of the fold change column in the stat dataframe"*: `log2fc_T4_over_T2`
>            - *"Name of the Corrected p-value column in the stat dataframe"*: `corrected_p_T4_over_T2`
>            - *"Name of the annotation column"*: `id`
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
