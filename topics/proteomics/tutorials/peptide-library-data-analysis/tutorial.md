---
layout: tutorial_hands_on

title: Peptide Library Data Analysis
level: Intermediate
questions:
  - How to utilize quantitative properties of amino acids and peptide sequence to analyse peptide data?
objectives:
  - Calculate descriptors
  - Quantitative analysis of peptide sequence properties
time_estimation: '20m'
contributors:
   - jaidevjoshi83
   - blankenberg
subtopic: property-prediction
tags: []
---

# Introduction



 Several computational methods have been proven very useful in the initial screening and prediction of peptides for various biological properties. These methods have emerged as effective alternatives to the lengthy and expensive traditional experimental approaches.  Properties associated with a group of peptide sequences such as overall charge, hydrophobicity profile, or k-mer composition can be utilized to compare peptide sequences and libraries.  In this tutorial, we will be discussing how peptide-based properties like charge, hydrophobicity, the composition of amino acids, etc. can be utilized to analyze the biological properties of peptides. Additionally, we will learn how to use different utilities of the Peptide Design and Analysis Under Galaxy (PDAUG) package to calculate various peptide-based descriptors, and use these descriptors and feature spaces to build informative plots.


### Easy access to tools, workflows and data from the docker image

An easy way to install and use the PDAUG toolset, and follow this tutorial is via a prebuilt docker image equipped with a PDAUG toolset, workflow, and data library. A prebuilds docker image can be downloaded and run by typing a simple command at the terminal after installing docker software on any operating system.

> <hands-on-title>Easy access of tools, workflows and data from docker image</hands-on-title>
>
> 1. Downloading the docker image from the docker hub using `docker pull jayadevjoshi12/galaxy_pdaug:latest` command.
> 2. Running the container with latest PDAUG tools `docker run -i -t -p 8080:80 jayadevjoshi12/galaxy_pdaug:latest`.
> 3. Workflow is available under the workflow section, use `admin` as username and `password` as a password to login as an administrator of your galaxy instance.
> 4. Use `admin` as username and `password` as a password to login galaxy instance, which is available at [localhost](http://127.0.0.1:8080) to access workflow and data.
>
{: .hands_on}


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


## Peptide Data

Several inbuilt data sets have been provided with the tool`PDAUG Peptide Data Access`. The antimicrobial peptides (AMPs) versus transmembrane peptides (TMPs) dataset was used as an example data set to understand the overall relation between features and biological properties of peptides. AMPs consist of an intersection of all activity annotations of the APD2 and CAMP databases, where gram-positive, gram-negative, and antifungal exact matches were observed. TMPs were extracted from alpha-helical transmembrane regions of proteins for classification.

> <hands-on-title>Fetching inbuild data</hands-on-title>
>
> 1. {% tool [PDAUG Peptide Data Access](toolshed.g2.bx.psu.edu/repos/jay/pdaug_peptide_data_access/pdaug_peptide_data_access/0.1.0) %} with the following parameters:
>    - *"Datasets"*: `AMPvsTMP`
>
{: .hands_on}


## Converting tabular data into fasta format

`PDAUG Peptide Data Access` tool returns data as a tabular file that contains sequences from both the classes. In order to utilize this data in the next steps, first we need to convert tabular data into fasta format. If data contains sequences from two different classes `PDAUG TSVtoFASTA` tool converts and splits data into two separate files for each of the class,  AMPs, and TMPs.  The reason behind converting and splitting the data is that all the downstream tools require two separate files if we are comparing two different peptide classes or calculating features.

> <hands-on-title>Converting tabular data into fasta formate</hands-on-title>
>
> 1. {% tool [PDAUG TSVtoFASTA](toolshed.g2.bx.psu.edu/repos/jay/pdaug_tsvtofasta/pdaug_tsvtofasta/0.1.0) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `PDAUG Peptide Data Access - AMPvsTMP (tabular)` (output of **PDAUG Peptide Data Access** {% icon tool %})
>    - *"Peptide Column"*:  `name`
>    - *"Method to convert data"*: `Split Data By Class Label`
>    - *"Column with the class label"*: `class label`
>
>
{: .hands_on}

## Analyzing peptide libraries (AMPs and TMPs) based on features and feature space

### Summary Plot for peptide libraries

In this step, we utilize `PDAUG Peptide Sequence Analysis` tool to compare peptide sequences based on hydrophobicity, hydrophobic movement, charge, amino acid fraction, and sequence length and create a summary plot.



> <hands-on-title>Generating a summary plot to assess peptide dataset</hands-on-title>
>
> 1. {% tool [PDAUG Peptide Sequence Analysis](toolshed.g2.bx.psu.edu/repos/jay/pdaug_peptide_sequence_analysis/pdaug_peptide_sequence_analysis/0.1.0) %} with the following parameters:
>    - *"Analysis options"*: `Plot Summary`
>        - {% icon param-file %} *"First input file"*: `PDAUG TSVtoFASTA on data 1 - first (fasta)` (first output of **PDAUG TSVtoFASTA** {% icon tool %})
>        - {% icon param-file %} *"Second input file"*: `PDAUG TSVtoFASTA on data 1 - second (fasta)` (second output of **PDAUG TSVtoFASTA** {% icon tool %})
>        - *"first input file"*: `TMPs`
>        - *"Second input file"*: `AMPs`
>
>   > <question-title></question-title>
>   > What can be concluded from the summary plot based on different properties?
>   >
>   >  > <solution-title></solution-title>
>   >  > The summary plot represents differences between two sets of peptides based on an amino acid fraction, global charge, sequence length, global hydrophobicity, glocal hydrophobic movement. Additionally, 3D scattered plot shows the clustering of peptides based on three features.
>   >  >
>   >  > 1. Leucine and Valine show relatively higher differences in terms of their fraction within both groups.
>   >  > 2. TMPs show a global charge in the range of 0-5 in comparison to AMPs which show a global charge in a range of 0-14.
>   >  > 3. AMPs show higher variability in terms of their length, global hydrophobic movement, and hydrophobicity in comparison to TMPs.
>   >  > 4. Hydrophobic properties are important in determining transmembrane properties of proteins and peptides which is evident with this summary plot.
>   >  > 5. Clustering of two different kinds of peptides can be observed with a 3D scattered plot based on their properties, however, we can also observe a few peptides with overlapping feature space.
>   >  {: .solution }
>   {: .question}
>
{: .hands_on}

![Summary Plot](../../images/SummaryPlot.png "Summary plot shows comparison between AMPs and TMPs")



### Assessing feature space distribution

In this tool, we have used `PDAUG Fisher's Plot` that compares two peptide libraries based on the feature space using the Fisher test.


> <hands-on-title>Generating a Fisher's plot to assess peptide dataset</hands-on-title>
>
> 1. {% tool [PDAUG Fisher's Plot](toolshed.g2.bx.psu.edu/repos/jay/pdaug_fishers_plot/pdaug_fishers_plot/0.1.0) %} with the following parameters:
>    - {% icon param-file %} *"First fasta file"*: `PDAUG TSVtoFASTA on data 1 - first (fasta)` (first output of **PDAUG TSVtoFASTA** {% icon tool %})
>    - {% icon param-file %} *"Second fasta file"*: `PDAUG TSVtoFASTA on data 1 - second (fasta)` (second output of **PDAUG TSVtoFASTA** {% icon tool %})
>    - *"Label for first population"*: `TMPs`
>    - *"Label for second population"*: `AMPs`
>
>   > <question-title></question-title>
>   > What does Fisher's plot represents?
>   >
>   >  > <solution-title></solution-title>
>   >  >
>   >  > Fisher's plot represents the difference between two groups of peptides based on their feature space. Each tiny square in this plot represents the feature space. Based on the sliding window Fisher's test was performed for each feature space to assess the presence of peptides from two different groups on each of the tiny squares.  The AMPs and TMPs in the feature space represented by their mean hydropathy and amino acid volume. Fisher's plot shows that the sequences with larger hydrophobic amino acids are more frequent in TMPs in comparison to AMPs.
>   >  {: .solution }
>   {: .question}
>
{: .hands_on}

![Fishers plot](../../images/FeatureSpace.png "TMPs peptides show amino acids with larger hydrophobic residues in compare to AMPs")

The AMPs and TMPs in the feature space represented by their mean hydropathy and amino acid volume. Fisher's plot shows that the sequences with larger hydrophobic amino acids are more frequent in TMPs in comparison to AMPs.

## Assessing the relation between peptide features by 3D scatter plot

### Calculating Sequence Property-Based Descriptors

In this step we will calculate Composition, Transition and Distribution (CTD) descriptos. Composition describptors are defined as the number of amino acids of a particular property divided by total number of amino acids.  Transition descriptors are representd as the number of transition from a particular property to different property divided by (total number of amino acids − 1). Distribution descriptors are derived by chain length and the amino acids of a particular property located on this length {% cite Govindan_Nair_2013 %}.


> <hands-on-title>Calculating descriptors for the peptide dataset</hands-on-title>
>
> 1. {% tool [PDAUG Sequence Property Based Descriptors](toolshed.g2.bx.psu.edu/repos/jay/pdaug_sequence_property_based_descriptors/pdaug_sequence_property_based_descriptors/0.1.0) %} with the following parameters:
>    - {% icon param-file %} *"Input fasta file"*: `PDAUG TSVtoFASTA on data 1 - first (fasta)` (first output of **PDAUG TSVtoFASTA** {% icon tool %})
>    - *"Descriptor Type"*: `CTD`
>
> 1. {% tool [PDAUG Sequence Property Based Descriptors](toolshed.g2.bx.psu.edu/repos/jay/pdaug_sequence_property_based_descriptors/pdaug_sequence_property_based_descriptors/0.1.0) %} with the following parameters:
>    - {% icon param-file %} *"Input fasta file"*: `PDAUG TSVtoFASTA on data 1 - second (fasta)` (second output of **PDAUG TSVtoFASTA** {% icon tool %})
>    - *"Descriptor Type"*: `CTD`
>
>
{: .hands_on}


### Adding the Class Label in both AMPs and TMPs

Class labels or target labels usually represents the class of peptides. Here in our data set, we have peptides, either as AMP or TMP. Since we have two classes we can represent these two classes with their actual labels AMPs and TMPs.


- **Adding Class Label (target labels) in AMPs and TMPs data**

> <hands-on-title>Adding Class Labels (target labels) to the tabular data</hands-on-title>
>
> 1. {% tool [PDAUG Add Class Label](toolshed.g2.bx.psu.edu/repos/jay/pdaug_addclasslabel/pdaug_addclasslabel/0.1.0) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `PDAUG Sequence Property Based Descriptors on data 2 - CTD (tabular)` (output of **PDAUG Sequence Property Based Descriptors** {% icon tool %})
>    - *"Class Label"*: `TMPs`
>
> 1. {% tool [PDAUG Add Class Label](toolshed.g2.bx.psu.edu/repos/jay/pdaug_addclasslabel/pdaug_addclasslabel/0.1.0) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `PDAUG Sequence Property Based Descriptors on data 3 - CTD (tabular)` (output of **PDAUG Sequence Property Based Descriptors** {% icon tool %})
>    - *"Class Label"*: `AMPs`
>
>
{: .hands_on}


### Merging the two tabular data files

We utilize `PDAUG Merge Dataframes` to merge two tabular data files.


> <hands-on-title>Merging two tabular data files</hands-on-title>
>
> 1. {% tool [PDAUG Merge Dataframes](toolshed.g2.bx.psu.edu/repos/jay/pdaug_merge_dataframes/pdaug_merge_dataframes/0.1.0) %} with the following parameters:
>    - {% icon param-files %} *"Input files"*: `PDAUG Add Class Label on data 6 - (tabular)` (output of **PDAUG Add Class Label** {% icon tool %}), `PDAUG Add Class Label on data 7 - (tabular)` (output of **PDAUG Add Class Label** {% icon tool %})
>    - *"Option to merg data"*: `Merge data without adding class label`
>
>
{: .hands_on}


### Plotting CTD descriptor data as Scatter plot

Tool `PDAUG Basic Plots` will be used to compare two peptide libraries based on three CTD descriptors `SecondaryStrD1100`, `SolventAccessibilityD2001`, and `NormalizedVDWVD3050` respectively. A 3D scatter plot will be generated.


> <hands-on-title>Generating a scatter plot to assess features</hands-on-title>
>
> 1. {% tool [PDAUG Basic Plots](toolshed.g2.bx.psu.edu/repos/jay/pdaug_basic_plots/pdaug_basic_plots/0.1.0) %} with the following parameters:
>    - *"Data plotting method"*: `Scatter Plot`
>        - {% icon param-file %} *"Input file"*: `PDAUG Merge Dataframes on data 9 and data 8 - (tabular)` (output of **PDAUG Merge Dataframes** {% icon tool %})
>        - *"Scatter Plot type"*: `3D`
>            - *"First feature"*: `_SecondaryStrD1100`
>            - *"Second feature"*: `_SolventAccessibilityD2001`
>            - *"Third feature"*: `_NormalizedVDWVD3050`
>        - *"Class label column"*: `Class_label`
>
>
{: .hands_on}

![3D Scatter plot ](../../images/3DScattered.png "3D scatter Plot shows relation between featues")

**Figure 3** Represent 3D scattered plot generated based on the CTD descriptors. Red dots represent TMPs and blue dots represent AMPs. Based on these 3 features, we can observe that both groups do not show any clear separation or cluster in the 3D space.


In this tutorial, we learned how to utilize inbuild data, calculate features, and utilize descriptors or features to assess biological properties. We also learned how to utilize various utilities of `PDAUG` to generate useful plots to include in our peptide research.

# Conclusion


In this tutorial, we learned an example flexible and extensible analysis of peptide data using PDAUG tools. We generated various plots based on the quantitative properties of amino acids and peptide sequences.

![Workflow](../../images/PeptideAnalysisWorkflow.png "Workflow used")

