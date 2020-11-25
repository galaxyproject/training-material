---
layout: tutorial_hands_on

title: Peptide Library Data Analysis
questions:
  - How to utilize quantitative properties of amino acids and peptide sequence to analyse peptide data?
objectives:
  - Calculate descriptors
  - Quantitative analysis of peptide sequence properties
time_estimation: '20m'
contributors:
   - jaidevjoshi83
   - blankenberg
---

# Introduction
{:.no_toc}

<!-- This is a comment. -->

Several computational methods have been proven very useful in the initial screening and prediction of peptides for various biological properties. These methods have emerged as effective alternatives to the lengthy and expensive traditional experimental approaches. In this tutorial, we will be discussing how to analyze peptide libraries on the basis of quantitative properties. In this tutorial, we will learn how to use different utilities of PDAUG to calculate various peptide-based features and utilize these features for various informative plots.


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


## Peptide Data

In this step, we will retrieve the inbuild dataset, which contains anti-microbial (AMPs) and transmembrane peptides (TMPs).


> ### {% icon hands_on %} Hands-on: Fetching inbuild data
> 
> 1. {% tool [PDAUG Peptide Data Access](https://toolshed.g2.bx.psu.edu/view/jay/pdaug_peptide_data_access/6c12ca9f5d88) %} with the following parameters:
>    - *"Datasets"*: `AMPvsTM` 
>
{: .hands_on}


## Converting tabular data into fasta format

In this step, we will be converting and splitting the tabular data file into fasta format. This tool splits data into two files based on their class labels.

> ### {% icon hands_on %} Hands-on: Converting tabular data into fasta formate
> 
> 1. {% tool [PDAUG TSVtoFASTA](https://toolshed.g2.bx.psu.edu/view/jay/pdaug_tsvtofasta/e272809a193f) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `PDAUG Peptide Data Access - AMPvsTM (tabular)` (output of **PDAUG Peptide Data Access** {% icon tool %})
>    - *"Data conversion"*: ` WithClassLabel `
>
>
{: .hands_on}

## Analyzing peptide libraries (AMPs and TMPs) based on features and feature space 

### Summary Plot for peptide libraries 

In this step, we utilize **PDAUG Peptide Sequence Analysis** tool to compare peptide sequences based on hydrophobicity, hydrophobic movement, charge, amino acid fraction, and sequence length and creates a summary plot.



> ### {% icon hands_on %} Hands-on: Generating a summary plot to assess peptide dataset
> 
> 1. {% tool [PDAUG Peptide Sequence Analysis](https://toolshed.g2.bx.psu.edu/view/jay/padaug_peptide_sequence_analysis/aa0416e63597) %} with the following parameters:
>    - *"Analysis options"*: `Plot Summary`
>        - {% icon param-file %} *"First input file"*: `PDAUG TSVtoFASTA on data 1 - WithClassLabel (tabular)` (first output of **PDAUG TSVtoFASTA** {% icon tool %})
>        - {% icon param-file %} *"Second input file"*: `PDAUG TSVtoFASTA on data 1 - WithClassLabel (tabular)` (second output of **PDAUG TSVtoFASTA** {% icon tool %})
>        - *"Second input file"*: `AMP`
>        - *"Second input file"*: `TM `
>
>
{: .hands_on}

![Summary Plot](../../images/SummaryPlot.png "Summary plot shows comparisoin between AMPs and TMPs")



### Assessing feature space distribution 

In this tool, we have used **PDAUG Fisher's Plot** that compare two peptide library based on the feature space using the Fisher test. 


> ### {% icon hands_on %} Hands-on: Generating a Fisher's plot to assess peptide dataset
> 
> 1. {% tool [PDAUG Fisher's Plot](https://toolshed.g2.bx.psu.edu/view/jay/pdaug_fishers_plot/548a9bddde6c) %} with the following parameters:
>    - {% icon param-file %} *"First fasta file"*: `PDAUG TSVtoFASTA on data 1 - WithClassLabel (tabular)` (first output of **PDAUG TSVtoFASTA** {% icon tool %})
>    - {% icon param-file %} *"Second fasta file"*: `PDAUG TSVtoFASTA on data 1 - WithClassLabel (tabular)` (second output of **PDAUG TSVtoFASTA** {% icon tool %})
>    - *"Label for first population"*: `AMP`
>    - *"Label for second population"*: `TM`
>
>
{: .hands_on}

![Fishers plot](../../images/FeatureSpace.png "TMPs peptides show amino acids with larger hydrophobic residues in compare to AMPs")

The AMPs and TMPs in the feature space represented by their mean hydropathy and amino acid volume. Fisher's plot shows that the sequences with larger hydrophobic amino acids are more frequent in TMPs in comparison to AMPs.

## Assessing the relation between peptide features by 3D scatter plot

### Calculating Sequence Property-Based Descriptors 

In this step, we will be utilizing the "PDAUG Sequence Property Based Descriptors" tool to calculate CTD (Composition Transition and Distribution) descriptor.


> ### {% icon hands_on %} Hands-on: Calculating descriptors for the peptide dataset
> 
> 1. {% tool [PDAUG Sequence Property Based Descriptors](https://toolshed.g2.bx.psu.edu/view/jay/pdaug_sequence_property_based_descriptors/6b6aec52638d) %} with the following parameters:
>    - {% icon param-file %} *"Input fasta file"*: `PDAUG TSVtoFASTA on data 1 - WithClassLabel (tabular)` (first output of **PDAUG TSVtoFASTA** {% icon tool %})
>    - *"Descriptor Type"*: `CTD`
>
> 1. {% tool [PDAUG Sequence Property Based Descriptors](https://toolshed.g2.bx.psu.edu/view/jay/pdaug_sequence_property_based_descriptors/6b6aec52638d) %} with the following parameters:
>    - {% icon param-file %} *"Input fasta file"*: `PDAUG TSVtoFASTA on data 1 - WithClassLabel (tabular)` (second output of **PDAUG TSVtoFASTA** {% icon tool %})
>    - *"Descriptor Type"*: `CTD`
>
>
{: .hands_on}



### Adding Class labels 

- **Adding class label in AMPs data**

> ### {% icon hands_on %} Hands-on: Adding class labels to the data frames
> 
> 1. {% tool [PDAUG Add Class Label](https://toolshed.g2.bx.psu.edu/view/jay/pdaug_addclasslabel/2df11ea23f10) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `PDAUG Sequence Property Based Descriptors on data 2 - CTD (tabular)` (output of **PDAUG Sequence Property Based Descriptors** {% icon tool %})
>    - *"Class Label"*: `1`
>
> 1. {% tool [PDAUG Add Class Label](https://toolshed.g2.bx.psu.edu/view/jay/pdaug_addclasslabel/2df11ea23f10) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `PDAUG Sequence Property Based Descriptors on data 3 - CTD (tabular)` (output of **PDAUG Sequence Property Based Descriptors** {% icon tool %})
>    - *"Class Label"*: `0`
>
>
{: .hands_on}


### Merging the two data frames 

We utilize **PDAUG Merge Dataframes** to merge two dataframes. 


> ### {% icon hands_on %} Hands-on: Merging two dataframes
> 
> 1. {% tool [PDAUG Merge Dataframes](https://toolshed.g2.bx.psu.edu/view/jay/pdaug_merge_dataframes/5bb52d4bf172) %} with the following parameters:
>    - {% icon param-files %} *"Input files"*: `PDAUG Add Class Label on data 6 - (tabular)` (output of **PDAUG Add Class Label** {% icon tool %}), `PDAUG Add Class Label on data 7 - (tabular)` (output of **PDAUG Add Class Label** {% icon tool %})
>
>
{: .hands_on}


### Plotting CTD descriptor data as Scatter plot

In this step, we utilize the **PDAUG Basic Plots** tool to compare two libraries based on three CTD descriptors SecondaryStrD1100, SolventAccessibilityD2001, and NormalizedVDWVD3050 respectively. A 3D scatter plot will be generated. 


> ### {% icon hands_on %} Hands-on: Generating a scatter plot to assess features
> 
> 1. {% tool [PDAUG Basic Plots](https://toolshed.g2.bx.psu.edu/view/jay/pdaug_basic_plots/7d247e27ff11) %} with the following parameters:
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

# Conclusion
{:.no_toc}

In this tutorial, we learned an example flexible and extensible analysis of peptide data using PDAUG tools. We generated various plots based on the quantitative properties of amino acids and peptide sequences.

![Workflow](../../images/PeptideAnalysisWorkflow.png "Workflow used")

