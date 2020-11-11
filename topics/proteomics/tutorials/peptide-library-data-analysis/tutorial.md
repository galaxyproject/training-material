---
layout: tutorial_hands_on

title: Peptide Library Data Analysis
questions:
  - How to utilize quantitative properties of amino acids and peptide sequence to analyse peptide data?
objectives:
  - Calculate descriptors
  - Qunatitative analysis of peptide sequence properties
time_estimation: ''
key_points:
  - The take-home messages
  - They will appear at the end of the tutorial
contributors:
   - jaidevjoshi83
   - blankenberg
---

# Introduction
{:.no_toc}

<!-- This is a comment. -->

Several computational methods have been proven very useful in the initial screening and prediction of peptides for various biological properties. These methods have emerged as effective alternatives to the lengthy and expensive traditional experimental approaches. In this tutorial, we will be discussing how to analyze peptide libraries on the basis of quantitative properties.   In this tutorial, we will learn how to use different utilities of PDAUG to calculate various peptide-based features and utilize these features for various informative plots.  


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


## Peptide Data

In this step, we will retrieve the inbuild dataset, which contains anti-microbial and transmembrane peptides.


> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PDAUG Peptide Data Access** {% icon tool %} with the following parameters:
>    - *"Datasets"*: `AMPvsTM` 
>
{: .hands_on}


## Converting tabular data into fasta formate

In this step, we will be converting and splitting the tabular data file into fasta format.  This tool splits data into two files based on their class labels. 

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PDAUG TSVtoFASTA** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `output1` (output of **PDAUG Peptide Data Access** {% icon tool %})
>    - *"Data conversion"*: ` WithClassLabel `
>
>
{: .hands_on}


## Calculating Sequence Property-Based Descriptors 

In this step, we will be utilizing the "PDAUG Sequence Property Based Descriptors" tool to calculate  CTD (Composition Transition and Distribution) descriptor.  


> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PDAUG Sequence Property Based Descriptors** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input fasta file"*: `OutFile1` (output of **PDAUG TSVtoFASTA** {% icon tool %})
>    - *"Descriptor Type"*: `CTD`
>
>
{: .hands_on}


## Summary Plot for peptide libraries 

In this step, we utilize **PDAUG Peptide Sequence Analysis** tool to compare peptide sequences based on hydrophobicity, hydrophobic movement, charge, amino acid fraction, and sequence length and creates a summary plot.



> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PDAUG Peptide Sequence Analysis** {% icon tool %} with the following parameters:
>    - *"Analysis options"*: `Plot Summary`
>        - {% icon param-file %} *"First input file"*: `OutFile1` (output of **PDAUG TSVtoFASTA** {% icon tool %})
>        - {% icon param-file %} *"Second input file"*: `OutFile2` (output of **PDAUG TSVtoFASTA** {% icon tool %})
>        - *"Second input file"*: `AMP`
>        - *"Second input file"*: `TM `
>
>
{: .hands_on}

![Alternative text](../../images/SummaryPlot.png " Summary plot shows comparisoin between AMPs and TMPs")


## Calculating Sequence Property-Based Descriptors

In this step, we will be utilizing the "PDAUG Sequence Property Based Descriptors" tool to calculate  CTD (Composition Transition and Distribution) descriptor.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PDAUG Sequence Property Based Descriptors** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input fasta file"*: `OutFile2` (output of **PDAUG TSVtoFASTA** {% icon tool %})
>    - *"Descriptor Type"*: `CTD`
>
>
{: .hands_on}


## Assessing feature space distribution 

In this tool, we have used **PDAUG Fisher's Plot** that compare two peptide library based on the feature space using the Fisher test. 


> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PDAUG Fisher's Plot** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"First fasta file"*: `OutFile1` (output of **PDAUG TSVtoFASTA** {% icon tool %})
>    - {% icon param-file %} *"Second fasta file"*: `OutFile2` (output of **PDAUG TSVtoFASTA** {% icon tool %})
>    - *"Lebel for first population"*: `AMP`
>    - *"Label for second population"*: `TM`
>
>
{: .hands_on}

![Alternative text](../../images/FeatureSpace.png " TMPs peptides show amino acides with larger hydrophobic residues in compare to AMPs ")


## Adding Class labels 

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PDAUG Add Class Label** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `output1` (output of **PDAUG Sequence Property Based Descriptors** {% icon tool %})
>
>
{: .hands_on}


> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PDAUG Add Class Label** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `output1` (output of **PDAUG Sequence Property Based Descriptors** {% icon tool %})
>    - *"Class Label"*: `1`
>
>
{: .hands_on}


## Merging the two data frames 

We utilize **PDAUG Merge Dataframes** to merge two dataframes. 


> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PDAUG Merge Dataframes** {% icon tool %} with the following parameters:
>    - {% icon param-files %} *"Input files"*: `OutFile1` (output of **PDAUG Add Class Label** {% icon tool %}), `OutFile1` (output of **PDAUG Add Class Label** {% icon tool %})
>
>
{: .hands_on}


##  Basic data plotting 

In this step, we utilize the **PDAUG Basic Plots** tool to compare two libraries based on three CTD descriptors SecondaryStrD1100, SolventAccessibilityD2001, and NormalizedVDWVD3050 respectively. A 3D scatter plot will be generated. 


> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PDAUG Basic Plots** {% icon tool %} with the following parameters:
>    - *"Data plotting method"*: `Scatter Plot`
>        - {% icon param-file %} *"Input file"*: `output1` (output of **PDAUG Merge Dataframes** {% icon tool %})
>        - *"Scatter Plot type"*: `3D`
>            - *"First feature"*: `_SecondaryStrD1100`
>            - *"Second feature"*: `_SolventAccessibilityD2001`
>            - *"Third feature"*: `_NormalizedVDWVD3050`
>        - *"Class label column"*: `Class_label`
>
>
{: .hands_on}

![Alternative text](../../images/3DScattered.png "3D scatter Plot")

# Conclusion
{:.no_toc}

In this tutorial, we learned the flexible and extensible analysis of the peptide data using PDAUG tools. We generated various plots based on the quantitative properties of amino acids and peptide sequences. 

![Alternative text](../../images/WorkFlow.png "Workflow used")

