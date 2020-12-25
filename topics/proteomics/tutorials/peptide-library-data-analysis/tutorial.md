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
---

# Introduction
{:.no_toc}


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
> 1. {% tool [PDAUG Peptide Data Access](toolshed.g2.bx.psu.edu/repos/jay/pdaug_peptide_data_access/pdaug_peptide_data_access/0.1.0) %} with the following parameters:
>    - *"Datasets"*: `AMPvsTM` 
>
{: .hands_on}


## Converting tabular data into fasta format

In this step, we will be converting and splitting the tabular data file into fasta format. This tool splits data into two files based on their class labels.

> ### {% icon hands_on %} Hands-on: Converting tabular data into fasta formate
> 
> 1. {% tool [PDAUG TSVtoFASTA](toolshed.g2.bx.psu.edu/repos/jay/pdaug_tsvtofasta/pdaug_tsvtofasta/0.1.0) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `PDAUG Peptide Data Access - AMPvsTM (tabular)` (output of **PDAUG Peptide Data Access** {% icon tool %})
>    - *"Data conversion"*: ` WithClassLabel `
>
>
{: .hands_on}

## Analyzing peptide libraries (AMPs and TMPs) based on features and feature space 

### Summary Plot for peptide libraries 

In this step, we utilize **PDAUG Peptide Sequence Analysis** tool to compare peptide sequences based on hydrophobicity, hydrophobic movement, charge, amino acid fraction, and sequence length and create a summary plot.



> ### {% icon hands_on %} Hands-on: Generating a summary plot to assess peptide dataset
> 
> 1. {% tool [PDAUG Peptide Sequence Analysis](toolshed.g2.bx.psu.edu/repos/jay/pdaug_peptide_sequence_analysis/pdaug_peptide_sequence_analysis/0.1.0) %} with the following parameters:
>    - *"Analysis options"*: `Plot Summary`
>        - {% icon param-file %} *"First input file"*: `PDAUG TSVtoFASTA on data 1 - first (fasta)` (first output of **PDAUG TSVtoFASTA** {% icon tool %})
>        - {% icon param-file %} *"Second input file"*: `PDAUG TSVtoFASTA on data 1 - second (fasta)` (second output of **PDAUG TSVtoFASTA** {% icon tool %})
>        - *"first input file"*: `AMP`
>        - *"Second input file"*: `TM `
>
>   > ### {% icon question %} Questions
>   > What can be concluded from the summary plot based on different properties?
>   >
>   >  > ### {% icon solution %} Solution
>   >  > The summary plot represents differences between two sets of peptides based on an amino acid fraction, global charge, sequence length, global hydrophobicity, glocal hydrophobic movement. Additionally, 3D scattered plot shows the clustering of peptides based on three features.
>   >  >  
>   >  > 1. Leucine and Valine show relatively higher differences in terms of their fraction within both groups. 
>   >  > 2. TMPs show a global charge in the range of 0-5 in comparison to AMPs which show a global charge in a range of 0-14.
>   >  > 3. AMPs show higher variability in terms of their length, global hydrophobic movement, and hydrophobicity in comparison to TMPs.
>   >  > 4. Hydrophobicity of peptides is important to determine their transmembrane properties which is evident with this summary plot.
>   >  > 5. Clustering of two different kinds of peptides can be observed with a 3D scattered plot based on their properties, however, we can also observe a few peptides with overlapping feature space.
>   >  {: .solution }
>   {: .question}
>
{: .hands_on}

![Summary Plot](../../images/SummaryPlot.png "Summary plot shows comparisoin between AMPs and TMPs")



### Assessing feature space distribution 

In this tool, we have used **PDAUG Fisher's Plot** that compare two peptide libraries based on the feature space using the Fisher test. 


> ### {% icon hands_on %} Hands-on: Generating a Fisher's plot to assess peptide dataset
> 
> 1. {% tool [PDAUG Fisher's Plot](toolshed.g2.bx.psu.edu/repos/jay/pdaug_fishers_plot/pdaug_fishers_plot/0.1.0) %} with the following parameters:
>    - {% icon param-file %} *"First fasta file"*: `PDAUG TSVtoFASTA on data 1 - first (fasta)` (first output of **PDAUG TSVtoFASTA** {% icon tool %})
>    - {% icon param-file %} *"Second fasta file"*: `PDAUG TSVtoFASTA on data 1 - second (fasta)` (second output of **PDAUG TSVtoFASTA** {% icon tool %})
>    - *"Label for first population"*: `AMP`
>    - *"Label for second population"*: `TM`
>
>   > ### {% icon question %} Questions
>   > What does Fisher's plot represents?
>   >
>   >  > ### {% icon solution %} Solution
>   >  > The summary plot represents differences between two sets of peptides based on an amino acid fraction, global charge, sequence length, global hydrophobicity, glocal hydrophobic movement. Additionally, 3D scattered plot shows the clustering of peptides based on three features.
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

In this step we will calculate CTD descriptos. Composition describptors are defined as the number of amino acids of a particular property divided by total number of amino acids.  Transition descriptors are representd as the number of transition from a particular property to different property divided by (total number of amino acids − 1). Distribution descriptors are derived by chain length and the amino acids of a particular property are located on this length{% cite Govindan_Nair_2013 %}.


> ### {% icon hands_on %} Hands-on: Calculating descriptors for the peptide dataset
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


### Adding Class labels in both AMPs and TMPs

Usually, class labels are represented by 0 or 1. If data has multi-class classification problems it can be represented by 0,1,2,3, etc. In addition to this, class labels can also be represented by a specific string such as "anticancer" and "non-anticancer" or "treated" and "untreated".

- **Adding class label in AMPs and TMPs data**

> ### {% icon hands_on %} Hands-on: Adding class labels to the tabular data
> 
> 1. {% tool [PDAUG Add Class Label](toolshed.g2.bx.psu.edu/repos/jay/pdaug_addclasslabel/pdaug_addclasslabel/0.1.0) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `PDAUG Sequence Property Based Descriptors on data 2 - CTD (tabular)` (output of **PDAUG Sequence Property Based Descriptors** {% icon tool %})
>    - *"Class Label"*: `1`
>
> 1. {% tool [PDAUG Add Class Label](toolshed.g2.bx.psu.edu/repos/jay/pdaug_addclasslabel/pdaug_addclasslabel/0.1.0) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `PDAUG Sequence Property Based Descriptors on data 3 - CTD (tabular)` (output of **PDAUG Sequence Property Based Descriptors** {% icon tool %})
>    - *"Class Label"*: `0`
>
>
{: .hands_on}


### Merging the two tabular data files 

We utilize **PDAUG Merge Dataframes** to merge two tabular data files. 


> ### {% icon hands_on %} Hands-on: Merging two tabular data files
> 
> 1. {% tool [PDAUG Merge Dataframes](toolshed.g2.bx.psu.edu/repos/jay/pdaug_merge_dataframes/pdaug_merge_dataframes/0.1.0) %} with the following parameters:
>    - {% icon param-files %} *"Input files"*: `PDAUG Add Class Label on data 6 - (tabular)` (output of **PDAUG Add Class Label** {% icon tool %}), `PDAUG Add Class Label on data 7 - (tabular)` (output of **PDAUG Add Class Label** {% icon tool %})
>
>
{: .hands_on}


### Plotting CTD descriptor data as Scatter plot

In this step, we utilize the **PDAUG Basic Plots** tool to compare two libraries based on three CTD descriptors SecondaryStrD1100, SolventAccessibilityD2001, and NormalizedVDWVD3050 respectively. A 3D scatter plot will be generated. 


> ### {% icon hands_on %} Hands-on: Generating a scatter plot to assess features
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

# Conclusion
{:.no_toc}

In this tutorial, we learned an example flexible and extensible analysis of peptide data using PDAUG tools. We generated various plots based on the quantitative properties of amino acids and peptide sequences.

![Workflow](../../images/PeptideAnalysisWorkflow.png "Workflow used")

