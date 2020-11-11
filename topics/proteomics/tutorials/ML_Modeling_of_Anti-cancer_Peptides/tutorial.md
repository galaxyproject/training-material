---
layout: tutorial_hands_on

title: ML Modeling of Anti-cancer Peptides
questions:
- Which ML algorithm is superior in predicting anti-cancer peptides?
objectives:
- Learn, how to calculate peptide descriptor 
- Learn, how to create training data set from features?
- Assessment of best ML algorithm in predicting anti-cancer peptide 
time_estimation: '1h'

contributors:
- contributor1
- contributor2

---

## Introduction
{:.no_toc}

<!-- This is a comment. -->

Several computational methods have been proven very useful in the initial screening and prediction of peptides for various biological properties. These methods have emerged as effective alternatives to the lengthy and expensive traditional experimental approaches. In this tutorial, we will be discussing how peptide-based properties like charge, hydrophobicity, the composition of amino acids, etc. can be utilized to predict the biological properties of peptides. In this tutorial, we will learn how to use different utilities of PDAUG to calculate various peptide-based descriptors and use these descriptors for machine learning modeling of peptides with known anti-cancer properties. We will use CTD (composition, transition, and distribution ) descriptor to define peptide sequences in the training set and will test 6 different machine learning algorithms.  We will also assess the effect of normalization on the accuracy of ML models.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Training data set  

A high-quality dataset was retrieved from a previously published work (Hajisharifi et al., 2014).  A balanced dataset, 138 ACPs (positive), and 138 non-ACPs (negative) sequences were obtained by randomly removing negative sequences. The length distribution of the positive dataset is somewhat different from the negative dataset.

Data can access in several ways.  

1. Files can be imported from the computer 
2. Using direct link. 


## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>      https://zenodo.org/record/4111092/files/ACPs.fasta
>      https://zenodo.org/record/4111092/files/non_ACPs.fasta
>    ```
>
>    {% include snippets/import_via_link.md %}
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

## Machine Learning Modeling of anti-cancer peptides 

Biological molecules such as proteins,  peptides, DNA, and RNA can be represented by the sequences, and based on these sequences various properties can be calculated.  These properties can be related to various biological properties of these sequences.  Properties associated with a peptide sequence such as overall charge, hydrophobicity profile, or k-mer composition can be utilized to build a machine learning model and predict the biological properties of unknown peptides. Finding anticancer peptides (ACPs) through wet-lab methods is costly and time-consuming; thus, the development of an efficient computational approach is useful to predict potential ACP peptides before wet-lab experimentation.  
![Alternative text](../../images/PDAUG_ML_1.png "ML algorithms use descriptors, a numerical representation of a sequence-based 
properties, for model building 
")

In this tutorial, we used 6 different machine learning algorithms for the prediction of ACPs using the peptide sequence-based CTD (Composition transition and distribution descriptors)descriptors or features. 
All the models were trained using our dataset that combines 138 Anticancer peptides and 138 non-anticancer peptides. We applied 10 fold cross validation on this data set.


## Calculating Peptide Descriptors 

A descriptor or feature is the quantitative or a qualitative measure of a property that is associated with a sequence. For example, a chemical compound can be described via its charge chemical formula, molecular weight, number of rotatable bonds, etc. Similarly, several properties can be deduced from the biological sequence that can be utilized to summarise a biological property such as anti-cancer activity.    In this example, we utilized CTD descriptors, describes as composition, transition, and distribution descriptors. 


- **Calculating descriptor for ACPs.** 

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Peptide Sequence Descriptors** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *""*: `output` (Input dataset)
>    - *"DesType"*: `CTD`
>
>
{: .hands_on}

- **Calculating descriptor for non-ACPs.**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Peptide Sequence Descriptors** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *""*: `output` (Input dataset)
>    - *"DesType"*: `CTD`
>
>
>
{: .hands_on}


## Adding class label in training data


The class label usually describes samples from two different groups, in our case ACPs and non-ACPs.   In the case of binary classification, usually, samples are labeled as "0" or 1. 


- **We marked ACPs as "1" or positive data.** 

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Add Class Label** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `output1` (output of **Peptide Sequence Descriptors** {% icon tool %})
>    - *"Class Label"*: `1`
>
>
>
{: .hands_on}

- **We marked non-ACPs as "0" or negative samples.** 

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Add Class Label** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `output1` (output of **Peptide Sequence Descriptors** {% icon tool %})
>
>
>
{: .hands_on}


## Merging data frame to combine Negative and Positive data set

In previous steps, we have calculated descriptors and labeled the data as positive and negative for ACPs and non-ACPs respectively, now we can use this dataset as a training dataset, however before to final step we have to merge both the dataset as one that represents our final training dataset with the class labels.  

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Merge dataframes** {% icon tool %} with the following parameters:
>    - {% icon param-files %} *"Input files"*: `OutFile1` (output of **Add Class Label** {% icon tool %}), `OutFile1` (output of **Add Class Label** {% icon tool %})
>
>
>
{: .hands_on}


## Applying 6 different ML algorithms on the training data set

In this step, we will apply six ML algorithms (LRC, RFC, GBC, DTC, SGDC & SVMC) with  10 fold cross-validation on the training data.   In cross-validation, positive and negative data are randomly divided into 10  parts each set has the 10th part of active as well as inactive peptides. The algorithm was trained on the 9 sets and the prediction was made on the remaining 10th set. This process was repeated for every set. Thus the final performance scores are calculated as a mean of all the folds. We used z-scaling to normalize the data before ML modeling. The entire workflow was applied to the four descriptor sets and accuracy was estimated based on accuracy, precision, recall, f1, and AUC.  


- **Linear regression Classifier (LRC)** 

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **ML Models** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `output1` (output of **Merge dataframes** {% icon tool %})
>    - *"Select Machine Learning algorithms"*: `LRC`
>        - *"Select advanced parameters"*: `No, use program defaults.`
>    - *"Choose the Test method"*: `Internal`
>
>
{: .hands_on}

- **Random Forest Classifier (RFC)**

> ### {% icon hands_on %} Hands-on: Task description
>
> 2. **ML Models** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `output1` (output of **Merge dataframes** {% icon tool %})
>    - *"Select Machine Learning algorithms"*: `RFC`
>        - *"Specify advanced parameters"*: `No, use program defaults.`
>    - *"Choose the Test method"*: `Internal`
>
>
{: .hands_on}


- **Gradient boosting Classifier (GBC)**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **ML Models** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `output1` (output of **Merge dataframes** {% icon tool %})
>    - *"Select Machine Learning algorithms"*: `GBC`
>        - *"Specify advanced parameters"*: `No, use program defaults.`
>    - *"Choose the Test method"*: `Internal`
>
>
>
{: .hands_on}


- **Decision tree (DTC)**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **ML Models** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `output1` (output of **Merge dataframes** {% icon tool %})
>    - *"Select Machine Learning algorithms"*: `DTC`
>        - *"Specify advanced parameters"*: `No, use program defaults.`
>    - *"Choose the Test method"*: `Internal`
>
{: .hands_on}


- **Stochastic gradient descent Classifier (SGDC)**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **ML Models** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `output1` (output of **Merge dataframes** {% icon tool %})
>    - *"Select Machine Learning algorithms"*: `SGDC`
>        - *"Specify advanced parameters"*: `No, use program defaults.`
>    - *"Choose the Test method"*: `Internal`
>
>
{: .hands_on}



- **Support vector machine Classifier (SVMC)**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **ML Models** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `output1` (output of **Merge dataframes** {% icon tool %})
>    - *"Select Machine Learning algorithms"*: `SVMC`
>        - *"Specify advanced parameters"*: `No, use program defaults.`
>    - *"Choose the Test method"*: `Internal`
>
{: .hands_on}


## Merging results in one file

In previous steps we have trained the machine learning models, these models return a TSV  that captures performance measures of these algorithms. We used the Marge Data Frame tool to combine these results as one file in this step.  

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Merge dataframes** {% icon tool %} with the following parameters:
>    - {% icon param-files %} *"Input files"*: `output1` (output of **ML Models** {% icon tool %}), `output1` (output of **ML Models** {% icon tool %}), `output1` (output of **ML Models** {% icon tool %}), `output1` (output of **ML Models** {% icon tool %}), `output1` (output of **ML Models** {% icon tool %})
>
>
{: .hands_on}


## Creating a final heat map to assess the results 

In the final step, a heat map will be generated which represents performance measures of various algorithms. We applied five different performance measures, accuracy, recall, F1-score, precision, and mean AUC (Area Under Curve) score.    


> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Heat Map** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *""*: `output1` (output of **Merge dataframes** {% icon tool %})
>
>
{: .hands_on}

![Alternative text](../../images/PDAUG_ML_2.png "ML algorithms use descriptors, a numerical representation of a sequence-based 
properties, for model building 
")

The brighter yellow color shows high-performance while the blue color shows a lower score.  Heart map suggests that algorithms GBC, LRC, and SVMC show high performs in comparison to the other three.   DTC shows an intermediate performance while RFC and SGDC performed poorly on this data set. 


# Conclusion
{:.no_toc}

In this tutorial, we learn how to utilize the quantitative properties of peptide sequences and apply the machine learning algorithms to predict the biological properties of the peptide sequence. 

![Alternative text](../../images/MLWorkFlow.png " Work flow ")