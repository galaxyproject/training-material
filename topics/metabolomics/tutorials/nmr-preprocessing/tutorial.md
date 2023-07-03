---
layout: tutorial_hands_on

title: Nuclear Magnetic Resonance: data preprocessing
zenodo_link: ''
questions:
- What are the main steps of untargeted 1H-NMR data preprocessing for metabolomic analyses?
- How to conduct 1H-NMR-based metabolomic data preprocessing using Galaxy?
objectives:
- To understand the steps necessary to perform preprocessing of untargeted 1H-NMR-based metabolomic data.
- To get familiar with the use of PEPS-NMR-based Galaxy modules dedicated to 1H-NMR data preprocessing.
- To evaluate the potential of a workflow approach when dealing with 1H-NMR data preprocessing.
time_estimation: 2H
key_points:
- To process untargeted 1H-NMR based metabolomic data preprocessing, various steps are needed.	
- Resources are available in Galaxy, but do not forget that you need appropriate knowledge to perform a relevant preprocessing.
contributors:
- mtremblayfr
- Cécile Canlet
- ManonMartin
---

# Introduction

<!-- This is a comment. -->

Metabolomics is an -omic science known for being one of the most closely related to phenotype. It focuses on studying the very small molecules which are called metabolites, to better understand matters linked to the metabolism. It involves the study of different types of matrices, such as blood, urine, tissues, in various organisms including plants and humans.
One of the three main technologies used to perform metabolomic analyses is Proton Nuclear Magnetic Resonance (1H NMR). Data analysis for this technology requires several steps, ranging from Fourier transform of Free Induction Decay (raw spectra) to statistical analysis and annotation. But, an inadequate preprocessing will never be compensated by a powerful data analysis. Some crucial steps must be carefully performed before the statistical analysis. Moreover, metabolomics studies can have several hundreds of samples. Manual preprocessing can be very time-consuming. 

To be able to perform semi-automatically a complete 1H NMR analysis with advanced methods (solvent suppression, baseline correction, etc.) in a single environment, the Wokflow4Metabolomics team provides Galaxy tools dedicated to metabolomics. This tutorial details the steps involved in the first part of untargeted 1H-NMR data processing: extracting information from FID data to obtain what is called a peak table. This step is commonly refered to as the preprocessing step. This tutorial will show you how to perform such a step using the Galaxy implementation of the PepsNMR R package (% cite Martin2018 %).

To illustrate this approach, we will use data from {% cite EscribanoVasquez2018 %}. One of the objectives of this work was to assess the influence of microbiota and high fat diet on the urinary metabolome. To analyze these data, we will then follow a Galaxy workflow developed by the Wokflow4metabolomics group ({% cite Giacomoni2014 %}, {% cite Guitton2017 %}).

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Overview

NMR data preprocessing is based on the PepsNMR R package (% cite Martin2018 %) (see slide 2). It covers several steps included in two tools (NMR_Read and NMR_Preprocessing, see slide 3).

In this tutorial we will focus on xx main steps:
- xx
- xx
- xx

# Data description

In this tutorial we will use the [AAP urine dataset](https://theses.hal.science/tel-02866073) generated in the lab of Claire Cherbuy (Micalis Institute, Jouy-en-Josas, France).
The intestinal microbiota is involved in the regulation of several metabolic pathways of the host, leading to important host-microbiota interaction axes such as those involved in metabolic signalling and immune/inflammatory responses. This microbial community has the enzymatic machinery necessary to metabolize nutriments coming from diet, and is a key factor of the host's energetic metabolism. Fat
consumption has increased considerably in recent decades. This high consumption of fat is harmful for health. On a high fat diet, the intestinal microbiota is in a dysbiotic state.
The objective of this study was to investigate the impact of a bacterial species, E. coli, which increases during fat consumption, on the metabolomic trajectory of mono-associated mice fed a standard and high fat diet. Two strains of E. coli were used in our study: the CEC15 strain we previously isolated from freshly pooled faecal samples of 15-day-old suckling rodents and the Nissle 1917 strain, as a representative member of the B2 E. coli phylogroup that has shown to increase in the human gut microbiote during the last decades.

To this aim, the authors collected samples from 59 male C57b mice divided in four microbiota groups (Germ Free / Germ Free+CEC 15 / Germ Free+Nissle1917 / Conventionally raised mice) and two diets (Normal / High fat) and performed 1-H NMR analysis.

# Get data

>    > <comment-title>On the supported files</comment-title>
>    >
>    > Only Bruker files are currently supported in the W4M plateform.
>    {: .comment}

The first step is to upload files into your Galaxy history. Bruker files have to be included in a zip file. Two directory structures are possible. FID can be 'structured' into sub-directories as illustrated in Figure 1, or not. You can upload data your own zip file from the 'Upload Data' button or you can used a shared zip file. 

![Figure 1: Example of FID's organization in directories and sub-directories](../../images/tutorial-nmr-read-sub-true.png)

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history and give it a name.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import your own `zip` file 
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> or create an history from a shared history
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md astype="" %} (https://workflow4metabolomics.usegalaxy.fr/u/mtremblayfranco/h/aap-urine-1)

> > <tip-title>Comment to W4M users</tip-title>
> > If you are a W4M user, please note that you can find at the following link a ready-to-start history:
> > [GTN_NMRpreprocessing](https://workflow4metabolomics.usegalaxy.fr/u/mtremblayfranco/h/gtnnmrpreprocessing).
> > We highly recommend to get started by importing this history.
> > this step corresponds to the dataset 1.
> {: .tip}

{: .hands_on}

You should have in your history a green zip file (`AAP_Urine`).

# Reading the zip file with the NMR_Read tool

The NMR_Read tool is used to read 'fid' files included in the zip file you previously added into your history. Parameters of this tool are linked with the possible options to organize fid files to preprocess (see slides 4 and 5).
This tool will return a FID data matrix saved as `NMR_Read_datamatrix` and a samplemetadata matrix saved as `NMR_Read_sampleMetadata`.

The dataMatrix file is a table containing intensities of NMR chemical shifts (variables, in rows) for every samples (in columns). The first column is for variables’ identifiers while the first row is for samples’ identifiers. 

The sampleMetadata file is a table containing information about the samples. The first column is for samples’ identifiers. Samples' identifiers of these two matrices are unique and identical which is manadatory for statistical analyses. 

You can add columns for analytical and biological information such as biological groups of interest for the supposed study (Sex, Exposition, Diet, Pathology...) to this table. If needed (further statistical analysis), you have to download the generated table in your local environment. For this, You can look at [the dedicated training material]({% link topics/galaxy-interface/tutorials/download-delete-data/tutorial.md %}).

> <hands-on-title> Read raw NMR fids </hands-on-title>
>
> 1. {% tool [NMR_Read](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Read/3.3.0) %} with the following parameters:
>    - *"Bruker FID file"*: `AAP_Urine`
>    - *"Presence of subdirectories?"*: ` TRUE `
>    - *"Use (sub)directories names as FID names?"*: ` TRUE `
>
>    > <comment-title> NMR_Read parameters </comment-title>
>    {% tool [NMR_Read](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    > If use of _title_ file and presence of sub-directories: set the FID Title line, `subdirs = TRUE`,  `dirs.names = FALSE`
>    > If use of _title_ file and no sub-directories: set the FID Title line, `subdirs = FALSE`,  `dirs.names = FALSE`
>    > If no use of _title_ file and presence of sub-directories: `subdirs = TRUE`,  `dirs.names = TRUE`
>    > If no use of _title_ file and no sub-directories: `subdirs = FALSE`,  `dirs.names = TRUE`
>    {: .comment}
>
>    > <comment-title> Comment to W4M users </comment-title>
>    >
>    > In the [AAP Urine](https://workflow4metabolomics.usegalaxy.fr/u/mtremblayfranco/h/gtnnmrpreprocessing) history,
>    > this step corresponds to the datasets number 3 to 6.
>    {: .comment}
>
{: .hands_on}

# Preprocessing fids with the NMR_Preprocessing tool

Your data are now ready for preprocessing. This step can be done with the NMR_Preprocessing tool. It will produce a datamatrix (`NMR_Preprocessing_dataMatrix`) including the preprocessed fid with row corresponding to chemical shifts and columns to samples, a variablemetadata file (`NMR_Preprocessing_variableMetadata`) including information on variables (NMR chemical shifts) and a graph with intermediate (depending on the graph option chosen for each sub-step) and final preprocessed spectra. These spectra will help you to have an idea about the effect of chosen parameters on the preprocessing step. 

The NMR_preprocessing tool includes several steps as described in Figure 2.

![Figure 2: Steps of NMR spectra preprocessin](../../images/tutorial-nmr-workflow.png)

## Group delay correction with **NMR_Preprocessing**

This step corresponds to the 1st order phase correction (see slides 6 and 7). 

## Solvent supression with **NMR_Preprocessing**

This step is explained in slides 8 and 9. `Smoothing parameter` determines how smooth is the solvent signal. Figure 3 below illustrates the effect of the `Smoothing parameter` spectra. Spectra have been obtained for a sample obtained with `Smoothing parameter` values of 1, 10^6 and 1^9.

![Figure 3: Effect of the `Smoothing parameter` in the Solvent suppression step](../../images/tutorial-nmr-workflow-solventsuppression.png)

> <hands-on-title> Effect of `Smoothing parameter` on signal intensity </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - {% icon param-file %} *"Data matrix of FIDs"*: the `dataMatrix` file
>    - {% icon param-file %} *"Sample metadata file"*: the `sampleMetadata` file
>    - In *"Solvent Suppression"*: the lambda smoother used to penalized the non-parametric estimation of the solvent signal
>        - *"Solvent Suppression: Smoothing parameter"*: `1.0`
>        - *"Display the FIDs after solvent suppression?"*: `yes`

> You can leave the other parameters with their default values.
>
> <question-title></question-title>
> 
> Based on explanations given in slide 8:
> 1. Which graph (Left/Middle/Right) corresponds to a value 1?
> 2. Which graph (Left/Middle/Right) corresponds to a value 10^6?
> 3. Which graph (Left/Middle/Right) corresponds to a value 10^9?
>
> > <solution-title></solution-title>
> > Explanations CCA
> > 1. 10^6
> > 2. 1
> > 3. 10^9
> >
> {: .solution}
>
{: .question}
>
>    > <comment-title> Comment to W4M users </comment-title>
>    >
>    > In the [AAP Urine](https://workflow4metabolomics.usegalaxy.fr/u/mtremblayfranco/h/gtnnmrpreprocessing) history, this step corresponds to the datasets number 11 to 14.
>    > You can also run this tool with default value (1000000.0, datasets 7 to 10) and 10^9 (datasets 15 to 18)
>    {: .comment}
>
{: .hands_on}

## Apodization **NMR_Preprocessing**

This step aims at improving the sensitivity by multiplying the FID by a factor called a weighting. Several classes of factors are available in W4M (Negative exponential, Gaussian, Hanning, Hamming, 
Cos2). You can refers to slides 8 and 9. 

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - {% icon param-file %} *"Data matrix of FIDs"*: the `dataMatrix` file
>    - {% icon param-file %} *"Sample metadata file"*: the `sampleMetadata` file
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>        - *"Apodization: Line broadening"*: `5`

> 2. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - {% icon param-file %} *"Data matrix of FIDs"*: the `dataMatrix` file
>    - {% icon param-file %} *"Sample metadata file"*: the `sampleMetadata` file
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>        - *"Apodization: Line broadening"*: `0.3`

> You can leave the other parameters with their default values.

> <question-title></Effect of the weighting function >
> Based on the "2016-61-UR-N4-CD" spectrum, which is the effet of the line broadening value on FID?
>
> > <solution-title></solution-title>
> >
> > CCA
> {: .solution}
>
{: .question}
> Several datasets are available in the history [AAP Urine](https://workflow4metabolomics.usegalaxy.fr/u/mtremblayfranco/h/gtnnmrpreprocessing), that corresponds to different "methods" and "Line broadening" values:
> - Method: negative exponential and Line broadening: 0.3 = datasets 19 - 22
> - Method: negative exponential and Line broadening: 5.0 = datasets 23 - 26
> - Method: cos2 and Phase: 0.0                           = datasets 27 - 30
> - Method: gauss and Line broadening: 5.0                = datasets 31 - 34
> - Method: hamming and Phase: 0.0                        = datasets 35 - 28
> 
{: .hands_on}
	
## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Solvent Suppression"*:
>        - *"Solvent Suppression: Smoothing parameter"*: `1000000000.0`
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `nearvalue`
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
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

## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>            - *"Line broadening"*: `0.3`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `nearvalue`
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
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

## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>            - *"Line broadening"*: `5.0`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `nearvalue`
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
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

## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `cos2`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `nearvalue`
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
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

## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `gauss`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `nearvalue`
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
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

## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `hamming`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `nearvalue`
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
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

## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `nearvalue`
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
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

## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: method"*: ` max `
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `nearvalue`
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
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

## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `nearvalue`
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
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

## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `all`
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
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

## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `window`
>            - In *"Search_zone"*:
>                - {% icon param-repeat %} *"Insert Search_zone"*
>                    - *"Search zone: left border"*: `2.0`
>                    - *"Search zone: right border"*: `-2.0`
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
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

## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `nearvalue`
>        - *"Shift Referencing: shiftHandling"*: ` cut `
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
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

## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `nearvalue`
>        - *"Shift Referencing: shiftHandling"*: ` circular `
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
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

## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `nearvalue`
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: smoothing parameter"*: `100.0`
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
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

## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `nearvalue`
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: smoothing parameter"*: `1000000000.0`
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
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

## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `nearvalue`
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: asymmetry parameter"*: `0.001`
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
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

## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `nearvalue`
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: asymmetry parameter"*: `0.5`
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
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

## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `nearvalue`
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: asymmetry parameter"*: `1.0`
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
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

## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `nearvalue`
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
>    - In *"Negative intensities to Zero"*:
>        - *"Set negative intensities to zero?"*: ` NO `
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
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           >
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

## Sub-step with **NMR_Preprocessing**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `nearvalue`
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: exclusion area(s)"*: ` NO `
>    - In *"Negative intensities to Zero"*:
>        - *"Set negative intensities to zero?"*: ` NO `
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
