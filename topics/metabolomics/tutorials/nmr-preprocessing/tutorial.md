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

Metabolomics is an -omic science known for being one of the most closely related to phenotypes. It focuses on studying the very small molecules which are called metabolites, to better understand matters linked to the metabolism. It involves the study of different types of matrices, such as blood, urine, tissues, in various organisms including plants and humans.
One of the three main technologies used to perform metabolomic analyses is Proton Nuclear Magnetic Resonance (1H NMR). Data analysis for this technology requires several steps, ranging from Fourier transform of Free Induction Decay (raw spectra) to statistical analysis and annotation. To be able to perform a complete 1H NMR analysis in a single environment, the Wokflow4Metabolomics team provides Galaxy tools dedicated to metabolomics. This tutorial details the steps involved in the first part of untargeted 1H-NMR data processing: extracting information from FID data to obtain what is called a peak table. This step is commonly refered to as the preprocessing step. This tutorial will show you how to perform such a step using the Galaxy implementation of the PEPSNMR R package (% cite Martin2018 %).
To illustrate this approach, we will use data from {% cite EscribanoVasquez2018 %}. One of the objectives of this work was to assess the influence of microbiota and high fat diet on the urinary metabolome. To analyze these data, we will then follow a Galaxy workflow developed by the Wokflow4metabolomics group ({% cite Giacomoni2014 %}, {% cite Guitton2017 %}).

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

# Overview

NMR data preprocessing covers several steps included in two tools (NMR_Read and NMR_Preprocessing). The NMR_preprocessing tool includs several steps as described in Figure 1. 
[![The full tutorial workflow](../../images/tutorial-nmr-workflow.png)](../../images/tutorial-nmr-workflow.png)

In this tutorial we will focus on xx main steps:
- xx
- xx
- xx

# Data description

In this tutorial we will use the [AAP urine dataset](https://theses.hal.science/tel-02866073) generated in the lab of Claire Cherbuy (Micalis Institute, Jouy-en-Josas, France).
The intestinal microbiota is involved in the regulation of several metabolic pathways of the host, leading to important host-microbiota interaction axes such as those involved in metabolic signalling and immune/inflammatory responses. This microbial community has the enzymatic machinery necessary to metabolize nutriments coming from diet, and is a key factor of the host's energetic metabolism. Fat
consumption has increased considerably in recent decades. This high consumption of fat is harmful for health. On a high fat diet, the intestinal microbiota is in a dysbiotic state.
The objective of this study was to investigate the impact of a bacterial species, E. coli, which increases during fat consumption, on the metabolomic trajectory of mono-associated mice fed a standard and high fat diet. Two strains of E. coli were used in our study: the CEC15 strain we previously isolated from freshly pooled faecal samples of 15-day-old suckling rodents and the Nissle 1917 strain, as a representative member of the B2 E. coli phylogroup
that has shown to increase in the human gut microbiote during the last decades.

To do so, the authors collected samples from 59 male C57b mice divided in four microbiota groups (Germ Free / Germ Free+CEC 15 / Germ Free+Nissle1917 / Conventionally raised mice) and two diets (Normal / High fat) and performed 1-H NMR analysis.

# Get data

>    > <comment-title>on the supported files</comment-title>
>    >
>    > Only Bruker files are currently supported in the W4M plateform.
>    {: .comment}

The first step is to upload files into Galazy. Bruker files have to be included in a zip file. Two directory structures are possible. FID can be 'structured' into subdirectories as illustrated in Figure 2 or not.

[!['Structured' files](../../images/tutorial-nmr-read-sub-true.png)](../../images/tutorial-nmr-read-sub-true.png)

> <hands-on-title>Data upload with <b>Get data</b></hands-on-title>
>
> 1. Create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the AAP urine `zip` file from a shared data library (https://workflow4metabolomics.usegalaxy.fr/u/mtremblayfranco/h/aap-urine)
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md astype="" %}
{: .hands_on}

You should have in your history a green zip file (`AAP2_Urine`).

# Reading the zip file with the NMR_Read tool

The NMR_Read tool is used 'fid' files included in the zip file you previously uploaded or in the shared history (https://workflow4metabolomics.usegalaxy.fr/u/mtremblayfranco/h/aap-urine). Parameters of this tool are linked with the possible options of the 'zip' file.
This tool will return a FID data matrix (including the raw spectra, saved as `NMR_Read_datamatrix`) and a samplemetadata matrix with information about these FIDs (NMR parameters acquisition, saved as `NMR_Read_sampleMetadata`).
  
> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [NMR_Read](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Read/3.3.0) %} with the following parameters:
>    - *"Presence of subdirectories?"*: ` TRUE `
>    - *"Use (sub)directories names as FID names?"*: ` TRUE `
>
>    > <comment-title> NMR_Read parameters </comment-title>
>    > If use of _title_ file and presence of sub-directories: set the FID Title line, `subdirs = TRUE`,  `dirs.names = FALSE`
>    > If use of _title_ file and no sub-directories: set the FID Title line, `subdirs = FALSE`,  `dirs.names = FALSE`
>    > If no use of _title_ file and presence of sub-directories: `subdirs = TRUE`,  `dirs.names = TRUE`
>    > If no use of _title_ file and no sub-directories: `subdirs = FALSE`,  `dirs.names = TRUE`
>    {: .comment}
>
{: .hands_on}



# Preprocessing fids with the NMR_Preprocessing tool

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
>    - In *"Solvent Suppression"*:
>        - *"Solvent Suppression: Smoothing parameter"*: `1.0`
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
