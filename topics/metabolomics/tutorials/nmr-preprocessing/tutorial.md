---
layout: tutorial_hands_on

title: 'Nuclear Magnetic Resonance Data preprocessing'
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
- To process untargeted 1H-NMR based metabolomic data, various steps are needed.	
- Resources are available in Galaxy, but do not forget that you need appropriate knowledge to perform a relevant preprocessing.
contributors:
- mtremblayfr
- ccanlet
- ManonMartin
---

# Introduction

<!-- This is a comment. -->

Metabolomics is an -omic science known for being one of the most closely related to phenotype. It focuses on studying the very small molecules which are called metabolites, to better understand matters linked to the metabolism. It involves the study of different types of matrices, such as blood, urine, tissues, in various organisms including plants and humans.
One of the three main technologies used to perform metabolomic analyses is Proton Nuclear Magnetic Resonance (1H NMR). Data analysis for this technology requires several steps, ranging from Fourier transform of Free Induction Decay (raw spectra) to statistical analysis and annotation. But, an inadequate preprocessing will never be compensated by a powerful data analysis. Some crucial steps must be carefully performed before the statistical analysis. Moreover, metabolomics studies can have several hundreds of samples. Manual preprocessing can be very time-consuming. 

To be able to perform a semi-automatically complete 1H NMR data preprocessing with advanced methods (solvent suppression, baseline correction, etc.) in a single environment, the Wokflow4Metabolomics team provides Galaxy tools dedicated to metabolomics. This tutorial details the steps involved in the first part of untargeted 1H-NMR data processing: extracting information from FID data to obtain what is called a peak table. This step is commonly refered to as the preprocessing step. This tutorial will show you how to perform such a step using the Galaxy implementation of the PepsNMR R package (% cite Martin2018 %).

To illustrate this approach, we will use data from {% cite EscribanoVasquez2018 %}. One of the objectives of this work was to assess the influence of microbiota and high fat diet on the urinary metabolome. To analyze these data, we will then follow a Galaxy workflow developed by the Wokflow4metabolomics group ({% cite Giacomoni2014 %}, {% cite Guitton2017 %}).

Since sometimes a couple of pictures is worth a thousand words, you will find in the following slides some material to help
you to understand how the NMR_Preprocessing tool works:
[link to slides](../../tutorials/nmr-preprocessing/slides.html).
This document is refered to as "Check the next X slides" in the present training material.
As an example, [check the 1st slide](../../tutorials/nmr-preprocessing/slides.html#pepsnmr_rpackage)
for complementary material about PepsNMR R package.


> <agenda-title></agenda-title>
> 
> In this tutorial, we will cover:
> 
> 1. TOC
> {:toc}
> 
{: .agenda}

# Overview

NMR data preprocessing is based on the PepsNMR R package (% cite Martin2018 %). It covers several steps included in two tools (NMR_Read and NMR_Preprocessing). You can see the NMR preprocessing workflow available in W4M [check the first slide](../../tutorials/nmr-preprocessing/slides.html#nmrpreprocessing_workflow).

In this tutorial we will focus on xx main steps:
- xx
- xx
- xx

# Data description

In this tutorial we will use the [GTN_NMRpreprocessing dataset](https://theses.hal.science/tel-02866073) generated in the lab of Claire Cherbuy (Micalis Institute, Jouy-en-Josas, France).
The intestinal microbiota is involved in the regulation of several metabolic pathways of the host, leading to important host-microbiota interaction axes such as those involved in metabolic signalling and immune/inflammatory responses. This microbial community has the enzymatic machinery necessary to metabolize nutriments coming from diet, and is a key factor of the host's energetic metabolism. Fat
consumption has increased considerably in recent decades. This high consumption of fat is harmful for health. On a high fat diet, the intestinal microbiota is in a dysbiotic state.
The objective of this study was to investigate the impact of a bacterial species, E. coli, which increases during fat consumption, on the metabolomic trajectory of mono-associated mice fed a standard and high fat diet. Two strains of E. coli were used in this study: the CEC15 strain previously isolated by C. Cherbuy's group from freshly pooled faecal samples of 15-day-old suckling rodents and the Nissle 1917 strain, as a representative member of the B2 E. coli phylogroup that has shown to increase in the human gut microbiote during the last decades.

To this aim, the authors collected samples from 59 male C57b mice divided in four microbiota groups (Germ Free / Germ Free+CEC 15 / Germ Free+Nissle1917 / Conventionally raised mice) and two diets (Normal / High fat) and performed 1-H NMR analysis.

# Get data

>    > <comment-title>On the supported files</comment-title>
>    > 
>    > Only Bruker files are currently supported in the W4M plateform.
>    {: .comment}

The first step is to upload files into your Galaxy history. Bruker files have to be included in a zip file. Two directory structures are possible. FID can be organized into sub-directories or not. You can upload data your own zip file from the 'Upload Data' button or you can used a shared zip file. 

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
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md astype="" %}
> 
> > <tip-title>Comment to W4M users</tip-title>
> > If you are a W4M user, please note that you can find at the following link a ready-to-start history:
> > [GTN_NMRpreprocessing](https://workflow4metabolomics.usegalaxy.fr/u/mtremblayfranco/h/gtnnmrpreprocessing).
> > We highly recommend to get started by importing this history.
> > this step corresponds to the dataset 1.
> {: .tip}
{: .hands_on}

You should have in your history a green zip file (`AAP_Urine`).

# Reading the zip file with the NMR_Read tool

The NMR_Read tool is used to read 'fid' files included in the zip file you previously added into your history. Parameters of this tool are linked with the possible options to organize fid files to preprocess ([check the next slide](../../tutorials/nmr-preprocessing/slides.html#nmrpreprocessing_workflow)).
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
>    > 
>    > If use of _title_ file and presence of sub-directories: set the FID Title line, `subdirs = TRUE`,  `dirs.names = FALSE`
>    > 
>    > If use of _title_ file and no sub-directories: set the FID Title line, `subdirs = FALSE`,  `dirs.names = FALSE`
>    > 
>    > If no use of _title_ file and presence of sub-directories: `subdirs = TRUE`,  `dirs.names = TRUE`
>    > 
>    > If no use of _title_ file and no sub-directories: `subdirs = FALSE`,  `dirs.names = TRUE`
>    {: .comment}
> 
>    > <comment-title> Comment to W4M users </comment-title>
>    > 
>    > In the [GTN_NMRpreprocessing](https://workflow4metabolomics.usegalaxy.fr/u/mtremblayfranco/h/gtnnmrpreprocessing) history,
>    > this step corresponds to the datasets number 3 to 6.
>    {: .comment}
> 
{: .hands_on}

# Preprocessing fids with the NMR_Preprocessing tool

Your data are now ready for preprocessing. This step can be done with the NMR_Preprocessing tool. It will produce a datamatrix (`NMR_Preprocessing_dataMatrix`) including the preprocessed fid with row corresponding to chemical shifts and columns to samples, a variablemetadata file (`NMR_Preprocessing_variableMetadata`) including information on variables (NMR chemical shifts) and a graph with intermediate (depending on the graph option chosen for each sub-step) and final preprocessed spectra. These spectra will help you to have an idea about the effect of chosen parameters on the preprocessing step. 

The NMR_preprocessing tool includes several steps as described in Figure 2. This tutorial ends up after the "Negative values zeroing step."

![Figure 2: Steps of NMR spectra preprocessing](../../images/tutorial-nmr-workflow.png)

## 1. Group delay correction

Phase correction ([check the next 3 slides](../../tutorials/nmr-preprocessing/slides.html#phase_shifts), [](../../tutorials/nmr-preprocessing/slides.html#group_delay_correction), [](../../tutorials/nmr-preprocessing/slides.html#group_delay_correction_illustration) and Figure 3) is a very important adjustment that needs to be made to a spectrum. Phase of a signal is related to the amount of signal observed above and below the baseline. Phase correction works to provide a signal in pure-absorption mode, which means a signal totally above and/or totally below the baseline.

![Figure 3: Illustration of the Group Delay recorded before the signal acquisition](../../images/tutorial-nmr-workflow-firstorderphasecorrection.png)

Phase correction involves adjusting both zero (ph0, see 5th step) and first-order (ph1) phases. This step corresponds to the 1st order phase correction. The first-order phase shift is due to the presence of a digital filter. The first tens of points in the FID are not part of the recorded signal and are called the group delay. Since the phase shift differs across signals, it introduces a first order phase shift linearly related to frequency.

> <hands-on-title> 1st order phase correction </hands-on-title>
> 
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - {% icon param-file %} *"Data matrix of FIDs"*: the `dataMatrix` file
>    - {% icon param-file %} *"Sample metadata file"*: the `sampleMetadata` file
>    - In *"Group delay correction"*: 
>        - *"Display the FIDs after 1st order phase correction?"*: `yes`
> 
> You can leave all parameters with their default values.
> 
{: .hands_on}

## 2. Solvent supression

[Check the next 2 slides](../../tutorials/nmr-preprocessing/slides.html#solvent_suppression), [](../../tutorials/nmr-preprocessing/slides.html#solvent_suppression_illustration) for explanations on solvent suppression. Smoothing parameter determines how smooth is the solvent signal after solvent suppression. Figure 4 below illustrates the effect of the Smoothing parameter on spectra. Spectra have been obtained with Smoothing parameter values of `1`, `10^6` and`1^9`.

![Figure 4: Effect of the Smoothing parameter in the Solvent suppression step](../../images/tutorial-nmr-workflow-solventsuppression.png)

> <hands-on-title> Effect of `Smoothing parameter` on signal intensity </hands-on-title>
> 
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - {% icon param-file %} *"Data matrix of FIDs"*: the `dataMatrix` file
>    - {% icon param-file %} *"Sample metadata file"*: the `sampleMetadata` file
>    - In *"Solvent Suppression"*: the lambda smoother is used to penalized the non-parametric estimation of the solvent signal
>        - *"Solvent Suppression: Smoothing parameter"*: `1.0`
>        - *"Display the FIDs after solvent suppression?"*: `no`
> 
> You can leave other parameters with their default values.
> 
> <question-title></question-title>
> 
> Based on explanations given in [](../../tutorials/nmr-preprocessing/slides.html#solvent_suppression_illustration), :
> 1. Which graph (Left/Middle/Right) corresponds to a value 1,
> 2. Which graph (Left/Middle/Right) corresponds to a value 10^6, 
> 3. Which graph (Left/Middle/Right) corresponds to a value 10^9 
> in the Figure 4?
> 
> > <solution-title></solution-title>
> > 
> > The smaller the lambda value, the smoother the signals. The value lambda = 1 corresponds to the second spectrum: the 
> > solvent signal is well suppressed, but this value is too small: all metabolite signals are diminished. A compromise 
> > must therefore be found between suppressing the solvent signal and smoothing the metabolite signals. The value lamba = 
> > 10^9 corresponds to spectrum n°3 and lambda = 10^6 to spectrum n°1. There is very little difference between these 2 
> > spectra. Default value seems to be a good compromise.
> > 
> {: .solution}
> 
{: .question}
> 
>    > <comment-title> Comment to W4M users </comment-title>
>    > 
>    > In the [GTN_NMRpreprocessing](https://workflow4metabolomics.usegalaxy.fr/u/mtremblayfranco/h/gtnnmrpreprocessing) history, this step corresponds to the datasets number 11 to 14.
>    > You can also run this tool with default value (1000000.0, datasets 7 to 10) and 10^9 (datasets 15 to 18)
>    {: .comment}
> 
{: .hands_on}

## 3. Apodization

This step aims at improving the sensitivity by multiplying the FID by a factor called a weighting. Several classes of factors are available in W4M (Negative exponential, Gaussian, Hanning, Hamming, 
Cos2). [Check the next 2 slides](../../tutorials/nmr-preprocessing/slides.html#apodization), [](../../tutorials/nmr-preprocessing/slides.html#apodization_illustration). 

> <hands-on-title> Task description </hands-on-title>
> 
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - {% icon param-file %} *"Data matrix of FIDs"*: the `dataMatrix` file
>    - {% icon param-file %} *"Sample metadata file"*: the `sampleMetadata` file
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>        - *"Apodization: Line broadening"*: `5`
>        - *"Display the FIDs after solvent suppression?"*: `no`
>        - *"Display the FIDs after Apodization?"*: `no`
> 
> 2. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - {% icon param-file %} *"Data matrix of FIDs"*: the `dataMatrix` file
>    - {% icon param-file %} *"Sample metadata file"*: the `sampleMetadata` file
>    - In *"Apodization"*:
>        - *"Apodization: method"*: `exp`
>        - *"Apodization: Line broadening"*: `0.3`
>        - *"Display the FIDs after Apodization?"*: `no`
> 
> You can leave other parameters with their default values.
>
> <question-title></Effect of the weighting function >
> 
> Based on the "2016-61-UR-N4-CD" spectrum, which is the effet of the line broadening value on FID?
> 
> > <solution-title></Effect of the weighting function>
> > 
> > CCA
> {: .solution}
> 
{: .question}
> 
>    > <comment-title> Comment to W4M users </comment-title>
>    > 
>    > In the [GTN_NMRpreprocessing](https://workflow4metabolomics.usegalaxy.fr/u/mtremblayfranco/h/gtnnmrpreprocessing) history, several datasets are available to evaluate effects of parameter on preprocessed spectra:
> - Method: negative exponential and Line broadening: 0.3 = datasets 19 - 22
> - Method: negative exponential and Line broadening: 5.0 = datasets 23 - 26
> - Method: cos2 and Phase: 0.0                           = datasets 27 - 30
> - Method: gauss and Line broadening: 5.0                = datasets 31 - 34
> - Method: hamming and Phase: 0.0                        = datasets 35 - 38
>    {: .comment}
{: .hands_on}
	
## 4. Fourier transform

Nest step corresponds to conversion of the signal in the time domain into a spectrum in the frequency domain, performed in the "Fourier transform" step ([Check the next  slide](../../tutorials/nmr-preprocessing/slides.html#fourier_transform)).

> <hands-on-title> Fourier transform </hands-on-title>
> 
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - {% icon param-file %} *"Data matrix of FIDs"*: the `dataMatrix` file
>    - {% icon param-file %} *"Sample metadata file"*: the `sampleMetadata` file
>    - In *"Fourier transform"*:
>        - *"Display the FIDs after solvent suppression?"*: `no`
> You can leave other parameters with their default values.
> {: .hands_on}

## 5. Zero order phase correction

In the fifth step, correction of the zero order phase is applied ([Check the next  slide](../../tutorials/nmr-preprocessing/slides.html#zero_order_phase)). The zero-order phase arises because the relative phase of the transmitter pulse and receiver are offset. This results in a mixing of the desired real part of the spectrum with a portion of the corresponding imaginary part, so one side of the base of each peak is observed to dip below the baseline). Zero-order phase correction undoes this mixing. The zero-order phase shift affects all frequencies in the same way.

> <hands-on-title> Zero order phase correction </hands-on-title>
> 
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - {% icon param-file %} *"Data matrix of FIDs"*: the `dataMatrix` file
>    - {% icon param-file %} *"Sample metadata file"*: the `sampleMetadata` file
>    - In *"Zero Order Phase Correction"*:
>        - *"Zero Order Phase Correction: method"*: ` RMS `
>        - *"Zero Order Phase Correction: exclusion area(s)"*: ` YES `
>        - In *"Zero Order Phase Correction: exclusion area(s)"*:
>           - *"Exclusion zone: left border"*: ` 5.1 `
>           - *"Exclusion zone: right border"*: ` 4.5 `
> 
> You can leave other parameters with their default values.
> 
>    > <comment-title> Comment to W4M users </comment-title>
>    > 
>    > In the [GTN_NMRpreprocessing](https://workflow4metabolomics.usegalaxy.fr/u/mtremblayfranco/h/gtnnmrpreprocessing) history, several datasets are available to evaluate effects of the method applied for phase correction:
> - Method: RMS = datasets 39 - 42
> - Method: MAX = datasets 43 - 46
>    {: .comment}
{: .hands_on}

## 6. Shift Referencing

A known standard (called internal reference compound), TMS or TSP, is usually added to the samples to refine the scale calibration. In this step, the chemical shift is defined relative to this reference and a ppm value is attributed to the reference peak (usually 0 ppm) and spectra are aligned to this peak. 

> <hands-on-title> Shift Referencing </hands-on-title>
> 
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - {% icon param-file %} *"Data matrix of FIDs"*: the `dataMatrix` file
>    - {% icon param-file %} *"Sample metadata file"*: the `sampleMetadata` file
>    - In *"Shift Referencing"*:
>        - *"Shift Referencing: definition of the search zone"*: `window`
>        - In *"Shift Referencing: definition of the search zone"*
>            - In *"Search zone"*
>               - *"Search zone: left border"*: -2.0
>               - *"Search zone: right border"*: -2.0
>        - *"Shift Referencing: shiftHandling"*: `zerofilling`
>        - *"Shift Referencing: the value of the reference peak in ppm"*: `0.0`
> 
> <question-title></Effect of Search zone parameter>
> 
> Run the NMR_Preprocessing tool with `Nearvalue` and `window` as values for the "Shift Referencing:   definition of the search zone" parameter, and `-2.0` and `2.0` respectively for left and right borders for the `window` search zone. What do you observe on spectra obtained for individual "X2016.61.UR.N4.CD"?
> 
> > <solution-title></solution-title>
> > 
> > Shift towards the right for the `window` value
> > 
> {: .solution}
> 
{: .question}
> 
>    > <comment-title> How does Shift referencing work in NMR_Preprocessing? </comment-title>
>    > 
>    > The algorithm proposes two ways to locate the reference compound peak in each spectrum within a range of intensities (search zone: `nearvalue`, `all`, (user-defined) `window`): it selects either the maximum intensity or the first peak in the search range higher than a predefined threshold. [Check the next  slide](../../tutorials/nmr-preprocessing/slides.html#zero_order_phase) to see the impact of the Search zone parameter.
>    > 
>    {: .comment}
> 
>    > <comment-title> Comment to W4M users </comment-title>
>    > Several datasets are available in the history [GTN_NMRpreprocessing](https://workflow4metabolomics.usegalaxy.fr/u/mtremblayfranco/h/gtnnmrpreprocessing), that corresponds to different search zone and shiftHandling values:
>    > 
>    > - definition of the search zone: `nearvalue`; shiftHandling: `zerofilling` = datasets 47 - 50
>    > - definition of the search zone: `all`; shiftHandling: `zerofilling`       = datasets 51 - 54
>    > - definition of the search zone: `window`; Search zone: left border: `-2.0`; Search zone: right border: `2.0`; shiftHandling: `zerofilling `                                      = datasets 55 - 58
>    >  - definition of the search zone: `nearvalue`; shiftHandling: `cut`        = datasets 59 - 62
>    > - definition of the search zone: `nearvalue`; shiftHandling: `circular`    = datasets 63 - 66
>    {: .comment}
> 
{: .hands_on}


## 7. Baseline correction

To ensure successful integration, baseline should be flat with no distortion. Baseline artefacts have to be removed. These artefacts result from multiple sources, such as the presence of macromolecules, a not entirely linear electronic detection process or calibration errors from the 180° pulse. Function used in this step estimates and removes the smoothed baseline from the spectra, based on two parameters `smoothing` and `asymmetry`.

> <hands-on-title> Baseline correction </hands-on-title>
> 
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - {% icon param-file %} *"Data matrix of FIDs"*: the `dataMatrix` file
>    - {% icon param-file %} *"Sample metadata file"*: the `sampleMetadata` file
>    - In *"Baseline Correction"*:
>        - *"Baseline Correction: smoothing parameter"*: `100000`
>        - *"Baseline Correction:asymmetry parameter"*: `0.05`
>        - *"Baseline Correction: exclusion area(s)"*: `NO`
>
>    > <comment-title> How does Baseline correction work in NMR_Preprocessing? </comment-title>
>    > 
>    > Algorithm uses uses asymmetric least squares with a roughness penalty:
>    > 
>    > - smoothing (default value=`1e7`): the larger it is, the smoother will be. With smoothing=`0`, the baseline will be equal to the signal and the corrected signal will be zero.
>    > 
>    > asymetry (default value=`0.05`): the smaller it is, the less the smoother will try to follow peaks when it is under the function and the more it will try to be under the function.
>    > 
>    {: .comment}
> 
>    > <comment-title> Comment to W4M users </comment-title>
>    > Several datasets are available in the history [GTN_NMRpreprocessing](https://workflow4metabolomics.usegalaxy.fr/u/mtremblayfranco/h/gtnnmrpreprocessing), that       >    > corresponds to different "smoothing" and "asymetry" values:
>    > - smoothing parameter: `100`; asymetry parameter: `0.05`   = datasets 67 - 70
>    > - smoothing parameter: `10^9`; asymetry parameter: `0.05`  = datasets 71 - 74
>    > - smoothing parameter: `10^5`; asymetry parameter: `0.001` = datasets 75 - 78
>    > - smoothing parameter: `10^5`; asymetry parameter: `0.5`   = datasets 79 - 82
>    > - smoothing parameter: `10^5`; asymetry parameter: `1.0`   = datasets 83 - 86
>    {: .comment}
> 
{: .hands_on}

> <question-title></Effect of the Asymetry parameter>
> 
> Run the NMR_Preprocessing tool with `0.001`, `0.5` and `1.0` as values for the Baseline correction: asymmetry parameter. What do you observe on spectra ofbatined for individual X2016.61.UR.N4.CD" (compare also with the default value?
> 
> > <solution-title></solution-title>
> > 
> > CCA
> >
> {: .solution}
> 
{: .question}

## 8. Negative values zeroing

Despite the application of baseline and phase corrections, spectra may still have negative intensities at specific frequency values. These cannot be properly interpreted and can have bad impacts on  statistical analyses. This filter simply sets them to zero.

> <hands-on-title> Negative values zeroing </hands-on-title>
> 
> 1. {% tool [NMR_Preprocessing](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/nmr_preprocessing/NMR_Preprocessing/3.3.0) %} with the following parameters:
>    - In *"Negative intensities to Zero?"*:
>        - *"Set negative intensities to zero?"*: `YES`
> 
>    > <comment-title> Comment to W4M users </comment-title>
>    > Several datasets are available in the history [GTN_NMRpreprocessing](https://workflow4metabolomics.usegalaxy.fr/u/mtremblayfranco/h/gtnnmrpreprocessing), that corresponds to `NO` and `YES` values:
>    > 
>    > - Set negative intensities to zero?: `YES` = datasets 89 - 90
>    > - Set negative intensities to zero?: `NO`  = datasets 91 - 94
>    {: .comment}
> 
{: .hands_on}


# Conclusion

Spectroscopic data pre-processing is a keystone procedure for enhanced information recovery in metabolomics, and the NMR_Preprocessing tool, based on the PepsNMR R package which is the only existing package able to comprehensively address a very detailed series of pre-processing steps designed for 1D 1 H NMR data. Use of this tool can increase spectral repeatability, suggesting an overall better recovery of information, and it can enhance predictive power in a classification context.
