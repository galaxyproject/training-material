---
layout: tutorial_hands_on

title: 'Mass spectrometry: LC-MS analysis'
zenodo_link: ''
questions:
- What are the main steps of untargetted LC-MS data processing for metabolomic analysis?
- How to conduct metabolomic data analysis from preprocessing to annotation using Galaxy?
objectives:
- To comprehend the diversity of LC-MS metabolomic data analysis.
- To get familiar with the main steps constituting a metabolomic workflow for untargetted LC-MS analysis.
- To evaluate the potential of a workflow approach when dealing with LC-MS metabolomic data.
time_estimation: '3h'
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- jfrancoismartin
- lecorguille
- melpetera
- yguitton

---


# Introduction
{:.no_toc}

**TODO** Explain why metabolomics, what do you want to do

To illustrate this approach, we will use data from {% cite Thvenot2015 %}. The objectives of this paper was to analyze
the inﬂuence of age, body mass index, and gender on the urine metabolome. To do so, the authors collected samples
from 183 employees from the French Alternative Energies and Atomic Energy Commission (CEA) and performed LC-HRMS LTQ-Orbitrap
(negative ionization mode) (**TODO** explain the terms).

Since the original dataset takes a few hours to be processed, we chose to take a limited subset of individuals for this tutorial.
This will allow you to perform an example of metabolomic workflow from pre-processing to annotation in a limited time, even though
the results obtained may not be reliable from a scientific point of view due to a sample size way too small. Nevertheless,
the chosen diversity of sample will allow you to explore the bases of a metabolomic workflow.

We chose a subset of 9 samples, composed of 6 biological samples and 3 quality-control pooled samples (QC pools - mix of all
biological samples).

To analyze these data, we will the follow a light version of the [LC-MS workflow](http://workflow4metabolomics.org/the-lc-ms-workflow),
developed by the [Wokflow4metabolomics group](http://workflow4metabolomics.org/), ({% cite Giacomoni2014 %}, {% cite Guitton2017 %}).
**TODO** Introduce with one or two sentence the workflow (explanation of the meaning of LC-MS, the big steps, etc).
This workflow takes as input **TODO** and perform several steps: pre-processing, statistics and annotation.


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Preprocessing with XCMS

The first step in the workflow is the pre-processing of the raw data with XCMS ({% cite Smith2006 %}).

XCMS is a free software dedicated to pre-processing of any types of mass spectrometry acquisition files from low to
high resolution, including FT-MS data coupled with different kind of chromatography (liquid or gaz). This software is
used worldwide by a huge community of specialists in metabolomics using mass spectrometry methods.

This software is based on different algorithms that have been published, and is provided and maintained using R software [5,6,7].

XCMS is able to read files with open format as mzXML and netCDF which are independent of the constructors' formats.

It is composed of R functions able to extract, filter, align and fill gap, with the possibility to annotate isotopes,
adducts and fragments using the R package CAMERA. This set of functions gives modularity, thus being particularly well
adapted to define workflows, one of the key points of Galaxy:

![Preprocessing of the raw data with XCMS (in blue)](../../images/tutorial-lcms-data-import-run-workflow.png)


## Uploading your data into Galaxy

In metabolomics studies, the number of samples can vary a lot (from a few ones to hundreds). Thus, extracting your
data from the raw files can be very fast as well as take quite a long time. To optimise as much as possible the
computing time, W4M core team chose to propose modules that can run single raw files for the first steps of
pre-processing, since the initial actions in the extraction process treat files independantly.

Since the first steps can be run on each file, the use of **Dataset collection** is recommanded in Galaxy to avoid
launching jobs manually for each sample. You can consider the Dataset collection option from the very beginning, while
uploading your data into Galaxy.

> ### {% icon hands_on %} Hands-on: Data upload with **Get data**
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
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

Any comment needed here?

## Data preparation before XCMS steps: **MSnbase readMSData**

This first step is only meant to prepare your data for XCMS. It takes as input your raw files and
prepares RData files for the first XCMS step.

> ### {% icon hands_on %} Hands-on: MSnbase readMSData
>
> 1. **MSnbase readMSData** {% icon tool %} with the following parameters:
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > There is only one parameter for this module, corresponding to the input file.
>    {: .comment}
>
{: .hands_on}


> ### {% icon question %} Question
>
> With this single input parameter, what actions should I take before clicking on 'Execute'?
>
> > ### {% icon solution %} Solution
> >
> > 1. Choosing the correct type of input (here **Dataset collection**)
> > 2. Selecting the correct input (here **mzML**)
> >
> {: .solution}
>
{: .question}

## First XCMS step: **peak picking**

***TODO*** step introduction

> ### {% icon hands_on %} Hands-on: xcms findChromPeaks (xcmsSet)
>
> 1. **xcms findChromPeaks (xcmsSet)** {% icon tool %} with the following parameters:
>    - *"Extraction method for peaks detection"*: `CentWave - chromatographic peak detection using the centWave method`
>        - *"Max tolerated ppm m/z deviation in consecutive scans in ppm"*: `3`
>        - *"Min,Max peak width in seconds"*: `5,20`
>        - In *"Advanced Options"*:
>            - *"Prefilter step for for the first analysis step (ROI detection)"*: `3,5000`
>            - *"Noise filter"*: `1000`
>    - In *"Resubmit your raw dataset or your zip file"*:
>        - *"Resubmit your dataset or your zip file"*: `no need`
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

## Sub-step with **xcms plot chromatogram**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **xcms plot chromatogram** {% icon tool %} with the following parameters:
>    - In *"Resubmit your raw dataset or your zip file"*:
>        - *"Resubmit your dataset or your zip file"*: `no need`
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

## Sub-step with **xcms findChromPeaks Merger**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **xcms findChromPeaks Merger** {% icon tool %} with the following parameters:
>    - In *"Resubmit your raw dataset or your zip file"*:
>        - *"Resubmit your dataset or your zip file"*: `no need`
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

## Sub-step with **xcms groupChromPeaks (group)**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **xcms groupChromPeaks (group)** {% icon tool %} with the following parameters:
>    - *"Method to use for grouping"*: `PeakDensity - peak grouping based on time dimension peak densities`
>        - *"Bandwidth"*: `5.0`
>        - *"Width of overlapping m/z slices"*: `0.01`
>    - *"Get the Peak List"*: `Yes`
>    - In *"Resubmit your raw dataset or your zip file"*:
>        - *"Resubmit your dataset or your zip file"*: `no need`
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

## Sub-step with **xcms adjustRtime (retcor)**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **xcms adjustRtime (retcor)** {% icon tool %} with the following parameters:
>    - *"Method to use for retention time correction"*: `PeakGroups - retention time correction based on aligment of features (peak groups) present in most/all samples.`
>        - *"Minimum required fraction of samples in which peaks for the peak group were identified"*: `0.8299`
>        - *"Smooth method"*: `loess - non-linear alignment`
>    - In *"Resubmit your raw dataset or your zip file"*:
>        - *"Resubmit your dataset or your zip file"*: `no need`
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

## Sub-step with **xcms plot chromatogram**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **xcms plot chromatogram** {% icon tool %} with the following parameters:
>    - In *"Resubmit your raw dataset or your zip file"*:
>        - *"Resubmit your dataset or your zip file"*: `no need`
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

## Sub-step with **xcms groupChromPeaks (group)**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **xcms groupChromPeaks (group)** {% icon tool %} with the following parameters:
>    - *"Method to use for grouping"*: `PeakDensity - peak grouping based on time dimension peak densities`
>        - *"Bandwidth"*: `5.0`
>        - *"Width of overlapping m/z slices"*: `0.01`
>    - *"Get the Peak List"*: `Yes`
>        - *"Convert retention time (seconds) into minutes"*: `Yes`
>        - *"Number of decimal places for retention time values reported in ions' identifiers."*: `2`
>        - *"Replace the remain NA by 0 in the dataMatrix"*: `Yes`
>    - In *"Resubmit your raw dataset or your zip file"*:
>        - *"Resubmit your dataset or your zip file"*: `no need`
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

## Sub-step with **xcms fillChromPeaks (fillPeaks)**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **xcms fillChromPeaks (fillPeaks)** {% icon tool %} with the following parameters:
>    - In *"Peak List"*:
>        - *"Convert retention time (seconds) into minutes"*: `Yes`
>        - *"Number of decimal places for retention time values reported in ions' identifiers."*: `2`
>    - In *"Resubmit your raw dataset or your zip file"*:
>        - *"Resubmit your dataset or your zip file"*: `no need`
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

## Sub-step with **CAMERA.annotate**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **CAMERA.annotate** {% icon tool %} with the following parameters:
>    - In *"Annotate Isotopes [findIsotopes]"*:
>        - *"Max. ion charge"*: `2`
>    - *"Mode"*: `Only groupFWHM and findIsotopes functions [quick]`
>    - In *"Statistics and results export: [diffreport]"*:
>        - *"Number of condition"*: `One condition`
>    - In *"Export options"*:
>        - *"Convert retention time (seconds) into minutes"*: `Yes`
>        - *"Number of decimal places for retention time values reported in ions' identifiers."*: `2`
>    - In *"Resubmit your raw dataset or your zip file"*:
>        - *"Resubmit your dataset or your zip file"*: `no need`
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


# Data processing: quality checks, normalisation, data filtering

In the previous step of LC-MS workflow, you saw how to extract features from your acquisition files. This data is
shaped in a format allowing the use of various standard statistical methods. However, being able to perform a
statistical analysis does not mean necessarily being able to highlight relevant information. Indeed, data are often affected
by various sources of unwanted variability. It can limit the effectiveness of statistical methods, leading sometimes to
difficulties in revealing investigated effects. Identifying such variability can help analysing your data at its full potential.

In this tutorial, we chose to limit the data processing to 3 steps:
 - overview of the variability in the data
 - signal drift correction
 - filtering of unreliable variables based on coefficients of variation

> ### {% icon comment %} Comments
> To get a little more information regarding data processing, do not hesitate to visit the usemetabo.oc plateform:
> [Link to LC-MS processing step](https://usemetabo.org/courses/w4mlc-ms-processing)
{: .comment}


## Step 1: global variability in the data

Next thing I will fill (MP)

> ### {% icon hands_on %} Hands-on: Using **Quality Metrics** to get an overview of your data
>
> 1. **Quality Metrics** {% icon tool %} with the following parameters:
>    - *"Coefficient of Variation"*: `no`
>    - *"Advanced parameters"*: `Use default`
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

## Step 2: handling the signal drift observed althrough the analytical sequence

> ### {% icon hands_on %} Hands-on: Data normalisation using the **Batch_correction** module
>
> 1. **Batch_correction** {% icon tool %} with the following parameters:
>    - *"Type of regression model "*: `linear`
>        - *"Factor of interest "*: `gender`
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

## Step 3: getting rid of unreliable variable using CV

> ### {% icon hands_on %} Hands-on: CV calculation using the **Quality Metrics** module
>
> 1. **Quality Metrics** {% icon tool %} with the following parameters:
>    - *"Coefficient of Variation"*: `no`
>    - *"Advanced parameters"*: `Use default`
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



> ### {% icon hands_on %} Hands-on: Data filtering using the **Generic_Filter** module
>
> 1. **Generic_Filter** {% icon tool %} with the following parameters:
>    - *"Deleting samples and/or variables according to Numerical values"*: `yes`
>        - In *"Identify the parameter to filter "*:
>            - {% icon param-repeat %} *"Insert Identify the parameter to filter "*
>                - *"On file"*: `Variable metadata`
>                - *"Name of the column to filter"*: `poolCV_over_sampleCV`
>                - *"Interval of values to remove"*: `upper`
>                    - *"Remove all values upper than"*: `1.0`
>            - {% icon param-repeat %} *"Insert Identify the parameter to filter "*
>                - *"On file"*: `Variable metadata`
>                - *"Name of the column to filter"*: `pool_CV`
>                - *"Interval of values to remove"*: `upper`
>                    - *"Remove all values upper than"*: `0.3`
>    - *"Deleting samples and/or variables according to Qualitative values"*: `yes`
>        - In *"Removing a level in factor"*:
>            - {% icon param-repeat %} *"Insert Removing a level in factor"*
>                - *"Name of the column to filter"*: `sampleType`
>                - *"Remove factor when"*: `pool`
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




# Statistical analysis to find variables of interest


## Computation of statistical indices

> ### {% icon hands_on %} Hands-on: Statistical analysis using the **Univariate** module
>
> 1. **Univariate** {% icon tool %} with the following parameters:
>    - *"Factor of interest"*: `bmi`
>    - *"Test"*: `Spearman correlation rank test (quantitative)`
>    - *"Method for multiple testing correction"*: `none`
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

## Reduction of the dataset to variables of interest (to do: rephrase this to something better)

> ### {% icon hands_on %} Hands-on: Variable filtering using the **Generic_Filter** module
>
> 1. **Generic_Filter** {% icon tool %} with the following parameters:
>    - *"Deleting samples and/or variables according to Numerical values"*: `yes`
>        - In *"Identify the parameter to filter "*:
>            - {% icon param-repeat %} *"Insert Identify the parameter to filter "*
>                - *"On file"*: `Variable metadata`
>                - *"Name of the column to filter"*: `bmi_spearman_none`
>                - *"Interval of values to remove"*: `upper`
>                    - *"Remove all values upper than"*: `0.05`
>            - {% icon param-repeat %} *"Insert Identify the parameter to filter "*
>                - *"On file"*: `Variable metadata`
>                - *"Name of the column to filter"*: `bmi_spearman_cor`
>                - *"Interval of values to remove"*: `between`
>                    - *"Remove all values between"*: `-0.9`
>                    - *"And"*: `0.9`
>    - *"Deleting samples and/or variables according to Qualitative values"*: `no`
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


# Annotation


> ### {% icon hands_on %} Hands-on: Annotating the data using the HMDB
>
> 1. **HMDB MS search** {% icon tool %} with the following parameters:
>    - *"Would you use a file "*: `YES`
>        - *"Do you have a header "*: `YES`
>        - *"Column of masses "*: `c3`
>    - *"Mass-to-charge ratio "*: `0.005`
>    - *"Number of maximum entries returned by the query "*: `3`
>    - *"Molecular Species "*: `Negatif Mode`
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
