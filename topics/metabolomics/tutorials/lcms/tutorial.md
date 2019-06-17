---
layout: tutorial_hands_on

title: 'Mass spectrometry: LC-MS analysis'
zenodo_link: 'http://doi.org/10.5281/zenodo.3244991'
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


# Preprocessing with xcms

The first step in the workflow is the pre-processing of the raw data with xcms ({% cite Smith2006 %}).

xcms is a free and open source software dedicated to pre-processing of any types of mass spectrometry acquisition files from low to
high resolution, including FT-MS data coupled with different kind of chromatography (liquid or gaz). This software is
used worldwide by a huge community of specialists in metabolomics using mass spectrometry methods.

This software is based on different algorithms that have been published, and is provided and maintained using R software [5,6,7].

xcms is able to read files with open format as mzXML, mzMl, mzData and netCDF which are independent of the constructors' formats.

It is composed of R functions able to extract, filter, align and fill gap, with the possibility to annotate isotopes,
adducts and fragments using the R package CAMERA. This set of functions gives modularity, thus being particularly well
adapted to define workflows, one of the key points of Galaxy:

![Preprocessing of the raw data with xcms (in blue)](../../images/tutorial-lcms-data-import-run-workflow.png)


## Importing the LC/MS data into Galaxy

In metabolomics studies, the number of samples can vary a lot (from a few ones to hundreds). Thus, extracting your
data from the raw files can be very fast as well as take quite a long time. To optimise as much as possible the
computing time, W4M core team chose to propose modules that can run single raw files for the first steps of
pre-processing, since the initial actions in the extraction process treat files independantly.

Since the first steps can be run on each file, the use of **Dataset collection** is recommanded in Galaxy to avoid
launching jobs manually for each sample. You can consider the Dataset collection option from the very beginning, while
uploading your data into Galaxy.

> ### {% icon hands_on %} Hands-on: Data upload the mzXML with **Get data**
>
> 1. Create a new history for this tutorial
> 2. Import the 9 mzXML files from [Zenodo](http://doi.org/10.5281/zenodo.3244991) or a shared data library inside a collection
>    - HU_neg_048.mzML
>    - HU_neg_090.mzML
>    - HU_neg_123.mzML
>    - HU_neg_157.mzML
>    - HU_neg_173.mzML
>    - HU_neg_192.mzML
>    - QC1_002.mzML
>    - QC1_008.mzML
>    - QC1_014.mzML
>
>    ```
>    https://zenodo.org/record/3244991/files/HU_neg_048.mzML
>    https://zenodo.org/record/3244991/files/HU_neg_090.mzML
>    https://zenodo.org/record/3244991/files/HU_neg_123.mzML
>    https://zenodo.org/record/3244991/files/HU_neg_157.mzML
>    https://zenodo.org/record/3244991/files/HU_neg_173.mzML
>    https://zenodo.org/record/3244991/files/HU_neg_192.mzML
>    https://zenodo.org/record/3244991/files/QC1_002.mzML
>    https://zenodo.org/record/3244991/files/QC1_008.mzML
>    https://zenodo.org/record/3244991/files/QC1_014.mzML
>    ```
>
>    {% include snippets/import_via_link.md collection=true collection_type="mzml" collection_name="sacurine"%}
>    {% include snippets/import_from_data_library.md %}
>
{: .hands_on}

You should have in your history a green Dataset collection (`sacurine`) with 9 datasets with as format mzml.

Their size can be checked in their information panel (i)

## Data preparation for xcms

This first step is only meant to read your mzXML and generate an object usable by xcms.

**MSnbase readMSData** takes as input your raw files and prepares RData files for the first xcms step.

> ### {% icon hands_on %} Hands-on: MSnbase readMSData
>
> 1. **MSnbase readMSData** {% icon tool %} with the following parameters:
>    - *"File(s) from your history containing your chromatograms"*:
>        - Click on the folder icon to select the Dataset collection: `sacurine`
>
{: .hands_on}


> ### {% icon question %} Question
>
> What do you get as output?
>
> > ### {% icon solution %} Solution
> >
> > 1. A **Dataset collection** containing 9 dataset
> > 2. The dataset are some RData object with the datatype **rdata.msnbase.raw**
> >
> {: .solution}
>
{: .question}

## Importing a sample metadata file

A sample metadata file contains for each of your raw files their metadata:
- class which will be used during the preprocessing steps
- number of batch which will be useful for a batch correction step
- different experimental conditions which can be used for the statistics

Note that you can either:
- upload one already filled
- use a template (because it can be painful to get the sample list without misspelling or omission)
  1. generate a template with the tool **xcms get a sampleMetadata file**
  2. fill it using your favorite table editor (Excel, Libre Office)
  3. upload it within Galaxy

> ### {% icon comment %} Comment
>
> The file have to be a `.tsv` (tab-separated values). Neither `.xlsx` nor `.odt` are supported.
{: .comment}

> ### {% icon hands_on %} Hands-on: Data upload the sampleMetada with **Get data**
>
> 1. Import the sampleMetadata_completed.tsv file from [Zenodo](http://doi.org/10.5281/zenodo.3244991) or from a shared data library
>    - sampleMetadata_completed.tsv
>
>    ```
>    https://zenodo.org/record/3244991/files/sampleMetadata_completed.tsv
>    ```
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
{: .hands_on}

> ### {% icon question %} Question
>
> 1. How many columns should I have in my sampleMetadata file?
> 2. What kind of class can I have?
>
> > ### {% icon solution %} Solution
> >
> > 1. At least 2, with the identify and the class. But as many as you need to describe the potencial variablity of your samples (ex: the person in charge of the sample preparation, the temperature, ...). The statistic analysis will expose the relevant parameters.
> > 2. Sample, QC, blank... The class (the 2nd column) is useful for the preprocessing step with xmcs to detect the metabolite across the samples. So it's important to separate the samples and the QC. If you don't have any specific class, just fill everywhere with `sample` or a dot `.`
> >
> {: .solution}
>
{: .question}

## First xcms step: **peak picking**

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


# W4M 3-tables format

To do: add a complement about the 3-table format, depending on the information already given in xcms part.
Inparticular, a focus of mandatory information for the present tutorial.

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

Commonly, LC-MS analysis generates a significant number of variables (hundreds to thousands). Getting a complete view of
such dataset may not be an easy task, but getting a glimpse of it is possible using some common unsupervised multivariate
analysis. One of the most commonly used method is the Principal Components Analysis (PCA). You can get a basic PCA along with
over useful information using the Quality Metrics tool available in the Quality Control section.

> ### {% icon hands_on %} Hands-on: Using **Quality Metrics** to get an overview of your data
>
> Execute **Quality Metrics** {% icon tool %} with the following parameters:
>    - *"Data matrix file"*: `The one from 'xcms fillChromPeaks' outputs`
>    - *"Sample metadata file"*: `Your original completed sampleMetadata file`
>    - *"Variable metadata file"*: `The one from 'xcms fillChromPeaks' or 'CAMERA.annotate' outputs`
>
>
>    > ### {% icon comment %} Comment
>    >
>    > You can leave other parameters with default values. 
>    {: .comment}
>
{: .hands_on}

For a first overview of your data, you can focus on the graphical output of this tool: **Quality_Metrics_figure.pdf**.
It provides a variety of useful information:
 - summary of the intensities in the dataMatrix file
 - view of these intensities with a color scale
 - 2-components PCA score plot to check for clusters or outliers
 - sum of intensities per sample according to injection order to check the presence of signal drift or batch effect
 - z-scores for intensity distribution and proportion of missing values
 - ions' sd and mean values

![Quality_Metrics_figure.pdf](../../images/QM_9samp_raw.png)

> ### {% icon question %} Question time: cross-referencing information
>
> Look at the 'Quality_Metrics_figure.pdf' file.
> 1. Look at the proportion of ions with a pool coefficient of variation <30%. Knowing that pools are identical samples, and that
values with a CV > 30% can not be considered stables, what can you conclude about your dataset regarding the possibility to
compare intensities between samples?
> 2. Look at the sum of intensities per sample according to injection order. What major information do you observe on the plot?
Can this observation help you understand the CV results you just looked at?
> 3. Now look at the PCA plot. Can you explain the first component (t1)? Use the previous plot to help you.
>
> > ### {% icon solution %} Solution
> >
> > 1. You can read on the plot that *pool CV < 30%: 26%*, meaning that the pool values are stable only for a quarter of the ions
in your dataset. If the pooled samples are not stable, this means that you can observe differences between samples even when
there are no biological differences. Thus, comparing samples becomes difficult and has high risk of being unreliable.
Consequently, with only a quarter of the ions being stable regarding pool intensities, performing statistical analyses on
this full dataset would probably lead to unreliable results.
> > 2. We can see on the figure that the global intensity of samples seems to decrease with the injection order. In particular,
the fact that the pooled samples' intensities decrease leads us to suspect a signal drift due to the clogging effect of successive
injection of samples.
> > This signal drift could be the reason why so many ions in the dataset leaded to high CV values for pools, since it prevents
at least part of the ions to be stable regarding pools' intensities.
> > 3. If we look closely at the samples' identifiers on the plot, it seems that the lowest numbers in IDs are at the right side
of the first component. Knowing that these numbers correspond to an order in the injection sequence, we can link it to the
previous picture's samples. Then, what we can observe is that the order of samples in the first component of PCA from right to left
corresponds approximately to the decreasing order of sums of intensities. Thus, we can conclude that the main variability in
the dataset may be due to the signal drift.
> >
> {: .solution}
>
{: .question}

## Step 2: handling the signal drift observed althrough the analytical sequence

It is known that when injecting successively a large number of samples, the system tends to get dirty, and this may cause a measure drift.
To prevent inability to catch signal anymore, in case of large injection series, the sequence is generally divided into several batches
and the source is cleaned between batches. Unfortunately, these signal drift and batch design can add significant variability in data,
making sample comparison complicated. In case data is impacted by these effects, it is highly recommanded to normalise the data to get
rid of these unwanted effects.

In our case study, we saw that the data seemed to be affected by signal drift. Thus, we will use the **Batch_correction** module to
get rid of it.

> ### {% icon hands_on %} Hands-on: Data normalisation using the **Batch_correction** module
>
> Execute **Batch_correction** {% icon tool %} with the following parameters:
>    - *"Data matrix file"*: `The one from 'xcms fillChromPeaks' outputs`
>    - *"Sample metadata file"*: `Your original completed sampleMetadata file`
>    - *"Variable metadata file"*: `The one from 'xcms fillChromPeaks' or 'CAMERA.annotate' outputs`
>    - *"Type of regression model "*: `linear`
>        - *"Factor of interest "*: `gender`
>
> You can leave other parameters with default values. 
>
>    > ### {% icon comment %} Comment
>    >
>    > The choice of the type of regression model to use depends on several parameters.
>    > In this case-study, since we only have 3 pools, there are only two possible choices: *linear* or *all loess sample*
>    > When possible, we recommend to use pools to correct the signal drift, that is why we chose to run the module with *linear*.
>    {: .comment}
>
{: .hands_on}

**What transformation have this module done to the ions' intensities?**

For each ion independently, the normalisation process works as described in the folowing picture:

![How this works](../../images/BC_theo.png)

The methodology is meant to correct for signal drift. In the module, it is combined with a correction for batch effect. Thus, if your
sequence is divided into several batches, the idea is to obtain something like the folowing:

![Before/after picture](../../images/BC_theo2.png)

In the case of *linear* regression model, the module performs some tests before applying the normalisation for quality purposes.
For some ions, if the normalisation process would have led to unconsistant results, the concerned ions are not corrected for signal drift.
This kind of quality checks depends on the type of regression model you use. Please refer to the module's help section for more information!


## Step 3: getting rid of unreliable variables using CV

Now that the data is corrected for signal drift, we expect to have stable intensities within pools. But is this always the case?
Truth is, even when correcting ions, we may not manage to get rid of analytical effect for 100% of ions. And even if we could,
LC-MS data may contain noise signal, or ions that are not reliable as they are too noisy. Thus, it is possible that the data still
contains unusable ions.

To filter the ions not reliable enough, we can consider CVs as a filtering indicator. The **Quality Metrics** module provides
different CV indicators depending on what is in your sample list. In particular, in the present case-study, it can compute pool CVs
as previously seen, but also a ratio between pool CVs and sample CVs. This is particularly of interest since we can expect that,
whatever the pool CV value, it will be lower than the corresponding sample CV value, since biological samples are supposed to be affected
by biological variability. Thus, we can filter the ions that do not respect this particular condition.

> ### {% icon hands_on %} Hands-on: CV calculation using the **Quality Metrics** module
>
> Execute **Quality Metrics** {% icon tool %} with the following parameters:
>    - *"Data matrix file"*: `The one from Batch_correction outputs`
>    - *"Sample metadata file"*: `Your original completed sampleMetadata file`
>    - *"Variable metadata file"*: `The one from Batch_correction outputs`
>
>
>    > ### {% icon comment %} Comment
>    >
>    > Do not forget you need to use this module again, since this time indicators will be computed on normalised intensities.
>    > What we are going to use this time is the tabular output, but while you are at it you can always check the pdf file.
>    {: .comment}
>
{: .hands_on}

The module provides a variableMetadata tabular output, containing all the computed CV values. You can then use these values to filter
your data using the **Generic_Filter** module, available in the **COMMON TOOLS > Data Handling** section. 


> ### {% icon hands_on %} Hands-on: Data filtering using the **Generic_Filter** module
>
> Execute **Generic_Filter** {% icon tool %} with the following parameters:
>    - *"Data matrix file"*: `The one from Batch_correction outputs`
>    - *"Sample metadata file"*: `Your original completed sampleMetadata file or the one from Quality_Metrics outputs`
>    - *"Variable metadata file"*: `The one from Quality_Metrics outputs`
>    - *"Deleting samples and/or variables according to Numerical values"*: `yes`
>        - {% icon param-repeat %} *"Identify the parameter to filter "*
>            - *"On file"*: `Variable metadata`
>            - *"Name of the column to filter"*: `poolCV_over_sampleCV`
>            - *"Interval of values to remove"*: `upper`
>                - *"Remove all values upper than"*: `1.0`
>        - {% icon param-repeat %} *"Insert Identify the parameter to filter "*
>            - *"On file"*: `Variable metadata`
>            - *"Name of the column to filter"*: `pool_CV`
>            - *"Interval of values to remove"*: `upper`
>                - *"Remove all values upper than"*: `0.3`
>    - *"Deleting samples and/or variables according to Qualitative values"*: `yes`
>        - {% icon param-repeat %} *"Removing a level in factor"*
>            - *"Name of the column to filter"*: `sampleType`
>            - *"Remove factor when"*: `pool`
>
>
>    > ### {% icon comment %} Comment
>    >
>    > You can see here that you can take the opportunity of this filtering step to get rid of the pools. Indeed, next step is
statistical analysis, and you do not need the pools anymore since they do not participate in the scientific design of the study. 
>    {: .comment}
>
{: .hands_on}


> ### {% icon question %} Questions
>
> 1. What does the *1.0* threshold mean in the hands-on exercice you just executed?
> 2. How many variables are left in your dataset? How many samples?
>
> > ### {% icon solution %} Solution
> >
> > 1. The *1.0* value corresponds to the maximum value kept in the dataset ('Interval of values to remove: *upper*') regarding the
*poolCV_over_sampleCV* column in your *Variable metadata* file. This means that any ion with a pool CV / sample CV ratio above 1
(*i.e.* a pool CV greater than the sample CV) is discarded from the dataset. 
> > 2. Filtering led to 2706 ions and 6 samples.
> >
> {: .solution}
>
{: .question}




# Statistical analysis to find variables of interest

The question of data filtering and correction must be addressed in all projects, even thought in some cases it may lead to 
the decision of no action on data. Once you applied your customed processing procedure, your tables are ready for the 
statistical analysis. 

There is a large variety of statistical analysis methods that you can apply on metabolomic LC-MS data. The most standard 
strategy is a combination of univariate analysis (such as applying a Mann-Whitney-Wilcoxon test on each ion indepedently)
and multivariate analysis (such as constructing a PLS model using all your ions at once). What you should keep in mind is
that the choice of your statistical analysis strategy depends on both your data characteristics (such as colinearity or
dataset size) and your study design. You should think carefully about what is appropriate for your own project. 

In this tutorial, we will take the example of univariate analysis, using the `bmi` column of the **sampleMetadata file** as 
our variable of interest (body mass index). Since this variable is quantitative, we will chose in this example to mesure 
the link between the BMI and the measured ions using **statistical correlation calculation**. For more examples of 
statistical analysis performed on LC-MS data, you can take a few minutes to watch the usemetabo.org open course video
[here](https://usemetabo.org/courses/w4mlc-ms-statistical-analysis). 

## Computation of statistical indices

First thing is to compute the correlation coefficients used to estimate the link between the variable of interest `bmi` 
and the ions that we have in our dataset. For this calculation we can use the **Univariate** module in the 
**Statistical Analysis** section. 

> ### {% icon hands_on %} Hands-on: Statistical analysis using the **Univariate** module
>
> Execute **Univariate** {% icon tool %} with the following parameters:
>    - *"Data matrix file"*: `The one from Generic_filter outputs`
>    - *"Sample metadata file"*: `The one from Generic_filter outputs`
>    - *"Variable metadata file"*: `The one from Generic_filter outputs`
>    - *"Factor of interest"*: `bmi`
>    - *"Test"*: `Spearman correlation rank test (quantitative)`
>    - *"Method for multiple testing correction"*: `none`
>
>    > ### {% icon comment %} Comment
>    >
>    > In this tutorial, we chose to perform the analysis without multiple testing correction. This choice is not 
based on a relevant statistical strategy (which would more likely be to *use* multiple testing correction). It is based on
the fact that with only 6 biological samples in a dataset of 2706 ions it is almost impossible to settle for correlation 
coefficients significantly different from zero. Consequently, to illustrate better the filtering step that will follow,
we chose not to apply the multiple testing correction, allowing us to obtain 'significant' results regarding statistical
indices. 
>    {: .comment}
>
{: .hands_on}

The module provides different types of output. You will find statistical indices in the *variableMetadata output* (such as
p-values and statistical indicators). 'Significant' results are illustrated by graphics in the *Univariate_figure.pdf* file.  

> ### {% icon question %} Questions
>
> How many *significant* variables were found?
>
> > ### {% icon solution %} Solution
> >
> > The module found that 61 variables have a correlation coefficient significantly different from 0. 
> >
> {: .solution}
>
{: .question}


## Reducing the dataset to keep ions of interest only

In untargeted metabolomics, statistical analysis is usually used as a filter to focus on a subset of variables with potential. 
This subset can be used to proceed to identification and thus biological interpretation. Hence, statistical indices are
generally associated with thresholds allowing us to determine which ions should be kept or discarded. 

In our example of correlation analysis, two indices can be used to filter the data.
 - p-values: it indicates whether it is likely for a given correlation coefficient not to be actually different from zero; 
considering a threshold of 0.05 generally corresponds to a misleading risk of 5%. 
 - correlation coefficient: it indicates if the correlation between a given ion and the variable of interest is strong or not;
it goes from -1 to 1, with 0 meaning no correlation; in our example we consider as a sufficiently strong link a coefficient with
absolute value above 0.9. 

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
