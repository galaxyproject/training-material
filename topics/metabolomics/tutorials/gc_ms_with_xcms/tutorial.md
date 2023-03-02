---
layout: tutorial_hands_on

title: 'Mass spectrometry: GC-MS analysis with XCMS, RAMClustR, and matchMS'
zenodo_link: 'https://zenodo.org/record/6878356'
questions:
- What are the main steps of GC-MS data processing for metabolomic analysis?
- What similarity metrics can be used to compare a pair of mass spectra and what are the differences between them?
objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentences describing what a learner should be able to do once they
  have completed the tutorial
- You can use Bloom's Taxonomy to write effective learning objectives
time_estimation: 2H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
requirements :
  - type: "internal"
    topic_name: metabolomics
    tutorials: 
      - lcms-preprocessing
contributors:
- xtrojak
- hechth

---


# Introduction

The study of metabolites in biological samples is routinely defined as metabolomics and provides the capability to investigate metabolism on a global and relatively unbiased scale in comparison to traditional targeted studies focused on specific pathways of metabolism and a small number of metabolites. The untargeted approach enables to detect thousands of metabolites in hypothesis-generating studies and to link previously unknown metabolites with biologically important roles. There are two major issues in contemporary metabolomics: the first is enormous loads of signal generated during the experiments, and the second is the fact that some metabolites in the studied samples may not be known to us. These obstacles make the task of processing and interpreting the metabolomics data a cumbersome and time-consuming process 

TBD

In this tutorial you will learn ......

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data preparation and prepocessing

## Import the data into Galaxy

> <hands-on-title> Upload data </hands-on-title>
>
> 1. Create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) into a collection:
>
>    ```
>    https://zenodo.org/record/6878356/files/8_qc_no_dil_milliq.raw
>    https://zenodo.org/record/6878356/files/21_qc_no_dil_milliq.raw
>    https://zenodo.org/record/6878356/files/29_qc_no_dil_milliq.raw
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md collection=true format="mzml" collection_name="input" renaming=false %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md astype="as a Collection" %}
>
> 3. Make sure your data is in a **collection**. You can always manually create the collection from separate files:
>
>    {% snippet faqs/galaxy/collections_build_list.md %}
>
> 4. Import the following extra files from Zenodo: 
>
>    ```
>    TBA
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
>    > <comment-title> The extra files </comment-title>
>    >
>    > The two additional files contain reference alkanes, reference spectral library, and sample metadata. The list of alkanes with retention time and carbon number or retention index is used to compute the retention index of the deconvoluted peaks. The alkanes should be measured ideally in the same batch as the input sample collection. A reference spectral library (`.msp`) is used for identification of spectra. The spectral library contains the recorded and annotated mass spectra of compounds which can be detected in the sample and confirmed via comparison with this library. The specific library is the in-house library of metabolite standards. Sample metadata corresponds to a table containing information about our samples. In particular, it contains sample name, type (QC, blank, sample, etc.), batch number, and injection order.
>    {: .comment}
>
{: .hands_on}

## Convert data to mzML

Our input data are in `.raw` format, which is not suitable for the downstream tools in this tutorial. We can use tool {% tool [msconvert](toolshed.g2.bx.psu.edu/repos/galaxyp/msconvert/msconvert/3.0.20287.2) %} to convert them to the appropriate format (`.mzML` in this case).

> <hands-on-title> Convert data to mzML </hands-on-title>
>
> 1. {% tool [msconvert](toolshed.g2.bx.psu.edu/repos/galaxyp/msconvert/msconvert/3.0.20287.2) %} with the following parameters:
>    - {% icon param-collection %} *"Input unrefined MS data"*: `output` (Input dataset collection)
>    - *"Do you agree to the vendor licenses?"*: `Yes`
>    - *"Output Type"*: `mzML`
>    - In *"Data Processing Filters"*:
>        - *"Apply peak picking?"*: `Yes`
>        - *"Apply m/z refinement with identification data?"*: `Yes`
>        - *"(Re-)calculate charge states?"*: `no`
>        - *"Filter m/z Window"*: `Yes`
>        - *"Filter out ETD precursor peaks?"*: `Yes`
>        - *"De-noise MS2 with moving window filter"*: `Yes`
>        - *"Demultiplex overlapping or MSX spectra"*: `Yes`
>    - In *"Scan Inclusion/Exclusion Filters"*:
>        - *"Filter MS Levels"*: `Yes`
>    - In *"General Options"*:
>        - *"Sum adjacent scans"*: `Yes`
>        - *"Output multiple runs per file"*: `Yes`
>
>    > <comment-title> Centroids </comment-title>
>    >
>    > msconvert with selected parameters computes centroids in the m/z domain. **TBD** exmplain why !!!
>    {: .comment}
>
{: .hands_on}

## Create XCMS object

The first part of data processing is using XCMS tool. This step is only meant to read our `.mzML` files and generate an object usable by XCMS tool. {% tool [MSnbase readMSData](toolshed.g2.bx.psu.edu/repos/lecorguille/msnbase_readmsdata/msnbase_readmsdata/2.16.1+galaxy0) %} takes as input our files and prepares `RData` files for the first XCMS step.

> <hands-on-title> Create XCMS object </hands-on-title>
>
> 1. {% tool [MSnbase readMSData](toolshed.g2.bx.psu.edu/repos/lecorguille/msnbase_readmsdata/msnbase_readmsdata/2.16.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"File(s) from your history containing your chromatograms"*: `output` (output of **msconvert** {% icon tool %})
>
>    {% snippet faqs/galaxy/tools_select_collection.md %}
{: .hands_on}

# Peak detection using XCMS

The first step in the workflow is to detect the peaks in our data using XCMS. This part, however, is covered by a [separate tutorial]({{ site.baseurl }}/topics/metabolomics/tutorials/lcms-preprocessing/tutorial.html). Although the tutorial is dedicated to LC-MS data, it can be followed also for our GC data. Therefore, in this section, we do not explain this part of the workflow in detials, but rather refer the reader to the dedicated tutorial. Please also pay attention to the parameter values for individual Galaxy tools, as these can differ from the refered tutorial and are adjusted to our dataset.

## Peak picking

The first step is to extract peaks from each of your data files independently. For this purpose, we use _centWave_ chromatographic peak detection algorithm implemented in {% tool [xcms findChromPeaks (xcmsSet)](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_xcmsset/abims_xcms_xcmsSet/3.12.0+galaxy0) %}.

> <hands-on-title> Peak picking </hands-on-title>
>
>  {% tool [xcms findChromPeaks (xcmsSet)](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_xcmsset/abims_xcms_xcmsSet/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"RData file"*: `xsetRData` (output of **MSnbase readMSData** {% icon tool %})
>    - *"Extraction method for peaks detection"*: `CentWave - chromatographic peak detection using the centWave method`
>        - *"Max tolerated ppm m/z deviation in consecutive scans in ppm"*: `3.0`
>        - *"Min,Max peak width in seconds"*: `1,15`
>        - In *"Advanced Options"*:
>            - *"Prefilter step for for the first analysis step (ROI detection)"*: `3,500`
>            - *"Noise filter"*: `1000`
>    - In *"Resubmit your raw dataset or your zip file"*:
>        - *"Resubmit your dataset or your zip file"*: `no need`
>  
>    You can leave the other parameters with their default values.
>
{: .hands_on}

## Determining shared ions

At this step, you obtain a dataset collection containing one `RData` file per sample, with independent lists of ions. Next, we want to identify the ions shared between samples. To do so, first you need to group your individual `RData` files into a single one.

> <hands-on-title> Merging files </hands-on-title>
>
>  {% tool [xcms findChromPeaks Merger](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_merge/xcms_merge/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"RData file"*: `xsetRData` (output of **xcms findChromPeaks (xcmsSet)** {% icon tool %})
>    - {% icon param-file %} *"Sample metadata file "*: `output` (Input dataset)
>    - In *"Resubmit your raw dataset or your zip file"*:
>        - *"Resubmit your dataset or your zip file"*: `no need`
>
>    You can leave the other parameters with their default values.
>
{: .hands_on}

Now we can proceed with the grouping and determining shared ions among samples. The aim of this step, called _grouping_, is to obtain a single matrix of ions’ intensities for all samples.

> <hands-on-title> Grouping peaks </hands-on-title>
>
>  {% tool [xcms groupChromPeaks (group)](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_group/abims_xcms_group/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"RData file"*: `xsetRData` (output of **xcms findChromPeaks Merger** {% icon tool %})
>    - *"Method to use for grouping"*: `PeakDensity - peak grouping based on time dimension peak densities`
>        - *"Bandwidth"*: `3.0`
>        - *"Minimum fraction of samples"*: `0.9`
>        - *"Width of overlapping m/z slices"*: `0.01`
>    - *"Get the Peak List"*: `Yes`
>    - In *"Resubmit your raw dataset or your zip file"*:
>        - *"Resubmit your dataset or your zip file"*: `no need`
>
>    You can leave the other parameters with their default values.
>
{: .hands_on}

## Retention time correction

A deviation in retention time occurs from a sample to another, especially when you inject large sequences of samples. This steps aims at correcting retention time drift for each peak among samples.

> <hands-on-title> Retention time correction </hands-on-title>
>
>  {% tool [xcms adjustRtime (retcor)](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_retcor/abims_xcms_retcor/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"RData file"*: `xsetRData` (output of **xcms groupChromPeaks (group)** {% icon tool %})
>    - *"Method to use for retention time correction"*: `PeakGroups - retention time correction based on aligment of features (peak groups) present in most/all samples.`
>        - *"Minimum required fraction of samples in which peaks for the peak group were identified"*: `0.7`
>        - *"Smooth method"*: `loess - non-linear alignment`
>    - In *"Resubmit your raw dataset or your zip file"*:
>        - *"Resubmit your dataset or your zip file"*: `no need`
>
>    You can leave the other parameters with their default values.
>
{: .hands_on}

## Second round of determining shared ions

By applying retention time correction, the used retention time values were modified. Consequently, applying this step on your data requires to complete it with an additional _grouping_ step.

> <hands-on-title> Grouping peaks </hands-on-title>
>
>  {% tool [xcms groupChromPeaks (group)](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_group/abims_xcms_group/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"RData file"*: `xsetRData` (output of **xcms adjustRtime (retcor)** {% icon tool %})
>    - *"Method to use for grouping"*: `PeakDensity - peak grouping based on time dimension peak densities`
>        - *"Bandwidth"*: `3.0`
>        - *"Minimum fraction of samples"*: `0.9`
>        - *"Width of overlapping m/z slices"*: `0.01`
>    - *"Get the Peak List"*: `Yes`
>        - *"Convert retention time (seconds) into minutes"*: `Yes`
>        - *"If NA values remain, replace them by 0 in the dataMatrix"*: `Yes`
>    - In *"Resubmit your raw dataset or your zip file"*:
>        - *"Resubmit your dataset or your zip file"*: `no need`
>
>    You can leave the other parameters with their default values.
>
{: .hands_on}

## Integrating areas of missing peaks

At this point, the peak list may contain `NA` values when peaks were not considered peaks in only some of the samples in the first peak picking step. In this step, we will integrate signal in the mz-rt area of an ion (chromatographic peak group) for samples in which no chromatographic peak for this ion was identified.

> <hands-on-title> Integrating areas of missing peaks </hands-on-title>
>
>  {% tool [xcms fillChromPeaks (fillPeaks)](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_fillpeaks/abims_xcms_fillPeaks/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"RData file"*: `xsetRData` (output of **xcms groupChromPeaks (group)** {% icon tool %})
>    - In *"Peak List"*:
>        - *"Convert retention time (seconds) into minutes"*: `Yes`
>    - In *"Resubmit your raw dataset or your zip file"*:
>        - *"Resubmit your dataset or your zip file"*: `no need`
>
>    You can leave the other parameters with their default values.
>
{: .hands_on}

# Peak deconvolution

The next step is deconvoluting the detected peaks in order to reconstruct the full spectra of the analysed compound. RAMClustR is used to group features based on correlations across samples in a hierarchy, focusing on consistency across samples.

## Sub-step with **RAMClustR**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [RAMClustR](toolshed.g2.bx.psu.edu/repos/recetox/ramclustr/ramclustr/1.2.4+galaxy2) %} with the following parameters:
>    - *"Choose input format:"*: `XCMS`
>        - In *"Input MS Data as XCMS"*:
>            - {% icon param-file %} *"Input XCMS"*: `xsetRData` (output of **xcms fillChromPeaks (fillPeaks)** {% icon tool %})
>        - In *"General parameters"*:
>            - *"Sigma r"*: `0.7`
>            - *"Maximum RT difference"*: `10.0`
>    - In *"Clustering"*:
>        - *"Minimal cluster size"*: `5`
>        - *"Maximal tree height"*: `0.9`
>    - In *"Normalisation"*:
>        - *"Normalisation method"*: `none`
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

# Retention index calculation

We developed a new package RIAssigner to compute retention indices for files in the `.msp` format using an indexed reference list in `.csv` or `.msp` format.

The output follows the same format as the input, but with added retention index values. These will be used at a later stage to improve compound identification with an additional filtering step. Multiple computation methods (piecewise-linear & cubic spline) are supported.

## Sub-step with **RIAssigner**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [RIAssigner](toolshed.g2.bx.psu.edu/repos/recetox/riassigner/riassigner/0.3.2+galaxy1) %} with the following parameters:
>    - In *"Query dataset"*:
>        - {% icon param-file %} *"Query compound list"*: `mass_spectra_merged` (output of **RAMClustR** {% icon tool %})
>    - In *"Reference dataset"*:
>        - {% icon param-file %} *"Reference compound list"*: `output` (Input dataset)
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

# Identification

The deconvoluted spectra are annotated for identification by comparing them with a reference spectral library. This library contains spectra of standards measured on the same instrument for optimal comparability. The matchms package is used for spectral matching. The cosine score with a greedy peak pairing heuristic was used to compute the number of matching ions with a given tolerance and the cosine scores for the matched peaks.

## Sub-step with **matchMS similarity**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [matchMS similarity](toolshed.g2.bx.psu.edu/repos/recetox/matchms/matchms/0.17.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Queries spectra"*: `output` (output of **RIAssigner** {% icon tool %})
>    - *"Symmetric"*: `Yes`
>    - In *"Algorithm Parameters"*:
>        - *"tolerance"*: `0.03`
>    - *"Apply RI filtering"*: `Yes`
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

## Sub-step with **matchms output formatter**

The output table contains the scores and number of matched ions of the deconvoluted spectra with the spectra in the reference library. The raw output is filtered to only contain the top matches (3 by default) and is then further filtered to contain only pairs with a score and number of matched ions larger than provided thresholds (0.65 & 3 by default).

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [matchms output formatter](toolshed.g2.bx.psu.edu/repos/recetox/matchms_formatter/matchms_formatter/0.1.4) %} with the following parameters:
>    - {% icon param-file %} *"Scores object"*: `similarity_scores` (output of **matchMS similarity** {% icon tool %})
>    - *"Formatting method"*: `Thresholding`
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