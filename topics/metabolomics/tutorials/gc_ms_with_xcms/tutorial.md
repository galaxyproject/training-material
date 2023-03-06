---
layout: tutorial_hands_on

title: 'Mass spectrometry: GC-MS data processing (with XCMS, RAMClustR, RIAssigner, and matchms)'
zenodo_link: 'https://zenodo.org/record/7701862'
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
- The processing of untargeted GC-MS metabolomic data can be done using open-source tools.
- TBA
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

The study of metabolites in biological samples is routinely defined as metabolomics and provides the capability to investigate metabolism on a global and relatively unbiased scale in comparison to traditional targeted studies focused on specific pathways of metabolism and a small number of metabolites. The untargeted approach enables to detect thousands of metabolites in hypothesis-generating studies and to link previously unknown metabolites with biologically important roles. There are two major issues in contemporary metabolomics: the first is enormous loads of signal generated during the experiments, and the second is the fact that some metabolites in the studied samples may not be known to us. These obstacles make the task of processing and interpreting the metabolomics data a cumbersome and time-consuming process {% cite Nash2019 %}.

A lot of packages are available for the analysis of GC-MS or LC-MS data. In this tutorial, we focus on open-source solutions integrated within Galaxy framework. In this tutorial, we will learn how to process the data samples and identify the present compounds. For demonstration, we use three GC-[EI+]-HRMS files from seminal plasma samples.

To process the data, we use several tools. **XCMS** ({% cite Smith2006 %}) is a general package for untargeted metabolomics profiling. It can be used for any type of mass spectrometry acquisition files from low to high resolution, including FT-MS data coupled with different kind of chromatography (liquid or gas). We use it to detect peaks within our samples. Once we have detected them, they need to be deconvoluted with focus on consistency across samples. For that, we use **RAMClustR** ({% cite broeckling2014ramclust %}) tool. To normalise the retention time of identified spectra in our sample, we compute retention index using **RIAssigner** ({% cite hecht2022riassigner %}) by comparing the data to well-defined list of alkanes. Finally, we identify detected spectra by aligning them with a database of known compounds. This can be achieved using **MatchMS** ({% cite Huber2020 %}), resulting into a table of identified compounds, weighted by a confidence score.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data preparation and prepocessing

Before we can start with the actual analysis pipeline, we first need to download and prepare our dataset. Many of the preprocessing steps can be run in parallel on individual samples. Therefore, we recommend using the Dataset collections in Galaxy. This can be achieved by using the dataset collection option from the beginning of your analysis, when uploading your data into Galaxy.

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
>    https://zenodo.org/record/7701862/files/8_qc_no_dil_milliq.raw
>    https://zenodo.org/record/7701862/files/21_qc_no_dil_milliq.raw
>    https://zenodo.org/record/7701862/files/29_qc_no_dil_milliq.raw
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md collection=true format="mzml" collection_name="input" renaming=false %}
>
> 3. Make sure your data is in a **collection**. You can always manually create the collection from separate files:
>
>    {% snippet faqs/galaxy/collections_build_list.md %}
>
> 4. Import the following extra files from [Zenodo]({{ page.zenodo_link }}):
>
>    ```
>    https://zenodo.org/record/7701862/files/reference_alkanes.csv
>    https://zenodo.org/record/7701862/files/reference_spectral_library.msp
>    https://zenodo.org/record/7701862/files/sample_metadata.tsv
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
>    > <comment-title> The extra files </comment-title>
>    >
>    > The two additional files contain **reference alkanes**, reference spectral library, and sample metadata. The list of alkanes with retention time and carbon number or retention index is used to compute the retention index of the deconvoluted peaks. The alkanes should be measured ideally in the same batch as the input sample collection.
>    > 
>    > The **reference spectral library** (`.msp`) is used for identification of spectra. It contains the recorded and annotated mass spectra of compounds which can be detected in the sample and confirmed via comparison with this library. The specific library is the in-house library of metabolite standards. **TODO**: what does it mean?
>    > 
>    > The **sample metadata** corresponds to a table containing information about our samples. In particular, it contains sample name, type (QC, blank, sample, etc.), batch number, and injection order.
>    {: .comment}
>
{: .hands_on}

As a result of this step, we should have in our history a green Dataset collection with three `.raw` files.

## Convert data to mzML

>    ***TODO***: *rename inputs and outputs in tools??*

Our input data are in `.raw` format, which is not suitable for the downstream tools in this tutorial. We use tool {% tool [msconvert](toolshed.g2.bx.psu.edu/repos/galaxyp/msconvert/msconvert/3.0.20287.2) %} to convert them to the appropriate format (`.mzML` in this case).

> <hands-on-title> Convert data to mzML </hands-on-title>
>
> 1. {% tool [msconvert](toolshed.g2.bx.psu.edu/repos/galaxyp/msconvert/msconvert/3.0.20287.2) %} with the following parameters:
>    - {% icon param-collection %} *"Input unrefined MS data"*: `input` (Input dataset collection)
>    - *"Do you agree to the vendor licenses?"*: `Yes`
>    - *"Output Type"*: `mzML`
>    - In *"Data Processing Filters"*:
>        - *"Apply peak picking?"*: `Yes`
>
>    > <comment-title> Centroids </comment-title>
>    >
>    > msconvert with selected parameters computes centroids in the m/z domain. **TBD** explain why !!!
>    {: .comment}
>
{: .hands_on}

## Create XCMS object

The first part of data processing is using **XCMS** tool. This step is only meant to read our `.mzML` files and generate an object usable by **XCMS** tool. {% tool [MSnbase readMSData](toolshed.g2.bx.psu.edu/repos/lecorguille/msnbase_readmsdata/msnbase_readmsdata/2.16.1+galaxy0) %} takes as input our files and prepares `RData` files for the first **XCMS** step.

> <hands-on-title> Create **XCMS** object </hands-on-title>
>
> 1. {% tool [MSnbase readMSData](toolshed.g2.bx.psu.edu/repos/lecorguille/msnbase_readmsdata/msnbase_readmsdata/2.16.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"File(s) from your history containing your chromatograms"*: `input.mzML` (output of **msconvert** {% icon tool %})
>
>    {% snippet faqs/galaxy/tools_select_collection.md %}
>
>    > <comment-title> Output - `input.raw.RData` </comment-title>
>    >
>    > Collection of `rdata.msnbase.raw` files. **TODO**
>    {: .comment}
{: .hands_on}

# Peak detection using XCMS

The first step in the workflow is to detect the peaks in our data using **XCMS**. This part, however, is covered by a [separate tutorial]({{ site.baseurl }}/topics/metabolomics/tutorials/lcms-preprocessing/tutorial.html). Although the tutorial is dedicated to LC-MS data, it can be followed also for our GC data. Therefore, in this section, we do not explain this part of the workflow in detials, but rather refer the reader to the dedicated tutorial. Please also pay attention to the parameter values for individual Galaxy tools, as these can differ from the refered tutorial and are adjusted to our dataset.

## Peak picking

The first step is to extract peaks from each of your data files independently. For this purpose, we use _centWave_ chromatographic peak detection algorithm implemented in {% tool [xcms findChromPeaks (xcmsSet)](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_xcmsset/abims_xcms_xcmsSet/3.12.0+galaxy0) %}.

> <hands-on-title> Peak picking </hands-on-title>
>
>  {% tool [xcms findChromPeaks (xcmsSet)](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_xcmsset/abims_xcms_xcmsSet/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"RData file"*: `input.raw.RData` (output of **MSnbase readMSData** {% icon tool %})
>    - *"Extraction method for peaks detection"*: `CentWave - chromatographic peak detection using the centWave method`
>        - *"Max tolerated ppm m/z deviation in consecutive scans in ppm"*: `3.0`
>        - *"Min,Max peak width in seconds"*: `1,15`
>        - In *"Advanced Options"*:
>            - *"Prefilter step for for the first analysis step (ROI detection)"*: `3,500`
>            - *"Noise filter"*: `1000`
>  
>    You can leave the other parameters with their default values.
>
{: .hands_on}

## Determining shared ions

At this step, you obtain a dataset collection containing one `RData` file per sample, with independent lists of ions. Next, we want to identify the ions shared between samples. To do so, first you need to group your individual `RData` files into a single one.

> <hands-on-title> Merging files </hands-on-title>
>
>  {% tool [xcms findChromPeaks Merger](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_merge/xcms_merge/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"RData file"*: `input.raw.xset.RData` (output of **xcms findChromPeaks (xcmsSet)** {% icon tool %})
>    - {% icon param-file %} *"Sample metadata file "*: `sample_metadata.tsv` (Input dataset)
>
>    You can leave the other parameters with their default values.
>
{: .hands_on}

Now we can proceed with the grouping and determining shared ions among samples. The aim of this step, called _grouping_, is to obtain a single matrix of ions’ intensities for all samples.

> <hands-on-title> Grouping peaks </hands-on-title>
>
>  {% tool [xcms groupChromPeaks (group)](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_group/abims_xcms_group/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"RData file"*: `xset.merged.RData` (output of **xcms findChromPeaks Merger** {% icon tool %})
>    - *"Method to use for grouping"*: `PeakDensity - peak grouping based on time dimension peak densities`
>        - *"Bandwidth"*: `3.0`
>        - *"Minimum fraction of samples"*: `0.9`
>        - *"Width of overlapping m/z slices"*: `0.01`
>    - *"Get the Peak List"*: `Yes`
>
>    You can leave the other parameters with their default values.
>
{: .hands_on}

## Retention time correction

A deviation in retention time occurs from a sample to another, especially when you inject large sequences of samples. This steps aims at correcting retention time drift for each peak among samples.

> <hands-on-title> Retention time correction </hands-on-title>
>
>  {% tool [xcms adjustRtime (retcor)](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_retcor/abims_xcms_retcor/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"RData file"*: `xset.merged.groupChromPeaks.RData` (output of **xcms groupChromPeaks (group)** {% icon tool %})
>    - *"Method to use for retention time correction"*: `PeakGroups - retention time correction based on aligment of features (peak groups) present in most/all samples.`
>        - *"Minimum required fraction of samples in which peaks for the peak group were identified"*: `0.7`
>        - *"Smooth method"*: `loess - non-linear alignment`
>
>    You can leave the other parameters with their default values.
>
{: .hands_on}

## Second round of determining shared ions

By applying retention time correction, the used retention time values were modified. Consequently, applying this step on your data requires to complete it with an additional _grouping_ step.

> <hands-on-title> Grouping peaks </hands-on-title>
>
>  {% tool [xcms groupChromPeaks (group)](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_group/abims_xcms_group/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"RData file"*: `xset.merged.groupChromPeaks.adjustRtime.RData` (output of **xcms adjustRtime (retcor)** {% icon tool %})
>    - *"Method to use for grouping"*: `PeakDensity - peak grouping based on time dimension peak densities`
>        - *"Bandwidth"*: `3.0`
>        - *"Minimum fraction of samples"*: `0.9`
>        - *"Width of overlapping m/z slices"*: `0.01`
>    - *"Get the Peak List"*: `Yes`
>        - *"Convert retention time (seconds) into minutes"*: `Yes`
>        - *"If NA values remain, replace them by 0 in the dataMatrix"*: `Yes`
>
>    You can leave the other parameters with their default values.
>
{: .hands_on}

## Integrating areas of missing peaks

At this point, the peak list may contain `NA` values when peaks were not considered peaks in only some of the samples in the first peak picking step. In this step, we will integrate signal in the mz-rt area of an ion (chromatographic peak group) for samples in which no chromatographic peak for this ion was identified.

> <hands-on-title> Integrating areas of missing peaks </hands-on-title>
>
>  {% tool [xcms fillChromPeaks (fillPeaks)](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_fillpeaks/abims_xcms_fillPeaks/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"RData file"*: `xset.merged.groupChromPeaks.adjustRtime.groupChromPeaks.RData` (output of **xcms groupChromPeaks (group)** {% icon tool %})
>    - In *"Peak List"*:
>        - *"Convert retention time (seconds) into minutes"*: `Yes`
>
>    You can leave the other parameters with their default values.
>
>    > <comment-title> Output </comment-title>
>    >
>    > TBA 
>    {: .comment}
>
{: .hands_on}

# Peak deconvolution

The next step is deconvoluting the detected peaks in order to reconstruct the full spectra of the analysed compound. {% tool [RAMClustR](toolshed.g2.bx.psu.edu/repos/recetox/ramclustr/ramclustr/1.2.4+galaxy2) %} is used to group features based on correlations across samples in a hierarchy, focusing on consistency across samples. While a feature typically is derived from a single compound, a spectrum of mass signals is more a more-accurate representation of the mass spectrometric signal for a given metabolite. **RAMClustR** uses a novel grouping method that operates in an unsupervised manner to group signals from MS data into spectra without relying on predictability of the in-source phenomenon.

**TODO** description based on https://pubs.acs.org/doi/10.1021/ac4019268

> <hands-on-title> Peak deconvolution </hands-on-title>
>
> 1. {% tool [RAMClustR](toolshed.g2.bx.psu.edu/repos/recetox/ramclustr/ramclustr/1.2.4+galaxy2) %} with the following parameters:
>    - *"Choose input format:"*: `XCMS`
>        - In *"Input MS Data as XCMS"*:
>            - {% icon param-file %} *"Input XCMS"*: `xset.merged.groupChromPeaks.adjustRtime.groupChromPeaks.fillChromPeaks.RData` (output of **xcms fillChromPeaks (fillPeaks)** {% icon tool %})
>        - In *"General parameters"*:
>            - *"Sigma r"*: `0.7`
>            - *"Maximum RT difference"*: `10.0`
>    - In *"Clustering"*:
>        - *"Minimal cluster size"*: `5`
>        - *"Maximal tree height"*: `0.9`
>
{: .hands_on}

The spectral data comes as `.msp` file, which is a text file structured according to the **NIST Search** spectra format. `.msp` is one of the generally accepted formats for mass spectra representations and it is compatible with lots of spectra processing programms (MS-DIAL, NIST MS Search, AMDIS, etc.). Because `.msp` files are text-based, they can be viewed as simple `txt` files. You can use any text editor that you have on your computer or use Galaxy built-in editor. In this tutorial we use the Galaxy editor to check the contents of the file:

> <hands-on-title> Data Exploration </hands-on-title>
>
> 1. Click *"Visualize this data"* {% icon galaxy-barchart %} icon next to the dataset in the Galaxy history
> 2. Select the *"Editor"* {% icon galaxy-eye %} tool. The contents of the file would look like this:
>
>    ```
>     NAME:C001
>     IONMODE:Negative
>     SPECTRUMTYPE:Centroid
>     RETENTIONTIME:383.27
>     Num Peaks:231
>     217.1073 64041926
>     243.0865 35597866
>     257.1134 31831229
>     224.061 27258239
>     258.11 24996353
>     241.0821 23957171
>     315.1188 13756744
>     ...
>
>     NAME:C002
>     IONMODE:Negative
>     SPECTRUMTYPE:Centroid
>     RETENTIONTIME:281.62
>     Num Peaks:165
>     307.1573 299174880
>     147.0654 298860831
>     149.0447 287809889
>     218.1066 118274758
>     189.076 112486871
>     364.1787 75134143
>     191.0916 52526567
>     308.1579 52057158
>     ...
>    ```
>
{: .hands_on}

`.msp` files can contain one or more mass spectra, these are split by an empty line. The individual spectra essentialy consist of two sections: metadata and peaks. The metadata consists of compound name, spectrum type (which is centroid in this case), ion mode, retention time, and number of *m/z* peaks. If the compound has not been identified as in our case, the **NAME** can be any arbitrary string. It is best however, if that string is unique within that `.msp` file. The metadata fields are usually unordered, so it is quite common for one `.msp` file to contain **NAME** as the first metadata key, and for another `.msp` to have this key somewhere in the middle. The keys themselves also aren't rigid and can have different names (e.g., **compound_name** instead of **NAME**) or store information not present in this `.msp`, such as ionization mode.

However, as we can observe, the metadata part is rather incomplete. We would like to gather more information about the detected spectra and identify the specific compounds corresponding to them.

**TODO** describe `Spec Abundance` file

# Retention index calculation

Retention index is a way how to convert equipment- and experiment-specific retention times into system-independent normalised constants. The retention index of a compound is computed from the retention time by interpolating between adjacent alkanes. This can be different for individual chromatographic system, but the derived retention indices are quite independent and allow comparing values measured by different analytical laboratories.

We use package {% tool [RIAssigner](toolshed.g2.bx.psu.edu/repos/recetox/riassigner/riassigner/0.3.2+galaxy1) %} to compute retention indices for files in the `.msp` format using an indexed reference list of alkanes in `.csv` or `.msp` format. The output follows the same format as the input, but with added retention index values. These can be used at a subsequent stage to improve compound identification. Multiple computation methods (e.g. piecewise-linear or cubic spline) are supported by the tool.

> <hands-on-title> Retention index calculation </hands-on-title>
>
> 1. {% tool [RIAssigner](toolshed.g2.bx.psu.edu/repos/recetox/riassigner/riassigner/0.3.2+galaxy1) %} with the following parameters:
>    - In *"Query dataset"*:
>        - {% icon param-file %} *"Query compound list"*: `Mass spectra from RAMClustR` (output of **RAMClustR** {% icon tool %})
>    - In *"Reference dataset"*:
>        - {% icon param-file %} *"Reference compound list"*: `reference_alkanes.csv` (Input dataset)
>
>
>    > <comment-title> Minutes vs. seconds </comment-title>
>    >
>    > You might notice that in the last **XCMS** step, we converted the retention times to minutes. However, in this step, we are using seconds again. The reason is that **RAMClustR** converted them internally again to seconds. Regarding the reference compound list, this database has already its retention times in seconds.
>    {: .comment}
>
>    > <comment-title> Kovats method </comment-title>
>    >
>    > TBD
>    {: .comment}
>
{: .hands_on}

# Identification

To identify and annotate the deconvoluted spectra, we compare them with a reference spectral library. This library contains spectra of standards measured on the same instrument for optimal comparability. The **matchMS** package is used for spectral matching. It was build to import and apply different similarity measures to compare large amounts of spectra. This includes common cosine scores, but can also easily be extended by custom measures.

## Compute similairity scores

> **TODO** how matchms is doing the matching/computing scores (normalised dot product) - peaks must be present in both reference and measured data + describe how it is diff from forward and reverse matching

We use the cosine score with a greedy peak pairing heuristic to compute the number of matching ions with a given tolerance and the cosine scores for the matched peaks.

> <hands-on-title> Compute similairity scores </hands-on-title>
>
> 1. {% tool [matchMS similarity](toolshed.g2.bx.psu.edu/repos/recetox/matchms/matchms/0.17.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Queries spectra"*: `RI using kovats of Mass spectra from RAMClustR` (output of **RIAssigner** {% icon tool %})
>    - *"Symmetric"*: `No` (if we were to query our spectra agains itself, we would select `Yes`)
>        - {% icon param-file %} *"Reference spectra"*: `reference_spectral_library.msp` (downloaded file from Zenodo)
>    - In *"Algorithm Parameters"*:
>        - *"tolerance"*: `0.03`
>    - *"Apply RI filtering"*: `Yes`
>
>    > <comment-title> Reference spectra </comment-title>
>    >
>    > More details about the reference spectra?
>    {: .comment}
>
{: .hands_on}

> <details-title> Overview of the spectral similarity scores </details-title>
> > <h5>Cosine Greedy</h5>
> > The cosine score, also known as the dot product, is based on representing the similarity of two spectra through the cosine of an  angle between the vectors that the spectra produce. Two peaks are considered as matching if their *m/z* values lie within the given   tolerance. Cosine greedy looks up matching peaks in a "greedy" way, which does not always lead to the most optimal alignments.
> >
> > This score was among the first to be used for looking up matching spectra in spectral libraries and to this day remains one of the most popular scoring methods for both library matching and molecular networking workflows.
> > <p>&nbsp;</p>
>
> > <h5>Cosine Hungarian</h5>
> > This method computes the similarities in the same way as the *Cosine Greedy* but with a difference in *m/z* peak alignment. The difference lies in that the Hungarian algorithm is used here to find matching peaks. This leads to the best peak pairs match, but can take significantly longer than the "greedy" algorithm.
> > <p>&nbsp;</p>
>
> > <h5>Modified Cosine</h5>
> > Modified Cosine is another, as its name states, representative of the family of cosine-based scores. This method aligns peaks by finding the best possible matches and consideres two peaks a match if their *m/z* values are within a tolerance before or after a mass-shift is applied. A mass shift is essentially a difference of precursor-*m/z* of two compared spectra. The similarity is then again expressed as a cosine of the angle between two vectors.
> > <p>&nbsp;</p>
>
> > <h5>Neutral Losses Cosine</h5>
> > Neutral Loss metric works similar to all described above with one major difference: instead of encoding the spectra as "intensity vs *m/z*" vector it encodes it to a "intensity vs *Δm/z*", where delta is computed as an *m/z* difference between precursor and a fragment *m/z*. This in theory, could better capture the underlying structural similarities between molecules.
> >
>
{: .details}

## Format the output

The output of the previous step is a `json` file. This format is very simple to read for our computers, but not so much for us. We can use {% tool [matchms output formatter](toolshed.g2.bx.psu.edu/repos/recetox/matchms_formatter/matchms_formatter/0.1.4) %} to convert the data to a tab-separated file with a scores matrix.

The output table contains the scores and number of matched ions of the deconvoluted spectra with the spectra in the reference library. The raw output is filtered to only contain the top matches (3 by default) and is then further filtered to contain only pairs with a score and number of matched ions larger than provided thresholds.

> <hands-on-title> Format the output </hands-on-title>
>
> 1. {% tool [matchms output formatter](toolshed.g2.bx.psu.edu/repos/recetox/matchms_formatter/matchms_formatter/0.1.4) %} with the following parameters:
>    - {% icon param-file %} *"Scores object"*: `CosineGreedy scores` (output of **matchMS similarity** {% icon tool %})
>
{: .hands_on}

> | query | reference | matches | scores |
> |-------|-----------|---------|--------|
> | C001  | Uridine_4TMS isomer 1   | 81       | 0.787    |
> | C004  | Asparagine_3TMS   | 56       | 0.909   |
> | C012  | Myo-inositol_6TMS | 29       | 0.688   |
> | ...   | ...       | ...     | ...    |
{: .matrix}

At this stage, all steps are complete: we have the list of identified spectra corresponding to a compound from the reference database. Since this process is dependent on the parameters we used along the way, the result should be considered dependent on them. To make enumerate this dependency, we each potential compound has assigned a confidence score.

# Conclusion

***TODO***