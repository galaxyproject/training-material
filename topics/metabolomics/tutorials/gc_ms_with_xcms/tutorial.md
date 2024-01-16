---
layout: tutorial_hands_on

title: 'Mass spectrometry: GC-MS data processing (with XCMS, RAMClustR, RIAssigner, and matchms)'
zenodo_link: 'https://zenodo.org/record/7890956'
level: Intermediate
questions:
- What are the main steps of gas chromatography-mass spectrometry (GC-MS) data processing for metabolomic analysis?
- What similarity metrics can be used to compare a pair of mass spectra and what are the differences between them?
- Do you know any alternative tools that can be used in place of the individual steps of this workflow?
objectives:
- To learn about the key steps in the preprocessing and analysis of untargeted GC-MS metabolomics data.
- To explore what open-source alternative tools can be used in the analysis of GC-MS data, learn about their possible parametrisations.
- To analyse authentic data samples and compare them with a data library of human metabolome, composed from a collection of mostly endogenous compounds.
time_estimation: 2H
key_points:
- The processing of untargeted GC-MS metabolomic data can be done using open-source tools.
- This processing workflow is complementary to CAMERA and metaMS.
requirements :
  - type: "internal"
    topic_name: metabolomics
    tutorials: 
      - lcms-preprocessing
contributions:
  authorship:
    - xtrojak
    - hechth
    - maximskorik
  editing:
    - hexylena

---

The study of metabolites in biological samples is routinely defined as metabolomics. Metabolomics' studies based on untargeted mass spectrometry provide the capability to investigate metabolism on a global and relatively unbiased scale in comparison to traditional targeted studies focused on specific pathways of metabolism and a small number of metabolites. The untargeted approach enables the detection of thousands of metabolites in hypothesis-generating studies and links previously unknown metabolites with biologically important roles {% cite Patti2012 %}. There are two major issues in contemporary mass spectrometry-based metabolomics: the first is enormous loads of signal generated during the experiments, and the second is the fact that some metabolites in the studied samples may not be known to us. These obstacles make the task of processing and interpreting the metabolomics data a cumbersome and time-consuming process {% cite Nash2019 %}.

Many packages are available for the analysis of GC-MS or LC-MS data - for more details see the reviews by {% cite Stanstrup2019 %} and {% cite Misra2021 %}. In this tutorial, we focus on open-source solutions integrated within the Galaxy framework, namely **XCMS**, **RAMClustR**, **RIAssigner**, and **matchms**. In this tutorial, we will learn how to (1) extract features from the raw data using **XCMS** ({% cite Smith2006 %}), (2) deconvolute the detected features into spectra with **RAMClustR** ({% cite broeckling2014ramclust %}), (3) compute retention indices with **RIAssigner** ({% cite hecht2022riassigner %}), and (4) identify the present compounds leveraging spectra and retention indices using **matchms** ({% cite Huber2020 %}). For demonstration, we use three GC-[EI+] high-resolution mass spectrometry data files generated from quality control seminal plasma samples.


> <details-title> Seminal plasma samples </details-title>
> 
> The seminal plasma samples were analyzed according to the standard operating procedure [(SOP) for metabolite profiling of seminal plasma via GC Orbitrap](https://zenodo.org/record/5734331).
> The 3 samples used in this training are pooled quality control (QC) samples coming from about 200 samples. The pooled samples were analyzed in dilution series to test the system suitability and the quality of the assay.
>
{: .details}

To process the data, we use several tools. **XCMS** ({% cite Smith2006 %}) is a general package for untargeted metabolomics profiling. It can be used for any type of mass spectrometry acquisition (centroid and profile) or resolution (from low to high resolution), including FT-MS data coupled with a different kind of chromatography (liquid or gas). We use it to detect chromatographic peaks within our samples. Once we have detected them, they need to be deconvoluted into spectra representing chemical compounds. For that, we use **RAMClustR** ({% cite broeckling2014ramclust %}) tool. To normalise the retention time of deconvoluted spectra in our sample, we compute the retention index using **RIAssigner** ({% cite hecht2022riassigner %}) by comparing the data to a well-defined list of reference compound (commonly alkanes) analyzed on the same GC column. Finally, we identify detected spectra by aligning them with a database of known compounds. This can be achieved using **matchms** ({% cite Huber2020 %}), resulting in a table of identified compounds weighted by a matching score.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data preparation and prepocessing

Before we can start with the actual analysis pipeline, we first need to download and prepare our dataset. Many of the preprocessing steps can be run in parallel on individual samples. Therefore, we recommend using the Dataset collections in Galaxy. This can be achieved by using the dataset collection option from the beginning of your analysis when uploading your data into Galaxy.

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
>    https://zenodo.org/record/7890956/files/8_qc_no_dil_milliq.raw
>    https://zenodo.org/record/7890956/files/21_qc_no_dil_milliq.raw
>    https://zenodo.org/record/7890956/files/29_qc_no_dil_milliq.raw
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md collection=true format="mzml" collection_name="input" renaming=false %}
>
> 3. Make sure your data is in a **collection**. You can always manually create the collection from separate files:
>
>    {% snippet faqs/galaxy/collections_build_list.md %}
>
>    In the further steps, this dataset collection will be referred to as `input` (and we recommend naming this collection like that to avoid confusion).
>
> 4. Import the following extra files from [Zenodo]({{ page.zenodo_link }}):
>
>    ```
>    https://zenodo.org/record/7890956/files/reference_alkanes.csv
>    https://zenodo.org/record/7890956/files/reference_spectral_library.msp
>    https://zenodo.org/record/7890956/files/sample_metadata.tsv
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
> 
>    Please pay attention to the format of all uploaded files, and make sure they were correctly imported.
>
>    {% snippet faqs/galaxy/datatypes_understanding_datatypes.md %}
>
>    > <comment-title> The extra files </comment-title>
>    >
>    > The three additional files contain the **reference alkanes**, the **reference spectral library**, and the **sample metadata**. Those files are auxiliary inputs used in the data processing and contain either extra information about the samples or serve as reference data for indexing and identification.
>    > 
>    > The **list of alkanes** (`.tsv` or `.csv`) with retention times and carbon number or retention index is used to compute the retention index of the deconvoluted peaks. The alkanes should be measured in the same batch as the input sample collection.
>    > 
>    > The **reference spectral library** (`.msp`) is used for the identification of spectra. It contains the recorded and annotated mass spectra of chemical standards, ideally from a similar instrument. The unknown spectra which can be detected in the sample can then be confirmed via comparison with this library. The specific library is an in-house library of metabolite standards measured at RECETOX using a GC Orbitrap.
>    > 
>    > The **sample metadata** (`.csv` or `.tsv`) is a table containing information about our samples. In particular, the tabular file contains for each sample its associated sample name, class (QC, blank, sample, etc.), batch number, and injection order. It is possible to add more columns to include additional details about the samples.
>    {: .comment}
>
{: .hands_on}

As a result of this step, we should have in our history a green Dataset collection with three `.raw` files as well as three separate files with reference alkanes, reference spectral library, and sample metadata.

## Convert the raw data to mzML

Our input data are in `.raw` format, which is not suitable for the downstream tools in this tutorial. We use the tool {% tool [msconvert](toolshed.g2.bx.psu.edu/repos/galaxyp/msconvert/msconvert/3.0.20287.2) %} to convert our samples to the appropriate format (`.mzML` in this case).

> <hands-on-title> Convert the raw data to mzML </hands-on-title>
>
> 1. {% tool [msconvert](toolshed.g2.bx.psu.edu/repos/galaxyp/msconvert/msconvert/3.0.20287.2) %} with the following parameters:
>    - {% icon param-collection %} *"Input unrefined MS data"*: `input` (Input dataset collection)
>    - *"Do you agree to the vendor licenses?"*: `Yes`
>    - *"Output Type"*: `mzML`
>    - In *"Data Processing Filters"*:
>        - *"Apply peak picking?"*: `Yes` (This option will compute centroids in the m/z domain.)
>
>    > <comment-title> Centroids </comment-title>
>    >
>    > `msconvert` with selected parameters computes centroids in the m/z domain. MS instruments continuously sample and record signals. Therefore, a mass peak for a single ion in one spectrum consists of multiple intensities at discrete m/z values. Centroiding is the process of reducing these mass peaks to a single representative signal, the centroid. This results in much smaller file sizes without losing too much information. In the further steps, **XCMS** uses the _centWave_ chromatographic peak detection algorithm, which was designed for centroided data. That is the reason we perform the centroiding prior to chromatographic peak detection in this step.
>    {: .comment}
>
{: .hands_on}

## Create the XCMS object

The first part of data processing is using the **XCMS** tool to detect peaks in the MS signal. For that, we first need to take the `.mzML` files and create a format usable by the **XCMS** tool. {% tool [MSnbase readMSData](toolshed.g2.bx.psu.edu/repos/lecorguille/msnbase_readmsdata/msnbase_readmsdata/2.16.1+galaxy0) %} ({% cite gatto2012msnbase %}. {% cite gatto2020msnbase %}) takes as input our files and prepares `RData` files for the first **XCMS** step.

> <hands-on-title> Create the XCMS object </hands-on-title>
>
> 1. {% tool [MSnbase readMSData](toolshed.g2.bx.psu.edu/repos/lecorguille/msnbase_readmsdata/msnbase_readmsdata/2.16.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"File(s) from your history containing your chromatograms"*: `input.mzML` (output of **msconvert** {% icon tool %})
>
>    {% snippet faqs/galaxy/tools_select_collection.md %}
>
>    > <comment-title> Output - `input.raw.RData` </comment-title>
>    >
>    > Collection of `rdata.msnbase.raw` files. `Rdata` file that is necessary in the next step of the workflow. These serve for an internal R representation of **XCMS** objects needed in the further steps.
>    {: .comment}
{: .hands_on}

# Peak detection using XCMS

The first step in the workflow is to detect the peaks in our data using **XCMS**. This part, however, is covered by a [separate tutorial]({{ site.baseurl }}/topics/metabolomics/tutorials/lcms-preprocessing/tutorial.html). Although the tutorial is dedicated to LC-MS data, it can also be followed for our GC-MS data. Therefore, in this section, we do not explain this part of the workflow in detail but rather refer the reader to the dedicated tutorial. Please also pay attention to the parameter values for individual Galaxy tools, as these can differ from the referred tutorial and are adjusted to our dataset.

> <details-title> Skip this step </details-title>
> 
> Since this step is already covered in a [separate tutorial]({{ site.baseurl }}/topics/metabolomics/tutorials/lcms-preprocessing/tutorial.html), it is possible to skip it. Instead, you can continue with [Peak deconvolution]({{ site.baseurl }}/topics/metabolomics/tutorials/gc_ms_with_xcms/tutorial.html#peak-deconvolution) step using a preprocessed **XCMS** object file prepared for you.
>
> > <hands-on-title> Upload data </hands-on-title>
> >
> > 1. Create a new history for this tutorial
> >
> >    {% snippet faqs/galaxy/histories_create_new.md %}
> >
> > 2. Import the following files from [Zenodo]({{ page.zenodo_link }}):
> >
> >    ```
> >    https://zenodo.org/record/7890956/files/XCMS_object.rdata.xcms.fillpeaks
> >    ```
> >
> >    {% snippet faqs/galaxy/datasets_import_via_link.md %}
> >
> >    The format of uploaded file containing **XCMS** object should be `rdata.xcms.fillpeaks`.
> >
> >    {% snippet faqs/galaxy/datatypes_understanding_datatypes.md %}
> >
> {: .hands_on}
>
{: .details}

## Peak picking

The first step is to extract peaks from each of your data files independently. For this purpose, we use the _centWave_ chromatographic peak detection algorithm implemented in {% tool [xcms findChromPeaks (xcmsSet)](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_xcmsset/abims_xcms_xcmsSet/3.12.0+galaxy0) %}.

> <hands-on-title> Peak picking </hands-on-title>
>
>  {% tool [xcms findChromPeaks (xcmsSet)](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_xcmsset/abims_xcms_xcmsSet/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"RData file"*: `input.raw.RData` (output of **MSnbase readMSData** {% icon tool %})
>    - *"Extraction method for peaks detection"*: `CentWave - chromatographic peak detection using the centWave method`
>        - *"Max tolerated ppm m/z deviation in consecutive scans in ppm"*: `3.0`
>        - *"Min,Max peak width in seconds"*: `1,30`
>        - In *"Advanced Options"*:
>            - *"Prefilter step for for the first analysis step (ROI detection)"*: `3,500`
>            - *"Noise filter"*: `1000`
>  
>    You can leave the other parameters with their default values.
>
{: .hands_on}

## Determining shared ions

At this step, you obtain a dataset collection containing one `RData` file per sample, with independent lists of ions. Next, we want to identify the ions shared between samples. To do so, first, you need to group your individual `RData` files into a single one.

> <hands-on-title> Merging files </hands-on-title>
>
>  {% tool [xcms findChromPeaks Merger](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_merge/xcms_merge/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"RData file"*: `input.raw.xset.RData` (output of **xcms findChromPeaks (xcmsSet)** {% icon tool %})
>    - {% icon param-file %} *"Sample metadata file"*: `sample_metadata.tsv` (Input dataset)
>
>    You can leave the other parameters with their default values.
>
{: .hands_on}

Now we can proceed with the grouping and determining shared ions among samples. The aim of this step, called _grouping_, is to obtain a single matrix of ions’ intensities across all samples.

> <hands-on-title> Grouping peaks </hands-on-title>
>
>  {% tool [xcms groupChromPeaks (group)](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_group/abims_xcms_group/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"RData file"*: `xset.merged.RData` (output of **xcms findChromPeaks Merger** {% icon tool %})
>    - *"Method to use for grouping"*: `PeakDensity - peak grouping based on time dimension peak densities`
>        - *"Bandwidth"*: `3.0`
>        - *"Minimum fraction of samples"*: `0.9`
>        - *"Width of overlapping m/z slices"*: `0.01`
>
>    You can leave the other parameters with their default values.
>
{: .hands_on}

## Retention time correction

A deviation in retention time occurs from one sample to another, especially when you inject large sequences of samples. This step aims at correcting retention time drift for each peak among samples.

> <hands-on-title> Retention time correction </hands-on-title>
>
>  {% tool [xcms adjustRtime (retcor)](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_retcor/abims_xcms_retcor/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"RData file"*: `xset.merged.groupChromPeaks.RData` (output of **xcms groupChromPeaks (group)** {% icon tool %})
>    - *"Method to use for retention time correction"*: `PeakGroups - retention time correction based on aligment of features (peak groups) present in most/all samples.`
>        - *"Minimum required fraction of samples in which peaks for the peak group were identified"*: `0.7`
>
>    You can leave the other parameters with their default values.
>
{: .hands_on}

## Second round of determining shared ions

By applying retention time correction, the used retention time values were modified. Consequently, applying this step to your data requires completing it with an additional _grouping_ step.

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
>        - *"If NA values remain, replace them by 0 in the dataMatrix"*: `No`
>
>    You can leave the other parameters with their default values.
>
{: .hands_on}

## Integrating areas of missing peaks

At this point, the peak list may contain `NA` values when peaks were not considered in only some of the samples in the first peak-picking step. In this step, we will integrate the signal in the m/z-rt area of an ion (chromatographic peak group) for samples in which no chromatographic peak for this ion was identified.

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
>    > After the `fillChromPeaks` step, you obtain your final intensity table. At this step, you have everything mandatory to begin analysing
>    > your data:
>    >  - A *sampleMetadata* file (if not done yet, to be completed with information about your samples)
>    >  - A *dataMatrix* file (with the intensities)
>    >  - A *variableMetadata* file (with information about ions such as retention times, m/z) 
>    {: .comment}
>
{: .hands_on}

# Peak deconvolution

The next step is deconvoluting the detected peaks in order to reconstruct the full spectra of the chemical compounds present in the sample, separated by the chromatography and ionized in the mass spectrometer. {% tool [RAMClustR](toolshed.g2.bx.psu.edu/repos/recetox/ramclustr/ramclustr/1.3.0+galaxy0) %} is used to group features based on correlations across samples in a hierarchy, focusing on consistency across samples. While each feature is typically derived from a single compound and represents a fragment, the whole mass spectrum can be used to more accurately identify the precursor compound or molecular ion. **RAMClustR** uses a novel grouping method that operates in an unsupervised manner to group peaks into spectra without relying on the predicting in-source phenomena (in-source fragments or adduct formation) or fragmentation mechanisms.

> <details-title> RAMClust method </details-title>
> 
> The RAMclust similarity scoring utilises a Gaussian function, allowing flexibility in tuning correlational and retention time similarity decay rates independently, based on the dataset and the acquisition instrumentation. The correlational relationship between two features can be described by either MS-MS, MS-idMS/MS, or idMS/MS-idMS/MS values, and we support several correlation methods (e.g. Pearson’s method) to calculate similarity. These are parametrised by `sigma_t` and `sigma_r`, Gaussian tuning parameters of retention time similarity and correlational score, respectively, between feature pairs.
> 
> Similarities are then converted to dissimilarities for clustering. The similarity matrix is then clustered using one of the available methods, such as average or complete linkage hierarchical clustering. The dendrogram is cut using the `cutreeDynamicTree` function from the package `dynamicTreeCut`. For this application, the minimum module size is set to 2, dictating that only clusters with two or more features are returned, as singletons are impossible to interpret intelligently.
> 
> Cluster membership, in conjunction with the abundance values from individual features in the input data, is used to create spectra. The mass to charge ratio is derived from the feature m/z, and the abundance for each m/z in the spectrum is derived from the weighted mean of the intensity values for that feature. These spectra are then exported as an `.msp` formatted file, which can be directly imported by **NIST MSsearch**, or used as input for **MassBank** or **NIST msPepSearch** batch searching.
>
{: .details}

> <hands-on-title> Peak deconvolution </hands-on-title>
>
> 1. {% tool [RAMClustR](toolshed.g2.bx.psu.edu/repos/recetox/ramclustr/ramclustr/1.3.0+galaxy0) %} with the following parameters:
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
>    You can leave the other parameters with their default values.
>
{: .hands_on}

The spectral data comes as an `.msp` file, which is a text file structured according to the **NIST MSSearch** spectra format. `.msp` is one of the generally accepted formats for mass spectral libraries (or collections of unidentified spectra, so called spectral archives), and it is compatible with lots of spectra processing programmes (MS-DIAL, NIST MS Search, AMDIS, etc.). Because `.msp` files are text-based, they can be viewed as simple `txt` files:

> <hands-on-title> Data Exploration </hands-on-title>
>
> Click *"View data"* {% icon galaxy-eye %} icon next to the dataset in the Galaxy history. The contents of the file would look like this:
> 
> 
> {% snippet faqs/galaxy/datasets_icons.md %}
> 
> 
> ```
> NAME:C001
> IONMODE:Negative
> SPECTRUMTYPE:Centroid
> RETENTIONTIME:383.27
> Num Peaks:231
> 217.1073 64041926
> 243.0865 35597866
> 257.1134 31831229
> 224.061 27258239
> 258.11 24996353
> 241.0821 23957171
> 315.1188 13756744
> ...
>
> NAME:C002
> IONMODE:Negative
> SPECTRUMTYPE:Centroid
> RETENTIONTIME:281.62
> Num Peaks:165
> 307.1573 299174880
> 147.0654 298860831
> 149.0447 287809889
> 218.1066 118274758
> 189.076 112486871
> 364.1787 75134143
> 191.0916 52526567
> 308.1579 52057158
> ...
> ```
>
> > <details-title> Negative ion mode </details-title>
> >
> > You might wonder how can the ionisation mode (_IONMODE_) for GC-MS data be negative when using electron impact (EI+) ionization. This is, of course, incorrect. This is actually just a default behaviour of **RAMClustR**. We can optionally change this by providing **RAMClustR** an experiment definition file. This file can be created manually or using the {% tool [RAMClustR define experiment](toolshed.g2.bx.psu.edu/repos/recetox/ramclustr_define_experiment/ramclustr_define_experiment/1.0.2) %} tool. There we can specify annotations such as what instrument we used or ionisation mode (which was EI+ in our case), and this will be transfered to the `.msp` file. Finally, we can provide such a file as an input to **RAMClustR** in the _Extras_ inputs section.
> >
> {: .details}
>
{: .hands_on}

`.msp` files can contain one or more mass spectra, these are split by an empty line. The individual spectra essentially consist of two sections: metadata and peaks. The metadata consists of compound name, spectrum type (which is centroid in this case), ion mode, retention time, and the number of *m/z* peaks. If the compound has not been identified, as in our case, the **NAME** can be any arbitrary string. It is best, however, if that string is unique within that `.msp` file. The metadata fields are usually unordered, so it is quite common for one `.msp` file to contain **NAME** as the first metadata key and for another `.msp` to have this key somewhere in the middle. The keys themselves also aren't strictly defined and can have different names (e.g., **compound_name** instead of **NAME**). The `msp` format is inherently not standardized and various flavours exist.

However, as we can observe, the metadata part is rather incomplete. We would like to gather more information about the detected spectra and identify the specific compounds corresponding to them.

The second output file is the so called **Spec Abundance** table, containing the expression of all deconvoluted spectra across all samples. It is the corresponding element to the peak intensity table obtained from **XCMS**, only for the whole deconvoluted spectrum. This file is used for further downstream processing and analysis of the data and eventual comparisons across different sample types or groups. For more information on this downstream data processing see this [tutorial]({{ site.baseurl }}/topics/metabolomics/tutorials/lcms-dataprocessing/tutorial.html).

# Retention index calculation

The retention index ({% cite van1963generalization %}) is a way how to convert equipment- and experiment-specific retention times into system-independent normalised constants for GC-based experiments. The retention index of a compound is computed from the retention time by interpolating between the retention times adjacent alkanes which are assigned a fixed retention index. This can be different for the individual chromatographic system, but the derived retention indices are independent and allow comparing values measured by different analytical laboratories.

We use the package {% tool [RIAssigner](toolshed.g2.bx.psu.edu/repos/recetox/riassigner/riassigner/0.3.4+galaxy1) %} to compute retention indices for files in the `.msp` format using an indexed reference list of alkanes in `.csv` or `.msp` format. The output follows the same format as the input but with added retention index values. These can be used at a subsequent stage to improve compound identification ({% cite kumari2011applying %}). Multiple computation methods (e.g. piecewise-linear or cubic spline) are supported by the tool.

> <hands-on-title> Retention index calculation </hands-on-title>
>
> 1. {% tool [RIAssigner](toolshed.g2.bx.psu.edu/repos/recetox/riassigner/riassigner/0.3.4+galaxy1) %} with the following parameters:
>    - In *"Query dataset"*:
>        - {% icon param-file %} *"Query compound list"*: `Mass spectra from RAMClustR` (output of **RAMClustR** {% icon tool %})
>    - In *"Reference dataset"*:
>        - {% icon param-file %} *"Reference compound list"*: `reference_alkanes.csv` (Input dataset)
>
>    You can leave the other parameters with their default values.
>
>    > <comment-title> Minutes vs. seconds </comment-title>
>    >
>    > You might notice that in the last **XCMS** step, we converted the retention times to minutes. However, in this step, we are using seconds again. The reason is that **RAMClustR** converted them internally again to seconds. Regarding the reference compound list, this database already has its retention times in seconds.
>    {: .comment}
>
>    > <details-title> Kováts method </details-title>
>    >
>    > The Kováts retention index (or Kováts index) of a compound is its retention time normalized to the retention times of adjacently eluting n-alkanes. It is based on the fact that the logarithm of the retention time converges to the number of carbon atoms in the alkane. For an isothermal chromatogram, you use the following equation to calculate the Kováts index: 
>    > 
>    > $$I = 100 \times \frac{\mathtt{log}~t_x - \mathtt{log}~t_n}{\mathtt{log}~t_{n+1} - \mathtt{log}~t_n}$$ 
>    > 
>    > where $$t_n$$ and $$t_{n+1}$$ are retention times of the reference n-alkane hydrocarbons eluting immediately before and after chemical compound *X* and $$t_x$$ is the retention time of compound *X*.
>    > 
>    > In practice, first, you run a chromatogram of a standard alkane mixture in the range of interest. Then you do a co-injection of your sample with the standard alkanes. The main advantage of this method is reproducibility.
>    >
>    {: .details}
>
{: .hands_on}

# Identification

To annotate and putatively identify the deconvoluted spectra, we compare them with a reference spectral library. This library contains spectra of standards measured on the same instrument for optimal comparability. The **matchms** package is used for spectral matching. It was built to import and apply different similarity measures to compare large numbers of spectra. This includes common cosine scores but can also easily be extended with custom similarity measures. For more information on how to use matchms from Python check out this [blog](https://blog.esciencecenter.nl/build-your-own-mass-spectrometry-analysis-pipeline-in-python-using-matchms-part-i-d96c718c68ee).

> <details-title> Reference spectral library </details-title>
> 
> Fragmentation patterns of ions (spectra) are widely used in the process of compound identification for EI+ ionization in complex mixtures. The number of identifiable compounds and associated data keeps growing with the increasing sensitivity and resolution of mass spectrometers. These data are clustered into collections of chemical structures and their spectra, also called spectral libraries ({% cite stein1995chemical %}, {% cite stein2012mass %}), which can be used for fast, reliable identifications for any compound whose fragmentation pattern is measured by the instrument. 
> 
> To identify the compounds after GC/MS analysis, the *library searching* is performed for every detected spectrum. This locates the most similar spectra in the reference library, providing a list of the potential identifications sorted by computed similarity. Indeed, a library search should also yield the confidence of compound identification, often expressed as a similarity score.
>
> Reference libraries' users usually expect their reliability to be quite high. However, there are many sources of error which we need to be aware of. A particularly error-prone aspect is annotation. Even an excellent spectrum becomes worthless for a wrongly annotated compound. A variety of computer-aided quality control measures can help to improve these issues, but libraries always require some level of expert human decision-making before a spectrum is added to the library. This is also important for redundant and unidentified spectra.
>
{: .details}


## Compute similarity scores

We use the cosine score with a greedy peak pairing heuristic to compute the number of matching ions with a given tolerance and the cosine scores for the matched peaks. Compared to more traditional methods of forward and reverse matching, **matchms** is performing matching on peaks present both in the measured sample and the reference data. This inherently yields higher scores than the forward and reverse scoring used by NIST.

> <details-title> Cosine score methods </details-title>
>
> One of the most common methods to compute the similarity score between two spectra is computing the cosine of the *angle* between them. Each spectrum can be represented as a vector with an axis along each of its *m/z* values. Then the cosine of the angle is basically a normalised *dot product* of these two vectors. The intensity and *m/z* weighting corrections have been used to optimise performance. For example, higher mass peaks can be given greater weight since they are more discriminating than lower mass peaks. *Reverse* scores assign no weight to peaks present only in the measured sample spectrum under the assumption that they arise from impurities.
>
{: .details}

> <hands-on-title> Compute similarity scores </hands-on-title>
>
> 1. {% tool [matchms similarity](toolshed.g2.bx.psu.edu/repos/recetox/matchms_similarity/matchms_similarity/0.20.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Queries spectra"*: `RI using kovats of Mass spectra from RAMClustR` (output of **RIAssigner** {% icon tool %})
>    - *"Symmetric"*: `No` (if we were to query our spectra agains itself, we would select `Yes`)
>        - {% icon param-file %} *"Reference spectra"*: `reference_spectral_library.msp` (downloaded file from Zenodo)
>    - In *"Algorithm Parameters"*:
>        - *"tolerance"*: `0.03`
>    - *"Apply RI filtering"*: `Yes` (we use this since our data is GC)
>
>    You can leave the other parameters with their default values.
>
>    > <comment-title> Used reference spectra </comment-title>
>    > The RECETOX Metabolome HR-[EI+]-MS library is a collection of mostly endogenous compounds from the MetaSci Human Metabolite Library. Analytes underwent methoximation/silylation prior to acquisition. Spectra were acquired at 70 eV on Thermo Fisher Q Exactive™ GC Orbitrap™ GC-MS/MS at 60000 resolving power.
>    {: .comment}
>
{: .hands_on}

> <details-title> Overview of the spectral similarity scores </details-title>
>
> ## Cosine Greedy
> The cosine score, also known as the dot product, is based on representing the similarity of two spectra through the cosine of an angle between the vectors that the spectra produce. Two peaks are considered as matching if their *m/z* values lie within the given tolerance. Cosine greedy looks up matching peaks in a "greedy" way, which does not always lead to the most optimal alignments.
>
> This score was among the first to be used for looking up matching spectra in spectral libraries and, to this day, remains one of the most popular scoring methods for both library matching and molecular networking workflows.
>
> ## Cosine Hungarian
> This method computes the similarities in the same way as the *Cosine Greedy* but with a difference in *m/z* peak alignment. The difference lies in that the Hungarian algorithm is used here to find matching peaks. This leads to the best peak pairs match but can take significantly longer than the "greedy" algorithm.
>
> ## Modified Cosine
> Modified Cosine is another, as its name states, representative of the family of cosine-based scores. This method aligns peaks by finding the best possible matches and considers two peaks a match if their *m/z* values are within tolerance before or after a mass-shift is applied. A mass shift is essentially a difference of precursor-*m/z* of two compared spectra. The similarity is then again expressed as a cosine of the angle between two vectors.
>
> ## Neutral Losses Cosine
> Neutral Loss metric works similarly to all described above with one major difference: instead of encoding the spectra as "intensity vs *m/z*" vector, it encodes it to an "intensity vs *Δm/z*", where delta is computed as an *m/z* difference between precursor and a fragment *m/z*. This, in theory, could better capture the underlying structural similarities between molecules.
>
{: .details}

## Format the output

The output of the previous step is a `json` file. This format is very simple to read for our computers, but not so much for us. We can use the {% tool [matchms scores formatter](toolshed.g2.bx.psu.edu/repos/recetox/matchms_formatter/matchms_formatter/0.20.0+galaxy0) %} to convert the data to a tab-separated file with a scores matrix.

The output table contains the scores and number of matched ions of the deconvoluted spectra with the spectra in the reference library. The raw output can be filtered to only contain the top matches (3 by default) and/or to contain only pairs with a score and number of matched ions larger than provided thresholds.

> <hands-on-title> Format the output </hands-on-title>
>
> 1. {% tool [matchms scores formatter](toolshed.g2.bx.psu.edu/repos/recetox/matchms_formatter/matchms_formatter/0.20.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Scores object"*: `CosineGreedy scores` (output of **matchms similarity** {% icon tool %})
>
>    You can leave the other parameters with their default values.
>
{: .hands_on}

> | query | reference | CosineGreedy_score | CosineGreedy_matches |
> |-------|-----------|---------|--------|
> | C044  | Guanine_3TMS    | 0.987    | 19    |
> | C068  | Tryptophan_3TMS | 0.981    | 9   |
> | C053  | Norleucine_2TMS | 0.980    | 12   |
> | ...   | ...             | ... | ...  |
{: .matrix}

At this stage, all steps are complete: we have the list of identified spectra corresponding to a compound from the reference database. Since this process is dependent on the parameters we used along the way, the result should be considered dependent on them. To enumerate this dependency, each potential compound has been assigned a confidence score.

# Conclusion

In this tutorial, we show how untargeted GC-MS data can be processed to obtain tentative identification of compounds contained in the analyzed samples.
Using **XCMS**, we obtain a peak and intensity table from our files which contains information about the detected ions (such as m/z and retention time) and their intensities across all samples.
Afterwards, we used **RAMClustR** to reconstruct the fragmentation spectra which we assume to originate from the analyzed compounds contained within our samples.
We use a correlation based approach which allows us to reconstruct also low-abundant features as long as their intensities are correlating across samples.
This method is different than the sample-wise deconvolution using **CAMERA** ({% cite Kuhl2012 %}).
To have more information available for compound identification, we computed the retention index using **RIAssigner** and added it to our deconvoluted spectra.
We then finally compare the deconvoluted spectra against the reference library using **matchms**, using a retention index threshold to pre-filter our matches and then computing the cosine similarity for the spectra within the retention index tolerance.

This workflow leverages the correlation of features across samples to improve the annotation performance, trying to capture even low abundant signals.
It is suitable for both high- and low-resolution GC-MS data acquired in either profile or centroid mode by adapting the individual tool parameters or in-/excluding individual steps which are specific to a certain data format.

Before further statistical analysis, the information from the scores table can be included in the library of deconvoluted spectra manually to add the putative identifications. The further processing can then leverage this information and the intensity information contained in the **SpecAbundance** output from **RAMclustR**.

Note that this tutorial doesn't cover all steps which are possible or compatible with this workflow. Other options include the {% tool [matchms filtering](toolshed.g2.bx.psu.edu/repos/recetox/matchms_filtering/matchms_filtering/0.17.0+galaxy0) %} tool with which you can normalize your spectra or remove low-abundant ions before matching.
Another option for normalization outside of the built-in normalization tools of **RAMClustR** is the {% tool [WaveICA](toolshed.g2.bx.psu.edu/repos/recetox/waveica/waveica/0.2.0+galaxy2) %} tool for correction of batch effects.
