---
layout: tutorial_hands_on

title: "Mass spectrometry imaging data preprocessing"
edam_ontology: "topic_0121"
zenodo_link: ""
questions:
  - "How to normalize data?"
  - "How to prepare the spectra for peak detection"
  - "How to do peak detection and follow up stepsn"

objectives:
  - "Understanding the different preprocessing steps that are necessary to prepare the MSI raw data for statistical follow up."
time_estimation: "30min"
key_points:
  - "Raw intensities need to be normalized"
  - "Raw spectra need to be preprocessed to improve peak detection. "
  - "Only real peaks are interesting to keep, noise and isotopic peaks should be removed. "
contributors:
  - foellmelanie
  - bgruening
---

# Introduction
{:.no_toc}

Preprocessing and peak detection are key steps for the analysis of mass spectrometry imaging data as they help to base the following statistical analysis on real peaks rather than on noise peaks. Several preprocessing steps have to be performed in order to prepare the mass spectra for peak picking, this includes spectra smoothing and baseline removal. Intensity values need to be normalized by transformation and spectra wide calibration. Then peak detection is performed, peak isotopes are removed and only peaks that occur in a decent number of spectra are kept. 


> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Preparing the spectra for peak detection

Several preprocessing steps are necessary to improve the performance of the peak detection step: 
  1. Transformation
  2. Smoothing
  3. Baseline removal
  4. SPectra wise calibration
  5. Warping/m/z calibration

TODO: provide test files; probably many m/z can be cut off to have a nice image with many pixels but not many m/z left?! Provide information only after QC - murine kidney etc...m/z are cut...Hauptfrage: Wie krüppelig sind unsere einzelnen Spektren in der Niere? 


> ### {% icon hands_on %} Hands-on: Preprocessing of MSI data
>
> 1. Create a new history and upload the MSI data and tabular file with m/z values of internal Calibrants as described here [MSI quality control]({{site.baseurl}}//galaxy-data-manipulation/tutorials/get-data/slides.html#18) 
>
> 2. Select the **MALDIquant preprocessing** tool
>    - Select `Read in only spectra of interest`
>    - Select 
>    - Set 
>    - Select
>    - Select
>
>    > ### {% icon tip %} Tip: Preprocessing of large files
>    > * MALDIquant might have trouble to process large files, but it has an option to `Read in only spectra of interest`. 
>    > * Alternatively, the **MSI filtering** tool could be used to split the data into several smaller datasets.
>    > * Another option, with better support for large files, is the **MSI preprocessing** tool from the Cardinal software suite.
>    {: .tip}
>
> 3. Maybe another QC tool? E.g. to check m/z shift? 
{: .hands_on}


# Detection of peaks in mass spectra from MSI data

Now we are ready for peak detection
  1. Transformation
  2. Smoothing
  3. Baseline removal
  4. SPectra wise calibration

> ### {% icon hands_on %} Hands-on: Peak detection for MSI data
> 1. Select the **MALDIquant peak detection** tool
>
> 2. Select 
> 2. Select 
> 3. Select 
> 4. Select `ppm range` 
> 5. Press **Execute**
{: .hands_on}


# Any follow up??? 

To read and extract information from MSI data files special tools are needed. The imzML exporter provides the option to export spectra, feature and intensity data from imzML or Analyze7.5
files into tabular files. This may be helpful to follow up on interesting aspects of the quality report with exact values. Standard statistical programmes and software to perform further plots are normally compatible with tabular files. Care has to be taken, that  In case the intensity matrix of an unprocessed file is exported the tabular file will be huge and it might happen that 


- spectra output: spectra in rows - for each spectrum: name, x and y coordinates,order, number of peaks (intensities > 0), total ion chromatogram (TIC), highest m/z feature per spectrum, optional count of input m/z per spectrum, optional spectrum annotation
- mz feature output: m/z in rows - for each m/z: name, m/z, mean, median, standard deviation (sd), standard error of the mean (sem), sum of all intensities per m/z, number of peaks (intensity > 0) per m/z

Bsp tables spectra, feature and int. matrix information???

> ### {% icon hands_on %} Hands-on: Export spectra and feature information from MSI data
>
> 1. Run **MSI imzML exporter** {% icon tool %} on the imzML file.
> 2. Select 'mz feature output' and 'pixel output' in `Multiple output files can be selected`.
> 3. Press **Execute**
>
>    > ### {% icon question %} Questions
>    >
>    > 1. How many spectra (pixel) does the dataset contain?
>    > 2. Which m/z feature is the highest in pixel X,Y OR which m/z feature is the one most often found (counting pixels?? hard)
>    > 3. How many numeric properties does the m/z feature output contain?
>    > 4. What is the m/z feature with the highest mean intensity? (Problematic in case of very many m/z where people start to scroll by hand???)
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. YYY - The pixel output file has XXX lines minus one header lines is YYY m/z features, as every row contains one spectrum.
>    > > 2. TTT - Run `Sort data on a column????` with pixelname xy_20_20
>    > > 3. 7 - The output file has 8???? columns and except for the name column, all columns contain numeric properties that can be used for further calculations.
>    > > 4. ZZZ - This can be either done by hand or with the tool `Sort data on a column????`
>    > {: .solution }
>    {: .question}
{: .hands_on}


# Concluding remarks
{:.no_toc}

- imzML (continuous or processed file type) and Analyze7.5 can be uploaded into Galaxy via the composite upload. 
- The MSI quality report should be applied several times in the analysis pipeline and especially to control that the preprocessing is going well. (Link preprocessing tutorial)
- The MSI imzML exporter can be used any time during the analysis when it is necessary to dig deeper into the data. Galaxy provides many text manipulation tools that can directly be applied on the exported tabular file to filter, sort, plot, ... the data.

????    This tutorial is based upon parts of the GalaxyP-101 tutorial (https://usegalaxyp.readthedocs.io/en/latest/sections/galaxyp_101.html).
