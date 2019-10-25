--- 
layout : tutorial_hands_on

title : 'Mass spectrometry : MSMS analysis with msPurity package'
level : Introductory
enable : true
zenodo_link : '' 
questions : 
    - Which biological questions are addressed by the tutorial? 
    - Which bioinformatics techniques are important to know for this type of data? 
objectives : 
    - The learning objectives are the goals of the tutorial 
    - They will be informed by your audience and will communicate to them and to yourself what you should focus on during the course 
    - They are single sentences describing what a learner should be able to do once they have done completed tutorial 
    - You can use Bloom's Taxonomy to write effective learning objectives 
time_estimation : 2H 
key_points : 
    - The take-home messages 
    - They will appear at the end of the tutorial 
requirements:
  - type: "internal"
    topic_name: metabolomics
    tutorials: 
      - lcms
contributors : 
    - jsaintvanne

--- 

# Introduction
{:.no_toc}

blabla 
 
> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Preprocessing with XCMS

The first step in the workflow is the pre-processing of the raw data with **XCMS** ({% cite Smith2006 %}).

**XCMS** {% icon tool %} is a free and open source software dedicated to pre-processing of any type of mass spectrometry acquisition files from low to
high resolution, including FT-MS data coupled with different kind of chromatography (liquid or gas). This software is
used worldwide by a huge community of specialists in metabolomics using mass spectrometry methods.

This software is based on different algorithms that have been published, and is provided and maintained using R software.

**MSnbase readMSData** {% icon tool %}, prior to **XCMS** {% icon tool %} is able to read files with open format as `mzXML`, `mzMl`, `mzData` and `netCDF`, which are independent of the constructors' formats. The **XCMS** {% icon tool %} package itself is composed of R functions able to extract, filter, align and fill gap, with the possibility to annotate isotopes,
adducts and fragments using the R package CAMERA ({% cite CAMERA %}). This set of functions gives modularity, and thus is particularly well
adapted to define workflows, one of the key points of Galaxy.

To be able to process your MS/MS datas, we need to start with the peakpicking of MS datas. One Galaxy Training already explains how to process with your MS datas. You should follow this link and complete this tutorial : [Mass spectrometry: LC-MS analysis](https://galaxyproject.github.io/training-material/topics/metabolomics/tutorials/lcms/tutorial.html). 

> ### {% icon hands_on %} Hands-on: Resume of the XCMS preprocessing
>
> 1. Import your datas into a Galaxy collection
> 
>    For a good workflow, you need to create a collection in Galaxy which contains **all** your datas. It can be MS only datas, MSMS with MS also files but it can **NOT** be only MSMS datas.
>
>    {% include snippets/create_dataset_collection.md %}
>    {% include snippets/import_via_link.md collection=true format="mzml" collection_name="sacurine" renaming=false %}
>    {% include snippets/import_from_data_library.md astype="as a Collection" %}
>
> 2. Prepare your MS datas with *MSnbase readMSData*
>
>    **MSnbase readMSData** {% icon tool %}, prior to **XCMS** {% icon tool %} is able to read files with open format as `mzXML`, `mzMl`, `mzData` and `netCDF`, which are independent of the constructors' formats.
>    - Input : your collection with all your files
>    - Output : a new collection with your files and their information after `readMSData` function. Format : `collection.raw.RData`.
>
> 3. First XCMS step : *peak picking*
>    
>    Open and run the **xcms findChromPeaks (xcmsSet)** {% icon tool %} tool and set your parameters.
>    - Input : collection of Rdata object(s) obtained just before (format : `collection.raw.RData`)
>    - Output : a collection with informations about peak picking for each files. Format : `collection.raw.xset.RData`.
>
>
> 4. *Merge samples* in one dataset with a sampleMetadata file
>
>    To merge your datas, you need to **input a sampleMetadata file** containing namefiles and their metadata informations like their class for example. If you don't add a sampleMetadata file, the **xcms findChromPeaks Merger** {% icon tool %} tool will **group all your files together**. You can also **create your sampleMetadata file** with W4M Galaxy tool **xcms get a sampleMetadata file** {% icon tool %} with the following parameters: *"RData file"* outputed from **MSnbase readMSData** {% icon tool %}.
>    - Input : collection of RData object(s) obtained during peak picking (format : `collection.raw.xset.RData`) and optionnaly a sampleMetadata file in `tabular` format.
>    - Output : one RData file named `xset.merged.RData` which contains a list of all your datas processed just before.
>
> 5. Second XCMS step : *determining shared ions across samples*
>
>    Open **xcms groupChromPeaks (group)** {% icon tool %} tool and run it with your parameters.
>    - Input : the RData file obtained just before with all datas and peaks informations for each files (format : `xset.merged.RData`).
>    - Output : one Rdata file named `xset.merged.groupChromPeaks.RData` containing all peaks grouped according to their mz and retention time.
>
> 6. Optional XCMS step : *retention time correction*
>
>    This step is optionnal, it aims to correct retention time drift for each peak among samples. You have to use the **xcms adjustRtime (retcor)** {% icon tool%} tool to correct this retention time.
>    - Input : the Rdata file named `xset.merged.groupChromPeaks.RData`
>    - Output : new RData file named `xset.merged.groupChromPeaks.adjustRtime.RData`
>    
>    > ### {% icon comment %} Important : Have to group again !
>    >
>    > Retention time have been modified. Consequently, it requires to complete it with an additionnal grouping step with **xcms groupChromPeaks (group)** {% icon tool%}. You will obtain a new Rdata file named `xset.merged.groupChromPeaks.adjustRtime.groupChromPeaks.RData`
>    {: .warning}
>
> 7. Final XCMS step *integrating areas of missing peaks*
>
>    This last step can be run after grouping your peaks and don't need the retention time correction. Run **xcms fillChromPeaks (fillPeaks)** {% icon tool%} with the following things.
>    - Input : you Rdata file `xset.merged.groupChromPeaks.*.RData`
>    - Output : a new RData file named `xset.merged.groupChromPeaks.*.fillChromPeaks.RData`
>
> > ### {% icon comment %} Important : Be careful of the file format
> >
> > During each step of preprocessing, your file has its format changed and can have also its name changed.
> > To be able to continue to MSMS processing, you need to have a RData object wich is **merged and grouped** (step 4 and 5) at least. It means that you should have a file named `xset.merged.groupChromPeaks.RData` (and maybe with some step more in it).
> {: .warning}
{: .hands_on}




# msPurity package

**msPurity** is a R package developped by Birmingham University team and published in 2017 ({% cite Lawson2017 %}).

This R package was developped to :
1. Assess the spectral quality of fragmentation spectra by evaluating the "precursor ion purity". 
2. Process fragmentation spectra. 
3. Perform spectral matching. 

**What is precursor ion purity?** 

What we call **"precursor ion purity"** is a measure of the contribution of a selected precursor peak in an isolation window used for fragmentation. The simple calculation involves dividing the intensity of the selected precursor peak by the total intensity of the isolation window. When assessing MS/MS spectra this calculation is done before and after the MS/MS scan of interest and the **purity is interpolated at the recorded time of the MS/MS acquisition**. Additionally, isotopic peaks can be removed, low abundance peaks are removed that are thought to have limited contribution to the resulting MS/MS spectra and the isolation efficiency of the mass spectrometer can be used to normalise the intensities used for the calculation.

To install this package, start R (version "3.6") and enter:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("msPurity") 
```

## purityA function

The `purityA` function is called to calculate the precursor purity of the fragmentation results. 

```R
mzMLpths <- c(file1, file2,...)
pa  <- purityA(mzMLpths)
```
In Galaxy instance, you need to run the **msPurity.purityA** {% icon tool %} tool with the collection containing all your mzML files. It is the collection with which we started the peak picking processus at the beginning of this tutorial.

  - Input : your collection with all your files (`mzML`, `mzXML`,...)
  - Output : one RData file containing all MS/MS datas from your files (`msPurity.purityA_on_yourdata.RData`)


## frag4feature function

The `frag4feature` function will link the fragmentation data back to the XCMS feature. 

```R
pa <- frag4feature(pa, xset)
```
You need to find the **msPurity.frag4feature** {% icon tool %} tool in the Galaxy instance to run this function. It will give you all your MS/MS datas which have a good precursor found during the XCMS processus.

  - Input : You need 2 inputs for this function
    1. Your `xset.merged.groupChromPeaks.*.RData` file containing all the XCMS processus results
    2. The Rdata file obtained just before with all your MS/MS datas (`msPurity.purityA_on_yourdata.RData`)
  - Output : one updated RData file named `msPurity.frag4feature_on_yourdata.RData`


## filterFragSpectra function

The fragmentation can be filtered prior to averaging using the `filterFragSpectra` function

```R
pa <- filterFragSpectra(pa)
```

## averageFragSpectra function

Averaging of the fragmentation spectra can be done with either `averageAllFragSpectra` or with `averageIntraFragSpectra` and `averageInterFragSpectra`. This will depend if the user wishes to treat the fragmentation spectra from within a file and between files. Another alternative is to ignore the averaging completely and just use the non-averaged fragmentation spectra for the spectral matching.

If the inter and intra fragmentation scans are to be treated differently the following should be followed:

```R
pa <- averageIntraFragSpectra(pa) # use parameters specific to intra spectra 
pa <- averageInterFragSpectra(pa) # use parameters specific to inter spectra
```

If the inter and intra fragmentation scans are to be treated the same the following workflow should be used.

```R
pa <- averageAllFragSpectra(pa) 
```

## creataDatabase function

An SQLite database is then created of the LC-MS/MS experiment. The SQLite schema of the spectral database can be detailed here.

```R
td <- tempdir()
q_dbPth <- createDatabase(pa, xset, outDir = td, dbName = 'lcmsms-processing.sqlite')
```

## spectralMatching function

blabla

## combineAnnotations function

blabla

## createMSP function

blabla

## flagRemove function

blabla


# Conclusion 
{:.no_toc} 

blabla 

>    > ### {% icon comment %} Comment
>    >
>    > Here we provided the sampleMetadata file so we know that the upload led to a 'tabular' file. But from experience we also know that it can happen that, when uploading a sampleMetadata table, user obtained other inappropriate types of data. This is generally due to the file not following all the requirements about the format (*e.g.* wrong separator, or lines with different numbers of columns).
>    > Thus, we highly recommend that you always take a second to check the data type after the upload. This way you can handle the problem right away if you appear to get one of these obvious issues.
>    {: .comment}