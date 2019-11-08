--- 
layout : tutorial_hands_on

title : 'Mass spectrometry : MSMS analysis with msPurity package'
level : Introductory
zenodo_link : 'https://zenodo.org/record/3244991' 
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
requirements :
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

The first step in the workflow is the pre-processing of the raw data with **XCMS** {% icon tool %} ({% cite Smith2006 %}).

**XCMS** {% icon tool %} is a free and open source software dedicated to pre-processing of any type of mass spectrometry acquisition files from low to
high resolution, including FT-MS data coupled with different kind of chromatography (liquid or gas). This software is
used worldwide by a huge community of specialists in metabolomics using mass spectrometry methods.

This software is based on different algorithms that have been published, and is provided and maintained using R software.

**MSnbase readMSData** {% icon tool %}, prior to **XCMS** {% icon tool %} is able to read files with open format as `mzXML`, `mzMl`, `mzData` and `netCDF`, which are independent of the constructors' formats. The **XCMS** {% icon tool %} package itself is composed of R functions able to extract, filter, align and fill gap, with the possibility to annotate isotopes,
adducts and fragments using the R package CAMERA ({% cite CAMERA %}). This set of functions gives modularity, and thus is particularly well
adapted to define workflows, one of the key points of Galaxy.

To be able to process your MS/MS datas, we need to **start with the peakpicking of MS datas**. One Galaxy Training already explains how to process with your MS datas. You should **follow this link and complete this tutorial** : [Mass spectrometry: LC-MS analysis](https://galaxyproject.github.io/training-material/topics/metabolomics/tutorials/lcms/tutorial.html). For MS/MS analysis you **don't really need to finish** this previous tutorial but for a better understanding of your datas, it is recommanded. In this tutorial, you **just have to compute** your datas with the **following steps** briefly describe in the *details* part below.

> ### {% icon details %} Some help : Resume of the XCMS preprocessing
>
> Here is a resume of the **XCMS** {% icon tool %} preprocessing. It has few steps that are an obligation to be able to obtain the final file. With this file, we can continue the workflow to process our MS/MS datas with **msPurity** {% icon tool %} package.
>
> > ### {% icon solution %} 1 - Import your datas into a Galaxy collection
> > 
> > For a good workflow, you need to create a collection in Galaxy which contains **all** your datas. It can be MS only datas, MSMS with MS also files but it can **NOT** be only MSMS datas.
> > 
> > > ### {% icon tip %} Tip: Create a dataset collection
> > > 1. Click on {% icon galaxy-selector %} icon (**Operation on multiple datasets**) on the top of the history
> > > 2. Select all the datasets for the collection
> > > 3. Expand **For all selected** menu
> > > 4. Select **Build dataset list**
> > > 5. Enter a name for the collection
> > > 6. Tick **Hide original elements?**
> > > 5. Click on **Create list** (and wait a bit)
> > {: .tip}
> >
> > > ### {% icon tip %} Tip: Importing data via links
> > >
> > > * Copy the link location
> > > * Open the Galaxy Upload Manager ({% icon galaxy-upload %} on the top-right of the tool panel)
> > > {% if include.collection %}
> > > * Click on **Collection** on the top
> > > {% endif %}
> > > {% if include.collection_type %}
> > > * Click on **Collection Type** and select `{{ include.collection_type }}`
> > > {% endif %}
> > > * Select **Paste/Fetch Data**
> > > * Paste the link into the text field
> > > {% if include.link %}
> > >   `{{ include.link }}`
> > > {% endif %}
> > > {% if include.link2 %}
> > >   `{{ include.link2 }}`
> > > {% endif %}
> > > {% if include.format %}
> > > * Change **Type** from "Auto-detect" to `{{ include.format }}`
> > > {% endif %}
> > > {% if include.genome %}
> > > * Change **Genome** to `{{ include.genome }}`
> > > {% endif %}
> > > * Press **Start**
> > > {% if include.collection %}
> > > * Click on **Build** when available
> > > {% if include.pairswaptext %}
> > > * Ensure that the forward and reverse reads are set to {{ include.pairswaptext }}, respectively.
> > >     * Click **Swap** otherwise
> > > {% endif %}
> > > * Enter a name for the collection
> > > {% if include.collection_name_convention %}
> > >     * A useful naming convention is to use {{ include.collection_name_convention }}
> > > {% endif %}
> > > {% if include.collection_name %}
> > >     * {{ include.collection_name }}
> > > {% endif %}
> > > * Click on **Create list** (and wait a bit)
> > > {% else %}
> > > * **Close** the window
> > > {% endif %}
> > > {% if include.renaming == undefined or include.renaming == true %}
> > > By default, Galaxy uses the URL as the name, so rename the files with a more useful name.
> > > {% endif %}
> > {: .tip}
> > > ### {% icon tip %} Tip: Importing data from a data library
> > >
> > > As an alternative to uploading the data from a URL or your computer, the files may also have been made available from a *shared data library*:
> > >
> > > * Go into **Shared data** (top panel) then **Data libraries**
> > > {% if include.path %}
> > > * {{ include.path }}
> > > {% else %}
> > > * Find the correct folder (ask your instructor)
> > > {% endif %}
> > > * Select the desired files
> > > * Click on the **To History** button near the top and select **{{ include.astype | default: "as Datasets" }}** from the dropdown menu
> > > * In the pop-up window, select the history you want to import the files to (or create a new one)
> > > * Click on **Import**
> > {: .tip}
> >
> {: .solution}
>
> > ### {% icon solution %} 2 - Prepare your MS datas with *MSnbase readMSData* {% icon tool %}
> >
> > **MSnbase readMSData** {% icon tool %}, prior to **XCMS** {% icon tool %} is able to read files with open format as `mzXML`, `mzMl`, `mzData` and `netCDF`, which are independent of the constructors' formats.
> >   - **Input** : your collection with all your files
> >   - **Output** : a new collection with your files and their information after `readMSData` function. Format : `collection.raw.RData`.
> {: .solution}
>
> > ### {% icon solution %} 3 - First XCMS step : *peak picking* {% icon tool %}
> >    
> >    Open and run the **xcms findChromPeaks (xcmsSet)** {% icon tool %} tool and set your parameters.
> >    - **Input** : collection of Rdata object(s) obtained just before (format : `collection.raw.RData`)
> >    - **Output** : a collection with informations about peak picking for each files. Format : `collection.raw.xset.RData`.
> {: .solution}
> 
> > ### {% icon solution %} 4 - *Merge samples* {% icon tool %} in one dataset with a sampleMetadata file
> >
> >    To merge your datas, you need to **input a sampleMetadata file** containing namefiles and their metadata informations like their class for example. If you don't add a sampleMetadata file, the **xcms findChromPeaks Merger** {% icon tool %} tool will **group all your files together**. You can also **create your sampleMetadata file** with W4M Galaxy tool **xcms get a sampleMetadata file** {% icon tool %} with the following parameters: *"RData file"* outputed from **MSnbase readMSData** {% icon tool %}.
> >    - **Input** : collection of RData object(s) obtained during peak picking (format : `collection.raw.xset.RData`) and optionnaly a sampleMetadata file in `tabular` format.
> >    - **Output** : one RData file named `xset.merged.RData` which contains a list of all your datas processed just before.
> {: .solution}
>
> > ### {% icon solution %} 5 - Second XCMS step : *determining shared ions across samples* {% icon tool %}
> >
> >    Open **xcms groupChromPeaks (group)** {% icon tool %} tool and run it with your parameters.
> >    - **Input** : the RData file obtained just before with all datas and peaks informations for each files (format : `xset.merged.RData`).
> >    - **Output** : one Rdata file named `xset.merged.groupChromPeaks.RData` containing all peaks grouped according to their mz and retention time.
> {: .solution}
>
> > ### {% icon solution %} 6 - Optional XCMS step : *retention time correction* {% icon tool %}
> >
> >    This step is optionnal, it aims to correct retention time drift for each peak among samples. You have to use the **xcms adjustRtime (retcor)** {% icon tool%} tool to correct this retention time.
> >    - **Input** : the Rdata file named `xset.merged.groupChromPeaks.RData`
> >    - **Output** : new RData file named `xset.merged.groupChromPeaks.adjustRtime.RData`
> >    
> > > ### {% icon comment %} Important : Have to group again !
> > >
> > > Retention time have been modified. Consequently, it requires to complete it with an additionnal grouping step with **xcms groupChromPeaks (group)** {% icon tool%}. You will obtain a new Rdata file named `xset.merged.groupChromPeaks.adjustRtime.groupChromPeaks.RData`
> >    {: .warning}
> {: .solution}
>
> > ### {% icon solution %} 7 - Final XCMS step *integrating areas of missing peaks* {% icon tool %}
> >
> >    This last step can be run after grouping your peaks and don't need the retention time correction. Run **xcms fillChromPeaks (fillPeaks)** {% icon tool%} with the following things.
> >    - **Input** : you Rdata file `xset.merged.groupChromPeaks.*.RData`
> >    - **Output** : a new RData file named `xset.merged.groupChromPeaks.*.fillChromPeaks.RData`
> {: .solution}
> > ### {% icon comment %} Important : Be careful of the file format
> >
> > During each step of preprocessing, your file has its format changed and can have also its name changed.
> > To be able to continue to MSMS processing, you need to have a RData object wich is **merged and grouped** (step 4 and 5) at least. It means that you should have a file named `xset.merged.groupChromPeaks.RData` (and maybe with some step more in it).
> {: .warning} 
{: .details}

When you have process **all or only needed** steps described before, you can continue with the MS/MS processing part with **msPurity** package. Don't forget to always check your files format ! For the next step you need to have this file `xset.merged.groupChromPeaks.*.RData` where * is the name of **optionnal** steps you could do during the pre-processing.


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

Given a vector of LC-MS/MS or DIMS/MS file paths the **precursor ion purity** of each MS/MS scan can be calculated. It then stores in the purityA S4 class object where **a dataframe of the purity results** can be accessed using the appropriate slot (*pa@puritydf*).

The **calculation** involves dividing the intensity of the selected precursor peak by the total intensity of the isolation window. It is performed before and after the MS/MS scan of interest and interpolated at the recorded time of the MS/MS acquisition.

Additionally, isotopic peaks can estimated and omitted from the calculation. Low abundance peaks are removed that are thought to have limited contribution to the resulting MS/MS spectra. The isolation efficiency of the mass spectrometer can be used to normalise the intensities used for the calculation.

Note that if there are any **mzML files that do not have MS/MS scans** a file ID is save but no assessments will be made.

```R
mzMLpths <- c(file1, file2,...)
pa  <- purityA(mzMLpths)
```

In Galaxy instance, you just need to run the **msPurity.purityA** {% icon tool %} tool with the collection containing all your mzML files. It is the collection with which we started the peak picking processus at the beginning of this tutorial. 

  - **Input** : your **collection** with all your files (`mzML`, `mzXML`,...)
  - **Output** : one RData file containing **pa** object with all MS/MS datas from your files (`msPurity.purityA_on_yourdata.RData`)

> ### {% icon solution %} Advanced parameters
>
> Different parameters are also available depending on the need of the user : 
>
> - **Use most intense peak within isolation window for precursor?** : you can set it to "yes" or "no". If you select "yes", the tool will ignore the recorded precursor within its own file and will use the most intense peak within the isolation window to calculate the precursor ion purity (*--mostIntense* parameter).
> - **Use nearest full scan to determine precursor?** : you can set it to "yes" or "no". If you select "yes", it will use the nearest full scan to the fragmentation scan to determine what is the m/z value of the precursor (*--nearest* parameter).
> - **Interpolation PPM** : have to enter a number. It correspond to the ppm tolerance for the precursor ion purity interpolation (ppm tolerance between the precursor ion found in the neighbour scans). The closest match within the window will be used for the interpolation. (*--ppmInterp* parameter).
> - **offsets** : two options are available here :
>   - **uses offsets determined in the file** : it will use the retention time offsets determined in the file directly for the precursor values
>   - **user supplied offset values** : with this option you have to enter 2 values. One corresponding to the left offsets (the *minoffset*) and the other one corresponding to the right offset (the *maxoffset*).
> - **Threshold to remove peaks below x % of the relative intensity of precursor of interest** : you have to enter a numeric value between 0 and 1. All peaks less than this percentage of the precursor ion of interest will be removed from the purity calculation. It is essentially a noise filter to remove peaks that are thought to have either none or very limited impact on the resulting fragmentation spectra (*--ilim* parameter)
> - **Normalisation for isolation efficiency** : you can choose 4 values for this parameter. If you set it to "none", the intensity of the isolation window stay like this. If you choose one of the following functions (*Gaussian*, *Raised Cosine* or *Calculated from Q-Exactive for +/-0.5 Da windows*), the intensity will be normalised by this function (*--iwNorm* parameter and *--iwNormFun* parameter).
> - **Handling of isotopic peaks** : you can choose 3 propositions (*--isotopes* parameter and *--im* parameter) :
>   - **Keep isotopes in precursor ion purity calculation** : with this option, nothing is done on isotopes in precursor ion purity calculation.
>   - **Exclude C12/C13 isotopes in precursor ion purity calculation** : with this option, you will exclude C13 isotopes (single, double and triple bonds). It is the default parameter.
>   - **Exclude a user supplied list of isotopes in purity calculation** : with this option, you can add a table as input. This table contains all your isotopes you want to excludes. Tabular files should be has the following : 
>
> *Add table*
{: .solution}


> ### {% icon solution %} Go further
> 
> - To verify your **pa** object and to learn what it contains, you can download your RData file from Galaxy (just clic on the *floppy disc* button) and open it in a R terminal or with RStudio application.
> ```R
> library(msPurity)
> load("./Downloads/msPurity.purityA_on_yourdata.RData")
> head(pa@puritydf)
```
> - During this first part of **msPurity** processing, only the *@puritydf* slot has been modified. It is now a table containing the following informations :
>   - pid: unique id for MS/MS scan
>   - fileid: unique id for mzML file
>   - seqNum: scan number
>   - precursorIntensity: precursor intensity value as defined in the mzML file
>   - precursorMZ: precursor m/z value as defined in the mzML file
>   - precursorRT: precursor RT value as defined in the mzML file
>   - precursorScanNum: precursor scan number value as defined in mzML file
>   - id: unique id (redundant)
>   - filename: mzML filename
>   - precursorNearest: MS1 scan nearest to the MS/MS scan purityA35
>   - aMz:  The m/z value in the "precursorNearest" MS1 scan which most closely matches theprecursorMZ value provided from the mzML file
>   - aPurity: The purity score for aMz
>   - apkNm: The number of peaks in the isolation window for aMz
>   - iMz:  The m/z value in the precursorNearest MS1 scan that is the most intense within theisolation window.
>   - iPurity: The purity score for iMz
>   - ipkNm: The number of peaks in the isolation window for iMz
>   - inPurity:  The interpolated purity score (the purity score is calculated at neighbouring MS1scans and interpolated at the point of the MS/MS acquisition)
>   - inpkNm: The interpolated number of peaks in the isolation window
> - The remaining slots for purityA class include :
>   - cores: The number of CPUs to be used for any further processing with this purityA object
>   - fileList: list of the files that have been processed
>   - mzRback: The backend library used by mzR to extract information from the mzML file (e.g. pwiz)
>   - grped_df: If frag4feature has been performed, a dataframe of the grouped XCMS features linked to the associated fragmentation spectra precursor * details is recorded here
>   - grped_MS/MS: If frag4feature has been performed, a list of fragmentation spectra assoicated with each grouped XCMS feature is recorded here
>   - f4f_link_type: If frag4feature has been performed, the linking method is recorded here
>   - av_spectra: if averageIntraFragSpectra, averageInterFragSpectra, or averageAllFragSpectra have been performed, the average spectra is recorded here
>   - av_intra_params: If averageIntraFragSpectra has been performed, the parameters are recorded here
>   - av_inter_params: if averageInterFragSpectra has been performed, the parameters are recorded here
>   - av_all_params: If averageAllFragSpectra has been performed, the parameters are recorded here
>   - db_path: If create_database has been performed, the resulting database path is recorded here
{: .solution}

## frag4feature function

The `frag4feature` function will link the fragmentation data back to the XCMS feature. This function matches precursor ions of MS/MS datas with features and store them into a new table.

```R
pa <- frag4feature(pa, xset)
```
You need to find the **msPurity.frag4feature** {% icon tool %} tool in the Galaxy instance to run this function. It will give you all your MS/MS datas which have a good precursor found during the XCMS processus.

  - Input : You need 2 inputs for this function
    1. Your `xset.merged.groupChromPeaks.*.RData` file containing all the XCMS processus results
    2. The Rdata file obtained just before with all your MS/MS datas (`msPurity.purityA_on_yourdata.RData`)
  - Output : one updated RData file named `msPurity.frag4feature_on_yourdata.RData`

> ### {% icon solution %} Go further
> The slot grped_df is a dataframe of the grouped XCMS features linked to a reference to any associated MS/MS scans in the region of the full width of the XCMS feature in each file. The dataframe contains the following columns.
>
>   - grpid: XCMS grouped feature id
>   - mz: derived from XCMS peaklist
>   - mzmin: derived from XCMS peaklist
>   - mzmax: derived from XCMS peaklist
>   - rt: derived from XCMS peaklist
>   - rtmin: derived from XCMS peaklist
>   - rtmax: derived from XCMS peaklist
>   - into: derived from XCMS peaklist
>   - intb: derived from XCMS peaklist
>   - maxo: derived from XCMS peaklist
>   - sn: derived from XCMS peaklist
>   - sample: derived from XCMS peaklist
>   - id: unique id of MS/MS scan
>   - precurMtchID: Associated nearest precursor scan id (file specific)
>   - precurMtchRT: Associated precursor scan RT
>   - precurMtchMZ: Associated precursor m/z
>   - precurMtchPPM: Associated precursor m/z parts per million (ppm) tolerance to XCMS feauture m/z
>   - inPurity: The interpolated purity score
{: .solution}

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

A query SQLite database can be matched against a library SQLite database with the spectralMatching function. The library spectral-database in most cases should contain the “known” spectra from either public or user generated resources. The library SQLite database by default contains data from MoNA including Massbank, HMDB, LipidBlast and GNPS. A larger database can be downloaded from here.

result <- spectralMatching(q_dbPth, q_xcmsGroups = c(12, 27), cores=1, l_accessions=c('CCMSLIB00000577898','CE000616'))

 Running msPurity spectral matching function for LC-MS(/MS) data

 Filter query dataset

 Filter library dataset

 aligning and matching

 Summarising LC feature annotations

 combineAnnotations function



## createMSP function

Create an MSP file for all the fragmentation spectra that has been linked to an XCMS feature via frag4feature. Can export all the associated scans individually or the averaged fragmentation spectra can be exported.

Additional metadata can be included in a dataframe (each column will be added to metadata of the MSP spectra). The dataframe must contain the column “grpid” corresponding to the XCMS grouped feature.

```
td <- tempdir()
createMSP(pa, msp_file_pth = file.path(td, 'out.msp'))
```

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
