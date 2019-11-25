--- 
layout : tutorial_hands_on

title : 'Mass spectrometry : MSMS analysis with msPurity package'
level : Introductory
zenodo_link : 'https://zenodo.org/record/3244991' 
questions : 
    - What are the main steps of MS/MS data processing for metabolomic analysis ? 
    - How te be able to annotate the maximum of spectra using Galaxy ? 
objectives : 
    - To be sure you have already comprehend the diversity of LC-MS analysis. 
    - To learn the principal functions of msPurity package through Galaxy.
    - To evaluate the potential of this new MS/MS workflow for MS/MS metabolomic analysis. 
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

**Mass spectrometry** (MS) is routinely used to *quantify*, *annotate*, and *identify* small molecules in complex biological matrices. An MS experimental workflow can consist of infusing a sample directly into a mass spectrometer without any prior chromatographic separation (direct infusion mass spectrometry; DIMS). But more often, the sample components are **spatially separated** via either gas (GC) or liquid chromatography (LC). 
?? A METTRE ?? The predictable mass fragmentation patterns observed from the electron ionization (EI) technique commonly used with gas chromatography allows for reliable **matching to mass spectral libraries** such as *NIST3* and the *Golm Metabolome Database*. ?? A METTRE ??
{: .text-justify}

#### Tandem Mass Spectrometry
{: .no_toc}

**Tandem mass spectrometry** (MS/MS or MS2) is a widely used approach for *structural annotation* and *identification of metabolites* in complex biological samples. The term *“tandem mass spectrometry”* is used when a single collision step is used, but product ions can be isolated for further collision to provide MSn spectra where n ≥ 3. A key component of any MS2 (or higher) technology is **the isolation of selected m/z windows** for gas-phase fragmentation and the **mapping back** of the fragmentation (product) spectrum to the selected m/z window.
{: .text-justify}

A metabolomics analytical and data analysis **workflow** that directly takes into consideration the purity of an isolation window has been developed, enabling deconvolution of MS2 spectra. The approach, demonstrated using an Agilent6520 Q-TOF instrument, requires sliding isolation windows to be acquired surrounding the precursor of interest. In these cases, simply **assessing the targeted precursor purity** can be useful in interpreting the MS2 spectra and aid in assessing the reliability of any subsequent annotation. What we call here *“precursor purity”* is calculated with a revised Michalski approach. Advances include that the metric is interpolated at the recorded time of the MS2 acquisition using bordering full scan spectra, the isolation efficiency of the mass spectrometer can be included within the calculation, and as per the PCI (Percent Chimera Intensity) approach, isotopic peaks of the targeted precursor can be removed and peaks with intensities **below a minimum percentage of the precursor peak intensity** are removed from the calculation. The software has been applied to investigate 12 DDA and one DIA metabolomics data sets for different biological samples retrieved from the data repositories [MetaboLights](https://www.ebi.ac.uk/metabolights/), [Metabolomics Workbench](https://www.metabolomicsworkbench.org/), and PRIMeData Resource of Plant Metabolomics ([DROP Met](http://prime.psc.riken.jp/)). We also detail how theoretical isolation windows can be assessed using MS1 data sets collected independent of MS2 acquisitions. The computational methods detailed in the msPurity paper ({% cite Lawson2017 %}) are available in the **R package *msPurity***. The package has been developed to **work as a standalone** or to be used **in conjunction with the metabolomics peak detection and processing R package *XCMS* 2.0**.
{: .text-justify}

#### What will we do ?
{: .no_toc}

To illustrate this approach, we will use data from [*Metabolights n°307*](https://www.ebi.ac.uk/metabolights/MTBLS307), extracted from {% cite VanDerHooft2016%}. The objective of this paper is to detect and visualize antihypertensive drug metabolites in untargeted metabolomics experiments based on the spectral similarity of their fragmentation spectra. To do so, the authors collected urines from a cohort of 26 patients and performed a Thermo Q-Exactive coupled to pHILIC chromatography using data dependent analysis (DDA) MS/MS gas-phase experiments. They found 165 separate drugs metabolites and during this tutorial we will try to find the best result out of these 165. For time reasons, we will just used a subset of the samples and should not find all these waiting metabolites.
{: .text-justify}

To analyze these data, we will start to follow a light version of the [LC-MS workflow](http://workflow4metabolomics.org/the-lc-ms-workflow),
developed by the [Wokflow4metabolomics group](http://workflow4metabolomics.org/) ({% cite Giacomoni2014 %}, {% cite Guitton2017 %}), then we can complet the MS/MS workflow developped by **msPurity** authors ({% cite Lawson2017%}). 
{: .text-justify}

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Preprocessing with XCMS

The first step of the workflow is the pre-processing of the raw data with **XCMS** ({% cite Smith2006 %}).
{: .text-justify}

**XCMS** {% icon tool %} is a free and open source software dedicated to pre-processing of any type of mass spectrometry acquisition files from low to high resolution, including FT-MS data coupled with different kind of chromatography (liquid or gas). This software is used worldwide by a huge community of specialists in metabolomics using mass spectrometry methods.
{: .text-justify}

This software is based on different algorithms that have been published, and is provided and maintained using R software.
{: .text-justify}

**MSnbase readMSData** {% icon tool %} function, prior to **XCMS**, is able to read files with open format as `mzXML`, `mzML`, `mzData` and `netCDF`, which are independent of the constructors' formats. The **XCMS** package itself is composed of R functions able to extract, filter, align and fill gap, with the possibility to annotate isotopes, adducts and fragments (using the R package CAMERA, {% cite CAMERA %}). This set of functions gives modularity, and thus is particularly well adapted to define workflows, one of the key points of Galaxy.
{: .text-justify}

First step of this tutorial is to download the data test. As describe in the introduction, we will use datas from {% cite VanDerHooft2016 %}. We will only process on a subset of these datas. So, you can **import your files directly in Galaxy by using the following URLs below** or download files into your computer (then upload them on Galaxy) : 
{: .text-justify}

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3244991.svg)](https://doi.org/10.5281/zenodo.3244991)
```
https://zenodo.org/record/3244991/files/HU_neg_048.mzML
https://zenodo.org/record/3244991/files/HU_neg_090.mzML
https://zenodo.org/record/3244991/files/HU_neg_123.mzML
https://zenodo.org/record/3244991/files/HU_neg_157.mzML
https://zenodo.org/record/3244991/files/HU_neg_173.mzML
https://zenodo.org/record/3244991/files/HU_neg_192.mzML
https://zenodo.org/record/3244991/files/QC1_002.mzML
https://zenodo.org/record/3244991/files/QC1_008.mzML
https://zenodo.org/record/3244991/files/QC1_014.mzML
```
This step is described as number 1 in details part below.

Then, to be able to process your MS/MS datas, we need to **start with the peakpicking of MS datas**. One Galaxy Training already explains how to process with your MS datas. You should **follow this link and complete this tutorial** : [Mass spectrometry: LC-MS analysis](https://galaxyproject.github.io/training-material/topics/metabolomics/tutorials/lcms/tutorial.html). For MS/MS analysis you **don't really need to finish** this previous tutorial but for a better understanding of your datas, it is recommanded. In this tutorial, you **just have to compute** your datas with the **following steps** briefly describe in the *details* part below (please follow aprameters values to have the same results during the training).
{: .text-justify}

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
> > **MSnbase readMSData** {% icon tool %}, prior to **XCMS** {% icon tool %} is able to read files with open format as `mzXML`, `mzML`, `mzData` and `netCDF`, which are independent of the constructors' formats.
> >   - **Input** : your collection with all your files
> >   - **Output** : a new collection with your files and their information after `readMSData` function. Format : `collection.raw.RData`.
> {: .solution}
>
> > ### {% icon solution %} 3 - First XCMS step : *peak picking* {% icon tool %}
> >    
> > Open and run the **xcms findChromPeaks (xcmsSet)** {% icon tool %} tool and set your parameters. For our training please enter the following parameters : 
> >   - **Extraction method for peaks detection** : there are 4 different possible options here. Please select `CentWave - chromatographic peak detection using the centWave method` for our training.
> >   - **Max tolerated ppm m/z deviation in consecutive scans in ppm** : it corresponds to the maximum deviation in ppm. Please set it to `5` for our training.
> >   - **Min, Max peak width in seconds** : it corresponds to the expected approximate peak width in chromatographic space. Please set it to `5,50` in our training to have a large example. 
> > 
> > Here, you have all the right parameters for our example. These are the two dataset collection you will have in your hsitory on the right of Galaxy instance : 
> >    - **Input** : collection of Rdata object(s) obtained just before (format : `collection.raw.RData`)
> >    - **Output** : a collection with informations about peak picking for each files. Format : `collection.raw.xset.RData`.
> {: .solution}
> 
> > ### {% icon solution %} 4 - *Merge samples* {% icon tool %} in one dataset with a sampleMetadata file
> >
> > To merge your datas, you need to **input a sampleMetadata file** containing namefiles and their metadata informations like their class for example. If you don't add a sampleMetadata file, the **xcms findChromPeaks Merger** {% icon tool %} tool will **group all your files together**. You can also **create your sampleMetadata file** with W4M Galaxy tool **xcms get a sampleMetadata file** {% icon tool %} with the following parameters: *"RData file"* outputed from **MSnbase readMSData** {% icon tool %}. Here is an example of the minimum expectations about a sampleMetadata file (important : don't write the format of the file, just their names) :
> > {: .text-justify}
> > 
> > | sample_name |  class  | 
> > |:-----------:|:-------:|
> > |    file1    |  homme  |
> > |-------------+---------|
> > |    file2    |  femme  |
> > |-------------+---------|
> > |    file3    |  homme  |
> > 
> > You have just to enter the following files in the tool and process it. You should now have these files in your history :
> >   - **Input** : 
> >     - collection of RData object(s) obtained during peak picking (format : `collection.raw.xset.RData`) 
> >     - optionnaly a sampleMetadata file in `tabular` format. For our example please start with `sampleMetadata_hommeVSfemme.csv` in tabular format.
> >   - **Output** : 
> >     - one RData file named `xset.merged.RData` which contains a list of all your datas processed just before.
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
{: .text-justify}

# Processing with msPurity package

**msPurity** is a R package developped by Birmingham University team and published in 2017 ({% cite Lawson2017 %}). It is also available via Galaxy instance through some [wrappers](https://github.com/computational-metabolomics/mspurity-galaxy/tree/master) that were developped by the same people. 
{: .text-justify}

This R package was developped to :
1. Assess the spectral quality of fragmentation spectra by evaluating the "precursor ion purity". 
2. Process fragmentation spectra. 
3. Perform spectral matching. 

**What is *precursor ion purity*?** 

What we call **"precursor ion purity"** is a measure of the contribution of a selected precursor peak in an isolation window used for fragmentation. The simple calculation involves dividing the intensity of the selected precursor peak by the total intensity of the isolation window. When assessing MS/MS spectra this calculation is done before and after the MS/MS scan of interest and the **purity is interpolated at the recorded time of the MS/MS acquisition**. Additionally, isotopic peaks can be removed, low abundance peaks are removed that are thought to have limited contribution to the resulting MS/MS spectra and the isolation efficiency of the mass spectrometer can be used to normalise the intensities used for the calculation.
{: .text-justify}

## 1 - Assessing the purity (*purityA function*)

The importance of **assessing the contribution of the precursor ion** within an isolation window for MS2 experiments has been previously detailed in proteomics, where precursor ion purity influences the quality and accuracy of matching to mass spectral libraries. But to date, there has been little attention to this data-processing technique in metabolomics. In **targeted and data-dependent acquisition** (DDA)-based experiments, where MS2 is performed on a dynamic list of precursor ions (as often determined by a preceding MS full scan), an isolation window is centered on the targeted precursor peak (m/z value). However, the **isolation window** can often contain **more than one distinct peak**. Fragmentation spectra resulting from these situations being termed *“chimeric”* and can be problematic for interpretation of the spectra and mass spectral library searching.
{: .text-justify}

The precursor purity metric is calculated as *“intensity of a selected precursor divided by the summed intensity of the isolation window”*. The impact of chimeric spectra on spectral matching and annotation depends on the purity of the isolation window fragmented (i.e., the ratio between the relative intensity of the precursor divided by the summed intensity of all ions within the isolation window). If the purity of the precursor ion is sufficiently low, it can often be difficult to determine the origin of the resulting product ion(s). This in turn can lead to erroneous spectral matching results or no spectral matches. Deconvolution of chimeric spectra however is possible and forms the basis of the data-analysis procedures applied to data independent acquisition (DIA) experiments.
{: .text-justify}

> ### {% icon comment %} How the purityA algorithm works ?
> 
> Given a vector of LC-MS/MS file paths the **precursor ion purity** of each MS/MS scan from each file can be calculated. It then stores in the purityA S4 class object where **a dataframe of the purity results** can be accessed using the appropriate slot (*pa@puritydf*). The **calculation** involves dividing the intensity of the selected precursor peak by the total intensity of the isolation window. It is performed before and after the MS/MS scan of interest and interpolated at the recorded time of the MS/MS acquisition. 
> {: .text-justify}
> 
> 
> Additionally, isotopic peaks can estimated and **omitted** from the calculation. Low abundance peaks are **removed** that are thought to have limited contribution to the resulting MS/MS spectra. The isolation efficiency of the mass spectrometer can be used to **normalise** the intensities used for the calculation. The **precursor ion purity** represents the measure of the contribution of a selected precursor peak in an isolation window used for fragmentation and can be used as away of assessing the spectral quality and level of "contamination" of fragmentation spectra. 
> {: .text-justify}
> 
> 
> Note that if there are any **files that do not have MS/MS scans** a file ID is save but no assessments will be made. 
> {: .text-justify}
> 
{: .comment}

Let's try assessing the purity of precursors with the Galaxy tool **msPurity.purityA {% icon tool %}**.

> ### {% icon hands_on %} Hands-on : msPurity.purityA {% icon tool %}
> 
> In Galaxy instance, you just need to run the **msPurity.purityA** {% icon tool %} tool with the collection containing all your files. It is the collection with which we started the peak picking processus at the beginning of this tutorial. 
>  - **\*.mzML file** : Please select the directory icon and then `select your dataset collection` you used at the begining of your XCMS preprocessing
>  - **Use most intense peak within isolation window for precursor?** : you can set it to "No" value if you want the process to keep the registered precursor from the file. If you set it to "Yes", it will selected the most intense peak within its isolation window. Please `set it to True`.
>  - **Use nearest full scan to determine precursor?** : you can set it to "No" value to keep the precursor scan already registered within the file. If you set it to "True", it will use the nearest full scan to determine what the m/z value is of the precursor. Please `set it to True`.
>  - **Interpolation PPM** : Enter the value you want as the ppm tolerance for the precursor ion purity interpolation (the closest match within the window will be used for the interpolation). Please `set it to 5`.
>
> You just have to adjust these 4 parameters and don't touch the others.
> 
> > ### {% icon comment %} Resume
> > To resume this tool, you will have these things in your history on the right side of Galaxy instance : 
> >  - **Input** : 
> >    - your **collection** with all your files (`mzML`, `mzXML`,...)
> >  - **Output** : 
> >    - one RData file containing **pa** object with all MS/MS datas from your files (`msPurity.purityA_on_yourdata.RData`)
> >    - one `tsv` file with all results from the assess of purity in the files you gave (*slot puritydf of pa object*). Each line correspond to one MS/MS scan (`msPurity.purityA_on_yourdata.tsv`). 
> {: .comment}
{: .hands_on}

> ### {% icon solution %} Advanced parameters
>
> Here, we list all parameters available on Galaxy tool, depending of the need of the user and not cited just before : 
>
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
> 
> | isotope_id | mass diff | abundance of isotope  | ppm tol for mz  | abundance buffer | charge | relative atomic mass (int) | xflag |
> |:----------:|:---------:|:---------------------:|:---------------:|:----------------:|:------:|:--------------------------:|:-----:|
> |     1      | 1.003355  |        1.07           |        4        |       0.1        |    1   |             12             |   1   |
> |------------+-----------+-----------------------+-----------------+------------------+--------+----------------------------+-------|
> |     2      | 1.9659026 |        0.2424         |        4        |       0.1        |    1   |             35             |   1   |
> |------------+-----------+-----------------------+-----------------+------------------+--------+----------------------------+-------|
> |     3      | 0.916291  |        0.4931         |        4        |       0.1        |    1   |             80             |   1   |
> {: .custom-class #table_isotopes}
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
>   - fileid: unique id for file
>   - seqNum: scan number
>   - precursorIntensity: precursor intensity value as defined in the file
>   - precursorMZ: precursor m/z value as defined in the file
>   - precursorRT: precursor RT value as defined in the file
>   - precursorScanNum: precursor scan number value as defined in file
>   - id: unique id (redundant)
>   - filename: filename
>   - precursorNearest: MS1 scan nearest to the MS/MS scan purityA35
>   - aMz:  The m/z value in the "precursorNearest" MS1 scan which most closely matches theprecursorMZ value provided from the file
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
>   - mzRback: The backend library used by mzR to extract information from the file (e.g. pwiz)
>   - grped_df: If frag4feature has been performed, a dataframe of the grouped XCMS features linked to the associated fragmentation spectra precursor * details is recorded here
>   - grped_MS/MS: If frag4feature has been performed, a list of fragmentation spectra assoicated with each grouped XCMS feature is recorded here
>   - f4f_link_type: If frag4feature has been performed, the linking method is recorded here
>   - av_spectra: if averageIntraFragSpectra, averageInterFragSpectra, or averageAllFragSpectra have been performed, the average spectra is recorded here
>   - av_intra_params: If averageIntraFragSpectra has been performed, the parameters are recorded here
>   - av_inter_params: if averageInterFragSpectra has been performed, the parameters are recorded here
>   - av_all_params: If averageAllFragSpectra has been performed, the parameters are recorded here
>   - db_path: If create_database has been performed, the resulting database path is recorded here
{: .solution}

The output is an RData file with the purityA S4 class object (referred to as pa for convenience throughout the manual). The object contains a slot (pa@puritydf) where the details of the purity assessments for each MS/MS scan. The purityA object can then be used for further processing including linking the fragmentation spectra to XCMS features, averaging fragmentation, database creation and spectral matching (from the created database).
There is also the additional output of the a tsv file of the pa@puritydf data frame.

## 2 - Match features with fragmentation spectra (*frag4feature function*)

The `frag4feature` function will link the fragmentation spectra (MS/MS) back to the XCMS features. This function matches precursor ions of MS/MS datas with features and store them into a new table. To be able to be matched, the associated acquisition time of the MS/MS event has to be **within the retention time window** defined for the individual peaks associated and the precursor m/z value also has to be **within the user ppm tolerance** to XCMS feature. 

> ### {% icon hands_on %} Hands-on : msPurity.frag4feature {% icon tool %}
>
> Just find the **msPurity.frag4feature** {% icon tool %} tool in the Galaxy instance to be able to run this function. As result, it will give you all your MS/MS datas which have a good precursor found during the XCMS processus. But before, you have to enter the right parameters : 
> {: .text-justify}
>  - **xcmsSet object** : the Rdata file you made during the XCMS preprocessing. Its name should be `xset.merged.groupChromPeaks.*.RData` (where \* is the name of optionnal tools you ran).
>  - **purityA object** : the RData file ouputed the step before containing all your MS/MS spectra. Should be named `msPurity.purityA_on_yourdata.RData`.
>  - **ppm error tolerance between precursor mz and XCMS feature mz** : fragmentation will be ignored if the precursor m/z value is not within the ppm tolerance to the XCMS feature m/z. Please **set it to 10** (default).
>  - **Should the most intense precursor be used within the isolation window?** : if you set it to "yes", the most intense precursor will be used. If it is on "no" the precursor closest the center of the isolation window will be used. Please **set it to "yes"**.
>  - **Was retention time correction used?** : just **precise if you used retention time correction** when you made the XCMS preprocessing. It is an optionnal step so maybe you have not made it.
>  - **For matching fragmentation to a feature, use the grouped feature range** : XCMS has already calculated individual peaks for each file, then it has grouped them together. If your MS/MS files have also MS datas, you can link MS2 spectra directly to a peak for each file. In this case, you can turn this parameter to "no". If you don't have any MS datas in each file, you can turn it to "yes". In this case, the full width of the grouped peaks will be used to link MS2 spectra with it. In our case, please **set it to "no"**, we have enough MS datas in our files.
> 
> You have an other parameter that is a threshold for the precursor ion purity. It is here to filter precursor with a high purity, but we will see the filter tool just after to be able to order our datas. 
> {: .text-justify}
>
> > ### {% icon comment %} Resume
> > To resume this tool, you will have these things in your history on the right side of Galaxy instance :
> >  - Input : You need 2 inputs for this function
> >    - Your `xset.merged.groupChromPeaks.*.RData` file containing all the XCMS processus results.
> >    - The Rdata file obtained just before with all your MS/MS datas (`msPurity.purityA_on_yourdata.RData`).
> >  - Output : 
> >    - one updated RData file named `msPurity.frag4feature_on_yourdata.RData`.
> >    - one `tsv` file containg all MS/MS spectra that matched with a XCMS feature (*slot grped_df of pa object*). Each line has informations about  the grouped XCMS features linked to the associated fragmentation spectra precursor details (`msPurity.frag4feature_on_yourdata.tsv`).
> {: .comment}
{: .hands_on}

> ### {% icon solution %} Go further
> The slot grped_df is a dataframe of the grouped XCMS features linked to a reference to any associated MS/MS scans in the region of the full width of the XCMS feature in each file. The dataframe contains the following columns.
> {: .text-justify}
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
>   - precurMtchPPM: Associated precursor m/z parts per million (ppm) tolerance to XCMS feature m/z
>   - inPurity: The interpolated purity score
>
> The slot grped_MS2 is a list of the associated fragmentation spectra for the grouped features.
> {: .text-justify}
{: .solution}

## 3 - (Optionnal) Filter your results (*filterFragSpectra function*)

The fragmentation result can be filtered prior to averaging using the `filterFragSpectra` function. You can filter your features according to different parameters : 
 - **precursor ion purity** : you can filter by enter the value that will be the threshold for the precursor ion purity value. Each MS/MS spectra with a precursor ion purity under this threshold will be tagged as *false* or will be removed.
 - **peak intensity** : you can fix a minimum for the peak intensity. All peaks under this value will be tagged or removed.
 - **relative abundance** : the relative abundance is calculated for each peak in the window. You can select a value that will be the minimum requiered as relative abundance. Every peak under this value will be removed or flagged. 
 - **signal-to-noise** : the ratio of a peak signal to the noise is calculated for each peak of a file. You can enter the minimum signa-to-noise ratio that will removed (or flagged) all peaks under this value. 

> ### {% icon hands_on %} Hands-on : msPurity.filterFragSpectra {% icon tool %}
>
> This function is used in **msPurity.filterFragSpectra {% icon tool %}** from Galaxy instance. It is an optionnal step of the process but we can make some tries to refine our MS/MS spectra and have a little bit less than thousands we had. Here is an example with the different parameters you have to complete to make a selection of your matched MS/MS spectra : 
> 
>  - **Miniumum precursor ion purity of the associated precursor for fragmentation spectra scan** : minimum precursor ion purity of the associated precursor for fragmentation spectra scan. For example, you can **set it to 0.9** to obtain MS/MS spectra with a good porecursor ion purity only.
>  - **Peak instensity threshold** : minimum intensity of a peak. 
>  - **Relative abundance threshold** : minimum relative abundance of a peak.
>  - **snr** : minimum signal-to-noise of a peak within each file. 
>  - **rmp** : select if you want to remved peak that do not meet the threshold criteria.
>  - **snmeth** : method to calculate signal to noise ration (median or mean).
>  - **allfrag** : choose to filter on all fragmentation spectra or only the fragmentation spectra grouped to XCMS features.
>
>
{: .hands_on}

## 4 - (Optionnal) Average your results (*averageFragSpectra function*)

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

## 5 - (Optionnal) Create a database (*creataDatabase function*)

An SQLite database is then created of the LC-MS/MS experiment. The SQLite schema of the spectral database can be detailed here.

```R
td <- tempdir()
q_dbPth <- createDatabase(pa, xset, outDir = td, dbName = 'lcmsms-processing.sqlite')
```

## 6 - (Optionnal) (*spectralMatching function*)

A query SQLite database can be matched against a library SQLite database with the spectralMatching function. The library spectral-database in most cases should contain the “known” spectra from either public or user generated resources. The library SQLite database by default contains data from MoNA including Massbank, HMDB, LipidBlast and GNPS. A larger database can be downloaded from here.

result <- spectralMatching(q_dbPth, q_xcmsGroups = c(12, 27), cores=1, l_accessions=c('CCMSLIB00000577898','CE000616'))

 Running msPurity spectral matching function for LC-MS(/MS) data

 Filter query dataset

 Filter library dataset

 aligning and matching

 Summarising LC feature annotations

 combineAnnotations function



## 7 - Create your MSP file (*createMSP function*)

Create an MSP file for all the fragmentation spectra that has been linked to an XCMS feature via frag4feature. Can export all the associated scans individually or the averaged fragmentation spectra can be exported.

Additional metadata can be included in a dataframe (each column will be added to metadata of the MSP spectra). The dataframe must contain the column “grpid” corresponding to the XCMS grouped feature.

```
td <- tempdir()
createMSP(pa, msp_file_pth = file.path(td, 'out.msp'))
```

## 8 - (Optionnal) (*flagRemove function*)

blabla


# Annotation

## With MetFrag

## With Sirius CSI-FingerID

## With .??


# Conclusion 
{:.no_toc} 

blabla 

>    > ### {% icon comment %} Comment
>    >
>    > Here we provided the sampleMetadata file so we know that the upload led to a 'tabular' file. But from experience we also know that it can happen that, when uploading a sampleMetadata table, user obtained other inappropriate types of data. This is generally due to the file not following all the requirements about the format (*e.g.* wrong separator, or lines with different numbers of columns).
>    > Thus, we highly recommend that you always take a second to check the data type after the upload. This way you can handle the problem right away if you appear to get one of these obvious issues.
>    {: .comment}
