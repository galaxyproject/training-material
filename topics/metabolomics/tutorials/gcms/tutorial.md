--- 
layout : tutorial_hands_on

title : 'Mass spectrometry : GC-MS analysis with metaMS package'
level : Introductory
zenodo_link : 'https://zenodo.org/record/3244991' 
questions : 
    - What are the main steps of GC-MS datas processing for metabolomic analysis ? 
    - How te be able to annotate the maximum of unknowns using Galaxy ? 
objectives : 
    - To be sure you have already comprehend the diversity of MS preprocessing analysis. 
    - To learn the principal functions of metaMS package through Galaxy.
    - To evaluate the potential of this new GC-MS workflow for GC-MS metabolomic analysis. 
time_estimation : 2H 
key_points : 
    - Have a good file containing all your peaks during the first stopover 
    - Find all your unknowns in your datas
    - Find your stanards if you have some 
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

Then, to be able to process your MS/MS datas, we need to **start with the peakpicking of MS datas**. One Galaxy Training already explains how to process with your MS datas. You should **follow this link and complete this tutorial** : [Mass spectrometry: LC-MS analysis](https://galaxyproject.github.io/training-material/topics/metabolomics/tutorials/lcms/tutorial.html). For MS/MS analysis you **don't really need to finish** this previous tutorial but for a better understanding of your datas, it is recommanded. In this tutorial, you **just have to compute** your datas with the **following steps** briefly describe in the *details* part below (please follow parameters values to have the same results during the training).
{: .text-justify}

> ### {% icon details %} Some help : Resume of the XCMS preprocessing
>
> Here is a resume of the **XCMS** {% icon tool %} preprocessing. Please follow instructions for this training tutorial. It has few steps that are an obligation to be able to obtain the final file. With this file, we can continue the workflow to process our MS/MS datas with **msPurity** {% icon tool %} package.
> {: .text-justify}
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
> > ### {% icon solution %} 3 - First XCMS step : *peak picking*
> >    
> > Open and run the **xcms findChromPeaks (xcmsSet)** {% icon tool %} tool and set your parameters. For our training please enter the following parameters : 
> >   - **Extraction method for peaks detection** : there are 4 different possible options here. Please select `CentWave - chromatographic peak detection using the centWave method` for our training.
> >   - **Max tolerated ppm m/z deviation in consecutive scans in ppm** : it corresponds to the maximum deviation in ppm. Please set it to `5` for our training.
> >   - **Min, Max peak width in seconds** : it corresponds to the expected approximate peak width in chromatographic space. Please set it to `5,20` in our training to have a large example. 
> > 
> > Here, you have all the right parameters for our example. These are the two dataset collection you will have in your hsitory on the right of Galaxy instance : 
> >    - **Input** : collection of Rdata object(s) obtained just before (format : `collection.raw.RData`)
> >    - **Output** : a collection with informations about peak picking for each files. Format : `collection.raw.xset.RData`.
> {: .solution}
> 
> > ### {% icon solution %} 4 - *Merge samples* in one dataset with a sampleMetadata file
> >
> > To merge your datas, you need to **input a sampleMetadata file** containing namefiles and their metadata informations like their class for example. If you don't add a sampleMetadata file, the **xcms findChromPeaks Merger** {% icon tool %} tool will **group all your files together**. You can also **create your sampleMetadata file** with W4M Galaxy tool **xcms get a sampleMetadata file** {% icon tool %} with the following parameters: *"RData file"* outputed from **MSnbase readMSData** {% icon tool %}. Here is an example of the minimum expectations about a sampleMetadata file (**important** : don't write the format of the file, just their names) :
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
> > ### {% icon solution %} 5 - Second XCMS step : *determining shared ions across samples*
> >
> > Open **xcms groupChromPeaks (group)** {% icon tool %} tool and run it with your parameters. This tool will group all your peaks found before together when they represent the same analyte across samples. It uses overlapping m/z bins and calculations of smoothed peak distributions in chromatographic time. In Galaxy instance please set the parameters as below : 
> >  - **Method to use for grouping** : please `keep the PeakDensity method` which will groups peak based on time dimension peak densities.
> >  - **Bandwidth** : it correspond to the standard deviation of the smoothing kernel to be used. Set it to `10` for our example.
> >  - **Minimum fraction of samples** : minimum fraction of samples in at least one sample group in which the peaks have to be present to be considered as a peak group. This value is the threshold and `don't change` it for our example.
> >  - **Minimum number of samples** : minimum number of sample(s) in which the peak have to be detected to be considered as a group. Please `don't change` it for our example.
> >  - **Width of overlapping m/z slices** : the m/z difference between two peaks to be able to be grouped. For our example please set it to `0.05`.
> > 
> > You don't need to change the other parameters, except if you want the peaklist as output. You should now have these files in your history : 
> >  - **Input** : 
> >    - the RData file obtained just before with all datas and peaks informations for each files (format : `xset.merged.RData`).
> >  - **Output** : 
> >    - one Rdata file named `xset.merged.groupChromPeaks.RData` containing all peaks grouped according to their mz and retention time.
> >    - one `pdf` file named `xset.merged.groupChromPeaks.plotChromPeakDensity.pdf`. It contains chromatographic density plots.
> >    - 2 `tsv` files. The first one named `xset.merged.group.dataMatrix.tsv` contains informations about ions intensities. The second one named `xset.merged.group.variableMetadata.tsv` contains an other table with informations about ions. These files are generated only if you set the parameter **Get peaklist** to `yes`. It is just more informations after the grouping step. 
> {: .solution}
>
> > ### {% icon solution %} 6 - Optional XCMS step : *retention time correction*
> >
> > This step is optionnal, it aims to correct retention time drift for each peak among samples. You have to use the **xcms adjustRtime (retcor)** {% icon tool%} tool to correct this retention time. We `don't use it in our example`. But if you will use it during an other processing, you should have these files in your history : 
> >  - **Input** : 
> >    - the Rdata file already with peaks grouped and named `xset.merged.groupChromPeaks.RData`
> >  - **Output** : 
> >    - new RData file named `xset.merged.groupChromPeaks.adjustRtime.RData` with correction of retention time
> >    - one `pdf` file named `xset.merged.groupChromPeaks_rawVSadjusted.adjustRtime.Rplots.pdf`  and containing a plot to show the difference between rt and adjusted rt for each file.
> >    
> > > ### {% icon comment %} Important : Have to group again !
> > >
> > > Retention time have been modified. Consequently, it requires to complete it with an additionnal grouping step with **xcms groupChromPeaks (group)** {% icon tool%}. You will obtain a new Rdata file named `xset.merged.groupChromPeaks.adjustRtime.groupChromPeaks.RData`
> >    {: .warning}
> {: .solution}
>
> > ### {% icon solution %} 7 - Final XCMS step *integrating areas of missing peaks*
> >
> > This last step can be run after grouping your peaks and don't need the retention time correction. Run **xcms fillChromPeaks (fillPeaks)** {% icon tool%} to identify, for each sample, peak groups where that sample is not represented. 
> >  - **Input** : 
> >    - your Rdata file `xset.merged.groupChromPeaks.*.RData`
> >  - **Output** : 
> >    - a new RData file named `xset.merged.groupChromPeaks.*.fillChromPeaks.RData`
> {: .solution}
> > ### {% icon comment %} Important : Be careful of the file format
> >
> > During each step of preprocessing, your file has its format changed and can have also its name changed.
> > To be able to continue to MSMS processing, you need to have a RData object wich is **merged and grouped** (step 4 and 5) at least. It means that you should have a file named `xset.merged.groupChromPeaks.RData` (and maybe with some step more in it).
> {: .warning} 
{: .details}


# Stopover : Verify your datas after the XCMS preprocessing

When you have process **all or only needed** steps described before, you can continue with the MS/MS processing part with **msPurity** package. Don't forget to always check your files format ! For the next step you need to have this file `xset.merged.groupChromPeaks.*.RData` where * is the name of **optionnal** steps you could do during the pre-processing. For our example, your file should be named `xset.merged.groupchromPeaks.RData`. 
{: .text-justify}

> ### {% icon comment %} Comment
> 
> The preprocessing part of this analysis can be **quite time-consuming**, and already corresponds to quite a few number of steps, depending of your analysis. We highly recommend, at this step of the MS/MS workflow, to split your analysis by beginning a new Galaxy history with **only the files you need** (final xset Rdata file and your data collection of mzML). This will help you in limiting selecting the wrong dataset in further analysis, and bring a little **tidiness** for future review of your MS/MS analysis process. You should also be able to make a better peakpicking in the future in the same history and it will not be polluated by MS/MS part of your process.
> {: .text-justify}
> 
> > ### {% icon tip %} Tip: Copy dataset to a new history
> >
> > 1. Click on the {% icon galaxy-gear %} icon (**History options**) on the top of the history panel
> > 2. Click on **Copy Dataset**
> > 3. Select the desired files
> > 4. Give a relevant name to the "New history"
> > 5. Click on the new history name in the green box that have just appear to switch to this history
> {: .tip}
> 
> To begin a new history with the files from your current history, you can **use the functionality ‘copy dataset’** and copy it into a new history (the option is hidden behind the notched wheel at the top right of the history).
> {: .text-justify}
> 
> You may have notice that the XCMS tools generate **output names that contain the different XCMS steps you used**, allowing easy traceability while browsing your history. Hence, we highly recommend you to rename it **with something short**, e.g. "xset", "XCMSSetObject", or anything not too long that you may find convenient.
> {: .text-justify}
> {% include snippets/rename_dataset.md %}
>
{: .comment}

Before the next step with msPurity package on MS/MS datas, here are some questions to be able to verify if your file is ready and if you have the same results as us. Please check these questions : 
{: .text-justify}

> ### {% icon question %} Question before MS/MS steps
> 
>  **1** - What are the steps of XCMS you made before your final file ?
> > ### {% icon solution %} Solution
> > 
> > Here are the different steps made for our example : 
> >  - **(Not with XCMS)** import your datas into Galaxy instance
> >  - **MSNbase readMSData** {% icon tool %} to read our MS datas
> >  - XCMS peakpicking with **xcms findChromPeaks (xcmsSet)** {% icon tool %} tool
> >  - (Not with XCMS but necessary) merge my datas into one file with **xcms findChromPeaks Merger** {% icon tool %} tool
> >  - XCMS grouping with **xcms groupChromPeaks (group)** {% icon tool %} tool
> >  - **(Not done)** XCMS retention time correction, then grouping again with xcms adjustRtime (retcor) {% icon tool %} tool
> >  - XCMS integration of missing peaks with **xcms fillChromPeaks (fillPeaks)** {% icon tool %} tool
> > 
> {: .solution}
> <br>
>  **2** - Concerning what we said before and the previous answer, what is the complete name of your final RData file ?
> > ### {% icon solution %} Solution
> > 
> > During each step of XCMS preprocessing, the name of the file which is processing is completed by the name of the step you were doing. So, finally your file should be name `xset.merged.groupChromPeaks.fillChromPeaks.RData`. That because (as seen in previous answer) you ran a grouping and the integration after merged datas.
> > {: .text-justify}
> > 
> {: .solution}
> <br>
> **3** - What is the size (in MB) of your final RData file ?
> > ### {% icon solution %} Solution
> > 
> > To be able to see the size of a file in your history, you just have to select it. It will deployed informations about it and you can see the size of yours. For our example, the size of the final file is **1.4 MB**.
> > {: .text-justify}
> > 
> {: .solution}
{: .question}


# Stopover : Verify your datas after the XCMS preprocessing



# Processing with msPurity package


## 1 - Assessing the purity (*purityA function*)


> ### {% icon hands_on %} Hands-on : msPurity.purityA {% icon tool %}
> 
> 
> > ### {% icon comment %} Resume
> > 
> {: .comment}
{: .hands_on}

> ### {% icon question %} Question
> 
> 
> > ### {% icon solution %} Solution
> > 
> > 
> {: .solution}
> 
{: .question}


> ### {% icon details %} Advanced parameters
> 
{: .details}

## 2 - Match features with fragmentation spectra (*frag4feature function*)



> ### {% icon hands_on %} Hands-on : msPurity.frag4feature {% icon tool %}
>
>
> > ### {% icon comment %} Resume
> > 
> {: .comment}
{: .hands_on}

> ### {% icon question %} Question
> 
> 
> 
> > ### {% icon solution %} Solution
> > 
> {: .solution}
> 
{: .question}

> ### {% icon solution %} Go further
> 
{: .solution}

## 3 - (Optionnal) Filter your results (*filterFragSpectra function*)


> ### {% icon hands_on %} Hands-on : msPurity.filterFragSpectra {% icon tool %}
> 
>
{: .hands_on}

## 4 - (Optionnal) Average your results (*averageFragSpectra function*)


> ### {% icon details %} How to process the averaging ?
> 
> 
{: .details}

## 5 - (Optionnal) Create a database (*creataDatabase function*)



## 6 - (Optionnal) (*spectralMatching function*)



## 7 - Create your MSP file (*createMSP function*)



## 8 - (Optionnal) (*flagRemove function*)

blabla

# Stopover : Verify your datas after the msPurity processing


# Annotation

## With MetFrag


## With Sirius CSI-FingerID



# Conclusion 
{:.no_toc} 

During this tutorial, the user learnt how to process with MS/MS datas. The peakpicking part is very very important for the rest of the annotation of MS/MS datas. 
