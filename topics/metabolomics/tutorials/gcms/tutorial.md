---
layout: tutorial_hands_on
draft: true

title: 'Mass spectrometry : GC-MS analysis with metaMS package'
level: Introductory
zenodo_link : 'https://zenodo.org/record/3244991' 
questions : 
    - What are the main steps of GC-MS datas processing for metabolomic analysis ? 
    - How te be able to annotate the maximum of unknowns using Galaxy ? 
objectives : 
    - To be sure you have already comprehend the diversity of MS pre-processing analysis. 
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


A lot of packages are available for the analysis of GC-MS or LC-MS data. Typically, hardware vendors provide software that is optimized for the instrument and allow a direct interaction of the lab scientist with the data. Some other open-source alternatives such as **XCMS** are also able to be integrated easily in web interfaces, allowing large numbers of files to be processed simultaneously. Because of the generality of packages like **XCMS**, several other packages have been developed to use the functionality of **XCMS** for optimal performance in a particular context. Package **metaMS** does so for the field of untargeted metabolomics, focuses on the GC-MS analysis during this tutorial. One of the goals **metaMS** was to set up a simple system with few user-settable parameters, capable of handling the vast majority of untargeted metabolomics experiments. 

During this tutorial, we will learn how to process easily a test dataset from raw files to the annotation using W4M Galaxy. Datas are from {% cite Dittami2012 %} and have been used as test dataset for the development of the Galaxy wrappers. 


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# First steps of pre-processing using a standard XCMS workflow (mandatory)

The first step of the workflow is the pre-processing of the raw data with **XCMS** ({% cite Smith2006 %}).
{: .text-justify}

**XCMS** {% icon tool %} is a free and open-source software dedicated to pre-processing any type of mass spectrometry acquisition files from low to high resolution, including FT-MS data coupled with different kind of chromatography (liquid or gas). This software is used worldwide by a huge community of specialists in metabolomics.
{: .text-justify}

This software is based on different algorithms that have been published, and is provided and maintained using R software.
{: .text-justify}

**MSnbase readMSData** {% icon tool %} function, prior to **XCMS**, is able to read files with open format as `mzXML`, `mzML`, `mzData` and `netCDF`, which are independent of the constructors' formats. The **XCMS** package itself is composed of R functions able to extract, filter, align and fill gap, with the possibility to annotate isotopes, adducts and fragments (using the R package CAMERA, {% cite CAMERA %}). This set of functions gives modularity, and thus is particularly well adapted to define workflows, one of the key points of Galaxy.
{: .text-justify}

First step of this tutorial is to download the data test. As describe in the introduction, we will use datas from {% cite Dittami2012 %}. We will only process on a subset of their data. 
So, you can **import your files directly in Galaxy from Zenodo (see hands-on below)** or download files into your computer using the following link then upload them on Galaxy: 
{: .text-justify}

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3631074.svg)](https://doi.org/10.5281/zenodo.3631074)


Then, to be able to pre-process our GC-MS data, we need to **start with the peakpicking of MS data**. 
One Galaxy Training material already explains how to act with MS data. We encourage you to **follow this link and complete the corresponding tutorial**: [Mass spectrometry: LC-MS preprocessing with XCMS]({% link topics/metabolomics/tutorials/lcms-preprocessing/tutorial.md %}). 
For GC-MS analysis you **don't really need to follow all of this previous tutorial** but for a better understanding of your data, it is recommanded to try it with their test dataset.
Concerning the current GC-MS tutorial, you **just have to compute the following steps and specific parameters** described in the hands-on part below (please follow parameters values to have the same results during the training).
{: .text-justify}

## 1 - Import the data into Galaxy 

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the following 6 `mzData` files into a collection named `mzData`
>    - Option 1: from a shared data library (ask your instructor)
>    - Option 2: from Zenodo using the URLs given below
>
>    ```
>    https://zenodo.org/record/3631074/files/alg11.mzData
>    https://zenodo.org/record/3631074/files/alg2.mzData
>    https://zenodo.org/record/3631074/files/alg3.mzData
>    https://zenodo.org/record/3631074/files/alg7.mzData
>    https://zenodo.org/record/3631074/files/alg8.mzData
>    https://zenodo.org/record/3631074/files/alg9.mzData
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md collection=true format="mzml" collection_name="mzData" renaming=false %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md astype="as a Collection" %}
>
> 3. Make sure your data is in a **collection**. Make sure it is named `mzData`
>    - If you forgot to select the collection option during import, you can create the collection now:
>
>    {% snippet faqs/galaxy/collections_build_list.md %}
>
> 4. Import the following 2 files from Zenodo or from a shared data library (ask your instructor). 
> Beware: these files must not be in a collection. 
>
>    ```
>    https://zenodo.org/record/3631074/files/sampleMetadata.tsv
>    https://zenodo.org/record/3631074/files/W4M0004_database_small.msp
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
>
{: .hands_on}

## 2 - First steps using XCMS

> <hands-on-title>First steps using a standard XCMS workflow</hands-on-title>
>
> 1. {% tool [MSnbase readMSData](toolshed.g2.bx.psu.edu/repos/lecorguille/msnbase_readmsdata/msnbase_readmsdata/2.16.1+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"File(s) from your history containing your chromatograms"*: `mzData` (Input dataset collection)
>
> 2. {% tool [xcms findChromPeaks (xcmsSet)](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_xcmsset/abims_xcms_xcmsSet/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"RData file"*: `mzData.raw.RData` (output of the **MSnbase readMSData** {% icon tool %} job)
>    - *"Extraction method for peaks detection"*: `MatchedFilter - peak detection in chromatographic space`
>        - *"Full width at half maximum of matched filtration gaussian model peak"*: `5`
>        - *"Step size to use for profile generation"*: `0.5`
>        - In *"Advanced Options"*:
>            - *"Maximum number of peaks that are expected/will be identified per slice"*: `500`
>            - *"Signal to Noise ratio cutoff"*: `2`
>            - *"Minimum difference in m/z for peaks with overlapping Retention Times"*: `0.5`
>
>
>    > <comment-title></comment-title>
>    >
>    > With GC-MS data in profile mode, we need to use the *MatchedFilter* algorithm instead of the *Centwave* one used in the LC-MS tutorial.
>    {: .comment}
>
> 3. {% tool [xcms findChromPeaks Merger](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_merge/xcms_merge/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"RData file"*: `mzData.raw.xset.RData` (output of the **xcms findChromPeaks (xcmsSet)** {% icon tool %} job)
>    - {% icon param-file %} *"Sample metadata file "*: `sampleMetadata.tsv` (One of the uploaded files from Zenodo)
>
>    > <comment-title></comment-title>
>    >
>    > To merge your data, you need to **input a sampleMetadata file** containing filenames and their metadata informations like their class for example.
If you don't add a sampleMetadata file, the **xcms findChromPeaks Merger** {% icon tool %} tool will **group all your files together**. 
You can also **create your sampleMetadata file** with W4M Galaxy tool **xcms get a sampleMetadata file** {% icon tool %} with the following parameters: *"RData file"* outputed from **MSnbase readMSData** {% icon tool %}. 
Here is an example of the minimum expectations about a sampleMetadata file (**important**: don't write the format of the file, just their names):
>    > {: .text-justify}
>    > 
>    > | sample_name |  class  | 
>    > |:-----------:|:-------:|
>    > |    file1    |   man   |
>    > |-------------+---------|
>    > |    file2    |  woman  |
>    > |-------------+---------|
>    > |    file3    |   man   |
>    > 
>    {: .comment}
>
{: .hands_on}

The output from **xcms findChromPeaks Merger** {% icon tool %} is a *.RData* file that is mandatory to proceed with the end of the extraction process.
There are two available options: 
- carrying on using the standard XCMS workflow similarly to the one used in the LC-MS tutorial
- using the metaMS strategy specifically designed for GC-MS data
The two options are illustrated in this tutorial. 

# Carrying on using the standard XCMS workflow (option 1)

This option follows the standard LC-MS workflow to obtain in the end a dataMatrix file and its corresponding variableMetadata file. 

> <hands-on-title>Example of end of extraction when using the standard XCMS workflow</hands-on-title>
>
> 1. {% tool [xcms groupChromPeaks (group)](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_group/abims_xcms_group/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"RData file"*: `xset.merged.RData` (output of the **xcms findChromPeaks Merger** {% icon tool %} job)
>    - *"Method to use for grouping"*: `PeakDensity - peak grouping based on time dimension peak densities`
>        - *"Bandwidth"*: `10.0`
>        - *"Width of overlapping m/z slices"*: `0.05`
>
> 2. {% tool [xcms fillChromPeaks (fillPeaks)](toolshed.g2.bx.psu.edu/repos/lecorguille/xcms_fillpeaks/abims_xcms_fillPeaks/3.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"RData file"*: `xset.merged.groupChromPeaks.RData` (output of the **xcms groupChromPeaks (group)** {% icon tool %} job)
>    - In *"Peak List"*:
>        - *"Convert retention time (seconds) into minutes"*: `Yes`
>        - *"Number of decimal places for retention time values reported in ions' identifiers."*: `2`
>        - *"Reported intensity values"*: `maxo`
>
>
{: .hands_on}

The outputs of this strategy are similar to the ones discribed in the LC-MS tutotial mentioned previously. 

> <comment-title>Important : Be careful of the file format</comment-title>
>
> During each step of pre-processing, your dataset has its format changed and can have also its name changed.
> To be able to continue to MSMS processing, you need to have a RData object wich is **merged and grouped** (from **xcms findChromPeaks Merger** {% icon tool %} and **xcms groupChromPeaks (group)** {% icon tool %}) at least. 
It means that you should have a file named `xset.merged.groupChromPeaks.RData` (and maybe with some step more in it).
{: .comment} 


# Stopover : Verify your data after the XCMS pre-processing

When you have processed **all or only needed** steps described before, you can continue with the MS/MS processing part with **msPurity** package. 
Don't forget to always check your files format! For the next step you need to have this file `xset.merged.groupChromPeaks.*.RData` where * is the name of **optionnal** steps you could do during the pre-processing. 
For our example, your file should be named `xset.merged.groupchromPeaks.RData`. 
{: .text-justify}

> <comment-title></comment-title>
> 
> The pre-processing part of this analysis can be **quite time-consuming**, and already corresponds to quite a few number of steps, depending of your analysis. We highly recommend, at this step of the MS/MS workflow, to split your analysis by beginning a new Galaxy history with **only the files you need** (final xset Rdata file and your data collection of mzML). This will help you in limiting selecting the wrong dataset in further analysis, and bring a little **tidiness** for future review of your MS/MS analysis process. You should also be able to make a better peakpicking in the future in the same history and it will not be polluated by MS/MS part of your process.
> {: .text-justify}
> 
> > <tip-title>Copy dataset to a new history</tip-title>
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
> {% snippet faqs/galaxy/datasets_rename.md %}
>
{: .comment}

Before the next step with msPurity package on MS/MS datas, here are some questions to be able to verify if your file is ready and if you have the same results as us. Please check these questions : 
{: .text-justify}

> <question-title>before MS/MS steps</question-title>
> 
>  **1** - What are the steps of XCMS you made before your final file ?
> > <solution-title></solution-title>
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
> > <solution-title></solution-title>
> > 
> > During each step of XCMS pre-processing, the name of the file which is processing is completed by the name of the step you were doing. So, finally your file should be name `xset.merged.groupChromPeaks.fillChromPeaks.RData`. That because (as seen in previous answer) you ran a grouping and the integration after merged datas.
> > {: .text-justify}
> > 
> {: .solution}
> <br>
> **3** - What is the size (in MB) of your final RData file ?
> > <solution-title></solution-title>
> > 
> > To be able to see the size of a file in your history, you just have to select it. It will deploy informations about it and you can see the size of yours. For our example, the size of the final file is **1.4 MB**.
> > {: .text-justify}
> > 
> {: .solution}
{: .question}


# Processing with metaMS part (option 2)

**metaMS** is a R package for MS-based metabolomics data. It can do basic peak picking and grouping using functions from **XCMS** and **CAMERA** packages. The main output of **metaMS** is a table of feature intensities in all samples which can be analyzed with multivariate methods immediately. The package also offers the possibility to create in-house databases of mass spectra (including retention information) of pure chemical compounds. These databases can then be used for annotation purposes. The most important functions of this package are *runGC* and *runLC* (and each one to create databases *createSTDdbGC* and *createSTDdbLC*).
{: .text-justify}
During this tutorial we are interested in GC-MS analysis, so we will use the *runGC* function of **metaMS** and described it in details to be able to understand this function. The standard workflow of **metaMS** for GC-MS data is the following : 

![Workflow metaMS](../../images/tuto_gcms_workflow_metaMS.png "Workflow of metaMS for GC datas")

The *runGC* function is implemented in **metaMS.runGC {% icon tool %} tool** in W4M Galaxy. It takes a vector of file names, corresponding to the samples, and a settings list as mandatory arguments. In addition, some extra arguments can be provided. In particular, a database of standards, as discussed later in the tutorial, can be provided for annotation purposes. This tool regroups all these steps that are described in the following parts to be able to understand all its functionalities and particularities. We will run the tool after we understand each of its steps because it is important to know what are the best parameters for our data and why each parameter is done. 
{: .text-justify}


## Peak picking

The peak picking is performed by the usual **XCMS** functions. A function has been written in **metaMS** to allow the individual parameters to be passed to the function as a settings list. The result is that the whole of the **XCMS** functionality is available, simply by changing the values of some settings, or by adding fields. 
 {: .text-justify}
Whereas the package is not up-to-date since the new version of **XCMS** (3.x). This new version brought a lot of new objects and transformed the peak picking process. To have the last version of this process, **metaMS** authorized **to start its function directly with the file containing all peak picking results**. 
{: .text-justify}
Due to this update, **we have already processed the peak picking during the first part** of this tutorial. So we can continue it with the file outputted from the peak picking part. This also allow us to make a good peak picking without the following step include in **metaMS** functions. So it takes less time of processing and we can verify our peaks with this cut between peak picking and the following steps of GC-MS analysis. 
{: .text-justify}

## Definition of pseudo-spectra

Rather than a feature-based analysis with individual peaks, as is the case with **XCMS**, **metaMS** performs a pseudospectrum-based analysis. So, the basic entity is a set of m/z values showing a chromatographic peak at the same retention time.
{: .text-justify}

![Pseudospectrum example](../../images/tuto_gcms_pseudospectrum_example.png "Pseudospectrum example from `msp` file")

This choice is motivated by several considerations. First of all, **in GC the amount of overlap is much less than in LC** : peaks are much narrower. This means that even a one- or two-second difference in retention time can be enough to separate the corresponding mass spectra. Secondly, fragmentation patterns for many compounds are **available in extensive libraries like the [NIST library](http://www.nist.gov/srd/nist1a.cfm "NIST library")**. In addition, the spectra are somewhat easier to interpret since adducts, such as found in LC, are not present. The main advantage of pseudo-spectra, however, is that their use allows the results to be interpreted directly as relative concentrations of chemical compounds : **a fingerprint in terms of chemical composition is obtained**, rather than a fingerprint in terms of hard-to-interpret features. The pseudo-spectra are obtained by simply clustering on retention time, using the *runCAMERA* function, which for GC data calls *groupFWHM*. All the usual parameters for the *groupFWHM* function are included in W4M Galaxy **metaMS.runGC {% icon tool %} tool**. The most important parameter is *perfwhm*, which determines the maximal retention time difference of features in one pseudospectrum. 
{: .text-justify}

The final step is to convert the **CAMERA** objects into easily handled lists, which are basically the R equivalent of the often-used `msp` format from the AMDIS software ({% cite Stein1999 %}). The `msp`  file is a nested list, with one entry for each sample, and each sample represented by a number of fields. The pseudo-spectra are three-column matrices, containing m/z, intensity and retention time information, respectively. They can be draw with the *plotPseudoSpectrum* function of **metaMS** package easily (Figure 2).
{: .text-justify}


## Annotation

Once we have identified our pseudo-spectra, we can start the annotation process. This is done by **comparing every pseudospectrum to a database of spectra**. As a similarity measure, we use the weighted dot product as it is fast, simple, and gives good results ({% cite Stein1994 %}). The first step in the comparison is based on retention, since a comparison of either retention time or retention index is much faster than a spectral comparison. The corresponding function is *matchSamples2DB*. Since the weighted dot product uses scaled mass spectra, the scaling of the database is done once, and then used in all comparisons.
{: .text-justify}

![Match spectra](../../images/tuto_gcms_match_spec.png "Best match between an experimental pseudospectrum (red) and a database entry (blue)")

This *matchSamples2DB* function returns a table where all patterns that have a match with a DB entry are shown in the first column, and the DB entry itself in the second column. If for a particular experimental pattern more than one match is found, the alternatives (with a lower match factor) are shown in the last column. To see the match for a particular pattern, one can use the function *matchExpSpec*, returning matchfactors (numbers between 0 and 1, where the latter means a perfect match) for all entries in the database (if the plotIt argument is TRUE, the best match is shown – see Figure 2). Samples may contain compounds that are not of any interest, such as plasticizers, internal standards, column material etc.... These can be filtered out before doing an annotation : **metaMS** allows certain categories of database entries (defined in slot *matchIrrelevants* of the settings object) to be removed before further annotation. If the spectra of these compounds are very specific (and they often are), the retention criterion may be bypassed by setting the maximal retention time difference to very high values, which leads to the removal of such spectra wherever they occur in the chromatogram.
{: .text-justify}


## Unknowns research

The most important aspect of untargeted metabolomics is the definition of unknowns, patterns that occur repeatedly in several samples, but for which no annotation has been found. In **metaMS** these unknowns are found by comparing all patterns within a certain retention time (or retention index) difference on their spectral characteristics. The same match function is used, but the threshold may be different from the threshold used to match with the database of standards. Likewise, the maximum retention time(index)difference may be different, too. 

In defining unknowns we have so far used settings that are more strict than when compared to a database : since all samples are typically measured in one single run, expected retention time differences are rather small. In addition, one would expect reproducible spectra for a single compound. A true unknown, or at least an interesting one, is also present in a significant fraction of the samples. All these parameters are gathered in thebetweenSampleselement of the settingsobject .Since the matching is done using scaled patterns, we need to create a scaled version of the experimental pseudo-spectra first.
{: .text-justify}

For large numbers of samples, this process can take quite some time (it scales quadratically), especiallyif the allowed difference in retention time is large. The result now is a list of two elements : the first is the annotation table that we also saw after the comparison with the database, and the second is a list of pseudo-spectra corresponding to unknowns. In the annotation table, negative indices correspond to the pseudo-spectra in this list.
{: .text-justify}


## Outputs and results

At this stage, all elements are complete : we have the list of pseudo-spectra with an annotation, either as a chemical standard from the database, or an unknown occurring in a sizeable fraction of the injections. The only things left to do is to calculate relative intensities for the pseudo-spectra, and to put the results in an easy-to-use table. This table consists of two parts. The first part is the information on the “features”, which here are the pseudo-spectra. The second part of the table contains the intensities of these features in the individual injections. 
{: .text-justify}

![Match spectra](../../images/tuto_gcms_finale_table.png "Final table with unknowns and compounds found during **metaMS** processus")

The first five lines are the standards, and the next ones are the unknowns that are identified by the pipeline.  In the manual interpretation of this kind of data, the intensities of one or two “highly specific” features are often used to achieve relative quantitation. In an automatic pipeline, this is a risky strategy: not only can the intensity of a peak vary quite dramatically (relative standard deviations of up to 30% are assumed acceptable in GC-MS, e.g. when SPME is applied), but these errors are all the more pronounced in high-intensity peaks (hence the common use of a relative standard deviation).

In addition, one is ignoring the information in the other peaks of the pseudospectrum. In **metaMS**, pseudospectrum intensity is expressed as a multiple of the corresponding reference pattern (either a database pattern or an unknown), where the intensity ratio is determined using robust regression to avoid one deviating feature to influence the results too much ({% cite Wehrens2014 %}). First, we define an object containing all relevant pseudo-spectra, and next the intensities are generated.
{: .text-justify}

In both cases, the result is a list containing a set of patterns corresponding with the compounds that have been found, either annotated or unknown, the relative intensities of these patterns in the individual annotations, and possibly the xcmsSetobject for further inspection. In practice, the *runGC* function is all that users need to use.
{: .text-justify}

That file can be used for database search online (as Golm ({% cite Kopka2005 %}) and MassBank ({% cite Horai2010 %})) or locally (NIST MSSEARCH) for NIST search a tutorial is available here.
{: .text-justify}

> <hands-on-title>metaMS.runGC</hands-on-title>
> 
> We now know each step of this *runGC* function. So, please open the **metaMS.runGC {% icon tool %} too** to run it. You should enter the following parameters for our tutorial : 
>   - **Rdata from xcms and merged** : here you have to select your file from **XCMS** where you made the peak picking, grouping and all the pre-processing. It should be named `xset.merged.groupdChromPeaks.RData`.
>   - **Settings** : you can keep it at *user_default* but to see all possible parameters please set it at `use_defnied`.
>     - **RT range option** : it ables to select a region of retention time. If you select to *show* it, you have to enter the window in minutes, separate by a coma (for example 5,20 to have results between 5 minutes and 20 minutes). For our tutorial, we `keep it to hide`. 
>     - **RT_Diff** : it is the allowed retention time difference in minutes between the same compound/unknown in different sample. For our tutorial, `keep it at 0.05` to have low differences between unknowns' retention times.
>     - **Min_Features** : this parameter is used during the comparison with database or unknowns. It corresponds to the minimal number of features required to have a valid pseudospectrum. For our tutorial, please `keep it to 5` to have really good compounds.
>     - **similarity_threshold** : this parameter is also used for comparison. It is the minimum similarity allowed between peaks mass spectra to be considers as equal. For our tutorial, please `keep it to 0.7`.
>     - **min.class.fract** : it corresponds to the minimal fraction of samples in which a pseudospectrum is present before it is regarded as an unknown. For the tutorial, please `keep it to 0.5`.
>     - **min.class.size** : it corresponds to the minimum number of samples in which a pseudospectrum should be present before it is regarded as an unknown. For our tutorial, please `set it to 2` because we have classes with only 2 samples. 
>   - **Use Personnal DataBase option** : you can compare your datas to a personnal database. If you want to do it start to choose `show` in this parameter. Then you will be able to select your file. If not, keep it to `hide` and you will only have unknowns as results.
>     - **DB file** : this parameter will appear if you choose to show it. You just have to `select your database in your file` to add this here. Be careful, this database has to respect some rules (please look at *?????????????????? part*).
>   - **Use RI option** : choose if you want to use the RI for standards.
>     - **RI file** : enter here `your RI file` which have to contains two columns : retention time and retention indices. 
>   - **Use RI as filter** : just to know if you want to use RI parameter as a filter.
>     - **RIshift** : if you want to use RI as filter, please precise here the RI shift. For our tutorial `keep the previous parameter to FALSE`. 
>
>
{: .hands_on}


# Take a look at your results after metaMS processing

We choose to separate our first W4M Galaxy tool into 2 parts: the processing of GC-MS data (**metaMS.runGC {% icon tool %}**) and the plotting results of these data (**metaMS.plot {% icon tool %}**). So we now have the first part describes just before and the second part we will describe just after. This part allows users to see the TIC (Total Ion Chromatogram), BPC (Base Peak Chromatogram), and also all EICs (Extracted Ion Chromatogram) you want, from our previous result outputted from **metaMS.runGC {% icon tool %} tool**. 
{: .text-justify}

If you separated your samples into different classes, this tool can constructs TICs and BPCs one class against one class, in a `pdf` file (Figure 5) : 
{: .text-justify}

![TIC](../../images/tuto_gcms_tic.png "TIC comparing 2 classes")

Concerning EICs, it is possible to choose for which compound you want to draw an EIC when you run the W4M Galaxy tool. According to your choice, you will obtain EICs for one compound in each sample you enter in the previous **metaMS** part. 
{: .text-justify}

![TIC](../../images/tuto_gcms_eic.png "Example of EIC of the 'Unknown 1' in sample 'alg2'")

> <hands-on-title>metaMS.plot</hands-on-title>
> 
> This tool is very easy to run. It is an obligation to process **metaMS.runGC {% icon tool %}** before this one. After that, you just have to choose if you want or not to draw your TIC, BPC or EIC : 
>   - **Rdata from new_metaMS_runGC** : the file you obtained with the **metaMS.runGC {% icon tool %}** tool. It should be named `runGC.RData`.
>   - **Do you want to process for TIC(s) ?** : if you select "yes" you will obtain the `pdf` file containing each TIC from each class against each others. 
>   - **Do you want to process for BPC(s) ?** : if you select "yes" you will obtain the `pdf` file containing each BPC from each class against each others. 
>   - **Do you want to process for EIC(s) ?** : if you select "yes" you will have to choose which compound(s) and unknown(s) you want to obtain its EIC.
>     - **EIC_Unknown** : here please choose which compound(s) or unknown(s) you want to obtain according to the `peaktable.tsv` file. For out tutorial it can be interesting to have a look at **all** the EICs. So put the `value to 0`. 
>
{: .hands_on}


# Conclusion 
 


