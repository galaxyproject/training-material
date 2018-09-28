---
layout: tutorial_hands_on

title: "Mass spectrometry imaging guide for beginners"
edam_ontology: "topic_0121"
zenodo_link: ""
questions:
  - "How to upload an imzML file into Galaxy"
  - "Getting an overview over the MSI data with the MSI quality report tool"
  - "How to export the underlaying data"

objectives:
  - "Basic understanding of handling imzML files in Galaxy and exploring the dataset with the quality report and imzML exporter tools"
time_estimation: "20min"
key_points:
  - "MSI datasets are uploaded via the composite function of Galaxy"
  - "The quality report tool provides many descriptive statistic plots that help to understand the properties of the data"
  - "The imzML exporter tools allowes The main data For analyzing cell culture or organic samples, search databases should include mycoplasma databases."
contributors:
  - foellmelanie
  - bgruening
---

# Introduction
{:.no_toc}

This tutorial introduces the handling of the mass spectrometry imaging (MSI) file types imzML in Galaxy and some first steps to explore the data. 
Mass spectrometry imaging is applied to measure the spatial distribution of biomolecules such as peptides, proteins, metabolites or chemical compounds. 
Depending on the analyte of interest, different ionization sources and mass analyzers are used in the mass spectrometer, for example MALDI-TOF or MALDI-FTICTR. 
Traditionally, the resulting MSI data is exported in the file format of the mass spectrometry vendor. With the introduction of the common standard imzML file format, more and more vendors provide an imzML export option and further software exists to convert proprietary files to an imzML file. 
The imzML file format was introduced to ease the exchange of MSI data between different instruments and data analysis software [Schramm et al., Journal of Proteomics, 2012](https://doi.org/10.1016/j.jprot.2012.07.026). 
Galaxy supports also the Analyze7.5 data format, that was originally used for MRI data, but ABSciex mass spectrometer make also use of this format. 
Independent of the file type all MSI data experiments contain meta data and for each measuring spot (pixel, with x and y coordinates) a full mass spectrum consisting of mass over charge (m/z) - intensity pairs. 
Often thousands of spectra with hundrets of m/z features are acquired leading to big and complex data. 
Before starting with the data analysis it is therefore helpful to visualize all levels of the data in different ways to better understand the data and to obtain an idea about it's quality and usefulness. 

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Uploading an imzML file into Galaxy

In case no direct imzML export is provided by the mass spectrometer software, proprietary files can be converted to imzML files with the tools provided on the following website: [ms-imaging.org](https://ms-imaging.org/wp/imzml/software-tools/).
The imzML file consists of two files: The first file contains the metadata in an XML file and has the extension .imzML. The second file contains the mass spectra data and is saved as binary file and its extension is .ibd. 
To be valid both files must have the same filename before the extension. 
Galaxy provides the 'composite' upload for files consisting of several components. 


TODO: provide test files; probably many m/z can be cut off to have a nice image with many pixels but not many m/z left?! Provide information only after QC - murine kidney etc...m/z are cut...Hauptfrage: Wie krüppelig sind unsere einzelnen Spektren in der Niere? 

TODO: comment: imzml can so far not be downloaded! 

> ### {% icon hands_on %} Hands-on: Uploading an imzML file
>
> 1. **Create a new history** and give it a name.
>
>    > ### {% icon tip %} Starting a new history
>    >
>    > * Click the **gear icon** at the top of the history panel
>    > * Select the option **Create New** from the menu
>    {: .tip}
>
> 2. **Upload imzML data** via the 'composite' option.
>    - Open the Galaxy **Upload Manager** next to the tool panel
>    - Select the **Composite** tab
>    - Set **Composite Type** to 'imzml'
>    - Press the first **Select** button (The imzML metadata component) and **choose local file**, select the .imzML file
>    - Press the second **Select** button (The mass spectral data component.) and **choose local file**, select the .ibd file
>    - Press **Start** and then **Close**
>
>    > ### {% icon tip %} Tip: FTP upload for large files
>    > * In case one subfile is larger than 2 GB the uploading needs to be done via ftp. 
>    > * The necessary steps are explained in this tutorial [Getting data into Galaxy]({{site.baseurl}}//galaxy-data-manipulation/tutorials/get-data/slides.html#18)
>    {: .tip}
>
> 3. **Rename dataset** If unhappy with the filename, rename the dataset.
>    - Click on the {% icon galaxy-pencil %} **pencil icon** for the dataset to edit its attributes.
>    - In the central panel, change the **Name** field to `my_first_imzML`
>    - Click the **Save** button>
>
>    > ### {% icon tip %} Tip: Uploading an Analyze7.5 file
>    > * Analyze7.5 files are also supported by Galaxy. 
>    > * The file consists of three components and is therefore uploaded via the 'composite' function, analogously to the imzML upload. 
>    > * The files to select in the 'composite' tab are the header file .hdr, the m/z values file .t2m and the spectra file .img.
>    {: .tip}
{: .hands_on}


# Exploring MSI data with the MSI quality report tool

Before starting any analysis, it is important to confirm that the quality of the acquired data is sufficient. 
Knowing the data's properties is important to choose the right preprocessing steps and parameters. 
The MSI quality report tool provides a fast way to obtain plenty of descriptive statistic plots for MSI data in imzML or Analyze7.5 format. 
The tool can use m/z values of internal calibrant or known sample features to provide intensity heatmap images and plots m/z accuracies.


> ### {% icon hands_on %} Hands-on: Running the MSI quality report tool
> 1. Upload the tabular file with the m/z values of the internal calibrants. 
>
>    > ### {% icon tip %} Tip: Importing data with copy and paste
>    >
>    > * Copy the table with m/z values
>    >```
>    >m/z       name
>    >805.42    tryptic_peptide
>    >1296.67   AngiotensinI
>    >```
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data** in the 'regular' tab
>    > * Paste the table into the field
>    > * Select the 'convert spaces to tabs' option which is available via the gear wheel above the input field
>    > * Press **Start** and then **Close**
>    {: .tip}
>
> 2. Select **MSI quality report** {% icon tool %} in the tool panel on the left side.
> 2. The imzML file will automatically be recognized and used as `input file`.
> 3. Select the tabular file with the m/z values `m/z of interest` drop down menu.
> 4. Set `ppm range` to at least the accuracy of the mass spectrometer, use '100' for this tutorial. 
> 5. Press **Execute**
{: .hands_on}

TODO: what plots are obtained, what do they say, what is good/bad quality
TODO: questions to number of spectra etc, how many m/z were valid? less because m/z range is cut off :D genial! 
TODO: questions: how far is Calibrant XY off (barplot?)
TODO: questions: why is number of calibrants everywhere same (not maximum probably): because not yet picked, some noise is everywhere
wuestions: number of peaks TIC correlation?! 

# Exporting MSI data to tabular files

To read and extract information from MSI data files special tools are needed. The imzML exporter provides the option to export spectra, feature and intensity data from imzML or Analyze7.5
files into tabular files. This may be helpful to follow up on interesting aspects of the quality report with exact values. Standard statistical programmes and software to perform further plots are normally compatible with tabular files. Care has to be taken, that  In case the intensity matrix of an unprocessed file is exported the tabular file will be huge and it might happen that 


- spectra output: spectra in rows - for each spectrum: name, x and y coordinates,order, number of peaks (intensities > 0), total ion chromatogram (TIC), highest m/z feature per spectrum, optional count of input m/z per spectrum, optional spectrum annotation
- mz feature output: m/z in rows - for each m/z: name, m/z, mean, median, standard deviation (sd), standard error of the mean (sem), sum of all intensities per m/z, number of peaks (intensity > 0) per m/z

Bsp tables spectra, feature and int. matrix information???

Von SGibb maldiquant: "This includes
checking the mass range, the length of each spectra and also visual exploration
of spectra to find and remove potentially defective measurements. a task often neglected" Cool wäre Datensatz mit ein paar zero TICs??? oder seehr niedrig? oder kaum peaks?


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
