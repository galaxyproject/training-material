---
layout: tutorial_hands_on

title: End-to-End Tissue Microarray Analysis with Galaxy-MTI
questions:
- What tools are available for pre-processing multiplex tissue images in Galaxy?
- What tools are available for downstream analysis of multiplex tissue images in Galaxy?
- How do I pre-process and analyze Tissue Microarray data?
- How can I visualize multiplex tissue images and associated data?
- How can I assign phenotypes to cells in an MTI dataset? 
objectives:
- Understand the tools available in Galaxy for multiplex tissue imaging analysis
- Analyze and visualize publicly available TMA data using Galaxy
time_estimation: 3H
key_points:
- Galaxy has tools and workflows that can be used to process and analyze multiplex tissue images
- Cell feature tables produced by the Galaxy TMA workflow can be used for downstream single-cell and spatial analyses
- There are powerful interactive visualization tools available in Galaxy that can combine the real images with associated data
- Tissue Microarray data can be analyzed using workflows that invoke MTI tools in batch
- Segmentation quality can vary significantly depending on features of the input image, tool used, and parameters
contributors:
- CameronFRWatson

---


# Introduction
{:.no_toc}

Multiplex tissue images are large, multi-channel images that contain intensity data for numerous biomarkers. The methods for generating multiplex tissue images are diverse, and each method can require specialized knowledge for downstream processing and analysis. The MCMICRO ({% cite Schapiro2021 %}) pipeline was developed to process multiplex images into single-cell data, and to have the range of tools to accomodate for different imaging methods. The tools used in the MCMICRO pipeline, in addition to tools for single-cell analysis, spatial analysis, and interactive visualization are available in Galaxy to facilitate comprehensive and accessible analyses of multiplex tissue images. The MCMICRO tools available in Galaxy are capable of processing Whole Slide Images (WSI) and Tissue Microarrays (TMA). WSIs are images in which a tissue section from a single sample occupies the entire microscope slide; whereas, TMAs multiplex smaller cores from multiple samples onto a single slide. This tutorial will demonstrate how to use the Galaxy multiplex imaging tools to process and analyze publicly available TMA test data provided by MCMICRO (Figure 1.).

Find a full example history [here](https://cancer.usegalaxy.org/u/watsocam/h/gtnexemplar002tma)

![exemplar viv](../../images/multiplex-tissue-imaging-TMA/ex2_combined_avivator.png "Fully registered image of the MCMICRO Exemplar-002 Tissue microarray. Exemplar-002 consists of four cores, each with a distinct tissue organization and expression of biomarkers. In the image, there are six biomarkers shown: DNA (white), CD163 (yellow), CD3D (blue), CD31 (red), VDAC1 (green), and Keratin (orange). This image is being viewed using Avivator, an interactive tool that allows the user to selectively view channels and adjust channel intensities.")


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Get data

Multiplex tissue images come in a variety of forms and file-types depending on the modality or platform used. For this tutorial, the Exemplar-002 data was imaged using Cyclic Immunofluorescence (CycIF) with a RareCyte slide scanner. Many of the steps in this workflow have platform-specific parameters, and the hands-on sections will show the best parameters for CycIF RareCyte images; however, notes will be made where critical differences may occur depending on the modality or platform throughout the tutorial.


The raw files for each round (10 in total) of the exemplar-002 data are available on [cancer.usegalaxy.org](https://cancer.usegalaxy.org) under **Data Libraries** (Figure 2.). Import the raw files into a new history as a **list collection**.


![getData](../../images/multiplex-tissue-imaging-TMA/ex2_getData.png "Finding the Exemplar-002 data on cancer.usegalaxy.org Data Libraries.")


> ### {% icon hands_on %} Hands-on: Data import to history
>
> 1. Create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the files for Exemplar-002 from the Shared Data Library to the new history
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
{: .hands_on}

> ### {% icon warning %} **Imaging platform differences**
> 
> The Exemplar-002 raw images are in *ome.tiff* format; however, commonly seen raw file-types are *ome.tiff*, *tiff*, *czi*, and *svs*. If your input images are not *ome.tiff* or *tiff*, you may have to edit the dataset attributes in Galaxy to allow tools to recognize them as viable inputs. 
>
{: .warning}


# Tile illumination correction with **BaSiC Illumination**

Commonly, raw MTI data will consist of one image per round of imaging. These individual round images are frequently captured in tiles, and there can be slight variations in how each tile was illuminated across the course of imaging. Prior to tile stitching and image registration, the tiles have to undergo illumination correction with **BaSiC Illumination** ({% cite Peng2017 %}) to account for this. Unlike many of the other tools in this workflow, BaSiC has no extra parameters to think about: Just input the collection of raw images and press *go*!

Two new list collections will appear in the history upon completion: 

  - BaSiC Illumination on Collection `X`: FFP (flat-field)
  - BaSiC Illumination on Collection `X`: DFP (deep-field)

> ### {% icon hands_on %} Hands-on: Illumination correction
>
> 1. {% tool [BaSiC Illumination](basic_illumination) %} with the following parameters:
>
>    - {% icon param-collection %} *"Raw Cycle Images: "*: List collection of raw images
>
{: .hands_on}


# Stitching and registration with **ASHLAR**

After illumination is corrected across round tiles, the tiles must be stitched together, and subsequently, each round mosaic must be registered together into a single pyramidal OME-TIFF file. **ASHLAR** from MCMICRO provides both of these functions. 

> ### {% icon comment %} Important detail: Marker File
>
> **ASHLAR** optionally reads a marker metadata file to name the channels in the output OME-TIFF image. This marker file will also be used in later steps. Make sure that the marker file is comma-separated and has the `marker_names` as the third column (Figure 3.). 
>
> ![markers](../../images/multiplex-tissue-imaging-TMA/ex2_markersFile.png "Markers file, used both in ASHLAR and downstream steps. Critically, the marker_names are in the third column.")
>
{: .comment}

> ### {% icon hands_on %} Hands-on: Image stitching and registration
>
> 1. {% tool [ASHLAR](ashlar) %} with the following parameters:
>
>    - {% icon param-collection %} *"Raw Images"*: List collection of raw images
>    - {% icon param-collection %} *"Deep Field Profile Images"*: List collection of DFP images produced by **BaSiC Illumination**
>    - {% icon param-collection %} *"Flat Field Profile Images"*: List collection of FFP images produced by **BaSiC Illumination**
>    - *"Flip X-axis"*: `No`
>    - *"Flip Y-axis"*: `No`
>    - *"Maximum allowed per-tile corrective shift"*: `30`
>    - *"Upgrade to BF6-Compliant OME-TIFF Pyramid"*: `Upgrade Pyramid`
>    - {% icon param-file %} *"Markers File (optional)"*: Comma-separated markers file with marker_names in third column
>
>    - In *"Advanced Options"*:
>        - *"Align Channel Number"*: `0` (Channel usually containing DAPI, hoescht, or other DNA marker)
>        - *"Sigma"*: `Not entered, left as default`
>        - *"Cyto mask channel"*: `Not entered, left as default`
>        - *"Flip output image horizontally"*: `No`
>        - *"Flip output image vertically"*: `No`
>        - *"Write output as a single pyramidal TIFF"*: `Yes`
>
{: .hands_on}

> ### {% icon warning %} **Imaging platform differences**
> 
> ASHLAR, among other tools in the MCMICRO and Galaxy-MTI pre-processing tools have some parameters that are specific to the 
imaging patform used. By default, ASHLAR is oriented to work with images from RareCyte scanners. AxioScan scanners render images
in a different orientation. Because of this, when using ASHLAR on AxioScan images, it is important to select the **Flip Y-Axis**
parameter to *Yes*
> 
> ASHLAR will work for most imaging modalities; however, certain modalities require different tools to be registered. For example,
multiplex immunohistochemistry (mIHC) images must use an aligner that registers each moving image to a reference Hematoxylin image. 
For this, Galaxy-MTI includes the alternative registration tool {% tool **PALOM** %}. 
>
{: .warning}


# TMA dearray with **UNetCoreograph**

Many downstream processing and analysis steps require each individual core from the TMA to be in a separate image file. To accomplish this from our registered ome.tiff image, we can use **UNetCoreograph** to detect and crop each core into separate files.

UNetCoreograph will output images (used for downstream steps), masks, and a preview image (Figure 4.).

![TMA](../../images/multiplex-tissue-imaging-TMA/ex2_dearray.png "Preview image from UNetCoreograph, outlines show detection of each individual core in the TMA.")

> ### {% icon hands_on %} Hands-on: TMA dearray
>
> 1. {% tool [UNetCoreograph](unet_coreograph) %} with the following parameters:
>
>    - {% icon param-file %} *"Registered TIFF"*: The output of **ASHLAR** (registered, pyramidal OME-TIFF file)
>    - *"Downsample factor"*: `5`
>    - *"Channel"*: `0`
>    - *"Buffer"*: `2.0`
>    - *"Sensitivity"*: `0.3`
>    - *"Cluster"*: `No`
>    - *"Tissue"*: `No`
>
>    > ### {% icon comment %} What about Whole Slide Images? 
>    >
>    > Whole slide images do not need to be dearrayed, so in most cases, this step can be skipped; however, UNetCoreograph has the *"Tissue"* option, which when selected, can act to separate the whole tissue from the background in a whole slide image which can be useful. In this case, it is important to toggle the *"Downsample factor"* as this often needs to be higher when extracting whole tissues.
>    {: .comment}
>
{: .hands_on}


# Nuclear segmentation with **Mesmer**

Cell segmentation is the basis for all downstream single-cell analyses. Different segmentation tools work highly variably depending on the imaging modality or platform used. Because of this, Galaxy-MTI has incorporated several cell segmentation tools so users may find the tool that works optimally for their data. 

Available segmentation tools in Galaxy-MTI:

  - Mesmer ({% cite Greenwald2021 %})
  - UnMicst and s3segmenter ({% cite Schapiro2021 %})
  - Cellpose ({% cite Stringer2020 %})
  - ilastik ({% cite Berg2019 %})

In this tutorial, we use **Mesmer** because it tends to perform generally well on a diverse range of image types, and has a limited number of parameters to understand. 

> ### {% icon comment %} Important detail: Running images in batches
>
> Now that each image has been split into individual core images, downstream tools must be run on the images separately. Luckily, Galaxy makes this easy by including the option to run each tool in batch across a collection of inputs. Next to the input for the tool, select {% icon param-collection %} (**Dataset collection**) as the input type, and pass the collection output by UNetCoreograph as input. 
>
{: .comment}

> ### {% icon hands_on %} Hands-on: Nuclear segmentation
>
> 1. {% tool [Mesmer](mesmer) %} with the following parameters:
>    - {% icon param-collection %} *"Image containing the nuclear marker(s) "*: Collection output of UNetCoreograph (images)
>    - *"The numerical index of the channel(s) from nuclear-image "*: `0`
>    - *"Compartment for segmentation prediction: "*: `Nuclear`
>    - *"Resolution of the image in microns-per-pixel"*: `0.65`
>    - *"Whether to np.squeeze the outputs before saving"*: `Yes`
>    - *"Segment with Cell Membrane"*: `No`
>
>    > ### {% icon comment %} np.squeeze
>    >
>    > The **np.squeeze** parameter is very important to select as `Yes` to make the output compatible with next steps
>    {: .comment}
>
{: .hands_on}

> ### {% icon warning %} **Imaging platform differences: Image resolution**
> 
> A crucial parameter for Mesmer and other segmentation tools is the **Image resolution**. This is reported in microns/pixel, and can vary depending on the imaging platform used and the settings at image acquisition. Mesmer accepts the resolution in microns/pixel; however, if using UNMICST, the resolution must be reported as a ratio of the resolution of UNMICST's training images (0.65). For example, when using UNMICST, if your images were captured at a resolution of 0.65, then the UNMICST value would be 1, but if your images were captured at 0.325 microns/pixel, then the value you would enter for UNMICST would be 0.5. 
>
{: .warning}


# Calculate single-cell features with **Quantification**

After generating a segmentation mask, the mask and the original registered image can be used to extract mean intensities for each marker in the panel, spatial coordinates, and morphological features for every cell. This step is performed by MCMICRO's **Quantification** module. 

Once again, as this is a TMA, we will be running this in batch mode for every core image and its segmentation mask. 

The quantification step will produce a CSV cell feature table for every image in the batch. 

> ### {% icon hands_on %} Hands-on: Quantification
>
> 1. {% tool [Quantification](quantification) %} with the following parameters:
>
>    - {% icon param-collection %} *"Registered TIFF "*: Collection output of UNetCoreograph (images)
>    - {% icon param-collection %} *"Primary Cell Mask "*: Collection output of Mesmer (or other segmentation tool)
>    - {% icon param-collection %} *"Additional Cell Masks "*: `Nothing Selected` (Other tools may produce multiple mask types)
>    - {% icon param-file %} *"Marker channels"*: Comma-separated markers file with marker_names in third column
>
>    > ### {% icon comment %} Mask metrics and Intensity metrics
>    >
>    > Leaving the *"mask metrics"* and *"intensity metrics"* blank will by default run all available metrics
>    >
>    {: .comment}
>
{: .hands_on}


# **Convert McMicro Output to Anndata**

Anndata ({% cite Virshup2021 %}) is a Python package and file format schema for working with annotated data matrices that has gained popularity in the single-cell analysis community. Many downstream analysis tools, including Scimap from MCMICRO, Scanpy ({% cite Wolf2018 %}), and Squidpy ({% cite Palla2022 %}) are built around anndata format files (h5ad). This tool splits the marker intensity data into a separate dataframe (`X`), and places all observational data (spatial coordinates, morphological features, etc.) in the cell feature table into a separate dataframe (`obs`) that shares the same indices as `X`. In downstream analyses, new categorical variables, such as phenotype assignments for each cell, are stored in the `obs` dataframe. 

Learn more about this file format at the [anndata documentation](https://anndata.readthedocs.io/en/latest/index.html).

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Convert McMicro Output to Anndata](scimap_mcmicro_to_anndata) %} with the following parameters:
>
>    - {% icon param-collection %} *"Select the input image or images"*: Collection output of Quantification (cellMaskQuant)
>    - In *"Advanced Options"*:
>        - *"Whether to remove the DNA channels from the final output"*: `No`
>        - *"Whether to log the data"*: `Yes`
>        - *"Name of the column that contains the CellID"*: `CellID` (Default)
>        - *"Whether to use unique name for cells/rows"*: `No`
>        - *"Column name to split the counts table and metadata"*: `X_centroid` (Default)
>
>    > ### {% icon warning %} Important parameter: Unique names for cells/rows
>    >
>    > Setting *"Whether to use unique name for cells/rows"* to `No` to ensures that downstream interactive visualizations will be able to map observational features to the mask CellIDs. 
>    {: .comment}
>
{: .hands_on}


# Scimap: **Single Cell Phenotyping**

There are several ways to classify cells available in Galaxy-MTI. Unsupervised approaches, such as Leiden clustering, can be performed on all cells and phenotypes can be manually annotated based on marker expression patterns observed by the user. This approach is time consuming, so here we will demonstrate automated phenotyping based on thresholds of specific lineage markers using MCMICRO's Scimap. Scimap phenotyping can either be provided a table of manual gate values for each marker of interest (which can be determined using the **GateFinder** tool in Galaxy-MTI), or by default, Scimap will fit a Gaussian Mixture Model (GMM) to the `log(intensity)` data for each marker to determine positive and negative populations for that marker. The marker intensity values are rescaled between (0,1) with 0.5 being the cut-off between negative and positive populations. Scimap uses a 'Phenotype workflow' to guide the classification of cells (Figure 5.). For more on how to construct a Scimap workflow, see the [Scimap documentation](https://scimap-doc.readthedocs.io/en/latest/tutorials/scimap-tutorial-cell-phenotyping/).


![phenoWF](../../images/multiplex-tissue-imaging-TMA/ex2_phenotypeWF.png "Example of a phenotype workflow compatible with Scimap. 'Pos' means that the marker must be positive to be classified as the respective phenotype. 'Anypos' means any, but not necessarily all, of the listed markers can be positive to call the respective phenotype.")

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Single Cell Phenotyping](scimap_phenotyping) %} with the following parameters:
>
>    - {% icon param-collection %} *"Select the input anndata"*: Output of **Convert MCMICRO output to Anndata**
>    - {% icon param-file %} *"Select the dataset containing manual gate information"*: (Optional) manually determined gates in CSV format. Gates will be determined automatically using a GMM for each marker if this file is not provided
>    - {% icon param-file %} *"Select the dataset containing gating workflow"*: CSV phenotype workflow (Figure 5.)
>    - *"Save the GMM gates plots If True"*: `Yes`
>
>
>    > ### {% icon comment %} Limitations of GMM automated phenotyping
>    >
>    > When manual gates are not provided, Scimap fits a GMM to determine a threshold between positive and negative cells. This automated gating works well when markers are highly abundant within the tissue, and the data shows a bimodal distribution (Figure 6A.). GMM gating can lead to spurious thresholds, however, when the data does not appear to be bimodal (Figure 6B.). This tends to happen when the marker is not highly abundant in the tissue, so there isn't a large positive population. Markers that have a highly continuous range of intensity, like certain functional markers, can also be problematic with GMM gating. It is recommended to always look at the GMM plots output by Scimap, and validate any potentially spurious gates manually. 
>    >
>    > ![GMM](../../images/multiplex-tissue-imaging-TMA/ex2_example_GMMs.png "Scimap automatic gating GMMs for two markers. (A) An example of a marker with a bimodal distribution and a reasonable looking gate. (B) An example of a marker with a unimodal distribution that is not ideal for fitting with a GMM, and would be a candidate for manual validation and gating.")
>    > 
>    {: .comment}
>
{: .hands_on}


# Interactive visualization of multiplex tissue images

Visual analysis is an important part of multiplex tissue imaging workflows. Galaxy-MTI has several tools that make interactive visualization easy, and can be used at various stages of analysis. 

## Converting UNetCoreograph images to OME-TIFF using the **Convert image** tool

UNetCoreograph outputs each individual core image in `tiff` format. Interactive visualization tools, such as **Vitessce** and **Avivator** require the images to be in `OME-TIFF` format to be viewed. Galaxy-MTI includes a conversion tool that can accomodate this, along with many other useful conversion functions. 

> ### {% icon hands_on %} Hands-on: Convert image
>
> 1. {% tool [Convert image](ip_convertimage) %} with the following parameters:
>    - *"Output data type"*: `OME TIFF`
>    - *"Extract series"*: `All series`
>    - *"Extract timepoint"*: `All timepoints`
>    - *"Extract channel"*: `All channels`
>    - *"Extract z-slice"*: `All z-slices`
>    - *"Extract range"*: `All images`
>    - *"Extract crop"*: `Full image`
>    - *"Tile image"*: `Tile image`
>    - *"Pyramid image"*: `Generate Pyramid`
>
{: .hands_on}


## **Rename OME-TIFF Channels**

Some tools can cause the channel names in an OME-TIFF image to be lost. To fix this, or to change the channel names to whatever the user prefers, the **Rename OME-TIFF Channels** tool can be invoked using a markers file similar to the one used in previous steps. 

> ### {% icon hands_on %} Hands-on: Rename channels
>
> 1. {% tool [Rename OME-TIFF Channels](toolshed.g2.bx.psu.edu/repos/watsocam/rename_tiff_channels/rename_tiff_channels/0.0.1.2) %} with the following parameters:
>
>    - {% icon param-file %} *"Input image in either tiff or OME-tiff format"*: File to be converted (or collection of files to run in batch)
>    - *"Format of input image"*: `ome.tiff`
>    - {% icon param-file %} *"Channel metadata CSV"*: Comma-separated markers file with marker_names in third column
>
{: .hands_on}


## Initial visualization with **Avivator**

For any `OME-TIFF` image in a Galaxy-MTI history, there will be an option to view the image using **Avivator**. This is a great way to perform an initial inspection of an image for QC purposes before continuing with downstream steps. The **Avivator** window can be launched by expanding the dataset information in the history panel and clicking the link (Figure 7.).

![vivhistory](../../images/multiplex-tissue-imaging-TMA/ex2_avivatorHistory.png "The highlighted link automatically appears for any OME-TIFF image (left) and, when clicked, launches an Avivator window to explore the image (right).")


## Generating an interactive visualization dashboard with **Vitessce**

**Vitessce** is a powerful visualization tool that creates interactive dashboards (Figure 8.) to look at a multiplex `OME-TIFF` images in conjunction with data generated during analysis and stored in an anndata file. The segmentation mask can be overlaid onto the image to qualitatively assess the segmentation performance. The mask can then be colored with associated observational data (Figure 9A.), such as `phenotype`, with the same colors appearing in barplots (Figure 9B.), UMAP representations, heatmaps, and marker intensity violin plots for comrehensive data exploration. 

![vitessce](../../images/multiplex-tissue-imaging-TMA/ex2_fullVitessce.png "A Full view of a vitesse dashboard for one core from Exemplar-002.")

![vitessce](../../images/multiplex-tissue-imaging-TMA/ex2_vitessce_zoomed.png "Each window in the dashboard can be resized to view the components in more detail. A closer look at the phenotype-labeled mask overlaid on the actual image (A), and the phenotype barplot (B).")

> ### {% icon hands_on %} Hands-on: Vitessce visualization
>
> 1. {% tool [Vitessce Visualization](vitessce_spatial) %} with the following parameters:
>
>    - {% icon param-file %} *"Select the OME Tiff image"*: OME-TIFF image to be viewed (or collection of files to run in batch)
>    - {% icon param-file %} *"Select masks for the OME Tiff image (Optional)"*: Output of Mesmer (or other segmentation tool)
>    - *"Whether to do phenotyping"*: `Yes`
>        - *"Select an embedding algorithm for scatterplot"*: `UMAP`
>        - *"Input phenotyping keys"*: `Multiple choices`
>            - *"Select the key(s)"*: `phenotype`
>
{: .hands_on}

# Next steps: Compositional and spatial analyses

Galaxy-MTI includes additional tools from **Scimap** and tools from the **Squidpy** package ({% cite Palla2022 %}) that can be used to perform a variety of downstream analyses. For example, once phenotypes have been assigned to individual cells, **Squidpy** has several methods for understanding the spatial organization of the tissue. Using **Squidpy**, a spatial neighborhood graph is first generated, from which the organization of specific phenotype groups and their interactions can be quantified. 

> ### {% icon hands_on %} Hands-on: Spatial analysis with **Squidpy**
>
> 1. {% tool **Squidpy Graph and Plotting** %} generate a spatial neighborhood graph with the following parameters:
>
>    - {% icon param-file %} *"Select the input anndata"*: Anndata file containing phenotype information (or other variable of interest)
>    - *"Select an analysis"*: `Spatial neighbors -- Create a graph from spatial coordinates`
>
> 2. {% tool **Squidpy Graph and Plotting** %} compute and plot a neighborhood enrichment analysis with the following parameters:
>
>    - {% icon param-file %} *"Select the input anndata"*: Output of step 1 (anndata file with spatial neighborhood graph)
>    - *"Select an analysis"*: `nhood_enrichment -- Compute neighborhood enrichment by permutation test`
>    - *"Key in anndata.AnnData.obs where clustering is stored"*: `phenotype`
>
>    > ### {% icon comment %} Neighborhood enrichment plot
>    >
>    > **Squidpy** was used to calculate neighborhood enrichments for each phenotype in core 2 of exemplar 2 (Figure 10.). This shows which phenotypes co-locate most frequently within the tissue. 
>    >
>    > ![enrichment](../../images/multiplex-tissue-imaging-TMA/ex2_squidpy_enrichment.png "The output of Squidpy's neighborhood enrichment on core 2 from Exemplar-002.")
>    >
>    {: .comment}
>
> 3. {% tool **Squidpy Graph and Plotting** %} calculate Ripley's L curves for each phenotype with the following parameters:
>
>    - {% icon param-file %} *"Select the input anndata"*: Output of step 1 (anndata file with spatial neighborhood graph)
>    - *"Select an analysis"*: `nhood_enrichment -- Compute neighborhood enrichment by permutation test`
>    - *"Key in anndata.AnnData.obs where clustering is stored"*: `phenotype`
>    - In *"Advanced Graph Options"*:
>        - *"Which Ripley's statistic to compute"*: `L`
>    - In *"Plotting Options"*:
>        - *"Ripley's statistic to be plotted"*: `L`
>
>    > ### {% icon comment %} Ripley's L plot
>    >
>    > **Squidpy** was used to calculate Ripley's L curves for each phenotype in core 2 of exemplar 2 (Figure 11.). This shows the overall organization of each phenotype in the tissue. If the curve for a given phenotype lies above the light grey null line (Example: Epithelial cells in Figure 11.), the phenotype is statistically significantly clustered. If the curve lies on the null line (Example: Myeloid lineage in Figure 11.), it's spatial distribution within the tissue is random. If the curve is underneath the null line (Example: T cells in Figure 11.), it's spatial distribution is statistically significantly dispersed. 
>    >
>    > ![enrichment](../../images/multiplex-tissue-imaging-TMA/ex2_squidpy_ripleys.png "The output of Squidpy's Ripley's L curve on core 2 from Exemplar-002.")
>    >
>    {: .comment}
>
{: .hands_on}


# Conclusion

In this tutorial, we demonstrated a complete multiplex tissue imaging analysis workflow performed entirely in a web browser using Galaxy-MTI. Using an example tissue microarray imaged with cylic immunofluoresence provided by MCMICRO, we...

  - Corrected illumination between imaging tiles
  - Stitched and registered input images to produce a single, pyramidal OME-TIFF image that is viewable in multiple built-in interactive viewing tools (Avivator, Vitessce)
  - Split the TMA into separate images for each core
  - Processed each core in parallel, beginning with nuclear segmentation
  - Quantified the mean marker intensities, morphological features, and spatial coordinates of each cell in each core
  - Converted the resulting tabular data to anndata format for convenient downstream anaylses and visualizations
  - Performed marker-based, automatically gated, phenotyping of cells
  - Prepared the dearrayed images and viewed them interactively in a dashboard combined with observational data

![galaxyWF](../../images/multiplex-tissue-imaging-TMA/ex2_galaxyWF.png "The entire workflow used in this tutorial.")