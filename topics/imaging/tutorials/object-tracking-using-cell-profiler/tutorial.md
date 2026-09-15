---
layout: tutorial_hands_on

title: Object tracking using CellProfiler
level: Intermediate
subtopic: analyses
zenodo_link: https://doi.org/10.5281/zenodo.4567084
requirements:
  -
    type: "internal"
    topic_name: galaxy-interface
    tutorials:
        - workflow-editor
questions:
- How to segment and track objects in fluorescence time-lapse microscopy images?
- Can the same workflow be reused on data from another discipline, such as Earth observation (tracking atmospheric rivers)?
objectives:
- Segment fluorescent objects using CellProfiler in Galaxy
- Track objects over multiple frames using CellProfiler in Galaxy
- Reuse the same CellProfiler tracking pipeline on an Earth-observation time series (cross-discipline reuse)
time_estimation: 1H
key_points:
- CellProfiler in Galaxy can be used to track objects in time-lapse microscopy images
- The same Galaxy CellProfiler pipeline can track objects across scientific domains — here a nucleus-tracking pipeline is reused, with only its parameters retargeted, to track atmospheric rivers in climate reanalysis data

contributions:
  authorship:
    - sunyi000
    - beatrizserrano
    - jkh1
    - annefou
  funding:
    - oscars

tags:
  - Isolated object tracking
  - Multi-channel image
  - Earth observation
---


Most biological processes are dynamic and observing them over time can provide valuable insights. Combining fluorescent markers with time-lapse imaging is a common approach to collect data on dynamic cellular processes such as cell division (e.g. {% cite Neumann2010 %}, {% cite Heriche2014%}). However, automated time-lapse imaging can produce large amounts of data that can be challenging to process. One of these challenges is the tracking of individual objects as it is often impossible to manually follow a large number of objects over many time points.

To demonstrate how automatic tracking can be applied in such situations, this tutorial will track dividing nuclei in a short time-lapse recording of one mitosis of a syncytial blastoderm stage _Drosophila_ embryo expressing a GFP-histone gene that labels chromatin.

Tracking is done by first segmenting objects then linking objects between consecutive frames. Linking is done by matching objects and several criteria or matching rules are available. Here we will link objects if they significantly overlap between the current and previous frames.


> <tip-title>Prerequisites</tip-title>
> It is expected that you are already familiar with the Galaxy interface and the workflow editor. If this is not the case, we recommend you to complete the requirements listed at the start of this tutorial.
{: .tip}

> <warning-title>Important information: CellProfiler in Galaxy</warning-title>
> The Galaxy {% tool [CellProfiler](toolshed.g2.bx.psu.edu/repos/bgruening/cp_cellprofiler/cp_cellprofiler/3.1.9+galaxy0) %}  tool takes two inputs: a CellProfiler pipeline and an image collection.
> Some pipelines created with stand-alone CellProfiler may not work with the Galaxy CellProfiler tool. Some reasons are:
>  * The pipeline was built with a different version of CellProfiler. The Galaxy tool currently uses CellProfiler 3.9.
>  * Modules used by the pipeline aren't available in Galaxy.
>  * Parameters for some CellProfiler modules are limited/constrained compared to the stand-alone version, most notably:
>    - Parameters require manual input from the user whereas, in the stand-alone version, some modules can inherit parameter values from other modules.
>    - Input and output file locations are set by Galaxy and can't be set by the user.
>    - Metadata extraction from file names is limited to a set of fixed patterns.
>
> It is recommended to build a CellProfiler pipeline using the Galaxy interface if the pipeline is to be run by Galaxy.
{: .warning}



> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


The **same** CellProfiler tracking pipeline can be applied to time series from very different scientific domains. This tutorial can be followed with a **bioimaging** dataset (dividing cell nuclei in a fluorescence microscopy recording) or an **Earth-observation** dataset (atmospheric rivers in a climate-reanalysis time series). The pipeline is the same in both paths — you build it once, exactly the same way; only the input data and a few names (the object name and the file-name marker) change. Choose the path that interests you most.

{% include _includes/cyoa-choices.html option1="Bioimaging" option2="Earth" default="Bioimaging" text="Pick a dataset: a fluorescence-microscopy recording of dividing cell nuclei, or a climate-reanalysis time series of atmospheric rivers. You build and run the *same* CellProfiler pipeline either way." %}

<div class="Earth" markdown="1">

> <comment-title> Tracking atmospheric rivers with a cell-tracking pipeline </comment-title>
>
> In the Earth-observation path we reuse this **same** CellProfiler object-tracking pipeline — built to follow dividing nuclei — to detect and track **atmospheric rivers**: long, narrow filaments of intense water-vapour transport in the atmosphere (the "rivers in the sky" that deliver much of the world's extreme rainfall). This is a cross-discipline experiment from the [OSCARS-FIESTA](https://fiesta.eu) project. Why does it work? An atmospheric river is a **bright, elongated filament on a dark background** in a map of vertically integrated vapour transport (IVT) — morphologically much like the bright nuclei the pipeline was built for. A river that intensifies, drifts and splits over time is the climate-science analogue of a dividing nucleus, so the same *segment-then-link-by-overlap* approach tracks it. The preprocessing notebooks, the full Galaxy run provenance, and a citable archive are in the [OSCARS-FIESTA example repository](https://github.com/annefou/fiesta-galaxy-cellprofiler-eo).
{: .comment}

</div>


# Get data

<div class="Bioimaging" markdown="1">

This tutorial will use a time-lapse recording of nuclei progressing through mitotic anaphase during early _Drosophila_ embryogenesis. The nuclei are labelled on chromatin with a GFP-histone marker and have been imaged every 7 seconds using a laser scanning confocal microscope with a 40X objective.
The images are saved as a zip archive on Zenodo and need to be uploaded to the Galaxy server before they can be used.

![Dividing nuclei](../../images/object-tracking-using-cell-profiler/Dividing_nuclei.gif "Time lapse recording of anaphase nuclei in a Drosophila embryo.")

> <hands-on-title>Data upload — bioimaging</hands-on-title>
>
> 1. Create a new history for this tutorial.
>    When you log in for the first time, an empty, unnamed history is created by default. You can simply rename it.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import {% icon galaxy-upload %} the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library.
>    - **Important:** If setting the type to 'Auto-detect', make sure that after upload, the datatype is set to zip.
>
>    ```
>    {{ page.zenodo_link }}/files/drosophila_sample.zip
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Rename {% icon galaxy-pencil %} the file to drosophila_embryo.zip
>
{: .hands_on}

</div>

<div class="Earth" markdown="1">

This path uses a time series of **atmospheric-river** activity over the North Pacific during an early-February 2017 event. Each frame is a map of the vertically integrated water-vapour transport (IVT) magnitude from the ERA5 climate reanalysis, rendered so that the atmospheric river appears as a bright filament on a dark background — one frame every 6 hours. The frames are saved as a zip archive on Zenodo and need to be uploaded to the Galaxy server before they can be used.

![An atmospheric river in an ERA5 IVT field](../../images/object-tracking-using-cell-profiler/input_atmospheric_river_eo.png "One frame of the ERA5 IVT magnitude over the North Pacific: a bright atmospheric-river filament transporting water vapour towards North America, on a dark background.")

> <hands-on-title>Data upload — Earth observation</hands-on-title>
>
> 1. Create a new history for this tutorial.
>    When you log in for the first time, an empty, unnamed history is created by default. You can simply rename it.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import {% icon galaxy-upload %} the atmospheric-river IVT frames from Zenodo.
>    - **Important:** If setting the type to 'Auto-detect', make sure that after upload, the datatype is set to zip.
>
>    ```
>    https://zenodo.org/records/20813832/files/npacific_ivt_frames.zip
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Rename {% icon galaxy-pencil %} the file to npacific_ivt_frames.zip
>
{: .hands_on}

> <comment-title> How the atmospheric-river frames were prepared </comment-title>
> The frames come from the ERA5 reanalysis {% cite Hersbach2020 %} accessed through the public, anonymous ARCO-ERA5 archive. For each 6-hourly time step we took the eastward and northward vertically integrated vapour-transport components and combined them into the IVT magnitude, IVT = √(u² + v²), over the North Pacific (2 – 11 February 2017). Each field was rescaled to an 8-bit image (bright river on dark background) and written as a three-channel PNG named `NPacific_IVT_0000.png`, `NPacific_IVT_0001.png`, … so the filename carries the time index — exactly the structure the *Metadata* step below expects. These frames are archived on Zenodo {% cite fouilloux2026npacificivt %}; the preprocessing notebooks, the full Galaxy run provenance, and a citable software archive are in the [OSCARS-FIESTA example repository](https://github.com/annefou/fiesta-galaxy-cellprofiler-eo).
{: .comment}

</div>


# Create a CellProfiler pipeline in Galaxy

In this section, we will build a CellProfiler pipeline from scratch in Galaxy.
We need to:
  - Read the images and the metadata
  - Convert the colour images to grayscale
  - Segment the objects (nuclei / atmospheric rivers)
  - Extract features from the segmented objects
  - Perform tracking
  - Produce some useful output files

A pipeline is built by chaining together Galaxy tools representing CellProfiler modules and must start with the {% tool [Starting modules](toolshed.g2.bx.psu.edu/repos/bgruening/cp_common/cp_common/3.1.9+galaxy1) %} tool and end with the {% tool [CellProfiler](toolshed.g2.bx.psu.edu/repos/bgruening/cp_cellprofiler/cp_cellprofiler/3.1.9+galaxy0) %}  tool.

![Image of the workflow](../../images/object-tracking-using-cell-profiler/CP_object_tracking_pipeline.png "Overview of the CellProfiler pipeline using Galaxy tools.")


> <details-title>More details about the pipeline steps</details-title>
>    - Metadata is needed to tell CellProfiler what a temporal sequence of images is and what the order of images is in the sequence.
>    - CellProfiler is designed to work primarily with grayscale images. Since we don't need the colour information, we convert colour images to grayscale type.
>    - Segmentation means identifying the objects in each image. In CellProfiler, this is done by thresholding the intensity level in each image.
>    - When we perform tracking we're usually interested in quantifying how some properties of the objects evolve over time. Also, sometimes we may want to do tracking by matching objects based on some property of the objects (e.g. a shape measurement). Therefore, for each segmented object, we compute some features, i.e. numerical descriptors of some properties of the object.
>    - Tracking will provide the information required to allow downstream data analysis tools to link the features into a multidimensional time series.
>
>
{: .details}

The pipeline you build is the **same** for both paths. Only two pieces of text differ, depending on which dataset you chose:

<div class="Bioimaging" markdown="1">

> <comment-title> Names to use in the bioimaging path </comment-title>
> - The **file-name marker** used to recognise and order the images is `GFPHistone`.
> - The **object name** we give the segmented objects is `Nuclei` (and the tracked image is `TrackedNuclei`).
{: .comment}

</div>

<div class="Earth" markdown="1">

> <comment-title> Names to use in the Earth-observation path </comment-title>
> - The **file-name marker** used to recognise and order the images is `IVT`.
> - The **object name** we give the segmented objects is `Rivers` (and the tracked image is `TrackedRivers`).
{: .comment}

</div>


## Create a new workflow

> <hands-on-title>Creating a new workflow</hands-on-title>
> 1. Create a new workflow
>
>    {% snippet faqs/galaxy/workflows_create_new.md %}
>
{: .hands_on }

The next steps will add new tools using the workflow editor.
Remember to save the workflow when done (or anytime) to not lose your input parameters.


## Read the images

> <hands-on-title>Reading images</hands-on-title>
>
> 1. Add Inputs → **Input Dataset** to the workflow
>
> 2. Add {% icon tool %} **Unzip** with the following parameter:
>   - {% icon param-select %} *"Extract single file"*: `All files`
>
>
{: .hands_on}


## Get the metadata

> <comment-title></comment-title>
>
> The Starting modules tool combines four CellProfiler modules that are always used together at the start of a pipeline. These modules are:
>  * Images
>  * Metadata
>  * NamesAndTypes
>  * Groups
{: .comment}

<div class="Bioimaging" markdown="1">

> <hands-on-title>Getting metadata — bioimaging</hands-on-title>
>
> 1. {% tool [Starting modules](toolshed.g2.bx.psu.edu/repos/bgruening/cp_common/cp_common/3.1.9+galaxy1) %} with the following parameters:
>
>     - Images
>        - *"Do you want to filter only the images?"*: `Select the images only`
>
>     - Metadata
>        - *"Do you want to extract the metadata?"*: `Yes, specify metadata`
>
>          - {% icon param-repeat %} Insert new metadata
>            - *"Metadata extraction method"*: `Extract from file/folder names`
>            - *"Metadata source"*: `File name`
>            - *"Select the pattern to extract metadata from the file name"*: `field1_field2_field3`
>            - *"Extract metadata from"*: `Images matching a rule`
>
>              - *"Match the following rules"*: `All`
>              - {% icon param-repeat %} Insert filtering rules:
>                - *"Select the filtering criteria"*: `File`
>                - *""operator*: `Does`
>                - *"contain"*: `Contain regular expression`
>                - *"match_value"*: `GFPHistone`
>
>     - NamesAndTypes
>       - *"Process 3D"*: `No, do not process 3D data`
>       - *"Assign a name to"*: `Give images one of several names depending on the metadata`
>         - {% icon param-repeat %} Insert another image:
>           - *"Match the following rules"*: `All`
>           - {% icon param-repeat %} Insert rule:
>
>             - *"Select rule criteria"*: `File`
>
>               - *"operator"*: `Does`
>               - *"contain"*: `Contain regular expression`
>               - *"match_value"*: `GFPHistone`
>
>           - *"Select the image type"*: `Color image`
>           - *"Name to assign these images"*: `OrigColor`
>           - *"Set intensity range from"*: `Image metadata`
>
>       - *"Image matching method"*: `Metadata`
>
>     - Groups
>       - *"Do you want to group your images?"*: `Yes, group the images`
>       - *"Metadata category"*: `field1`
>
>
{: .hands_on}

> <question-title></question-title>
>
> How are we capturing metadata and what type of metadata are we getting?
>
> > <solution-title></solution-title>
> >
> > Metadata is obtained from the filenames by extracting three text fields separated by an underscore. The metadata we get is "DrosophilaEmbryo", "GFPHistone" and numbers with leading zeros. These three fields represent respectively the sample identifier, the marker visualized in the image and the index in the time series.
> >
> {: .solution}
>
{: .question}

</div>

<div class="Earth" markdown="1">

> <hands-on-title>Getting metadata — Earth observation</hands-on-title>
>
> 1. {% tool [Starting modules](toolshed.g2.bx.psu.edu/repos/bgruening/cp_common/cp_common/3.1.9+galaxy1) %} with the following parameters:
>
>     - Images
>        - *"Do you want to filter only the images?"*: `Select the images only`
>
>     - Metadata
>        - *"Do you want to extract the metadata?"*: `Yes, specify metadata`
>
>          - {% icon param-repeat %} Insert new metadata
>            - *"Metadata extraction method"*: `Extract from file/folder names`
>            - *"Metadata source"*: `File name`
>            - *"Select the pattern to extract metadata from the file name"*: `field1_field2_field3`
>            - *"Extract metadata from"*: `Images matching a rule`
>
>              - *"Match the following rules"*: `All`
>              - {% icon param-repeat %} Insert filtering rules:
>                - *"Select the filtering criteria"*: `File`
>                - *""operator*: `Does`
>                - *"contain"*: `Contain regular expression`
>                - *"match_value"*: `IVT`
>
>     - NamesAndTypes
>       - *"Process 3D"*: `No, do not process 3D data`
>       - *"Assign a name to"*: `Give images one of several names depending on the metadata`
>         - {% icon param-repeat %} Insert another image:
>           - *"Match the following rules"*: `All`
>           - {% icon param-repeat %} Insert rule:
>
>             - *"Select rule criteria"*: `File`
>
>               - *"operator"*: `Does`
>               - *"contain"*: `Contain regular expression`
>               - *"match_value"*: `IVT`
>
>           - *"Select the image type"*: `Color image`
>           - *"Name to assign these images"*: `OrigColor`
>           - *"Set intensity range from"*: `Image metadata`
>
>       - *"Image matching method"*: `Metadata`
>
>     - Groups
>       - *"Do you want to group your images?"*: `Yes, group the images`
>       - *"Metadata category"*: `field1`
>
>
{: .hands_on}

> <question-title></question-title>
>
> How are we capturing metadata and what type of metadata are we getting?
>
> > <solution-title></solution-title>
> >
> > Metadata is obtained from the filenames by extracting three text fields separated by an underscore. The metadata we get is "NPacific", "IVT" and numbers with leading zeros. These three fields represent respectively the region identifier, the variable visualized in the image (the integrated vapour transport) and the index in the time series.
> >
> {: .solution}
>
{: .question}

</div>

> <question-title></question-title>
>
> How will CellProfiler form a temporal sequence?
>
> > <solution-title></solution-title>
> >
> > Images will be grouped based on field1 which here is the sample/region identifier and ordered alpha-numerically (by default) which will order them by field3 (the time series index) since fields 1 and 2 are constant.
> >
> {: .solution}
>
{: .question}


## Convert the images to grayscale

> <hands-on-title>Colour to grayscale conversion</hands-on-title>
>
> 1. {% tool [ColorToGray](toolshed.g2.bx.psu.edu/repos/bgruening/cp_color_to_gray/cp_color_to_gray/3.1.9+galaxy0) %} with the following parameters:
>    - *"Select the input CellProfiler pipeline"*: Connect output of **Starting Modules** {% icon tool %} to input of {% tool ColorToGray %}
>    - *"Enter the name of the input image"*: `OrigColor`
>    - *"Conversion method"*: `Combine`
>        - *"Name the output image"*: `OrigGray`
>        - *"Image type"*: `RGB`
>          - *"Relative weight of the red channel"*: `1`
>          - *"Relative weight of the green channel"*: `1`
>          - *"Relative weight of the blue channel"*: `1`
>
>
{: .hands_on}


## Segmentation

The first step to track objects starts with the identification of those objects on the images.

<div class="Bioimaging" markdown="1">

> <hands-on-title>Nuclei segmentation</hands-on-title>
>
> 1. {% tool [IdentifyPrimaryObjects](toolshed.g2.bx.psu.edu/repos/bgruening/cp_identify_primary_objects/cp_identify_primary_objects/3.1.9+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: output of **ColorToGray** {% icon tool %}
>    - *"Use advanced settings?"*: `Yes, use advanced settings`
>        - *"Enter the name of the input image (from NamesAndTypes)"*: `OrigGray`
>        - *"Enter the name of the primary objects to be identified"*: `Nuclei`
>        - *"Typical minimum diameter of objects, in pixel units (Min)"*: `30`
>        - *"Typical maximum diameter of objects, in pixel units (Max)"*: `9999`
>        - *"Discard objects outside diameter range?"*: `Yes`
>        - *"Discard objects touching the border of the image?"*: `Yes`
>        - *"Threshold strategy"*: `Global`
>            - *"Thresholding method"*: `Otsu`
>                - *"Two-class or three-class thresholding?"*: `Three classes`
>                    - *"Assign pixels in the middle intensity class to the foreground or the background?"*: `Background`
>                - *Threshold smoothing scale"*: `1.3488`
>                - *Threshold correction factor"*: `1`
>                - *Lower bound on threshold"*: `0.01`
>                - *Upper bound on threshold"*: `1`
>        - *"Method to distinguish clumped objects"*: `Intensity`
>            - *"Method to draw dividing lines between clumped objects"*: `Intensity`
>                - *"Automatically calculate size of smoothing filter for declumping?"*: `Yes`
>                - *"Automatically calculate minimum allowed distance between local maxima?"*: `Yes`
>                - *"Speed up by using lower-resolution mage to find local maxima?"*: `Yes`
>        - *"Fill holes in identified objects"*: `After both thresholding and declumping`
>        - *"Handling of objects if excessive number of objects identified"*: `Continue`
>
{: .hands_on}

</div>

<div class="Earth" markdown="1">

> <hands-on-title>Atmospheric-river segmentation</hands-on-title>
>
> 1. {% tool [IdentifyPrimaryObjects](toolshed.g2.bx.psu.edu/repos/bgruening/cp_identify_primary_objects/cp_identify_primary_objects/3.1.9+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: output of **ColorToGray** {% icon tool %}
>    - *"Use advanced settings?"*: `Yes, use advanced settings`
>        - *"Enter the name of the input image (from NamesAndTypes)"*: `OrigGray`
>        - *"Enter the name of the primary objects to be identified"*: `Rivers`
>        - *"Typical minimum diameter of objects, in pixel units (Min)"*: `30`
>        - *"Typical maximum diameter of objects, in pixel units (Max)"*: `9999`
>        - *"Discard objects outside diameter range?"*: `Yes`
>        - *"Discard objects touching the border of the image?"*: `Yes`
>        - *"Threshold strategy"*: `Global`
>            - *"Thresholding method"*: `Otsu`
>                - *"Two-class or three-class thresholding?"*: `Three classes`
>                    - *"Assign pixels in the middle intensity class to the foreground or the background?"*: `Background`
>                - *Threshold smoothing scale"*: `1.3488`
>                - *Threshold correction factor"*: `1`
>                - *Lower bound on threshold"*: `0.01`
>                - *Upper bound on threshold"*: `1`
>        - *"Method to distinguish clumped objects"*: `Intensity`
>            - *"Method to draw dividing lines between clumped objects"*: `Intensity`
>                - *"Automatically calculate size of smoothing filter for declumping?"*: `Yes`
>                - *"Automatically calculate minimum allowed distance between local maxima?"*: `Yes`
>                - *"Speed up by using lower-resolution mage to find local maxima?"*: `Yes`
>        - *"Fill holes in identified objects"*: `After both thresholding and declumping`
>        - *"Handling of objects if excessive number of objects identified"*: `Continue`
>
{: .hands_on}

</div>


## Feature extraction

Once the objects of interest are identified, we extract features, i.e. numerical descriptors of object properties. We do this because we may be interested in analysing the evolution of these properties over time or want to use them in the tracking procedure to match objects over time.

<div class="Bioimaging" markdown="1">

> <hands-on-title>Shape features</hands-on-title>
>
> 1. {% tool [MeasureObjectSizeShape](toolshed.g2.bx.psu.edu/repos/bgruening/cp_measure_object_size_shape/cp_measure_object_size_shape/3.1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: output of **IdentifyPrimaryObjects** {% icon tool %}
>    - In *"new object"*:
>        - {% icon param-repeat %} *"Insert new object"*
>            - *"Enter the name of the object to measure"*: `Nuclei`
>    - *"Calculate the Zernike features?"*: `No`
>
{: .hands_on}

> <hands-on-title>Intensity features</hands-on-title>
>
> 1. {% tool [MeasureObjectIntensity](toolshed.g2.bx.psu.edu/repos/bgruening/cp_measure_object_intensity/cp_measure_object_intensity/3.1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: output of **MeasureObjectSizeShape** {% icon tool %}
>    - In *"new image"*:
>        - {% icon param-repeat %} *"Insert new image"*
>            - *"Enter the name of an image to measure"*: `OrigGray`
>    - In *"new object"*:
>        - {% icon param-repeat %} *"Insert new object"*
>            - *"Enter the name of the objects to measure"*: `Nuclei`
>
{: .hands_on}

</div>

<div class="Earth" markdown="1">

> <hands-on-title>Shape features</hands-on-title>
>
> 1. {% tool [MeasureObjectSizeShape](toolshed.g2.bx.psu.edu/repos/bgruening/cp_measure_object_size_shape/cp_measure_object_size_shape/3.1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: output of **IdentifyPrimaryObjects** {% icon tool %}
>    - In *"new object"*:
>        - {% icon param-repeat %} *"Insert new object"*
>            - *"Enter the name of the object to measure"*: `Rivers`
>    - *"Calculate the Zernike features?"*: `No`
>
{: .hands_on}

> <hands-on-title>Intensity features</hands-on-title>
>
> 1. {% tool [MeasureObjectIntensity](toolshed.g2.bx.psu.edu/repos/bgruening/cp_measure_object_intensity/cp_measure_object_intensity/3.1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: output of **MeasureObjectSizeShape** {% icon tool %}
>    - In *"new image"*:
>        - {% icon param-repeat %} *"Insert new image"*
>            - *"Enter the name of an image to measure"*: `OrigGray`
>    - In *"new object"*:
>        - {% icon param-repeat %} *"Insert new object"*
>            - *"Enter the name of the objects to measure"*: `Rivers`
>
{: .hands_on}

</div>


## Track the objects

With the objects and the relevant features measured, we are now ready to start the tracking step!

<div class="Bioimaging" markdown="1">

> <hands-on-title>Object tracking</hands-on-title>
>
> 1. {% tool [TrackObjects](toolshed.g2.bx.psu.edu/repos/bgruening/cp_track_objects/cp_track_objects/3.1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: output of **MeasureObjectIntensity** {% icon tool %}
>    - *"Enter the name of the objects to track"*: `Nuclei`
>    - *"Choose a tracking method"*: `Overlap`
>        - *"Maximum pixel distance to consider matches"*: `50`
>        - *"Filter objects by lifetime?"*: `No`
>        - *"Select display option?"*: `Color and Number`
>        - *"Save color-coded image?"*: `Yes`
>            - *"Name the output image"*: `TrackedNuclei`
>
{: .hands_on}

</div>

<div class="Earth" markdown="1">

> <hands-on-title>Object tracking</hands-on-title>
>
> 1. {% tool [TrackObjects](toolshed.g2.bx.psu.edu/repos/bgruening/cp_track_objects/cp_track_objects/3.1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: output of **MeasureObjectIntensity** {% icon tool %}
>    - *"Enter the name of the objects to track"*: `Rivers`
>    - *"Choose a tracking method"*: `Overlap`
>        - *"Maximum pixel distance to consider matches"*: `50`
>        - *"Filter objects by lifetime?"*: `No`
>        - *"Select display option?"*: `Color and Number`
>        - *"Save color-coded image?"*: `Yes`
>            - *"Name the output image"*: `TrackedRivers`
>
{: .hands_on}

</div>


## Visualize results

To make sure that the tracking has gone as expected, we will have a look at the original images together with the results of the segmentation step. And we will visualize them together (in tiles) for easier comparison.

<div class="Bioimaging" markdown="1">

> <hands-on-title>Visualize segmentation outcome</hands-on-title>
>
> 1. {% tool [OverlayOutlines](toolshed.g2.bx.psu.edu/repos/bgruening/cp_overlay_outlines/cp_overlay_outlines/3.1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: output of **TrackObjects** {% icon tool %}
>    - *"Display outlines on a blank image?"*: `No`
>        - *"Enter the name of image on which to display outlines"*: `OrigGray`
>        - *"Outline display mode"*: `Color`
>            - In *"Outline"*:
>                - {% icon param-repeat %} *"Insert Outline"*
>                    - *"Enter the name of the objects to display"*: `Nuclei`
>                    - *"Select outline color"*: {% color_picker #ff0000 %} (red)
>    - *"Name the output image"*: `OutlineImage`
>    - *"How to outline"*: `Inner`
>
{: .hands_on}

> <hands-on-title>Tiling images</hands-on-title>
>
> 1. {% tool [Tile](toolshed.g2.bx.psu.edu/repos/bgruening/cp_tile/cp_tile/3.1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: output of **OverlayOutlines** {% icon tool %}
>    - *"Enter the name of an input image"*: `OrigColor`
>    - *"Name the output image"*: `TiledImages`
>    - *"Tile assembly method"*: `Within cycles`
>        - *"Automatically calculate number of rows?"*: `No`
>            - *"Final number of rows"*: `1`
>        - *"Automatically calculate number of columns?"*: `Yes`
>        - *"Image corner to begin tiling"*: `top left`
>        - *"Direction to begin tiling"*: `row`
>        - *"Use meander mode?"*: `No`
>        - In *"Another image"*:
>            - {% icon param-repeat %} *"Insert Another image"*
>                - *"Enter the name of an additional image to tile"*: `OutlineImage`
>            - {% icon param-repeat %} *"Insert Another image"*
>                - *"Enter the name of an additional image to tile"*: `TrackedNuclei`
>
{: .hands_on}

</div>

<div class="Earth" markdown="1">

> <hands-on-title>Visualize segmentation outcome</hands-on-title>
>
> 1. {% tool [OverlayOutlines](toolshed.g2.bx.psu.edu/repos/bgruening/cp_overlay_outlines/cp_overlay_outlines/3.1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: output of **TrackObjects** {% icon tool %}
>    - *"Display outlines on a blank image?"*: `No`
>        - *"Enter the name of image on which to display outlines"*: `OrigGray`
>        - *"Outline display mode"*: `Color`
>            - In *"Outline"*:
>                - {% icon param-repeat %} *"Insert Outline"*
>                    - *"Enter the name of the objects to display"*: `Rivers`
>                    - *"Select outline color"*: {% color_picker #ff0000 %} (red)
>    - *"Name the output image"*: `OutlineImage`
>    - *"How to outline"*: `Inner`
>
{: .hands_on}

> <hands-on-title>Tiling images</hands-on-title>
>
> 1. {% tool [Tile](toolshed.g2.bx.psu.edu/repos/bgruening/cp_tile/cp_tile/3.1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: output of **OverlayOutlines** {% icon tool %}
>    - *"Enter the name of an input image"*: `OrigColor`
>    - *"Name the output image"*: `TiledImages`
>    - *"Tile assembly method"*: `Within cycles`
>        - *"Automatically calculate number of rows?"*: `No`
>            - *"Final number of rows"*: `1`
>        - *"Automatically calculate number of columns?"*: `Yes`
>        - *"Image corner to begin tiling"*: `top left`
>        - *"Direction to begin tiling"*: `row`
>        - *"Use meander mode?"*: `No`
>        - In *"Another image"*:
>            - {% icon param-repeat %} *"Insert Another image"*
>                - *"Enter the name of an additional image to tile"*: `OutlineImage`
>            - {% icon param-repeat %} *"Insert Another image"*
>                - *"Enter the name of an additional image to tile"*: `TrackedRivers`
>
{: .hands_on}

</div>


## Save the images and features

The tiled images and the features computed in previous steps are now exported to the Galaxy history ready to be analyzed.
> <hands-on-title>Save the images</hands-on-title>
>
> 1. {% tool [SaveImages](toolshed.g2.bx.psu.edu/repos/bgruening/cp_save_images/cp_save_images/3.1.9+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: output of **Tile** {% icon tool %}
>    - *"Select the type of image to save"*: `Image`
>        - *"Saved the format to save the image(s)"*: `png`
>    - *"Enter the name of the image to save"*: `TiledImages`
>    - *"Select method for constructing file names"*: `From image filename`
>        - *"Enter the image file name"*: `OrigColor`
>        - *"Append a suffix to the image file name"*: `Yes`
>            - *"Text to append to the image file name"*: `_tile`
>    - *"Overwrite existing files without warning?"*: `Yes`
>    - *"When to save"*: `Every cycle`
>    - *"Record the file and path information to the saved image?"*: `No`
>
>
{: .hands_on}

> <hands-on-title>Export tabular data to character-delimited text files</hands-on-title>
>
> 1. {% tool [ExportToSpreadsheet](toolshed.g2.bx.psu.edu/repos/bgruening/cp_export_to_spreadsheet/cp_export_to_spreadsheet/3.1.9+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: output of **SaveImages** {% icon tool %}
>    - *"Select the column delimiter"*: `Comma (",")`
>    - *"Add a prefix to file names?"*: `Do not add prefix to the file name`
>    - *"Overwrite existing files without warning"*: `Yes`
>    - *"Add image metadata columns to your object data file"*: `Yes`
>    - *"Representation of Nan/Inf"*: `NaN`
>    - *"Calculate the per-image mean values for object measurements?"*: `Yes`
>    - *"Calculate the per-image median values for object measurements?"*: `Yes`
>    - *"Calculate the per-image standard deviation values for object measurements?"*: `Yes`
>    - *"Create a GenePattern GCT file?"*: `No`
>
>    > <comment-title></comment-title>
>    >
>    > By default all measurement types are exported.
>    > This is in contrast to the stand-alone version of CellProfiler which allows choosing which measurement types to save.
>    {: .comment}
>
{: .hands_on}


# Run the pipeline with **CellProfiler**

As mentioned in the introduction, the tool {% tool [CellProfiler](toolshed.g2.bx.psu.edu/repos/bgruening/cp_cellprofiler/cp_cellprofiler/3.1.9+galaxy0) %}  takes as input the pipeline that we have been assembling over this tutorial and the input data. We are now ready to finally run CellProfiler!
> <hands-on-title>Running the pipeline with CellProfiler</hands-on-title>
>
> 1. {% tool [CellProfiler](toolshed.g2.bx.psu.edu/repos/bgruening/cp_cellprofiler/cp_cellprofiler/3.1.9+galaxy0) %}  with the following parameters:
>    - {% icon param-file %} *"Pipeline file"*: `output of **ExportToSpreadsheet** {% icon tool %}
>    - *"Are the input images packed into a tar archive?"*: `No`
>        - {% icon param-collection %} *"Images"*: output of **Unzip** {% icon tool %}
>    - *"Detailed logging file?"*: `Yes`
>
> 2. Run the workflow.
>
>    {% snippet faqs/galaxy/workflows_run.md %}
>
>    > <comment-title></comment-title>
>    >
>    > If the pipeline fails, inspect the CellProfiler log file for clues about errors.
>    > A common source of errors is not providing the correct information to a module such as pointing a CellProfiler module at a non-existent dataset (e.g. with a typo in the name).
>    {: .comment}
>
{: .hands_on}


<div class="Bioimaging" markdown="1">

![Sample output](../../images/object-tracking-using-cell-profiler/CP_object_tracking_sample_output.png "Sample of images produced by the pipeline.")

> <question-title></question-title>
>
> How many nuclei have been tracked?
>
> > <solution-title></solution-title>
> > Five objects were detected but only four of these were nuclei.
> >
> {: .solution}
>
{: .question}

> <question-title></question-title>
>
> How could we get rid of tracks that don't correspond to nuclei?
>
> > <solution-title></solution-title>
> >
> > Objects that are not nuclei appear to be transient so one solution would be to require the tracks to have a minimal length. This can be done using the *"Filter objects by lifetime"* parameter of the TrackObjects module.
> >
> {: .solution}
>
{: .question}

> <question-title></question-title>
>
> Where is the data needed to plot change in nuclei eccentricity over time?
>
> > <solution-title></solution-title>
> > The data is in the csv file called Nuclei in the CellProfiler pipeline output data set.
> >
> {: .solution}
>
{: .question}

</div>

<div class="Earth" markdown="1">

![Atmospheric-river tracking output](../../images/object-tracking-using-cell-profiler/output_atmospheric_river_tracking_eo.png "Sample output for the Earth-observation path, for one frame: the input IVT field (left), the identified atmospheric rivers outlined in red (middle), and the tracked rivers colour-coded and numbered (right). Each colour/number is followed across the time series.")

The pipeline outlines the atmospheric-river filaments and links them across the 6-hourly frames, assigning each a track number — exactly as it links dividing nuclei across microscopy frames.

> <question-title></question-title>
>
> Where is the data needed to plot how an atmospheric river's length changes over time?
>
> > <solution-title></solution-title>
> > The data is in the csv file called Rivers in the CellProfiler pipeline output data set: each row is one river in one frame, with its shape measurements (e.g. major axis length) and its TrackObjects label, so you can follow a single track over time.
> >
> {: .solution}
>
{: .question}

> <comment-title> Is this a scientifically valid atmospheric-river tracker? </comment-title>
> This tutorial shows that the *method* transfers: the same segment-then-link-by-overlap pipeline detects and tracks atmospheric rivers. It uses CellProfiler's generic intensity thresholding, **not** the established meteorological criteria {% cite Guan2015 %} (e.g. an IVT threshold plus a length > 2000 km and length-to-width > 2 geometry test). In the [OSCARS-FIESTA study](https://github.com/annefou/fiesta-galaxy-cellprofiler-eo), the same overlap-tracking approach with those criteria applied recovered five distinct atmospheric-river tracks over this ten-day window and followed one persistent river for about five and a half days. That repository (with a citable archive and a full provenance trail) is the place to go for a quantitative atmospheric-river analysis.
{: .comment}

</div>


# Conclusion


We've run a CellProfiler pipeline on Galaxy to segment and track objects in a time series, exported images to visually inspect the outcome, and saved tables of computed object features in comma-separated text files for future analysis.

You also saw that the **same pipeline** can be applied to data from a completely different discipline: a pipeline built to track dividing nuclei in fluorescence microscopy was reused — only its object name and file-name marker retargeted — to detect and track atmospheric rivers in a climate-reanalysis time series. A "bright objects that move, merge and split on a dark background" problem looks the same to CellProfiler whether the objects are nuclei under a microscope or rivers of water vapour seen from space. This kind of cross-discipline reuse is exactly what FAIR, well-described Galaxy workflows make possible.
