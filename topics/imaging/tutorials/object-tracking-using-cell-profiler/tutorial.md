---
layout: tutorial_hands_on

title: Object tracking using CellProfiler
zenodo_link: http://doi.org/10.5281/zenodo.4567084
questions:
- How to segment and track objects in fluorescence time-lapse microscopy images?
objectives:
- Segment fluorescent objects using CellProfiler in Galaxy
- Track objects over multiple frames using CellProfiler in Galaxy
time_estimation: 1H
key_points:
- CellProfiler in Galaxy can be used to track objects in time-lapse microscopy images
contributors:
- sunyi000
- beatrizserrano
- jkh1

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

Most biological processes are dynamic and observing them over time can provide valuable insights. Combining fluorescent markers with time-lapse imaging is a common approach to collect data on dynamic cellular processes such as cell division (e.g. {% cite Neumann2010 %}, {% cite Heriche2014%}). However, automated time-lapse imaging can produce large amount of data that can be challenging to process. One of these challenges is the tracking of individual objects as it is often impossible to manually follow a large number of objects over many time points. 
To demonstrate how automatic tracking can be applied in such situations, this tutorial will track dividing nuclei in a short time lapse recording of one mitosis of a syncytial blastoderm stage Drosophila embryo expressing a GFP-histone gene that labels chromatin.
Tracking is done by first segmenting objects then linking objects between consecutive frames. Linking is done by matching objects and several criteria or matching rules are available. Here we will link objects if they significantly overlap between the current and previous frames.



> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# CellProfiler in Galaxy

Give some background about what the trainees will be doing in the section.
Remember that many people reading your materials will likely be novices,
so make sure to explain all the relevant concepts.

# Get data

This tutorial will use a time-lapse recording of nuclei progressing through mitotic anaphase during early Drosophila embryogenesis. The nuclei are labelled on chromatin with a GFP-histone marker and have been imaged every 7 seconds using a laser scanning confocal microscope with a 40X objective.
The images are saved as a zip archive on Zenodo and need to be uploaded to the Galaxy server before they can be used.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial.  
>    When you log in for the first time, an empty, unnamed history is created by default. You can simply rename it.
>    {% snippet snippets/create_new_history.md %}
>
> 2. Import {% icon galaxy-upload %} the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library.
>    - **Important:** If setting the type to 'Auto-detect', make sure that after upload, the datatype is set to zip.
>
>    ```
>    https://zenodo.org/api/files/e5d1bd5c-60a0-42e4-8f0d-a2ebc863c5d9/drosophila_sample.zip
>    ```
>
>    {% snippet snippets/import_via_link.md %}
>    {% snippet snippets/import_from_data_library.md %}
>
> 3. Rename {% icon galaxy-pencil %} the file to drosophila_sample.zip
>
> 4. Unzip the file with the **Unzip file** {% icon tool %} tool with the following parameters:
>    - {% icon param-file %} *"input_file"*: `drosophila_sample.zip`
>
> 5. Rename {% icon galaxy-pencil %} the resulting collection to `anaphase nuclei`. You may need to refresh the history for the new name to appear.
>    {% snippet snippets/rename_collection.md %}
{: .hands_on}



# Option 1: Create a CellProfiler pipeline in Galaxy

In this section, we will build a CellProfiler pipeline from scratch in Galaxy.
We need to:  
  - Read the images and the metadata  
  - Convert the colour images to grayscale  
  - Segment the nuclei  
  - Extract features from the segmented nuclei  
  - Perform tracking  
  - Produce some useful output files  

A pipeline is built by chaining together Galaxy tools representing CellProfiler modules and must start with the 'Starting modules'{% icon tool %} tool and end with the 'CellProfiler'{% icon tool %} tool.

![Image of the pipeline](../../images/image_name "The pipeline")


> ### {% icon details %} More details about the pipeline steps
>    - Metadata is needed to tell CellProfiler what a temporal sequence of images is and what the order of images is in the sequence.
>    - CellProfiler is designed to work primarily with grayscale images. Since we don't need the colour information, we convert colour images to grayscale type.  
>    - Segmentation means identifying the nuclei in each image. In CellProfiler this is done by thresholding the intensity level in each image.  
>    - We often perform tracking because we're interested in quantifying how some properties of the objects evolve over time. Therefore for each segmented object we compute some features, i.e. numerical descriptors of some properties of the object.
>    - Tracking will provide the information required to allow downstream data analysis tools to link the features into a multidimensional time series 
> 
>
{: .details}


## Read the images and their metadata

> ### {% icon hands_on %} Hands-on: Reading images
>
> 1. {% tool [Starting modules](https://usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/bgruening/cp_common/cp_common/3.1.9+galaxy1) %} with the following parameters:
>    
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

## Change the images to gray color

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [ColorToGray](toolshed.g2.bx.psu.edu/repos/bgruening/cp_color_to_gray/cp_color_to_gray/3.1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **Starting Modules** {% icon tool %})
>    - *"Enter the name of the input image"*: `OrigColor`
>    - *"Conversion method"*: `Combine`
>        - *"Image type"*: `RGB`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Segmentation

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [IdentifyPrimaryObjects](toolshed.g2.bx.psu.edu/repos/bgruening/cp_identify_primary_objects/cp_identify_primary_objects/3.1.9+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **ColorToGray** {% icon tool %})
>    - *"Use advanced settings?"*: `Yes, use advanced settings`
>        - *"Enter the name of the input image (from NamesAndTypes)"*: `OrigGray`
>        - *"Enter the name of the primary objects to be identified"*: `Embryos`
>        - *"Typical minimum diameter of objects, in pixel units (Min)"*: `30`
>        - *"Typical maximum diameter of objects, in pixel units (Max)"*: `9999`
>        - *"Threshold strategy"*: `Global`
>            - *"Thresholding method"*: `Otsu`
>                - *"Two-class or three-class thresholding?"*: `Three classes`
>                    - *"Assign pixels in the middle intensity class to the foreground or the background?"*: `Background`
>                - *"Lower bound on threshold"*: `0.01`
>        - *"Method to distinguish clumped objects"*: `Intensity`
>            - *"Method to draw dividing lines between clumped objects"*: `Intensity`
>                - *"Automatically calculate size of smoothing filter for declumping?"*: `Yes`
>                - *"Automatically calculate minimum allowed distance between local maxima?"*: `Yes`
>        - *"Handling of objects if excessive number of objects identified"*: `Continue`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Measure nuclei size and shape

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [MeasureObjectSizeShape](toolshed.g2.bx.psu.edu/repos/bgruening/cp_measure_object_size_shape/cp_measure_object_size_shape/3.1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **IdentifyPrimaryObjects** {% icon tool %})
>    - In *"new object"*:
>        - {% icon param-repeat %} *"Insert new object"*
>            - *"Enter the name of the object to measure"*: `Embryos`
>    - *"Calculate the Zernike features?"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Measure the nuclei intensity

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [MeasureObjectIntensity](toolshed.g2.bx.psu.edu/repos/bgruening/cp_measure_object_intensity/cp_measure_object_intensity/3.1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **MeasureObjectSizeShape** {% icon tool %})
>    - In *"new image"*:
>        - {% icon param-repeat %} *"Insert new image"*
>            - *"Enter the name of an image to measure"*: `OrigGray`
>    - In *"new object"*:
>        - {% icon param-repeat %} *"Insert new object"*
>            - *"Enter the name of the objects to measure"*: `Embryos`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Track nuclei

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [TrackObjects](toolshed.g2.bx.psu.edu/repos/bgruening/cp_track_objects/cp_track_objects/3.1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **MeasureObjectIntensity** {% icon tool %})
>    - *"Enter the name of the objects to track"*: `Embryos`
>    - *"Choose a tracking method"*: `Overlap`
>        - *"Filter objects by lifetime?"*: `No`
>        - *"Save color-coded image?"*: `Yes`
>            - *"Name the output image"*: `TrackedCells`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **OverlayOutlines**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [OverlayOutlines](toolshed.g2.bx.psu.edu/repos/bgruening/cp_overlay_outlines/cp_overlay_outlines/3.1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **TrackObjects** {% icon tool %})
>    - *"Display outlines on a blank image?"*: `No`
>        - *"Enter the name of image on which to display outlines"*: `OrigGray`
>        - *"Outline display mode"*: `Color`
>            - In *"Outline"*:
>                - {% icon param-repeat %} *"Insert Outline"*
>                    - *"Enter the name of the objects to display"*: `Embryos`
>                    - *"Select outline color"*: `#ff0000`
>    - *"Name the output image"*: `OutlineImage`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Tile**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Tile](toolshed.g2.bx.psu.edu/repos/bgruening/cp_tile/cp_tile/3.1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **OverlayOutlines** {% icon tool %})
>    - *"Enter the name of an input image"*: `OrigColor`
>    - *"Name the output image"*: `AdjacentImage`
>    - *"Tile assembly method"*: `Within cycles`
>        - *"Automatically calculate number of rows?"*: `No`
>            - *"Final number of rows"*: `1`
>        - *"Automatically calculate number of columns?"*: `Yes`
>        - *"Use meander mode?"*: `No`
>        - In *"Another image"*:
>            - {% icon param-repeat %} *"Insert Another image"*
>                - *"Enter the name of an additional image to tile"*: `OutlineImage`
>            - {% icon param-repeat %} *"Insert Another image"*
>                - *"Enter the name of an additional image to tile"*: `TrackedCells`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SaveImages**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [SaveImages](toolshed.g2.bx.psu.edu/repos/bgruening/cp_save_images/cp_save_images/3.1.9+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **Tile** {% icon tool %})
>    - *"Select the type of image to save"*: `Image`
>        - *"Saved the format to save the image(s)"*: `png`
>    - *"Enter the name of the image to save"*: `AdjacentImage`
>    - *"Select method for constructing file names"*: `Single name`
>        - *"Enter single file name"*: `Specimen-FrameNumber`
>    - *"Overwrite existing files without warning?"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **ExportToSpreadsheet**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [ExportToSpreadsheet](toolshed.g2.bx.psu.edu/repos/bgruening/cp_export_to_spreadsheet/cp_export_to_spreadsheet/3.1.9+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **SaveImages** {% icon tool %})
>    - *"Select the column delimiter"*: ``
>    - *"Add a prefix to file names?"*: `Do not add prefix to the file name`
>    - *"Calculate the per-image median values for object measurements?"*: `No`
>    - *"Calculate the per-image standard deviation values for object measurements?"*: `No`
>    - *"Create a GenePattern GCT file?"*: `No`
>    - *"Export all measurement types?"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **CellProfiler**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [CellProfiler](toolshed.g2.bx.psu.edu/repos/bgruening/cp_cellprofiler/cp_cellprofiler/3.1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Pipeline file"*: `output_pipeline` (output of **ExportToSpreadsheet** {% icon tool %})
>    - *"Are the input images packed into a tar archive?"*: `No`
>        - {% icon param-collection %} *"Images"*: `output` (Input dataset collection)
>    - *"Detailed logging file?"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}


# Option 2: Upload a CellProfiler pipeline and run it in Galaxy

> ### {% icon hands_on %} Hands-on: Pipeline upload
>
> 1. Import {% icon galaxy-upload %} the cppipe file from your local computer.
>    > ### {% icon comment %} Comment
>    > Make sure the pipeline version is compatible with the CellProfiler version.
>    {: .comment}
>
> 2. Run {% tool [CellProfiler](toolshed.g2.bx.psu.edu/repos/bgruening/cp_cellprofiler/cp_cellprofiler/3.1.9) %} with the following parameters:
>    - {% icon param-file %} *"Pipeline file"*: `uploaded_pipeline` (output of **ExportToSpreadsheet** {% icon tool %})
>    - *"Are the input images packed into a tar archive?"*: `No`
>        - {% icon param-collection %} *"Images"*: `output` (Input dataset collection)
>    - *"Detailed logging file?"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}


# Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*



# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
