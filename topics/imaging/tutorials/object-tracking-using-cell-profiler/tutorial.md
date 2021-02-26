---
layout: tutorial_hands_on

title: Object tracking using CellProfiler
zenodo_link: http://doi.org/10.5281/zenodo.4567084
questions:
- Which biological questions are addressed by the tutorial?
- Which bioinformatics techniques are important to know for this type of data?
objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentences describing what a learner should be able to do once they
  have completed the tutorial
- You can use Bloom's Taxonomy to write effective learning objectives
time_estimation: 1H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- sunyi000
- beatrizserrano
- jkh1

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

General introduction about the topic and then an introduction of the
tutorial (the questions and the objectives). It is nice also to have a
scheme to sum up the pipeline used during the tutorial. The idea is to
give to trainees insight into the content of the tutorial and the (theoretical
and technical) key concepts they will learn.

You may want to cite some publications; this can be done by adding citations to the
bibliography file (`tutorial.bib` file next to your `tutorial.md` file). These citations
must be in bibtex format. If you have the DOI for the paper you wish to cite, you can
get the corresponding bibtex entry using [doi2bib.org](https://doi2bib.org).

With the example you will find in the `tutorial.bib` file, you can add a citation to
this article here in your tutorial like this:
{% raw %} `{% cite Batut2018 %}`{% endraw %}.
This will be rendered like this: {% cite Batut2018 %}, and links to a
[bibliography section](#bibliography) which will automatically be created at the end of the
tutorial.


**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**


***TODO***: Explain the biology

***TODO***: Explain object tracking

***TODO***: Explain the dataset (Drosophila embryos) from [CellProfiler examples](https://cellprofiler.org/examples). We can insert here a few of the images so that the 'movement' is shown.



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

***TODO***: Maybe here we can say that CP in Galaxy has some minor changes compared to the UI. Or maybe not, not sure about this section.

## Title for a subsection
Section and subsection titles will be displayed in the tutorial index on the left side of
the page, so try to make them informative and concise!

# Hands-on Sections
Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    https://zenodo.org/api/files/e5d1bd5c-60a0-42e4-8f0d-a2ebc863c5d9/drosophila_sample.zip
>    ```
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
{: .hands_on}


***TODO***: We have a zip in zenodo, need to use the tool unzip to get all the images into a collection. 



# Option 1: Create your CellProfiler pipeline in Galaxy

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


## Sub-step with **ColorToGray**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [ColorToGray](toolshed.g2.bx.psu.edu/repos/bgruening/cp_color_to_gray/cp_color_to_gray/3.1.9) %} with the following parameters:
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

## Sub-step with **IdentifyPrimaryObjects**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [IdentifyPrimaryObjects](toolshed.g2.bx.psu.edu/repos/bgruening/cp_identify_primary_objects/cp_identify_primary_objects/3.1.9) %} with the following parameters:
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

## Sub-step with **MeasureObjectSizeShape**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [MeasureObjectSizeShape](toolshed.g2.bx.psu.edu/repos/bgruening/cp_measure_object_size_shape/cp_measure_object_size_shape/3.1.9) %} with the following parameters:
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

## Sub-step with **MeasureObjectIntensity**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [MeasureObjectIntensity](toolshed.g2.bx.psu.edu/repos/bgruening/cp_measure_object_intensity/cp_measure_object_intensity/3.1.9) %} with the following parameters:
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

## Sub-step with **TrackObjects**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [TrackObjects](toolshed.g2.bx.psu.edu/repos/bgruening/cp_track_objects/cp_track_objects/3.1.9) %} with the following parameters:
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
> 1. {% tool [OverlayOutlines](toolshed.g2.bx.psu.edu/repos/bgruening/cp_overlay_outlines/cp_overlay_outlines/3.1.9) %} with the following parameters:
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
> 1. {% tool [Tile](toolshed.g2.bx.psu.edu/repos/bgruening/cp_tile/cp_tile/3.1.9) %} with the following parameters:
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
> 1. {% tool [SaveImages](toolshed.g2.bx.psu.edu/repos/bgruening/cp_save_images/cp_save_images/3.1.9) %} with the following parameters:
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
> 1. {% tool [ExportToSpreadsheet](toolshed.g2.bx.psu.edu/repos/bgruening/cp_export_to_spreadsheet/cp_export_to_spreadsheet/3.1.9) %} with the following parameters:
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
> 1. {% tool [CellProfiler](toolshed.g2.bx.psu.edu/repos/bgruening/cp_cellprofiler/cp_cellprofiler/3.1.9) %} with the following parameters:
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


# Option 2: Upload your CellProfiler pipeline and run it in Galaxy

## Sub-step 1

***TODO***: How to upload a pipeline file (take the snippet upload data sas example) and call it `uploaded_pipeline`.

## Sub-step with **CellProfiler**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [CellProfiler](toolshed.g2.bx.psu.edu/repos/bgruening/cp_cellprofiler/cp_cellprofiler/3.1.9) %} with the following parameters:
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