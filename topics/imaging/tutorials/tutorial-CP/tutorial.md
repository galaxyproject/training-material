---
layout: tutorial_hands_on

title: Nucleoli segmentation and feature extraction using CellProfiler
zenodo_link: ''
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
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- contributor1
- contributor2

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

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Title for your first section

Give some background about what the trainees will be doing in the section.
Remember that many people reading your materials will likely be novices,
so make sure to explain all the relevant concepts.

## Title for a subsection
Section and subsection titles will be displayed in the tutorial index on the left side of
the page, so try to make them informative and concise!

# Hands-on Sections
Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

# Title of the section usually corresponding to a big step in the analysis

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


## Sub-step with **IdentifyPrimaryObjects**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **IdentifyPrimaryObjects** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **Starting Modules** {% icon tool %})
>    - *"Use advanced settings?"*: `Yes, use advanced settings`
>        - *"Enter the name of the input image (from NamesAndTypes)"*: `DNA`
>        - *"Enter the name of the primary objects to be identified"*: `Nuclei`
>        - *"Typical minimum diameter of objects, in pixel units (Min)"*: `15`
>        - *"Typical maximum diameter of objects, in pixel units (Max)"*: `200`
>        - *"Threshold strategy"*: `Global`
>            - *"Thresholding method"*: `Otsu`
>                - *"Two-class or three-class thresholding?"*: `Two classes`
>                - *"Threshold correction factor"*: `0.9`
>        - *"Method to distinguish clumped objects"*: `Shape`
>            - *"Method to draw dividing lines between clumped objects"*: `Shape`
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

## Sub-step with **ConvertObjectsToImage**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **ConvertObjectsToImage** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **IdentifyPrimaryObjects** {% icon tool %})
>    - *"Enter the name of the input objects you want to convert to an image"*: `Nuclei`
>    - *"Enter the name of the resulting image"*: `MaskNuclei`
>    - *"Select the color format"*: `Binary (black & white)`
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

## Sub-step with **DisplayDataOnImage**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **DisplayDataOnImage** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **ConvertObjectsToImage** {% icon tool %})
>    - *"Display object or image measurements?"*: `Object`
>        - *"Enter the name of the input objects"*: `Nuclei`
>        - *"Measurement category"*: `Number`
>    - *"Enter the name of the image on which to display the measurements"*: `DNA`
>    - *"Display mode"*: `Text`
>        - *"Text color"*: `#ff0000`
>        - *"Number of decimals"*: `0`
>    - *"Name the output image that has the measurements displayed"*: `ImageDisplay`
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
> 1. **SaveImages** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **DisplayDataOnImage** {% icon tool %})
>    - *"Select the type of image to save"*: `Image`
>        - *"Saved the format to save the image(s)"*: `tiff`
>    - *"Enter the name of the image to save"*: `ImageDisplay`
>    - *"Select method for constructing file names"*: `From image filename`
>        - *"Enter the image name (from NamesAndTypes) to be used as file prefix"*: `DNA`
>        - *"Append a suffix to the image file name?"*: `Yes`
>            - *"Text to append to the image name"*: `_nucleiNumbers`
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

## Sub-step with **EnhanceOrSuppressFeatures**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **EnhanceOrSuppressFeatures** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **SaveImages** {% icon tool %})
>    - *"Enter the name of the input image"*: `DNA`
>    - *"Enter a name for the resulting image"*: `DNAdarkholes`
>    - *"Select the operation"*: `Enhance`
>        - *"Feature type"*: `Dark holes`
>            - *"Maximum hole size"*: `15`
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

## Sub-step with **MaskImage**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **MaskImage** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **EnhanceOrSuppressFeatures** {% icon tool %})
>    - *"Enter the name of the input image"*: `DNAdarkholes`
>    - *"Enter the name of the resulting image"*: `MaskDNAdarkholes`
>    - *"Use objects or an image as a mask?"*: `Objects`
>        - *"Enter the name objects to mask the input image"*: `Nuclei`
>    - *"Invert the mask?"*: `No`
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
> 1. **IdentifyPrimaryObjects** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **MaskImage** {% icon tool %})
>    - *"Use advanced settings?"*: `Yes, use advanced settings`
>        - *"Enter the name of the input image (from NamesAndTypes)"*: `MaskDNAdarkholes`
>        - *"Enter the name of the primary objects to be identified"*: `Nucleoli`
>        - *"Typical minimum diameter of objects, in pixel units (Min)"*: `2`
>        - *"Typical maximum diameter of objects, in pixel units (Max)"*: `15`
>        - *"Discard objects touching the border of the image?"*: `No`
>        - *"Threshold strategy"*: `Global`
>            - *"Thresholding method"*: `Otsu`
>                - *"Two-class or three-class thresholding?"*: `Two classes`
>                - *"Threshold correction factor"*: `0.9`
>        - *"Method to distinguish clumped objects"*: `Shape`
>            - *"Method to draw dividing lines between clumped objects"*: `Shape`
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

## Sub-step with **ConvertObjectsToImage**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **ConvertObjectsToImage** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **IdentifyPrimaryObjects** {% icon tool %})
>    - *"Enter the name of the input objects you want to convert to an image"*: `Nucleoli`
>    - *"Enter the name of the resulting image"*: `MaskNucleoli`
>    - *"Select the color format"*: `Binary (black & white)`
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

## Sub-step with **GrayToColor**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **GrayToColor** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **ConvertObjectsToImage** {% icon tool %})
>    - *"Enter the name of the resulting image"*: `CombinedMask`
>    - *"Select a color scheme"*: `RGB`
>        - *"Enter the name of the image to be colored red"*: `MaskNucleoli`
>        - *"Relative weight for the red image"*: `0.8`
>        - *"Enter the name of the image to be colored blue"*: `MaskNuclei`
>        - *"Relative weight for the blue image"*: `0.5`
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
> 1. **SaveImages** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **GrayToColor** {% icon tool %})
>    - *"Select the type of image to save"*: `Image`
>        - *"Saved the format to save the image(s)"*: `tiff`
>    - *"Enter the name of the image to save"*: `CombinedMask`
>    - *"Select method for constructing file names"*: `From image filename`
>        - *"Enter the image name (from NamesAndTypes) to be used as file prefix"*: `DNA`
>        - *"Append a suffix to the image file name?"*: `Yes`
>            - *"Text to append to the image name"*: `_combinedMask`
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

## Sub-step with **IdentifyPrimaryObjects**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **IdentifyPrimaryObjects** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **SaveImages** {% icon tool %})
>    - *"Use advanced settings?"*: `Yes, use advanced settings`
>        - *"Enter the name of the input image (from NamesAndTypes)"*: `DNA`
>        - *"Enter the name of the primary objects to be identified"*: `NucleiIncludingTouchingBorders`
>        - *"Typical minimum diameter of objects, in pixel units (Min)"*: `15`
>        - *"Typical maximum diameter of objects, in pixel units (Max)"*: `200`
>        - *"Discard objects outside the diameter range?"*: `No`
>        - *"Discard objects touching the border of the image?"*: `No`
>        - *"Threshold strategy"*: `Global`
>            - *"Thresholding method"*: `Otsu`
>                - *"Two-class or three-class thresholding?"*: `Two classes`
>                - *"Threshold correction factor"*: `0.9`
>        - *"Method to distinguish clumped objects"*: `Shape`
>            - *"Method to draw dividing lines between clumped objects"*: `Shape`
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

## Sub-step with **ConvertObjectsToImage**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **ConvertObjectsToImage** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **IdentifyPrimaryObjects** {% icon tool %})
>    - *"Enter the name of the input objects you want to convert to an image"*: `NucleiIncludingTouchingBorders`
>    - *"Enter the name of the resulting image"*: `Image_NucleiIncludingTouchingBorders`
>    - *"Select the color format"*: `Binary (black & white)`
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

## Sub-step with **ImageMath**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **ImageMath** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **ConvertObjectsToImage** {% icon tool %})
>    - *"Enter a name for the resulting image"*: `BG`
>    - *"Operation"*: `Subtract`
>        - In *"First Image"*:
>            - *"Image or measurement?"*: `Image`
>                - *"Enter the name of the first image"*: `DNA`
>        - In *"Second Image"*:
>            - *"Image or measurement?"*: `Image`
>                - *"Enter the name of the second image"*: `Image_NucleiIncludingTouchingBorders`
>    - *"Ignore the image masks?"*: `No`
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

## Sub-step with **MeasureGranularity**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **MeasureGranularity** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **ImageMath** {% icon tool %})
>    - In *"new image"*:
>        - {% icon param-repeat %} *"Insert new image"*
>            - *"Enter the name of a greyscale image to measure"*: `DNA`
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

## Sub-step with **MeasureTexture**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **MeasureTexture** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **MeasureGranularity** {% icon tool %})
>    - In *"new image"*:
>        - {% icon param-repeat %} *"Insert new image"*
>            - *"Enter the name of an image to measure"*: `DNA`
>    - *"Measure images or objects?"*: `Objects`
>        - In *"new object"*:
>            - {% icon param-repeat %} *"Insert new object"*
>                - *"Enter the names of the objects to measure"*: `Nuclei`
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
> 1. **MeasureObjectIntensity** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **MeasureTexture** {% icon tool %})
>    - In *"new image"*:
>        - {% icon param-repeat %} *"Insert new image"*
>            - *"Enter the name of an image to measure"*: `DNA`
>    - In *"new object"*:
>        - {% icon param-repeat %} *"Insert new object"*
>            - *"Enter the name of the objects to measure"*: `Nuclei`
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
> 1. **MeasureObjectSizeShape** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **MeasureObjectIntensity** {% icon tool %})
>    - In *"new object"*:
>        - {% icon param-repeat %} *"Insert new object"*
>            - *"Enter the name of the object to measure"*: `Nuclei`
>        - {% icon param-repeat %} *"Insert new object"*
>            - *"Enter the name of the object to measure"*: `Nucleoli`
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

## Sub-step with **RelateObjects**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **RelateObjects** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **MeasureObjectSizeShape** {% icon tool %})
>    - *"Parent objects"*: `Nuclei`
>    - *"Child objects"*: `Nucleoli`
>    - *"Calculate child-parent distances?"*: `Both`
>        - *"Calculate distances to other parents?"*: `No`
>    - *"Do you want to save the children with parents as a new object set?"*: `Yes`
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

## Sub-step with **MeasureImageQuality**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **MeasureImageQuality** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **RelateObjects** {% icon tool %})
>    - *"Calculate blur metrics?"*: `Yes`
>    - *"Calculate thresholds?"*: `Yes`
>        - *"Use all thresholding methods?"*: `No`
>            - In *"new threshold method"*:
>                - {% icon param-repeat %} *"Insert new threshold method"*
>                    - *"Select a thresholding method"*: `Otsu`
>                        - *"Two-class or three-class thresholding?"*: `Two classes`
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

## Sub-step with **MeasureImageAreaOccupied**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **MeasureImageAreaOccupied** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **MeasureImageQuality** {% icon tool %})
>    - In *"new area"*:
>        - {% icon param-repeat %} *"Insert new area"*
>            - *"Measure the area occupied in a binary image, or in objects?"*: `Objects`
>                - *"Enter the name of the objects to measure"*: `Nuclei`
>        - {% icon param-repeat %} *"Insert new area"*
>            - *"Measure the area occupied in a binary image, or in objects?"*: `Objects`
>                - *"Enter the name of the objects to measure"*: `Nucleoli`
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

## Sub-step with **MeasureImageIntensity**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **MeasureImageIntensity** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **MeasureImageAreaOccupied** {% icon tool %})
>    - In *"new image"*:
>        - {% icon param-repeat %} *"Insert new image"*
>            - *"Enter the name of the image to measure"*: `DNA`
>            - *"Measure the intensity only from areas enclosed by objects?"*: `No`
>        - {% icon param-repeat %} *"Insert new image"*
>            - *"Enter the name of the image to measure"*: `DNA`
>            - *"Measure the intensity only from areas enclosed by objects?"*: `Yes`
>                - *"Enter the name of the objects that the intensity will be aggregated within"*: `Nuclei`
>        - {% icon param-repeat %} *"Insert new image"*
>            - *"Enter the name of the image to measure"*: `BG`
>            - *"Measure the intensity only from areas enclosed by objects?"*: `No`
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
> 1. **ExportToSpreadsheet** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select the input CellProfiler pipeline"*: `output_pipeline` (output of **MeasureImageIntensity** {% icon tool %})
>    - *"Select the column delimiter"*: `Tab`
>    - *"Add a prefix to file names?"*: `Do not add prefix to the file name`
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
> 1. **CellProfiler** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Pipeline file"*: `output_pipeline` (output of **ExportToSpreadsheet** {% icon tool %})
>    - *"Are the input images packed into a tar archive?"*: `Yes`
>        - {% icon param-file %} *"A tarball of images"*: `output_tar` (output of **IDR Download** {% icon tool %})
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


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.