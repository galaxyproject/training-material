---
layout: tutorial_hands_on

title: Voronoi Segmentation
zenodo_link: 'https://doi.org/10.5281/zenodo.5494629'
questions:
  - How to use Galaxy for Voronoi Segmentation?
  - How should images be prepared before applying Voronoi segmentation?
  - How can Voronoi segmentation be used to analyze spatial relationships and divide an image into distinct regions based on proximity?
objectives:
  - "What Galaxy tools can I use to perform Voronoi Segmentation in Galaxy."
time_estimation: 3H
key_points:
  - Learn how to prepare images for Voronoi segmentation.
  - Learn to use Voronoi Segmentation to identify different regions in an image
contributors:
- annefou

---


# Introduction

<!-- This is a comment. -->

Voronoi segmentation is a technique used to divide an image or space into regions
based on the proximity to a set of defined points, called seeds or sites. Each 
region, known as a Voronoi cell, contains all locations that are closer to its 
seed than to any other. This approach is especially useful when analyzing spatial 
relationships, as it reveals how different areas relate in terms of distance and 
distribution. Voronoi segmentation is widely applicable for tasks where it's 
important to understand the proximity or neighborhood structure of points, such 
as organizing space, studying clustering patterns, or identifying regions of 
influence around each point in various types of data.


## Voronoi Segmentation for bioimage analysis

In bioimage analysis, Voronoi segmentation is a valuable tool for studying the 
spatial organization of cells, tissues, or other biological structures within an 
image. By dividing an image into regions around each identified cell or structure, 
Voronoi segmentation enables researchers to analyze how different cell types are 
distributed, measure distances between cells, and examine clustering patterns. This 
can provide insights into cellular interactions, tissue organization, and functional 
relationships within biological samples, such as identifying the proximity of immune 
cells to tumor cells or mapping neuron distributions within brain tissue.

## Voronoi Segmentation for Earth Observation

In Earth observation, Voronoi segmentation is used to analyze spatial patterns and distributions in satellite or aerial images. By creating regions based on proximity to specific points, such as cities, vegetation clusters, or monitoring stations, Voronoi segmentation helps in studying how features are organized across a landscape. This method is particularly useful for mapping resource distribution, analyzing urban growth, monitoring vegetation patterns, or assessing land use changes. For instance, it can help divide an area into regions of influence around weather stations or identify how different land cover types interact spatially, aiding in environmental monitoring and planning.


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get data


## Bioimage data

This tutorial will use an image dataset from the [BioImage archive](https://www.ebi.ac.uk/bioimage-archive/). This dataset is specifically prepared for training nuclear segmentation. 

The images are saved in the BioImage archive and can be uploaded to the Galaxy server with the corresponding BioImage Archive retrieval tool.

![S-BIAD634:IM276 dataset](../../images/voronoi-segmentation/BIAD634_IM276.png "An annotated fluorescence image dataset for training nuclear segmentation methods")

> <hands-on-title>Data upload with Bioimage Archive Tool</hands-on-title>
>
> 1. Create a new history for this tutorial.
>    When you log in for the first time, an empty, unnamed history is created by default. You can simply rename it.
> 
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. {% tool [FTP Link for BioImage Archive](toolshed.g2.bx.psu.edu/repos/bgruening/bia_download/bia_download/0.1.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Storage mode"*: `fire` (Storage mode is always fire)
>    - *"The path of accession"*: `S-BIAD/634/S-BIAD634`
>    > <comment-title> BioImage Archive </comment-title>
>    >
>    > This tool will upload all the files into your Galaxy history which can be very inconvenient when you have large dataset. 
>    > In that case, you can delete data files you do not plan to use for your analysis.
>    {: .comment}
>
> 3. Rename {% icon galaxy-pencil %} the file `Neuroblastoma_0.tif` to `input_image.tif`  
> 4. Check that the datatype
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to `input`
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
>
{: .hands_on}

## Earth Observation (EO) data

![Datasets of SH's AI4ER MRes Project](https://edsbook.org/_images/7d3b3ce159046d8da12d413a00c69137e4a073dcf1ee27d7cd4e33af6d93d526.png "a top-down RGB image of forest, captured by drone, aircraft or satellite.")

> <hands-on-title> EO Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial.
>    When you log in for the first time, an empty, unnamed history is created by default. You can simply rename it.
> 
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    - **Important:** If setting the type to 'Auto-detect', make sure that after upload, the datatype is set to tiff.
>
>    ```
>    https://zenodo.org/records/5494629/files/Sep_2014_RGB_602500_646500.tif
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename {% icon galaxy-pencil %} the file `Sep_2014_RGB_602500_646500.tif` to `input_image.tif`
> 4. Check that the datatype
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to `input`
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

# Data preparation

## Sub-step with **Convert image format**

> <hands-on-title> Select channel for Voronoi Segmentation </hands-on-title>
>
> 1. {% tool [Convert image format](toolshed.g2.bx.psu.edu/repos/imgteam/bfconvert/ip_convertimage/6.7.0+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Input Image"*: `output` (Input dataset)
>    - *"Extract series"*: `All series`
>    - *"Extract timepoint"*: `All timepoints`
>    - *"Extract channel"*: `Extract channel`
>        - *"Channel id"*: `{'id': 2, 'output_name': 'output'}`
>    - *"Extract z-slice"*: `All z-slices`
>    - *"Extract range"*: `All images`
>    - *"Extract crop"*: `Full image`
>    - *"Tile image"*: `No tiling`
>    - *"Pyramid image"*: `No Pyramid`
>
>    > <comment-title> Why do we need to select a single channel? </comment-title>
>    >
>    > Select a single channel from the input image. Note that some tools number channels starting from 1, while others start from 0.
>    {: .comment}
>
{: .hands_on}

## Sub-step with **Convert image format**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Convert image format](toolshed.g2.bx.psu.edu/repos/imgteam/bfconvert/ip_convertimage/6.7.0+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Input Image"*: `output` (Input dataset)
>    - *"Extract series"*: `All series`
>    - *"Extract timepoint"*: `All timepoints`
>    - *"Extract channel"*: `Extract channel`
>        - *"Channel id"*: `{'id': 2, 'output_name': 'output'}`
>    - *"Extract z-slice"*: `All z-slices`
>    - *"Extract range"*: `All images`
>    - *"Extract crop"*: `Full image`
>    - *"Tile image"*: `No tiling`
>    - *"Pyramid image"*: `No Pyramid`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Convert binary image to label map**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Convert binary image to label map](toolshed.g2.bx.psu.edu/repos/imgteam/binary2labelimage/ip_binary_to_labelimage/0.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Binary image"*: `output` (Input dataset)
>    - *"Mode"*: `Connected component analysis`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Filter 2-D image**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Filter 2-D image](toolshed.g2.bx.psu.edu/repos/imgteam/2d_simple_filter/ip_filter_standard/1.12.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input image"*: `output` (output of **Convert image format** {% icon tool %})
>    - *"Filter type"*: `Gaussian`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Compute Voronoi tessellation**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Compute Voronoi tessellation](toolshed.g2.bx.psu.edu/repos/imgteam/voronoi_tesselation/voronoi_tessellation/0.22.0+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input image"*: `output` (output of **Convert binary image to label map** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Threshold image**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Threshold image](toolshed.g2.bx.psu.edu/repos/imgteam/2d_auto_threshold/ip_threshold/0.18.1+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Input image"*: `output` (output of **Filter 2-D image** {% icon tool %})
>    - *"Thresholding method"*: `Manual`
>        - *"Threshold value"*: `{'id': 3, 'output_name': 'output'}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Process images using arithmetic expressions**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Process images using arithmetic expressions](toolshed.g2.bx.psu.edu/repos/imgteam/image_math/image_math/1.26.4+galaxy1) %} with the following parameters:
>    - *"Expression"*: `tessellation * (mask / 255) * (1 - seeds / 255)`
>    - In *"Input images"*:
>        - {% icon param-repeat %} *"Insert Input images"*
>            - {% icon param-file %} *"Image"*: `result` (output of **Compute Voronoi tessellation** {% icon tool %})
>            - *"Variable for representation of the image within the expression"*: `tessellation`
>        - {% icon param-repeat %} *"Insert Input images"*
>            - {% icon param-file %} *"Image"*: `output` (Input dataset)
>            - *"Variable for representation of the image within the expression"*: `seeds`
>        - {% icon param-repeat %} *"Insert Input images"*
>            - {% icon param-file %} *"Image"*: `output` (output of **Threshold image** {% icon tool %})
>            - *"Variable for representation of the image within the expression"*: `mask`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Colorize label map**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Colorize label map](toolshed.g2.bx.psu.edu/repos/imgteam/colorize_labels/colorize_labels/3.2.1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input image (label map)"*: `result` (output of **Process images using arithmetic expressions** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
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

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.


