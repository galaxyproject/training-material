---
layout: tutorial_hands_on

title: "Introduction to Image Analysis using Galaxy"
zenodo_link: https://zenodo.org/record/3362976
level: Introductory
questions:
  - "How do I use Galaxy with imaging data?"
  - "How do I convert images using Galaxy?"
  - "How do I display images in Galaxy?"
  - "How do I filter images in Galaxy?"
  - "How do I segment simple images in Galaxy?"
objectives:
  - "How to handle images in Galaxy."
  - "How to perform basic image processing in Galaxy."
key_points:
- The **Image Info** tool can provide valuable metadata information of an image.
- TIFF files cannot viewed directly in most web browser, so a visualization plugin must be used.
- For visualization, images with a bit-depth more than 8-bit have to be histogram equalized.
time_estimation: "1H"
follow_up_training:
  -
    type: "internal"
    topic_name: imaging
    tutorials:
      - hela-screen-analysis
contributions:
  authorship:
    - thomaswollmann
    - shiltemann
    - kostrykin
  funding:
    - elixir-europe
tags:
  - HeLa

---


Image analysis is the extraction of meaningful information from images by means of digital image processing techniques. Imaging is an important component in a wide range of scientific fields of study, such as astronomy, medicine, physics, biology, geography, chemistry, robotics, and industrial manufacturing.

This tutorial shows how to use Galaxy to perform basic image analysis tasks such as format conversion, image enhancement, segmentation, and feature extraction.

> <agenda-title></agenda-title>
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Getting Data

The dataset required for this tutorial is available from [Zenodo](https://zenodo.org/record/3362976) and
contains a screen of [DAPI](https://en.wikipedia.org/wiki/DAPI) stained [HeLa](https://en.wikipedia.org/wiki/HeLa) nuclei ([more information](https://zenodo.org/record/3362976)). We will use a sample image from this dataset for training basic image processing skills in Galaxy.

Our objective is to automatically count the number of cells contained in this image. In order to achieve this, we will enhance the quality of the image, automatically detect the nuclei and segment the nuclei and count them.


> <hands-on-title>Data upload</hands-on-title>
>
> 1. If you are logged in, create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the following dataset from [Zenodo](https://zenodo.org/record/3362976) or from the data library (ask your instructor).
>    - **Important:** Choose the type of data as `zip`.
>
>    ```
>    https://zenodo.org/record/3362976/files/B2.zip
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. {% tool [Unzip](toolshed.g2.bx.psu.edu/repos/imgteam/unzip/unzip/6.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"input_file"*: `B2.zip`
>    - *"Extract single file"*: `Single file`
>    - *"Filepath"*: `B2--W00026--P00001--Z00000--T00000--dapi.tif`
>
> 4. Rename {% icon galaxy-pencil %} the dataset to `input.tif`
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
{: .hands_on}


# Image Metadata Extraction

Now, we can extract metadata from an image.

> <hands-on-title>Extract Image Metadata</hands-on-title>
>
> 1. {% tool [Show image info](toolshed.g2.bx.psu.edu/repos/imgteam/image_info/ip_imageinfo/5.7.1+galaxy1) %} with the following parameters to extract metadata from the image:
>    - {% icon param-file %} *"Input Image"*: `input.tif` file (output of the previous step)
> 2. Click on the {% icon galaxy-eye %} (eye) icon next to the file name, to look at the file content and search for image acquisition information
>
>    > <question-title></question-title>
>    >
>    > 1. What is the datatype?
>    > 2. What are the pixel dimentions?
>    > 3. How many bits per pixel are used?
>    >
>    > > <solution-title></solution-title>
>    > > 1. TIFF
>    > > 2. 1344x1024
>    > > 3. 16
>    > {: .solution }
>    {: .question}
{: .hands_on}

# Visual Inspection of TIFF Images

Not all tools can handle all image formats. Especially proprietary microscope image formats should be converted to TIFF ([supported formats](https://docs.openmicroscopy.org/bio-formats/5.7.1/supported-formats.html)). However, TIFF cannot be displayed in most web browsers directly. Therefore, for visual inspection of TIFF images, we use a TIFF visualization plugin in Galaxy.

> <hands-on-title>Visual Inspection of TIFF Images</hands-on-title>
>
> 1. Click on the title of your file to see the row of small icons for saving, linking, etc.:
> ![Screenshot of Galaxy icons. Seven small blue icons are shown on a green background. From left to right they are: floppy disk, link, information, redo, bar chart, flow chart and a question mark.](../../images/imaging-introduction/LittleJobIcons.png)
> 2. Click on the **visualise icon** {% icon galaxy-visualise %} and then select the **Tiff Viewer** visualization plugin.
>
{: .hands_on}

Your image should look something like this:

![raw input image](../../images/imaging-introduction/viz_input.png){: width="75%"}

> <question-title></question-title>
>
> You can observe that the image content is barely visible. Why?
>
> > <solution-title></solution-title>
> > The original image is 16-bit and the intensity values are spread over a larger range than the
> > display can render. Therefore, for improved visibility the intensity histogram of the image can
> > be normalized first.
> {: .solution }
{: .question}


Next we will normalize the histogram to improve the contrast. We do this using a [Contrast Limited Adaptive Histogram Equalization (CLAHE)](https://en.wikipedia.org/wiki/Adaptive_histogram_equalization) approach.

> <hands-on-title>Normalize Histogram</hands-on-title>
>
> 1. {% tool [Perform histogram equalization](toolshed.g2.bx.psu.edu/repos/imgteam/2d_histogram_equalization/ip_histogram_equalization/0.18.1+galaxy0) %} with the following parameters to normalize the histogram of the image:
>    - {% icon param-file %} *"Input image"*: `input.tif` file
>    - *"Histogram equalization algorithm"*: `CLAHE`
> 2. Rename {% icon galaxy-pencil %} the generated file to `input_normalized`.
> 3. Click on the **visualise icon** {% icon galaxy-visualise %} of the file to visually inspect the image using the **Tiff Viewer** visualization plugin.
{: .hands_on}

Your image should now look something like this:

![viz_normalized output image](../../images/imaging-introduction/viz_normalized.png){: width="75%"}

We can now clearly make out the presence of the stained nuclei. Next we will automatically detect these features and segment the image.

# Image Filtering

Specific features of interest (e.g., edges, noise) can be enhanced or suppressed by using an image filter.

> <hands-on-title>Filter image</hands-on-title>
>
> 1. {% tool [Filter 2-D image](toolshed.g2.bx.psu.edu/repos/imgteam/2d_simple_filter/ip_filter_standard/1.12.0+galaxy1) %} with the following parameters to smooth the image:
>    - *"Filter type"*: `Gaussian`
>    - *"Sigma"*: `3`
>    - {% icon param-file %} *"Source file"*: `input.tif` file
> 2. Rename {% icon galaxy-pencil %} the generated file to `input_smoothed`.
> 3. {% tool [Perform histogram equalization](toolshed.g2.bx.psu.edu/repos/imgteam/2d_histogram_equalization/ip_histogram_equalization/0.18.1+galaxy0) %} with the following parameters to normalize the histogram of the image:
>    - {% icon param-file %} *"Input image"*: `input_smoothed` file (output of {% tool [Filter 2-D image](toolshed.g2.bx.psu.edu/repos/imgteam/2d_simple_filter/ip_filter_standard/1.12.0+galaxy1) %})
>    - *"Histogram equalization algorithm"*: `CLAHE`
> 4. Rename {% icon galaxy-pencil %} the generated file to `input_smoothed_normalized`.
> 5. Click on the **visualise icon** {% icon galaxy-visualise %} of the file to visually inspect the image and compare the result with `input_normalized`. You can observe that `input_smoothed_normalized` has significantly reduced noise.
{: .hands_on}

Your image should now look something like this:

![viz_smoothed_normalized output image](../../images/imaging-introduction/viz_smoothed_normalized.png){: width="75%"}


# Segmentation

Objects of interest like nuclei can be segmented by using a smoothed image and thresholding. Moreover, the results can be overlayed with the original image.

> <hands-on-title>Segment image</hands-on-title>
>
> 1. {% tool [Threshold image](toolshed.g2.bx.psu.edu/repos/imgteam/2d_auto_threshold/ip_threshold/0.18.1+galaxy3) %} with the following parameters to segment the image:
>    - {% icon param-file %} *"Input image"*: `input_smoothed` file (output of {% tool [Filter 2-D image](toolshed.g2.bx.psu.edu/repos/imgteam/2d_simple_filter/ip_filter_standard/1.12.0+galaxy1) %})
>    - *"Thresholding method"*: `Globally adaptive / Otsu`
> 2. Rename {% icon galaxy-pencil %} the generated file to `input_segmented`.
> 3. {% tool [Convert binary image to label map](toolshed.g2.bx.psu.edu/repos/imgteam/binary2labelimage/ip_binary_to_labelimage/0.5+galaxy0) %} with the following parameters to segment the image:
>    - {% icon param-file %} *"Binary image"*: `input_segmented` file (output of {% tool [Threshold image](toolshed.g2.bx.psu.edu/repos/imgteam/2d_auto_threshold/ip_threshold/0.18.1+galaxy3) %})
> 4. Rename {% icon galaxy-pencil %} the generated file to `input_segmented_labeled`
>
>    > <question-title></question-title>
>    >
>    > 1. What does {% tool [Convert binary image to label map](toolshed.g2.bx.psu.edu/repos/imgteam/binary2labelimage/ip_binary_to_labelimage/0.5+galaxy0) %} do? (Hint: check the tool help section)
>    > 2. View the `input_segmented_labeled` image from the last step, what do you see? Can you explain this result?
>    > 3. Exercise: Try to make the information in this image better visible (Hint: use {% tool [Perform histogram equalization](toolshed.g2.bx.psu.edu/repos/imgteam/2d_histogram_equalization/ip_histogram_equalization/0.18.1+galaxy0) %})
>    >
>    > > <solution-title></solution-title>
>    > > 1. The tool assigns each connected component (e.g., segmented cell) in the image a unique object ID called *label* and stores it as the intensity value.
>    > > 2. The image looks completely black.
>    > >    The object labels generated by {% tool [Convert binary image to label map](toolshed.g2.bx.psu.edu/repos/imgteam/binary2labelimage/ip_binary_to_labelimage/0.5+galaxy0) %} are relatively low. Since the labels are stored as intensity values, these are too low to be visible in this case. Nevertheless, there is more information in this image than meets the eye.
>    > > 3. To make the labeled objects visible, the values have to be stretched to a larger range of visible intensity values. We
>    > >    can do that by equalizing the histogram again. To this end, use {% tool [Perform histogram equalization](toolshed.g2.bx.psu.edu/repos/imgteam/2d_histogram_equalization/ip_histogram_equalization/0.18.1+galaxy0) %} with the following parameters to normalize the intensity values:
>    > >    - {% icon param-file %} *"Input image"*: `input_segmented_labeled` file (output of {% tool [Convert binary image to label map](toolshed.g2.bx.psu.edu/repos/imgteam/binary2labelimage/ip_binary_to_labelimage/0.5+galaxy0) %})
>    > >    - *"Histogram equalization algorithm"*: `CLAHE`
>    > >
>    > >    The information contained in the original image has now become visible to the human eye:
>    > >    ![normalized viz_segmented file](../../images/imaging-introduction/viz_segmented.png)
>    > >
>    > {: .solution }
>    {: .question}
>
>
> 7. {% tool [Overlay images](toolshed.g2.bx.psu.edu/repos/imgteam/overlay_images/ip_overlay_images/0.0.4+galaxy4) %} with the following parameters to convert the image to PNG:
>    - *"Type of the overlay"*: `Segmentation contours over image`
>    - {% icon param-file %} *"Intensity image"*: `input_normalized` file
>    - {% icon param-file %} *"Label map"*: `input_segmented_labeled` file (output of {% tool [Convert binary image to label map](toolshed.g2.bx.psu.edu/repos/imgteam/binary2labelimage/ip_binary_to_labelimage/0.5+galaxy0) %})
>    - *"Contour thickness"*: `2`
>    - *"Contour color"*: `red`
>    - *"Show labels"*: `yes`
>    - *"Label color"*: `yellow`
> 8. Click on the {% icon galaxy-eye %} (eye) icon next to the file name, to look at the file content and assess the segmentation performance.
> 9. {% tool [Count objects in label map](toolshed.g2.bx.psu.edu/repos/imgteam/count_objects/ip_count_objects/0.0.5-2) %} with the following parameters to count the segmented objects in the image:
>    - {% icon param-file %} *"Source file"*: `input_segmented_labeled` file (output of {% tool [Convert binary image to label map](toolshed.g2.bx.psu.edu/repos/imgteam/binary2labelimage/ip_binary_to_labelimage/0.5+galaxy0) %})
>
>    > <question-title></question-title>
>    >
>    > How many objects were segmented?
>    >
>    > > <solution-title></solution-title>
>    > >  The {% tool [Count objects in label map](toolshed.g2.bx.psu.edu/repos/imgteam/count_objects/ip_count_objects/0.0.5-2) %} tool counted 425 objects.
>    > {: .solution }
>    {: .question}
{: .hands_on}

The resulting image should look something like this:

![segmentation mask output image](../../images/imaging-introduction/viz_segmentation_mask.png){: width="75%"}

We see the segmentation mask overlayed; each detected object (nucleus) is labeled with its ID value.

We see that with the help of just a few simple steps, we were able to detect the locations of the stained nuclei, and count them.

# Conclusion


In this exercise you imported images into Galaxy, extracted meta information from an image, learned how to visualize microscopy images, filtered the image, and segmented cells using Galaxy.
