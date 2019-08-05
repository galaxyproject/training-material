---
layout: tutorial_hands_on

title: "Introduction to image analysis using Galaxy"
zenodo_link: https://zenodo.org/record/3360236
questions:
  - "How do I use Galaxy with imaging data?"
  - "How do I convert images using Galaxy?"
  - "How do I display images in Galaxy?"
  - "How do I filter images in Galaxy?"
  - "How do I segment simple images in Galaxy?"
objectives:
  - "How to handle images in Galaxy."
  - "How to perform basic image processing in Galaxy."
requirements:
time_estimation: "30M"

contributors:
  - thomaswollmann

---

# Introduction
{:.no_toc}

This tutorial shows how to use Galaxy for image analysis. The data used in this tutorial is available at [Zenodo](https://zenodo.org/record/3360236).

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Getting data

The dataset required for this tutorial contains a screen of DAPI stained HeLa nuclei ([more information](https://zenodo.org/record/3360236)). We will use a sample image from this dataset for training basic image processing skills in Galaxy.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. If you are logged in, create a new history for this tutorial
>
>    {% include snippets/create_new_history.md %}
>
> 2. Import the following dataset and choose the type of data as `zip`.
>
>    ```
>    https://zenodo.org/record/3360236/files/data.zip
>    ```
>
>    {% include snippets/import_via_link.md %}
>
> 3. "Use the unzip tool to extract a single image from the zipped screen": `unzip`
>    - {% icon param-file %} *"input_file"*: `Zipped ` input file
>    - *"Extract single file"*: `Single file`
>    - *"Filepath"*: `B2--W00026--P00001--Z00000--T00000--dapi.tif`
>
> 4. Rename dataset to `input.tif`
>
>    {% include snippets/rename_dataset.md %}
{: .hands_on}


# Image Metadata Extraction

Now, we can extract metadata from an image.

> ### {% icon hands_on %} Hands-on: Extract Image Metadata
>
> 1. **Image Info** {% icon tool %} with the following parameters to extract metadata from the image:
>    - {% icon param-file %} *"Input Image"*: `input.tif` file (output of the previous step)
> 2. Click on the {% icon galaxy-eye %} (eye) icon next to the file name, to look at the file content and search for image acquisition information
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What is the datatype?
>    > 2. What are the pixel dimentions?
>    > 3. How many bits per pixel are used?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. TIFF
>    > > 2. 1344x1024
>    > > 3. 16
>    > {: .solution }
>    {: .question}
{: .hands_on}

# Image Conversion

Not all tools can handle all image formats. Especially proprietary microscope image formats should be converted to TIFF ([supported formats](https://docs.openmicroscopy.org/bio-formats/5.7.1/supported-formats.html)). However, TIFF can not be displayed in the browser. Therefore, we convert `input.tif` to a PNG for visualization.

> ### {% icon hands_on %} Hands-on: Convert Image
>
> 1. **Convert image** {% icon tool %} with the following parameters to convert the image to PNG:
>    - {% icon param-file %} *"Input Image"*: `input.tif` file (output of the previous step)
>    - *"Output data type"*: `PNG`
> 2. Rename the generated file to `viz_input`
> 3. Click on the {% icon galaxy-eye %} (eye) icon next to the file name, to look at the file content
>
>    > ### {% icon question %} Questions
>    >
>    > 1. You can observe that the image content is barely visible. Why?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. The original image is 16-bit and the intensity values are spreaded over a larger range than the display can render. Therefore, for improved visibility the intensity histogram of the image can be normalized first.
>    > {: .solution }
>    {: .question}
{: .hands_on}

> ### {% icon hands_on %} Hands-on: Normalize Histogram and Convert Image
>
> 1. **Histogram equalization** {% icon tool %} with the following parameters to normalize the histogram of the image:
>    - {% icon param-file %} *"Source file"*: `input.tif` file (output of the previous step)
>    - *"Histogram Equalization Algorithm"*: `CLAHE`
> 2. Rename the generated file to `input_normalized`
> 3. **Convert image** {% icon tool %} with the following parameters to convert the image to PNG:
>    - {% icon param-file %} *"Input Image"*: `input_normalized` file (output of the histogram equalization)
>    - *"Output data type"*: `PNG`
> 4. Rename the generated file to `viz_normalized`
> 5. Click on the {% icon galaxy-eye %} (eye) icon next to the file name, to look at the file content
{: .hands_on}

# Image Filtering

Specific features of interest (e.g., edges, noise) can be enhanced or suppressed by using an image filter.

> ### {% icon hands_on %} Hands-on: Filter image
>
> 1. **Filter Image** {% icon tool %} with the following parameters to smooth the image:
>    - *"Image type"*: `Gaussian Blur`
>    - *"Radius/Sigma"*: `3`
>    - {% icon param-file %} *"Source file"*: `input.tif` file (output of the previous step)
> 2. Rename the generated file to `input_smoothed`
> 3. **Histogram equalization** {% icon tool %} with the following parameters to normalize the histogram of the image:
>    - {% icon param-file %} *"Source file"*: `input_smoothed` file (output of filter image)
>    - *"Histogram Equalization Algorithm"*: `CLAHE`
> 4. Rename the generated file to `input_smoothed_normalized`
> 5. **Convert image** {% icon tool %} with the following parameters to convert the image to PNG:
>    - {% icon param-file %} *"Input Image"*: `input_smoothed_normalized` file (output of the histogram equalization)
>    - *"Output data type"*: `PNG`
> 6. Rename the generated file to `viz_smoothed_normalized_normalized`
> 7. Click on the {% icon galaxy-eye %} (eye) icon next to the file name, to look at the file content and compare the result with `viz_normalized`.
{: .hands_on}

# Segmentation

Objects of interest like nuclei can be segmented by using a smoothed image and thresholding. Moreover, the results can be overlayed with the original image.

> ### {% icon hands_on %} Hands-on: Segment image
>
> 1. **Auto Threshold** {% icon tool %} with the following parameters to segment the image:
>    - {% icon param-file %} *"Source file"*: `input_smoothed` file (output of the previous step)
>    - *"Threshold Algorithm"*: `Otsu`
>    - *"Dark Background"*: `Yes`
> 2. Rename the generated file to `input_segmented`
> 3. **Convert image** {% icon tool %} with the following parameters to convert the image to PNG:
>    - {% icon param-file %} *"Input Image"*: `input_segmented` file (output of auto threshold)
>    - *"Output data type"*: `PNG`
> 4. Rename the generated file to `viz_segmented`
> 5. **Overlay Segmentation Mask** {% icon tool %} with the following parameters to convert the image to PNG:
>    - {% icon param-file %} *"Image Source File"*: `viz_normalized` file (output of previous step)
>    - {% icon param-file %} *"Mask Source File"*: `viz_segmented` file (output of convert image)
>    - *"Image Is Greyscale"*: `Yes`
>    - *"Thickness"*: `3`
>    - *"Stroke Color"*: `red`
>    - *"Plot Labels"*: `yes`
>    - *"Label Color"*: `yellow`
> 6. Click on the {% icon galaxy-eye %} (eye) icon next to the file name, to look at the file content and assess the segmentation performance
{: .hands_on}


# Conclusion
{:.no_toc}

In this exercise you imported images into Galaxy, extracted meta information from an image, converted between file formats, learned how to visualize microscopy images, filtered the image, and segmented cells using Galaxy.
