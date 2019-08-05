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

The dataset required for this tutorial contains a screen of [DAPI](https://en.wikipedia.org/wiki/DAPI) stained [HeLa](https://en.wikipedia.org/wiki/HeLa) nuclei ([more information](https://zenodo.org/record/3360236)). We will use a sample image from this dataset for training basic image processing skills in Galaxy.


> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. If you are logged in, create a new history for this tutorial
>
>    {% include snippets/create_new_history.md %}
>
> 2. Import the following dataset
>    - **Important:** Choose the type of data as `zip`.
>
>    ```
>    https://zenodo.org/record/3360236/files/data.zip
>    ```
>
>    {% include snippets/import_via_link.md %}
>
> 3. **Unzip file** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"input_file"*: `Zipped ` input file
>    - *"Extract single file"*: `Single file`
>    - *"Filepath"*: `B2--W00026--P00001--Z00000--T00000--dapi.tif`
>
> 4. Rename {% icon galaxy-pencil %} the dataset to `input.tif`
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
>    - {% icon param-file %} *"Input Image"*: `input.tif` file
>    - *"Output data type"*: `PNG`
> 2. Rename {% icon galaxy-pencil %} the generated file to `viz_input`
> 3. Click on the {% icon galaxy-eye %} (eye) icon next to the file name to look at the file content
>
{: .hands_on}

Your image should look something like this:

![raw input image](../../images/imaging-introduction/viz_input.png){: width="75%"}

> ### {% icon question %} Questions
>
> You can observe that the image content is barely visible. Why?
>
> > ### {% icon solution %} Solution
> > The original image is 16-bit and the intensity values are spread over a larger range than the
> > display can render. Therefore, for improved visibility the intensity histogram of the image can
> > be normalized first.
> {: .solution }
{: .question}


Next we will normalize the histogram.

> ### {% icon hands_on %} Hands-on: Normalize Histogram and Convert Image
>
> 1. **Histogram equalization** {% icon tool %} with the following parameters to normalize the histogram of the image:
>    - {% icon param-file %} *"Source file"*: `input.tif` file
>    - *"Histogram Equalization Algorithm"*: `CLAHE`
> 2. Rename {% icon galaxy-pencil %} the generated file to `input_normalized`
> 3. **Convert image** {% icon tool %} with the following parameters to convert the image to PNG:
>    - {% icon param-file %} *"Input Image"*: `input_normalized` file (output of **Histogram equalization** {% icon tool %})
>    - *"Output data type"*: `PNG`
> 4. Rename {% icon galaxy-pencil %} the generated file to `viz_normalized`
> 5. Click on the {% icon galaxy-eye %} (eye) icon next to the file name, to look at the file content
{: .hands_on}

Your image should now look something like this:

![viz_normalized output image](../../images/imaging-introduction/viz_normalized.png){: width="75%"}

We can now clearly make out the presence of the stained nuclei. Next we will automatically detect these features and segment the image.

# Image Filtering

Specific features of interest (e.g., edges, noise) can be enhanced or suppressed by using an image filter.

> ### {% icon hands_on %} Hands-on: Filter image
>
> 1. **Filter Image** {% icon tool %} with the following parameters to smooth the image:
>    - *"Image type"*: `Gaussian Blur`
>    - *"Radius/Sigma"*: `3`
>    - {% icon param-file %} *"Source file"*: `input.tif` file
> 2. Rename {% icon galaxy-pencil %} the generated file to `input_smoothed`
> 3. **Histogram equalization** {% icon tool %} with the following parameters to normalize the histogram of the image:
>    - {% icon param-file %} *"Source file"*: `input_smoothed` file (output of **Filter image** {% icon tool %})
>    - *"Histogram Equalization Algorithm"*: `CLAHE`
> 4. Rename {% icon galaxy-pencil %} the generated file to `input_smoothed_normalized`
> 5. **Convert image** {% icon tool %} with the following parameters to convert the image to PNG:
>    - {% icon param-file %} *"Input Image"*: `input_smoothed_normalized` file (output of **Histogram equalization** {% icon tool %})
>    - *"Output data type"*: `PNG`
> 6. Rename {% icon galaxy-pencil %} the generated file to `viz_smoothed_normalized`
> 7. Click on the {% icon galaxy-eye %} (eye) icon next to the file name, to look at the file content and compare the result with `viz_normalized`. You can observe that `viz_smoothed_normalized` has significant reduced noise.
{: .hands_on}

Your image should now look something like this:

![viz_smoothed_normalized output image](../../images/imaging-introduction/viz_smoothed_normalized.png){: width="75%"}


# Segmentation

Objects of interest like nuclei can be segmented by using a smoothed image and thresholding. Moreover, the results can be overlayed with the original image.

> ### {% icon hands_on %} Hands-on: Segment image
>
> 1. **Auto Threshold** {% icon tool %} with the following parameters to segment the image:
>    - {% icon param-file %} *"Source file"*: `input_smoothed` file (output of **Filter image** {% icon tool %})
>    - *"Threshold Algorithm"*: `Otsu`
>    - *"Dark Background"*: `Yes`
> 2. Rename {% icon galaxy-pencil %} the generated file to `input_segmented`
> 3. **Binary 2 Label** {% icon tool %} with the following parameters to segment the image:
>    - {% icon param-file %} *"Binary Image File"*: `input_segmented` file (output of **Auto Threshold** {% icon tool %})
> 4. Rename {% icon galaxy-pencil %} the generated file to `input_segmented_labeled`
> 5. **Convert image** {% icon tool %} with the following parameters to convert the image to PNG:
>    - {% icon param-file %} *"Input Image"*: `input_segmented_labeled` file (output of **Binary 2 Label** {% icon tool %})
>    - *"Output data type"*: `PNG`
> 6. Rename {% icon galaxy-pencil %} the generated file to `viz_segmented`
> 7. **Overlay Segmentation Mask** {% icon tool %} with the following parameters to convert the image to PNG:
>    - {% icon param-file %} *"Image Source File"*: `viz_normalized` file
>    - {% icon param-file %} *"Mask Source File"*: `viz_segmented` file
>    - *"Image Is Greyscale"*: `Yes`
>    - *"Thickness"*: `3`
>    - *"Stroke Color"*: `red`
>    - *"Plot Labels"*: `yes`
>    - *"Label Color"*: `yellow`
> 8. Click on the {% icon galaxy-eye %} (eye) icon next to the file name, to look at the file content and assess the segmentation performance
{: .hands_on}

The resulting images should look something like this:

![viz_segmented output image](../../images/imaging-introduction/viz_segmented.png){: width="40%"}
![segmentation mask output image](../../images/imaging-introduction/viz_segmentation_mask.png){: width="40%"}

We see that we were able to detect the locations of the nuclei and segment the cells in the image.

# Conclusion
{:.no_toc}

In this exercise you imported images into Galaxy, extracted meta information from an image, converted between file formats, learned how to visualize microscopy images, filtered the image, and segmented cells using Galaxy.
