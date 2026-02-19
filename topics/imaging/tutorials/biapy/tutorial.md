---
layout: tutorial_hands_on
title: "Execute a BiaPy workflow in Galaxy"
zenodo_link: https://zenodo.org/records/10973241
level: Intermediate
subtopic: advanced
questions:
- "What is BiaPy and how does it streamline deep learning workflows for bioimage analysis?"
- "How can we make Deep Learning (DL) models accessible to a broader audience?"
- "How can I execute a BiaPy pipeline directly within the Galaxy platform?"
- "How do I utilize pre-trained models from the BioImage.IO repository to perform inference on image data?"
objectives:
- Learn to configure and run a BiaPy workflow by editing a YAML file to define hardware settings, data paths, and model selection.
- Execute an inference workflow in Galaxy using two different pre-trained models sourced from BioImage.IO.
time_estimation: 2H
key_points:
- BiaPy is an open-source tool designed to lower the technical barriers for using DL in bioimage analysis.
- In Galaxy, BiaPy can run BioImage.IO pre-trained models and provides task-aware pre/post-processing (e.g. instance segmentation decoding) and summary statistics—beyond raw model predictions.
- The BiaPy pipeline can be controlled via a YAML configuration file, which specifies the task type and model source.
contributions:
  authorship:
    - rmassei
    - danifranco
  reviewing:
    - beatrizserrano
    - kostrykin
    - arrmunoz
tags: 
  - Image segmentation
  - Image annotation
  - Deep learning
  - Conversion
  - Overlay
  - 3D image
  - Volume rendering
---

The application of supervised and unsupervised **Deep Learning (DL)** methods in bioimage analysis have been constantly increasing in biomedical research in the last decades ({% cite esteva2021deep %}). 
DL algorithms allow automatically classifying complex biological structures by learning complex patterns and features directly from large-scale imaging data, medical scans, or high-throughput biological datasets ({% cite franco2025biapy %}). Furthermore, trained models can be easily
shared on online repositories (e.g., [BioImage.IO](https://bioimage.io/#/models)) to be reused by other scientists and support open science. 

However, running DL models often requires high-level programming skills which can often be a barrier to general audience especially the 
one without a proper computational background. Additionally, many DL models require GPU acceleration, which is not always accessible to all researchers. 
Such obstacles might limit the practical and routine adoption of DL models in bioimaging. 

So, how to make DL models accessible to a larger audience? Well, [BiaPy](https://biapyx.github.io/) is an open source framework that streamlines the use of common deep-learning workflows for a large variety of bioimage analysis tasks, including 2D and 3D semantic segmentation, instance segmentation, object detection, image denoising, single image super-resolution, self-supervised learning (for model pretraining), image classification and image-to-image translation. 

In this training, you will learn how to execute a BiaPy workflow directly in Galaxy by running [inference](https://en.wikipedia.org/wiki/Deep_learning) on a set of images using two pre-trained models from BioImage.IO defined in a 
BiaPy YAML configuration file. In particular, we will execute the CartoCell pipeline a high-content pipeline for 3D image analysis, unveils cell morphology patterns in epithelia ({% cite andres2023cartocell %}).

![example-yaml.png](../../images/biapy/example-yaml.png "Example of a BiaPy YAML file where the model with ID in BioImage.IO 'merry-water-buffalo' is defined (red box). A BiaPy YAML configuration file includes information about the hardware to be used, such as the number of CPUs or GPUs, the specific image analysis task, the model to be used, optional hyperparameters, the optimizer, and the paths for loading and storing data.")

You will perform a comparative analysis of the segmentation performance of two models from BioImage.IO, namely
[venomous-swan](https://bioimage.io/#/artifacts/venomous-swan) and [merry-water-buffalo](https://bioimage.io/#/artifacts/merry-water-buffalo).
Both model can perform cyst segmentation for fluorescence microscopy images and have the same 3D [U-Net](https://en.wikipedia.org/wiki/U-Net) + Residual Blocks base architecture. 
However, venomous-swan enhances the 3D Residual U-Net with Squeeze-and-Excitation (SE) blocks.

Our goal is to check if whether the inclusion of SE blocks in the venomous-swan model leads to improved segmentation accuracy 
compared to the merry-water-buffalo model. This comparison will help us understand the impact of the SE blocks on the model's ability to segment 
cysts in fluorescence microscopy images.

Finally, we will assess the models using various metrics as well as a qualitative assessments of the segmentation masks will be conducted to 
visually inspect the differences in segmentation quality between the two models.

Let's start with BiaPy!

> <agenda-title></agenda-title>
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Getting test data and the BiaPy YAML configuration file

The dataset required for this tutorial is available from [Zenodo]({{ page.zenodo_link }}). The CartoCell dataset contains whole epithelial cysts acquired at low resolution with minimal human intervention ([more information]({{ page.zenodo_link }})). The dataset is divided into *test*, *train* and *validation* data, each folder containing images and associated segmentation masks.

In order to simplify the upload, we already prepared the test images and YAML files in the **Data Library** that you can access on the left panel in Galaxy.

{% snippet faqs/galaxy/datasets_import_from_data_library.md %}

After importing the data from the Data Library, you should have the following files in your history: 
- `01_raw_image.tiff`
- `01_raw_mask.tiff`
- `02_raw_image.tiff`
- `02_raw_mask.tiff`
- `conf_cartocell_swam.yaml`
- `conf_cartocell_buffalo.yaml`

## Run inference using the BioImage.IO pre-trained model

Now we can set up the BiaPy tool with the ['venomous-swam' model](https://bioimage.io/#/artifacts/venomous-swan) which is defined in `conf_cartocell_swam.yaml`

> <hands-on-title>Configure the BiaPy Tool with 'venomous-swam'</hands-on-title>
>
> 1. {% tool [Build a workflow with BiaPy](toolshed.g2.bx.psu.edu/repos/iuc/biapy/biapy/3.6.5+galaxy0) %} with the following parameters to extract metadata from the image:
>
>- *Do you have a configuration file?* : `Yes, I have one and I want to run BiaPy directly`
>
>- *Select a configuration file*: `conf_cartocell_swam.yaml`
>
>- *Select the model checkpoint (if needed)* : Leave it blank. A checkpoint is a local file containing the trained model weights (e.g. .safetensors/.pth). In this tutorial we load a pre-trained model from the BioImage Model Zoo (BioImage.IO), so no local checkpoint file is required.
>
>- In the test data section, select the raw images to run predictions on and the ground truth/target images to evaluate those predictions. If no target data is provided, evaluation metrics will not be computed. **Make sure the files are in the same order so each raw image is paired with its corresponding target image**.
>
>     - *Specify the test raw images*: `01_raw_image.tiff` and `02_raw_image.tiff`
>
>     - *Specify the test ground truth/target images*: `01_raw_mask.tiff` and `02_raw_mask.tiff`
>
>- On *Select output* check the boxes:
>     - {% icon param-check %} `Test predictions (if exist)`
>     - {% icon param-check %} `Post-processed test predictions (if exist)` 
>     - {% icon param-check %} `Evaluation metrics (if exist, on test data)`
{: .hands_on}

Once the tool finishes running, you will have three different datasets in your history.

**1. Test predictions**: Full-size output images produced by the model on the test set. Because the model predicts small, overlapping patches, these patch outputs are merged back together to form one prediction per original image.

![test-prediction.png](../../images/biapy/test-prediction.png "Test predition on 22th z-stack of 02_raw_image.tiff"){: style="width:50%;" } 

**2. Post-Processed Test Prediction**: Test predictions after automatic “clean-up” steps defined in the configuration. These steps can refine the raw output (for example, removing small spurious regions or separating touching objects). 
In the YAML file definition, [Voronoi tessellation](https://en.wikipedia.org/wiki/Voronoi_diagram) is automatically applied to ensure that all instances touch each other.

![post-test-prediction.png](../../images/biapy/post-test-prediction.png "Post test predition on 22th z-stack of 02_raw_image.tiff"){: style="width:50%;" } 

**3. Test metrics:** Numerical scores that measure how well the predictions match the ground truth (if provided). 
In instance segmentation, the report typically includes:

- Intersection Over Union (IoU) per output channel (how well pixel regions overlap). This metric, also referred as the [Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index), is essentially a method to quantify the percent of overlap between the target mask and the prediction output.

- Matching metrics (how well individual predicted objects match true objects), shown for raw predictions and post-processed predictions.

...but you can find more info on the test metrics in [BiaPy documentation](https://biapy.readthedocs.io/en/latest/)!

## Visualize the results

As first step, we can visualize one slice of the segmentation on the original image. We will work with `02_raw_image.tiff`

> <hands-on-title>Extract 2D results from the BiaPy output</hands-on-title>
>
> 1. {% tool [Extract Dataset](__EXTRACT_DATASET__) %} with the following parameters:
>- Input List: '"Build a workflow with BiaPy on dataset 2, 3, and others: Post-processed test predictions"'
>- How should a dataset be selected?:
>     - Select by Index
>         - Element index: `1`
>
> 2. Rename {% icon galaxy-pencil %} the dataset to `biapy_prediction_swam.tiff`
>
> 3. {% tool [Convert image format](toolshed.g2.bx.psu.edu/repos/imgteam/bfconvert/ip_convertimage/6.7.0+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Input Image"*: `biapy_prediction_swam.tiff`
>    - *"Extract series"*: `All series`
>    - *"Extract timepoint"*: `All timepoints`
>    - *"Extract channel"*: `All channels`
>    - *"Extract z-slice"*: `Extract z-slice`
>        - *"Z-slice id"* `22`
>    - *"Extract range"*: `All images`
>    - *"Extract crop"*: `Full image`
>    - *"Tile image"*: `No tiling`
>    - *"Pyramid image"*: `No Pyramid`
>
> 4. Rename {% icon galaxy-pencil %} the dataset to `stack_biapy_prediction_swam.tiff`
>
> 5. {% tool [Convert image format](toolshed.g2.bx.psu.edu/repos/imgteam/bfconvert/ip_convertimage/6.7.0+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Input Image"*: `01_raw_image.tiff`
>    - *"Extract series"*: `All series`
>    - *"Extract timepoint"*: `All timepoints`
>    - *"Extract channel"*: `All channels`
>    - *"Extract z-slice"*: `Extract z-slice`
>        - *"Z-slice id"* `22`
>    - *"Extract range"*: `All images`
>    - *"Extract crop"*: `Full image`
>    - *"Tile image"*: `No tiling`
>    - *"Pyramid image"*: `No Pyramid`
>
> 6. Rename {% icon galaxy-pencil %} the dataset to `stack_raw_image.tiff`
>
> 7. {% tool [Overlay images](toolshed.g2.bx.psu.edu/repos/imgteam/overlay_images/ip_overlay_images/0.0.4+galaxy4) %} with the following parameters to convert the image to PNG:
>    - *"Type of the overlay"*: `Segmentation contours over image`
>    - {% icon param-file %} *"Intensity image"*: `stack_raw_image.tiff` file
>    - {% icon param-file %} *"Label map"*: `stack_biapy_prediction_swam.tiff` file (output of {% tool [Convert binary image to label map](toolshed.g2.bx.psu.edu/repos/imgteam/binary2labelimage/ip_binary_to_labelimage/0.5+galaxy0) %})
>    - *"Contour thickness"*: `1`
>    - *"Contour color"*: `green`
>    - *"Show labels"*: `no`
>8. Rename {% icon galaxy-pencil %} the dataset to `2D-overlay_biapy_swam.tiff`
{: .hands_on}

The segmentation results for the 22th z-stack are shown below:

![comparison_2D_swam.png](../../images/biapy/comparison_2D_swam.png "Segmentation results obtained with the 'venomous-swan' model. From left to right: the raw input image, predicted segmentation labels, and an overlay of the labels on the raw image."){: width="50%"}

We can also do better and visualize the full 3D segmentation using the [LibCarna](https://github.com/kostrykin/LibCarna) tool in Galaxy!

> <hands-on-title>Visual 3D</hands-on-title>
>
> 1. {% tool [Render 3-D image data](toolshed.g2.bx.psu.edu/repos/imgteam/libcarna_render/libcarna_render/0.2.0+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Input image (3-D)"*: `01_raw_image.tiff` 
>    - *"Unit of the intensity values"*: `No unit`
>    - *"Coordinate system"*: `Point Y to the Top`
>    - *"Rendering mode"*: `Maximum Intensity Projection (MIP)`
>    - *"Color map"*: `gist_gray`
>    - *"Camera parameters"*:
>      - *"Distance"*: `100`
>    - *"Render mask overlay"*:
>      - {% icon param-file %} *"Mask overlay (3-D) "*: `biapy_prediction_swam.tiff`
>    - *"Video parameters"*:
>      - *"Frames"*: `400`
{: .hands_on}

<video loop="true" autoplay="autoplay" muted width="75%">
    <source src="../../images/biapy/carto_segm.mp4" type="video/mp4"/>
</video>

Pretty cool, eh? 

We can do the same for also for `02_raw_image.tiff`:

<video loop="true" autoplay="autoplay" muted width="75%">
    <source src="../../images/biapy/carto_segm_2.mp4" type="video/mp4"/>
</video>

## Compare different pre-trained models 

Let's now run the BiaPy tool again but this time with the ['merry-water-buffalo'](https://bioimage.io/#/artifacts/merry-water-buffalo) model:

> <hands-on-title>Configure the BiaPy Tool for 'merry-water-buffalo'</hands-on-title>
>
> 1. {% tool [Build a workflow with BiaPy](toolshed.g2.bx.psu.edu/repos/iuc/biapy/biapy/3.6.5+galaxy0) %} with the following parameters to extract metadata from the image:
>
>- *Do you have a configuration file?* : `Yes, I have one and I want to run BiaPy directly`
>
>- *Select a configuration file*: `conf_cartocell_buffalo.yaml`
>
>- *Select the model checkpoint (if needed)* : Leave it blank. We will load the pre-trained model directly from the BioImage.IO, so no checkpoint file is required.
>
>- In the test data section, select the raw images to run predictions on and the ground truth/target images to evaluate those predictions. If no target data is provided, evaluation metrics will not be computed. **Make sure the files are in the same order so each raw image is paired with its corresponding target image**.
>
>     - *Specify the test raw images*: `01_raw_image.tiff` and `02_raw_image.tiff`
>
>     - *Specify the test ground truth/target images*: `01_raw_mask.tiff` and `02_raw_mask.tiff`
>
>- On *Select output* check the boxes:
>     - {% icon param-check %} `Test predictions (if exist)`
>     - {% icon param-check %} `Post-processed test predictions (if exist)` 
>     - {% icon param-check %} `Evaluation metrics (if exist, on test data)`
{: .hands_on}

We can visualize again the results using the previous approach:

> <hands-on-title>Extract the results from the BiaPy output</hands-on-title>
>
> 1. {% tool [Extract Dataset](__EXTRACT_DATASET__) %} with the following parameters:
>- Input List: '"Build a workflow with BiaPy on dataset 2, 3, and others: Post-processed test predictions"'
>- How should a dataset be selected?:
>     - Select by Index
>         - Element index: `1`
>
> 2. Rename {% icon galaxy-pencil %} the dataset to `biapy_prediction_buffalo.tiff`
>
> 3. {% tool [Convert image format](toolshed.g2.bx.psu.edu/repos/imgteam/bfconvert/ip_convertimage/6.7.0+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Input Image"*: `biapy_prediction_buffalo.tiff`
>    - *"Extract series"*: `All series`
>    - *"Extract timepoint"*: `All timepoints`
>    - *"Extract channel"*: `All channels`
>    - *"Extract z-slice"*: `Extract z-slice`
>        - *"Z-slice id"* `22`
>    - *"Extract range"*: `All images`
>    - *"Extract crop"*: `Full image`
>    - *"Tile image"*: `No tiling`
>    - *"Pyramid image"*: `No Pyramid`
>
> 4. Rename {% icon galaxy-pencil %} the dataset to `stack_biapy_prediction_buffalo.tiff`
>
> 5. {% tool [Convert image format](toolshed.g2.bx.psu.edu/repos/imgteam/bfconvert/ip_convertimage/6.7.0+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Input Image"*: `01_raw_image.tiff`
>    - *"Extract series"*: `All series`
>    - *"Extract timepoint"*: `All timepoints`
>    - *"Extract channel"*: `All channels`
>    - *"Extract z-slice"*: `Extract z-slice`
>        - *"Z-slice id"* `22`
>    - *"Extract range"*: `All images`
>    - *"Extract crop"*: `Full image`
>    - *"Tile image"*: `No tiling`
>    - *"Pyramid image"*: `No Pyramid`
>
> 6. Rename {% icon galaxy-pencil %} the dataset to `stack_raw_image.tiff`
>
> 7. {% tool [Overlay images](toolshed.g2.bx.psu.edu/repos/imgteam/overlay_images/ip_overlay_images/0.0.4+galaxy4) %} with the following parameters to convert the image to PNG:
>    - *"Type of the overlay"*: `Segmentation contours over image`
>    - {% icon param-file %} *"Intensity image"*: `stack_raw_image.tiff` file
>    - {% icon param-file %} *"Label map"*: `stack_biapy_prediction_buffalo.tiff` file (output of {% tool [Convert binary image to label map](toolshed.g2.bx.psu.edu/repos/imgteam/binary2labelimage/ip_binary_to_labelimage/0.5+galaxy0) %})
>    - *"Contour thickness"*: `1`
>    - *"Contour color"*: `green`
>    - *"Show labels"*: `no`
>8. Rename {% icon galaxy-pencil %} the dataset to `2D-overlay_biapy_buffalo.tiff`
{: .hands_on}

Results will look like this:

![comparison_2D_buffalo.png](../../images/biapy/comparison_2D_buffalo.png "Segmentation results with the 'merry-water-buffalo' model from BioImage.IO"){: width="50%"}

From a visual inspection of the results, the **'venomous-swam'** model appears to produce sharper contours and more clearly defined cells, whereas **'merry-water-buffalo'** seems better at capturing cells with fewer merges. However, its segmentation is slightly noisier.

It is hard to say that a prediction is better than other by just looking at a slice when working in 3D! 

The **Test metrics** will give us a better overview on which model is performing better!

**'venomous-swam'** struggles mainly with missing objects (low recall)

Let's evaluate segmentation at a IoU ≥ 0.5, a moderate and not-too-strict matching threshold:

    Precision: 0.328
    Recall: 0.171
    F1: 0.222
    PQ: 0.157

**water-buffalo** is clearly stronger at instance detection and segmentation

At a commonly used matching threshold (IoU ≥ 0.5), averaged across the 2 test images:

    Precision: 0.761
    Recall: 0.631
    F1: 0.689
    PQ: 0.473

So **Water-buffalo** is better both in:

- finding objects (higher recall),
- keeping predictions correct (higher precision).

## Conclusions

In this tutorial, you executed BiaPy inference pipelines directly in Galaxy using YAML configuration files, and compared two pre-trained BioImage.IO models on the same 3D dataset.

You learned how to:

- Run BiaPy in **“configuration-file driven”** mode, which makes analyses easier to reproduce and share
- Load **pre-trained BioImage.IO models** without providing a local checkpoint
- Inspect results both as a **2D slice overlay** and as a **3D rendering**
- Use the **evaluation metrics** to compare models objectively rather than relying on visual inspection alone

In our example, the two models showed different strengths: one produced cleaner contours in a slice view, while the other achieved stronger quantitative performance (higher object-level matching metrics at common IoU thresholds). This illustrates why combining **qualitative visualization** with **quantitative scoring** is important when selecting a model.

### Where to go next

- Edit the YAML to test different post-processing settings (e.g., instance separation parameters)
- Run inference on additional images or your own data (keeping image/mask pairing consistent)
- Try other BioImage.IO models and compare them using the same workflow and plots

# Optional: Extract complete training workflow from history
As an optional step, you can extract a complete workflow from your Galaxy history. This allows you to save and reuse the entire training process as a reproducible and shareable workflow.

{% snippet faqs/galaxy/workflows_extract_from_history.md %}

