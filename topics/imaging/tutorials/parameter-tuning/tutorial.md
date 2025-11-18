---
layout: tutorial_hands_on

title: "Parameter tuning and optimization - Evaluating nuclei segmentation with Galaxy"
level: Introductory
questions:
  - What are the challenges of using the same settings for every biological image, and how does parameter tuning address these challenges?
  - How can we choose the best filters for preprocessing images for nuclei segmentation?
objectives:
  - Understand the importance of parameter tuning in bioimage analysis for achieving accurate results
  - Learn how to perform parameter tuning for segmentation in Galaxy
key_points:
  - Parameter tuning helps achieve reliable and repeatable results
time_estimation: "1H"
follow_up_training:
  -
    type: "internal"
    topic_name: imaging
    tutorials:
      - hela-screen-analysis
contributions:
  authorship:
    - rmassei
  reviewing:
    - kostrykin
    - beatrizserrano
  funding:
    - nfdi4bioimage
    - dfg
tags:
  - bioimaging
---


Parameter tuning is super important in image analysis. 
When you're doing image analysis — like segmentation, quantification, or feature extraction — the settings you choose make a big difference in how accurate your results are. 
But biological images can vary strongly in terms of quality, contrast, and structure. 
So you can't just use the same settings for every image. 

That's where parameter tuning comes in. 
By adjusting these settings carefully, researchers can make sure their tools are capturing the right biological info without adding noise or messing up the images. 
This helps make sure their work is repeatable and reliable, which is key for getting accurate scientific results. 
Basically, it's all about finding the perfect balance to get the most useful and accurate data possible.

In this tutorial, we will show how to perform parameter tuning for segmentation. In particular, we will fine-tune a parameter of Gaussian filters for nuclei segmentation by
quantitatively comparing ground-truth segmentation with the ones processed in Galaxy. As a dataset,
we will use GFP-GOWT1 mouse stem cell images, which were acquired from {% cite Brtov2011 %}.
The training dataset contains a collection of images plus their ground truth ({% cite Maka2014 %}). The images are
publicly available at the [ISBI Cell Tracking Challenge](https://celltrackingchallenge.net/2d-datasets/) repository.

![graph_comparison.png](../../images/parameter-tuning/graph_comparison.png){: width="50%"}

So, let's proceed!

## Get the image data

> <hands-on-title>Data Upload</hands-on-title>
>
> 1. Create a new history for this tutorial.
>
> 2. Download the following image and import it into your Galaxy history.
>    - [`example_image.tiff`](workflows/test-data/example_image.tiff)
>    - [`example_image_ground_truth.tiff`](workflows/test-data/example_image_ground_truth.tiff)
>    
>    If you are importing the image via URL:
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    If you are importing the image from the shared data library:
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets appropriately if needed (e.g. `"original_image"`, `"ground_truth"`)
>
> 4. Confirm the datatypes are correct (`tiff` for both images)
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
> 
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
{: .hands_on}


## Define a series of values to test

The parameter of Gaussian filters that you are going to tune in this tutorial is *sigma*, that is the size of the filters. You can now input a series of sigma values you want to test. Values need to be imported as a tabular file as
shown in the figure:

![fetch_data.jpg](../../images/parameter-tuning/fetch_data.jpg)

Next, copy and paste or manually add the desired values. Additionally, change the file format to `tabular`. Then, rename the file to `sigma_values` and click `Start`.
The dataset will appear now in your history.


You now have everything you need to build the parameter tuning workflow!

## Evaluate nuclei segmentation - Full workflow

> <hands-on-title>Create a workflow to evaluate nuclei segmentation</hands-on-title>
>
> 1. Create a new workflow in the workflow editor.
>
>    {% snippet faqs/galaxy/workflows_create_new.md %}
>
> 
> 2. Add three {% icon tool %} **Input dataset**:
> - {% icon param-file %} **1: Input Dataset**, {% icon param-file %} **2: Input Dataset** and {% icon param-file %} **3: Input Dataset** will appear in your workflow. 
> Change the "Label" of these inputs to *Dataset to perform segmentation*, *Input ground truth* and *Sigma values to test*
> 2. Add {% tool [Split by group](toolshed.g2.bx.psu.edu/repos/bgruening/split_file_on_column/tp_split_on_column/0.6) %} 
> 3. Add {% tool [Parse parameter value](param_value_from_file) %} from the list of tools with the following recommended parameters:
>   - **Select type of parameter to parse**: ***Float***
> 4. Add {% tool [Filter 2-D image](toolshed.g2.bx.psu.edu/repos/imgteam/2d_simple_filter/ip_filter_standard/1.12.0+galaxy1) %} from the list of tools with the following recommended parameters: 
>    - **Filter Type**: ***Gaussian*** 
>    - **Sigma**: Click the **"Add Connection to Module"** symbol
> 5. Add {% tool [Perform histogram equalization](toolshed.g2.bx.psu.edu/repos/imgteam/2d_histogram_equalization/ip_histogram_equalization/0.18.1+galaxy0) %} from the list of tools with the following recommended parameters:
>    - **Histogram equalization algorithm**: ***CLAHE***
> 6. Add {% tool [Threshold image](toolshed.g2.bx.psu.edu/repos/imgteam/2d_auto_threshold/ip_threshold/0.18.1+galaxy3) %} from the list of tools with the following recommended parameters:
>    - **Thresholding method** : ***Globally Adaptive/Otsu***
>    - **Offset**: 0.0
> 7. Add {% tool [Compute image segmentation and object detection performance measures](toolshed.g2.bx.psu.edu/repos/imgteam/segmetrics/ip_segmetrics/1.4.0-2) %} from the list of tools with the following recommended parameters:
>    - **Unzip**: ***No***
>    - **Segmentation is uniquely labeled**: ***No***
>    - **Ground truth is uniquely labeled**: ***No***
>    - Performance measure(s):
>      - ***Region-based / Dice***
>      - ***Region-based / Jaccard Coefficent***
>      - ***Region-based / Rand Index***
> 8. Add {% tool [Collapse Collection](toolshed.g2.bx.psu.edu/repos/nml/collapse_collections/collapse_dataset/5.1.0) %} from the list of tools with the following recommended parameters:
>    - **Prepend File name**: ***Yes***
>    - **Where to add dataset name**: ***Same line and only once per dataset***
> 9. Connect the following inputs:
>     - Connect the output of {% icon param-file %} **3: Sigma values to test** to the {% icon param-file %} *"File to split"*
>     input of {% icon tool %} **4: Split by group**.
>     - Connect the output of {% icon param-file %} **4: Split by group** to the {% icon param-file %} *"Object ID"*
>     input of {% icon tool %} **5: Input file containing parameter to parse out of**. 
>     - Connect the output of {% icon param-file %} **5: Input file containing parameter to parse out of** to the {% icon param-file %} *"Sigma"*
>     input of {% icon tool %} **6: Filter 2-D image**. 
>     - Connect the output of {% icon param-file %} **2: Dataset to perform segmentation** to the {% icon param-file %} *"Input Image"*
>     input of {% icon tool %} **6: Filter 2-D image**. 
>     - Connect the output of {% icon param-file %} **6: Filter 2-D image** to the {% icon param-file %} *"Input Image"*
>     input of {% icon tool %} **7: Perform histogram equalization**. 
>     - Connect the output of {% icon param-file %} **7: Perform histogram equalization** to the {% icon param-file %} *"Input Image"*
>     input of {% icon tool %} **8: Threshold image**. 
>     - Connect the output of {% icon param-file %} **8: Threshold image** to the {% icon param-file %} *"Segmented images"*
>     input of {% icon tool %} **9: Compute image segmentation and object detection performance measures**. 
>     - Connect the output of {% icon param-file %} **2: Input Ground Truth** to the {% icon param-file %} *"Ground truth images"*
>     input of {% icon tool %} **9: Compute image segmentation and object detection performance measures**. 
>     - Connect the output of {% icon param-file %} **9: Compute image segmentation and object detection performance measures** to the {% icon param-file %} *"Collection of files to collapse into single dataset"*
>     input of {% icon tool %} **10: Collapse Collection**.
{: .hands_on}


This is how the workflow should look like!

![workflow_overview.png](../../images/parameter-tuning/workflow_overview.png)

In this workflow, the {% tool [Parse parameter value](param_value_from_file) %} tool is passing *sigma* values one by one to the workflow downstream. In this way, a collection
of datasets is created for each sigma value and each one is processed individually.
Such strategy can be adapted for any parameter of the workflow and, in general, for any tool.

Coming back to the results, here are the results from segmentation with four different sigmas...

![sigma_comparisons.png](../../images/parameter-tuning/sigma_comparisons.png){: width="50%"}


...and below we can see the output of the analysis:

![output_workflow.png](../../images/parameter-tuning/output_workflow.png)

We are comparing three different performance measures region-based segmentation [Dice](https://en.wikipedia.org/wiki/Dice-S%C3%B8rensen_coefficient), 
[Jaccard](https://en.wikipedia.org/wiki/Jaccard_index) and [Rand](https://en.wikipedia.org/wiki/Rand_index).
All three indexes measure the agreement of the predicted segmentation image foreground and
the ground truth image foreground


Overall, the results suggest that as the *sigma* value increases too much,
the coefficient values decrease. Lower *sigma* values result in segmentations 
that more closely match the ground truth, with the highest Dice coefficient of 
0.74 achieved at a *sigma* value of 5.

## Conclusions

In this tutorial, we demonstrate how to test different *sigma* values to 
evaluate segmentation with a ground truth image. 

Overall, the {% tool [Parse parameter value](param_value_from_file) %} tool 
can be easily reused for testing other parameters and applied to different workflows.