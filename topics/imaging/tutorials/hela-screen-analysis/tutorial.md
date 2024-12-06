---
layout: tutorial_hands_on

title: "Analyse HeLa fluorescence siRNA screen"
zenodo_link: https://zenodo.org/record/3362976
level: Intermediate
questions:
  - "How do I analyze a HeLa fluorescence siRNA screen?"
  - "How do I segment cell nuclei?"
  - "How do I extract features from segmentations?"
  - "How do I filter segmentations by morphological features?"
  - "How do I apply a feature extraction workflow to a screen?"
  - "How do I visualize feature extraction results?"
objectives:
  - "How to segment cell nuclei in Galaxy."
  - "How to extract features from segmentations in Galaxy."
  - "How to filter segmentations by morphological features in Galaxy."
  - "How to extract features from an imaging screen in Galaxy."
  - "How to analyse extracted features from an imaging screen in Galaxy."
key_points:
- Galaxy workflows can be used to scale image analysis pipelines to whole screens.
- Segmented objects can be filtered using the **Filter label map by rules** tool.
- Galaxy charts can be used to compare features extracted from screens showing cells with different treatments.
requirements:
  -
    type: "internal"
    topic_name: imaging
    tutorials:
      - imaging-introduction
follow_up_training:
  -
    type: "internal"
    topic_name: statistics
    tutorials:
      - machinelearning
time_estimation: "1H"
contributions:
  authorship:
    - thomaswollmann
    - kostrykin
  funding:
    - elixir-europe
tags:
  - HeLa

---


This tutorial shows how to segment and extract features from cell nuclei Galaxy for image analysis. As example use case, this tutorial shows you how to compare the phenotypes of PLK1 threated cells in comparison to a control. The data used in this tutorial is available at [Zenodo](https://zenodo.org/record/3362976).

RNA interference (RNAi) is used in the example use case for silencing genes by way of mRNA degradation. Gene knockdown by this method is achieved by introducing small double-stranded interfering RNAs (siRNA) into the cytoplasm. Small interfering RNAs can originate from inside the cell or can be exogenously introduced into the cell. Once introduced into the cell, exogenous siRNAs are processed by the RNA-induced silencing complex (RISC).The siRNA is complementary to the target mRNA to be silenced, and the RISC uses the siRNA as a template for locating the target mRNA. After the RISC localizes to the target mRNA, the RNA is cleaved by a ribonuclease. RNAi is widely used as a laboratory technique for genetic functional analysis. RNAi in organisms such as C. elegans and Drosophila melanogaster provides a quick and inexpensive means of investigating gene function. Insights gained from experimental RNAi use may be useful in identifying potential therapeutic targets, drug development, or other applications. RNA interference is a very useful research tool, allowing investigators to carry out large genetic screens in an effort to identify targets for further research related to a particular pathway, drug, or phenotype.

The example used in this tutorial deals with PLK1 knocked down cells. PLK1 is an early trigger for G2/M transition. PLK1 supports the functional maturation of the centrosome in late G2/early prophase and establishment of the bipolar spindle. PLK1 is being studied as a target for cancer drugs. Many colon and lung cancers are caused by K-RAS mutations. These cancers are dependent on PLK1.

> <agenda-title></agenda-title>
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Getting data

The dataset required for this tutorial contains a screen of DAPI stained HeLa nuclei ([more information](https://zenodo.org/record/3360236)). We will use a sample image from this dataset for training basic image processing skills in Galaxy.

> <hands-on-title>Data upload</hands-on-title>
>
> 1. If you are logged in, create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import {% icon galaxy-upload %} the following dataset from [Zenodo]( https://zenodo.org/record/3362976) or from the data library (ask your instructor).
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
> 4. Rename {% icon galaxy-pencil %} the dataset to `testinput.tif`
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
> 5. {% tool [Unzip](toolshed.g2.bx.psu.edu/repos/imgteam/unzip/unzip/6.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"input_file"*: `B2.zip`
>    - *"Extract single file"*: `All files`
>
> 6. Rename {% icon galaxy-pencil %} the resulting collection to `control`
>
>    {% snippet faqs/galaxy/collections_rename.md %}
>
> 7. Import {% icon galaxy-upload %} the following dataset from [Zenodo]( https://zenodo.org/record/3362976) or from the data library (ask your instructor).
>    - **Important:** Choose the type of data as `zip`.
>    ```
>    https://zenodo.org/record/3362976/files/B3.zip
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 8. {% tool [Unzip](toolshed.g2.bx.psu.edu/repos/imgteam/unzip/unzip/6.0+galaxy0) %} to extract the zipped screen:
>    - {% icon param-file %} *"input_file"*: `B3.zip`
>    - *"Extract single file"*: `All files`
>
> 9. Rename {% icon galaxy-pencil %} the collection to `PLK1`
> 9. Upload {% icon galaxy-upload %} the following segmentation filter rules as a new pasted file (format: `tabular`):
>    ```
>    	area	eccentricity
>    min	500	0.
>    max	100000	0.5
>    ```
>
>    {% snippet faqs/galaxy/datasets_create_new_file.md format="tabular" %}
>
> 9. Rename {% icon galaxy-pencil %} dataset to `rules`
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
{: .hands_on}


# Create feature extraction workflow
First, we will create and test a workflow which extracts mean DAPI intensity, area, and major axis length of cell nuclei from an image.

> <hands-on-title>Create feature extraction workflow</hands-on-title>
>
> 1. {% tool [Filter 2-D image](toolshed.g2.bx.psu.edu/repos/imgteam/2d_simple_filter/ip_filter_standard/1.12.0+galaxy1) %} with the following parameters to smooth the image:
>    - {% icon param-file %} *"Input image"*: `testinput.tif` file
>    - *"Filter type"*: `Gaussian`
>    - *"Sigma"*: `3`
> 2. {% tool [Threshold image](toolshed.g2.bx.psu.edu/repos/imgteam/2d_auto_threshold/ip_threshold/0.18.1+galaxy3) %} with the following parameters to segment the image:
>    - {% icon param-file %} *"Input image"*: output of {% tool [Filter 2-D image](toolshed.g2.bx.psu.edu/repos/imgteam/2d_simple_filter/ip_filter_standard/1.12.0+galaxy1) %}
>    - *"Thresholding method"*: `Globally adaptive / Otsu`
> 3. {% tool [Convert binary image to label map](toolshed.g2.bx.psu.edu/repos/imgteam/binary2labelimage/ip_binary_to_labelimage/0.5+galaxy0) %} with the following parameters to split touching objects:
>    - {% icon param-file %} *"Binary image"*: output of {% tool [Threshold image](toolshed.g2.bx.psu.edu/repos/imgteam/2d_auto_threshold/ip_threshold/0.18.1+galaxy3) %}
>    - *"Mode":* `Watershed transform`
>    - *"Minimum distance between two objects"*: `20`
> 4. {% tool [Extract image features](toolshed.g2.bx.psu.edu/repos/imgteam/2d_feature_extraction/ip_2d_feature_extraction/0.18.1+galaxy0) %} with the following parameters to extract features from the segmented objects:
>    - {% icon param-file %} *"Label map"*: output of {% tool [Convert binary image to label map](toolshed.g2.bx.psu.edu/repos/imgteam/binary2labelimage/ip_binary_to_labelimage/0.5+galaxy0) %}
>    - *"Use the intensity image to compute additional features"*: `No intensity image`
>    - *"Select features to compute"*: `Select features`
>    - *"Available features"*:
>        - {% icon param-check %} `Label from the label map`
>        - {% icon param-check %} `Area`
>        - {% icon param-check %} `Eccentricity`
>        - {% icon param-check %} `Major axis length`
> 5. {% tool [Filter label map by rules](toolshed.g2.bx.psu.edu/repos/imgteam/2d_filter_segmentation_by_features/ip_2d_filter_segmentation_by_features/0.0.1-4) %} with the following parameters to filter the label map from 3. with the extracted features and a set of rules:
>    - {% icon param-file %} *"Label map"*: output of {% tool [Convert binary image to label map](toolshed.g2.bx.psu.edu/repos/imgteam/binary2labelimage/ip_binary_to_labelimage/0.5+galaxy0) %}
>    - {% icon param-file %} *"Features"*: output of {% tool [Extract image features](toolshed.g2.bx.psu.edu/repos/imgteam/2d_feature_extraction/ip_2d_feature_extraction/0.18.1+galaxy0) %}
>    - {% icon param-file %} *"Rules"*: `rules` file
> 6. {% tool [Extract image features](toolshed.g2.bx.psu.edu/repos/imgteam/2d_feature_extraction/ip_2d_feature_extraction/0.18.1+galaxy0) %} with the following parameters to extract features the final readout from the segmented objects:
>    - {% icon param-file %} *"Label map"*: output of {% tool [Filter label map by rules](toolshed.g2.bx.psu.edu/repos/imgteam/2d_filter_segmentation_by_features/ip_2d_filter_segmentation_by_features/0.0.1-4) %}
>    - *"Use the intensity image to compute additional features"*: `Use intensity image`
>    - {% icon param-file %} *"Intensity image"*: `testinput.tif` file
>    - *"Select features to compute"*: `Select features`
>    - *"Available features"*:
>      - {% icon param-check %} `Mean Intensity (requires original image)`
>      - {% icon param-check %} `Area`
>      - {% icon param-check %} `Major axis length`
> 7. Now we can extract the workflow for batch processing:
>    - Name it "feature_extraction".
>    - Don't treat `B2.zip` and `B3.zip` as inputs (the workflow is supposed to be applied to the images directly).
>    - Exclude {% tool [Unzip](toolshed.g2.bx.psu.edu/repos/imgteam/unzip/unzip/6.0+galaxy0) %} by unchecking the tool (3 times).
>
>    {% snippet faqs/galaxy/workflows_extract_from_history.md %}
>
> 8. Edit the workflow you just created:
>    - Select "Input dataset" from the list of tools. The step {% icon param-file %} **8: Input Dataset** appears.
>    - Change the "Label" of {% icon param-file %} **8: Input Dataset** to `input image`.
>    - Change the "Label" of {% icon param-file %} **1: rules** to `filter rules`.
>    - Connect the output of {% icon param-file %} **8: input image** to the input of {% icon tool %} **2: Filter 2-D image**.
>    - Connect the output of {% icon param-file %} **8: input image** to the "Intensity image" input of {% icon tool %} **7: Extract image features**.
>    - Mark the results of {% icon tool %} **6: Filter label map by rules** and {% icon tool %} **7: Extract image features** as the primary outputs of the workflow (by clicking on the checkboxes of the outputs).
>
{: .hands_on}

The resulting workflow should look something like this:

![feature extraction workflow](../../images/hela-screen-analysis/feature_extraction_workflow.png "Feature extraction workflow.")

# Apply workflow to screen

Now we want to apply our extracted workflow to a series of images and merge the results. For this purpose, we create a workflow which uses the previously created workflow as a sub-workflow.

> <hands-on-title>Create screen analysis workflow</hands-on-title>
>
> 1. Create a new workflow in the workflow editor.
>
>    {% snippet faqs/galaxy/workflows_create_new.md %}
>
> 2. Select "Input dataset collection" from the list of tools. The step {% icon param-collection %} **1: Input Dataset Collection** appears in your workflow. Change the "Label" of this step to `input images`.
> 3. Add the input dataset {% icon param-file %} **2: rules** to your workflow (select "Input dataset" from the list of tools and set the "Label" of the newly created step to `rules`).
> 4. Add the {% icon workflow %} **feature_extraction** workflow as a sub-workflow:
>    - Expand the "Workflows" section in the list of tools and select "feature_extraction" to add it to the workflow.
>    - Connect the output of {% icon param-file %} **1: input images** to the "input image" input of {% icon workflow %} **3: feature_extraction**.
>    - Connect the output of {% icon param-file %} **2: rules** to the "filter rules" input of {% icon workflow %} **3: feature_extraction**.
> 5. Create the step {% icon tool %} **4: Collapse Collection** in the workflow (by choosing "Collapse Collection" from the list of tools).
>    - Connect the output "output (tabular)" of {% icon workflow %} **3: feature_extraction** to {% icon tool %} **4: Collapse Collection**.
>    - Set *"Keep one header line"* of {% icon tool %} **4: Collapse Collection**: `Yes`
>    - Set *"Prepend File name"* of {% icon tool %} **4: Collapse Collection**: `No`
>    - Mark the output of {% icon tool %} **4: Collapse Collection** as the primary workflow output.
> 6. Save your workflow and name it `analyze_screen`.
{: .hands_on}

The resulting workflow should look something like this:

![screen analysis workflow](../../images/hela-screen-analysis/analyze_screen_workflow.png "Full screen analysis workflow.")


> <hands-on-title>Run screen analysis workflow</hands-on-title>
>
> 1. Run the "analyze_screen" workflow on the `control` screen using the `rules` file.
>
>    {% snippet faqs/galaxy/workflows_run.md %}
>
> 2. Run the "analyze_screen" workflow on the `PLK1` screen using the `rules` file.
>
{: .hands_on}

# Plot feature extraction results

Finally, we want to plot the results for better interpretation.

> <hands-on-title>Plot feature extraction results</hands-on-title>
>
> 1. Click on the **Visualize** {% icon galaxy-barchart %} icon of the {% icon tool %} **4: Collapse Collection** results.
> 2. Run **Box plot (jqPlot)** with the following parameters:
>    - *"Provide a title"*: `Screen features`
>    - *"X-Axis label"*:
>    - *"Y-Axis label"*:
>    - *"1: Data series"*:
>        - *"Provide a label"*: `Mean intensity`
>        - *"Observations"*: `Column 1`
>    - *"2: Data series"*:
>        - *"Provide a label"*: `Area`
>        - *"Observations"*: `Column 2`
>    - *"3: Data series"*:
>        - *"Provide a label"*: `Major axis length`
>        - *"Observations"*: `Column 3`
>
>    > <question-title></question-title>
>    >
>    > Plot the feature distribution of PLK1 and control. What differences do you observe between the screens?
>    >
>    > > <solution-title></solution-title>
>    > > The phenotype of PLK1 threated cells show a higher mean intensity and a shorter major axis in comparison to the control.
>    > {: .solution }
>    {: .question}
{: .hands_on}


One of the resulting plots should look something like this:

![feature extraction results box plot](../../images/hela-screen-analysis/result_boxplot.png){: width="100%"}

# Conclusion


In this exercise you imported images into Galaxy, segmented cell nuclei, filtered segmentations by morphological features, extracted features from segmentations, scaled your workflow to a whole screen, and plotted the feature extraction results using Galaxy.
