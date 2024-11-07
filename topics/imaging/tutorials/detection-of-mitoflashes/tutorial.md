---
layout: tutorial_hands_on

title: Tracking of mitochondria and capturing mitoflashes
zenodo_link: ''
questions:
- What are mitoflashes, and why are they relevant for understanding mitochondrial function?
- How can bioinformatics and image analysis tools help in tracking and analyzing mitochondrial dynamics?
objectives:
- Understand the biological significance of mitoflashes and their implications for cellular health.
- Learn to track mitochondrial movements in live-cell imaging data.
- Use image analysis tools to detect and quantify mitoflashes.
time_estimation: 1H
key_points:
- Mitoflashes are superoxide burst events in mitochondria that signal functional states and stress responses.
- Bioinformatics tools can facilitate the visualization, tracking, and quantitative analysis of these events.
- Proper tracking and analysis of mitoflashes contribute to understanding mitochondrial behavior and health in cells.
contributors:
- Diana Chiang Jurado
- Leonid Kostrykin

---

# Introduction

Mitochondria act like fuel stations for the cell, supplying the energy needed to keep it functioning and healthy. During certain activities, they can produce bursts of superoxide, known as 'mitoflashes.' These short, intense events occur in individual mitochondria and can be observed in live heart cells and other cell types using a confocal microscope. Commonly detected with a mitochondria-targeted circularly permuted fluorescent protein (mt-cpYFP), mitoflashes provide real-time insights into mitochondrial respiration function in situ and act as a biosensor for superoxide levels, reflecting the activity of the mitochondrial electron transport chain.

Mitoflashes are triggered by an increase in basal reactive oxygen species (ROS) levels in mitochondria. This burst event activates the mitochondrial permeability transition pore (mPTP), causing a depolarization of the mitochondrial membrane potential (ΔΨm) and a mild alkalization of the matrix. Although mitoflashes can indicate changes in ROS, they do not precisely quantify ROS levels. These events can reveal active mitochondria from static ones, becoming mitoflashes great biomarkers for tracking mitochondrial health and function.

The frequency and kinetics of mitoflashes hold significant physiological and pathophysiological implications. They correlate with key processes such as muscle contraction, cell differentiation, neuron development, wound healing, and lifespan prediction. The frequency and characteristics of mitoflashes vary by cell type; for example, adult cardiomyocytes display approximately 3.8 ± 0.5 mitoflashes, while primary cultured hippocampal neurons show around 31 ± 4 per cell.

In this tutorial, you will learn to track mitochondria in live-cell imaging data and detect mitoflashes using specialized Galaxy tools for image analysis. You will work with time-lapse microscopy data, often stored as TIFF files with image stacks, to observe and analyze mitochondrial events across multiple time points. By identifying these events and quantifying their frequency through intensity measurements fitted to a curve, you'll gain insights into mitochondrial behavior and activity over time.

> In this tutorial, we will cover:
>
> 1. Understanding the concept and biological relevance of mitoflashes.
> 2. Processing live-cell imaging data for mitochondrial tracking.
> 3. Detecting and analyzing mitoflashes with Galaxy tools.
> {: .agenda}

# Preparing Your Data

First, we need to upload the data we'll work with. Ensure that you have the necessary files prepared from Zenodo or a shared data library within Galaxy.

## Data Upload

> <hands-on-title>Data Upload</hands-on-title>
>
> 1. Create a new history for this tutorial in Galaxy.
> 2. Import the mitoflash imaging data from [Zenodo]({{https://zenodo.org/records/14051280/files/Mitoflash_CM-8bit.tif}}) or from the shared data library:
>
>    ```
>    File1: Mitoflash_CM-8bit.tif 
>    ```
>
>    **Tip**: If any extra files are included, remove them to keep the workspace clean.
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets to keep track of them easily, e.g., "MitoFlash."
> 4. Confirm that the datatypes are correct for each file:
>    - `Mitoflash_CM-8bit.tif` should be an image file format.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
> 5. Tag each dataset with a label such as "mitoflash" for easy identification later if you are working multiple images.
>
{: .hands_on}

# Analyzing Mitochondrial Movement

In this section, we will focus on identifying mitochondrial regions in a time-lapse sequence by detecting spots based on local intensity maxima. This step provides insights into mitochondrial energy production and cellular health dynamics.

## Step 1: Detecting Mitochondrial Regions

> <hands-on-title>Detecting Mitochondrial Regions</hands-on-title>
>
> 1. {% tool [Spot Detection](toolshed.g2.bx.psu.edu/repos/imgteam/spot_detection/ip_spot_detection/1.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Image input"*: `Mitochondrial_Flashes.tiff`
>    - Set parameters such as "Thresholding" to optimize mitochondrial detection.
>
>    This tool detects mitochondria as spots or 'mitoflashes' based on intensity maxima, identifying them as distinct regions. Adjusting the threshold parameter can improve the accuracy of spot detection, allowing for better identification of mitoflash events.
>
{: .hands_on}

## Step 2: Tracking Mitochondrial Movement

> <hands-on-title>Tracking Mitochondrial Movement</hands-on-title>
>
> 1. {% tool [Perform linking in time series (nearest neighbors)](toolshed.g2.bx.psu.edu/view/imgteam/points_association_nn) %} with the following parameters:
>    - {% icon param-file %} *"Coordinates (and intensities) of input points"*: Output from the **Spot Detection** tool.
>    - Set "Neighborhood size" and other tracking parameters to optimize linking across.
>
>    **Comment**:  This tool links the detected mitoflash events (spots) across consecutive frames by finding the nearest neighboring spot in each frame, effectively tracking mitochondrial movement and mitoflash events over time. Observing these tracks provides insight into mitochondrial dynamics, helping to analyze how mitochondria move, change, and interact within the cell, which can reflect their functional states and roles in cellular processes.
>
{: .hands_on}

# Detecting Mitoflashes

Mitoflashes are identified based on sudden changes in fluorescence intensity in mitochondria, signifying superoxide bursts. We will use a curve-fitting tool to capture these fluctuations and measure their characteristics.

## Step 3: Curve Fitting for Mitoflash Detection

> <hands-on-title>Mitoflash Detection</hands-on-title>
>
> 1. {% tool [Perform curve fitting](toolshed.g2.bx.psu.edu/repos/imgteam/curve_fitting/ip_curve_fitting/0.0.3-2) %} with the following parameters:
>    - {% icon param-file %} *"File name of input data points (xlsx)"*: Select `spots_linked.xlsx`, the output file from the **Perform linking in time series (nearest neighbors)** tool.
>    - **Degree of the polynomial function**: Set to `2nd degree` to capture intensity peaks that characterize mitoflash events.
>    - **Penalty**: Choose *Least absolute deviations (LAD)* for robust fitting to intensity fluctuations.
>    - **Alpha**: Set a significance level, such as `0.01`, to generate assistive curves if needed.
>
>    This step applies polynomial curve fitting to model intensity fluctuations over time that represent mitoflash events. After running the tool, key parameters to analyze from the output file include:
>    - **Amplitude (F/F0)**: The peak value in the **CURVE** or **CURVE_A** column, representing the magnitude of the superoxide burst or intensity change.
>    - **Tpk (Time to Peak)**: The frame number at which the **CURVE** or **CURVE_A** reaches its maximum value, indicating the speed of the burst.
>    - **T50 (Duration)**: The number of frames where the **CURVE** or **CURVE_A** values are at least half of the peak value, representing the duration of the mitoflash.
>
>    These derived parameters provide insights into mitochondrial function, particularly in response to cellular stress or dynamic changes in activity levels.
>
{: .hands_on}


# Analysis and Interpretation

Once mitoflashes have been detected, we can analyze their frequency, duration, and intensity. This information is crucial for understanding the physiological significance of these events in processes like muscle contraction, neuron development, and wound healing.

> <details-title>More details about mitoflash dynamics</details-title>
>
> Mitoflashes represent brief but significant superoxide burst events. They provide insight into mitochondrial responses to energy demands and stress, acting as indicators of mitochondrial health. The kinetics of mitoflashes, including frequency and parameters like amplitude, Tpk, and T50, are influenced by cellular conditions, such as bioenergetics and starvation.
>
{: .details}

# Conclusion

In this tutorial, we have covered key techniques for tracking mitochondria and detecting mitoflashes. By processing and analyzing live-cell imaging data, you can gain valuable insights into mitochondrial behavior and health in cells. Tracking mitoflashes can contribute to studies in metabolism, aging, and diseases linked to mitochondrial dysfunction.
