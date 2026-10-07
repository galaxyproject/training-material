---
layout: tutorial_hands_on

title: "Use Pycytominer in Galaxy for processing high dimensional image-based readouts"
level: Intermediate
subtopic: analyses
questions:
  - "How do I process high-dimensional readouts using Galaxy?"
  - "How can I process CellProfiler or DeepProfiler feature readouts with Pycytominer in Galaxy?"
objectives:
  - "Learn how to execute the five main Pycytominer steps"
  - "Create a complete workflow for processing data with Pycytominer"
key_points:
- CellProfiler readouts can be easily processed in Galaxy using the Pycytominer CLI
- Reproducible workflows can be created to process readouts consistently

requirements:
  -
    type: "internal"
    topic_name: imaging
    tutorials:
      - imaging-introduction
  -
    type: "internal"
    topic_name: imaging
    tutorials:
      - tutorial-CP

time_estimation: "1H"
contributions:
  authorship:
    - rmassei
  funding:
    - nfdi4bioimage
tags:
  - Object feature extraction
  - high-throughput
  - data cleaning
  - cellprofiler
  - deepprofiler

---

High-content imaging screens generate thousands of microscopy images capturing how cells respond to genetic or chemical perturbations. In particular, [Cell Painting assays](https://en.wikipedia.org/wiki/Cell_painting) are a high-content/high-throughput imaging methods designed to reveal a broad range of cellular phenotypes. Image analysis makes it possible to extract numerical information from cell shape, intensity, texture, granularity, producing what are known as "morphological profiles". These features combined form morphological profiles offering a window into a variety of biological processes, such as how cells react to genetic modifications, drug exposure, and shifts in their environment ({% cite seal2025cell %}). Extracting biological meaning from these images requires transforming raw features and readouts into clean, comparable profiles. Because of the sheer number of values involved, these results can be extremely large, and specific frameworks are needed to process such morphological profiles correctly.

In this context, **Pycytominer** is a Python toolkit for processing high dimensional readouts from high-throughput image-based profiling experiments ({% cite serrano2025reproducible %}).

![pycytominer-logo.png](../../images/pycytominer/pycytominer-logo.png){: width="50%"}

In this tutorial, you will learn how to run a Pycytominer pipeline using Galaxy. We will follow the different steps explained in the [Pycytominer documentation](https://pycytominer.readthedocs.io/en/stable/tutorials/introduction_to_pycytominer.html). If you want a more comprehensive explanation of each step, please feel free to visit the main [Pycytominer main documentation page](https://pycytominer.readthedocs.io/en/stable/tutorials/introduction_to_pycytominer.html) or the [GitHub repository](https://github.com/cytomining/pycytominer)!
each step, please feel free to visit the main [Pycytominer main documentation page](https://pycytominer.readthedocs.io/en/stable/tutorials/introduction_to_pycytominer.html) or the [GitHub repository](https://github.com/cytomining/pycytominer)!

> <agenda-title></agenda-title>
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Getting data

A synthetic dataset necessary for this tutorial can be created following the instructions of the [Pycytominer documentation](https://pycytominer.readthedocs.io/en/stable/tutorials/introduction_to_pycytominer.html#Tutorial-Data).

**Experimental design**:

| **Property** | **Value** |
|:-----|:------:|
| Plates (biological replicates)    | 1     |
| Wells per plate   | 6 (2 × DMSO vehicle control, 2 × Compound A, 2 × Compound B)      |
| Cells per well    | ~100      |
| Total single-cell measurements  | ~600    |
| Morphological features   | 11 (across three compartments)      |
|:-----|:------:|


For simplicity, we provide the generated files for this tutorial.

> <hands-on-title>Data Upload</hands-on-title>
>
> 1. If you are logged in, create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Download the following image-based profiles and import them into your Galaxy history.
>    - [`01_platemap.tsv`](workflows/test-data/01_platemap.tsv)
>    - [`01_single_cells.tsv`](workflows/test-data/01_single_cells.tsv)
>    
>    If you are importing the files via URL:
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    If you are importing the files from the shared data library:
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Confirm the datatypes are correct (`tabular` for both profiles)
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
> 
{: .hands_on}

## Step 1: Aggregate — From Cells to Wells

**Aggregation** collapses single-cell measurements into a single profile per well or per sample by computing a summary statistic (such as the median) across all cells.

> <hands-on-title>Aggregate plate readouts with Pycytominer</hands-on-title>
>
> 1. {% tool [Aggregate readouts](toolshed.g2.bx.psu.edu/repos/imgteam/pycytominer_aggregate/pycytominer_aggregate/1.6.1+galaxy0) %} with the following parameters to aggregate readouts:
>    - {% icon param-file %} *"Input feature-readouts table"*: `01_single_cells.tsv` file
>    - *"Aggregation Column"*: Select "c1:Metadata_Plate" and "c2:Metadata_Well"
>    - *"Aggregation function"*: `Mean`
>
> 2. Rename {% icon galaxy-pencil %} the generated file to `01_output_aggregate.tsv`.
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
> 3. Click on the **visualise icon** {% icon galaxy-visualise %} of the file to visually inspect the image-based profiles using the **Tabulator** visualisation plugin.
{: .hands_on}

The 600 single-cell measurements are now aggregated into 6 profiles, one per well: for each feature, the values of the ~100 cells in a well are summarized into a single value (here, the mean).

![01-aggregate.png](../../images/pycytominer/01-aggregate.png)

## Step 2: Annotate — Adding Experimental Context

**Annotation** merges these profiles with experimental metadata (i.e. plate and well identifiers, and other conditions) so each profile is linked to what was done to the cells.

> <hands-on-title> Annotate readouts with metadata using Pycytominer</hands-on-title>
>
> 1. {% tool [Annotate readouts with metadata](toolshed.g2.bx.psu.edu/repos/imgteam/pycytominer_annotate/pycytominer_annotate/1.6.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input feature-readouts table"*: `01_output_aggregate.tsv` file
>    - *"Column describing the wells in the feature-readouts table"*: Select "c2:Metadata_Well"
>    - {% icon param-file %} *"Input platemap table"*: `01_platemap.tsv` file
>    - *"Column describing the wells in the platemap"*: Select "c1:well_position"
> 2. Rename {% icon galaxy-pencil %} the generated file to `02_output_annotated.tsv`.
> 3. Click on the **visualise icon** {% icon galaxy-visualise %} of the file to visually inspect the image-based profiles using the **Tabulator** visualisation plugin.
{: .hands_on}

Three additional columns are now added to the table: "Metadata_treatment", "Metadata_cell_line" and "Metadata_concentration_um". Pycytominer adds the Metadata_ prefix to the plate map columns (treatment → Metadata_treatment) to distinguish them from the morphological features. All this information is important to give more context to the data. Metadata_treatment allows us to identify the DMSO control wells used for normalization.

![02-annotate.png](../../images/pycytominer/02-annotate.png)

## Step 3: Normalize by Removing Technical Variation

**Normalization** rescales features so they can be compared with each other. Without it, features with large absolute values (e.g. cell area) would dominate any downstream distance calculation, regardless of whether they carry biological signal. A common approach is to standardize each feature against control samples: here, every well is expressed relative to the DMSO control wells. In experiments with several plates or batches, normalizing each plate against its own controls also corrects technical variation, such as differences in staining, imaging conditions or cell density.

> <hands-on-title>Normalize readouts with Pycytominer</hands-on-title>
>
> 1. {% tool [Normalize readouts](toolshed.g2.bx.psu.edu/repos/imgteam/pycytominer_normalize/pycytominer_normalize/1.6.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input feature-readouts table"*: `02_output_annotated.tsv` file
>    - *"Column with normalization values"*: Select "c1:Metadata_treatment"
>    - *"Value"*: Type "DMSO"
>    - *"Normalization method"*: Select "Standardize"
> 2. Rename {% icon galaxy-pencil %} the generated file to `03_output_normalized.tsv`.
> 3. Click on the **visualise icon** {% icon galaxy-visualise %} of the file to visually inspect the image-based profiles using the **Tabulator** visualization plugin.
{: .hands_on}

With the standardize normalization method, feature becomes a z-score based on the mean and standard deviation of the DMSO wells. So DMSO wells end up around 0, and treated wells show how many standard deviations they differ from the control. All features are now expressed in the same unit (standard deviations from the DMSO control), so they can be compared with each other.

![03-normalize.png](../../images/pycytominer/03-normalize.png)

## Step 4: Feature Selection — Keeping Only Informative Features

**Feature selection** removes uninformative or redundant features, such as those with low variance, high correlation with other features, or missing values, yielding a compact and reliable feature set.

> <hands-on-title>Select informative features with Pycytominer</hands-on-title>
>
> 1. {% tool [Select informative features](toolshed.g2.bx.psu.edu/repos/imgteam/pycytominer_feature_select/pycytominer_feature_select/1.6.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input feature-readouts table"*: '03_output_normalized.tsv` file
>    - *"Operations"*: Select "Variance Threshold" and "Blocklist"
> 2. Rename {% icon galaxy-pencil %} the generated file to `04_output_features.tsv`.
> 3. Click on the **visualise icon** {% icon galaxy-visualise %} of the file to visually inspect the image-based profiles using the **Tabulator** visualization plugin.
{: .hands_on}

 Variance Threshold removes features that barely vary across samples, and Blocklist removes features that are known to be noisy or uninformative in image-based profiling. Thanks to the Variance Threshold operation, the table now has 15 columns instead of 16: Cells_AreaShape_EulerNumber was removed because it has the same value in every well, so it carries no information. 

![04-features.png](../../images/pycitominer/04-features.png)

## Step 5: Consensus — Collapsing Replicates

The **Compute Consensus** tool collapses replicate profiles into one consensus profile per treatment group by computing the mean across all replicates.

> <hands-on-title>Compute consensus profiles with Pycytominer</hands-on-title>
>
> 1. {% tool [Compute consensus](toolshed.g2.bx.psu.edu/repos/imgteam/pycytominer_consensus/pycytominer_consensus/1.6.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input feature-readouts table"*: 04_output_features.tsv` file
>    - *"Column with unique condition"*: Select "c1:Metadata_treatment", "c2:Metadata_cell_line" and "c3:Metadata_concentration_um"
>    - *"Reduction operation"*: Select "Mean"
> 2. Rename {% icon galaxy-pencil %} the generated file to `05_output_consensus.tsv`.
> 3. Click on the **visualize icon** {% icon galaxy-visualise %} of the file to visually inspect the image-based profiles using the **Tabulator** visualization plugin.
{: .hands_on}

The table now has 3 profiles, one per treatment (DMSO, Compound A and Compound B).

![05-consensus.png](../../images/pycytominer/05-consensus.png)

## A full workflow for table readouts processing

You can now create a workflow from the different Pycytominer steps in your history:

> <hands-on-title> Extract Pycytominer workflow from history  </hands-on-title>
> 1. Now we can extract the workflow for batch processing:
>
>    {% snippet faqs/galaxy/workflows_extract_from_history.md %}
>
>    - Name it "pycytominer-full-steps".
>    - Uncheck `01_platemap.tsv` and `01_single_cells.tsv` as inputs (the workflow is supposed to be applied to the image-based profiles directly).
> 2. Edit the workflow you just created:
>    - Select "Input dataset" from the list of tools. The step {% icon param-file %} **8: Input Dataset** appears.
>    - Select "Input dataset" from the list of tools. The step {% icon param-file %} **9: Input Dataset** appears.
>    - Change the "Label" of {% icon param-file %} **8: Input Dataset** to `input table readouts`.
>    - Change the "Label" of {% icon param-file %} **9: Input Dataset** to `input plate metadata`.
>    - Connect the output of {% icon param-file %} **8: input table readouts** to the input of {% icon tool %} **3: Aggregate readouts**.
>    - Connect the output of {% icon param-file %} **9: input plate metadata** to the "Input plate Table" input of {% icon tool %} **4: Annotate readouts with metadata**.
>    - Mark the results of {% icon tool %} **7: Compute Consensus Profile** as the primary outputs of the workflow (by clicking on the checkboxes of the outputs).
>
{: .hands_on}

You have now a Pycytominer automatized workflow in Galaxy! 

![06-final-workflow.png](../../images/pycitominer/06-final-workflow.png)

# Conclusion

The following tutorial uses high-content imaging analysis as an example, but the same Pycytominer tools can be used
in many other contexts for data wrangling, normalization, and annotation... Find your own solution!