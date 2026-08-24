---
layout: tutorial_hands_on

title: "Use Pycytominer in Galaxy for processing high dimensional readouts from high-throughput image-based profiling experiments"
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

High-content imaging screens generate thousands of microscopy images capturing how cells respond to genetic or chemical perturbations. In particular, [cell painting assays](https://en.wikipedia.org/wiki/Cell_painting) are a high-content/high-throughput imaging methods designed to reveal a broad range of cellular phenotypes. Image analysis makes it possible to extract numerical information from cell shape, producing what are known as "morphological profiles". These morphological profiles offer a window into a variety of biological processes, such as how cells react to genetic modifications, drug exposure, and shifts in their environment ({% cite seal2025cell %}). Extracting biological meaning from these images requires transforming raw features and readouts into clean, comparable profiles. Because of the sheer number of values involved, these results can be extremely large, and specific frameworks are needed to process such morphological profiles correctly.

In this context, **Pycitominer** is a Python toolkit for processing high dimensional readouts from high-throughput image-based profiling experiments ({% cite serrano2025reproducible %}).

![pycitominer-logo.png](../../images/pycitominer/pycitominer-logo.png){: width="50%"}

In this tutorial, you will learn how to run a Pycitominer pipeline using Galaxy. We will follow the different steps explained in the [Pycitominer documentation](https://pycytominer.readthedocs.io/en/stable/tutorials/introduction_to_pycytominer.html). If you want a more comphresince explanation of
each step, please fell free to visit the main [Pycitominer main documentation page](https://pycytominer.readthedocs.io/en/stable/tutorials/introduction_to_pycytominer.html) or the [GitHub repository](https://github.com/cytomining/pycytominer)!

> <agenda-title></agenda-title>
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Getting data

The data necessary for this tutorial can be created following the instruction of the Pycitominer documentation. However, for simplicity, we already made them available for you here in the training!

> <hands-on-title>Data Upload</hands-on-title>
>
> 1. Create a new history for this tutorial.
>
> 2. Download the following image and import it into your Galaxy history.
>    - [`01_platemap.tsv`](workflows/test-data/01_platemap.tsv)
>    - [`01_single_cells.tsv`](workflows/test-data/01_single_cells.tsv)
>    
>    If you are importing the image via URL:
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    If you are importing the image from the shared data library:
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Confirm the datatypes are correct (`tabular` for both images)
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
> 
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
{: .hands_on}

# Step 1: Aggregate — From Cells to Wells

**Aggregation** collapses single-cell measurements into a single profile per well or per sample by computing a summary statistic (such as the median) across all cells.

# Step 2: Annotate — Adding Experimental Context

**Annotation** merges these profiles with experimental metadata (i.e. plate and well identifiers, and other conditions) so each profile is linked to what was done to the cells.

# Step3: Normalize — Removing Technical Variation

**Normalization** rescales features to make them comparable across plates and batches, commonly by standardizing each feature against control samples to correct for plate-to-plate variation.

# Step 4: Feature Selection — Keeping Only Informative Features

**Feature selection** removes uninformative or redundant features, such as those with low variance, high correlation with other features, or missing values, yielding a compact and reliable feature set.

# Step 5: Consensus — Collapsing Replicates

**Compute Consensus** replicate profiles into one consensus profile per treatment group by computing the median across all replicates.

# A full workflow for data treatment

# Conclusion



