---
layout: tutorial_hands_on

title: Creating FAIR Quality assessment reports and draft of Data Papers from EML
  metadata with MetaShRIMPS
zenodo_link: https://doi.org/10.5281/zenodo.8130567
questions:
- How to improve the FAIR quality of an EML metadata ?
- How to use metadatas for machine actionnable processes ?
- What is the point of having a FAIR metadata ?
objectives:
- Learn how to use the interactive tool Metashrimps
- Understand the challenges MetaShRIMPS is trying to respond to
- Learn how to create a FAIR Quality assessment report from a metadata using EML standard
- Understand the concept of Data Paper and learn how to produce it
- Explain the necessity of using such tools when producing ecological metadata
time_estimation: 30mn
key_points:
- The take-home messages
- They will appear at the end of the tutorial
tags:
  - Metadata
  - EML
  - FAIR
  - Data Paper
contributors:
- TanguyGen
- yvanlebras

---


# Introduction

This tutorial aims to teach how to use the interactive tool MetaShRIMPS, available on Galaxy Ecology,
to produce Data Papers drafts and FAIR quality assessment reports from metadata using EML
standard.
This tutorial purpose is also to explain the purpose of this tool which is to improve the overall FAIR quality
of metadata.
MetaShRIMPS is a tool developped by the PNDB (French Biodiversity Data Hub), a research structure initiated by the 
French Natural History Museum (MNHN).

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Why do we need FAIR Data ?
<img src="./Images/FAIR_data_principles.jpg" alt="FAIR Data Principles" width="500"/>
# How can this tool help in improving metadata quality
Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

## Get data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

# Title of the section usually corresponding to a big step in the analysis

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> <details-title> More details about the theory </details-title>
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:



## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
