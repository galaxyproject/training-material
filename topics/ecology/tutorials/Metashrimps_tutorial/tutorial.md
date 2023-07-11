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

> <comment-title>Why do we need FAIR Data and Metadata?</comment-title>
> [FAIR](https://www.go-fair.org/fair-principles/) stand for **Findable, Accessible, Interoperable, Reusable**. 
>
>
><img src="./Images/FAIR_data_principles.jpg" alt="FAIR Data Principles" width="500"/>
>This principles were [officialy instaured in 2016](https://www.go-fair.org/fair-principles/) to improve and the access and usabiliy of data by the machine and to help making data reusable and shareable for users.
>Metadata is the data used to describe and explain all the context behind the production of data. It is necessary to produce a rich and FAIR metadata in order 
>to permit external users to understand and reuse data for their own studies.
>The purpose of this tool is to help the user improve their metadata quality in order to increase its value to the scientific community and to help highlighting
>the work of all the producers of the data.
{:  .comment}
> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Why do we need FAIR Data and Metadata?


## How can this tool improve the quality metadata?

This tool aims to give an easy access to a quality assessment report of metadatas which could guide the producers of
data/metadata to the production of a FAIR metadata.


# How to use MetaShRIMPS?
## Get data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2.Import this file from [Zenodo]({{ https://doi.org/10.5281/zenodo.8130567 }}) to test it
>     -> `{{ Training Data for "Creating Quality FAIR assessment reports and draft of Data Papers from EML metadata with MetaShRIMPS" }}`):
>    ```
>    https://zenodo.org/record/8130567/files/Kakila_database_marine_mammal.xml
>    ```
>
>    {% snippet topics/ecology/tutorials/Metashrimps_tutorial/import_files.md %}

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
