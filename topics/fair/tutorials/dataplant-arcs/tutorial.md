---
layout: tutorial_hands_on
title: DataPLANT ARCs
abbreviations:
  FAIR: Findable, Accessible, Interoperable, Reusable
  RDM: Research Data Management
  DMP: Data Management Plan
  PID: Persistent Identifier

zenodo_link: 'https://zenodo.org/records/13935061'
questions:
- What are ARCs?
- Why should I create ARCs?
- How can I create my ARC?
time_estimation: "2H"
key_points:
- Annotated Research Contexts (ARCs) ..
- ARCitect is a useful tool to help you create your ARC
tags:
- plants
- fair
- data stewardship
priority: 2
contributions:
  authorship:
    - Brilator
    - CMR248
    - Freymaurer
    - Martin-Kuhl
    - SabrinaZander
    - StellaEggels
  editing:
    - shiltemann
  funding:
    - nfdi4plants

subtopic: dataplant

requirements:
  - type: "internal"
    topic_name: fair
    tutorials:
      - fair-intro

extra:
  arcitect_version: "0.48"
---

intro

> <comment-title> ARCitect version </comment-title>
>
> This tutorial uses version {{page.extra.arcitect_version}} of ARCitect. Since this tool is still
> under rapid development, a newer version of ARCitect may already be available and screenshots may be outdated (despite our best
> efforts to keep this tutorial updated)
>
{: .comment}

## Background

background about DataPLANT, RO-crates, FAIR, RDM, etc



## Get Set up for the tutorial

For this tutorial you need three things:
 1. An account on DataHUB
 2. The ARCitect tool installed
 3. Example data to fill your ARC with

Below we will walk you through each of these 3 steps

### Create a DataHUB account

DataHUB is a GitLab instance designed to hosts ARCs. Here you can work collaboratively to create your ARC, and once you are ready to publish your ARC, you can do this from here as well.

TODO: discuss multiple DataHUBs?

> <hands-on-title> Create a DataHUB account </hands-on-title>
>
> 1. Already have an account? Please [log in](https://git.nfdi4plants.org/])
> 2. New to DataHUB? Create an account: [https://git.nfdi4plants.org/](https://git.nfdi4plants.org/)
> 3. If you are not familiar with Git or GitLab that is ok, ARCitect will take care of syncing to DataHub
>    - TODO: link to section/slides with more info about DataHUB
>
>
{: .hands_on}


### Download demo data


TODO: a few sentences about the demo data (what kind of data, where did it come from, what was the study?)


> <hands-on-title> Download example ARC data </hands-on-title>
>
> 1. Download the zip file from Zenodo (<2MB) to your machine
>
>    ```
>    https://zenodo.org/records/13935061/files/arcitect-demo-data.zip
>    ```
>
> 2. Unzip the file
>    - You can put the files wherever is most convenient for you
>
> 3. Have a quick look at the files and folders
>    - This example data mimics you might have your project data organised on your computer
>    - In the next part we will show how to organise this data in an ARC-compatible way
>
{: .hands_on}


### Install ARCitect

Now we will install the ARCitect tool.

TODO: make this a choose-your-own adventure section, instructions per OS


> <hands-on-title> Install the ARCitect tool on Linux  </hands-on-title>
>
> 1. Follow the installation instructions on the [DataPLANT knowledgebase](https://nfdi4plants.org/nfdi4plants.knowledgebase/) for your OS:
>    - [Windows Installation Instructions](https://nfdi4plants.org/nfdi4plants.knowledgebase/docs/ARCitect-Manual/arcitect_installation_windows.html)
>    - [macOS Installation Instructions](https://nfdi4plants.org/nfdi4plants.knowledgebase/docs/ARCitect-Manual/arcitect_installation_macos.html)
>    - [Linux Installation Instructions](https://nfdi4plants.org/nfdi4plants.knowledgebase/docs/ARCitect-Manual/arcitect_installation_linux.html)
>
> 2. Start ARCitect to verify that everything is installed correctly. You should see a screen like this:
>
>    ![screenshot of the ARCitect interface after startup](images/arcitect-home.png)
{: .hands_on}


## ARCitect: Structure your data

Now that you have everything set up for the course, we can start creating our ARC.

First we will organize our example data into the ARC structure. ARCs build on the [ISA Abstract Model](https://isa-specs.readthedocs.io/en/latest/isamodel.html)
for metadata annotation. The ISA model comes with a hierarchy (ISA: Investigation - Study - Assay)
that aligns well with most projects in (plant) biology labs. It allows grouping multiple assays (measurements) to one study,
and multiple studies to one investigation.

![Overview of the ISA model](images/isa-model.png "Image source left panel: [isa-tools.org](https://isa-tools.org/format/specification.html)")

Your ARC has one isa.investigation.xlsx workbook at its root, where metadata about the investigataion is captured. Similarly, each study or assay that you add to your ARC contains one isa.study.xlsx or isa.assay.xlsx, respectively.


## ARCitect: Add your metadata (SWATE)


## ARCitect: Syncing to DataHUB



## Publishing your ARC


