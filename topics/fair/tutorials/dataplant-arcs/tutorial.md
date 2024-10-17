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

In this tutorial we will guide you through the process of creating your ARC (Annotated Research Context)


TODO: bit of background/intro about DataPLANT, ARCs etc

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


## Get Set up for the tutorial

For this tutorial you need three things:
 1. An account on DataHUB
 2. The ARCitect tool installed
 3. Example data to fill your ARC with

Below we will walk you through each of these 3 steps

### Create a DataHUB account

DataHUB is a GitLab instance designed to hosts ARCs. Here you can work collaboratively to create your ARC, and once you are ready to publish your ARC, you can do this from here as well.

TODO: details box discussing multiple DataHUBs?

> <hands-on-title> Create a DataHUB account </hands-on-title>
>
> 1. Already have an account? Please [log in](https://git.nfdi4plants.org/])
> 2. New to DataHUB? Create an account: [https://git.nfdi4plants.org/](https://git.nfdi4plants.org/)
>    - If you are not familiar with Git or GitLab that is ok, ARCitect will take care of syncing to DataHub behind the scenes
>
{: .hands_on}


### Download demo data

In this tutorial we will create an ARC for an example investigation. In this example investigation we grew Talinum plants under
drought conditions, with watered plants as control. Under these conditions the plants switch their type of photosynthesis
(water: C3, drought: CAM). To study this change, we profiled gene expression (Transciptomics) using RNASeq, and metabolites
(Metabolomics) via GC-MS (Gas-Chromatography Mass Spectrometry)


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
>    TODO: screenshot of data?
>
{: .hands_on}


### Install ARCitect

Now we will install the ARCitect tool. This tool will help you create and fill your ARC.


> <hands-on-title> Install the ARCitect tool on Linux  </hands-on-title>
>
> 1. Follow the installation instructions on the [DataPLANT knowledgebase](https://nfdi4plants.org/nfdi4plants.knowledgebase/) for your OS:
>    - [Windows Installation Instructions](https://nfdi4plants.org/nfdi4plants.knowledgebase/docs/ARCitect-Manual/arcitect_installation_windows.html)
>    - [macOS Installation Instructions](https://nfdi4plants.org/nfdi4plants.knowledgebase/docs/ARCitect-Manual/arcitect_installation_macos.html)
>    - [Linux Installation Instructions](https://nfdi4plants.org/nfdi4plants.knowledgebase/docs/ARCitect-Manual/arcitect_installation_linux.html)
>
>    > <comment-title> ARCitect version </comment-title>
>    >
>    > This tutorial uses version {{page.extra.arcitect_version}} of ARCitect. Since this tool is still
>    > under rapid development, a newer version of ARCitect may already be available. Feel free to install
>    > the latest version of ARCitect, but be aware that some screenshots and instructions may be outdated.
>    >
>    {: .comment}
>
>
>
> 2. Start ARCitect to verify that everything is installed correctly. You should see a screen like this:
>
>    ![screenshot of the ARCitect interface after startup](images/arcitect-home.png)
>
{: .hands_on}


## ARCitect: Initialize your ARC structure

Now that you have everything set up for the course, we can start creating our ARC. First we will organize our example data into the ARC structure, then we will add structured metadata about our experiment.


### ISA model

ARCs build on the [ISA Abstract Model](https://isa-specs.readthedocs.io/en/latest/isamodel.html)
for metadata annotation.


The ISA model comes with a hierarchy (ISA: Investigation - Study - Assay)
that aligns well with most projects in (plant) biology labs. It allows grouping multiple assays (measurements) to one study,
and multiple studies to one investigation.

![Overview of the ISA model](images/isa-model.png "Image source left panel: [isa-tools.org](https://isa-tools.org/format/specification.html)")

Your ARC has one isa.investigation.xlsx workbook at its root, where metadata about the investigataion is captured. Similarly, each study or assay that you add to your ARC contains one isa.study.xlsx or isa.assay.xlsx, respectively.


> <comment-title> Before you start </comment-title>
>
> Before creating your own ARC, invest some time to think about the following questions.
>
> - What is my **investigation**?
> - What is my **study**?
> - Which **assay** did I perform?
> - What is my (raw) **dataset**?
> - What **protocols** did I use?
>
{: .comment}


### The ARCitect interface

The ARCitect window consists of 4 main parts; on the left is the menu panel, next to it is a panel displaying the ARC folder structure, then the main panel where we will configure our ARC, and on the righ is (optionally) a help panel.

![](images/ARCitect-GUI.png "Overview of ARCitect interface. A) menu panel, B) ARC folder structure, C) Main panel, D) Help panel (off by default).")

Let's have a closer look at the menu panel:

![screenshot of the ARCitect menu panel](images/arcitect-sidebar.png)

This menu bar allows you to:
1. **Login**: Log in to DataHUB (required to sync your ARC to DataHUB)
2. **New ARC**: Create a new ARC
3. **Open ARC**: Open an existing ARC from your machine
4. **Download ARC**: From DataHUB to your machine (e.g. from a colleague or published)
5. **Save ARC**: Save your changes to your machine
6. **Explorer**: Open the local ARC folder on your computer
7. **Commit**: prepare your ARC changes for sync to DataHUB
8. **DataHUB sync**: push your comitted changes to DataHUB
9. **History**: Show previous versions of your ARC
10. **Validation**: Configure validation (checks) for our ARC
11. **Services**: Check on the status of ARC services (e.g. is DataHUB down for maintenance?)
12. **Settings**: For example turn on tool tips or the help panel



### Creating a new ARC

Let's start by creating a new ARC.

> <hands-on-title> Create a new ARC in ARCitect </hands-on-title>
>
> 1. In the left-hand panel of ARCitect, select **New ARC**
>    - Navigate to a location on your machine where you would like to create the ARC folder
>    - Choose a descriptive name for your ARC, e.g. `TalinumPhotosynthesis`
>
>    > <comment-title>What's in a name?</comment-title>
>    >
>    > By default, the name give your ARC here will be used for:
>    >  - the ARC folder on your machine
>    >  - the ARC repository on DataHUB at `https://git.nfdi4plants.org/<YourUserName>/<YourARC>` (see next steps)
>    >  - the identifier for your investigation
>    >
>    > **Note:** Make sure that no ARC exists at `https://git.nfdi4plants.org/<YourUserName>/<YourARC>`. Otherwise you will sync to that ARC.
>    {: .comment}
>
> 2. If all went well, you should see the folder structure that was created for your ARC:
>    ![](images/new-arc.png)
{: .hands_on}

### Add investigation-level metadata

We will start by adding some basic data about our ARC. We do this in the investigation-level metadata.

> <hands-on-title>Adding investigation information</hands-on-title>
>
> 1. In the ARC folder panel (second panel from the left), click on the investigation
>    - This is the root of the ARC folder structure (e.g. *TalinumPhotosynthesis*)
>    - This will bring up the investigation information in the main panel
>    - Notice that the Identifier has already been filled in for you (but you can still change it)
>
>    ![screenshot of the investigation metadata screen](images/investigation-empty.png)
>
> 2. **Add a title** for your investigation
>    - This should be a concise but descriptive title
>    - Tip: think of what you would name your publication based on this investigation
>    - For this tutorial, you could your title to `Effect of Drought on Talinum Photosynthesis`
>
> 3. **Add a description** for your investigation
>   - Here you can describe your main research questions and methods in a few sentences or paragraphs.
>   - Tip: think of this as a short abstract like you would provide in a publication
>   - For this tutorial you can add any text here
>
> 4. **Add an ARC Contact** (you)
>    - Click on the **+ (plus) icon** under Contacts
>    - Click on the **⌄ (dropdown) icon** on the right of the contact to expand it
>    - Fill in at least your First name, Last name, and email.
>    - **Tip:** if you have an [ORCiD](https://orcid.org) ID, fill that in and hit "Search" to autofill these fields
>
> 5. Your investigation tab should look something like this
>
>    ![screenshot of investigation metadata filled in](images/investigation-filled.png)
>
{: .hands_on}


### Saving your ARC

It is a good idea to frequently save your work. There are 2 levels of saving your ARC; you can:
1. Save it locally (to your computer)
2. Save it online to DataHUB

By default your ARC on DataHUB is private to you, but you have the option of sharing it with your collaborators.

> <comment-title>Why and how often to sync to DataHUB?</comment-title>
> Syncing your ARC to DataHUB frequently has several advantages:
>  - It makes it easy to collaborate with colleagues on your ARC
>  - It is a backup of your ARC in case anything happens to your computer
>  - It is version controlled. This means you can always go back to any previous version your ARC if needed
>
> We advise you to save your ARC to your local machine as often as you can, and sync to DataHUB every time you have
> completed a unit of work (e.g. "added data files", "added study metadata", "added protocol document" etc), this
> makes it easy for others (and you) to see what you have changed.
>
{: .comment}

Since we have just completed our first unit of work (initialized our ARC and added basic investigation metadata), let's save and sync our ARC to DataHUB now.

> <hands-on-title>Save and Sync our ARC</hands-on-title>
>
> 1. Save locally simply by hitting the **Save ARC** button on the left
>    - That's it! Too easy!
>    - Do this as often as possible to avoid losing work
>
> 2. Click **Login** at the top of the menu panel
>    - Select `git.nfdi4plant.org` as your host
>    - Click the **Login** button
>    - Since you already logged into DataHub via your browser earlier, ARCitect should now automatically log you in
>
> 3. To check that you are successfully logged in, the Login button should now be replace with your DataHub username
>    ![screenshot of the menu panel as it looks when you are logged in](images/arcitect-loggedin.png)
>
> 4. Click on **Commit** in the menu
>    - This screen will show you all files that have chanced since your last sync
>    - You can choose whether to sync all files, or only some
>    - Choosing which files to sync is called "committing" those file.
>    - In our case only one file has changed, `isa.investigation.xlsx`. By default it is already selected to be committed
>
> 5. Add a **commit message**
>    - in a few words, describe your changes
>    - this will help you and your colleagues
>    - For example, for this changes you can write `added investigation metadata`
>    - You can leave the rest of the fields to the default values
>
>    ![screenshot of the commit menu](images/arcitect-commit.png)
>
> 6. Click the **Commit** button
>    - There will be a confirmation screen like the one below
>    - Click the **Ok** button
>
>    ![](images/arcitect-commit-confirmation.png)
>
> 7. Click on **DataHUB Sync** in the menu
>
>    ![Screenshot of the DataHUB sync menu](images/arcitect-sync.png)
>
> 8. Click on **+ (plus)** to add a 'remote' (DataHUB location to save your ARC)
>    - you only have to do this once
>    - you can keep all the default settings
>    - Click **Add**
>
>    ![](images/arcitect-add-remote.png)
>
> 9. Click **Push** to sync from your computer to DataHUB
>    - **Pull** will do the reverse, it will 'pull' in changes from DataHUB into your local copy of the ARC
>
> 10. Check that the push was successfull
>     - Navigate to `https://git.nfdi4plants.org/YourUsername/TalinumPhotosynthesis`
>     - Replace 'YourUsername' with your DataHUB username in the URL
>     - If you used a different Identifier for your ARC, replace 'TalinumPhotosynthesis' with your ARC Identifier in the URL
>     - You should see your ARC on DataHUB:
>
>     ![screenshot of your new ARC on DataHUB](images/datahub-firstpush.png)
>
{: .hands_on}

Congrats! You have initialized, saved, and synced your ARC! Now let's continue with the good part, filling your ARC.

### Create your assays and studies




## ARCitect: Adding files to your ARC

add example files to various folders (protocols, assays datasets, workflows, etc?)


## ARCitect: Add your experimental metadata (SWATE)


## ARCitect: Syncing to DataHUB

### ARC validation on DataHUB

## Publishing your ARC

### Enabling Invenio testing

how to enable the invenio testing, how to troubleshoot errors


### Publishing Process

Once all of the test pass, you can start the submission process for publishing your ARC



## Further Reading

Link to knowledgebase

Link to quickstart videos https://nfdi4plants.org/nfdi4plants.knowledgebase/docs/guides/arcitect_QuickStart_Videos.html

Link to follow up tutorials
 - Tutorial for CWL/Workflows part of the ARC?
 - DataHUb tutorial (sharing, groups, etc?)
 - Creating of metadata templates?






