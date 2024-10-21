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
time_estimation: 2H
contributions:
  authorship:
    - shiltemann
    - Brilator
    - SabrinaZander
    - CMR248
    - Freymaurer
    - Martin-Kuhl
    - StellaEggels
  #editing:
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
 1. An account on [DataHUB](https://git.nfdi4plants.org)
 2. The [ARCitect](https://www.nfdi4plants.de/nfdi4plants.knowledgebase/docs/ARCitect-Manual/index.html) tool installed
 3. Example data to fill your ARC with

Below we will walk you through each of these 3 setup steps.

### Create a DataHUB account

DataHUB is a GitLab instance designed to hosts ARCs. Here you can work collaboratively to create your ARC, and once you are ready to publish your ARC, you can do this from here as well.

> <hands-on-title> Create a DataHUB account </hands-on-title>
>
> 1. Already have an account? Please [log in](https://git.nfdi4plants.org/])
> 2. New to DataHUB? Create an account: [https://git.nfdi4plants.org/](https://git.nfdi4plants.org/)
>    - If you are not familiar with Git or GitLab that is ok, ARCitect will take care of syncing to DataHub behind the scenes
>
{: .hands_on}


### Download demo data

In this tutorial we will create an ARC for an example research data. In this example investigation we grew Talinum plants under
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
>    - This example data mimics how you might have your project data organised on your computer
>    - In the next part we will show how to organise this data in an ARC-compatible way
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

Now that you have everything set up for the course, we can start creating our ARC. First we will organize our example data into the ARC structure, then we will add structured metadata describing our experiments and data.


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
4. {% icon download-cloud %} **Download ARC**: From DataHUB to your machine (e.g. from a colleague or published)
5. {% icon save %} **Save ARC**: Save your changes to your machine
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
>
> 3. Click on the **Explorer** button in the menu to see the files created on your machine
>
{: .hands_on}

The basic ARC structure is initialized for you, and looks something like this:

```bash
my-arc
├── studies
│   ├── study1              # subfolder for each study
│   │   ├── protocols       # folder for protocol documents
│   │   ├── resources       # folder for any other supporting information and files
│   │   isa.study.xlsx      # study-level metadata (e.g. sample sheet, culturing/harvesting conditions)
│   │   README.md           # Free-text description of your study
│   ├── assay2
│   │   ├── ...
├── assays
│   ├── assay1              # subfolder for each assay
│   │   ├── dataset         # assay result files (e.g. sequencing files)
│   │   ├── protocols       # assay protocol files
│   │   isa.assay.xlsx      # assay-level metadata (e.g. instrument)
│   │   README.md
│   ├── assay2
│   │   ├── ...
├── runs                    # analysis results
│   ├── ...
├── workflows               # analysis workflows (CWL, or any other scripts/code used)
│   ├── ...
├── isa.investigation.xlsx  # investigation-level metadata
└── README.md
```


### Add investigation-level metadata

We will start by adding some basic data about our ARC. We do this in the investigation-level metadata.

> <hands-on-title>Adding investigation information</hands-on-title>
>
> 1. In the ARC folder panel (second panel from the left), click on the investigation
>    - This is the root of the ARC folder structure (e.g. *TalinumPhotosynthesis*)
>    - This will bring up the investigation information in the main panel
>    - Notice that the Identifier has already been filled in for you (but you can still change it)
>    - Behind the scenes this information will be saved as an Excel sheet `isa.investigation.xlsx`
>
>    ![screenshot of the investigation metadata screen](images/investigation-empty.png)
>
> 2. **Add a title** for your investigation
>    - This should be a concise but descriptive title
>    - Tip: think of what you would name your publication based on this investigation
>    - For this tutorial, you could set your title to `Effect of Drought on Talinum Photosynthesis`
>
> 3. **Add a description** for your investigation
>   - Here you can describe your main research questions and methods in a few sentences or paragraphs.
>   - Tip: think of this as a short abstract like you would provide in a publication
>   - For this tutorial you can add any text here
>
> 4. **Add an ARC Contact** (you)
>    - Click on the **{% icon plus %} (plus) icon** under Contacts
>    - Click on the **{% icon dropdown %} (dropdown) icon** on the right of the contact to expand it
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
2. Save it online to DataHUB (syncing)

By default your ARC on DataHUB is private to you, but you have the option of sharing it with your collaborators.

> <comment-title>Why and how often to sync to DataHUB?</comment-title>
>
> Syncing your ARC to DataHUB frequently has several advantages:
>  - It makes it easy to collaborate with colleagues on your ARC
>  - It is a backup of your ARC in case anything happens to your computer
>  - It is version controlled. This means you can always go back to any previous version your ARC if needed
>
> We advise you to save your ARC to your local machine as often as you can, and sync to DataHUB every time you have
> completed a unit of work (e.g. "added data files", "added culturing metadata", "added protocol documents" etc), this
> makes it easy for others (and future you) to see what you have changed.
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
>    ![Screenshot of the commit confirmation dialogue window](images/arcitect-commit-confirmation.png)
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
>    ![screenshot of the remote adding dialogue window](images/arcitect-add-remote.png)
>
> 9. Click **Push** to sync from your computer to DataHUB
>    - **Pull** will do the reverse, it will 'pull' in changes from DataHUB into your local copy of the ARC,
>      for example if a collaborator made changes, or you worked on it from a different machine.
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

Congrats! You have initialized, saved, and synced your ARC! Now let's continue with the good part, actually filling your ARC.


## Structure your data in the ARC

### Think about your data

Next we will divide our research into studies and assays. As a rule of thumb, create a study for every research question, and an assay for every measurement you performed (e.g. sequencing, mass-spec, imaging). An investigation can have one study, or several. A study can have one or more assays, and each assay can be associated with one or more studies.


> <hands-on-title> Evaluate the example data </hands-on-title>
>
> 1. Recall our example research:
>
>    *In this example investigation we grew Talinum plants under drought conditions, with watered plants as control.
>    Under these conditions the plants switch their type of photosynthesis (water: C3, drought: CAM). To study this
>    change, we profiled gene expression (Transciptomics) using RNASeq, and metabolites (Metabolomics) via GC-MS
>    (Gas-Chromatography Mass Spectrometry)*
>
> 2. Have a look at the demo data we downloaded
>
> > <question-title> Which studies and assays do we have? </question-title>
> >
> > 1. Which studies do we have?
> > 2. Which assays did we perform?
> >
> > > <solution-title></solution-title>
> > > 1. We had one research question, the effect of drought on photosynthesis, so we will create 1 study
> > > 2. We performed 2 measurements, RNA sequencing and GC-MS, so we will create 2 assays in our ARC
> > {: .solution}
> {: .question}
{: .hands_on}

### Create a study

Let's create our study now.

> <hands-on-title>Create an ARC study</hands-on-title>
>
> 1. Click on the **{% icon plus %} (plus) icon** next to the studies folder
> 2. Provide a study identifier, e.g. `talinum_drought`
>    ![screenshot of the study creation dialogue window](images/arcitect-new-study.png)
> 3. Click **{% icon plus %} New Study** button to confirm
> 4. **Expand** the new study folder that was created in your ARC
>    - you should see the protocols and resources folders created
>      ![screenshot of the ARC structure after adding the new study](images/arcitect-new-study-files.png)
{: .hands_on}

As we did for the investigation, we can also add some basic study information

> <hands-on-title>Add basic study metadata</hands-on-title>
>
> 1. Click on your study `talinum_drought`
> 2. The main panel will now show study metadata, here you can add basic information about your study
>    such as description, contacts, associated publications.
>
>    ![screenshot of the study metadata screen](images/arcitect-study-metadata.png)
>
> 3. **Add a description** for the study
> 4. **Add yourself** as a contact
>    - Remember that you can use your ORCiD ID for a fast way to fill this data
> 5. **Add a publication**
>    - Normally you would add the paper you published based on the data in this ARC
>    - For this tutorial, add any publication (e.g. one of yours, or one you find interesting)
>    - The easiest way to do this, is to add the **DOI link** or **Pubmed ID**, then the information will be autofilled
>    - If you don't have a publication in mind, try one of the following papers:
>      - Pubmed ID: `7825135`
>      - DOI: `https://doi.org/10.1136/bmj.331.7531.1498`
>
>    ![Example of study metadata filled](images/arcitect-study-metadata-filled.png)
{: .hands_on}

Since we have once again completed a unit of work, let's save and sync our ARC again!

> <hands-on-title>Save & Sync</hands-on-title>
>
> 1. Save your ARC locally
>    - **Save ARC** button in menu
> 2. Sync your ARC to DataHUB
>    - **DataHub Sync** button in menu
>    - Commit all changed files
>    - Provide a good commit message (e.g. `added drought study`)
{: .hands_on}


### Create an assay

Next, we will create our 2 assays. Recall that we performed two measurements/assays, RNA sequencing, and metabolomics.

> <hands-on-title>Create the RNA-seq assay</hands-on-title>
>
> 1. Click on the **{% icon plus %} (plus) icon** next to the assays folder
> 2. Give your assay a name, e.g. `rnaseq`
> 3. In **measurement type** field, choose `mRNA Sequencing`
>    - this field is **ontology based**
>    - as soon as you start typing, suggetions for known terms will pop up
>    - Notice the TSR and TAN fields are automatically filled in once you select a term. These fields help make your ARC machine readable.
>
> 4. Similarly, set the  **Technology Type** to `Next-generation Sequencing`
>
>    > <question-title></question-title>
>    > Look through our demo data, can you find out which sequencing machine was used?
>    > > <solution-title></solution-title>
>    > > In the example data, we've got a **methods** folder where we've stored notes on our assays,
>    > > In the file called **illumina_libraries** we see that the sequencer we used was a `Illumina HiSeq 2000`
>    > {: .solution}
>    {: .question}
>
> 5. Under **Technology platform**, fill in `Illumina HiSeq 2000`
>
>    ![screenshot of assay metadata filled](images/arcitect-assay.png)
>
{: .hands_on}


Awesome, now let's repeat the process for our metabolomics assay.

> <hands-on-title>Save & Sync</hands-on-title>
>
> 1. Create a new assay called `metabolomics`
> 2. Under **measurement type** choose `Gas Chromatography Mass Spectrometry`
> 3. We can leave the Technology type/platform fields empty for now
>
>   ![screenshot of metabolomics assay](images/arcitect-assay2.png)
{: .hands_on}


And once again, let's save our ARC and sync it to DataHUB.

> <hands-on-title>Save & Sync</hands-on-title>
>
> 1. Save your ARC locally
>    - **Save ARC** button in menu
> 2. Sync your ARC to DataHUB
>    - **DataHub Sync** button in menu
>    - Commit all changed files
>    - Provide a good commit message (e.g. `added rnaseq and metabolomics assays`)
{: .hands_on}



## Adding files to your ARC

As our next step, let's add the files we have in our demo folder to our ARC.


> <hands-on-title>Adding protocols</hands-on-title>
>
> 1. Open the folder with demo data on your computer
> 2. Look through the **methods** folder
>    - This contains files with notes about methods we used in our research
>    - In our ARC, these usually go in **protocols** folders
>    - Studies and assays both have protocols folders
>
>    > <question-title></question-title>
>    >
>    > Where should each of the 4 files in the methods folder go in our ARC?
>    >
>    > > <solution-title></solution-title>
>    > > - `plant_material.txt`: contains a general sample sheet, and should go in the study-level protocol folder
>    > > - `RNA_extraction.txt` and `Illumina_libraries.txt` are specific to the RNA sequencing assay, so should go in the protocols folder for the rnaseq assay
>    > > - `metabolite_extraction.txt` is specific to our metabolomics assay, so should go in the metabolomics protocol folder
>    > {: .solution}
>    {: .question}
>
> 3. There are two ways to add files to your ARC, we will do both now, in the future you can choose whichever option you prefer
>
> 4. **Adding files via ARCitect**
>    - Expand the `studies` folder
>    - Expand the `talinum_drought` folder
>    - Right-click on the `protocols` folder
>    - Choose `Import files` option
>    - Browse to the demo data and choose the `methods/plant_material.txt` file
>
>    ![screenshot of menu to add documents](images/arcitect-import-files.png)
>
> 5. **Adding files via file explorer**
>    - Click on `Explorer` in the left-hand menu
>    - Open the `assays/metabolomics/protocols` folder of your ARC
>    - In another window, open the demo data folder
>    - Copy or drag the `metabolite_extraction.txt` file from the demo data folder to your ARC folder
>
> 6. **Add the two remaining protocol files** to the ARC in your preferred way
>    - Files `illumina_libraries.txt` and `RNA_extraction.txt`
>    - Both should go to `assays/rnaseq/protocols/` folder in your ARC
>
{: .hands_on}

Note that adding files via ARCitect import option will create a copy of the file, rather than moving it, so if you have large data files you might want to move the files via the explorer option.

Next, let's add the RNA-seq data we have to the ARC structure.

> <hands-on-title>Add RNA-seq data</hands-on-title>
>
> 1. Look in the `rnaseq_data` folder of the demo data
>    - Files with the `.fastq.gz` suffix are raw sequencing files
>    - Raw data goes in the `assays/rnaseq/dataset/` folder
>    - `NGS_samplesheet.xlsx` can go into the assay's protocol folder
>
> 2. Move the files to the ARC in your preferred way
>
>    ![screenshot of arc explorer panel after adding files](images/arcitect-assay-files-rnaseq.png)
>
{: .hands_on}

Finally, let's do the same for our metabolomics data

> <hands-on-title>Add Metabolomics data</hands-on-title>
>
> 1. Look in the `metabolomics_data` folder of the demo data
>    - Folders with the `.D` suffix are raw data files
>    - The rest of the files can go to the protocols folder
>
> 2. Move all the files to the ARC metabolomics assay folder
>    - Note that if you move them via ARCitect, you should use the **Import Directories** option for the raw data folders
>
>    ![screenshot of arcitect files structure after adding metabolomics files](images/arcitect-assay-files-metabolomics.png)
>
{: .hands_on}



Congrats! you have added data to your ARC! As always, let's make sure to save & sync our work before we continue.

> <hands-on-title>Save & Sync</hands-on-title>
>
> 1. Save your ARC locally
> 2. Sync your ARC to DataHUB
>    - Provide a good commit message (e.g. `added protocols and raw data files`)
{: .hands_on}


## Add your experimental metadata (SWATE)

Now that we have all our data in the ARC, we need to add metadata about or research in a structured way.
We already put a lot of information about our samples and experiments and assays in the protocols folders,
but these are just free-text notes; the information is not structured yet.

The protocols folder is a create place

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






