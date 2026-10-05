---
layout: tutorial_hands_on
title: Introduction to the Dataverse Integration in Galaxy
level: Introductory
questions:
  - What is Dataverse?
  - How can you browse existing Dataverse repositories from within Galaxy?
  - How can you integrate a private Dataverse repository?
objectives:
  - Import datasets to Galaxy using the Dataverse integration in Galaxy
time_estimation: 1H
key_points:
  - With the Dataverse integration in Galaxy, you can directly browse and import public or private datasets, eliminating manual download cycles.
contributions:
  authorship:
    - anuprulez
    - dadrasarmin
  editing:
    - Sch-Da
---

## What is Dataverse?

At its core, **[Dataverse](https://dataverse.org/)** is an open-source web application designed to share, preserve, cite, explore, and analyze research data. Developed and led by the [Institute for Quantitative Social Science (IQSS)](https://www.iq.harvard.edu/) at [Harvard University](https://www.harvard.edu/), it serves as a robust repository platform used by universities, research institutions, and laboratories worldwide to host data across all scientific disciplines from genomics and structural biology to the social sciences.

For researchers, Dataverse acts as a digital archive that ensures data is not just stored, but made visible and accessible to the global scientific community. It automates the creation of standard academic citations and permanent identifiers, making data sharing an easy part of the research lifecycle.

## How Dataverse Organizes Data

To effectively use Dataverse, it helps to understand its hierarchical structure. Think of it as a nested system that moves from institutional hosting down to individual files:

* **[Dataverse Collections (or "Dataverses")](https://guides.dataverse.org/en/6.10.1/user/dataverse-management.html):** A container that can hold other dataverses or datasets. An institution, department, laboratory, or individual researcher can have their own dedicated Dataverse collection, which they can brand, manage, and configure independently.
* **[Datasets](https://guides.dataverse.org/en/6.10.1/user/dataset-management.html):** A dataset in Dataverse is a logical grouping of research data. It contains the actual data files, documentation, code, and the crucial **metadata** that describes the data.
* **[Files](https://guides.dataverse.org/en/6.10.1/user/tabulardataingest/index.html):** These are the actual building blocks uploaded by the researcher, such as text files, CSV spreadsheets, FASTQ sequence files, imaging data, or R/Python analysis scripts.

This flexible hierarchy means a single institution can manage its entire data output under one parent Dataverse collection, while creating sub-collections for individual projects, grants, or researchers.

## Why Use Dataverse?

Dataverse is built from the ground up around the **[FAIR Data Principles** (Findable, Accessible, Interoperable, and Reusable)](https://www.go-fair.org/fair-principles/). It achieves this through several important features:

* **[Persistent Identifiers](https://en.wikipedia.org/wiki/Persistent_identifier):** Every time a dataset is published, Dataverse automatically assigns it a [Digital Object Identifier (DOI)](https://en.wikipedia.org/wiki/Digital_object_identifier). This ensures the dataset can be permanently cited in academic journals, giving authors proper academic credit.
* **Rich Metadata Standards:** Dataverse enforces standard metadata schemas (like [Dublin Core](https://en.wikipedia.org/wiki/Dublin_Core) or [DDI](https://ddialliance.org/ddi-codebook)) and offers domain-specific metadata fields (such as life sciences or geospatial data). This makes your data highly searchable and discoverable by global search engines like Google Dataset Search.
* **Rigorous [Version Control](https://en.wikipedia.org/wiki/Version_control):** Research is iterative. Dataverse allows you to update datasets while maintaining a clear, public history of previous versions. If you add new samples or fix a data error, users can see exactly what changed between version 1.0 and 2.0.
* **Customizable [Access Control](https://en.wikipedia.org/wiki/Computer_access_control):** Depositing data doesn't mean losing control. Dataverse allows you to keep datasets restricted during a peer-review process, requires users to sign data use agreements, or release data fully into the public domain via [Creative Commons](https://en.wikipedia.org/wiki/Creative_Commons) (e.g., CC0) waivers.

## Galaxy and Dataverse: Complementary Roles

In computational workflows, data generation, analysis, and archiving go hand in hand. While **Galaxy** is a platform for executing complex workflows and managing interactive data analysis, it is not designed to act as a permanent data repository.

By connecting Galaxy to Dataverse, researchers can directly export their finalized analysis histories, processed datasets, and verified workflows into a formal Dataverse repository. This creates an end-to-end, reproducible research pipeline: you import raw data, analyze it transparently in Galaxy, and deposit the verified results directly into Dataverse, ready for publication and community reuse. This data could even be accessed again with Galaxy for further analysis at a later stage.

## How can you browse existing Dataverse repositories from within Galaxy?

**Please note:** the Galaxy Dataverse integration is currently only available on [usegalaxy.eu](https://usegalaxy.eu).
Galaxy has integrated Dataverse directly into its **File Sources (Remote Files)** framework. This allows researchers to search, browse, and import remote repository datasets directly into an active Galaxy history without downloading them locally first.

### Method 1: Browsing Pre-configured Public Dataverses

Galaxy instances frequently maintain a list of pre-configured public repositories for major data platforms.

> <hands-on-title>Browse Pre-configured Dataverses</hands-on-title>
>
> 1. Click the **Upload Data** icon at the top of the Galaxy tool panel.
> 2. Inside the upload dialog box, select the **Choose from repository** button at the bottom.
> 3. In the search/filter box at the top of the repository browser, type `Dataverse`.
> 4. A list of globally available, pre-configured Dataverse instances will appear. Click on an instance to explore its public directories, datasets, and files. It will look similar to this:
>
> ![Pre-Configured Dataverse](https://galaxyproject.org/images/news/2026-01-15-dataverse/preconfigured_dataverse.png)
>
{: .hands_on}

From here, you can select your instance and browse or search for existing files. **Note**: the search is case-sensitive. 
If you have already found the file you want to work with, scroll down to "Importing Discovered Files into Your History".

### Method 2: Adding and Browsing a Custom Dataverse Repository

If you need to connect to an institutional Dataverse instance not listed by default, or if you need to access your own private datasets, you must first link your account using an [API](https://en.wikipedia.org/wiki/API) token.

> <hands-on-title>Generate an API Token in Dataverse</hands-on-title>
>
> 1. Log into the target Dataverse instance (e.g., [Harvard Dataverse](https://dataverse.harvard.edu/)).
> 2. Click on your account name and select **API Token**.
> 3. Click **Create Token** and **Copy to Clipboard** to copy the token string.
{: .hands_on}

> <hands-on-title>Link the Repository in Galaxy</hands-on-title>
>
> 1. In Galaxy, navigate to the top menu and **click on your user name** and then **Preferences**.
> 2. Select **Manage your Repositories** and click the **Create** button.
> 3. Select **Dataverse**.
> 4. Complete the configuration form:
>  * **Name:** A descriptive title (e.g., *My Institutional Dataverse*).
>  * **Dataverse instance endpoint:** The full base URL of the repository (e.g., `https://dataverse.harvard.edu`).
>  * **API Token:** Paste the token generated in Step 1.
>  * **Allow Galaxy to export data:** Set to **Yes** if you want the flexibility to push finalized Galaxy histories back to Dataverse later.
> 5. Click **Create**.
>
{: .hands_on}

The remote repository configuration menu looks like this:

![Remote Repositories](https://galaxyproject.org/images/news/2026-01-15-dataverse/bsc_dv_1.png)

Your completed form could look like this (but with information of the desired Dataverse):

![Remote Repositories filled form](https://galaxyproject.org/images/news/2026-01-15-dataverse/bsc_dv_2.png)

After saving, your new entry appears in **My Repositories**.

![My Repositories](https://galaxyproject.org/images/news/2026-01-15-dataverse/bsc_dv_3.png)

> <hands-on-title>Importing Discovered Files into Your History</hands-on-title>
>
> Once the connection is established (via Method 1 or Method 2), navigating and uploading the files follows a standardized routine:
>
> 1. Open the **Upload Data** menu and click **Choose from repository**.
> 2. Select your newly configured or pre-configured Dataverse source.
> 3. **Browse or Search:** You can navigate through the hierarchical Dataverse collections or use the search filter for dataset names.
> 4. **Select Files:** Once you find the correct dataset container, click into it, tick the checkbox next to the desired files and hit **Select**.
> 5. **Execute Import:** The selected items will populate the standard Galaxy upload queue. Click **Start** to run the import process.
>
> Once the progress bar turns green, the files will appear in the right-hand history panel as standard Galaxy datasets, ready to be routed into any analytical tool or workflow.
>
{: .hands_on}
