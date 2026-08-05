---
layout: tutorial_hands_on
title: Accessing IIIF from within Galaxy
level: Introductory
questions:
  - What is IIIF?
  - How can you add existing IIIF repositories to Galaxy?
objectives:
  - Import datasets to Galaxy using the IIIF integration in Galaxy
requirements:
  - type: "internal"
    topic_name: digital-humanities
    tutorials:
      - introduction_to_dh
time_estimation: 15M
key_points:
  - With the IIIF integration in Galaxy, you can directly browse IIIF from within Galaxy, eliminating manual download cycles.
contributions:
  authorship:
    - Sch-Da
  editing:
    - gautschr
  funding:
    - deKCD
    - mwk
    - deNBI
---

## What is IIIF?

IIIF stands for Image Interoperability Framework. According to [their website](https://iiif.io/), "IIIF (pronounced triple eye eff) is a set of open standards for delivering high-quality, attributed digital objects online at scale."
It comes with an international community and is backed by cultural institutions. 
IIIF allows you to access high-quality images quickly and flexibly.

By connecting Galaxy to IIIF, researchers can directly browse IIIF resources and import images into their Galaxy history for further analysis.

## Adding IIIF repositories to Galaxy

**Please note:** the Galaxy IIIF integration is currently only available on [usegalaxy.eu](https://usegalaxy.eu).

Galaxy integrates IIIF through its File Sources (Remote Files) framework, enabling researchers to search, browse, and import images from remote IIIF repositories directly into an active Galaxy history without first downloading them locally.

### Adding and Browsing your IIIF repository

> <hands-on-title>Add your IIIF manifest to Galaxy</hands-on-title>
>
> 1. In Galaxy, navigate to the top menu, **click on your user name** and then **Preferences**.
> 2. Select **Manage your Repositories** and click the **Create** button.
> 3. Scroll down and select **IIIF Manifest**.
> 4. Complete the configuration form:
>  * **Name:** A descriptive title (e.g., *Bodleian library*).
>  * If you want, you can add an optional description.
>  * **IIIF Manifest URL:** The full base URL of the IIIF manifest (e.g., `https://iiif.bodleian.ox.ac.uk/iiif/collection/top`).
> 5. Click **Create**.
>
{: .hands_on}

After saving, your new entry appears in **My Repositories**.

![My Repositories] {% link topics/galaxy-interface/tutorials/iiif/bodleian.png %}

Once you have added your IIIF Manifest, you can access its content and load images into your history for further processing.

> <hands-on-title>Importing Discovered Files into Your History</hands-on-title>
>
> Once the connection to your IIIF manifest is established, navigating and uploading the files follows a standardised routine:
>
> 1. Open the **Upload Data** menu and click **Choose from repository**.
> 2. Select your newly configured IIIF source.
> 3. **Browse or Search:** You can navigate the hierarchical collections or use the dataset-name search filter.
> 4. **Select Files:** Once you find the correct folder, click into it, tick the checkbox next to the desired files and hit **Select**.
> 5. **Execute Import:** The selected items will populate the standard Galaxy upload queue. Click **Start** to run the import process.
>
> Once the progress bar turns green, the files will appear in the right-hand history panel as standard Galaxy datasets, ready to be routed into any analytical tool or workflow.
>
{: .hands_on}

Depending on your needs, you can crop or preprocess the uploaded images in Galaxy, segment them or use them with other tools, for example, for optical character recognition. 
Why not try adding your favourite IIIF images now?
