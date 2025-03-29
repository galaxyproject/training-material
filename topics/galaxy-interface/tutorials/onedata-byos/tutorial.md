---
layout: tutorial_hands_on

title: "Onedata user-owned storage"
subtopic: upload
tags:
  - storage
contributors:
  - lopiola
  - bwalkowi

time_estimation: "15m"
level: Introductory

requirements:
 - type: "internal"
   topic_name: galaxy-interface
   tutorials:
     - onedata-getting-started

questions:
- How to use Onedata as a Storage Location for Galaxy datasets?
- What is the difference between Remote File Source and Storage Location?
- What permissions are required for Onedata Storage Location?
objectives:
- Configure Onedata as a Storage Location for Galaxy datasets.
- Learn how to manage Storage Location preferences.
- Understand the requirements and implications of using Onedata storage.
key_points:
- Storage Location stores Galaxy datasets directly in your Onedata Space.
- Write access to the Onedata Space is required.
- Storage Location is different from Remote File Source - it acts as an Object Store.
- You can use multiple Storage Locations and set preferences between them.
---


> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Prerequisites

1. This tutorial assumes that you have basic knowledge about Onedata and access 
   to a Onedata ecosystem. If needed, follow 
   [this tutorial]({% link topics/galaxy-interface/tutorials/onedata-getting-started/tutorial.html %})
   first!
2. To use Onedata as user-owned storage, you need the **domain** of the
   **Onezone service**, a **Space** name, and a suitable **access token**. Here
   is the relevant 
   [guide on how to get them]({% link topics/galaxy-interface/tutorials/onedata-getting-started/tutorial.html#accessing-onedata-services %}). 
3. The Galaxy server must be properly configured by the admins for the Onedata
   Storage Location (BYOS) templates to be available. Here is the 
   [corresponding tutorial]({% link topics/admin/tutorials/onedata-configuration/tutorial.html %}).


# Introduction

Your account in Onedata can be used to store your Galaxy datasets, effectively
extending your quota. Onedata then serves as a so-called **Storage Location**,
acting as an **Object Store** for the Galaxy server. Below hands-on tutorial
will help you configure and use a Onedata Storage Location using the so-called
Bring Your Own Storage (BYOS) approach.


> <tip-title>Disambiguation</tip-title>
> While Onedata can be used for both, a **Storage Location** is not the same as
> a **Remote File Source**. In this tutorial, you will be setting up a
> Onedata-based **Storage Location**, which allows storing your Galaxy datasets
> directly in a Onedata Space in a transparent way. If you are looking to use a
> Remote File Source, refer to 
> [this tutorial]({% link topics/galaxy-interface/tutorials/onedata-remote-import/tutorial.html %}).
{: .tip}


# Configuration

Follow these steps:

1. Log in to your Galaxy account.
2. Go to the **Manage Your Storage Locations** section of the **Preferences** menu.
3. Click **Create** in the top right corner to create a new Storage Location.
   ![Create a Remote File Source](../../images/onedata-byos/configure-create-sl.png)
4. Choose the **Onedata Storage** template. If there is no such template, the 
   Galaxy server is not configured to support it. Consider contacting its admins.
5. Fill in the information, following the hints visible on the form.
   ![Configure Onedata Storage Location](../../images/onedata-byos/configure-onedata-sl.png)
6. Click **Create** to finalize.

You have now configured a new Storage Location, but you still need to tell
Galaxy your preferences so that it is effectively used among different available
options:

1. Go to the **Preferred Storage Location** section of the **Preferences** menu.
2. Choose the newly added Storage Location:
   ![Prefer Onedata Storage Location](../../images/onedata-byos/preferred-storage-location.png)


# See it in action

Upload some new data to Galaxy, or run a workflow to produce results. Then
navigate to your Onedata account (e.g. [https://datahub.egi.eu](https://datahub.egi.eu))
and open the Space (and the path) that you have put down in the config. You
should see the Galaxy data:

![Data in Onedata](../../images/onedata-byos/data-in-onedata.png)


# Troubleshooting

If you experience any problems, take a look at the 
[troubleshooting]({% link topics/galaxy-interface/tutorials/onedata-getting-started/tutorial.html#troubleshooting %}) 
guide.
