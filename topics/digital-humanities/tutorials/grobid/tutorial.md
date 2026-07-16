---
layout: tutorial_hands_on

title: Text mining articles with GROBID
level: Introductory
questions:
  - What is GROBID
  - How to get started with GROBID in Galaxy?
objectives:
  - Log in to Galaxy
  - Upload files to the platform
  - Use GROBID
time_estimation: 1H
key_points:
  - You can use GROBID from within Galaxy
tags:
- text mining
contributions:
  authorship:
    - Sch-Da
  funding:
    - deKCD
    - mwk
    - deNBI
---

This tutorial shows you how you can use GROBID on Galaxy.
The first couple of steps derive from [A short introduction to Galaxy]({% link topics/introduction/tutorials/galaxy-intro-short/tutorial.html %}).

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# What is GROBID?
GROBID stands for **G**ene**R**ation of **BI**bliographic **D**ata. It is a text mining library for large-scale extraction of bibliographical metadata.


# Get started in Galaxy

## Create an account on Galaxy
To use Galaxy's full potential, you must register and create an account. You can skip this step if you already have a Galaxy account.

{% snippet faqs/galaxy/account_create.md %}

Alternatively, you can log in using a single sign-on of your choice, for example, from [IAM4NFDI](https://iam.services.base4nfdi.de/faq_ENG/) on [Galaxy Europe](https://usegalaxy.eu/).

 ![Screenshot of Galaxy Europe register window with the IAM4NFDI login button highlighted](../../images/iam4nfdi.png)

## Log in to Galaxy

> <hands-on-title>Log in to Galaxy</hands-on-title>
> 1. Open your favourite browser (Chrome, Safari, Edge or Firefox as your browser, not Internet Explorer!)
> 2. Browse to your Galaxy instance, for example [Galaxy Europe](https://usegalaxy.eu/)
> 3. Log in with your credentials
>
> ![Screenshot of Galaxy Europe with the register or login button highlighted](../../images/Galaxyeulogin.png)
>
>   > <comment-title>Different Galaxy servers</comment-title>
>   >  This is an image of Galaxy Europe, located at [usegalaxy.eu](https://usegalaxy.eu)
>   >
>   > The particular Galaxy server you are using may look slightly different and have a different web address.
>   >
>   > You can also find more possible Galaxy servers at the top of this tutorial in **Available on these Galaxies**
>   {: .comment}
{: .hands_on}

The Galaxy homepage is divided into four sections (panels):
* The Activity Bar on the left: _This is where you will navigate to the resources in Galaxy (Tools {% icon tool %}, Workflows {% icon galaxy-workflows-activity %}, Histories {% icon galaxy-history-storage-choice %}, etc.)_
* Currently active "Activity Panel" on the left: _By default, the {% icon tool %} **Tools** activity will be active and its panel will be expanded_
* Viewing panel in the middle: _The main area for context for your analysis_
* History of analysis and files on the right: _Shows your "current" history; i.e.: Where any new files for your analysis will be stored_

![Screenshot of the Galaxy interface with aforementioned structure]({% link topics/introduction/images/galaxy_interface.png %})

The first time you use Galaxy, your history panel is empty.

## Name your current history

Your "History" is on the panel on the right. It is a record of the actions you have taken.

> <hands-on-title>Name history</hands-on-title>
> 1. Go to the **History** panel (on the right)
> 2. Click {% icon galaxy-pencil %} (**Edit**) next to the history name (which by default is "Unnamed history")
>
>    ![Screenshot of the galaxy interface with the history name being edited, it currently reads "Unnamed history", the default value. An input box is below it.]({% link shared/images/rename_history.png %}){:width="250px"}
>
>    > <comment-title></comment-title>
>    >
>    > In some previous versions of Galaxy, you will need to click the history name to rename it as shown here:
>    > ![Screenshot of the galaxy interface with the history name being edited, it currently reads "Unnamed history", the default value.](../../../../shared/images/rename_history_old.png){:width="320px"}
>    {: .comment}
>
> 3. Type in a new name, for example, "Testing GROBID"
> 4. Click **Save**
>
> > <comment-title>Renaming not an option?</comment-title>
> > If renaming does not work, you may not be logged in, so try logging in to Galaxy first. Anonymous users can only have one history, and they cannot rename it.
> {: .comment}
>
{: .hands_on}

## Upload a file to Galaxy

There are many different ways of getting data to Galaxy.
The "Activity Bar" is located on the leftmost part of the interface. It shows you various options.

For this tutorial, we suggest you use a scientific article of your choice. Make sure that the article provider allows text mining before you take the next steps.

> <hands-on-title>Upload a file</hands-on-title>
> 1. At the top of the **Activity Bar**, click the {% icon galaxy-upload %} **Upload** activity
>
>    ![upload data button shown in the galaxy interface]({% link topics/introduction/images/upload-data.png %})
>
>    This brings up a box:
>
>    ![the Galaxy upload dialogue, the 'regular' tab is active with a large textarea to paste subsequent URL]({% link topics/introduction/images/upload-box.png %})
>
> 3. Click **Paste/Fetch data**
> 4. Paste in the address of the files you want to upload in here.
> 5. Click **Start**
> 6. Click **Close**
{: .hands_on}


Your uploaded file is now in your current history.
When the file has been uploaded to Galaxy, it will turn green.

> <comment-title></comment-title>
> After this, you will see your first history item (called a "dataset") in Galaxy's right panel. It will go through
> the grey (preparing/queued) and yellow (running) states to become green (success).
>
{: .comment}

The contents of the file will be displayed in the central Galaxy panel. If the dataset is large, you will see a warning message which explains that only the first megabyte is shown.

> <hands-on-title>View the text files content</hands-on-title>
> 1. Click the {% icon galaxy-eye %} (eye) icon next to the dataset name, to look at the file content
>
>    ![galaxy history view showing a single dataset mutant_r1.fastq. Display link is being hovered.](../../images/eye-icon.png){:width="520px"}
>
> 2. Check the datatype - is it **PDF**? Then you are all set. Otherwise, adapt the datatype.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}


If you want to work with multiple articles at once, we suggest uploading all of them at once and creating a **collection** out of them.

When does it make sense for you to create a collection?

{% snippet faqs/galaxy/histories_datasets_vs_collections.md box_type="none"%}

How can you create collections from single datasets?

{% snippet faqs/galaxy/collections_autobuild_list.md %}

Your articles are now on Galaxy, great!
The next steps now depend on your article: was it born-digital? Then you can jump ahead to using GROBID.

Is it a scanned text? In that case, your text will not be machine-readable or searchable.
In this case, you need to first perform optical character recognition (OCR) to properly use GROBID with your article.
This is possible in Galaxy in various ways, you could, for example, use {% tool [Tesseract](toolshed.g2.bx.psu.edu/repos/iuc/tesseract/tesseract/5.5.2+galaxy4) %} on Galaxy an make your PDF machine readable.

Once your article is ready, we can run GROBID on it.



