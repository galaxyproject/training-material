---
layout: tutorial_hands_on

title: Text-Mining articles with GROBID
level: Introductory
draft: true
questions:
  - What is GROBID?
  - How to get started with GROBID in Galaxy?
objectives:
  - Log in to Galaxy
  - Upload files to the platform
  - Use GROBID to extract information from publications
time_estimation: 1H
key_points:
  - You can use GROBID from within Galaxy to mine text.
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
[GROBID](https://github.com/grobidOrg/grobid) stands for **G**ene**R**ation **O**f **BI**bliographic **D**ata. It is a text mining library for large-scale extraction of bibliographical metadata and is now available on Galaxy.

# Why use GROBID in Galaxy?
[Galaxy](https://galaxyproject.org/) is an open-source data analysis platform. 
Using GROBID on Galaxy gives you access to 4000+ tools to work with - all without programming skills. 
Moreover, Galaxy allows you to leverage high-performance computing (HPC) resources from within your browser.
This is particularly useful if you are working on large amounts of text.

If your texts are not machine-readable, you could use Tesseract directly in Galaxy, or perform further text-mining tasks or visualise your results after running GROBID. All directly in Galaxy. 
You can chain different tools into [workflows](https://gxy.io/GTN:T00151) that you can rerun easily, allowing you to automate your analysis steps and make your tasks shareable.

Additionally, [FAIR](https://galaxyproject.org/fair/) and reproducible research are at the core of what we do. 
Galaxy supports best practices in [Research Data Management](https://gxy.io/GTN:T00575) (RDM) across all steps of the research data life cycle and makes reproducibility a built-in default to support your research. Learn more on Galaxy in our [Introduction to Digital Humanities](https://gxy.io/GTN:T00557) and more on the platform's RDM features in the [Introduction to Galaxy as an RDM platform](https://gxy.io/GTN:T00575). 

# Get started in Galaxy

## Create an account on Galaxy
To use Galaxy's full potential, you must register and create an account. You can skip this step if you already have a Galaxy account.

{% snippet faqs/galaxy/account_create.md %}

Alternatively, you can access Galaxy using a single sign-on of your choice, for example, from [IAM4NFDI](https://iam.services.base4nfdi.de/faq_ENG/) on [Galaxy Europe](https://usegalaxy.eu/). The interface could look like this:

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
> > If renaming does not work, you may not be logged in, so try logging in to Galaxy first. Anonymous users can have only one history and cannot rename it.
> {: .comment}
>
{: .hands_on}

## Upload a file to Galaxy

There are many ways to get data into Galaxy.
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
> 4. Paste in the address of the files you want to upload here.
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
> 1. Click the {% icon galaxy-eye %} (eye) icon next to the dataset name to look at the file content.
>
>    ![galaxy history view showing a single dataset. Display link is being hovered.]({% link topics/digital-humanities/images/History-view.png %})
>
> 2. Check the datatype - is it **PDF**? Then you are all set. Otherwise, adapt the datatype.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}


If you want to work with multiple articles at once, we suggest uploading them all at once and creating a **collection** from them.

When does it make sense for you to create a collection?

{% snippet faqs/galaxy/histories_datasets_vs_collections.md %}

How can you create collections from single datasets?

{% snippet faqs/galaxy/collections_autobuild_list.md %}

Your articles are now on Galaxy, great!
The next steps now depend on your article: was it born-digital? Then you can jump ahead to using GROBID.

Is it a scanned text? In that case, your text will not be machine-readable or searchable.
In this case, you need to first perform optical character recognition (OCR) to properly use GROBID with your article.
This is possible in Galaxy in various ways; you could, for example, use {% tool [Tesseract](toolshed.g2.bx.psu.edu/repos/iuc/tesseract/tesseract/5.5.2+galaxy4) %} on Galaxy and make your PDF machine-readable.

Once your article is ready, we can run GROBID on it.

# Run GROBID on Galaxy

Navigate to the activity bar on the left-hand side and click on {% icon tool %} tools. Enter **GROBID** in the search bar.
Click on the search result {% tool [GROBID extract](toolshed.g2.bx.psu.edu/repos/iuc/grobid_grobid/grobid_grobid/0.9.0+galaxy0) %} (or this link) to see the tool.

Now, you can select how to process your PDF or PDFs. 

> <hands-on-title>Processing files with GROBID</hands-on-title>
> 1. Click on the dropdown menu to select **Scientific article PDF to TEI**
>
> 2. Scientific article PDF: Select the article or articles you want to process. You can select single datasets, multiple datasets and a Dataset collection.
>     You can use the dropdown menu or drag and drop them from your history on the right-hand side.
>     ![Screenshot of GROBID with tool parameters filled]({% link topics/digital-humanities/images/grobid.png %})
>
> 3. Optionally, you can click on **Advanced options** to further specify your output.
>    
>
> 4. Once you are done, click on the {% icon workflow-run %} **Run Tool** button at the top or bottom of the tool to run it.
>
{: .hands_on}

The job is now running. 
It will create a new dataset in your history, on top of your upload. 
The dataset will go through the grey (preparing/queued) and yellow (running) states to become green (success). 

Once it is green, you can click on the  {% icon galaxy-eye %} (eye) icon next to the dataset name to look at the file content.
In this example, we chose to convert the article to TEI.

> <comment-title>What is TEI?</comment-title>
> TEI, or the **T**ext **E**ncoding **I**nitiative is an international project to develop guidelines for "machine-actionable cultural heritage texts."
> For more details and guidelines, visit the [TEI website](https://tei-c.org/).
>
{: .comment}

The output of your file is an annotated TEI file.

Should you notice that you made a mistake or selected the wrong input file or output format, do not worry.
It is easy to re-run jobs in Galaxy:

{% snippet faqs/galaxy/tools_rerun.md %}

This allows you to see exactly which parameters you used when you ran the job before and adapt them if needed.
Make your changes, then click **Run Tool**.
This creates a new dataset with updated parameters.

In the next step, we want to extract the datasets mentioned in the scientific articles. 
We will use GROBID DataStet to achieve this:
Navigate to the activity bar on the left-hand side once more and click on {% icon tool %} tools. Enter **GROBID** in the search bar.
Click on the search result {% tool [GROBID DataStet](toolshed.g2.bx.psu.edu/repos/iuc/grobid_dataset_annotate/grobid_dataset_annotate/0.8.1+galaxy0) %} (or this link).

> <hands-on-title>Using GROBID DataStet</hands-on-title>
> 1. Scientific article PDF or TEI XML: Select the TEI output you created earlier. 
>    If you have not renamed it, it will be called **GROBID on dataset** or something similar.
>
> 2. Output format: Select **Enriched TEI XML** in the dropdown menu.
>     ![Screenshot of GROBID  Datastet with tool parameters filled]({% link topics/digital-humanities/images/datastet.png %}) 
>
> 3. Click on the {% icon workflow-run %} **Run Tool** button at the top or bottom of the tool to run it.
>
{: .hands_on}

**Beware**: If you want a JSON output, you need to run it on the whole PDF. 
In this case, selecting both TEI and JSON will result in an error.
For more information, see the tool help.

{% snippet faqs/galaxy/tools_help.md %}

The job starts running. Once it turns green, it is ready for you to inspect.
Click on the {% icon galaxy-eye %} (eye) icon next to the dataset name to look at the new file's content.
The result is an enriched TEI XML file.
Depending on your needs, you may already be finished or want to use additional tools to achieve your goals.
You could now extract specific lines from this (and all the other created TEI files) by using {% tool [Select lines that match an expression](Grep1) %} or use different tools to achieve your goal.

Give it a try and see how Grobid can help you mine your articles!
