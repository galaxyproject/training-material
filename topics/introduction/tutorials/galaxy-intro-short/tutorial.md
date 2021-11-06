---
layout: tutorial_hands_on

title: "A short introduction to Galaxy"
zenodo_link: "https://doi.org/10.5281/zenodo.582600"
level: Introductory
tags:
  - español
questions:
  - "How to get started in Galaxy"
objectives:
  - "Learn how to upload a file"
  - "Learn how to use a tool"
  - "Learn how to view results"
  - "Learn how to view histories"
  - "Learn how to extract and run a workflow"
time_estimation: "40m"
key_points:
  - "The Galaxy interface has tools on the left, viewing pane in the middle, and a history of your data analysis on the right."
  - "You can create a new history for each analysis. All your histories are saved."
  - "To get data into Galaxy, you can upload a file by pasting in a web address. There are other ways to get data into Galaxy (not covered in this tutorial): you can upload a file from your computer, and you can import an entire history."
  - "Choose a tool and change any settings for your analysis."
  - "Run the tool. The output files will be saved at the top of your history."
  - "View the output files by clicking on the eye icon."
  - "View all your histories and move files between them. Switch to a different history."
  - "Log out of your Galaxy server. When you log back in (to the same server), your histories will all be there."
subtopic: core
translations:
  - es
contributors:
  - annasyme
  - nsoranzo

---

# Overview
{:.no_toc}

* This is a short introduction to the Galaxy user interface - the web page that you interact with.
* We will cover key tasks in Galaxy: uploading files, using tools, viewing histories, and running workflows.

> ### Agenda
> 1. TOC
> {:toc}
>
{: .agenda}

## What does Galaxy look like?


> ### {% icon hands_on %} Hands-on: Log in to Galaxy
> 1. Open your favorite browser (Chrome, Safari or Firefox as your browser, not Internet Explorer!)
> 2. Browse to your Galaxy instance
> 3. Log in or register
>
> ![Screenshot of Galaxy Australia with the register or login button highlighted](../../images/galaxy-login.png)
>
>   > ### {% icon comment %} Different Galaxy servers
>   >  This is an image of Galaxy Australia, located at [usegalaxy.org.au](https://usegalaxy.org.au/)
>   >
>   > The particular Galaxy server that you are using may look slightly different and have a different web address:
>   > - The main Galaxy server is [usegalaxy.org](https://usegalaxy.org/)
>   > - The European Galaxy server is [usegalaxy.eu](https://usegalaxy.eu/)
>   >
>   > You can also find more possible Galaxy servers at the top of this tutorial in **Available on these Galaxies**
>   {: .comment}
{: .hands_on}

The Galaxy homepage is divided into three panels:
* Tools on the left
* Viewing panel in the middle
* History of analysis and files on the right

![Screenshot of the Galaxy interface, the tools panel is on the left, the main panel is in the center, and the history is on the right.](../../images/galaxy_interface.png)

The first time you use Galaxy, there will be no files in your history panel.

# Key Galaxy actions

## Name your current history

Your "History" is in the panel at the right.

> ### {% icon hands_on %} Hands-on: Name history
> 1. Go to the **History** panel (on the right)
> 2. Click on the history name (which by default is "Unnamed history")
>
>    ![Screenshot of the galaxy interface with the history name being edited, it currently reads "Unnamed history", the default value.](../../../../shared/images/rename_history.png){:width="320px"}
>
> 3. Type in a new name, for example, "My Analysis"
> 4. Press <kbd>Enter</kbd> on your keyboard to save it
>
> > ### {% icon comment %} Renaming not an option?
> > If renaming does not work, it is possible you aren't logged in, so try logging in to Galaxy first. Anonymous users are only permitted to have one history, and they cannot rename it.
> {: .comment}
>
{: .hands_on}

## Upload a file

Your "Tools" are in the panel at the left.

> ### {% icon hands_on %} Hands-on: Upload a file from URL
> 1. At the top of the **Tools** panel (on the left), click {% icon galaxy-upload %} **Upload**
>
>    ![upload data button](../../images/upload-data.png)
>
>    This brings up a box:
>
>    ![filebox](../../images/upload-box.png){:width="500px"}
>
> 3. Click **Paste/Fetch data**
> 4. Paste in the address of a file:
>
>    ```
>    https://zenodo.org/record/582600/files/mutant_R1.fastq
>    ````
>
> 5. Click **Start**
> 6. Click **Close**
>
{: .hands_on}

Your uploaded file is now in your current history.
When the file has uploaded to Galaxy, it will turn green.

> ### {% icon comment %} Comment
> After this you will see your first history item (called a "dataset") in Galaxy's right panel. It will go through
> the gray (preparing/queued) and yellow (running) states to become green (success).
>
{: .comment}

What is this file?

> ### {% icon hands_on %} Hands-on: View the dataset content
> 1. Click on the {% icon galaxy-eye %} (eye) icon next to the dataset name, to look at the file content
>
>    ![eye icon with dropdown titled View Data](../../images/eye-icon.png){:width="320px"}
{: .hands_on}

The contents of the file will be displayed in the central Galaxy panel.

This file contains DNA sequencing reads from a bacteria, in FASTQ format:

   ![fastq](../../images/fastq.png){:width="620px"}

## Use a tool

Let's look at the quality of the reads in this file.

> ### {% icon hands_on %} Hands-on: Use a tool
> 1. Type **FastQC** in the tools panel search box (top)
> 2. Click on the {% tool [FastQC](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.73+galaxy0) %} tool
>
>    The tool will be displayed in the central Galaxy panel.
>
> 3. Select the following parameters:
>    - {% icon param-file %} *"Raw read data from your current history"*: the FASTQ dataset that we uploaded
>    - No change in the other parameters
> 4. Click **Execute**
>
{: .hands_on}

This tool will run and two new output datasets will appear at the top of your history panel.

## View results

We will now look at the output dataset called *FastQC on data 1: Webpage*.

> ### {% icon comment %} Comment
> * Note that Galaxy has given this dataset a name according to both the tool name ("FastQC") and the input ("data 1") that it used.
> * The name "data 1" means the dataset number 1 in Galaxy's current history (our FASTQ file).
>
{: .comment}


> ### {% icon hands_on %} Hands-on: View results
> * Once it's green, click on the {% icon galaxy-eye %} (eye) icon next to the "Webpage" output dataset.
>
>    The information is displayed in the central panel
>
>    ![fastqc-out](../../images/fastqc-out.png){:width="620px"}
{: .hands_on}

This tool has summarised information about all of the reads in our FASTQ file.

> ### {% icon question %} Questions
>
> 1. What was the length of the reads in the input FASTQ file?
> 2. Do these reads have higher quality scores in the centre or at the ends?
>
>   > ### {% icon solution %} Solutions
>   > 1. 150 bp
>   > 2. In the center
>   {: .solution}
{: .question}


## Run another tool

Let's run a tool to filter out lower-quality reads from our FASTQ file.


> ### {% icon hands_on %} Hands-on: Run another tool
> 1. Type **Filter by quality** in the tools panel search box (top)
> 2. Click on the tool {% tool [Filter by quality](toolshed.g2.bx.psu.edu/repos/devteam/fastq_quality_filter/cshl_fastq_quality_filter/1.0.2+galaxy0) %}
> 3. Set the following parameters:
>    - {% icon param-file %} *"Input FASTQ file"*: our initial FASTQ dataset
>    - *"Quality cut-off value"*: 35
>    - *"Percent of bases in sequence that must have quality equal to / higher than cut-off value"*: 80
> 4. Click **Execute**
{: .hands_on}

After the tool has run, its output dataset will appear at the top of your History panel.
* This dataset will be called "Filter by quality on data 1".
* Remember that Galaxy has named this file according to the tool it used ("Filter by quality") and the input dataset ("data 1").
* The actual numbers in front of the datasets in the history are not important.

What are the results from this filtering tool?

We could click on the eye icon to view the contents of this output file, but it will not be very informative - we will just see a list of reads.

> ### {% icon hands_on %} Hands-on: Get metadata about a file
> 1. Click on the output dataset name in the History panel.
>
>    This expands the information about the file.
>
>    ![filter1 with 3 arrows click on filename to expand, click filter settings and how many reads were filtered out](../../images/filter-fastq1.png)
>
{: .hands_on}

> ### {% icon question %} Questions
>
> How many read has been discarded
>
>   > ### {% icon solution %} Solutions
>   > 1786 low-quality reads were discarded
>   {: .solution}
{: .question}

## Re-run that tool with changed settings

We can now try to filter our input reads to an even higher standard, and see how this changes the resulting output (an exploratory analysis). We will change the filter settings and re-run the tool.

> ### {% icon hands_on %} Hands-on: Re-run the tool
> 1. Click on the {% icon galaxy-refresh %} icon (**Run this job again**) for the output dataset of **Filter by quality** {% icon tool %}
>
>    ![rerun the job](../../images/rerun.png)
>
>    This brings up the tool interface in the central panel with the parameters set to the values used previously to generate this dataset.
>
> 2. Change the settings to something even stricter
>
>    For example, you might decide you want 80 percent of bases to have a quality of 36 or higher, instead of 35.
>
> 3. Click **Execute**
> 4. View the results: Click on the output dataset name to expand the information. (*Note*: not the {% icon galaxy-eye %} (eye) icon.)
{: .hands_on}

> ### {% icon question %} Questions
>
> How many reads were discarded under these new filtering conditions?
>
{: .question}

You can re-run a tool many times with different settings. Each time you re-run the tool, its new output datasets will appear at the top of your current history.


## Convert your analysis history into a workflow

When you look carefully at your history, you can see that it contains all the steps of our analysis, from the beginning to the end. By building this history we have actually built a complete record of our analysis with Galaxy preserving all parameter settings applied at every step. But when you need to analyze new data, it would be tedious to do each step over again. Wouldn't it be nice to just convert this history into a workflow that we will be able to execute again and again?

Galaxy makes this very easy with the `Extract workflow` option. This means any time you want to build a workflow, you can just perform the steps once manually, and then convert it to a workflow, so that next time it will be a lot less work to do the same analysis.

> ### {% icon hands_on %} Hands-on: Extract workflow
>
> 1. **Clean up** your history: remove any failed (red) jobs from your history by clicking on the {% icon galaxy-cross %} button.
>
>    This will make the creation of the workflow easier.
>
> 2. Click on {% icon galaxy-gear %} (**History options**) at the top of your history panel and select **Extract workflow**.
>
>    ![`Extract Workflow` entry in the history options menu](../../images/history_menu_extract_workflow.png)
>
>    The central panel will show the content of the history in reverse order (oldest on top), and you will be able to choose which steps to include in the workflow.
>
>    ![Selection of steps for `Extract Workflow` from history](../../images/intro_short_workflow_extract.png)
>
> 3. Replace the **Workflow name** to something more descriptive, for example: `QC and filtering`.
>
> 4. **Rename** the workflow input in the box at the top of second column to: `FASTQ reads`
>
> 5. If there are any steps that shouldn't be included in the workflow, you can **uncheck** them in the first column of boxes. In this case, uncheck the second **Filter by quality** tool at the bottom, where we used a too high quality cut-off.
>
> 6. Click on the **Create Workflow** button near the top.
>
>    You will get a message that the workflow was created.
>
{: .hands_on}

In a minute we will see how to find the extracted workflow and how to use it.

## Create a new history

Let's create a new history.

> ### {% icon hands_on %} Hands-on: New history
> 1. Create a new history
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Rename your history, *e.g.* "Next Analysis"
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
{: .hands_on}

This new history does not have any datasets in it yet.

## Look at all your histories

Where is your first history, called "My Analysis"?

> ### {% icon hands_on %} Hands-on: View histories
> 1. Click on the **View all histories** ({% icon galaxy-columns %} icon) at the top right of your history
>
>    ![view all histories](../../images/galaxy_interface_history_switch.png){:width="320px"}
>
>    A new page will appear with all your histories displayed here.
>
> 2. Copy a dataset into your new history
>    1. Click on the FASTQ dataset in "My Analysis" history
>    2. Drag it into the "Next Analysis" history
>
>    This makes a copy of the dataset in the new history (without actually using additional disk space).
>
> 3. Click on the Home icon {% icon galaxy-home %} (or **Analyze Data** on older versions of Galaxy) in the top panel to go back to your analysis window
>
> ![Copy a dataset between histories](../../images/copy-dataset.gif "Copy a dataset between histories by dragging it")
>
{: .hands_on}

Your main Galaxy window will now show "Next Analysis" as the current history, and it will have one dataset in it.

At any time, you can go back into the "View all histories" page and "Switch to" a different history.


## Run workflow in the new history

Now that we have built our workflow, let's use it to re-create our small analysis in a single step. The same workflow could also be used on some new FASTQ data to quickly repeat the same analysis on different inputs.

> ### {% icon hands_on %} Hands-on: Run workflow
>
> 1. Click on **Workflow** in the top menu bar of Galaxy.
>    - Here you have a list of all your workflows.
>    - Your newly created workflow should be listed at the top:
>
>    ![`Your workflows` list](../../images/intro_short_workflow_list.png)
>
>    If you click on a workflow name, you can see all available actions for the workflow, e.g. edit, copy, rename, delete.
>
> 2. Click on the {% icon workflow-run %} (*Run workflow*) button next to your workflow.
>    - The central panel will change to allow you to configure and launch the workflow.
>
>    ![Run workflow form](../../images/intro_short_run_workflow.png)
>
> 3. Check that the *"FASTQ reads"* input is set to the FASTQ dataset we have copied to the new history.
     - In this page we could change any parameter for the tools composing the workflow as we would do when running them one by one.
>
> 4. Click the **Run Workflow** button at the top-right of the screen.
>
> 5. You should see a message that the workflow was successfully invoked. Then jobs will start to run and datasets appear in your "Next Analysis" history, replicating the steps of your previous history.
>
{: .hands_on}


# Conclusion
{:.no_toc}

{% icon trophy %} Well done! You have completed the short introduction to Galaxy, where you named the history, uploaded a file, used a tool, viewed results and run a workflow. Additional tutorials are available for a more in-depth introduction to Galaxy's features.
