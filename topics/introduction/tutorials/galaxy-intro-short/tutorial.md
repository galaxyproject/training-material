---
layout: tutorial_hands_on

title: "A short and cool introduction to Galaxy"
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
  - "Learn how to share a history"
time_estimation: "40m"
key_points:
  - "The Galaxy interface has an activity bar on the left, a tool (or other activated) panel next to it (if expanded), viewing pane in the middle, and a history of your data analysis on the right."
  - "You can create a new history for each analysis. All your histories are saved."
  - "To get data into Galaxy, you can upload a file by pasting in a web address. There are other ways to get data into Galaxy (not covered in this tutorial): you can upload a file from your computer, and you can import an entire history."
  - "Choose a tool and change any settings for your analysis."
  - "Run the tool. The output files will be saved at the top of your history."
  - "View the output files by clicking on the eye icon."
  - "View all your histories and move files between them. Switch to a different history."
  - "Log out of your Galaxy server. When you log back in (to the same server), your histories will all be there."
subtopic: first
priority: 2
translations:
  - es
contributions:
  authorship:
  - annasyme
  - nsoranzo
  editing:
  - bebatut
  - ahmedhamidawan
  funding:
  - AustralianBioCommons
  - elixir-europe
  - uni-freiburg
---

# Overview


* This is a short introduction to the Galaxy user interface - the web page that you interact with.
* We will cover key tasks in Galaxy: uploading files, using tools, viewing histories, and running workflows.

> <agenda-title></agenda-title>
> 1. TOC
> {:toc}
>
{: .agenda}

## What does Galaxy look like?


> <hands-on-title>Log in to Galaxy</hands-on-title>
> 1. Open your favorite browser (Chrome, Safari or Firefox as your browser, not Internet Explorer!)
> 2. Browse to your Galaxy instance
> 3. Log in or register
>
> ![Screenshot of Galaxy Australia with the register or login button highlighted](../../images/galaxy-login.png)
>
>   > <comment-title>Different Galaxy servers</comment-title>
>   >  This is an image of Galaxy Australia, located at [usegalaxy.org.au](https://usegalaxy.org.au/)
>   >
>   > The particular Galaxy server that you are using may look slightly different and have a different web address.
>   >
>   > You can also find more possible Galaxy servers at the top of this tutorial in **Available on these Galaxies**
>   {: .comment}
{: .hands_on}

The Galaxy homepage is divided into four sections (panels):
* The Activity Bar on the left: _This is where you will navigate to the resources in Galaxy (Tools, Workflows, Histories etc.)_
* Currently active "Activity Panel" on the left: _By default, the {% icon tool %} **Tools** activity will be active and its panel will be expanded_
* Viewing panel in the middle: _The main area for context for your analysis_
* History of analysis and files on the right: _Shows your "current" history; i.e.: Where any new files for your analysis will be stored_

![Screenshot of the Galaxy interface with aforementioned structure](../../images/galaxy_interface.png)

The first time you use Galaxy, there will be no files in your history panel.

# Key Galaxy actions

## Name your current history

Your "History" is in the panel at the right.

> <hands-on-title>Name history</hands-on-title>
> 1. Go to the **History** panel (on the right)
> 2. Click on {% icon galaxy-pencil %} (**Edit**) next to the history name (which by default is "Unnamed history")
>
>    ![Screenshot of the galaxy interface with the history name being edited, it currently reads "Unnamed history", the default value. An input box is below it.]({% link shared/images/rename_history.png %}){:width="250px"}
>
>    > <comment-title></comment-title>
>    >
>    > In some previous versions of Galaxy, you will need to click on the history name to rename it as shown here:
>    > ![Screenshot of the galaxy interface with the history name being edited, it currently reads "Unnamed history", the default value.](../../../../shared/images/rename_history_old.png){:width="320px"}
>    {: .comment}
>
> 3. Type in a new name, for example, "My Analysis"
> 4. Click on **Save**
>
> > <comment-title>Renaming not an option?</comment-title>
> > If renaming does not work, it is possible you aren't logged in, so try logging in to Galaxy first. Anonymous users are only permitted to have one history, and they cannot rename it.
> {: .comment}
>
{: .hands_on}

## Upload a file

The "Activity Bar" can be seen on the left-most part of the interface.

> <hands-on-title>Upload a file from URL</hands-on-title>
> 1. At the top of the **Activity Bar**, click on the {% icon galaxy-upload %} **Upload** activity
>
>    ![upload data button shown in the galaxy interface](../../images/upload-data.png)
>
>    This brings up a box:
>
>    ![the complicated galaxy upload dialog, the 'regular' tab is active with a large textarea to paste subsequent URL](../../images/upload-box.png)
>
> 3. Click **Paste/Fetch data**
> 4. Paste in the address of a file:
>
>    ```
>    https://zenodo.org/record/582600/files/mutant_R1.fastq
>    ```
>
> 5. Click **Start**
> 6. Click **Close**
>
{: .hands_on}

Your uploaded file is now in your current history.
When the file has uploaded to Galaxy, it will turn green.

> <comment-title></comment-title>
> After this you will see your first history item (called a "dataset") in Galaxy's right panel. It will go through
> the gray (preparing/queued) and yellow (running) states to become green (success).
>
{: .comment}

> <details-title>Is this step taking a while?</details-title>
>   Sometimes during courses, data upload gets a little slow. You can also import data through a history link.
>
> 1. Import history from: [example input history](https://humancellatlas.usegalaxy.eu/u/wendi.bacon.training/h/short-introduction-to-galaxy--input)
>
>    {% snippet faqs/galaxy/histories_import.md %}
> 2. **Rename** {% icon galaxy-pencil %} the the history to your name of choice.
>
{: .details}

What is this file?

> <hands-on-title>View the dataset content</hands-on-title>
> 1. Click on the {% icon galaxy-eye %} (eye) icon next to the dataset name, to look at the file content
>
>    ![galaxy history view showing a single dataset mutant_r1.fastq. Display link is being hovered.](../../images/eye-icon.png){:width="320px"}
{: .hands_on}

The contents of the file will be displayed in the central Galaxy panel.

This file contains DNA sequencing reads from a bacteria, in FASTQ format:

   ![preview of a fastq file showing the 4 line structure described in fig caption. 3 reads are shown.](../../images/fastq.png "A FastQ file of course has four lines per record: the record identifier (`@mutant-no_snps.gff-24960/`), the sequence (`AATG…`), the plus character (`+`), and then the quality scores for the sequence (`5??A…`)."){:width="620px"}

## Use a tool

Let's look at the quality of the reads in this file.

> <hands-on-title>Use a tool</hands-on-title>
> 1. Type **FastQC** in the tools panel search box (top)
> 2. Click on the {% tool [FastQC](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.74+galaxy1) %} tool
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

{% snippet faqs/galaxy/tutorial_mode.md %}

## View results

We will now look at the output dataset called *FastQC on data 1: Webpage*.

> <comment-title></comment-title>
> * Note that Galaxy has given this dataset a name according to both the tool name ("FastQC") and the input ("data 1") that it used.
> * The name "data 1" means the dataset number 1 in Galaxy's current history (our FASTQ file).
>
{: .comment}


> <hands-on-title>View results</hands-on-title>
> * Once it's green, click on the {% icon galaxy-eye %} (eye) icon next to the "Webpage" output dataset.
>
>    The information is displayed in the central panel
>
>    ![Graph from fastqc's report. fastqc's images themselves are inaccessible, but this graph shows overall mostly green (good) sequences scores across the length of the read.](../../images/fastqc-out.png){:width="620px"}
{: .hands_on}

This tool has summarised information about all of the reads in our FASTQ file.

> <question-title></question-title>
>
> 1. What was the length of the reads in the input FASTQ file?
> 2. Do these reads have higher quality scores in the centre or at the ends?
>
>   > <solution-title></solution-title>
>   > 1. 150 bp
>   > 2. In the center
>   {: .solution}
{: .question}


## Run another tool

Let's run a tool to filter out lower-quality reads from our FASTQ file.


> <hands-on-title>Run another tool</hands-on-title>
> 1. Type **Filter by quality** in the tools panel search box (top)
> 2. Click on the tool {% tool [Filter by quality](toolshed.g2.bx.psu.edu/repos/devteam/fastq_quality_filter/cshl_fastq_quality_filter/1.0.2+galaxy2) %}
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

> <hands-on-title>Get metadata about a file</hands-on-title>
> 1. Click on the dataset (around where its name is) in the History panel.
>
>    This expands the information about the file. _(By default all datasets are collapsed)_
>
>    ![Diagram of how to locate the information. As above, clicking on the name expands the dataset, and an info section is shown with the filter settings. 1786 (14%) discarded.](../../images/filter-fastq1.png)
>
{: .hands_on}

> <question-title></question-title>
>
> How many reads have been discarded?
>
>   > <solution-title></solution-title>
>   > 1786 low-quality reads were discarded
>   {: .solution}
{: .question}

## Re-run that tool with changed settings

We can now try to filter our input reads to an even higher standard, and see how this changes the resulting output (an exploratory analysis). We will change the filter settings and re-run the tool.

> <hands-on-title>Re-run the tool</hands-on-title>
> 1. Click on the {% icon galaxy-refresh %} icon (**Run Job Again**) for the output dataset of **Filter by quality** {% icon tool %}
>
>    ![A dataset is expanded showing the Run Job Again button highlighted](../../images/rerun.png)
>
>    This brings up the tool interface in the central panel with the parameters set to the values used previously to generate this dataset.
>
> 2. Change the settings to something even stricter
>
>    For example, you might decide you want 80 percent of bases to have a quality of 36 or higher, instead of 35.
>
> 3. Click **Execute**
> 4. View the results: Click on the output dataset name to expand the information
>
>    > <comment-title></comment-title>
>    > Not the {% icon galaxy-eye %} (eye) icon.
>    {: .comment}
>
{: .hands_on}

> <question-title></question-title>
>
> How many reads were discarded under these new filtering conditions?
>
>   > <solution-title></solution-title>
>   > If you selected 80% of bases with 36 as quality cut-off, then 11517 reads (92%) should have been discarded, which indicates that we have gone too far with the filtering in this case.
>   {: .solution}
{: .question}

You can re-run a tool many times with different settings. Each time you re-run the tool, its new output datasets will appear at the top of your current history.

## Share your history

Finally, let's imagine that you had a problem in your analysis and you want to ask for help. The easiest way to ask for help is to share your history. Try and create a link for your history and share it with...yourself!

{% snippet faqs/galaxy/histories_sharing.md %}


## Convert your analysis history into a workflow

When you look carefully at your history, you can see that it contains all the steps of our analysis, from the beginning (at the bottom) to the end (on top). The history in Galaxy records details of every tool you run and preserves all parameter settings applied at each step. But when you need to analyze new data, it would be tedious to do each step one-by-one again. Wouldn't it be nice to just convert this history into a workflow that we will be able to execute again and again?

Galaxy makes this very easy with the `Extract workflow` option. This means any time you want to build a workflow, you can just perform the steps once manually, and then convert it to a workflow, so that next time it will be a lot less work to do the same analysis.

> <hands-on-title>Extract workflow</hands-on-title>
>
> 1. **Clean up** your history: remove any failed (red) jobs from your history by clicking on the {% icon galaxy-delete %} button.
>
>    This will make the creation of the workflow easier.
>
> 2. Click on {% icon galaxy-history-options %} (**History options**) at the top of your history panel and select **Extract workflow**.
>
>    !['Extract Workflow' entry in the history options menu](../../images/history_menu_extract_workflow.png)
>
>    The central panel will show the content of the history in reverse order (oldest on top), and you will be able to choose which steps to include in the workflow.
>
>    ![Selection of steps for Extract Workflow from history. All three of fastqc, filter by quality, and the second filter by quality are selected.](../../images/intro_short_workflow_extract.png)
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

> <hands-on-title>New history</hands-on-title>
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

## Look at multiple histories

Where is your first history, called "My Analysis"?

> <hands-on-title>View histories in <i>History Multiview</i></hands-on-title>
> One of the ways to view multiple histories at once in Galaxy is through the **History Multiview**
> 1. There are multiple ways to get to the multview:
>    - Click on the {% icon galaxy-multihistory %} **History Multiview** activity in the activity bar
>    - Or, click on {% icon galaxy-history-options %} (**History options**) and then click on the {% icon galaxy-columns %} **Show Histories side-by-side** option
>    ![History options menu dropdown showing you have 163 histories, and a show histories side-by-side button.](../../images/galaxy_interface_history_switch.png){:width="320px"}
>
>    A new page will appear with your histories displayed side-by-side here.
>
> 2. Copy a dataset into your new history
>    1. Click on the FASTQ dataset in "My Analysis" history
>    2. Drag it into the "Next Analysis" history
>
>    ![Gif of copying datasets between histories in the side-by-side history view. For now this feature is not keyboard accessible, it is a known issue.](../../images/copy-dataset.gif "Copy a dataset between histories by dragging it")
>
>    This makes a copy of the dataset in the new history (without actually using additional disk space).
>
> 3. Click on the Home icon {% icon galaxy-home %} (or **Analyze Data** on older versions of Galaxy) in the top panel to go back to your analysis window
>
{: .hands_on}

Your main Galaxy window will now show "Next Analysis" as the current history, and it will have one dataset in it.

> <comment-title></comment-title>
>   This is not the only way to view your histories in Galaxy:
>
> 1. An exhaustive list (table) of all your histories is available in the **My Histories** tab in the *Histories List* accessible via clicking the {% icon galaxy-histories-activity %} **Histories** activity
> 2. You can quickly switch to another history by clicking on the {% icon switch-histories %} **Switch to history** button next to {% icon galaxy-history-options %} **History options**:
>
>    {% snippet faqs/galaxy/histories_switch.md %}
>
> _Here is a tutorial on [Understanding Galaxy history system]( {% link topics/galaxy-interface/tutorials/history/tutorial.md %} )_
{: .comment}

## Run workflow in the new history

Now that we have built our workflow, let's use it to re-create our small analysis in a single step. The same workflow could also be used on some new FASTQ data to quickly repeat the same analysis on different inputs.

> <hands-on-title>Run workflow</hands-on-title>
>
> 1. Click on the {% icon galaxy-workflows-activity %} **Workflows** activity in the activity bar.
>
>    Here you have a list of all your workflows _(the **My Workflows** tab is active by default)_.
>
>    Your newly created workflow should be listed at the top:
>
>    ![Workflow list page showing a single workflow named QC and filtering.](../../images/intro_short_workflow_list.png)
>
>    You can see all available actions for the workflow on the workflow card, e.g. edit, copy, rename, share etc. Any other options (e.g.: delete, export etc.) are available by clicking on the {% icon galaxy-dropdown %} **Workflow actions** button on the top right of the card.
>
> 2. Click on the {% icon workflow-run %} (*Run workflow*) button on the bottom right of the workflow card.
>
>    The central panel will change to allow you to configure and launch the workflow.
>
>    ![Run workflow form with a single input: FASTQ reads. mutant_r1.fastq is selected as the input dataset for that parameter.](../../images/intro_short_run_workflow.png)
>
> 3. Check that the *"FASTQ reads"* input is set to the FASTQ dataset we have copied to the new history.
>
>    In this page we could change any parameter for the tools composing the workflow as we would do when running them one by one.
>
> 4. Click the **Run Workflow** button at the top-right of the screen.
>
>    You should see a message that the workflow was successfully invoked. Then jobs will start to run and datasets appear in your "Next Analysis" history, replicating the steps of your previous history.
>
{: .hands_on}


# Conclusion

Well done! You have completed the short introduction to Galaxy, where you named the history, uploaded a file, used a tool, viewed results and run a workflow. Additional tutorials are available for a more in-depth introduction to Galaxy's features.
