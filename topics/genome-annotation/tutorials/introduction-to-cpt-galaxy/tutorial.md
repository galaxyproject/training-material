---
layout: tutorial_hands_on
topic_name: genome-annotation
tutorial_name: introduction-to-cpt-galaxy
---

Galaxy is a platform for doing reproducible bioinformatics research. It provides a friendly interface to vast number of complex command line tools, and it encourages consistent science by using identical software and interfaces across all Galaxy instances. Galaxy aims to make computational biology accessible to research scientists that may not have computer programming or system administration experience. The [Center for Phage Technology (CPT)](https://cpt.tamu.edu/) is dependent on Galaxy for all computer-based analysis. Long-running jobs can be launched while a scientist returns to their lab work; meanwhile, Galaxy keeps track of the progress of the analysis and automatically saves the work done up until that point.

The Galaxy interface consists of three panels; on the left is a column containing all available tools. In the center is where analysis will occur and resulted viewed. On the right is the current history and ability to access other available histories.

> ### {% icon tip %} Note that…
> If the user is lost within the program, clicking **Analyze Data** within the blue bar at the top will return the user to the Galaxy home page. Refreshing the page will also return one to the home page without losing any work done up until that point; Galaxy automatically saves all progress.
{: .tip}

> ### Agenda
>
> In this tutorial, the following Galaxy characteristics will be reviewed:
>
> 1. Tools
> 2. Histories
> 3. Workflows
>
{: .agenda}

# Tools

> ### {% icon comment %} What is a Galaxy Tool?
> A tool is something that generates/transforms data. Within Galaxy, tools are simple interfaces to the complex software behind them. Some examples of a tool would be Gene Caller, which reads the genome and yields a list of possible gene locations, or BLAST, which searches protein sequences in the genome against a database.
{: .comment}

Nearly all Galaxy tools process input files and produce output files. A new user ay be wondering how to upload their data as input. In the top right-hand column of the Tools column is an upload symbol.

<!-- INSERT FRESH IMAGE -->

Clicking on this will bring up an upload menu that will allow the import of data into Galaxy. Files can be dragged and dropped into this box; alternatively, clicking on **Choose local file** in the bottom menu will allow the user to browse files on their local device for upload.

<!-- INSERT FRESH IMAGE -->

> ### {% icon comment %} Advanced Users…
> * **Choose FTP file** allows one to select a file that has been previously uploaded via FTP. *This is required for files >2GB*.
> * **Paste/Fetch data** allows one to paste in a bit of text or a URL. Galaxy will import that into the history panel on the right.
> * [This tutorial](https://galaxyproject.github.io/training-material/topics/introduction/tutorials/galaxy-intro-ngs-data-managment/tutorial.html) offers many different examples of various means of uploading data into Galaxy.
{: .comment}

Once the file has been detected by Galaxy, it will appear in the upload window.

<!-- INSERT FRESH IMAGE -->

> ### {% icon tip %} File Format Issues?
> If Galaxy does not detect the file type properly, the user can set the file type. Although it is a rare occurrence, be sure to double check the file is formatted properly before overriding Galaxy.
{: .tip}

When all of the files desired to be uploaded have been selected, click **Start** in the bottom right of the upload menu.

<!-- INSERT FRESH IMAGE -->

The dataset will indicate to you that it is uploaded in the upload window by yielding a 100% status; at this time, the window can be closed. In the history column, the freshly uploaded dataset will turn yellow…

<!-- INSERT FRESH IMAGE -->

… followed by green when it is ready.

<!-- INSERT FRESH IMAGE -->

There is now data in Galaxy that is ready to be processed by one of many available tools! At the top of the tool panel is a search bar. Alternatively, clicking one of the bold, underlined selections will reveal multiple tools of a certain type which the user can choose from. When selecting the tool, be sure to read the *What it does* text at the bottom of the tool interface page that appears; it will give the user important information regarding running the tool and the consequences of it.

<!-- INSERT FRESH IMAGES -->

Review the options in the tool interface. Keep in mind that many options are set to default values. Uploaded files or other datasets in the current history can be used as inputs.

<!-- INSERT FRESH IMAGES -->

When configuration of the tool is complete, **execute it**, and it will appear as a set of output files in the history on the right.

> ### {% icon tip %} Organizing 
>In order to keep track of work done, tools implemented, and workflows executed, it is advised that the user edits the name of the first dataset that appears

# Histories

