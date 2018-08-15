---
layout: tutorial_hands_on
topic_name: introduction
tutorial_name: galaxy-intro-peaks2genes
---

# Introduction
{:.no_toc}

We stumbled upon a paper [Li et al., Cell Stem Cell 2012](https://www.ncbi.nlm.nih.gov/pubmed/22862943) that contains the analysis of possible target genes of an interesting protein in mice. The targets were obtained by ChIP-seq and the raw data is available through [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37268).
The list of genes however is neither in the supplement of the paper nor part of the GEO submission.
The closest thing we could find is a list of the regions where the signal is significantly enriched (so called *peaks*):

1 | 3660676 | 3661050 | 375 | 210 | 62.0876250438913 | -2.00329386666667
1 | 3661326 | 3661500 | 175 | 102 | 28.2950833625942 | -0.695557142857143
1 | 3661976 | 3662325 | 350 | 275 | 48.3062708406486 | -1.29391285714286
1 | 3984926 | 3985075 | 150 | 93 | 34.1879823073944 | -0.816992
1 | 4424801 | 4424900 | 100 | 70 | 26.8023246007435 | -0.66282

**Table 1** Subsample of the available file

The goal of this exercise is to **turn this list of genomic regions into a list of possible target genes**.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Pretreatments

> ### {% icon hands_on %} Hands-on: Open Galaxy
> 1. Browse to a Galaxy instance: the one recommended by your instructor or one in the list **Galaxy instance** on the head of this page
> 2. Log in or register (top panel)
>
>    ![Login or Register on the top panel](../../images/login_register.png)
{: .hands_on}

The Galaxy interface consist of three main parts. The available tools are listed on the left, your analysis history is recorded on the right, and the middle pane will show the tools and datasets.

![Galaxy interface](../../images/galaxy_interface.png "The Galaxy interface")

Let's start with a fresh history.

> ### {% icon hands_on %} Hands-on: Create history
>
> 1. Make sure you have an empty analysis history.
>
>    > ### {% icon tip %} Starting a new history
>    >
>    > * Click the {% icon galaxy-gear %} (gear) icon (**History options**)  at the top of the history panel
>    > * Select the option **Create New** from the menu
>    {: .tip}
>
> 2. **Rename your history** to make it easy to recognize
>
>    > ### {% icon tip %} Rename a history
>    >
>    > * Click on the title of the history (by default the title is `Unnamed history`)
>    >
>    >   ![Renaming history](../../../../shared/images/rename_history.png)
>    >
>    > * Type `Galaxy Introduction` as the name
>    >
>    {: .tip}
>
{: .hands_on}

## Data upload

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Download the list of peak regions (the file [`GSE37268_mof3.out.hpeak.txt.gz`](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE37268&format=file&file=GSE37268%5Fmof3%2Eout%2Ehpeak%2Etxt%2Egz)) from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37268) to your computer
> 2. Click on the upload button in the upper left ot the interface
>
>    ![Upload icon](../../../galaxy-data-manipulation/images/upload_button.png)
>
> 3. Press **Choose local file** and search for your file on your computer
> 4. Select `interval` as **Type**
> 5. Press **Start**
> 6. Press **Close**
> 7. Wait for the upload to finish
>
{: .hands_on}

Galaxy will automatically unpack the file.

> ### {% icon comment %} Comment
> After this you will see your first history item in Galaxy’s right pane. It will go through
> the gray (preparing/queued) and yellow (running) states to become green (success):
>
> ![History section](../../images/intro_01.png)
{: .comment}


> ### {% icon tip %} Tip: Importing data via links
>
> You could also import this file via its link:
>
> * Copy the link location
> * Open the Galaxy Upload Manager
> * Select **Paste/Fetch Data**
> * Paste the link into the text field
> * Select `interval` as **Type**
> * Press **Start**
{: .tip}

> ### {% icon hands_on %} Hands-on: Inspect and edit attribute of a file
>
> 1. Click on the file in the history panel
>
>    Some meta-information (e.g. format, reference database) about the file and the header of the file are then displayed:
>
>    ![Expanded file in history](../../images/input_expanded_file.png)
>
> 2. Click on the {% icon galaxy-eye %} (eye) icon (**View data**) in your dataset in the history
>
>    The content of the file is displayed in the middle panel
>
> 3. Click on the {% icon galaxy-pencil %} (pencil) icon (**Edit attributes**) in your dataset in the history
>
>    A form to edit dataset attributes is displayed in the middle panel
>
> 4. Search for `mm9` in **Database/Build** attribute and select `Mouse July 2007 (NCBI37/mm9)`
> 5. Click on **Save** on the top
>
{: .hands_on}

In order to find the related genes to these peak regions,
we also need a list of genes in mice, which we can obtain from UCSC.

> ### {% icon hands_on %} Hands-on: Data upload from UCSC
>
> 1. Search for `UCSC Main` in the tool search bar (top left)
>
>     ![UCSC Main tool in tools section](../../images/101_01.png)
>
> 2. Click on `UCSC Main` {% icon tool %}
>
>     You will be taken to the **UCSC table browser**, which looks something like this:
>
>     ![UCSC table browser interface](../../images/intro_02.png)
>
> 3. Set the following options:
>     - *"clade"*: `Mammal`
>     - *"genome"*: `Mouse`
>     - *"assembly"*: `July 2007 (NCBI37/mm9)`
>     - *"group"*: `Genes and Gene Predictions`
>     - *"track"*: `RefSeq Genes`
>     - *"table"*: `refGene`
>     - *"region"*: `genome`
>     - *"output format"*: `BED - browser extensible data`
>     - *"Send output to"*: `Galaxy` (only)
>
> 4. Click on the **get output** button
>
>    You will see the next screen:
>
>    ![Output settings](../../images/intro_03.png)
>
> 5. Make sure that *"Create one BED record per"* is set to `Whole Gene`
> 6. Click on the **Send Query to Galaxy** button
> 7. Wait for the upload to finish
> 8. Rename our dataset ({% icon galaxy-pencil %} (pencil) icon) to something more recognizable (`Genes`)
>
{: .hands_on}

> ### {% icon comment %} BED file format
> The **BED - Browser Extensible Data** format provides a flexible way to encode gene regions. BED lines have three required fields:
> - chromosome ID
> - start position (0-based)
> - end position (end-exclusive)
>
> There can be up to and nine additional optional fields, but the number of fields per line must be consistent throughout any single set of data.
>
> You can find more information about it at [UCSC](https://genome.ucsc.edu/FAQ/FAQformat#format1) including a description of the optional fields.
{: .comment}

Now we collected all the data we need to start our analysis.

# Part 1: Naive approach

## File preparation

Let's have a look at our files to see what we actually have here.

> ### {% icon hands_on %} Hands-on: View file content
>
> 1. Click on the {% icon galaxy-eye %} (eye) icon (**View data**) of the peak file to view the content of it
>
>    It should look like this:
>
>    ![Contents of the peak file](../../images/intro_04.png)
>
> 2. View the content of the regions of the genes from UCSC
>
>    ![Contents of UCSC file](../../images/intro_05.png)
>
{: .hands_on}

> ### {% icon question %} Questions
>
> While the file from UCSC has labels for the columns, the peak file does not. Can you guess what the columns stand for?
>
{: .question}

This peak file is not in any standard format and just by looking at it, we cannot find out what the numbers in the different columns mean. In the paper the authors mention that they used the peak caller [HPeak](https://www.ncbi.nlm.nih.gov/pubmed/20598134).

By looking at the HPeak manual we can find out that the columns contain the following information:

 - chromosome name by number
 - start coordinate
 - end coordinate
 - length
 - location within the peak that has the highest hypothetical DNA fragment coverage (summit)
 - not relevant
 - not relevant

In order to compare the two files, we have to make sure that the chromosome names follow the same format.
As we directly see, the peak file lacks `chr` before any chromosome number. But what happens with chromosome 20 and 21? Will it be X and Y instead? Let's check:

> ### {% icon hands_on %} Hands-on: View end of file
>
> 1. Search for **Select last** {% icon tool %} tool and run **Select last lines from a dataset (tail)** with the following settings:
>     - *"Text file"*: our peak file `GSE37268_mof3.out.hpeak.txt`
>     - *"Operation"*: `Keep last lines`
>     - *"Number of lines"*: Choose a value, e.g. `100`
> 2. Click **Execute**
> 3. Wait for the job to finish
> 4. Inspect the file through the {% icon galaxy-eye %} (eye) icon (**View data**)
>
>    > ### {% icon question %} Questions
>    >
>    > 1. How are the chromosomes named?
>    > 2. How are the chromosomes X and Y named?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. The chromosomes are just given by their number. In the gene file from UCSC, they started with `chr`
>    > > 2. The chromosomes X and Y are names 20 and 21
>    > {: .solution }
>    {: .question}
{: .hands_on}

In order to convert the chromosome names we have therefore two things to do:

1. add `chr`
2. change 20 and 21 to X and Y

> ### {% icon hands_on %} Hands-on: Adjust chromosome names
>
> 1. **Replace Text** {% icon tool %}: Run **Replace Text in a specific column** with the following settings:
>     - *"File to process"*: our peak file `GSE37268_mof3.out.hpeak.txt`
>     - *"in column"*: `Column:1`
>     - *"Find pattern"*: `[0-9]+`
>
>         This will look for numerical digits
>
>     - *"Replace with"*: `chr&`
>
>         `&` is a placeholder for the find result of the pattern search
>
> 2. **Replace Text** {% icon tool %}: Let's rerun the tool with
>    - *"File to process"*: the output from the last run, e.g. something like `Replace Text on data ...`
>    - *"in column"*: `Column:1`
>    - *"Find pattern"*: `chr20`
>    - *"Replace with"*: `chrX`
>
>    > ### {% icon tip %} Tip: Rerunning a tool
>    >
>    > * Expand the dataset information
>    > * Press the {% icon galaxy-refresh %} icon (**Run this job again**)
>    {: .tip}
>
> 3. **Replace Text** {% icon tool %}: Rerun this tool to do the same for chromosome Y
>    - *"File to process"*: the output from the **last** run
>    - *"in column"*: `Column:1`
>    - *"Find pattern"*: `chr21`
>    - *"Replace with"*: `chrY`
>
> 4. Inspect the latest file through the {% icon galaxy-eye %} (eye) icon
>
>    Have we been successful?
>
>    We have quite some files now and should take care that we don't loose track.
>
> 5. Rename the file to something more handy, e.g. `Peak regions`
{: .hands_on}

## Analysis

Our goal is still to compare the 2 region files (the genes file and the peak file from the publication)
to know which peaks are related to which genes. If you really only want to know which peaks are located **inside** genes you
can skip the next step. Otherwise, it might be reasonable to include the promoter region into the comparison, e.g. because
you want to include Transcriptions factors in ChIP-seq experiments.

> ### {% icon hands_on %} Hands-on: Add promoter region to gene records
>
> 1. **Get Flanks** {% icon tool %}: Run **Get flanks returns flanking region/s for every gene** with the following settings:
>     - *"Select data"*: file from UCSC
>     - *"Region"*: `Around Start`
>     - *"Location of the flanking region/s"*: `Upstream`
>     - *"Offset"*: `10000`
>     - *"Length of the flanking region(s)"*: `12000`
>
>     This tool returns flanking regions for every gene
>
> 2. Compare the rows of the resulting BED file with the input to find out how the start and end positions changed
>
>    > ### {% icon tip %} Tip: Inspecting several files using the scratchbook
>    >
>    > * Click **Enable/Disable Scratchbook** on the top panel
>    >
>    >    ![Enable/Disable Scratchbook](../../images/intro_scratchbook_enable.png)
>    >
>    > * Click on the {% icon galaxy-eye %} (eye) icon of the files to inspect
>    > * Click on **Show/Hide Scratchbook**
>    >
>    >    ![Show/Hide Scratchbook](../../images/intro_scratchbook_show_hide.png)
>    {: .tip}
>
> 3. Rename your dataset to reflect your findings
{: .hands_on}

You might have noticed that the UCSC file is in `BED` format and has a database associated to it. That's what we want for our peak file as well

> ### {% icon hands_on %} Hands-on: Change format and database
>
> 1. Click on the {% icon galaxy-pencil %} (pencil) icon in the history entry of our peak region file
> 2. Switch to the **Convert** tab
> 3. Select `Convert Genomic Intervals to BED`
> 4. Press **Convert datatype**
> 5. Check that the "Database/Build" is `mm9` (the database build for mice used in the paper)
{: .hands_on}

It's time to find the overlapping intervals (finally!). To do that, we want to extract the genes which overlap/intersect with our peaks.

> ### {% icon hands_on %} Hands-on: Find Overlaps
>
> 1. **Intersect** {% icon tool %}: Run **Intersect the intervals of two datasets** with the following settings:
>     - *"Return"*: `Overlapping Intervals`
>     - *"of"*: the UCSC file with promoter regions
>     - *"that intersect"*: our converted peak region file
>     - *"for at least"*: `1`
>
>    > ### {% icon comment %} Comments
>    > The order of the inputs is important! We want to end up with a list of genes, so the corresponding dataset needs to be the first input.
>    {: .comment}
{: .hands_on}

We now have the list of genes (column 4) overlapping with the peak regions.
To get a better overview of the genes we obtained, we want to look at their distribution across the different chromosomes.
We will regroup the table by chromosome and count the number of genes with peaks on each chromosome

> ### {% icon hands_on %} Hands-on: Count genes on different chromosomes
>
> 1. **Group** {% icon tool %}: Run **Group data by a column and perform aggregate operation on other columns** with the following settings:
>     - *"Select data"* to the result of the intersection
>     - *"Group by column"*:`Column 1`
>     - Press **Insert Operation** and choose:
>         - *"Type"*: `Count`
>         - *"On column"*: `Column 1`
>         - *"Round result to nearest integer?"*: `No`
>
>    > ### {% icon question %} Questions
>    >
>    > Which chromosome contained the highest number of target genes?
>    >
>    > > ### {% icon solution %} Solution
>    > > The result varies with different settings. If you followed step by step, it should be chromosome 7 with 1675 genes.
>    > {: .solution }
>    {: .question}
>
{: .hands_on}

## Visualization

Since we have some nice data, let's draw a barchart out of it!

> ### {% icon hands_on %} Hands-on: Draw barchart
>
> 1. Click on {% icon galaxy-barchart %} (visualize) icon at the latest history item
> 2. Select `Bar diagram`
> 3. Choose a title at **Provide a title**, e.g. `Gene counts per chromosome`
> 4. Switch to the **Select data** tab and play around with the settings
>
{: .hands_on}

## Extracting workflow

When you look carefully at your history, you can see that it contains all steps of our analysis, from the beginning to the end. By building this history we have actually built a complete record of our analysis with Galaxy preserving all parameter settings applied at every step.
Wouldn't it be nice to just convert this history into a workflow that we'll be able to execute again and again?

Galaxy makes this very simple with the `Extract workflow` option. This means that any time you want to build a workflow, you can just perform it manually once, and then convert it to a workflow, so that next time it will be a lot less work to do the same analysis. It also allows you to easily share or publish your analysis.

> ### {% icon hands_on %} Hands-on: Extract workflow
>
> 1. **Clean up** your history
>
>    If you had any failed jobs (red), please remove those datasets from your history by clicking on the {% icon galaxy-cross %} (cross) button. This will make the creation of a workflow easier.
>
> 2. Go to the history **History options** ({% icon galaxy-gear %} (gear) icon)
> 3. Select the `Extract Workflow` option.
>
>    ![Extracting workflow in history menu](../../images/history_menu_extract_workflow.png)
>
>    The center panel will change and you will be able to choose which steps to include/exclude and how to name the newly created workflow.
>
> 4. **Uncheck** any steps that shouldn't be included in the workflow
>
>    Since we did some steps which where specific to our custom peak file, we might want to exclude:
>    - **Select last** {% icon tool %}
>    - all **Replace Text** {% icon tool %} steps
>    - **Convert Genomic Intervals to strict BED** {% icon tool %}
>    - **Get flanks** {% icon tool %}
>
> 5. Rename the workflow to something descriptive, e.g. `From peaks to genes`
>
> 6. Click on the **Create Workflow** button near the top.
>
>    You will get a message that the workflow was created. But where did it go?
>
> 7. Click on **Workflow** in the top menu of Galaxy
>
>    Here you have a list of all your workflows
>
> 8. Select the newly generated workflow and click on **Edit**
>
>    You should see something similar to this:
>
>    ![Editing workflow interface](../../images/intro_06.png)
>
>    > ### {% icon comment %} The workflow editor
>    > We can examine the workflow in Galaxy's workflow editor. Here you can view/change the parameter settings of each step, add and remove tools, and connect an output from one tool to the input of another, all in an easy and graphical manner. You can also use this editor to build workflows from scratch.
>    {: .comment}
>
>     Although we have our two inputs in the workflow they are missing their connection to the first tool (**Intersect** {% icon tool %}), because we didn't carry over some of the intermediate steps.
>
> 9. Connect each input dataset to the **Intersect** {% icon tool %} tool by dragging the arrow pointing outwards on the right of its box (which denotes an output) to an arrow on the left of the **Intersect** box pointing inwards (which denotes an input)
> 10. Rename the input datasets to `Reference regions` and `Peak regions`
> 11. Click on the {% icon galaxy-gear %} (gear) icon at the top right
> 12. Press **Auto Re-layout** to clean up our view
>     ![Auto re-layouting](../../images/intro_07.png)
> 13. Click on the {% icon galaxy-gear %} (gear) icon at the top right
> 14. Press **Save** to save your changes
>
>    > ### {% icon tip %} Tip: Hiding intermediate steps
>    > When a workflow is executed, the user is usually primarily interested in the final product and not in all intermediate steps. By default all the outputs of a workflow will be shown, but we can explicitly tell Galaxy which output to show and which to hide for a given workflow. This behaviour is controlled by the little asterisk next to every output dataset:
>    >
>    > ![Workflow editor mark output](../../../../shared/images/workflow_editor_mark_output.png)
>    >
>    > If you click on this asterisk for any of the output datasets, then *only* files with an asterisk will be shown, and all outputs without an asterisk will be hidden (Note that clicking *all* outputs has the same effect as clicking *none* of the outputs, in both cases all the datasets will be shown).
>    {: .tip}
>
{: .hands_on}

Now it's time to reuse our workflow for a more sophisticated approach.

# Part 2: More sophisticated approach

In part 1 we used an overlap definition of 1 bp (default setting). In order to get a more meaningful definition, we now want to use the information of the position of the peak summit and check for overlap of the summits with genes.

## Preparation

Create a new history and name it. If you forgot how to do that, you can have a look at the beginning of this tutorial.
The history is now empty, but we need our peak file again. Before we upload it twice, we can copy it from our former history:

> ### {% icon hands_on %} Hands-on: Copy history items
>
> 1. Click on the **View all histories** ({% icon galaxy-columns %} icon) at the top right of your history
>
>       You should see both of your histories side-by-side now
>
> 2. Use drag-and-drop with your mouse to copy the edited peak file (after the replace steps) but still in interval format, which contains the summit information, to your new history.
> 3. Click on **Analyze Data** in the top panel to go back to your analysis window
>
{: .hands_on}

## Create peak summit file

We need to generate a new BED file from the original peak file that contains the positions of the peak summits. The start of the summit is the start of the peak (column 2) plus the location within the peak that has the highest hypothetical DNA fragment coverage (column 5). As the end we simply define `start + 1`.

> ### {% icon hands_on %} Hands-on: Create peak summit file
>
> 1. **Compute** {% icon tool %}: Run **Compute an expression on every row** with the following settings:
>   - *"Add expression"*: `c2+c5`
>   - *"as a new column to"*: our peak file
>   - *"Round result?"*: `YES`
> 2. **Compute an expression on every row** {% icon tool %}: rerun this tool on the last result with:
>   - *"Add expression"*: `c8+1`
>   - *"as a new column to"*: the result from step 1
>   - *"Round result?"*: `YES`
>
{: .hands_on}

Now we cut out just the chromosome plus the start and end of the summit:

> ### {% icon hands_on %} Hands-on: Cut out columns
> 1. **Cut** {% icon tool %}: Run **Cut columns from a table** with the following settings:
>   - *"Cut columns"*: `c1,c8,c9`
>   - *"Delimited by Tab"*: `Tab`
>   - *"From"*: our latest history item
>
>    The output from **Cut** will be in `tabular` format.
>
> 2. Change the format to `interval` since that's what the tool **Intersect** expects.
{: .hands_on}

## Get gene names

The RefSeq genes we downloaded from UCSC did only contain the RefSeq identifiers, but not the gene names. To get a list of gene names in the end, we use another BED file from the Data Libraries.

> ### {% icon comment %} Comments
> There are several ways to get the gene names in, if you need to do it yourself. One way is to retrieve a mapping through Biomart and then join the two files (**Join two Datasets side by side on a specified field** {% icon tool %}). Another is to get the full RefSeq table from UCSC and manually convert it to BED format.
{: .comment}

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Import from [Zenodo](https://zenodo.org/record/1025586) or from the data library the file `mm9.RefSeq_genes_from_UCSC.bed`
>
>      ```
>      https://zenodo.org/record/1025586/files/mm9.RefSeq_genes_from_UCSC.bed
>      ```
>
>    > ### {% icon tip %} Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * Press **Start**
>    {: .tip}
>
>    > ### {% icon tip %} Tip: Importing data from a data library
>    >
>    > * Go into **Shared data** (top panel) then **Data libraries**
>    > * Click on "Training data" and then "Introduction - From peaks to genes"
>    > * Select interesting file
>    > * Click on "Import selected datasets into history"
>    > * Import in a new history
>    {: .tip}
>
>    As default, Galaxy takes the link as name, so rename them.
>
> 2. Inspect the file content to check if it contains gene names
>
{: .hands_on}

## Repeat workflow

It's time to reuse the workflow we created earlier.

> ### {% icon hands_on %} Hands-on: Run a workflow
> 1. Open the workflow menu (top menu bar)
> 2. Find the workflow you made in the previous section, and select the option **Run**
> 3. Choose as inputs our imported gene BED file and the result of the **Cut** tool
> 4. Click **Run workflow**
>
>    The outputs should appear in the history but it might take some time until they are finished.
>
{: .hands_on}

We used our workflow to rerun our analysis with the peak summits. The **Group** tool again produced a list containing the amount of genes found in each chromosome.
But woudln't it be more interesting to know about the amount of peaks in each unique gene? Let's rerun the workflow with different settings!

> ### {% icon hands_on %} Hands-on: Run a workflow with changed settings
> 1. Open the workflow menu (top menu bar)
> 2. Find the workflow you made in the previous section, and select the option **Run**
> 3. Choose as inputs our imported gene BED file and the result of the **Cut** tool
> 4. Click on the title of the Group tool to expand the options.
> 5. Change the following settings by clicking at the {% icon galaxy-pencil %} (pencil) icon on the left:
>     - **Group by column**: `7`
>     - **Operation -> On column**: `7`
> 6. Click **Run workflow**
{: .hands_on}

Congratulations! You should have a file with all the unique gene names and a count on how many peaks they contained.

> ### {% icon question %} Questions
>
> The list of unique genes is not sorted. Try to sort it on your own!
>
> > ### {% icon solution %} Solution
> > You can use the tool "Sort data in ascending or descending order" on column 2 and a numerical sort.
> {: .solution }
{: .question}


# Share your work

One of the most important features of Galaxy comes at the end of an analysis. When you have published striking findings, it is important that other researchers are able to reproduce your in-silico experiment. Galaxy enables users to easily share their workflows and histories with others.

To share a history, click on the {% icon galaxy-gear %} (gear) symbol in the history pane and select `Share or Publish`. On this page you can do 3 things:

1. **Make accessible via Link**

    This generates a link that you can give out to others. Anybody with this link will be able to view your history.

2. **Publish History**

    This will not only create a link, but will also publish your history. This means your history will be listed under `Shared Data → Published Histories` in the top menu.

3. **Share with Individual Users**

    This will share the history only with specific users on the Galaxy instance.


> ### {% icon hands_on %} Hands-on: Share history and workflow
>
> 1. Share one of your histories with your neighbour
> 2. See if you can do the same with your workflow!
> 3. Find the history and/or workflow shared by your neighbour
>
>    Histories shared with specific users can be accessed by those users in their history menu ({% icon galaxy-gear %} (gear) icon) under `Histories shared with me`.
>
{: .hands_on}

# Conclusion
{:.no_toc}

{% icon trophy %} You have just performed your first analysis in Galaxy. You also created a workflow from your analysis so you can easily repeat the exact same analysis on other datasets. Additionally you shared your results and methods with others.
