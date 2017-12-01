---
layout: tutorial_hands_on
topic_name: introduction
tutorial_name: galaxy-intro-history
---

When data is uploaded from your computer or analysis is done on existing data using Galaxy, each output from those steps
generates a dataset. These datasets (and the output datasets from later analysis on them) are stored by Galaxy in
**Histories**.

## The Current History

All users have one 'current' history, which can be thought of as **a workspace or a current working directory** in
bioinformatics terms. You current history displayed in the right hand side of the main 'Analyze Data' Galaxy page in
what is called the history panel.

The history panel displays output datasets in the order in which they were created with the oldest/first shown on the
bottom. As new analyses are done and new output datasets are generated, the newest datasets are added to the top of the
the history panel. In this way, the history panel displays the history of your **analysis over time**.

**Users that have registered an account and logged in can have many histories** and the history panel allows switching
between them and creating new ones. This can be useful to organize different analyses.

**Anonymous users** (if your Galaxy allows them) are users that have not registered an account. Anonymous users are
only allowed one history. On our main, public Galaxy server, users are encouraged to register and log in with the
benefit that they can work on many histories and switch between them.

The histories of anonymous users are only associated through your browser's session. **If you close the browser or
clear you sessions - that history will be lost!** We can not recover it for you if it is.
{: .alert .alert-warning}

### Current history controls

![Current history buttons](../../images/current-history-buttons.png)

Above the current history panel are three buttons: the refresh, history options, and 'view all histories' button.

The refresh button will entirely reload the history being viewed. This can be helpful if you believe the history
interface needs to be updated or isn't updating properly.

The history options button opens the history options menu which allows you to perform history-related tasks.

The 'view all histories' button sends you to the interface for managing multiple histories

## History Information

Histories also store information apart from the datasets they contain. They can be named/re-named, tagged, and
annotated.

### Renaming a history

All histories begin with the name 'Unnamed history'. Non-anonymous users can rename the history as they see fit:

1. Click the existing name. A text input field will appear with the current name.
2. Entered a new name or edit the existing the way you'd like.
3. Press 'Enter' to save the new name. The input field will disappear and the new name display.
4. To cancel renaming, press 'Esc' or click outside the input field.

![Renaiming history](../../images/renaming.png)

### Tagging a history

Tags are short pieces of text used to describe the thing they're attached to and many things in Galaxy can be tagged.
Each item can have many tags and you can add new tags or remove them at any time. Tags can be another useful way to
organize and search your data. For instance, you might tag a history with the type of analysis you did in it: 'assembly'
or 'variants'. Or you may tag them according to data sources or some other metadata: 'long-term-care-facility' or
'yellowstone park:2014'.

To tag a history:

1. Click the tag button at the top of the history panel. An input field showing existing tags (if any) will appear.
2. Begin typing your new tag in the field. Any tags that you've used previously will show below your partial entry -
  allowing you to use this 'autocomplete' data to re-use your previous tags without typing them in full.
3. Press enter or select one of the previous tags with your arrow keys or mouse.
4. To remove an existing tag, click the small 'X' on the tag or use the backspace key while in the input field.

![Using tags](../../images/tags.png)

### Annotating a history

Sometimes tags and names are not enough to describe the work done within a history. Galaxy allows you to create history
annotations: longer text entries that allow more formatting options. Newlines and whitespace are preserved. Later, if
you publish or share the history, the annotation will be displayed automatically - allowing you to share additional
notes about the analysis.

To annotate a history:

1. Click the annotation button at the top of the history panel. A larger text section will appear displaying any
  existing annotation (or, if there's none, italic text saying you can click on the control to create an annotation).
2. Click the annotation section. The a larger input field will appear.
3. Add any annotations you desire. 'Return' will create a line break and white space is preserved. (Tabs cannot be
  entered since the 'Tab' button is used to switch between controls on the page - tabs can be pasted in however).
4. To save the annotation, click the 'Done' button.

![Annotations to history](../../images/annotations.png)

### History size

As datasets are added to a history, Galaxy will store them in files on its file system. The total size of these files
for all the datasets in a history is displayed underneath the history name. For example, if a history has 200 megabytes
of dataset data on Galaxy's filesystem, '200 MB' will be displayed underneath the name.

If your Galaxy server uses quotas, the total combined size of all your histories will be compared to your quota and
if you're using more than the quota allows Galaxy will prevent you from running any new jobs until you've deleted some
datasets and brought that total below the quota.


## History Panel Datasets

Datasets in the history panel show the state of the job that has generated or will generate the data.

There are several different 'states' a dataset can be in:

1. When you first upload a file or run a tool, the dataset will be in the **queued** state. This indicates that the
  job that will create this dataset has not yet started and is in line to begin.
1. When the job starts the dataset will be in the **running** state. The job that created these datasets is now
  running on Galaxy's cluster.
1. When the job has completed successfully, the datasets will be in the **ok** state.
1. If there's been an error while running the tool, the datasets will be in the **error** state.
1. If a previously running or queued job has been paused by Galaxy, the dataset will be in the **paused** state.
  You can re-start/resume paused jobs using the options menu above the history panel and selecting 'Resume Paused Jobs'.

![States in history](../../images/states.png)

Datasets in the panel are initially shown in a 'summary' view, that only displays:

1. a **number** indicating in what order (or what step) this dataset was created
1. the dataset **name**
1. a **view** button: click this to view the dataset contents in raw format in the browser
1. an **edit** button: click this to edit dataset properties
1. a **delete** button: click this to delete the dataset from the history (*don't worry*, you can undo this action)


**Note:** some of the buttons above may be disabled if the dataset is in a state that shouldn't allow doing the
action. For example, the 'edit' button is disabled for datasets that are still queued or running.
{: .alert .alert-warning}

![Datasets as summary in history](../../images/summary.png)

You can click the dataset name and the view will expand to show more details:

1. a short description of the data
1. the file **format**
1. the reference sequence (or **database**) for the data
1. (optionally) some information/output from the job that produced this dataset
1. a row of buttons that allow further actions on the dataset
1. a **peek** of the data: a couple of rows of data with the column headers (if available)

![Details](../../images/details.png)


**Note:** many of these details are only displayed if the dataset has finished running, is in the 'ok' state, and
is not deleted. Otherwise, you may only see a shorter message describing the dataset's state (e.g. 'this dataset
is waiting to run')
{: .alert .alert-warning}

## Managing Datasets Individually

### Hiding and unhiding datasets

Some procedures in Galaxy such as workflows will often **hide** history datasets in order to simplify the history
and hide intermediate steps of an automated analysis. These hidden datasets won't normally appear in the history panel
but their still mentioned in the history subtitle (the smaller, grey text that appears below the history name). If you
history has hidden datasets, the number will appear there (e.g. '3 hidden') as a clickable link. If you click this link,
the hidden datasets are shown. Each hidden dataset has a link in the top of the summary view that allows you to unhide
it. You can click that link again (which will now be 'hide hidden') to make them not shown again.

![Hiding dataset](../../images/hide.png)

### Deleting and undeleting datasets

You can **delete** any dataset in your history by clicking the delete button. This does not immediately remove the
dataset's data from Galaxy and **it is reversible**. When you delete a dataset from the history, it will be removed
from the panel but (like hidden datasets) the total number of deleted datasets is shown in the history subtitle as a
link. Clicking this link (e.g. '3 deleted') will make the deleted datasets visible and each deleted dataset will have a
link for manually undeleting it above its title. You can click that link again (which will now be 'hide deleted') to
make them not shown again.

![Deleting dataset](../../images/delete.png)

### Purging datasets and removing them permanently from Galaxy

If you are showing deleted datasets and *your Galaxy allows users to purge datasets*, you will see an additional link
in the top of each deleted dataset titled **'Permanently remove it from disk**'. Clicking this will remove the file
that contains that dataset's data and will decrease the disk space used by the history. **This action is not reversible
and cannot be undone**.

If your Galaxy doesn't allow users to purge their datasets, you will not see that link.

### Admins may purge your deleted datasets

Depending on the policy of your Galaxy server, administrators will often run scripts that search for and purge the
datasets you've marked as deleted. Often deleted datasets and histories are purged based on the age of the deletion
(e.g. datasets that have been marked as deleted for 90 days or more). Check with the administrators of your instance to
find out the policy used.

## Managing Multiple Datasets Easily

### Muti-selection

You can also hide, delete, and purge multiple datasets at once by **multi-selecting datasets**:

1. Click the multiselect button containing the checkbox to the right of the history size.
1. Checkboxes will appear inside each dataset in the history.
1. Scroll and click the checkboxes next to the datasets you want to manage. You can also 'shift+click' to select a
  range of datasets.
1. Click the 'All' button to select all datasets. Click the 'None' button to clear/remove all selections.
1. Click the 'For all selected...' to choose the action. The action will be performed on all selected datasets. (If
  an action doesn't apply to a selected dataset - like deleting a deleted dataset - nothing will happen.)
1. You can click the multiselect button again to hide the checkboxes again.

![Multiselecting datasets](../../images/multiselect.png)

### Searching for datasets

You can filter what datasets are shown and search for datasets using the search bar at the top of the panel. Enter
any text that a dataset you'd be looking for would contain, including:

* the name or part of the name
* any text (or partial text) from the description and info fields
* the file format or reference database
* any text or partial text from the annotation or tags of a dataset

For example:

* To find all vcf files you might enter: `vcf` alone.
* To find all files whose names contain data 1, you can enter: `data 1`
* To search for a VCF file named 'VCF filter on data 1' and tagged with 'experiment 1', you could enter:
  `vcf filter on data 1 experiment 1`


**Note:** searches are case-insensitive. For example, `VCF` and `vcf` are equivalent.
{: .alert .alert-warning}

![Searches for dataset](../../images/basic-search.png)

### Clearing a search

You can clear a search and show all visible datasets by clicking the round 'X' button in the right of the search bar
or - while entering text in the search bar - hitting the escape key ('esc').

### Advanced searching

You can also specify what dataset properties you're searching using keyword search. These are the property names
followed by '=' then the value. When using these only the property named is searched for that value:

* *name*: for example, `name=my-dataset` would show all datasets whose names contain 'my-dataset' - *but* not
  search other properties such as annotation, etc.
* *format*: to search for a particular format: `format=vcf`
* *database*: to search for all datasets with a particular reference: `database=hg19`
* *annotation*, *description*, *info*: the annotation, summary description, or job info respectively
* *tag*: to show only datasets with that tag (or partial tag): `tag=experiment1`. Note: you can search for
  datasets that have multiple tags by re-using the tag keyword search: `tag=experiment1 tag=to_publish`


**Note** that keyword searches can be combined: `database=mm10 annotation=successful`<br>
You can enclose text and include spaces using double quotes: `name="My Dataset" annotation="First run"`.<br>
 All keyword searches are AND'd - all searches should be true in order for the dataset to be shown.
{: .alert .alert-warning}

If you find normal searching is showing too many datasets, and not what you're looking for, try the advanced keyword
search.

![Advanced search](../../images/adv-search.png)

### Search and multiselect

It's often useful to combine search and multiselect. Multiselections will persist over searches and the All/None buttons
will only apply to those datasets that are currently shown with the given search.

So, for example, to select and manipulate two different sets of datasets: a) fastqsanger files tagged with 'low quality'
and b) hg19 reference BAM files whose names contain "Output":

1. In the search bar, enter: `format=fastqsanger tag="low quality"` and hit enter.
1. Click the multiselect button to show the checkboxes.
1. Click the 'All' button to select all the fastqsanger files
1. In the search bar again, enter: `database=hg19 format=BAM name=output` and hit enter.
1. Click the 'All' button again to select all the BAM files.
1. Clear the search using the clear button. All datasets are shown and both fastqsanger and BAM files are still
  selected.
1. You can now perform some action on those two sets of datasets.

### Undeleting ...deleted histories

If you have not purged a history and only deleted it, it is possible to 'undelete' it and reverse or undo the deletion.
Since one of the purposes of deleting hisories is to remove them from view, we'll use the interface to specifically
search for deleted histories and then to undelete the one we're interested in.

There are two ways to do this currently: via the multi-history panel and through the saved histories page. The multi-
history method is presented here:

Click the multi-history icon at the top right of the 'Analyze Data' (home) page. Note: you must be logged in to
see the icon and use the multi-history page. You should see all the (non-deleted) histories that you've created.
![Multi-history button](../../images/undelete.multihistory-button.png)

Click the '...' icon button in the grey header at the top of the page. You should see a dialog that presents some options for viewing the histories. Click the 'include deleted histories' option.
![Multihistory options](../../images/undelete.multihistory-options.png)

The page should reload and now both non-deleted and deleted histories will be displayed. Deleted histories will
have a small message under their titles stating 'This history has been deleted'.

![Deleted and non-deleted histories](../../images/undelete.thishasbeendeleted.png)
Now click the small button with the down arrow just above the deleted history you want to undelete. Then click
the 'Undelete' option there. Your history should now be undeleted.
![Undelete button](../../images/undelete.undelete-button.png)

Click the 'Switch to' button at the top of that history and then click 'done' at the very top left to return to
the 'Analyze Data' page.
![Switch to history](../../images/undelete.switchto.png)

## Dataset Collections

When you have multiple datasets that will be sent through the same analysis it can often be useful to place those
datasets in a dataset collection. When collections are used as input when running a tool, you're telling Galaxy to run
that tool on each  dataset in the collection using the same settings. This happens automatically and there's no need to
fill in the tool form more than once.

It may be helpful to metaphorically think of dataset collections as containers (or directories) in which you place
datasets. To create a dataset collection from datasets in a history:

1. Use multiselect to select the datasets you'd like to make into a collection.
1. From the drop down menu labeled 'For all selected...', choose one of the collections types.
1. Depending on the type you chose, some configuration and options may be displayed, and when you've chosen those,
  the collection is created and will show in the current history panel.
1. Note: that your datasets are still shown in the history alongside the new collection. You can hide them if you like
  using the multiselect menu.

### Viewing a dataset collection

You can view what datasets were inside a collection by clicking on the collection title. The history panel will be
replaced by a list of the collection contents and each are expandable as a normal dataset in a history would be. You
can click the '< Back to (you history name)' link at the top to return to the history view (see below for examples).

They current collection types are:

### Dataset pairs

A common pattern of dataset files are pairs of read files - often some form of fastq files - where one file contains
the forward reads and one file contains the reverse reads. Many bioinformatic tools accept these pairs and Galaxy can
further simplify this by placing both files into on 'Dataset Pair' collection. Only two files will be added to the
collection: forward and reverse.

![Dataset pairs](../../images/pair.png)

### Dataset list

Choose 'Dataset List' when you have a set of files that are of the same type and will be run through some similar
analysis. The datasets in a dataset list must have unique names (e.g. you cannot have two datasets in a dataset list
with the name '1.bed').

![Dataset lists](../../images/list.png)

### List of dataset pairs

Think of this as a collection of collections: multiple dataset pairs contained in a dataset list. The interface used
to create this is currently the most flexible and potentially most complicated. It will attempt to automatically pair
datasets sent to the interface based on the dataset names. You are free to select your own pairs, however, and change
the order of the collection. Click the help text at the top of the interface to see more information.

![List of dataset pairs](../../images/list-pairs.png)
