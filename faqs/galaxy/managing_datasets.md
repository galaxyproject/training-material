---
title: Managing Datasets
area: Datasets
box_type: tip
layout: faq
contributors: [mruqqie]
---

**What are Datasets?**
- Datasets are the inputs and outputs of each step in an analysis project in Galaxy.

**Getting Datasets**

This can be done by:
- Clicking the {% galaxy-upload %} icon on the left, below 'search tools'.
- FTP upload (incase of large files).

**Dataset Icons and Text**
- In the right corner:
  - {% galaxy-eye %}: To display data in browser
  - {% galaxy-pencil %}: To edit attributes
  - {% galaxy-cross %}: The delete icon
  
- In the lower right corner:
  - Edit dataset tags
  - Edit dataset annotation
  
- In the upper left corner:
  - Dataset name
  - Dataset size/number of lines (actual or estimated)
  - Format datatypes
  - Database
  - {% galaxy-info %} Info (optional)
  
- In the lower left corner
  - Download
  - View Details
  - Run this job again {% galaxy-refresh %}
  - Display in trackster (optional)
  - Display at UCSC main (optional)
  - Display at Ensembl Current (optional)
  
**Data size and disk Quotas**
- The size limit for a file loaded using FTP is 50G.
- The size limit for a job's output is:
  - 50G on the Test server
  - 200G on the Main server

**Format**
- The format of a dataset is ideally defined by the assigned datatype attribute.
- To initially assign a dataset's datatype attribute, the uploaded/imported file can be specified with some import tools or be named with the appropriate file extension.
- To specify, modify or correct a dataset's datatype attribute after upload:
  - Click on the "pencil" icon {% galaxy-pencil %} in the right corner of the dataset's box to reach the "Edit Attributes" form.
  - Click on the "Change data type" section of the form to make changes.
  - Click on Save.
- To transform a dataset format (original → new datatype attribute), use one of the many tools in the Convert Formats group.

**Visualize**
- Click on the eye icon {% galaxy-eye %} to display data in browser. This works for most datatypes except compressed datatypes such as BAM.
- Direct links to view a dataset within a browser may include:
  - Trackster "Galaxy Track Browser (GTB)"
  - GeneTrack
  - UCSC
  - Ensembl

**Copy**
- To copy the datasets from one history to another history:
  - From the right history pane's top Options menu, select Copy Datasets.
  - On the form in the center pane, specify the From and To history/histories.
- To Copy a Hidden dataset:
  - In the From histories right pane, click on the gear icon {% galaxy-gear %} to unhide hidden datasets.
  - Use the toogle at the top of the history panel (directly below the history name) to view them.
  
**Clone**
- Cloning a history creates the exact copy of that history. The cloned history's default name is the original history's name prefixed by *Clone of*.
- Options are:
  - Clone all history items, including deleted items.
  - Clone only items that are not deleted.

**Delete vs Delete Permanently**
- Deleted datasets and histories can be recovered by users as they are retained in Galaxy for a time period (several months) set by the instance administrator.
- Permanently deleted datasets and histories cannot be recovered by the user or administrator.
- To review or adjust an individual dataset, click on the name to expand it.
  - If it is only deleted, but not permenently deleted (purged), you'll see a message with links to recover or to purge.
    - Click on Undelete it to recover the dataset.
    - Click on Permenently remove it from disk to purge the dataset and remove it from the account quota calculation.
- To review or adjust multiple datasets in batch:
  - Click on the "checked box" icon {% galaxy-selector %} near the top right of the history panel to switch into "Operations on Mulitple Datasets" mode.
  - Several options to show, hide, delete, undelete, purge, and group datasets are available. A selection box will be available for each individual dataset.
  - Check the datasets you want to modify and chose your option.
- To avoid losing important data by accident is to clearly name all important histories and datasets.
  - To name a dataset:
    - Click on the pencil icon {% galaxy-pencil %} on the top right of a Dataset to reach the Edit Attributes form.
    - A dataset's primary Name, Info , Annotation, and Notes can be adjusted on the first tab, and other metadata attributes can be adjusted on the other tabs.
  - To name a history:
    - Click near the top of the right history pane where the default text Unnamed history is located. Enter the new name and and press "enter/return".

**Searching Datasets**
- A search bar is at the top of every History.
- Type text into the field and press "enter/return".
- Clear the search by:
  - Removing the text in the bar and pressing enter.
  - By pressing the ESC key while the text is highlighted.
  - Clicking on the "clear search" button on the right side of the bar.
- To search for terms that include whitespaces (e.g. a dataset named 'My Dataset' which includes a space):
  - Enclose a search term with double quotes (eg. "My Data") and this will match the dataset in the example but exclude any dataset that may have matched an un-enclosed term (such as '~MyInterval' or 'Some Data').
- Search terms are applied to many of the 'fields' that describe a dataset and, by default, check all fields against each term.
  - To apply a search term to only one field, use the field name followed by an equals sign and then followed by the term (or double quote enclosed term). E.g. name="My Dataset" format=interval.
  - Here is a list of field names that can be used:
    - name or title: the name or partial name of the dataset(s)
    - database or genome_build: the name or partial name of the genome database (e.g. hg19)
    - format or file_ext: the file type/datatype of the dataset(s)
    - description: the text or partial text of description of the dataset(s) (the text just below the dataset title shown when the dataset is expanded)
    - info: the text or partial text of the info box of the dataset(s) (the text just below the format and database and above the download and info buttons shown when the dataset is expanded).
    - annotation: the text or partial text of your annotation on the dataset(s)
    - tag: the text or partial text of any applied you to the dataset(s)
