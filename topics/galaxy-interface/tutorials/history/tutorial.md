---
layout: tutorial_hands_on
redirect_from:
  - /topics/galaxy-ui/tutorials/history/tutorial

title: "Understanding Galaxy history system"
zenodo_link: ""
level: Introductory
questions:
  - "How do Galaxy histories work?"
objectives:
  - "Gain understanding on navigating and manipulating histories"
time_estimation: "30m"
key_points:
  - "In Galaxy, data is stored in datasets and organized in histories"
  - "The right pane shows the current history. It displays the name of the history, the tags and annotations assigned to it, and the size of the data in it"
  - "The datasets are shown differently depending on their status: queued, running, ok, error, or paused"
  - "The summary view for each dataset can be expanded to show additional information"
  - "Datasets deleted are hidden from history. They can be restored to the history, or permanently deleted from the server, or purged. Purged files cannot be restored"
  - "Histories can also be deleted and purged, similar to datasets"
  - "Datasets can be organized into collections within a history with search and multi-select. Jobs applied to collections use the same settings for each dataset in the collection"
contributions:
  authorship:
  - nekrut
  editing:
  - polkhe
  - martenson
  - gallardoalba
  - shiltemann
  - bebatut
  - bgruening
  - hexylena
  - ahmedhamidawan
subtopic: histories

recordings:
- youtube_id: C3_HSgbbDWc
  speakers:
  - shiltemann
  captioners:
  - shiltemann
  length: 18M
  date: '2023-05-19'
  galaxy_version: '23.01'

---

{% snippet faqs/gtn/galaxy_tested_with.md version="24.1" %}

When data is uploaded from your computer or analysis is done on existing data using Galaxy, each output from those steps
generates a dataset. These datasets (and the output datasets from later analysis on them) are stored by Galaxy in
**Histories**.

# The history panel

All users have one 'current' history, which can be thought of as **a workspace** or **a current working directory** in
bioinformatics terms. Your current history is displayed in the right hand side of the main 'Analyze Data' Galaxy page in
what is called the history panel.

![Screenshot of the Galaxy UI with the toolbox, center panel, and history visible](../../images/history.png "Galaxy History is simply the right panel of the interface. It lists all datasets uploaded or produced during an analysis. Galaxy was designed around history with the idea of <i>reproducibility</i> in mind!")

The history panel displays output datasets in the order in which they were created, with the oldest/first shown at the
bottom. As new analyses are done and new output datasets are generated, the newest datasets are added to the top of the
the history panel. In this way, the history panel displays the history of your **analysis over time**.

**Users that have registered an account and logged in can have many histories** and the history panel allows switching
between them and creating new ones. This can be useful to organize different analyses.

**Anonymous users** (if your Galaxy allows them) are users that have not registered an account. Anonymous users are
only allowed one history. Users are encouraged to register and log in with the benefit that they can work on many histories and switch between them.

> <warning-title>Anonymous Users: Beware</warning-title>
> The histories of anonymous users are only associated through your browser's session. **If you close the browser or
> clear you sessions - that history will be lost!** We can not recover it for you if it is.
{: .warning}

## Global history controls

{% snippet faqs/galaxy/histories_top_control_buttons.md box_type="none"%}

### {% icon new-history %} "**Create new history**"

{% snippet faqs/galaxy/histories_create_new.md box_type="none"%}

### {% icon switch-histories %} "**Switch to history**"

{% snippet faqs/galaxy/histories_switch.md box_type="none"%}

### {% icon galaxy-history-options %} "**History options**" 

{% snippet faqs/galaxy/histories_options.md box_type="none"%}

# History manipulation

## Renaming a history

All histories begin with the name 'Unnamed history'. Non-anonymous users can rename the history as they see fit:

{% snippet faqs/galaxy/histories_rename.md box_type="none"%}

## Tagging a history

{% snippet faqs/galaxy/histories_tagging.md box_type="none"%}

## Annotating a history

{% snippet faqs/galaxy/histories_annotation.md box_type="none"%}

## Undeleting ... deleted histories

{% snippet faqs/galaxy/histories_undelete.md box_type="none"%}




## History size, storage selection, and views

The lower part of the history header contains a number of buttons:

![Buttons for figuring out history size, selecting storage location, and controlling views](../../../../shared/images/history_size_storage_views.png "Using these buttons one can obtain an information about history's disk footpring, change storage options, and control dataset views.")

> <warning-title>Some buttons are instance-specific</warning-title>
> The list of buttons shown above may vary depending on which Galaxy instance you are using. For example, at the time of writing the **Preferred storage** {% icon galaxy-history-storage-choice %} button is only available on https://usegalaxy.org.
{: .warning}

- {% icon galaxy-history-size %} **History size** - shows history storage overview in the central pane of the interface.
- {% icon galaxy-history-storage-choice %} **Preferred storage location** - allows users to specify where history datasets will be stored. This button is only available on Galaxy instances with scratch storage such as usegalaxy.org. Scratch storage allows users to much larger storage allocation for a limited amount time.
- {% icon galaxy-show-active %} **Show active** - shows active (non-deleted and non-hidden) datasets in the history.
- {% icon galaxy-delete %} **Include deleted** - include deleted datasets into the history view. If you delete a dataset is does not disappear unless you explicitly purge it.
- {% icon galaxy-show-hidden %} **Include hidden** - include hidden datasets into the history view. Any dataset in history can be hidden. For example, workflow executions frequently hide intermediate datasets so that they do now complicate history view. It is a way to hide non-importnat datasets from the view. 

# History datasets

So far we only discussed functions and controls affecting the *entrire* history. Yet history is a collection of datasets. Now it is time to discuss interfacse elements of individual datasets.  

## Datasets can be individual or bundled into collections

A history dataset can exist by itself, as an independent entity, or as a part of a **collection**. Collections make it possible to analyze datasets with hundreds of thousands of samples.

{% snippet faqs/galaxy/histories_datasets_vs_collections.md box_type="none"%}

## United Colors of Galaxy: Dataset states

{% snippet faqs/galaxy/histories_dataset_colors.md box_type="none"%}

## Dataset snippet in detail

{% snippet faqs/galaxy/histories_dataset_item.md box_type="none"%}

## Managing Datasets Individually

### Hiding and unhiding datasets

Datasets in Galaxy history can be hidden. This is useful for reducing complexity of history. For example, some intermediate datasets generating during an analysis of workflow execution are not important and there is no need to see them. 

#### Hiding datasets

{% snippet faqs/galaxy/datasets_hidden.md box_type="none"%}

#### Unhiding datasets

{% snippet faqs/galaxy/datasets_unhidden.md box_type="none"%}

### Deleting and undeleting datasets

You can **delete** any dataset in your history. Unless you explicitly tell Galaxy to delete a dataset permanently (see below) this does not immediately remove the dataset from Galaxy" **it is reversible**. When you delete a dataset from the history, it will be removed from the panel but (just like hidden datasets).

#### Deleting datasets

{% snippet faqs/galaxy/datasets_deleting.md box_type="none"%}

#### Undeleting datasets

{% snippet faqs/galaxy/datasets_undelete.md box_type="none"%}

> <comment-title>Admins may purge your deleted datasets</comment-title>
> Depending on the policy of your Galaxy server, administrators may run scripts that search for and purge the
> datasets you've marked as deleted. Often, deleted datasets and histories are purged based on the age of the deletion
> (e.g. datasets that have been marked as deleted for 90 days or more). Check with the administrators of your Galaxy instance to
> find out the policy used.
{: .comment}

## Tagging datasets

{% snippet faqs/galaxy/datasets_add_tag.md box_type="none"%}

## Managing Multiple Datasets Easily

{% snippet faqs/galaxy/datasets_multiple.md box_type="none"%}

<!-- Including search to deduplicate -->
{% include topics/galaxy-interface/tutorials/search/search.md %}

