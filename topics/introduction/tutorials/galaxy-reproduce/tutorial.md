---
layout: tutorial_hands_on

title: "How to reproduce published Galaxy analyses"
zenodo_link: https://zenodo.org/record/1319069/files/iris.csv
level: Introductory
questions:
  - "How to reproduce published Galaxy results (workflows and histories)"
objectives:
  - "Learn how to load published data into Galaxy"
  - "Learn how to run a published Galaxy workflow"
  - "Learn how histories can be used, published, and re-used."
time_estimation: "1H"
key_points:
  - "TODO"
contributors:
  - foellmelanie
  - annefou
---


# Introduction
{:.no_toc}

This training will demonstrate how to reproduce analyses performed in the Galaxy framework. Before we start with the hands-on part, we would like to give you some information about Galaxy. 

Galaxy is a scientific workflow, data integration and data and analysis persistence and publishing platform. Galaxy is an open-source platform for accessible, reproducible, and transparent computational research. While Galaxy was started to allow non-bioinformaticians to analyze DNA sequencing data, it nowadays enables analysis tasks of many different domains including typical omics-type of analyses in biology, machine learning, ecology and climate science. Galaxy is easy to use because it is accessible via a web-browser and provides a graphical user interface which enables access to pre-installed tools and large computational resources. In Galaxy, all analyses are stored in so-called histories. The history keeps track of all the tools, tool versions and parameters that were used in the analysis. From such a history, a workflow can be extracted; this workflow can be used to easily repeat the analysis on different data. Both, histories and workflows can either be shared privately with colleagues or publicly, for example as part of a published manuscript. 

For more background information about Galaxy, have a look into the [Galaxy publication](https://academic.oup.com/nar/article/46/W1/W537/5001157). The very technical details about technologies that enable reproducible analyses within Galaxy are described in this [publication](https://www.sciencedirect.com/science/article/pii/S2405471218301406). 


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


## What does Galaxy look like?

Many different Galaxy servers exist. Some are public, some are private, some focus on a specific topic and others like the usegalaxy.* servers are more general. To reproduce published results it is highly recommended to use the same Galaxy server that was used in the original study. In case this was a private server that is not accessible to you, you might want to use one of the main Galaxy servers: usegalaxy.org, usegalaxy.eu, usegalaxy.org.au. To learn more about the different Galaxy servers visit the slides options for using Galaxy (https://training.galaxyproject.org/training-material/topics/introduction/tutorials/options-for-using-galaxy/slides.html#1). The particular Galaxy server that you are using may look slightly different than the one shown in this training. Galaxy instance administrators can choose the exact version of Galaxy they would like to offer and can customize its look and feel to some extent. The basic functionality will be rather similar across instances, so don’t worry! In this training we will use the European Galaxy server on which the original analysis was performed and shared. 


> ### {% icon hands_on %} Hands-on: Log in or register
> 1. Open your favorite browser (Chrome/Chromium, Safari or Firefox, but not Internet Explorer/Edge!)
> 2. Browse to a [Galaxy instance](https://galaxyproject.org/use/) of your choice
> 3. Choose *Login or Register* from the navigation bar at the top of the page
> 4. If you have previously registered an account with this particular instance of Galaxy (user accounts are *not* shared between public servers!), proceed by logging in with your registered *public name*, or email address, and your password.
>
>    If you need to create a new account, click on *Register here* instead.
>
{: .hands_on}


The Galaxy interface consists of three main parts:

1. The available tools are listed on the left
2. Your analysis history is recorded on the right
3. The central panel will let you run analyses and view outputs

![Galaxy ecosystem]({% link shared/images/galaxy_interface.png %})


# Create a history and load data into it

Each analysis in Galaxy starts by creating a new analysis history and loading data into it. Galaxy supports a huge variety of data types and data sources. Different ways of bringing data into Galaxy are explained here (https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/get-data/slides.html#1 ). To reproduce published results, the data needs to be loaded from the public repository where the authors have deposited the data. This is most often done by importing data via a web link. 


> ### {% icon hands_on %} Hands-on: Create history
>
> 1. Make sure you start from an empty analysis history.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. **Rename your history** to be meaningful and easy to find. For instance, you can choose **Reproduction of published Galaxy results** as the name of your new history.
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
{: .hands_on}


> ### {% icon comment %} Background about the dataset
> The Iris flower data set, also known as Fisher’s or Anderson's Iris data set, is a multivariate dataset introduced by the British statistician and biologist Ronald Fisher in his 1936 paper ({% cite Fisher1936 %}).
> Each row of the table represents an iris flower sample, describing its species and the dimensions in centimeters of its botanical parts, the sepals and petals.
> You can find more detailed information about this dataset on its dedicated [Wikipedia page](https://en.wikipedia.org/wiki/Iris_flower_data_set).
{: .comment}


> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. {% tool [Import](upload1) %} the file `iris.csv` from [Zenodo](https://zenodo.org/record/1319069/files/iris.csv) or from the data library (ask your instructor)
>
>    ```
>    https://zenodo.org/record/1319069/files/iris.csv
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
>
> 2. **Rename** {% icon galaxy-pencil %} the dataset to `iris`
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
> 3. Check the **datatype**
>    - Click on the history item to expand it to get more information.
>    - The datatype of the iris dataset should be `csv`.
>    - **Change** {% icon galaxy-pencil %} the datatype *if* it is different than `csv`.
>      - Option 1: Datatypes can be **autodetected**
>      - Option 2: Datatypes can be **manually set**
>
>    {% snippet faqs/galaxy/datasets_detect_datatype.md datatype="datatypes" %}
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="csv" %}
>
> 4. Add an `#iris` tag {% icon galaxy-tags %} to the dataset
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
>    Make sure the tag starts with a hash symbol (`#`), which will make the tag stick not only to this dataset, but also to any results derived from it.
>    This will help you make sense of your history.
>
{: .hands_on}

> ### {% icon comment %} Different types of datasets
> Some input datasets might need more specialized treatment than explained here. 
> Collections contain several single dataset tied together. In case a workflow input requires a collection, you’ll need to build a collection out of your files after uploading them.
> A specialized training explains how to use collections (https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/collections/tutorial.html ). 
>
> {% snippet faqs/galaxy/collections_build_list.md %}
>
> A few data files contain more than one subfile. These are uploaded via the composite data function. 
> TODO: FAQ: Uploading composite data
> In case you want to run a published Galaxy workflow on your own data, you can find explanations about the options to upload your own data here (https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/get-data/slides.html#1). 
{: .comment}


# Import and run a Galaxy workflow

Galaxy workflows may be published either directly via the Galaxy server or on public workflow repositories such as WorkflowHub (https://workflowhub.eu/ ). The workflow may be present one of the three ways: 1) as a .ga file or url link, which needs to be imported into Galaxy, 2) as a link from a personal Galaxy server account that needs to be added to the own Galaxy account, 3) as a link that directly starts running the workflow in a specific Galaxy server. 


> ### {% icon hands_on %} Import and run workflow available as .ga file or link
>
> 1. **Import** the workflow either via url directly from [Zenodo](TODO) or by uploading the .ga file
>
>    ```
>    TODO Zenodo link
>    ```
>    {% snippet faqs/galaxy/workflows_import.md %}
>
> 2. **Start** the workflow by clicking on the run symbol on the last column in the workflow overview list
>
> 3. Select the `iris` dataset as the input dataset.
>
> 4. **Run** the workflow
>
>    > ### {% icon question %} Question
>    >
>    > How many history items do you have after running the workflow?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > 12, out of which 4 are shown and 8 hidden (at top of history right ander the history name)
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

TODO: box with information about hidden datasets
TODO: box with information about option 2 and 3

By starting the workflow all jobs are sent to the Galaxy cluster for analysis. Sometimes it can take a bit until the datasets show up in your history. The jobs are processed one after the other or in parallel if the same input is used for several steps. Grey means waiting to run, yellow means running and green means finished. Red means there was an error.

{% snippet faqs/galaxy/analysis_troubleshooting.md %}

> ### {% icon comment %} The tool stays grey
> This scenario will likely not happen with this training analysis, but might happen with a real workflow. The tool runs are scheduled on the computing infrastructures according to their consumption of cores and memory. Thus, tools that are assigned lots of cores and/or memory need to wait until an appropriate computing spot is available. Several Galaxy server display the current computational load which gives you an idea how busy it currently is. 
{: .comment}

Each dataset that turns green can already be inspected even though later datasets are still running. The second part of the training will focus on how to inspect datasets in a history. 


# Inspecting the analysis history

Each history item represents one dataset, except when collections are used. History items are numbered, duplicates are not possible because any type of dataset manipulation will automatically generate a separate dataset in the history. Some tools produce several output files and therefore the number of history items can be larger than the number of steps in a workflow. Dataset names in the analysis history are not relevant for the tool run, therefore they can be adjusted in order to make the history easier to understand. The default name of a history item is composed of the tool name that was run and the history item number(s) of the input file(s), e.g. 'Unique on data 5'

> ### {% icon hands_on %} Inspect history
>
> 1. **Vizualize** the two scatter plots by clicking on their eye icons (view data)
>
>    {% snippet faqs/galaxy/features_scratchbook.md %}
>
> 2. **Show** all datasets by clicking on `hidden` on top of your history, right below the history name
>
> 3. **Compare** the `Convert csv to tabular` file with the `Datamash` file side by side
>
> 4. **Track** how the `Datamash` results where obtained by clicking on the on the `Datamash` item in the history and then on its i icon (view details). The performed operations can be found in the section `Tool parameters`
>
{: .hands_on}


> ### {% icon question %} Questions
>
> 1. What are the different Iris species?
> 2. How many lines has the `Convert csv to tabular` file?
> 3. Which column of the  `Remove beginning` file contains sepal length and which petal length?
>
> > ### {% icon solution %} Solution
> >
> > 1.  The 3 different Iris species are:
> >     - setosa
> >     - versicolor
> >     - virginica
> > 
> > 2. 151 lines (by clicking on the file one can see the line count under its name)
> > 
> > 3. Column 1 and 3 (the dataset was generated by removing the header line from data 2, thus the content of the columns is the same as in data file 2)
>    {: .question}


# Manipulating the analysis

Maybe you are interested in changing some of the original tool parameters and see how this changes the result. The parameter changes can be either done inside the workflow editor, which makes it easy to change many parameters at once (TODO link workflow training) or directly in the history. To repeat a single analysis step with new parameters this can be done directly in the Galaxy history with the `re-run` option. 

> ### {% icon hands_on %} Manipulate the analysis steps
>
> 1. **Rerun** the Scatterplot to plot Sepal length vs. Petal length
>
>    {% snippet faqs/galaxy/tools_rerun.md %}
>
>>   - *"Column to plot on x-axis"*: `1`
>    - *"Column to plot on y-axis"*: `3`
>
{: .hands_on}


# Additional Exercise: Import a published analysis history and explore it

Often not only workflows and raw data are published but also the full Galaxy histories. These histories can either be inspected via their provided link, or imported in order to enable manipulating them in your own Galaxy account. 

> ### {% icon hands_on %} Hands-on: Importing a history
>
>
> 1. {% tool [Import](upload1) %} a published history shared via an EU Server account
>
>    ```
>    https://usegalaxy.eu/u/annefou/h/galaxy101-for-everyone-diamond-dataset
>    ```
>
>    {% snippet faqs/galaxy/histories_import.md %}
> 
{: .hands_on}

TODO: Questions/tasks for this history

# Conclusion
{:.no_toc}

{% icon trophy %} Well done! You have just reproduced your first analysis in Galaxy.

