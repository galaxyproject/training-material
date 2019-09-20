---
layout: tutorial_hands_on

title: "Galaxy 101 for everyone"
zenodo_link: https://zenodo.org/record/1319069/files/iris.csv
level: Introductory
questions:
  - "What are the differences between the Iris species?"
objectives:
  - "Familiarize yourself with the basics of Galaxy"
  - "Learn how to obtain data from external sources"
  - "Learn how to tag datasets"
  - "Learn how to run tools"
  - "Learn how histories work"
  - "Learn how to share your work"
time_estimation: "1H30M"
key_points:
  - "Galaxy provides an easy-to-use graphical user interface for often complex command-line tools"
  - "Galaxy keeps a full record of your analysis in a history"
  - "Workflows enable you to repeat your analysis on different data"
  - "Galaxy can connect to external sources for data import and visualization purposes"
  - "Galaxy provides ways to share your results and methods with others"
contributors:
  - annefou
  - nagoue
  - chrisbarnettster
  - michelemaroni89
  - olanag1
  - tnabtaf
---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

This practical aims to familiarize you with the Galaxy user interface. 
It will teach you how to perform basic tasks such as importing data, running tools, working with histories, creating workflows, and sharing your work.
Not everyone has the same background and that's ok! 

**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

> ### {% icon comment %} Background
> The Iris flower data set or Fisher’s Iris data set is a multivariate data set introduced by the British statistician and biologist Ronald Fisher in his 1936 paper ({% cite Fisher1936 %}). 
> Each row of the table represents an iris flower, including its species and dimensions of its botanical parts, sepal and petal, in centimeters.
> For more history of this dataset read here [Wikipedia](https://en.wikipedia.org/wiki/Iris_flower_data_set).
{: .comment}

## What does Galaxy look like?

> ### {% icon hands_on %} Hands-on: Log in or register
> Browse to your [Galaxy instance](https://galaxyproject.org/use/) and log in or register.
> 1. Open your favorite browser (Chrome, Safari or Firefox as your browser, not Internet Explorer!)
> 2. Browse to your Galaxy instance
> 3. Log in or register
>
>   > ### {% icon comment %} Different Galaxy servers
>   > The particular Galaxy server that you are using may look slightly different than the one shown in this training. Don't worry!
>   {: .comment}
>
{: .hands_on}


The Galaxy interface consists of three main parts:

1. The available tools are listed on the left
2. Your analysis history is recorded on the right
3. The central panel will let you run analyses and view outputs

![Galaxy ecosystem]({{ site.baseurl }}{% link shared/images/galaxy_interface.png %})


# Create a history

Galaxy allows you to create histories. Overall a history represents an experimental lab book, or a recipe very much like a cooking recipe with a list of ingredients (datasets) and a set of instructions 
(pipeline of operations) that describes how to prepare or make something (such as a plot, or even a new dataset).
The order of operations is important as very often the next operation takes as input the result of the previous operations. For instance, when baking
a cake, you would first sift flour and then mix it with eggs as it would be impossible to sift flour afterwards.
That is what we call a pipeline.

Then the finalized pipeline can be serialized as a workflow.A workflow is the concatenation of one or multiple histories as a series of building blocks 
for replicating an experimental result or a recipe. If we use cooking as an analogy, a workflow could represent an entire menu with all the recipes for each meal.
In other words, using a workflow makes it possible to apply the same procedure to a different dataset, just by changing the input. 

> ### {% icon hands_on %} Hands-on: Create history
>
> Make sure you start from an empty analysis history.
>
>    {% include snippets/create_new_history.md %}
>![Rename the history]({{ site.baseurl }}{% link shared/images/rename_history.png %})

> **Rename your history** to be meaningful and easy to find. For instance, you can choose **Galaxy 101 for everyone** as the name of your new history. 
>    {% include snippets/rename_history.md %}
{: .hands_on}

## Upload Iris dataset

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Import `iris.csv` from [Zenodo](https://zenodo.org/record/1319069/files/iris.csv) or from the data library (ask your instructor)
>
>    ```
>    https://zenodo.org/record/1319069/files/iris.csv
>    ```
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
>    As default, Galaxy takes the link as name, so rename them.
>
> 2. Rename the dataset to `iris`
>
>    {% include snippets/rename_dataset.md %}
>
> 3. Check the datatype. The datatype of the iris dataset is `csv`. Change datatype
> if it is different then `csv`.
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 4. Add the tag `iris` to the dataset
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

# Pre-processing

A pre-processing step can be required to proceed analysis. In this case, format convertion and header removal have to be processed.


## Convert dataset **csv_to_tabular**

> ### {% icon hands_on %} Hands-on: Converting dataset format 
>    
>    * Click on the {% icon galaxy-pencil %} **pencil icon** for the dataset to edit its attributes
>    * In the central panel, click on the {% icon galaxy-gear %} **Convert** tab on the top
>    * Select `Convert CSV to tabular`
>    * Click the **Convert datatype** button
>
>    > ### {% icon comment %} Comment
>    > A New dataset is created in the history panel. You can click on the dataset name to get an overview
>    {: .comment}
>
> 2. Rename the dataset to `iris tabular`
>
>    {% include snippets/rename_dataset.md %}
> 3. Add the tag `preprocessing` to the dataset
>
>    {% include snippets/add_tag.md %}
> 4. Inspect the generated file by clicking on the {% icon galaxy-eye %} (eye) icon (**View data**)
{: .hands_on}


## Remove header

> ### {% icon comment %} Search bar
>
> Different Galaxy servers may have tools available under different sections, therefore it is often useful to use the **search bar** at the top of the tool panel to find your tool.
>
> Additionally different servers may have multiple, similarly named tools which accomplish similar functions. For these tutorials, you should select precisely the one that is described. However, in your real analyses, you'll need to search among the various options to find the one that works for you.
>
{: .comment}

> ### {% icon hands_on %} Hands-on: Removing header
>
> 1. **Remove beginning** {% icon tool %} with the following parameters:
>
>    - *Remove first*: `1` to remove the first line only.
>    - *from*: {% icon param-file %}: select **iris tabular**
>    - **Execute**
>
>    > ### {% icon comment %} Comment
>    >
>    > Use the **tools search box** to find **Removing beginning** {% icon tool %}. 
>    {: .comment}
>
>    ![Settings for the `Remove beginning` tool](../../images/101_foreveryone_remove_beginning.png)

> 2. Rename the dataset to `iris clean`
>
>    {% include snippets/rename_dataset.md %}
> 3. Add the tag `clean` to the dataset
>
>    {% include snippets/add_tag.md %}
> 4. Inspect the generated file by clicking on the {% icon galaxy-eye %} (eye) icon (**View data**)
{: .hands_on}



# Data Analysis: What does the dataset contain?

Now we are going to inspect the dataset using simple tools in order to get used to galaxy interface and answer basic questions.

## How many different species are in the dataset?

> ### {% icon hands_on %} Hands-on: Filtering dataset 
>
> 1.   **Cut** {% icon tool %} with the following parameters:
>    - *"Cut columns"* should be changed to `c5`
>    - *"Delimited by"* should be kept to `Tab`
>    - *"From"* select `iris clean` dataset
>
> 2.    Click **Execute**
>
> 3. Wait for the job to finish
>
> 4. View the resulting file (with the {% icon galaxy-eye %} (eye) icon). Use the output file as input of the second tool.
>
> 5. Add the tag `analysis` to the output dataset
> 
> 6.   **Unique** {% icon tool %} with the following parameters:
>    - *"File to scan for unique values"* select your last output file
>
> 7. Click **Execute**
>
> 8. Wait for the job to finish
>
> 9. View the resulting file (with the {% icon galaxy-eye %} (eye) icon).
>
> 10. Add the tag `analysis` to the output dataset
>
> 11. Examine the output file
>
{: .hands_on}

> ### {% icon question %} Questions
> 
> 1. How many different species are in the dataset?
> 2. What are the different Iris species?
> 
> > ### {% icon solution %} Solution
> >
> > 1. There are 3 species.
> > 2. The 3 different Iris species are:
> > - setosa
> > - versicolor
> > - virginica
> {: .solution}
{: .question}


> ### {% icon hands_on %} Hands-on: Grouping dataset
>
> Another way round to answer this question with only one tool:
>
> 1.   **Group** {% icon tool %} with the following parameters:
>    - *"Select data"* select `iris clean` dataset
>    - *"Group by column"* should be changed to `Column: 5`
>
> 2. Click **Execute**
>
> 3. Wait for the job to finish
>
> 4. Add the tag `analysis` to the output dataset
>
> 5. View the resulting file (with the {% icon galaxy-eye %} (eye) icon).
>
{: .hands_on}


> ### {% icon question %} Question
> 1. How many different species are in the dataset?
> 2. What are the different Iris species?
>
> > ### {% icon solution %} Solution
> >
> > 1. There are 3 species.
> > 2. The 3 different Iris species are:
> > - setosa
> > - versicolor
> > - virginica
> {: .solution}
{: .question}

## How many samples by species?

> ### {% icon hands_on %} Hands-on: Grouping dataset and adding information
>
> 1.   **Group** {% icon tool %} with the following parameters:
>    - *"Select data"* select `iris clean` dataset
>    - *"Group by column"* should be changed to `Column: 5`
>    - *"Insert operation"* click on the icon.
>    - *"Type"* select `Count`
>    - *"On column"* select `Column: 1`
> 2. Click **Execute**
> 3. Wait for the job to finish
> 4. View the resulting file (with the {% icon galaxy-eye %} (eye) icon).
> 5. Add the tag `analysis` to the output dataset
>
{: .hands_on}

> ### {% icon question %} Question
> 1. How many samples by species are in the dataset?
>
> > ### {% icon solution %} Solution
> > 
> > 1. We have 50 samples per species:
> > ```
> > 1 | 2
> > ---- | ----------
> > setosa | 50
> > versicolor | 50
> > virginica | 50
> > ```
> {: .solution}
{: .question}

# Analysis: How to differentiate the different Iris species?

Our objective is to find what differentiate the different Iris species. We know that we have **3** species of iris flowers, with
**50** samples for each:
- setosa,
- versicolor,
- virginica.

These species look very much alike as shown on the figure below.

![3 species of Iris flowers](../../images/iris_flowers.png "3 species of Iris flowers")

And our objective is to find out whether the features we have been given for each species can help us to highlight the differences between the 3 species.

In our dataset, we have the following features measured for each sample:
- Petal length
- Petal width
- Sepal length
- Sepal width

> ### {% icon comment %} petal and sepal
> The image below shows you what is a sepal and petal.
> ![Sepal and petal](../../images/iris_sepal_petal.png "Sepal and petal of Iris flowers")
{: .comment}

## Summary and descriptive statistics with **Datamash**

> ### {% icon hands_on %} Hands-on: Get the mean and sample standard deviation of Iris flower features
>
> 1. **Datamash** {% icon tool %} with the following parameters:
>    - *Input tabular dataset*: {% icon param-file %}: select **iris tabular**
>    - *"Group by fields"*: `5`
>    - *"Input file has a header line"*: `Yes`
>    - *"Print header line"*: `No`
>    - *"Sort input"*: `Yes`
>    - "Print all fields from input file": `No`
>    - *"Ignore case when grouping"*: `Yes`
>    - In *"Operation to perform on each group"*:
>        - {% icon param-repeat %} *"Insert Operation to perform on each group"*
>            - *"Type"*: `Mean`
>            - *"On column"*: `c1`
>        - {% icon param-repeat %} *"Insert Operation to perform on each group"*
>            - *"Type"*: `Sample Standard deviation`
>            - *"On column"*: `c1`
>        - {% icon param-repeat %} *"Insert Operation to perform on each group"*
>            - *"Type"*: `Mean`
>            - *"On column"*: `c2`
>        - {% icon param-repeat %} *"Insert Operation to perform on each group"*
>            - *"Type"*: `Sample Standard deviation`
>            - *"On column"*: `c2`
>        - {% icon param-repeat %} *"Insert Operation to perform on each group"*
>            - *"Type"*: `Mean`
>            - *"On column"*: `c3`
>        - {% icon param-repeat %} *"Insert Operation to perform on each group"*
>            - *"Type"*: `Sample Standard deviation`
>            - *"On column"*: `c3`
>        - {% icon param-repeat %} *"Insert Operation to perform on each group"*
>            - *"Type"*: `Mean`
>            - *"On column"*: `c4`
>        - {% icon param-repeat %} *"Insert Operation to perform on each group"*
>            - *"Type"*: `Sample Standard deviation`
>            - *"On column"*: `c4`
>
> 2. Rename the dataset to `iris summary and statistics`
>
>    {% include snippets/rename_dataset.md %}
>
> 3. Add the tag `analysis` to the dataset
>
>    {% include snippets/add_tag.md %}
> 4. Inspect the generated file by clicking on the {% icon galaxy-eye %} (eye) icon (**View data**)
>
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. Can we differentiate the different Iris flower species?
>
> > ### {% icon solution %} Solution
> >
> > 1. From the results, we can see that the average Setosa petal length is lower than 1.5 with a relatively small standard deviation (<0.2).
> > The same can be observed for Setosa petal widths. These numbers are much smaller (width and length) then versicolor and virginica petals.
> > We can then use these characteristics to differentiate Iris Setosa from the two other species (versicolor and virginica). On the other hand,
> > we cannot easily differentiate Iris Versicolor from Iris Virginica. Further analysis is necessary.
> >
> {: .solution}
>
{: .question}

## Visualize Iris dataset with **Scatterplot w ggplot2**

Let's visualize the Iris dataset to see how the features depend on each other, and 
check whether we can spot any immediate patterns.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Scatterplot w ggplot2** {% icon tool %} with the following parameters:
>    - *Input tabular dataset*: {% icon param-file %}: select **iris clean**
>    - *"Column to plot on x-axis"*: `1`
>    - *"Column to plot on y-axis"*: `2`
>    - *"Plot title"*: `Sepal length as a function of sepal width`
>    - *"Label for x axis"*: `Sepal length`
>    - *"Label for y axis"*: `Sepal width`
>    - In *"Advanced Options"*:
>        - *"Data point options"*: `User defined point options`
>            - *"relative size of points"*: `2.0`
>        - *"Plotting multiple groups"*: `Plot multiple groups of data on one plot`
>            - *"column differentiating the different groups"*: `5`
>            - *"Color schemes to differentiate your groups"*: `Set 2 - predefined color pallete (discrete, max=8 colors)`
>        - *"Axis title options"*: `Default`
>        - *"Axis text options"*: `Default`
>        - *"Plot title options"*: `Default`
>        - *"Axis scaling"*: `Automatic axis scaling`
>
> 2. Click **Execute** to perform the graph. Your new output dataset will look something like this:
>
>    ![Contents of the `Group` output dataset](../../images/101_foreveryone_scatter.png)
>
> 3. Add the tag `plot` to the dataset
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. What does this scatter plot tell us about Iris species?
> 2. Make a new scatter plot with Petal length versus Petal width. Can we differentiate the three Iris species?
>
> > ### {% icon solution %} Solution
> >
> > 1. We get similar results than with Summary and statistics: Iris Setosa can clearly be distingished from Iris versicolor and
> > iris virginica. We can also see that Sepal width and length are not sufficient features to differentiate Iris versicolor from Iris
> > virginica.
> > 2. **Scatterplot w ggplot2** {% icon tool %} with the following parameters:
> >  - *Input tabular dataset*: {% icon param-file %}: select **iris clean**
> >  - *"Column to plot on x-axis"*: `3`
> >  - *"Column to plot on y-axis"*: `4`
> >  - *"Plot title"*: `Petal length as a function of petal width`
> >  - *"Label for x axis"*: `Petal length`
> >  - *"Label for y axis"*: `Petal width`
> >  - In *"Advanced Options"*:
> >      - *"Data point options"*: `User defined point options`
> >          - *"relative size of points"*: `2.0`
> >      - *"Plotting multiple groups"*: `Plot multiple groups of data on one plot`
> >          - *"column differentiating the different groups"*: `5`
> >          - *"Color schemes to differentiate your groups"*: `Set 2 - predefined color pallete (discrete, max=8 colors)`
> >      - *"Axis title options"*: `Default`
> >      - *"Axis text options"*: `Default`
> >      - *"Plot title options"*: `Default`
> >      - *"Axis scaling"*: `Automatic axis scaling`
> >
> > Click **Execute** to perform the graph. Your new output dataset will look something like this:
> >
> >   ![Contents of the `Group` output dataset](../../images/101_foreveryone_scatter_petal.png)
> > 3. Add the tag `plot` to the dataset
> >
> > We can better differentiate the 3 Iris species but for some samples the petal length versus width is still insufficient
> > to differentiate Iris versicolor from Iris virginica. And as before, Iris setosa can easily distinguish from the two other species.
> {: .solution}
{: .question}

# Share your work

One of the most important features of Galaxy comes at the end of an analysis. When you have published striking findings, it is important that other researchers are able to reproduce your in-silico experiment. Galaxy enables users to easily share their workflows and histories with others.

To share a history, click on the {% icon galaxy-gear %} icon in the history panel and select `Share or Publish`. On this page you can do 3 things:

1. **Make History Accessible via Link**. This generates a link that you can give out to others. Anybody with this link will be able to view your history.
2. **Make History Accessible and Publish**. This will not only create a link, but will also publish your history. This means your history will be listed under `Shared Data → Histories` in the top menu.
3. **Share with a user**. This will share the history only with specific users on the Galaxy instance.

> ### {% icon comment %} Permissions
> Different servers have different default permission settings. Some servers create all of your datasets completely private to you, while others make them accessible if you know the secret ID.
>
> Be sure to select **Also make all objects within the History accessible** whenever you make a history accessible via link, otherwise whomever you send your link to might not be able to see your history.
{: .comment}

> ### {% icon hands_on %} Hands-on: Share history and workflow
>
> 1. Share one of your histories with your neighbour.
> 2. See if you can do the same with your workflow!
> 3. Find the history and/or workflow shared by your neighbour. Histories shared with specific users can be accessed by those users in their {% icon galaxy-gear %} history menu under `Histories shared with me`.
{: .hands_on}

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
