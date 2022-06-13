---
layout: tutorial_hands_on

title: Data Manipulation Olympics
tags:
- workflows
zenodo_link: 'https://zenodo.org/record/6638036'
questions:
- How can I do basic data manipulation in Galaxy?
- Which tools are available to convert, reformat, filter, sort etc my text-based data?
objectives:
- Familiarize yourself with data manipulation tools in Galaxy
- Perform basic text manipulation tasks in Galaxy
- Become comfortable converting text-based files in a variety of ways.
time_estimation: 30m
key_points:
- There are a lot of tools available in Galaxy for performing basic data manipulation tasks
- Bacis data manipulation is often needed between steps in a larger scientific analysis in order to connect outputs from one tool to input of another.
- There are often multiple ways/tools to achieve the same end result
- Having a basic understanding of data manipulation tools will make it easier to do exploratory data analysis
contributions:
  authorship:
    - shiltemann
  editing:
    - hexylena
  funding:
    - erasmusplus
level: Introductory

subtopic: next-steps
---


# Introduction
{:.no_toc}

Scientific analyses often consist of a number of tools that run one after the other, in order to go from the raw data to scientific insight. Between these specialized tools, simple data manipulation steps are often needed as a kind of "glue" between tools. For example, the output of tool A may produce a file that contains all the information needed as input for tool B, but tool B expects the columns in a different order. Or in genomic data analysis, some tools expect chromosome X to be listed as `chrX`, while others simply expect `X`. In these situations, extra data manipulation steps are needed to prepare files for input to analysis tools.

Galaxy has a large collection of tools to perform such basic data manipulation tasks, and becoming familiar with these operations will allow to perform your analysis more easily in Galaxy (and outside).


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Background

In this tutorial, we will use as our dataset a table with results from the Olympics, from the games in Athens in 1896 until Rio 2016. The objective is to familiarize you with a large number of the most important data manipulation tools in Galaxy. Much like the Olympics, there are many different disciplines (types of operations), and for each operation there are often multiple techniques (tools) available to athletes (data analysts, you) that are great for achieving the goal.


![image of olympic rings, logo and two atheletes around the words "Data Analysis Olympics"](./images/cover.jpg)


We will show you many of these commonly needed data manipulation operations, and some examples of how to perform them in Galaxy. We also provide many exercises so that you can train your skills and become a data manipulation Olympian!


# Upload Data

Before we can do any manipulation, we will need some data. Let's upload our table with olympics results now.

> ### {% icon hands_on %} Hands-on: Get data
>
> 1. {% tool [Import](upload1) %} the file `olympics.tsv` via link
>
>    ```
>    https://zenodo.org/record/6638036/files/olympics.tsv
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 2. **Expand** on the item in your history to see some metadata and a short preview of the contents.
>
> 3. **View** {% icon galaxy-eye %} the dataset by clicking on the eye icon.
>
>    > ### {% icon question %} Question
>    >
>    > 1. What is the format of the file?
>    > 2. What does each row represent?
>    > 3. How many lines are in the file?
>    > 4. How many columns?
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. When you expand the dataset, you will see `format: tabular`, this is another term for a tsv file.
>    > > 2. Each row represents an athlete's participation in an event. If an athlete competes in multiple events, there is a line for each event.
>    > > 3. Look at the expanded view in the history, there are ~270,000 rows in the dataset
>    > > 4. There are 14 columns in this file. There are multiple ways to find this answer:
>    > >    - Count the columns (only doable for small files)
>    > >    - In the expanded view, scroll sideways on the dataset preview, at the top the columns are numbered
>    > >    - Click on the {% icon galaxy-info %} i icon on the dataset, here you will find more detailed information about the file and the job that created it.
>    > >      At the bottom is also a preview (peek) of the dataset, and numbered columns
>    > >
>    > > ![a screenshot of the expanded view of the dataset in the history, it shows the datatype, number of lines in the file, and a preview
>    > > of the dataset with numbered columns](./images/columns-number.png)
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}



# Choose your adventure!

This tutorial is structured a bit differently than most. **You do not have to do the steps in the order they are presented below.** Every section in this tutorial uses the dataset you just uploaded (the `olympics.tsv` file) as input, so you can jump to any section in this tutorial right now if you have a particular data manipulation operation in mind you want to learn more about.




# File Format Conversion

The file we uploaded is a `.tsv` file. This stands for *tab-separated values*. This means that this is a file containing rows and columns, where a TAB character is used to signify a column ends and a new one begins. Galaxy is great at understanding tab-separated files files, and most of the data manipulation tools are designed to work with such files.

A similar format you may come across a lot in data science, is the `.csv` file, or *comma-separated values* file. This is the same as `.tsv`, but uses comma (`,`) characters to indicate new columns, instead of TAB (`\t`) characters.

Galaxy can convert these two formats into each other.



> ### {% icon hands_on %} Hands-on: Convert TSV to CSV
>
> 1. Convert the TSV file into a CSV file
>
>    {% snippet faqs/galaxy/datasets_convert_datatype.md conversion="csv (using 'Convert tabular to CSV')" %}
>
> 2. {% icon galaxy-eye %} **View the converted dataset**
>
>    > ### {% icon question %} Question
>    >
>    > 1. Where are the commas?
>    > 2. Why are some values in quotes?
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. Galaxy understands this fomat to be a table with rows and columns, so when you view the file, Galaxy hides the commas separating the columns, to make the file easier to read. If you download the file and view it in a text editor, you will see that there are 13 commas on each line, separating the 14 columns of the file.
>    > >
>    > > 2. If the data in a column contains a comma (e.g. in this file we have events such as `swimming 5,000 meters`), we put the value in quotes to signifiy that that comma is part of the data, not a column delimiter.
>    > >
>    > {: .solution}
>    {: .question}
>
>
{: .hands_on}


## Converting vs Changing the datatype

The **file format conversion** step changed the dataset itself (by changing TABs to commas), and therefore created another dataset in your history. It is also possible to **change the datatype**; this does not change the underlying file/data, but just tells Galaxy what type of file it is. Galaxy uses this information in a variety of ways, such as:
  - how to display the data when you click the {% icon galaxy-eye %} **View** button
  - which visualization options to offer
  - filtering the possible inputs in a tool form to only those of the correct datatype
  - which file format conversions to make available
  - etc


When you upload a file to Galaxy, by default it will attempt to auto-detect the file format. This is very often correct, but some times you may need to manually set the datatype to the correct value.


> ### {% icon hands_on %} Hands-on: Change datatype from CSV to TXT
>
> 1. **Change** the `csv` file into a `txt` file
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="txt" %}
>
> 2. {% icon galaxy-eye %} **View the dataset** again
>
>    > ### {% icon question %} Question
>    >
>    > 1. What do you see?
>    > 2. Why didn't this step create a new item in your history?
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. Since Galaxy now no longer know that this is a tabular file with rows and columns,
>    > > it displays the data as-is, no longer hiding the commas.
>    > >
>    > > 2. The file itself did not change, only the metadata. We simply told Galaxy that this file is not a `csv` file, but just a file containing text.
>    > >
>    > {: .solution}
>    {: .question}
>
> 3. **Change** the file back from `txt` to `tabular`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular" %}
>
{: .hands_on}




# Summary Statistics



# Sort by column

We have a lot of data in this file, but it isn't really ordered in any logical way. Let's change that.


> ### {% icon hands_on %} Hands-on: Sort table based on a column
>
> 1. We will sort the file in chronological order based on the year of the Olympic games
>
>    > ### {% icon question %} Questions
>    >
>    > 1. Which column contains the year?
>    > 2. Do we want ascending or descending order if we want the oldest games at the top?
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. Column 9
>    > > 2. The file should be sorted in ascending (increasing) order
>    > >
>    > {: .solution}
>    {: .question}
>
>
> 2. {% tool [Sort](sort1) %} with the following parameters:
>    - {% icon param-file %} *"Sort dataset"*: `olympics.tsv`
>    - {% icon param-select %} *"on column"*: `Column 9`
>    - {% icon param-select %} *"with flavor"*: `Numerical sort`
>    - {% icon param-select %} *"everything in"*: `Ascending order`
>    - {% icon param-text %} *"Number of header lines to skip"*: `1`
>
> 3. {% icon galaxy-eye %} **View** the sorted file.
>
>    > ### {% icon question %} Question
>    >
>    > Which athlete is listed at the top of the file now?
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. Giuseppe Rivabella
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

This is great, but maybe it would make more sense to sort alphabetically by athlete name *within each year*.

## Sort on multiple columns at once

So we want to sort twice, first by year, an then within each year, we sort again alphabetically by name. The sort tool can do this!


> ### {% icon hands_on %} Hands-on: Sort table based on a column
>
> 1. We will sort the file in chronological order based on the year of the Olympic games
>
>    > ### {% icon question %} Questions
>    >
>    > 1. Which column contains athlete names?
>    > 2. Do we want ascending or descending order if we want to sort alphabetically with A first?
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. Column 2
>    > > 2. The file should be sorted in ascending (increasing) order
>    > >
>    > {: .solution}
>    {: .question}
>
>
> 2. {% icon rerun %} **Rerun** the sort tool with the following parameters:
>    - All parameter from the first step should already be set for you, and should remain the same
>    - {% icon param-repeat %} Insert Column Selection
>      - {% icon param-select %} *"on column"*: `Column 2`
>      - {% icon param-select %} *"with flavor"*: `Alphabetical sort`
>      - {% icon param-select %} *"everything in"*: `Ascending order`
>
> 3. {% icon galaxy-eye %} **View** the sorted file.
>
>    > ### {% icon question %} Question
>    >
>    > Which athlete is listed at the top now? Which discipline (sport) did they compete in?
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. A. Tryfiatis-Trypiapis, who competed in the [Cycling Men's 12-Hours Race](https://en.wikipedia.org/wiki/Cycling_at_the_1896_Summer_Olympics_%E2%80%93_Men%27s_12_hour_race)
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

## Exercise

Ok, time to train! let's see if you can use the sort tool to answer the following questions:

> ### {% icon question %} Exercise: Reverse the sort
>
> Which athlete comes *last by alphabet*, in the *most recent* Olympics?
>
> > ### {% icon solution %} Answer
> >
> >  . We do this by repeating the previous sort (on year and then name), but changing the order to *descending* for both, to get the answer to the top of the file.
> >
> {: .solution}
{: .question}


> ### {% icon question %} Exercise: Sorting by age
>
> 1. What is the oldest age of a competing athelete?
> 2. Of all athletes of this age, which comes first alphabeticall?
>
> > ### {% icon solution %} Answer
> >
> >  . We do this by repeating the previous sort (on year and then name), but changing the order to *descending* for both, to get the answer to the top of the file.
> >
> {: .solution}
{: .question}



# Filter by column

show only winter olympics


# Find and Replace




Tip: Star your favorite tools

# Transpose


# Join columns from two different files


> ### {% icon hands_on %} Hands-on: Get data
>
> 1. {% tool [Import](upload1) %} the file `country-information.tsv` via link
>
>    ```
>    https://zenodo.org/record/6638036/files/country-information.tsv
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>>
{: .hands_on}

IOC ID to country name

Add Continent column

# Concatenating files

add years 2016 and onwards?

# Rearrange columns

# Remove columns

# Group on column

# Search in file

grep?

# Split file

Winter and Summer

# Concatenate file

Join winter and summer again

# Compute expressions

datamash, group

count number of countries

get number of athletes per country

average number of medals per athlete per country

country with highest numer of medals

sport with the tallest athletes

# Visualisation

# Exercises (with answers)

This section provides a number of exercises that require you to combine one or more of the techniques you learned in this tutorial. This is a great place to start if you want to practice your data manipulation skills!
