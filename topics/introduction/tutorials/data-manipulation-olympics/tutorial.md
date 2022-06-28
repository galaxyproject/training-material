---
layout: tutorial_hands_on

title: Data Manipulation Olympics
zenodo_link: 'https://zenodo.org/record/6772556'
tags:
- cyoa
questions:
- How can I do basic data manipulation in Galaxy?
- Which tools are available to convert, reformat, filter, sort etc my text-based data?
objectives:
- Familiarize yourself with data manipulation tools in Galaxy
- Perform basic text manipulation tasks in Galaxy
- Become comfortable converting text-based files in a variety of ways.
time_estimation: 1h
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


<!-- Note to contributors: feel free to add sections here to include additional visualisation options. Make sure each section is independent of each other, i.e. each section should start with the olympics.tsv file. Also make sure to include many exercises (with answers) for your section -->


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


# Cheatsheet

Here is an overview table of the different data manipulations in this tutorial, with links to the tools in Galaxy.

If you've opened this tutorial via the {% icon level %} icon in Galaxy (top menu bar), you can click on the tool names in the last column to quickly open them in Galaxy and start using them on your own!


| Operation              | Description            | Galaxy Tool    |
|------------------------|------------------------|----------------|
| Convert format         | Change the file format | {% icon galaxy-pencil%} Edit attributes |
| Word count             | Count the number of lines, words and characters in a file | {% tool [Line/Word/Character count](wc_gnu) %} |
| Sort on a column       | Change the order of the rows based on values in one or more columns | {% tool [Sort](sort1) %} |
| Filter                 | Remove rows from a file based on values in one or more columns | {% tool [Filter](Filter1) %} |
| Replace Text | in a specific column | {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/1.1.3) %}|

TIP: If you find yourself frequently using the same tool often but struggle to find it in the long list of tools, you can **star** your favourite tools in Galaxy!

{% snippet faqs/galaxy/tools_favorite.md %}


In this tutorial, these tools are explained in more detail, and we provide some exercises for you to practice.


# Background

In this tutorial, we will use as our dataset a table with results from the Olympics, from the games in Athens in 1896 until Tokyo in 2020. The objective is to familiarize you with a large number of the most important data manipulation tools in Galaxy. Much like the Olympics, there are many different disciplines (types of operations), and for each operation there are often multiple techniques (tools) available to athletes (data analysts, you) that are great for achieving the goal.


![image of olympic rings, logo and two atheletes around the words "Data Analysis Olympics"](./images/cover.jpg)


We will show you many of these commonly needed data manipulation operations, and some examples of how to perform them in Galaxy. We also provide many exercises so that you can train your skills and become a data manipulation Olympian!


> ### {% icon tip %} Also check out the Data Visualisation Olympics tutorial!
>
> There is a related tutorial that uses this dataset to perform a range of different visualisations in Galaxy.
>
> Interested? Check it out [here]({% link topics/introduction/tutorials/data-visualisation-olympics/tutorial.md %})
>
{: .comment}


# Upload Data

Before we can do any manipulation, we will need some data. Let's upload our table with olympics results now.

> ### {% icon hands_on %} Hands-on: Get data
>
> 1. {% tool [Import](upload1) %} the file `olympics.tsv` via link
>
>    ```
>    {{page.zenodo_link}}/files/olympics.tsv
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
>    > 3. How many lines are in the file? (Hint: {% tool [Line/Word/Character count](wc_gnu) %})
>    > 4. How many columns?
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. When you expand the dataset, you will see `format: tabular`, this is another term for a tab-separated (`tsv`) file.
>    > > 2. Each row represents an athlete's participation in an event. If an athlete competes in multiple events, there is a line for each event.
>    > > 3. 234,523. Look at the expanded view in the history, this tells us there are ~250,000 rows in the dataset. We can get the exact number using the  {% tool [Line/Word/Character count](wc_gnu) %} tool.
>    > > 4. There are 17 columns in this file. There are multiple ways to find this answer:
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


## About this dataset

The data was [obtained](https://github.com/UOSCS/Olympic_Athletes) from [Olympedia](https://www.olympedia.org/). The file `olympics.tsv` contains
234,522 rows and 17 columns. Each row corresponds to an individual athlete competing in an individual Olympic event. The columns are:

- **athlete_id** - Unique number for each athlete
- **name** - Athlete's name
- **sex** - M or F
- **birth_year** - 4-digit number
- **birth_day** - e.g. 24 July
- **birth_place** - town and/or country
- **height** - In centimeters
- **weight** - In kilograms
- **team** - Team name
- **noc** - National Olympic Committee 3-letter code
- **games** - Year and season
- **year** - Integer
- **season** - Summer or Winter
- **city** - Host city
- **sport** - Sport
- **event** - Event
- **medal** - Gold, Silver, Bronze, or NA

We will use this dataset to practice our data manipulation skills in Galaxy.


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
>    > 1. What do you notice?
>    > 2. Why are some values in quotes?
>    > 3. Scroll down a bit. Why are some rows displayed differently than others?
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. Galaxy does not display the table as nicely as before.
>    > >    This is because Galaxy is optimized to work with `tsv` files. For most rows you now see commas separating the different columns.
>    > > 2. If the data in a column contains a comma (e.g. in this file we have events such as `swimming 5,000 meters`), we put the value in quotes to signifiy that that comma is part of the data, not a column delimiter.
>    > > 3. This is a bit of a quirk in Galaxy, but when no columns contain a `,` as part of the value, Galaxy displays it as a table, but when columns contain a comma (and are in quotation marks), it displays it as text. Galaxy is optimized for `tsv` files, so if you have `csv` data, it is best to convert the datatype to `tabular` before you start your analysis.
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
>    > > it displays the data as-is, no longer hiding any commas.
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

We have a lot of data in this file, but it is ordered by the athlete ID number, which is a somewhat arbitrary and meaningless number. But we can sort the rows in this file to something more convenient, for example alphabetically by name of the athelete, or chronologically by year of the Olympics.

> ### {% icon hands_on %} Hands-on: Sort table based on a column
>
> 1. We will sort the file in chronological order based on the year of the Olympic games
>
>    > ### {% icon question %} Questions
>    >
>    > 1. Which column contains the year?
>    > 2. Do we want ascending or descending order if we want the oldest games at the top?
>    > 3. What should we do with the very first row (the one with the header names?)
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. Column 12
>    > > 2. The file should be sorted in ascending (increasing) order
>    > > 3. The header line should always stay on top, so we want to ignore it when we sort the file.
>    > >
>    > {: .solution}
>    {: .question}
>
>
> 2. {% tool [Sort](sort1) %} with the following parameters:
>    - {% icon param-file %} *"Sort dataset"*: `olympics.tsv`
>    - {% icon param-select %} *"on column"*: `Column 12`
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
>    > > 1. Karakatsanis. Who competed in a Shooting event 1896 Summer Olympics in Athens.
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
> 2. {% icon galaxy-refresh %} **Rerun** the sort tool with the following parameters:
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
>    > > 1. A. Grigoriadis. He competed in the 500 meters freestyle swiming event.
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

## Exercises

Ok, time to train! let's see if you can use the sort tool to answer the following questions:

> ### {% icon question %} Exercise: Reverse the sort
>
> Which athlete comes *last by alphabet*, in the *most recent* Olympics?
>
> > ### {% icon solution %} Answer
> >
> > `Žolt Peto` who competed in tabletennis at the 2020 Summer Olympics in Tokyo.
> > We do this by repeating the previous sort (on year and then name), but changing the order to *descending* for both, to get the answer to the top of the file.
> >
> {: .solution}
{: .question}


> ### {% icon question %} Exercise: sort by height
>
> 1. What is the height of the tallest competing athelete? Which athlete(s) are of this height?
> 2. What is the shortest?
> 3. Who was the tallest athlete from the most recent Olympics? How tall were they?
>
> > ### {% icon solution %} Hints
> >
> > 1. Height is listed in column 7. Since this is a number, we need *numerical sort*, and becaue we want the tallest on top, we will need to sort in *descending* (decreasing) order.
> > 2. Rerun the tool for step 1, but change the order to *ascending*
> > 3. First sort by year (descending), then by height (descending)
> >
> {: .solution}
>
> > ### {% icon solution %} Answers
> >
> >  1. Adam Sandurski from poland is the tallest athlete in the file, at 214 cm tall.
> >  2. Here we don't get a clear answer. That is because not all athletes have height data available, and blank fields are being sorted to the top. We can filter out rows with empty values to get our anwer (see Filter section to learn how to do this). For now we cannot answer this question with just the sort tool. Some times multiple tools are required to perform such tasks. The [exercise section at the end of this tutorial](#exercises-with-answers) has many exercises that require a combination of different tools.
> >  3. Gennaro Di Mauro, 210 cm.
> >
> {: .solution}
{: .question}



# Filter by column

This file contains a lot of data, we may only be interested in a subset of this data. For example, we may only want to look at one particular Olympics, or one particular sport. In such cases we can filter the dataset. This will remove any rows that are not of interest to us (i.e. that don't meet certain criteria).


> ### {% icon hands_on %} Hands-on: Filter table based on a column
>
> We will filter the file to show only winter Olympics
>
> 1. Look at the `olympics.tsv` file and answer the following questions
>
>    > ### {% icon question %} Questions
>    >
>    > 1. Which column contains this information?
>    > 2. Which values can this column have? (make sure to notice capitalisation, 'Winter' is not the same as 'winter' to these tools)
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. Column 13, the column with the *season* header
>    > > 2. The values can be `Summer` or `Winter`
>    > >
>    > {: .solution}
>    {: .question}
>
>
> 2. {% tool [**Filter** data on any column using simple expressions](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `olympics.tsv`
>    - {% icon param-select %} *"With the following condition"*: `c13=='Winter'`
>    - {% icon param-text %} *"Number of header lines to skip"*: `1`
>
> 3. {% icon galaxy-eye %} **View** the filtered file.
>
>    > ### {% icon question %} Question
>    >
>    > How many lines are in this file? (Hint: expand the dataset in your history or use {% tool [Line/Word/Character count](wc_gnu) %} )
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. 44,681 (this is including the header line)
>    > >
>    > {: .solution}
>    {: .question}
>
> 4. Repeat the step for the Summer olympics
>
>    > ### {% icon question %} Question
>    >
>    > 1. How many lines do you expect in the this file?
>    > 2. How many lines are in this file? Were you right?
>    >
>    > > ### {% icon solution %} Hints
>    > >
>    > > 1. Use the {% tool [Line/Word/Character count](wc_gnu) %} to find the number of lines in the `olympics.tsv` file and subtract the number of rows in the Winter olympics file
>    > > 2. Be careful to consider whether these counts include the header line of the file or not
>    > >
>    > {: .solution}
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. The original file has 234,523 lines, and the Winter olympics had 44,681 lines. So we would expect 234,523 - 44,681 = 189,842 rows of data. Since we have subtracted the header line in this equation as well, we expect the Summer olypmics file to have 1 more line that this, so 189,843 total lines.
>    > > 2. 189,843. If you were off by one or two lines, it may have been that you counted the header lines double
>    > >
>    > > It is always useful to take a moment to think about the expected outcome, this makes it easier to spot mistakes and will save you time in the long run.
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

## Exercises

Ok, time to train! let's see if you can use the filter tool to answer the following questions:





# Unique



# Find and Replace


# Transpose


# Join columns from two different files


> ### {% icon hands_on %} Hands-on: Get data
>
> 1. {% tool [Import](upload1) %} the file `country-information.tsv` via link
>
>    ```
>    {{page.zenodo_link}}/files/country-information.tsv
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

- TODO: shortest athlete (couldnt ansswer with only sort tool)


> ### {% icon question %} Exercise: sort by height
>
> 1.
> 2.
> 3.
>
> > ### {% icon solution %} Hints
> >
> > 1.
> > 2.
> > 3.
> >
> {: .solution}
>
> > ### {% icon solution %} Answers
> >
> >  1.
> >  2.
> >  3.
> {: .solution}
{: .question}

