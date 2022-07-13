---
layout: tutorial_hands_on

title: Data Manipulation Olympics
zenodo_link: 'https://zenodo.org/record/6803028'
tags:
- cyoa
questions:
- How can I do basic data manipulation in Galaxy?
- Which tools are available to convert, reformat, filter, sort etc my text-based data?
objectives:
- Familiarize yourself with data manipulation tools in Galaxy
- Perform basic text manipulation tasks in Galaxy
- Become comfortable converting text-based files in a variety of ways.
- Reason about the expected outcome of tools
time_estimation: 1h
key_points:
- There are a lot of tools available in Galaxy for performing basic data manipulation tasks
- Bacis data manipulation is often needed between steps in a larger scientific analysis in order to connect outputs from one tool to input of another.
- There are often multiple ways/tools to achieve the same end result
- Having a basic understanding of data manipulation tools will make it easier to do exploratory data analysis
- Always read the help text of the tool before using it to get a full understanding of its workings
- Always try to formulate the output you are expecting from a tool. This will make it easier to spot mistakes as soon as possible.
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


<!--
Note to contributors: feel free to add sections here to include additional data manipulation options.
Make sure each section is independent of each other, i.e. each section should start with the olympics.tsv file.
Also make sure to include many exercises (with answers) for your section!
-->


<!-- set up some variables to easily update tool versions throughout tutorial, since most tools are used many times %} -->
{% assign version_datamash="toolshed.g2.bx.psu.edu/repos/iuc/datamash_ops/datamash_ops/1.1.0+galaxy2" %}
{% assign version_column_maker="toolshed.g2.bx.psu.edu/repos/devteam/column_maker/Add_a_column1/1.6" %}
{% assign version_replace_text_column="toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regexColumn1/1.0.1" %}
{% assign version_replace_text_line="toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.1" %}}
{% assign version_cat="cat1" %}
{% assign version_remove_beginning="Remove beginning1" %}
{% assign version_join="join1" %}
{% assign version_remove_columns_by_header="toolshed.g2.bx.psu.edu/repos/iuc/column_remove_by_header/column_remove_by_header/0.0.1" %}
{% assign version_cut_columns="Cut1" %}
{% assign version_paste="Paste1" %}
{% assign version_split="toolshed.g2.bx.psu.edu/repos/bgruening/split_file_on_column/tp_split_on_column/0.4" %}


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


| Operation              | Description                        | Galaxy Tool    |
|------------------------|------------------------------------|----------------|
| Convert format         | Change the file format             | {% icon galaxy-pencil%} Edit attributes |
| Word count             | Count the number of lines, words and characters in a file | {% tool [Line/Word/Character count](wc_gnu) %} |
| Sort on a column       | Change the order of the rows based on values in one or more columns | {% tool [Sort](sort1) %} |
| Filter                 | Remove rows based on values in one or more columns | {% tool [Filter](Filter1) %}|
| Counting               | Count occurrences of values in a column   | {% tool [**Count**](Count1) %}, {% tool [Datamash]({{version_datamash}}) %}|
| Group on a column      | And perform simple operations (count, mean, min, max etc) | {% tool [**Group**](Grouping1) %}, {% tool [Datamash]({{version_datamash}}) %} |
| Compute an expression  | Over each row                      | {% tool [Compute]({{version_column_maker}}) %} |
| Find and Replace       | in a specific column               | {% tool [Column Regex Find and Replace]({{version_replace_text_column}}) %}|
| Find and Replace       | on every line                      | {% tool [Regex Find and Replace]({{version_replace_text_line}}) %}|
| Join two Datasets      | side by side on a specified field  | {% tool [Join two Datasets]({{version_join}}) %} |
| Concatenate datasets   | one after the other                | {% tool [Concatenate datasets]({{version_cat}}) %} |
| Remove Beginning       | Good for removing header lines     | {% tool [Remove beginning of a file]({{version_remove_beginning}}) %}
| Cut Columns            | By header name                     | {% tool [Remove columns by heading]({{version_remove_columns_by_header}}) %}|
| Cut Columns            | By column number                   | {% tool [Cut columns from a table]({{version_cut_columns}}) %}|
| Paste                  | Two files side by side             | {% tool [Paste]({{version_paste}}) %} |
| Split file             | Based on values of a column        | {% tool [Split]({{version_split}}) %} |

**TIP: Adding tools to your Favourites:** If you find yourself frequently using the same tool often but struggle to find it back in the long list of tools, you can **star** your favourite tools in Galaxy!

{% snippet faqs/galaxy/tools_favorite.md %}


In this tutorial, these tools are explained in more detail, and we provide some exercises for you to practice.


# Background

In this tutorial, we will use as our dataset a table with results from the Olympics, from the games in Athens in 1896 until Tokyo in 2020. The objective is to familiarize you with a large number of the most important data manipulation tools in Galaxy. Much like the Olympics, there are many different disciplines (types of operations), and for each operation there are often multiple techniques (tools) available to athletes (data analysts, you) that are great for achieving the goal.


![image of olympic rings, logo and two athletes around the words "Data Analysis Olympics"](./images/cover.jpg)


We will show you many of these commonly needed data manipulation operations, and some examples of how to perform them in Galaxy. We also provide many exercises so that you can train your skills and become a data manipulation Olympian!


# Upload Data

Before we can do any manipulation, we will need some data. Let's upload our table with Olympics results now.

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
- **height** - In centimeters (or `NA` if data not known)
- **weight** - In kilograms (or `NA` if data not known)
- **team** - Team name
- **noc** - National Olympic Committee 3-letter code
- **games** - Year and season
- **year** - Integer
- **season** - Summer or Winter
- **city** - Host city
- **sport** - Sport
- **event** - Event
- **medal** - Gold, Silver, Bronze (or `NA` if no medal was won)

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



# Sorting

We have a lot of data in this file, but it is ordered by the athlete ID number, which is a somewhat arbitrary and meaningless number. But we can sort the rows in this file to something more convenient, for example alphabetically by name of the athlete, or chronologically by year of the Olympics.

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

Ok, time to train! Let's see if you can use the sort tool to answer the following questions:

> ### {% icon question %} Exercise: Reverse the sort
>
> Which athlete comes *last by alphabet*, in the *most recent* Olympics?
>
> > ### {% icon solution %} Answer
> >
> > `Žolt Peto` who competed in table tennis at the 2020 Summer Olympics in Tokyo.
> > We do this by repeating the previous sort (on year and then name), but changing the order to *descending* for both, to get the answer to the top of the file.
> >
> {: .solution}
{: .question}


> ### {% icon question %} Exercise: sort by height
>
> 1. What is the height of the tallest competing athlete? Which athlete(s) are of this height?
> 2. What is the shortest?
> 3. Who was the tallest athlete from the most recent Olympics? How tall were they?
>
> > ### {% icon solution %} Hints
> >
> > 1. Height is listed in column 7. Since this is a number, we need *numerical sort*, and because we want the tallest on top, we will need to sort in *descending* (decreasing) order.
> > 2. Rerun the tool for step 1, but change the order to *ascending*
> > 3. First sort by year (descending), then by height (descending)
> >
> {: .solution}
>
> > ### {% icon solution %} Answers
> >
> >  1. Adam Sandurski from Poland is the tallest athlete in the file, at 214 cm tall.
> >  2. Here we don't get a clear answer. That is because not all athletes have height data available, and blank fields are being sorted to the top. We can filter out rows with empty values to get our answer (see Filter section to learn how to do this). For now we cannot answer this question with just the sort tool. Some times multiple tools are required to perform such tasks. The [exercise section at the end of this tutorial](#exercises-with-answers) has many exercises that require a combination of different tools.
> >  3. Gennaro Di Mauro, 210 cm.
> >
> {: .solution}
{: .question}


# Counting

A common operation we might want to perform on tables of data, is simple counting. How many times does a certain value appear? For our dataset for instance, we might want to know how many countries participated in each Olympics, how many women, etc; any column that has categorical data that we can count. The tool {% tool [**Count** occurrences of each record](Count1) %} does exactly this.


> ### {% icon hands_on %} Hands-on: Count occurrences of values in columns
>
> Let's start by simply counting how many different Olympic Games we have in our dataset, and how many times it appears (so how many participations there were each year)
>
> 1. {% tool [**Count** occurrences of each record](Count1)%} with the following parameters
>    - *"from dataset"*: `olympics.tsv`
>    - *"Count occurrences of values in column(s)"*: `Column 11`
>
> 2. {% icon galaxy-eye %} **View** the results.
>
>    > ### {% icon question %} Question
>    >
>    > 1. How many different Olympic games are in our file?
>    > 2. Which Olympic games had the most participations? (Tip: set the parameter *"How should the results be sorted?"* to `most common values first`)
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. 52. there are 53 lines in the resulting file, with one line containing the value of the column header (`games`).
>    > > 2. 1996 Summer Olympics. (10501 participations)
>    > >
>    > > The resulting file looks something like:
>    > >
>    > > ```
>    > > 615	1896 Summer Olympics
>    > > 2503	1900 Summer Olympics
>    > > 2643	1904 Summer Olympics
>    > > 3213	1908 Summer Olympics
>    > > 4610	1912 Summer Olympics
>    > > 3448	1920 Summer Olympics
>    > > 5242	1924 Summer Olympics
>    > > 358	1924 Winter Olympics
>    > > 4493	1928 Summer Olympics
>    > > ...
>    > > ```
>    > {: .solution}
>    {: .question}
>
> You may have noticed that we could have selected multiple columns in this tool. This lets us count on combinations of columns.
> Let's try counting the number of men and women in each olympic games.
>
> 3. {% tool [**Count** occurrences of each record](Count1)%} with the following parameters
>    - *"from dataset"*: `olympics.tsv`
>    - *"Count occurrences of values in column(s)"*: `Column 11, Column 3`
>
> 4. {% icon galaxy-eye %} **View** the results.
>
>    > ### {% icon question %} Question
>    >
>    > You see the resulting file has a line for every combination of the two columns (games and sex), providing the count for each.
>    >
>    > 1. How many women were in the first Olympic games?
>    > 2. Which Olympic games had the most women participants? (Tip: set the parameter *"How should the results be sorted?"* to `most common values first`)
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. 2. (note that we cannot be sure if this is two different women, or 1 woman participating twice).
>    > > 2. 2020 Summer Olympics (4652)
>    > >
>    > > The resulting file looks somtthing like:
>    > >
>    > > ```
>    > > 615	1896 Summer Olympics
>    > > 2503	1900 Summer Olympics
>    > > 2643	1904 Summer Olympics
>    > > 3213	1908 Summer Olympics
>    > > 4610	1912 Summer Olympics
>    > > 3448	1920 Summer Olympics
>    > > 5242	1924 Summer Olympics
>    > > 358	1924 Winter Olympics
>    > > 4493	1928 Summer Olympics
>    > > ...
>    > > ```
>    > {: .solution}
>    {: .question}
>
>
{: .hands_on}


Let's say we wanted to know how many different sports there were in each Olympics. If we used the counting tool above, we would get a results file for each combination of sport and olympics, with the number of lines (participations) of each. But we don't really care about the number of lines that have this combination, just the total number of unique sports in each games.

To answer these types of questions we can use a slightly more advanced tool, called {% tool [Datamash]({{version_datamash}}) %}. This tool can do a lot more than counting, but here we will show you how to use it to count in this way.


> ### {% icon hands_on %} Hands-on: Counting number of different sports per Olympic Games
>
> We will now determine how many different sport there were in each of the different Olympics
>
> 1. {% tool [Datamash]({{version_datamash}}) %} with the following parameters:
>    - *"Input tabular dataset"*: `olympics.tsv`
>    - *"Group by fields"*: `11` (the `games` column)
>    - *"Input file has a header line"*: `yes`
>    - *"Sort Input"*: `yes`
>    - *"Operation to perform on each group"*: `Count Unique values`
>      - *"On Column"*: `Column: 15` (the `sport` column)
>
> 2. {% icon galaxy-eye %} **View** the results.
>
>    > ### {% icon question %} Question
>    >
>    > 1. Why did we tell the tool to sort the file?
>    > 2. How many sport were in the first Olympics? How many in the latest?
>    > 3. Which Olympics had the most different sports?
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. This is for technical reasons, the tool assumes a sorted file, so it will restart the counter every time it encounters a value in a column that is different from the previous one.
>    > > 2. 10 and 38.
>    > > 3. The 2020 Summer Olympics had the most different sports (38)
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}


Datamash can do a lot more than counting,  and we will showcase some of these other operatins in the [Grouping](#grouping) section.

## Exercises

Ok, let's practice!

> ### {% icon question %} Exercise: Number of participations per country
>
> 1. Which country has had the most participations in the Olympics?
> 2. How many countries participated in the first Olympics? How many in the last?
>
> > ### {% icon solution %} Hints
> >
> > 1. Since we are counting participations (rows), we can use the simple {% tool [Count](Count1)%} tool
> > 2. Since we are counting a bit more complex question, we need the {% tool [Datamash]({{version_datamash}}) %} tool
> >
> {: .solution}
>
> > ### {% icon solution %} Answers
> >
> >  1. The United States with 17,286 participations
> >  2. 15 and 250.
> >
> {: .solution}
>
> > ### {% icon solution %} Full Solutions
> >
> >  1. {% tool [Count](Count1) %} with the following parameters:
> >     - *"from dataset"*: `olympics.tsv`
> >     - *"Count occurrences of values in column(s)"*: `Column 9` (the `team` column)
> >     - *"How should the results be sorted?"*: `With the most common values first`
> >
> >     This gives an output like:
> >     ```
> >     17286	United States
> >     11700	France
> >     10230	Great Britain
> >     8898	Italy
> >     7988	Canada
> >     ```
> >
> >  2. {% tool [Datamash]({{version_datamash}}) %} with the following parameters:
> >   - *"Input tabular dataset"*: `olympics.tsv`
> >   - *"Group by fields"*: `11` (the `games` column)
> >   - *"Input file has a header line"*: `yes`
> >   - *"Sort Input"*: `yes`
> >   - *"Operation to perform on each group"*: `Count Unique values`
> >     - *"On Column"*: `Column: 9` (the `sport` column)
> >
> >   This gives an output like:
> >   ```
> >   1896 Summer Olympics	15
> >   1900 Summer Olympics	29
> >   1904 Summer Olympics	13
> >   ...
> >   2016 Summer Olympics	280
> >   2018 Winter Olympics	108
> >   2020 Summer Olympics	250
> >   ```
> >
> {: .solution}
{: .question}



# Filtering

This file contains a lot of data, but we may only be interested in a subset of this data. For example, we may only want to look at one particular Olympics, or one particular sport. In such cases we can filter the dataset. This will create a new dataset, removing any rows that are not of interest to us (i.e. that don't meet the criteria we provide).


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
>    > > ### {% icon solution %} Answers
>    > >
>    > > 1. Column 13, the column with the *season* header
>    > > 2. The values can be `Summer` or `Winter`
>    > >
>    > {: .solution}
>    {: .question}
>
> 2. Open the tool {% tool [**Filter** data on any column using simple expressions](Filter1) %}
>    - read the help text at the bottom of the tool carefully
>    - get familiar with how to indicate columns and write expressions for filtering
>
>    > ### {% icon question %} Questions
>    >
>    > 1. How would you write the expressions for the following conditions:
>    >    1. column 6 must be contain 'Yes'
>    >    2. column 13 must be smaller than 75
>    >    3. column 7 cannot be 'NA'
>    >    4. column 2 cannot be empty
>    >
>    > 2. It is also possible to combine multiple conditions, using `and`, `or`, `not` and parentheses
>    >    How would you write expressions for the following filtering conditions:
>    >    1. column 5 is larger than 2 or smaller than -2
>    >    2. column 5 is larger than 2 and smaller than 10
>    >    3. the sum of columns 4 and 6 is greater than or equal to 25,000
>    >
>    > > ### {% icon solution %} Answers
>    > >
>    > > 1. The answers are:
>    > >    1. `c6=='Yes'`
>    > >    2. `c13<75`
>    > >    3. `c7!='NA'`
>    > >    4. `c2!=''`
>    > >
>    > > 2. The answers are:
>    > >    1. `c5>2 or c5<-2`
>    > >    2. `c5>2 and c5<10`
>    > >    3. `c4+c6 >= 25000`
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
>    > > 44,681 (this is including the header line)
>    > >
>    > {: .solution}
>    {: .question}
>
> 4. Repeat the step for the Summer Olympics
>
>    > ### {% icon question %} Questions
>    >
>    > 1. How many lines do you expect in the this file?
>    > 2. How many lines are in this file? Were you right?
>    >
>    > > ### {% icon solution %} Hints
>    > >
>    > > 1. Use the {% tool [Line/Word/Character count](wc_gnu) %} to find the number of lines in the `olympics.tsv` file and subtract the number of rows in the Winter Olympics file
>    > > 2. Be careful to consider whether these counts include the header line of the file or not
>    > >
>    > {: .solution}
>    >
>    > > ### {% icon solution %} Answers
>    > >
>    > > 1. The original file has 234,523 lines, and the Winter Olympics had 44,681 lines. So we would expect 234,523 - 44,681 = 189,842 rows of data. Since we have subtracted the header line in this equation as well, we expect the Summer Olympics file to have 1 more line that this, so 189,843 total lines.
>    > > 2. 189,843. If you were off by one or two lines, it may have been that you counted the header lines double
>    > > <br>
>    > > It is always useful to take a moment to think about the expected outcome, this makes it easier to spot mistakes and will save you time in the long run.
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

## Exercises

Ok, time to train! let's see if you can use the filter tool to answer the following questions:


> ### {% icon question %} Exercise: Medal winners
>
> 1. How many gold medals were handed out?
> 2. How many total medals?
> 3. How many medals were handed out during the 2018 Olympics?
> 4. How many medals were won by individuals with a height between 170 and 180 cm?
> 5. How many gold medals were won by individuals shorter than 160cm or taller than 190?
>
> > ### {% icon solution %} Hints
> >
> > - Column 17 contains information about medals
> > - The possible values are `Gold`, `Silver`, `Bronze`, and `` (empty).
> > - Expand the output or use the tool {% tool [Line/Word/Character count](wc_gnu) %} to see the number of lines in the file
> > - Don't forget that the output (and line count) may include the header line
> > - Do not use quotes on number columns (e.g. year)
> > - You may need parentheses for complex conditions
> >
> {: .solution}
>
> > ### {% icon solution %} Answers
> >
> >  1. 8,110.  Expression: `c17=='Gold'`
> >  2. 24,633. Expression: `c17=='Gold' or c17=='Silver' or c17=='Bronze'`, or `c17!=''`
> >  3. 131.    Expression: `c17=='Gold' and c12==2018` (note: do not use quotes around `2018`, as it is a numerical value)
> >  4. 8,086   Expression: `c17!='' and c7>=170 and c7<=180`
> >  5. 812     Expression: `c17=='Gold' and (c7<160 or c7>190)` (note: parentheses are important here)
> >
> {: .solution}
{: .question}





# Grouping

Often we may want to group rows based on a value in a column, and perform some operation on the resulting rows. For example we would like to group the olympics data by one value (e.g. year, country, sport), and determine some value for each group (e.g. number of medals won, average age of athletes, etc).

In the [counting](#counting) section of this tutorial we show how to get answers that require a count (e.g. number of medals won), but sometimes we want to do something more complex, like calculating the average height of athletes in a group, say per country or per sport. This section will show some example of these types of questions.

We can use the {% tool [Datamash]({{version_datamash}}) %} tool for this purpose.

> ### {% icon hands_on %} Hands-on: Number of athletes per Olympics
>
> We would like to answer the following question: *How tall was the tallest athlete of each sport?*
>
> 1. Open the {% tool [Datamash]({{version_datamash}}) %} tool and read the help section at the bottom
>
>    > ### {% icon question %} Question
>    >
>    > Which settings do you think we need to provide to answer our question?
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > - *"Group by fields"*: We want to group by team (Column 15).
>    > > - *"Sort"*: `Yes`. This may not be obvious, but because our file is currently not sorted by our chosen group (sport), we need to tell the tool to do this.
>    > > - *"Skip NA or NaN values"*: since we do have NA values for athletes for whom height data is unknown, we should set this to `Yes`
>    > > - Our file has a header line, so we should indicate this as well
>    > > - *"Operation"*: we want to determine the **Maximum** height (Column 7)
>    > >
>    > {: .solution}
>    {: .question}
>
> 2. {% tool [Datamash]({{version_datamash}}) %} with the following parameters:
>    - *"Input tabular dataset"*: `olympics.tsv`
>    - *"Group by fields"*: `15`
>    - *"Sort input"*: `yes`
>    - *"Input file has a header line"*: `yes`
>    - *"Skip NA or NaN values"*: `yes`
>    - *"Operation to perform on each group"*: `Maximum`
>      - *"On Column"*: `Column: 7`
>
> 3. {% icon galaxy-eye %} **View** the results.
>
>    > ### {% icon question %} Question
>    >
>    > 1. How tall was the tallest athlete in basketball? And what about karate?
>    > 2. Why do some sports have a value of `inf`?
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. Basketball's tallest athlete was 192cm. For Karate it is 163.
>    > > 2. Our dataset had quite a number of `NA` (unknown) values in the height column, especially for the earlier Olympics. For sports that had only NA values, there is no maximum so the tool outputs `inf` instead.
>    > {: .solution}
>    {: .question}
>
{: .hands_on}


## Grouping on multiple columns

You may have noticed that we could also provide multiple columns to group on. If we do this, we can compute values for combinations of groups, such as sex and sport, to find e.g. the tallest woman in basketball or the shortest man per Olympics. There are also many more options for the computation we perform, so perhaps we are more interested not in the tallest athlete, but the average height. Let's perform some of these slightly more advanced queries now.


> ### {% icon hands_on %} Hands-on: Tallest man and woman per sport
>
> The question we would like to answer here, is what is the average height for men and women per sport?
>
> 1.  Open the {% tool [Datamash]({{version_datamash}}) %} tool
>     - Which parameters do you think we need?
>     - Refer to the help section at the bottom of the tool page if you need more information
>
> 2. {% tool [Datamash]({{version_datamash}}) %} with the following parameters:
>    - *"Input tabular dataset"*: `olympics.tsv`
>    - *"Group by fields"*: `15,3` (The sports column, and the sex column)
>    - *"Sort input"*: `yes`
>    - *"Input file has a header line"*: `yes`
>    - *"Print header line"*: `yes`
>    - *"Skip NA or NaN values"*: `yes`
>    - *"Operation to perform on each group"*: `Mean`
>      - *"On Column"*: `Column: 7`
>
> 3. {% icon galaxy-eye %} **View** the results.
>    - Notice the header line in this output that we requested with the *"Print header line parameter"*. Adding this line will help you remember which columns you grouped on and which computation you performed. In our case it was obvious, but if you have a dataset with multiple columns with similar values, this can be useful
>    - See if you can answer the following questions based on the output file.
>
>    > ### {% icon question %} Question
>    >
>    > 1. What is the average height of women participating in archery?
>    > 2. What is the average height of men participating in [ballooning](https://en.wikipedia.org/wiki/Ballooning_at_the_1900_Summer_Olympics)?
>    > 3. Why do some values have `nan` instead of a height?
>    > 4. Why do some sports not have a value for one of the sexes?
>    > 5. Can you find a sport where women were taller than the men? (Hint: it starts with the letter A)
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. 167.25677031093 cm
>    > > 2. 170 cm
>    > > 3. If none of the rows in the group had height data available, it will output `nan` (not a number) instead. This is most common for sports that were only featured a few times in the early years of the Olympics.
>    > > 4. Sports such as artistic swimming only exist for women, so no M appears in the data for that group, so there simply is no row for the mean height of men doing artistic swimming in our output.
>    > > 5. [Art Competitions](https://en.wikipedia.org/wiki/Art_competitions_at_the_Summer_Olympics)
>    > >
>    > > If all went well, your output file should look something like:
>    > >
>    > > ```
>    > > GroupBy(sport)	     GroupBy(sex)  mean(height)
>    > > Aeronautics         M             nan
>    > > Alpine Skiing       F             167.38324708926
>    > > Alpine Skiing       M             178.18747142204
>    > > Alpinism            M             nan
>    > > Archery             F             167.25677031093
>    > > Archery             M             178.5865470852
>    > > Art Competitions    F             175.33333333333
>    > > Art Competitions    M             173.97260273973
>    > > Artistic Gymnastics F             156.15316901408
>    > > ```
>    > {: .solution}
>    {: .question}
>
{: .hands_on}


## Exercises

> ### {% icon question %} Exercise: Grouping and computing
>
> 1. How tall is the shortest woman Badminton player to win a gold medal?
> 2. What is the average weight of athletes from Denmark (DEN) in the 2020 Olympics?
>
> > ### {% icon solution %} Hints
> >
> > 1. We need to group on 3 columns: medals, sport and sex (note: the order you provide the columns determines the order they are listed in in the output)
> > 2. We need to group on 2 columns: country and year, then compute the average (mean) over column 8 (weight)
> >
> {: .solution}
>
> > ### {% icon solution %} Answers
> >
> >  1. 161 cm.
> >  2.
> >
> {: .solution}
>
> > ### {% icon solution %} Full Solutions
> >
> > 1. {% tool [Datamash]({{version_datamash}}) %} with the following parameters:
> >    - *"Input tabular dataset"*: `olympics.tsv`
> >    - *"Group by fields"*: `17,15,3` (The medal, sports, and the sex columns)
> >    - *"Sort input"*: `yes`
> >    - *"Input file has a header line"*: `yes`
> >    - *"Print header line"*: `yes`
> >    - *"Skip NA or NaN values"*: `yes`
> >    - *"Operation to perform on each group"*: `Minimum`
> >      - *"On Column"*: `Column: 7`
> >
> >    This will give an output like below, scroll down to find the gold medalists, then badminton, then F (if you used a different order in the Group by fields parameter, this file may look a bit different, but will still provide the same information)
> >
> >    ```
> >    GroupBy(medal)	GroupBy(sport)	GroupBy(sex)	min(height)
> >    Bronze	Alpine Skiing	F	156
> >    Bronze	Alpine Skiing	M	167
> >    Bronze	Archery	F	155
> >    Bronze	Archery	M	166
> >    Bronze	Art Competitions	F	-inf
> >    Bronze	Art Competitions	M	172
> >    Bronze	Artistic Gymnastics	F	136
> >    ```
> >
> > 2. {% tool [Datamash]({{version_datamash}}) %} with the following parameters:
> >    - *"Input tabular dataset"*: `olympics.tsv`
> >    - *"Group by fields"*: `12,9` (year and team columns)
> >    - *"Sort input"*: `yes`
> >    - *"Input file has a header line"*: `yes`
> >    - *"Print header line"*: `yes`
> >    - *"Skip NA or NaN values"*: `yes`
> >    - *"Operation to perform on each group"*: `Minimum`
> >      - *"On Column"*: `Column: 8`
> >
> >    This will give an output like below:
> >
> >    ```
> >    GroupBy(medal)	GroupBy(sport)	GroupBy(sex)	min(height)
> >    Bronze	Alpine Skiing	F	156
> >    Bronze	Alpine Skiing	M	167
> >    Bronze	Archery	F	155
> >    Bronze	Archery	M	166
> >    Bronze	Art Competitions	F	-inf
> >    Bronze	Art Competitions	M	172
> >    Bronze	Artistic Gymnastics	F	136
> >    ```
> {: .solution}
>
{: .question}


# Computing

Sometimes we want to use the data in our column to compute a new value, and add that to the table. For instance, for our dataset we could caluclate athtletes BMI (using height and weight columns), or their age at time of participations (from year of birth and year of the Olymics). By adding these computed values as a new colum to our datset, we can more easily query the dataset for these values. We can do these types of operations using the {% tool [Compute - an expression on every row]({{version_column_maker}}) %} tool.

As an example, let's calculate the age of each athlete at the time of participation, and add this as a new column to our dataset.


> ### {% icon hands_on %} Hands-on: Compute age of athletes
>
> 1. Open the {% tool [Compute an expression on every row]({{version_column_maker}}) %} tool.
>    - read the help text at the bottom of the tool
>    - what parameters do you think we need to use?
>
> 2. {% tool [Compute an expression on every row]({{version_column_maker}}) %} with the following parameters:.
>    - *"Add expression"*: `c12-c4` (year column minus the year_of_birth column)
>    - *"As a new column to"*: `olympics.tsv`
>    - *"Round result?"*: `Yes`
>    - *"Input has a header line with column names?"*: `Yes`
>    - *"The new column name"*: `age`
>
> 3. {% icon galaxy-eye %} **View** the results.
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What changed in our file?
>    > 2. How old was Arnaud Boetsch during his Olympic tennis participation?
>    >
>    > > ### {% icon solution %} Answers
>    > >
>    > > 1. A new `age` column was added to the end of our file, the value is the age of the athlete in years at time of the olympics.
>    > > 2. Arnaud Boetsch is listed on the first two lines, they turned 27 the year of their olympics.
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

This was a simple computation, but much more complex mathematical expressions can be computed with this tool. See the help section at the bottom of the tool for a list of all supported operations. In the exercise below, we will compute the BMI for each athlete as an example.


## Exercises

BMI stands for Body Mass Index, is a metric to provide a very crude measure of how healthy your weight is. The formula to compute BMI is:

$$ BMI = weight / (height^2) $$

(with weight in kilograms and height in meters).


Let's use the {% tool [Compute]({{version_column_maker}}) %} tool to compute this data for all athletes and add it as a new column.




> ### {% icon question %} Exercise: Calculating BMI
>
> 1. How would you express this calculation in the tool?
>    - Remember that our height is in cm, and the formula expects height in meters
>    - For each tabular file, Galaxy will try to determine whether a row is numerical or not. While height and weight are numbers, we also have a lot of "NA" values here, which Galaxy sees as a word, and may not automatically interpret the column as numerical, so to be safe, we will manually tell the tool the column are numeric,
>      - do this using e.g. `int(c3)` (`int` stands for integer, if you had a column with decimal numbers, you would say `float(c3)`)
>
> 2. What is the BMI for Arnaud Boetsch?
>
> 3. Why does the output have fewer lines than the input?
>
> > ### {% icon solution %} Hints
> >
> > - division is `/` and multiplication is ` * ` .
> > - The tool does not recognize the `^` operation as exponentiation. You can use `height * height` or `pow(height,2)`
> > - Parentheses may be required.
> > - Use `int(column)` to tell the tool the columns contain numbers
> >   - e.g. for column 3 + column you would use `int(c3) + int(c4)` in the expression
> >   - this is only needed because some of our rows have "NA" in this columns, so Galaxy is not sure if it is a number column, a string (word) column, or a mixed column.
> >
> {: .solution}
>
> > ### {% icon solution %} Answers
> >
> >  1. `int(c8)/(int(c7)*int(c7))*10000` (other variations are possible)
> >  2. 22.69
> >  3. The tool only outputs lines for which it was able to perform the computation, so any lines which had `NA` in the height or weight column are skipped. You could use the [join operation](#joining-files) to re-obtain the missing lines, see also one of the exercises at the end of this tutorial.
> {: .solution}
>
{: .question}



# Find and Replace

Often you may need to change the contents of a file a bit to fit the expectations of an analysis tool. For instance, our file uses `NA` for missing values, but other conventions included leaving the cell empty, or using `NaN` (not a number) instead. Or, when working with chromosomal data, you may need to add or remove the `chr` prefix from a column before using it as input to a certain tool. In such situations, we can find all occurrences of a certain pattern in our file, and replace it with another value.

If we want to perform such a replacement on a single column in our data, we can use the {% tool [Column Regex Find and Replace]({{version_replace_text_column}}) %} tool, or if we want to do a replacement on a line by line basis (e.g. if our data isn't tabular and therefore doesn't have columns), we can use the {% tool [Regex Find and Replace]({{version_replace_text_line}}) %} tool. Both of these tools are similar, in that they use *regular expressions* (or *regex*) to define a pattern to search for. Regular expressions are a standardized way to describe patterns, and while they can seem quite complex at first, with just a few of the basics and a bit of practice, you can perform powerful operations on large datasets with ease.

A few of the basics of regular expression, plus some links to further resources are given in the box below:

{% snippet faqs/galaxy/analysis_regular_expressions.md %}

Let's start with a simple example:

> ### {% icon hands_on %} Hands-on: Find and Replace
>
> Our file uses a mix of `Athina` and `Athens` to indicate the Capital City of Greece in the `city` column.
> Let's standardize this by replacing occurrences of `Athina` with `Athens`.
>
> 1. Open {% tool [Column Regex Find and Replace]({{version_replace_text_column}}) %}
>    - Read the help text at the bottom, what settings do you think we need to use?
>    - Read the Regular expressions 101 FAQ below
>
>      {% snippet faqs/galaxy/analysis_regular_expressions.md %}
>
> 2. {% tool [Column Regex Find and Replace]({{version_replace_text_column}}) %} with the following parameters:
>    - {% icon param-file %} *"Select cells from"*: `olympics.tsv`
>    - {% icon param-select %}*"using column"*: `Column 14`
>    - {% icon param-repeat %} *"Check"*
>      - *"Find Regex"*: `Athina`
>      - *"Replacement"*: `Athens`
>
> 3. {% icon galaxy-eye %} **View** the results.
>    - Look at the file before and after. Athlete 7 (Patrick Chila) near the top of the `olympics.tsv` file, had a value of Athina in the city column. Verify that it has been changed to Athens.
>
>    > ### {% icon question %} Question
>    >
>    > Why did we use the column replace tool, and not the general replace tool?
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > It is safer to use the column replace tool in our case. We know the city name only occurs in one of the columns.
>    > > If we had used the line replace tool, it would have replaced all occurrences of `Athina`, which may have unforeseen consequences (e.g. maybe somebody was named Athina, in that case we don't want to replace it)
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

This was rather simple example, so let's try a few more examples with slightly more complex expressions.


## Exercises

You may have noticed that our file has a lot of missing data. Especially for the earlier years, things like height, weight and birthday of athletes was not registered, or simply not known. In some columns you see these missing values have been replaced with an `NA` (not available) value. In other columns (for example birth place), the cells have simply been left empty.

Different tools may expect different ways of handling missing data. So you may have to change your missing data from empty to `NA`, `NaN`, or something else, between analysis steps


> ### {% icon hands_on %} Hands-on: Fill empty cells
>
> We will now replace empty cells in the `birth_place` column, do use `NA` instead
>
> 1. Open {% tool [Column Regex Find and Replace]({{version_replace_text_column}}) %}
>    - Read the help text at the bottom, what settings do you think we need to use?
>    - Read the Regular expressions 101 FAQ below
>
>      {% snippet faqs/galaxy/analysis_regular_expressions.md %}
>
>    > ### {% icon question %} Question
>    >
>    > 1. Should we use the column or line replace tool?
>    > 2. What is the expression for an empty line? (hint: see regular expressions 101 box above)
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. Since we only want to replace in one column, we will use the column replace tool
>    > > 2. `^$` indicates an empty line (`^` indicated the beginning, and `$` the end). The value in the column is treated as the line, since we are looking only in 1 column.
>    > >
>    > {: .solution}
>    {: .question}
>
> 2. {% tool [Column Regex Find and Replace]({{version_replace_text_column}}) %} with the following parameters:
>    - {% icon param-file %} *"Select cells from"*: `olympics.tsv`
>    - {% icon param-select %}*"using column"*: `Column 6`
>    - {% icon param-repeat %} *"Check"*
>      - *"Find Regex"*: `^$`
>      - *"Replacement"*: `NA`
>
> 3. {% icon galaxy-eye %} **View** the results.
>    - Look at the file before replacement, find a line without a value in the `birth_place` column. Verify that it has been replaced in the resulting file.
>
{: .hands_on}

Let's do another example, this one using capture groups.

Look at the `birth_day` column. It has values in a format like `12 December`. Suppose we have a tool that expects this data to be in the reverse format, `December 12`. The file is too big to change this manually in every column. But with regular expression tools we can make this replacement easily

> ### {% icon hands_on %} Hands-on: Reverse Birthday Format
>
> We will now change the format in birthday column from `day month` to `month day`
>
> 1. Open {% tool [Column Regex Find and Replace]({{version_replace_text_column}}) %}
>    - Read the help text at the bottom, what settings do you think we need to use?
>    - Read the Regular expressions 101 FAQ below
>
>      {% snippet faqs/galaxy/analysis_regular_expressions.md %}
>
>    > ### {% icon question %} Questions
>    >
>    > 1. How do we match on the birthday format? How strict/exact shoule we be here?
>    > 2. How do we captures both the day and the month?
>    > 3. How do we refer to the values we captured (for the replacement value)
>    >
>    > > ### {% icon solution %} Hints
>    > >
>    > > 1. Birthday is one or more digits, followed by a space, followed by one or more letters.
>    > > 2. Remember that you can capture values using parentheses `(..)`
>    > > 3. `\1` will be the variable containing the value in the first capture, `\2` the second, etc
>    > >
>    > {: .solution}
>    >
>    > > ### {% icon solution %} Answers
>    > >
>    > > 1. There are multiple solutions here, depending on how strict you want to be
>    > >    - `\d+ ([a-zA-Z]+)` (not strict, would also match on `142 Septober`
>    > >    - `[0-9]{1,2} (January|February|March|April|May|June|July|August|September|October|November|December)` (this is much more strict, only matches on 1 or 2 digits, followed by a space, followed by one of the months. But this would still match on 42 April, and may miss it if one month names didn't start with a capital letter)
>    > >    - `[123]?[0-9] [a-zA-Z]+` this will only allow dates below 40
>    > >
>    > >    There are different ways to express this, and there is no one perfect solution. If you know your data is clean, and you do not have to worry about values like 42 Septober,
>    > >    then you can be less strict. If your data is less clean, you may be more worried about capturing things that aren't valid birthdays and might want to be a bit stricter.
>    > >    It all depends on your data and use case.
>    > >
>    > > 2. `([\d]{1,2}) ([a-zA-Z]+)` captures both the day and the month
>    > > 3. We can use `\1` to refer to the day we captured, and `\2` for the month in our replacement.
>    > >
>    > {: .solution}
>    {: .question}
>
> 2. {% tool [Column Regex Finda and Replace]({{version_replace_text_column}}) %} with the following parameters:
>    - {% icon param-file %} *"Select cells from"*: `olympics.tsv`
>    - {% icon param-select %}*"using column"*: `Column 5`
>    - {% icon param-repeat %} *"Check"*
>      - *"Find Regex"*: `(\d{1,2}) ([a-zA-Z]+)` (or your own variation)
>      - *"Replacement"*: `\2 \1` (first month, then space, then day)
>
> 3. {% icon galaxy-eye %} **View** the results.
>    - Did it work? If not, no shame, it often takes some trial and error to get your expression right. Click the rerun {% icon galaxy-refresh %} button on the tool,
>      and tweak your expression until it works.
>
{: .hands_on}




# Removing Columns

We can remove columns from a table using either {% tool [Remove columns by heading]({{version_remove_columns_by_header}}) %} if your table has a header line, or {% tool [Cut columns from a table]({{version_cut_columns}}) %} if it does not (in this case we just indicate columns by their number). These tools can also be used to change the order of columns in your file

To do the reverse, adding one or more columns, we can use the {% tool [Paste]({{version_paste}}) %} tool. This assumes we have the same number of rows in both files, already in the correct order. It is a very "dumb" tool that simple combines the two files next to each other.


> ### {% icon hands_on %} Hands-on: Remove columns
>
> Suppose we want to simplify our file a bit. All we want is file with 4 columns: athlete name, sport, olympic games, and medals.
>
> 1. Open the {% tool [Remove columns by heading]({{version_remove_columns_by_header}}) %} tool and read the help text at the bottom.
>    - Which settings do you think we need?
>
> 2. {% tool [Remove columns by heading]({{version_remove_columns_by_header}}) %} with the following parameters:
>    - *"Tabular File"*: `olympics.tsv`
>    - *"Header name"*: `name`
>    - {% icon param-repeat %} *"Header name"*: `sport`
>    - {% icon param-repeat %} *"Header name"*: `games`
>    - {% icon param-repeat %} *"Header name"*: `medals`
>    - *"Keep named columns"*: `Yes`
>
> 3. {% icon galaxy-eye %} **View** the results.
>
>    > ### {% icon question %} Question
>    >
>    > 1. How many rows and columns do you expect the output to have?
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. We expect the same number of rows as the original dataset, but now only the 4 columns we requested to keep.
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}


Notice that during this step, we also changed the order of the columns. This tool can also be used to re-arrange columns, if you supply all column names but in a different order.



## Exercises


> ### {% icon question %} Exercise: Removing Columns
>
> 1. Create a file with 4 columns: name, sport, games, medal (same task as the previous hands-on), but use the {% tool [Cut columns from a table]({{version_cut_columns}}) %} tool instead.
> 2. Create a file that is similar to `olympics.tsv`, but without the first column (athlete_id column)
> 3. You want to keep all the columns of `olympics.tsv`, but change the order so that `sport` and `event` come right after the athlete's name.
>
> > ### {% icon solution %} Hints
> >
> > 1. We need to determine the column numbers and provide these to the tool, rather than the column headers
> > 2. Think about what the *"Keep named colums"* parameter does to simplify the settings here
> > 3. Which of the two tools would be easier to use? (this can be a personal preference, but think about whether you would rather provide column names or numbers)
> >
> {: .solution}
>
> > ### {% icon solution %} Full Solutions
> >
> > 1. {% tool [Cut columns from a table]({{version_cut_columns}}) %} using the following parameters:
> >    - *"Cut Columns"*: `c2,c15,c11,c17`
> >    - *"Delimited By"*: `TAB`
> >    - *"From"*: `olympics.tsv`
> >
> > 2. {% tool [Remove columns by heading]({{version_remove_columns_by_header}}) %} with the following parameters:
> >    - *"Tabular File"*: `olympics.tsv`
> >    - *"Header name"*: `athlete_id`
> >    - *"Keep named columns"*: `No`
> >
> >    **Note:** the *"Keep named columns"* parameter determines wheter we keep or remove the columns we specified.
> >    You could have obtained the same result by supplying all column names except the first one, and selecting
> >    *"Keep named columns"*: `No`, but that would have been a lot more work.
> >
> > 3. {% tool [Cut columns from a table]({{version_cut_columns}}) %} using the following parameters:
> >    - *"Cut Columns"*: `c1,c2,c15,c16,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c17`
> >    - *"Delimited By"*: `TAB`
> >    - *"From"*: `olympics.tsv`
> >
> >    **Note:** you can also use the {% tool [Remove columns by heading]({{version_remove_columns_by_header}}) %}, it just requires a bit more typing,
> >    but on the other hand it is also a bit less error-prone (i.e. it is easier to mix up column numbers than column names).
> >
> {: .solution}
>
{: .question}




# Joining Files

This file contains a lot of information, but we may want to add more information. For example if we had a file with information about each country (population, capital city, etc), we could join that information with our olympics data, to get a single file with all information in every row.

For example, if we would like to be able to group by continent, to e.g. count athletes, medals etc per continent, we will have to add a `continent` column to our file. To do this we would need a second file that maps each country to the corresponding continent. This is what we will do in the next hands-on section.

We obtained country information data from [DataHub](https://datahub.io/core/country-codes). More information about this file can be found in the description there. It has 56 columns with a wide variety of data about each country (from country codes, to capital city, languages spoken, etc)

> ### {% icon hands_on %} Hands-on: Get data
>
> 1. {% tool [Import](upload1) %} the file `country-information.tsv` via link
>
>    ```
>    {{page.zenodo_link}}/files/country-information.tsv
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 2. {% icon galaxy-eye %} **View** file.
>
>    > ### {% icon question %} Question
>    >
>    > 1. How many columns does this file have?
>    > 2. Which column(s) in this file are the same as in the `olympics.tsv` file?
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. The country information file has 56 columns (see [DataHub](https://datahub.io/core/country-codes) for more details).
>    > > 2. Both files have a `NOC` column with the 3-letter country code (`NOC` stands for National Olympic Committee). We can use this column to join the appropriate country data to each row in our `olympics.tsv` file.
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

We would now like to take our Olympics dataset as the basis, and add columns to every row of this file with some information about the country. In order to join, we will need to have one column that is shared between the two files, on which we can match. The `NOC` column is perfect for this because it is a defined standard. Both files also contain a column with the country name in it, which is also a possible candidate to use for joining, but because it is less standardised, it is safer to use the NOC column. For example, if one file uses "Netherlands", while the other uses "The Netherlands" to indicate the same country, the joining will fail for these rows. So always make sure the columns you join on are compatible!

> ### {% icon hands_on %} Hands-on: Joining country information to the Olympics data.
>
> 1. Open {% tool [Join two Datasets side by side on a specified field]({{version_join}}) %}, and read the help text at the bottom
>    - Which settings do you think we need to use?
>
> 2. {% tool [Join two Datasets side by side on a specified field]({{version_join}}) %} with the following parameters:
>    - *"Join"*: `olympics.tsv`
>    - *"using column"*: `Column 10` (the `noc` column)
>    - *"with"*: `country-information.tsv`
>    - *"and column"*: `Column 2` (the `NOC` column)
>    - *"Keep the header lines?"*: `Yes`
>
> 3. {% icon galaxy-eye %} **View** the results.
>
>    > ### {% icon question %} Question
>    >
>    > 1. What do you expect the output to look like? Were you right?
>    > 2. How many columns are in the resulting file?
>    > 3. What is a possible downside to this approach?
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. All the columns from the country information file are added to the end of each row of our olympics dataset
>    > > 2. Our olympics datset had 17 columns, the country information file has 56 columns. Therefore we have 17+56=73 columns columns in our resulting file.
>    > > 3. There is a lot of data duplication in this file now. The exact same country information is added to every line of every athlete from a certain country.
>    > >    This means much larger file size, and more usage of your quota.
>    > >    If you do not need all these columns, it could save you a lot of space to remove unneeded columns from the `country-information.tsv` file, before joining.
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}



# Concatenating

Concatenation of two files simple means putting the contents of the two files together, one after the other. Our dataset was created in 2021, but since then we've had another Olympic event, the 2022 Winter Olympics in Beijing. If we have the same data for this latest Olympics, we could simply add the rows from the 2022 games to our current file with data, in order to create a single file with all data from 1896 to 2022.

First, let's get this data for the 2022 Olympics

> ### {% icon hands_on %} Hands-on: Get 2022 Olympics data
>
> 1. {% tool [Import](upload1) %} the file `country-information.tsv` via link
>
>    ```
>    {{page.zenodo_link}}/files/olympics_2022.tsv
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 2. View {% icon galaxy-eye %} the new dataset, does it have the same structure as our original `olympics.tsv` file?
>
>    > ### {% icon question %} Question
>    >
>    > 1. Does the new file have the same structure?
>    > 2. Can we simply add the lines of the new files to the end of our existing olympics dataset?
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. Yes, this file has all the same columns, in the same order, so concatenation should be relatively straightforward.
>    > > 2. If we simply put the contents of this file after our existing dataset, we will have a second header line in the middle of our data rows. It is best to remove the header line from the second dataset after we have verified it is compatible. This way we will only add the real data rows to our dataset.
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}


Since this new dataset has the exact same structure (number and order of columns), we can simple add the lines from this file to the end of our existing `olympics.tsv` file.


> ### {% icon hands_on %} Hands-on: Concatenate the two files
>
> First we need to remove the header line from the 2022 file
>
> 1. {% tool [Remove beginning of a file]({{version_remove_beginning}}) with the following parameters:
>    - *"Remove first"*: `1`
>    - *"from"*: `olympics_2022.tsv`
>
> Now we can perform the concatenation:
>
> 2. {% tool [Concatenate datasets tail-to-head]({{version_cat}}) %} with the following parameters:
>    - *"Concatenate Datasets"*: `olympics.tsv` (this file will be first)
>    - {% icon param-repeat%} *"Insert Dataset"*: `olympics_2022.tsv`
>
> 3. {% icon galaxy-eye %} **View** the results.
>
>    > ### {% icon question %} Question
>    >
>    > 1. How many lines do you expect in the new file? Were you correct? (Hint: use {% tool [Line/Word/Character count](wc_gnu) %} to count lines)
>    > 2. Where are the lines of the 2022 Olympics?
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > 1. The `olympics.tsv` file had 234,523 lines, and the `olympics_2022.tsv` file had 4076 lines. Both of these numbers include a header line, which we removed for the second file, so we expect our concatenated file to contain 234,523 + 4076 - 1 = 238,598 lines (which it does).
>    > > 2. The new file has the entire contents of `olympics.tsv` at the beginning of the file, followed by the contents of the `olympics_2022.tsv` file at the end.
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

Now this only works so simply because our two datasets had the same structure; the same number of columns, in the same order. If your data comes from different sources, you may have to do some additional data manipulation before you can concatenate, e.g. to make sure the columns match, or how each file deals with missing data (empty cells, `NA`, `NaN` or something else).


# Splitting Files

This dataset contains a lot of data, we might want to split the data into one file per Olympic games, or have one file for all the winter games, and another file for all the summer games. In these situtations, where we want to create new, smaller files, based on the values in a column, we can use the tool {% tool [Split file according to the values of a column]({{version_split}}) %}.

Tip: use this tool only if you want a separate file for **all values** in a column. For example, if you want a file for each Olympic games, use the split tool, but if you just want a file for one particular Olympics (say Tokyo 2020), you could simply use the [Filter](#filtering) tools to extract the data from the full dataset, without making files for all the other Olympics as well.


> ### {% icon hands_on %} Hands-on: Splitting the dataset by Olympics
>
> We would like to create a separate file for each Olympic event.
>
> 1. Open the tool {% tool [Split file according to the values of a column]({{version_split}}) %}, and read the help text at the bottom
>    - which settings do you think we need to use?
>
> 2. {% tool [Split file according to the values of a column]({{version_split}}) %} with the following parameters:
>    - *"File to select"*: `olympics.tsv`
>    - *"on column"*: `Column 11` (the `games` column)
>    - *"Include the header in all splitted files?"*: `Yes`
>
> 3. {% icon galaxy-eye %} **View** the results.
>
>    > ### {% icon question %} Questions
>    >
>    > 1. How many different files do we get?
>    > 2. How many lines are in the file for the 1908 Summer Olympics?
>    >
>    > > ### {% icon solution %} Answers
>    > >
>    > > 1. We get 1 collection with 52 files in it, one per Olympic games:
>    > >    ![screenshot of the output collection of the split tool](./images/split-file.png)
>    > > 2. 3214 lines
>    > {: .solution}
>    {: .question}
>
{: .hands_on}


## Exercises

Let's practice this a bit more, see if you can use the split file to answer the following questions:

> ### {% icon question %} Exercise: sort by height
>
> 1. Create two files, one for Summer Olympics, one for Winter Olympics. Which has more lines?
> 2. Split the file by sport, how many sports have there been at the Olympics?
> 3. Split the file by medal, would you expect the output files to be equal sizes?
>
> > ### {% icon solution %} Hints
> >
> > 1. Split on the `season` column (Column 13)
> > 2. Split on the `sport` column (Column 15)
> > 3. Split on the `medal` column (Column 17)
> >
> {: .solution}
>
> > ### {% icon solution %} Answers
> >
> >  1. The summer olympics file has more lines (~180,000) than the winter olympics file (~45,000) (the first winter Olympics wasn't held until 1924)
> >  2. Each file in the resulting collection represents a sport, there are 91 datasets in the collection, representing 91 different sports in the history of the Olympics.
> >  3. The files for gold silver and bronze are roughly equal, but not exactly (think of situations like shared second place, then 2 silver medals are handed out and no bronze medal)
> {: .solution}
>
{: .question}



# Conclusion

These tools and operations covered in the tutorial are just a few examples of some of the most common operations. There are many more tools available in Galaxy that perform other data manipulations. We encourage you to look around the `General Text Tools` section in the Galaxy toolbox (left panel) to find and explore more data manipulation tools. The more comfortable you are performing these kinds of steps, the more you can get out of Galaxy!


# Exercises: Putting it all together!

This section provides a number of exercises that require you to combine two or more of the techniques you learned in this tutorial. This is a great way to practice your data manipulation skills. Full solutions are provided for every exercise (i.e. all tools and settings), but for many of these exercises there will be multiple solutions, so if you obtained the same results in a different way, that is correct too!

- TODO: shortest athlete (couldnt answer with only sort tool, also need to filter out empty lines)
-  TODO: number of unique athletes who received silver medals during the 2018 olympics

- TODO: first find and replace in weight column to replace values such as "68-70" with just one number (or median) orso (or use compute to add new weight column), then can use datamash to find mean weight of athletes per country or sport orso

- TODO: calculate BMI, lines with NA will be missing from output, join back with original file.

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
>
> > ### {% icon solution %} Full solution
> >
> >  1.
> >  2.
> >  3.
> {: .solution}
>
{: .question}



> ### {% icon hands_on %} Hands-on: Sort table based on a column
>
>
> 1. {% icon galaxy-eye %} **View** the results.
>
>    > ### {% icon question %} Question
>    >
>    >
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > >
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

