---
layout: tutorial_hands_on

title: Data Manipulation Olympics
tags:
- workflows
zenodo_link: 'https://zenodo.org/record/6627896'
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

In this tutorial, we will use as our dataset a table with results from the Olympics, from the games in Athens in 1986 until Rio 2016. The objective is to familiarize you with a large number of the most important data manipulation tools in Galaxy. Much like the Olympics, there are many different disciplines (types of operations), and for each operation there are often multiple techniques (tools) available to athletes (data analysts, you) that are great for achieving the goal.


![image of olympic rings, logo and two atheletes around the words "Data Analysis Olympics"](./images/cover.jpg)


We will show you many of these commonly needed data manipulation operations, and some examples of how to perform them in Galaxy. We also provide many exercises so that you can train your skills and become a data manipulation Olympian!


# Upload Data

Before we can do any manipulation, we will need some data. Let's upload our table with olympics results now.

> ### {% icon hands_on %} Hands-on: Get data
>
> 1. {% tool [Import](upload1) %} the file `olympics_sports_almanac.csv` via link
>
>    ```
>    https://zenodo.org/record/6627896/files/olympics-sports-almanac.csv
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 2. **View** {% icon galaxy-eye %} the dataset
>    - We uploaded a `csv` file, a comma-separated file. Galaxy displays this as a table, hiding the commas for convenience.
>
>    > ### {% icon question %} Question
>    >
>    > Why are most of the values in quotation marks?
>    >
>    > > ### {% icon solution %} Answer
>    > >
>    > > This is a CSV file, that means a comma-separated file. It is assumed that each comma in the file signifies the start of a new column.
>    > >
>    > > If the data in a column contains a comma (e.g. in this file we have events such as `swimming 5,000 meters`), we put the value in quotes to signifiy that that comma is part of the data, not a column delimiter.
>    > >
>    > > Often, csv files put all non-numeric values in quotes for this reason.
>    > {: .solution}
>    {: .question}
>
{: .hands_on}



# Choose your adventure!

This tutorial is structured a bit differently than most. **You do not have to do the steps in the order they are presented below.** Every section in this tutorial uses the dataset you just uploaded (the `olympics-sports-almanactsv` file) as input, so you can jump to any section in this tutorial right now if you have a particular data manipulation operation in mind you want to learn more about.




# File Format Conversion

Galaxy understands the `csv` format just fine, but most tools work best with tab-separated (tabular) files. This format uses TAB characters to signify where a new column begins instead of commas. Galaxy can do the conversion between these formats for you.


Hands-on: CSV to TSV




# Find and Replace

An added benefit of using tab-seperated files instead of comma-separated files is that it is much less likely that we need a tab character as part of the value in a column than a comma. Certainly in our dataset none of the columns contain a tab character as part of the value, so we can remove the quotation marks from our file



remove quotes from TSV file


Tip: Star your favorite tools

# Summary Statistics



# Sort by column

sort by date


# Filter by column

show only winter olympics


# Join columns from two different files


> ### {% icon hands_on %} Hands-on: Get data
>
> 1. {% tool [Import](upload1) %} the file `country-information.tsv` via link
>
>    ```
>    https://zenodo.org/record/6627896/files/country-information.tsv
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>>
{: .hands_on}

IOC ID to country name

Add Continent column

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

