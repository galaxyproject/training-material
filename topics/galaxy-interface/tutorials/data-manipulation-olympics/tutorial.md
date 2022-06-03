---
layout: tutorial_hands_on

title: Data Manipulation Olympics
tags:
- workflows
zenodo_link: ''
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
level: Beginner

subtopic: analyse
---


# Introduction
{:.no_toc}

Scientific analyses often consist of a number of tools that run one after the other, to go from raw data to biological insight. Between these specialized tools, simple data manipulation steps are often needed as a kind of "glue" between tools. For example, the output of tool A may produce a file that contains all the information needed as input for tool B, but tool B expects the columns in a different order.

Galaxy has a large collection of tools to perform such basic data manipulation tasks, and becoming familiar with these will allow to perform your analysis more easily.

In this tutorial, we will use historical dataset about the Olympics to familiarize you with a large number of the most important data manipulation tools. Much like the Olympics, there are many different disciplines (types of operations), and for each operation there are often multiple ways to accomplish the goal. We will show you all of these different operations and some examples on how to perform them, and provide many exercises so you can train your skills and become a data manipulation Olympian!


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}



# Upload Data

> ### {% icon hands_on %} Hands-on: Get data
>
> 1. {% tool [Import](upload1) %} the file `olympics_sports_almanac.csv` via link
>
>    ```
>    https://zenodo.org/record/6611769/files/olympics-sports-almanac.csv
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
{: .hands_on}


# File Format Conversion

CSV to TSV

# Find and Replace

remove quotes from TSV file

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
>    https://zenodo.org/record/6611769/files/country-information.tsv
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
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

