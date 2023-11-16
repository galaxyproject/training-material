---
layout: tutorial_hands_on

title: 'Data Manipulation Olympics - JQ'
zenodo_link: 'https://zenodo.org/record/6803028'
tags:
- cyoa
- json
questions:
- How can I do basic data manipulation with JQ?
- Which functions are available to convert, reformat, filter, sort etc my JSON-based data?
objectives:
- Familiarize yourself with data manipulation filters in JQ
- Perform basic data manipulation tasks in JQ
- Become comfortable converting JSON-based files in a variety of ways.
- Reason about the expected outcome of tools
time_estimation: 1h
key_points:
- There are a lot of functions available in JQ for performing basic data manipulation tasks
- Basic data manipulation is often needed between steps in a larger scientific analysis in order to connect outputs from one tool to input of another.
- There are often multiple ways/tools to achieve the same end result
- Having a basic understanding of data manipulation tools will make it easier to do exploratory data analysis
- Always read the help text of the tool before using it to get a full understanding of its workings
- Always try to formulate the output you are expecting from a tool. This will make it easier to spot mistakes as soon as possible.
contributions:
  authorship:
    - shiltemann
    - hexylena
  funding:
    - gallantries
level: Introductory

subtopic: olympics
---



Scientific analyses often consist of a number of tools that run one after the other, in order to go from the raw data to scientific insight. Between these specialized tools, simple data manipulation steps are often needed as a kind of "glue" between tools. For example, the output of tool A may produce a file that contains all the information needed as input for tool B, but tool B expects the columns in a different order. Or in genomic data analysis, some tools expect chromosome X to be listed as `chrX`, while others simply expect `X`. In these situations, extra data manipulation steps are needed to prepare files for input to analysis tools.

<!--
Note to contributors: feel free to add sections here to include additional data manipulation options.
Make sure each section is independent of each other, i.e. each section should start with the olympics.tsv file.
Also make sure to include many exercises (with answers) for your section!
-->

Galaxy has a large collection of tools to perform such basic data manipulation tasks, and becoming familiar with these operations will allow to perform your analysis more easily in Galaxy (and outside).


> <agenda-title></agenda-title>
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

Operation                | Description                                                         | Galaxy Tool
------------------------ | ------------------------------------                                | ----------------
Compute an expression    | Over each row, add it as a new column                               | `.key1 * .key2`
Concatenate datasets     | one after the other                                                 | `cat *.json | jq`
Convert format           | Change the file format                                              | Often xsv, other tools.
Counting                 | Count the number of objects in an array                             | `length`
Cut Columns              | By header name                                                      | `.key`
Filter                   | Remove objects based on values                                      | `select(.key == "value")`
Find and Replace         | in a specific column                                                | `gsub("(?<a>...)", "xyz\\(.a)abc")`
Group on a column        | And perform simple operations (count, mean, min, max etc)           | `group_by(.key1, -.key2)`
Join two Datasets        | side by side on a specified field                                   | `[JOIN(INDEX(input[]; .right); .[]; .left; add)]`
Paste                    | Two files side by side                                              | `-s '[.[][]]'`
Remove Beginning         | Good for removing header lines                                      | `last`, `index`
Select First lines       | Good for finding top 10s or saving header lines                     | `first`, `index`
Sort on a column         | Change the order of the rows based on values in one or more columns | `sort`, `sort_by(.key)`
Split file               | Based on values of a column                                         | n/a (requires bash/awk/etc.)
Unique                   | Remove duplicate rows                                               | `unique`

In this tutorial, these functions are explained in more detail, and we provide some exercises for you to practice.

# Background

In this tutorial, we will use as our dataset a table with results from the Olympics, from the games in Athens in 1896 until Tokyo in 2020. The objective is to familiarize you with a large number of the most important data manipulation tools in Galaxy. Much like the Olympics, there are many different disciplines (types of operations), and for each operation there are often multiple techniques (tools) available to athletes (data analysts, you) that are great for achieving the goal.


![image of olympic rings, logo and two athletes around the words "Data Analysis Olympics"]({% link topics/introduction/tutorials/data-manipulation-olympics/images/cover.jpg %})


We will show you many of these commonly needed data manipulation operations, and some examples of how to perform them in Galaxy. We also provide many exercises so that you can train your skills and become a data manipulation Olympian!


# Download Data

Before we can do any manipulation, we will need some data. Let's download our table with Olympics results now.

We can use `jq` just by itself to pretty print the results:

```bash
cat olympics.json | jq -S . | head -n 20
```

> <question-title></question-title>
>
> 1. What is the format of the file?
> 2. How is it structured?
> 3. How many lines are in the file? (Hint: `wc`)
> 4. How many columns?
>
> > <solution-title></solution-title>
> >
> > 1. It is a JSON file, or a Javascript Object Notation
> > 2. It is structured as an array (we can tell by the leading `[`) with objects representing an athlete's participation in an event. If an athlete competes in multiple events, there is an entry for each event.
> > 3. `1`. The entire JSON file is on a single line, which makes it rather unpleasant to access via normal text editors. However it contains 234,522 records (`cat olympics.json | jq -c '.[]' | wc -l`)
> > 4. There are 17 columns in this file. There are multiple ways to find this answer:
> >    - Count the fields manually
> >    - Advanced JQ: `cat olympics.json| jq '.[0] | keys | length'`
> {: .solution}
{: .question}


## About this dataset

The data was [obtained](https://github.com/UOSCS/Olympic_Athletes) from [Olympedia](https://www.olympedia.org/). The file `olympics.json` contains
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

This tutorial is structured a bit differently than most. **You do not have to do the steps in the order they are presented below.** Every section in this tutorial uses the dataset you just uploaded (the `olympics.json` file) as input, so you can jump to any section in this tutorial right now if you have a particular data manipulation operation in mind you want to learn more about.

# File Format Conversion

The file we uploaded is a `.json` file. This stands for *JavaScript Object Notation*. This means that this is a file containing arrays (`[]`) and objects (`{}`). 

First we need to understand a bit about file structure, and how JQ can convert back and forth between different structures:

The `.` operator is the identity operator, it gives you back exactly the same structure

```bash
echo '[{"a": 1}, {"a": 2}]' | jq '.'
```

You can use `.[]` to loop over the items in an array:

```bash
echo '[{"a": 1}, {"a": 2}]' | jq '.[]'
```

And supply a number to access an individual array item:

```bash
echo '[{"a": 1}, {"a": 2}]' | jq '.[0]'
```

And `.<key>` to access a specific key

```bash
echo '[{"a": 1}, {"a": 2}]' | jq '.[].a'
```

But that's not all! Sometimes you may want the resulting output to be back in an array, rather than by line:

```bash
echo '[{"a": 1}, {"a": 2}]' | jq '[.[].a]'
```

By wrapping the entire expression in `[]` we say we want the output of that statement in an array. We can actually do the same to get each entry as a new object:

```bash
echo '[{"a": 1}, {"a": 2}]' | jq '{"results": .[].a}'
```

Or an object with an array:

```bash
echo '[{"a": 1}, {"a": 2}]' | jq '{"results": [.[].a]}'
```

Or an array of objects:

```bash
echo '[{"a": 1}, {"a": 2}]' | jq '[{"results": .[].a}]'
```

We'll need to use these reshaping abilities as we go, to get data in the right shape to be worked with. Sometimes you'll want an array of objects, or objects with items in their keys.

## Pipes

Lastly, pipes are incredibly important in JQ, we'll use them to combine multiple filters, or make our pipelines more readable:

```bash
echo '[{"a": 1}, {"a": 2}]' | jq '.'
```

Here we pipe the entire incoming object (JQ operates on streams of JSON objects, so, there can be multiple), to the identity operator again but this time with `[]` to loop over every object in that item. (Really this is the same as `.[]` but we use a pipe for demonstration.)

```bash
echo '[{"a": 1}, {"a": 2}]' | jq '. | .[]'
```

We can then pass that stream of individual items (no longer a single JSON object) on to another filter, one extracting an `.a`:

```bash
echo '[{"a": 1}, {"a": 2}]' | jq '. | .[] | .a'
```


# Sorting

We have a lot of data in this file, but it is ordered by the athlete ID number, which is a somewhat arbitrary and meaningless number. But we can sort the rows in this file to something more convenient, for example alphabetically by name of the athlete, or chronologically by year of the Olympics.

In `jq` we can use the `sort` or `sort_by` function. From the documentation:

```bash
echo '[8,3,null,6]' | jq '.sort'
```

`sort_by` however lets us sort multiple entries in an array, by a child key:

```bash
echo '[{"foo":4, "bar":10}, {"foo":3, "bar":100}, {"foo":2, "bar":1}]' | jq '.sort_by(.foo)'
echo '[{"foo":4, "bar":10}, {"foo":3, "bar":100}, {"foo":2, "bar":1}]' | jq '.sort_by(.bar)'
```

So let's sort our file in chronological order, based on the year of the Olympic games:

> <question-title></question-title>
>
> 1. Which key contains the year?
> 2. Do we want ascending or descending order if we want the oldest games at the top?
>
> > <solution-title></solution-title>
> >
> > 1. `.year`
> > 2. The file should be sorted in ascending (increasing) order
> >
> {: .solution}
{: .question}

```bash
cat olympics.json | jq -c 'sort_by(.year) | .[]'
```

> <tip-title>What does `.[]` do?</tip-title>
> We use this to first take the identity operator, taking exactly the array that was passed to us, and then `[]` lets us loop over the array, returning a stream of multiple individual objects. The `-c` flag lets us see that more easily, with one object per line.
{: .tip}

If we wanted to do it in reverse, we could just use the `reverse` function, or just use `-` prefixed to our sort factor.

```bash
cat olympics.json | jq -c 'sort_by(.year) | reverse | .[]'
cat olympics.json | jq -c 'sort_by(-.year) | .[]'
```

> <question-title></question-title>
>
> 1. Write a JQ command to access only the first entry.
> 2. Which athlete is listed at the top of the file now?
>
> > <solution-title></solution-title>
> > 1. We can simply access the 0th element:
> >    ```bash
> >    cat olympics.json | jq -c 'sort_by(.year) | .[0]'
> >    ```
> > 1. J. Defert. Who competed in a Tennis event 1896 Summer Olympics in Athens.
> >
> {: .solution}
{: .question}

**Rename** the resulting file to somthing meaningful (e.g. `olympics-chronological.json`)

This is great, but maybe it would make more sense to sort alphabetically by athlete name *within each year*.

## Sort on multiple columns at once

So we want to sort twice, first by year, an then within each year, we sort again alphabetically by name. The `sort_by` function can do this too!

We will sort the file in chronological order based on the year of the Olympic games

```bash
cat olympics.json | jq -c 'sort_by(.year, .name) | .[]'
```

> <question-title></question-title>
>
> Which athlete is listed at the top now? Which discipline (sport) did they compete in?
>
> > <solution-title></solution-title>
> >
> > 1. A. Grigoriadis. He competed in the 500 meters freestyle swimming event.
> >
> {: .solution}
{: .question}


## Exercises

Ok, time to train! Let's see if you can use the sort tool to answer the following questions:

> <question-title noprefix>Exercise: Reverse the sort</question-title>
>
> Which athlete comes *last by alphabet*, in the *most recent* Olympics?
>
> > <solution-title></solution-title>
> >
> > `Žolt Peto` who competed in table tennis at the 2020 Summer Olympics in Tokyo.
> > <br>
> > We do this by repeating the previous sort (on year and then name), but changing the order to *descending* for both, to get the answer to the top of the file.
> >
> {: .solution}
{: .question}


> <question-title noprefix>Exercise: sort by height</question-title>
>
> 1. What is the height of the tallest competing athlete? Which athlete(s) are of this height?
> 2. What is the shortest?
> 3. Who was the tallest athlete from the most recent Olympics? How tall were they?
>
> > <tip-title>Removing null values</tip-title>
> > This will be covered more during the Filtering section, but for now simply use this filter:
> > ```bash
> > cat olympics.json | jq '[.[] | select(.height != null)] | ...your sort_by here...'
> > ```
> {: .tip}
>
> > <solution-title noprefix>Hints</solution-title>
> >
> > 1. We can use `.height`, and because we want the tallest on top, we will need to sort in *descending* (decreasing) order. Unfortunately you might discover there are null values.
> > 2. Rerun the same query as step 1, but change the order to *ascending*
> > 3. First sort by year (descending), then by height (descending)
> >
> {: .solution}
>
> > <solution-title noprefix>Answers</solution-title>
> >
> >  1. Adam Sandurski from Poland is the tallest athlete in the file, at 214 cm tall.
> >  2. Lyton Mphande from Seol is the shortest at 127 cm.
> >  3. Gennaro Di Mauro, 210 cm. (2020 Summer Olympics in Tokyo)
> >
> {: .solution}
>
> > <solution-title noprefix>Full Solutions</solution-title>
> > 1. `jq -c '[.[] | select(.height != null)] | sort_by(-.height) | .[0]' < olympics.json`
> > 2. `jq -c '[.[] | select(.height != null)] | sort_by(.height) | .[0]' < olympics.json`
> > 3. `jq -c '[.[] | select(.height != null)] | sort_by(-.year, -.height) | .[0]' < olympics.json`
> {: .solution}
>
{: .question}


# Filtering

This file contains a lot of data, but we may only be interested in a subset of this data. For example, we may only want to look at one particular Olympics, or one particular sport. In such cases we can filter the dataset. This will create a new dataset, removing any rows that are not of interest to us (i.e. that don't meet the criteria we provide).


We will filter the file to show only winter Olympics
Look at the `olympics.json` file and answer the following questions

> <question-title></question-title>
>
> 1. Which key contains this information?
> 2. Which values can this column have? (make sure to notice capitalisation, 'Winter' is not the same as 'winter' to these tools)
>
> > <solution-title></solution-title>
> >
> > 1. `.season`
> > 2. The values can be `Summer` or `Winter` (`jq -c '[.[].season] | sort | unique' < olympics.json`)
> >
> {: .solution}
{: .question}

We'll be using the `select()` filter to select entries matching specific conditions:

> <question-title></question-title>
>
> 1. How would you write the expressions for the following conditions:
>    1. `enrolled` must be 'Yes'
>    2. `age` must be smaller than 75
>    3. `height` cannot be null
>    4. `birthplace` cannot be empty
>
> 2. It is also possible to combine multiple conditions, using `and`, `or`, `not` and parentheses
>    How would you write expressions for the following filtering conditions:
>    1. `height` is larger than 200 or smaller than 160
>    2. `height` is larger than 200 and smaller than 210
>
> > <solution-title></solution-title>
> >
> > 1. The answers are:
> >    1. `select(.enrolled == "Yes")`
> >    2. `select(.age < 75)`
> >    3. `select(.height != null)``
> >    4. `select(.birthplace != "")`
> >
> > 2. The answers are:
> >    1. `select(.height > 200 or .height < 160)`
> >    2. `select(.height > 200 and .height < 210)`
> >
> {: .solution}
{: .question}



Ok, great, now that you've got the hang of writing expressions for this tool, let's create a file with only Winter Olympics. Make sure it is contained in an array, in case we want to do further sorting.

```bash
jq -c '[.[] | select(.season == "Winter")]' < olympics.json > winter.json
```

> <question-title></question-title>
>
> How many lines are in this file? (Hint: use the `length` filter.)
>
> > <solution-title></solution-title>
> >
> > 44,680
> >
> {: .solution}
{: .question}


**Repeat** the step for the Summer Olympics

```bash
jq -c '[.[] | select(.season == "Summer")]' < olympics.json > summer.json
```


> <question-title></question-title>
>
> 1. How many lines do you expect in the this file?
> 2. How many lines are in this file? Were you right?
>
> > <solution-title noprefix>Hints</solution-title>
> >
> > 1. Use the `length` filter, assuming your data is in an array, per the original prompt.
> > 2. Be careful to consider whether these counts include the header line of the file or not
> >
> {: .solution}
>
> > <solution-title noprefix>Answers</solution-title>
> >
> > 1. The original file has 234,522 lines, and the Winter Olympics had 44,680 lines. So we would expect 234,522 - 44,680 = 189,842 rows of data.
> > It is always useful to take a moment to think about the expected outcome, this makes it easier to spot mistakes and will save you time in the long run.
> >
> {: .solution}
{: .question}

## Exercises

Ok, time to train! let's see if you can use the `select` filter to answer the following questions:


> <question-title noprefix>Exercise: Medal winners</question-title>
>
> 1. How many gold medals were handed out?
> 2. How many total medals?
> 3. How many medals were handed out during the 2018 Olympics?
> 4. How many medals were won by individuals with a height between 170 and 180 cm? (inclusive)
> 5. How many gold medals were won by individuals shorter than 160cm or taller than 190?
>
> > <solution-title noprefix>Hints</solution-title>
> >
> > - Column 17 contains information about medals
> > - The possible values are `Gold`, `Silver`, `Bronze`, and `` (empty).
> > - Don't forget that the output (and line count) may include the header line
> > - Do not use quotes on number columns (e.g. year)
> > - You may need parentheses for complex conditions
> > - You may want to `[ ... query ...] | length` to get a nice length
> >
> {: .solution}
>
> > <solution-title noprefix>Answers</solution-title>
> >
> >  1. 8,110   (Expression: `select(.medal == "Gold")`)
> >  2. 24,633  (Expression: `select(.medal == "Gold" or .medal == "Silver" or .medal == "Bronze")`, or `select(.medal != null)`)
> >  3. 131     (Expression: `select(.medal == "Gold" and .year == 2018)` (note: do not use quotes around `2018`, as it is a numerical value))
> >  4. 8,086   (Expression: `select(.medal != null and .height >= 170 and .height <=180)`)
> >  5. 2,333   (Expression: `select(.medal != null and (.height < 160 or .height > 190))` (note: parentheses are important here))
> >
> > Note: these numbers are found by determining the number of lines in the file after each filtering step, and subtracting 1 for the header line.
> >
> {: .solution}
{: .question}


# Counting

A common operation we might want to perform on tables of data, is simple counting. How many times does a certain value appear? For our dataset for instance, we might want to know how many countries participated in each Olympics, how many women, etc; any column that has categorical data that we can count.

Let's start by simply counting how many different Olympic Games we have in our dataset, and how many times it appears (so how many participations there were each year)

We'll need to use the `group_by` filter which takes a key, and then emits arrays with objects with those matching keys:

```bash
echo '[{"foo":1, "bar":10}, {"foo":3, "bar":100}, {"foo":1, "bar":1}]' | jq 'group_by(.foo)
```

So let's try that with ours:

```bash
jq '. | group_by(.year)' < olympics.json
```

We now have a structure like

```
[
    [
        {"year" 1896, ...},
        {"year" 1896, ...},
        {"year" 1896, ...}
    ],
    [
        {"year" 1900, ...},
        {"year" 1900, ...},
        {"year" 1900, ...}
    ],
    ...
]
```

We'll re-structure to one "year group" per line:

```bash
jq '. | group_by(.year) | .[]' < olympics.json
```

So we can again use the `length` function to determine how many items are in each via

```bash
jq '. | group_by(.year) | .[] | [. | length]' < olympics.json
```

And all that's left is to include the game name:

```bash
jq '. | group_by(.year) | .[] | [. | length, .[0].games]' < olympics.json
```

But the output format is not so nice, so let's introduce the `@tsv` filter to turn the resulting data into a table:

```bash
jq '. | group_by(.year) | .[] | [. | length, .[0].games] | @tsv' < olympics.json
```

Where did the " and `\t` come from? We don't want those! That is a result of the resulting data being in a JSON string, so we can tell jq to use `raw` mode with `-r`

```bash
jq -r '. | group_by(.year) | .[] | [. | length, .[0].games] | @tsv' < olympics.json
```

> <question-title></question-title>
>
> 1. How many different Olympic games are in our file?
> 2. Which Olympic games had the most participations? (Tip: set the parameter *"How should the results be sorted?"* to `most common values first`)
>
> > <solution-title></solution-title>
> >
> > 1. 52 games (`[. | group_by(.games) | .[] | [. | length, .[0].games]] | length`)
> >
> >    The resulting file looks something like:
> >
> >     ```
> >     615	1896 Summer Olympics
> >     2503	1900 Summer Olympics
> >     2643	1904 Summer Olympics
> >     3213	1908 Summer Olympics
> >     4610	1912 Summer Olympics
> >     3448	1920 Summer Olympics
> >     5242	1924 Summer Olympics
> >     358	1924 Winter Olympics
> >     4493	1928 Summer Olympics
> >     ...
> >     ```
> >
> > 2. 1996 Summer Olympics. (10501 participations)
> >
> {: .solution}
{: .question}

You may have guessed that like `sort_by`, that we could have selected multiple columns in the `group_by` step. This lets us count on combinations of columns.

Let's try counting the number of men and women in each olympic games.

```bash
jq '. | group_by(.games, .sex) | .[] | [. | length, .[0].games, .[0].sex] | @tsv' -r < olympics.json
```


> <question-title></question-title>
>
> You see the resulting file has a line for every combination of the two columns (games and sex), providing the count for each.
>
> 1. How many women were in the first Olympic games?
>
> 2. Which Olympic games had the most women participants?
>
> > <solution-title></solution-title>
> >
> > 1. 2 women participated in the 1896 Olympics. (note that we cannot be sure if this is two different women, or 1 woman participating twice).
> >    The file looks something like this:
> >    ```
> >    2	F	1896 Summer Olympics
> >    43	F	1900 Summer Olympics
> >    17	F	1904 Summer Olympics
> >    55	F	1908 Summer Olympics
> >    97	F	1912 Summer Olympics
> >    132	F	1920 Summer Olympics
> >    269	F	1924 Summer Olympics
> >    ```
> >
> > 2. 2020 Summer Olympics (4652)
> >
> {: .solution}
{: .question}


Let's say we wanted to know how many different sports there were in each Olympics. If we used the counting tool above, we would get a results file for each combination of sport and olympics, with the number of lines (participations) of each. But we don't really care about the number of lines that have this combination, just the total number of unique sports in each games.

We can use the `unique` filter in our pipeline to discover this. First let's do our group by and iterate over each resulting group:

```bash
jq '. | group_by(.games) | .[]' < olympics.json
```

For fun this time let's organise the output as an object. We know one value we want is the games:

```bash
jq '. | group_by(.games) | .[] | {"games": .[0].games}' < olympics.json
```

But the other we want is the list of sports. We can use our `.[]` to loop over every sport item in that group

```bash
jq '. | group_by(.games) | .[] | {"games": .[0].games, "sports": [.[].sport] }' < olympics.json
```

We're almost there! We just need to unique this list and determine it's length:

```bash
jq '. | group_by(.games) | .[] | {"games": .[0].games, "sports": [.[].sport] | unique | length }' < olympics.json
```

Or if you preferred as a table. But here it seems we can wrap it in `()` to make the parsing clearer:

```bash
jq '. | group_by(.games) | .[] | [.[0].games, ([.[].sport] | unique | length)] | @tsv' -r < olympics.json
```

> <question-title></question-title>
>
> 2. How many sport were in the first Olympics? How many in the latest?
> 3. Which Olympics had the most different sports?
>
> > <solution-title></solution-title>
> >
> > 2. 10 and 38.
> > 3. The 2020 Summer Olympics had the most different sports (38)
> >
> {: .solution}
{: .question}

Save the output as something descriptive.

## Exercises

Ok, let's practice!

> <question-title noprefix>Exercise: Number of participations per country</question-title>
>
> 1. Which country has had the most participations in the Olympics?
> 2. How many countries participated in the first Olympics? How many in the last?
>
> > <solution-title noprefix>Hints</solution-title>
> >
> > 1. Since we are counting instances of a key, we can use `group_by(.team)` and then loop over that to print out the length, and the team name of each of those items.
> > 2. This is basically the same question as "how many women" participated, try modifying that query.
> >
> {: .solution}
>
> > <solution-title noprefix>Answers</solution-title>
> >
> >  1. The United States with 17,286 participations (`cat olympics.json | jq '. | group_by(.team) | .[] |  [. | length, .[0].team] | @tsv' -r`)
> >  2. 15 and 250. (`cat olympics.json | jq '. | group_by(.games) | .[] | [.[0].games, ([.[].team] | unique | length)] | @tsv' -r`)
> >
> {: .solution}
{: .question}

# Computing

Sometimes we want to use the data in our column to compute a new value, and add that to the table. For instance, for our dataset we could caluclate athtletes BMI (using height and weight columns), or their age at time of participation (from year of birth and year of the Olymics). By adding these computed values as a new colum to our datset, we can more easily query the dataset for these values. We can do these types of operations using the `+=` operation on a given object.

As an example, let's calculate the age of each athlete at the time of participation, and add this as a new column to our dataset.

First let's structure our query to loop over every object in the array:

```bash
jq '.[]' < olympics.json
```

Then let's try adding a static value to it:

```bash
jq '.[] | . += {"marco": "polo"}' < olympics.json
```

Do you see the resulting key in the new file? Ok, let's replace it with a computed value:

```bash
jq '.[] | . += {"age": .year - .birth_year}' < olympics.json
```

This errors, but JQ is a bit sensitive to syntax, so we can add some `()` to make things more explicit

```bash
jq '.[] | . += {"age": (.year - .birth_year)}' < olympics.json
```

Ah and then we discover null values!

```bash
jq '.[] | select(.birth_year != null) | . += {"age": (.year - .birth_year)}' < olympics.json
```

And there we go! But this is structured as a stream of objects, maybe we want it again as an array `[]` to ensure it remains in the same shape as our input dataset.

```bash
jq '[.[] | select(.birth_year != null) | . += {"age": (.year - .birth_year)}]' < olympics.json > olympics-with-age.json
```

> <question-title></question-title>
>
> 1. How old was Arnaud Boetsch during his Olympic tennis participation?
>
> > <solution-title></solution-title>
> >
> > 2. Arnaud Boetsch is listed on the first two lines, who turned 27 the year of their Olympics.
> >
> {: .solution}
{: .question}

This was a simple computation, but much more complex mathematical expressions can be computed with this tool. See the help section at the bottom of the tool for a list of all supported operations. In the exercise below, we will compute the BMI for each athlete as an example.


## Exercises

BMI stands for Body Mass Index, is a metric to provide a very crude measure of how healthy your weight is. The formula to compute BMI is:

$$ BMI = weight / (height^2) $$

(with weight in kilograms and height in meters).

> <tip-title>BMI is problematic</tip-title>
> [Adolphe Quetelet](https://en.wikipedia.org/wiki/Adolphe_Quetelet), who invented [BMI](https://en.wikipedia.org/wiki/Body_mass_index), was not a doctor of medicine, instead he was a Belgian astronomer, mathematician, statistician, and sociologist. His data consisted probably entirely of cisgender, white, european men and women {% cite Eknoyan_2007 %}.
> In his defense, he did argue it should not be used at individual levels, however it came to be used as such due to simplicity. (See [Wikipedia's History section](https://en.wikipedia.org/wiki/Body_mass_index#History)). BMI is a poor measure of health, especially for populations with high muscle content, and cannot be simply re-used as-is with anyone other than Anglo Saxons {% cite Caleyachetty_2021 %}, [Wikipedia:Limitations](https://en.wikipedia.org/wiki/Body_mass_index#Limitations).
>
> However as it is an easy to calculate metric, we include the calculation here.
{:.tip}


Let's compute this data for all athletes and add it as a new column!

> <question-title noprefix>Exercise: Calculating BMI</question-title>
>
> 1. How would you express this calculation in jq?
>    - Remember that our height is in cm, and the formula expects height in meters
>    - And that we have null values we must filter out
>
> 2. What is the BMI for Arnaud Boetsch?
>
> > <solution-title noprefix>Hints</solution-title>
> >
> > - division is `/` and multiplication is ` * ` .
> > - Generally we cannot use `^`
> > - Parentheses may be required.
> > - use `select(.value != null)` to remove nulls.
> > - If that isn't eough, you can check that the type is a number with `(.value|type) == "number"` to ensure it really is a number and not e.g. a string.
> > - remember to wrap everything in `[...]` to retain the data shape
> >
> {: .solution}
>
> > <solution-title noprefix>Answers</solution-title>
> > 1. other variations are possible:
> >    ```bash
> >    jq '[.[] | select((.height|type) == "number" and (.weight|type) == "number") | . + {"bmi": (.weight / (.height/100 * .height/100))}]' < olympics.json > olympics-bmi.json
> >    ```
> > 2. 22.69
> >    ```bash
> >    jq '.[] | select(.name=="Arnaud Boetsch") | .bmi' < olympics-bmi.json
> >    ```
> {: .solution}
{: .question}


# Find and Replace

Often you may need to change the contents of a file a bit to fit the expectations of an analysis tool. For instance, our file uses `NA` for missing values, but other conventions included leaving the cell empty, or using `NaN` (not a number) instead. Or, when working with chromosomal data, you may need to add or remove the `chr` prefix from a column before using it as input to a certain tool. In such situations, we can find all occurrences of a certain pattern in our file, and replace it with another value.

If we want to perform such a replacement on a single column in our data, we can use the `gsub` function.

A few of the basics of regular expression, plus some links to further resources are given in the box below:

{% snippet faqs/galaxy/analysis_regular_expressions.md %}

Let's start with a simple example:
Our file uses a mix of `Athina` and `Athens` to indicate the Capital City of Greece in the `city` column.
Let's standardize this by replacing occurrences of `Athina` with `Athens`.

Let's start by filtering out the old spelling:

```bash
jq '.[] | select(.city == "Athina")' < olympics.json
```

We can access just the .city and try out our gsub next:

```bash
jq '.[] | select(.city == "Athina") | .city | gsub("Athina"; "Athens")' < olympics.json
```

Looks good! Let's use the technique we learned in [Computing](#computing) to update the key with a new key.

```bash
jq '.[] | select(.city == "Athina") | . += {"city": (.city | gsub("Athina"; "Athens"))}' < olympics.json
```

Also looks good, ok, that's it!

```bash
jq '[.[] | . += {"city": (.city | gsub("Athina"; "Athens"))}]' < olympics.json
```


Look at the file before and after. Athlete 7 (Patrick Chila) near the top of the `olympics.tsv` file, had a value of Athina in the city column. Verify that it has been changed to Athens.

This was rather simple example, so let's try a few more examples with slightly more complex expressions.


## Exercises

You may have noticed that our file has a lot of missing data. Especially for the earlier years, things like height, weight and birthday of athletes was not registered, or simply not known. In some columns you see these missing values have been replaced with an `NA` (not available) value. In other columns (for example birth place), the cells have simply been left empty.

Different tools may expect different ways of handling missing data. So you may have to change your missing data from empty to `NA`, `NaN`, or something else, between analysis steps

> <hands-on-title>Fill empty cells</hands-on-title>
>
> We will now replace empty cells in the `birth_place` column, to use `null` instead.
> Read the [documentation for PCRE and gsub](https://stedolan.github.io/jq/manual/#RegularexpressionsPCRE)
>
> > <question-title></question-title>
> >
> > 1. Should we use gsub?
> > 2. What is the expression for a line with just whitespace?
> >
> > > <solution-title noprefix>Answer</solution-title>
> > >
> > > 1. Yes, definitely
> > > 2. `^ *$` indicates an empty line, potentially with spaces. (`^` indicated the beginning, zero or more spaces, and `$` the end). The value in the column is treated as the line, since we are looking only in 1 column.
> > >
> > {: .solution}
> >
> > > <solution-title noprefix>Hints</solution-title>
> > > Note that you cannot use `gsub` to replace something with null. Instead you'll need to use `(if ... then "value" else "other" end)`.
> > > You might want to use the gsub instead to ensure lines with single blank characters, and empty lines, all map to the same null result.
> > {: .solution}
> >
> > > <solution-title noprefix>Full Solution</solution-title>
> > > other variations are possible:
> > > ```bash
> > > cat olympics.json | jq '.[] | (if (.birth_place|gsub("^ *$"; "")) == "" then null else .birth_place end)'
> > > ```
> > {: .solution}
> {: .question}
>
{: .hands_on}

Let's do another example, this one using capture groups.

Look at the `birth_day` column. It has values in a format like `12 December`. Suppose we have a tool that expects this data to be in the reverse format, `December 12`. The file is too big to change this manually in every column. But with regular expression tools we can make this replacement easily

We will now change the format in birthday column from `day month` to `month day`, as our boss is American and requested the silly format.

Read the [documentation for PCRE and gsub](https://stedolan.github.io/jq/manual/#RegularexpressionsPCRE)
- Read the text, what settings do you think we need to use?
- Read the Regular expressions 101 FAQ below

  {% snippet faqs/galaxy/analysis_regular_expressions.md %}

> <question-title></question-title>
>
> 1. How do we match on the birthday format? How strict/exact shoule we be here?
> 2. How do we captures both the day and the month?
> 3. How do we refer to the values we captured (for the replacement value)
>
> > <solution-title noprefix>Hints</solution-title>
> >
> > 1. Birthday is one or more digits, followed by a space, followed by one or more letters.
> > 2. Remember that you can capture values using parentheses `(..)`
> > 3. We can name our capturing groups for easier referencing, by providing a name. `(?<day>...)` can be later referred to as `\(.day)`.
> >
> {: .solution}
>
> > <solution-title noprefix>Answers</solution-title>
> >
> > 1. There are multiple solutions here, depending on how strict you want to be
> >    - `\d+ ([a-zA-Z]+)` (not strict, would also match on `142 Septober`
> >    - `[0-9]{1,2} (January|February|March|April|May|June|July|August|September|October|November|December)` (this is much more strict, only matches on 1 or 2 digits, followed by a space, followed by one of the months. But this would still match on 42 April, and may miss it if one month names didn't start with a capital letter)
> >    - `[123]?[0-9] [a-zA-Z]+` this will only allow dates below 40
> >
> >    There are different ways to express this, and there is no one perfect solution. If you know your data is clean, and you do not have to worry about values like 42 Septober,
> >    then you can be less strict. If your data is less clean, you may be more worried about capturing things that aren't valid birthdays and might want to be a bit stricter.
> >    It all depends on your data and use case.
> >
> > 2. `(?<day>[\d]{1,2}) (?<mon>[a-zA-Z]+)` captures both the day and the month
> > 3. We can use `\(.day)` to refer to the day we captured, and `\(.mon)` for the month in our replacement.
> >
> > If you'd like to be fancier than just using `gsub`, you can also try the other capturing functions. E.g. `capture` returns an object with the named matches, which you can add back directly in to the original object:
> >
> > ```bash
> > jq '.[] | . += (.birth_day | capture("(?<birth_day>[0-9]+) (?<birth_month>[A-Z][a-z]+)\\s*"))' < olympics.json
> > ```
> >
> {: .solution}
{: .question}


# Removing Columns

We can remove columns from a table using either `del(.key)` if we just want to remove a few values, or otherwise by selecting, if we just want to keep a few values.

Suppose we want to simplify our file a bit. All we want is file with 4 columns: athlete name, sport, olympic games, and medals.

We can most easily do this by returning an array:

```bash
jq '.[] | [.name, .sport, .games, .medal]' < olympics.json
```

But perhaps we want something a bit more structured:

```bash
jq '.[] | {"name": .name, "sport": .sport, "games": .games, "medal": .medal}' < olympics.json
```

And of course we'll keep our list shape:

```bash
jq '[.[] | {"name": .name, "sport": .sport, "games": .games, "medal": .medal}]' < olympics.json > subset.json
```

Notice that during this step, we also changed the order of the keys. This tool can also be used to re-arrange keys, if you supply all key names but in a different order.

## Exercises

> <question-title noprefix>Exercise: Removing Columns</question-title>
>
> 1. Create a file that is similar to `olympics.json`, but without the first column (athlete_id column)
> 2. Which of the two methods would be easier to use and when? (this can be a personal preference, but think about whether you would subtract or add)
>
> > <solution-title noprefix>Hints</solution-title>
> >
> > 1. Try the `del(.key)` function
> >
> {: .solution}
>
> > <solution-title noprefix>Full Solutions</solution-title>
> > 1. `jq '[.[] | del(.athlete_id)]' < olympics.json`
> > 2. The `del()` function is significantly easier if you are just removing keys, but may be cumbersome if you have a lot of key and only want a few. In that case it may be faster to type those out instead.
> {: .solution}
{: .question}


# Unique

Sometimes, in the course of our data manipulations, we may end up with a file that has duplicate values. In order to filter these out, we can use the `unique` filter.

Let's say we would like to create a list of all unique athletes (id and name).

First we will cut just the `athlete_id` and `name` columns from our dataset

```bash
jq '.[] | {"name": .name, "athlete_id": .athlete_id}' < olympics.json
```

> <question-title></question-title>
>
> 1. Do you see duplication? Why is that?
>
> > <solution-title></solution-title>
> >
> > 1. Yes. For all athletes who participated more than once, the row will be identical.
> >
> {: .solution}
{: .question}

Now let's remove those duplicates. For this filter to work, it needs to be passed an array, not individual values in a stream:

```bash
jq '[.[] | {"name": .name, "athlete_id": .athlete_id}] | unique' < olympics.json
```

> <question-title></question-title>
>
> How many unique athletes do we have?
>
> > <solution-title></solution-title>
> > 94,733
> >
> > ```bash
> > jq '[.[] | {"name": .name, "athlete_id": .athlete_id}] | unique | length' < olympics.json
> > ```
> >
> {: .solution}
{: .question}


# Joining Files

This file contains a lot of information, but we may want to add more information. For example if we had a file with information about each country (population, capital city, etc), we could join that information with our Olympics data, to get a single file with all information in every row.

For example, if we would like to be able to group by continent, to e.g. count athletes, medals etc per continent, we will have to add a `continent` column to our file. To do this we would need a second file that maps each country to the corresponding continent. This is what we will do in the next hands-on section.

We obtained country information data from [DataHub](https://datahub.io/core/country-codes). More information about this file can be found in the description there. It has 56 columns with a wide variety of data about each country (from country codes, to capital city, languages spoken, etc)

Download

```
{{page.zenodo_link}}/files/country-information.json
```

> <question-title></question-title>
>
> 1. How many keys does this file have?
> 2. Which keys(s) in this file are the same as in the `olympics.tsv` file?
>
> > <solution-title></solution-title>
> >
> > 1. The country information file has 56 columns (see [DataHub](https://datahub.io/core/country-codes) for more details).
> > 2. Both files have a `NOC` column with the 3-letter country code (`NOC` stands for National Olympic Committee). However, one is lowercase.
> {: .solution}
{: .question}

We would now like to take our Olympics dataset as the basis, and add columns to every row of this file with some information about the country. In order to join, we will need to have one column that is shared between the two files, on which we can match. The `NOC` column is perfect for this because it is a defined standard. Both files also contain a column with the country name in it, which is also a possible candidate to use for joining, but because it is less standardised, it is safer to use the NOC column. For example, if one file uses "Netherlands", while the other uses "The Netherlands" to indicate the same country, the joining will fail for these rows. So always make sure the columns you join on are compatible!

We can use the join and index commands, per [SQL Style Operators](https://stedolan.github.io/jq/manual/#SQL-StyleOperators) documentation and [this stack overflow answer](https://stackoverflow.com/a/72609410)

```bash
jq '[JOIN(INDEX(input[]; .NOC); .[]; .noc; add)]' olympics.json country-information.json
```

> <question-title></question-title>
>
> 1. What do you expect the output to look like? Were you right?
> 2. How many columns are in the resulting file? What about the NOC column?
> 3. What is a possible downside to this approach?
>
> > <solution-title></solution-title>
> >
> > 1. All the columns from the country information file are added to the end of each row of our olympics dataset
> > 2. Our olympics datset had 17 columns, the country information file has 56 columns. Therefore we have 17+56=73 columns columns in our resulting file. This also means the NOC column
> >    we joined on appears twice in our output.
> > 3. There is a lot of data duplication in this file now. The exact same country information is added to every line of every athlete from a certain country.
> >    This means much larger file size.
> >    If you do not need all these columns, it could save you a lot of space to remove unneeded columns from the `country-information.json` file, before joining.
> >
> {: .solution}
{: .question}

# Concatenating

Concatenation of two files simple means putting the contents of the two files together, one after the other. Our dataset was created in 2021, but since then we've had another Olympic event, the 2022 Winter Olympics in Beijing. If we have the same data for this latest Olympics, we could simply add the rows from the 2022 games to our current file with data, in order to create a single file with all data from 1896 to 2022.

First, let's get this data for the 2022 Olympics

```
wget {{page.zenodo_link}}/files/olympics_2022.json
```

View the new dataset, does it have the same structure as our original `olympics.json` file?

> <question-title></question-title>
>
> 1. Does the new file have the same structure?
> 2. Can we simply add the lines of the new files to the end of our existing olympics dataset?
>
> > <solution-title></solution-title>
> >
> > 1. Yes, this file has all the same columns, in the same order, so concatenation should be relatively straightforward.
> > 2. If we simply put the contents of this file after our existing dataset, we will have a second header line in the middle of our data rows. It is best to remove the header line from the second dataset after we have verified it is compatible. This way we will only add the real data rows to our dataset.
> >
> {: .solution}
{: .question}

Since this new dataset has the exact same structure (number and order of columns), we can simple add the lines from this file to the end of our existing `olympics.json` file.
For this, we'll need to use a new argument to `jq`, namely `-s` or `--slurp` which allows us to access all files as if they were one big indexed array.

```bash
cat olympics.json olympics_2022.json | jq -s '.[0] + .[1]' > olympics-combined.json
```

Alternatively we could loop over every array in the parent, and combine them. This would support more than two files:

```bash
cat olympics.json olympics_2022.json | jq -s '[.[][]]' > olympics-combined.json
```

Now this only works so simply because our two datasets had the same structure; two arrays. If your data comes from different sources, you may have to do some additional data manipulation before you can concatenate, e.g. to make sure the columns match, or how each file deals with missing data (empty cells, `NA`, `NaN` or something else).


# Splitting Files

This dataset contains a lot of data, we might want to split the data into one file per Olympic games, or have one file for all the winter games, and another file for all the summer games. In these situtations, where we want to create new, smaller files, based on the values in a column, we can use the tool.

> <tip-title>Not generically possible with only JQ</tip-title>
> It is not generically possible to e.g. split a file by value, with only jq. We need to use other tools, e.g. for loops, in addition to accomplish this.
{: .tip}

We would like to create a separate file for each Olympic event. First let's get a list of every event.

```bash
jq '[.[] | .games] | unique[]' -r < olympics.json
```

Easily done. Now we'll construct a loop over these results

```bash
jq '[.[] | .games] | unique[]' -r < olympics.json | while read -r game; do
    echo "Game: $game";
done
```

Here we can then use this to create a new file per the game identifier:

```bash
jq '[.[] | .games] | unique[]' -r < olympics.json | while read -r game; do
    echo "Selecting: $game";
    jq '[.[] | select(.games == "'$game'")]' < olympics.json > olympics-"$game".json
done
```


> <question-title></question-title>
>
> 1. How many different files do we get?
> 2. How many lines are in the file for the 1908 Summer Olympics?
>
> > <solution-title></solution-title>
> >
> > 1. We get 1 collection with 52 files in it, one per Olympic games:
> > 2. 61049 lines, 3213 entries.
> {: .solution}
{: .question}


## Exercises

Let's practice this a bit more, see if you can use the split file to answer the following questions:

> <question-title noprefix>Exercise: sort by height</question-title>
>
> 1. Create two files, one for Summer Olympics, one for Winter Olympics. Which has more lines?
> 2. Split the file by sport, how many sports have there been at the Olympics?
> 3. Split the file by medal, would you expect the output files to be equal sizes?
>
> > <solution-title noprefix>Answers</solution-title>
> >
> >  1. The summer olympics file has more lines than the winter olympics file (the first winter Olympics wasn't held until 1924)
> >  2. Each file in the resulting collection represents a sport, there are 91 datasets, representing 91 different sports in the history of the Olympics.
> >  3. The files for gold silver and bronze are roughly equal, but not exactly (think of situations like shared second place, then 2 silver medals are handed out and no bronze medal)
> {: .solution}
>
{: .question}



# Conclusion

These tools and operations covered in the tutorial are just a few examples of some of the most common operations. There are many more tools available in Galaxy that perform other data manipulations. We encourage you to look around the `General Text Tools` section in the Galaxy toolbox (left panel) to find and explore more data manipulation tools. The more comfortable you are performing these kinds of steps, the more you can get out of Galaxy!


# Exercises: Putting it all together!

This section provides a number of exercises that require you to combine two or more of the techniques you learned in this tutorial. This is a great way to practice your data manipulation skills. Full solutions are provided for every exercise (i.e. all tools and settings), but for many of these exercises there will be multiple solutions, so if you obtained the same results in a different way, that is correct too!

> <question-title noprefix>Exercise 1: Finding shortest/lightest athlete</question-title>
>
> If you have done exercises in the [sorting](#sorting) section, you noticed that finding the shortest athlete ever to compete was not easy,
> because all the rows with missing height data (`NA`) in the column were sorted to the top. We need to filter out these values first, then
> perform our sort, so that our answer is on top.
>
> 1. Find the shortest athlete ever to compete in the Olympics
> 2. Find the shortest athlete of the Winter Olympics
> 2. Find the lightest athlete of the *most recent* Summer Olympics
>
> > <solution-title noprefix>Hints</solution-title>
> >
> > 1. You will need to filter out the columns with (`NA`) in the height column first
> > 2. You will need to filter by season as well
> > 3. You will need to filter out missing data in the weight column, filter out Summer Olympics, then sort (by 2 columns)
> > 4. Does the order in which you perform these steps matter?
> >
> {: .solution}
>
> > <solution-title noprefix>Answers</solution-title>
> >
> >  1. Lyton Mphande and  Rosario Briones were both 127 cm tall, competing in boxing and gymnastics respectively
> >  2. Carolyn Krau was a 137 cm tall figure skater.
> >  3. Flávia Saraiva was the lightest athlete (31kg), she was a Artistic Gymnast from Brazil.
> >
> {: .solution}
>
> > <solution-title noprefix>Full solution</solution-title>
> >
> > 1. First we filter out the NA values from the height column:
> >
> >    ```bash
> >    cat olympics.json | jq '[.[] | select(.height != null)]'
> >    ```
> >
> >    Then we can sort by height, in ascending order to get the shortest athletes on top:
> >
> >    ```bash
> >    cat olympics.json | jq '[.[] | select(.height != null)] | sort_by(.height) | .[0]'
> >    ```
> >
> > 2. We can take the output from the first exercise, and filter for only Winter Olympics:
> >
> >    ```bash
> >    cat olympics.json | jq '[.[] | select(.height != null and .season == "Winter")] | sort_by(.height) | .[0]'
> >    ```
> >
> > 3. First we filter out the NA values from the weight column:
> >
> >    ```bash
> >    cat olympics.json | jq '[.[] | select(.weight != null and .games == "2020 Summer Olympics")] | sort_by(.weight) | .[0]'
> >    ```
> >
> > 4. In this case, changing the order of these operations will still have led to the same answer. However, sorting takes more time than filtering,
> >    so it is usually faster to filter first, then do the sorting on a smaller file. Especially noticeable if you have large dataset.
> >
> {: .solution}
>
{: .question}

Ok, let's try another exercise! You might have noticed that some of the values in the weight column were classes, e.g. `63-78` (kg).
In this exercise, we will first clean up this weight column to avoid ranges (by taking the lower number in the range for every
row a weight range is used), and then answering the question a couple of questions around the weight of athletes.

> <question-title noprefix>Exercise 2: Data cleaning and computations of the weight column</question-title>
>
> 1. Get a list of all the values that occur in the weight value, take note of all the values that are not a single number or `NA`; anything else should be cleaned up
> 2. Clean up the weight value so that we only have single numbers; weight classes (e.g. `63-78`) should be replace by the lower bound (`63`) of that class
> 3. How heavy was the lightest woman competing in the Biathlon? And the heaviest?
>
> > <solution-title noprefix>Hints</solution-title>
> >
> > 1. You can use ` | sort | uniq -c` (the unix tools) to simplify some of the counting.
> > 2. You can use the `gsub` filter to clean up the weight.
> >    This tool is a bit complex, if you haven't done so yet, we strongly recommend you work through the [Find and Replace section](#find-and-replace) first.
> >    To check if your cleaning worked, re-run the counting from step 1.
> >
> {: .solution}
>
> > <solution-title noprefix>Answers</solution-title>
> >
> > 1. We see values such as `100-105`, but also `77,5`, `100, 104` (notice the space!) and even `76, 77, 79` to indicate weight ranges or other unexpected weight values.
> >    All these variations must be converted to single numbers. For this example we will simply convert these ranges to their lower number, and remove the rest when we
> >    clean the column (taking the upper or median value would be equally valid,  but for our solution we will remove everything after the first number.
> >
> > 2. There are many possible regular expression that will clean this column. You can use a single check to catch all the exceptions, or use multiple, there is no one right answer.
> >    Check if your cleaning is successful using the select/sort -u/uniq -c pipeline again. Make sure you only have single number values in the weight column before
> >    proceeding to question 3 in this exercise. Also make sure the weights you have make sense (e.g. if you end up with columns having a weight of `0`, something went wrong somewhere!)
> >
> > 3. lightest woman weighed 45 kilograms, the heaviest 82.
> >
> {: .solution}
>
> > <solution-title noprefix>Full Solutions</solution-title>
> >
> > 1. `cat olympics.json | jq '.[] | select(.weight != null) | .weight' | sort | uniq -c`
> >
> >    This should output a file like:
> >
> >    ```
> >    2 "100, 104"
> >    2 "100, 107"
> >    2 "100-105"
> >    1 "100-106"
> >    ...
> >    ```
> >
> >    We will want to change weight ranges to single numbers in the next step
> >
> > 2. We can use the `tonumber` filter here to convert the output of `gsub` into a proper number.
> >
> >    ```bash
> >    cat olympics.json | jq '[.[] | select(.weight != null) | . + {"weight_fixed": (if (.weight|type) == "string" then (.weight |gsub("[-,].*"; "") | tonumber) else .weight end)}]'
> >    ```
> >
> >
> >    **NOTE:** there will be a lot of valid answers here, you could do it all in one check like the solution above, or in multiple checks (e.g. one for each different kind of way to
> >    denote a weight range you found in step 1). As long as a re-run of the counting step shows only single numbers in the
> >    column after your find & replace step, then your answer was correct! Also make sure the weights make sense (e.g. if your resulting column has weights of `0` or other small numbers, this
> >    doesn't make sense and you will have to tweak your expression. Don't feel bad if this takes a bit of trial and error, that is expected!
> >
> > 3. Doing this in two steps is a bit easier to read/cleaner:
> >    
> >    ```bash
> >    cat olympics.json | jq '[.[] | select(.weight != null) | . + {"weight_fixed": (if (.weight|type) == "string" then (.weight |gsub("[-,].*"; "") | tonumber) else .weight end)}]' > olympics-fixed-weights.json
> >    cat olympics-fixed-weights.json | jq '. | group_by(.sport, .sex) | .[] | {"sport": .[0].sport, "sex": .[0].sex, "weight_min": ([.[] | .weight_fixed] | min), "weight_max": ([.[] | .weight_fixed] | max)}' -c
> >    ```
> >
> >    Your resulting file should look something like:
> >
> >    ```
> >    {"sport":"Alpine Skiing","sex":"F","weight_min":45,"weight_max":90}
> >    {"sport":"Alpine Skiing","sex":"M","weight_min":50,"weight_max":138}
> >    {"sport":"Archery","sex":"F","weight_min":42,"weight_max":95}
> >    {"sport":"Archery","sex":"M","weight_min":46,"weight_max":130}
> >    {"sport":"Art Competitions","sex":"M","weight_min":59,"weight_max":93}
> >    ```
> >
> {: .solution}
{: .question}

Congratulations! You have now mastered the basics of data manipulation! There are a lot more data manipulation operations available in Galaxy that you may need.
Please explore the tools for yourself, and check back with this tutorial often, we plan to add more sections and exercises over time!
