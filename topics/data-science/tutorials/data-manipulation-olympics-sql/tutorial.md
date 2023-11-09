---
layout: tutorial_hands_on

title: 'Data Manipulation Olympics - SQL'
zenodo_link: 'https://zenodo.org/record/6803028'
tags:
- cyoa
- sql
questions:
- How can I do basic data manipulation in SQL?
- Which functions are available to convert, reformat, filter, sort etc my data stored in a database?
objectives:
- Familiarize yourself with data manipulation in SQL
- Perform basic SQL query tasks in Galaxy
- Reason about the expected outcome of tools
time_estimation: 1h
key_points:
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

notebook:
    language: sql

abbreviations:
    SQL: "Structured Query Language"

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

| Operation                | Description                                                         | Galaxy Tool                                   |
| ------------------------ | ------------------------------------------------------------------- | --------------------------------------------- |
| Compute on rows          | to derive new column values from existing ones                      | `SELECT x * 2 FROM y`                         |
| Concatenate datasets     | one after the other                                                 | `SELECT * FROM x; union all; SELECT * FROM y` |
| Counting                 | Count occurrences of values in a column                             | `SELECT count(x) FROM y where x = 'value'`    |
| Cut Columns              | By header name                                                      | `SELECT x, y, z FROM a`                       |
| Filter                   | Remove rows based on values in one or more columns                  | `... WHERE x = 'value'`                       |
| Find and Replace         | in a specific column                                                | `REPLACE()`, `regexp_replace` in postgresql   |
| Group on a column        | And perform simple operations (count, mean, min, max etc)           | `... GROUP BY x ...`                          |
| Join two Datasets        | side by side on a specified field                                   | `SELECT * FROM x, y JOIN x.id = y.id`         |
| Select First lines       | Good for finding top 10s or saving header lines                     | `... LIMIT 10`                                |
| Sort on a column         | Change the order of the rows based on values in one or more columns | `... ORDER BY x ASC`                          |
| Unique                   | Remove duplicate rows                                               | `SELECT DISTINCT x FROM y`                    |


In this tutorial, these functions are explained in more detail, and we provide some exercises for you to practice.

# Background

In this tutorial, we will use as our dataset a table with results from the Olympics, from the games in Athens in 1896 until Tokyo in 2020. The objective is to familiarize you with a large number of the most important data manipulation tools in Galaxy. Much like the Olympics, there are many different disciplines (types of operations), and for each operation there are often multiple techniques (tools) available to athletes (data analysts, you) that are great for achieving the goal.


![image of olympic rings, logo and two athletes around the words "Data Analysis Olympics"]({% link topics/introduction/tutorials/data-manipulation-olympics/images/cover.jpg %})


We will show you many of these commonly needed data manipulation operations, and some examples of how to perform them in Galaxy. We also provide many exercises so that you can train your skills and become a data manipulation Olympian!

# Preamble
```sql
# This preamble sets up the sql "magic" for jupyter. Use %%sql in your cells to write sql!
!python3 -m pip install ipython-sql sqlalchemy
!wget -c {{page.zenodo_link}}/files/olympics.db
import sqlalchemy
engine = sqlalchemy.create_engine("sqlite:///olympics.db")
%load_ext sql
%sql sqlite:///olympics.db
%config SqlMagic.displaycon=False
```

# Download Data

Before we can do any manipulation, we will need some data. Let's download our table with Olympics results now.

```sql
SELECT * FROM olympics LIMIT 10;
```

And now we can start querying the database:

```sql
SELECT
    name
FROM
    sqlite_schema
```

> <question-title></question-title>
>
> 1. What tables are available?
> 2. How are they structured?
>
> > <solution-title></solution-title>
> >
> > 1. `countries`, `olympics`, `olympics_2022`
> > 2. Each are tables with 10 or more columns.
> {: .solution}
{: .question}


## About this dataset

The data was [obtained](https://github.com/UOSCS/Olympic_Athletes) from [Olympedia](https://www.olympedia.org/). The `olympics` table contains
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

This tutorial is structured a bit differently than most. **You do not have to do the steps in the order they are presented below.** Every section in this tutorial uses the dataset you just uploaded (the `olympics.db` file) as input, so you can jump to any section in this tutorial right now if you have a particular data manipulation operation in mind you want to learn more about.

# Sorting

We have a lot of data in this file, but it is ordered by the athlete ID number, which is a somewhat arbitrary and meaningless number. But we can sort the rows in this file to something more convenient, for example alphabetically by name of the athlete, or chronologically by year of the Olympics.

In `SQL` we can use the `ORDER BY` clause. We'll start by limiting our results, as every table in this dataset is quite large.

```sql
SELECT NOC, `CLDR display name` FROM countries LIMIT 30;
```

You can use `ORDER BY column-name ASC` or `ORDER BY column-name DESC` to sort the data ascending or descending.

```sql
SELECT NOC, `CLDR display name` FROM countries ORDER BY NOC ASC LIMIT 30;
```

```sql
SELECT NOC, `CLDR display name` FROM countries ORDER BY NOC DESC LIMIT 30;
```

So let's sort our file in chronological order, based on the year of the Olympic games:

> <question-title></question-title>
>
> 1. Which column contains the year?
>
> > <solution-title></solution-title>
> >
> > 1. `year`
> >
> {: .solution}
{: .question}

```sql
SELECT * FROM olympics ORDER BY year LIMIT 30;
```

If we wanted to do it in reverse, we could just use `order by year desc`

```sql
SELECT * FROM olympics ORDER BY year DESC LIMIT 30;
```

> <question-title></question-title>
>
> 1. Write a query to access only the first entry.
> 2. Which athlete is listed at the top of the file now?
>
> > <solution-title></solution-title>
> > 1. We can use `LIMIT` for this.
> >    ```sql
> >    select * from olympics order by year limit 1;
> >    ```
> > 1. J. Defert. Who competed in a Tennis event 1896 Summer Olympics in Athens.
> >
> {: .solution}
{: .question}

This is great, but maybe it would make more sense to sort alphabetically by athlete name *within each year*.

## Sort on multiple columns at once

So we want to sort twice, first by year, an then within each year, we sort again alphabetically by name.

We will sort the file in chronological order based on the year of the Olympic games

```sql
SELECT * FROM olympics ORDER BY year, name LIMIT 30;
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
> >
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
> > ```sql
> > SELECT * from olympics where height is not null ... ;
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
> > 1. `select * from olympics  where height is not null order by height desc limit 1;`
> > 2. `select * from olympics  where height is not null order by height asc limit 1`
> > 3. `select * from olympics  where height is not null order by year desc, height desc limit 1`
> {: .solution}
>
{: .question}


# Filtering

This file contains a lot of data, but we may only be interested in a subset of this data. For example, we may only want to look at one particular Olympics, or one particular sport. In such cases we can filter the dataset. This will create a new dataset, removing any rows that are not of interest to us (i.e. that don't meet the criteria we provide).

We will filter the file to show only winter Olympics
Look at the `olympics` table and answer the following questions

> <question-title></question-title>
>
> 1. Which key contains this information?
> 2. Which values can this column have? (make sure to notice capitalisation, 'Winter' is not the same as 'winter' to these tools)
>
> > <solution-title></solution-title>
> >
> > 1. `season`
> > 2. The values can be `Summer` or `Winter` (`select distinct season from olympics`)
> >
> {: .solution}
{: .question}

We'll be using the `WHERE` filter to select entries matching specific conditions:

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
> >    1. `select * from olympics where enrolled = 'Yes'`
> >    2. `select * from olympics where age < 75`
> >    3. `select * from olympics where height is not null``
> >    4. `select * from olympics where birthplace != ""`
> >
> > 2. The answers are:
> >    1. `select * from olympics where height > 200 or height < 160`
> >    2. `select * from olympics where height > 200 and height < 210`
> >
> {: .solution}
{: .question}



Ok, great, now that you've got the hang of writing expressions for this tool, let's create a file with only Winter Olympics. Make sure it is contained in an array, in case we want to do further sorting.

```sql
CREATE TABLE winter AS SELECT * FROM olympics WHERE season = 'Winter'
```

> <question-title></question-title>
>
> How many entries are in this file? (Hint: use `count(*)`)
>
> > <solution-title></solution-title>
> >
> > 44,680
> >
> {: .solution}
{: .question}


**Repeat** the step for the Summer Olympics

```sql
CREATE TABLE summer AS SELECT * FROM olympics WHERE season = 'Summer'
```

> <question-title></question-title>
>
> 1. How many lines do you expect in the this file?
> 2. How many lines are in this file? Were you right?
>
> > <solution-title noprefix>Hints</solution-title>
> >
> > 1. Use the `count(*)` select
> > 2. Be careful to consider whether these counts include the header line of the file or not
> >
> {: .solution}
>
> > <solution-title noprefix>Answers</solution-title>
> >
> > 1. The original file has 234,522 entries, and the Winter Olympics had 44,680 entries. So we would expect 234,522 - 44,680 = 189,842 rows of data.
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
> >
> {: .solution}
>
> > <solution-title noprefix>Answers</solution-title>
> >
> >  1. 8,110   (Expression: `SELECT count(*) FROM olympics WHERE medal == "Gold"`)
> >  2. 24,633  (Expression: `SELECT count(*) FROM olympics WHERE medal == "Gold" or medal == "Silver" or medal == "Bronze")`, or `medal != null`)
> >  3. 131     (Expression: `SELECT count(*) FROM olympics WHERE medal == "Gold" and year == 2018` (note: do not use quotes around `2018`, as it is a numerical value))
> >  4. 8,086   (Expression: `SELECT count(*) FROM olympics WHERE medal is not null and height >= 170 and height <=180`)
> >  5. 2,333   (Expression: `SELECT count(*) FROM olympics WHERE medal is not null and (height < 160 or height > 190)` (note: parentheses are important here))
> >
> > Note: these numbers are found by determining the number of lines in the file after each filtering step, and subtracting 1 for the header line.
> >
> {: .solution}
{: .question}


# Counting

A common operation we might want to perform on tables of data, is simple counting. How many times does a certain value appear? For our dataset for instance, we might want to know how many countries participated in each Olympics, how many women, etc; any column that has categorical data that we can count.

Let's start by simply counting how many different Olympic Games we have in our dataset, and how many times it appears (so how many participations there were each year)

We'll need to use the `group by` syntax which takes a key, and then groups by those values.

> <question-title></question-title>
>
> 1. How many different Olympic games are in our file?
> 2. Which Olympic games had the most participations? (Tip: use order by)
>
> > <solution-title></solution-title>
> >
> > 1. 52 games (`select count(*), games from olympics group by games`)
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

You may have guessed that like `order by`, that we could have selected multiple columns in the `group by` step. This lets us count on combinations of columns.

Let's try counting the number of men and women in each olympic games.

```sql
select count(*), games, sex from olympics group by games, sex
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
> > 1. 2 women participated in the 1896 Olympics. (note that we cannot be sure if this is two different women, or 1 woman participating twice, in this query. Do you know any way we could query that? Try it out!)
> >    The results looks something like this:
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

Let's say we wanted to know how many different sports there were in each Olympics. If we used the counting query above, we would get resultsfor each combination of sport and olympics, with the number of lines (participations) of each. But we don't really care about the number of lines that have this combination, just the total number of unique sports in each games.

We can use the `distinct` filter in our pipeline to discover this. First let's do our group by and iterate over each resulting group:

```sql
select games, sport from olympics group by games;
```

And let's count all of their appearances

```sql
select games, count(sport) as sports from olympics group by games;
```

But those results still aren't distinct, those numbers are far too high. So let's use distinct:

```sql
select games, count(distinct sport) as sports from olympics group by games;
```

We're almost there! Let's sort this

```sql
select games, count(distinct sport) as sports from olympics group by games order by sports asc;
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
> > 1. Since we are counting instances of a key, we can use `group by team` and then loop over that to print out the length, and the team name of each of those items.
> > 2. This is basically the same question as "how many women" participated, try modifying that query.
> >
> {: .solution}
>
> > <solution-title noprefix>Answers</solution-title>
> >
> >  1. The United States with 17,286 participations (`select team, count(team) as count from olympics group by team order by count desc;`)
> >  2. 15 and 250. (`select games, count(distinct team) as teams from olympics group by games;`)
> >
> {: .solution}
{: .question}


# Grouping

Often we may want to group rows based on a value in a column, and perform some operation on the resulting rows. For example we would like to group the olympics data by one value (e.g. year, country, sport), and determine some value for each group (e.g. number of medals won, average age of athletes, etc).

In the [counting](#counting) section of this tutorial we show how to get answers that require a count (e.g. number of medals won), but sometimes we want to do something more complex, like calculating the average height of athletes in a group, say per country or per sport. This section will show some example of these types of questions.

We can use continue to use group by for this, but now we'll need the max and min aggregate operations. Essentially every time we use `group by` we need to use an aggregation like finding the maximum, minimum, or counting the number of results.

> <hands-on-title>Tallest athlete per sport</hands-on-title>
>
> We would like to answer the following question: *How tall was the tallest athlete of each sport?*
>
> > <question-title></question-title>
> >
> > 1. How tall was the tallest athlete in basketball? And what about karate?
> > 2. Why do some sports have null values?
> >
> > > <solution-title></solution-title>
> > > ```sql
> > > select max(height), min(height),sport from olympics group by sport
> > > ```
> > >
> > > 1. Basketball's tallest athlete was 192cm. For Karate it is 163.
> > > 2. Our dataset had quite a number of `null` (unknown) values in the height column, especially for the earlier Olympics. These are preserved in the outputs.
> > {: .solution}
> {: .question}
{: .hands_on}


## Grouping on multiple columns

You may have noticed that we could also provide multiple columns to group on. If we do this, we can compute values for combinations of groups, such as sex and sport, to find e.g. the tallest woman in basketball or the shortest man per Olympics. There are also many more options for the computation we perform, so perhaps we are more interested not in the tallest athlete, but the average height. Let's perform some of these slightly more advanced queries now.


> <hands-on-title>Average height of men and women per sport</hands-on-title>
>
> The question we would like to answer here, is what is the average height for men and women per sport?
>
> ```sql
> select avg(height), sex, sport from olympics group by sex, sport;
> ```
>
> See if you can answer the following questions based on the output file.
>
> > <question-title></question-title>
> >
> > 1. What is the average height of women participating in archery?
> > 2. What is the average height of men participating in [ballooning](https://en.wikipedia.org/wiki/Ballooning_at_the_1900_Summer_Olympics)?
> > 3. Why do some values have `null` instead of a height?
> > 4. Why do some sports not have a value for one of the sexes?
> > 5. Can you find a sport where women were taller than the men? (Hint: it starts with the letter A)
> >
> > > <solution-title></solution-title>
> > >
> > > 1. 167.25677031093 cm
> > > 2. 170 cm
> > > 3. If none of the rows in the group had height data available, it will output `nan` (not a number) instead. This is most common for sports that were only featured a few times in the early years of the Olympics.
> > > 4. Sports such as artistic swimming only exist for women, so no M appears in the data for that group, so there simply is no row for the mean height of men doing artistic swimming in our output.
> > > 5. [Art Competitions](https://en.wikipedia.org/wiki/Art_competitions_at_the_Summer_Olympics)
> > >
> > > If all went well, your output file should look something like:
> > >
> > > ```
> > > GroupBy(sport)	     GroupBy(sex)  mean(height)
> > > Aeronautics         M             nan
> > > Alpine Skiing       F             167.38324708926
> > > Alpine Skiing       M             178.18747142204
> > > Alpinism            M             nan
> > > Archery             F             167.25677031093
> > > Archery             M             178.5865470852
> > > Art Competitions    F             175.33333333333
> > > Art Competitions    M             173.97260273973
> > > Artistic Gymnastics F             156.15316901408
> > > ```
> > {: .solution}
> {: .question}
>
{: .hands_on}


## Exercises

> <question-title noprefix>Exercise: Grouping and computing</question-title>
>
> 1. How tall is the shortest woman Badminton player to win a gold medal?
> 2. What is the average height of athletes from team Denmark in the 1964 Olympics? (Note: 1964 has summer and winter olympics)
>
> > <solution-title noprefix>Hints</solution-title>
> >
> > 1. We need to group on 3 columns: medal, sport and sex, and then select the `min`.
> > 2. We need to group on 2 columns: country (team) and year, then compute the average over height.
> {: .solution}
>
> > <solution-title noprefix>Answers</solution-title>
> >
> >  1. 161 cm.
> >  2. mean height: 175.91304347826, standard deviation: 7.0335410308672`
> {: .solution}
>
> > <solution-title noprefix>Full Solutions</solution-title>
> > ```sql
> > select min(height), medal, sport, sex from olympics group by medal, sport, sex;
> > select avg(height), team, games from olympics where team = "Denmark" group by team, games;
> > ```
> >
> {: .solution}
>
{: .question}


# Computing

Sometimes we want to use the data in our column to compute a new value, and add that to the table. For instance, for our dataset we could caluclate athtletes BMI (using height and weight columns), or their age at time of participation (from year of birth and year of the Olymics). By adding these computed values as a new colum to our datset, we can more easily query the dataset for these values. We can do these types of operations on the fly, and then if we like, store them as (temporary) tables.

As an example, let's calculate the age of each athlete at the time of participation, and add this as a new column to our dataset.

```sql
select year - birth_year as age, games from olympics LIMIT 30;
```

If we want to save that result to make it easier to query, then we have a couple options.

1. Don't store it, calculate on demand. Very storage efficient.

   ```sql
   select noc, name, ..., year - birth_year as age from olympics
   ```

1. Create a temporary table with this data, great if we only need it temporarily

   ```sql
   create temporary table olympics_ages as select *, year - birth_year as age from olympics
   ```

2. Create a new permanent table with this data

   ```sql
   create table olympics_ages as select *, year - birth_year as age from olympics
   ```

3. Update the existing table by adding it as a new column

   ```sql
   alter table olympics add column age int;
   update olympics set age = year - birth_year;
   ```

> <question-title></question-title>
>
> 1. How old was Arnaud Boetsch during his Olympic tennis participation?
>
> > <solution-title></solution-title>
> >
> > 1. Arnaud Boetsch is listed on the first two lines, who turned 27 the year of their Olympics.
> >
> {: .solution}
{: .question}

This was a simple computation, but much more complex mathematical expressions can be computed with this tool. In the exercise below, we will compute the BMI for each athlete as an example.


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
> 1. How would you express this calculation in SQL?
>    - Remember that our height is in cm, and the formula expects height in meters
>    - And that we have null values, but we can ignore those as they will make the final calculation null as well
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
> >
> >    ```sql
> >    create temporary table olympics_bmi as select *, weight / (height / 100 * height / 100) as bmi from olympics;
> >    ```
> >
> > 2. 22.69
> >
> >    ```sql
> >    select * from olympics_bmi where name like 'Arnaud Boetsch';
> >    ```
> {: .solution}
{: .question}


# Find and Replace

Often you may need to change the contents of a file a bit to fit the expectations of an analysis tool. For instance, our database uses `null` for missing values, but other conventions included leaving the cell empty instead. Or, when working with chromosomal data, you may need to add or remove the `chr` prefix from a column before using it as input to a certain tool. In such situations, we can find all occurrences of a certain pattern in our file, and replace it with another value.

If we want to perform such a replacement on a single column in our data, we can use an update statement.

A few of the basics of regular expression, plus some links to further resources are given in the box below:

{% snippet faqs/galaxy/analysis_regular_expressions.md %}

Let's start with a simple example:
Our file uses a mix of `Athina` and `Athens` to indicate the Capital City of Greece in the `city` column.
Let's standardize this by replacing occurrences of `Athina` with `Athens`.

Let's start by filtering out the old spelling:

```sql
select * from olympics where city = 'Athina' limit 30;
```

Let's try replacing it:

```sql
update olympics
set city = 'Athens'
where city = 'Athina';
```

Look at the file before and after. Athlete 7 (Patrick Chila) near the top of the `olympics.tsv` file, had a value of Athina in the city column. Verify that it has been changed to Athens.

This was rather simple example, so let's try a few more examples with slightly more complex expressions.

## Exercises

You may have noticed that our file has a lot of missing data. Especially for the earlier years, things like height, weight and birthday of athletes was not registered, or simply not known. In some columns you see these missing values have been replaced with an `NA` (not available) value. In other columns (for example birth place), the cells have simply been left empty.

Different tools may expect different ways of handling missing data. So you may have to change your missing data from empty to `NA`, `NaN`, or something else, between analysis steps

> <hands-on-title>Fill empty cells</hands-on-title>
>
> We will now replace empty cells in the `birth_place` column, to use `null` instead.
>
> > <solution-title noprefix>Hints</solution-title>
> > Remember that comparison to nulls is done with `is` instead of `=`
> {: .solution}
>
> > <solution-title noprefix>Answers</solution-title>
> > other variations are possible:
> > ```sql
> > update olympics set birth_place = null where birth_place = '';
> > ```
> {: .solution}
>
{: .hands_on}

Let's do another example, this one splitting and re-constructing strings.

Look at the `birth_day` column. It has values in a format like `12 December`. Suppose we have a tool that expects this data to be in the reverse format, `December 12`. We would not want to do this manually, but with sql we can make this replacement easily

We will now change the format in birthday column from `day month` to `month day`, as our boss is American and requested the silly format.

First we need to understand that sqlite does not ship a regex engine, thus we cannot use familiar regular expressions. Instead we can make use of `instr(string, search)` to find the location of a substring like ` `, identifying where the day stops and the month starts. Then we can use `substr(string, start)` and `substr(string, start, end)` to chop up our date string.

In SQL, concatenation is done with `||`.

> <question-title></question-title>
>
> 1. How do we captures both the day and the month?
>
> > <solution-title noprefix>Hints</solution-title>
> >
> > 1. We should use something like `substr(birth_day, instr(birth_day, ' '))`
> >
> {: .solution}
>
> > <solution-title noprefix>Answers</solution-title>
> >
> > ```sql
> > select
> >   birth_day,
> >   instr(birth_day, ' ') as idx,
> >   substr(birth_day, instr(birth_day, ' ')) as month,
> >   substr(birth_day, 0, instr(birth_day, ' ')) as day
> > from olympics
> > limit 10
> > ```
> > 
> > We can make our final query:
> > 
> > ```sql
> > select
> >   substr(birth_day, instr(birth_day, ' ')) || substr(birth_day, 0, instr(birth_day, ' ')) as birth_day_new
> > from olympics
> > limit 10
> > ```
> > And then store this as a new column (using `update`) or as a new table (using `create [temporary] table`)
> >
> {: .solution}
{: .question}


# Removing Columns

In sqlite you cannot remove columns, it is not supported. Proper databases like Postgres and MySQL support this operation.
Instead for sqlite you could select the columns you want to keep, create a new table from that, and then delete the original table.

Other databases:

```
alter table delete column NAME from olympics
```

# Unique

Sometimes, in the course of our data manipulations, we may end up with a file that has duplicate values. In order to filter these out, we can use the `distinct` filter.

Let's say we would like to create a list of all unique athletes (id and name).

First we will just select the `athlete_id` and `name` columns from our dataset

```sql
select name, athlete_id from olympics LIMIT 30
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

Now let's remove those duplicates.

```sql
select distinct name, athlete_id from olympics limit 30
```

> <question-title></question-title>
>
> How many unique athletes do we have? Note that you cannot count multiple columns, so choose one that is correct.
>
> > <solution-title></solution-title>
> > 94,733
> >
> > ```sql
> > select count(distinct athlete_id) from olympics;
> > ```
> >
> {: .solution}
{: .question}


# Joining Datasets

This database contains a lot of information, but we may want to add more information. For example if we had a file with information about each country (population, capital city, etc), we could join that information with our Olympics data, to get a single result with all information in every row.

For example, if we would like to be able to group by continent, to e.g. count athletes, medals etc per continent, we will have to add a `continent` column to our file. To do this we would need a second file that maps each country to the corresponding continent. This is what we will do in the next hands-on section.

We obtained country information data from [DataHub](https://datahub.io/core/country-codes). More information about this file can be found in the description there. It has 56 columns with a wide variety of data about each country (from country codes, to capital city, languages spoken, etc)

It is available in the countries table.

> <question-title></question-title>
>
> 2. Which keys(s) in this file are the same as in the `olympics.tsv` file?
>
> > <solution-title></solution-title>
> >
> > 2. Both files have a `NOC` column with the 3-letter country code (`NOC` stands for National Olympic Committee). However, one is lowercase.
> {: .solution}
{: .question}

We would now like to take our Olympics dataset as the basis, and add columns to every row of this file with some information about the country. In order to join, we will need to have one column that is shared between the two files, on which we can match. The `NOC` column is perfect for this because it is a defined standard. Both files also contain a column with the country name in it, which is also a possible candidate to use for joining, but because it is less standardised, it is safer to use the NOC column. For example, if one file uses "Netherlands", while the other uses "The Netherlands" to indicate the same country, the joining will fail for these rows. So always make sure the columns you join on are compatible!

We can use the join commands.

```sql
select * from olympics left join countries on olympics.noc = countries.NOC LIMIT 30
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
> > 3. There is a lot of data duplication in the output now. The exact same country information is added to every line of every athlete from a certain country.
> >    This means much larger response size.
> >    If you do not need all these columns, it could save you a lot of space to select only specific columns that you require.
> >
> {: .solution}
{: .question}

# Concatenating

Concatenation of two files simple means putting the contents of the two files together, one after the other. Our dataset was created in 2021, but since then we've had another Olympic event, the 2022 Winter Olympics in Beijing. If we have the same data for this latest Olympics, we could simply add the rows from the 2022 games to our current file with data, in order to create a single file with all data from 1896 to 2022.

View the table `olympics_2022`, does it have the same structure as our original `olympics` table?

> <question-title></question-title>
>
> 1. Does the new table have the same structure?
> 2. Can we simply add the lines of the new table to the end of our existing olympics dataset?
>
> > <solution-title></solution-title>
> > 1. Yes, this file has all the same columns, in the same order, so concatenation should be relatively straightforward.
> > 2. Yes.
> >
> {: .solution}
{: .question}

Since this new dataset has the exact same structure (number and order of columns), we can simple add the lines from this file to the end of our existing `olympic` table.
For this, we'll need to use the `union all` which takes two separate sql queries and unifies the results.

```sql
select * from olympics
union all
select * from olympics_2022
LIMIT 30;
```

(Note: We are limiting the outputs to ensure your browser does not crash loading all of the data.)

Now this only works so simply because our two datasets had the same structure. If your data comes from different sources, you may have to do some additional data manipulation before you can union, e.g. to make sure the columns match, or how each file deals with missing data (empty cells, `NA`, `NaN` or something else).

# Conclusion

These operations covered in the tutorial are just a few examples of some of the most common operations. There are many more available. We encourage you to look around the documentation of sqlite or your database. The more comfortable you are performing these kinds of steps, the more you can get out of SQL!

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
> >    ```sql
> >    select * from olympics where height is not null
> >    ```
> >
> >    Then we can sort by height, in ascending order to get the shortest athletes on top:
> >
> >    ```sql
> >    select * from olympics where height is not null order by height desc limit 10;
> >    ```
> >
> > 2. We can take the output from the first exercise, and filter for only Winter Olympics:
> >
> >    ```sql
> >    select * from olympics where height is not null and season = 'Winter' order by height desc limit 10;
> >    ```
> >
> > 3. First we filter out the NA values from the weight column:
> >
> >    ```sql
> >    select * from olympics where weight is not null and games = '2020 Summer Olympics' order by weight asc limit 10;
> >    ```
> >
> {: .solution}
>
{: .question}

Congratulations! You have now mastered the basics of data manipulation! There are a lot more data manipulation operations available that you may need. Please explore the tools for yourself, and check back with this tutorial often, we plan to add more sections and exercises over time!
