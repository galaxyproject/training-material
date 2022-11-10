---
layout: tutorial_hands_on
title: Advanced SQL
level: Introductory
zenodo_link:
requirements:
- type: "internal"
  topic_name: data-science
  tutorials:
      - sql-basic
follow_up_training:
- type: "internal"
  topic_name: data-science
  tutorials:
      - sql-python

questions:
- "How can I calculate sums, averages, and other summary values?"
- "How can I combine data from multiple tables?"
- "How should I format data in a database, and why?"
- "How can I create, modify, and delete tables and data?"
- "How can I access databases from programs written in Python?"
objectives:
- "Define aggregation and give examples of its use."
- "Write queries that compute aggregated values."
- "Trace the execution of a query that performs aggregation."
- "Explain how missing data is handled during aggregation."
- "Explain the operation of a query that joins two tables."
- "Explain how to restrict the output of a query containing a join to only include meaningful combinations of values."
- "Write queries that join tables on equal keys."
- "Explain what primary and foreign keys are, and why they are useful."
- "Explain what an atomic value is."
- "Distinguish between atomic and non-atomic values."
- "Explain why every value in a database should be atomic."
- "Explain what a primary key is and why every record should have one."
- "Identify primary keys in database tables."
- "Explain why database entries should not contain redundant information."
- "Identify redundant information in databases."
- "Write statements that create tables."
- "Write statements to insert, modify, and delete records."
- "Write short programs that execute SQL queries."
- "Trace the execution of a program that contains an SQL query."
- "Explain why most database applications are written in a general-purpose language rather than in SQL."
time_estimation:  3H
key_points:
- "Use aggregation functions to combine multiple values."
- "Aggregation functions ignore `null` values."
- "Aggregation happens after filtering."
- "Use GROUP BY to combine subsets separately."
- "If no aggregation function is specified for a field, the query may return an arbitrary value for that field."
- "Use JOIN to combine data from two tables."
- "Use table.field notation to refer to fields when doing joins."
- "Every fact should be represented in a database exactly once."
- "A join produces all combinations of records from one table with records from another."
- "A primary key is a field (or set of fields) whose values uniquely identify the records in a table."
- "A foreign key is a field (or set of fields) in one table whose values are a primary key in another table."
- "We can eliminate meaningless combinations of records by matching primary keys and foreign keys between tables."
- "The most common join condition is matching keys."
- "Every value in a database should be atomic."
- "Every record should have a unique primary key."
- "A database should not contain redundant information."
- "Units and similar metadata should be stored with the data."
- "Use CREATE and DROP to create and delete tables."
- "Use INSERT to add data."
- "Use UPDATE to modify existing data."
- "Use DELETE to remove data."
- "It is simpler and safer to modify data when every record has a unique primary key."
- "Do not create dangling references by deleting records that other records refer to."
- "General-purpose languages have libraries for accessing databases."
- "To connect to a database, a program must use a library specific to that database manager."
- "These libraries use a connection-and-cursor model."
- "Programs can read query results in batches or all at once."
- "Queries should be written using parameter substitution, not string formatting."

contributors:
- carpentries
- hexylena
- avans-atgm

subtopic: sql

notebook:
    language: sql

abbreviations:
    SQL: "Structured Query Language"

tags:
- SQL
---

{% include _includes/quiz.html id="recap.yml" %}

> <comment-title></comment-title>
>
> This tutorial is **significantly** based on [the Carpentries](https://carpentries.org) [Databases and SQL](https://github.com/swcarpentry/sql-novice-survey/) lesson, which is licensed CC-BY 4.0.
>
> Abigail Cabunoc and Sheldon McKay (eds): "Software Carpentry: Using Databases and SQL."  Version 2017.08, August 2017,
> [github.com/swcarpentry/sql-novice-survey](https://github.com/swcarpentry/sql-novice-survey), [https://doi.org/10.5281/zenodo.838776](https://doi.org/10.5281/zenodo.838776)
>
> Adaptations have been made to make this work better in a GTN/Galaxy environment.
{: .comment}


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


```sql
# This preamble sets up the sql "magic" for jupyter. Use %%sql in your cells to write sql!
!python3 -m pip install ipython-sql sqlalchemy
!wget -c http://swcarpentry.github.io/sql-novice-survey/files/survey.db
import sqlalchemy
engine = sqlalchemy.create_engine("sqlite:///survey.db")
%load_ext sql
%sql sqlite:///survey.db
%config SqlMagic.displaycon=False
```

# Aggregation

We now want to calculate ranges and averages for our data.
We know how to select all of the dates from the `Visited` table:

```sql
SELECT dated FROM Visited;
```


but to combine them,
we must use an aggregation function
such as `min` or `max`.
Each of these functions takes a set of records as input,
and produces a single record as output:

```sql
SELECT min(dated) FROM Visited;
```

![SQL Aggregation](../../images/carpentries-sql/sql-aggregation.svg)

```sql
SELECT max(dated) FROM Visited;
```

`min` and `max` are just two of
the aggregation functions built into SQL.
Three others are `avg`,
`count`,
and `sum`:

```sql
SELECT avg(reading) FROM Survey WHERE quant = 'sal';
```

```sql
SELECT count(reading) FROM Survey WHERE quant = 'sal';
```

```sql
SELECT sum(reading) FROM Survey WHERE quant = 'sal';
```

We used `count(reading)` here,
but we could just as easily have counted `quant`
or any other field in the table,
or even used `count(*)`,
since the function doesn't care about the values themselves,
just how many values there are.

SQL lets us do several aggregations at once.
We can,
for example,
find the range of sensible salinity measurements:

```sql
SELECT min(reading), max(reading) FROM Survey WHERE quant = 'sal' AND reading <= 1.0;
```

We can also combine aggregated results with raw results,
although the output might surprise you:

```sql
SELECT person, count(*) FROM Survey WHERE quant = 'sal' AND reading <= 1.0;
```

Why does Lake's name appear rather than Roerich's or Dyer's?
The answer is that when it has to aggregate a field,
but isn't told how to,
the database manager chooses an actual value from the input set.
It might use the first one processed,
the last one,
or something else entirely.

Another important fact is that when there are no values to aggregate ---
for example, where there are no rows satisfying the `WHERE` clause ---
aggregation's result is "don't know"
rather than zero or some other arbitrary value:

```sql
SELECT person, max(reading), sum(reading) FROM Survey WHERE quant = 'missing';
```

One final important feature of aggregation functions is that
they are inconsistent with the rest of SQL in a very useful way.
If we add two values,
and one of them is null,
the result is null.
By extension,
if we use `sum` to add all the values in a set,
and any of those values are null,
the result should also be null.
It's much more useful,
though,
for aggregation functions to ignore null values
and only combine those that are non-null.
This behavior lets us write our queries as:

```sql
SELECT min(dated) FROM Visited;
```

instead of always having to filter explicitly:

```sql
SELECT min(dated) FROM Visited WHERE dated IS NOT NULL;
```

Aggregating all records at once doesn't always make sense.
For example,
suppose we suspect that there is a systematic bias in our data,
and that some scientists' radiation readings are higher than others.
We know that this doesn't work:

```sql
SELECT person, count(reading), round(avg(reading), 2)
FROM  Survey
WHERE quant = 'rad';
```

because the database manager selects a single arbitrary scientist's name
rather than aggregating separately for each scientist.
Since there are only five scientists,
we could write five queries of the form:

```sql
SELECT person, count(reading), round(avg(reading), 2)
FROM  Survey
WHERE quant = 'rad'
AND   person = 'dyer';
```

but this would be tedious,
and if we ever had a data set with fifty or five hundred scientists,
the chances of us getting all of those queries right is small.

What we need to do is
tell the database manager to aggregate the hours for each scientist separately
using a `GROUP BY` clause:

```sql
SELECT   person, count(reading), round(avg(reading), 2)
FROM     Survey
WHERE    quant = 'rad'
GROUP BY person;
```

`GROUP BY` does exactly what its name implies:
groups all the records with the same value for the specified field together
so that aggregation can process each batch separately.
Since all the records in each batch have the same value for `person`,
it no longer matters that the database manager
is picking an arbitrary one to display
alongside the aggregated `reading` values.

> <tip-title>Know Excel? It's just a pivot table.</tip-title>
> `GROUP BY` is basically just a pivot table for Excel users, it lets you build
> nice summary tables which aggregate your results.
>
> And if you didn't already know the Excel equivalent, now you know what to
> look for when you need it!
{: .tip}

Just as we can sort by multiple criteria at once,
we can also group by multiple criteria.
To get the average reading by scientist and quantity measured,
for example,
we just add another field to the `GROUP BY` clause:

```sql
SELECT   person, quant, count(reading), round(avg(reading), 2)
FROM     Survey
GROUP BY person, quant;
```

Note that we have added `quant` to the list of fields displayed,
since the results wouldn't make much sense otherwise.

Let's go one step further and remove all the entries
where we don't know who took the measurement:

```sql
SELECT   person, quant, count(reading), round(avg(reading), 2)
FROM     Survey
WHERE    person IS NOT NULL
GROUP BY person, quant
ORDER BY person, quant;
```


Looking more closely,
this query:

1.  selected records from the `Survey` table where the `person` field was not null;
2.  grouped those records into subsets so that the `person` and `quant` values in each subset were the same;
3.  ordered those subsets first by `person`, and then within each sub-group by `quant`; and
4.  counted the number of records in each subset, calculated the average `reading` in each, and chose a `person` and `quant` value from each (it doesn't matter which ones, since they're all equal).

> <question-title>Counting Temperature Readings</question-title>
>
> How many temperature readings did Frank Pabodie record,
> and what was their average value?
>
> > <solution-title></solution-title>
> >
> > ```
> > SELECT count(reading), avg(reading) FROM Survey WHERE quant = 'temp' AND person = 'pb';
> > ```
> >
> > |count(reading)|avg(reading)|
> > |--------------|------------|
> > |2             |-20.0       |
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

> ## Averaging with NULL
>
> The average of a set of values is the sum of the values
> divided by the number of values.
> Does this mean that the `avg` function returns 2.0 or 3.0
> when given the values 1.0, `null`, and 5.0?
>
> > <solution-title></solution-title>
> > The answer is 3.0.
> > `NULL` is not a value; it is the absence of a value.
> > As such it is not included in the calculation.
> >
> > You can confirm this, by executing this code:
> >
> > ```
> > SELECT AVG(a) FROM (
> >     SELECT 1 AS a
> >     UNION ALL SELECT NULL
> >     UNION ALL SELECT 5);
> > ```
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

> <question-title>What Does This Query Do?</question-title>
>
> We want to calculate the difference between
> each individual radiation reading
> and the average of all the radiation readings.
> We write the query:
>
> ```
> SELECT reading - avg(reading) FROM Survey WHERE quant = 'rad';
> ```
>
> What does this actually produce, and can you think of why?
>
> > <solution-title></solution-title>
> > The query produces only one row of results when we what we really want is a result for each of the readings.
> > The `avg()` function produces only a single value, and because it is run first, the table is reduced to a single row.
> > The `reading` value is simply an arbitrary one.
> >
> > To achieve what we wanted, we would have to run two queries:
> >
> > ```
> > SELECT avg(reading) FROM Survey WHERE quant='rad';
> > ```
> >
> > This produces the average value (6.5625), which we can then insert into a second query:
> >
> > ```
> > SELECT reading - 6.5625 FROM Survey WHERE quant = 'rad';
> > ```
> >
> > This produces what we want, but we can combine this into a single query using subqueries.
> >
> > ```
> > SELECT reading - (SELECT avg(reading) FROM Survey WHERE quant='rad') FROM Survey WHERE quant = 'rad';
> > ```
> >
> > This way we don't have execute two queries.
> >
> > In summary what we have done is to replace `avg(reading)` with `(SELECT avg(reading) FROM Survey WHERE quant='rad')` in the original query.
> >
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

> <question-title>Ordering When Concatenating</question-title>
>
> The function `group_concat(field, separator)`
> concatenates all the values in a field
> using the specified separator character
> (or ',' if the separator isn't specified).
> Use this to produce a one-line list of scientists' names,
> such as:
>
> ```
> William Dyer, Frank Pabodie, Anderson Lake, Valentina Roerich, Frank Danforth
> ```
>
> Can you find a way to order the list by surname?
{: .question}

```sql
-- Try solutions here!
```


# Combining Data

In order to submit our data to a web site
that aggregates historical meteorological data,
we might need to format it as
latitude, longitude, date, quantity, and reading.
However,
our latitudes and longitudes are in the `Site` table,
while the dates of measurements are in the `Visited` table
and the readings themselves are in the `Survey` table.
We need to combine these tables somehow.

This figure shows the relations between the tables:

![Survey Database Structure](../../images/carpentries-sql/sql-join-structure.svg)

The SQL command to do this is `JOIN`.
To see how it works,
let's start by joining the `Site` and `Visited` tables:

```sql
SELECT * FROM Site JOIN Visited;
```

`JOIN` creates
the cross product
of two tables,
i.e.,
it joins each record of one table with each record of the other table
to give all possible combinations.
Since there are three records in `Site`
and eight in `Visited`,
the join's output has 24 records (3 * 8 = 24) .
And since each table has three fields,
the output has six fields (3 + 3 = 6).

What the join *hasn't* done is
figure out if the records being joined have anything to do with each other.
It has no way of knowing whether they do or not until we tell it how.
To do that,
we add a clause specifying that
we're only interested in combinations that have the same site name,
thus we need to use a filter:

```sql
SELECT * FROM Site JOIN Visited ON Site.name = Visited.site;
```

`ON` is very similar to `WHERE`,
and for all the queries in this lesson you can use them interchangeably.
There are differences in how they affect [outer joins][outer],
but that's beyond the scope of this lesson.
Once we add this to our query,
the database manager throws away records
that combined information about two different sites,
leaving us with just the ones we want.

Notice that we used `Table.field` to specify field names
in the output of the join.
We do this because tables can have fields with the same name,
and we need to be specific which ones we're talking about.
For example,
if we joined the `Person` and `Visited` tables,
the result would inherit a field called `id`
from each of the original tables.

We can now use the same dotted notation
to select the three columns we actually want
out of our join:

```sql
SELECT Site.lat, Site.long, Visited.dated
FROM   Site JOIN Visited
ON     Site.name = Visited.site;
```


If joining two tables is good,
joining many tables must be better.
In fact,
we can join any number of tables
simply by adding more `JOIN` clauses to our query,
and more `ON` tests to filter out combinations of records
that don't make sense:

```sql
SELECT Site.lat, Site.long, Visited.dated, Survey.quant, Survey.reading
FROM   Site JOIN Visited JOIN Survey
ON     Site.name = Visited.site
AND    Visited.id = Survey.taken
AND    Visited.dated IS NOT NULL;
```

We can tell which records from `Site`, `Visited`, and `Survey`
correspond with each other
because those tables contain
primary keys
and foreign keys.
A primary key is a value,
or combination of values,
that uniquely identifies each record in a table.
A foreign key is a value (or combination of values) from one table
that identifies a unique record in another table.
Another way of saying this is that
a foreign key is the primary key of one table
that appears in some other table.
In our database,
`Person.id` is the primary key in the `Person` table,
while `Survey.person` is a foreign key
relating the `Survey` table's entries
to entries in `Person`.

Most database designers believe that
every table should have a well-defined primary key.
They also believe that this key should be separate from the data itself,
so that if we ever need to change the data,
we only need to make one change in one place.
One easy way to do this is
to create an arbitrary, unique ID for each record
as we add it to the database.
This is actually very common:
those IDs have names like "student numbers" and "patient numbers",
and they almost always turn out to have originally been
a unique record identifier in some database system or other.
As the query below demonstrates,
SQLite [automatically numbers records][rowid] as they're added to tables,
and we can use those record numbers in queries:

```sql
SELECT rowid, * FROM Person;
```

> <question-title>Listing Radiation Readings</question-title>
>
> Write a query that lists all radiation readings from the DR-1 site.
> > <solution-title></solution-title>
> >
> > ```
> > SELECT Survey.reading
> > FROM Site JOIN Visited JOIN Survey
> > ON Site.name = Visited.site
> > AND Visited.id = Survey.taken
> > WHERE Site.name = 'DR-1'
> > AND Survey.quant = 'rad';
> > ```
> >
> > |reading   |
> > |----------|
> > |9.82      |
> > |7.8       |
> > |11.25     |
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

> <question-title>Where's Frank?</question-title>
>
> Write a query that lists all sites visited by people named "Frank".
> > <solution-title></solution-title>
> >
> > ```
> > SELECT DISTINCT Site.name
> > FROM Site JOIN Visited JOIN Survey JOIN Person
> > ON Site.name = Visited.site
> > AND Visited.id = Survey.taken
> > AND Survey.person = Person.id
> > WHERE Person.personal = 'Frank';
> > ```
> >
> > |name   |
> > |-------|
> > |DR-3   |
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

> <question-title>Reading Queries</question-title>
>
> Describe in your own words what the following query produces:
>
> ```
> SELECT Site.name FROM Site JOIN Visited
> ON Site.lat < -49.0 AND Site.name = Visited.site AND Visited.dated >= '1932-01-01';
> ```
{: .question}

```sql
-- Try solutions here!
```

> <question-title>Who Has Been Where?</question-title>
>
> Write a query that shows each site with exact location (lat, long) ordered by visited date,
> followed by personal name and family name of the person who visited the site
> and the type of measurement taken and its reading. Please avoid all null values.
> Tip: you should get 15 records with 8 fields.
> > <solution-title></solution-title>
> >
> > ```
> > SELECT Site.name, Site.lat, Site.long, Person.personal, Person.family, Survey.quant, Survey.reading, Visited.dated
> > FROM Site JOIN Visited JOIN Survey JOIN Person
> > ON Site.name = Visited.site
> > AND Visited.id = Survey.taken
> > AND Survey.person = Person.id
> > WHERE Survey.person IS NOT NULL
> > AND Visited.dated IS NOT NULL
> > ORDER BY Visited.dated;
> > ```
> >
> > name   |  lat        |  long       |  personal   | family   | quant     | reading   |     dated
> >--------|-------------|-------------|-------------|----------|-----------|-----------|-----------
> >DR-1    |    -49.85   |   -128.57   |  William    | Dyer     |   rad     |    9.82   |   1927-02-08
> >DR-1    |    -49.85   |   -128.57   |  William    | Dyer     |   sal     |    0.13   |   1927-02-08
> >DR-1    |    -49.85   |   -128.57   |  William    | Dyer     |   rad     |    7.8    |   1927-02-10
> >DR-1    |    -49.85   |   -128.57   |  William    | Dyer     |   sal     |    0.09   |   1927-02-10
> >DR-3    |    -47.15   |   -126.72   |  Anderson   | Lake     |   sal     |    0.05   |   1930-01-07
> >DR-3    |    -47.15   |   -126.72   |  Frank      | Pabodie  |   rad     |    8.41   |   1930-01-07
> >DR-3    |    -47.15   |   -126.72   |  Frank      | Pabodie  |   temp    |    -21.5  |   1930-01-07
> >DR-3    |    -47.15   |   -126.72   |  Frank      | Pabodie  |   rad     |    7.22   |   1930-01-12
> >DR-3    |    -47.15   |   -126.72   |  Anderson   | Lake     |   sal     |    0.1    |   1930-02-26
> >DR-3    |    -47.15   |   -126.72   |  Frank      | Pabodie  |   rad     |    4.35   |   1930-02-26
> >DR-3    |    -47.15   |   -126.72   |  Frank      | Pabodie  |   temp    |    -18.5  |   1930-02-26
> >MSK-4   |    -48.87   |   -123.4    |  Anderson   | Lake     |   rad     |    1.46   |   1932-01-14
> >MSK-4   |    -48.87   |   -123.4    |  Anderson   | Lake     |   sal     |    0.21   |   1932-01-14
> >MSK-4   |    -48.87   |   -123.4    |  Valentina  | Roerich  |   sal     |    22.5   |   1932-01-14
> >DR-1    |    -49.85   |   -128.57   |  Valentina  | Roerich  |   rad     |    11.25  |   1932-03-22
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

A good visual explanation of joins can be found [in the SQL Join Visualizer][joinref]

[outer]: https://en.wikipedia.org/wiki/Join_%28SQL%29#Outer_join
[rowid]: https://www.sqlite.org/lang_createtable.html#rowid
[joinref]: https://sql-joins.leopard.in.ua/

# Data Hygiene

Now that we have seen how joins work, we can see why the relational
model is so useful and how best to use it.  The first rule is that
every value should be atomic, i.e., not
contain parts that we might want to work with separately.  We store
personal and family names in separate columns instead of putting the
entire name in one column so that we don't have to use substring
operations to get the name's components.  More importantly, we store
the two parts of the name separately because splitting on spaces is
unreliable: just think of a name like "Eloise St. Cyr" or "Jan Mikkel
Steubart".

The second rule is that every record should have a unique primary key.
This can be a serial number that has no intrinsic meaning,
one of the values in the record (like the `id` field in the `Person` table),
or even a combination of values:
the triple `(taken, person, quant)` from the `Survey` table uniquely identifies every measurement.

The third rule is that there should be no redundant information.
For example,
we could get rid of the `Site` table and rewrite the `Visited` table like this:

|id   |lat   |long   |dated      |
|-----|------|-------|-----------|
|619  |-49.85|-128.57| 1927-02-08|
|622  |-49.85|-128.57| 1927-02-10|
|734  |-47.15|-126.72| 1930-01-07|
|735  |-47.15|-126.72| 1930-01-12|
|751  |-47.15|-126.72| 1930-02-26|
|752  |-47.15|-126.72| None    |
|837  |-48.87|-123.40| 1932-01-14|
|844  |-49.85|-128.57| 1932-03-22|

In fact,
we could use a single table that recorded all the information about each reading in each row,
just as a spreadsheet would.
The problem is that it's very hard to keep data organized this way consistent:
if we realize that the date of a particular visit to a particular site is wrong,
we have to change multiple records in the database.
What's worse,
we may have to guess which records to change,
since other sites may also have been visited on that date.

The fourth rule is that the units for every value should be stored explicitly.
Our database doesn't do this,
and that's a problem:
Roerich's salinity measurements are several orders of magnitude larger than anyone else's,
but we don't know if that means she was using parts per million instead of parts per thousand,
or whether there actually was a saline anomaly at that site in 1932.

Stepping back,
data and the tools used to store it have a symbiotic relationship:
we use tables and joins because it's efficient,
provided our data is organized a certain way,
but organize our data that way because we have tools to manipulate it efficiently.
As anthropologists say,
the tool shapes the hand that shapes the tool.

> <question-title>Identifying Atomic Values</question-title>
>
> Which of the following are atomic values? Which are not? Why?
>
> *   New Zealand
> *   87 Turing Avenue
> *   January 25, 1971
> *   the XY coordinate (0.5, 3.3)
>
> > <solution-title></solution-title>
> > New Zealand is the only clear-cut atomic value.
> >
> > The address and the XY coordinate contain more than one piece of information
> > which should be stored separately:
> > - House number, street name
> > - X coordinate, Y coordinate
> >
> > The date entry is less clear cut, because it contains month, day, and year elements.
> > However, there is a `DATE` datatype in SQL, and dates should be stored using this format.
> > If we need to work with the month, day, or year separately, we can use the SQL functions available for our database software
> > (for example [`EXTRACT`](https://docs.oracle.com/cd/B19306_01/server.102/b14200/functions050.htm) or [`STRFTIME`](http://www.sqlite.org/lang_datefunc.html) for SQLite).
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

> <question-title>Identifying a Primary Key</question-title>
>
> What is the primary key in this table?
> I.e., what value or combination of values uniquely identifies a record?
>
> |latitude|longitude|date      |temperature|
> |--------|---------|----------|-----------|
> |57.3    |-22.5    |2015-01-09|-14.2      |
>
> > <solution-title></solution-title>
> > Latitude, longitude, and date are all required to uniquely identify the temperature record.
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

# Creating and Modifying Data

So far we have only looked at how to get information out of a database,
both because that is more frequent than adding information,
and because most other operations only make sense
once queries are understood.
If we want to create and modify data,
we need to know two other sets of commands.

The first pair are [`CREATE TABLE`][create-table] and [`DROP TABLE`][drop-table].
While they are written as two words,
they are actually single commands.
The first one creates a new table;
its arguments are the names and types of the table's columns.
For example,
the following statements create the four tables in our survey database:

```sql
CREATE TABLE Person(id text, personal text, family text);
CREATE TABLE Site(name text, lat real, long real);
CREATE TABLE Visited(id integer, site text, dated text);
CREATE TABLE Survey(taken integer, person text, quant text, reading real);
```

We can get rid of one of our tables using:

```sql
DROP TABLE Survey;
```

Be very careful when doing this:
if you drop the wrong table, hope that the person maintaining the database has a backup,
but it's better not to have to rely on it.

Different database systems support different data types for table columns,
but most provide the following:

|data type|  use                                       |
|---------|  ----------------------------------------- |
|INTEGER  |  a signed integer                          |
|REAL     |  a floating point number                   |
|TEXT     |  a character string                        |
|BLOB     |  a "binary large object", such as an image |

Most databases also support Booleans and date/time values;
SQLite uses the integers 0 and 1 for the former,
and represents the latter as text or numeric fields.

An increasing number of databases also support geographic data types,
such as latitude and longitude.
Keeping track of what particular systems do or do not offer,
and what names they give different data types,
is an unending portability headache.

> <tip-title>Which database should I use?</tip-title>
> SQLite is fantastic for small databases or embedded into applications where
> you want to be able to use SQL to query and process data.
>
> However for any real analysis PostgreSQL is usually the best choice, it
> scales incredibly well and can meet a wide range of use cases. It has good
> data type support.
{: .tip}

> <tip-title>Do you have geographic data?</tip-title>
> Use Postgres. The [PostGIS](https://postgis.net/) library is fantastic and industry standard for storing geographic data in a database.
{: .tip}

When we create a table,
we can specify several kinds of constraints on its columns.
For example,
a better definition for the `Survey` table would be:

```sql
CREATE TABLE Survey(
    taken   integer not null, -- where reading taken
    person  text,             -- may not know who took it
    quant   text not null,    -- the quantity measured
    reading real not null,    -- the actual reading
    primary key(taken, quant),
    foreign key(taken) references Visited(id),
    foreign key(person) references Person(id)
);
```

Once again,
exactly what constraints are available
and what they're called
depends on which database manager we are using.

Once tables have been created,
we can add, change, and remove records using our other set of commands,
`INSERT`, `UPDATE`, and `DELETE`.

Here is an example of inserting rows into the `Site` table:

```sql
INSERT INTO Site (name, lat, long) VALUES ('DR-1', -49.85, -128.57);
INSERT INTO Site (name, lat, long) VALUES ('DR-3', -47.15, -126.72);
INSERT INTO Site (name, lat, long) VALUES ('MSK-4', -48.87, -123.40);
```

We can also insert values into one table directly from another:

```sql
CREATE TABLE JustLatLong(lat real, long real);
INSERT INTO JustLatLong SELECT lat, long FROM Site;
```

Modifying existing records is done using the `UPDATE` statement.
To do this we tell the database which table we want to update,
what we want to change the values to for any or all of the fields,
and under what conditions we should update the values.

For example, if we made a mistake when entering the lat and long values
of the last `INSERT` statement above, we can correct it with an update:

```sql
UPDATE Site SET lat = -47.87, long = -122.40 WHERE name = 'MSK-4';
```

Be careful to not forget the `WHERE` clause or the update statement will
modify *all* of the records in the database.

Deleting records can be a bit trickier,
because we have to ensure that the database remains internally consistent.
If all we care about is a single table,
we can use the `DELETE` command with a `WHERE` clause
that matches the records we want to discard.
For example,
once we realize that Frank Danforth didn't take any measurements,
we can remove him from the `Person` table like this:

```sql
DELETE FROM Person WHERE id = 'danforth';
```

But what if we removed Anderson Lake instead?
Our `Survey` table would still contain seven records
of measurements he'd taken,
but that's never supposed to happen:
`Survey.person` is a foreign key into the `Person` table,
and all our queries assume there will be a row in the latter
matching every value in the former.

This problem is called referential integrity:
we need to ensure that all references between tables can always be resolved correctly.
One way to do this is to delete all the records
that use `'lake'` as a foreign key
before deleting the record that uses it as a primary key.
If our database manager supports it,
we can automate this
using cascading delete.
However,
this technique is outside the scope of this chapter.

> <tip-title>Hybrid Storage Models</tip-title>
>
> Many applications use a hybrid storage model
> instead of putting everything into a database:
> the actual data (such as astronomical images) is stored in files,
> while the database stores the files' names,
> their modification dates,
> the region of the sky they cover,
> their spectral characteristics,
> and so on.
> This is also how most music player software is built:
> the database inside the application keeps track of the MP3 files,
> but the files themselves live on disk.
{: .tip}

> <question-title>Replacing NULL</question-title>
>
> Write an SQL statement to replace all uses of `null` in
> `Survey.person` with the string `'unknown'`.
>
> > <solution-title></solution-title>
> > ```
> > UPDATE Survey SET person = 'unknown' WHERE person IS NULL;
> > ```
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

> <question-title>Backing Up with SQL</question-title>
>
> SQLite has several administrative commands that aren't part of the
> SQL standard.  One of them is `.dump`, which prints the SQL commands
> needed to re-create the database.  Another is `.read`, which reads a
> file created by `.dump` and restores the database.  A colleague of
> yours thinks that storing dump files (which are text) in version
> control is a good way to track and manage changes to the database.
> What are the pros and cons of this approach?  (Hint: records aren't
> stored in any particular order.)
>
> > <solution-title></solution-title>
> > #### Advantages
> > - A version control system will be able to show differences between versions
> > of the dump file; something it can't do for binary files like databases
> > - A VCS only saves changes between versions, rather than a complete copy of
> > each version (save disk space)
> > - The version control log will explain the reason for the changes in each version
> > of the database
> >
> > #### Disadvantages
> > - Artificial differences between commits because records don't have a fixed order
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

[create-table]: https://www.sqlite.org/lang_createtable.html
[drop-table]: https://www.sqlite.org/lang_droptable.html
