---
layout: tutorial_hands_on
title: Introduction to SQL
level: Introductory
zenodo_link:
requirements:  []
follow_up_training:
- type: "internal"
  topic_name: data-science
  tutorials:
      - sql-advanced

questions:
- "How can I get data from a database?"
- "How can I sort a query's results?"
- "How can I remove duplicate values from a query's results?"
- "How can I select subsets of data?"
- "How can I calculate new values on the fly?"
- "How do databases represent missing information?"
- "What special handling does missing information require?"

objectives:
- "Explain the difference between a table, a record, and a field."
- "Explain the difference between a database and a database manager."
- "Write a query to select all values for specific fields from a single table."
- "Write queries that display results in a particular order."
- "Write queries that eliminate duplicate values from data."
- "Write queries that select records that satisfy user-specified conditions."
- "Explain the order in which the clauses in a query are executed."
- "Write queries that calculate new values for each selected record."
- "Explain how databases represent missing information."
- "Explain the three-valued logic databases use when manipulating missing information."
- "Write queries that handle missing information correctly."

time_estimation:  3H
key_points:
- "A relational database stores information in tables, each of which has a fixed set of columns and a variable number of records."
- "A database manager is a program that manipulates information stored in a database."
- "We write queries in a specialized language called SQL to extract information from databases."
- "Use SELECT... FROM... to get values from a database table."
- "SQL is case-insensitive (but data is case-sensitive)."
- "The records in a database table are not intrinsically ordered: if we want to display them in some order, we must specify that explicitly with ORDER BY."
- "The values in a database are not guaranteed to be unique: if we want to eliminate duplicates, we must specify that explicitly as well using DISTINCT."
- "Use WHERE to specify conditions that records must meet in order to be included in a query's results."
- "Use AND, OR, and NOT to combine tests."
- "Filtering is done on whole records, so conditions can use fields that are not actually displayed."
- "Write queries incrementally."
- "Queries can do the usual arithmetic operations on values."
- "Use UNION to combine the results of two or more queries."
- "Databases use a special value called NULL to represent missing information."
- "Almost all operations on NULL produce NULL."
- "Queries can test for NULLs using IS NULL and IS NOT NULL."

contributors:
- carpentries
- hexylena
- dirowa
- bazante1
- avans-atgm

subtopic: sql

notebook:
    language: sql

abbreviations:
    SQL: "Structured Query Language"

tags:
- SQL
---

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

# Selecting Data

A relational database
is a way to store and manipulate information.
Databases are arranged as table.
Each table has columns (also known as fields) that describe the data,
and rows (also known as records) which contain the data.

When we are using a spreadsheet,
we put formulas into cells to calculate new values based on old ones.
When we are using a database,
we send commands
(usually called queries)
to a database manager:
a program that manipulates the database for us.
The database manager does whatever lookups and calculations the query specifies,
returning the results in a tabular form
that we can then use as a starting point for further queries.

Queries are written in a language called {SQL},
SQL provides hundreds of different ways to analyze and recombine data.
We will only look at a handful of queries,
but that handful accounts for most of what scientists do.

> <tip-title>Changing Database Managers</tip-title>
>
> Many database managers --- Oracle,
> IBM DB2, PostgreSQL, MySQL, Microsoft Access, and SQLite ---  understand
> SQL but each stores data in a different way,
> so a database created with one cannot be used directly by another.
> However, every database manager
> can import and export data in a variety of formats like .csv, SQL,
> so it *is* possible to move information from one to another.
{: .tip}

Before we get into using {SQL} to select the data, let's take a look at the tables of the database we will use in our examples:

**Person**: people who took readings.

|id      |personal |family
|--------|---------|----------
|dyer    |William  |Dyer
|pb      |Frank    |Pabodie
|lake    |Anderson |Lake
|roe     |Valentina|Roerich
|danforth|Frank    |Danforth

**Site**: locations where readings were taken.

|name |lat   |long   |
|-----|------|-------|
|DR-1 |-49.85|-128.57|
|DR-3 |-47.15|-126.72|
|MSK-4|-48.87|-123.4 |

**Visited**: when readings were taken at specific sites.

|id   |site |dated     |
|-----|-----|----------|
|619  |DR-1 |1927-02-08|
|622  |DR-1 |1927-02-10|
|734  |DR-3 |1930-01-07|
|735  |DR-3 |1930-01-12|
|751  |DR-3 |1930-02-26|
|752  |DR-3 |None      |
|837  |MSK-4|1932-01-14|
|844  |DR-1 |1932-03-22|

**Survey**: the actual readings.  The field `quant` is short for quantitative and indicates what is being measured.  Values are `rad`, `sal`, and `temp` referring to 'radiation', 'salinity' and 'temperature', respectively.

|taken|person|quant|reading|
|-----|------|-----|-------|
|619  |dyer  |rad  |9.82   |
|619  |dyer  |sal  |0.13   |
|622  |dyer  |rad  |7.8    |
|622  |dyer  |sal  |0.09   |
|734  |pb    |rad  |8.41   |
|734  |lake  |sal  |0.05   |
|734  |pb    |temp |-21.5  |
|735  |pb    |rad  |7.22   |
|735  |None  |sal  |0.06   |
|735  |None  |temp |-26.0  |
|751  |pb    |rad  |4.35   |
|751  |pb    |temp |-18.5  |
|751  |lake  |sal  |0.1    |
|752  |lake  |rad  |2.19   |
|752  |lake  |sal  |0.09   |
|752  |lake  |temp |-16.0  |
|752  |roe   |sal  |41.6   |
|837  |lake  |rad  |1.46   |
|837  |lake  |sal  |0.21   |
|837  |roe   |sal  |22.5   |
|844  |roe   |rad  |11.25  |

Notice that three entries --- one in the `Visited` table,
and two in the `Survey` table --- don't contain any actual
data, but instead have a special `None` entry:
we'll return to these missing values.

For now,
let's write an SQL query that displays scientists' names.
We do this using the SQL command `SELECT`,
giving it the names of the columns we want and the table we want them from.
Our query and its output look like this:

```sql
SELECT family, personal FROM Person;
```

The semicolon at the end of the query
tells the database manager that the query is complete and ready to run.
We have written our commands in upper case and the names for the table and columns
in lower case,
but we don't have to:
as the example below shows,
SQL is case insensitive.

```sql
SeLeCt FaMiLy, PeRsOnAl FrOm PeRsOn;
```

You can use SQL's case insensitivity to your advantage. For instance,
some people choose to write SQL keywords (such as `SELECT` and `FROM`)
in capital letters and **field** and **table** names in lower
case. This can make it easier to locate parts of an SQL statement. For
instance, you can scan the statement, quickly locate the prominent
`FROM` keyword and know the table name follows.  Whatever casing
convention you choose, please be consistent: complex queries are hard
enough to read without the extra cognitive load of random
capitalization.  One convention is to use UPPER CASE for SQL
statements, to distinguish them from tables and column names. This is
the convention that we will use for this lesson.

> <question-title>Is a personal and family name column a good design?</question-title>
> If you were tasked with designing a database to store this same data, is storing the name data in
> this way the best way to do it? Why or why not?
>
> Can you think of any names that would be difficult to enter in such a schema?
>
> > <solution-title></solution-title>
> > No, it is generally not. There are a lot of [falsehoods that programmers believe about names](https://shinesolutions.com/2018/01/08/falsehoods-programmers-believe-about-names-with-examples/).
> > The situation is much more complex as you can read in that article, but names vary wildly and
> > generally placing constraints on how names are entered is only likely to frustrate you or your
> > users later on when they need to enter data into that database.
> >
> > In general you should consider using a single text field for the name and allowing users to
> > specify them as whatever they like (if it is a system with registration), or asking what they
> > wish to be recorded (if you are doing this sort of data collection).
> >
> > If you are doing scientific research, you might know that names are generally very poor
> > identifiers of a single human, and in that case consider recording their
> > [ORCiD](https://orcid.org/) which will help you reference that individual when you are
> > publishing later.
> >
> > This is also a good time to consider what data you really *need* to collect. If you are working
> > in the EU under GDPR, do you really need their full legal name? Is that necessary? Do you have a
> > plan for ensuring that data is correct when publishing, if any part of their name has changed
> > since?
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

While we are on the topic of SQL's syntax, one aspect of SQL's syntax
that can frustrate novices and experts alike is forgetting to finish a
command with `;` (semicolon).  When you press enter for a command
without adding the `;` to the end, it can look something like this:

~~~
SELECT id FROM Person
...>
...>
~~~

This is SQL's prompt, where it is waiting for additional commands or
for a `;` to let SQL know to finish.  This is easy to fix!  Just type
`;` and press enter!

Now, going back to our query,
it's important to understand that
the rows and columns in a database table aren't actually stored in any particular order.
They will always be *displayed* in some order,
but we can control that in various ways.
For example,
we could swap the columns in the output by writing our query as:

```sql
SELECT personal, family FROM Person;
```

or even repeat columns:

```sql
SELECT id, id, id FROM Person;
```

As a shortcut,
we can select all of the columns in a table using `*`:

```sql
SELECT * FROM Person;
```

> <question-title>Selecting Site Names</question-title>
>
> Write a query that selects only the `name` column from the `Site` table.
>
> > <solution-title></solution-title>
> >
> > ~~~
> > SELECT name FROM Site;
> > ~~~
> >
> > |name      |
> > |----------|
> > |DR-1      |
> > |DR-3      |
> > |MSK-4     |
> {: .solution}
{: .question}


```sql
-- Try solutions here!
```

> <question-title>Query Style</question-title>
>
> Many people format queries as:
>
> ~~~
> SELECT personal, family FROM person;
> ~~~
>
> or as:
>
> ~~~
> select Personal, Family from PERSON;
> ~~~
>
> What style do you find easiest to read, and why?
{: .question}

```sql
-- Try solutions here!
```


# Sorting and Removing Duplicates

In beginning our examination of the Antarctic data, we want to know:

* what kind of quantity measurements were taken at each site;
* which scientists took measurements on the expedition;

To determine which measurements were taken at each site,
we can examine the `Survey` table.
Data is often redundant,
so queries often return redundant information.
For example,
if we select the quantities that have been measured
from the `Survey` table,
we get this:

```sql
SELECT quant FROM Survey;
```


This result makes it difficult to see all of the different types of
`quant` in the Survey table.  We can eliminate the redundant output to
make the result more readable by adding the `DISTINCT` keyword to our
query:

```sql
SELECT DISTINCT quant FROM Survey;
```

If we want to determine which visit (stored in the `taken` column)
have which `quant` measurement,
we can use the `DISTINCT` keyword on multiple columns.
If we select more than one column,
distinct *sets* of values are returned
(in this case *pairs*, because we are selecting two columns):

```sql
SELECT DISTINCT taken, quant FROM Survey;
```

Notice in both cases that duplicates are removed
even if the rows they come from didn't appear to be adjacent in the database table.


Our next task is to identify the scientists on the expedition by looking at the `Person` table.
As we mentioned earlier,
database records are not stored in any particular order.
This means that query results aren't necessarily sorted,
and even if they are,
we often want to sort them in a different way,
e.g., by their identifier instead of by their personal name.
We can do this in SQL by adding an `ORDER BY` clause to our query:

```sql
SELECT * FROM Person ORDER BY id;
```

|id     |personal |family  |
|-------|---------|--------|
|danfort|Frank    |Danforth|
|dyer   |William  |Dyer    |
|lake   |Anderson |Lake    |
|pb     |Frank    |Pabodie |
|roe    |Valentina|Roerich |

By default, when we use `ORDER BY`,
results are sorted in ascending order of the column we specify
(i.e.,
from least to greatest).

We can sort in the opposite order using `DESC` (for "descending"):

> <tip-title>A note on ordering</tip-title>
>
> While it may look that the records are consistent every time we ask for them in this lesson, that is because no one has changed or modified any of the data so far. Remember to use `ORDER BY` if you want the rows returned to have any sort of consistent or predictable order.
{: .tip}
```sql
SELECT * FROM person ORDER BY id DESC;
```

(And if we want to make it clear that we're sorting in ascending order,
we can use `ASC` instead of `DESC`.)


In order to look at which scientist measured quantities during each visit,
we can look again at the `Survey` table.
We can also sort on several fields at once.
For example,
this query sorts results first in ascending order by `taken`,
and then in descending order by `person`
within each group of equal `taken` values:

```sql
SELECT taken, person, quant FROM Survey ORDER BY taken ASC, person DESC;
```

This query gives us a good idea of which scientist was involved in which visit,
and what measurements they performed during the visit.

Looking at the table, it seems like some scientists specialized in
certain kinds of measurements.  We can examine which scientists
performed which measurements by selecting the appropriate columns and
removing duplicates.

```sql
SELECT DISTINCT quant, person FROM Survey ORDER BY quant ASC;
```

> <question-title>Finding Distinct Dates</question-title>
>
> Write a query that selects distinct dates from the `Visited` table.
>
> > <solution-title></solution-title>
> >
> > ~~~
> > SELECT DISTINCT dated FROM Visited;
> > ~~~
> >
> > |dated     |
> > |----------|
> > |1927-02-08|
> > |1927-02-10|
> > |1930-01-07|
> > |1930-01-12|
> > |1930-02-26|
> > |&nbsp;    |
> > |1932-01-14|
> > |1932-03-22|
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

> <question-title>Displaying Full Names</question-title>
>
> Write a query that displays the full names of the scientists in the `Person` table,
> ordered by family name.
>
> > <solution-title></solution-title>
> >
> > ~~~
> > SELECT personal, family FROM Person ORDER BY family ASC;
> > ~~~
> >
> > |personal  |family    |
> > |----------|----------|
> > |Frank     |Danforth  |
> > |William   |Dyer      |
> > |Anderson  |Lake      |
> > |Frank     |Pabodie   |
> > |Valentina |Roerich   |
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

> <tip-title>Is sorting names useful?</tip-title>
> If you are someone with a name which falls at the end of the alphabet, you've likely been
> penalised for this your entire life. Alphabetically sorting names should always be looked at
> critically and through a lens to whether you are fairly reflecting everyone's contributions,
> rather than just the default sort order.
>
> There are many options, either by some metric of contribution that everyone could agree on, or
> better, consider random sorting, like the GTN uses with our [Hall of Fame]({% link hall-of-fame.md %})
> page where we intentionally order randomly to tell contributors that no one persons
> contributions matter more than anothers.
>
> > The evidence provided in a variety of studies leaves no doubt that an
> > alphabetical author ordering norm disadvantages researchers with
> > last names toward the end of the alphabet. There is furthermore con-
> > vincing evidence that researchers are aware of this and that they
> > react strategically to such alphabetical discrimination, for example
> > with their choices of who to collaborate with. See {% cite Weber_2018 %} for more.
> {: .quote}
{: .tip}

> <tip-title>Name collation</tip-title>
> When you are sorting things in SQL, you need to be aware of something called collation which can
> affect your results if you have values that are not the letters A-Z. Collating is the process of
> sorting values, and this affects many human languages when storing data in a database.
>
> Here is a Dutch example. In the old days their alphabet contained a `ÿ` which was later replaced
> with `ĳ`, a digraph of two characters squished together. This is commonly rendered as `ij`
> however, two separate characters, due to the internet and widespread use of keyboards featuring
> mainly ascii characters. However, it is still the 25th letter of their alphabet.
>
> ```
> sqlite> create table nl(value text);
> sqlite> insert into nl values ('appel'), ('beer'), ('index'), ('ijs'), ('jammer'), ('winkel'), ('zon');
> sqlite> select * from nl order by value;
> appel
> beer
> index
> ijs
> jammer
> winkel
> zon
> ```
>
> Find a dutch friend and ask them if this is the correct order for this list. Unfortunately it
> isn't. Even though it is `ij` as two separate characters, it should be sorted as if it was `ĳ` or
> `ÿ`, before `z`. Like so: appel, beer, index, jammer, winkel, ijs, zon
>
> While there is not much you can do about it now (you're just beginning!) it is something you
> should be aware of. When you later need to know about this, you will find the term 'collation'
> useful, and you'll find the procedure is different for every database engine.
{: .tip}

# Filtering

One of the most powerful features of a database is
the ability to filter data,
i.e.,
to select only those records that match certain criteria.
For example,
suppose we want to see when a particular site was visited.
We can select these records from the `Visited` table
by using a `WHERE` clause in our query:

```sql
SELECT * FROM Visited WHERE site = 'DR-1';
```

The database manager executes this query in two stages.
First,
it checks at each row in the `Visited` table
to see which ones satisfy the `WHERE`.
It then uses the column names following the `SELECT` keyword
to determine which columns to display.

This processing order means that
we can filter records using `WHERE`
based on values in columns that aren't then displayed:

```sql
SELECT id FROM Visited WHERE site = 'DR-1';
```

![SQL Filtering in Action]`(../../images/carpentries-sql/sql-filter.svg)

We can use many other Boolean operators to filter our data.
For example,
we can ask for all information from the DR-1 site collected before 1930:

```sql
SELECT * FROM Visited WHERE site = 'DR-1' AND dated < '1930-01-01';
```

> <tip-title>Date Types</tip-title>
>
> Most database managers have a special data type for dates.
> In fact, many have two:
> one for dates,
> such as "May 31, 1971",
> and one for durations,
> such as "31 days".
> SQLite doesn't:
> instead,
> it stores dates as either text
> (in the ISO-8601 standard format "YYYY-MM-DD HH:MM:SS.SSSS"),
> real numbers
> ([Julian days](https://en.wikipedia.org/wiki/Julian_day), the number of days since November 24, 4714 BCE),
> or integers
> ([Unix time](https://en.wikipedia.org/wiki/Unix_time), the number of seconds since midnight, January 1, 1970).
> If this sounds complicated,
> it is,
> but not nearly as complicated as figuring out
> [historical dates in Sweden](https://en.wikipedia.org/wiki/Swedish_calendar).
{: .tip}

> <tip-title>'30: 1930 or 2030?</tip-title>
> Storing the year as the last two digits causes problems in databases, and is part of what caused
> [Y2K](https://en.wikipedia.org/wiki/Year_2000_problem). Be sure to use the databases' built in
> format for storing dates, if it is available as that will generally avoid any major issues.
>
> Similarly there is a ["Year 2038 problem"](https://en.wikipedia.org/wiki/Year_2000_problem#Year_2038_problem),
> as the timestamps mentioned above that count seconds since Jan 1, 1970 were running out of space
> on 32-bit machines. Many systems have since migrated to work around this with 64-bit timestamps.
{: .tip}

If we want to find out what measurements were taken by either Lake or Roerich,
we can combine the tests on their names using `OR`:

```sql
SELECT * FROM Survey WHERE person = 'lake' OR person = 'roe';
```

Alternatively,
we can use `IN` to see if a value is in a specific set:

```sql
SELECT * FROM Survey WHERE person IN ('lake', 'roe');
```

We can combine `AND` with `OR`,
but we need to be careful about which operator is executed first.
If we *don't* use parentheses,
we get this:

```sql
SELECT * FROM Survey WHERE quant = 'sal' AND person = 'lake' OR person = 'roe';
```

which is salinity measurements by Lake,
and *any* measurement by Roerich.
We probably want this instead:

```sql
SELECT * FROM Survey WHERE quant = 'sal' AND (person = 'lake' OR person = 'roe');
```

We can also filter by partial matches.  For example, if we want to
know something just about the site names beginning with "DR" we can
use the `LIKE` keyword.  The percent symbol acts as a
wildcard, matching any characters in that
place.  It can be used at the beginning, middle, or end of the string
[See this page on wildcards](https://www.w3schools.com/sql/sql_wildcards.asp) for more information:

```sql
SELECT * FROM Visited WHERE site LIKE 'DR%';
```


Finally,
we can use `DISTINCT` with `WHERE`
to give a second level of filtering:

```sql
SELECT DISTINCT person, quant FROM Survey WHERE person = 'lake' OR person = 'roe';
```

But remember:
`DISTINCT` is applied to the values displayed in the chosen columns,
not to the entire rows as they are being processed.

> <tip-title>Growing Queries</tip-title>
>
> What we have just done is how most people "grow" their {SQL} queries.
> We started with something simple that did part of what we wanted,
> then added more clauses one by one,
> testing their effects as we went.
> This is a good strategy --- in fact,
> for complex queries it's often the *only* strategy --- but
> it depends on quick turnaround,
> and on us recognizing the right answer when we get it.
>
> The best way to achieve quick turnaround is often
> to put a subset of data in a temporary database
> and run our queries against that,
> or to fill a small database with synthesized records.
> For example,
> instead of trying our queries against an actual database of 20 million Australians,
> we could run it against a sample of ten thousand,
> or write a small program to generate ten thousand random (but plausible) records
> and use that.
{: .tip}

> <question-title>Fix This Query</question-title>
>
> Suppose we want to select all sites that lie within 48 degrees of the equator.
> Our first query is:
>
> ~~~
> SELECT * FROM Site WHERE (lat > -48) OR (lat < 48);
> ~~~
>
> Explain why this is wrong,
> and rewrite the query so that it is correct.
>
> > <solution-title></solution-title>
> >
> > Because we used `OR`, a site on the South Pole for example will still meet
> > the second criteria and thus be included. Instead, we want to restrict this
> > to sites that meet _both_ criteria:
> >
> > ~~~
> > SELECT * FROM Site WHERE (lat > -48) AND (lat < 48);
> > ~~~
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

> <question-title>Finding Outliers</question-title>
>
> Normalized salinity readings are supposed to be between 0.0 and 1.0.
> Write a query that selects all records from `Survey`
> with salinity values outside this range.
>
> > <solution-title></solution-title>
> >
> > ~~~
> > SELECT * FROM Survey WHERE quant = 'sal' AND ((reading > 1.0) OR (reading < 0.0));
> > ~~~
> >
> > |taken     |person    |quant     |reading   |
> > |----------|----------|----------|----------|
> > |752       |roe       |sal       |41.6      |
> > |837       |roe       |sal       |22.5      |
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

> <question-title> Matching Patterns</question-title>
>
> Which of these expressions are true?
>
> 1. `'a' LIKE 'a'`
> 2. `'a' LIKE '%a'`
> 3. `'beta' LIKE '%a'`
> 4. `'alpha' LIKE 'a%%'`
> 5. `'alpha' LIKE 'a%p%'`
>
> > <solution-title></solution-title>
> >
> > 1. True because these are the same character.
> > 2. True because the wildcard can match _zero_ or more characters.
> > 3. True because the `%` matches `bet` and the `a` matches the `a`.
> > 4. True because the first wildcard matches `lpha` and the second wildcard matches zero characters (or vice versa).
> > 5. True because the first wildcard matches `l` and the second wildcard matches `ha`.
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

> <tip-title>Case-insensitive matching</tip-title>
> But what about if you don't care about if it's `ALPHA` or `alpha` in the database, and you are
> using a language that has a notion of case (unlike e.g. Chinese, Japenese)?
>
> Then you can use the `ILIKE` operator for 'case Insensitive LIKE'.
> for example the following are true:
>
> - `'a' ILIKE 'A'`
> - `'AlPhA' ILIKE '%lpha'`
{: .tip}


# Calculating New Values

After carefully re-reading the expedition logs,
we realize that the radiation measurements they report
may need to be corrected upward by 5%.
Rather than modifying the stored data,
we can do this calculation on the fly
as part of our query:

```sql
SELECT 1.05 * reading FROM Survey WHERE quant = 'rad';
```

When we run the query,
the expression `1.05 * reading` is evaluated for each row.
Expressions can use any of the fields,
all of usual arithmetic operators,
and a variety of common functions.
(Exactly which ones depends on which database manager is being used.)
For example,
we can convert temperature readings from Fahrenheit to Celsius
and round to two decimal places:

```sql
SELECT taken, round(5 * (reading - 32) / 9, 2) FROM Survey WHERE quant = 'temp';
```


As you can see from this example, though, the string describing our
new field (generated from the equation) can become quite unwieldy. {SQL}
allows us to rename our fields, any field for that matter, whether it
was calculated or one of the existing fields in our database, for
succinctness and clarity. For example, we could write the previous
query as:

```sql
SELECT taken, round(5 * (reading - 32) / 9, 2) as Celsius FROM Survey WHERE quant = 'temp';
```

We can also combine values from different fields,
for example by using the string concatenation operator `||`:

```sql
SELECT personal || ' ' || family FROM Person;
```

But of course that can also be solved by simply having a single name field which avoids other
issues.

> <question-title> Fixing Salinity Readings</question-title>
>
> After further reading,
> we realize that Valentina Roerich
> was reporting salinity as percentages.
> Write a query that returns all of her salinity measurements
> from the `Survey` table
> with the values divided by 100.
>
> > <solution-title></solution-title>
> >
> > ~~~
> > SELECT taken, reading / 100 FROM Survey WHERE person = 'roe' AND quant = 'sal';
> > ~~~
> >
> > |taken     |reading / 100|
> > |----------|-------------|
> > |752       |0.416        |
> > |837       |0.225        |
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

> <question-title>Unions</question-title>
>
> The `UNION` operator combines the results of two queries:
>
> ~~~
> SELECT * FROM Person WHERE id = 'dyer' UNION SELECT * FROM Person WHERE id = 'roe';
> ~~~
>
> |id  |personal |family |
> |----|-------- |-------|
> |dyer|William  |Dyer   |
> |roe |Valentina|Roerich|
>
> The `UNION ALL` command is equivalent to the `UNION` operator,
> except that `UNION ALL` will select all values.
> The difference is that `UNION ALL` will not eliminate duplicate rows.
> Instead, `UNION ALL` pulls all rows from the query
> specifics and combines them into a table.
> The `UNION` command does a `SELECT DISTINCT` on the results set.
> If all the records to be returned are unique from your union,
> use `UNION ALL` instead, it gives faster results since it skips the `DISTINCT` step.
> For this section, we shall use UNION.
>
> Use `UNION` to create a consolidated list of salinity measurements
> in which Valentina Roerich's, and only Valentina's,
> have been corrected as described in the previous challenge.
> The output should be something like:
>
> |taken|reading|
> |-----|-------|
> |619  |0.13   |
> |622  |0.09   |
> |734  |0.05   |
> |751  |0.1    |
> |752  |0.09   |
> |752  |0.416  |
> |837  |0.21   |
> |837  |0.225  |
>
> > <solution-title></solution-title>
> >
> > ~~~
> > SELECT taken, reading FROM Survey WHERE person != 'roe' AND quant = 'sal' UNION SELECT taken, reading / 100 FROM Survey WHERE person = 'roe' AND quant = 'sal' ORDER BY taken ASC;
> > ~~~
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

> <question-title>Selecting Major Site Identifiers</question-title>
>
> The site identifiers in the `Visited` table have two parts
> separated by a '-':
>
> ~~~
> SELECT DISTINCT site FROM Visited;
> ~~~
>
> |site |
> |-----|
> |DR-1 |
> |DR-3 |
> |MSK-4|
>
> Some major site identifiers (i.e. the letter codes) are two letters long and some are three.
> The "in string" function `instr(X, Y)`
> returns the 1-based index of the first occurrence of string Y in string X,
> or 0 if Y does not exist in X.
> The substring function `substr(X, I, [L])`
> returns the substring of X starting at index I, with an optional length L.
> Use these two functions to produce a list of unique major site identifiers.
> (For this data,
> the list should contain only "DR" and "MSK").
>
> > <solution-title></solution-title>
> > ```
> > SELECT DISTINCT substr(site, 1, instr(site, '-') - 1) AS MajorSite FROM Visited;
> > ```
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

# Missing Data

Real-world data is never complete --- there are always holes.
Databases represent these holes using a special value called `null`.
`null` is not zero, `False`, or the empty string;
it is a one-of-a-kind value that means "nothing here".
Dealing with `null` requires a few special tricks
and some careful thinking.

By default, the Python SQL interface does not display NULL values in its output, instead it shows `None`.

To start,
let's have a look at the `Visited` table.
There are eight records,
but #752 doesn't have a date --- or rather,
its date is null:

```sql
SELECT * FROM Visited;
```

Null doesn't behave like other values.
If we select the records that come before 1930:

```sql
SELECT * FROM Visited WHERE dated < '1930-01-01';
```

we get two results,
and if we select the ones that come during or after 1930:

```sql
SELECT * FROM Visited WHERE dated >= '1930-01-01';
```

we get five,
but record #752 isn't in either set of results.
The reason is that
`null<'1930-01-01'`
is neither true nor false:
null means, "We don't know,"
and if we don't know the value on the left side of a comparison,
we don't know whether the comparison is true or false.
Since databases represent "don't know" as null,
the value of `null<'1930-01-01'`
is actually `null`.
`null>='1930-01-01'` is also null
because we can't answer to that question either.
And since the only records kept by a `WHERE`
are those for which the test is true,
record #752 isn't included in either set of results.

Comparisons aren't the only operations that behave this way with nulls.
`1+null` is `null`,
`5*null` is `null`,
`log(null)` is `null`,
and so on.
In particular,
comparing things to null with = and != produces null:

```sql
SELECT * FROM Visited WHERE dated = NULL;
```

produces no output, and neither does:

```sql
SELECT * FROM Visited WHERE dated != NULL;
```

To check whether a value is `null` or not,
we must use a special test `IS NULL`:

```sql
SELECT * FROM Visited WHERE dated IS NULL;
```

or its inverse `IS NOT NULL`:

```sql
SELECT * FROM Visited WHERE dated IS NOT NULL;
```

Null values can cause headaches wherever they appear.
For example,
suppose we want to find all the salinity measurements
that weren't taken by Lake.
It's natural to write the query like this:

```sql
SELECT * FROM Survey WHERE quant = 'sal' AND person != 'lake';
```

but this query filters omits the records
where we don't know who took the measurement.
Once again,
the reason is that when `person` is `null`,
the `!=` comparison produces `null`,
so the record isn't kept in our results.
If we want to keep these records
we need to add an explicit check:

```sql
SELECT * FROM Survey WHERE quant = 'sal' AND (person != 'lake' OR person IS NULL);
```


We still have to decide whether this is the right thing to do or not.
If we want to be absolutely sure that
we aren't including any measurements by Lake in our results,
we need to exclude all the records for which we don't know who did the work.

In contrast to arithmetic or Boolean operators, aggregation functions
that combine multiple values, such as `min`, `max` or `avg`, *ignore*
`null` values. In the majority of cases, this is a desirable output:
for example, unknown values are thus not affecting our data when we
are averaging it. Aggregation functions will be addressed in more
detail in [the next section](#).

> <question-title>Sorting by Known Date</question-title>
>
> Write a query that sorts the records in `Visited` by date,
> omitting entries for which the date is not known
> (i.e., is null).
>
> > <solution-title></solution-title>
> >
> > ~~~
> > SELECT * FROM Visited WHERE dated IS NOT NULL ORDER BY dated ASC;
> > ~~~
> >
> > |id        |site      |dated     |
> > |----------|----------|----------|
> > |619       |DR-1      |1927-02-08|
> > |622       |DR-1      |1927-02-10|
> > |734       |DR-3      |1930-01-07|
> > |735       |DR-3      |1930-01-12|
> > |751       |DR-3      |1930-02-26|
> > |837       |MSK-4     |1932-01-14|
> > |844       |DR-1      |1932-03-22|
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

> <question-title>NULL in a Set</question-title>
>
> What do you expect the following query to produce?
>
> ~~~
> SELECT * FROM Visited WHERE dated IN ('1927-02-08', NULL);
> ~~~
>
> What does it actually produce?
> > <solution-title></solution-title>
> >
> > You might expect the above query to return rows where dated is either '1927-02-08' or NULL.
> > Instead it only returns rows where dated is '1927-02-08', the same as you would get from this
> > simpler query:
> >
> > ~~~
> > SELECT * FROM Visited WHERE dated IN ('1927-02-08');
> > ~~~
> >
> > The reason is that the `IN` operator works with a set of *values*, but NULL is by definition
> > not a value and is therefore simply ignored.
> >
> > If we wanted to actually include NULL, we would have to rewrite the query to use the IS NULL condition:
> >
> > ~~~
> > SELECT * FROM Visited WHERE dated = '1927-02-08' OR dated IS NULL;
> > ~~~
> {: .solution}
{: .question}

```sql
-- Try solutions here!
```

> <question-title>Pros and Cons of Sentinels</question-title>
>
> Some database designers prefer to use
> a sentinel value
> to mark missing data rather than `null`.
> For example,
> they will use the date "0000-00-00" to mark a missing date,
> or -1.0 to mark a missing salinity or radiation reading
> (since actual readings cannot be negative).
> What does this simplify?
> What burdens or risks does it introduce?
{: .question}

```sql
-- Try solutions here!
```
