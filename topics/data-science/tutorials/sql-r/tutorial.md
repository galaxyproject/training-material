---
layout: tutorial_hands_on
title: SQL with R
level: Intermediate
zenodo_link:
requirements:
- type: "internal"
  topic_name: data-science
  tutorials:
      - r-basics
      - r-advanced
- type: "internal"
  topic_name: data-science
  tutorials:
      - sql-advanced
follow_up_training:  []

questions:
- "How can I access databases from programs written in R?"
objectives:
- "Write short programs that execute SQL queries."
- "Trace the execution of a program that contains an SQL query."
- "Explain why most database applications are written in a general-purpose language rather than in SQL."
time_estimation:  45M
key_points:
- "Data analysis languages have libraries for accessing databases."
- "To connect to a database, a program must use a library specific to that database manager."
- "R's libraries can be used to directly query or read from a database."
- "Programs can read query results in batches or all at once."
- "Queries should be written using parameter substitution, not string formatting."
- "R has multiple helper functions to make working with databases easier."

contributors:
- carpentries
- hexylena
- avans-atgm

subtopic: sql

notebook:
    language: r

tags:
- SQL
- R
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

For this tutorial we need to download a database that we will use for the queries.

```r
download.file("http://swcarpentry.github.io/sql-novice-survey/files/survey.db", destfile="survey.db")
```

# Programming with Databases - R

Let's have a look at how to access a database from
a data analysis language like R.
Other languages use almost exactly the same model:
library and function names may differ,
but the concepts are the same.

Here's a short R program that selects latitudes and longitudes
from an SQLite database stored in a file called `survey.db`:

```r
library(RSQLite)
connection <- dbConnect(SQLite(), "survey.db")
results <- dbGetQuery(connection, "SELECT Site.lat, Site.long FROM Site;")
print(results)
dbDisconnect(connection)
```

The program starts by importing the `RSQLite` library.
If we were connecting to MySQL, DB2, or some other database,
we would import a different library,
but all of them provide the same functions,
so that the rest of our program does not have to change
(at least, not much)
if we switch from one database to another.

Line 2 establishes a connection to the database.
Since we're using SQLite,
all we need to specify is the name of the database file.
Other systems may require us to provide a username and password as well.

On line 3, we retrieve the results from an SQL query.
It's our job to make sure that SQL is properly formatted;
if it isn't,
or if something goes wrong when it is being executed,
the database will report an error.
This result is a dataframe with one row for each entry and one column for each column in the database.

Finally, the last line closes our connection,
since the database can only keep a limited number of these open at one time.
Since establishing a connection takes time,
though,
we shouldn't open a connection,
do one operation,
then close the connection,
only to reopen it a few microseconds later to do another operation.
Instead,
it's normal to create one connection that stays open for the lifetime of the program.

Queries in real applications will often depend on values provided by users.
For example,
this function takes a user's ID as a parameter and returns their name:

```r
library(RSQLite)

connection <- dbConnect(SQLite(), "survey.db")

getName <- function(personID) {
  query <- paste0("SELECT personal || ' ' || family FROM Person WHERE id =='",
                  personID, "';")
  return(dbGetQuery(connection, query))
}

print(paste("full name for dyer:", getName('dyer')))

dbDisconnect(connection)
```

We use string concatenation on the first line of this function
to construct a query containing the user ID we have been given.
This seems simple enough,
but what happens if someone gives us this string as input?

```
dyer'; DROP TABLE Survey; SELECT '
```

It looks like there's garbage after the user's ID,
but it is very carefully chosen garbage.
If we insert this string into our query,
the result is:

```
SELECT personal || ' ' || family FROM Person WHERE id='dyer'; DROP TABLE Survey; SELECT '';
```

If we execute this,
it will erase one of the tables in our database.

This is called an SQL injection attack,
and it has been used to attack thousands of programs over the years.
In particular,
many web sites that take data from users insert values directly into queries
without checking them carefully first.
A very [relevant XKCD](https://xkcd.com/327/) that explains the
dangers of using raw input in queries a little more succinctly:

![A 4 panel comic, in the first panel a person is shown answering the phone, hearing that their son's school has some computer trouble. In panel 2 they apologises asking if their child broke something. In panel 3, the unseen person on the other end of the phone call asks if they really named their son Robert'); Drop table students;--? They respond saying 'oh yes. little bobby tables we call him.' In the 4th panel the caller says 'well we have lost this years student records, I hope you're happy.' They respond 'And I hope you've learned to sanitize your database inputs'.](../../images/xkcd/exploits_of_a_mom.png)

Since an unscrupulous parent might try to smuggle commands into our queries in many different ways,
the safest way to deal with this threat is
to replace characters like quotes with their escaped equivalents,
so that we can safely put whatever the user gives us inside a string.
We can do this by using a prepared statement
instead of formatting our statements as strings.
Here's what our example program looks like if we do this:

```r
library(RSQLite)
connection <- dbConnect(SQLite(), "survey.db")

getName <- function(personID) {
  query <- "SELECT personal || ' ' || family FROM Person WHERE id == ?"
  return(dbGetPreparedQuery(connection, query, data.frame(personID)))
}

print(paste("full name for dyer:", getName('dyer')))

dbDisconnect(connection)
```

The key changes are in the query string and the `dbGetQuery` call (we use dbGetPreparedQuery instead).
Instead of formatting the query ourselves,
we put question marks in the query template where we want to insert values.
When we call `dbGetPreparedQuery`,
we provide a dataframe
that contains as many values as there are question marks in the query.
The library matches values to question marks in order,
and translates any special characters in the values
into their escaped equivalents
so that they are safe to use.

> <question-title>Filling a Table vs. Printing Values</question-title>
>
> Write an R program that creates a new database in a file called
> `original.db` containing a single table called `Pressure`, with a
> single field called `reading`, and inserts 100,000 random numbers
> between 10.0 and 25.0.  How long does it take this program to run?
> How long does it take to run a program that simply writes those
> random numbers to a file?
{: .question}

> <question-title>Filtering in SQL vs. Filtering in R</question-title>
>
> Write an R program that creates a new database called
> `backup.db` with the same structure as `original.db` and copies all
> the values greater than 20.0 from `original.db` to `backup.db`.
> Which is faster: filtering values in the query, or reading
> everything into memory and filtering in R?
{: .question}

## Database helper functions in R

R's database interface packages (like `RSQLite`) all share
a common set of helper functions useful for exploring databases and
reading/writing entire tables at once.

To view all tables in a database, we can use `dbListTables()`:

```r
connection <- dbConnect(SQLite(), "survey.db")
dbListTables(connection)
```


To view all column names of a table, use `dbListFields()`:

```r
dbListFields(connection, "Survey")
```


To read an entire table as a dataframe, use `dbReadTable()`:

```r
dbReadTable(connection, "Person")
```


Finally to write an entire table to a database, you can use `dbWriteTable()`.
Note that we will always want to use the `row.names = FALSE` argument or R
will write the row names as a separate column.
In this example we will write R's built-in `iris` dataset as a table in `survey.db`.

```r
dbWriteTable(connection, "iris", iris, row.names = FALSE)
head(dbReadTable(connection, "iris"))
```

And as always, remember to close the database connection when done!

```r
dbDisconnect(connection)
```
