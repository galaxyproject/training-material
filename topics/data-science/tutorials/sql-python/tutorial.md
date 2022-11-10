---
layout: tutorial_hands_on
title: SQL with Python
level: Intermediate
zenodo_link:
requirements:
- type: "internal"
  topic_name: data-science
  tutorials:
      - sql-advanced
      #- python-basic
follow_up_training:  []

questions:
- "How can I access databases from programs written in Python?"
objectives:
- "Write short programs that execute SQL queries."
- "Trace the execution of a program that contains an SQL query."
- "Explain why most database applications are written in a general-purpose language rather than in SQL."
time_estimation:  45M
key_points:
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
    language: python

tags:
- SQL
- Python
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

```python
!wget -c http://swcarpentry.github.io/sql-novice-survey/files/survey.db
```


# Programming with Databases - Python

Let's have a look at how to access a database from
a general-purpose programming language like Python.
Other languages use almost exactly the same model:
library and function names may differ,
but the concepts are the same.

Here's a short Python program that selects latitudes and longitudes
from an SQLite database stored in a file called `survey.db`:

```python
import sqlite3

connection = sqlite3.connect("survey.db")
cursor = connection.cursor()
cursor.execute("SELECT Site.lat, Site.long FROM Site;")
results = cursor.fetchall()
for r in results:
    print(r)
cursor.close()
connection.close()
```

The program starts by importing the `sqlite3` library.
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
Line 3 then uses this connection to create a cursor.
Just like the cursor in an editor,
its role is to keep track of where we are in the database.

On line 4, we use that cursor to ask the database to execute a query for us.
The query is written in SQL,
and passed to `cursor.execute` as a string.
It's our job to make sure that SQL is properly formatted;
if it isn't,
or if something goes wrong when it is being executed,
the database will report an error.

The database returns the results of the query to us
in response to the `cursor.fetchall` call on line 5.
This result is a list with one entry for each record in the result set;
if we loop over that list (line 6) and print those list entries (line 7),
we can see that each one is a tuple
with one element for each field we asked for.

Finally, lines 8 and 9 close our cursor and our connection,
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

```python
import sqlite3

def get_name(database_file, person_id):
    query = "SELECT personal || ' ' || family FROM Person WHERE id='" + person_id + "';"

    connection = sqlite3.connect(database_file)
    cursor = connection.cursor()
    cursor.execute(query)
    results = cursor.fetchall()
    cursor.close()
    connection.close()

    return results[0][0]

print("Full name for dyer:", get_name('survey.db', 'dyer'))
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

Since a villain might try to smuggle commands into our queries in many different ways,
the safest way to deal with this threat is
to replace characters like quotes with their escaped equivalents,
so that we can safely put whatever the user gives us inside a string.
We can do this by using a prepared statement
instead of formatting our statements as strings.
Here's what our example program looks like if we do this:

```python
import sqlite3

def get_name(database_file, person_id):
    query = "SELECT personal || ' ' || family FROM Person WHERE id=?;"

    connection = sqlite3.connect(database_file)
    cursor = connection.cursor()
    cursor.execute(query, [person_id])
    results = cursor.fetchall()
    cursor.close()
    connection.close()

    return results[0][0]

print("Full name for dyer:", get_name('survey.db', 'dyer'))
```

The key changes are in the query string and the `execute` call.
Instead of formatting the query ourselves,
we put question marks in the query template where we want to insert values.
When we call `execute`,
we provide a list
that contains as many values as there are question marks in the query.
The library matches values to question marks in order,
and translates any special characters in the values
into their escaped equivalents
so that they are safe to use.

We can also use `sqlite3`'s cursor to make changes to our database,
such as inserting a new name.
For instance, we can define a new function called `add_name` like so:

```python
import sqlite3

def add_name(database_file, new_person):
    query = "INSERT INTO Person (id, personal, family) VALUES (?, ?, ?);"

    connection = sqlite3.connect(database_file)
    cursor = connection.cursor()
    cursor.execute(query, list(new_person))
    cursor.close()
    connection.close()


def get_name(database_file, person_id):
    query = "SELECT personal || ' ' || family FROM Person WHERE id=?;"

    connection = sqlite3.connect(database_file)
    cursor = connection.cursor()
    cursor.execute(query, [person_id])
    results = cursor.fetchall()
    cursor.close()
    connection.close()

    return results[0][0]

# Insert a new name
add_name('survey.db', ('barrett', 'Mary', 'Barrett'))
# Check it exists
print("Full name for barrett:", get_name('survey.db', 'barrett'))
```


Note that in versions of sqlite3 >= 2.5, the `get_name` function described
above will fail with an `IndexError: list index out of range`,
even though we added Mary's
entry into the table using `add_name`.
This is because we must perform a `connection.commit()` before closing
the connection, in order to save our changes to the database.

```python
import sqlite3

def add_name(database_file, new_person):
    query = "INSERT INTO Person (id, personal, family) VALUES (?, ?, ?);"

    connection = sqlite3.connect(database_file)
    cursor = connection.cursor()
    cursor.execute(query, list(new_person))
    cursor.close()
    connection.commit()
    connection.close()


def get_name(database_file, person_id):
    query = "SELECT personal || ' ' || family FROM Person WHERE id=?;"

    connection = sqlite3.connect(database_file)
    cursor = connection.cursor()
    cursor.execute(query, [person_id])
    results = cursor.fetchall()
    cursor.close()
    connection.close()

    return results[0][0]

# Insert a new name
add_name('survey.db', ('barrett', 'Mary', 'Barrett'))
# Check it exists
print("Full name for barrett:", get_name('survey.db', 'barrett'))
```


> <question-title>Filling a Table vs. Printing Values</question-title>
>
> Write a Python program that creates a new database in a file called
> `original.db` containing a single table called `Pressure`, with a
> single field called `reading`, and inserts 100,000 random numbers
> between 10.0 and 25.0.  How long does it take this program to run?
> How long does it take to run a program that simply writes those
> random numbers to a file?
>
> > <solution-title></solution-title>
> > ```
> > import sqlite3
> > # import random number generator
> > from numpy.random import uniform
> >
> > random_numbers = uniform(low=10.0, high=25.0, size=100000)
> >
> > connection = sqlite3.connect("original.db")
> > cursor = connection.cursor()
> > cursor.execute("CREATE TABLE Pressure (reading float not null)")
> > query = "INSERT INTO Pressure (reading) VALUES (?);"
> >
> > for number in random_numbers:
> >     cursor.execute(query, [number])
> >
> > cursor.close()
> > connection.commit() # save changes to file for next exercise
> > connection.close()
> > ```
> > {: .python}
> >
> > For comparison, the following program writes the random numbers
> > into the file `random_numbers.txt`:
> >
> > ```
> > from numpy.random import uniform
> >
> > random_numbers = uniform(low=10.0, high=25.0, size=100000)
> > with open('random_numbers.txt', 'w') as outfile:
> >     for number in random_numbers:
> >         # need to add linebreak \n
> >         outfile.write("{}\n".format(number))
> > ```
> > {: .python}
> {: .solution}
{: .question}

> <question-title>Filtering in SQL vs. Filtering in Python</question-title>
>
> Write a Python program that creates a new database called
> `backup.db` with the same structure as `original.db` and copies all
> the values greater than 20.0 from `original.db` to `backup.db`.
> Which is faster: filtering values in the query, or reading
> everything into memory and filtering in Python?
>
> > <solution-title></solution-title>
> > The first example reads all the data into memory and filters the
> > numbers using the if statement in Python.
> >
> > ```
> > import sqlite3
> >
> > connection_original = sqlite3.connect("original.db")
> > cursor_original = connection_original.cursor()
> > cursor_original.execute("SELECT * FROM Pressure;")
> > results = cursor_original.fetchall()
> > cursor_original.close()
> > connection_original.close()
> >
> > connection_backup = sqlite3.connect("backup.db")
> > cursor_backup = connection_backup.cursor()
> > cursor_backup.execute("CREATE TABLE Pressure (reading float not null)")
> > query = "INSERT INTO Pressure (reading) VALUES (?);"
> >
> > for entry in results:
> >     # number is saved in first column of the table
> >     if entry[0] > 20.0:
> >         cursor_backup.execute(query, entry)
> >
> > cursor_backup.close()
> > connection_backup.commit()
> > connection_backup.close()
> > ```
> >
> > In contrast the following example uses the conditional ``SELECT`` statement
> > to filter the numbers in SQL.
> > The only lines that changed are in line 5, where the values are fetched
> > from `original.db` and the for loop starting in line 15 used to insert
> > the numbers into `backup.db`.
> > Note how this version does not require the use of Python's if statement.
> >
> > ```
> > import sqlite3
> >
> > connection_original = sqlite3.connect("original.db")
> > cursor_original = connection_original.cursor()
> > cursor_original.execute("SELECT * FROM Pressure WHERE reading > 20.0;")
> > results = cursor_original.fetchall()
> > cursor_original.close()
> > connection_original.close()
> >
> > connection_backup = sqlite3.connect("backup.db")
> > cursor_backup = connection_backup.cursor()
> > cursor_backup.execute("CREATE TABLE Pressure (reading float not null)")
> > query = "INSERT INTO Pressure (reading) VALUES (?);"
> >
> > for entry in results:
> >     cursor_backup.execute(query, entry)
> >
> > cursor_backup.close()
> > connection_backup.commit()
> > connection_backup.close()
> > ```
> >
> {: .solution}
{: .question}

> <question-title>Generating Insert Statements</question-title>
>
> One of our colleagues has sent us a
> CSV
> file containing
> temperature readings by Robert Olmstead, which is formatted like
> this:
>
> ```
> Taken,Temp
> 619,-21.5
> 622,-15.5
> ```
> {: .output}
>
> Write a small Python program that reads this file in and prints out
> the SQL `INSERT` statements needed to add these records to the
> survey database.  Note: you will need to add an entry for Olmstead
> to the `Person` table.  If you are testing your program repeatedly,
> you may want to investigate SQL's `INSERT or REPLACE` command.
{: .question}
