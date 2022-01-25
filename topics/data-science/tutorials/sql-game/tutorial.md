---
layout: tutorial_hands_on

title: "SQL Educational Game - Murder Mystery"
level: Intermediate
zenodo_link: ""
requirements:
- type: "internal"
  topic_name: data-science
  tutorials:
      - sql-basic
      - sql-advanced
follow_up_training: []

questions:
- "Who did the crime?"
objectives:
- Explore SQL City and discover who committed the murder
- Reinforce your experiences with SQL such as querying, filtering, and joining data.
time_estimation: 2H
key_points:
- Learning SQL can be fun!
notebook:
  language: python
subtopic: sql
contributors:
- hexylena
- NUKnightLab
- erasmusplus
- avans-atgm
tags:
- game
- SQL
- jupyter-notebook
---

This is not a tutorial like most GTN content but a fun exercise for you to play around and learn a bit about SQL in a more 'practical', and hopefully re-inforce the skills you covered in Basic and Advanced SQL skills. It makes use of the [NUKnightLab/sql-mysteries](https://github.com/NUKnightLab/sql-mysteries) SQL murder mystery project and released under open licenses:

> Original code for NUKnightLab/sql-mysteries is released under the MIT License.
> Original text and other content is released under Creative Commons CC BY-SA 4.0.
{: .quote}

Download the database and connector:

```python
# This preamble sets up the sql "magic" for jupyter. Use %%sql in your cells to write sql!
!python3 -m pip install ipython-sql sqlalchemy
!wget -c https://github.com/NUKnightLab/sql-mysteries/raw/master/sql-murder-mystery.db
```

Setup the database connection:

```python
import sqlalchemy
engine = sqlalchemy.create_engine("sqlite:///sql-murder-mystery.db")
%load_ext sql
%sql sqlite:///sql-murder-mystery.db
%config SqlMagic.displaycon=False
```

## Tables

Which tables are available to you? What columns do they contain? Here's a handy reference for you:

```python
import pandas as pd
from sqlalchemy import MetaData
m = MetaData()
m.reflect(engine)
results = []
for table in m.tables.values():
    results.append([table.name, ', '.join([c.name for c in table.c])])
pd.set_option('display.max_colwidth', None)
pd.DataFrame(results, columns=["Table", "Columns"])
```


## The search!

Begin your search for the truth

> A crime has taken place and the police are both useless *and* corrupt and it
> is up to you and your community to solve the mystery. They failed to secure
> their database, and now their crime scene reports are public, it is time to
> figure out who the murderer was.
>
> You know that the crime was a **murder** that occurred sometime on **Jan.15,
> 2018** and that it took place in **SQL City**.
>
> All the clues to this mystery are buried in a huge database, and you need to
> use SQL to navigate through this vast network of information. Your first step
> to solving the mystery is to retrieve the corresponding crime scene report
> from the police department’s database. From there, you can use your SQL
> skills to find the murderer.
{: .quote}

```python
%%sql
select * from crime_scene_report limit 8;
```

Try using the 'Insert → Cell Below' functionality to keep track of important query results as you go!

## Solution

Write the following queries in your SQL environment to check whether you've found the right murderer:

```python
%%sql
INSERT INTO solution VALUES (1, "Insert the name of the person you found here");
SELECT value FROM solution;
```
