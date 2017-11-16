---
layout: tutorial_hands_on
topic_name: admin
tutorial_name: database-schema
---

Galaxy Database Schema
======================


# Requirements

For the hands-on examples you need access to a Galaxy server and access to its PostgreSQL database. You can set-up this yourself, or use the Galaxy Docker Image provided by Björn Grüning (https://github.com/bgruening/docker-galaxy-stable). During this tutorial, we will work with the Galaxy Docker Image.

Setting up Docker and using the Galaxy Docker Image:
(please do this before the tutorial, preferably when you are still at home using a fast internet connection)

Follow the instruction to install the Docker engine: https://docs.docker.com/engine/installation/   
Execute:  docker run -d -p 8080:80 bgruening/galaxy-stable
(this will download the Galaxy Docker Image, when executed the first time, and
start it)
Test the Docker image in your web browser: http://localhost:8080
(see first paragraph of Björn’s introduction, for special cases when using the
Docker Toolbox on Mac and Windows). Quit by: docker kill NAME (you get the name with:  docker ps )

# Introduction

## Database versus Object Model

The session description is database centric and we’ll be focusing on the relational
database that backs Galaxy servers.  But that’s only half the picture of the this data.
The other is the object model which is the object-oriented view of this same data.
The object model is used by the code to manipulate and access the database.
The translation between the two worlds is handled by an object-relational mapping implemented with SQLAlchemy (https://www.sqlalchemy.org).

Today we are covering the database and how to access it with SQL.  We aren’t going to cover the corresponding object model or object relational mapping.

## Database Platform

The default out-of-the-box Galaxy installation uses SQLite (https://www.sqlite.org).
SQLite is a lightweight database management system (DBMS) that can be packaged inside Galaxy and does not require any additional steps at initial setup time.

However, SQLite is not the recommended DBMS for running a Galaxy server. The recommended production DMBS for
Galaxy is PostgreSQL (https://www.postgresql.org). PostgreSQL offers a full set of DBMS features and robust support for multiple simultaneous users.

This workshop will be entirely based in PostgreSQL (also referred to as Postgres).

## What is in (and not in) the Galaxy database?

The Galaxy database contains management information about your server.  The database tracks users, groups, jobs, histories, datasets, workflows and so on.  

What’s not in the database is the data. Datasets are stored outside the database. The database does keep metadata – information about the datasets such as data type. The tools themselves are not stored in the database either.

## Understanding the Database Schema

#### ER diagrams and SchemaSpy

Entity-relationship diagrams are a way to understand tables and the relationships between them inside a relational database.  SchemaSpy (http://schemaspy.sourceforge.net/) is a free (and remarkable tool) for generating ER diagrams.  We’ve used it generate a description of the database backing the server in this container.  See

    <!-- TODO: following link is broken, update with new link once fixed on Hub -->
    https://galaxyproject.org/schema/SchemaSpy/index.html

The “Tables” tab is a good place to start learning the structure of the database.  Each table represents a different type of thing, and often that thing is itself a relationship. For example, each record in the dataset table has information about a specific dataset, while records in the history_dataset_association table have information about what histories that dataset is in.

Each SchemaSpy table’s page shows the attributes in that table, as well as any constraints on those attributes, and the relationships between that table and other tables.

Also see the “Run SchemaSpy in this container” section below for how to install and then run SchemaSpy yourself.

#### Database conventions
The Galaxy database uses a number of naming and design conventions.  Understanding these can make navigating the database much easier.

#### id attributes
Every table has an id column that uniquely identifies each row.  (The id column is the primary key in database terminology.) Beyond uniquely identifying a row in the table, ID values have no meaning.  ID values are unique within a table, but not across the database.

#### Relationships between tables, and `_id` columns
Relationships between tables are implemented by exporting id columns from one table into another.  Imported ids are called foreign keys in database nomenclature, and are uniformly named
   table_the_id_came_from_id

There are a few notable exceptions to this rule.  If the ID is from a table that is prefixed with galaxy_, for example, galaxy_user or galaxy_session, the  galaxy_ will be dropped from the column name.  For example, galaxy_user.id becomes user_id in the over 50 tables it is imported into

#### Relationship tables

As mentioned previously, some tables, such as history_dataset_association represent relationships between things, rather than things themselves.  In this case history_dataset_association describes relationships between datasets and histories.

Relationship table names typically contain the names of tables they are relating, suffixed with `_association`.

Why are nulls allowed in almost every column?
We have no idea.  In practice, they aren’t nulls in most of those columns.

Why aren’t there comments, on anything?
PostgreSQL supports comments to table definitions, but there are none shown in the SchemaSpy report. Why? The table definitions are actually generated by SQLAlchemy, the object-relational mapping software used by Galaxy, and SQLAlchemy does not support it.

There is nothing in the database that results from direct manipulation of the table definitions through DDL.  Everything comes in through SQLAlchemy.

## Start Docker and Galaxy

> ### {% icon hands_on %} ***Hands on!***
>
>   1. Start the Galaxy Docker Image -  this time as an interactive session
>
>    ```sh
>       docker run -i -t -p 8080:80 bgruening/galaxy-stable /bin/bash
>    ```
>
>   2. Start Galaxy and its PostgreSQL server
>
>    ```sh
>       startup > log 2>&1 &
>    ```
>
>   3. Follow the startup process
>
>    ```sh
>       tail -f log
>    ```


## Important tables

> ### {% icon hands_on %} ***Hands on!***
>
>   1. Connect to the PostgreSQL database (change to user galaxy first)
>
>    ```sh 
>        su galaxy
>        psql -d galaxy -U galaxy
>    ```
>
>   2. List all tables
>
>    ```sql
>       \dt
>    ```


Enter `q` to exit the view results page, and space to see the next results page.

### Table “galaxy_user”

> ### {% icon hands_on %} Hands-on
>
>    ```sql
>        select * from galaxy_user;
>    ```

As described in Björn’s introduction, an Admin user is already pre-set (email: ‘admin@galaxy.org’, password: ‘admin’). Now let’s add (i.e. register) a new user via the Galaxy website. And check the database:

> ### {% icon hands_on %} Hands-on
>
>    ```sql
>        select * from galaxy_user;
>    ```

### Table “job”

> ### {% icon hands_on %} Hands-on
>
>    ```sql
>        select * from job;
>    ```

Run a few jobs on the galaxy website (e.g _upload file_ a simple table and _add column_ with ‘Iterate’ no and yes) and check the database again:

> ### {% icon hands_on %} Hands-on
>
>    ```sql
>        select * from job \x\g\x
>    ```

### Table “job_parameter”

> ### {% icon hands_on %} Hands-on
>
>   ```sql
>       select * from job_parameter;
>   ```


### Table “history”

> ### {% icon hands_on %} Hands-on
>
>   ```sql
>       select * from history;
>   ```

Give your current history a name and check the database again.


### Table “dataset”

> ### {% icon hands_on %} Hands-on
>
>   ```sql
>       select * from dataset;
>   ```

### Table “history_dataset_association”

> ### {% icon hands_on %} Hands-on
>
>   ```sql
>       select * from history_dataset_association;
>   ```


## More (hands-on) Examples, not covered by the reports app

Have a look at the reports up (which is also provided in the Docker Image):

http://admin:admin@localhost:8080/reports/

Depending on your local needs, some queries are missing, like:

## Jobs per tool per year  /  jobs per tool since 2015

You can add the numbers per month from the reports, or:


>   ```sql
>       select j.id, j.create_time from job j limit 5;
>   ```
>
>   ```sql
>        select j.id, j.create_time from job j
>            where j.create_time >= '2015-12-31'
>            and j.create_time < '2016-12-31';
>    ```
>
>   ```sql
>        select j.id,j.create_time from job j
>           where EXTRACT(year FROM j.create_time) = 2016
>           and j.tool_id='upload1';`
>   ```
>
>...and now include the user
>
>   ```sql
>        select count(j.id) from job j, galaxy_user u
>           where j.user_id = u.id
>           and u.email = 'hansrudolf.hotz@fmi.ch'
>           and EXTRACT(year FROM j.create_time) = 2016
>           and j.tool_id='upload1';
>   ```
>
>   ```sql
>        select u.email, count(*) from job j, galaxy_user u
>           where j.user_id = u.id
>           and EXTRACT(year FROM j.create_time) = 2016
>           and j.tool_id='upload1'
>           GROUP BY u.email;
>   ```

### Jobs per tool of a certain version

Imagine the current version of a tool is working fine, however a previous version had a bug: now you wanna warn
all the users who have used the broken version, without alerting users who never used the broken one.

The following example is from the development server at the FMI

>   ```sql
>       select distinct(j.tool_version) from job j
>           where j.tool_id = 'qAlign';
>   ```
>
>   ```sql
>        select j.user_id from job j
>           where j.tool_id = 'qAlign'
>           and j.tool_version = '1.0.4quasr';
>   ```
>
>   ```sql
>       select u.email, j.create_time from job j, galaxy_user u
>           where j.user_id = u.id
>           and j.tool_id = 'qAlign'
>           and j.tool_version = '1.0.4quasr';
>   ```


## All users running a job using a certain parameter

> ### {% icon hands_on %} ***Hands on!***
>
>   ```sql
>       select jp.name, jp.value  from job_parameter jp
>           where name = 'iterate'`
>   ```
>
>   ```sql
>       select u.email, jp.name, jp.value
>           from job_parameter jp, job j, galaxy_user u
>           where jp.name = 'iterate'
>           and j.tool_id = 'addValue'
>           and jp.job_id = j.id
>           and j.user_id = u.id;
>   ```
>

## Close PostgreSQL client and quit docker

>   Close the PostgreSQL client
>
>   ```sql
>       \q
>   ```
>
>   Quit the interactive docker (change back to root first)
>
>   ```sh
>       exit
>       exit
>   ```

## Other Topics

#### How to move from MySQL to PostgreSQL

Slides: https://docs.google.com/presentation/d/1N3BDNQT3s7eQEO3BO89TQbTYwKp92fxHDRWQwC-T1kA

https://wiki.galaxyproject.org/Community/Log/2015/MySQL2PostgreSQL

#### Is there ever a need to manually change the contents of a table

Slides:
https://docs.google.com/presentation/d/1l4DD0IaJjuvk1zAT1Sjv26bLyrSOg3VUm7rD-TQl_Zs

#### Run SchemaSpy in this container

To run SchemaSpy in your container you’ll need to get it, and also install some required software packages.

>   ```sh
>   wget http://downloads.sourceforge.net/project/schemaspy/schemaspy/SchemaSpy%205.0.0/schemaSpy_5.0.0.jar
>   apt-get update
>   apt-get install libpostgresql-jdbc-java
>   apt-get install graphviz
>   ```
>
>    To run SchemaSpy:
>   ```
>   java -jar schemaSpy_5.0.0.jar -t pgsql -db galaxy -u galaxy -host localhost -s public -dp /usr/share/java/postgresql-jdbc4-9.2.jar -o SpyOut
>   ```
>

The SpyOut directory will contain the generated reports and diagrams, anchored at index.html.



# Conclusion

There is a lot of information stored in the Galaxy database. Use this information for trouble shooting when necessary and use it as a source for extendend user statistics.


# :clap: Thank you
