---
layout: tutorial_hands_on
topic_name: admin
tutorial_name: dev-to-production
---

Move from dev instance to production instance
=============================================

:grey_question: ***Questions***

- *Why do I need to configure Galaxy for production ?*
- *How to run a galaxy in a production environment ?*

:dart: ***Objectives***

- *Learn to install a Galaxy server.*
- *Upgrade a basic galaxy installation to a production environment*.
- *Get a basic understanding of entry points in the main galaxy configuration file.*

:heavy_check_mark: ***Requirements***

- *[Galaxy Server Administration]({{site.baseurl}}/topics/admin)*

:hourglass: ***Time estimation*** *TODO*

# Introduction

In this tutorial you will learn to install and configure a galaxy instance.
The basic installation instructions are suitable for developpement only.
For setting up a Galaxy for a multi-user production environment additonal steps are necessary.
After a basic installation, this tutorial present the main steps for moving from a basic installation to a production environment.

<img src="../../images/scheme-dev_to_production.png" alt="Flowchart of the processes of configuring galaxy for a production environment. After the basic installation the developer settings should be disabled. The next steps are to switch to a database server and use proxy server. A compute cluster should be used. The following would be to clean up database and rotate log files. Check the local data and enable upload via FTP and the galaxy instance is ready for the production." style="width: 600px;"/>


# Basic Installation

## Clone it and run it  

{% icon hands_on %} ***Hands on!***

Galaxy's source code is hosted in a [GitHub](https://github.com/galaxyproject/galaxy) repository.

To get the code run:

	git clone -b release_16.07 https://github.com/galaxyproject/galaxy.git

To start the galaxy go to the galaxy directory and run:

	sh run.sh

This will start up the server on localhost and port 8080, so Galaxy can be accessed from your web browser at http://localhost:8080 .

Galaxy's server will start printing its output to your terminal. To stop the Galaxy server, just hit Ctrl-c in the terminal from which Galaxy is running.

## Access Galaxy over the network
In the basic installation Galaxy is bind to the loopback interface.
To bind Galaxy to any avalaible network interface edit the config/galaxy.ini file and change the host setting to:

	host = 0.0.0.0

## What did you just installed ?
The galaxy you have just installed is configured with the following:

- [SQLite](https://www.sqlite.org/): a servless database.
- A built-in HTTP server, written in Python.

The tools are run locally and the galaxy server itself run in a single process.

## Usefull things to know ...

- The command: git clone -b release_16.07 https://github.com/galaxyproject/galaxy.git

*We use the "-b" option to clone a specific branch of the repository.
In the example above we clone the branch release_16.07 because the lastest release is 16.07.
A new version of Galaxy is released every 3 months.
It is strongly encouraged to run this tutorial with the latest release.
See the [GitHub](https://github.com/galaxyproject/galaxy) repository to get the latest one.*

- Galaxy configuration filenames

*In the config/ directory all the files are suffixed with ".sample".
Sample files are default files provided by Galaxy.
When you change a configuration it is strongly encouraged to work directly on a filename without the 'sample' suffix.
So, go ahead save the 'galaxy.ini.sample' and create a 'galaxy.ini'*

## What's next ?
Next, you will configure your Galaxy server for production.
The first thing to do is to disable the developer settings.

# Disable developer settings.
As previously mentionned the basic installation is for development only.
A set of option are usefull for development but become irrelevant for production:

- debug:

	debug = False

Setting debug to False disable middleware that loads the entire response in memory for displaying debugging information in the page.
If left enabled, the proxy server may timeout waiting for a response or your Galaxy process may run out of memory if it's serving large files.

- use_interactive:

	use_interactive = False

Setting use_interactive to false disable displaying and live debugging of tracebacks via the web.
Leaving it enabled will expose your configuration (database password, id_secret, etc.).

- filter-with:

## Subpart 1

Short introduction about this subpart.

{% icon hands_on %} ***Hands on!***

1. First step
2. Second step
3. Third step

## Subpart 2

Short introduction about this subpart.

{% icon hands_on %} ***Hands on!***

1. First step
2. Second step
3. Third step

Some blabla

{% icon hands_on %} ***Hands on!***

1. First step
2. Second step
3. Third step

# Part 2

Short introduction about this subpart.

{% icon hands_on %} ***Hands on!***

1. First step
2. Second step
3. Third step

## Subpart 2

Short introduction about this subpart.

{% icon hands_on %} ***Hands on!***

1. First step
2. Second step
3. Third step

Some blabla

{% icon hands_on %} ***Hands on!***

1. First step
2. Second step
3. Third step

# Conclusion

Conclusion about the technical key points. And then relation between the technics and the biological question to end with a global view.

:grey_exclamation: ***Key Points***

- *Simple sentence to sum up the first key point of the tutorial (Take home message)*
- *Second key point*
- *Third key point*
- *...*

# :clap: Thank you
