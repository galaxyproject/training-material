---
layout: tutorial_hands_on

title: "Galaxy Monitoring with gxadmin"
zenodo_link: ""
questions:
  - What is gxadmin
  - What can it do?
  - How to write a query?
objectives:
  - Learn gxadmin basics
  - See some queries and learn how they help debug production issues
time_estimation: "30m"
key_points:
  - gxadmin is a tool to run common database queries useful for Galaxy admins
  - new queries are welcome and easy to contribute
contributors:
  - hexylena
subtopic: monitoring
tags:
  - monitoring
  - ansible
  - git-gat
---

We will just briefly cover the features available in `gxadmin`, there are lots of queries that may or may not be useful for your Galaxy instance and you will have to read the documentation before using them.

It started life as a small shell script that Helena wrote because she couldn't remember what [Gravity](https://github.com/galaxyproject/gravity) was called or where it could be found. Some of the functions needed for things like swapping zerglings are still included in gxadmin but are highly specific to UseGalaxy.eu and not generally useful.

Since then it became the home for "all of the SQL queries we [galaxy admins] run regularly." @hexylena and @natefoo often shared SQL queries with each other in private chats, but this wasn't helpful to the admin community at large, so they decided to put them all in `gxadmin` and make it as easy to install as possible. We are continually trying to make this tool more generic and generally useful, if you notice something that's missing or broken, or have a new query you want to run, just [let us know](https://github.com/usegalaxy-eu/gxadmin/issues/new).

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="gxadmin" %}

## Installing gxadmin

It's simple to install gxadmin. Here's how you do it, if you haven't done it already.

> <hands-on-title>Installing gxadmin with Ansible</hands-on-title>
>
> 1. Edit your `requirements.yml` and add the following:
>
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -28,3 +28,5 @@
>       version: 0.1.0
>     - src: galaxyproject.pulsar
>       version: 1.0.8
>    +- src: galaxyproject.gxadmin
>    +  version: 0.0.8
>    {% endraw %}
>    ```
>    {: data-commit="Add requirement"}
>
>    {% snippet topics/admin/faqs/diffs.md %}
>
> 2. Install the role with:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-galaxy install -p roles -r requirements.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 3. Add the role to your playbook:
>
>    {% raw %}
>    ```diff
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -33,3 +33,4 @@
>         - galaxyproject.nginx
>         - galaxyproject.tusd
>         - galaxyproject.cvmfs
>    +    - galaxyproject.gxadmin
>    {% endraw %}
>    ```
>    {: data-commit="Add the gxadmin role"}
>
> 4. Run the playbook
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
{: .hands_on}

With that, `gxadmin` should be installed! Now, test it out:

> <hands-on-title>Test out gxadmin</hands-on-title>
>
> 1. Run `gxadmin` as the galaxy user and list recently registered users:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > sudo -u galaxy gxadmin query latest-users
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
>    > <code-in-title>Output</code-in-title>
>    > ```bash
>    >  id |          create_time          | disk_usage | username |       email        | groups | active
>    > ----+-------------------------------+------------+----------+--------------------+--------+--------
>    >   1 | 2021-06-09 12:25:59.299651+00 | 218 kB     | admin    | admin@example.org  |        | f
>    > (1 rows)
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
{: .hands_on}

> ```bash
> 1.sh
> ```
> {: data-test="true"}
{: .hidden}

## Configuration

If `psql` runs without any additional arguments, and permits you to access your galaxy database then you do not need to do any more configuration for gxadmin.
Otherwise, you may need to set some of the [PostgreSQL environment variables](https://github.com/usegalaxy-eu/gxadmin#postgres)

## Overview

`gxadmin` has several categories of commands, each with different focuses. This is not a technically meaningful separation, it is just done to make the interface easier for end users.

Category      | Keyword             | Purpose
---           | --                  | --
Configuration | `config`            | Commands relating to galaxy's configuration files like XML validation.
Filters       | `filter`            | Transforming streams of text.
Galaxy Admin  | `galaxy`            | Miscellaneous galaxy related commands like a cleanup wrapper.
uWSGI         | `uwsgi`             | If you're using [systemd for Galaxy](https://github.com/usegalaxy-eu/ansible-galaxy-systemd/) and a handler/zergling setup, then this lets you manage your handlers and zerglings.
DB Queries    | `{csv,tsv,i,}query` | Queries against the database which return tabular output.
Report        | `report`            | Queries which return more complex and structured markdown reports.
Mutations     | `mutate`            | These are like queries, except they mutate the database. All other queries are read-only.
Meta          | `meta`              | More miscellaneous commands, and a built-in updating function.

## Admin Favourite Queries

**@slugger70's favourite**: `gxadmin query old-histories`. He contributed this function to find old histories, as their instance has a 90 day limit on histories, anything older than that might be automatically removed. This helps their group identify any histories that can be purged in order to save space. Running this on UseGalaxy.eu, we have some truly ancient histories, and maybe could benefit from a similar policy.

id  |        update-time         | user-id | email |           name           | published | deleted | purged | hid-counter
--- | -------------------------- | ------- | ----- | ------------------------ | --------- | ------- | ------ | ------------
361 | 2013-02-24 16:27:29.197572 |     xxx | xxxx  | Unnamed history          | f         | f       | f      |           6
362 | 2013-02-24 15:31:05.804747 |     xxx | xxxx  | Unnamed history          | f         | f       | f      |           1
347 | 2013-02-22 15:59:12.044169 |     xxx | xxxx  | Unnamed history          | f         | f       | f      |          19
324 | 2013-02-22 15:57:54.500637 |     xxx | xxxx  | Exercise 5               | f         | f       | f      |          64
315 | 2013-02-22 15:50:51.398894 |     xxx | xxxx  | day5 practical           | f         | f       | f      |          90
314 | 2013-02-22 15:45:47.75967  |     xxx | xxxx  | 5. Tag Galaxy-Kurs       | f         | f       | f      |          78

**@natefoo's favourite**: `gxadmin query job-inputs`. He contributed this function which helps him debug jobs which are not running and should be.

hda-id   | hda-state | hda-deleted | hda-purged |  d-id   | d-state | d-deleted | d-purged | object-store-id
-------- | --------- | ----------- | ---------- | ------- | ------- | --------- | -------- | ----------------
8638197  |           | f           | f          | 8246854 | running | f         | f        | files9
8638195  |           | f           | f          | 8246852 | running | f         | f        | files9
8638195  |           | f           | f          | 8246852 | running | f         | f        | files9

**@bgruening's favourite**: `gxadmin query latest-users` let's us see who has recently joined our server. We sometimes notice that people are running a training on our infrastructure and they haven't registered for [training infrastructure as a service](https://galaxyproject.eu/tiaas) which helps us coordinate infrastructure for them so they don't have bad experiences.

id    |        create_time         | disk_usage | username | email | groups | active
----- | -------------------------- | ---------- | -------- | ----- | ------ | -------
3937  | 2019-01-27 14:11:12.636399 | 291 MB     | xxxx     | xxxx  |        | t
3936  | 2019-01-27 10:41:07.76126  | 1416 MB    | xxxx     | xxxx  |        | t
3935  | 2019-01-27 10:13:01.499094 | 2072 kB    | xxxx     | xxxx  |        | t
3934  | 2019-01-27 10:06:40.973938 | 0 bytes    | xxxx     | xxxx  |        | f
3933  | 2019-01-27 10:01:22.562782 |            | xxxx     | xxxx  |        | f

**@hexylena's favourite** `gxadmin report job-info`. This command gives more information than you probably need on the execution of a specific job, formatted as markdown for easy sharing with fellow administrators.

```text
# Galaxy Job 5132146

Property      | Value
------------- | -----
         Tool | toolshed.g2.bx.psu.edu/repos/bgruening/canu/canu/1.7
        State | running
      Handler | handler_main_2
      Created | 2019-04-20 11:04:40.854975+02 (3 days 05:49:30.451719 ago)
Job Runner/ID | condor / 568537
        Owner | e08d6c893f5

## Destination Parameters

Key                     |   Value
---                     |   ---
description             |   `canu`
priority                |   `-128`
request_cpus            |   `20`
request_memory          |   `64G`
requirements            |   `GalaxyGroup == "compute"`
tmp_dir                 |   `True`

## Dependencies

Name   |   Version   |   Dependency Type   |   Cacheable   |   Exact   |   Environment Path                          |   Model Class
---    |   ---       |   ---               |   ---         |   ---     |   ---                                       |   ---
canu   |   1.7       |   conda             |   false       |   true    |   /usr/local/tools/_conda/envs/__canu@1.7   |   MergedCondaDependency

## Tool Parameters

Name                 |   Settings
---------            |   ------------------------------------
minOverlapLength     |   500
chromInfo            |   /opt/galaxy/tool-data/shared/ucsc/chrom/?.len
stage                |   all
contigFilter         |   {lowCovDepth: 5, lowCovSpan: 0.5, minLength: 0, minReads: 2, singleReadSpan: 1.0}
s                    |   null
mode                 |   -nanopore-raw
dbkey                |   ?
genomeSize           |   300000
corOutCoverage       |   40
rawErrorRate         |
minReadLength        |   1000
correctedErrorRate   |

## Inputs

Job ID    |   Name                |   Extension     |   hda-id    |   hda-state   |   hda-deleted   |   hda-purged   |   ds-id     |   ds-state   |   ds-deleted   |   ds-purged   |   Size
----      |   ----                |   ----          |   ----      |   ----        |   ----          |   ----         |   ----      |   ----       |   ----         |   ----        |   ----
4975404   |   Osur_record.fastq   |   fastqsanger   |   9517188   |               |   t             |   f            |   9015329   |   ok         |   f            |   f           |   3272 MB
4975404   |   Osur_record.fastq   |   fastqsanger   |   9517188   |               |   t             |   f            |   9015329   |   ok         |   f            |   f           |   3272 MB

## Outputs

Name                                          |   Extension   |   hda-id    |   hda-state   |   hda-deleted   |   hda-purged   |   ds-id     |   ds-state   |   ds-deleted   |   ds-purged   |   Size
----                                          |   ----        |   ----      |   ----        |   ----          |   ----         |   ----      |   ----       |   ----         |   ----        |   ----
Canu assembler on data 41 (trimmed reads)     |   fasta.gz    |   9520369   |               |   f             |   f            |   9018510   |   running    |   f            |   f           |
Canu assembler on data 41 (corrected reads)   |   fasta.gz    |   9520368   |               |   f             |   f            |   9018509   |   running    |   f            |   f           |
Canu assembler on data 41 (unitigs)           |   fasta       |   9520367   |               |   f             |   f            |   9018508   |   running    |   f            |   f           |
Canu assembler on data 41 (unassembled)       |   fasta       |   9520366   |               |   f             |   f            |   9018507   |   running    |   f            |   f           |
Canu assembler on data 41 (contigs)           |   fasta       |   9520365   |               |   f             |   f            |   9018506   |   running    |   f            |   f           |
```

# `gxadmin` for Monitoring

`gxadmin` already supported `query`, `csvquery`, and `tsvquery` for requesting data from the Galaxy database in tables, CSV, or TSV formats, but we recently implemented `influx` queries which output data in a format that [Telegraf](https://github.com/influxdata/telegraf) can consume.

So running `gxadmin query queue-overview` normally shows something like:

tool_id                                                                                |  tool_version  |    destination_id    |     handler     |  state  | job_runner_name | count
-------------------------------------------------------------------------------------- | -------------- | -------------------- | --------------- | ------- | --------------- | ------
toolshed.g2.bx.psu.edu/repos/iuc/unicycler/unicycler/0.4.6.0                           | 0.4.6.0        | 12cores_180G_special | handler_main_4  | running | condor          |     1
toolshed.g2.bx.psu.edu/repos/iuc/unicycler/unicycler/0.4.6.0                           | 0.4.6.0        | 12cores_180G_special | handler_main_5  | running | condor          |     1
toolshed.g2.bx.psu.edu/repos/devteam/freebayes/freebayes/1.1.0.46-0                    | 1.1.0.46-0     | 12cores_12G          | handler_main_3  | running | condor          |     2
toolshed.g2.bx.psu.edu/repos/iuc/qiime_extract_barcodes/qiime_extract_barcodes/1.9.1.0 | 1.9.1.0        | 4G_memory            | handler_main_1  | running | condor          |     1
toolshed.g2.bx.psu.edu/repos/iuc/hisat2/hisat2/2.1.0+galaxy3                           | 2.1.0+galaxy3  | 8cores_20G           | handler_main_11 | running | condor          |     1
toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.72                                | 0.72           | 20G_memory           | handler_main_11 | running | condor          |     4
ebi_sra_main                                                                           | 1.0.1          | 4G_memory            | handler_main_3  | running | condor          |     2
ebi_sra_main                                                                           | 1.0.1          | 4G_memory            | handler_main_4  | running | condor          |     3


The `gxadmin iquery queue-overview` is run by our Telegraf monitor on a regular basis, allowing us to consume the data:

```
queue-overview,tool_id=toolshed.g2.bx.psu.edu/repos/iuc/unicycler/unicycler/0.4.6.0,tool_version=0.4.6.0,state=running,handler=handler_main_4,destination_id=12cores_180G_special,job_runner_name=condor count=1
queue-overview,tool_id=toolshed.g2.bx.psu.edu/repos/iuc/unicycler/unicycler/0.4.6.0,tool_version=0.4.6.0,state=running,handler=handler_main_5,destination_id=12cores_180G_special,job_runner_name=condor count=1
queue-overview,tool_id=toolshed.g2.bx.psu.edu/repos/devteam/freebayes/freebayes/1.1.0.46-0,tool_version=1.1.0.46-0,state=running,handler=handler_main_3,destination_id=12cores_12G,job_runner_name=condor count=1
queue-overview,tool_id=toolshed.g2.bx.psu.edu/repos/devteam/vcffilter/vcffilter2/1.0.0_rc1+galaxy1,tool_version=1.0.0_rc1+galaxy1,state=queued,handler=handler_main_11,destination_id=4G_memory,job_runner_name=condor count=1
queue-overview,tool_id=toolshed.g2.bx.psu.edu/repos/iuc/hisat2/hisat2/2.1.0+galaxy3,tool_version=2.1.0+galaxy3,state=running,handler=handler_main_11,destination_id=8cores_20G,job_runner_name=condor count=1
queue-overview,tool_id=toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.72,tool_version=0.72,state=running,handler=handler_main_11,destination_id=20G_memory,job_runner_name=condor count=4
queue-overview,tool_id=ebi_sra_main,tool_version=1.0.1,state=running,handler=handler_main_3,destination_id=4G_memory,job_runner_name=condor count=2
queue-overview,tool_id=ebi_sra_main,tool_version=1.0.1,state=running,handler=handler_main_4,destination_id=4G_memory,job_runner_name=condor count=3
```

And produce [some nice graphs](https://stats.galaxyproject.eu/) from it.

You can use an influx configuration like:

```toml
[[inputs.exec]]
    commands = ["/usr/bin/galaxy-queue-size"]
    timeout = "10s"
    data_format = "influx"
    interval = "1m"
```

This often requires a wrapper script, because you'll need to pass environment variables to the `gxadmin` invocation, e.g.:

```bash
#!/bin/bash
export PGUSER=galaxy
export PGHOST=dbhost
gxadmin iquery queue-overview --short-tool-id
gxadmin iquery workflow-invocation-status
```


> <tip-title>Which queries support iquery?</tip-title>
> This data is not currently exposed, so, just try the queries. But it's easy to add influx support when missing! [Here is an example](https://github.com/usegalaxy-eu/gxadmin/blob/e0ec0174ebbdce1acd8c40c7431308934981aa0c/parts/22-query.sh#L54), we set the variables in a function:
>
> ```
> fields="count=1"
> tags="tool_id=0"
> ```
>
> This means: column 0 is a tag named tool_id, and column 1 is a field (real value) named count.
> [Here is an example](https://github.com/usegalaxy-eu/gxadmin/blob/e0ec0174ebbdce1acd8c40c7431308934981aa0c/parts/22-query.sh#L1987) that has multiple fields that are stored.
>
{: .tip}

# Implementing a Query

Queries are really easy to implement! All you have to do is add your SQL, with a small bash function to wrap it. `gxadmin` supports 'local' functions, which you can add locally without contributing back. We strongly encourage you to contribute your functions back to `gxadmin` though, you'll never know who wants to know the same thing about their db.

`gxadmin` will look for local functions in `~/.config/gxadmin-local.sh`

## A basic function

> <hands-on-title>Implementing a basic function</hands-on-title>
>
> 1. If `~/.config/` does not exist, create that directory with `mkdir -p ~/.config/`
>
> 2. Open `~/.config/gxadmin-local.sh` in a text editor.
>
> 3. Add the following to the file and save it.
>
>    ```bash
>    local_hello() { ## : Says hi
>    	echo "hi!"
>    }
>    ```
>
> 4. Run `gxadmin local`
>
>    > <question-title></question-title>
>    >
>    > What do you see?
>    >
>    > > <solution-title></solution-title>
>    > > It should look like:
>    > >
>    > > ```console
>    > > gxadmin usage:
>    > >
>    > > Local: (These can be configured in /home/hxr/.config/gxadmin-local.sh)
>    > >
>    > >     local hello          Says hi
>    > >
>    > > help / -h / --help : this message. Invoke '--help' on any subcommand for help specific to that subcommand
>    > > ```
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
> 5. Run `gxadmin local hello`
>
>    > <question-title></question-title>
>    >
>    > What do you see?
>    >
>    > > <solution-title></solution-title>
>    > > It should output 'hi!'
>    > {: .solution}
>    >
>    {: .question}
{: .hands_on}

This is the simplest possible function you can add, and is pretty limited in its functionality. This can offer you a nice place to put all of your existing bash scripts, and have autogenerated help for them.

## Adding help

Every function is improved by documentation! Let's add that now:

> <hands-on-title>Adding help</hands-on-title>
>
> 1. Open `~/.config/gxadmin-local.sh` in a text editor.
>
> 2. Update your function to add the `handle_help` call:
>
>    ```bash
>    local_hello() { ## : Says hi
>    	handle_help "$@" <<-EOF
>    		Greets you
>    	EOF
>    	echo "hi!"
>    }
>    ```
>
> 4. Run `gxadmin local hello --help`
>
>    > <question-title></question-title>
>    >
>    > What do you see?
>    >
>    > > <solution-title></solution-title>
>    > > It should look like:
>    > >
>    > > ```console
>    > > $ gxadmin local hello --help
>    > > **NAME**
>    > >
>    > > local hello - Says hi
>    > >
>    > > **SYNOPSIS**
>    > >
>    > > gxadmin local hello
>    > >
>    > > **NOTES**
>    > >
>    > > Greets you
>    > > ```
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
{: .hands_on}


## Adding a query

The bulk of gxadmin is not functions calling shell commands though, it's mostly SQL queries. So let's find the N most recent jobs

> <hands-on-title>Adding a query</hands-on-title>
>
> 1. Open `~/.config/gxadmin-local.sh` in a text editor.
>
> 2. Add a new function:
>
>    ```bash
>    local_query-latest() { ## [jobs|10]: Queries latest N jobs (default to 10)
>    	handle_help "$@" <<-EOF
>    		Find information about the latest jobs on your server.
>    	EOF
>
>    	# Value of first argument, or 10 if isn't supplied
>    	job_limit=${1:-10}
>
>    	# Here we store the query in a bash variable named QUERY
>    	read -r -d '' QUERY <<-EOF
>    		SELECT
>    			id, tool_id, tool_version, state
>    		FROM job
>    		ORDER BY id desc
>    		LIMIT ${job_limit}
>    	EOF
>    }
>    ```
>
> 4. Run `gxadmin local query-latest 5` to select
>
>    > <question-title></question-title>
>    >
>    > What do you see?
>    >
>    > > <solution-title></solution-title>
>    > > It should similar to the following, assuming you've run tools in your Galaxy
>    > >
>    > > ```console
>    > > $ gxadmin local query-latest 5
>    > >  id  |         tool_id          | tool_version | state
>    > > -----+--------------------------+--------------+-------
>    > >  103 | upload1                  | 1.1.6        | ok
>    > >  102 | upload1                  | 1.1.6        | ok
>    > >  101 | upload1                  | 1.1.6        | error
>    > >  100 | circos                   | 0.91         | ok
>    > >   99 | circos                   | 0.91         | ok
>    > > (5 rows)
>    > > ```
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
{: .hands_on}

# Summary

There are a lot of queries, all tailored to specific use cases, some of these may be interesting for you, some may not. These are [all documented](https://github.com/usegalaxy-eu/gxadmin#commands) with example inputs and outputs in the gxadmin readme, and help is likewise available from the command line.

{% snippet topics/admin/faqs/missed-something.md step=10 %}
