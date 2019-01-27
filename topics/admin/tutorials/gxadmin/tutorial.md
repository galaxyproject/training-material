---
layout: tutorial_hands_on

title: "gxadmin"
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
  - gxadmin formalizes + shares many queries that admins have been running for a long time in private
  - new queries are welcome and easy to contribute
contributors:
  - erasche
---

# Overview
{:.no_toc}

We will just briefly cover the features available in `gxadmin`, there are lots of queries that may or may not be useful for your Galaxy instance and you will have to read the documentation before using them.

It started life as a small shell script that Helena wrote because she couldn't remember what [Gravity](https://github.com/galaxyproject/gravity) was called or where it could be found. Some of the functions needed for things like swapping zerglings are still included in gxadmin but are highly specific to UseGalaxy.eu and not generally useful.

Since then it became the home for "all of the SQL queries we [galaxy admins] run regularly." @erasche and @natefoo often shared SQL queries with each other in private chats, but this wasn't helpful to the admin community at large, so they decided to put them all in `gxadmin` and make it as easy to install as possible. We are continually trying to make this tool more generic and generally useful, if you notice something that's missing or broken, or have a new query you want to run, just [let us know](https://github.com/usegalaxy-eu/gxadmin/issues/new).

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Configuration

If `psql` runs without any additional arguments, and permits you to access your galaxy database then you do not need to do any more configuration for gxadmin.
Otherwise, you may need to set some of the [PostgreSQL environment variables](https://github.com/usegalaxy-eu/gxadmin#query-setup)

# Queries

**@slugger70's favourite**: `gxadmin query old-histories`. He contributed this function to find old histories, as their instance has a 90 day limit on histories, anything older than that might be automatically removed. This helps their group identify any histories that can be purged in order to save space. Running this on UseGalaxy.eu, we have some truly ancient histories, and maybe could benefit from a similar policy.

id  |        update-time         | user-id | email |           name           | published | deleted | purged | hid-counter
--- | -------------------------- | ------- | ----- | ------------------------ | --------- | ------- | ------ | ------------
361 | 2013-02-24 16:27:29.197572 |     xxx | xxxx  | Unnamed history          | f         | f       | f      |           6
362 | 2013-02-24 15:31:05.804747 |     xxx | xxxx  | Unnamed history          | f         | f       | f      |           1
347 | 2013-02-22 15:59:12.044169 |     xxx | xxxx  | Unnamed history          | f         | f       | f      |          19
324 | 2013-02-22 15:57:54.500637 |     xxx | xxxx  | Exercise 5               | f         | f       | f      |          64
315 | 2013-02-22 15:50:51.398894 |     xxx | xxxx  | day5 practical           | f         | f       | f      |          90
314 | 2013-02-22 15:45:47.75967  |     xxx | xxxx  | 5. Tag Galaxy-Kurs       | f         | f       | f      |          78
330 | 2013-02-22 15:44:11.099138 |     xxx | xxxx  | Day 4 practise           | f         | f       | f      |          92
320 | 2013-02-22 15:29:07.426646 |     xxx | xxxx  | Day 5                    | f         | f       | f      |          55
327 | 2013-02-22 15:24:32.155946 |     xxx | xxxx  | Galaxy course day 5 - 22 | f         | f       | f      |          41
345 | 2013-02-22 15:22:04.899407 |     xxx | xxxx  | Day 5 - Question 2       | f         | f       | f      |           9
346 | 2013-02-22 15:16:38.934597 |     xxx | xxxx  | Galaxy Course Day 5 - Q5 | f         | f       | f      |          46
326 | 2013-02-22 14:43:10.233117 |     xxx | xxxx  | Day 5                    | f         | f       | f      |          22
343 | 2013-02-22 14:10:29.035706 |     xxx | xxxx  | Unnamed history          | f         | f       | f      |           1
246 | 2013-02-22 13:57:37.332747 |     xxx | xxxx  | Galaxy Course Day3       | f         | f       | f      |          35

**@natefoo's favourite**: `gxadmin query job-inputs`. He contributed this function which helps him debug jobs which are not running and should be. The query can

hda-id   | hda-state | hda-deleted | hda-purged |  d-id   | d-state | d-deleted | d-purged | object-store-id
-------- | --------- | ----------- | ---------- | ------- | ------- | --------- | -------- | ----------------
8638197  |           | f           | f          | 8246854 | running | f         | f        | files9
8638195  |           | f           | f          | 8246852 | running | f         | f        | files9
8638195  |           | f           | f          | 8246852 | running | f         | f        | files9

**@bgruening's favourite**: `gxamdin query latest-users` let's us see who has recently joined our server. We sometimes notice that people are running a training on our infrastructure and they haven't registered for [training infrastructure as a service](https://galaxyproject.eu/tiaas) which helps us coordinate infrastructure for them so they don't have bad experiences.

id    |        create_time         | disk_usage | username | email | groups | active
----- | -------------------------- | ---------- | -------- | ----- | ------ | -------
3937  | 2019-01-27 14:11:12.636399 | 291 MB     | xxxx     | xxxx  |        | t
3936  | 2019-01-27 10:41:07.76126  | 1416 MB    | xxxx     | xxxx  |        | t
3935  | 2019-01-27 10:13:01.499094 | 2072 kB    | xxxx     | xxxx  |        | t
3934  | 2019-01-27 10:06:40.973938 | 0 bytes    | xxxx     | xxxx  |        | f
3933  | 2019-01-27 10:01:22.562782 |            | xxxx     | xxxx  |        | f
3932  | 2019-01-27 09:37:34.11496  | 8299 kB    | xxxx     | xxxx  |        | t
3931  | 2019-01-27 05:53:33.072407 |            | xxxx     | xxxx  |        | t
3930  | 2019-01-27 05:05:07.323792 | 85 GB      | xxxx     | xxxx  |        | t
3929  | 2019-01-27 02:58:08.580324 | 795 MB     | xxxx     | xxxx  |        | t
3928  | 2019-01-26 21:07:20.483847 |            | xxxx     | xxxx  |        | t
3927  | 2019-01-26 20:24:20.016207 | 0 bytes    | xxxx     | xxxx  |        | t
3926  | 2019-01-26 07:51:06.98045  | 331 GB     | xxxx     | xxxx  |        | t
3925  | 2019-01-26 07:47:29.132961 | 0 bytes    | xxxx     | xxxx  |        | t
3924  | 2019-01-26 06:48:04.093405 | 0 bytes    | xxxx     | xxxx  |        | f
3923  | 2019-01-26 04:51:03.47129  | 15 GB      | xxxx     | xxxx  |        | t
3922  | 2019-01-25 20:04:40.934584 | 7886 MB    | xxxx     | xxxx  |        | t

**@erasche's favourite** `gxadmin iquery queue-overview`. `gxadmin` already supported `query`, `csvquery`, and `tsvquery` for requesting data from the Galaxy database in tables, CSV, or TSV formats, but we recently implemented `influx` queries which output data in a format that [Telegraf](https://github.com/influxdata/telegraf) can consume.

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

# Implementing a Query

We will not do a proper hands-on, but instead show [@slugger70's PR to find old histories](https://github.com/usegalaxy-eu/gxadmin/pull/5/files#diff-f905e55928aae903b7e13cc72842af3c), he implemented a function, provided some help output in a formatted manner, and then wrote his SQL query. If you don't feel comfortable writing bash, just share any SQL you've written and we will help you add it.

# Summary

There are a lot of queries, all tailored to specific use cases, some of these may be interesting for you, some may not. These are [all documented](https://github.com/usegalaxy-eu/gxadmin#commands) with example inputs and outputs in the gxadmin readme, and help is likewise available from the command line.
