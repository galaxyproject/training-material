---
layout: tutorial_hands_on

title: "Galaxy Monitoring with Telegraf and Grafana"
zenodo_link: ""
questions:
  - How to monitor Galaxy with Telegraf
  - How do I set up InfluxDB
  - How can I make graphs in Grafana?
  - How can I best alert on important metrics?
objectives:
  - Setup InfluxDB
  - Setup Telegraf
  - Setup Grafana
  - Create several charts
time_estimation: "2h"
tags:
  - monitoring
subtopic: features
key_points:
  - Galaxy supports pluggable monitoring extensions.
  - Use grafana or the reports webapp to monitor your service.
contributors:
  - erasche
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
---


# Overview
{:.no_toc}

Monitoring is an incredibly important part of server monitoring and maintenance. Being able to observe trends and identify hot spots by collecting metrics gives you a significan ability to respond to any issues that arise in production. Monitoring is quite easy to get started with, it can be as simple as writing a quick shell script in order to start collecting metrics.


> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

This tutorial explicitly assumes you are starting from a setup like what was created in the *Galaxy installation with Ansible* tutorial.

# InfluxDB

[InfluxDB](https://www.influxdata.com/) provides the data storage for monitoring. It is a <abbr title="Time Series Database">TSDB</abbr>, so it has been designed specifically for storing time-series data like monitoring and metrics.

> ### {% icon hands_on %} Hands-on: Setting up InfluxDB
>
> 1. Edit your `requirements.yml` and add the following:
>
>    ```yaml
>    - src: mtchavez.influxdb
>      version: v6.0.0
>    ```
>
> 2. `ansible-galaxy install -p roles -r requirements.yml`
>
> 3. Create a new playbook, `monitoring.yml` with the following:
>
>    ```yaml
>    ---
>    - hosts: monitoring
>      become: true
>      roles:
>        - mtchavez.influxdb
>    ```
>
>    During this tutorial we will install everything on the same host, but often one keeps the monitoring infrastructure (Grafana, InfluxDB) on a separate host.
>
> 4. Edit the inventory file (`hosts`) an add a group for monitoring like:
>
>    ```ini
>    [monitoring]
>    training-0.example.org ansible_connection=local
>    ```
>
>    Ensure that the hostname is the full hostname of your machine.
>
> 4. Run the playbook:
>
>    ```
>    ansible-playbook -i hosts monitoring.yml
>    ```
>
{: .hands_on}

This will setup an InfluxDB server listening on port `:8086`. The service is currently unauthenticated but it is only listening on `localhost` so it is less of a concern. The service can be authenticated and SSL configured quite easily but that is outside the scope of this tutorial.

You can access the Influx service by running the command `influx`.

```
$ influx
Connected to http://localhost:8086 version 1.7.7
InfluxDB shell version: 1.7.7
> show databases
name: databases
name
----
_internal
> use _internal
Using database _internal
> show measurements
name: measurements
name
----
cq
database
httpd
queryExecutor
runtime
shard
subscriber
tsm1_cache
tsm1_engine
tsm1_filestore
tsm1_wal
write
```

This provides commands like `show databases` and others, but we will not use this interface very often.

# Grafana

[Grafana](https://grafana.com/) provides a visual interface to our metrics. It provides a nice query builder that provides a uniform experience across multiple backend databases, and many attractive graphing options. Another benefit is that many of the UseGalaxy.\* servers share their publish their dashboards publicly, and you can easily copy these and use them on your own server.

> ### {% icon hands_on %} Hands-on: Setting up Grafana
>
> 1. Edit your `requirements.yml` and add the following:
>
>    ```yaml
>    - src: cloudalchemy.grafana
>      version: 0.14.2
>    ```
>
> 2. `ansible-galaxy install -p roles -r requirements.yml`
>
> 3. Add `cloudalchemy.grafana` to your `monitoring.yml` playbook
>
> 4. Edit the file `group_vars/galaxyservers.yml` and set the following variables:
>
>    ```yaml
>    ---
>    grafana_url: "https://{{ inventory_hostname }}/grafana/"
>
>    grafana_security:
>        # Feel free to choose any other values here too
>        admin_user: admin
>        admin_password: password
>
>    # These datasources will be automatically included into Grafana
>    grafana_datasources:
>     - name: Galaxy
>       type: influxdb
>       access: proxy
>       url: http://127.0.0.1:8086
>       isDefault: true
>       version: 1
>       editable: false
>       database: telegraf
>    ```
>
> 5. Run the playbook:
>
>    ```
>    ansible-playbook -i hosts monitoring.yml
>    ```
>
> 5. Update the nginx configuration in `templates/nginx/galaxy.j2` to include the following at the end, before the last curly brace
>
>    ```nginx
>        ...
>        location /grafana/ {
>            proxy_pass http://127.0.0.1:3000/;
>        }
>    ```
>
> 5. Run the Galaxy playbook which includes nginx:
>
>    ```
>    ansible-playbook -i hosts galaxy.yml
>    ```
>
{: .hands_on}

This has now deployed Grafana on your domain under `/grafana/`, with the username and password you set. The datasource, from which Grafana obtains data, is preconfigured in Grafana.

# Telegraf

We use [Telegraf](https://github.com/influxdata/telegraf) for monitoring as it is incredibly easy to get started with, and it natively integrates with InfluxDB. They have extension documentation on how to configure different types of monitoring, and Telegraf supports [a huge array of inputs](https://github.com/influxdata/telegraf#input-plugins).

> ### {% icon hands_on %} Hands-on: Dependencies
>
> 1. Edit your `requirements.yml` and add the following:
>
>    ```yaml
>    - src: dj-wasabi.telegraf
>      version: 0.12.0
>    ```
>
> 2. Install the requirements
>
>    ```
>    ansible-galaxy install -p roles -r requirements.yml
>    ```
>
> 3. Add an entry to the end of your `galaxy.yml` playbook under `roles:`
>
>    ```yaml
>    - dj-wasabi.telegraf
>    ```
>
> 4. Open your group variables file, and add the following variables:
>
>    ```yaml
>    telegraf_agent_output:
>      - type: influxdb
>        config:
>        - urls = ["http://127.0.0.1:8086"]
>        - database = "telegraf"
>
>    telegraf_plugins_default:
>      - plugin: cpu
>      - plugin: disk
>      - plugin: kernel
>      - plugin: processes
>      - plugin: io
>      - plugin: mem
>      - plugin: system
>      - plugin: swap
>      - plugin: net
>      - plugin: netstat
>
>    telegraf_plugins_extra:
>      listen_galaxy_routes:
>        plugin: "statsd"
>        config:
>          - service_address = ":8125"
>          - metric_separator = "."
>          - allowed_pending_messages = 10000
>    ```
>
>    This configures telegraf to output to the configured influxdb server in the `telegraf` database. A number of plugins are enabled as `defaults`, this is useful if you have telegraf configured for multiple machines; you can have a base configuration that applies to all machines (perhaps in `group_vars/all.yml`), and then `extra` configuration that is per-machine.
>
>    We have configured the `statsd` plugin for telegraf, as we will use it to receive Galaxy timing data.
>
> 5. We need to enable Galaxy to send data to Telegraf:
>
>    In your group variables file, edit the Galaxy configuration in` galaxy_config > galaxy` and add `statsd_host: localhost` and `statsd_influxdb: true`. It should look like:
>
>    ```yaml
>    galaxy_config:
>      galaxy:
>        ...
>        statsd_host: localhost
>        statsd_influxdb: true
>    ```
>
> 6. Run the `galaxy.yml` playbook
>
{: .hands_on}

# Monitoring with Grafana

The stats have been collecting in InfluxDB for a few minutes, so now we will now configure Grafana with dashboards and alerting rules.

## Importing a dashboard

For any public Grafana dashboard, you can copy the dashboard for your own use. This is a nice features of Grafana that has really helped it spread in the Galaxy community, any cool thing one of builds, everyone else can copy and build upon.

> ### {% icon hands_on %} Hands-on: Import a dashboard
>
> 1. [Visit UseGalaxy.eu's Node Detail dashboard](https://stats.galaxyproject.eu/d/000000023/node-detail-infrastructure?orgId=1)
>
> 2. Look for the sharing icon at the top and click it
>
> 3. Under the "Export" tab, click "Save to file"
>
> 4. On your own Grafana server, on the home page, hover over the `+` icon and use "Import" from the menu.
>
> 5. Click "Upload .json file" and select the json dashboard you downloaded
>
> 6. Click "Import".
>
{: .hands_on}

With this, your first dashboard should be live! You should see some data from your Galaxy instance, like CPU/load/memory/etc. This can give you a nice `htop` like view into your systems, all collected in one easy dashboard. At the top you will see a box labelled "Host" with a dropdown. If you have more systems, you can click here to select between different machines.

## Setting up a Galaxy dashboard

Importing dashboards is a good start, but it's more interesting to create our own that's personalise to our needs.

> ### {% icon hands_on %} Hands-on: Create a dashboard
>
> 1. Again find the `+` icon in Grafana and create a dashboard. This will bring you to a new screen
>
> 2. Click *Add Query*, and you will be dropped into the Grafana query builder
>
>    ![Grafana query builder interface](../../images/grafana-query-builder.png)
>
>    This is the query builder interface. The interface somewhat resembles a SQL query, selecting data *from* a database, *where* it meets some condition, *select*ing some specific data, and *grouping by* time period. If this isn't immediately clear how it behaves, hopefully it will become more clear once you have built some queries.
>
> 3. Let's build a query:
>    - From:
>      - *"select measurement"*: `galaxy.`
>    - Select:
>      - *"field(value)"*: `field(mean)`
>    - Group by:
>      - *"fill(null)": `fill(none)`
>      - add new (+): `tag(path)`
>    - Alias by: `[[tag_path]]`
>
> 4. At the top of the page it probably says "Last 6 hours", click this to change to "Last 30 minutes"
>
> 5. Remember to save the dashboard, and give it a name like "Galaxy"
>
{: .hands_on}

## Styling

## Monitoring

