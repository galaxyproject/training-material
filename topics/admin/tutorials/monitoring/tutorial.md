---
layout: tutorial_hands_on

title: "Galaxy Monitoring"
zenodo_link: ""
questions:
  - How to monitor a Galaxy service with the Reports?
objectives:
  - Setup and start the Galaxy reports app.
time_estimation: "30m"
key_points:
  - Galaxy supports pluggable monitoring extensions.
  - The Reports webapp is one option to monitor your system.
contributors:
  - erasche
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

# Telegraf

We use [Telegraf](https://github.com/influxdata/telegraf) for monitoring as it is incredibly easy to get started with.

> ### {% icon hands_on %} Hands-on: Dependencies
>
> 1. Edit your `requirements.yml` and add the following:
>
>    ```yaml
>    - dj-wasabi.telegraf
>    ```
>
> 2. `ansible-galaxy install -p roles -r requirements.yml`
>
{: .hands_on}

This role has many convenient variables for us to set. We'll do that now:

> ### {% icon hands_on %} Hands-on: Telegraf Configuration
>
> 1. Open your group variables file, and add the following variables:
>
>    ```yaml
>    telegraf_agent_version: 1.9.3
>
>    telegraf_agent_output:
>      - type: influxdb
>        config:
>        - urls = ["YOUR-INFLUX-URL"]
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
> 2. We need to enable Galaxy to send data to Telegraf:
>
>    In your group variables file, edit the Galaxy configuration in` galaxy_config > galaxy` and add `statsd_host: localhost` and `statsd_influxdb: true`. It should look like:
>
>    ```yaml
>    galaxy_config:
>      galaxy:
>        statsd_host: localhost
>        statsd_influxdb: true
>    ```
>
> 3. Run the playbook
>
{: .hands_on}

The stats should show up in Grafana now, allowing you to see node detail information and route timings. If you need dashboards, you can [copy the EU dashboards](http://influx.training.galaxyproject.eu:8080/d/000000023/node-detail?orgId=1)
