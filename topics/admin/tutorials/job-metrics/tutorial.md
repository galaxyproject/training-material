---
layout: tutorial_hands_on

title: "Recording Job Metrics"
zenodo_link: ""
questions:
objectives:
  - What are job metrics?
  - What sort of information can I collect?
  - Where can I find this information?
time_estimation: "15m"
key_points:
contributors:
  - hexylena
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
      - connect-to-compute-cluster
---

# Overview
{:.no_toc}

Job metrics record properties of the jobs that are executed, information that can help you plan for trainings or plan capacity for further expansions of your Galaxy server.

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Metrics

Galaxy includes a built-in framework to collect job metrics and store these in its database. [Some work was done](https://github.com/galaxyproject/galaxy/blob/dev/lib/galaxy/config/sample/job_metrics_conf.xml.sample) to try and analyse job runtime metrics to optimise cluster allocation based on job inputs, and enhance job submission ({% cite Tyryshkina_2019 %}). More work will be done in this area.

> ### {% icon comment %} Note
>
> Job metrics are only visible to Galaxy *admin users*, unless you set `expose_potentially_sensitive_job_metrics: true`, like UseGalaxy.eu does. EU's intention with this is to empower users and make everything as transparent as possible.
>
{: .comment}

## Setting up Galaxy

By default, Galaxy enables the `core` metrics:

![screenshot of galaxy metrics](../../images/job-metrics-basic.png)

These include very basic submission parameters. We want more information!

> ### {% icon hands_on %} Hands-on: Setting up the job metrics file
>
> 1. Create the file `files/galaxy/config/job_metrics_conf.xml` with the following contents:
>
>    ```xml
>    <?xml version="1.0"?>
>    <job_metrics>
>      <core />
>      <cpuinfo />
>      <meminfo />
>      <uname />
>      <env />
>      <cgroup />
>     <hostname />
>    </job_metrics>
>    ```
>
>    You can see the [sample file](https://github.com/galaxyproject/galaxy/blob/dev/lib/galaxy/config/sample/job_metrics_conf.xml.sample) for further options regarding metrics.
>
> 2. Edit your playbook to install the package named `cgroup-tools` in a pre-task (with git/make/etc). This package is required to use `cgget` which is used in metrics collection.
>
> 3. Edit the group variables file, `group_vars/galaxyservers.yml`:
>
>    You'll need to make two edits:
>    - Setting the `job_metrics_config_file`, to tell Galaxy where to look for the job metrics configuration.
>    - Adding the file to the list of `galaxy_config_files` to deploy it to the server:
>
>    {% raw %}
>    ```diff
>    --- galaxyservers.yml.old
>    +++ galaxyservers.yml
>
>    + galaxy_job_metrics_config_file: "{{ galaxy_config_dir }}/job_metrics_conf.xml"
>
>      galaxy_config:
>        galaxy:
>    +     job_metrics_config_file: "{{ galaxy_job_metrics_config_file }}"
>          brand: "My Galaxy"
>          admin_users: admin@example.org
>          database_connection: "postgresql:///galaxy?host=/var/run/postgresql"
>    @@ -120,6 +121,8 @@ gie_proxy_setup_service: systemd
>      gie_proxy_sessions_path: "{{ galaxy_mutable_data_dir }}/interactivetools_map.sqlite"
> 
>      galaxy_config_files:
>    +   - src: files/galaxy/config/job_metrics_conf.xml
>    +     dest: "{{ galaxy_job_metrics_config_file }}"
>        - src: files/galaxy/config/tool_conf_interactive.xml
>          dest: "{{ galaxy_config_dir }}/tool_conf_interactive.xml"
>        - src: files/galaxy/config/job_conf.xml
>    ```
>    {% endraw %}
>
> 4. Run the playbook
>
>    ```
>    ansible-playbook galaxy.yml
>    ```
{: .hands_on}


## Generating Metrics

With this, the job metrics tracking should be set up. Now when you run a job, you will see many more metrics:

> ### {% icon hands_on %} Hands-on: Generate some metrics
>
> 1. Run a job (any tool is fine, even upload)
>
> 2. View the information of the output dataset ({% icon galaxy-info %})
>
{: .hands_on}

![advanced metrics](../../images/job-metrics-advanced.png)


## What should I collect?

There is not a good rule we can tell you, just choose what you think is useful or will be. Numeric parameters are "cheaper" than the text parameters like uname to store, eventually you may find yourself wanting to remove old job metrics if you decide to collect the environment variables or similar.
