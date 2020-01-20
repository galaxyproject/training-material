---
layout: tutorial_hands_on

title: "Galaxy Monitoring with Reports"
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
  - natefoo
  - bgruening
  - slugger70
  - hexylena
subtopic: features
tags:
  - monitoring
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
---


# Overview
{:.no_toc}

The reports application gives some pre-configured analytics screens.

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Setting up Reports

The reports application is included with the Galaxy codebase and this tutorial assumes you've already done all of the setup required for Galaxy, Supervisord, uWSGI, and NGINX.

> ### {% icon hands_on %} Hands-on: Setup Reports
>
>
> 1. First we add a basic configuration of the Reports app to the playbook templates. Create `templates/galaxy/config/` folder, if it doesn't exist, and create `templates/galaxy/config/reports.yml` with the following contents:
>
>    {% raw %}
>    ```yml
>    uwsgi:
>        socket: 127.0.0.1:9001
>        buffer-size: 16384
>        processes: 1
>        threads: 4
>        offload-threads: 2
>        static-map: /static/style={{ galaxy_server_dir }}/static/style/blue
>        static-map: /static={{ galaxy_server_dir }}/static
>        static-map: /favicon.ico=static/favicon.ico
>        master: true
>        virtualenv: {{ galaxy_venv_dir }}
>        pythonpath: {{ galaxy_server_dir }}/lib
>        mount: /reports=galaxy.webapps.reports.buildapp:uwsgi_app()
>        manage-script-name: true
>        thunder-lock: false
>        die-on-term: true
>        hook-master-start: unix_signal:2 gracefully_kill_them_all
>        hook-master-start: unix_signal:15 gracefully_kill_them_all
>        py-call-osafterfork: true
>        enable-threads: true
>    reports:
>        cookie-path: /reports
>        database_connection: "postgresql:///galaxy?host=/var/run/postgresql"
>        file_path: /data
>        filter-with: proxy-prefix
>        template_cache_path: "{{ galaxy_mutable_data_dir }}/compiled_templates"
>    ```
>    {% endraw %}
>
> 2. In your `galaxyservers` group variables file, tell the playbook to deploy the reports configuration file:
>
>    {% raw %}
>    ```yml
>    galaxy_config_templates:
>    ...
>    - src: templates/galaxy/config/reports.yml
>      dest: "{{ galaxy_config_dir }}/reports.yml"
>    ```
>    {% endraw %}
>
>
> 3. Similar to Galaxy we will again use Supervisor to manage the Reports process. In the same file (`galaxyservers` group variables file) edit the `supervisor_programs` section - add an entry for the Reports webapp:
>
>    {% raw %}
>    ```yml
>    supervisor_programs:
>      ....
>      - name: reports
>        state: present
>        command: uwsgi --yaml {{ galaxy_config_dir }}/reports.yml
>        configuration: |
>          autostart=true
>          autorestart=true
>          startretries=1
>          startsecs=10
>          user=galaxy
>          umask=022
>          directory={{ galaxy_server_dir }}
>          environment=HOME={{ galaxy_mutable_data_dir }},VIRTUAL_ENV={{ galaxy_venv_dir }},PATH={{ galaxy_venv_dir }}/bin:%(ENV_PATH)s
>    ```
>    {% endraw %}
>
> 4. Then we need to tell NGINX it should serve our Reports app under `<server_url>/reports` url. Edit your `galaxyservers` group variables file, and under the NGINX configuraiton, add a block for proxying the reports application. It should look like:
>
>    ```nginx
>    location /reports {
>        uwsgi_pass           127.0.0.1:9001;
>        uwsgi_param          UWSGI_SCHEME $scheme;
>        include              uwsgi_params;
>    }
>    ```
>
> 5. Run the playbook
>
> 6. The reports application should be available, under `<server_url>/reports/`. Note the trailing slash.
>
{: .hands_on}
