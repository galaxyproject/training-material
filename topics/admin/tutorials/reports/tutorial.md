---
layout: tutorial_hands_on

title: "Galaxy Monitoring with Reports"
zenodo_link: ""
questions:
  - How to monitor a Galaxy service with the Reports application?
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
subtopic: monitoring
tags:
  - ansible
  - monitoring
  - git-gat
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
---

The reports application gives some pre-configured analytics screens. These are very easy to setup and can help with debugging issues in Galaxy.

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="reports" %}

# Setting up Reports

The reports application is included with the Galaxy codebase and this tutorial assumes you've already done all of the setup required for Galaxy, systemd, uWSGI, and NGINX.

> <hands-on-title>Setup Reports</hands-on-title>
>
>
> 1. First we add a basic configuration of the Reports app to the playbook templates. Create `templates/galaxy/config/` folder, if it doesn't exist, and create `templates/galaxy/config/reports.yml` with the following contents:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/templates/galaxy/config/reports.yml
>    @@ -0,0 +1,26 @@
>    +uwsgi:
>    +    socket: 127.0.0.1:9001
>    +    buffer-size: 16384
>    +    processes: 1
>    +    threads: 4
>    +    offload-threads: 2
>    +    static-map: /static/style={{ galaxy_server_dir }}/static/style/blue
>    +    static-map: /static={{ galaxy_server_dir }}/static
>    +    static-map: /favicon.ico=static/favicon.ico
>    +    master: true
>    +    virtualenv: {{ galaxy_venv_dir }}
>    +    pythonpath: {{ galaxy_server_dir }}/lib
>    +    mount: /reports=galaxy.webapps.reports.buildapp:uwsgi_app()
>    +    manage-script-name: true
>    +    thunder-lock: false
>    +    die-on-term: true
>    +    hook-master-start: unix_signal:2 gracefully_kill_them_all
>    +    hook-master-start: unix_signal:15 gracefully_kill_them_all
>    +    py-call-osafterfork: true
>    +    enable-threads: true
>    +reports:
>    +    cookie-path: /reports
>    +    database_connection: "{{ galaxy_config.galaxy.database_connection }}"
>    +    file_path: /data
>    +    filter-with: proxy-prefix
>    +    template_cache_path: "{{ galaxy_mutable_data_dir }}/compiled_templates"
>    {% endraw %}
>    ```
>    {: data-commit="Setup reports config file"}
>
>    {% snippet topics/admin/faqs/diffs.md %}
>
> 2. In your `galaxyservers` group variables file, tell the playbook to deploy the reports configuration file:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -51,6 +51,7 @@ galaxy_root: /srv/galaxy
>     galaxy_user: {name: galaxy, shell: /bin/bash}
>     galaxy_commit_id: release_22.05
>     galaxy_force_checkout: true
>    +galaxy_reports_path: "{{ galaxy_config_dir }}/reports.yml"
>     miniconda_prefix: "{{ galaxy_tool_dependency_dir }}/_conda"
>     miniconda_version: 4.7.12
>     miniconda_manage_dependencies: false
>    @@ -131,6 +132,8 @@ galaxy_config_templates:
>         dest: "{{ galaxy_config.galaxy.dependency_resolvers_config_file }}"
>       - src: templates/galaxy/config/tool_destinations.yml
>         dest: "{{ galaxy_config.galaxy.tool_destinations_config_file }}"
>    +  - src: templates/galaxy/config/reports.yml
>    +    dest: "{{ galaxy_reports_path }}"
>     
>     galaxy_local_tools:
>     - testing.xml
>    {% endraw %}
>    ```
>    {: data-commit="Deploy reports config to the config directory"}
>
>
> 3. Similar to Galaxy we will again use systemd to manage the Reports process.
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -144,6 +144,7 @@ galaxy_dynamic_job_rules:
>     
>     # systemd
>     galaxy_manage_systemd: true
>    +galaxy_manage_systemd_reports: yes
>     galaxy_systemd_env: [DRMAA_LIBRARY_PATH="/usr/lib/slurm-drmaa/lib/libdrmaa.so.1"]
>     
>     # Certbot
>    {% endraw %}
>    ```
>    {: data-commit="Enable the reports systemd unit"}
>
> 4. Then we need to tell NGINX it should serve our Reports app under `<server_url>/reports` url. Edit your `templates/nginx/galaxy.j2` file, and within the server block, add a block for proxying the reports application. It should look like:
>
>    {% raw %}
>    ```diff
>    --- a/templates/nginx/galaxy.j2
>    +++ b/templates/nginx/galaxy.j2
>    @@ -91,4 +91,10 @@ server {
>         }
>     
>         {{ tiaas_nginx_routes }}
>    +
>    +    location /reports/ {
>    +        uwsgi_pass           127.0.0.1:9001;
>    +        uwsgi_param          UWSGI_SCHEME $scheme;
>    +        include              uwsgi_params;
>    +    }
>     }
>    {% endraw %}
>    ```
>    {: data-commit="Configure a location block for Reports in NGINX"}
>
> 5. Run the playbook:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 6. The reports application should be available, under `<server_url>/reports/`.>
{: .hands_on}

> ```bash
> 1.sh
> ```
> {: data-test="true"}
{: .hidden}

> <comment-title>Insecure!</comment-title>
> But notice that your Reports server is not secured! Check out the [External Authentication]({% link topics/admin/tutorials/external-auth/tutorial.md %}) tutorial for information on securing Reports.
{: .comment}

{% snippet topics/admin/faqs/missed-something.md step=13 %}
