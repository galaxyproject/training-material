---
layout: tutorial_hands_on

title: "Galaxy Reports"
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
  - erasche
---


# Overview
{:.no_toc}



> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}


> ### {% icon hands_on %} Hands-on: Setup Reports
>
> 1. Edit your `galaxyservers` group variables file, and under the NGINX configuraiton, add a block for proxying the reports application. It should look like:
>
>    ```nginx
>    location /reports {
>        uwsgi_pass           127.0.0.1:9001;
>        uwsgi_param          UWSGI_SCHEME $scheme;
>        include              uwsgi_params;
>    }
>    ```
>
> 2. In your `supervisor_programs`, add an entry for the reports webapp:
>
>    {% raw %}
>    ```yaml
>    supervisor_programs:
>      ....
>      - name: reports
>        state: present
>        command: uwsgi --yaml {{ galaxy_config_dir }}/reports.yml
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
> 3. Create `files/galaxy/config/reports.yml` with the following contents:
>
>    ```yml
>    ```
>
{: .hands_on}




# Monitoring and maintenance

## Running the Reports Application

### Section 1 - Configure reports

Begin by making a copy of the reports config file to your config directory, and editing it:

```console
$ sudo -u galaxy cp /srv/galaxy/server/config/reports.ini.sample /srv/galaxy/config/reports.ini
$ sudo -u galaxy -e /srv/galaxy/config/reports.ini
```

Since we serve Galaxy at the root of our webserver, we'll need to serve Reports from a subdirectory: `/reports`. This is the default if we enable the `proxy-prefix` filter, all we need to do is uncomment the `proxy-prefix` setting. We also need to point the reports application at Galaxy's PostgreSQL database:

```ini
filter-with = proxy-prefix
cookie_path = /reports
database_connection = postgresql:///galaxy?host=/var/run/postgresql
file_path = /srv/galaxy/data
```

